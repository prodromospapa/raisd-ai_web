#!/usr/bin/env python3
import os
import uuid
import shlex
import shutil
import subprocess
from io import StringIO
from typing import Optional, Tuple
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless servers
import json
import time
import threading
from concurrent.futures import ThreadPoolExecutor
import matplotlib.pyplot as plt
import pysam

from collections import Counter
from math import lgamma
from types import SimpleNamespace


from flask import (
    Flask, render_template, request, redirect, url_for,
    send_from_directory, flash
)
from werkzeug.utils import secure_filename

# ------------------------------ Config ------------------------------
BASE_DIR      = os.path.abspath(os.path.dirname(__file__))
DATA_DIR      = os.path.join(BASE_DIR, "data")
RUNS_DIR      = os.path.join(BASE_DIR, "runs")
ALLOWED_EXT   = {"vcf", "gz", "bcf"}
MAX_CONTENT   = 2 * 1024 * 1024 * 1024  # 2 GB
PLOTS_FILE    = "sweep_TR_plot.png"
RAISD_REPORT  = "RAiSD_Report.results"
RAISD_GRIDDIR = "RAiSD_Grid.results"
# How long (seconds) to keep an unfinished run after the client last ping.
# If the client leaves the page (stops polling / heartbeat), the background
# cleaner will remove the run directory after this timeout. Default: 5 minutes.
RUN_DIR_TIMEOUT_SECONDS = int(os.environ.get('RUN_DIR_TIMEOUT_SECONDS', '300'))

app = Flask(__name__)
app.config.update(MAX_CONTENT_LENGTH=MAX_CONTENT, SECRET_KEY=os.urandom(16))
os.makedirs(RUNS_DIR, exist_ok=True)

# In-memory map of run_id -> last_seen (epoch seconds). This is used by the
# background cleaner to determine which run directories can be deleted when the
# client abandons the page. We store timestamps here as a best-effort signal.
# Note: on restart this map is empty; the cleanup thread will still remove
# stale directories based on mtime if needed.
_run_last_seen = {}


def touch_run(run_id: str) -> None:
    """Record that a client is still interested in the given run_id."""
    _run_last_seen[run_id] = int(time.time())


def _cleanup_worker():
    """Background thread that removes run directories not seen recently.

    It iterates RUNS_DIR every minute and removes any run directory whose
    last-seen timestamp is older than RUN_DIR_TIMEOUT_SECONDS. If a run_id is
    not present in the in-memory map, the directory's mtimes are used as a
    fallback to estimate last activity.
    """
    while True:
        try:
            now = int(time.time())
            for name in os.listdir(RUNS_DIR):
                run_dir = os.path.join(RUNS_DIR, name)
                if not os.path.isdir(run_dir):
                    continue
                last_seen = _run_last_seen.get(name)
                if last_seen is None:
                    # fallback: use latest mtime inside the directory
                    try:
                        mt = max(os.path.getmtime(os.path.join(run_dir, p))
                                 for p in os.listdir(run_dir))
                        last_seen = int(mt)
                    except Exception:
                        # if we can't stat contents, use dir mtime
                        last_seen = int(os.path.getmtime(run_dir))

                # If the run finished (RAiSD report present) don't clean up
                # here; removal is left to the client or explicit cleanup.
                if os.path.exists(os.path.join(run_dir, RAISD_REPORT)):
                    continue

                if (now - last_seen) > RUN_DIR_TIMEOUT_SECONDS:
                    try:
                        shutil.rmtree(run_dir)
                        _run_last_seen.pop(name, None)
                        app.logger.info(f"Cleaned stale run dir: {run_dir}")
                    except Exception as e:
                        app.logger.warning(f"Failed to remove stale run dir {run_dir}: {e}")
        except Exception as e:
            app.logger.exception(f"Run dir cleanup thread error: {e}")
        time.sleep(60)


# Start cleanup thread as a daemon so it doesn't block process exit.
try:
    t = threading.Thread(target=_cleanup_worker, daemon=True)
    t.start()
except Exception:
    app.logger.warning("Failed to start run-dir cleanup thread")

# Executor for running analysis jobs concurrently. Default is number of CPUs
# or 2 if unavailable. Can be overridden with MAX_CONCURRENT_ANALYSES env var.
MAX_CONCURRENT = int(os.environ.get('MAX_CONCURRENT_ANALYSES', str(max(1, (os.cpu_count() or 2)))))
_executor = ThreadPoolExecutor(max_workers=MAX_CONCURRENT)
# Map run_id -> Future so we can track running jobs (best-effort)
_run_futures = {}
_futures_lock = threading.Lock()

# Tools are expected inside the **conda env** PATH
BCFTOOLS = shutil.which("bcftools") or "bcftools"
RAISD_AI = shutil.which("RAiSD-AI") or os.environ.get("RAISD_AI", "RAiSD-AI")
RAISD_AI_ZLIB = shutil.which("RAiSD-AI-ZLIB") or os.environ.get("RAISD_AI_ZLIB", "RAiSD-AI-ZLIB")

# ---------------------- bcftools plugin helper ----------------------

def _bcftools_env():
    env = os.environ.copy()
    if "BCFTOOLS_PLUGINS" not in env:
        # conda installs plugins here typically
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            guess = os.path.join(conda_prefix, "libexec", "bcftools")
            if os.path.isdir(guess):
                env["BCFTOOLS_PLUGINS"] = guess
    return env

# --------------------------- Small helpers --------------------------

def _exe_available(prog: str) -> bool:
    return bool(shutil.which(prog) or (os.path.isabs(prog) and os.access(prog, os.X_OK)))


def choose_raisd_exe(vcf_path: str) -> str:
    p = vcf_path.lower()
    if p.endswith('.bcf') or p.endswith('.vcf.gz'):
        return RAISD_AI_ZLIB
    return RAISD_AI

def allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT


def append_log(run_dir: str, msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        with open(os.path.join(run_dir, "run.log"), "a", encoding="utf-8") as f:
            f.write(f"[{ts}] {msg}\n")
    except Exception:
        # Logging should never break the flow
        pass


def _needs_index(path: str) -> Tuple[bool, Optional[str]]:
    """
    Determine whether an index is expected for the given HTS file and which
    pathname we'd expect to exist.

    Returns (needs_index, expected_index_path or None).
    """
    p = path.lower()
    if p.endswith(".bcf"):
        idx = path + ".csi"
        return (not os.path.exists(idx), idx)
    if p.endswith(".vcf.gz"):
        # prefer CSI; accept TBI as well
        csi = path + ".csi"
        tbi = path + ".tbi"
        return (not (os.path.exists(csi) or os.path.exists(tbi)), csi)
    return (False, None)


def ensure_hts_index(path: str, bcftools: str = BCFTOOLS) -> Optional[str]:
    """
    Ensure a CSI/TBI index exists for .bcf or .vcf.gz. Returns the index path
    if created or already present; otherwise None.
    """
    needs, idx = _needs_index(path)
    if not needs:
        return idx
    # Create CSI index to be safe for large coordinates
    cmd = [bcftools, "index", "-f", "-c", path]
    proc = subprocess.run(cmd, env=_bcftools_env(), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to index file with bcftools: {proc.stderr}\nCMD: {' '.join(cmd)}")
    if idx and not os.path.exists(idx):
        # bcftools might have produced a .tbi for vcf.gz if -c was ignored; check it
        tbi = path + ".tbi"
        if os.path.exists(tbi):
            return tbi
        raise FileNotFoundError(f"Indexing reported success but index not found for {path}")
    return idx


def norm(x):
    x = np.asarray(x, float)
    s = x.sum()
    return x / s if s > 0 else x


def jsd(p, q, eps=1e-12):
    p = norm(np.clip(p, eps, None))
    q = norm(np.clip(q, eps, None))
    m = 0.5 * (p + q)
    return float(0.5 * ((p * np.log(p / m)).sum() + (q * np.log(q / m)).sum()) / np.log(2.0))


def load_sfs_matrix(path):
    df = pd.read_csv(path, index_col=0).astype(float)
    cols = sorted(df.columns, key=lambda c: int(c))   # ensure numeric order
    df = df[cols].div(df.sum(axis=1), axis=0)         # row-normalize to probabilities
    return df


def infer_samples_and_ploidy(vcf_path, max_records=200):
    vf = pysam.VariantFile(vcf_path)
    n_ind = len(vf.header.samples)

    counts = Counter()
    for i, rec in enumerate(vf.fetch()):
        if i >= max_records:
            break
        for smp in rec.samples.values():
            gt = smp.get("GT")
            if not gt:
                continue
            L = sum(1 for a in gt if a is not None)
            if L > 0:
                counts[L] += 1

    ploidy = counts.most_common(1)[0][0] if counts else 2
    mixed = (len(counts) > 1)
    return n_ind * ploidy, ploidy, mixed, dict(counts)


def sfs_from_bcftools(vcf_path, samples, bcftools=BCFTOOLS):
    # Requires bcftools; sums AC over multiallelic ALTs; produces SFS[0..samples]
    cmd = (
        f"{bcftools} +fill-tags {shlex.quote(vcf_path)} -Ou -- -t AC | "
        f"{bcftools} query -f '%AC\\n' -"
    )
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1, env=_bcftools_env())
    assert p.stdout is not None  # for type checkers; PIPE ensures stdout is set
    sfs = np.zeros(samples + 1, np.int64)
    for line in p.stdout:
        parts = [x for x in line.strip().split(",") if x and x != "."]
        if parts:
            ac = sum(int(x) for x in parts)
            if 0 <= ac <= samples:
                sfs[ac] += 1
    p.wait()
    if p.returncode != 0:
        raise RuntimeError("bcftools pipeline failed; check bcftools and plugins path")
    return sfs

# ------------------ hypergeometric projection (expected) ------------------

def _pre_lgamma(N):
    G = np.zeros(N + 2, float)
    for m in range(1, N + 2):
        G[m] = lgamma(m)
    return G


def _logC(n, k, G):
    return G[n + 1] - G[k + 1] - G[n - k + 1]


def _pmf_i(i, nproj, N, G):
    jmin = max(0, nproj + i - N)
    jmax = min(i, nproj)
    js = np.arange(jmin, jmax + 1, dtype=int)
    denom = _logC(N, i, G)
    pmf = np.exp(_logC(nproj, js, G) + _logC(N - nproj, i - js, G) - denom)
    s = pmf.sum()
    return jmin, jmax, (pmf / s if s > 0 else pmf)


def project_expected_df(df, target):
    N = max(int(c) for c in df.columns)          # reference chromosome count from table
    if target > N:
        raise ValueError(f"target ({target}) > reference N ({N}); projection is downsampling only.")
    G, cache = _pre_lgamma(N), {}
    rows = []
    for _, row in df.iterrows():
        p_full = np.zeros(N + 1, float)          # bins 0..N
        for c, v in row.items():
            p_full[int(c)] = float(v)
        z = np.zeros(target + 1, float)          # bins 0..target
        for i, w in enumerate(p_full):
            if w <= 0:
                continue
            if i not in cache:
                cache[i] = _pmf_i(i, target, N, G)
            jmin, jmax, pmf = cache[i]
            z[jmin:jmax + 1] += w * pmf
        rows.append(norm(z[1:]))                 # drop bin 0 -> bins 1..target
    return pd.DataFrame(rows, index=df.index, columns=[str(k) for k in range(1, target + 1)])

# ----------------------- RAiSD parsing & plotting ------------------------

HEADER = [r"$\mu_{Var}$", r"$\mu_{SFS}$", r"$\mu_{LD}$", r"$\mu$", r"$sweep_{TR}$", r"$\mu_{var}^{sweep_{TR}}$"]


def read_raisd_report(filename: str) -> pd.DataFrame:
    with open(filename, 'r') as f:
        content = f.read()
    chunks = [chunk.strip() for chunk in content.split('//') if chunk.strip()]
    dfs = []
    for chunk in chunks:
        raw_lines = [ln for ln in chunk.split('\n') if ln.strip()]
        data_lines = [
            ln for ln in raw_lines
            if ln and ln.strip()[0].isdigit() and len(ln.split()) >= 9
        ]
        if not data_lines:
            continue
        df = pd.read_csv(StringIO('\n'.join(data_lines)), header=None, sep=r'\s+', engine='python', index_col=0)
        if df.shape[1] > len(HEADER):
            df = df.iloc[:, -len(HEADER):]
        elif df.shape[1] < len(HEADER):
            continue
        df.columns = HEADER
        dfs.append(df)
    if not dfs:
        raise ValueError("No parsable data blocks found in RAiSD_Report.results")
    return pd.concat(dfs)


def plot_sweep_tr(df: pd.DataFrame, out_png: str):
    plt.figure(figsize=(10, 6))
    plt.plot(df.index, df[r"$\mu_{var}^{sweep_{TR}}$"], label=r"$\mu_{var}^{sweep_{TR}}$")
    plt.xlabel('Position')
    plt.ylabel(r"$\mu_{var}^{sweep_{TR}}$")
    plt.title(r"Plot of $\mu_{var}^{sweep_{TR}}$")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def _sanitize_metric_name(m: str) -> str:
    # make a filesystem- and url-safe short name for a metric header
    # remove non-alphanumeric and collapse spaces/brace characters
    import re
    s = re.sub(r'[^0-9A-Za-z]+', '_', m)
    s = re.sub(r'_+', '_', s).strip('_')
    return s or 'metric'


def plot_metric(df: pd.DataFrame, metric: str, out_png: str):
    # Generic plotter for a single metric series in the RAiSD report
    plt.figure(figsize=(10, 6))
    if metric not in df.columns:
        # nothing to plot
        plt.text(0.5, 0.5, f'No data for {metric}', horizontalalignment='center', verticalalignment='center')
    else:
        plt.plot(df.index, df[metric], label=metric)
        plt.xlabel('Position')
        plt.ylabel(metric)
        plt.title(metric)
        plt.legend()
        plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

# ------------------------------ Pipeline ------------------------------

def compute_best_match(species: str, vcf_path: str):
    sfs_table = os.path.join(DATA_DIR, species, 'sfs.csv')
    if not os.path.exists(sfs_table):
        raise FileNotFoundError(f"Missing SFS table: {sfs_table}")
    df = load_sfs_matrix(sfs_table)
    target, ploidy, mixed, counts = infer_samples_and_ploidy(vcf_path)
    proj = project_expected_df(df, target)
    sfs = sfs_from_bcftools(vcf_path, target)
    q = norm(sfs[1:])
    best_pop, best_jsd = min(((pop, jsd(proj.loc[pop].values, q)) for pop in proj.index), key=lambda x: x[1])
    demographic_model, population = best_pop.rsplit('_', 1)
    return {
        'target': target,
        'ploidy': ploidy,
        'mixed_ploidy': mixed,
        'ploidy_counts': counts,
        'demographic_model': demographic_model,
        'population': population,
        'best_jsd': best_jsd,
    }


def run_raisd(run_dir: str, model_path: str, vcf_path: str, chromosome: int, grid=300, ra_exe: Optional[str] = None):
    # Clean old grid dir if RAiSD creates it
    grid_dir = os.path.join(run_dir, RAISD_GRIDDIR)
    if os.path.isdir(grid_dir):
        shutil.rmtree(grid_dir, ignore_errors=True)

    ra_bin = ra_exe or choose_raisd_exe(vcf_path)
    # Pass raw path strings as separate argv items. Do NOT use shlex.quote here
    # because subprocess.run is called with a list (no shell) and quoting
    # would insert literal quote characters into the argument values.
    cmd = [
        ra_bin, '-n', 'results',
        '-mdl', model_path,
        '-f',
        '-op', 'SWP-SCN',
        '-I', vcf_path,
        '-frm',
        '-G', str(grid),
        '-pci', '1', '1',
        '-R'
    ]
    append_log(run_dir, f"Running RAiSD: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True)
    # Persist stdout/stderr for inspection
    try:
        if proc.stdout:
            with open(os.path.join(run_dir, 'raisd_stdout.log'), 'w', encoding='utf-8') as f:
                f.write(proc.stdout)
        if proc.stderr:
            with open(os.path.join(run_dir, 'raisd_stderr.log'), 'w', encoding='utf-8') as f:
                f.write(proc.stderr)
    except Exception:
        pass
    append_log(run_dir, f"RAiSD finished with code {proc.returncode}")
    if proc.returncode != 0:
        raise RuntimeError(f"RAiSD failed: {proc.stderr}\n\nCMD: {' '.join(cmd)}")

    report_path = os.path.join(run_dir, RAISD_REPORT)
    if not os.path.exists(report_path):
        raise FileNotFoundError(f"Expected RAiSD report not found at {report_path}")
    return report_path

# ------------------------------- Routes -------------------------------

@app.route('/')
def index():
    species_options = sorted([d for d in os.listdir(DATA_DIR) if os.path.isdir(os.path.join(DATA_DIR, d))]) or ["HomSap"]
    return render_template('index.html', species_options=species_options)


@app.route('/analyze', methods=['POST'])
def analyze():
    # allow analysis to reference a previously-uploaded staged file using run_id
    run_id = request.form.get('run_id', '').strip()
    uploaded_filename = request.form.get('uploaded_filename', '').strip()

    species    = request.form.get('species', 'HomSap').strip()
    chromosome = request.form.get('chromosome', '').strip()
    grid_str   = request.form.get('grid', '300').strip()
    file       = request.files.get('vcf')

    if not species:
        flash("Species is required.")
        return redirect(url_for('index'))

    try:
        chromosome = int(chromosome)
    except Exception:
        flash("Chromosome must be an integer.")
        return redirect(url_for('index'))

    # Validate grid parameter (RAiSD -G)
    try:
        grid = int(grid_str)
        if grid < 1:
            raise ValueError('grid must be >= 1')
    except Exception:
        flash("Grid must be a positive integer.")
        return redirect(url_for('index'))

    if not file or file.filename == '':
        # if a staged upload exists, prefer that file
        if run_id and uploaded_filename:
            run_dir = os.path.join(RUNS_DIR, run_id)
            vcf_path = os.path.join(run_dir, uploaded_filename)
            if not os.path.exists(vcf_path):
                flash("Staged upload not found. Please upload again.")
                return redirect(url_for('index'))
        else:
            flash("Please upload a VCF/BCF file.")
            return redirect(url_for('index'))

    # determine filename from uploaded file or staged upload
    if file and getattr(file, 'filename', ''):
        fname = file.filename
    else:
        fname = uploaded_filename

    if not fname or not allowed_file(fname):
        flash("Unsupported file type. Use .vcf, .vcf.gz, or .bcf")
        return redirect(url_for('index'))

    # If the client provided a staged run_id and no fresh file, reuse that run dir
    if (not file or not getattr(file, 'filename', '')) and run_id and uploaded_filename:
        run_dir = os.path.join(RUNS_DIR, run_id)
        if not os.path.isdir(run_dir):
            flash("Staged run directory not found. Please upload again.")
            return redirect(url_for('index'))
        filename = uploaded_filename
        vcf_path = os.path.join(run_dir, filename)
        append_log(run_dir, f"Using staged upload: {filename}")
        touch_run(run_id)
    else:
        # create a new run dir and save the uploaded file
        run_id = uuid.uuid4().hex[:12]
        run_dir = os.path.join(RUNS_DIR, run_id)
        os.makedirs(run_dir, exist_ok=True)
        filename = secure_filename(fname)
        vcf_path = os.path.join(run_dir, filename)
        if file and getattr(file, 'filename', ''):
            file.save(vcf_path)
            append_log(run_dir, f"Uploaded file saved: {filename} ({os.path.getsize(vcf_path)} bytes)")
            # record run seen right after creation
            touch_run(run_id)
        else:
            # This shouldn't normally happen because earlier checks ensured a file is present
            flash("No file available for analysis.")
            return redirect(url_for('index'))

    # Start the long-running analysis in a background thread and return
    # immediately. This avoids HTTP timeouts (reverse proxies or browsers)
    # when processing large uploads which can take many minutes.
    def _analyze_worker(run_dir, filename, vcf_path, species, chromosome, grid):
        meta = {
            'run_id': os.path.basename(run_dir),
            'species': species,
            'chromosome': chromosome,
            'grid': grid,
            'uploaded_name': filename,
            'uploaded_url': None,
            'plot_url': None,
            'report_url': None,
            'match': None,
            'error': None,
        }
        try:
            if shutil.which(BCFTOOLS) is None and BCFTOOLS == "bcftools":
                raise EnvironmentError("bcftools not found in PATH (activate your conda env)")
            ra_exe = choose_raisd_exe(vcf_path)
            append_log(run_dir, f"Selected RAiSD executable: {ra_exe}")
            if not _exe_available(ra_exe):
                exe_name = os.path.basename(ra_exe)
                if exe_name.endswith("RAiSD-AI-ZLIB"):
                    raise EnvironmentError("RAiSD-AI-ZLIB not found; build it and add to PATH or set RAISD_AI_ZLIB env var")
                else:
                    raise EnvironmentError("RAiSD-AI not found; build it and add to PATH or set RAISD_AI env var")

            # Auto-index for BCF/VCF.GZ
            try:
                append_log(run_dir, "Checking/creating index for uploaded file if needed...")
                idx_path = ensure_hts_index(vcf_path)
                if idx_path:
                    append_log(run_dir, f"Index present at: {idx_path}")
                else:
                    append_log(run_dir, "No index required for this file type.")
            except Exception as ie:
                if _needs_index(vcf_path)[0]:
                    raise
                append_log(run_dir, f"Indexing skipped (plain VCF) or warning: {ie}")

            match = compute_best_match(species, vcf_path)
            meta['match'] = match
            append_log(run_dir, (
                f"Best match -> model: {match['demographic_model']}, population: {match['population']}, "
                f"JSD: {match['best_jsd']:.6f}, target chromosomes: {match['target']}, "
                f"ploidy: {match['ploidy']}, mixed_ploidy: {match['mixed_ploidy']}"
            ))
            model_path = os.path.join(
                DATA_DIR, species, match['demographic_model'], match['population'], str(chromosome), 'RAiSD_Model.model'
            )
            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Model file missing: {model_path}")

            report_path = run_raisd(run_dir, model_path, vcf_path, chromosome, grid=grid, ra_exe=ra_exe)
            df = read_raisd_report(report_path)
            # Generate one plot per metric in HEADER and record their paths
            plots = []
            for h in HEADER:
                short = _sanitize_metric_name(h)
                fname = f"plot_{short}.png"
                plot_path = os.path.join(run_dir, fname)
                try:
                    # prefer the specific metric plotter
                    plot_metric(df, h, plot_path)
                except Exception:
                    try:
                        plot_sweep_tr(df, plot_path)
                    except Exception:
                        # if plotting fails, create a placeholder PNG
                        try:
                            plt.figure(figsize=(6, 3))
                            plt.text(0.5, 0.5, 'Plot failed', horizontalalignment='center', verticalalignment='center')
                            plt.tight_layout()
                            plt.savefig(plot_path, dpi=100)
                            plt.close()
                        except Exception:
                            pass
                plots.append({'metric': h, 'filename': fname})
            append_log(run_dir, f"Generated plots and parsed report: {RAISD_REPORT}")

            grid_dir = os.path.join(run_dir, RAISD_GRIDDIR)
            if os.path.isdir(grid_dir):
                shutil.rmtree(grid_dir, ignore_errors=True)
                append_log(run_dir, f"Cleaned grid directory: {RAISD_GRIDDIR}")

            # Build URLs that the result page will use
            # url_for needs an application context when called from a
            # background thread; push one here so URL generation doesn't
            # fail when multiple analyses run concurrently.
            try:
                with app.app_context():
                    meta['plot_url'] = url_for('runs_file', run_id=meta['run_id'], filename=PLOTS_FILE)
                    meta['report_url'] = url_for('runs_file', run_id=meta['run_id'], filename=RAISD_REPORT)
                    # build per-metric plot URLs
                    meta['plots'] = []
                    try:
                        for p in plots:
                            fname = p.get('filename')
                            if not fname:
                                continue
                            meta['plots'].append({'metric': p.get('metric'), 'url': url_for('runs_file', run_id=meta['run_id'], filename=fname)})
                    except Exception:
                        # if URL building fails here, final endpoint will synthesize
                        pass
            except Exception:
                # fallback to None; final endpoint can synthesize URLs later
                meta['plot_url'] = None
                meta['report_url'] = None

            # Remove some large temporary files to save space
            try:
                input_path = os.path.join(run_dir, filename)
                for p in [input_path, input_path + '.csi', input_path + '.tbi']:
                    try:
                        if os.path.exists(p):
                            os.remove(p)
                    except Exception:
                        pass
                for p in ['raisd_stdout.log', 'raisd_stderr.log', 'RAiSD_Info.results']:
                    try:
                        fp = os.path.join(run_dir, p)
                        if os.path.exists(fp):
                            os.remove(fp)
                    except Exception:
                        pass
            except Exception:
                pass

        except Exception as e:
            append_log(run_dir, f"Error: {e}")
            meta['error'] = str(e)
        finally:
            # persist meta so the final-render endpoint can read it
            try:
                with open(os.path.join(run_dir, 'result_meta.json'), 'w', encoding='utf-8') as mf:
                    mf.write(json.dumps(meta))
            except Exception:
                pass

    # submit the analysis to the shared executor so multiple analyses can
    # run concurrently without spawning unbounded threads.
    try:
        fut = _executor.submit(_analyze_worker, run_dir, filename, vcf_path, species, chromosome, grid)
        with _futures_lock:
            _run_futures[run_id] = fut
    except Exception as e:
        append_log(run_dir, f"Failed to submit analysis job: {e}")
        app.logger.exception("Failed to submit analysis job")
        # Prefer returning JSON so client-side JS gets a structured error.
        return json.dumps({'error': f'Failed to start analysis: {e}'}), 500, {'Content-Type': 'application/json'}

    # Return JSON immediately with the run_id; the client will poll /runs/<run_id>/tail
    try:
        return json.dumps({'run_id': run_id}), 200, {'Content-Type': 'application/json'}
    except Exception as e:
        app.logger.exception("Failed to return analyze response")
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}


@app.route('/runs/<run_id>/<path:filename>')
def runs_file(run_id, filename):
    return send_from_directory(os.path.join(RUNS_DIR, run_id), filename, as_attachment=False)


@app.route('/runs/<run_id>/tail')
def runs_tail(run_id):
    """Return the tail of the run log or RAiSD_Info.results as JSON along with a 'done' flag.

    done is True when the RAiSD report file exists in the run directory (analysis finished).
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}

    # mark as seen (client is polling the tail endpoint)
    touch_run(run_id)

    # Prefer RAiSD_Info.results if present (it contains detailed pipeline output), otherwise run.log
    info_path = os.path.join(run_dir, 'RAiSD_Info.results')
    log_path = os.path.join(run_dir, 'run.log')
    source = None
    if os.path.exists(info_path):
        source = info_path
    elif os.path.exists(log_path):
        source = log_path

    lines = ''
    if source:
        try:
            with open(source, 'r', encoding='utf-8', errors='replace') as f:
                data = f.read()
            # return last 1000 characters to keep payload small
            lines = data[-5000:]
        except Exception as e:
            lines = f"(failed to read log: {e})"
    else:
        lines = '(no log available yet)'

    report_path = os.path.join(run_dir, RAISD_REPORT)
    done = os.path.exists(report_path)
    return json.dumps({'lines': lines, 'done': bool(done)}), 200, {'Content-Type': 'application/json'}


@app.route('/runs/<run_id>/cleanup', methods=['POST'])
def runs_cleanup(run_id):
    """Remove the temporary run directory for the given run_id.

    This is intended to be called by the client when the user navigates away
    from the results page. The endpoint is deliberately simple: if the
    directory exists it will be removed; otherwise a 404 is returned.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    try:
        shutil.rmtree(run_dir)
        _run_last_seen.pop(run_id, None)
        return json.dumps({'status': 'removed'}), 200, {'Content-Type': 'application/json'}
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}


@app.route('/runs/<run_id>/final')
def runs_final(run_id):
    """Return the final rendered result page when analysis completes.

    If the analysis is still running, return 202 Accepted with a short
    JSON payload so the client can retry. When ready, render `result.html`.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}

    # Prefer result_meta.json written by the background worker
    meta_path = os.path.join(run_dir, 'result_meta.json')
    report_path = os.path.join(run_dir, RAISD_REPORT)

    # mark as seen
    touch_run(run_id)

    meta = {}
    if os.path.exists(meta_path):
        try:
            with open(meta_path, 'r', encoding='utf-8') as mf:
                meta = json.load(mf)
        except Exception:
            meta = {}

    # If report exists (analysis completed successfully) or meta has error/result, render final page
    if os.path.exists(report_path) or meta.get('report_url') or meta.get('error'):
        # If result_meta.json didn't include match/urls, try to synthesize
        match = meta.get('match') or {}
        # Templates expect attribute access (match.best_jsd) — convert dict to
        # a SimpleNamespace so dot-notation works when meta['match'] is a dict.
        if isinstance(match, dict):
            try:
                match = SimpleNamespace(**match)
            except Exception:
                # Fallback: leave as dict if conversion fails
                pass
        uploaded_name = meta.get('uploaded_name') or ''
        plot_url = meta.get('plot_url') or url_for('runs_file', run_id=run_id, filename=PLOTS_FILE)
        report_url = meta.get('report_url') or url_for('runs_file', run_id=run_id, filename=RAISD_REPORT)
        try:
            return render_template(
                'result.html',
                run_id=run_id,
                species=meta.get('species', ''),
                chromosome=meta.get('chromosome', ''),
                grid=meta.get('grid', ''),
                match=match or {},
                plot_url=plot_url,
                report_url=report_url,
                plots=meta.get('plots', []),
                log_url=None,
                uploaded_name=uploaded_name,
                uploaded_url=None,
            )
        except Exception as e:
            # Log and return a structured JSON error so client-side JS can
            # report the problem instead of getting an HTML 500 that is
            # harder to surface to the user.
            app.logger.exception(f"Failed to render final result for run {run_id}: {e}")
            return json.dumps({'error': f'Failed to render result: {e}'}), 500, {'Content-Type': 'application/json'}

    # Still running
    return json.dumps({'status': 'running'}), 202, {'Content-Type': 'application/json'}



@app.route('/runs/<run_id>/heartbeat', methods=['POST'])
def runs_heartbeat(run_id):
    """Client should POST to this endpoint periodically while the results
    page is open. That keeps the temporary run directory alive. If the
    client leaves the page, the background cleaner will remove the run
    after RUN_DIR_TIMEOUT_SECONDS.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    touch_run(run_id)
    return json.dumps({'status': 'ok'}), 200, {'Content-Type': 'application/json'}


@app.route('/healthz')
def healthz():
    return {"status": "ok"}


@app.route('/chromosomes')
def get_chromosomes():
    species = request.args.get('species', 'HomSap').strip()
    species_dir = os.path.join(DATA_DIR, species)
    chrom_file = os.path.join(species_dir, 'chromosomes.txt')
    if not os.path.isfile(chrom_file):
        return json.dumps({'error': 'chromosomes.txt not found'}), 404, {'Content-Type': 'application/json'}
    try:
        with open(chrom_file, 'r') as f:
            chromosomes = [line.strip() for line in f if line.strip()]
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}
    return json.dumps({
        'species': species,
        'chromosomes': chromosomes
    }), 200, {'Content-Type': 'application/json'}


@app.route('/upload', methods=['POST'])
def upload():
    """Stage the uploaded file immediately and return a run_id + filename for later analysis."""
    file = request.files.get('vcf')
    if not file or not getattr(file, 'filename', ''):
        return json.dumps({'error': 'No file uploaded'}), 400, {'Content-Type': 'application/json'}

    fname = secure_filename(getattr(file, 'filename', '') or '')
    run_id = uuid.uuid4().hex[:12]
    run_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(run_dir, exist_ok=True)
    vcf_path = os.path.join(run_dir, fname)
    file.save(vcf_path)
    append_log(run_dir, f"Staged upload saved: {fname} ({os.path.getsize(vcf_path)} bytes)")
    # record that the run was just created / seen
    touch_run(run_id)
    return json.dumps({'run_id': run_id, 'filename': fname}), 200, {'Content-Type': 'application/json'}


if __name__ == '__main__':
    # When running via `python app.py` we enable threaded mode and read
    # configuration from environment variables so the development entry
    # point is less likely to drop long-running background work.
    # For production use, prefer running under gunicorn/uwsgi with
    # an appropriate timeout and worker configuration (see README).
    host = os.environ.get('HOST', '0.0.0.0')
    port = int(os.environ.get('PORT', '5000'))
    debug = bool(os.environ.get('FLASK_DEBUG', '0') == '1')
    threaded = bool(os.environ.get('FLASK_THREADED', '1') == '1')
    # disable reloader to avoid accidental duplicate background threads
    app.run(host=host, port=port, debug=debug, threaded=threaded, use_reloader=False)