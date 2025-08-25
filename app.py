#!/usr/bin/env python3
import os
import uuid
import shlex
import shutil
import subprocess
from io import StringIO, BytesIO
from typing import Optional, Tuple
from functools import lru_cache
from datetime import datetime
import glob
import gzip

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

from flask import (
    Flask, render_template, request, redirect, url_for,
    send_from_directory, flash, Response
)
from werkzeug.utils import secure_filename

# ------------------------------ Config ------------------------------
BASE_DIR      = os.path.abspath(os.path.dirname(__file__))
DATA_DIR      = os.path.join(BASE_DIR, "data")
RUNS_DIR      = os.path.join(BASE_DIR, "runs")
ALLOWED_EXT   = {"vcf", "gz", "bcf"}
MAX_CONTENT   = 2 * 1024 * 1024 * 1024  # 2 GB
RAISD_REPORT  = "RAiSD_Report.results"
RAISD_GRIDDIR = "RAiSD_Grid.results"

RUN_DIR_TIMEOUT_SECONDS = int(os.environ.get('RUN_DIR_TIMEOUT_SECONDS', '300'))
FINAL_WAIT_TIMEOUT_SECONDS = int(os.environ.get('FINAL_WAIT_TIMEOUT_SECONDS', '7200'))  # 2h
FINAL_WAIT_POLL_INTERVAL = float(os.environ.get('FINAL_WAIT_POLL_INTERVAL', '1.0'))     # s

app = Flask(__name__)
app.config.update(MAX_CONTENT_LENGTH=MAX_CONTENT, SECRET_KEY=os.urandom(16))
os.makedirs(RUNS_DIR, exist_ok=True)

_run_last_seen = {}

def touch_run(run_id: str) -> None:
    _run_last_seen[run_id] = int(time.time())

def _cleanup_worker():
    while True:
        try:
            now = int(time.time())
            for name in os.listdir(RUNS_DIR):
                run_dir = os.path.join(RUNS_DIR, name)
                if not os.path.isdir(run_dir):
                    continue
                last_seen = _run_last_seen.get(name)
                if last_seen is None:
                    try:
                        mt = max(os.path.getmtime(os.path.join(run_dir, p))
                                 for p in os.listdir(run_dir))
                        last_seen = int(mt)
                    except Exception:
                        last_seen = int(os.path.getmtime(run_dir))
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

try:
    t = threading.Thread(target=_cleanup_worker, daemon=True)
    t.start()
except Exception:
    app.logger.warning("Failed to start run-dir cleanup thread")

MAX_CONCURRENT = int(os.environ.get('MAX_CONCURRENT_ANALYSES', str(max(1, (os.cpu_count() or 2)))))
_executor = ThreadPoolExecutor(max_workers=MAX_CONCURRENT)
_run_futures = {}
_futures_lock = threading.Lock()

# Tool discovery
BCFTOOLS = shutil.which("bcftools") or "bcftools"
RAISD_AI = shutil.which("RAiSD-AI") or os.environ.get("RAISD_AI", "RAiSD-AI")
RAISD_AI_ZLIB = shutil.which("RAiSD-AI-ZLIB") or os.environ.get("RAISD_AI_ZLIB", "RAiSD-AI-ZLIB")

# ---------------------- bcftools plugin helper ----------------------
def _bcftools_env():
    env = os.environ.copy()
    if "BCFTOOLS_PLUGINS" not in env:
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
    needs_z = p.endswith('.bcf') or p.endswith('.vcf.gz')
    cand = []
    env_z = os.environ.get("RAISD_AI_ZLIB")
    env_ai = os.environ.get("RAISD_AI")
    if needs_z:
        cand += [env_z, "RAiSD-AI-ZLIB"]
    cand += [env_ai, "RAiSD-AI", "RAiSD"]
    for exe in cand:
        if not exe:
            continue
        if os.path.isabs(exe) and os.access(exe, os.X_OK):
            return exe
        w = shutil.which(exe)
        if w:
            return w
    raise EnvironmentError(
        "No RAiSD executable found. Install one of RAiSD-AI-ZLIB, RAiSD-AI, or RAiSD "
        "and ensure it’s on PATH or set RAISD_AI[_ZLIB]."
    )

def allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT

def append_log(run_dir: str, msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        with open(os.path.join(run_dir, "run.log"), "a", encoding="utf-8") as f:
            f.write(f"[{ts}] {msg}\n")
    except Exception:
        pass

def _needs_index(path: str) -> Tuple[bool, Optional[str]]:
    p = path.lower()
    if p.endswith(".bcf"):
        idx = path + ".csi"
        return (not os.path.exists(idx), idx)
    if p.endswith(".vcf.gz"):
        csi = path + ".csi"
        tbi = path + ".tbi"
        return (not (os.path.exists(csi) or os.path.exists(tbi)), csi)
    return (False, None)

def ensure_hts_index(path: str, bcftools: str = BCFTOOLS) -> Optional[str]:
    needs, idx = _needs_index(path)
    if not needs:
        return idx
    cmd = [bcftools, "index", "-f", "-c", path]
    proc = subprocess.run(cmd, env=_bcftools_env(), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to index file with bcftools: {proc.stderr}\nCMD: {' '.join(cmd)}")
    if idx and not os.path.exists(idx):
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
    cols = sorted(df.columns, key=lambda c: int(c))
    df = df[cols].div(df.sum(axis=1), axis=0)
    return df

def infer_samples_and_ploidy(vcf_path, max_records=200):
    vf = pysam.VariantFile(vcf_path)
    n_ind = len(vf.header.samples)
    from collections import Counter
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
    cmd = (
        f"{bcftools} +fill-tags {shlex.quote(vcf_path)} -Ou -- -t AC | "
        f"{bcftools} query -f '%AC\\n' -"
    )
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1, env=_bcftools_env())
    assert p.stdout is not None
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

# ---------------------- Gene annotation helpers ----------------------
ANNOTATION_DIR_TMPL = os.path.join(DATA_DIR, "{species}", "annotation", "{chr}.tsv")

@lru_cache(maxsize=64)
def _load_annotation(species: str, chromosome: str) -> pd.DataFrame:
    path = ANNOTATION_DIR_TMPL.format(species=species, chr=str(chromosome))
    if not os.path.exists(path):
        raise FileNotFoundError(f"Annotation file not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype={"chromosome": str})

    def _first_non_null(series):
        for v in series:
            if pd.notna(v):
                return v
        return np.nan

    gene_level = (
        df.sort_values(["gene_id", "start"]).
          groupby(["chromosome", "gene_id"], as_index=False).
          agg(start=("start", "min"),
               end=("end", "max"),
               gene_name=("gene_name", _first_non_null),
               biotype=("biotype", _first_non_null))
    )
    gene_level["label"] = gene_level["gene_name"].fillna(gene_level["gene_id"]).astype(str)
    gene_level["start"] = gene_level["start"].astype(int)
    gene_level["end"] = gene_level["end"].astype(int)
    return gene_level

def _genes_in_window(species: str, chromosome: str, start: int, end: int, biotype: Optional[str] = None) -> pd.DataFrame:
    lo, hi = (int(start), int(end)) if start <= end else (int(end), int(start))
    genes = _load_annotation(species, str(chromosome))
    genes_chr = genes[genes["chromosome"].astype(str) == str(chromosome)].copy()
    hits = genes_chr[(genes_chr["start"] <= hi) & (genes_chr["end"] >= lo)].copy()
    if biotype:
        hits = hits[hits["biotype"].astype(str).str.lower() == str(biotype).lower()]
    hits["overlap_bp"] = (np.minimum(hits["end"], hi) - np.maximum(hits["start"], lo) + 1).clip(lower=0)
    hits = hits.sort_values(["start", "end", "label"]).reset_index(drop=True)
    # Strand removed from output
    return hits[["label", "start", "end", "biotype", "overlap_bp"]]

# ------------------ projection helpers ------------------
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
    N = max(int(c) for c in df.columns)
    if target > N:
        raise ValueError(f"target ({target}) > reference N ({N}); projection is downsampling only.")
    G, cache = _pre_lgamma(N), {}
    rows = []
    for _, row in df.iterrows():
        p_full = np.zeros(N + 1, float)
        for c, v in row.items():
            p_full[int(c)] = float(v)
        z = np.zeros(target + 1, float)
        for i, w in enumerate(p_full):
            if w <= 0:
                continue
            if i not in cache:
                cache[i] = _pmf_i(i, target, N, G)
            jmin, jmax, pmf = cache[i]
            z[jmin:jmax + 1] += w * pmf
        rows.append(norm(z[1:]))
    return pd.DataFrame(rows, index=df.index, columns=[str(k) for k in range(1, target + 1)])

# ----------------------- RAiSD report + plotting ------------------------
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

def _sanitize_metric_name(m: str) -> str:
    import re
    s = re.sub(r'[^0-9A-Za-z]+', '_', m)
    s = re.sub(r'_+', '_', s).strip('_')
    return s or 'metric'

def plot_metric(df: pd.DataFrame, metric: str, out_png: str):
    plt.figure(figsize=(10, 6))
    if metric not in df.columns:
        plt.text(0.5, 0.5, f'No data for {metric}', ha='center', va='center')
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

def plot_metric_bytes(df: pd.DataFrame, metric: str, xmin: Optional[float] = None, xmax: Optional[float] = None) -> bytes:
    """Render a single-metric plot to PNG bytes; optional x-range can be supplied."""
    plt.figure(figsize=(10, 6))
    try:
        if metric not in df.columns:
            plt.text(0.5, 0.5, f'No data for {metric}', ha='center', va='center')
        else:
            plt.plot(df.index, df[metric], label=metric)
            plt.xlabel('Position')
            plt.ylabel(metric)
            plt.title(metric)
            plt.legend()
            plt.grid(True)
            if xmin is not None or xmax is not None:
                lo = xmin if xmin is not None else df.index.min()
                hi = xmax if xmax is not None else df.index.max()
                plt.xlim(lo, hi)
        plt.tight_layout()
        buf = BytesIO()
        plt.savefig(buf, dpi=150, format='png')
        buf.seek(0)
        return buf.read()
    finally:
        plt.close()

# ----------------- Cached report loader for interactive JSON plotting ---------
@lru_cache(maxsize=256)
def _cached_report(report_path: str) -> pd.DataFrame:
    return read_raisd_report(report_path)


@app.route('/runs/<run_id>/plot')
def runs_plot_dynamic(run_id):
    """Serve a single-metric plot PNG for a run with optional ?metric=...&xmin=...&xmax=..."""
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}

    metric = request.args.get('metric') or ''
    xmin = request.args.get('xmin')
    xmax = request.args.get('xmax')

    report_path = os.path.join(run_dir, RAISD_REPORT)
    if not os.path.exists(report_path):
        found = _find_and_normalize_report(run_dir, prefix='results')
        if found:
            report_path = found
    if not os.path.exists(report_path):
        return json.dumps({'error': 'RAiSD report not found'}), 404, {'Content-Type': 'application/json'}

    try:
        df = read_raisd_report(report_path)
    except Exception as e:
        app.logger.exception('Failed to read RAiSD report for dynamic plot')
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

    # normalize/resolve metric names: allow sanitized names (from filenames) to match report columns
    metric_param = (metric or '').strip()
    chosen_metric = None
    # precompute sanitized map
    san_map = { _sanitize_metric_name(col): col for col in df.columns }
    # direct match
    if metric_param in df.columns:
        chosen_metric = metric_param
    else:
        # try exact sanitized match of the param
        mnorm = _sanitize_metric_name(metric_param)
        if mnorm in san_map:
            chosen_metric = san_map[mnorm]
        else:
            # handle common filename forms like 'plot_mu_var' or 'plot_mu_var.png'
            s = metric_param
            if s.endswith('.png'):
                s = s.rsplit('/', 1)[-1]
            if s.startswith('plot_'):
                s2 = s[len('plot_'):]
            else:
                s2 = s
            s2norm = _sanitize_metric_name(s2)
            if s2norm in san_map:
                chosen_metric = san_map[s2norm]
    if chosen_metric is None:
        app.logger.debug(f"Dynamic plot: metric resolution failed. requested={metric_param}. available={list(df.columns)}. sanitized={list(san_map.keys())}")
    try:
        xmin_f = float(xmin) if xmin is not None and xmin != '' else None
        xmax_f = float(xmax) if xmax is not None and xmax != '' else None
    except Exception:
        return json.dumps({'error': 'xmin/xmax must be numeric if provided'}), 400, {'Content-Type': 'application/json'}

    try:
        # If no column matched, pass the original metric (plot_metric_bytes will render a 'No data' image)
        png = plot_metric_bytes(df, chosen_metric or metric_param, xmin=xmin_f, xmax=xmax_f)
        return (png, 200, {'Content-Type': 'image/png'})
    except Exception as e:
        app.logger.exception('Dynamic plot generation failed')
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/metric_data')
def runs_metric_data(run_id):
    """Return JSON series data for a metric: positions (index) and values.
    Optional query params: metric=, xmin=, xmax=, maxpoints= (downsample if large).
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    report_path = os.path.join(run_dir, RAISD_REPORT)
    if not os.path.exists(report_path):
        found = _find_and_normalize_report(run_dir, prefix='results')
        if found:
            report_path = found
    if not os.path.exists(report_path):
        return json.dumps({'error': 'report not found'}), 404, {'Content-Type': 'application/json'}
    metric = (request.args.get('metric') or '').strip()
    xmin = request.args.get('xmin')
    xmax = request.args.get('xmax')
    maxpoints = request.args.get('maxpoints')
    try:
        df = _cached_report(report_path)
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

    # resolve metric similar to runs_plot_dynamic
    metric_param = metric
    chosen_metric = None
    san_map = { _sanitize_metric_name(col): col for col in df.columns }
    if metric_param in df.columns:
        chosen_metric = metric_param
    else:
        mnorm = _sanitize_metric_name(metric_param)
        if mnorm in san_map:
            chosen_metric = san_map[mnorm]
        else:
            s = metric_param
            if s.endswith('.png'):
                s = s.rsplit('/', 1)[-1]
            if s.startswith('plot_'):
                s2 = s[len('plot_'):]
            else:
                s2 = s
            s2norm = _sanitize_metric_name(s2)
            if s2norm in san_map:
                chosen_metric = san_map[s2norm]
    if chosen_metric is None:
        # default to first column
        chosen_metric = df.columns[0]

    # slice by xmin/xmax if provided
    sub = df
    try:
        if xmin is not None and xmin != '':
            xmin_f = float(xmin)
            sub = sub[sub.index >= xmin_f]
        if xmax is not None and xmax != '':
            xmax_f = float(xmax)
            sub = sub[sub.index <= xmax_f]
    except Exception:
        pass

    xs = sub.index.values
    ys = sub[chosen_metric].values
    # Downsample if needed
    try:
        mp = int(maxpoints) if maxpoints else None
    except Exception:
        mp = None
    if mp and mp > 0 and xs.shape[0] > mp:
        # simple stride subsampling
        stride = max(1, xs.shape[0] // mp)
        xs = xs[::stride]
        ys = ys[::stride]

    payload = {
        'metric_requested': metric_param,
        'metric': chosen_metric,
        'sanitized_metric': _sanitize_metric_name(chosen_metric),
        'count': int(xs.shape[0]),
        'xmin': float(sub.index.min()) if sub.shape[0] else None,
        'xmax': float(sub.index.max()) if sub.shape[0] else None,
        'global_xmin': float(df.index.min()),
        'global_xmax': float(df.index.max()),
        'positions': xs.tolist(),
        'values': [float(v) for v in ys.tolist()],
    }
    return json.dumps(payload), 200, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/report_csv')
def runs_report_csv(run_id):
    """Serve the parsed RAiSD report (the DataFrame used for plotting) as a CSV file.
    Columns are the original metric headers; index is genomic position.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    report_path = os.path.join(run_dir, RAISD_REPORT)
    if not os.path.exists(report_path):
        found = _find_and_normalize_report(run_dir, prefix='results')
        if found:
            report_path = found
    if not os.path.exists(report_path):
        return json.dumps({'error': 'report not found'}), 404, {'Content-Type': 'application/json'}
    try:
        df = _cached_report(report_path)
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

    # Prepare CSV: include index label
    buf = StringIO()
    try:
        df.to_csv(buf, index_label='Position')
    except Exception as e:
        return json.dumps({'error': f'failed to serialize csv: {e}'}), 500, {'Content-Type': 'application/json'}
    csv_data = buf.getvalue()
    return Response(csv_data, mimetype='text/csv', headers={'Content-Disposition': 'attachment; filename=RAiSD_Report.csv'})

# Expected plots (we wait for all of them)
def _expected_plot_filenames():
    return [f"plot_{_sanitize_metric_name(h)}.png" for h in HEADER]

def _plots_exist(run_dir: str) -> bool:
    for fn in _expected_plot_filenames():
        if not os.path.exists(os.path.join(run_dir, fn)):
            return False
    return True

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

def _find_and_normalize_report(run_dir: str, prefix: str = "results") -> Optional[str]:
    patterns = [
        os.path.join(run_dir, f"RAiSD_Report.{prefix}"),
        os.path.join(run_dir, f"RAiSD_Report.{prefix}.gz"),
        os.path.join(run_dir, "RAiSD_Report.*"),
        os.path.join(run_dir, "**", f"RAiSD_Report.{prefix}"),
        os.path.join(run_dir, "**", f"RAiSD_Report.{prefix}.gz"),
        os.path.join(run_dir, "**", "RAiSD_Report.*"),
    ]
    candidates = []
    for pat in patterns:
        for path in glob.glob(pat, recursive=True):
            if os.path.basename(path).startswith("RAiSD_Report"):
                try:
                    candidates.append((os.path.getmtime(path), path))
                except Exception:
                    continue
    if not candidates:
        return None
    candidates.sort(reverse=True)
    src = candidates[0][1]
    dst = os.path.join(run_dir, RAISD_REPORT)
    if src.endswith(".gz"):
        with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:  # type: ignore[assignment]
            shutil.copyfileobj(f_in, f_out)  # type: ignore[arg-type]
        return dst
    if os.path.abspath(src) != os.path.abspath(dst):
        shutil.copyfile(src, dst)
    return dst

def run_raisd(run_dir: str, model_path: str, vcf_path: str, chromosome: int, grid=300, ra_exe: Optional[str] = None):
    grid_dir = os.path.join(run_dir, RAISD_GRIDDIR)
    if os.path.isdir(grid_dir):
        shutil.rmtree(grid_dir, ignore_errors=True)
    ra_bin = ra_exe or choose_raisd_exe(vcf_path)
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

    report_path = os.path.join(run_dir, RAISD_REPORT)
    if not os.path.exists(report_path):
        found = _find_and_normalize_report(run_dir, prefix="results")
        if found:
            append_log(run_dir, f"Normalized report to: {found}")
            report_path = found

    if proc.returncode != 0:
        err_tail = (proc.stderr or "").splitlines()[-30:]
        raise RuntimeError(
            "RAiSD failed with non-zero exit code.\n"
            f"CMD: {' '.join(cmd)}\n"
            f"stderr (tail):\n" + "\n".join(err_tail)
        )
    if not os.path.exists(report_path):
        std_tail = (proc.stdout or "").splitlines()[-30:]
        err_tail = (proc.stderr or "").splitlines()[-30:]
        raise FileNotFoundError(
            "RAiSD finished but no report was produced.\n"
            f"Searched under: {run_dir}\n"
            f"CMD: {' '.join(cmd)}\n"
            "Possible causes: incompatible flags, empty scan region, or a build that writes elsewhere.\n"
            "stdout (tail):\n" + "\n".join(std_tail) + "\n\n"
            "stderr (tail):\n" + "\n".join(err_tail)
        )
    return report_path

# ------------------------------- Routes -------------------------------
@app.route('/')
def index():
    species_options = sorted([d for d in os.listdir(DATA_DIR) if os.path.isdir(os.path.join(DATA_DIR, d))]) or ["HomSap"]
    return render_template('index.html', species_options=species_options)

@app.route('/analyze', methods=['POST'])
def analyze():
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
    try:
        grid = int(grid_str)
        if grid < 1:
            raise ValueError('grid must be >= 1')
    except Exception:
        flash("Grid must be a positive integer.")
        return redirect(url_for('index'))

    if not file or file.filename == '':
        if run_id and uploaded_filename:
            run_dir = os.path.join(RUNS_DIR, run_id)
            vcf_path = os.path.join(run_dir, uploaded_filename)
            if not os.path.exists(vcf_path):
                flash("Staged upload not found. Please upload again.")
                return redirect(url_for('index'))
        else:
            flash("Please upload a VCF/BCF file.")
            return redirect(url_for('index'))

    if file and getattr(file, 'filename', ''):
        fname = file.filename
    else:
        fname = uploaded_filename
    if not fname or not allowed_file(fname):
        flash("Unsupported file type. Use .vcf, .vcf.gz, or .bcf")
        return redirect(url_for('index'))

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
        run_id = uuid.uuid4().hex[:12]
        run_dir = os.path.join(RUNS_DIR, run_id)
        os.makedirs(run_dir, exist_ok=True)
        filename = secure_filename(fname)
        vcf_path = os.path.join(run_dir, filename)
        if file and getattr(file, 'filename', ''):
            file.save(vcf_path)
            append_log(run_dir, f"Uploaded file saved: {filename} ({os.path.getsize(vcf_path)} bytes)")
            touch_run(run_id)
        else:
            flash("No file available for analysis.")
            return redirect(url_for('index'))

    def _analyze_worker(run_dir, filename, vcf_path, species, chromosome, grid):
        meta = {
            'run_id': os.path.basename(run_dir),
            'species': species,
            'chromosome': chromosome,
            'grid': grid,
            'uploaded_name': filename,
            'uploaded_url': None,
            'report_url': None,
            'match': None,
            'model_path': None,
            'error': None,
            'plots_ready': False,
            'plots': [],
            'plot_url': None,  # will point to the currently selected/first plot
        }
        try:
            if not _exe_available(BCFTOOLS):
                raise EnvironmentError(f"bcftools not found or not executable: {BCFTOOLS} (activate your conda env)")
            ra_exe = choose_raisd_exe(vcf_path)
            append_log(run_dir, f"Selected RAiSD executable: {ra_exe}")
            if not _exe_available(ra_exe):
                raise EnvironmentError(f"RAiSD executable not usable: {ra_exe}")

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
            meta['model_path'] = model_path

            report_path = run_raisd(run_dir, model_path, vcf_path, chromosome, grid=grid, ra_exe=ra_exe)
            df = read_raisd_report(report_path)

            # Generate one plot per metric (NO separate sweep_TR summary plot)
            plots = []
            for h in HEADER:
                short = _sanitize_metric_name(h)
                fname = f"plot_{short}.png"
                plot_path = os.path.join(run_dir, fname)
                try:
                    plot_metric(df, h, plot_path)
                except Exception:
                    # minimal placeholder
                    try:
                        plt.figure(figsize=(6, 3))
                        plt.text(0.5, 0.5, 'Plot failed', ha='center', va='center')
                        plt.tight_layout()
                        plt.savefig(plot_path, dpi=100)
                        plt.close()
                    except Exception:
                        pass
                plots.append({'metric': h, 'filename': fname})

            # URLs for template
            try:
                with app.app_context():
                    meta['report_url'] = url_for('runs_file', run_id=meta['run_id'], filename=RAISD_REPORT)
                    meta['plots'] = []
                    for p in plots:
                        fname = p.get('filename')
                        if not fname:
                            continue
                        meta['plots'].append({
                            'metric': p.get('metric'),
                            'filename': fname,
                            'url': url_for('runs_file', run_id=meta['run_id'], filename=fname)
                        })
                    # default selection: put sweep_TR (sanitized) first if present
                    preferred = _sanitize_metric_name(r"$\mu_{var}^{sweep_{TR}}$")
                    i = next((k for k, pl in enumerate(meta['plots']) if _sanitize_metric_name(pl['metric']) == preferred), 0)
                    if i != 0:
                        meta['plots'].insert(0, meta['plots'].pop(i))
                    meta['plot_url'] = meta['plots'][0]['url'] if meta['plots'] else None
            except Exception:
                meta['report_url'] = None

            # mark plots ready
            meta['plots_ready'] = True
            append_log(run_dir, f"Generated plots and parsed report: {RAISD_REPORT}")

            # clean big temporaries
            try:
                input_path = os.path.join(run_dir, filename)
                for pth in [input_path, input_path + '.csi', input_path + '.tbi']:
                    try:
                        if os.path.exists(pth):
                            os.remove(pth)
                    except Exception:
                        pass
                for pth in ['raisd_stdout.log', 'raisd_stderr.log', 'RAiSD_Info.results']:
                    try:
                        fp = os.path.join(run_dir, pth)
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
            try:
                with open(os.path.join(run_dir, 'result_meta.json'), 'w', encoding='utf-8') as mf:
                    mf.write(json.dumps(meta))
            except Exception:
                pass

    try:
        fut = _executor.submit(_analyze_worker, run_dir, filename, vcf_path, species, chromosome, grid)
        with _futures_lock:
            _run_futures[run_id] = fut
    except Exception as e:
        append_log(run_dir, f"Failed to submit analysis job: {e}")
        app.logger.exception("Failed to submit analysis job")
        return json.dumps({'error': f'Failed to start analysis: {e}'}), 500, {'Content-Type': 'application/json'}

    try:
        return json.dumps({
            'run_id': run_id,
            'final_url': url_for('runs_final', run_id=run_id) + '?nowait=1',
            'wait_url': url_for('runs_final', run_id=run_id)
        }), 200, {'Content-Type': 'application/json'}
    except Exception as e:
        app.logger.exception("Failed to return analyze response")
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/<path:filename>')
def runs_file(run_id, filename):
    return send_from_directory(os.path.join(RUNS_DIR, run_id), filename, as_attachment=False)

@app.route('/runs/<run_id>/tail')
def runs_tail(run_id):
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    touch_run(run_id)
    info_path = os.path.join(run_dir, 'RAiSD_Info.results')
    log_path = os.path.join(run_dir, 'run.log')
    source = info_path if os.path.exists(info_path) else (log_path if os.path.exists(log_path) else None)
    lines = ''
    if source:
        try:
            with open(source, 'r', encoding='utf-8', errors='replace') as f:
                data = f.read()
            lines = data[-5000:]
        except Exception as e:
            lines = f"(failed to read log: {e})"
    else:
        lines = '(no log available yet)'
    done = os.path.exists(os.path.join(run_dir, RAISD_REPORT)) and _plots_exist(run_dir)
    return json.dumps({'lines': lines, 'done': bool(done)}), 200, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/cleanup', methods=['POST'])
def runs_cleanup(run_id):
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    try:
        shutil.rmtree(run_dir)
        _run_last_seen.pop(run_id, None)
        return json.dumps({'status': 'removed'}), 200, {'Content-Type': 'application/json'}
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

def _ready_or_error(run_dir: str) -> Tuple[bool, dict]:
    meta_path = os.path.join(run_dir, 'result_meta.json')
    report_path = os.path.join(run_dir, RAISD_REPORT)
    meta = {}
    if os.path.exists(meta_path):
        try:
            with open(meta_path, 'r', encoding='utf-8') as mf:
                meta = json.load(mf)
        except Exception:
            meta = {}
    ready = (os.path.exists(report_path) and (_plots_exist(run_dir) or meta.get('plots_ready')))
    if ready or meta.get('error'):
        return True, meta
    return False, meta

@app.route('/runs/<run_id>/final')
def runs_final(run_id):
    """
    Default: BLOCK until analysis is ready (report + ALL plots).
    Add ?nowait=1 for legacy non-blocking JSON while running.
    Optional: ?plot=<metric|sanitized|filename> to select a plot as primary.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    touch_run(run_id)

    if request.args.get('nowait', '0') == '1':
        report_path = os.path.join(run_dir, RAISD_REPORT)
        if not (os.path.exists(report_path) and _plots_exist(run_dir)):
            return json.dumps({'status': 'running'}), 202, {'Content-Type': 'application/json'}

    t0 = time.time()
    while True:
        ready, meta = _ready_or_error(run_dir)
        if ready:
            break
        if time.time() - t0 > FINAL_WAIT_TIMEOUT_SECONDS:
            return json.dumps({'error': 'timeout waiting for results'}), 504, {'Content-Type': 'application/json'}
        time.sleep(FINAL_WAIT_POLL_INTERVAL)
        touch_run(run_id)

    meta_path = os.path.join(run_dir, 'result_meta.json')
    report_path = os.path.join(run_dir, RAISD_REPORT)
    if not os.path.exists(report_path):
        found = _find_and_normalize_report(run_dir, prefix="results")
        if found:
            report_path = found

    def _load_meta(path: str, attempts: int = 3, delay: float = 0.05):
        m = {}
        for i in range(attempts):
            if not os.path.exists(path):
                break
            try:
                with open(path, 'r', encoding='utf-8') as mf:
                    txt = mf.read().strip()
                if not txt:
                    raise ValueError('empty meta file')
                m = json.loads(txt)
                # basic sanity: species + grid keys expected
                if 'species' in m and 'grid' in m:
                    return m
            except Exception as e:
                if i == attempts - 1:
                    append_log(run_dir, f"Meta load failed (final attempt): {e}")
                time.sleep(delay)
        return m

    meta = _load_meta(meta_path)
    if not meta:
        append_log(run_dir, 'Meta data missing or unreadable; proceeding with defaults')

    raw_match = meta.get('match') or {}

    def _get_val(src, key):
        def _norm(k):
            return str(k).replace('_', '').replace(' ', '').lower()
        if isinstance(src, dict):
            if key in src:
                return src[key]
            for k in src.keys():
                if _norm(k) == _norm(key):
                    return src[k]
            return None
        try:
            if hasattr(src, key):
                return getattr(src, key)
        except Exception:
            pass
        for attr in dir(src):
            if attr.startswith('_'):
                continue
            try:
                if _norm(attr) == _norm(key):
                    return getattr(src, attr)
            except Exception:
                continue
        return None

    expected_keys = [
        ('target', 0), ('ploidy', 2), ('mixed_ploidy', False), ('ploidy_counts', {}),
        ('demographic_model', ''), ('population', ''), ('best_jsd', 0.0)
    ]
    match = {}
    for k, default in expected_keys:
        v = _get_val(raw_match, k)
        match[k] = v if v is not None else default

    # Build/synthesize plots list
    plots = meta.get('plots') or []
    if not plots:
        # synthesize from files
        try:
            for fn in sorted(os.listdir(run_dir)):
                if fn.startswith('plot_') and fn.lower().endswith('.png'):
                    plots.append({
                        'metric': fn.replace('plot_', '').replace('.png', ''),
                        'filename': fn,
                        'url': url_for('runs_file', run_id=run_id, filename=fn)
                    })
        except Exception:
            pass

    # Selection via query param ?plot=
    desired = (request.args.get('plot') or '').strip()
    def _matches(p, needle: str) -> bool:
        if not needle:
            return False
        if needle == p.get('filename'):
            return True
        if needle == p.get('url'):
            return True
        m = p.get('metric') or ''
        if needle == m or needle == _sanitize_metric_name(m):
            return True
        return False

    sel_idx = next((i for i, p in enumerate(plots) if _matches(p, desired)), None)
    if sel_idx is not None and sel_idx != 0:
        plots.insert(0, plots.pop(sel_idx))

    plot_url = plots[0]['url'] if plots else None

    if meta.get('error') and not os.path.exists(report_path):
        return json.dumps({'error': meta.get('error')}), 500, {'Content-Type': 'application/json'}

    # Build ordered summary fields to guarantee visibility in template
    summary_fields = []
    s = str(meta.get('species') or '').replace('_', ' ').strip()
    if s:
        summary_fields.append(('Species', s))
    chrom = meta.get('chromosome')
    if chrom not in (None, ''):
        summary_fields.append(('Chromosome', chrom))
    grid = meta.get('grid')
    if grid not in (None, ''):
        summary_fields.append(('Grid', grid))
    target = match.get('target')
    if target not in (None, ''):
        summary_fields.append(('Target copies', target))
    dm = match.get('demographic_model')
    if dm:
        summary_fields.append(('Demographic model', dm))
    pop = match.get('population')
    if pop:
        summary_fields.append(('Population', pop))
    jsd_v = match.get('best_jsd')
    if isinstance(jsd_v, (int, float)):
        summary_fields.append(('JSD', f"{jsd_v:.6f}"))

    try:
        return render_template(
            'result.html',
            run_id=run_id,
            species=meta.get('species', ''),
            chromosome=meta.get('chromosome', ''),
            grid=meta.get('grid', ''),
            match=match or {},
            model_path=meta.get('model_path'),
            plot_url=plot_url,   # main plot is the first in `plots`
            report_url=meta.get('report_url') or url_for('runs_file', run_id=run_id, filename=RAISD_REPORT),
            plots=plots,
            log_url=None,
            uploaded_name=meta.get('uploaded_name') or '',
            uploaded_url=None,
            summary_fields=summary_fields,
        )
    except Exception as e:
        app.logger.exception(f"Failed to render final result for run {run_id}: {e}")
        return json.dumps({'error': f'Failed to render result: {e}'}), 500, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/heartbeat', methods=['POST'])
def runs_heartbeat(run_id):
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
    return json.dumps({'species': species, 'chromosomes': chromosomes}), 200, {'Content-Type': 'application/json'}


@app.route('/genes', methods=['GET'])
def genes_endpoint():
    """Genes overlapping window.
    Supports multi-biotype filtering via biotypes=bt1,bt2 (legacy 'biotype' kept for backward compatibility).
    Returns list of available_biotypes (unfiltered) for UI checkbox population.
    """
    species = (request.args.get('species') or '').strip()
    chromosome = (request.args.get('chromosome') or '').strip()
    start = request.args.get('start')
    end   = request.args.get('end')
    # Accept plural param first, else fallback to single
    biotypes_param = request.args.get('biotypes')
    if biotypes_param:
        biotypes = [b.strip() for b in biotypes_param.split(',') if b.strip()]
    else:
        bt_single = (request.args.get('biotype') or '').strip()
        biotypes = [bt_single] if bt_single else []

    if not species or not chromosome or start is None or end is None:
        return json.dumps({'error': 'species, chromosome, start, end are required'}), 400, {'Content-Type': 'application/json'}

    try:
        start_i = int(float(start))
        end_i   = int(float(end))
    except Exception:
        return json.dumps({'error': 'start and end must be numbers'}), 400, {'Content-Type': 'application/json'}

    try:
        # Load all hits first to derive available biotypes, then filter locally
        all_hits = _genes_in_window(species, chromosome, start_i, end_i, biotype=None)
        available_biotypes = sorted(set(str(b) for b in all_hits['biotype'] if pd.notna(b)))
        if biotypes:
            hits = all_hits[all_hits['biotype'].isin(biotypes)].copy()
        else:
            hits = all_hits
    except FileNotFoundError as e:
        return json.dumps({'error': str(e)}), 404, {'Content-Type': 'application/json'}
    except Exception as e:
        app.logger.exception("genes_endpoint error")
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

    payload = {
        'species': species,
        'chromosome': chromosome,
        'start': min(start_i, end_i),
        'end': max(start_i, end_i),
        'biotypes': biotypes,
        'available_biotypes': available_biotypes,
        'count': int(hits.shape[0]),
        'genes': hits.to_dict(orient='records'),
    }
    return json.dumps(payload), 200, {'Content-Type': 'application/json'}

@app.route('/upload', methods=['POST'])
def upload():
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
    touch_run(run_id)
    return json.dumps({'run_id': run_id, 'filename': fname}), 200, {'Content-Type': 'application/json'}

if __name__ == '__main__':
    host = os.environ.get('HOST', '0.0.0.0')
    port = int(os.environ.get('PORT', '5000'))
    debug = bool(os.environ.get('FLASK_DEBUG', '0') == '1')
    threaded = bool(os.environ.get('FLASK_THREADED', '1') == '1')
    app.run(host=host, port=port, debug=debug, threaded=threaded, use_reloader=False)
