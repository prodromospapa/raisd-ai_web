#!/usr/bin/env python3
import os
import uuid
import shlex
import shutil
import subprocess
from io import StringIO, BytesIO
from typing import Optional, Tuple, Dict
from functools import lru_cache
from datetime import datetime
import glob
import gzip
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless servers
import json
import time
import threading
from concurrent.futures import ThreadPoolExecutor, Future
import matplotlib.pyplot as plt
import pysam
import signal

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
# Retain completed result folders for this long since last activity; then auto-clean to save disk
RESULTS_RETENTION_SECONDS = int(os.environ.get('RESULTS_RETENTION_SECONDS', '7200'))  # 2h
# Short grace window before deleting a run on explicit /cleanup, to survive accidental refreshes
CLEANUP_GRACE_SECONDS = int(os.environ.get('CLEANUP_GRACE_SECONDS', '30'))

app = Flask(__name__)
app.config.update(MAX_CONTENT_LENGTH=MAX_CONTENT, SECRET_KEY=os.urandom(16))
os.makedirs(RUNS_DIR, exist_ok=True)

_run_last_seen = {}
_scheduled_cleanup = {}

def touch_run(run_id: str) -> None:
    _run_last_seen[run_id] = int(time.time())
    # Cancel any scheduled cleanup when activity is observed
    try:
        t = _scheduled_cleanup.pop(run_id, None)
        if t:
            try:
                t.cancel()
            except Exception:
                pass
    except Exception:
        pass

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
                report_present = os.path.exists(os.path.join(run_dir, RAISD_REPORT))
                # Choose timeout based on whether results exist
                timeout = RESULTS_RETENTION_SECONDS if report_present else RUN_DIR_TIMEOUT_SECONDS
                if timeout is not None and timeout >= 0:
                    if (now - last_seen) > timeout:
                        # Skip deletion if a rescan job is currently running for this run
                        try:
                            active_rescan = False
                            # Access trackers safely even if not yet defined at module import time
                            g = globals()
                            rf = g.get('_rescan_futures') or {}
                            rm = g.get('_rescan_meta') or {}
                            for jid, fut in list(rf.items()):
                                try:
                                    meta = rm.get(jid) if isinstance(rm, dict) else None
                                    if meta and meta.get('run_id') == name and fut and (not fut.done()):
                                        active_rescan = True
                                        break
                                except Exception:
                                    continue
                            if active_rescan:
                                continue
                        except Exception:
                            # Be conservative: if we can't determine, skip this cycle
                            continue
                        try:
                            shutil.rmtree(run_dir)
                            _run_last_seen.pop(name, None)
                            app.logger.info(f"Cleaned stale run dir (report={'yes' if report_present else 'no'}): {_shorten_path(run_dir)}")
                        except Exception as e:
                            app.logger.warning(f"Failed to remove stale run dir {_shorten_path(run_dir)}: {e}")
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

# GPU handling: allow starting the web app with --gpu to request GPU-accelerated
# RAiSD-AI scanning. We detect GPUs at startup and record whether the platform
# appears to have one. If the user requests --gpu but no GPU is present, we
# will inform them and not pass the -gpu argument to the RAiSD executable.
USE_GPU_FOR_RAISD = False
_GPU_CHECKED = False
_GPU_AVAILABLE = False

def _detect_gpu() -> bool:
    """Return True if a GPU device appears available on this host.

    Detection is conservative and uses a few heuristics: nvidia-smi presence
    and usable, /dev/nvidia0 device, or rocm-smi for ROCm systems. This avoids
    importing heavy GPU libs and works without additional deps.
    """
    try:
        # Fast check: nvidia-smi
        nvsmi = shutil.which('nvidia-smi')
        if nvsmi:
            proc = subprocess.run([nvsmi, '-L'], capture_output=True, text=True, timeout=3)
            if proc.returncode == 0 and proc.stdout and 'GPU' in proc.stdout:
                return True
        # Check for device node
        if os.path.exists('/dev/nvidia0'):
            return True
        # ROCm / AMD: rocm-smi may exist
        rmsmi = shutil.which('rocm-smi')
        if rmsmi:
            proc = subprocess.run([rmsmi, '--showproductname'], capture_output=True, text=True, timeout=3)
            if proc.returncode == 0 and proc.stdout and proc.stdout.strip():
                return True
        # As a last resort, look for obvious PCI vendor strings (best-effort, may be slow)
        lspci = shutil.which('lspci')
        if lspci:
            proc = subprocess.run([lspci], capture_output=True, text=True, timeout=3)
            if proc.returncode == 0 and proc.stdout:
                out = proc.stdout.lower()
                if 'nvidia' in out or 'amd' in out or 'advanced micro devices' in out:
                    return True
    except Exception:
        pass
    return False

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

def _any_report_exists(run_dir: str) -> bool:
    try:
        if os.path.exists(os.path.join(run_dir, RAISD_REPORT)):
            return True
        pats = [
            os.path.join(run_dir, 'RAiSD_Report.*'),
            os.path.join(run_dir, '**', 'RAiSD_Report.*'),
        ]
        for pat in pats:
            if glob.glob(pat, recursive=True):
                return True
    except Exception:
        pass
    return False

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
        # Sanitize any absolute filesystem paths in messages before writing logs
        safe_msg = _sanitize_log_message(msg)
        with open(os.path.join(run_dir, "run.log"), "a", encoding="utf-8") as f:
            f.write(f"[{ts}] {safe_msg}\n")
    except Exception:
        pass


def _shorten_path(p: str) -> str:
    """Return a short/redacted representation of a filesystem path.

    Rules:
    - If path is empty or not an absolute path, return as-is (or basename for files).
    - For absolute paths, return '<REDACTED>/(basename or last two components)'.
    """
    try:
        if not p:
            return p
        # If p looks like multiple paths joined (e.g. CMD), don't attempt to split; handle single path
        if '\n' in p or '\t' in p:
            return p
        if os.path.isabs(p):
            # Keep last two path components for context where possible
            parts = p.rstrip(os.path.sep).split(os.path.sep)
            if len(parts) >= 3:
                return os.path.join('<REDACTED>', parts[-2], parts[-1])
            return os.path.join('<REDACTED>', parts[-1])
        # Not absolute: return basename for safety
        return os.path.basename(p) or p
    except Exception:
        return '<REDACTED_PATH>'


def _sanitize_log_message(msg: str) -> str:
    """Sanitize a free-form log message by replacing absolute paths with shortened forms.

    This is a conservative sanitizer: it replaces sequences that look like absolute
    Unix paths (starting with '/') with a shortened representation.
    """
    try:
        # Replace any /... sequences with shortened version while preserving small tokens like '/dev/nvidia0'
        def repl(match):
            p = match.group(0)
            return _shorten_path(p)

        # A simple regex to match absolute Unix paths (space/newline/beginning-delimited)
        return re.sub(r'(?<!\S)/[A-Za-z0-9_./\\-]+', repl, str(msg))
    except Exception:
        return '<REDACTED_MSG>'

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
    # Accept a few common species folder name variants so callers may pass
    # either 'Homo sapiens' or 'Homo_sapiens'. Try the provided species
    # value first, then fall back to underscore/space variants.
    candidates = []
    s = str(species or '')
    candidates.append(s)
    # replace spaces with underscores and vice-versa
    if ' ' in s:
        candidates.append(s.replace(' ', '_'))
    if '_' in s:
        candidates.append(s.replace('_', ' '))
    # also try a collapsed-underscore form (remove duplicate underscores)
    candidates.append(s.replace('  ', ' ').replace('__', '_'))
    # normalize case-preserving duplicates removed
    seen = set()
    tried_paths = []
    for cand in candidates:
        if not cand or cand in seen:
            continue
        seen.add(cand)
        path = ANNOTATION_DIR_TMPL.format(species=cand, chr=str(chromosome))
        tried_paths.append(path)
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, sep="\t", dtype={"chromosome": str})
                app.logger.info(f"Loaded annotation file: {_shorten_path(path)}")
            except Exception as read_e:
                app.logger.exception(f"Failed to read annotation file {_shorten_path(path)}: {read_e}")
                raise
            break
    else:
        tried_short = [ _shorten_path(p) for p in tried_paths ]
        app.logger.warning(f"_load_annotation: no annotation file found for species variants {candidates!r} chromosome={chromosome}. Tried: {tried_short}")
        raise FileNotFoundError(f"Annotation file not found. Tried: {', '.join(tried_paths)}")

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
    # Sort so that the gene with the largest overlap is first (descending), then by start for stability.
    # This supports UI requirement: "highest is in the first row" meaning highest overlap shown at top.
    try:
        hits = hits.sort_values(["overlap_bp", "start", "end", "label"], ascending=[False, True, True, True])
    except Exception:
        # Fallback to original ordering if sort fails unexpectedly
        hits = hits.sort_values(["start", "end", "label"])
    hits = hits.reset_index(drop=True)
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

def normalize_sfs_df(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame containing only the polymorphic bins by dropping
    the smallest and largest integer-labeled columns (commonly 0 and n).

    The function accepts DataFrames whose column labels are int-like strings
    or ints. It will preserve the input index and return floating probabilities
    normalized per-row across the remaining bins.
    """
    # Map to integer column labels
    cols_int = sorted(int(c) for c in df.columns)
    if len(cols_int) <= 2:
        # Nothing meaningful remains after dropping first/last
        raise ValueError("SFS table must have more than two bins to normalize (need polymorphic bins)")
    desired_ints = cols_int[1:-1]
    desired_cols = [str(i) for i in desired_ints]
    # Ensure columns exist; missing ones filled with zeros
    out = pd.DataFrame(index=df.index, columns=desired_cols, dtype=float)
    for c in desired_cols:
        if c in df.columns:
            out[c] = df[c].astype(float)
        else:
            out[c] = 0.0
    # Renormalize rows to sum to 1 (polymorphic-only probabilities)
    row_sums = out.sum(axis=1).replace(0, 1.0)
    out = out.div(row_sums, axis=0)
    return out

def project_expected_df(df, target):
    # Determine integer column labels and use the maximum label as reference N.
    cols = sorted(int(c) for c in df.columns)
    if len(cols) == 0:
        raise ValueError("DataFrame has no columns")
    min_col, max_col = cols[0], cols[-1]
    # Use max label as N. This supports inputs that omit the 0 bin (columns start at 1).
    N = max_col
    if target > N:
        raise ValueError(f"target ({target}) > reference N ({N}); projection is downsampling only.")
    G, cache = _pre_lgamma(N), {}
    rows = []
    for _, row in df.iterrows():
        # Build full-length p vector indexed 0..N
        p_full = np.zeros(N + 1, float)
        for c, v in row.items():
            idx = int(c)
            if 0 <= idx <= N:
                p_full[idx] = float(v)
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
# RAiSD native output order (do NOT change this list ordering; it must match tool output):
RAISD_HEADER_RAW = [
    r"$\mu_{Var}$",
    r"$\mu_{SFS}$",
    r"$\mu_{LD}$",
    r"$\mu$",
    r"$sweep_{TR}$",
    r"$\mu_{var}^{sweep_{TR}}$",
]

# Desired display / report order requested: muVar, muSFS, muLD, mu, mu_var_sweep_TR, sweepTR
# i.e. swap the last two compared to raw tool output.
DISPLAY_HEADER = [
    r"$\mu_{Var}$",
    r"$\mu_{SFS}$",
    r"$\mu_{LD}$",
    r"$\mu$",
    r"$\mu_{var}^{sweep_{TR}}$",
    r"$sweep_{TR}$",
]

# UI requirement: only expose these metrics in the dropdown (plots for others may still be generated internally)
# Keep internal/raw RAiSD header strings here; provide a small display map so the UI
# shows nicer names (μ and sweepTR) while the code continues to use the original
# RAiSD column names for data lookups.
EXPOSED_METRICS = {r"$\mu$", r"$sweep_{TR}$"}
EXPOSED_METRIC_DISPLAY = {
    r"$\mu$": "μ",
    r"$sweep_{TR}$": "sweepTR",
}

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
        if df.shape[1] > len(RAISD_HEADER_RAW):
            df = df.iloc[:, -len(RAISD_HEADER_RAW):]
        elif df.shape[1] < len(RAISD_HEADER_RAW):
            continue
        # First assign raw header names matching the native RAiSD order
        df.columns = RAISD_HEADER_RAW
        # Reorder to desired display order without altering values-to-metric mapping
        df = df[DISPLAY_HEADER]
        dfs.append(df)
    if not dfs:
        raise ValueError("No parsable data blocks found in RAiSD_Report.results")
    return pd.concat(dfs)

def _sanitize_metric_name(m: str) -> str:
    import re
    s = re.sub(r'[^0-9A-Za-z]+', '_', m)
    s = re.sub(r'_+', '_', s).strip('_')
    return s or 'metric'

def plot_metric_bytes(df: pd.DataFrame, metric: str, xmin: Optional[float] = None, xmax: Optional[float] = None) -> bytes:
    """Render a single-metric plot to PNG bytes; optional x-range can be supplied."""
    plt.figure(figsize=(10, 6))
    try:
        display_metric = EXPOSED_METRIC_DISPLAY.get(metric, metric)
        if metric not in df.columns:
            plt.text(0.5, 0.5, f'No data for {display_metric}', ha='center', va='center')
        else:
            plt.plot(df.index, df[metric], label=display_metric)
            plt.xlabel('Position')
            plt.ylabel(display_metric)
            plt.title(display_metric)
            plt.legend()
            plt.grid(True)
            if xmin is not None or xmax is not None:
                lo = xmin if xmin is not None else df.index.min()
                hi = xmax if xmax is not None else df.index.max()
                plt.xlim(lo, hi)
            # Default y-axis range for sweep_TR metric should be 0..1
            try:
                # Compare sanitized names so either raw LaTeX header or sanitized forms match
                if _sanitize_metric_name(metric) == _sanitize_metric_name(r"$sweep_{TR}$"):
                    # Make it a bit taller than [0, 1] so values at 1.0 are clearly visible
                    plt.ylim(0.0, 1.02)
            except Exception:
                # If anything goes wrong, don't prevent plotting
                pass
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
    # Mark activity so cleanup doesn't remove an actively viewed run
    touch_run(run_id)

    metric = request.args.get('metric') or ''
    xmin = request.args.get('xmin')
    xmax = request.args.get('xmax')

    # Prefer persisted per-(dm,pop) report to avoid using the generic RAiSD_Report.results
    report_path = _choose_persisted_report(run_dir) or os.path.join(run_dir, RAISD_REPORT)
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
    touch_run(run_id)
    # Prefer persisted per-(dm,pop) report for CSV export too
    report_path = _choose_persisted_report(run_dir) or os.path.join(run_dir, RAISD_REPORT)
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
    touch_run(run_id)
    # Prefer persisted per-(dm,pop) report so CSV reflects the selected model/population
    report_path = _choose_persisted_report(run_dir) or os.path.join(run_dir, RAISD_REPORT)
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
    # Name the CSV after the underlying report when available (e.g., RAiSD_Report.<dm>__<pop>.csv)
    try:
        base = os.path.basename(report_path)
        if base.endswith('.results'):
            csv_name = base[:-8] + '.csv'
        else:
            csv_name = 'RAiSD_Report.csv'
    except Exception:
        csv_name = 'RAiSD_Report.csv'
    return Response(csv_data, mimetype='text/csv', headers={'Content-Disposition': f'attachment; filename={csv_name}'})

# Static PNG plots removed; readiness hinges solely on report presence
def _plots_exist(run_dir: str) -> bool:  # backward compatibility shim
    return True

# ------------------------------ Pipeline ------------------------------
def compute_best_match(species: str, vcf_path: str):
    sfs_table = os.path.join(DATA_DIR, species, 'sfs.csv')
    if not os.path.exists(sfs_table):
        raise FileNotFoundError(f"Missing SFS table: {sfs_table}")
    df = load_sfs_matrix(sfs_table)
    # Normalize/convert table to polymorphic-only bins (drop 0 and n)
    try:
        df = normalize_sfs_df(df)
    except Exception:
        # If normalization fails, fall back to the raw df but continue
        pass
    target, ploidy, mixed, counts = infer_samples_and_ploidy(vcf_path)
    # Attempt to project the expected SFS to the observed target. If the
    # observed sample (target) is larger than the reference N in the
    # persisted SFS table, instead project the observed SFS down to the
    # reference N and compare in that reduced space. This allows users to
    # upload larger input files without forcing `sfs.csv` to be re-generated.
    try:
        proj = project_expected_df(df, target)
        sfs = sfs_from_bcftools(vcf_path, target)
        q = norm(sfs[1:])  # observed SFS (skip invariant)
    except ValueError as e:
        # Detect the specific 'target > reference N' error from project_expected_df
        msg = str(e)
        if 'target' in msg and 'reference N' in msg:
            # Reference N (max column) from the SFS table
            try:
                N_ref = max(int(c) for c in df.columns)
            except Exception:
                # Re-raise original if we can't determine a reference
                raise
            # Project expected SFS rows to the reference N so both sides use the
            # same dimensionality for comparison.
            proj = project_expected_df(df, int(N_ref))
            # Read observed SFS (full length) and project it down to N_ref.
            sfs = sfs_from_bcftools(vcf_path, target)
            obs_full = np.asarray(sfs, float)
            total = obs_full.sum()
            if total > 0:
                obs_full = obs_full / total
            # Build a one-row DataFrame with columns '0'..'<M>' so we can reuse
            # project_expected_df machinery to downsample the observed SFS.
            cols = [str(i) for i in range(0, len(obs_full))]
            obs_df = pd.DataFrame([obs_full], index=['_obs_'], columns=cols)
            obs_proj = project_expected_df(obs_df, int(N_ref)).loc['_obs_'].values
            q = norm(obs_proj)
        else:
            # Unknown ValueError - propagate
            raise
    pop_jsd_pairs = []  # list of (full_name, jsd_val)
    for pop in proj.index:
        try:
            d = jsd(proj.loc[pop].values, q)
        except Exception:
            d = float('inf')
        pop_jsd_pairs.append((pop, d))
    pop_jsd_pairs.sort(key=lambda x: x[1])
    best_pop, best_jsd = pop_jsd_pairs[0]
    demographic_model, population = best_pop.split("=")
    best_expected = proj.loc[best_pop].values.tolist()
    # Build structured lists
    top_matches = []
    for name, val in pop_jsd_pairs[:6]:
        try:
            dm, popn = name.split('=')
        except ValueError:
            dm, popn = name, ''
        # Attach the projected expected SFS for the same target so UI can draw without re-fetch
        try:
            expected_vec = proj.loc[name].values.tolist()
        except Exception:
            expected_vec = None
        item = {'demographic_model': dm, 'population': popn, 'jsd': float(val)}
        if expected_vec is not None:
            item['expected_sfs'] = expected_vec
        top_matches.append(item)
    # Full list (may be large). We omit expected SFS for non-top matches to keep size smaller.
    all_jsd = []
    for name, val in pop_jsd_pairs:
        try:
            dm, popn = name.split('=')
        except ValueError:
            dm, popn = name, ''
        item = {'demographic_model': dm, 'population': popn, 'jsd': float(val)}
        # Attach expected_sfs if available so UI can render picked traces without extra projection
        try:
            item_expected = proj.loc[name].values.tolist()
            item['expected_sfs'] = item_expected
        except Exception:
            pass
        all_jsd.append(item)
    return {
        'target': target,
        'ploidy': ploidy,
        'mixed_ploidy': mixed,
        'ploidy_counts': counts,
        'demographic_model': demographic_model,
        'population': population,
        'best_jsd': best_jsd,
        'input_sfs': q.tolist(),
        'best_expected_sfs': best_expected,
        'top_matches': top_matches,
        'all_jsd': all_jsd,
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


def _choose_persisted_report(run_dir: str) -> Optional[str]:
    """Return the best persisted per-(dm,pop) report for a run.
    Preference order:
      1. `result_meta.json` -> last_rescan.report if exists
      2. any files matching `RAiSD_Report.<dm>__<pop>.results` (newest)
      3. if only RAiSD_Report.results exists and meta.match present: copy it to a persisted name and return that
      4. otherwise None
    """
    meta_path = os.path.join(run_dir, 'result_meta.json')
    meta = {}
    try:
        if os.path.exists(meta_path):
            with open(meta_path, 'r', encoding='utf-8') as mf:
                meta = json.load(mf)
    except Exception:
        meta = {}

    # 1. last_rescan.report
    last = meta.get('last_rescan') or {}
    rpt = last.get('report')
    if rpt:
        p = os.path.join(run_dir, rpt)
        if os.path.exists(p):
            return p

    # 2. look for named reports produced by RAiSD (-n <name> -> RAiSD_Report.<name>)
    # Prefer the newest matching the "<dm>__<pop>" or "<dm>" patterns; fall back to any RAiSD_Report.*
    pattern = os.path.join(run_dir, 'RAiSD_Report.*')
    candidates = []
    for path in glob.glob(pattern):
        try:
            candidates.append((os.path.getmtime(path), path))
        except Exception:
            continue
    if candidates:
        candidates.sort(reverse=True)
        return candidates[0][1]

    # 3. fallback: generic report if present
    generic = os.path.join(run_dir, RAISD_REPORT)
    if os.path.exists(generic):
        return generic
    return None

def run_raisd(run_dir: str, model_path: str, vcf_path: str, chromosome: int, grid=300, ra_exe: Optional[str] = None):
    grid_dir = os.path.join(run_dir, RAISD_GRIDDIR)
    if os.path.isdir(grid_dir):
        shutil.rmtree(grid_dir, ignore_errors=True)
    ra_bin = ra_exe or choose_raisd_exe(vcf_path)
    # Derive demographic model and population from the model_path (expected layout: .../<species>/<dm>/<pop>/<chrom>/RAiSD_Model.model)
    try:
        parts = os.path.normpath(model_path).split(os.sep)
        # last is RAiSD_Model.model, parent is chrom, grandparent is pop, great-grandparent is dm
        if len(parts) >= 4:
            dm = parts[-4]
            pop = parts[-3]
        else:
            dm = 'results'
            pop = ''
    except Exception:
        dm = 'results'
        pop = ''

    def _sanitize_name(s: str) -> str:
        return re.sub(r'[^A-Za-z0-9_\-]+', '_', (s or '').strip()).strip('_')

    # RAiSD -n should be the base name WITHOUT the RAiSD_Report. prefix.
    # We'll pass <dm>__<pop> (or <dm>) so RAiSD creates RAiSD_Report.<dm>__<pop>
    # and RAiSD_Info.<dm>__<pop>.
    name_base = 'results'
    try:
        if dm:
            if pop:
                name_base = f"{_sanitize_name(dm)}__{_sanitize_name(pop)}"
            else:
                name_base = f"{_sanitize_name(dm)}"
    except Exception:
        name_base = 'results'

    cmd = [
        ra_bin, '-n', name_base,
        '-mdl', model_path,
        '-f',
        '-op', 'SWP-SCN',
        '-I', vcf_path,
        '-frm',
        '-G', str(grid),
        '-pci', '1', '1',
        '-R'
    ]
    # Append -gpu if requested and a GPU is available. The USE_GPU_FOR_RAISD
    # global is set at process startup when the user passed --gpu and a GPU
    # was detected. If --gpu was requested but no GPU was found, USE_GPU_FOR_RAISD
    # will be False and we will not add the flag.
    try:
        if USE_GPU_FOR_RAISD:
            cmd.append('-gpu')
    except Exception:
        pass
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

    # Determine the expected RAiSD-produced report name and copy to generic filename for compatibility
    expected_report_named = os.path.join(run_dir, f"RAiSD_Report.{name_base}")
    report_path = os.path.join(run_dir, RAISD_REPORT)
    if os.path.exists(expected_report_named):
        try:
            if os.path.abspath(expected_report_named) != os.path.abspath(report_path):
                shutil.copyfile(expected_report_named, report_path)
                append_log(run_dir, f"Linked primary report: {os.path.basename(expected_report_named)} -> {RAISD_REPORT}")
        except Exception as ce:
            append_log(run_dir, f"Warning: failed to link named report: {ce}")
    elif not os.path.exists(report_path):
        # Fallback to previous search-based normalization
        found = _find_and_normalize_report(run_dir, prefix="results")
        if found:
            append_log(run_dir, f"Normalized report to: {os.path.basename(found)}")
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
            f"Searched under: {_shorten_path(run_dir)}\n"
            f"CMD: {' '.join(cmd)}\n"
            "Possible causes: incompatible flags, empty scan region, or a build that writes elsewhere.\n"
            "stdout (tail):\n" + "\n".join(std_tail) + "\n\n"
            "stderr (tail):\n" + "\n".join(err_tail)
        )
    # Successful run: remove grid directory and RAiSD_Info.<name> to save space
    try:
        if os.path.isdir(grid_dir):
            shutil.rmtree(grid_dir, ignore_errors=True)
            append_log(run_dir, f"Removed grid directory: {RAISD_GRIDDIR}")
        # Delete RAiSD_Info.<name> (not needed)
        info_named = os.path.join(run_dir, f"RAiSD_Info.{name_base}")
        if os.path.exists(info_named):
            try:
                os.remove(info_named)
                append_log(run_dir, f"Removed RAiSD info file: {os.path.basename(info_named)}")
            except Exception:
                pass
    except Exception as ce:
        append_log(run_dir, f"Warning: failed to remove grid directory {RAISD_GRIDDIR}: {ce}")
    return report_path

# ------------------------------- Routes -------------------------------
# In-memory tracker for rescan jobs (per-process). Also persist basic status to disk for UI polling.
_rescan_futures: Dict[str, Future] = {}
_rescan_meta: Dict[str, dict] = {}

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
    user_chrom_input = chromosome  # keep original text
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
            # If this is a staged upload (no file in this POST), require explicit confirmation
            # to actually start analysis. Clients should pass start_now=1 to trigger processing.
            start_now = request.form.get('start_now') == '1'
            if not start_now:
                # Touch the run dir so it isn't immediately cleaned up, and record staged metadata.
                touch_run(run_id)
                try:
                    staged_meta = {
                        'run_id': run_id,
                        'uploaded_name': uploaded_filename,
                        'species': request.form.get('species', '').strip(),
                        'chromosome': request.form.get('chromosome', '').strip(),
                        'grid': request.form.get('grid', '').strip(),
                        'status': 'staged'
                    }
                    with open(os.path.join(run_dir, 'staged_upload.json'), 'w', encoding='utf-8') as sf:
                        sf.write(json.dumps(staged_meta))
                except Exception:
                    pass
                # Return a lightweight JSON response indicating the upload is staged and won't start yet.
                return json.dumps({'run_id': run_id, 'status': 'staged'}), 200, {'Content-Type': 'application/json'}
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
        # Ensure minimal meta exists early so the UI can read species/chromosome
        try:
            meta_path = os.path.join(run_dir, 'result_meta.json')
            meta_early = {
                'run_id': run_id,
                'species': species,
                'chromosome': chromosome,
                'grid': grid_str if isinstance(grid_str, str) else str(grid_str),
                'uploaded_name': filename,
                'plots_ready': False,
                'plots': []
            }
            with open(meta_path, 'w', encoding='utf-8') as mf:
                mf.write(json.dumps(meta_early))
        except Exception:
            pass
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
            # Write minimal meta immediately so page loads have metadata available
            try:
                meta_path = os.path.join(run_dir, 'result_meta.json')
                meta_early = {
                    'run_id': run_id,
                    'species': species,
                    'chromosome': chromosome,
                    'grid': grid_str if isinstance(grid_str, str) else str(grid_str),
                    'uploaded_name': filename,
                    'plots_ready': False,
                    'plots': []
                }
                with open(meta_path, 'w', encoding='utf-8') as mf:
                    mf.write(json.dumps(meta_early))
            except Exception:
                pass
        else:
            flash("No file available for analysis.")
            return redirect(url_for('index'))

    # --------------------------------------------------------------
    # Chromosome inference / validation from VCF header
    # If the uploaded VCF header lists contigs and the user-provided chromosome
    # (integer) is not directly represented, attempt smart matching:
    #  - Accept exact numeric match
    #  - Accept 'chrN' vs 'N' differences
    #  - If single contig present, override chromosome with that contig
    #  - If cannot reconcile, show an error listing available contigs
    # --------------------------------------------------------------
    try:
        contigs = []
        try:
            vf_tmp = pysam.VariantFile(vcf_path)
            contigs = [str(k) for k in vf_tmp.header.contigs.keys()]
        except Exception:
            contigs = []
        if contigs:
            contig_set = set(contigs)
            numeric_str = str(chromosome)
            chr_pref = f"chr{numeric_str}"
            chosen_contig = None
            if numeric_str in contig_set:
                chosen_contig = numeric_str
            elif chr_pref in contig_set:
                chosen_contig = chr_pref
            else:
                # Attempt to find a single distinct contig if only one present
                if len(contigs) == 1:
                    chosen_contig = contigs[0]
                    append_log(run_dir, f"Auto-overriding chromosome {chromosome} -> {chosen_contig} (single contig in VCF)")
                else:
                    # Try to strip 'chr' prefixes and compare
                    stripped_map = {c[3:]: c for c in contigs if isinstance(c, str) and c.lower().startswith('chr') and len(c) > 3}
                    if numeric_str in stripped_map:
                        chosen_contig = stripped_map[numeric_str]
            if chosen_contig is None:
                # Could not reconcile; present helpful error if mismatch likely
                # (Only if user selection not among available options and not single contig case already handled)
                if str(chromosome) not in contig_set and f"chr{chromosome}" not in contig_set:
                    shown = [str(x) for x in contigs[:20]]
                    flash(f"Selected chromosome {chromosome} not found in VCF. Available: {', '.join(shown)}{'...' if len(contigs)>20 else ''}")
                    # Clean up created run_dir so user can retry without clutter
                    try:
                        append_log(run_dir, "Removing run directory due to chromosome mismatch")
                        shutil.rmtree(run_dir, ignore_errors=True)
                    except Exception:
                        pass
                    return redirect(url_for('index'))
            else:
                # Normalize chromosome directory name (drop leading 'chr' if present) for model path usage later
                chosen_contig_str = str(chosen_contig)
                model_chr = re.sub(r'^chr', '', chosen_contig_str, flags=re.IGNORECASE)
                if model_chr.isdigit():
                    chromosome = int(model_chr)
                append_log(run_dir, f"Chromosome validated/inferred: user={user_chrom_input} -> using_contig={chosen_contig} -> model_dir_component={chromosome}")
    except Exception as ce:
        append_log(run_dir, f"Chromosome inference warning: {ce}")

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
            'analysis_input': None,  # basename of VCF/BCF actually used for analysis (post-conversion)
        }
        try:
            # If input is a BCF, convert to compressed VCF (.vcf.gz) so the rest of the
            # pipeline uniformly works with VCF.GZ (e.g. RAiSD-AI-ZLIB selection, indexing).
            # This also future-proofs any code paths that assume text VCF when parsing.
            converted_path = None
            if vcf_path.lower().endswith('.bcf'):
                base = os.path.splitext(vcf_path)[0]  # drop .bcf
                converted_path = base + '.vcf.gz'
                append_log(run_dir, f"Converting BCF to VCF.GZ: {os.path.basename(vcf_path)} -> {os.path.basename(converted_path)}")
                conv_cmd = [BCFTOOLS, 'view', '-Oz', '-o', converted_path, vcf_path]
                proc_conv = subprocess.run(conv_cmd, capture_output=True, text=True, env=_bcftools_env())
                if proc_conv.returncode != 0 or not os.path.exists(converted_path):
                    raise RuntimeError(
                        "Failed to convert BCF to VCF.GZ.\n" +
                        f"CMD: {' '.join(conv_cmd)}\nstdout:\n{proc_conv.stdout}\nstderr:\n{proc_conv.stderr}"
                    )
                # Prefer working on converted file from here on.
                vcf_path = converted_path
                append_log(run_dir, "BCF conversion successful; proceeding with .vcf.gz")
            else:
                append_log(run_dir, "Input not BCF; no conversion needed")
            meta['analysis_input'] = os.path.basename(vcf_path)
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
            # Persist detailed match/JSD data for later hover usage
            try:
                details = {
                    'observed_sfs': match.get('input_sfs'),
                    'best_expected_sfs': match.get('best_expected_sfs'),
                    'top_matches': match.get('top_matches'),  # may include expected_sfs per top item
                    'all_jsd': match.get('all_jsd'),  # now may include expected_sfs for all entries
                }
                with open(os.path.join(run_dir, 'match_details.json'), 'w', encoding='utf-8') as md:
                    md.write(json.dumps(details))
            except Exception as det_e:
                append_log(run_dir, f"Warning: failed to write match_details.json: {det_e}")
            model_path = os.path.join(
                DATA_DIR, species, match['demographic_model'], match['population'], str(chromosome), 'RAiSD_Model.model'
            )
            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Model file missing: {model_path}")
            meta['model_path'] = model_path

            report_path = run_raisd(run_dir, model_path, vcf_path, chromosome, grid=grid, ra_exe=ra_exe)
            df = read_raisd_report(report_path)
            try:
                _cached_report.cache_clear()  # type: ignore[attr-defined]
            except Exception:
                pass

            # Do NOT generate static plot PNGs; we only serve dynamic data via /metric_data
            # Build plots list as dicts including a human-friendly display label for the UI
            plots = []
            for h in DISPLAY_HEADER:
                if h in EXPOSED_METRICS:
                    display = EXPOSED_METRIC_DISPLAY.get(h, h)
                    plots.append({'metric': h, 'display': display})

            # URLs for template
            try:
                with app.app_context():
                    meta['report_url'] = url_for('runs_file', run_id=meta['run_id'], filename=RAISD_REPORT)
                    # Build metric list for UI (dynamic plotting only)
                    preference_order = [r"$\mu$", r"$sweep_{TR}$"]
                    # Persist plot dicts (metric + display) for the template
                    meta['plots'] = [{'metric': p['metric'], 'display': p.get('display', p['metric'])} for p in plots]
                    meta['plots'].sort(key=lambda pl: (preference_order.index(pl['metric']) if pl['metric'] in preference_order else 999, pl['metric']))
                    meta['plot_url'] = None  # no static image
            except Exception:
                meta['report_url'] = None

            # mark plots ready (logical readiness of dynamic metrics)
            meta['plots_ready'] = True
            append_log(run_dir, f"Parsed report (no static plots created): {RAISD_REPORT}")

            # No extra persisted copies; rely on RAiSD -n naming and generic aliasing handled in run_raisd

            # clean big temporaries
            try:
                # Remove ONLY the original uploaded file if it differs from the analysis input.
                # Keep the actual analysis_input (and its index) so we can recompute match details if needed later.
                analysis_input = meta.get('analysis_input')
                input_path = os.path.join(run_dir, filename)
                if os.path.basename(input_path) != analysis_input:
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
    raisd_stdout = os.path.join(run_dir, 'raisd_stdout.log')
    raisd_stderr = os.path.join(run_dir, 'raisd_stderr.log')
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
    done = _any_report_exists(run_dir)
    # also include RAiSD stdout/stderr if present (trimmed)
    rs_out = ''
    rs_err = ''
    try:
        if os.path.exists(raisd_stdout):
            with open(raisd_stdout, 'r', encoding='utf-8', errors='replace') as rf:
                rs_out = rf.read()[-5000:]
    except Exception:
        rs_out = '(failed to read RAiSD stdout)'
    try:
        if os.path.exists(raisd_stderr):
            with open(raisd_stderr, 'r', encoding='utf-8', errors='replace') as rf:
                rs_err = rf.read()[-5000:]
    except Exception:
        rs_err = '(failed to read RAiSD stderr)'

    return json.dumps({'lines': lines, 'raisd_stdout': rs_out, 'raisd_stderr': rs_err, 'done': bool(done)}), 200, {'Content-Type': 'application/json'}


@app.route('/runs/<run_id>/cancel', methods=['POST'])
def runs_cancel(run_id):
    """Best-effort cancellation: try to cancel the submitted future, and kill any
    processes whose current working directory is the run directory. Returns JSON
    describing what was attempted."""
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}

    # Try to cancel the future if present
    fut = None
    with _futures_lock:
        fut = _run_futures.get(run_id)
    cancelled = False
    try:
        if fut:
            cancelled = fut.cancel()
    except Exception:
        cancelled = False

    # Best-effort: kill processes whose cwd is inside the run dir
    killed = []
    try:
        for pid_str in os.listdir('/proc'):
            if not pid_str.isdigit():
                continue
            pid = int(pid_str)
            try:
                # readlink may raise if process exits
                proc_cwd = os.readlink(f'/proc/{pid}/cwd')
                # match exact run_dir or a subpath
                try:
                    # normalization
                    if os.path.commonpath([proc_cwd, run_dir]) == run_dir:
                        try:
                            os.kill(pid, signal.SIGTERM)
                            killed.append(pid)
                        except Exception:
                            pass
                except Exception:
                    # fallback to simple prefix check
                    if proc_cwd.startswith(run_dir):
                        try:
                            os.kill(pid, signal.SIGTERM)
                            killed.append(pid)
                        except Exception:
                            pass
            except Exception:
                continue
        # give processes a moment, then force kill remaining
        time.sleep(0.5)
        for pid in list(killed):
            try:
                os.kill(pid, 0)
                try:
                    os.kill(pid, signal.SIGKILL)
                except Exception:
                    pass
            except Exception:
                # process already gone
                try:
                    killed.remove(pid)
                except Exception:
                    pass
    except Exception:
        # ignore errors in best-effort cleanup
        pass

    append_log(run_dir, f"Cancellation requested: fut_cancelled={cancelled}, pids_killed={killed}")
    return json.dumps({'cancelled': cancelled, 'pids_killed': killed}), 200, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/rescan', methods=['POST'])
def runs_rescan(run_id):
    """
    Start a background rescan using a specific demographic_model/population for the same run
    and same chromosome/species/input file. Returns a job id to poll for completion.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    try:
        payload = request.get_json(silent=True) or {}
        dm = (payload.get('demographic_model') or '').strip()
        pop = (payload.get('population') or '').strip()
        # optional grid override from client
        grid_override = payload.get('grid')
        if not dm or not pop:
            return json.dumps({'error': 'demographic_model and population are required'}), 400, {'Content-Type': 'application/json'}
    except Exception as e:
        return json.dumps({'error': f'invalid json: {e}'}), 400, {'Content-Type': 'application/json'}

    # Prevent starting another rescan if one is in-flight for this run
    for jid, fut in list(_rescan_futures.items()):
        meta = _rescan_meta.get(jid) or {}
        if meta.get('run_id') == run_id and not fut.done():
            return json.dumps({'error': 'a rescan is already running', 'job_id': jid}), 409, {'Content-Type': 'application/json'}

    meta_path = os.path.join(run_dir, 'result_meta.json')
    if not os.path.exists(meta_path):
        return json.dumps({'error': 'run metadata missing'}), 500, {'Content-Type': 'application/json'}
    try:
        with open(meta_path, 'r', encoding='utf-8') as mf:
            meta = json.load(mf)
    except Exception as e:
        return json.dumps({'error': f'failed to read meta: {e}'}), 500, {'Content-Type': 'application/json'}
    species = meta.get('species')
    chromosome = meta.get('chromosome')
    grid = meta.get('grid') or 300
    # If client supplied a grid override, validate and use it
    try:
        if grid_override is not None and str(grid_override).strip() != '':
            grid = int(grid_override)
            if grid < 1:
                raise ValueError('grid must be >= 1')
    except Exception:
        # ignore invalid override and keep meta/grid default
        grid = int(grid) if grid else 300
    analysis_input = meta.get('analysis_input')
    if analysis_input:
        vcf_path = os.path.join(run_dir, analysis_input)
    else:
        # Fallback to uploaded_name if analysis_input missing
        vcf_path = os.path.join(run_dir, meta.get('uploaded_name') or '')
    if not (species and isinstance(chromosome, int) and os.path.exists(vcf_path)):
        return json.dumps({'error': 'incomplete run metadata or missing input'}), 400, {'Content-Type': 'application/json'}
    # Normalize grid to int
    try:
        grid = int(grid)
    except Exception:
        grid = 300

    model_path = os.path.join(DATA_DIR, species, dm, pop, str(chromosome), 'RAiSD_Model.model')
    if not os.path.exists(model_path):
        return json.dumps({'error': f'model not found for {dm}/{pop} chr{chromosome}'}), 404, {'Content-Type': 'application/json'}

    def _sanitize(s: str) -> str:
        return re.sub(r'[^A-Za-z0-9_\-]+', '_', s or '').strip('_')

    # Expected RAiSD-produced name components from -n
    base_name = f"{_sanitize(dm)}__{_sanitize(pop)}"
    persisted_name = f"RAiSD_Report.{base_name}"
    persisted_path = os.path.join(run_dir, persisted_name)

    job_id = uuid.uuid4().hex[:12]

    # If a persisted report already exists for this model/pop, reuse it immediately
    if os.path.exists(persisted_path):
        append_log(run_dir, f"Reusing existing report for {dm}/{pop}: {persisted_name}")
        # update meta and write status file
        try:
            try:
                with open(meta_path, 'r', encoding='utf-8') as mf:
                    cur = json.load(mf)
            except Exception:
                cur = {}
            cur['model_path'] = model_path
            cur['plots_ready'] = True
            cur['last_rescan'] = {
                'job_id': job_id,
                'demographic_model': dm,
                'population': pop,
                'grid': grid,
                'timestamp': time.time(),
                'report': persisted_name,
                'reused': True
            }
            with open(meta_path, 'w', encoding='utf-8') as mf:
                mf.write(json.dumps(cur))
        except Exception:
            pass
        status_path = os.path.join(run_dir, f'rescan_{job_id}.json')
        try:
            with open(status_path, 'w', encoding='utf-8') as sf:
                sf.write(json.dumps({
                    'job_id': job_id,
                    'status': 'done',
                    'reused': True,
                    'demographic_model': dm,
                    'population': pop,
                    'grid': grid,
                    'report': persisted_name
                }))
        except Exception:
            pass
        # create a resolved Future so polling endpoints see it as done
        fut = Future()
        fut.set_result(None)
        _rescan_futures[job_id] = fut
        _rescan_meta[job_id] = {'run_id': run_id, 'demographic_model': dm, 'population': pop}
        return json.dumps({'job_id': job_id, 'status': 'done', 'reused': True}), 200, {'Content-Type': 'application/json'}

    def _rescan_worker():
        status_path = os.path.join(run_dir, f'rescan_{job_id}.json')
        def write_status(st: str, err: Optional[str] = None):
            try:
                payload = {'job_id': job_id, 'status': st}
                if err is not None:
                    payload['error'] = err
                # Always include dm/pop in status for clients to update UI
                payload['demographic_model'] = dm
                payload['population'] = pop
                payload['grid'] = str(grid)
                with open(status_path, 'w', encoding='utf-8') as sf:
                    sf.write(json.dumps(payload))
            except Exception:
                pass
        try:
            write_status('running')
            append_log(run_dir, f"Rescan start: {dm} / {pop} (job {job_id})")
            # Ensure index exists
            try:
                ensure_hts_index(vcf_path)
            except Exception as ie:
                # If not required, continue; else fail
                if _needs_index(vcf_path)[0]:
                    raise
                append_log(run_dir, f"Indexing skipped/warning: {ie}")
            ra_exe = choose_raisd_exe(vcf_path)
            report_path = run_raisd(run_dir, model_path, vcf_path, chromosome, grid=grid, ra_exe=ra_exe)
            # Ensure we have the named RAiSD report path (copy may have created generic alias already)
            try:
                if not os.path.exists(persisted_path) and os.path.exists(report_path):
                    # If report_path is generic, prefer the named one if present
                    named_guess = os.path.join(run_dir, f"RAiSD_Report.{base_name}")
                    src = named_guess if os.path.exists(named_guess) else report_path
                    if os.path.abspath(src) != os.path.abspath(persisted_path):
                        shutil.copyfile(src, persisted_path)
            except Exception:
                pass
            # Remove RAiSD_Info.<base_name> if present
            try:
                info_named = os.path.join(run_dir, f"RAiSD_Info.{base_name}")
                if os.path.exists(info_named):
                    os.remove(info_named)
            except Exception:
                pass
            # Invalidate cached DataFrame so subsequent metric_data/report_csv reads the new report
            try:
                _cached_report.cache_clear()  # type: ignore[attr-defined]
            except Exception:
                pass
            # Update meta with selected model path and last_rescan details
            try:
                with open(meta_path, 'r', encoding='utf-8') as mf:
                    cur = json.load(mf)
            except Exception:
                cur = {}
            cur['model_path'] = model_path
            cur['plots_ready'] = True
            cur['last_rescan'] = {
                'job_id': job_id,
                'demographic_model': dm,
                'population': pop,
                'grid': grid,
                'timestamp': time.time(),
                'report': persisted_name
            }
            try:
                with open(meta_path, 'w', encoding='utf-8') as mf:
                    mf.write(json.dumps(cur))
            except Exception:
                pass
            # Try to compute and persist the projected expected SFS for this (dm,pop)
            # so clients can render the picked-model trace without a separate fetch.
            try:
                # Determine target sample size
                targ = None
                try:
                    if isinstance(cur.get('match'), dict):
                        tval = cur['match'].get('target')
                        if isinstance(tval, int) and tval > 0:
                            targ = int(tval)
                except Exception:
                    targ = None
                if not targ:
                    # fallback to analysis_input in meta
                    ai = cur.get('analysis_input') or cur.get('uploaded_name')
                    if ai:
                        vcf_candidate = os.path.join(run_dir, ai)
                        if os.path.exists(vcf_candidate):
                            try:
                                targ, _, _, _ = infer_samples_and_ploidy(vcf_candidate)
                            except Exception:
                                targ = None
                if targ and cur.get('species'):
                    species_name = str(cur.get('species'))
                    sfs_table = os.path.join(DATA_DIR, species_name, 'sfs.csv')
                    if os.path.exists(sfs_table):
                        try:
                            df = load_sfs_matrix(sfs_table)
                            try:
                                dfn = normalize_sfs_df(df)
                            except Exception:
                                dfn = df
                            proj = project_expected_df(dfn, int(targ))
                            key = f"{dm}={pop}"
                            if key not in proj.index:
                                # forgiving match
                                def _normk(s):
                                    return re.sub(r'[^A-Za-z0-9]+', '', (s or '').lower())
                                wanted = _normk(key)
                                found = None
                                for idx in proj.index:
                                    if _normk(idx) == wanted:
                                        found = idx
                                        break
                                if found is not None:
                                    key = found
                            if key in proj.index:
                                arr = proj.loc[key].values.tolist()
                                # update match_details.json if present
                                try:
                                    md_path = os.path.join(run_dir, 'match_details.json')
                                    if os.path.exists(md_path):
                                        with open(md_path, 'r', encoding='utf-8') as mdf:
                                            md = json.load(mdf)
                                    else:
                                        md = {}
                                    # update top_matches and all_jsd entries
                                    for list_key in ('top_matches', 'all_jsd'):
                                        if isinstance(md.get(list_key), list):
                                            for ent in md[list_key]:
                                                try:
                                                    if (str(ent.get('demographic_model') or '')).strip().lower() == (str(dm) or '').strip().lower() and (str(ent.get('population') or '')).strip().lower() == (str(pop) or '').strip().lower():
                                                        ent['expected_sfs'] = arr
                                                except Exception:
                                                    continue
                                    # persist updated match_details.json
                                    try:
                                        with open(md_path, 'w', encoding='utf-8') as mdf:
                                            mdf.write(json.dumps(md))
                                    except Exception:
                                        pass
                                except Exception:
                                    pass
                                # Also update result_meta.json's match.all_jsd if present
                                try:
                                    if isinstance(cur.get('match'), dict) and isinstance(cur['match'].get('all_jsd'), list):
                                        updated = False
                                        for ent in cur['match']['all_jsd']:
                                            try:
                                                if (str(ent.get('demographic_model') or '')).strip().lower() == (str(dm) or '').strip().lower() and (str(ent.get('population') or '')).strip().lower() == (str(pop) or '').strip().lower():
                                                    ent['expected_sfs'] = arr
                                                    updated = True
                                            except Exception:
                                                continue
                                        if updated:
                                            try:
                                                with open(meta_path, 'w', encoding='utf-8') as mf:
                                                    mf.write(json.dumps(cur))
                                            except Exception:
                                                pass
                                except Exception:
                                    pass
                        except Exception:
                            pass
            except Exception:
                # non-fatal: best-effort persistence only
                pass
            append_log(run_dir, f"Rescan complete: {dm}/{pop} (job {job_id})")
            # Final status with report filename for client
            try:
                with open(status_path, 'w', encoding='utf-8') as sf:
                    sf.write(json.dumps({
                        'job_id': job_id,
                        'status': 'done',
                        'demographic_model': dm,
                        'population': pop,
                        'report': persisted_name
                    }))
            except Exception:
                pass
        except Exception as e:
            append_log(run_dir, f"Rescan error ({job_id}): {e}")
            write_status('error', str(e))

    try:
        fut = _executor.submit(_rescan_worker)
        _rescan_futures[job_id] = fut
        _rescan_meta[job_id] = {'run_id': run_id, 'demographic_model': dm, 'population': pop}
        return json.dumps({'job_id': job_id, 'status': 'started'}), 202, {'Content-Type': 'application/json'}
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/rescan_status')
def runs_rescan_status(run_id):
    job_id = request.args.get('job') or ''
    if not job_id:
        return json.dumps({'error': 'job id required'}), 400, {'Content-Type': 'application/json'}
    meta = _rescan_meta.get(job_id)
    if not meta or meta.get('run_id') != run_id:
        return json.dumps({'error': 'unknown job'}), 404, {'Content-Type': 'application/json'}
    run_dir = os.path.join(RUNS_DIR, run_id)
    status_path = os.path.join(run_dir, f'rescan_{job_id}.json')
    st = {'job_id': job_id, 'status': 'unknown'}
    try:
        if os.path.exists(status_path):
            with open(status_path, 'r', encoding='utf-8') as sf:
                st = json.load(sf)
        else:
            fut = _rescan_futures.get(job_id)
            if fut is None:
                st = {'job_id': job_id, 'status': 'unknown'}
            elif fut.done():
                st = {'job_id': job_id, 'status': 'done'}
            else:
                st = {'job_id': job_id, 'status': 'running'}
    except Exception as e:
        st = {'job_id': job_id, 'status': 'error', 'error': str(e)}
    return json.dumps(st), 200, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/cleanup', methods=['POST'])
def runs_cleanup(run_id):
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    # Schedule deletion with a short grace period so accidental refresh won't lose the run
    def _do_delete():
        try:
            # If there was activity since scheduling, skip delete
            last = _run_last_seen.get(run_id)
            if last and (int(time.time()) - last) < CLEANUP_GRACE_SECONDS:
                return
            shutil.rmtree(run_dir)
            _run_last_seen.pop(run_id, None)
        except Exception:
            pass
        finally:
            try:
                _scheduled_cleanup.pop(run_id, None)
            except Exception:
                pass
    try:
        # If there is an existing timer, cancel and reschedule
        old = _scheduled_cleanup.get(run_id)
        if old:
            try:
                old.cancel()
            except Exception:
                pass
        t = threading.Timer(max(0, CLEANUP_GRACE_SECONDS), _do_delete)
        t.daemon = True
        _scheduled_cleanup[run_id] = t
        t.start()
        return json.dumps({'status': 'scheduled', 'grace_seconds': CLEANUP_GRACE_SECONDS}), 200, {'Content-Type': 'application/json'}
    except Exception as e:
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}

@app.route('/runs/<run_id>/expected_sfs')
def runs_expected_sfs(run_id):
    """Return expected SFS vector for a specific (dm,pop) projected to this run's target.
    Query params: dm=, pop=
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    touch_run(run_id)
    dm = (request.args.get('dm') or '').strip()
    pop = (request.args.get('pop') or '').strip()
    if not dm or not pop:
        return json.dumps({'error': 'dm and pop are required'}), 400, {'Content-Type': 'application/json'}
    meta_path = os.path.join(run_dir, 'result_meta.json')
    try:
        with open(meta_path, 'r', encoding='utf-8') as mf:
            meta = json.load(mf)
    except Exception as e:
        return json.dumps({'error': f'failed to read meta: {e}'}), 500, {'Content-Type': 'application/json'}
    species = meta.get('species') or ''
    # Determine target number of copies (n) from prior computation; fallback to recompute from analysis_input
    target = None
    try:
        m = meta.get('match') or {}
        if isinstance(m, dict):
            t = m.get('target')
            if isinstance(t, int) and t > 0:
                target = t
    except Exception:
        pass
    if not target:
        try:
            ai = meta.get('analysis_input')
            if ai:
                vcf_path = os.path.join(run_dir, ai)
                if os.path.exists(vcf_path):
                    target, _, _, _ = infer_samples_and_ploidy(vcf_path)
        except Exception:
            pass
    if not target:
        return json.dumps({'error': 'could not determine target sample size'}), 500, {'Content-Type': 'application/json'}
    # Build expected SFS for this (dm,pop)
    sfs_table = os.path.join(DATA_DIR, species, 'sfs.csv')
    if not os.path.exists(sfs_table):
        return json.dumps({'error': 'SFS table missing'}), 404, {'Content-Type': 'application/json'}
    try:
        df = load_sfs_matrix(sfs_table)
        try:
            df = normalize_sfs_df(df)
        except Exception:
            pass
        proj = project_expected_df(df, int(target))
    except Exception as e:
        return json.dumps({'error': f'failed to compute expected SFS: {e}'}), 500, {'Content-Type': 'application/json'}
    key = f"{dm}={pop}"
    if key not in proj.index:
        # Attempt a forgiving match ignoring spaces/underscores and case
        def norm(s: str) -> str:
            return re.sub(r'[^A-Za-z0-9]+', '', (s or '').lower())
        wanted = norm(key)
        found_key = None
        for idx in proj.index:
            if norm(str(idx)) == wanted:
                found_key = idx
                break
        if found_key is None:
            return json.dumps({'error': 'model/pop not found in SFS table'}), 404, {'Content-Type': 'application/json'}
        key = found_key
    arr = proj.loc[key].values.tolist()
    return json.dumps({'demographic_model': dm, 'population': pop, 'target': int(target), 'expected_sfs': [float(x) for x in arr]}), 200, {'Content-Type': 'application/json'}

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
    # Consider the run ready only if a parsed report exists AND the
    # metadata explicitly marks plots_ready. This avoids returning 'ready'
    # when RAiSD has created a partial report file but parsing/plot readiness
    # hasn't completed.
    ready = _any_report_exists(run_dir) and bool(meta.get('plots_ready', False))
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

    # If user hit the page normally and no report is present yet, show processing page (non-blocking)
    if request.args.get('nowait', '0') != '1' and not _any_report_exists(run_dir):
        try:
            return render_template('processing.html', run_id=run_id)
        except Exception:
            # Fallback: minimal JSON status if template missing
            return json.dumps({'status': 'running'}), 202, {'Content-Type': 'application/json'}

    if request.args.get('nowait', '0') == '1':
        # Non-blocking status check for processing page: use _ready_or_error so
        # we consider both the presence of a parsed report and meta['plots_ready'].
        ready, meta = _ready_or_error(run_dir)
        if not ready:
            return json.dumps({'status': 'running'}), 202, {'Content-Type': 'application/json'}
        return json.dumps({'status': 'ready'}), 200, {'Content-Type': 'application/json'}

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
    # When rendering final results, prefer persisted per-(dm,pop) report so UI points at it
    report_path = _choose_persisted_report(run_dir) or os.path.join(run_dir, RAISD_REPORT)
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

    # Early fallback: if match info absent (e.g. user hit endpoint before worker wrote meta), parse run.log
    if not raw_match or not raw_match.get('demographic_model'):
        try:
            log_path = os.path.join(run_dir, 'run.log')
            if os.path.exists(log_path):
                with open(log_path, 'r', encoding='utf-8', errors='replace') as lf:
                    lines = lf.readlines()[-500:]  # tail only
                import re
                # Parse best match line
                best_re = re.compile(r"Best match -> model: (?P<model>[^,]+), population: (?P<pop>[^,]+), JSD: (?P<jsd>[0-9.eE+-]+), target chromosomes: (?P<target>\d+), ploidy: (?P<ploidy>\d+), mixed_ploidy: (?P<mixed>True|False)")
                for line in reversed(lines):  # search from end
                    m = best_re.search(line)
                    if m:
                        parsed = {
                            'demographic_model': m.group('model').strip(),
                            'population': m.group('pop').strip(),
                            'best_jsd': float(m.group('jsd')),
                            'target': int(m.group('target')),
                            'ploidy': int(m.group('ploidy')),
                            'mixed_ploidy': (m.group('mixed') == 'True')
                        }
                        raw_match.update({k: v for k, v in parsed.items() if v is not None})
                        append_log(run_dir, 'Hydrated best match from log tail')
                        break
                # Parse species from RAiSD command line if still missing
                if not meta.get('species'):
                    cmd_re = re.compile(r"Running RAiSD: .*? -mdl (.*?)/data/([^/]+)/")
                    for line in reversed(lines):
                        m2 = cmd_re.search(line)
                        if m2:
                            species_guess = m2.group(2)
                            meta['species'] = species_guess
                            append_log(run_dir, f'Hydrated species from log: {species_guess}')
                            break
                if raw_match and not meta.get('match'):
                    meta['match'] = raw_match
                    # Persist hydrated meta so subsequent calls are fast
                    try:
                        with open(meta_path, 'w', encoding='utf-8') as mf:
                            mf.write(json.dumps(meta))
                    except Exception:
                        pass
        except Exception as log_e:
            append_log(run_dir, f'Log hydration failed: {log_e}')

    # Recovery: Ensure best match fields are populated (user requirement to always show true values, not n/a)
    try:
        needed_keys = ['demographic_model', 'population', 'best_jsd', 'target', 'ploidy', 'mixed_ploidy']
        missing_or_placeholder = False
        for k in needed_keys:
            v = raw_match.get(k)
            if v in (None, '', 'n/a'):
                missing_or_placeholder = True
                break
        if missing_or_placeholder:
            # Attempt recomputation if we still have the analysis input file
            analysis_input = meta.get('analysis_input')
            if analysis_input:
                vcf_candidate = os.path.join(run_dir, analysis_input)
                if os.path.exists(vcf_candidate):
                    append_log(run_dir, 'Recomputing best match metadata (previous values incomplete)')
                    try:
                        recomputed = compute_best_match(meta.get('species') or '', vcf_candidate)
                        raw_match.update(recomputed)
                        meta['match'] = raw_match
                        # Persist updated meta for future calls
                        with open(meta_path, 'w', encoding='utf-8') as mf:
                            mf.write(json.dumps(meta))
                        append_log(run_dir, 'Best match recomputation successful')
                    except Exception as re:
                        append_log(run_dir, f'Recomputation of best match failed: {re}')
                else:
                    append_log(run_dir, f'Analysis input file missing; cannot recompute match: {vcf_candidate}')
            else:
                append_log(run_dir, 'No analysis_input recorded; cannot recompute match')
            # Secondary fallback: scan run directory for any VCF/BCF that can be used to compute match
            try:
                if any(raw_match.get(k) in (None, '', 'n/a') for k in ['demographic_model', 'population']) and meta.get('species'):
                    patterns = ['*.vcf.gz', '*.vcf', '*.bcf']
                    candidates = []
                    for pat in patterns:
                        candidates.extend(glob.glob(os.path.join(run_dir, pat)))
                    # Prefer compressed VCF if multiple
                    if candidates:
                        candidates.sort(key=lambda p: (0 if p.endswith('.vcf.gz') else (1 if p.endswith('.vcf') else 2), len(p)))
                        for cand in candidates:
                            try:
                                append_log(run_dir, f'Attempting fallback best-match computation using {os.path.basename(cand)}')
                                recomputed2 = compute_best_match(str(meta.get('species') or ''), cand)
                                raw_match.update(recomputed2)
                                meta['match'] = raw_match
                                with open(meta_path, 'w', encoding='utf-8') as mf:
                                    mf.write(json.dumps(meta))
                                append_log(run_dir, 'Fallback best match computation successful')
                                break
                            except Exception as fe:
                                append_log(run_dir, f'Fallback best match attempt failed for {cand}: {fe}')
            except Exception as fb_e:
                append_log(run_dir, f'Unexpected error during fallback best match scan: {fb_e}')
    except Exception as rec_e:
        append_log(run_dir, f'Unexpected error during best match recovery: {rec_e}')

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

    # Preserve extra arrays if present in raw_match
    for extra_key in ('input_sfs', 'best_expected_sfs', 'top_matches', 'all_jsd'):
        if extra_key in raw_match and extra_key not in match:
            match[extra_key] = raw_match[extra_key]

    # Fallback: load match_details.json if arrays missing
    if ('input_sfs' not in match or not match.get('input_sfs')) or ('top_matches' not in match or not match.get('top_matches')):
        try:
            md_path = os.path.join(run_dir, 'match_details.json')
            if os.path.exists(md_path):
                with open(md_path, 'r', encoding='utf-8') as mdf:
                    det = json.load(mdf)
                for k in ('observed_sfs', 'best_expected_sfs', 'top_matches', 'all_jsd'):
                    if k in det:
                        if k == 'observed_sfs':
                            match.setdefault('input_sfs', det[k])
                        else:
                            match.setdefault(k, det[k])
        except Exception as _md_e:
            append_log(run_dir, f"Warning: failed to load match_details.json: {_md_e}")

    # Build/synthesize plots list
    plots = meta.get('plots') or []
    if not plots:
        # fallback: infer exposed metrics directly
        plots = [{'metric': m, 'display': EXPOSED_METRIC_DISPLAY.get(m, m)} for m in DISPLAY_HEADER if m in EXPOSED_METRICS]

    # Apply exposure filter and ensure display label present for each plot dict
    filtered = []
    for p in plots:
        m = p.get('metric') or ''
        if m in EXPOSED_METRICS:
            filtered.append({'metric': m, 'display': p.get('display') or EXPOSED_METRIC_DISPLAY.get(m, m)})
    plots = filtered

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

    plot_url = None  # no static image

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
    # Always include best-match fields with placeholders if unavailable so UI always shows them
    dm = match.get('demographic_model')
    summary_fields.append(('Demographic model', dm if dm else ''))
    pop = match.get('population')
    summary_fields.append(('Population', pop if pop else ''))
    jsd_v = match.get('best_jsd')
    if isinstance(jsd_v, (int, float)):
        summary_fields.append(('JSD', f"{jsd_v:.6f}"))
    else:
        summary_fields.append(('JSD', ''))

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
    # Normalize for UI: strip leading 'chr' so dropdowns show consistent values like '1', 'X', etc.
    norm_chroms = [c.replace('chr', '') if isinstance(c, str) and c.lower().startswith('chr') else c for c in chromosomes]
    # Deduplicate while preserving order
    seen = set(); ordered = []
    for c in norm_chroms:
        if c not in seen:
            seen.add(c); ordered.append(c)
    return json.dumps({'species': species, 'chromosomes': ordered}), 200, {'Content-Type': 'application/json'}


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
        # Available biotypes for UI (as native strings)
        available_biotypes = sorted({str(b) for b in all_hits['biotype'] if pd.notna(b)})
        if biotypes:
            # Robust membership test: compare normalized strings to handle numpy.str_ etc.
            biotypes_norm = [str(b) for b in biotypes]
            hits = all_hits[all_hits['biotype'].astype(str).isin(biotypes_norm)].copy()
        else:
            hits = all_hits
        # Ensure JSON-serializable native Python types (int/str/None) for downstream json.dumps
        # Convert numeric columns to Python ints and replace NaN with None for safety
        if not hits.empty:
            # Force integer types for start/end/overlap_bp where possible
            for col in ('start', 'end', 'overlap_bp'):
                if col in hits.columns:
                    # convert via astype(int) may fail on NaN; coerce then convert
                    hits[col] = hits[col].where(pd.notna(hits[col]), None)
                    hits[col] = hits[col].apply(lambda v: int(v) if v is not None else None)
            # Ensure label and biotype are plain Python strings or None
            for col in ('label', 'biotype'):
                if col in hits.columns:
                    hits[col] = hits[col].where(pd.notna(hits[col]), None)
                    hits[col] = hits[col].apply(lambda v: str(v) if v is not None else None)
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
    # Inspect VCF header contigs immediately and return them so the UI can react
    contigs = []
    suggested = None
    try:
        try:
            vf = pysam.VariantFile(vcf_path)
            contigs = [str(k) for k in vf.header.contigs.keys()]
        except Exception:
            contigs = []
        if contigs:
            if len(contigs) == 1:
                # normalize suggestion to stripped form
                suggested = contigs[0].replace('chr', '') if isinstance(contigs[0], str) and contigs[0].lower().startswith('chr') else contigs[0]
            else:
                # prefer plain numeric contig '1' over 'chr1' when present
                # normalize contig list for suggestion checks
                contig_norm = [c.replace('chr', '') if isinstance(c, str) and c.lower().startswith('chr') else c for c in contigs]
                if '1' in contig_norm:
                    suggested = '1'
                else:
                    # As fallback, pick the first normalized contig
                    suggested = contig_norm[0] if contig_norm else None
    except Exception as e:
        append_log(run_dir, f"Upload-time contig inspection failed: {e}")

    # Return contigs normalized (strip 'chr') for UI convenience but keep original filename
    contigs_norm = [c.replace('chr', '') if isinstance(c, str) and c.lower().startswith('chr') else c for c in contigs]
    payload = {'run_id': run_id, 'filename': fname, 'contigs': contigs_norm}
    if suggested:
        payload['suggested_chromosome'] = suggested
    return json.dumps(payload), 200, {'Content-Type': 'application/json'}

if __name__ == '__main__':
    import argparse
    host = os.environ.get('HOST', '0.0.0.0')
    port = int(os.environ.get('PORT', '5000'))
    debug = bool(os.environ.get('FLASK_DEBUG', '0') == '1')
    threaded = bool(os.environ.get('FLASK_THREADED', '1') == '1')

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--gpu', dest='gpu', action='store_true', help='Request GPU-accelerated RAiSD-AI scanning when available')
    # parse known args so we don't interfere with Flask's own argv handling
    args, _others = parser.parse_known_args()
    # Check if user requested GPU; if so, try to detect one. If none found, warn and ignore.
    if args.gpu:
        try:
            available = _detect_gpu()
        except Exception:
            available = False
        if available:
            USE_GPU_FOR_RAISD = True
            app.logger.info('GPU requested and appears available: RAiSD-AI will be invoked with -gpu')
        else:
            USE_GPU_FOR_RAISD = False
            # Inform user immediately on stderr/stdout as well as app logger
            msg = 'No GPU detected on this host; ignoring --gpu and running RAiSD on CPU.'
            try:
                print(msg)
            except Exception:
                pass
            app.logger.warning(msg)
    # Validate host name/service before attempting to bind. Some environments set HOST
    # to values that cannot be resolved which leads to an immediate socket error like
    # "Name or service not known" when Flask attempts to start. Try a pre-check and
    # fall back to a safe default (0.0.0.0) if resolution fails.
    try:
        import socket
        try:
            # getaddrinfo will raise socket.gaierror on unresolved names
            socket.getaddrinfo(host, port)
        except Exception as gai_e:
            app.logger.warning(f"HOST value '{host}' failed to resolve: {gai_e}; falling back to '0.0.0.0'")
            host = '0.0.0.0'
    except Exception:
        # If socket import or check fails for any reason, proceed with given host
        pass

    app.run(host=host, port=port, debug=debug, threaded=threaded, use_reloader=False)
