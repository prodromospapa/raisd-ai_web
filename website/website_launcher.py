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
import copy
from concurrent.futures import ThreadPoolExecutor, Future
import matplotlib.pyplot as plt
import pysam
import signal

from collections import Counter
from math import lgamma
from scipy.special import gammaln
import math

from flask import (
    Flask, render_template, request, redirect, url_for,
    send_from_directory, flash, Response
)
from werkzeug.utils import secure_filename

# ------------------------------ Config ------------------------------
BASE_DIR      = os.path.abspath(os.path.dirname(__file__))
DATA_DIR      = os.path.join(BASE_DIR, "../data")
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

# Server-side combined precompute options (opt-in via env vars)
PRECOMPUTE_COMBINED_ENABLED = str(os.environ.get('SERVER_PRECOMPUTE_COMBINED', '')).lower() in ('1', 'true', 'yes')
_PRECOMPUTE_COMBINED_METRICS = os.environ.get('SERVER_COMBINED_METRICS')
PRECOMPUTE_COMBINED_METRICS = None
if _PRECOMPUTE_COMBINED_METRICS:
    PRECOMPUTE_COMBINED_METRICS = [s.strip() for s in _PRECOMPUTE_COMBINED_METRICS.split(',') if s.strip()]


@app.route('/ext_gene_link')
def ext_gene_link():
    """Redirect helper: given a gene label and optional gene_id and species,
    choose a safe external page and redirect the browser there. This centralizes
    sanitization and avoids exposing problematic query parameters to clients.

    Params: label, gene_id, species
    """
    label = (request.args.get('label') or '').strip()
    gene_id = (request.args.get('gene_id') or '').strip()
    species = (request.args.get('species') or '').strip()
    # Sanitize: keep only alnum, whitespace, underscore, hyphen
    def clean(s):
        return re.sub(r'[^A-Za-z0-9_\-\s]+', ' ', s or '').strip()

    cleaned = clean(label) or clean(gene_id) or 'gene'
    q = re.sub(r'\s+', ' ', cleaned)
    q_enc = q.replace(' ', '+')
    # Prefer GeneCards search by default
    target = f'https://www.genecards.org/Search/Keyword?query={q_enc}'
    # If the gene_id looks like an Ensembl id, prefer Ensembl search/summary
    if gene_id and re.match(r'^ENS[A-Z0-9]+', gene_id, re.IGNORECASE):
        sp = species.replace(' ', '_') if species else ''
        # Ensembl search results URL or direct summary when gene_id known
        if sp:
            target = f'https://www.ensembl.org/{sp}/Gene/Summary?g={gene_id}'
        else:
            target = f'https://www.ensembl.org/Search/Results?q={gene_id}'
    # Fallback to NCBI gene search if label was extremely unusual
    if not re.search(r'[A-Za-z0-9]', cleaned):
        target = f'https://www.ncbi.nlm.nih.gov/gene/?term={q_enc}'
    return redirect(target, code=302)

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

# Additional distance metrics (match names used in UI)
def _prep_pq(p, q, eps=1e-12):
    p = np.asarray(p, float); q = np.asarray(q, float)
    p = np.clip(p, 0, None); q = np.clip(q, 0, None)
    ps = p.sum(); qs = q.sum()
    if ps <= 0 or not np.isfinite(ps): p = np.ones_like(p); ps = p.sum()
    if qs <= 0 or not np.isfinite(qs): q = np.ones_like(q); qs = q.sum()
    p /= ps; q /= qs
    p = np.clip(p, eps, 1.0); q = np.clip(q, eps, 1.0)
    return p, q

def l1(p, q):
    p, q = _prep_pq(p, q)
    return float(np.abs(p - q).sum())

def l2(p, q):
    p, q = _prep_pq(p, q)
    return float(np.sqrt(((p - q) ** 2).sum()))

def cosine_distance(p, q, eps=1e-12):
    p = np.asarray(p, float); q = np.asarray(q, float)
    p = np.clip(p, 0, None); q = np.clip(q, 0, None)
    num = float((p * q).sum())
    den = math.sqrt(float((p * p).sum())) * math.sqrt(float((q * q).sum()))
    if den <= eps: return 1.0
    return 1.0 - (num / den)

def hellinger(p, q):
    p, q = _prep_pq(p, q)
    return float(np.linalg.norm(np.sqrt(p) - np.sqrt(q)) / math.sqrt(2.0))

def bhattacharyya(p, q, eps=1e-12):
    p, q = _prep_pq(p, q, eps)
    bc = float(np.sum(np.sqrt(p * q)))
    bc = min(max(bc, eps), 1.0)
    return float(-math.log(bc))

def bray_curtis(p, q, eps=1e-12):
    p, q = _prep_pq(p, q, eps)
    num = float(np.abs(p - q).sum())
    den = float((p + q).sum())
    if den <= eps: return 0.0
    return num / den

def canberra(p, q, eps=1e-12):
    p, q = _prep_pq(p, q, eps)
    num = np.abs(p - q)
    den = np.abs(p) + np.abs(q) + eps
    return float((num / den).sum())

def chi_square(p, q, eps=1e-12):
    p, q = _prep_pq(p, q, eps)
    return float(0.5 * (((p - q) ** 2) / (p + q + eps)).sum())

def corr_distance(p, q, eps=1e-12):
    p = np.asarray(p, float); q = np.asarray(q, float)
    if np.allclose(p, p.mean()) or np.allclose(q, q.mean()):
        return 1.0
    r = np.corrcoef(p, q)[0, 1]
    if not np.isfinite(r): return 1.0
    return float(1.0 - r)

def ks_distance(p, q):
    p, q = _prep_pq(p, q)
    P = np.cumsum(p); Q = np.cumsum(q)
    return float(np.max(np.abs(P - Q)))

def aitchison(p, q, eps=1e-12):
    p, q = _prep_pq(p, q, eps)
    lp = np.log(p); lq = np.log(q)
    cp = lp - lp.mean(); cq = lq - lq.mean()
    return float(np.linalg.norm(cp - cq))

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

@lru_cache(maxsize=256)
def _load_annotation(species: str, chromosome: str, ensembl_version: Optional[str] = None) -> pd.DataFrame:
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
    # If an explicit ensembl_version was requested, prefer files under
    # annotation/ensembl_<version>/<chr>.tsv, falling back to legacy paths.
    ver_variants = []
    if ensembl_version:
        ver_variants.append(os.path.join('{species}', 'annotation', f'ensembl_{ensembl_version}', '{chr}.tsv'))
    # Always try the legacy single annotation directory patterns too
    ver_variants.append(os.path.join('{species}', 'annotation', '{chr}.tsv'))

    df = None
    for cand in candidates:
        if not cand or cand in seen:
            continue
        seen.add(cand)
        found = False
        # Try versioned then legacy variants
        for vt in ver_variants:
            path = os.path.join(DATA_DIR, vt.format(species=cand, chr=str(chromosome)))
            tried_paths.append(path)
            if os.path.exists(path):
                try:
                    df = pd.read_csv(path, sep="\t", dtype={"chromosome": str})
                    app.logger.info(f"Loaded annotation file: {_shorten_path(path)}")
                except Exception as read_e:
                    app.logger.exception(f"Failed to read annotation file {_shorten_path(path)}: {read_e}")
                    raise
                found = True
                break
        if found:
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

    # Safety: ensure df was loaded (static analyzers may not follow the earlier raise)
    if df is None:
        raise FileNotFoundError(f"Annotation file not found for species={species} chromosome={chromosome}")

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

def _genes_in_window(species: str, chromosome: str, start: int, end: int, biotype: Optional[str] = None, ensembl_version: Optional[str] = None) -> pd.DataFrame:
    lo, hi = (int(start), int(end)) if start <= end else (int(end), int(start))
    # Default behaviour uses no explicit ensembl version. Callers may pass
    # species values like 'Homo_sapiens' or include an explicit version as
    # a suffix (handled upstream). This function remains simple and callers
    # should request specific versions via the genes endpoint.
    genes = _load_annotation(species, str(chromosome), ensembl_version=ensembl_version)
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
def _logC(n, k):
    """Stable log binomial coefficient using gammaln.

    Accepts scalars or numpy arrays for n and k.
    """
    return gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)


def _pmf_i(i, nproj, N):
    """Return (jmin, jmax, pmf_array) for projecting a derived-count i from
    sample size N down to nproj haplotypes.

    The pmf_array has length nproj+1 and is normalized; entries outside the
    support jmin..jmax are zero. Computation is done in log-space with a
    log-sum-exp stabilization to avoid underflow for large N.
    """
    jmin = max(0, nproj + i - N)
    jmax = min(i, nproj)
    if jmax < jmin:
        return jmin, jmax, np.zeros(nproj + 1, float)
    js = np.arange(jmin, jmax + 1, dtype=int)
    denom = _logC(N, i)
    # compute log-probabilities for supported js
    logp = _logC(nproj, js) + _logC(N - nproj, i - js) - denom
    mx = np.max(logp)
    if not np.isfinite(mx):
        return jmin, jmax, np.zeros(nproj + 1, float)
    pmf_vals = np.exp(logp - mx)
    out = np.zeros(nproj + 1, float)
    out[js] = pmf_vals
    s = out.sum()
    if s > 0:
        out /= s
    return jmin, jmax, out

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
    cache = {}
    rows = []
    for _, row in df.iterrows():
        # Build full-length p vector indexed 0..N
        p_full = np.zeros(N + 1, float)
        for c, v in row.items():
            idx = int(c)
            if 0 <= idx <= N:
                p_full[idx] = float(v)
        # We will build a full-length vector indexed 0..N and then project into
        # target+1 bins (0..target). However the SFS table and UI use polymorphic
        # bins (1..n-1). After projection we will drop the invariant bins (0 and
        # target) and return columns 1..(target-1).
        z = np.zeros(target + 1, float)
        for i, w in enumerate(p_full):
            if w <= 0:
                continue
            if i not in cache:
                cache[i] = _pmf_i(i, target, N)
            jmin, jmax, pmf = cache[i]
            z += w * pmf
        # drop invariant bins (0 and target) to produce polymorphic-only vector
        if target <= 1:
            polymorphic = np.asarray([], float)
        else:
            polymorphic = z[1:target]
        rows.append(norm(polymorphic))
    # columns correspond to bins 1..(target-1)
    cols_out = [str(k) for k in range(1, target)]
    return pd.DataFrame(rows, index=df.index, columns=cols_out)

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

# Desired display / report order required by RAiSD AI report UI:
# r"$\\mu_{Var}$, r"$\\mu_{SFS}$, r"$\\mu_{LD}$, r"$\\mu$", r"$sweep_{TR}$", r"$\\mu_{var}^{sweep_{TR}}$"
# Keep order exactly as requested by user so display/reports follow this sequence.
DISPLAY_HEADER = [
    r"$\mu_{Var}$",
    r"$\mu_{SFS}$",
    r"$\mu_{LD}$",
    r"$\mu$",
    r"$sweep_{TR}$",
    r"$\mu_{var}^{sweep_{TR}}$",
]

# UI requirement: only expose a subset of metrics in the dropdown (plots for others
# may still be generated internally). Use the raw RAiSD header LaTeX strings here.
# Exposed metrics should follow the same logical ordering (μ then sweep_TR).
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
    # Allow clients to request a specific persisted report via query params
    # (report=<basename>) or via dm/pop shorthand (dm=..., pop=...). If the
    # requested persisted report is missing, fall back to the best available
    # persisted or generic report instead of returning 404 to avoid broken
    # links when the UI selects a model that hasn't been rescanned yet.
    req_report = (request.args.get('report') or '').strip()
    req_dm = (request.args.get('dm') or '').strip()
    req_pop = (request.args.get('pop') or '').strip()

    report_path = None
    # 1) honor explicit `report=` (basename only)
    if req_report:
        # Normalize to avoid directory traversal
        rp = os.path.basename(req_report)
        cand = os.path.join(run_dir, rp)
        if os.path.exists(cand):
            report_path = cand
        else:
            app.logger.debug(f"Requested report {rp} not found for run {run_id}; will fallback")

    # 2) honor dm/pop pair -> RAiSD_Report.<dm>_<pop>
    if report_path is None and req_dm and req_pop:
        safe_dm = re.sub(r'[^A-Za-z0-9_\-]+', '_', req_dm).strip('_')
        safe_pop = re.sub(r'[^A-Za-z0-9_\-]+', '_', req_pop).strip('_')
        persisted_name = f"RAiSD_Report.{safe_dm}_{safe_pop}"
        cand = os.path.join(run_dir, persisted_name)
        if os.path.exists(cand):
            report_path = cand
        else:
            app.logger.debug(f"Requested dm/pop persisted report {persisted_name} not found; will fallback")

    # 3) prefer any persisted report (best match) or generic alias
    if report_path is None:
        report_path = _choose_persisted_report(run_dir) or os.path.join(run_dir, RAISD_REPORT)
        if not os.path.exists(report_path):
            found = _find_and_normalize_report(run_dir, prefix='results')
            if found:
                report_path = found
    # If still missing, return a friendly error (avoid raw 404 HTML)
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
    # Name the CSV after the underlying report when available (e.g., RAiSD_Report.<dm>_<pop>.csv)
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
def _compute_combined_scores_server(all_distances, metrics_to_use=None):
    """Compute normalized average-rank combined scores for a list of all_distances.

    all_distances: list of {'demographic_model', 'population', 'distances': {metric:val}}
    metrics_to_use: list of metric names to include (default: all keys found)
    Returns: list of normalized combined scores aligned to all_distances (index order)
    """
    try:
        if not all_distances:
            return []
        N = len(all_distances)
        # Determine metrics
        sample = all_distances[0].get('distances', {}) if all_distances else {}
        if metrics_to_use is None:
            metrics = [m for m in sample.keys()]
        else:
            metrics = [m for m in metrics_to_use if m in sample]
        if not metrics:
            return [None] * N
        # Build per-metric ranks (1..N)
        ranks = {m: [0] * N for m in metrics}
        for m in metrics:
            arr = []
            for idx, entry in enumerate(all_distances):
                v = entry.get('distances', {}).get(m, None)
                arr.append((idx, float('inf') if v is None else v))
            arr.sort(key=lambda x: (x[1], x[0]))
            for rnk, (idx, _) in enumerate(arr, start=1):
                ranks[m][idx] = rnk
        # Average ranks and normalize
        avg = []
        for i in range(N):
            s = 0
            for m in metrics:
                s += ranks[m][i] if ranks[m][i] else N
            avg_rank = s / len(metrics)
            avg.append(avg_rank / N)
        return avg
    except Exception:
        return [None] * (len(all_distances) if all_distances else 0)


def compute_best_match(species: str, vcf_path: str, precompute_combined: bool = False, combined_metrics: Optional[list] = None, primary_metric: Optional[str] = None):
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
        # observed SFS should match the polymorphic bins in the SFS table:
        # use bins 1..(target-1). If target is small, fall back to available polymorphic bins.
        if target <= 1:
            # no polymorphic bins; normalize empty array to avoid errors
            q = np.asarray([], float)
        else:
            # sfs is length target+1 (0..target); take 1..target-1
            q = norm(sfs[1:target])
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
            # project_expected_df now returns polymorphic bins (1..N_ref-1), so obs_proj
            # corresponds to those bins already; normalize to be safe
            q = norm(obs_proj)
        else:
            # Unknown ValueError - propagate
            raise
    # First compute JSD per projected row so we can still return JSD values
    pop_jsd_pairs = []  # list of (full_name, jsd_val)
    for pop in proj.index:
        try:
            d = jsd(proj.loc[pop].values, q)
        except Exception:
            d = float('inf')
        pop_jsd_pairs.append((pop, d))
    pop_jsd_pairs.sort(key=lambda x: x[1])
    best_jsd = pop_jsd_pairs[0][1]

    # Compute a richer per-population distance dictionary for the UI to use
    metrics_available = [
        "JSD", "L1", "L2", "Cosine", "Hellinger", "Bhattacharyya",
        "BrayCurtis", "Canberra", "ChiSquare", "CorrDist", "KS", "Aitchison",
    ]
    all_distances = []
    all_jsd = []
    for name, jsd_val in pop_jsd_pairs:
        try:
            dm, popn = name.split('=')
        except Exception:
            dm, popn = name, ''
        distances = {}
        # compute metric values with best-effort exception handling
        try:
            distances['JSD'] = jsd(proj.loc[name].values, q)
        except Exception:
            distances['JSD'] = None
        try:
            distances['L1'] = l1(proj.loc[name].values, q)
        except Exception:
            distances['L1'] = None
        try:
            distances['L2'] = l2(proj.loc[name].values, q)
        except Exception:
            distances['L2'] = None
        try:
            distances['Cosine'] = cosine_distance(proj.loc[name].values, q)
        except Exception:
            distances['Cosine'] = None
        try:
            distances['Hellinger'] = hellinger(proj.loc[name].values, q)
        except Exception:
            distances['Hellinger'] = None
        try:
            distances['Bhattacharyya'] = bhattacharyya(proj.loc[name].values, q)
        except Exception:
            distances['Bhattacharyya'] = None
        try:
            distances['BrayCurtis'] = bray_curtis(proj.loc[name].values, q)
        except Exception:
            distances['BrayCurtis'] = None
        try:
            distances['Canberra'] = canberra(proj.loc[name].values, q)
        except Exception:
            distances['Canberra'] = None
        try:
            distances['ChiSquare'] = chi_square(proj.loc[name].values, q)
        except Exception:
            distances['ChiSquare'] = None
        try:
            distances['CorrDist'] = corr_distance(proj.loc[name].values, q)
        except Exception:
            distances['CorrDist'] = None
        try:
            distances['KS'] = ks_distance(proj.loc[name].values, q)
        except Exception:
            distances['KS'] = None
        try:
            distances['Aitchison'] = aitchison(proj.loc[name].values, q)
        except Exception:
            distances['Aitchison'] = None

        # Build legacy all_jsd entry (keep expected_sfs when available)
        entry = {'demographic_model': dm, 'population': popn, 'jsd': float(jsd_val)}
        try:
            entry_expected = proj.loc[name].values.tolist()
            entry['expected_sfs'] = entry_expected
        except Exception:
            pass
        all_jsd.append(entry)

        all_distances.append({'demographic_model': dm, 'population': popn, 'distances': distances})

    # If requested, compute combined normalized average-rank scores server-side
    if precompute_combined:
        try:
            combined_scores = _compute_combined_scores_server(all_distances, combined_metrics)
            if combined_scores and len(combined_scores) == len(all_distances):
                for i, val in enumerate(combined_scores):
                    try:
                        # Ensure distances dict exists
                        if 'distances' not in all_distances[i] or all_distances[i]['distances'] is None:
                            all_distances[i]['distances'] = {}
                        all_distances[i]['distances']['Combined'] = (None if val is None else float(val))
                    except Exception:
                        # Non-fatal: continue for other entries
                        continue
                if 'Combined' not in metrics_available:
                    metrics_available.append('Combined')
        except Exception:
            # best-effort: ignore failures and return without Combined
            pass
    # Determine which metric should be used to pick the primary/top matches
    sel_metric = (primary_metric or 'JSD').strip()
    # If Combined requested but not present, attempt to compute it now using provided combined_metrics
    if sel_metric == 'Combined' and 'Combined' not in metrics_available:
        try:
            combined_scores = _compute_combined_scores_server(all_distances, combined_metrics)
            if combined_scores and len(combined_scores) == len(all_distances):
                for i, val in enumerate(combined_scores):
                    try:
                        if 'distances' not in all_distances[i] or all_distances[i]['distances'] is None:
                            all_distances[i]['distances'] = {}
                        all_distances[i]['distances']['Combined'] = (None if val is None else float(val))
                    except Exception:
                        continue
                if 'Combined' not in metrics_available:
                    metrics_available.append('Combined')
        except Exception:
            pass

    # Build an ordering by the selected metric (fall back to JSD when metric missing)
    order = []  # list of (name_key, metric_value)
    for entry in all_distances:
        dm = entry.get('demographic_model')
        popn = entry.get('population')
        name_key = f"{dm}={popn}" if popn else dm
        val = None
        try:
            val = entry.get('distances', {}).get(sel_metric)
        except Exception:
            val = None
        if val is None:
            # fallback to JSD value from all_jsd array
            try:
                # find matching all_jsd entry
                found = next((a for a in all_jsd if a.get('demographic_model') == dm and a.get('population') == popn), None)
                val = found.get('jsd') if found else float('inf')
            except Exception:
                val = float('inf')
        order.append((name_key, float('inf') if val is None else float(val)))
    order.sort(key=lambda x: (x[1], x[0]))

    # Select top (best) by chosen metric
    demographic_model = ''
    population = ''
    best_expected = None
    best_primary_value = None
    if order:
        best_name = order[0][0]
        try:
            best_primary_value = float(order[0][1]) if order[0][1] is not None else None
        except Exception:
            best_primary_value = None
        if '=' in best_name:
            demographic_model, population = best_name.split('=', 1)
        else:
            demographic_model, population = best_name, ''
        try:
            best_expected = proj.loc[best_name].values.tolist()
        except Exception:
            # try alternative lookup by splitting
            try:
                key = f"{demographic_model}={population}" if population else demographic_model
                best_expected = proj.loc[key].values.tolist()
            except Exception:
                best_expected = None
    # Determine best_jsd by looking up the JSD value for the selected-best entry
    best_jsd = None
    try:
        found_jsd = next((a.get('jsd') for a in all_jsd if a.get('demographic_model') == demographic_model and a.get('population') == population), None)
        best_jsd = float(found_jsd) if found_jsd is not None else None
    except Exception:
        best_jsd = None

    # Build top_matches according to selected metric ordering, also include JSD for each entry
    top_matches = []
    for name_key, val in order[:6]:
        try:
            if '=' in name_key:
                dm, popn = name_key.split('=', 1)
            else:
                dm, popn = name_key, ''
        except Exception:
            dm, popn = name_key, ''
        try:
            expected_vec = proj.loc[name_key].values.tolist()
        except Exception:
            expected_vec = None
        item = {'demographic_model': dm, 'population': popn}
        # include the metric value under a predictable key
        item_key = sel_metric.lower().replace(' ', '_') if sel_metric else 'distance'
        try:
            item[item_key] = float(val)
        except Exception:
            item[item_key] = None
        # attach jsd if available
        try:
            jsd_val = next((a.get('jsd') for a in all_jsd if a.get('demographic_model') == dm and a.get('population') == popn), None)
            item['jsd'] = float(jsd_val) if jsd_val is not None else None
        except Exception:
            item['jsd'] = None
        if expected_vec is not None:
            item['expected_sfs'] = expected_vec
        top_matches.append(item)

    return {
        'target': target,
        'ploidy': ploidy,
        'mixed_ploidy': mixed,
        'ploidy_counts': counts,
        'demographic_model': demographic_model,
        'population': population,
        'best_jsd': best_jsd,
        'selected_metric': sel_metric,
        'best_selected_value': best_primary_value,
        'input_sfs': q.tolist(),
        'best_expected_sfs': best_expected,
        'top_matches': top_matches,
        'all_jsd': all_jsd,
        'all_distances': all_distances,
        'metrics_available': metrics_available,
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
    2. any files matching `RAiSD_Report.<dm>_<pop>.results` (newest)
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
    # Prefer the newest matching the "<dm>_<pop>" or "<dm>" patterns; fall back to any RAiSD_Report.*
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
    # We'll pass <dm>_<pop> (or <dm>) so RAiSD creates RAiSD_Report.<dm>_<pop>
    # and RAiSD_Info.<dm>_<pop>.
    name_base = 'results'
    try:
        if dm:
            if pop:
                # Use single underscore to join demography and population
                name_base = f"{_sanitize_name(dm)}_{_sanitize_name(pop)}"
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
    # If RAiSD didn't produce a named report (RAiSD_Report.<name_base>),
    # ensure we persist one by copying the generic RAiSD_Report.results to
    # RAiSD_Report.<name_base>. This guarantees downstream code can look
    # for per-(dm,pop) persisted reports even when the RAiSD executable
    # didn't honor the -n flag or produced only the generic filename.
    try:
        if not os.path.exists(expected_report_named) and os.path.exists(report_path):
            try:
                shutil.copyfile(report_path, expected_report_named)
                append_log(run_dir, f"Created persisted named report: {os.path.basename(expected_report_named)} from {os.path.basename(report_path)}")
            except Exception as ce:
                append_log(run_dir, f"Warning: failed to create persisted named report: {ce}")
    except Exception:
        # Non-fatal: continue even if persistence attempt fails
        pass
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
        # Also remove any RAiSD_Grid.<name_base> artifact (file or directory)
        grid_named = os.path.join(run_dir, f"RAiSD_Grid.{name_base}")
        try:
            if os.path.isdir(grid_named):
                shutil.rmtree(grid_named, ignore_errors=True)
                append_log(run_dir, f"Removed named grid directory: {os.path.basename(grid_named)}")
            elif os.path.exists(grid_named):
                try:
                    os.remove(grid_named)
                    append_log(run_dir, f"Removed named grid file: {os.path.basename(grid_named)}")
                except Exception:
                    pass
        except Exception as ce:
            append_log(run_dir, f"Warning: failed to remove named RAiSD_Grid.{name_base}: {ce}")
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
    # Provide a list of available distance metrics for the initial form
    metrics_available = [
        "JSD", "L1", "L2", "Cosine", "Hellinger", "Bhattacharyya",
        "BrayCurtis", "Canberra", "ChiSquare", "CorrDist", "KS", "Aitchison",
    ]
    return render_template('index.html', species_options=species_options, metrics_available=metrics_available)

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

    # Capture metric selections from the form (optional)
    distance_metric = (request.form.get('distance_metric') or '').strip()
    # combined_metrics may be a single string or multiple; ensure list
    combined_metrics_raw = request.form.getlist('combined_metrics') if hasattr(request.form, 'getlist') else []
    if not combined_metrics_raw:
        cm_single = (request.form.get('combined_metrics') or '').strip()
        combined_metrics_raw = [s.strip() for s in cm_single.split(',') if s.strip()] if cm_single else []

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
                            'distance_metric': distance_metric,
                            'combined_metrics': combined_metrics_raw,
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
            # Try to detect single contig now to prefill the first page
            detected_chr = None
            detected_contig = None
            try:
                if _exe_available(BCFTOOLS):
                    p1 = subprocess.Popen([BCFTOOLS, 'query', '-f', '%CHROM\n', vcf_path], stdout=subprocess.PIPE, env=_bcftools_env(), text=True)
                    assert p1.stdout is not None
                    out = subprocess.check_output(['sort', '-u'], stdin=p1.stdout, text=True)
                    p1.wait()
                    qlines = [ln.strip() for ln in out.splitlines() if ln.strip()]
                    if len(qlines) == 1:
                        detected_contig = qlines[0]
                        # set detected_chr only if it's numeric or chr+numeric
                        m = re.sub(r'^chr', '', detected_contig, flags=re.IGNORECASE)
                        if m.isdigit():
                            detected_chr = int(m)
            except Exception:
                detected_chr = None
                detected_contig = None

            meta_early = {
                'run_id': run_id,
                'species': species,
                'chromosome': detected_chr if detected_chr is not None else chromosome,
                'using_contig': detected_contig,
                'grid': grid_str if isinstance(grid_str, str) else str(grid_str),
                'uploaded_name': filename,
                'distance_metric': distance_metric,
                'combined_metrics': combined_metrics_raw,
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
                # Try to detect single contig now to prefill the first page
                detected_chr = None
                detected_contig = None
                try:
                    if _exe_available(BCFTOOLS):
                        p1 = subprocess.Popen([BCFTOOLS, 'query', '-f', '%CHROM\n', vcf_path], stdout=subprocess.PIPE, env=_bcftools_env(), text=True)
                        assert p1.stdout is not None
                        out = subprocess.check_output(['sort', '-u'], stdin=p1.stdout, text=True)
                        p1.wait()
                        qlines = [ln.strip() for ln in out.splitlines() if ln.strip()]
                        if len(qlines) == 1:
                            detected_contig = qlines[0]
                            m = re.sub(r'^chr', '', detected_contig, flags=re.IGNORECASE)
                            if m.isdigit():
                                detected_chr = int(m)
                except Exception:
                    detected_chr = None
                    detected_contig = None

                meta_early = {
                    'run_id': run_id,
                    'species': species,
                    'chromosome': detected_chr if detected_chr is not None else chromosome,
                    'using_contig': detected_contig,
                    'grid': grid_str if isinstance(grid_str, str) else str(grid_str),
                    'uploaded_name': filename,
                    'distance_metric': distance_metric,
                    'combined_metrics': combined_metrics_raw,
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
        # Strict detection using bcftools + sort -u only (per request).
        contigs = []
        single_contig_auto_ok = None

        if not _exe_available(BCFTOOLS):
            flash(f"Chromosome detection requires 'bcftools' to be installed and available. Please install/activate bcftools and try again.")
            try:
                append_log(run_dir, "Removing run directory due to missing bcftools for chromosome detection")
                shutil.rmtree(run_dir, ignore_errors=True)
            except Exception:
                pass
            return redirect(url_for('index'))

        # Run: bcftools query -f '%CHROM\n' <vcf> | sort -u
        try:
            p1 = subprocess.Popen([BCFTOOLS, 'query', '-f', '%CHROM\n', vcf_path], stdout=subprocess.PIPE, env=_bcftools_env(), text=True)
            assert p1.stdout is not None
            out = subprocess.check_output(['sort', '-u'], stdin=p1.stdout, text=True)
            p1.wait()
            if p1.returncode not in (0, None):
                raise RuntimeError('bcftools query failed')
            qlines = [ln.strip() for ln in out.splitlines() if ln.strip()]
            uniq = qlines
            if len(uniq) == 0:
                flash("Uploaded VCF contains no variant records; provide a VCF with variants.")
                try:
                    append_log(run_dir, "Removing run directory due to empty VCF (no variants)")
                    shutil.rmtree(run_dir, ignore_errors=True)
                except Exception:
                    pass
                return redirect(url_for('index'))
            if len(uniq) > 1:
                shown = uniq[:20]
                flash(f"This tool requires VCF data from a single chromosome. Uploaded VCF contains multiple chromosomes: {', '.join(shown)}{'...' if len(uniq)>20 else ''}")
                try:
                    append_log(run_dir, "Removing run directory due to multi-chromosome VCF")
                    shutil.rmtree(run_dir, ignore_errors=True)
                except Exception:
                    pass
                return redirect(url_for('index'))
            # Single contig found via bcftools+sort: use it and continue
            contigs = [uniq[0]]
            # Decide whether it's safe to auto-override a user-supplied numeric chromosome.
            try:
                c = uniq[0]
                if re.match(r'^(?:chr)?\d+$', c, flags=re.IGNORECASE) or re.match(r'^(?:chr)?(?:X|Y|MT|M)$', c, flags=re.IGNORECASE):
                    single_contig_auto_ok = True
                    # Note detected contig; do NOT overwrite user input here.
                    append_log(run_dir, f"Detected single contig via bcftools: {c} (will not override user selection here)")
                else:
                    single_contig_auto_ok = False
                    flash(f"Uploaded VCF has a single contig named '{c}', which is non-standard. Please ensure you select the intended chromosome in the form; no auto-selection was performed.")
            except Exception:
                single_contig_auto_ok = None
        except subprocess.TimeoutExpired:
            flash("Chromosome detection timed out (bcftools). Try again or ensure bcftools is functional.")
            try:
                append_log(run_dir, "Removing run directory due to bcftools timeout")
                shutil.rmtree(run_dir, ignore_errors=True)
            except Exception:
                pass
            return redirect(url_for('index'))
        except Exception as e:
            flash(f"Chromosome detection failed: {e}")
            try:
                append_log(run_dir, f"Removing run directory due to bcftools detection error: {e}")
                shutil.rmtree(run_dir, ignore_errors=True)
            except Exception:
                pass
            return redirect(url_for('index'))

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
                # Attempt to find a single distinct contig if only one present and auto-override allowed
                if len(contigs) == 1 and (single_contig_auto_ok is None or single_contig_auto_ok):
                    chosen_contig = contigs[0]
                    # Record that a single contig was detected; do not log an assertion of overriding.
                    append_log(run_dir, f"Single-contig VCF detected: using_contig={chosen_contig}")
                else:
                    # Try to strip 'chr' prefixes and compare
                    stripped_map = {c[3:]: c for c in contigs if isinstance(c, str) and c.lower().startswith('chr') and len(c) > 3}
                    if numeric_str in stripped_map:
                        chosen_contig = stripped_map[numeric_str]
                    # If still unresolved, see if first record CHROM matches any contig entry
                    if chosen_contig is None and uniq:
                        fr = uniq[0]
                        if fr in contig_set:
                            chosen_contig = fr
                        elif isinstance(fr, str) and fr.lower().startswith('chr') and fr[3:] in stripped_map:
                            chosen_contig = stripped_map[fr[3:]]
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
                # Only override the user's provided chromosome if they left the field empty
                if model_chr.isdigit():
                    if not user_chrom_input:
                        chromosome = int(model_chr)
                        append_log(run_dir, f"Auto-filled chromosome from VCF: using_contig={chosen_contig} -> model_dir_component={chromosome}")
                    else:
                        # Only log the 'not overriding' decision when the detected contig
                        # differs from the user's explicit selection; avoid noisy duplicate logs
                        try:
                            u_norm = str(user_chrom_input).lower() if user_chrom_input is not None else ''
                            c_norm = str(chosen_contig).lower() if chosen_contig is not None else ''
                        except Exception:
                            u_norm = str(user_chrom_input) if user_chrom_input is not None else ''
                            c_norm = str(chosen_contig) if chosen_contig is not None else ''
                        if u_norm != c_norm:
                            append_log(run_dir, f"Detected contig {chosen_contig}, but user selected {user_chrom_input}; not overriding user choice")
                append_log(run_dir, f"Chromosome validated/inferred: user={user_chrom_input} -> using_contig={chosen_contig} -> model_dir_component={chromosome}")
                # Update early metadata so the UI can reflect the detected chromosome
                try:
                    meta_path = os.path.join(run_dir, 'result_meta.json')
                    if os.path.exists(meta_path):
                        try:
                            with open(meta_path, 'r', encoding='utf-8') as mf:
                                meta_json = json.load(mf)
                        except Exception:
                            meta_json = {}
                        meta_json['chromosome'] = chromosome
                        meta_json['using_contig'] = str(chosen_contig)
                        with open(meta_path, 'w', encoding='utf-8') as mf:
                            json.dump(meta_json, mf)
                except Exception:
                    pass
                # Do not flash a user-visible message here; UI will read result_meta.json to show the detected chromosome.
    except Exception as ce:
        append_log(run_dir, f"Chromosome inference warning: {ce}")

    def _analyze_worker(run_dir, filename, vcf_path, species, chromosome, grid, primary_metric=None, combined_metrics=None):
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
            # Preserve the user's metric preferences so downstream readers (runs_final) can
            # use them when recomputing or rendering without relying on the earlier temp file.
            'distance_metric': primary_metric,
            'combined_metrics': combined_metrics or [],
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

            match = compute_best_match(
                species, vcf_path,
                precompute_combined=PRECOMPUTE_COMBINED_ENABLED,
                combined_metrics=combined_metrics if combined_metrics is not None else PRECOMPUTE_COMBINED_METRICS,
                primary_metric=primary_metric
            )
            meta['match'] = match
            # Log both the canonical JSD value (always computed) and the user's
            # selected primary metric (which may be 'Combined' or another distance).
            try:
                sel_metric_log = (match.get('selected_metric') or 'JSD')
                sel_val_raw = match.get('best_selected_value')
                if isinstance(sel_val_raw, (int, float)):
                    sel_val_log = f"{float(sel_val_raw):.6f}"
                else:
                    sel_val_log = ''
            except Exception:
                sel_metric_log = 'JSD'
                sel_val_log = ''
            # Prepare a safe JSD string for logging
            try:
                jsd_raw = match.get('best_jsd')
                if isinstance(jsd_raw, (int, float)):
                    jsd_log = f"{float(jsd_raw):.6f}"
                else:
                    jsd_log = 'n/a'
            except Exception:
                jsd_log = 'n/a'
            append_log(run_dir, (
                f"Best match -> model: {match.get('demographic_model')}, population: {match.get('population')}, "
                f"selected_metric: {sel_metric_log}{('='+sel_val_log) if sel_val_log else ''}, "
                f"JSD: {jsd_log}, target chromosomes: {match.get('target')}, "
                f"ploidy: {match.get('ploidy')}, mixed_ploidy: {match.get('mixed_ploidy')}"
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
        fut = _executor.submit(_analyze_worker, run_dir, filename, vcf_path, species, chromosome, grid, distance_metric or None, combined_metrics_raw or None)
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


@app.route('/runs/<run_id>/match_heavy')
def runs_match_heavy(run_id):
    """Return the heavier match payload for a run (all_distances, all_jsd, match_details).

    This endpoint is intended for asynchronous lazy-loading by the client so that
    the initial results page can render quickly with a lightweight match object.
    The response is gzipped to reduce network transfer size when large arrays are present.
    """
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return json.dumps({'error': 'run not found'}), 404, {'Content-Type': 'application/json'}
    touch_run(run_id)
    meta_path = os.path.join(run_dir, 'result_meta.json')
    match = None
    try:
        if os.path.exists(meta_path):
            with open(meta_path, 'r', encoding='utf-8') as mf:
                meta = json.load(mf)
                match = meta.get('match')
    except Exception:
        match = None

    # Read separately persisted match_details.json if present
    details = None
    details_path = os.path.join(run_dir, 'match_details.json')
    try:
        if os.path.exists(details_path):
            with open(details_path, 'r', encoding='utf-8') as df:
                details = json.load(df)
    except Exception:
        details = None

    payload = {
        'match': match,
        'match_details': details,
    }
    try:
        raw = json.dumps(payload)
    except Exception as e:
        app.logger.exception(f"Failed to serialize match_heavy for run {run_id}: {e}")
        return json.dumps({'error': 'serialization failed'}), 500, {'Content-Type': 'application/json'}

    try:
        gz = gzip.compress(raw.encode('utf-8'))
        headers = {'Content-Type': 'application/json', 'Content-Encoding': 'gzip'}
        return (gz, 200, headers)
    except Exception:
        # Fallback to plain JSON if gzip fails
        return raw, 200, {'Content-Type': 'application/json'}


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

    # Expected RAiSD-produced name components from -n (use single underscore)
    base_name = f"{_sanitize(dm)}_{_sanitize(pop)}"
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
                # Accept either the legacy form or the new form that includes an
                # optional selected_metric field: "selected_metric: NAME[=VALUE], "
                best_re = re.compile(
                    r"Best match -> model: (?P<model>[^,]+), population: (?P<pop>[^,]+), "
                    r"(?:selected_metric:\s*(?P<sel_metric>[^,=]+)(?:=(?P<sel_val>[0-9.eE+-]+))?,\s*)?"
                    r"JSD:\s*(?P<jsd>[0-9.eE+-]+|n/a), target chromosomes: (?P<target>\d+), "
                    r"ploidy: (?P<ploidy>\d+), mixed_ploidy: (?P<mixed>True|False)"
                )
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
                        recomputed = compute_best_match(
                            meta.get('species') or '', vcf_candidate,
                            precompute_combined=PRECOMPUTE_COMBINED_ENABLED,
                            combined_metrics=meta.get('combined_metrics') or PRECOMPUTE_COMBINED_METRICS,
                            primary_metric=(meta.get('distance_metric') or 'JSD')
                        )
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
                                recomputed2 = compute_best_match(
                                    str(meta.get('species') or ''), cand,
                                    precompute_combined=PRECOMPUTE_COMBINED_ENABLED,
                                    combined_metrics=meta.get('combined_metrics') or PRECOMPUTE_COMBINED_METRICS,
                                    primary_metric=meta.get('distance_metric')
                                )
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
    for extra_key in ('input_sfs', 'best_expected_sfs', 'top_matches', 'all_jsd', 'all_distances', 'metrics_available', 'selected_metric', 'best_selected_value'):
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

    # If the user originally requested a non-JSD metric on the main page, ensure
    # the best-match is computed using that metric rather than defaulting to JSD.
    try:
        sel_metric = (meta.get('distance_metric') or '')
        # Only recompute if user explicitly requested a non-JSD metric AND
        # the recorded match was not already computed with that metric (or lacks the selected value).
        if sel_metric and str(sel_metric).strip() and str(sel_metric).strip().upper() != 'JSD':
            current_match = meta.get('match') or {}
            recorded_sel = (current_match.get('selected_metric') or '') if isinstance(current_match, dict) else ''
            recorded_best = current_match.get('best_selected_value') if isinstance(current_match, dict) else None
            need_recompute = False
            if not current_match:
                need_recompute = True
            else:
                # If the recorded selected metric doesn't match what the user asked for,
                # or the recorded value is missing, then recompute. Otherwise skip.
                if str(recorded_sel).strip() != str(sel_metric).strip() or recorded_best in (None, ''):
                    need_recompute = True

            if need_recompute:
                analysis_input = meta.get('analysis_input') or meta.get('uploaded_name')
                if analysis_input:
                    candidate = os.path.join(run_dir, analysis_input)
                    if os.path.exists(candidate):
                        try:
                            append_log(run_dir, f"Recomputing best match using selected metric: {sel_metric}")
                            recomputed = compute_best_match(
                                meta.get('species') or '', candidate,
                                precompute_combined=PRECOMPUTE_COMBINED_ENABLED,
                                combined_metrics=meta.get('combined_metrics') or PRECOMPUTE_COMBINED_METRICS,
                                primary_metric=sel_metric
                            )
                            # Update match and persist to meta so UI and future calls reflect the user's choice
                            match = recomputed
                            meta['match'] = recomputed
                            try:
                                with open(meta_path, 'w', encoding='utf-8') as mf:
                                    mf.write(json.dumps(meta))
                            except Exception:
                                pass
                        except Exception as e:
                            append_log(run_dir, f"Selected-metric recompute failed: {e}")
    except Exception:
        # non-fatal
        pass

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
    # Show only the user's selected primary metric in the Best Match summary.
    try:
        sel_metric = (match.get('selected_metric') or meta.get('distance_metric') or '').strip()
        sel_val = match.get('best_selected_value')
        if not sel_metric:
            # default to JSD when no selection available
            sel_metric = 'JSD'
            sel_val = match.get('best_jsd')
        if isinstance(sel_val, (int, float)):
            summary_fields.append((sel_metric, f"{sel_val:.6f}"))
        else:
            summary_fields.append((sel_metric, ''))
    except Exception:
        # fallback: always show JSD
        try:
            jsd_v = match.get('best_jsd')
            summary_fields.append(('JSD', f"{jsd_v:.6f}" if isinstance(jsd_v, (int, float)) else ''))
        except Exception:
            summary_fields.append(('JSD', ''))

    try:
        # Build a pruned copy of match for embedding in the template. We must avoid
        # including large numeric arrays (expected_sfs etc.) and also avoid passing
        # any non-serializable objects (functions, methods). Use deepcopy to avoid
        # mutating the persisted meta in-place.
        pruned_match = {}
        try:
            pruned_match = copy.deepcopy(match) if isinstance(match, dict) else {}
        except Exception:
            # Fallback: shallow copy
            pruned_match = dict(match) if isinstance(match, dict) else {}

        # Remove unexpectedly large arrays to avoid huge payloads, but keep
        # small SFS vectors so the client can render the SFS comparison modal
        # and mini-plot without an extra fetch. Only drop if length is absurd.
        for key in ('best_expected_sfs', 'input_sfs'):
            try:
                val = pruned_match.get(key)
                if isinstance(val, (list, tuple)) and len(val) > 1024:
                    pruned_match.pop(key, None)
                # otherwise keep small arrays (observed/expected SFS) for client-side plotting
            except Exception:
                pruned_match.pop(key, None)

        # Prune expected_sfs inside all_jsd list entries
        if isinstance(pruned_match.get('all_jsd'), list):
            for ent in (pruned_match.get('all_jsd') or []):
                try:
                    if isinstance(ent, dict):
                        ent.pop('expected_sfs', None)
                except Exception:
                    pass

        # For all_distances keep only demographic_model, population and distances dict
        if isinstance(pruned_match.get('all_distances'), list):
            new_all_dist = []
            for ent in (pruned_match.get('all_distances') or []):
                try:
                    if isinstance(ent, dict):
                        d = {'demographic_model': ent.get('demographic_model'), 'population': ent.get('population'), 'distances': {}}
                        # copy distances but ensure no nested arrays or callables
                        rawd = ent.get('distances') or {}
                        for k2, v2 in (rawd.items() if isinstance(rawd, dict) else []):
                            # avoid callables or big arrays
                            if callable(v2):
                                continue
                            if isinstance(v2, (list, tuple)) and len(v2) > 1024:
                                # skip unexpectedly huge lists
                                continue
                            d['distances'][k2] = v2
                        new_all_dist.append(d)
                except Exception:
                    continue
            pruned_match['all_distances'] = new_all_dist
        # Ensure compatibility: if we have all_distances but no all_jsd, synthesize
        # a lightweight all_jsd list (with JSD values) so legacy templates work.
        try:
            if isinstance(pruned_match.get('all_distances'), list) and not isinstance(pruned_match.get('all_jsd'), list):
                synth = []
                for ent in (pruned_match.get('all_distances') or []):
                    try:
                        dm = ent.get('demographic_model')
                        popn = ent.get('population')
                        jsdv = None
                        try:
                            jsdv = (ent.get('distances') or {}).get('JSD')
                        except Exception:
                            jsdv = None
                        synth.append({'demographic_model': dm, 'population': popn, 'jsd': (None if jsdv is None else float(jsdv))})
                    except Exception:
                        continue
                pruned_match['all_jsd'] = synth
        except Exception:
            pass

        # Conversely, if we have all_jsd but not all_distances, build an
        # all_distances list so the new metric selector can use it (with only JSD present).
        try:
            if isinstance(pruned_match.get('all_jsd'), list) and not isinstance(pruned_match.get('all_distances'), list):
                synth2 = []
                for ent in (pruned_match.get('all_jsd') or []):
                    try:
                        dm = ent.get('demographic_model')
                        popn = ent.get('population')
                        jsdv = ent.get('jsd')
                        synth2.append({'demographic_model': dm, 'population': popn, 'distances': {'JSD': (None if jsdv is None else float(jsdv))}})
                    except Exception:
                        continue
                pruned_match['all_distances'] = synth2
        except Exception:
            pass

        # Ensure metrics_available exists when possible
        try:
            if 'metrics_available' not in pruned_match:
                if isinstance(pruned_match.get('all_distances'), list) and pruned_match['all_distances']:
                    first = pruned_match['all_distances'][0]
                    if isinstance(first.get('distances'), dict):
                        pruned_match['metrics_available'] = list(first.get('distances').keys())
        except Exception:
            pass
        return render_template(
            'result.html',
            run_id=run_id,
            species=meta.get('species', ''),
            chromosome=meta.get('chromosome', ''),
            grid=meta.get('grid', ''),
            match=pruned_match,
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
        # Optional ensembl_version param to select which annotation set to use
        ensembl_version = (request.args.get('ensembl_version') or '').strip() or None
        # Load all hits first to derive available biotypes, then filter locally
        all_hits = _genes_in_window(species, chromosome, start_i, end_i, biotype=None, ensembl_version=ensembl_version)
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


@app.route('/annotation_versions', methods=['GET'])
def annotation_versions():
    """Return available Ensembl versions for a species.

    Looks for data/<species>/annotation/versions.json and returns its contents if present.
    Otherwise scans subdirectories under data/<species>/annotation/ for folders named ensembl_<v>.
    """
    species = (request.args.get('species') or '').strip()
    if not species:
        return json.dumps({'error': 'species required'}), 400, {'Content-Type': 'application/json'}
    base = os.path.join(DATA_DIR, species, 'annotation')
    versions = []
    try:
        # Prefer scanning the annotation directory for subfolders named
        # `ensembl_<digits>` rather than relying on a manifest file. This
        # avoids the need to update versions.json and keeps the website in
        # sync with the filesystem.
        if os.path.isdir(base):
            for name in os.listdir(base):
                m = re.match(r'^ensembl_(\d+)$', name)
                if m:
                    versions.append(m.group(1))
        # Sort numerically when possible so the latest release is last
        def _num_key(s):
            try:
                return int(s)
            except Exception:
                return float('inf')
        versions = sorted([str(v) for v in versions], key=_num_key)
    except Exception as e:
        app.logger.exception(f"annotation_versions error for species={species}: {e}")
        return json.dumps({'error': str(e)}), 500, {'Content-Type': 'application/json'}
    return json.dumps({'species': species, 'ensembl_versions': versions}), 200, {'Content-Type': 'application/json'}

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
    # Inspect VCF chromosomes using the bcftools pipeline the UI expects
    # Preferred command: bcftools query -f '%CHROM\n' <vcf>  (we enumerate unique values in Python)
    contigs = []
    suggested = None
    try:
        if _exe_available(BCFTOOLS):
            try:
                proc = subprocess.run([BCFTOOLS, 'query', '-f', '%CHROM\\n', vcf_path], capture_output=True, text=True, env=_bcftools_env(), timeout=30)
                if proc.returncode == 0:
                    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
                    # sort-unique semantics: stable sort of unique set
                    contigs = sorted(list(dict.fromkeys(lines)))
                else:
                    append_log(run_dir, f"bcftools query failed during upload inspection: rc={proc.returncode} stderr={proc.stderr}")
                    contigs = []
            except Exception as be:
                append_log(run_dir, f"bcftools query raised during upload inspection: {be}")
                contigs = []
        # Fallback to header contigs via pysam if bcftools unavailable or produced no output
        if not contigs:
            try:
                vf = pysam.VariantFile(vcf_path)
                contigs = [str(k) for k in vf.header.contigs.keys()]
            except Exception:
                contigs = []

        # Decide on a suggested chromosome only when there's exactly one clear contig
        if contigs:
            if len(contigs) == 1:
                c0 = contigs[0]
                # Normalize suggested form for UI: strip leading 'chr' if present
                if isinstance(c0, str) and c0.lower().startswith('chr'):
                    suggested = c0[3:]
                else:
                    suggested = c0
            else:
                # Prefer a plain numeric '1' suggestion if present; otherwise, no automatic suggestion
                contig_norm = [c[3:] if isinstance(c, str) and c.lower().startswith('chr') else c for c in contigs]
                if '1' in contig_norm:
                    suggested = '1'
                else:
                    suggested = None
    except Exception as e:
        append_log(run_dir, f"Upload-time contig inspection failed: {e}")

    # Return contigs normalized (strip leading 'chr') for UI convenience but keep original filename
    contigs_norm = [c[3:] if isinstance(c, str) and c.lower().startswith('chr') else c for c in contigs]
    payload = {'run_id': run_id, 'filename': fname, 'contigs': contigs_norm}
    if suggested is not None:
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
