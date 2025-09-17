import os
import sys
import io
import math
import tempfile
import shutil
import time
import random
import re
import atexit
import signal
import subprocess
import threading
import argparse
import shlex
import gzip
import concurrent.futures
import stdpopsim as sps
import errno


# Structured exception raised when RAM limit enforcement aborts the run.
class ExceedsMaxRam(Exception):
    """Raised when total system RAM usage exceeds configured threshold."""
    pass

# optional progress bar (not required)
try:
    from tqdm.auto import tqdm
except Exception:
    try:
        from tqdm import tqdm
    except Exception:
        tqdm = None

# Try to import demes.ms.to_ms once at module import time for reliability
try:
    from demes.ms import to_ms as _demes_to_ms
except Exception:
    _demes_to_ms = None



def _ensure_mutations_on_ts(ts, mu_fallback=1e-8):
    """Ensure a tree sequence has at least one mutation by applying sim_mutations.

    Returns the (possibly modified) tree sequence. Non-fatal on error.
    """
    try:
        if ts is None:
            return ts
        nsites = int(getattr(ts, 'num_sites', 0))
        if nsites > 0:
            return ts
        # Try to add mutations using msprime if available
        import msprime as _msp
        dg = bool(getattr(ts, 'discrete_genome', False)) if hasattr(ts, 'discrete_genome') else False
        return _msp.sim_mutations(ts, rate=mu_fallback, discrete_genome=dg)
    except Exception:
        return ts


def present_size_of(deme):
    """Return a present-day size estimate for a deme object.

    This is a minimal compatibility helper used by builders. It prefers
    attributes commonly found on demes/tskit objects and falls back to 1e4.
    """
    try:
        if hasattr(deme, 'present_size') and deme.present_size is not None:
            return float(deme.present_size)
        if hasattr(deme, 'initial_size') and deme.initial_size is not None:
            return float(deme.initial_size)
        if hasattr(deme, 'start_size') and deme.start_size is not None:
            return float(deme.start_size)
    except Exception:
        pass
    return 1e4


def reorder_put_first(user_order, first):
    """Return a copy of user_order with `first` moved to the front if present."""
    if not user_order:
        return user_order
    out = [x for x in user_order if x != first]
    if first in user_order:
        return [first] + out
    return out


def harmonic_number(n):
    """Return the nth harmonic number H_n = sum_{k=1..n} 1/k (float)."""
    if n <= 0:
        return 0.0
    s = 0.0
    for k in range(1, int(n) + 1):
        s += 1.0 / k
    return s


def to_ms(graph, N0=1.0, samples=None):
    """Produce a minimal ms-like demography string for the given demes.Graph.

    This is a conservative placeholder that emits no special -e* flags and
    leaves demography details minimal so builders can proceed.
    """
    # Prefer using demes.ms.to_ms if demes is available so migration and
    # other demography flags are emitted properly. Fall back to a simple -I
    # block otherwise.
    if _demes_to_ms is not None:
        try:
            # demes.ms.to_ms in recent versions requires N0 as a keyword-only arg.
            # Try the most complete signature first, then fallback to older variants.
            if samples is not None:
                try:
                    return _demes_to_ms(graph, N0=N0, samples=samples)
                except TypeError:
                    # signature may not accept samples or N0 together; try other combos
                    pass
            try:
                return _demes_to_ms(graph, N0=N0)
            except TypeError:
                # Try older signatures: samples only, then bare graph
                if samples is not None:
                    try:
                        return _demes_to_ms(graph, samples=samples)
                    except TypeError:
                        return _demes_to_ms(graph)
                return _demes_to_ms(graph)
        except Exception:
            # fall through to simple fallback on any runtime error
            pass
    if samples is None:
        samples = []
    # create a trivial ms command fragment with no demographic events
    npop = max(1, len(samples))
    base = f"-I {npop} " + ' '.join(str(int(s)) for s in samples)
    return base
def strip_I_block(ms_cmd):
    base = f"-I {npop} " + ' '.join(str(int(s)) for s in samples)
    return base
def strip_I_block(ms_cmd):
    toks = ms_cmd.split()
    if '-I' not in toks:
        return ms_cmd
    i = toks.index('-I')
    if i+1 >= len(toks):
        return ms_cmd
    try:
        npop = int(toks[i+1])
    except ValueError:
        return ms_cmd
    j = i + 2 + npop
    return ' '.join(toks[:i] + toks[j:])

def strip_all_growth(ms_cmd):
    toks = ms_cmd.split()
    out = []
    i = 0
    while i < len(toks):
        if toks[i] in ('-G','-g'):
            i += 2
            continue
        out.append(toks[i]); i += 1
    return ' '.join(out)

def remap_indices_ms_1based(ms_cmd, mapping):
    toks = ms_cmd.split(); out = []; i = 0
    while i < len(toks):
        t = toks[i]
        if t == '-n' and i+2 < len(toks):
            try:
                idx = int(float(toks[i+1]))
                out += ['-n', str(mapping.get(idx, idx)), toks[i+2]]
                i += 3; continue
            except Exception:
                pass
        if t == '-en' and i+3 < len(toks):
            try:
                idx = int(float(toks[i+2]))
                out += ['-en', toks[i+1], str(mapping.get(idx, idx)), toks[i+3]]
                i += 4; continue
            except Exception:
                pass
        if t == '-es' and i+3 < len(toks):
            try:
                idx = int(float(toks[i+2]))
                out += ['-es', toks[i+1], str(mapping.get(idx, idx)), toks[i+3]]
                i += 4; continue
            except Exception:
                pass
        if t == '-m' and i+3 < len(toks):
            try:
                ii = int(float(toks[i+1])); jj = int(float(toks[i+2]))
                out += ['-m', str(mapping.get(ii, ii)), str(mapping.get(jj, jj)), toks[i+3]]
                i += 4; continue
            except Exception:
                pass
        if t == '-em' and i+4 < len(toks):
            try:
                ii = int(float(toks[i+2])); jj = int(float(toks[i+3]))
                out += ['-em', toks[i+1], str(mapping.get(ii, ii)), str(mapping.get(jj, jj)), toks[i+4]]
                i += 5; continue
            except Exception:
                pass
        if t == '-ej' and i+3 < len(toks):
            try:
                ii = int(float(toks[i+2])); jj = int(float(toks[i+3]))
                out += ['-ej', toks[i+1], str(mapping.get(ii, ii)), str(mapping.get(jj, jj))]
                i += 4; continue
            except Exception:
                pass
        out.append(t); i += 1
    return ' '.join(out)

def map_ej_to_ed(ms_cmd):
    toks = ms_cmd.split(); out = []; i = 0
    while i < len(toks):
        if toks[i] == '-ej' and i + 3 < len(toks):
            out += ['-ed', toks[i+1], toks[i+2], toks[i+3]]
            i += 4; continue
        out.append(toks[i]); i += 1
    return ' '.join(out)

def shift_to_0_based(ms_cmd):
    toks = ms_cmd.split(); out = []; i = 0
    while i < len(toks):
        t = toks[i]
        if t == '-en' and i + 3 < len(toks):
            try:
                out += ['-en', toks[i+1], str(int(float(toks[i+2]))-1), toks[i+3]]
                i += 4; continue
            except ValueError:
                pass
        if t == '-em' and i + 4 < len(toks):
            try:
                out += ['-em', toks[i+1], str(int(float(toks[i+2]))-1), str(int(float(toks[i+3]))-1), toks[i+4]]
                i += 5; continue
            except ValueError:
                pass
        # Attempt to stop any spawned watcher process (best-effort)
        try:
            if _TMP_WATCHER_PID:
                try:
                    os.kill(_TMP_WATCHER_PID, signal.SIGTERM)
                except Exception:
                    pass
        except Exception:
            pass
        if t == '-m' and i + 3 < len(toks):
            # migration rate: -m i j rate  (convert population indices to 0-based)
            try:
                out += ['-m', str(int(float(toks[i+1]))-1), str(int(float(toks[i+2]))-1), toks[i+3]]
                i += 4; continue
            except ValueError:
                pass
        if t == '-ed' and i + 3 < len(toks):
            try:
                out += ['-ed', toks[i+1], str(int(float(toks[i+2]))-1), str(int(float(toks[i+3]))-1)]
                i += 4; continue
            except ValueError:
                pass
        if t == '-es' and i + 3 < len(toks):
            try:
                out += ['-es', toks[i+1], str(int(float(toks[i+2]))-1), toks[i+3]]
                i += 4; continue
            except ValueError:
                pass
        if t == '-n' and i + 2 < len(toks):
            try:
                out += ['-n', str(int(float(toks[i+1]))-1), toks[i+2]]
                i += 3; continue
            except ValueError:
                pass
        out.append(t); i += 1
    return ' '.join(out)

def stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step=1.05):
    en = []
    graph_order = [d.name for d in graph.demes]
    name_to_idx = {name: i+1 for i, name in enumerate(graph_order)}
    for deme in graph.demes:
        i1 = name_to_idx[deme.name]
        for ep in deme.epochs:
            t0, t1 = ep.start_time, ep.end_time
            n0, n1 = ep.start_size, ep.end_size
            if ep.size_function == 'exponential' and t0 > t1:
                fold = n0 / n1 if n1 > 0 else float('inf')
                fold = max(fold, 1.0 + 1e-9)
                steps = max(1, math.ceil(math.log(fold) / math.log(max_fold_per_step)))
                for s in range(1, steps+1):
                    frac = s/steps
                    n_t = n0 * (n1/n0)**frac
                    t_gen = t0 - frac*(t0 - t1)
                    t_ms = t_gen / (4.0 * N0)
                    en.append((t_ms, i1, n_t / N0))
    en.sort(key=lambda x: -x[0])
    return en

def sanitize_ms_demography(ms_demog_str):
    toks = ms_demog_str.split(); out = []; i = 0
    # include common ms-style -e* flags with their expected argument counts so
    # sanitize preserves their arguments instead of dropping them. Missing
    # entries here caused flags like '-es' to be emitted without their args,
    # which breaks downstream MS-like parsers (see CmdLineParseException).
    # Flags with fixed argument counts
    arg_counts = {'-n':2,'-en':3,'-m':3,'-em':4,'-ej':3, '-ed':3, '-es':3, '-eg':4}
    def is_flag(tok):
        return tok.startswith('-')

    while i < len(toks):
        t = toks[i]
        if t in arg_counts:
            need = arg_counts[t]
            # Special-case '-n': it may be followed by multiple size values for the same index
            # (demes.ms.to_ms can emit multiple sizes in sequence). Convert sequences like
            # '-n IDX SIZE1 SIZE2' into '-n IDX SIZE1 -n IDX SIZE2' so downstream ms-like
            # parsers don't see orphan numeric tokens.
            if t == '-n':
                # ensure index token exists
                if i + 1 >= len(toks):
                    i += 1
                    continue
                idx_tok = toks[i+1]
                # collect following non-flag tokens as sizes
                j = i + 2
                sizes = []
                while j < len(toks) and not is_flag(toks[j]):
                    sizes.append(toks[j]); j += 1
                # if no sizes found, just emit the original -n and index (defensive)
                if not sizes:
                    out.extend(toks[i:min(len(toks), i+1+need)])
                    i = min(len(toks), i+1+need)
                    continue
                # emit one '-n idx size' for each size token
                for sz in sizes:
                    out.extend(['-n', idx_tok, sz])
                i = j
                continue
            # default handling for flags with fixed argument counts
            take_end = min(len(toks), i + 1 + need)
            out.extend(toks[i:take_end])
            i = take_end
            continue
        if is_flag(t):
            # Unknown flag: keep the flag and also keep following tokens until next flag
            out.append(t)
            j = i + 1
            while j < len(toks) and not is_flag(toks[j]):
                out.append(toks[j]); j += 1
            i = j
            continue
        # bare token (shouldn't usually happen) - keep it
        out.append(t); i += 1
    return ' '.join(out)


def _primary_simulation_is_neutral(meta):
    """Return True when the primary simulation should be considered neutral.

    Rules (mirrors the inline logic used elsewhere):
    - Engines 'ms', 'scrm', 'msprime' are considered neutral-only.
    - For other engines, absence of selection metadata (sel_args, sweep_pos, sel_2Ns)
      means the primary run is neutral.
    """
    try:
        if not isinstance(meta, dict):
            return False
        eng0 = meta.get('engine')
        if eng0 in ('ms', 'scrm', 'msprime'):
            return True
        sel_present = bool(meta.get('sel_args')) or bool(meta.get('sweep_pos')) or bool(meta.get('sel_2Ns'))
        return not sel_present
    except Exception:
        return False

def _make_affinity_preexec(cpu_list):
    """Return a preexec_fn that sets CPU affinity for the child process (best-effort)."""
    def _set_affinity():
        try:
            os.sched_setaffinity(0, cpu_list)
        except Exception:
            pass
    return _set_affinity


# ----- child process tracking for safer cleanup -----
# Keep a thread-safe set of child PIDs started by this process so we can
# attempt to terminate them if the parent receives SIGINT/SIGTERM.
_child_pids_lock = threading.Lock()
_child_pids = set()

def _register_child_pid(pid: int):
    try:
        with _child_pids_lock:
            _child_pids.add(int(pid))
    except Exception:
        pass

def _unregister_child_pid(pid: int):
    try:
        with _child_pids_lock:
            _child_pids.discard(int(pid))
    except Exception:
        pass

def _terminate_children(timeout=2.0):
    """Best-effort: send SIGTERM to known child PIDs, then SIGKILL if needed."""
    try:
        with _child_pids_lock:
            pids = list(_child_pids)
    except Exception:
        pids = []
    if not pids:
        return
    # Attempt graceful termination
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except Exception:
            pass
        # Also attempt to terminate the child's process group if possible
        try:
            pg = os.getpgid(pid)
            if pg:
                try:
                    os.kill(-pg, signal.SIGTERM)
                except Exception:
                    pass
        except Exception:
            pass
    # brief wait
    try:
        end = time.time() + float(timeout)
        while time.time() < end:
            alive = []
            for pid in pids:
                try:
                    # 0 signal checks existence
                    os.kill(pid, 0)
                    alive.append(pid)
                except Exception:
                    pass
            if not alive:
                break
            time.sleep(0.05)
            pids = alive
    except Exception:
        pass
    # Force kill any remaining
    for pid in pids:
        try:
            os.kill(pid, signal.SIGKILL)
        except Exception:
            pass
        try:
            pg = os.getpgid(pid)
            if pg:
                try:
                    os.kill(-pg, signal.SIGKILL)
                except Exception:
                    pass
        except Exception:
            pass


def _run_and_register(*popen_args, **popen_kwargs):
    """Lightweight wrapper around subprocess.Popen that registers child PID.

    Mirrors the common subset of subprocess.run(...) used in this file and
    returns a subprocess.CompletedProcess-like object with attributes
    stdout, stderr, returncode. If 'check=True' is passed and the returncode
    is non-zero, raise subprocess.CalledProcessError like subprocess.run.
    """
    # Map signature to Popen constructor
    stdin = popen_kwargs.pop('stdin', None)
    capture_output = popen_kwargs.pop('capture_output', False)
    text_mode = popen_kwargs.pop('text', False)
    check = popen_kwargs.pop('check', False)
    env = popen_kwargs.pop('env', None)
    pre = popen_kwargs.pop('preexec_fn', None)
    # Build Popen args
    try:
        proc = subprocess.Popen(*popen_args, stdin=stdin, stdout=subprocess.PIPE if capture_output else None,
                                stderr=subprocess.PIPE if capture_output else None, env=env, preexec_fn=pre, text=text_mode)
    except TypeError:
        # fallback for older Python versions without text in Popen
        proc = subprocess.Popen(*popen_args, stdin=stdin, stdout=subprocess.PIPE if capture_output else None,
                                stderr=subprocess.PIPE if capture_output else None, env=env, preexec_fn=pre)
    # register pid
    try:
        _register_child_pid(proc.pid)
    except Exception:
        pass
    try:
        out, err = proc.communicate()
        rc = proc.returncode
    finally:
        # unregister regardless of outcome
        try:
            _unregister_child_pid(proc.pid)
        except Exception:
            pass
    class _CP:
        def __init__(self, stdout, stderr, returncode):
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode
    cp = _CP(out, err, rc)
    if check and cp.returncode and cp.returncode != 0:
        raise subprocess.CalledProcessError(cp.returncode, popen_args[0] if popen_args else None, output=cp.stdout, stderr=cp.stderr)
    return cp


# ---------- local temp-dir helpers ----------
# Use a per-run temporary directory placed next to this script so callers
# can inspect files if needed and so cleanup can remove the entire folder.
_LOCAL_TMP_ROOT = None
_TMP_WATCHER_PID = None
# Lock to make local tmpdir creation thread-safe (avoid races when multiple
# workers try to create/cleanup the same repo-local tmpdir simultaneously).
try:
    _LOCAL_TMP_LOCK = threading.Lock()
except Exception:
    _LOCAL_TMP_LOCK = None

# Module-level friendly message populated when RAM threshold is exceeded
# so signal handlers can access it even if running in a different frame.
_LAST_RAM_EXCEEDED_MESSAGE = None


def _is_process_alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except Exception:
        return False


def _spawn_tmp_watcher(tmp_path, parent_pid):
    """Spawn a detached watcher that removes tmp_path if parent_pid dies.

    This forks a child process that detaches and monitors the parent.
    It's best-effort: it will remove the directory unless
    SIMULATOR_PRESERVE_TMP is set.
    """
    global _TMP_WATCHER_PID
    if not tmp_path or not os.path.isdir(tmp_path):
        return
    # If a watcher is already recorded and alive, avoid spawning another one.
    try:
        global _TMP_WATCHER_PID
        if _TMP_WATCHER_PID and _is_process_alive(_TMP_WATCHER_PID):
            return
    except Exception:
        pass
    try:
        pid = os.fork()
    except OSError:
        return
    if pid > 0:
        # parent: record watcher pid
        _TMP_WATCHER_PID = pid
        return
    # child: detach
    try:
        os.setsid()
    except Exception:
        pass
    try:
        # loop until parent disappears
        while True:
            if not _is_process_alive(parent_pid) or os.getppid() == 1:
                break
            time.sleep(0.5)
    except Exception:
        pass
    # attempt removal: only remove the specific run tmp_path. Do NOT remove
    # the shared parent '.tmp' directory — that can race with other runs and
    # cause unrelated run folders to be deleted by older watchers.
    try:
        if os.path.isdir(tmp_path) and not os.environ.get('SIMULATOR_PRESERVE_TMP'):
            try:
                shutil.rmtree(tmp_path, ignore_errors=True)
            except Exception:
                pass
    except Exception:
        pass
    try:
        os._exit(0)
    except Exception:
        os._exit(0)

def _ensure_local_tmpdir():
    """Create and return a local tmp directory inside the simulator script dir.

    The directory is created lazily on first call and returned for reuse. Files
    should be created inside this directory. Cleanup is performed by
    _cleanup_local_tmpdir().
    """
    global _LOCAL_TMP_ROOT
    # If a worker-specific tmpdir is requested via env, prefer and return it
    worker_tmp_env = os.environ.get('SIMULATOR_WORKER_TMP')
    if worker_tmp_env:
        try:
            os.makedirs(worker_tmp_env, exist_ok=True)
            return worker_tmp_env
        except Exception:
            # fallback to normal behaviour
            pass
    # Fast path + guard: if already created, return it. Use a lock to avoid
    # races between concurrent callers creating the same repo-local tmpdir.
    try:
        if _LOCAL_TMP_LOCK:
            _LOCAL_TMP_LOCK.acquire()
        if _LOCAL_TMP_ROOT and os.path.isdir(_LOCAL_TMP_ROOT):
            return _LOCAL_TMP_ROOT
    finally:
        try:
            if _LOCAL_TMP_LOCK:
                _LOCAL_TMP_LOCK.release()
        except Exception:
            pass
    try:
        base = os.path.dirname(__file__) or os.getcwd()
    except Exception:
        base = os.getcwd()
    # Repo-local parent (hidden) named '.tmp' to group all runs; per-run subfolders inside it
    repo_parent = os.path.join(base, '.tmp')
    run_name = f"simulator_run_{os.getpid()}_{int(time.time())}_{random.randrange(10**6)}"
    path = os.path.join(repo_parent, run_name)
    # Prefer using system tmpdir if it's a tmpfs (in-RAM) for speed.
    try:
        system_tmp = tempfile.gettempdir()
    except Exception:
        system_tmp = None

    def _is_tmpfs(p):
        try:
            # find best mount match from /proc/mounts
            mounts = []
            with open('/proc/mounts', 'rt') as m:
                for ln in m:
                    parts = ln.split()
                    if len(parts) >= 3:
                        mounts.append((parts[1], parts[2]))
            # choose longest mountpoint prefix that matches p
            best = ('', '')
            for mp, fstype in mounts:
                if p.startswith(mp) and len(mp) > len(best[0]):
                    best = (mp, fstype)
            return best[1] in ('tmpfs', 'ramfs')
        except Exception:
            return False

    try:
        # Allow overriding via environment: force repo-local tmp if set
        force_repo = os.environ.get('SIMULATOR_USE_REPO_TMP', '')
        force_system = os.environ.get('SIMULATOR_USE_SYSTEM_TMP', '')
        if force_repo:
            use_system = False
        elif force_system:
            use_system = True
        else:
            use_system = bool(system_tmp and _is_tmpfs(system_tmp))

        if use_system and system_tmp:
            # create a hidden subdir under system tmp (tmpfs) for speed
            name2 = f"simulator_tmp_{os.getpid()}_{int(time.time())}_{random.randrange(10**6)}"
            path2 = os.path.join(system_tmp, name2)
            try:
                os.makedirs(path2, exist_ok=True)
                _LOCAL_TMP_ROOT = path2
                # debug print
                if os.environ.get('SIMULATOR_DEBUG_TMP'):
                    sys.stderr.write(f"DEBUG: using system tmpdir: {_LOCAL_TMP_ROOT}\n")
                    sys.stderr.flush()
                try:
                    _spawn_tmp_watcher(_LOCAL_TMP_ROOT, os.getpid())
                except Exception:
                    pass
                return _LOCAL_TMP_ROOT
            except Exception:
                pass

        # fallback to repo-local hidden '.tmp' parent with per-run subfolder
        try:
            # Use a few retries to tolerate transient filesystem races (e.g. NFS or concurrent processes)
            created = False
            for attempt in range(3):
                try:
                    os.makedirs(repo_parent, exist_ok=True)
                    os.makedirs(path, exist_ok=True)
                    _LOCAL_TMP_ROOT = path
                    created = True
                    break
                except Exception:
                    time.sleep(0.05 + 0.05 * attempt)
            if not created:
                raise
            if os.environ.get('SIMULATOR_DEBUG_TMP'):
                sys.stderr.write(f"DEBUG: using repo-local tmpdir: {_LOCAL_TMP_ROOT}\n")
                sys.stderr.flush()
            try:
                _spawn_tmp_watcher(_LOCAL_TMP_ROOT, os.getpid())
            except Exception:
                pass
        except Exception:
            # final fallback to system tempdir
            _LOCAL_TMP_ROOT = tempfile.gettempdir()
    except Exception:
        _LOCAL_TMP_ROOT = tempfile.gettempdir()
    return _LOCAL_TMP_ROOT


def _cleanup_local_tmpdir():
    """Remove the local tmp directory and its contents if it exists."""
    global _LOCAL_TMP_ROOT
    if not _LOCAL_TMP_ROOT:
        return
    try:
        # only remove directories we created (hidden prefix) to be safe
        base_name = os.path.basename(_LOCAL_TMP_ROOT or '')
        parent = os.path.dirname(_LOCAL_TMP_ROOT or '')
        # If this is a per-run folder inside repo '.tmp', remove run folder and remove parent
        if parent and os.path.basename(parent) == '.tmp' and os.path.isdir(_LOCAL_TMP_ROOT):
            # allow preserving tmp for inspection when requested
            if not os.environ.get('SIMULATOR_PRESERVE_TMP'):
                try:
                    # remove the entire run folder
                    shutil.rmtree(_LOCAL_TMP_ROOT, ignore_errors=True)
                except Exception:
                    pass
                if os.environ.get('SIMULATOR_DEBUG_TMP'):
                    try:
                        sys.stderr.write(f"DEBUG: removed repo-local tmpdir: {_LOCAL_TMP_ROOT}\n")
                        sys.stderr.flush()
                    except Exception:
                        pass
                # Do NOT remove the repo-level '.tmp' parent directory here.
                # Removing the shared parent can race with other concurrent runs
                # and cause unrelated run folders to be deleted. Keep only the
                # per-run folder removal above.
                if os.environ.get('SIMULATOR_DEBUG_TMP'):
                    try:
                        sys.stderr.write(f"DEBUG: preserved .tmp parent directory: {parent}\n")
                        sys.stderr.flush()
                    except Exception:
                        pass
        # keep previous safety: support older naming convention
        elif base_name.startswith('.simulator_tmp_') and os.path.isdir(_LOCAL_TMP_ROOT):
            shutil.rmtree(_LOCAL_TMP_ROOT, ignore_errors=True)
            if os.environ.get('SIMULATOR_DEBUG_TMP'):
                try:
                    sys.stderr.write(f"DEBUG: removed old-style tmpdir: {_LOCAL_TMP_ROOT}\n")
                    sys.stderr.flush()
                except Exception:
                    pass
        # Also support system tmp names created by this helper
        elif os.path.basename(_LOCAL_TMP_ROOT).startswith('simulator_tmp_') and os.path.isdir(_LOCAL_TMP_ROOT):
            shutil.rmtree(_LOCAL_TMP_ROOT, ignore_errors=True)
            if os.environ.get('SIMULATOR_DEBUG_TMP'):
                try:
                    sys.stderr.write(f"DEBUG: removed system tmpdir: {_LOCAL_TMP_ROOT}\n")
                    sys.stderr.flush()
                except Exception:
                    pass
    except Exception:
        pass
    finally:
        # ensure watcher pid reference is cleared on cleanup
        try:
            global _TMP_WATCHER_PID
            _TMP_WATCHER_PID = None
        except Exception:
            pass
        _LOCAL_TMP_ROOT = None


# Register cleanup handlers so tmp dir is removed on process exit or common termination signals
def _register_cleanup_handlers():
    try:
        atexit.register(_cleanup_local_tmpdir)
    except Exception:
        pass
    # Module-level place to store a friendly RAM-exceeded message populated by the monitor
    try:
        _LAST_RAM_EXCEEDED_MESSAGE
    except NameError:
        _LAST_RAM_EXCEEDED_MESSAGE = None

    def _signal_handler(signum, frame):
        # If the RAM monitor prepared a friendly message, print it first
        try:
            msg = globals().get('_LAST_RAM_EXCEEDED_MESSAGE')
            if msg:
                try:
                    sys.stderr.write(str(msg).rstrip() + "\n")
                    try:
                        sys.stderr.flush()
                    except Exception:
                        pass
                except Exception:
                    pass
        except Exception:
            pass
        # Cleanup local temporary directory
        try:
            _cleanup_local_tmpdir()
        except Exception:
            pass
        # Revert handler to default first to avoid recursive invocation
        try:
            signal.signal(signum, signal.SIG_DFL)
        except Exception:
            pass
        # Attempt to terminate any known child PIDs
        try:
            _terminate_children()
        except Exception:
            pass
        # As a best-effort, try killing the whole process group so
        # subprocesses sharing our group receive the signal too. Using
        # the default handler above prevents our custom handler re-entry.
        try:
            pgid = os.getpgid(os.getpid())
            os.kill(-pgid, signum)
        except Exception:
            pass
    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        try:
            signal.signal(sig, _signal_handler)
        except Exception:
            pass

# call registration at import time
_register_cleanup_handlers()

# ---------- conversion helpers ----------

def ms_like_to_vcf(ms_text: str, length: int, chrom: str = 'chr1', ploidy: int = 1) -> str:
    """Convert a single-replicate ms-like block to a simple VCF string.

    Expects a block starting with 'segsites:' and a following 'positions:' line, then haplotypes.
    """
    if not ms_text:
        # minimal empty VCF
        header = [
            '##fileformat=VCFv4.2',
            '##source=stdpopsim2ms.py',
            f'##contig=<ID={chrom},length={int(length) if length else 0}>',
            '##INFO=<ID=SIM,Number=0,Type=Flag,Description="Placeholder">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        ]
        return '\n'.join(header) + '\n'
    lines = [l.rstrip('\n') for l in ms_text.splitlines()]
    segsites = 0
    pos_vals = []
    hap_lines = []
    i = 0
    # find first segsites (allow leading whitespace)
    while i < len(lines):
        ln = lines[i].lstrip()
        if ln.startswith('segsites:'):
            try:
                segsites = int(ln.split()[1])
            except Exception:
                segsites = 0
            # positions line should follow within a few lines
            j = i + 1
            pos_line_idx = None
            while j < len(lines) and j <= i + 10:
                if lines[j].lstrip().startswith('positions:'):
                    pos_line_idx = j
                    break
                if lines[j].lstrip().startswith('segsites:'):
                    break
                j += 1
            if pos_line_idx is None:
                # consider no variants
                segsites = 0
                pos_vals = []
                hap_lines = []
                break
            try:
                pos_tokens = lines[pos_line_idx].split(':', 1)[1].strip().split()
            except Exception:
                pos_tokens = []
            try:
                pos_vals = [float(p) for p in pos_tokens]
            except Exception:
                pos_vals = []
            # collect haplotypes after positions
            k = pos_line_idx + 1
            raw_hap_lines = []
            while k < len(lines):
                s = lines[k].rstrip('\n')
                s_strip = s.strip()
                if not s_strip or s_strip.lstrip().startswith('//') or s_strip.lstrip().startswith('segsites:') or s_strip.lstrip().startswith('positions:'):
                    break
                raw_hap_lines.append(s)
                k += 1

            # Normalize collected lines into pure 0/1 haplotype strings of length `segsites`.
            # Be permissive: accept lines that contain spaces, other separators, or multiple
            # concatenated haplotypes on the same physical line. This handles msms output
            # variants that may include spacing or markers when -SI/-Smark/-SFC are used.
            hap_lines = []
            if segsites > 0:
                for line in raw_hap_lines:
                    ls = line.strip()
                    # Fast-path: line already a contiguous 0/1 string
                    if ls and set(ls) <= {'0', '1'} and len(ls) == segsites:
                        hap_lines.append(ls)
                        continue

                    # Remove any non 0/1 characters to collapse tokens like '0 1 0' -> '010'
                    only01 = ''.join(ch for ch in line if ch in ('0', '1'))
                    if len(only01) == segsites:
                        hap_lines.append(only01)
                        continue

                    # If the line contains multiple concatenated haplotypes (e.g. '010...010...'),
                    # split into chunks of segsites
                    if len(only01) > segsites and (len(only01) % segsites) == 0:
                        for i0 in range(0, len(only01), segsites):
                            chunk = only01[i0:i0+segsites]
                            if len(chunk) == segsites:
                                hap_lines.append(chunk)
                        continue

                    # Finally, consider whitespace-separated tokens that are individual alleles
                    toks = [t for t in ls.split() if t in ('0', '1')]
                    if len(toks) == segsites:
                        hap_lines.append(''.join(toks))
                        continue

                # Fallback: if parsing produced very few haplines, attempt to recover
                # by concatenating all 0/1 characters present in the raw haplotype
                # region and splitting into fixed-length chunks of `segsites`.
                if len(hap_lines) < 2:
                    combined = ''.join(ch for ch in ''.join(raw_hap_lines) if ch in ('0', '1'))
                    if segsites > 0 and len(combined) >= segsites and (len(combined) % segsites) == 0:
                        recovered = [combined[i0:i0+segsites] for i0 in range(0, len(combined), segsites)]
                        if len(recovered) > len(hap_lines):
                            hap_lines = recovered

                # Recovery: some ms/msms variants emit haplotype lines that are
                # uniformly shorter than the declared segsites (for example
                # segsites-1). In that case we can pad each cleaned hap string
                # with '0' (ancestral) at the end to restore the expected
                # length and proceed with conversion.
                if raw_hap_lines and segsites > 0:
                    cleaned = [''.join(ch for ch in l if ch in ('0', '1')) for l in raw_hap_lines]
                    # If a large fraction of cleaned lines already match segsites,
                    # prefer those as the haplotypes (handles cases where only a
                    # subset matched earlier fast-paths).
                    matching = [c for c in cleaned if len(c) == segsites]
                    if len(matching) >= max(2, int(0.5 * len(cleaned))):
                        hap_lines = matching
                    else:
                        # Choose the most common cleaned length and, if it
                        # represents a majority and is within a small delta of
                        # segsites, pad those lines to recover haplotypes. Also
                        # handle the mixed case where most lines are short and a
                        # few are already the correct length: pad the shorter
                        # ones and include the correct ones in original order.
                        from collections import Counter
                        cnt = Counter(len(c) for c in cleaned)
                        most_len, most_count = cnt.most_common(1)[0]
                        # Mixed-length case: some lines already match segsites,
                        # others are uniformly a bit shorter -> pad the short
                        # ones and keep the correct ones so total count remains.
                        correct_count = sum(1 for c in cleaned if len(c) == segsites)
                        if correct_count >= 1 and most_len < segsites and (segsites - most_len) <= 2 and (most_count + correct_count) == len(cleaned):
                            # pad short ones and include correct-length ones in original order
                            hap_lines = [ (c if len(c) == segsites else c + ('0' * (segsites - len(c)))) for c in cleaned ]
                        elif most_count >= max(2, int(0.5 * len(cleaned))) and segsites - most_len <= 2:
                            # pad all cleaned lines that have the most common length
                            hap_lines = [c + ('0' * (segsites - len(c))) for c in cleaned if len(c) == most_len]
                        else:
                            # Fallback: if all lines are same length and close to segsites,
                            # pad them all (handles uniform shortness)
                            lens = set(len(c) for c in cleaned)
                            if len(lens) == 1:
                                only_len = next(iter(lens))
                                if only_len == segsites:
                                    hap_lines = cleaned
                                elif only_len < segsites and segsites - only_len <= 2:
                                    hap_lines = [c + ('0' * (segsites - len(c))) for c in cleaned]
            break
        i += 1
    if segsites == 0 or not pos_vals or not hap_lines:
        # header-only VCF for empty replicate
        header = [
            '##fileformat=VCFv4.2',
            '##source=stdpopsim2ms.py',
            f'##contig=<ID={chrom},length={int(length) if length else 0}>',
            '##INFO=<ID=SIM,Number=0,Type=Flag,Description="Placeholder">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        ]
        return '\n'.join(header) + '\n'
    # normalize positions to integer coordinates
    used = set()
    pos_ints = []
    L = int(length) if length else 0
    for f in pos_vals:
        try:
            p = int(round(float(f) * float(L)))
        except Exception:
            p = 1
        if p < 1:
            p = 1
        if L and p > L:
            p = L
        while p in used and (not L or p < L):
            p += 1
        used.add(p)
        pos_ints.append(p)
    n_hap = len(hap_lines)
    if ploidy < 1:
        ploidy = 1
    if n_hap % ploidy != 0:
        raise SystemExit(f'ERROR: Total haplotypes {n_hap} not divisible by ploidy {ploidy} for VCF conversion.')
    n_ind = n_hap // ploidy
    header = [
        '##fileformat=VCFv4.2',
        '##source=stdpopsim2ms.py',
        f'##contig=<ID={chrom},length={L}>',
        '##INFO=<ID=SIM,Number=0,Type=Flag,Description="Placeholder">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    ]
    sample_names = [f'Ind{i+1}' for i in range(n_ind)]
    header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sample_names))
    body = []
    for vidx, p in enumerate(pos_ints, start=1):
        ref = 'A'; alt = 'T'
        genos = []
        for i in range(0, n_hap, ploidy):
            alleles = ['1' if hap_lines[i + a][vidx - 1] == '1' else '0' for a in range(ploidy)]
            genos.append('|'.join(alleles))
        body.append(f'{chrom}\t{p}\tsnp{vidx}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t' + '\t'.join(genos))
    return '\n'.join(header + body) + '\n'

def _msprime_worker_run(payload, env=None):
    """Run a small msprime worker in-process using stdpopsim.

    Accepts a payload dict similar to the external msprime_worker.py JSON and
    returns (stdout_text, stderr_text, returncode).
    """
    out_parts = []
    err_parts = []
    exit_code = 0
    # If an env mapping was provided, temporarily apply it so msprime and
    # underlying native libraries see the intended OpenMP/BLAS thread caps.
    # Restore the original environment after the run to avoid global side-effects.
    old_env_vals = {}
    applied_env = False
    if env:
        try:
            for k, v in env.items():
                # Only set string values in os.environ
                old_env_vals[k] = os.environ.get(k)
                os.environ[k] = str(v)
            applied_env = True
        except Exception:
            # If we fail to set env, continue without failing the worker
            applied_env = False
    # Optional per-worker peak sampler: when SIMULATOR_PEAK_DIR is set in the
    # worker env (or in the process env), spawn a small background thread that
    # samples this process RSS and writes peak_<pid>.txt in that directory.
    peak_thread = None
    peak_stop_event = None
    try:
        peak_dir = None
        try:
            if env and isinstance(env, dict) and env.get('SIMULATOR_PEAK_DIR'):
                peak_dir = env.get('SIMULATOR_PEAK_DIR')
            else:
                peak_dir = os.environ.get('SIMULATOR_PEAK_DIR')
        except Exception:
            peak_dir = None

        peak_interval = None
        try:
            if env and isinstance(env, dict) and env.get('SIMULATOR_PEAK_INTERVAL'):
                peak_interval = float(env.get('SIMULATOR_PEAK_INTERVAL'))
            else:
                peak_interval = float(os.environ.get('SIMULATOR_PEAK_INTERVAL', 0.02))
        except Exception:
            peak_interval = 0.02

        if peak_dir:
            try:
                os.makedirs(peak_dir, exist_ok=True)
            except Exception:
                pass

            def _sampler_current_process(peak_dir_local, stop_event, interval):
                try:
                    import psutil
                except Exception:
                    return
                try:
                    p = psutil.Process(os.getpid())
                except Exception:
                    return
                peak = 0
                try:
                    while not stop_event.is_set():
                        try:
                            rss = int(getattr(p.memory_info(), 'rss', 0) or 0)
                            if rss > peak:
                                peak = rss
                        except Exception:
                            break
                        stop_event.wait(interval)
                    # final sample
                    try:
                        rss = int(getattr(p.memory_info(), 'rss', 0) or 0)
                        if rss > peak:
                            peak = rss
                    except Exception:
                        pass
                finally:
                    try:
                        path = os.path.join(peak_dir_local, f"peak_{os.getpid()}.txt")
                        with open(path, 'w') as pf:
                            pf.write(str(peak))
                    except Exception:
                        pass

            peak_stop_event = threading.Event()
            peak_thread = threading.Thread(target=_sampler_current_process, args=(peak_dir, peak_stop_event, peak_interval), daemon=True)
            try:
                peak_thread.start()
            except Exception:
                peak_thread = None
                peak_stop_event = None

    except Exception:
        # If sampler setup fails, continue without breaking the worker
        peak_thread = None
        peak_stop_event = None

    try:
        species = payload.get('species')
        model_id = payload.get('model_id') or payload.get('model')
        chromosome = payload.get('chromosome')
        length = int(payload.get('length')) if payload.get('length') is not None else None
        pop_order = payload.get('pop_order') or []
        counts = payload.get('counts') or payload.get('counts_disc') or []
        ploidy = int(payload.get('ploidy') or 1)
        rep_count = int(payload.get('rep_count') or 1)
        seed = payload.get('seed')

        if species is None or model_id is None:
            return ('', '# ERROR: payload must include species and model_id\n', 4)

        try:
            sp = sps.get_species(species)
        except Exception as e:
            return ('', f"# ERROR: failed to resolve species '{species}': {e}\n", 5)

        try:
            model_obj = sp.get_demographic_model(model_id)
        except Exception as e:
            return ('', f"# ERROR: failed to get demographic model '{model_id}': {e}\n", 6)

        chosen_contig = None
        if chromosome:
            try:
                if length:
                    try:
                        chosen_contig = sp.get_contig(chromosome, mutation_rate=getattr(model_obj, 'mutation_rate', None), left=0, right=max(0, int(length) - 1))
                    except Exception:
                        chosen_contig = sp.get_contig(chromosome)
                else:
                    chosen_contig = sp.get_contig(chromosome)
            except Exception as e:
                return ('', f"# ERROR: failed to resolve contig '{chromosome}': {e}\n", 7)
        else:
            try:
                genome = getattr(sp, 'genome', None)
                contigs = getattr(genome, 'contigs', []) if genome is not None else []
                if contigs:
                    chosen_contig = contigs[0]
                else:
                    chosen_contig = None
            except Exception:
                chosen_contig = None

        if chosen_contig is None:
            return ('', '# ERROR: could not determine contig for msprime worker; provide --chromosome.\n', 8)

        
        # Resolve a robust mutation rate for msprime mutations if needed
        mu_local = None
        try:
            # 1) contig-specific mutation rate
            if hasattr(chosen_contig, 'mutation_rate') and chosen_contig.mutation_rate is not None:
                mu_local = float(chosen_contig.mutation_rate)
        except Exception:
            mu_local = None
        if mu_local is None:
            try:
                # 2) model-level mutation rate
                model_mu = getattr(model_obj, 'mutation_rate', None)
                if model_mu is not None:
                    mu_local = float(model_mu)
            except Exception:
                mu_local = None
        if mu_local is None:
            try:
                # 3) species genome default
                genome = getattr(sp, 'genome', None)
                for attr in ('per_site_mutation_rate','mutation_rate','mu'):
                    if genome is not None and hasattr(genome, attr):
                        val = getattr(genome, attr)
                        if val is not None:
                            mu_local = float(val)
                            break
            except Exception:
                mu_local = None
        if mu_local is None:
            mu_local = 1e-8
        samples_per_pop = {}
        try:
            for i, pop in enumerate(pop_order):
                hapc = int(counts[i]) if i < len(counts) else 0
                samples_per_pop[pop] = int(hapc // ploidy)
        except Exception:
            samples_per_pop = {p: 0 for p in pop_order}

        try:
            if seed is not None:
                import random as _rnd
                _rnd.seed(int(seed))
                try:
                    import numpy as _np_local
                    _np_local.random.seed(int(seed))
                except Exception:
                    pass
        except Exception:
            pass

        try:
            eng = sps.get_engine('msprime')
        except Exception as e:
            return ('', f"# ERROR: failed to get msprime engine: {e}\n", 9)

        # local converter from tree sequence to ms-like (simple)
        
        def _ts_to_ms_like(ts, length_val, pl=1, suppress_sfs=False):
            # Stream haplotypes and compute per-replicate SFS without building a genotype matrix.
            if ts is None:
                return 'segsites: 0\n'
            try:
                n_hap = int(ts.num_samples)
            except Exception:
                n_hap = 0
            if n_hap <= 0:
                return 'segsites: 0\n'
            # Prepare per-haplotype buffers and per-site derived counts
            hap_buffers = [bytearray() for _ in range(n_hap)]
            counts = []
            pos_vals = []
            try:
                it = ts.variants()
            except Exception:
                it = []
            for v in it:
                try:
                    g = v.genotypes
                    # append one column to each hap buffer
                    for i in range(n_hap):
                        hap_buffers[i].append(49 if g[i] == 1 else 48)  # '1' or '0'
                    # store per-site derived count
                    try:
                        counts.append(int(g.sum()))
                    except Exception:
                        # fallback if sum not available
                        c = 0
                        for i in range(n_hap):
                            if g[i] == 1:
                                c += 1
                        counts.append(int(c))
                    # positions: normalized to [0,1] if length provided
                    p = getattr(v, 'position', 0.0)
                    if length_val:
                        try:
                            denom = float(length_val)
                        except Exception:
                            denom = None
                        if denom and denom > 0:
                            pos_vals.append(float(p) / denom)
                        else:
                            pos_vals.append(float(p))
                    else:
                        pos_vals.append(float(p))
                except Exception:
                    continue
            segsites = len(counts)
            if segsites == 0:
                return 'segsites: 0\n'
            # Build ms-like text
            out_lines = []
            out_lines.append(f"segsites: {segsites}")
            out_lines.append('positions: ' + ' '.join(str(p) for p in pos_vals))
            out_lines.extend(b.decode('ascii') for b in hap_buffers)
            # Optionally append a per-replicate SFS line to make downstream SFS robust
            if not suppress_sfs:
                try:
                    import numpy as _np_local
                    binc = _np_local.bincount(_np_local.asarray(counts, dtype=_np_local.int32), minlength=n_hap+1)
                    sfs_vec = binc[1:n_hap]
                    # store raw counts (un-normalized) so the writer can choose normalization
                    out_lines.append('#SFS_PER_REP\t' + '\t'.join(str(int(x)) for x in sfs_vec.tolist()))
                except Exception:
                    pass
            return '\n'.join(out_lines) + '\n'

        for ridx in range(rep_count):
            try:
                # Prefer passing an explicit seed to the engine when provided
                if seed is not None:
                    try:
                        ts = eng.simulate(model_obj, chosen_contig, samples_per_pop, seed=int(seed))
                    except TypeError:
                        # Some engine wrappers do not accept a seed kwarg; fall back
                        ts = eng.simulate(model_obj, chosen_contig, samples_per_pop)
                else:
                    ts = eng.simulate(model_obj, chosen_contig, samples_per_pop)

                ts = _ensure_mutations_on_ts(ts)
                try:
                    import msprime as _msp
                    _num_sites = int(getattr(ts, 'num_sites', 0))
                    if _num_sites == 0:
                        _dg = bool(getattr(ts, 'discrete_genome', False)) if hasattr(ts, 'discrete_genome') else False
                        _mu = mu_local if isinstance(mu_local, (int, float)) and mu_local > 0 else 1e-8
                        ts = _msp.sim_mutations(ts, rate=_mu, discrete_genome=_dg)
                except Exception:
                    pass

                # Prefer writing VCF directly from the tree sequence when available
                try:
                    only_sfs_worker = bool(payload.get('only_sfs', False)) if isinstance(payload, dict) else False
                    if only_sfs_worker:
                        # Compute SFS but avoid producing the main ms/vcf text to save memory.
                        try:
                            sfs_af = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
                            sfs_af_line = '#SFS_AF\t' + '\t'.join(str(float(x)) for x in (sfs_af.tolist() if hasattr(sfs_af, 'tolist') else list(sfs_af))) + '\n'
                            err_parts.append(sfs_af_line)
                        except Exception:
                            # non-fatal: fall back to emitting zero SFS marker
                            try:
                                err_parts.append('#SFS_AF\t' + '\t'.join('0' for _ in range(max(0, int(getattr(ts, "num_samples", 0)) - 0))))
                            except Exception:
                                pass
                    else:
                        buf = io.StringIO()
                        if hasattr(ts, 'write_vcf'):
                            cid = str(chromosome) if chromosome is not None else '1'
                            try:
                                ts.write_vcf(buf, contig_id=cid)
                            except TypeError:
                                ts.write_vcf(buf)
                            vcf_text = buf.getvalue()
                            out_parts.append(vcf_text)
                        else:
                            suppress_flag = bool(payload.get('suppress_worker_sfs', False)) if isinstance(payload, dict) else False
                            ms_block = _ts_to_ms_like(ts, length or (getattr(chosen_contig, 'length', None) or 0), ploidy, suppress_sfs=suppress_flag)
                            out_parts.append(ms_block)
                except Exception:
                    # Fallback to ms-like conversion if VCF writing fails (or SFS-only fails)
                    try:
                        if not bool(payload.get('only_sfs', False)):
                            suppress_flag = bool(payload.get('suppress_worker_sfs', False)) if isinstance(payload, dict) else False
                            ms_block = _ts_to_ms_like(ts, length or (getattr(chosen_contig, 'length', None) or 0), ploidy, suppress_sfs=suppress_flag)
                            out_parts.append(ms_block)
                        else:
                            # If only_sfs, ensure an SFS marker is emitted
                            try:
                                sfs_af = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
                                sfs_af_line = '#SFS_AF\t' + '\t'.join(str(float(x)) for x in (sfs_af.tolist() if hasattr(sfs_af, 'tolist') else list(sfs_af))) + '\n'
                                err_parts.append(sfs_af_line)
                            except Exception:
                                pass
                    except Exception:
                        err_parts.append(f"# ERROR: failed to produce output for replicate {ridx}\n")

                # Compute and emit allele-frequency-spectrum on stderr as a separate marker
                try:
                    payload['suppress_worker_sfs'] = bool(payload.get('suppress_worker_sfs', False)) if isinstance(payload, dict) else False
                except Exception:
                    payload['suppress_worker_sfs'] = False

                if not payload.get('suppress_worker_sfs', False):
                    # emit a small debug marker for number of sites so the controller can see why SFS may be zero
                    try:
                        nsites_debug = int(getattr(ts, 'num_sites', 0))
                        err_parts.append(f"#DEBUG_NUM_SITES\t{nsites_debug}\n")
                    except Exception:
                        pass
                    try:
                        sfs_af = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
                        sfs_af_line = '#SFS_AF\t' + '\t'.join(str(float(x)) for x in (sfs_af.tolist() if hasattr(sfs_af, 'tolist') else list(sfs_af))) + '\n'
                        err_parts.append(sfs_af_line)
                    except Exception:
                        # non-fatal: don't crash worker if SFS calc fails
                        pass

            except Exception as e:
                import traceback as _tb
                err_parts.append(f"# ERROR: msprime simulate failed for replicate {ridx}: {e}\n{_tb.format_exc()}\n")
                exit_code = 10
                break

    except Exception as e:
        import traceback as _tb
        return ('', f"# ERROR: msprime worker internal failure: {e}\n{_tb.format_exc()}\n", 11)
    finally:
        # restore environment to previous state
        if applied_env:
            try:
                for k, old in old_env_vals.items():
                    if old is None:
                        try:
                            del os.environ[k]
                        except Exception:
                            pass
                    else:
                        os.environ[k] = old
            except Exception:
                pass

    # Stop peak sampler thread if it was started
    try:
        if peak_stop_event is not None:
            try:
                peak_stop_event.set()
            except Exception:
                pass
        if peak_thread is not None:
            try:
                peak_thread.join(timeout=1.0)
            except Exception:
                pass
    except Exception:
        pass

    # Attempt to cleanup any worker-specific tmpdir referenced in env
    try:
        wk = None
        try:
            wk = os.environ.get('SIMULATOR_WORKER_TMP')
        except Exception:
            wk = None
        if wk and os.path.isdir(wk):
            try:
                shutil.rmtree(wk, ignore_errors=True)
                if os.environ.get('SIMULATOR_DEBUG_TMP'):
                    try:
                        sys.stderr.write(f"DEBUG: removed worker tmpdir at process exit: {wk}\n")
                        sys.stderr.flush()
                    except Exception:
                        pass
            except Exception:
                pass
    except Exception:
        pass

    return (''.join(out_parts), ''.join(err_parts), exit_code)

def vcf_to_ms_like(vcf_text, length, ploidy=1):
    """Convert a single-replicate VCF string into a minimal ms-like string.

    Expects a single VCF (one replicate). Positions are normalized to [0,1]
    by dividing POS by length. Supports phased (|) or unphased (/).
    """
    lines = vcf_text.splitlines()
    # find header line with samples
    header_idx = None
    for i,l in enumerate(lines):
        if l.startswith('#CHROM'):
            header_idx = i; break
    if header_idx is None:
        raise SystemExit('ERROR: VCF header (#CHROM) not found for conversion to ms-like.')
    samples = lines[header_idx].strip().split('\t')[9:]
    if not samples:
        # no samples -> empty ms-like
        return 'segsites: 0\n'
    n_ind = len(samples)
    n_hap = n_ind * (ploidy or 1)
    # collect variant lines
    var_lines = [l for l in lines[header_idx+1:] if l and not l.startswith('#')]
    if not var_lines:
        return 'segsites: 0\n'
    pos_list = []
    hap_mat = [[] for _ in range(n_hap)]
    for vl in var_lines:
        parts = vl.split('\t')
        try:
            pos = int(parts[1])
        except Exception:
            raise SystemExit('ERROR: Failed to parse POS in VCF for conversion to ms-like.')
        # normalize to [0,1]
        posf = float(pos) / float(length) if length and length > 0 else 0.0
        pos_list.append(posf)
        fmt = parts[8].split(':') if len(parts) > 8 else ['GT']
        # genotype fields
        gts = [p.split(':')[0] if len(p.split(':'))>0 else '.' for p in parts[9:9+n_ind]]
        # for each sample, decompose into alleles
        hidx = 0
        for s_gt in gts:
            if s_gt == '.' or s_gt == './.' or s_gt == '.|.':
                # unknown -> treat as 0
                alleles = ['0'] * (ploidy or 1)
            else:
                sep = '|' if '|' in s_gt else ('/' if '/' in s_gt else None)
                if sep is None:
                    # single value -> assume haploid
                    alleles = [s_gt]
                else:
                    alleles = s_gt.split(sep)
            # pad/truncate to ploidy
            if len(alleles) < (ploidy or 1):
                alleles = alleles + ['0'] * ((ploidy or 1) - len(alleles))
            for a in alleles[:(ploidy or 1)]:
                # convert to 0/1 char
                ch = '1' if a.strip() == '1' else '0'
                hap_mat[hidx].append(ch)
                hidx += 1
    # build ms-like text
    out_lines = []
    out_lines.append(f"segsites: {len(pos_list)}")
    out_lines.append('positions: ' + ' '.join(str(p) for p in pos_list))
    # each haplotype line
    for h in hap_mat:
        out_lines.append(''.join(h))
    return '\n'.join(out_lines) + '\n'

# ---------- compression helpers ----------

def _write_bgzip_text(path: str, text: str):
    """Write BGZF-compressed text suitable for bcftools/tabix indexing.

    Strategy order:
      1. Use external 'bgzip -c'.
      2. Use pysam.tabix_compress (if pysam available).
      3. Fail with clear error (no silent raw gzip fallback) to avoid unusable files.

    After writing, verify BGZF header (FEXTRA flag + 'BC' subfield) else raise.
    """
    tmp_path = None
    bgzip_exe = shutil.which('bgzip')
    try:
        tmpdir = _ensure_local_tmpdir()
        with tempfile.NamedTemporaryFile('w', delete=False, dir=tmpdir) as tf:
            tf.write(text)
            tmp_path = tf.name
        if bgzip_exe:
            with open(path, 'wb') as outfh:
                proc = subprocess.run([bgzip_exe, '-c', tmp_path], stdout=outfh, stderr=subprocess.PIPE)
            if proc.returncode != 0:
                raise RuntimeError(f"bgzip failed: {proc.stderr.decode().strip()}")
        else:
            try:
                import pysam  # type: ignore
                # pysam.tabix_compress creates bgzip file; force overwrite
                pysam.tabix_compress(tmp_path, path, force=True)
            except Exception as e:
                raise SystemExit('# ERROR: Cannot produce BGZF output (.vcf.gz). Install htslib bgzip or pysam. Details: ' + str(e))
        # Verify BGZF: read first 18 bytes
        with open(path, 'rb') as fh:
            head = fh.read(200)
        if len(head) < 18 or head[0:2] != b'\x1f\x8b' or head[2] != 0x08:
            raise SystemExit('# ERROR: Output not gzip format for BGZF.')
        flg = head[3]
        if not (flg & 0x04):  # FEXTRA bit
            raise SystemExit('# ERROR: gzip file missing FEXTRA; not BGZF (install bgzip).')
        # BGZF extra subfield contains 'BC' marker; search in first 200 bytes
        if b'BC' not in head:
            raise SystemExit('# ERROR: gzip extra field lacks BGZF BC marker; indexing would fail.')
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try: os.unlink(tmp_path)
            except Exception: pass

def _vcf_text_to_bcf(path_vcf: str, path_bcf: str):
    """Convert VCF (plain or gz) text file at path_vcf into BCF at path_bcf using bcftools.

    This function requires bcftools available in PATH. It will overwrite path_bcf if exists.
    """
    bcftools_exe = shutil.which('bcftools') or 'bcftools'
    # Use bcftools view -Ob -o out.bcf in a subprocess
    # Read the input VCF and normalize the #CHROM header line to use tabs
    # (bcftools requires tab-delimited header). If normalization is needed,
    # write to a temporary file and call bcftools on that file.
    tmp_vcf_to_use = path_vcf
    tmp_created = None
    try:
        try:
            with open(path_vcf, 'r', encoding='utf-8', errors='replace') as fh:
                lines = fh.readlines()
        except Exception:
            lines = None
        if lines:
            new_lines = []
            header_changed = False
            for ln in lines:
                # Accept a CHROM header even if indented; bcftools requires the
                # '#CHROM' header to start at column 0 with tab-separated fields.
                if ln.lstrip().startswith('#CHROM'):
                    parts = ln.lstrip().strip().split()
                    if len(parts) >= 9:
                        ln_fixed = '\t'.join(parts) + '\n'
                        if ln_fixed != ln:
                            header_changed = True
                            ln = ln_fixed
                new_lines.append(ln)
            if header_changed:
                tmpdir = _ensure_local_tmpdir()
                fd, tmp_created = tempfile.mkstemp(suffix='.vcf', dir=tmpdir)
                os.close(fd)
                with open(tmp_created, 'w', encoding='utf-8') as ofh:
                    ofh.writelines(new_lines)
                tmp_vcf_to_use = tmp_created
        proc = subprocess.run([bcftools_exe, 'view', '-Ob', '-o', path_bcf, tmp_vcf_to_use], capture_output=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr.decode().strip())
    except Exception as e:
        raise SystemExit(f'# ERROR: failed to convert VCF to BCF via bcftools: {e}')
    finally:
        try:
            if tmp_created and os.path.exists(tmp_created):
                os.unlink(tmp_created)
        except Exception:
            pass

def _ensure_vcf_has_samples(vcf_text: str, meta: dict) -> str:
    """Ensure VCF header '#CHROM' line contains sample columns.

    If header lacks sample columns (only 9 fields), derive sample count from
    metadata (hap_counts or counts_disc and ploidy) and append sample names
    Ind1..IndN to the header line.
    """
    if not vcf_text:
        return vcf_text
    lines = vcf_text.splitlines()
    # find header index
    hidx = None
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith('#CHROM'):
            hidx = i
            break
    if hidx is None:
        return vcf_text
    header_line = lines[hidx].rstrip('\n')
    parts = re.split(r'\t|\s+', header_line.strip())
    if len(parts) >= 10:
        return vcf_text
    # try to derive sample count from metadata
    total_hap = None
    try:
        if meta and isinstance(meta, dict):
            hap_counts = meta.get('hap_counts') or meta.get('counts_disc') or []
            if hap_counts:
                # hap_counts might be list of ints per pop
                total_hap = sum(int(x) for x in hap_counts)
    except Exception:
        total_hap = None
    try:
        ploidy = int(meta.get('ploidy', 1)) if meta and isinstance(meta, dict) else 1
    except Exception:
        ploidy = 1
    n_ind = None
    if total_hap and ploidy and total_hap >= ploidy:
        try:
            n_ind = int(total_hap // ploidy)
        except Exception:
            n_ind = None
    if not n_ind or n_ind <= 0:
        return vcf_text
    sample_names = [f'Ind{i+1}' for i in range(n_ind)]
    new_header = header_line.strip() + '\t' + '\t'.join(sample_names)
    lines[hidx] = new_header
    return '\n'.join(lines) + '\n'

def _adjust_vcf_contig_length(vcf_text: str, chrom_id=None, new_length=None) -> str:
    """Replace contig length for contig with ID==chrom_id in a VCF text block.

    If chrom_id is None, replace the first contig line's length. Returns modified text.
    """
    if not vcf_text:
        return vcf_text
    if new_length is None:
        return vcf_text
    try:
        new_length_int = int(new_length)
    except Exception:
        return vcf_text
    def repl(match):
        id_val = match.group(1)
        if chrom_id is None or str(id_val) == str(chrom_id):
            return f'##contig=<ID={id_val},length={new_length_int}>'
        return match.group(0)
    # Match lines like ##contig=<ID=1,length=248956422>
    out_text, nsub = re.subn(r'##contig=<ID=([^,>]+),length=\d+>', repl, vcf_text)
    if nsub == 0 and chrom_id is None:
        # nothing replaced, try a looser match for any contig and change first occurrence
        out_text = re.sub(r'(##contig=<ID=[^,>]+,length=)\d+', r"\1" + str(new_length_int), vcf_text, count=1)
    return out_text

# ---------- SFS utilities ----------

def _parse_ms_like_replicates(ms_text):
    """Yield haplotype matrix (list of strings) per replicate from ms-like output."""
    lines = ms_text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        ln_stripped = ln.lstrip()
        # accept comment/rep separator even if indented
        if ln_stripped.startswith('//'):
            i += 1
            continue
        if ln_stripped.startswith('segsites:'):
            # parse segsites from stripped content
            try:
                segsites = int(ln_stripped.split()[1])
            except Exception:
                # skip malformed replicate
                i += 1
                continue
            pos_line_idx = None
            # look ahead for positions line (allow leading whitespace)
            j = i + 1
            while j < len(lines) and j <= i + 10:
                if lines[j].lstrip().startswith('positions:'):
                    pos_line_idx = j
                    break
                if lines[j].lstrip().startswith('segsites:'):
                    break
                j += 1
            if pos_line_idx is None:
                i += 1
                continue
            hap_start = pos_line_idx + 1
            hap_lines = []
            k = hap_start
            while k < len(lines):
                hl_raw = lines[k]
                hl = hl_raw.strip()
                if not hl:
                    break
                if hl.lstrip().startswith('//') or hl.lstrip().startswith('segsites:') or hl.lstrip().startswith('positions:'):
                    break
                # accept only contiguous 0/1 strings (strip to remove whitespace)
                if set(hl) <= {'0', '1'}:
                    hap_lines.append(hl)
                k += 1
            if hap_lines:
                yield hap_lines
            i = k
        else:
            i += 1

def compute_sfs(ms_text, normalized=False, mode='mean'):
    """Compute unfolded SFS from ms-like output.

    Returns (header_line, data_lines) where data_lines is list of strings.
    mode: 'mean' averages across replicates; 'per-rep' outputs each replicate separately.
    """
    replicates = list(_parse_ms_like_replicates(ms_text))
    if not replicates:
        return ('#SFS', ['# No replicates parsed'])
    # assume all replicates have same haplotype count
    n_hap = len(replicates[0])
    for haps in replicates:
        if len(haps) != n_hap:
            raise SystemExit('ERROR: Inconsistent haplotype counts across replicates for SFS.')
    bins = n_hap - 1  # k=1..n_hap-1
    if mode == 'per-rep':
        lines = []
        for ridx, haps in enumerate(replicates, start=1):
            if not haps:
                continue
            L = len(haps[0])
            # verify sequences length
            for hp in haps:
                if len(hp) != L:
                    raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
            sfs = [0]*(n_hap+1)
            segsites = L
            for col in range(L):
                k = sum(1 for hp in haps if hp[col] == '1')
                if 0 < k < n_hap:
                    sfs[k] += 1
            if segsites == 0:
                norm = [0]*(n_hap-1)
            else:
                if normalized:
                    norm = [sfs[k]/segsites for k in range(1,n_hap)]
                else:
                    norm = [sfs[k] for k in range(1,n_hap)]
            lines.append('rep'+str(ridx)+'\t'+'\t'.join(str(x) for x in norm))
        header = '#SFS\t' + '\t'.join(str(k) for k in range(1,n_hap))
        return header, lines
    else:  # mean
        accum = [0.0]*(n_hap-1)
        count_used = 0
        for haps in replicates:
            if not haps:
                continue
            L = len(haps[0])
            for hp in haps:
                if len(hp) != L:
                    raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
            sfs = [0]*(n_hap+1)
            segsites = L
            if segsites == 0:
                continue
            for col in range(L):
                k = sum(1 for hp in haps if hp[col] == '1')
                if 0 < k < n_hap:
                    sfs[k] += 1
            if normalized:
                row = [sfs[k]/segsites for k in range(1,n_hap)]
            else:
                row = [sfs[k] for k in range(1,n_hap)]
            accum = [a + r for a,r in zip(accum,row)]
            count_used += 1
        if count_used == 0:
            mean_vals = [0]*(n_hap-1)
        else:
            mean_vals = [a / count_used for a in accum]
        header = '#SFS\t' + '\t'.join(str(k) for k in range(1,n_hap))
        line = 'mean\t' + '\t'.join(str(x) for x in mean_vals)
        return header, [line]

def compute_sfs_fast(ms_text, normalized=False, mode='mean'):
    """Accelerated SFS computation using NumPy where available.

    Falls back to compute_sfs if NumPy is not installed. Interface-compatible.
    """
    try:
        import numpy as _np  # local import to avoid hard dependency
    except Exception:
        return compute_sfs(ms_text, normalized=normalized, mode=mode)

    replicates = list(_parse_ms_like_replicates(ms_text))
    if not replicates:
        return ('#SFS', ['# No replicates parsed'])
    n_hap = len(replicates[0])
    for haps in replicates:
        if len(haps) != n_hap:
            raise SystemExit('ERROR: Inconsistent haplotype counts across replicates for SFS.')
    header = '#SFS\t' + '\t'.join(str(k) for k in range(1, n_hap))
    # Memory-efficient helper that avoids forming an n_hap x L matrix.
    def _sfs_vals_from_haps(haps):
        if not haps:
            return _np.zeros(n_hap-1, dtype=(_np.float32 if normalized else float))
        L = len(haps[0])
        for hp in haps:
            if len(hp) != L:
                raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
        if L == 0:
            return _np.zeros(n_hap-1, dtype=(_np.float32 if normalized else float))
        # Choose the smallest unsigned dtype that can hold counts up to n_hap
        if n_hap <= 255:
            _dtype_counts = _np.uint8
        elif n_hap <= 65535:
            _dtype_counts = _np.uint16
        else:
            _dtype_counts = _np.uint32
        counts_per_site = _np.zeros(L, dtype=_dtype_counts)
        # Accumulate per-site derived counts streaming one haplotype at a time
        for hp in haps:
            arr = _np.frombuffer(hp.encode('ascii'), dtype='S1') == b'1'
            # Add directly; bool will upcast to chosen unsigned dtype
            counts_per_site += arr
        binc = _np.bincount(counts_per_site, minlength=n_hap+1)
        raw = binc[1:n_hap]
        if normalized:
            denom = raw.sum()
            return ((raw / denom).astype(_np.float32) if denom > 0 else _np.zeros_like(raw, dtype=_np.float32))
        # Non-normalized: prefer integer vector to minimize memory; callers may upcast if needed
        return raw.astype(_np.int32, copy=False)

    if mode == 'per-rep':
        out_lines = []
        for ridx, haps in enumerate(replicates, start=1):
            vals = _sfs_vals_from_haps(haps)
            out_lines.append(f"rep{ridx}\t" + '\t'.join(str(v) for v in vals.tolist()))
        return header, out_lines
    else:  # mean mode
        accum = _np.zeros(n_hap-1, dtype=(_np.float32 if normalized else float))
        used = 0
        for haps in replicates:
            vals = _sfs_vals_from_haps(haps)
            if vals.size == 0:
                continue
            accum += vals
            used += 1
        if used == 0:
            mean_vals = [0]*(n_hap-1)
        else:
            # Convert accumulator to float64 for final averaging to reduce rounding error
            try:
                mean_vals = (accum.astype(_np.float64) / float(used)).tolist()
            except Exception:
                mean_vals = (accum / float(used)).tolist()
        # If normalized requested, ensure the mean vector sums to 1 (numerical guard)
        if normalized:
            try:
                s = float(sum(mean_vals))
                if s > 0.0:
                    mean_vals = [v / s for v in mean_vals]
            except Exception:
                pass
        return header, ['mean\t' + '\t'.join(str(v) for v in mean_vals)]

def _stream_sfs_mean_from_file(path, n_hap_expected=None, normalized=False):
    """Stream-parse an ms-like file and compute the mean SFS across replicates.

    - Reads line-by-line to avoid holding entire data in RAM.
    - For each replicate, allocates a per-site counter vector of minimal dtype and
      accumulates haplotype rows one-by-one.
    - Returns (header, ["mean\t..."]).
    """
    import numpy as _np
    used = 0
    accum = None  # float32 when normalized; otherwise int32 for compactness
    header = None

    # Select minimal dtype for per-site counts based on expected n_hap
    def _counts_dtype(nh):
        if nh is None:
            return _np.uint16
        if nh <= 255:
            return _np.uint8
        if nh <= 65535:
            return _np.uint16
        return _np.uint32

    with open(path, 'rt') as fh:
        line = fh.readline()
        while line:
            # find next replicate
            if not line.lstrip().startswith('segsites:'):
                line = fh.readline();
                continue
            # segsites line
            try:
                segsites = int(line.split()[1])
            except Exception:
                segsites = 0
            # find positions: line
            pos_line = None
            for _ in range(10):
                line = fh.readline()
                if not line:
                    break
                s = line.strip()
                if s.startswith('positions:'):
                    pos_line = s
                    break
                if s.startswith('segsites:'):
                    # malformed, step back by one as outer while will see it
                    break
            if pos_line is None:
                # skip malformed replicate
                continue
            # parse number of positions to infer L
            try:
                L = len(pos_line.split(':', 1)[1].strip().split())
            except Exception:
                L = 0
            # if no sites, skip replicate for mean (match previous behavior)
            if L == 0:
                # advance to next replicate boundary
                while True:
                    pos = fh.tell()
                    line = fh.readline()
                    if not line:
                        break
                    st = line.strip()
                    if not st or st.startswith('//') or st.startswith('segsites:'):
                        # rewind one line to let outer loop see 'segsites:'
                        try:
                            fh.seek(pos)
                        except Exception:
                            pass
                        break
                # no contribution
                continue
            # Prepare per-site counts
            counts_dtype = _counts_dtype(n_hap_expected)
            counts = _np.zeros(L, dtype=counts_dtype)
            n_haps_seen = 0
            # Read hap lines
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line:
                    break
                s = line.strip()
                if (not s) or s.startswith('//') or s.startswith('segsites:') or s.startswith('positions:'):
                    # end of replicate block
                    try:
                        fh.seek(pos)
                    except Exception:
                        pass
                    break
                # accept only 0/1 hap strings
                if set(s) <= {'0','1'} and len(s) == L:
                    arr = _np.frombuffer(s.encode('ascii'), dtype='S1') == b'1'
                    counts += arr
                    n_haps_seen += 1
                else:
                    # skip malformed line in this context
                    continue
            # Validate n_hap across replicates if expected provided
            if n_hap_expected is not None and n_haps_seen != int(n_hap_expected):
                raise SystemExit('ERROR: Inconsistent haplotype counts in streaming SFS computation.')
            n_hap = int(n_hap_expected) if n_hap_expected is not None else int(n_haps_seen)
            if header is None:
                header = '#SFS\t' + '\t'.join(str(k) for k in range(1, n_hap))
            binc = _np.bincount(counts, minlength=n_hap+1)
            raw = binc[1:n_hap]
            if normalized:
                denom = raw.sum()
                row = (raw/denom).astype(_np.float32) if denom > 0 else _np.zeros_like(raw, dtype=_np.float32)
            else:
                # keep integers until final mean to minimize memory
                row = raw.astype(_np.int32, copy=False)
            if accum is None:
                # initialize accum matching the selected dtype for row
                accum = row.copy()
            else:
                # ensure shapes match in case of unexpected variability
                if accum.shape != row.shape:
                    raise SystemExit('ERROR: SFS length mismatch across replicates in streaming computation.')
                accum += row
            used += 1
            # continue scanning: read next line for outer loop
            line = fh.readline()

    if used == 0 or accum is None:
        if header is None and n_hap_expected:
            header = '#SFS\t' + '\t'.join(str(k) for k in range(1, int(n_hap_expected)))
        mean_vals = [0]*( (int(n_hap_expected)-1) if n_hap_expected else 0 )
    else:
        # Convert accumulator to float64 for final averaging to reduce rounding error
        try:
            mean_vals = (accum.astype(_np.float64) / float(used)).tolist()
        except Exception:
            mean_vals = (accum / float(used)).tolist()
    # If normalized, renormalize final mean to guard against tiny numerical drift
    if normalized:
        try:
            s = float(sum(mean_vals))
            if s > 0.0:
                mean_vals = [v / s for v in mean_vals]
        except Exception:
            pass
    return header or '#SFS', ['mean\t' + '\t'.join(str(v) for v in mean_vals)]

def compute_sfs_stream_mean(ms_text, n_hap_expected, normalized=False, use_temp_file=True):
    """Compute mean SFS from ms-like text in a memory-efficient, streaming fashion.

    If use_temp_file is True, spill text to a temporary file first so we can drop
    the in-memory copy early (callers can del the string). Otherwise, wrap with
    StringIO; temp file is preferred for large inputs.
    """
    import tempfile, os
    if use_temp_file:
        tmp = None
        tmpdir = _ensure_local_tmpdir()
        try:
            with tempfile.NamedTemporaryFile('w', delete=False, dir=tmpdir) as tf:
                tf.write(ms_text)
                tmp = tf.name
            return _stream_sfs_mean_from_file(tmp, n_hap_expected=n_hap_expected, normalized=normalized)
        finally:
            if tmp and os.path.exists(tmp):
                try: os.unlink(tmp)
                except Exception: pass
    else:
        # Fallback without disk (still stream per line, but holds ms_text in memory)
        import io, numpy as _np
        used = 0
        accum = None
        header = '#SFS\t' + '\t'.join(str(k) for k in range(1, int(n_hap_expected)))
        # Reuse file-based implementation by first writing to disk is better; keep a simple fallback
        buf = io.StringIO(ms_text)
        # Write to a temp file anyway to reuse tested code path
        return compute_sfs_stream_mean(ms_text, n_hap_expected, normalized, use_temp_file=True)

# ---------- discoal builder ----------

def build_discoal_command(*, species, model_id, user_order, individual_counts, discoal_pop0,
                          reps, length, max_fold_per_step, sweep_time, x, s=None, min_snps=None, chr_name=None,
                          disable_en_ladder=False, seed=None, time_units='gens'):
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    # Prefer mutation rate specified on the demographic model if present
    try:
        model_mu = getattr(model, 'mutation_rate', None)
        if model_mu is not None:
            model_mu = float(model_mu)
        else:
            model_mu = None
    except Exception:
        model_mu = None
    # Prefer mutation rate specified on the demographic model if present
    try:
        model_mu = getattr(model, 'mutation_rate', None)
        if model_mu is not None:
            model_mu = float(model_mu)
        else:
            model_mu = None
    except Exception:
        model_mu = None
    raw_demog = model.model
    try:
        graph = raw_demog.to_demes()
    except AttributeError:
        graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')

    graph_names = [d.name for d in graph.demes]
    ploidy = sp.ploidy
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if discoal_pop0 is None:
        discoal_pop0 = user_order[0]
    desired_disc_order = reorder_put_first(user_order, discoal_pop0)
    N0 = deme_present_sizes[discoal_pop0]

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    # use model's mutation rate if available
    if model_mu is not None:
        mu = model_mu
    # use model's mutation rate if available
    if model_mu is not None:
        mu = model_mu
    if chr_name:
        try:
            chosen_contig = sp.get_contig(chr_name)
        except Exception:
            chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig, 'mutation_rate'):
                    mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig, 'recombination_map') and hasattr(chosen_contig.recombination_map, 'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    counts_disc = [name_to_hap[n] for n in desired_disc_order]

    total_hap = sum(counts_disc)
    H = harmonic_number(total_hap - 1)
    if length is None:
        if min_snps is not None:
            # discoal theta later set as theta = 2 * N0 * mu * ploidy * L
            # Expected segregating sites S ≈ theta * H (same harmonic expectation as ms) =>
            # S ≈ (2 * N0 * mu * ploidy * L) * H
            # Solve L = min_snps / (2 * N0 * mu * ploidy * H)
            theta_factor = 2.0 * N0 * mu * ploidy
            if theta_factor * H <= 0:
                raise SystemExit('ERROR: Non-positive theta factor * H, cannot derive length.')
            safety = 1.15  # bias upward so realized S more likely >= target
            length = int(math.ceil(min_snps / (theta_factor * H) * safety))
        else:
            # New fallback: use full chromosome length if chromosome specified
            if chr_name and chosen_contig is not None and hasattr(chosen_contig, 'length'):
                try:
                    length = int(getattr(chosen_contig, 'length'))
                except Exception:
                    raise SystemExit('ERROR: Failed to obtain chromosome length; provide --length or --target-snps.')
            else:
                raise SystemExit('ERROR: Provide --length or --target-snps (or specify chromosome for automatic full length).')

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n, 0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_no_growth = strip_all_growth(ms_no_I)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_disc_idx = {n: i+1 for i,n in enumerate(desired_disc_order)}
    mapping = {name_to_graph_idx[n]: name_to_disc_idx[n] for n in desired_disc_order}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_remap = map_ej_to_ed(demog_remap)
    demog_0b = shift_to_0_based(demog_remap)

    if chr_name:
        theta = 2.0 * N0 * mu * ploidy * length
        rho   = 2.0 * N0 * rrate * ploidy * length
    else:
        theta = 2.0 * ploidy * N0 * mu * length
        rho   = 2.0 * ploidy * N0 * rrate * length

    sel_args = []
    # discoal -ws t x a : t interpreted here as ORIGIN time (when allele became beneficial) per user spec.
    # If user supplied per-generation s, compute 2Ns using N0
    if s is not None:
        try:
            a = 2.0 * N0 * float(s)
        except Exception:
            raise SystemExit('ERROR: Failed to convert --sel-s to 2Ns using N0.')
    else:
        a = None
    # Interpret sweep_time according to time_units
    sweep_time_user = sweep_time
    sweep_time_4N = None
    if sweep_time_user is not None:
        if time_units == 'gens':
            try:
                sweep_time_4N = float(sweep_time_user) / (4.0 * float(N0))
            except Exception:
                raise SystemExit('ERROR: Failed to convert --sweep-time (generations) to discoal 4N units using N0.')
        else:  # time_units == '4N'
            sweep_time_4N = float(sweep_time_user)

    if (a is not None) and (sweep_time_4N is not None):
        if x is None:
            x = 0.5
        sel_args = ['-ws', str(sweep_time_4N), '-x', str(x), '-a', str(a)]

    npop_disc = len(desired_disc_order)
    total_hap_disc = sum(counts_disc)
    seed_part = f" -seed {seed}" if seed is not None else ''
    discoal_cmd = (
        f"# pop0 (sweep-capable) deme: {discoal_pop0}; indices 0-based; order: {desired_disc_order}\n"
        f"discoal {total_hap_disc} {reps} {length} -t {theta} -r {rho} -p {npop_disc} {' '.join(map(str,counts_disc))}{seed_part} "
        f"{' '.join(sel_args)} {demog_0b}".strip()
    )
    meta = {
        'N0': N0, 'mu': mu, 'rrate': rrate, 'counts_disc': counts_disc, 'npop': npop_disc,
        'demog': demog_0b, 'pop_order': desired_disc_order, 'pop0': discoal_pop0,
        'ploidy': ploidy, 'sel_args': sel_args, 'length': length,
        'chromosome': chr_name if chr_name else None, 'species_id': species,
        'sweep_pos': x if sel_args else None,
        'sel_2Ns': a if sel_args else None,
        'sweep_time': sweep_time if sel_args else None,
        'fixation_time': None,
        'growth_discretization': {
            'requested_disable': bool(disable_en_ladder),
            'max_fold_per_step': float(max_fold_per_step),
            'applied': (not disable_en_ladder and len(en_ladder) > 0),
            'steps': len(en_ladder),
            'reason': ('disabled_by_flag' if disable_en_ladder else ('stepwise_en_ladder' if len(en_ladder) > 0 else 'no_exponential_epochs')),
        },
        'base_comments': [f"# pop0 (sweep-capable) deme: {discoal_pop0}; indices 0-based; order: {desired_disc_order}",
                          f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    if sel_args and x is not None:
        try:
            sweep_bp = int(round(float(x) * length)) if length else None
        except Exception:
            sweep_bp = None
        loc_line = f"# sweep-pos={x}"
        if sweep_bp is not None:
            loc_line += f" sweep-bp={sweep_bp}"
        if sweep_time is not None:
            loc_line += f" sweep-time={sweep_time}"
        meta['base_comments'].append(loc_line)
    # no debug printing
    if seed is not None:
        meta['seed'] = seed
    return discoal_cmd, meta

# ---------- ms (neutral) builder ----------

def build_ms_command(*, species, model_id, user_order, individual_counts, pop0,
                     reps, length, max_fold_per_step, chr_name=None,
                     min_snps=None, disable_en_ladder=False, seed=None):
    engine = 'ms'
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try: graph = raw_demog.to_demes()
    except AttributeError: graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')
    ploidy = sp.ploidy
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if pop0 is None: pop0 = user_order[0]
    if pop0 not in deme_present_sizes:
        raise SystemExit(f"ERROR: pop0 '{pop0}' not in model demes {list(deme_present_sizes)}")
    N0 = deme_present_sizes[pop0]

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try: chosen_contig = sp.get_contig(chr_name)
        except Exception: chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig,'mutation_rate'): mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig,'recombination_map') and hasattr(chosen_contig.recombination_map,'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    total_hap_placeholder = sum(ploidy * c for c in individual_counts)
    if length is None:
        if min_snps is not None:
            H = harmonic_number(total_hap_placeholder - 1)
            theta_per_bp = 4.0 * N0 * mu  # theta = 4 N0 mu L
            if theta_per_bp * H <= 0:
                raise SystemExit('ERROR: Cannot derive length for ms (theta_per_bp*H<=0).')
            safety = 1.10
            length = int(math.ceil(min_snps / (theta_per_bp * H) * safety))
        elif chr_name and chosen_contig is not None:
            length = int(getattr(chosen_contig,'length'))
        else:
            # New fallback: if chromosome specified but length not retrievable -> error message adjusted
            raise SystemExit('ERROR: Provide --length or --target-snps (or specify valid chromosome) for engine=ms/msms.')

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    # Ensure consistency with discoal: put the sweep-capable population first
    # so msms -I ordering and -SI initial-frequency vectors align when
    # --sweep-pop is provided. If pop0 is None, default to the first user_order.
    desired_order = reorder_put_first(user_order, pop0)
    counts_order = [name_to_hap[n] for n in desired_order]
    total_hap = sum(counts_order)

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n,0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    ms_no_growth = strip_all_growth(ms_no_I)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_user_idx = {n: i+1 for i,n in enumerate(desired_order)}
    mapping = {name_to_graph_idx[n]: name_to_user_idx[n] for n in desired_order if n in name_to_graph_idx}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_final = sanitize_ms_demography(demog_remap)

    theta = 4.0 * N0 * mu * length
    rho = 4.0 * N0 * rrate * length
    npop = len(desired_order)

    base_parts = [engine, str(total_hap), str(reps),
                  '-t', str(theta), '-r', str(rho), str(length),
                  '-I', str(npop)] + list(map(str, counts_order))
    # Intentionally NOT adding -seeds to avoid external seed file creation per user request.
    base_cmd_line = ' '.join(base_parts) + ' ' + demog_final
    full_cmd = (f"# engine={engine}\n# pop0={pop0}\n# order={desired_order}\n# N0={N0} mu={mu} r={rrate} length={length}\n{base_cmd_line}")
    meta = {
        'engine': engine, 'N0': N0, 'mu': mu, 'rrate': rrate,
        'counts_disc': counts_order, 'npop': npop, 'demog': demog_final,
        'pop_order': desired_order, 'pop0': pop0, 'ploidy': ploidy,
        'sel_args': [], 'length': length,
        'chromosome': chr_name if chr_name else None,
        'species_id': species,
        'sweep_pos': None,
        'sel_2Ns': None,
        'sweep_time': None,
        'fixation_time': None,
        'growth_discretization': {
            'requested_disable': bool(disable_en_ladder),
            'max_fold_per_step': float(max_fold_per_step),
            'applied': (not disable_en_ladder and bool(en_ladder)),
            'steps': len(en_ladder),
            'reason': ('disabled_by_flag' if disable_en_ladder else ('stepwise_en_ladder' if en_ladder else 'no_exponential_epochs')),
        },
    'base_comments': ['# engine='+engine, f"# pop0={pop0}", f"# order={desired_order}",
              f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    # no debug printing
    if seed is not None:
        meta['seed'] = seed
    return full_cmd, meta

# ---------- scrm (SMC approximation) builder ----------
def build_scrm_command(*, species, model_id, user_order, individual_counts, pop0,
                       reps, length, max_fold_per_step, chr_name=None,
                       min_snps=None, disable_en_ladder=False, seed=None):
    """Build a scrm command string mirroring the ms builder but emitting 'scrm' as the executable.

    scrm is a commonly-used SMC approximation simulator with ms-like CLI. This builder
    mirrors the ms command construction so downstream handling (parsing -t/-r/-I etc.) works
    unchanged while allowing users to request the approximation engine via --engine scrm.
    """
    # reuse most logic from build_ms_command but swap engine name
    engine = 'scrm'
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try: graph = raw_demog.to_demes()
    except AttributeError: graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')
    ploidy = sp.ploidy
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if pop0 is None: pop0 = user_order[0]
    if pop0 not in deme_present_sizes:
        raise SystemExit(f"ERROR: pop0 '{pop0}' not in model demes {list(deme_present_sizes)}")
    N0 = deme_present_sizes[pop0]

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try: chosen_contig = sp.get_contig(chr_name)
        except Exception: chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig,'mutation_rate'): mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig,'recombination_map') and hasattr(chosen_contig.recombination_map,'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    total_hap_placeholder = sum(ploidy * c for c in individual_counts)
    if length is None:
        if min_snps is not None:
            H = harmonic_number(total_hap_placeholder - 1)
            theta_per_bp = 4.0 * N0 * mu  # theta = 4 N0 mu L
            if theta_per_bp * H <= 0:
                raise SystemExit('ERROR: Cannot derive length for scrm (theta_per_bp*H<=0).')
            safety = 1.10
            length = int(math.ceil(min_snps / (theta_per_bp * H) * safety))
        elif chr_name and chosen_contig is not None:
            length = int(getattr(chosen_contig,'length'))
        else:
            raise SystemExit('ERROR: Provide --length or --target-snps (or specify valid chromosome) for engine=scrm.')

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    # Put the sweep population (pop0) first in the sampling order to match
    # discoal behaviour and ensure msms -SI init-frequency vectors align.
    desired_order = reorder_put_first(user_order, pop0)
    counts_order = [name_to_hap[n] for n in desired_order]
    total_hap = sum(counts_order)

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n,0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    ms_no_growth = strip_all_growth(ms_no_I)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_user_idx = {n: i+1 for i,n in enumerate(desired_order)}
    mapping = {name_to_graph_idx[n]: name_to_user_idx[n] for n in desired_order if n in name_to_graph_idx}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_final = sanitize_ms_demography(demog_remap)

    theta = 4.0 * N0 * mu * length
    rho = 4.0 * N0 * rrate * length
    npop = len(desired_order)

    # scrm uses an ms-like CLI; emit 'scrm' as executable and keep ms-like flags for downstream parsing
    base_parts = [engine, str(total_hap), str(reps),
                  '-t', str(theta), '-r', str(rho), str(length),
                  '-I', str(npop)] + list(map(str, counts_order))
    base_cmd_line = ' '.join(base_parts) + ' ' + demog_final
    full_cmd = (f"# engine={engine}\n# pop0={pop0}\n# order={desired_order}\n# N0={N0} mu={mu} r={rrate} length={length}\n{base_cmd_line}")
    meta = {
        'engine': engine, 'N0': N0, 'mu': mu, 'rrate': rrate,
        'counts_disc': counts_order, 'npop': npop, 'demog': demog_final,
        'pop_order': desired_order, 'pop0': pop0, 'ploidy': ploidy,
        'sel_args': [], 'length': length,
        'chromosome': chr_name if chr_name else None,
        'species_id': species,
        'sweep_pos': None,
        'sel_2Ns': None,
        'sweep_time': None,
        'fixation_time': None,
        'growth_discretization': {
            'requested_disable': bool(disable_en_ladder),
            'max_fold_per_step': float(max_fold_per_step),
            'applied': (not disable_en_ladder and bool(en_ladder)),
            'steps': len(en_ladder),
            'reason': ('disabled_by_flag' if disable_en_ladder else ('stepwise_en_ladder' if en_ladder else 'no_exponential_epochs')),
        },
    'base_comments': ['# engine='+engine, f"# pop0={pop0}", f"# order={desired_order}",
              f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    # no debug printing
    if seed is not None:
        meta['seed'] = seed
    return full_cmd, meta

# ---------- msprime (in-process stdpopsim engine) builder ----------
def build_msprime_command(*, species, model_id, user_order, individual_counts, pop0,
                          reps, length, max_fold_per_step, chr_name=None,
                          min_snps=None, disable_en_ladder=False, seed=None):
    """Prepare metadata for running stdpopsim.engines._MsprimeEngine in-process.

    This returns a small command-like string (for debug) and a meta dict compatible with other builders.
    """
    engine = 'msprime'
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try: graph = raw_demog.to_demes()
    except AttributeError: graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')
    ploidy = sp.ploidy
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if pop0 is None: pop0 = user_order[0]
    if pop0 not in deme_present_sizes:
        raise SystemExit(f"ERROR: pop0 '{pop0}' not in model demes {list(deme_present_sizes)}")
    N0 = deme_present_sizes[pop0]

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try: chosen_contig = sp.get_contig(chr_name)
        except Exception: chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig,'mutation_rate'): mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig,'recombination_map') and hasattr(chosen_contig.recombination_map,'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    # Put sweep population first so msms -I/-SI vectors align with discoal.
    desired_order = reorder_put_first(user_order, pop0)
    counts_order = [name_to_hap[n] for n in desired_order]

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n,0) for n in graph_order]

    # For debug and compatibility, produce a short pseudo-command line
    total_hap = sum(counts_order)
    # Provide an ms-like engine line so the runner can find token[2]=reps
    engine_line = f"{engine} {total_hap} {reps} -L {length}"
    pseudo_cmd = (f"# engine={engine}\n# pop0={pop0}\n# order={desired_order}\n# N0={N0} mu={mu} r={rrate} length={length}\n"
                  f"{engine_line}")
    meta = {
        'engine': engine, 'N0': N0, 'mu': mu, 'rrate': rrate,
        'counts_disc': counts_order, 'npop': len(desired_order), 'demog': None,
        'pop_order': desired_order, 'pop0': pop0, 'ploidy': ploidy,
        'length': length, 'reps': reps, 'chromosome': chr_name if chr_name else None,
        'species_id': species, 'graph': graph, 'model_obj': model,
        'samples_graph_order': samples_graph_order,
        'base_comments': ['# engine=msprime', f"# pop0={pop0}", f"# order={desired_order}"],
        # Growth discretization is a no-op for msprime (continuous-size support)
        'growth_discretization': {
            'requested_disable': bool(disable_en_ladder),
            'max_fold_per_step': float(max_fold_per_step),
            'applied': False,
            'steps': 0,
            'reason': 'msprime_supports_continuous_growth',
        }
    }
    # no debug printing
    if seed is not None:
        meta['seed'] = seed
    return pseudo_cmd, meta

# ---------- msms (selection) builder with SFC + Smark ----------

def build_msms_command(*, species, model_id, user_order, individual_counts, pop0,
                       reps, length, max_fold_per_step, chr_name=None,
                       min_snps=None, disable_en_ladder=False, a=None, s=None, x=None,
                       fixation_time=None, sweep_time=None, seed=None, time_units='gens'):
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try: graph = raw_demog.to_demes()
    except AttributeError: graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')

    ploidy = sp.ploidy
    if pop0 is None:
        pop0 = user_order[0]
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if pop0 not in deme_present_sizes:
        raise SystemExit(f"ERROR: pop0 '{pop0}' not in model demes {list(deme_present_sizes)}")
    N0 = float(deme_present_sizes[pop0])

    # If user provided per-generation s, convert to 2Ns and prefer it over a
    if s is not None:
        try:
            a = 2.0 * N0 * float(s)
        except Exception:
            raise SystemExit('ERROR: Failed to convert --sel-s to 2Ns using N0.')

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try: chosen_contig = sp.get_contig(chr_name)
        except Exception: chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig,'mutation_rate'): mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig,'recombination_map') and hasattr(chosen_contig.recombination_map,'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    total_hap_placeholder = sum(ploidy * c for c in individual_counts)
    if length is None:
        if min_snps is not None:
            H = harmonic_number(total_hap_placeholder - 1)
            theta_per_bp = 4.0 * N0 * mu  # neutral expectation
            if theta_per_bp * H <= 0:
                raise SystemExit('ERROR: Cannot derive length for msms (theta_per_bp*H<=0).')
            # Under selective sweeps segregating sites are reduced; inflate length more aggressively
            base_len = min_snps / (theta_per_bp * H)
            safety = 1.25 if (a is not None and a > 0) else 1.10
            length = int(math.ceil(base_len * safety))
        elif chr_name and chosen_contig is not None:
            length = int(getattr(chosen_contig,'length'))
        else:
            raise SystemExit('ERROR: Provide --length or --target-snps (or specify chromosome) for engine=msms.')

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    # Put the sweep population first in msms sampling order to align -I and -SI vectors.
    #
    # Rationale: msms expects positional vectors (-I counts, -SI init freqs) where
    # each position corresponds to a population index. discoal already reorders
    # the sample list so the sweep-capable population comes first; to keep
    # semantics consistent across builders we do the same for msms. This ensures
    # that specifying --sweep-pop=<name> applies selection to the named deme
    # regardless of the order the user passed to --pop-order.
    # Note: neutral engines (ms, scrm, msprime) ignore selection flags and do
    # not require a sweep-pop argument.
    desired_order = reorder_put_first(user_order, pop0)
    counts_order = [name_to_hap[n] for n in desired_order]
    total_hap = sum(counts_order)

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n,0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    ms_no_growth = strip_all_growth(ms_no_I)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_user_idx = {n: i+1 for i,n in enumerate(desired_order)}
    mapping = {name_to_graph_idx[n]: name_to_user_idx[n] for n in desired_order if n in name_to_graph_idx}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_final = sanitize_ms_demography(demog_remap)

    # time-invariant detection (for completeness; we will prefer SFC anyway)
    time_invariant = not any(flag in (' ' + demog_final + ' ')
                             for flag in (' -en ', ' -em ', ' -ej ', ' -eg ', ' -eG ', ' -es ', ' -ema '))

    theta = 4.0 * N0 * mu * length
    rho   = 4.0 * N0 * rrate * length
    npop  = len(desired_order)

    parts = ['msms', str(total_hap), str(reps),
             '-t', str(theta), '-r', str(rho), str(length),
             '-N', str(int(round(N0))),
             '-I', str(npop), *map(str, counts_order)]

    sel_args = []
    if (a is not None) or (x is not None):
        if a is None:
            # should not happen because we convert s->a earlier, but keep defensive check
            raise SystemExit('ERROR: engine=msms sweep requested but --sel-2Ns or --sel-s missing.')
        # ensure numeric
        try:
            a = float(a)
        except Exception:
            raise SystemExit('ERROR: --sel-2Ns/--sel-s could not be interpreted as a number.')
        # genic (h≈0.5)
        sel_args += ['-SaA', str(a), '-SAA', str(2.0 * a)]
        # selected site position + mark it in output
        if x is None:
            x = 0.5
        sel_args += ['-Sp', str(x), '-Smark']
    # Derive origin time (backwards) and handle fixation-time / sweep_time semantics.
    # Approx fixation duration in coalescent units (duration from origin to fixation):
    # T_duration ≈ ln(a)/a for a>1 else ln(1+a)/a. (a = 2Ns)
    # If fixation_time > 0 (meaning the allele was fixed in the past at time "fixation_time"), we will
    # instruct msms to set the allele frequency to ~1.0 at that time via `-SI <t> ...` and NOT use `-SFC`.
    # If fixation_time == 0 (fixation at present), we keep the existing behavior: initialize at a small
    # frequency at tstart and use `-SFC` to condition on fixation at present.
        if a <= 0:
            raise SystemExit('ERROR: --sel-2Ns/--sel-s must be > 0 for msms selection.')
        if a > 1.0:
            t_duration = math.log(a) / a
        else:
            t_duration = math.log(1.0 + a) / a
    # interpret fixation_time or sweep_time according to time_units, convert to 4N units when needed
        fixation_time_gen = None
        fixation_time_4N = 0.0
        if fixation_time is not None:
            if time_units == 'gens':
                fixation_time_gen = fixation_time
                try:
                    if not (N0 and float(N0) > 0):
                        raise Exception('Invalid N0 for conversion')
                    fixation_time_4N = float(fixation_time_gen) / (4.0 * float(N0))
                except Exception:
                    raise SystemExit('ERROR: Failed to convert --fixation-time (generations) to 4N units using N0.')
            else:
                fixation_time_4N = float(fixation_time)
        # sweep_time represents the origin time when the allele became beneficial (same semantics as discoal)
        sweep_time_gen = None
        sweep_time_4N = None
        if fixation_time is not None:
            if time_units == 'gens':
                fixation_time_gen = fixation_time
                try:
                    if not (N0 and float(N0) > 0):
                        raise Exception('Invalid N0 for conversion')
                    fixation_time_4N = float(fixation_time_gen) / (4.0 * float(N0))
                except Exception:
                    raise SystemExit('ERROR: Failed to convert --fixation-time (generations) to 4N units using N0.')
            else:
                fixation_time_4N = float(fixation_time)
        if sweep_time is not None:
            if time_units == 'gens':
                sweep_time_gen = sweep_time
                try:
                    if not (N0 and float(N0) > 0):
                        raise Exception('Invalid N0 for conversion')
                    sweep_time_4N = float(sweep_time_gen) / (4.0 * float(N0))
                except Exception:
                    raise SystemExit('ERROR: Failed to convert --sweep-time (generations) to msms 4N units using N0.')
            else:
                sweep_time_4N = float(sweep_time)

        # t_duration is the approximate time (in 4N units) from origin to fixation
        tstart_duration = t_duration

    # Priority: if sweep_time provided, treat it as origin time where selected allele was at low freq
    # and will sweep to fixation at present (we initialize at a small freq at sweep_time and use -SFC to
    # condition on fixation at present). If sweep_time is not provided but fixation_time is provided
    # and >0, that means the allele was already fixed at that past time; we model this by setting the
    # allele frequency ~1 at fixation_time via -SI and NOT using -SFC.
        # initialize at a small frequency at origin time (t_duration + fixation_time) and use -SFC to
        # condition on fixation at present.
        if fixation_time is None or float(fixation_time) == 0.0:
            tstart = tstart_duration + (fixation_time_4N or 0.0)
            f0 = max(1.0 / (2.0 * N0), 1e-6)
            init_freqs = ['0'] * npop
            sel_deme_index = desired_order.index(pop0)  # 0-based index into our order list
            init_freqs[sel_deme_index] = str(f0)
            sel_args += ['-SI', str(tstart), str(npop), *init_freqs, '-SFC']
        else:
            # fixation_time > 0: user wants the allele to be fixed at time fixation_time in the past.
            # Emit -SI <fixation_time_4N> ... with frequency ~1.0 for the sweep population and DO NOT add -SFC.
            # Use a numerically safe frequency slightly below 1.0: 1 - 1/(2N0) (haploid sample scaling)
            try:
                safe_one = 1.0 - (1.0 / (2.0 * float(N0)))
                if safe_one <= 0.0:
                    safe_one = 0.999999
            except Exception:
                safe_one = 0.999999
            init_freqs = ['0'] * npop
            sel_deme_index = desired_order.index(pop0)
            init_freqs[sel_deme_index] = str(safe_one)
            sel_args += ['-SI', str(fixation_time_4N), str(npop), *init_freqs]
        if sweep_time is not None:
            # initialize at small frequency at sweep origin (sweep_time_4N) and use -SFC to condition on fixation
            tstart = tstart_duration + (sweep_time_4N or 0.0)
            f0 = max(1.0 / (2.0 * N0), 1e-6)
            init_freqs = ['0'] * npop
            sel_deme_index = desired_order.index(pop0)
            init_freqs[sel_deme_index] = str(f0)
            sel_args += ['-SI', str(tstart), str(npop), *init_freqs, '-SFC']
        elif fixation_time is None or float(fixation_time) == 0.0:
            # default / present-fixation behavior: initialize a small freq at tstart and use -SFC
            tstart = tstart_duration + (fixation_time_4N or 0.0)
            f0 = max(1.0 / (2.0 * N0), 1e-6)
            init_freqs = ['0'] * npop
            sel_deme_index = desired_order.index(pop0)  # 0-based index into our order list
            init_freqs[sel_deme_index] = str(f0)
            sel_args += ['-SI', str(tstart), str(npop), *init_freqs, '-SFC']
        else:
            # fixation_time > 0: user wants the allele to be fixed at time fixation_time in the past.
            # Emit -SI <fixation_time_4N> ... with frequency ~1.0 for the sweep population and DO NOT add -SFC.
            try:
                safe_one = 1.0 - (1.0 / (2.0 * float(N0)))
                if safe_one <= 0.0:
                    safe_one = 0.999999
            except Exception:
                safe_one = 0.999999
            init_freqs = ['0'] * npop
            sel_deme_index = desired_order.index(pop0)
            init_freqs[sel_deme_index] = str(safe_one)
            sel_args += ['-SI', str(fixation_time_4N), str(npop), *init_freqs]

    if sel_args:
        parts.extend(sel_args)
    parts.extend(demog_final.split())

    if seed is not None:
        parts.extend(['-seed', str(seed)])
    full_cmd = (f"# engine=msms\n# pop0={pop0}\n# order={desired_order}\n# N0={N0} mu={mu} r={rrate} length={length}\n{' '.join(parts)}")
    meta = {
        'engine': 'msms', 'N0': N0, 'mu': mu, 'rrate': rrate,
        'counts_disc': counts_order, 'npop': npop, 'demog': demog_final,
        'pop_order': desired_order, 'pop0': pop0, 'ploidy': ploidy,
        'sel_args': sel_args, 'length': length,
        'chromosome': chr_name if chr_name else None,
        'species_id': species,
        'sweep_pos': x if sel_args else None,
        'sel_2Ns': a if sel_args else None,
        'sel_s': (s if sel_args else None),
        'sweep_time': None,
        'fixation_time': fixation_time if sel_args else None,
        'growth_discretization': {
            'requested_disable': bool(disable_en_ladder),
            'max_fold_per_step': float(max_fold_per_step),
            'applied': (not disable_en_ladder and bool(en_ladder)),
            'steps': len(en_ladder),
            'reason': ('disabled_by_flag' if disable_en_ladder else ('stepwise_en_ladder' if en_ladder else 'no_exponential_epochs')),
        },
        'base_comments': ['# engine=msms', f"# pop0={pop0}", f"# order={desired_order}",
              f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    if sel_args and x is not None:
        try:
            sweep_bp = int(round(float(x) * length)) if length else None
        except Exception:
            sweep_bp = None
        loc_line = f"# sweep-pos={x}"
        if sweep_bp is not None:
            loc_line += f" sweep-bp={sweep_bp}"
        if fixation_time:
            loc_line += f" fixation-time={fixation_time}"
        meta['base_comments'].append(loc_line)
    # no debug printing
    if seed is not None:
        meta['seed'] = seed
    return full_cmd, meta

# ---------- CLI ----------

class _RequiredAnnotatingFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """HelpFormatter that appends '(required)' to help text for required options."""
    def add_argument(self, action):
        if action.required:
            if action.help:
                if '(required)' not in action.help:
                    action.help += ' (required)'
            else:
                action.help = '(required)'
        super().add_argument(action)

def parse_args():
    """Unified parser: always show selection arguments; validate post parse based on --engine."""
    ap = argparse.ArgumentParser(
        description=(
            "Simulate genomic data using stdpopsim demographies and write ms-like or VCF/BCF outputs.\n"
            "Engines: discoal (sweep origin), ms (neutral), msms (sweep fixation), scrm (SMC), msprime (in-process).\n"
            "Outputs: ms/ms.gz (raw ms-like) or vcf/vcf.gz/bcf (converted). For VCF/BCF with multiple replicates, one file per replicate is written inside a folder named after the output stem.\n"
            "Extras: paired-neutral baseline (engine defaults to --engine), unfolded SFS computation (VCF auto-converted internally), target-snps length derivation and per-replicate enforcement with optional overshoot shrink (see --target-snps-tol), show-command preview for external engines.\n"
            "Species may be given by identifier or full common/scientific name (e.g. Homo sapiens)."
        ),
        formatter_class=_RequiredAnnotatingFormatter,
    )
    g = ap.add_argument_group('Core')
    g.add_argument('--engine', required=True, choices=['discoal','ms','msms','scrm','msprime'],
                   help="Simulator engine: discoal (selection origin), ms (neutral), msms (selection fixation), scrm (SMC approximation), msprime (stdpopsim/msprime in-process).")
    g.add_argument('--species-id', dest='species', required=True,
                   help="Species identifier or name as accepted by stdpopsim (e.g. 'Homo sapiens').")
    g.add_argument('--model-id', dest='model', required=True,
                   help="Demographic model id for the chosen species (see stdpopsim models).")
    g.add_argument('--pop-order', dest='demes', required=True,
                   help="Comma list of deme/population names specifying the sample order (e.g. 'popA,popB').")
    g.add_argument('--sample-individuals', dest='samples_config', required=True,
                   help="Comma integers: individuals sampled per deme in same order as --pop-order.")
    g.add_argument('--length', type=int, default=None,
                   help="Sequence length (bp). If omitted, may be derived from --target-snps or contig info.")
    g.add_argument('--target-snps', dest='min_snps', type=int, default=None,
                   help="Target/minimum segregating sites; derives --length when not provided.")
    g.add_argument('--target-snps-tol', dest='min_snps_tol', type=float, default=0.0,
             help=("Tolerance for --target-snps overshoot, specified as a percent in [0,100]. "
                 "Internally the value is converted to a fractional tolerance (percent/100). Set 0 to disable shrink (default 0). "
                 "Example: 2 means ~2% overshoot tolerated before shrink attempts."))
    g.add_argument('--chromosome', dest='chromosome', default=None,
                   help="Chromosome/contig name for mutation/recombination rates/length (optional).")
    g.add_argument('--replicates', dest='reps', type=int, default=1,
                   help="Number of simulation replicates to run.")
    g.add_argument('--output', dest='out', default=None,
                   help=("Output file path. If omitted, modelid_pop0[_chrom][_sweeppop].EXT is written (EXT per --output-format). "
                         "Use '-' for stdout."))
    g.add_argument('--temp', dest='temp_loc', choices=['local','system'], default='local',
                   help=("Where to place intermediate temporary files: 'local' (repo-local hidden folder)"
                         " or 'system' (use Python's tempfile.gettempdir()). Default: local."))
    g.add_argument('--output-format', dest='format', default='ms', choices=['ms','ms.gz','vcf','vcf.gz','bcf'],
                   help="Output format. 'ms'/'ms.gz' write raw ms-like (optionally gzipped); 'vcf'/'vcf.gz' convert single replicate to VCF.")
    g.add_argument('--parallel', dest='parallel', type=int, default=1,
                   help="Parallel workers for replicate chunks.")
    g.add_argument('--sims-per-work', dest='sims_per_work', type=int, default=None,
             help=("Number of simulation replicates to group per worker invocation (per core/thread) for external engines. "
                 "If omitted the runner will split total replicates into exactly --parallel parts (parts need not be equal); "
                 "the base size is floor(replicates/parallel) and any remainder is added to the final part. "
                 "This flag is ignored for engine=msprime (internal in-process batching). "
                 "Minimum 1; values greater than total replicates are capped to total replicates."))
    g.add_argument('--max-ram-percent', dest='max_ram_percent', type=float, default=80.0,
             help=("Numeric threshold (percent 0..100) used by the RAM watchdog. Interpretation depends on"
                 " --max-ram-cap (system|run). Default 80.0."))
    g.add_argument('--max-ram-cap', dest='max_ram_cap', choices=['system', 'run'], default='system',
             help=("How to interpret --max-ram-percent: 'system' (default) compares to total system RAM usage percent;"
                 " 'run' compares the current process RSS as a percentage of total system RAM (i.e. the run's share)."
                 " Note: if system RAM usage reaches 99%% the run will be terminated regardless of this setting."))
    # Show-only command flag
    g.add_argument('--show-command', action='store_true', help="Print the external simulator command that would run and exit (no execution).")
    g.add_argument('--progress', action='store_true', help="Show progress bar (requires tqdm).")
    g.add_argument('--sfs', nargs='?', const=True, default=False,
                   help="Compute unfolded SFS. Optional value sets output basename; default name modelid_pop0[_chrom][_sweeppop].sfs.")
    # per-chunk output deletion is enabled by default when --sfs is requested
    g.add_argument('--sfs-normalized', action='store_true', help='Normalize SFS counts to frequencies.')
    # Deprecated flags --force-resim / --no-force-resim removed: worker AF markers are authoritative.
    g.add_argument('--sfs-mode', choices=['mean','per-rep'], default='mean', help="SFS aggregation mode.")

    pn = ap.add_argument_group('Paired neutral baseline')
    pn.add_argument('--paired-neutral', nargs='?', const=True, default=False,
                    help=(
                        'Also run a neutral simulation with identical parameters (selection removed). '
                        'Optionally provide a name/path for the neutral artifact, e.g., --paired-neutral neutral.ms. '
                        'If omitted, the neutral filename is the primary name with _neutral before the extension. '
                        'In SFS-only mode (no --output), the provided value serves as the neutral SFS basename.'
                    ))
    pn.add_argument('--neutral-engine', choices=['discoal','ms','msms','scrm','msprime'], help='Engine to use for neutral run (default: same as --engine).')

    gd = ap.add_argument_group('Demography')
    gd.add_argument('--disable-growth-discretization', dest='no_en_ladder', action='store_true', help='Disable exponential growth discretization.')
    gd.add_argument('--growth-max-fold', dest='max_fold_per_step', type=float, default=1.05, help='Max fold change per discretized growth step.')

    sg = ap.add_argument_group('Selection (engine-specific), required for discoal/msms')
    sg.add_argument('--sweep-pop', dest='discoal_pop0', help='Sweep population (discoal/msms). Ignored for ms.')
    sg.add_argument('--sweep-pos', dest='x', type=float, help='Selected site position as percent in [0,100] (required for discoal/msms). Ignored for ms,scrm and msprime.')
    sg.add_argument('--sel-s', dest='s', type=float, help='Selection coefficient s per generation (required for discoal/msms). Converted to 2Ns using present-size N for the sweep population. Ignored for ms,scrm and msprime.')
    sg.add_argument('--sweep-time', dest='sweep_time', type=float, help='Origin time for the allele (units set by --time-units; default gens).')
    sg.add_argument('--fixation-time', dest='fix_time', type=float, default=None, help='Fixation time back (units set by --time-units; gens=generations ago). If omitted, no explicit past-fixation time is assumed.')
    sg.add_argument('--time-units', dest='time_units', choices=['gens','4N'], default='gens',
                   help="Units for --sweep-time and --fixation-time: 'gens' = generations ago, '4N' = coalescent units (4N generations). Default: gens.")

    args = ap.parse_args()
    # No legacy debug flag; use --show-command for previews.
    engine = args.engine

    # Start optional RAM monitor thread if requested (default 80%). It checks
    # total system RAM usage (percent) and sets args._ram_exceeded when the
    # threshold is crossed so the runner can abort safely. Uses psutil when
    # available; otherwise monitoring is disabled with an informational note.
    # Validate max_ram_percent: must be a float in [0,100]
    try:
        max_ram_pct_raw = getattr(args, 'max_ram_percent', None)
        if max_ram_pct_raw is None:
            max_ram_pct = None
        else:
            max_ram_pct = float(max_ram_pct_raw)
            if not (0.0 <= max_ram_pct <= 100.0):
                raise ValueError('--max-ram-percent must be between 0 and 100')
    except Exception as e:
        # If validation fails, stop early with a clear message
        sys.exit(f"ERROR: invalid --max-ram-percent: {e}")

    # cap mode: 'system' (compare vm.percent) or 'run' (compare process RSS/total RAM)
    try:
        max_ram_cap = getattr(args, 'max_ram_cap', 'system')
        if max_ram_cap not in (None, 'system', 'run'):
            raise ValueError('--max-ram-cap must be "system" or "run"')
    except Exception as e:
        sys.exit(f"ERROR: invalid --max-ram-cap: {e}")

    # Normalize target-snps tolerance: require percent (0..100) only; convert to internal fraction
    try:
        raw_tol = getattr(args, 'min_snps_tol', 0.0)
        if raw_tol is None:
            tol_pct = 0.0
        else:
            tol_pct = float(raw_tol)
        # Enforce percent-only input
        if tol_pct < 0.0 or tol_pct > 100.0:
            raise ValueError('--target-snps-tol must be between 0 and 100 (percent)')
        # Store internal fraction (0..1)
        setattr(args, 'min_snps_tol', float(tol_pct) / 100.0)
    except Exception as e:
        sys.exit(f"ERROR: invalid --target-snps-tol: {e}")

    # Normalize sweep position: accept percent in [0,100] and convert to internal fraction [0,1]
    try:
        raw_x = getattr(args, 'x', None)
        if raw_x is None:
            x_pct = None
        else:
            x_pct = float(raw_x)
        if x_pct is not None:
            if x_pct < 0.0 or x_pct > 100.0:
                raise ValueError('--sweep-pos must be between 0 and 100 (percent)')
            # convert percent to fraction
            setattr(args, 'x', float(x_pct) / 100.0)
    except Exception as e:
        sys.exit(f"ERROR: invalid --sweep-pos: {e}")

    def _ram_monitor_loop(stop_event, threshold_pct, args_obj):
        try:
            import psutil
        except Exception:
            # psutil not available; do not monitor
            try:
                sys.stderr.write('# INFO: psutil not available; --max-ram-percent/--max-ram-cap monitoring disabled.\n')
            except Exception:
                pass
            return
        try:
            pfunc = psutil.virtual_memory
        except Exception:
            try:
                sys.stderr.write('# INFO: psutil.virtual_memory not available; RAM monitoring disabled.\n')
            except Exception:
                pass
            return
        # Poll every 0.5s (tunable); when threshold exceeded, set flag on args and exit
        while not stop_event.is_set():
            try:
                vm = pfunc()
                system_used = float(getattr(vm, 'percent', 0.0))

                # Safety cutoff: if system usage reaches or exceeds 99% always trigger
                if system_used >= 99.0:
                    trigger_reason = 'system_critical'
                    trigger_value = system_used
                else:
                    trigger_reason = None
                    trigger_value = None

                # Only evaluate configured threshold if provided and not already critical
                if trigger_reason is None and threshold_pct is not None:
                    cap_mode = getattr(args_obj, 'max_ram_cap', 'system')
                    if cap_mode == 'system' or cap_mode is None:
                        # compare system percent
                        if system_used > float(threshold_pct):
                            trigger_reason = 'system'
                            trigger_value = system_used
                    else:
                        # run mode: compare our process RSS as percent of total RAM
                        try:
                            proc = psutil.Process(os.getpid())
                            rss = float(proc.memory_info().rss)
                            total = float(getattr(vm, 'total', 0) or 0)
                            if total > 0:
                                run_pct = rss / total * 100.0
                                if run_pct > float(threshold_pct):
                                    trigger_reason = 'run'
                                    trigger_value = run_pct
                        except Exception:
                            # If any issue reading process mem, fall back to system check
                            if system_used > float(threshold_pct):
                                trigger_reason = 'system'
                                trigger_value = system_used

                if trigger_reason is not None:
                    # mark on args so the rest of the runner can react
                    try:
                        # Close any active progress bar to avoid tqdm re-rendering
                        try:
                            pb = getattr(args_obj, '_progress_bar', None)
                            if pb is not None:
                                try:
                                    pb.close()
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        try:
                            suggestion_parts = [
                                'reduce the number of parallel workers with --parallel',
                                "if not using msprime, lower --sims-per-work",
                                'reduce --replicates or --sample-individuals',
                                'or rerun with a larger --max-ram-percent / different --max-ram-cap'
                            ]
                            suggestion_text = '; '.join(suggestion_parts)
                            if trigger_reason == 'system' or trigger_reason == 'system_critical':
                                friendly = (
                                    f"System memory usage {trigger_value:.1f}% exceeded the configured system limit of {threshold_pct}% (or hit the 99% safety cutoff).\n"
                                    f"Suggested actions: {suggestion_text}.\n"
                                )
                            else:
                                # run mode
                                try:
                                    # convert threshold to GB for clarity
                                    total = float(getattr(vm, 'total', 0) or 0)
                                    limit_bytes = (float(threshold_pct) / 100.0) * total if total > 0 else None
                                    if limit_bytes is not None:
                                        limit_gb = limit_bytes / (1024**3)
                                        friendly = (
                                            f"Run memory {trigger_value:.1f}% of total RAM exceeded configured limit {threshold_pct}% (~{limit_gb:.2f} GB).\n"
                                            f"Suggested actions: {suggestion_text}.\n"
                                        )
                                    else:
                                        friendly = (f"Run memory {trigger_value:.1f}% exceeded configured limit of {threshold_pct}%.\n")
                                except Exception:
                                    friendly = (f"Run memory {trigger_value:.1f}% exceeded configured limit of {threshold_pct}%.\n")
                        except Exception:
                            friendly = (f"Memory threshold exceeded ({trigger_value:.1f}%).\n")
                        try:
                            setattr(args_obj, '_ram_exceeded', True)
                            setattr(args_obj, '_ram_exceeded_message', friendly)
                            # Also store a module-level message so the signal handler
                            # (which may run in the main thread) can print it.
                            try:
                                globals()['_LAST_RAM_EXCEEDED_MESSAGE'] = friendly
                            except Exception:
                                pass
                            # Best-effort: request our own process to terminate so the
                            # registered _signal_handler runs and can print the message.
                            try:
                                os.kill(os.getpid(), signal.SIGTERM)
                            except Exception:
                                pass
                        except Exception:
                            pass
                        # Best-effort: try to terminate child processes we started
                        try:
                            _terminate_children()
                        except Exception:
                            pass
                        # Do not print the friendly message here; let the top-level
                        # handler print it once to avoid duplicate messages.
                    except Exception:
                        pass
                    # stop monitoring and exit
                    break
            except Exception:
                pass
            stop_event.wait(0.5)

    # Launch monitor if a numeric threshold is present
    try:
        if max_ram_pct is not None:
            ram_ev = threading.Event()
            ram_th = threading.Thread(target=_ram_monitor_loop, args=(ram_ev, max_ram_pct, args), daemon=True)
            try:
                setattr(args, '_ram_monitor_event', ram_ev)
                setattr(args, '_ram_monitor_thread', ram_th)
                setattr(args, '_ram_exceeded', False)
                ram_th.start()
            except Exception:
                try:
                    sys.stderr.write('# INFO: Failed to start RAM monitor thread; continuing without enforcement.\n')
                except Exception:
                    pass
    except Exception:
        pass
    argv_full = sys.argv[1:]
    # Misuse checks
    # Allow msms to accept --sweep-time (origin time) as an alternative to --fixation-time.
    # Required selection params
    if engine == 'discoal':
        # Neutral discoal run permitted when all selection params omitted.
        # If the user requests selection we require --sel-s and --sweep-pos
        # and at least one of --sweep-time or --fixation-time.
        user_attempted_selection = (getattr(args, 's', None) is not None) or (args.x is not None) or (args.sweep_time is not None) or (getattr(args, 'fix_time', None) is not None)
        if user_attempted_selection:
            # require sel-s and sweep-pos
            if getattr(args, 's', None) is None or args.x is None:
                raise SystemExit('ERROR: discoal selection requires --sel-s and --sweep-pos.')
            # require at least one of sweep-time or fixation-time
            has_sweep_time = args.sweep_time is not None
            has_fix_time = getattr(args, 'fix_time', None) is not None
            if not (has_sweep_time or has_fix_time):
                raise SystemExit('ERROR: discoal selection requires --sweep-time or --fixation-time.')
            # Validate non-negativity of provided time value(s)
            if has_sweep_time and getattr(args, 'sweep_time', None) is not None and args.sweep_time < 0:
                sys.exit('ERROR: --sweep-time must be >= 0.')
            if has_fix_time and getattr(args, 'fix_time', None) is not None and args.fix_time < 0:
                sys.exit('ERROR: --fixation-time must be >= 0.')
    elif engine == 'msms':
        # msms selection: if user attempts selection require --sel-s and --sweep-pos
        # and exactly one of --sweep-time or --fixation-time.
        user_provided_fix = any(a.startswith('--fixation-time') for a in sys.argv[1:])
        user_attempted_selection = (getattr(args, 's', None) is not None) or (args.x is not None) or user_provided_fix
        if user_attempted_selection:
            if getattr(args, 's', None) is None or args.x is None:
                sys.exit('ERROR: msms selection requires both --sel-s and --sweep-pos.')
            has_sweep_time = args.sweep_time is not None
            has_fix_time = getattr(args, 'fix_time', None) is not None
            if not (has_sweep_time ^ has_fix_time):
                sys.exit('ERROR: msms selection requires exactly one of --sweep-time or --fixation-time (they are mutually exclusive).')
            if getattr(args, 'fix_time', None) is not None and getattr(args, 'fix_time', None) < 0:
                sys.exit('ERROR: --fixation-time must be >= 0 (time back in gens or 4N units as specified).')
            # Also validate sweep_time when provided for msms
            if args.sweep_time is not None and args.sweep_time < 0:
                sys.exit('ERROR: --sweep-time must be >= 0.')
    else:  # ms or other neutral-like engines
        sel_flags = ['--sweep-pos', '--sel-s', '--sweep-time', '--fixation-time', '--sweep-pop']
        provided = [fl for fl in sel_flags if any(a == fl or a.startswith(fl + '=') for a in argv_full)]
        if provided:
            # some neutral engines (ms, scrm) ignore selection flags
            sys.stderr.write('# INFO: Ignoring selection flags ' + ', '.join(provided) + " because engine in ('ms','scrm') is neutral only. Use --engine discoal or --engine msms for selection.\n")
    # target-snps tolerance dependency (no hard error: just ignore if target absent)
    tol_flag_provided = any(a == '--target-snps-tol' or a.startswith('--target-snps-tol=') for a in argv_full)
    if args.min_snps is None and tol_flag_provided:
        # Inform user and neutralize tolerance value
        sys.stderr.write('# INFO: Ignoring --target-snps-tol (no --target-snps provided). It only applies when deriving/ enforcing target SNP count.\n')
        setattr(args, 'min_snps_tol', 0.0)
    # SFS dependency (convert prior hard error into informational ignore)
    have_sfs = bool(getattr(args, 'sfs', False))
    sfs_norm_flag = any(a == '--sfs-normalized' or a.startswith('--sfs-normalized=') for a in argv_full)
    sfs_mode_flag = any(a == '--sfs-mode' or a.startswith('--sfs-mode=') for a in argv_full)
    if not have_sfs and (sfs_norm_flag or sfs_mode_flag):
        sys.stderr.write('# INFO: Ignoring --sfs-normalized/--sfs-mode (require --sfs). Provide --sfs to enable SFS calculation.\n')
        # Leave defaults; ensure normalization flag off
        if hasattr(args, 'sfs_normalized'):
            setattr(args, 'sfs_normalized', False)
    # Validate chromosome argument (if provided) against stdpopsim for the given species
    chr_name = getattr(args, 'chromosome', None)
    if chr_name:
        try:
            sp = sps.get_species(args.species)
        except Exception:
            # Try resolving a full species name (e.g. 'Homo sapiens') to its ID
            resolved = None
            try:
                name_to_id = {sp_obj.name: sp_obj.id for sp_obj in sps.all_species()}
                # Direct match first
                if args.species in name_to_id:
                    resolved = name_to_id[args.species]
                else:
                    # Case-insensitive match
                    wanted = args.species.strip().lower()
                    for full_name, sp_id in name_to_id.items():
                        if full_name.lower() == wanted:
                            resolved = sp_id
                            break
                if resolved:
                    args.species = resolved
                    sp = sps.get_species(args.species)
                else:
                    raise SystemExit(f"ERROR: Cannot validate --chromosome because species '{args.species}' not found in stdpopsim.")
            except SystemExit:
                raise
            except Exception:
                raise SystemExit(f"ERROR: Cannot validate --chromosome because species '{args.species}' not found in stdpopsim.")
        # Try to resolve contig using stdpopsim API; if not found, attempt to list available contig names
        contig_ok = False
        try:
            c = sp.get_contig(chr_name)
            if c is not None:
                contig_ok = True
        except Exception:
            contig_ok = False
        if not contig_ok:
            # Attempt to enumerate contigs to provide helpful feedback
            contig_names = []
            try:
                genome = getattr(sp, 'genome', None)
                if genome is not None and hasattr(genome, 'contigs'):
                    contig_names = [getattr(cc, 'name', str(cc)) for cc in genome.contigs]
            except Exception:
                contig_names = []
            if contig_names:
                sample = ', '.join(contig_names[:50])
                more = '...' if len(contig_names) > 50 else ''
                raise SystemExit(f"ERROR: Chromosome '{chr_name}' not found for species '{args.species}'. Available contigs: {sample}{more}")
            else:
                raise SystemExit(f"ERROR: Chromosome '{chr_name}' not found for species '{args.species}' in stdpopsim.")
    return args

# ----------------------- Simulation Orchestration -----------------------
def _pre_run_setup(args):
    """Pre-run hook: set temp dir preference from args.temp_loc.

    Behaviour:
    - default 'local' to avoid RAM overflow
    - 'system' uses Python's system tempdir (e.g. /tmp via tempfile.gettempdir())
    """
    try:
        pref = getattr(args, 'temp_loc', 'local')
        if pref == 'local':
            os.environ['SIMULATOR_USE_REPO_TMP'] = '1'
        else:  # 'system'
            os.environ.pop('SIMULATOR_USE_REPO_TMP', None)
    # Do not set debug or preserve env flags; temporary files are not preserved
    except Exception:
        pass
    return None

def _post_run_teardown(args, result):
    """Post-run hook: perform cleanup of local temporary directory."""
    try:
        _cleanup_local_tmpdir()
    except Exception:
        pass
    # Stop any background RAM monitor thread started in _run_simulation_core
    try:
        ev = getattr(args, '_ram_monitor_event', None)
        th = getattr(args, '_ram_monitor_thread', None)
        if ev is not None:
            try:
                ev.set()
            except Exception:
                pass
        if th is not None and isinstance(th, threading.Thread):
            try:
                th.join(timeout=1.0)
            except Exception:
                pass
    except Exception:
        pass
    return result

def run_simulation_from_args(args):
    """Thin orchestrator delegating to the original core implementation."""
    _pre_run_setup(args)
    result = _run_simulation_core(args)
    return _post_run_teardown(args, result)

def _run_simulation_core(args):
    # Track whether user explicitly provided --output (before we potentially assign a default)
    setattr(args, '_user_out_provided', args.out is not None and args.out != '')
    user_order = [d.strip() for d in args.demes.split(',') if d.strip()]

    # Automatic primary seed (used for reproducibility and optionally neutral run)
    # Primary seed retained only for optional reuse; individual parallel chunks get unique seeds later.
    try:
        import random as _rnd
        primary_seed = _rnd.randint(1, 2**31 - 1)
    except Exception:
        primary_seed = None

    # Resolve species name -> id if needed
    try:
        _species_map = {sp.name: sp.id for sp in sps.all_species()}
    except Exception:
        _species_map = {}
    if _species_map:
        parts = args.species.strip().split()
        if len(parts) == 2:
            args.species = parts[0].title() + ' ' + parts[1].lower()
        if args.species in _species_map.values():
            pass
        elif args.species in _species_map:
            args.species = _species_map[args.species]
        else:
            sys.exit(f"ERROR: Unknown species '{args.species}'.")

    # parallel worker count: cap differs by engine preference
    try:
        max_cpu = os.cpu_count() or 1
    except Exception:
        max_cpu = 1
    # If using msprime prefer physical cores as the parallel limit; else use logical CPUs
    cap_cpu = max_cpu
    try:
        if getattr(args, 'engine', None) == 'msprime':
            try:
                import psutil
                phys = psutil.cpu_count(logical=False) or max_cpu
            except Exception:
                phys = max_cpu
            cap_cpu = phys
    except Exception:
        cap_cpu = max_cpu
    args.parallel = min(max(1, args.parallel), cap_cpu)
    # Do not allow more parallel workers than total replicates requested
    try:
        args.parallel = min(int(args.parallel), max(1, int(getattr(args, 'reps', 1))))
    except Exception:
        pass
    # If the user specified a chromosome/contig but DID NOT provide an explicit
    # --length or --target-snps, prefer using the contig's native length. If the
    # user has explicitly supplied --length or --target-snps we respect that and
    # do NOT force the full contig length.
    if getattr(args, 'chromosome', None) is not None and (getattr(args, 'length', None) is None and getattr(args, 'min_snps', None) is None):
        try:
            sp_tmp = sps.get_species(args.species)
            try:
                cont = sp_tmp.get_contig(args.chromosome)
            except Exception:
                cont = None
            if cont is not None and hasattr(cont, 'length'):
                try:
                    args.length = int(getattr(cont, 'length'))
                except Exception:
                    # if we cannot read contig length, leave args.length unchanged and let builders error naturally
                    pass
        except Exception:
            # best-effort only; do not fatal here
            pass

    # validate model & complete demes list
    try:
        species_obj = sps.get_species(args.species)
    except Exception as e:
        sys.exit(f"ERROR: Unknown species '{args.species}': {e}")
    valid_models = [m.id for m in species_obj.demographic_models]
    if args.model not in valid_models:
        sys.exit(f"ERROR: Unknown model '{args.model}'. Available models: {','.join(valid_models)}")
    try:
        model_obj = species_obj.get_demographic_model(args.model)
        raw_demog_m = model_obj.model
        try: graph_m = raw_demog_m.to_demes()
        except AttributeError: graph_m = raw_demog_m
        model_demes = [d.name for d in getattr(graph_m, 'demes', [])]
    except Exception:
        model_demes = []
    missing = []
    if model_demes:
        unknown_user = [d for d in user_order if d not in model_demes]
        if unknown_user:
            sys.exit(f"ERROR: Unknown deme(s) {unknown_user}. Available: {','.join(model_demes)}")
        missing = [d for d in model_demes if d not in user_order]
        if missing:
            user_order = user_order + missing

    # parse samples
    raw_sample_tokens = [t.strip() for t in args.samples_config.split(',') if t.strip()!='']
    try:
        individual_counts = [int(x) for x in raw_sample_tokens]
    except ValueError:
        sys.exit('ERROR: --sample-individuals must be a comma list of integers.')
    if len(individual_counts) > len(user_order):
        sys.exit('ERROR: More sample counts provided than demes.')
    padded = False
    if len(individual_counts) < len(user_order):
        individual_counts += [0] * (len(user_order) - len(individual_counts))
        padded = True

    engine = args.engine

    # sweep flag validation (engine specific)
    if engine == 'discoal':
        # Treat as selection only if user supplied at least one of the sweep arguments.
        # Accept either --sweep-time (origin time) or --fixation-time as a valid time input.
        user_attempted_selection = (
            (getattr(args, 's', None) is not None) or
            (getattr(args, 'x', None) is not None) or
            (getattr(args, 'sweep_time', None) is not None) or
            (getattr(args, 'fix_time', None) is not None)
        )
        if user_attempted_selection:
            # require sel-s and sweep-pos
            if getattr(args, 's', None) is None or args.x is None:
                sys.exit('ERROR: discoal sweep requires --sel-s and --sweep-pos.')
            # require at least one of sweep-time or fixation-time
            has_sweep_time = getattr(args, 'sweep_time', None) is not None
            has_fix_time = getattr(args, 'fix_time', None) is not None
            if not (has_sweep_time or has_fix_time):
                sys.exit('ERROR: discoal sweep requires --sweep-time or --fixation-time.')
            # validate non-negativity for whichever time was provided
            if has_sweep_time and getattr(args, 'sweep_time', None) is not None and args.sweep_time < 0:
                sys.exit('ERROR: --sweep-time must be >= 0.')
            if has_fix_time and getattr(args, 'fix_time', None) is not None and args.fix_time < 0:
                sys.exit('ERROR: --fixation-time must be >= 0.')
        else:
            sys.stderr.write('# INFO: No sweep parameters supplied; running neutral discoal simulation.\n')
    elif engine == 'msms':
        # Because fix_time has a default (0.0), do not treat its mere presence as a selection attempt.
        user_provided_fix = any(a.startswith('--fixation-time') for a in sys.argv[1:])
        user_attempted_selection = (getattr(args, 's', None) is not None) or (args.x is not None) or user_provided_fix
        if user_attempted_selection:
            if getattr(args, 's', None) is None or args.x is None:
                sys.exit('ERROR: msms sweep requires both --sel-s and --sweep-pos.')
            if args.fix_time is not None and args.fix_time < 0:
                sys.exit('ERROR: --fixation-time must be >= 0 (time back in 4N units).')
        else:
            # engines that are neutral-only (ms, scrm) — selection flags are ignored
            sys.stderr.write('# INFO: No sweep parameters supplied; running neutral simulation (engine neutral-only).\n')

    # Adaptive target-S segregating sites refinement if user gave --target-snps without explicit --length
    def _pilot_build(engine_name, length_val, reps_override=1, min_snps_forward=None):
        if engine_name == 'discoal':
            cmd, meta_local = build_discoal_command(
                    species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    discoal_pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                    max_fold_per_step=args.max_fold_per_step, sweep_time=(getattr(args,'sweep_time', None) if getattr(args,'sweep_time', None) is not None else getattr(args,'fix_time', None)),
                    x=getattr(args,'x',None), s=getattr(args,'s',None), min_snps=min_snps_forward, chr_name=args.chromosome,
                    disable_en_ladder=args.no_en_ladder, seed=None, time_units=getattr(args,'time_units','gens'),
                )
        elif engine_name == 'msms':
            cmd, meta_local = build_msms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=min_snps_forward, disable_en_ladder=args.no_en_ladder,
                a=getattr(args,'a',None), s=getattr(args,'s',None), x=getattr(args,'x',None), fixation_time=(getattr(args,'fix_time', None) if getattr(args,'fix_time', None) is not None else None), sweep_time=(getattr(args,'sweep_time', None) if getattr(args,'sweep_time', None) is not None else getattr(args,'fix_time', None)), seed=None, time_units=getattr(args,'time_units','gens'),
            )
        elif engine_name == 'scrm':
            cmd, meta_local = build_scrm_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=min_snps_forward, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        elif engine_name == 'msprime':
            cmd, meta_local = build_msprime_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=min_snps_forward, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        else:  # ms
            cmd, meta_local = build_ms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=min_snps_forward, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        return cmd, meta_local

    # Adaptive target-S refinement involves executing pilot simulations. Skip when only showing command.
    if args.length is None and args.min_snps is not None and not getattr(args, 'show_command', False):
        # Initial builder-derived length (using existing formula logic in builders)
        init_cmd, init_meta = _pilot_build(engine, None, reps_override=1, min_snps_forward=args.min_snps)
        candidate_len = init_meta['length']
        target_snps = args.min_snps
        max_iter = 6
        last_obs = None
        for it in range(max_iter):
            pilot_cmd, _m = _pilot_build(engine, candidate_len, reps_override=1, min_snps_forward=None)
            raw_line = None
            for l in pilot_cmd.splitlines():
                if l.startswith(engine+' '):
                    raw_line = l; break
            if raw_line is None:
                break
            toks = shlex.split(raw_line)
            try:
                # Pilot runs should be limited to a small CPU set to avoid oversubscription
                total_cpus = os.cpu_count() or 1
                # assign the pilot run to CPUs 0..min(1,total_cpus-1)
                pilot_cpus = list(range(min(1, total_cpus))) if total_cpus > 0 else [0]
                pre = _make_affinity_preexec(pilot_cpus)
                r = _run_and_register(toks, capture_output=True, text=True, check=True, preexec_fn=pre)
                out_txt = r.stdout
            except Exception:
                break
            m = re.search(r'^segsites:\s*(\d+)', out_txt, re.MULTILINE)
            segs = int(m.group(1)) if m else 0
            last_obs = segs
            if segs >= target_snps:
                break
            if segs == 0:
                candidate_len = int(candidate_len * 2)
            else:
                factor = target_snps / float(segs)
                inflation = 1.10 if engine in ('ms','scrm') else (1.35 if engine == 'msms' else 1.50)
                candidate_len = int(math.ceil(candidate_len * factor * inflation))
        # Optional shrink if overshoot beyond internal tolerance
        tol = getattr(args,'min_snps_tol', 0.0)
        if tol > 0 and last_obs is not None and last_obs > target_snps * (1.0 + tol):
            shrink_iters = 3
            for _ in range(shrink_iters):
                if last_obs <= target_snps * (1.0 + tol):
                    break
                ratio = target_snps / float(last_obs)
                candidate_len = max(100, int(candidate_len * ratio * 0.98))
                pilot_cmd, _m = _pilot_build(engine, candidate_len, reps_override=1, min_snps_forward=None)
                raw_line = None
                for l in pilot_cmd.splitlines():
                    if l.startswith(engine+' '):
                        raw_line = l; break
                if raw_line is None:
                    break
                toks = shlex.split(raw_line)
                try:
                    # Pilot shrink runs also limited to small CPU set
                    total_cpus = os.cpu_count() or 1
                    pilot_cpus = list(range(min(1, total_cpus))) if total_cpus > 0 else [0]
                    pre = _make_affinity_preexec(pilot_cpus)
                    r = _run_and_register(toks, capture_output=True, text=True, check=True, preexec_fn=pre)
                    out_txt = r.stdout
                except Exception:
                    break
                m = re.search(r'^segsites:\s*(\d+)', out_txt, re.MULTILINE)
                last_obs = int(m.group(1)) if m else 0
                if last_obs < target_snps:
                    candidate_len = int(math.ceil(candidate_len * (target_snps / max(1,last_obs)) * 1.05))
                    break
        refined_len = candidate_len
    # refined length computed; proceed
        args.length = refined_len

    # build command (final, using possibly refined args.length unless in debug mode which skips refinement)
    # Inform about growth-discretization flags with msprime (no effect across the pipeline)
    try:
        if engine == 'msprime' and (getattr(args, 'no_en_ladder', False) or (float(getattr(args, 'max_fold_per_step', 1.05)) != 1.05)):
            sys.stderr.write('# INFO: --disable-growth-discretization/--growth-max-fold have no effect for engine=msprime (continuous growth simulated natively).\n')
    except Exception:
        pass
    if engine == 'discoal':
        sim_cmd, meta = build_discoal_command(
            species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
            discoal_pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
            max_fold_per_step=args.max_fold_per_step, sweep_time=(getattr(args,'sweep_time',None) if getattr(args,'sweep_time',None) is not None else getattr(args,'fix_time', None)),
            x=getattr(args,'x',None), s=getattr(args,'s',None), min_snps=args.min_snps, chr_name=args.chromosome,
            disable_en_ladder=args.no_en_ladder, seed=None, time_units=getattr(args,'time_units','gens'),
        )
    elif engine == 'msms':
        sim_cmd, meta = build_msms_command(
            species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
            pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
            max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
            min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder,
            a=getattr(args,'a',None), s=getattr(args,'s',None), x=getattr(args,'x',None), fixation_time=(getattr(args,'fix_time', None) if getattr(args,'fix_time', None) is not None else None), sweep_time=(getattr(args,'sweep_time', None) if getattr(args,'sweep_time', None) is not None else getattr(args,'fix_time', None)), seed=None, time_units=getattr(args,'time_units','gens'),
        )
    else:
        # ms neutral only: selection flags were already warned & ignored in parse_args.
        if engine == 'scrm':
            sim_cmd, meta = build_scrm_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=args.length,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        elif engine == 'msprime':
            sim_cmd, meta = build_msprime_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=args.length,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        else:
            sim_cmd, meta = build_ms_command(
            species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
            # For engine=ms there is no --sweep-pop/--sweep-pos group; pop0 falls back to first deme.
            pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=args.length,
            max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
            min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, seed=None,
        )

    # Show-only mode: print the external command and exit (do not execute). For msprime, no external command.
    if getattr(args, 'show_command', False):
        if engine == 'msprime':
            sys.stderr.write('# ERROR: engine=msprime is an internal stdpopsim/msprime engine; there is no external shell command to show.\n')
            sys.exit(1)
        raw_line = None
        for ln in sim_cmd.splitlines():
            if ln.strip().startswith(engine + ' '):
                raw_line = ln.strip(); break
        if raw_line is None:
            raw_line = sim_cmd.strip()
        sys.stdout.write(raw_line + ('\n' if not raw_line.endswith('\n') else ''))
        return

    # ---------------- Simulation runner + enforcement helpers ----------------
    def _simulate_current_length(length_val, show_progress=True):
        """Build and execute simulation at given length. Returns meta, raw_cmd_line, stdout, stderr."""
        # rebuild command for current length (meta updated)
        if engine == 'discoal':
            sim_cmd_loc, meta_loc = build_discoal_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                discoal_pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, sweep_time=(getattr(args, 'sweep_time', None) if getattr(args,'sweep_time', None) is not None else getattr(args,'fix_time', None)),
                x=getattr(args, 'x', None), s=getattr(args, 's', None), min_snps=args.min_snps, chr_name=args.chromosome,
                disable_en_ladder=args.no_en_ladder, seed=None,
            )
        elif engine == 'msms':
            sim_cmd_loc, meta_loc = build_msms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder,
                a=getattr(args, 'a', None), s=getattr(args, 's', None), x=getattr(args, 'x', None), fixation_time=(getattr(args, 'fix_time', None) if getattr(args,'fix_time', None) is not None else None), sweep_time=(getattr(args,'sweep_time', None) if getattr(args,'sweep_time', None) is not None else getattr(args,'fix_time', None)), seed=None,
            )
        elif engine == 'scrm':
            sim_cmd_loc, meta_loc = build_scrm_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        elif engine == 'msprime':
            sim_cmd_loc, meta_loc = build_msprime_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        else:
            sim_cmd_loc, meta_loc = build_ms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, seed=None,
            )
        raw_cmd_line_loc = None
        for ln in sim_cmd_loc.splitlines():
            if ln.startswith(engine + ' '):
                raw_cmd_line_loc = ln.strip(); break
        if raw_cmd_line_loc is None:
            sys.exit('ERROR: Failed to locate simulator command line (builder).')
        # No legacy debug bypass; always execute unless --show-command was used earlier.
        cmd_tokens_loc = shlex.split(raw_cmd_line_loc)
        total_reps_loc = int(cmd_tokens_loc[2])
        # progress handling (optionally suppressed for enforcement iterations)
        # Honor the user-requested parallel worker count (bounded by total replicates).
        # This controls the number of concurrent Python worker threads used to
        # run replicate chunks. Use the CLI --parallel as-is for all engines.
        threads_loc = max(1, min(int(getattr(args, 'parallel', 1)), total_reps_loc))
        per_rep_mode_loc = bool(args.progress and show_progress and total_reps_loc > 1)

        # Ignore --sims-per-work for msprime (internal engine handles batching itself)
        if engine == 'msprime' and getattr(args, 'sims_per_work', None) is not None:
            try:
                sys.stderr.write('# INFO: --sims-per-work ignored for engine=msprime (internal batching).\n')
            except Exception:
                pass
            args.sims_per_work = None

        # When --sims-per-work is not provided, we will split total replicates
        # into exactly `threads_loc` parts (the requested parallel workers).
        # The parts need not be equal: we use floor division and add the
        # remainder to the last chunk. For msprime, --sims-per-work is ignored
        # and chunks are per-replicate.
        max_allowed_sims_per_work = total_reps_loc  # absolute upper bound
        if engine != 'msprime' and getattr(args, 'sims_per_work', None) is None:
            try:
                base_spw = max(1, total_reps_loc // threads_loc)
            except Exception:
                base_spw = 1
            try:
                rem = total_reps_loc - base_spw * threads_loc
                if rem > 0:
                    sys.stderr.write(f"# INFO: --sims-per-work not provided; splitting {total_reps_loc} replicates into {threads_loc} parts (base {base_spw}, remainder {rem} added to last part).\n")
                else:
                    sys.stderr.write(f"# INFO: --sims-per-work not provided; splitting {total_reps_loc} replicates into {threads_loc} equal parts of {base_spw}.\n")
            except Exception:
                pass
        # If user provided a value, validate and cap it (msprime remains None)
        if engine != 'msprime' and getattr(args, 'sims_per_work', None) is not None:
            try:
                user_spw = int(args.sims_per_work)
                if user_spw < 1:
                    sys.stderr.write(f"# INFO: --sims-per-work {user_spw} < 1; using 1 instead.\n")
                    user_spw = 1
                if user_spw > max_allowed_sims_per_work:
                    sys.stderr.write(f"# INFO: --sims-per-work {user_spw} > total replicates {max_allowed_sims_per_work}; capping to {max_allowed_sims_per_work}.\n")
                    user_spw = max_allowed_sims_per_work
                args.sims_per_work = user_spw
            except Exception:
                # If parsing/capping fails, leave as None so the default
                # splitting behaviour (into --parallel parts) is used.
                args.sims_per_work = None

        # Now produce the final chunks list. If per-rep mode, keep per-replicate
        # chunks. Otherwise split total_reps_loc into pieces of size
        # args.sims_per_work (last chunk may be smaller). This allows more
        # chunks than --parallel; the executor will limit concurrency.
        # Chunking: do not force per-replicate chunks just because progress is enabled.
        # Engine msprime always uses per-replicate chunks because it runs in-process.
        if engine == 'msprime':
            chunks_loc = [1] * total_reps_loc
        else:
                # If sims_per_work is None (user omitted), split into threads_loc parts
                # with remainder added to the last part. Otherwise use the provided size.
                spw_val = getattr(args, 'sims_per_work', None)
                if spw_val is None:
                    try:
                        base = max(1, total_reps_loc // threads_loc)
                    except Exception:
                        base = 1
                    rem = total_reps_loc - base * threads_loc
                    chunks_loc = [base] * threads_loc
                    if rem > 0:
                        # Add remainder to the last chunk
                        chunks_loc[-1] = chunks_loc[-1] + rem
                else:
                    spw = max(1, int(spw_val))
                    if spw >= total_reps_loc:
                        chunks_loc = [total_reps_loc]
                    else:
                        full = total_reps_loc // spw
                        rem = total_reps_loc % spw
                        chunks_loc = [spw] * full
                        if rem:
                            chunks_loc.append(rem)
        # Prepare unique seeds per chunk for engines supporting explicit seeding
        chunk_process_count = None
        # Engines that support per-replicate seeding/mode. Include msprime so it can run multiple
        # replicates with distinct seeds in per-rep mode like other in-process engines.
        per_rep_engine = engine in ('discoal', 'msms', 'msprime') and per_rep_mode_loc
        # For in-process msprime runs, respect the thread caps that will be
        # assigned per worker below; do not force a global msprime internal
        # threads override here (we set caps per worker when launching them).
        # Generate unique seeds per chunk for external engines and msprime; for ms we will use
        # the Hudson ms "-seeds" triplet to ensure independence across parallel chunks.
        if engine in ('discoal', 'msms', 'msprime', 'ms'):
            if per_rep_mode_loc:
                chunk_process_count = len(chunks_loc)  # each chunk = 1 replicate
            else:
                chunk_process_count = len(chunks_loc)
            try:
                unique_seeds = set()
                seeds_for_chunks = []
                import random as _rnd2
                while len(seeds_for_chunks) < chunk_process_count:
                    s = _rnd2.randint(1, 2**31 - 1)
                    if s in unique_seeds:
                        continue
                    unique_seeds.add(s)
                    seeds_for_chunks.append(s)
            except Exception:
                seeds_for_chunks = [None] * chunk_process_count
        else:
            seeds_for_chunks = [None] * len(chunks_loc)

        def run_chunk_loc(rc, replicate_index=None, chunk_idx=None):
            toks = cmd_tokens_loc[:]; toks[2] = str(rc)
            # Prepare a default environment for subprocesses early so all branches
            # can reference `env` safely. We default BLAS/OpenMP thread caps to 1
            # here to avoid accidental oversubscription for external engines.
            env = os.environ.copy()
            try:
                for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
                    env[k] = env.get(k, '1')
            except Exception:
                pass
            # Ensure the current Python interpreter's bin directory is on PATH so
            # executables installed into the same conda env (e.g., msms, scrm)
            # are discoverable by subprocess.run().
            try:
                exec_dir = os.path.dirname(sys.executable)
                if exec_dir and exec_dir not in (env.get('PATH') or ''):
                    env['PATH'] = exec_dir + os.pathsep + (env.get('PATH') or '')
            except Exception:
                pass
            # Ensure preexec_fn variable exists for all branches (may be set later).
            pre = None
            # If the engine supports explicit per-chunk seeding, attach the seed for reproducibility.
            seed_val = None
            if seeds_for_chunks and chunk_idx is not None and chunk_idx < len(seeds_for_chunks):
                seed_val = seeds_for_chunks[chunk_idx]
                # Build contig and samples mapping
                try:
                    meta_lp = meta_loc
                    model_obj = meta_lp.get('model_obj')
                    species_id = meta_lp.get('species_id')
                    sp = sps.get_species(species_id)
                    chr_name = meta_lp.get('chromosome') or getattr(args, 'chromosome', None)
                    # Fallback: if no chromosome provided, pick the first available contig name
                    if not chr_name:
                        try:
                            genome_obj = getattr(sp, 'genome', None)
                            contigs_list = getattr(genome_obj, 'contigs', []) if genome_obj is not None else []
                            if contigs_list:
                                # Prefer .name, else .id, else string repr
                                cand = getattr(contigs_list[0], 'name', None) or getattr(contigs_list[0], 'id', None) or str(contigs_list[0])
                                if isinstance(cand, str) and cand.strip():
                                    chr_name = cand
                        except Exception:
                            # leave chr_name as None; subsequent get_contig will likely fail and be reported below
                            pass
                    # Build a contig spanning [0, length_val) with model's mutation_rate where possible
                    try:
                        chosen_contig = sp.get_contig(chr_name, mutation_rate=getattr(model_obj, 'mutation_rate', None), left=0, right=max(0, int(length_val) - 1))
                    except Exception:
                        # fallback to requesting contig without bounds
                        try:
                            chosen_contig = sp.get_contig(chr_name)
                        except Exception:
                            chosen_contig = None
                    # If still None, report a clear error early
                    if chosen_contig is None:
                        avail_names = []
                        try:
                            genome_obj = getattr(sp, 'genome', None)
                            contigs_list = getattr(genome_obj, 'contigs', []) if genome_obj is not None else []
                            avail_names = [getattr(c, 'name', getattr(c, 'id', str(c))) for c in contigs_list]
                        except Exception:
                            pass
                        hint = ("; available: " + ", ".join(map(str, avail_names[:20])) + ("..." if len(avail_names) > 20 else "")) if avail_names else ""
                        return (rc, '', f"# ERROR: msprime engine could not resolve a contig. Provide --chromosome or choose one{hint}.\n", 1, replicate_index)
                    # no debug printing
                    # if model provides a mutation rate and chosen_contig exists, try to set it to avoid stdpopsim warning
                    try:
                        if chosen_contig is not None and getattr(model_obj, 'mutation_rate', None) is not None:
                            try:
                                chosen_contig.mutation_rate = float(getattr(model_obj, 'mutation_rate'))
                            except Exception:
                                pass
                    except Exception:
                        pass
                    # Prepare samples per-population as individuals (use counts from meta if available)
                    samples_per_pop = {}
                    pop_order = meta_lp.get('pop_order') or []
                    counts = meta_lp.get('counts_disc') or []
                    pl = meta_lp.get('ploidy', 1) or 1
                    # counts_disc are haplotypes; convert to individuals
                    for i, pop in enumerate(pop_order):
                        try:
                            hapc = counts[i]
                            samples_per_pop[pop] = int(hapc // pl)
                        except Exception:
                            samples_per_pop[pop] = 0
                except Exception as e:
                    return (rc, '', f'# ERROR: failed to prepare msprime inputs: {e}\n', 1, replicate_index)
                # Run rc replicates, collect ms-like output blocks
                out_acc = []
                err_acc = []
                # Ensure we respect the requested replicate count. Prefer the rc passed into
                # run_chunk_loc (which is the number of replicates for this chunk). Fall back
                # to the builder-provided 'reps' metadata if rc is None or invalid.
                try:
                    rep_count = int(rc) if rc is not None else int(meta_lp.get('reps', 1))
                except Exception:
                    rep_count = int(meta_lp.get('reps', 1) or 1)
                if engine == 'msprime':
                    # Build JSON payload for worker and run worker.py in a subprocess.
                    payload = {
                        'species': meta_loc.get('species_id') or args.species,
                        'model_id': meta_loc.get('model_id') or args.model,
                        'chromosome': meta_loc.get('chromosome') or args.chromosome,
                        'length': int(meta_loc.get('length') or length_val or 0),
                        'pop_order': meta_loc.get('pop_order') or [],
                        'counts': meta_loc.get('counts_disc') or [],
                        'ploidy': meta_loc.get('ploidy') or 1,
                        'rep_count': rep_count,
                        'seed': seed_val,
                    }
                    try:
                        payload['suppress_worker_sfs'] = (not bool(getattr(args, 'sfs', False)))
                    except Exception:
                        payload['suppress_worker_sfs'] = False
                    # Allow caller to request that workers do not emit per-rep SFS markers
                    # Workers should produce allele-frequency-spectrum markers (#SFS_AF);
                    # do not ask workers to suppress SFS emission.
                    payload['suppress_worker_sfs'] = False
                    try:
                        # Per-worker threading: allocate native threads based on
                        # physical CPU cores (not logical/threaded cores). This
                        # avoids hyperthreading oversubscription when running
                        # multiple msprime workers in parallel. We compute a
                        # fair share of physical cores per concurrently-running
                        # worker and set OpenMP/BLAS caps accordingly.
                        try:
                            import psutil
                            phys_total = psutil.cpu_count(logical=False) or (os.cpu_count() or 1)
                        except Exception:
                            phys_total = os.cpu_count() or 1
                        req_workers = max(1, int(getattr(args, 'parallel', 1)))
                        workers = max(1, min(req_workers, (chunk_process_count or 1)))
                        base = max(1, int(phys_total) // workers)
                        rem = int(phys_total) % workers
                        if chunk_idx is not None:
                            per_worker_threads = base + (1 if (chunk_idx % workers) < rem else 0)
                        else:
                            per_worker_threads = base
                        env_worker = os.environ.copy()
                        for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
                            env_worker[k] = str(per_worker_threads)
                        try:
                            out_txt, err_txt, rc = _msprime_worker_run(payload, env=env_worker)
                            out_acc.append(out_txt)
                            if rc != 0:
                                err_acc.append(err_txt or f'# ERROR: msprime worker exit code {rc}\n')
                        except Exception as e:
                            err_acc.append(f'# ERROR: failed to run msprime worker: {e}\n')
                    except Exception as e:
                        err_acc.append(f'# ERROR: msprime worker setup failed: {e}\n')
                else:
                    # Non-msprime engines: run the external process and stream its output
                    # to a per-worker output file to avoid returning large strings to the parent.
                    try:
                        # Determine worker tmpdir
                        try:
                            only_sfs_local = bool(getattr(args, '_suppress_primary_output', False) or (getattr(args, 'sfs', False) and not getattr(args, '_user_out_provided', False)))
                        except Exception:
                            only_sfs_local = False
                        # Do not preserve per-worker tmpdirs automatically; they are removed after use.
                        worker_tmp = env.get('SIMULATOR_WORKER_TMP') or _ensure_local_tmpdir()
                        try:
                            os.makedirs(worker_tmp, exist_ok=True)
                        except Exception:
                            pass
                        out_path = None
                        try:
                            # create a unique output file per chunk
                            out_name = f"worker_out_{os.getpid()}_{int(time.time())}_{random.randrange(10**6)}.out"
                            out_path = os.path.join(worker_tmp, out_name)
                            # Redirect child stdout/stderr directly to files to avoid
                            # any parent-side buffering. This is the most memory-safe
                            # approach: the parent never receives the child's stdout
                            # bytes over a pipe.
                            err_path = None
                            try:
                                err_name = f"worker_err_{os.getpid()}_{int(time.time())}_{random.randrange(10**6)}.err"
                                err_path = os.path.join(worker_tmp, err_name)
                                with open(out_path, 'wb') as ofh, open(err_path, 'wb') as efh:
                                    # Ensure child process runs in the per-worker tmp directory so
                                    # any auxiliary files created by external engines (e.g. Hudson
                                    # 'seedms' from ms) are placed inside worker_tmp and removed
                                    # during cleanup. This avoids leaving seed files in the repo
                                    # when engine == 'ms'.
                                    try:
                                        proc = subprocess.Popen(toks, stdout=ofh, stderr=efh, env=env, preexec_fn=pre, cwd=worker_tmp)
                                    except Exception:
                                        # Fallback: if worker_tmp is invalid for cwd, run without cwd
                                        proc = subprocess.Popen(toks, stdout=ofh, stderr=efh, env=env, preexec_fn=pre)
                                    _register_child_pid(proc.pid)
                                    # Optional external-worker peak sampler: monitor child pid RSS
                                    ext_peak_thread = None
                                    ext_peak_stop = None
                                    try:
                                        peak_dir_ext = os.environ.get('SIMULATOR_PEAK_DIR') if os.environ.get('SIMULATOR_PEAK_DIR') else (env.get('SIMULATOR_PEAK_DIR') if env and isinstance(env, dict) else None)
                                    except Exception:
                                        peak_dir_ext = None
                                    try:
                                        peak_interval_ext = float(os.environ.get('SIMULATOR_PEAK_INTERVAL', 0.02))
                                    except Exception:
                                        peak_interval_ext = 0.02
                                    if peak_dir_ext:
                                        try:
                                            os.makedirs(peak_dir_ext, exist_ok=True)
                                        except Exception:
                                            pass
                                        def _sampler_child_pid(child_pid, peak_dir_local, stop_event, interval):
                                            try:
                                                import psutil
                                            except Exception:
                                                return
                                            try:
                                                p = psutil.Process(child_pid)
                                            except Exception:
                                                return
                                            peak = 0
                                            try:
                                                while not stop_event.is_set():
                                                    try:
                                                        rss = int(getattr(p.memory_info(), 'rss', 0) or 0)
                                                        if rss > peak:
                                                            peak = rss
                                                    except Exception:
                                                        break
                                                    stop_event.wait(interval)
                                                # final sample
                                                try:
                                                    rss = int(getattr(p.memory_info(), 'rss', 0) or 0)
                                                    if rss > peak:
                                                        peak = rss
                                                except Exception:
                                                    pass
                                            finally:
                                                try:
                                                    path = os.path.join(peak_dir_local, f"peak_{child_pid}.txt")
                                                    with open(path, 'w') as pf:
                                                        pf.write(str(peak))
                                                except Exception:
                                                    pass
                                        try:
                                            ext_peak_stop = threading.Event()
                                            ext_peak_thread = threading.Thread(target=_sampler_child_pid, args=(proc.pid, peak_dir_ext, ext_peak_stop, peak_interval_ext), daemon=True)
                                            ext_peak_thread.start()
                                        except Exception:
                                            ext_peak_thread = None
                                            ext_peak_stop = None
                                    try:
                                        rc_proc = proc.wait()
                                    finally:
                                        try: _unregister_child_pid(proc.pid)
                                        except Exception: pass
                                        # signal sampler to stop and join
                                        try:
                                            if ext_peak_stop is not None:
                                                ext_peak_stop.set()
                                            if ext_peak_thread is not None:
                                                ext_peak_thread.join(timeout=1.0)
                                        except Exception:
                                            pass
                                # record output file path
                                out_acc.append(out_path)
                                # record stderr only if the file has content (avoid treating empty stderr as an error)
                                try:
                                    if err_path and os.path.isfile(err_path) and os.path.getsize(err_path) > 0:
                                        err_acc.append(f'FILE:{err_path}')
                                except Exception:
                                    # If anything goes wrong checking stderr, fall back to recording the marker
                                    try:
                                        if err_path:
                                            err_acc.append(f'FILE:{err_path}')
                                    except Exception:
                                        pass
                                # If the process exit code indicates failure, add an explicit error marker
                                if rc_proc and rc_proc != 0:
                                    err_acc.append(f'# ERROR: engine process failed with code {rc_proc}\n')
                            except Exception as e:
                                err_acc.append(f'# ERROR: running engine process to file failed: {e}\n')
                        except Exception as e:
                            err_acc.append(f'# ERROR: running engine process to file failed: {e}\n')
                    except Exception as e:
                        err_acc.append(f'# ERROR: failed to run engine process: {e}\n')
                # Return: rc, stdout-accumulated, stderr-accumulated, error-flag (1 if errors), replicate_index
                # Attempt to cleanup per-chunk tmpdir if it was created
                try:
                    wk = None
                    try:
                        wk = env.get('SIMULATOR_WORKER_TMP')
                    except Exception:
                        wk = None
                    # Determine whether to preserve the worker tmpdir.
                    preserve_local = False
                    try:
                        if env.get('SIMULATOR_PRESERVE_WORKER_TMP'):
                            preserve_local = True
                    except Exception:
                        pass
                    try:
                        if os.environ.get('SIMULATOR_PRESERVE_TMP'):
                            preserve_local = True
                    except Exception:
                        pass
                    # If the caller requested primary output (args._user_out_provided),
                    # keep worker tmpdirs until parent streams them; otherwise remove.
                    try:
                        user_wants_out = bool(getattr(args, '_user_out_provided', False))
                    except Exception:
                        user_wants_out = False
                    if wk and os.path.isdir(wk) and not preserve_local and not user_wants_out:
                        try:
                            shutil.rmtree(wk, ignore_errors=True)
                            if os.environ.get('SIMULATOR_DEBUG_TMP'):
                                try:
                                    sys.stderr.write(f"DEBUG: removed worker tmpdir at chunk end: {wk}\n")
                                    sys.stderr.flush()
                                except Exception:
                                    pass
                        except Exception:
                            pass
                except Exception:
                    pass
                return (rc, out_acc, err_acc, (1 if err_acc else None), replicate_index)
        results_loc = []
        bar_loc = None
        progress_done_loc = 0
        if args.progress and show_progress:
            if tqdm is not None:
                bar_loc = tqdm(total=total_reps_loc, desc='replicates', unit='rep')
                try:
                    setattr(args, '_progress_bar', bar_loc)
                except Exception:
                    pass
            else:
                if total_reps_loc == 1:
                    sys.stderr.write('# progress: 0/1\n')
                else:
                    sys.stderr.write('# WARNING: --progress requested but tqdm not installed; using textual updates.\n')

        # local helper to abort cleanly on RAM exceed: close progress bar then exit
        def _abort_due_to_ram(msg=None):
            try:
                if bar_loc is not None:
                    try:
                        bar_loc.close()
                    except Exception:
                        pass
            except Exception:
                pass
            # If caller provided an explicit message, use it; otherwise raise a
            # structured ExceedsMaxRam exception containing actionable suggestions.
            # Raise a module-level ExceedsMaxRam so callers can catch it.
            suggestion_parts = []
            try:
                suggestion_parts.append('Reduce --parallel')
                suggestion_parts.append('lower --sims-per-work (for heavy per-worker memory)')
                suggestion_parts.append('or lower --replicates or --sample-individuals')
            except Exception:
                pass
            suggestion = '; '.join(p for p in suggestion_parts if p)
            try:
                arg_msg = getattr(args, '_ram_exceeded_message', None)
            except Exception:
                arg_msg = None
            final_msg = msg or arg_msg or (f"Aborting run due to exceeding max RAM. Suggestion: {suggestion}")
            raise ExceedsMaxRam(final_msg)

        # helper to update the progress bar only when RAM hasn't been exceeded
        def _maybe_update_bar(delta=1):
            """Safely update the progress UI.

            If a tqdm bar is active, call its update method. Otherwise, emit
            a textual progress line to stderr. Respect a RAM-exceeded flag
            by no-op'ing when set.
            """
            try:
                if getattr(args, '_ram_exceeded', False):
                    return
                if bar_loc is not None:
                    try:
                        bar_loc.update(delta)
                    except Exception:
                        # Best-effort: ignore progress update failures
                        pass
            except Exception:
                pass
        # Special-case msprime: run each chunk in its own process (ProcessPool)
        # so msprime uses physical cores rather than Python threads.
        if engine == 'msprime':
            payloads = []
            envs = []
            try:
                import psutil
                phys_total = psutil.cpu_count(logical=False) or (os.cpu_count() or 1)
            except Exception:
                phys_total = os.cpu_count() or 1
            req_workers = max(1, int(getattr(args, 'parallel', 1)))
            workers = max(1, min(req_workers, (chunk_process_count or 1)))
            base = max(1, int(phys_total) // workers)
            rem = int(phys_total) % workers
            for idx, chunk_size in enumerate(chunks_loc):
                # seed for this chunk if available
                seed_val = None
                if seeds_for_chunks and idx < len(seeds_for_chunks):
                    seed_val = seeds_for_chunks[idx]
                # build msprime payload (JSON-serializable)
                payload = {
                    'species': meta_loc.get('species_id') or args.species,
                    'model_id': meta_loc.get('model_id') or args.model,
                    'chromosome': meta_loc.get('chromosome') or getattr(args, 'chromosome', None),
                    'length': int(meta_loc.get('length') or length_val or 0),
                    'pop_order': meta_loc.get('pop_order') or [],
                    'counts': meta_loc.get('counts_disc') or [],
                    'ploidy': meta_loc.get('ploidy') or 1,
                    'rep_count': int(chunk_size),
                    'seed': seed_val,
                }
                try:
                    payload['suppress_worker_sfs'] = (not bool(getattr(args, 'sfs', False)))
                except Exception:
                    payload['suppress_worker_sfs'] = False
                # Workers should provide SFS_AF markers; do not suppress SFS emission.
                payload['suppress_worker_sfs'] = False
                # If the top-level invocation requested only SFS (no primary output), ask workers
                # to avoid producing the full ms/vcf text to save memory/time.
                try:
                    only_sfs_top = bool(getattr(args, '_suppress_primary_output', False) or (getattr(args, 'sfs', False) and not getattr(args, '_user_out_provided', False)))
                except Exception:
                    only_sfs_top = False
                payload['only_sfs'] = bool(only_sfs_top)
                # Preserve per-worker tmpdirs automatically when the user requested
                # a primary output (user provided --out) and this is not a
                # SFS-only invocation; workers' outputs are needed for combining.
                # Never preserve per-worker tmpdirs; always clean up after workers.
                payload['preserve_worker_tmp'] = False
                per_worker_threads = base + (1 if (idx % workers) < rem else 0)
                env_worker = os.environ.copy()
                for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
                    env_worker[k] = str(per_worker_threads)
                # allocate a per-worker tmpdir beneath the chosen root so workers don't share one global dir
                try:
                    root_tmp = _ensure_local_tmpdir()
                    # child workers may run on other hosts/processes; use deterministic subdir name per idx
                    worker_name = f"worker_{idx}_{os.getpid()}_{int(time.time())}_{random.randrange(10**6)}"
                    worker_tmp = os.path.join(root_tmp, worker_name) if root_tmp else None
                    if worker_tmp:
                        try:
                            os.makedirs(worker_tmp, exist_ok=True)
                            env_worker['SIMULATOR_WORKER_TMP'] = worker_tmp
                            # ask external subprocesses to use this tmpdir too
                            env_worker['TMPDIR'] = worker_tmp
                            env_worker['TMP'] = worker_tmp
                            env_worker['TEMP'] = worker_tmp
                        except Exception:
                            pass
                except Exception:
                    pass
                payloads.append(payload)
                envs.append(env_worker)

            # Run payloads in separate processes (capped by 'workers' so
            # concurrency honors --parallel and physical-core sharing)
            futs = []
            # Decide how many concurrent msprime worker processes to spawn.
            # Decide max concurrent msprime worker processes. Cap to the
            # requested parallelism, the number of payloads, and the number
            # of physical CPU cores to avoid oversubscription.
            max_proc = min(len(payloads), req_workers, phys_total)

            if getattr(args, '_ram_exceeded', False):
                _abort_due_to_ram()
            with concurrent.futures.ProcessPoolExecutor(max_workers=max_proc) as pex:
                for i, (pld, envw) in enumerate(zip(payloads, envs)):
                    futs.append((i, pex.submit(_msprime_worker_run, pld, envw)))
                for idx, fut in futs:
                    try:
                        if getattr(args, '_ram_exceeded', False):
                            _abort_due_to_ram()
                        out_txt, err_txt, rc = fut.result()
                    except Exception as e:
                        out_txt, err_txt, rc = '', f'# ERROR: msprime worker failed in process: {e}\n', 1
                    results_loc.append((int(payloads[idx].get('rep_count',1)), out_txt, err_txt, (1 if rc else None), idx))
                    # Update progress
                    if bar_loc is not None:
                        _maybe_update_bar(payloads[idx].get('rep_count', 1))
                    elif args.progress and show_progress:
                        progress_done_loc += payloads[idx].get('rep_count', 1)
                        sys.stderr.write(f'# progress: {progress_done_loc}/{total_reps_loc}\n')

                    # Attempt immediate cleanup of per-worker tmpdir created for this worker
                    try:
                        envw_local = envs[idx] if idx < len(envs) else None
                        wk = None
                        if envw_local:
                            wk = envw_local.get('SIMULATOR_WORKER_TMP') or envw_local.get('TMPDIR') or None
                        # Determine preservation: payload-level or env-level or global preserve
                        preserve = False
                        try:
                            if payloads[idx].get('preserve_worker_tmp', False):
                                preserve = True
                        except Exception:
                            pass
                        try:
                            if envw_local and envw_local.get('SIMULATOR_PRESERVE_WORKER_TMP'):
                                preserve = True
                        except Exception:
                            pass
                        if os.environ.get('SIMULATOR_PRESERVE_TMP'):
                            preserve = True
                        if wk and os.path.isdir(wk) and not preserve:
                            try:
                                shutil.rmtree(wk, ignore_errors=True)
                                if os.environ.get('SIMULATOR_DEBUG_TMP'):
                                    try:
                                        sys.stderr.write(f"DEBUG: removed worker tmpdir after worker completion: {wk}\n")
                                        sys.stderr.flush()
                                    except Exception:
                                        pass
                            except Exception:
                                pass
                    except Exception:
                        pass
            results_loc.sort(key=lambda x: x[4])
        else:
            if per_rep_mode_loc:
                # Abort early if RAM monitor signalled an exceedance
                if getattr(args, '_ram_exceeded', False):
                    _abort_due_to_ram()
                # Limit concurrent chunk runners to threads_loc (honours --parallel).
                # Submit all chunk tasks but allow executor to queue excess tasks so
                # only up to `threads_loc` run at once and remaining chunks wait.
                with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(len(chunks_loc), threads_loc))) as ex:
                    futs = []
                    rep_counter = 0
                    for idx, c in enumerate(chunks_loc):
                        futs.append(ex.submit(run_chunk_loc, c, rep_counter, idx))
                        rep_counter += 1
                    for fut in concurrent.futures.as_completed(futs):
                        if getattr(args, '_ram_exceeded', False):
                            # stop waiting for remaining futures and abort
                            _abort_due_to_ram()
                        rc, out, err, code, rep_idx = fut.result()
                        results_loc.append((rc, out, err, code, rep_idx))
                        # Attempt immediate cleanup of per-chunk tmpdir if set in env during the chunk
                        try:
                            # run_chunk_loc sets SIMULATOR_WORKER_TMP in its local env var if used for subprocesses
                            # try to get from current process env as a best-effort (may not be accurate)
                            wk = None
                            try:
                                wk = os.environ.get('SIMULATOR_WORKER_TMP')
                            except Exception:
                                wk = None
                            # If the run_chunk returned a rep_idx, payloads may contain preserve flags
                            preserve = False
                            try:
                                # if results list contains env info embed later, skip; use global preserve if set
                                if os.environ.get('SIMULATOR_PRESERVE_TMP'):
                                    preserve = True
                            except Exception:
                                pass
                            if wk and os.path.isdir(wk) and not preserve:
                                try:
                                    shutil.rmtree(wk, ignore_errors=True)
                                    if os.environ.get('SIMULATOR_DEBUG_TMP'):
                                        try:
                                            sys.stderr.write(f"DEBUG: removed chunk tmpdir after completion: {wk}\n")
                                            sys.stderr.flush()
                                        except Exception:
                                            pass
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        if bar_loc is not None:
                            # rc is the number of replicates executed in this chunk
                            try:
                                _maybe_update_bar(int(rc))
                            except Exception:
                                _maybe_update_bar(1)
                        elif args.progress and show_progress:
                            try:
                                progress_done_loc += int(rc)
                            except Exception:
                                progress_done_loc += 1
                            sys.stderr.write(f'# progress: {progress_done_loc}/{total_reps_loc}\n')
                results_loc.sort(key=lambda x: (x[4] if len(x) > 4 else 0))
            else:
                if len(chunks_loc) == 1:
                    rc, out, err, code, rep_idx = run_chunk_loc(chunks_loc[0], 0, 0)
                    results_loc.append((rc, out, err, code, rep_idx))
                    if bar_loc is not None:
                        _maybe_update_bar(total_reps_loc)
                    elif args.progress and show_progress and total_reps_loc == 1:
                        sys.stderr.write('# progress: 1/1\n')
                else:
                    # Abort early if RAM monitor signalled an exceedance
                    if getattr(args, '_ram_exceeded', False):
                        _abort_due_to_ram()
                    # When chunking into multiple chunks we still want to cap
                    # concurrent runners to `threads_loc` (derived from --parallel)
                    # so remaining chunks are queued and executed as workers free up.
                    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(len(chunks_loc), threads_loc))) as ex:
                        futs = [ex.submit(run_chunk_loc, c, i, i) for i, c in enumerate(chunks_loc)]
                        for fut in concurrent.futures.as_completed(futs):
                            if getattr(args, '_ram_exceeded', False):
                                _abort_due_to_ram()
                            rc, out, err, code, rep_idx = fut.result()
                            results_loc.append((rc, out, err, code, rep_idx))
                            if bar_loc is not None:
                                try:
                                    _maybe_update_bar(int(rc))
                                except Exception:
                                    _maybe_update_bar(1)
                            elif args.progress and show_progress:
                                try:
                                    progress_done_loc += int(rc)
                                except Exception:
                                    progress_done_loc += 1
                                sys.stderr.write(f'# progress: {progress_done_loc}/{total_reps_loc}\n')
        if bar_loc is not None:
            bar_loc.close()
        def _expand_marker_text(s):
            try:
                if not s:
                    return ''
                # marker form: FILE:/abs/path
                if isinstance(s, str) and s.startswith('FILE:'):
                    p = s[5:]
                    try:
                        if os.path.isfile(p):
                            with open(p, 'r', errors='replace') as fh:
                                return fh.read()
                    except Exception:
                        return s
                # If it's a path to a file, read it
                if isinstance(s, str) and os.path.isfile(s):
                    try:
                        with open(s, 'r', errors='replace') as fh:
                            return fh.read()
                    except Exception:
                        return s
                return str(s)
            except Exception:
                try:
                    return str(s)
                except Exception:
                    return ''

        for tpl in results_loc:
            if len(tpl) == 5:
                rc, out, err, code, _r = tpl
            else:
                rc, out, err, code = tpl
            if code is not None:
                sys.stderr.write(f'# ERROR: {engine} exited with code {code}.\n')
                if out:
                    sys.stderr.write(_expand_marker_text(out))
                if err:
                    sys.stderr.write(_expand_marker_text(err))
                sys.exit(code)
        # If the caller only requested SFS (no primary output) and this
        # is not the msprime in-process branch, avoid concatenating the
        # potentially large stdout blobs into memory. Instead stream them
        # to a temp file and compute the mean SFS into a compact textual
        # representation that we return in place of the full ms/vcf text.
        try:
            want_only_sfs = bool(getattr(args, 'sfs', False) and not getattr(args, '_user_out_provided', False))
        except Exception:
            want_only_sfs = False

        # Whether user requested deletion of per-chunk primary outputs after
        # extracting SFS. This lets users keep disk/RAM low by removing heavy
        # ms/vcf files as soon as their SFS contribution is captured.
        try:
            # Discard per-chunk primary outputs by default only when --sfs is requested
            # and the user did NOT provide an explicit --output. If the user provided
            # --output we should keep the primary ms/vcf files and not replace
            # primary output with the compact SFS text.
            discard_requested = bool(getattr(args, 'sfs', False) and not getattr(args, '_user_out_provided', False))
        except Exception:
            discard_requested = False
        if discard_requested:
            try:
                sys.stderr.write('# INFO: per-chunk primary outputs will be removed after SFS extraction (SFS-only mode default). Use --output to keep primary outputs.\n')
            except Exception:
                pass

        # Helper: read stdout content from a results_loc entry which may be
        # either a string or a filepath produced by workers.
        def _read_stdout_entry(e):
            try:
                stdout_val = e[1] if len(e) > 1 else None
                # If worker returned a list of file paths, read each file and
                # concatenate. If it's a single filepath string, read it.
                if isinstance(stdout_val, (list, tuple)):
                    parts = []
                    for item in stdout_val:
                        try:
                            if isinstance(item, str) and os.path.isfile(item):
                                with open(item, 'r') as rfh:
                                    parts.append(rfh.read())
                            else:
                                parts.append(str(item))
                        except Exception:
                            # ignore individual read errors
                            pass
                    return ''.join(parts)
                if isinstance(stdout_val, str) and os.path.isfile(stdout_val):
                    try:
                        with open(stdout_val, 'r') as rfh:
                            return rfh.read()
                    except Exception:
                        return ''
                return stdout_val or ''
            except Exception:
                return ''

        def _expand_err_entry(e):
            try:
                v = e[2] if len(e) > 2 else None
                if not v:
                    return ''
                if isinstance(v, str) and v.startswith('FILE:'):
                    p = v[5:]
                    try:
                        if os.path.isfile(p):
                            with open(p, 'r', errors='replace') as fh:
                                return fh.read()
                    except Exception:
                        return v
                if isinstance(v, str) and os.path.isfile(v):
                    try:
                        with open(v, 'r', errors='replace') as fh:
                            return fh.read()
                    except Exception:
                        return v
                return str(v)
            except Exception:
                return ''

        def _read_stdout_and_maybe_delete(e):
            # Similar to _read_stdout_entry but removes file-backed outputs
            # when the user requested discard behavior.
            try:
                stdout_val = e[1] if len(e) > 1 else None
                if isinstance(stdout_val, (list, tuple)):
                    parts = []
                    for item in stdout_val:
                        try:
                            if isinstance(item, str) and os.path.isfile(item):
                                try:
                                    with open(item, 'r') as rfh:
                                        parts.append(rfh.read())
                                except Exception:
                                    pass
                                if discard_requested:
                                    try:
                                        os.unlink(item)
                                    except Exception:
                                        pass
                            else:
                                parts.append(str(item))
                        except Exception:
                            pass
                    return ''.join(parts)
                if isinstance(stdout_val, str) and os.path.isfile(stdout_val):
                    try:
                        with open(stdout_val, 'r') as rfh:
                            data = rfh.read()
                    except Exception:
                        data = ''
                    if discard_requested:
                        try:
                            os.unlink(stdout_val)
                        except Exception:
                            pass
                    return data
                return stdout_val or ''
            except Exception:
                return ''

        # SFS-only fast path for non-msprime engines: stream worker outputs
        # (including file-backed outputs) into a single temporary file and
        # compute the mean SFS from that file.
        concatenated_stdout_loc = ''
        concatenated_stderr_loc = ''
        if want_only_sfs and engine != 'msprime':
            import tempfile as _tempfile, os as _os
            tmp_path = None
            try:
                tmpdir = _ensure_local_tmpdir()
                discard_files = bool(want_only_sfs)
                with _tempfile.NamedTemporaryFile('w', delete=False, dir=tmpdir) as tf:
                    for t in results_loc:
                        try:
                            # stream file-backed outputs or in-memory strings
                            stdout_val = t[1] if len(t) > 1 else None
                            if isinstance(stdout_val, (list, tuple)):
                                for item in stdout_val:
                                    try:
                                        if isinstance(item, str) and os.path.isfile(item):
                                            with open(item, 'r') as rfh:
                                                chunk = rfh.read()
                                            if chunk:
                                                tf.write(chunk)
                                            if discard_files:
                                                try:
                                                    os.unlink(item)
                                                except Exception:
                                                    pass
                                        else:
                                            s = str(item)
                                            if s:
                                                tf.write(s)
                                    except Exception:
                                        # ignore individual read/remove errors
                                        pass
                            elif isinstance(stdout_val, str) and os.path.isfile(stdout_val):
                                try:
                                    with open(stdout_val, 'r') as rfh:
                                        chunk = rfh.read()
                                    if chunk:
                                        tf.write(chunk)
                                    if discard_files:
                                        try:
                                            os.unlink(stdout_val)
                                        except Exception:
                                            pass
                                except Exception:
                                    pass
                            else:
                                stdout_chunk = _read_stdout_entry(t)
                                if stdout_chunk:
                                    tf.write(stdout_chunk)
                        except Exception:
                            pass
                    tmp_path = tf.name
                # Compute SFS from the streamed file and return a compact SFS text
                # Compute SFS from the streamed file and return a compact SFS text
                try:
                    header_sfs, lines_sfs = _STREAM_SFS_PLACEHOLDER(...)
                    _sfs_text_tmp = header_sfs + '\n' + '\n'.join(lines_sfs) + '\n'
                    if want_only_sfs:
                        concatenated_stdout_loc = _sfs_text_tmp
                    else:
                        try:
                            with open(tmp_path, 'r', errors='replace') as _rf:
                                concatenated_stdout_loc = _rf.read()
                        except Exception:
                            concatenated_stdout_loc = ''
                except Exception:
                    try:
                        with open(tmp_path, 'r', errors='replace') as _rf:
                            concatenated_stdout_loc = _rf.read()
                    except Exception:
                        concatenated_stdout_loc = ''
                        v = e[2] if len(e) > 2 else None
                        if not v:
                            return ''
                        if isinstance(v, str) and v.startswith('FILE:'):
                            p = v[5:]
                            try:
                                if os.path.isfile(p):
                                    with open(p, 'r', errors='replace') as fh:
                                        return fh.read()
                            except Exception:
                                return v
                        if isinstance(v, str) and os.path.isfile(v):
                            try:
                                with open(v, 'r', errors='replace') as fh:
                                    return fh.read()
                            except Exception:
                                return v
                        return str(v)
                    except Exception:
                        return ''

                concatenated_stderr_loc = ''.join((_expand_err_entry(t) for t in results_loc))
            finally:
                try:
                    if tmp_path and _os.path.exists(tmp_path):
                        _os.unlink(tmp_path)
                except Exception:
                    pass
        else:
            # Default concatenation path: either stream-and-delete to save RAM
            # (when discard_requested and we don't need the full ms output for
            # target-snps enforcement), or concatenate into memory as before.
            if want_only_sfs and discard_requested and getattr(args, 'min_snps', None) is None:
                # Stream worker outputs to a temporary file, compute SFS from it,
                # and avoid retaining full ms/vcf text in memory.
                import tempfile as _tempfile, os as _os
                tmp_path = None
                try:
                    tmpdir = _ensure_local_tmpdir()
                    with _tempfile.NamedTemporaryFile('w', delete=False, dir=tmpdir) as tf:
                        for t in results_loc:
                            try:
                                stdout_val = t[1] if len(t) > 1 else None
                                if isinstance(stdout_val, (list, tuple)):
                                    for item in stdout_val:
                                        try:
                                            if isinstance(item, str) and os.path.isfile(item):
                                                with open(item, 'r') as rfh:
                                                    chunk = rfh.read()
                                                if chunk:
                                                    tf.write(chunk)
                                                try:
                                                    os.unlink(item)
                                                except Exception:
                                                    pass
                                            else:
                                                s = str(item)
                                                if s:
                                                    tf.write(s)
                                        except Exception:
                                            pass
                                elif isinstance(stdout_val, str) and os.path.isfile(stdout_val):
                                    try:
                                        with open(stdout_val, 'r') as rfh:
                                            chunk = rfh.read()
                                        if chunk:
                                            tf.write(chunk)
                                        try:
                                            os.unlink(stdout_val)
                                        except Exception:
                                            pass
                                    except Exception:
                                        pass
                                else:
                                    stdout_chunk = _read_stdout_entry(t)
                                    if stdout_chunk:
                                        tf.write(stdout_chunk)
                            except Exception:
                                pass
                        tmp_path = tf.name
                    # Compute SFS from the streamed file and return compact SFS text
                    # Compute SFS from the streamed file and return compact SFS text
                    try:
                        header_sfs, lines_sfs = _STREAM_SFS_PLACEHOLDER(...)
                        _sfs_text_tmp = header_sfs + '\n' + '\n'.join(lines_sfs) + '\n'
                        if want_only_sfs:
                            concatenated_stdout_loc = _sfs_text_tmp
                        else:
                            try:
                                with open(tmp_path, 'r', errors='replace') as _rf:
                                    concatenated_stdout_loc = _rf.read()
                            except Exception:
                                concatenated_stdout_loc = ''
                    except Exception:
                        try:
                            with open(tmp_path, 'r', errors='replace') as _rf:
                                concatenated_stdout_loc = _rf.read()
                        except Exception:
                            concatenated_stdout_loc = ''

                    concatenated_stderr_loc = ''.join((_expand_err_entry(t) for t in results_loc))
                finally:
                    try:
                        if tmp_path and _os.path.exists(tmp_path):
                            _os.unlink(tmp_path)
                    except Exception:
                        pass
            else:
                # Default behavior: concatenate into memory
                parts_out = []
                for t in results_loc:
                    if discard_requested:
                        parts_out.append(_read_stdout_and_maybe_delete(t))
                    else:
                        parts_out.append(_read_stdout_entry(t))
                concatenated_stdout_loc = ''.join(parts_out)
            def _expand_err_entry(e):
                try:
                    v = e[2] if len(e) > 2 else None
                    if not v:
                        return ''
                    if isinstance(v, str) and v.startswith('FILE:'):
                        p = v[5:]
                        try:
                            if os.path.isfile(p):
                                with open(p, 'r', errors='replace') as fh:
                                    return fh.read()
                        except Exception:
                            return v
                    if isinstance(v, str) and os.path.isfile(v):
                        try:
                            with open(v, 'r', errors='replace') as fh:
                                return fh.read()
                        except Exception:
                            return v
                    return str(v)
                except Exception:
                    return ''

            concatenated_stderr_loc = ''.join((_expand_err_entry(t) for t in results_loc))
            # Cleanup per-worker output files now that we've read them
            try:
                for t in results_loc:
                    try:
                        stdout_val = t[1] if len(t) > 1 else None
                        if isinstance(stdout_val, str) and os.path.isfile(stdout_val):
                            try:
                                os.unlink(stdout_val)
                            except Exception:
                                pass
                    except Exception:
                        pass
            except Exception:
                pass
        # If msprime workers wrote VCF directly (start with VCF header), mark produced_format
        try:
            if engine == 'msprime' and isinstance(meta_loc, dict):
                sniff_blob = (concatenated_stdout_loc or '') + '\n' + (concatenated_stderr_loc or '')
                if '#CHROM' in sniff_blob or '##fileformat' in sniff_blob:
                    meta_loc['produced_format'] = 'vcf'
                else:
                    # default to ms-like if not detected
                    meta_loc.setdefault('produced_format', 'ms')
        except Exception:
            pass
        # Persist chunk seeding metadata for potential neutral pairing reuse
        try:
            meta_loc.setdefault('runtime', {})
            meta_loc['runtime']['chunk_sizes'] = chunks_loc
            # Persist chunk seeds for engines that support explicit seeding. Include ms and msprime.
            if engine in ('discoal', 'msms', 'msprime', 'ms'):
                meta_loc['runtime']['chunk_seeds'] = seeds_for_chunks
            meta_loc['runtime']['per_rep_mode'] = per_rep_mode_loc
        except Exception:
            pass
        return meta_loc, raw_cmd_line_loc, concatenated_stdout_loc, concatenated_stderr_loc

    # Initial simulation at refined/derived length
    meta, raw_cmd_line, concatenated_stdout, concatenated_stderr = _simulate_current_length(args.length, show_progress=True)
    # Debug-only dump removed per cleanup

    # Enforcement: ensure ALL replicates meet/exceed target and optionally shrink overshoot
    if args.min_snps is not None:
        target_snps = args.min_snps
        tol_internal = getattr(args,'min_snps_tol', 0.0)
        max_enforce_iters = 5
        enforce_iter = 0  # count of inflation iterations
        length_current = args.length
        pattern_segs = re.compile(r'^segsites:\s*(\d+)', re.MULTILINE)
        while enforce_iter < max_enforce_iters:
            seg_list = [int(x) for x in pattern_segs.findall(concatenated_stdout)]
            if not seg_list:
                break
            min_segs = min(seg_list)
            if min_segs < target_snps:
                factor = target_snps / max(1, min_segs)
                inflation = 1.05 if engine in ('ms','scrm') else (1.15 if engine == 'msms' else 1.20)
                new_length = int(math.ceil(length_current * factor * inflation))
                if new_length <= length_current:
                    new_length = length_current + 1
                # no debug printing
                length_current = new_length
                meta, raw_cmd_line, concatenated_stdout, concatenated_stderr = _simulate_current_length(length_current, show_progress=False)
                enforce_iter += 1
                continue
            break
        args.length = length_current
    # After enforcement loop, optionally perform overshoot shrink and record stats
    if args.min_snps is not None:
        seg_list_final = [int(x) for x in pattern_segs.findall(concatenated_stdout)]
        # no debug printing
        overshoot_status = None
        shrink_bracket_attempts = 0
        shrink_extra_reductions = 0
        shrink_bs_iters = 0
        if tol_internal > 0:
            seg_list_curr = [int(x) for x in pattern_segs.findall(concatenated_stdout)]
            if seg_list_curr:
                hi_len = args.length
                hi_stdout = concatenated_stdout
                hi_meta = meta
                # Track the corresponding command line for the current 'hi' candidate.
                # raw_cmd_line at this point reflects the last simulation before shrink attempts.
                hi_raw_cmd_line = raw_cmd_line
                hi_seg_list = seg_list_curr
                lo_len = None
                lo_seg_list = None
                attempts = 0
                while attempts < 6:
                    max_segs_curr = max(hi_seg_list)
                    target_cap = target_snps * (1.0 + tol_internal)
                    scale_cap = target_cap / float(max_segs_curr)
                    trial_len = int(math.ceil(hi_len * max(0.60, min(0.95, scale_cap))))
                    if trial_len >= hi_len or trial_len < 100:
                        break
                    meta_trial, raw_cmd_line_trial, stdout_trial, stderr_trial = _simulate_current_length(trial_len, show_progress=False)
                    seg_list_trial = [int(x) for x in pattern_segs.findall(stdout_trial)]
                    if not seg_list_trial:
                        break
                    if min(seg_list_trial) >= target_snps:
                        hi_len, hi_stdout, hi_meta, hi_seg_list = trial_len, stdout_trial, meta_trial, seg_list_trial
                        hi_raw_cmd_line = raw_cmd_line_trial
                        attempts += 1
                        shrink_bracket_attempts = attempts
                        continue
                    else:
                        lo_len = trial_len
                        lo_seg_list = seg_list_trial
                        break
                if lo_len is None:
                    probe_len = int(hi_len * 0.85)
                    extra = 0
                    while extra < 6 and probe_len >= 100:
                        meta_trial, raw_cmd_line_trial, stdout_trial, stderr_trial = _simulate_current_length(probe_len, show_progress=False)
                        seg_list_trial = [int(x) for x in pattern_segs.findall(stdout_trial)]
                        if not seg_list_trial:
                            break
                        if min(seg_list_trial) >= target_snps:
                            hi_len, hi_stdout, hi_meta, hi_seg_list = probe_len, stdout_trial, meta_trial, seg_list_trial
                            hi_raw_cmd_line = raw_cmd_line_trial
                            probe_len = int(probe_len * 0.85)
                            shrink_extra_reductions = extra + 1
                        else:
                            lo_len = probe_len
                            lo_seg_list = seg_list_trial
                            break
                        extra += 1
                if lo_len is not None:
                    bs_iter = 0
                    while bs_iter < 12 and lo_len + 1 < hi_len:
                        mid = (lo_len + hi_len) // 2
                        if mid == hi_len:
                            break
                        meta_mid, raw_cmd_line_mid, stdout_mid, stderr_mid = _simulate_current_length(mid, show_progress=False)
                        seg_list_mid = [int(x) for x in pattern_segs.findall(stdout_mid)]
                        if not seg_list_mid:
                            break
                        if min(seg_list_mid) >= target_snps:
                            hi_len, hi_stdout, hi_meta, hi_seg_list = mid, stdout_mid, meta_mid, seg_list_mid
                            hi_raw_cmd_line = raw_cmd_line_mid
                        else:
                            lo_len = mid
                        bs_iter += 1
                        shrink_bs_iters = bs_iter
                args.length = hi_len
                concatenated_stdout = hi_stdout
                meta = hi_meta
                # Update the raw command line to remain consistent with the final (possibly shrunken) length.
                raw_cmd_line = hi_raw_cmd_line
                max_final = max(hi_seg_list)
                if max_final <= target_snps * (1.0 + tol_internal):
                    overshoot_status = 'met'
                else:
                    overshoot_status = 'unmet_minimal_feasible'
                    # no debug printing
                # no debug printing
            else:
                overshoot_status = 'no_data'
        final_seg_list = [int(x) for x in pattern_segs.findall(concatenated_stdout)] if 'pattern_segs' in locals() else []
        if final_seg_list:
            meta.setdefault('enforcement_stats', {})
            meta['enforcement_stats'].update({
                'target_snps': target_snps,
                'tol': tol_internal,
                'segsites_min': min(final_seg_list),
                'segsites_max': max(final_seg_list),
                'overshoot': overshoot_status or ('n/a' if tol_internal == 0 else 'unknown'),
                'enforce_iters': enforce_iter,
                'shrink_bracket_attempts': shrink_bracket_attempts,
                'shrink_extra_reductions': shrink_extra_reductions,
                'shrink_bs_iters': shrink_bs_iters,
            })
    # Decide default naming pieces
    model_id = getattr(args, 'model', 'model')
    pop0 = None
    if engine == 'discoal':
        pop0 = getattr(args, 'discoal_pop0', None)
    elif engine == 'msms':
        pop0 = getattr(args, 'discoal_pop0', None)
    else:  # ms
        pop0 = getattr(args, 'discoal_pop0', None)
    if not pop0:
        # fallback to first deme in user_order
        pop0 = user_order[0] if user_order else 'pop0'
    chrom = getattr(args, 'chromosome', None)
    sweeppop = None
    if engine in ('discoal','msms'):
        sweeppop = getattr(args, 'discoal_pop0', None)
    # Base stem
    parts_name = [model_id, pop0]
    if chrom:
        parts_name.append(str(chrom))
    if sweeppop and sweeppop != pop0:
        parts_name.append(sweeppop)
    default_stem = '_'.join(parts_name)
    # If paired-neutral was requested but the primary run is neutral, ignore as redundant.
    try:
        if getattr(args, 'paired_neutral', False):
            if _primary_simulation_is_neutral(meta):
                eng0 = meta.get('engine')
                reason = 'engine is neutral-only' if eng0 in ('ms','scrm','msprime') else 'no sweep parameters provided'
                try:
                    sys.stderr.write(f"# INFO: primary simulation is neutral ({reason}); ignoring --paired-neutral.\n")
                except Exception:
                    pass
                setattr(args, 'paired_neutral', False)
    except Exception:
        pass

    if args.sfs:
        # track whether we've already written the SFS file to avoid duplicate writes/messages
        sfs_written = False
        # defer actual file creation until the very end of the run
        sfs_pending = None
        # If engine is msprime, prefer in-process allele_frequency_spectrum() to compute SFS
        if engine == 'msprime':
            try:
                sys.stderr.write('#DEBUG: entering msprime in-process SFS branch\n')
            except Exception:
                pass
            # Use a simple, serial in-process msprime SFS approach (like junk/1.sfs.py):
            # - no windows, no multithreading
            # - run reps sequentially, compute ts.allele_frequency_spectrum
            # - produce textual SFS (#SFS header and mean/per-rep lines)
            try:
                eng = sps.get_engine('msprime')
            except Exception as e:
                # fallback to text-based parsing if msprime engine unavailable
                if args.sfs_mode == 'mean':
                    header_sfs, lines_sfs = compute_sfs_stream_mean(concatenated_stdout, n_hap_expected=None, normalized=args.sfs_normalized)
                else:
                    header_sfs, lines_sfs = compute_sfs_fast(concatenated_stdout, normalized=args.sfs_normalized, mode=args.sfs_mode)
                sfs_text = header_sfs + '\n' + '\n'.join(lines_sfs) + '\n'
            else:
                # Prepare species/model/samples
                sp = None
                model_obj = None
                try:
                    sp = sps.get_species(getattr(args, 'species', None))
                except Exception:
                    sp = None
                if isinstance(meta, dict) and meta.get('model_obj') is not None:
                    model_obj = meta.get('model_obj')
                else:
                    try:
                        if sp is not None:
                            model_obj = sp.get_demographic_model(getattr(args, 'model', None))
                    except Exception:
                        model_obj = None

                # Determine ploidy and counts (counts_disc are hap counts)
                pl_use = 1
                if sp is not None:
                    try: pl_use = int(getattr(sp, 'ploidy', pl_use) or pl_use)
                    except Exception: pl_use = pl_use
                desired_order = None
                counts_disc = None
                if isinstance(meta, dict):
                    desired_order = meta.get('pop_order') or None
                    counts_disc = meta.get('counts_disc') or None
                if not desired_order:
                    desired_order = user_order
                if not counts_disc:
                    try:
                        counts_disc = [int(pl_use * int(x)) for x in individual_counts]
                    except Exception:
                        counts_disc = [0] * len(desired_order)

                # Build samples_per_pop as individuals
                samples_per_pop = {}
                for i, pop in enumerate(desired_order):
                    hapc = int(counts_disc[i]) if i < len(counts_disc) else 0
                    samples_per_pop[pop] = int(hapc // pl_use)

                # Choose contig if possible; prefer contig with model mutation_rate to avoid warnings
                chr_name = getattr(args, 'chromosome', None) or (meta.get('chromosome') if isinstance(meta, dict) else None)
                chosen_contig = None
                try:
                    model_rate = getattr(model_obj, 'mutation_rate', None) if model_obj is not None else None
                    if sp is not None:
                        try:
                            if model_rate is not None:
                                chosen_contig = sp.get_contig(chr_name, mutation_rate=float(model_rate), left=0, right=(int(args.length) if args.length else None))
                            else:
                                chosen_contig = sp.get_contig(chr_name)
                        except Exception:
                            try:
                                chosen_contig = sp.get_contig(chr_name)
                            except Exception:
                                chosen_contig = None
                except Exception:
                    chosen_contig = None

                # Ensure chosen_contig has mutation_rate matching model when available
                try:
                    if chosen_contig is not None and model_obj is not None and getattr(model_obj, 'mutation_rate', None) is not None:
                        try: chosen_contig.mutation_rate = float(getattr(model_obj, 'mutation_rate'))
                        except Exception: pass
                except Exception:
                    pass

                # Compute total haplotypes and header
                total_hap = sum(int(c) for c in counts_disc)
                n_hap = int(total_hap) if total_hap > 0 else 1
                header = '#SFS\t' + '\t'.join(str(k) for k in range(1, n_hap))

                # Determine replicate count
                try:
                    reps = int(meta.get('reps', 1) if isinstance(meta, dict) else int(getattr(args, 'reps', 1)))
                except Exception:
                    reps = int(getattr(args, 'reps', 1))

                per_rep_lines = []
                accum = [0.0] * (n_hap - 1)
                used = 0

                # Debug: print msprime SFS context
                try:
                    sys.stderr.write(f"#DEBUG msprime SFS context: counts_disc={counts_disc} total_hap={total_hap} n_hap={n_hap} mu_local={mu_local}\n")
                except Exception:
                    pass

                # Check for per-rep SFS stored from the primary simulate loop
                stored_sfs = None
                try:
                    stored_sfs = meta.get('runtime', {}).get('per_rep_sfs') if isinstance(meta, dict) else None
                except Exception:
                    stored_sfs = None

                # If stored per-rep SFS arrays are available, compute the requested output directly
                if stored_sfs:
                    try:
                        # build header
                        total_hap = sum(int(c) for c in counts_disc)
                        n_hap = int(total_hap) if total_hap > 0 else 1
                        header = '#SFS\t' + '\t'.join(str(k) for k in range(1, n_hap))

                        if args.sfs_mode == 'per-rep':
                            per_rep_lines = []
                            for ridx, sfs_r in enumerate(stored_sfs):
                                if sfs_r is None:
                                    per_rep_lines.append(f"rep{ridx+1}\t" + '\t'.join('0' for _ in range(n_hap-1)))
                                    continue
                                arr = list(map(float, sfs_r))
                                if len(arr) < n_hap + 1:
                                    arr = arr + [0.0] * (n_hap + 1 - len(arr))
                                vals = arr[1:n_hap]
                                if args.sfs_normalized:
                                    denom = sum(vals)
                                    if denom > 0:
                                        vals_out = [v/denom for v in vals]
                                    else:
                                        vals_out = [0.0] * len(vals)
                                else:
                                    vals_out = [int(v) for v in vals]
                                per_rep_lines.append(f"rep{ridx+1}\t" + '\t'.join(str(x) for x in vals_out))
                            sfs_text = header + '\n' + '\n'.join(per_rep_lines) + '\n'
                            # write and finish
                            if isinstance(args.sfs, str) and args.sfs is not True and args.sfs != 'True':
                                sfs_path = args.sfs if args.sfs.endswith('.sfs') else args.sfs + '.sfs'
                            else:
                                if getattr(args, 'out', None) and args.out not in (None, '-'):
                                    base = args.out
                                    for ext_candidate in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                                        if base.endswith(ext_candidate):
                                            base = base[:-len(ext_candidate)]; break
                                    sfs_path = base + '.sfs'
                                else:
                                    sfs_path = default_stem + '.sfs'
                            # Defer creating the SFS file until the run fully completes.
                            sfs_pending = (sfs_path, sfs_text)
                            skip_sfs_loop = True
                        else:
                            # compute mean across available stored_sfs
                            accum = None
                            used = 0
                            for sfs_r in stored_sfs:
                                if sfs_r is None:
                                    continue
                                arr = list(map(float, sfs_r))
                                if accum is None:
                                    # ensure proper length
                                    if len(arr) < n_hap + 1:
                                        arr = arr + [0.0] * (n_hap + 1 - len(arr))
                                    accum = [float(x) for x in arr[1:n_hap]]
                                    used = 1
                                else:
                                    if len(arr) < n_hap + 1:
                                        arr = arr + [0.0] * (n_hap + 1 - len(arr))
                                    vals = [float(x) for x in arr[1:n_hap]]
                                    accum = [a + v for a, v in zip(accum, vals)]
                                    used += 1
                            if accum is None or used == 0:
                                mean_vals = [0.0] * (n_hap - 1)
                            else:
                                mean_vals = [a / used for a in accum]
                            if args.sfs_normalized:
                                denom_m = sum(mean_vals)
                                if denom_m > 0:
                                    mean_vals = [v/denom_m for v in mean_vals]
                                else:
                                    mean_vals = [0.0] * len(mean_vals)
                            sfs_text = ('#SFS\t' + '\t'.join(str(k) for k in range(1, n_hap)) + '\n') + 'mean\t' + '\t'.join(str(x) for x in mean_vals) + '\n'
                            # write and finish
                            if isinstance(args.sfs, str) and args.sfs is not True and args.sfs != 'True':
                                sfs_path = args.sfs if args.sfs.endswith('.sfs') else args.sfs + '.sfs'
                            else:
                                if getattr(args, 'out', None) and args.out not in (None, '-'):
                                    base = args.out
                                    for ext_candidate in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                                        if base.endswith(ext_candidate):
                                            base = base[:-len(ext_candidate)]; break
                                    sfs_path = base + '.sfs'
                                else:
                                    sfs_path = default_stem + '.sfs'
                            # Defer creating the SFS file until the run fully completes.
                            sfs_pending = (sfs_path, sfs_text)
                            skip_sfs_loop = True
                    except Exception:
                        # fall through to other handling
                        pass

                # If stored_sfs is not available but we have captured simulator output
                # (concatenated_stdout), try to extract any per-replicate SFS markers
                # produced by in-process msprime sim ("#SFS_PER_REP\t...") and use
                # them directly. If none are present, fall back to parsing haplotype
                # blocks via the streaming/text parsers.
                # Preference order: worker-provided #SFS_AF (authoritative) ->
                # worker #SFS_PER_REP -> text-based parsers -> resimulate via msprime
                # (resimulate only when no worker SFS markers or text parsing succeed).
                force_resim_flag = False
                skip_sfs_loop = False
                if stored_sfs is None and (concatenated_stdout or concatenated_stderr):
                    try:
                        # First, look for explicit per-rep SFS markers written by _ts_to_ms_like
                        try:
                            import re as _re
                            # Prefer allele-frequency-spectrum markers (exact from msprime)
                            combined_out = (concatenated_stdout or '') + '\n' + (concatenated_stderr or '')
                            sfs_af_lines = _re.findall(r'(?m)^#SFS_AF\s+(.+)$', combined_out)
                            sfs_lines = _re.findall(r'(?m)^#SFS_PER_REP\s+(.+)$', combined_out)
                        except Exception:
                            sfs_lines = []
                            sfs_af_lines = []
                        # If we found AF markers, treat them as stored_sfs with float values
                        if sfs_af_lines:
                            parsed_af = []
                            for ln in sfs_af_lines:
                                try:
                                    vals = [float(x) for x in ln.strip().split() if x is not None and x != '']
                                except Exception:
                                    vals = []
                                parsed_af.append(vals)
                            stored_sfs = parsed_af
                            # If workers provided exact allele-frequency-spectrum
                            # values, treat them as authoritative and do not
                            # re-simulate, even if force-resim was requested.
                            try:
                                force_resim_flag = False
                            except Exception:
                                pass
                        elif sfs_lines:
                            parsed = []
                            for ln in sfs_lines:
                                try:
                                    vals = [float(x) for x in ln.strip().split() if x is not None and x != '']
                                except Exception:
                                    vals = []
                                parsed.append(vals)
                            stored_sfs = parsed
                        

                        if stored_sfs and not force_resim_flag:
                            # If we obtained stored_sfs from the stdout markers, reuse the
                            # existing stored_sfs-handling above by letting the code fall
                            # through (it will detect stored_sfs and format output).
                            pass
                        else:
                            # No explicit markers found: fall back to the text-based SFS parsers
                            if args.sfs_mode == 'mean':
                                # If available, derive n_hap from meta to optimize dtype; else let it infer
                                try:
                                    _pl = int(meta.get('ploidy') or 1)
                                    _hap_list = meta.get('counts_disc') or []
                                    _n_hap = int(sum(int(h) for h in _hap_list)) if _hap_list else None
                                except Exception:
                                    _n_hap = None
                                header_sfs, lines_sfs = compute_sfs_stream_mean(concatenated_stdout, n_hap_expected=_n_hap, normalized=args.sfs_normalized)
                            else:
                                header_sfs, lines_sfs = compute_sfs_fast(concatenated_stdout, normalized=args.sfs_normalized, mode=args.sfs_mode)
                            sfs_text = header_sfs + '\n' + '\n'.join(lines_sfs) + '\n'
                            # write file and exit early from --sfs branch
                            if isinstance(args.sfs, str) and args.sfs is not True and args.sfs != 'True':
                                if args.sfs.endswith('.sfs'):
                                    sfs_path = args.sfs
                                else:
                                    sfs_path = args.sfs + '.sfs'
                            else:
                                if getattr(args, 'out', None) and args.out not in (None, '-'):
                                    base = args.out
                                    for ext_candidate in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                                        if base.endswith(ext_candidate):
                                            base = base[:-len(ext_candidate)]; break
                                    sfs_path = base + '.sfs'
                                else:
                                    sfs_path = default_stem + '.sfs'
                            try:
                                with open(sfs_path, 'w') as f:
                                    f.write('# SFS output\n')
                                    f.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                                    f.write(sfs_text)
                                    try:
                                        f.flush()
                                        os.fsync(f.fileno())
                                    except Exception:
                                        pass
                                if not sfs_written:
                                    sys.stderr.write(f"# INFO: wrote SFS to {sfs_path}\n")
                                    sfs_written = True
                            except Exception as e:
                                sys.stderr.write(f"# ERROR: failed to write SFS file {sfs_path}: {e}\n")
                                # Per user preference, never dump the full SFS to terminal stdout.
                                # Do not print sfs_text; only emit an error summary.
                                if not args.out or args.out == '-':
                                    sys.stderr.write(f"# ERROR: SFS not written to file and will not be printed to terminal.\n")
                            # skip the rest of msprime in-process SFS logic
                            # Only skip when not forcing resimulation
                            if not force_resim_flag:
                                skip_sfs_loop = True
                            stored_sfs = []
                    except Exception:
                        # fallback to re-simulate below if parsing fails
                        stored_sfs = None

                if not skip_sfs_loop:
                    for ridx in range(reps):
                        ts = None
                        sfs_raw = None
                        # If we have stored per-rep SFS from the primary loop, use it
                        if stored_sfs is not None and ridx < len(stored_sfs):
                            sfs_raw = stored_sfs[ridx]
                            if sfs_raw is None:
                                # previously failed to compute SFS for this rep
                                if args.sfs_mode == 'per-rep':
                                    per_rep_lines.append(f"rep{ridx+1}\t" + '\t'.join('0' for _ in range(n_hap-1)))
                                continue
                        else:
                            # No stored SFS available: run simulate to obtain ts and compute SFS
                            try:
                                ts = eng.simulate(model_obj, chosen_contig, samples_per_pop)
                                ts = _ensure_mutations_on_ts(ts)
                                try:
                                    import msprime as _msp
                                    _num_sites = int(getattr(ts, 'num_sites', 0))
                                    if _num_sites == 0:
                                        _dg = bool(getattr(ts, 'discrete_genome', False)) if hasattr(ts, 'discrete_genome') else False
                                        _mu = mu_local if isinstance(mu_local, (int,float)) and mu_local>0 else 1e-8
                                        ts = _msp.sim_mutations(ts, rate=_mu, discrete_genome=_dg)
                                except Exception:
                                    pass
                            except Exception as e:
                                if args.sfs_mode == 'per-rep':
                                    per_rep_lines.append(f"rep{ridx+1}\t" + '\t'.join('0' for _ in range(n_hap-1)))
                                continue
                            try:
                                # Announce when we compute via allele_frequency_spectrum
                                if force_resim_flag:
                                    try:
                                        sys.stderr.write(f"#DEBUG: force-resim computing allele_frequency_spectrum for rep {ridx+1}\n")
                                    except Exception:
                                        pass
                                sfs_raw = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
                            except Exception as e:
                                if args.sfs_mode == 'per-rep':
                                    per_rep_lines.append(f"rep{ridx+1}\t" + '\t'.join('0' for _ in range(n_hap-1)))
                                continue

                        # Process this replicate's SFS values (ensure this runs inside the loop)
                        try:
                            arr = list(map(float, sfs_raw))
                        except Exception:
                            arr = []
                        if len(arr) < n_hap + 1:
                            arr = arr + [0.0] * (n_hap + 1 - len(arr))
                        vals = arr[1:n_hap]
                        if args.sfs_normalized:
                            denom = sum(vals)
                            if denom > 0:
                                vals_out = [v/denom for v in vals]
                            else:
                                vals_out = [0.0] * len(vals)
                        else:
                            # non-normalized per-rep: present integer counts like other engines
                            vals_out = [int(v) for v in vals]

                        if args.sfs_mode == 'per-rep':
                            per_rep_lines.append(f"rep{ridx+1}\t" + '\t'.join(str(x) for x in vals_out))
                        else:
                            accum = [a + float(b) for a, b in zip(accum, vals_out)]
                            used += 1

                # finalize text
                if args.sfs_mode == 'per-rep':
                    sfs_text = header + '\n' + '\n'.join(per_rep_lines) + '\n'
                else:
                    if used == 0:
                        mean_vals = [0.0] * (n_hap - 1)
                    else:
                        mean_vals = [a / used for a in accum]
                    # if normalized was requested, ensure mean sums to 1 (re-normalize mean of frequencies if needed)
                    if args.sfs_normalized:
                        denom_m = sum(mean_vals)
                        if denom_m > 0:
                            mean_vals = [v/denom_m for v in mean_vals]
                        else:
                            mean_vals = [0.0] * len(mean_vals)
                    sfs_text = header + '\n' + 'mean\t' + '\t'.join(str(x) for x in mean_vals) + '\n'
        else:
            # Use accelerated implementation if NumPy is available; falls back automatically.
            if args.sfs_mode == 'mean':
                try:
                    _hap_list = meta.get('counts_disc') or []
                    _n_hap = int(sum(int(h) for h in _hap_list)) if _hap_list else None
                except Exception:
                    _n_hap = None
                header_sfs, lines_sfs = compute_sfs_stream_mean(concatenated_stdout, n_hap_expected=_n_hap, normalized=args.sfs_normalized)
            else:
                header_sfs, lines_sfs = compute_sfs_fast(concatenated_stdout, normalized=args.sfs_normalized, mode=args.sfs_mode)
            sfs_text = header_sfs + '\n' + '\n'.join(lines_sfs) + '\n'
        # Determine SFS file path
        if isinstance(args.sfs, str) and args.sfs is not True and args.sfs != 'True':
            # User provided custom base/path; honor it (append .sfs if missing)
            if args.sfs.endswith('.sfs'):
                sfs_path = args.sfs
            else:
                sfs_path = args.sfs + '.sfs'
        else:
            # No custom SFS name: if --output supplied (and not '-') derive from it; else use default stem
            if getattr(args, 'out', None) and args.out not in (None, '-'):
                # Strip known output extensions then add .sfs
                base = args.out
                for ext_candidate in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                    if base.endswith(ext_candidate):
                        base = base[:-len(ext_candidate)]
                        break
                sfs_path = base + '.sfs'
            else:
                sfs_path = default_stem + '.sfs'
        try:
            with open(sfs_path, 'w') as f:
                f.write('# SFS output\n')
                f.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                f.write(sfs_text)
                try:
                    f.flush()
                    os.fsync(f.fileno())
                except Exception:
                    pass
            if not sfs_written:
                sys.stderr.write(f"# INFO: wrote SFS to {sfs_path}\n")
                sfs_written = True
        except Exception as e:
            sys.stderr.write(f"# ERROR: failed to write SFS file {sfs_path}: {e}\n")
            # Do not print the SFS to terminal stdout per user request.
            if not args.out or args.out == '-':
                sys.stderr.write(f"# ERROR: SFS not written to file and will not be printed to terminal.\n")
        # If user requested only SFS (no explicit --output provided), exit now
        # EXCEPT when a paired neutral run is requested, in which case we continue
        # so we can also generate the neutral paired output/SFS.
        # Special case: if the user explicitly provided --output-format but did not
        # provide --output, create a primary simulation file named after the SFS
        # basename (if given) or the default stem. This allows "--sfs --output-format vcf"
        # to produce both the SFS and the simulation VCF named consistently.
        argv_full_main = sys.argv[1:]
        output_format_flag = any(a == '--output-format' or a.startswith('--output-format=') for a in argv_full_main)
        if output_format_flag and not getattr(args, '_user_out_provided', False):
            # Prefer using the --sfs basename for primary output if the user
            # supplied a string value for --sfs (argparse sets it to True when
            # the flag is provided without a value). If --sfs ends with '.sfs',
            # strip that and replace with the desired output extension. Otherwise
            # fall back to the default stem.
            ext_map_tmp = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
            desired_ext_tmp = ext_map_tmp.get(args.format, '.ms')
            chosen_base = None
            try:
                if isinstance(args.sfs, str) and args.sfs:
                    chosen_base = args.sfs
                    # If user gave a trailing .sfs, remove it so we can append the
                    # correct extension. Do not strip other extensions to avoid
                    # surprising behavior when the user intentionally provided them.
                    if chosen_base.endswith('.sfs'):
                        chosen_base = chosen_base[:-4]
            except Exception:
                chosen_base = None
            if chosen_base:
                args.out = chosen_base + desired_ext_tmp
            else:
                args.out = default_stem + desired_ext_tmp
            # primary output is written below.
            setattr(args, '_user_out_provided', True)
            # Ensure we won't suppress primary output later
            try:
                setattr(args, '_suppress_primary_output', False)
            except Exception:
                pass
            try:
                sys.stderr.write(f"# INFO: --output not provided but --output-format detected; writing primary output to {args.out}\n")
            except Exception:
                pass

        if (
            getattr(args, 'sfs', False)
            and not getattr(args, '_user_out_provided', False)
            and not getattr(args, 'paired_neutral', False)
        ):
            try:
                sys.stderr.flush()
                sys.stdout.flush()
            except Exception:
                pass
    # Decide if we should suppress primary output (SFS-only scenario)
    # If user did not supply --output but requested --sfs, treat as SFS-only and do not create primary ms/vcf file.
    write_primary_output = True
    if (not getattr(args, '_user_out_provided', False)) and getattr(args, 'sfs', False):
        write_primary_output = False
        setattr(args, '_suppress_primary_output', True)
    # But if --output-format was explicitly set, we intend to write the primary output too.
    try:
        _argv_all = sys.argv[1:]
        _oflag = any(a == '--output-format' or a.startswith('--output-format=') for a in _argv_all)
        if _oflag and not getattr(args, '_user_out_provided', False):
            write_primary_output = True
            setattr(args, '_suppress_primary_output', False)
    except Exception:
        pass
    # Handle primary output: if args.out is None -> write file with default stem; if '-' -> stdout
    if args.out is None and write_primary_output:
        fmt_ext = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
        ext = fmt_ext.get(args.format, '.ms')
        args.out = default_stem + ext
        try:
            # Only inform when the user did not explicitly provide --output
            if not getattr(args, '_user_out_provided', False) and args.out and args.out != '-':
                sys.stderr.write(f"# INFO: no --output provided; primary sweep output will be written to {args.out}\n")
        except Exception:
            pass
    if args.out == '-':
        if concatenated_stdout:
            sys.stdout.write(concatenated_stdout)
    elif not args.sfs and not concatenated_stdout and args.format in ('ms','ms.gz'):
        # nothing to write
        pass
    if concatenated_stderr:
        # Filter out worker-emitted SFS/debug marker lines so we don't print
        # the full SFS or debug site counts to the terminal. Keep other stderr.
        try:
            filtered_lines = []
            for ln in concatenated_stderr.splitlines(True):
                stripped = ln.lstrip()
                if stripped.startswith('#SFS_AF') or stripped.startswith('#SFS_PER_REP') or stripped.startswith('#DEBUG_NUM_SITES'):
                    # swallow marker line
                    continue
                filtered_lines.append(ln)
            if filtered_lines:
                sys.stderr.write(''.join(filtered_lines))
        except Exception:
            # on any issue, fall back to writing the original stderr
            sys.stderr.write(concatenated_stderr)
    if write_primary_output and args.out and args.out != '-':
        # Normalize output filename extension based on requested format if user omitted one.
        ext_map = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
        desired_ext = ext_map.get(args.format, '')
        if desired_ext and not args.out.endswith(desired_ext):
            # Strip any recognized extension (handle multi-part .vcf.gz first) then append desired
            for kext in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                if args.out.endswith(kext):
                    args.out = args.out[:-len(kext)]
                    break
            args.out = args.out + desired_ext
        fmt = args.format
        # Safety: if the engine produced ms-like stdout and user requested ms output,
        # ensure the output filename uses the correct .ms / .ms.gz extension.
        try:
            if concatenated_stdout and concatenated_stdout.lstrip().startswith('ms ') and fmt in ('ms', 'ms.gz'):
                desired_ext = '.ms.gz' if fmt == 'ms.gz' else '.ms'
                if args.out and not args.out.endswith(desired_ext):
                    for kext in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                        if args.out.endswith(kext):
                            args.out = args.out[:-len(kext)]
                            break
                    args.out = args.out + desired_ext
        except Exception:
            pass
        # Write primary output first
        if fmt == 'bcf':
            # Produce per-replicate VCF texts then convert each to BCF
            def _build_vcf_meta_lines_bcf():
                info_pairs = []
                info_pairs.append(('engine', meta.get('engine', engine)))
                info_pairs.append(('species', meta.get('species_id')))
                info_pairs.append(('model', model_id))
                info_pairs.append(('pop0', meta.get('pop0')))
                if meta.get('pop_order'):
                    info_pairs.append(('order', '|'.join(map(str, meta.get('pop_order')))))
                if meta.get('counts_disc'):
                    hap_list = meta.get('counts_disc') or []
                    pl = meta.get('ploidy', 1) or 1
                    try:
                        sample_list = [int(h // pl) for h in hap_list]
                    except Exception:
                        sample_list = hap_list
                    info_pairs.append(('sample_counts', '|'.join(map(str, sample_list))))
                    info_pairs.append(('hap_counts', '|'.join(map(str, hap_list))))
                theta_val = rho_val = None
                try:
                    toks_tmp = raw_cmd_line.split()
                    if '-t' in toks_tmp:
                        theta_val = toks_tmp[toks_tmp.index('-t')+1]
                    if '-r' in toks_tmp:
                        rho_val = toks_tmp[toks_tmp.index('-r')+1]
                except Exception:
                    pass
                if theta_val: info_pairs.append(('theta', theta_val))
                if rho_val: info_pairs.append(('rho', rho_val))
                info_pairs.append(('length', meta.get('length')))
                info_pairs.append(('N0', meta.get('N0')))
                info_pairs.append(('mu', meta.get('mu')))
                info_pairs.append(('r', meta.get('rrate')))
                info_pairs.append(('ploidy', meta.get('ploidy')))
                if meta.get('chromosome') is not None:
                    info_pairs.append(('chromosome', meta.get('chromosome')))
                if meta.get('sweep_pos') is not None:
                    sp = meta.get('sweep_pos')
                    info_pairs.append(('sweep_pos', sp))
                    try:
                        if sp is not None and meta.get('length'):
                            sweep_bp = int(round(float(sp) * int(meta.get('length'))))
                            info_pairs.append(('sweep_bp', sweep_bp))
                    except Exception:
                        pass
                if meta.get('sel_2Ns') is not None:
                    info_pairs.append(('sel_2Ns', meta.get('sel_2Ns')))
                if meta.get('sweep_time') is not None:
                    info_pairs.append(('sweep_time', meta.get('sweep_time')))
                if meta.get('fixation_time') is not None:
                    info_pairs.append(('fixation_time', meta.get('fixation_time')))
                kv_dict = {k:v for k,v in info_pairs if v is not None}
                core_keys = ['engine','species','model','pop0','order','sample_counts','hap_counts']
                param_keys = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                sel_keys = ['sweep_pos','sweep_bp','sel_2Ns','sweep_time','fixation_time']
                def build_group_line(prefix, keys):
                    vals = [f"{k}={kv_dict[k]}" for k in keys if k in kv_dict]
                    if not vals:
                        return None
                    return '# ' + prefix + ': ' + ', '.join(vals)
                meta_comment_lines_local = []
                for grp, keys in (('core', core_keys), ('params', param_keys), ('selection', sel_keys)):
                    ln = build_group_line(grp, keys)
                    if ln:
                        meta_comment_lines_local.append(ln)
                return meta_comment_lines_local

            per_rep_vcfs = []
            if meta.get('engine') == 'msprime' and meta.get('produced_format') == 'vcf' and concatenated_stdout:
                vtxt = concatenated_stdout
                blocks = []
                if '##fileformat' in vtxt:
                    parts = vtxt.split('##fileformat')
                    for p in parts:
                        p = p.strip()
                        if not p:
                            continue
                        blocks.append('##fileformat' + p if not p.startswith('##fileformat') else p)
                else:
                    blocks = [vtxt]
                for vb in blocks:
                    try:
                        vb = _adjust_vcf_contig_length(vb, meta.get('chromosome'), meta.get('length'))
                    except Exception:
                        pass
                    vcf_lines = vb.splitlines()
                    meta_lines = _build_vcf_meta_lines_bcf()
                    try:
                        header_idx = next(i for i,l in enumerate(vcf_lines) if l.startswith('#CHROM'))
                    except StopIteration:
                        header_idx = len(vcf_lines)
                    vcf_meta_lines = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines]
                    vcf_lines_final = vcf_lines[:header_idx] + vcf_meta_lines + vcf_lines[header_idx:]
                    per_rep_vcfs.append('\n'.join(vcf_lines_final) + ('\n' if vcf_lines_final and not vcf_lines_final[-1].endswith('\n') else ''))
            else:
                # Split ms-like stdout into replicate blocks
                text = concatenated_stdout or ''
                blocks = []
                cur = []
                for ln in (text.splitlines(True)):
                    if ln.lstrip().startswith('//'):
                        if cur:
                            blocks.append(''.join(cur))
                            cur = []
                        continue
                    if ln.lstrip().startswith('ms ') or ln.lstrip().startswith('msms ') or ln.lstrip().startswith('discoal ') or ln.lstrip().startswith('scrm '):
                        continue
                    cur.append(ln)
                if cur:
                    blocks.append(''.join(cur))
                seg_hdr_pat = re.compile(r'^\s*segsites\s*:', re.IGNORECASE | re.MULTILINE)
                for blk in blocks:
                    # Skip blocks that don't contain any segsites header
                    if not seg_hdr_pat.search(blk or ''):
                        continue
                    vtxt = ms_like_to_vcf(blk, meta['length'], chrom=meta.get('chromosome') or 'chr1', ploidy=meta.get('ploidy',1))
                    vcf_lines = vtxt.splitlines()
                    meta_lines = _build_vcf_meta_lines_bcf()
                    try:
                        header_idx = next(i for i,l in enumerate(vcf_lines) if l.startswith('#CHROM'))
                    except StopIteration:
                        header_idx = len(vcf_lines)
                    vcf_meta_lines = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines]
                    vcf_lines_final = vcf_lines[:header_idx] + vcf_meta_lines + vcf_lines[header_idx:]
                    per_rep_vcfs.append('\n'.join(vcf_lines_final) + ('\n' if vcf_lines_final and not vcf_lines_final[-1].endswith('\n') else ''))

            # Write BCFs
            if args.reps and args.reps > 1 and args.out and args.out != '-':
                parent = os.path.dirname(args.out)
                base = os.path.basename(args.out)
                # derive folder name by stripping known extensions
                folder_base = base
                if base.endswith('.vcf.gz'):
                    folder_base = base[:-7]
                elif base.endswith('.vcf') or base.endswith('.bcf'):
                    folder_base = base[:-4]
                out_dir = os.path.join(parent, folder_base) if parent else folder_base
                try:
                    os.makedirs(out_dir, exist_ok=True)
                except Exception:
                    pass
                root = folder_base
                for i, vtxt in enumerate(per_rep_vcfs, start=1):
                    path_out = os.path.join(out_dir, f"{root}_{i}.bcf")
                    tmp_vcf = None
                    try:
                        tmpdir = _ensure_local_tmpdir()
                        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.vcf', dir=tmpdir) as tf:
                            tf.write(vtxt)
                            tmp_vcf = tf.name
                        _vcf_text_to_bcf(tmp_vcf, path_out)
                    finally:
                        try:
                            if tmp_vcf is not None and os.path.exists(tmp_vcf):
                                os.unlink(tmp_vcf)
                        except Exception:
                            pass
            else:
                tmp_vcf = None
                try:
                    tmpdir = _ensure_local_tmpdir()
                    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.vcf', dir=tmpdir) as tf:
                        tf.write(per_rep_vcfs[0] if per_rep_vcfs else '')
                        tmp_vcf = tf.name
                    _vcf_text_to_bcf(tmp_vcf, args.out)
                finally:
                    try:
                        if tmp_vcf is not None and os.path.exists(tmp_vcf):
                            os.unlink(tmp_vcf)
                    except Exception:
                        pass
        elif fmt.startswith('vcf'):
            # Build metadata lines for header injection
            def _build_vcf_meta_lines():
                info_pairs = []
                info_pairs.append(('engine', meta.get('engine', engine)))
                info_pairs.append(('species', meta.get('species_id')))
                info_pairs.append(('model', model_id))
                info_pairs.append(('pop0', meta.get('pop0')))
                if meta.get('pop_order'):
                    info_pairs.append(('order', '|'.join(map(str, meta.get('pop_order')))))
                if meta.get('counts_disc'):
                    hap_list = meta.get('counts_disc') or []
                    pl = meta.get('ploidy', 1) or 1
                    try:
                        sample_list = [int(h // pl) for h in hap_list]
                    except Exception:
                        sample_list = hap_list
                    info_pairs.append(('sample_counts', '|'.join(map(str, sample_list))))
                    info_pairs.append(('hap_counts', '|'.join(map(str, hap_list))))
                theta_val = rho_val = None
                try:
                    toks_tmp = raw_cmd_line.split()
                    if '-t' in toks_tmp:
                        theta_val = toks_tmp[toks_tmp.index('-t')+1]
                    if '-r' in toks_tmp:
                        rho_val = toks_tmp[toks_tmp.index('-r')+1]
                except Exception:
                    pass
                if theta_val: info_pairs.append(('theta', theta_val))
                if rho_val: info_pairs.append(('rho', rho_val))
                info_pairs.append(('length', meta.get('length')))
                info_pairs.append(('N0', meta.get('N0')))
                info_pairs.append(('mu', meta.get('mu')))
                info_pairs.append(('r', meta.get('rrate')))
                info_pairs.append(('ploidy', meta.get('ploidy')))
                if meta.get('chromosome') is not None:
                    info_pairs.append(('chromosome', meta.get('chromosome')))
                if meta.get('sweep_pos') is not None:
                    sp = meta.get('sweep_pos')
                    info_pairs.append(('sweep_pos', sp))
                    try:
                        if sp is not None and meta.get('length'):
                            sweep_bp = int(round(float(sp) * int(meta.get('length'))))
                            info_pairs.append(('sweep_bp', sweep_bp))
                    except Exception:
                        pass
                if meta.get('sel_2Ns') is not None:
                    info_pairs.append(('sel_2Ns', meta.get('sel_2Ns')))
                if meta.get('sweep_time') is not None:
                    info_pairs.append(('sweep_time', meta.get('sweep_time')))
                if meta.get('fixation_time') is not None:
                    info_pairs.append(('fixation_time', meta.get('fixation_time')))
                kv_dict = {k:v for k,v in info_pairs if v is not None}
                core_keys = ['engine','species','model','pop0','order','sample_counts','hap_counts']
                param_keys = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                sel_keys = ['sweep_pos','sweep_bp','sel_2Ns','sweep_time','fixation_time']
                def build_group_line(prefix, keys):
                    vals = [f"{k}={kv_dict[k]}" for k in keys if k in kv_dict]
                    if not vals:
                        return None
                    return '# ' + prefix + ': ' + ', '.join(vals)
                meta_comment_lines_local = []
                for grp, keys in (('core', core_keys), ('params', param_keys), ('selection', sel_keys)):
                    ln = build_group_line(grp, keys)
                    if ln:
                        meta_comment_lines_local.append(ln)
                return meta_comment_lines_local

            # Prepare per-replicate VCF texts
            per_rep_vcfs = []
            if meta.get('engine') == 'msprime' and meta.get('produced_format') == 'vcf' and concatenated_stdout:
                vtxt = concatenated_stdout
                blocks = []
                if '##fileformat' in vtxt:
                    parts = vtxt.split('##fileformat')
                    for p in parts:
                        p = p.strip()
                        if not p:
                            continue
                        blocks.append('##fileformat' + p if not p.startswith('##fileformat') else p)
                else:
                    blocks = [vtxt]
                for vb in blocks:
                    try:
                        vb = _adjust_vcf_contig_length(vb, meta.get('chromosome'), meta.get('length'))
                    except Exception:
                        pass
                    vcf_lines = vb.splitlines()
                    meta_lines = _build_vcf_meta_lines()
                    try:
                        header_idx = next(i for i,l in enumerate(vcf_lines) if l.startswith('#CHROM'))
                    except StopIteration:
                        header_idx = len(vcf_lines)
                    vcf_meta_lines = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines]
                    vcf_lines_final = vcf_lines[:header_idx] + vcf_meta_lines + vcf_lines[header_idx:]
                    per_rep_vcfs.append('\n'.join(vcf_lines_final) + ('\n' if vcf_lines_final and not vcf_lines_final[-1].endswith('\n') else ''))
            else:
                # Split ms-like stdout into replicate blocks by '//'
                text = concatenated_stdout or ''
                blocks = []
                cur = []
                for ln in (text.splitlines(True)):
                    if ln.lstrip().startswith('//'):
                        if cur:
                            blocks.append(''.join(cur))
                            cur = []
                        continue
                    # Skip echoed engine lines
                    if ln.lstrip().startswith('ms ') or ln.lstrip().startswith('msms ') or ln.lstrip().startswith('discoal ') or ln.lstrip().startswith('scrm '):
                        continue
                    cur.append(ln)
                if cur:
                    blocks.append(''.join(cur))
                seg_hdr_pat = re.compile(r'^\s*segsites\s*:', re.IGNORECASE | re.MULTILINE)
                for blk in blocks:
                    # Skip blocks that don't contain any segsites header
                    if not seg_hdr_pat.search(blk or ''):
                        continue
                    vtxt = ms_like_to_vcf(blk, meta['length'], chrom=meta.get('chromosome') or 'chr1', ploidy=meta.get('ploidy',1))
                    vcf_lines = vtxt.splitlines()
                    meta_lines = _build_vcf_meta_lines()
                    try:
                        header_idx = next(i for i,l in enumerate(vcf_lines) if l.startswith('#CHROM'))
                    except StopIteration:
                        header_idx = len(vcf_lines)
                    vcf_meta_lines = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines]
                    vcf_lines_final = vcf_lines[:header_idx] + vcf_meta_lines + vcf_lines[header_idx:]
                    per_rep_vcfs.append('\n'.join(vcf_lines_final) + ('\n' if vcf_lines_final and not vcf_lines_final[-1].endswith('\n') else ''))

            # If multiple replicates, write each to its own file under a folder named after --output
            if args.reps and args.reps > 1 and args.out and args.out != '-':
                parent = os.path.dirname(args.out)
                base = os.path.basename(args.out)
                # derive folder name by stripping known extensions
                folder_base = base
                if fmt == 'vcf.gz' and base.endswith('.vcf.gz'):
                    folder_base = base[:-7]
                elif fmt == 'vcf' and base.endswith('.vcf'):
                    folder_base = base[:-4]
                elif fmt == 'bcf' and base.endswith('.bcf'):
                    folder_base = base[:-4]
                out_dir = os.path.join(parent, folder_base) if parent else folder_base
                try:
                    os.makedirs(out_dir, exist_ok=True)
                except Exception:
                    pass
                root = folder_base
                if root.endswith('.vcf.gz'):
                    root = root[:-7]
                    extw = '.vcf.gz'
                elif root.endswith('.vcf'):
                    root = root[:-4]
                    extw = '.vcf'
                else:
                    root = root[:-4] if root.endswith('.bcf') else root
                    extw = '.bcf' if fmt == 'bcf' else ('.vcf.gz' if fmt=='vcf.gz' else '.vcf')
                for i, vtxt in enumerate(per_rep_vcfs, start=1):
                    fname = f"{root}_{i}{extw}"
                    path_out = os.path.join(out_dir, fname)
                    if fmt == 'vcf.gz':
                        _write_bgzip_text(path_out, vtxt)
                    elif fmt == 'bcf':
                        tmp_vcf = None
                        try:
                            tmpdir = _ensure_local_tmpdir()
                            with tempfile.NamedTemporaryFile('w', delete=False, suffix='.vcf', dir=tmpdir) as tf:
                                tf.write(vtxt)
                                tmp_vcf = tf.name
                            _vcf_text_to_bcf(tmp_vcf, path_out)
                        finally:
                            try:
                                if tmp_vcf is not None and os.path.exists(tmp_vcf):
                                    os.unlink(tmp_vcf)
                            except Exception:
                                pass
                    else:
                        with open(path_out, 'w') as f:
                            f.write(vtxt)
            else:
                # Single output file
                vcf_text = per_rep_vcfs[0] if per_rep_vcfs else ''
                if fmt == 'vcf.gz':
                    _write_bgzip_text(args.out, vcf_text)
                else:
                    with open(args.out, 'w') as f:
                        f.write(vcf_text)
        else:
                f_open = gzip.open if fmt == 'ms.gz' else open
                mode = 'wt' if fmt == 'ms.gz' else 'w'
                # Keep file open for the entire write sequence so subsequent writes don't target a closed file.
                with f_open(args.out, mode) as f:
                    f.write(raw_cmd_line + '\n')
                    # Build unified comma-separated info line
                    info_pairs = []
                    info_pairs.append(('engine', meta.get('engine', engine)))
                    info_pairs.append(('species', meta.get('species_id')))
                    info_pairs.append(('model', model_id))
                    info_pairs.append(('pop0', meta.get('pop0')))
                    if meta.get('pop_order'):
                        info_pairs.append(('order', '|'.join(map(str, meta.get('pop_order')))))
                    if meta.get('counts_disc'):
                        hap_list = meta.get('counts_disc') or []
                        pl = meta.get('ploidy', 1) or 1
                        try:
                            sample_list = [int(h // pl) for h in hap_list]
                        except Exception:
                            sample_list = hap_list
                        info_pairs.append(('sample_counts', '|'.join(map(str, sample_list))))
                        info_pairs.append(('hap_counts', '|'.join(map(str, hap_list))))
                    # parse theta/rho from command
                    theta_val = rho_val = None
                    try:
                        toks_tmp = raw_cmd_line.split()
                        if '-t' in toks_tmp:
                            theta_val = toks_tmp[toks_tmp.index('-t')+1]
                        if '-r' in toks_tmp:
                            rho_val = toks_tmp[toks_tmp.index('-r')+1]
                    except Exception:
                        pass
                    if theta_val: info_pairs.append(('theta', theta_val))
                    if rho_val: info_pairs.append(('rho', rho_val))
                    info_pairs.append(('length', meta.get('length')))
                    info_pairs.append(('N0', meta.get('N0')))
                    info_pairs.append(('mu', meta.get('mu')))
                    info_pairs.append(('r', meta.get('rrate')))
                    info_pairs.append(('ploidy', meta.get('ploidy')))
                    if meta.get('chromosome') is not None:
                        info_pairs.append(('chromosome', meta.get('chromosome')))
                    if meta.get('sweep_pos') is not None:
                        sp = meta.get('sweep_pos')
                        info_pairs.append(('sweep_pos', sp))
                        try:
                            if sp is not None and meta.get('length'):
                                sweep_bp = int(round(float(sp) * int(meta.get('length'))))
                                info_pairs.append(('sweep_bp', sweep_bp))
                        except Exception: pass
                    if meta.get('sel_2Ns') is not None:
                        info_pairs.append(('sel_2Ns', meta.get('sel_2Ns')))
                    if meta.get('sweep_time') is not None:
                        info_pairs.append(('sweep_time', meta.get('sweep_time')))
                    if meta.get('fixation_time') is not None:
                        info_pairs.append(('fixation_time', meta.get('fixation_time')))
                    # demog intentionally omitted from metadata output per user request
                    # Group related keys: core, params, selection (enforcement line removed per user request)
                    core_keys = ['engine','species','model','pop0','order','sample_counts','hap_counts']
                    param_keys = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                    sel_keys = ['sweep_pos','sweep_bp','sel_2Ns','sweep_time','fixation_time']
                    # Enforcement diagnostics retained internally but not emitted.
                    kv_dict = {k:v for k,v in info_pairs if v is not None}
                    # Inject enforcement stats if present
                    enf = meta.get('enforcement_stats') or {}
                    for k,v in enf.items():
                        kv_dict[k] = v
                    def build_line(group_name, keys):
                        vals = [f"{k}={kv_dict[k]}" for k in keys if k in kv_dict]
                        if not vals:
                            return None
                        return '# ' + group_name + ': ' + ', '.join(vals)
                    lines_info = []
                    for grp_name, keys in (
                        ('core', core_keys),
                        ('params', param_keys),
                        ('selection', sel_keys),
                    ):
                        ln = build_line(grp_name, keys)
                        if ln:
                            lines_info.append(ln)
                    for ln in lines_info:
                        f.write(ln + '\n')
                    # Write body: if msprime produced VCF, convert to ms-like per replicate; else, write captured ms-like output.
                    if meta.get('engine') == 'msprime' and meta.get('produced_format') == 'vcf' and concatenated_stdout:
                        vcf_text_all = concatenated_stdout
                        # Build per-rep seed list from runtime chunk seeds
                        runtime = meta.get('runtime', {}) or {}
                        chunk_seeds = runtime.get('chunk_seeds') or []
                        chunk_sizes = runtime.get('chunk_sizes') or []
                        seed_lines = []
                        for ci, base_seed in enumerate(chunk_seeds):
                            if base_seed is None:
                                continue
                            reps_in_chunk = chunk_sizes[ci] if ci < len(chunk_sizes) else 1
                            for r in range(reps_in_chunk):
                                seed_lines.append(str(int(base_seed) + int(r)))
                        # Split concatenated VCFs on '##fileformat'
                        blocks = []
                        if '##fileformat' in vcf_text_all:
                            parts = vcf_text_all.split('##fileformat')
                            for p in parts:
                                p = p.strip()
                                if not p:
                                    continue
                                blocks.append('##fileformat' + ('\n' if not p.startswith('\n') else '') + p)
                        else:
                            blocks = [vcf_text_all]
                        # Convert each block to ms-like and write with replicate separators and seed lines
                        for ib, vb in enumerate(blocks):
                            vb_adj = vb
                            try:
                                vb_adj = _adjust_vcf_contig_length(vb_adj, meta.get('chromosome'), meta.get('length'))
                            except Exception:
                                pass
                            # Prepend per-rep seed (plain integer on its own line), if available
                            if ib < len(seed_lines):
                                f.write(seed_lines[ib] + '\n')
                            # Replicate separator
                            f.write('//\n')
                            try:
                                ms_blk = vcf_to_ms_like(vb_adj, meta.get('length'), ploidy=meta.get('ploidy', 1))
                            except Exception:
                                ms_blk = 'segsites: 0\n'
                            if not ms_blk.endswith('\n'):
                                ms_blk += '\n'
                            f.write(ms_blk)
                    else:
                        # seed lines not written in header; per-replicate seeds are emitted before each replicate block
                        # format line removed per user request
                        lines = (concatenated_stdout or '').splitlines(True)
                        for ln in lines:
                            stripped_out = ln.lstrip()
                            if stripped_out.startswith('#SFS') or stripped_out.startswith('#SFS_AF') or stripped_out.startswith('#SFS_PER_REP') or stripped_out.startswith('#DEBUG_NUM_SITES'):
                                continue
                            # Avoid duplicating echoed command lines. For ms and discoal we skip the engine line.
                            # For msms, some versions may echo both an 'msms ...' line and a compatibility 'ms ...' line;
                            # keep only the 'msms' line we already wrote (raw_cmd_line) and drop any stray 'ms ' line.
                            if ln.startswith(engine+' '):
                                continue
                            if engine == 'msms' and ln.startswith('ms '):
                                continue
                            f.write(ln)
    # If a paired neutral run is requested, release large primary output buffers now to free RAM before neutral simulation.
    if getattr(args, 'paired_neutral', False):
        try:
            # Keep raw_cmd_line and meta (needed for deriving neutral command), drop bulky stdout/stderr strings.
            concatenated_stdout = None  # type: ignore
            concatenated_stderr = None  # type: ignore
            # Also drop any intermediate oversized variables if present.
            if 'hi_stdout' in locals():
                hi_stdout = None  # type: ignore
            if 'lo_seg_list' in locals():
                lo_seg_list = None  # type: ignore
            if 'hi_seg_list' in locals():
                hi_seg_list = None  # type: ignore
            if 'final_seg_list' in locals():
                final_seg_list = None  # type: ignore
            import gc as _gc
            _gc.collect()
            # no debug printing
        except Exception:
            pass
    # Optional paired neutral run (clean implementation)
    if getattr(args, 'paired_neutral', False):
        fmt = args.format
        neut_engine = args.neutral_engine or engine
        # Always generate an independent seed for the neutral run
        try:
            import random as _rseed
            neut_seed = _rseed.randint(1,2**31-1)
        except Exception:
            neut_seed = None
        # Build / derive neutral command
        if neut_engine == engine and neut_engine in ('discoal','msms'):
            base_line = raw_cmd_line
            toks = shlex.split(base_line)

            def _strip_discoal(ts):
                """Remove discoal selection options: -ws <t>, -x <pos>, -a <2Ns>."""
                out = []
                i = 0
                flags = {'-ws': 1, '-x': 1, '-a': 1}
                while i < len(ts):
                    t = ts[i]
                    if t in flags:
                        i += 1 + flags[t]
                        continue
                    out.append(t)
                    i += 1
                return out

            def _strip_msms(ts):
                """Remove msms selection options: -SaA <a> -SAA <2a> -Sp <x> -Smark -SI <t> <npop> <freqs...> -SFC."""
                out = []
                i = 0
                # flags with 1 argument
                one_arg = {'-SaA', '-SAA', '-Sp'}
                # flags with variable args
                while i < len(ts):
                    t = ts[i]
                    if t in one_arg:
                        i += 2  # skip flag and its single value
                        continue
                    if t == '-Smark' or t == '-SFC':
                        i += 1
                        continue
                    if t == '-SI':
                        # structure: -SI <tstart> <npop> <f1> ... <fN>
                        j = i + 1
                        # need at least two more tokens for tstart and npop
                        if j + 1 < len(ts):
                            # skip tstart and npop
                            j += 2
                            try:
                                npop_val = int(ts[i + 2])
                            except Exception:
                                npop_val = 0
                            j += max(0, npop_val)
                        i = j
                        continue
                    out.append(t)
                    i += 1
                return out

            # Build a neutral command by stripping selection args from the primary command.
            try:
                if neut_engine == 'discoal':
                    stripped = _strip_discoal(toks)
                else:  # msms
                    stripped = _strip_msms(toks)
                neut_cmd = ' '.join(stripped)
            except Exception:
                neut_cmd = base_line
            # Create a neutral metadata dict from primary meta and clear selection fields.
            try:
                if isinstance(meta, dict):
                    neut_meta = dict(meta)
                else:
                    neut_meta = {}
            except Exception:
                neut_meta = {}
            neut_meta['engine'] = neut_engine
            # Clear selection-specific metadata for the neutral run
            for k in ('sel_args','sweep_pos','sel_2Ns','sel_s','sweep_time','fixation_time'):
                if k in neut_meta:
                    neut_meta[k] = None if not isinstance(neut_meta[k], list) else []
            # Defer creating the SFS file until the run fully completes (only if available).
            if 'sfs_path' in locals() and 'sfs_text' in locals():
                sfs_pending = (sfs_path, sfs_text)
        else:
            if neut_engine == 'discoal':
                neut_cmd, neut_meta = build_discoal_command(
                    species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    discoal_pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length, max_fold_per_step=args.max_fold_per_step,
                    sweep_time=(getattr(args,'sweep_time',None) if getattr(args,'sweep_time',None) is not None else getattr(args,'fix_time', None)), x=None, s=None, min_snps=None, chr_name=args.chromosome, disable_en_ladder=args.no_en_ladder,
                    seed=neut_seed, time_units=getattr(args,'time_units','gens')
                )
            elif neut_engine=='msms':
                neut_cmd, neut_meta = build_msms_command(species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length, max_fold_per_step=args.max_fold_per_step,
                    chr_name=args.chromosome, min_snps=None, disable_en_ladder=args.no_en_ladder, a=None, s=None, x=None, fixation_time=0.0, sweep_time=None,
                    seed=neut_seed, time_units=getattr(args,'time_units','gens'))
            elif neut_engine=='ms':
                neut_cmd, neut_meta = build_ms_command(species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length, max_fold_per_step=args.max_fold_per_step,
                    chr_name=args.chromosome, min_snps=None, disable_en_ladder=args.no_en_ladder,
                    seed=neut_seed)
            elif neut_engine=='scrm':
                neut_cmd, neut_meta = build_scrm_command(species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length, max_fold_per_step=args.max_fold_per_step,
                    chr_name=args.chromosome, min_snps=None, disable_en_ladder=args.no_en_ladder,
                    seed=neut_seed)
            elif neut_engine=='msprime':
                neut_cmd, neut_meta = build_msprime_command(species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length, max_fold_per_step=args.max_fold_per_step,
                    chr_name=args.chromosome, min_snps=None, disable_en_ladder=args.no_en_ladder,
                    seed=neut_seed)
            else:
                # Fallback to ms neutral if an unknown engine is somehow specified
                neut_cmd, neut_meta = build_ms_command(species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length, max_fold_per_step=args.max_fold_per_step,
                    chr_name=args.chromosome, min_snps=None, disable_en_ladder=args.no_en_ladder,
                    seed=neut_seed)
        # Extract raw neutral line
        raw_neut=None
        for ln in neut_cmd.splitlines():
            if ln.startswith(neut_engine+' '):
                raw_neut=ln.strip(); break
        if raw_neut is None:
            sys.stderr.write('# ERROR: failed to extract neutral command line; skipping neutral run.\n')
        else:
            # Run neutral simulation with its own progress bar if requested
            try:
                neut_tokens = shlex.split(raw_neut)
                total_reps_neut = int(neut_tokens[2])
            except Exception:
                # Fallback: treat as single run
                total_reps_neut = 1
            neut_threads = max(1, min(int(getattr(args, 'parallel', 1)), total_reps_neut))
            # Mirror primary run chunking: respect --parallel and --sims-per-work
            try:
                neut_threads = max(1, min(int(getattr(args, 'parallel', 1)), total_reps_neut))
            except Exception:
                neut_threads = max(1, int(getattr(args, 'parallel', 1) or 1))
            # Per-rep progress mode similar to primary: show progress when requested and multiple reps
            per_rep_mode_neut = bool(getattr(args, 'progress', False) and total_reps_neut > 1)

            if neut_engine == 'msprime':
                # msprime runs in-process per-replicate
                neut_chunks = [1] * total_reps_neut
            else:
                # If user provided --sims-per-work, use it (cap to total_reps_neut). Otherwise split into neut_threads parts
                spw_val = getattr(args, 'sims_per_work', None)
                if spw_val is None:
                    try:
                        base = max(1, total_reps_neut // neut_threads)
                    except Exception:
                        base = 1
                    rem = total_reps_neut - base * neut_threads
                    neut_chunks = [base] * neut_threads
                    if rem > 0:
                        neut_chunks[-1] = neut_chunks[-1] + rem
                else:
                    try:
                        spw = max(1, int(spw_val))
                    except Exception:
                        spw = 1
                    if spw >= total_reps_neut:
                        neut_chunks = [total_reps_neut]
                    else:
                        full = total_reps_neut // spw
                        rem = total_reps_neut % spw
                        neut_chunks = [spw] * full
                        if rem:
                            neut_chunks.append(rem)

            # Unique seeds per chunk for engines that support explicit seeding
            if neut_engine in ('discoal','msms','ms','msprime'):
                try:
                    import random as _rndN
                    unique_seeds_n = set(); seeds_for_chunks_n = []
                    while len(seeds_for_chunks_n) < len(neut_chunks):
                        sN = _rndN.randint(1, 2**31 - 1)
                        if sN in unique_seeds_n:
                            continue
                        unique_seeds_n.add(sN)
                        seeds_for_chunks_n.append(sN)
                except Exception:
                    seeds_for_chunks_n = [None] * len(neut_chunks)
            else:
                seeds_for_chunks_n = [None] * len(neut_chunks)
            # Prepare progress bar
            neut_bar = None
            neut_progress_done = 0
            if getattr(args, 'progress', False):
                if tqdm is not None:
                    neut_bar = tqdm(total=total_reps_neut, desc='neutral', unit='rep')
                else:
                    sys.stderr.write('# neutral progress: 0/' + str(total_reps_neut) + '\n')
            def run_neut_chunk(rc, chunk_idx=None, rep_index=None):
                toks_loc = neut_tokens[:]; toks_loc[2] = str(rc)
                # Force a unique seed for this replicate
                if neut_engine in ('discoal','msms','ms','msprime'):
                    seed_val_n = seeds_for_chunks_n[chunk_idx] if chunk_idx is not None and chunk_idx < len(seeds_for_chunks_n) else None
                    if neut_engine in ('discoal','msms'):
                        # Remove any pre-existing -seed to avoid duplicates, then add
                        if '-seed' in toks_loc:
                            try:
                                idxs = [i for i,t in enumerate(toks_loc) if t=='-seed']
                                for idx in reversed(idxs):
                                    del toks_loc[idx:idx+2]
                            except Exception:
                                pass
                        if seed_val_n is not None:
                            toks_loc.extend(['-seed', str(int(seed_val_n))])
                    elif neut_engine == 'ms':
                        # Intentionally avoid using ms -seeds here to prevent the
                        # Hudson 'seedms' external seed file from being created.
                        # We keep seeds_for_chunks_n for bookkeeping but do not
                        # append '-seeds' to the command line. ms will use its
                        # default RNG behavior (or the user can supply explicit
                        # seeding elsewhere).
                        pass
                    elif neut_engine == 'msprime':
                        # handled below in-process; nothing to append to toks_loc
                        pass
                # Special in-process handling for msprime neutral runs
                if neut_engine == 'msprime':
                    try:
                        eng = sps.get_engine('msprime')
                    except Exception as e:
                        return (rc, '', f'# ERROR: neutral msprime: failed to get engine: {e}\n', 1, rep_index)
                    # Prepare inputs from neut_meta
                    try:
                        sp = sps.get_species(neut_meta.get('species_id'))
                        model_obj = neut_meta.get('model_obj')
                        chr_name = neut_meta.get('chromosome') or getattr(args, 'chromosome', None)
                        length_n = int(neut_meta.get('length') or 0)
                        # Choose contig; try bounded first
                        try:
                            chosen_contig = sp.get_contig(chr_name, mutation_rate=getattr(model_obj, 'mutation_rate', None), left=0, right=max(0, length_n - 1))
                        except Exception:
                            try:
                                chosen_contig = sp.get_contig(chr_name)
                            except Exception:
                                chosen_contig = None
                        if chosen_contig is None:
                            return (rc, '', '# ERROR: neutral msprime: could not resolve contig; provide --chromosome.\n', 1, rep_index)
                        # Samples per pop from hap counts and ploidy
                        samples_per_pop = {}
                        pop_order = neut_meta.get('pop_order') or []
                        counts_disc = neut_meta.get('counts_disc') or []
                        pl = neut_meta.get('ploidy', 1) or 1
                        for i,pop in enumerate(pop_order):
                            try:
                                samples_per_pop[pop] = int((counts_disc[i] if i < len(counts_disc) else 0) // pl)
                            except Exception:
                                samples_per_pop[pop] = 0
                        # Seed kwargs
                        kwargs = {}
                        if 'seed_val_n' in locals() and seed_val_n is not None:
                            kwargs['seed'] = int(seed_val_n)
                        # Prefer model mutation rate if available
                        try:
                            mr = getattr(model_obj, 'mutation_rate', None)
                            if mr is not None:
                                try:
                                    chosen_contig.mutation_rate = float(mr)
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        ts = eng.simulate(model_obj, chosen_contig, samples_per_pop, **kwargs)
                        ts = _ensure_mutations_on_ts(ts)
                        buf = io.StringIO()
                        if hasattr(ts, 'write_ms'):
                            try:
                                import tskit as _tskit
                                _tskit.write_ms(ts, buf)
                            except Exception:
                                try:
                                    ts.write_ms(buf)
                                except Exception:
                                    return (rc, '', '# ERROR: neutral msprime: failed to write ms text from tree sequence.\n', 1, rep_index)
                            neut_meta['produced_format'] = 'ms'
                        elif hasattr(ts, 'write_vcf'):
                            ts.write_vcf(buf)
                            neut_meta['produced_format'] = 'vcf'
                        else:
                            return (rc, '', '# ERROR: neutral msprime: no writer available on tree sequence.\n', 1, rep_index)
                        return (rc, buf.getvalue(), '', None, rep_index)
                    except Exception as e:
                        tb = traceback.format_exc()
                        try:
                            sys.stderr.write(tb + '\n')
                        except Exception:
                            pass
                        return (rc, '', f'# ERROR: neutral msprime simulate failed: {e}\n', 1, rep_index)
                try:
                    # Prevent simulator subprocesses from spawning many OpenMP/BLAS threads
                    env_n = os.environ.copy()
                    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
                        env_n[k] = env_n.get(k, '1')
                    total_cpus_n = os.cpu_count() or 1
                    workers_n = max(1, int(getattr(args, 'parallel', 1)))
                    block_size_n = max(1, total_cpus_n // workers_n)
                    # map the neutral chunk index into the worker slots to avoid out-of-range CPU blocks
                    worker_slot_n = (chunk_idx or 0) % workers_n
                    start_cpu_n = worker_slot_n * block_size_n
                    end_cpu_n = min(total_cpus_n, start_cpu_n + block_size_n)
                    if start_cpu_n >= total_cpus_n:
                        cpu_list_n = list(range(total_cpus_n))
                    else:
                        cpu_list_n = list(range(start_cpu_n, end_cpu_n)) if end_cpu_n > start_cpu_n else [start_cpu_n % total_cpus_n]
                    pre_n = _make_affinity_preexec(cpu_list_n)
                    rloc = _run_and_register(toks_loc, capture_output=True, text=True, check=True, env=env_n, preexec_fn=pre_n)
                    return (rc, rloc.stdout, rloc.stderr, None, rep_index)
                except subprocess.CalledProcessError as e_c:
                    return (rc, e_c.stdout, e_c.stderr, e_c.returncode, rep_index)
            neut_results = []
            if len(neut_chunks) == 1 and not per_rep_mode_neut:
                rc, out, err, code, rep_idx = run_neut_chunk(neut_chunks[0], 0, 0)
                neut_results.append((rc, out, err, code, rep_idx))
                if neut_bar is not None:
                    neut_bar.update(total_reps_neut)
                elif getattr(args, 'progress', False) and total_reps_neut == 1 and tqdm is None:
                    sys.stderr.write('# neutral progress: 1/1\n')
            else:
                # Parallel execution
                with concurrent.futures.ThreadPoolExecutor(max_workers=len(neut_chunks) if not per_rep_mode_neut else neut_threads) as ex_n:
                    futs_n = []
                    if per_rep_mode_neut:
                        rep_counter_n = 0
                        for idxN, cN in enumerate(neut_chunks):
                            futs_n.append(ex_n.submit(run_neut_chunk, cN, idxN, rep_counter_n))
                            rep_counter_n += 1
                    else:
                        for idxN, cN in enumerate(neut_chunks):
                            futs_n.append(ex_n.submit(run_neut_chunk, cN, idxN, idxN))
                    for fut_n in concurrent.futures.as_completed(futs_n):
                        rcN, outN, errN, codeN, rep_idxN = fut_n.result()
                        neut_results.append((rcN, outN, errN, codeN, rep_idxN))
                        if neut_bar is not None:
                            try:
                                neut_bar.update(int(rcN if not per_rep_mode_neut else 1))
                            except Exception:
                                neut_bar.update(1)
                        elif getattr(args, 'progress', False):
                            neut_progress_done += (rcN if not per_rep_mode_neut else 1)
                            sys.stderr.write(f'# neutral progress: {neut_progress_done}/{total_reps_neut}\n')
            if neut_bar is not None:
                neut_bar.close()
            # Sort and concatenate stdout in replicate order
            try:
                neut_results.sort(key=lambda x: (x[4] if x[4] is not None else 0))
            except Exception:
                pass
            # Aggregate stdout / stderr
            neut_out_txt = ''.join([r[1] for r in neut_results])
            neut_err_txt = ''.join([r[2] for r in neut_results if r[2]])
            # Optional debug dump: write a short summary of neut_results and
            # the aggregated stdout to a debug file when RAISD_DEBUG_NEUT is set.
            try:
                if os.getenv('RAISD_DEBUG_NEUT'):
                    dbg_path = os.path.join(os.getcwd(), 'neutral_debug.log')
                    with open(dbg_path, 'w') as dbgf:
                        dbgf.write(f"# neutral_results entries: {len(neut_results)}\n")
                        for idx, r in enumerate(neut_results):
                            rc_val, out_val, err_val, code_val, rep_idx_val = r
                            dbgf.write(f"-- entry {idx}: rc={rc_val} code={code_val} rep_idx={rep_idx_val} out_len={len(out_val)} err_len={len(err_val)}\n")
                            # write first 1000 chars of stdout for inspection
                            excerpt = out_val[:1000].replace('\r', '')
                            dbgf.write(excerpt + ("\n...[truncated]\n" if len(out_val) > 1000 else "\n"))
                        dbgf.write('\n# aggregated neut_out_txt first lines:\n')
                        for ln in neut_out_txt.splitlines()[:200]:
                            dbgf.write(ln + '\n')
            except Exception:
                pass
            # Report any non-zero return codes
            for r in neut_results:
                if r[3]:
                    sys.stderr.write(f"# ERROR: neutral replicate returned code {r[3]} (engine={neut_engine}).\n")
            # Determine neutral output path or basename from --paired-neutral value (if provided)
            paired_val = getattr(args, 'paired_neutral', False)
            neut_out_path = None
            user_specified_neut = False
            if isinstance(paired_val, str) and paired_val.strip():
                # If the user explicitly provided a --paired-neutral name, use it
                # exactly as given (do not auto-append '_neutral'). This preserves
                # user intent for custom naming.
                neut_out_path = paired_val.strip()
                user_specified_neut = True
            # If no paired-neutral name was provided, derive one from the primary
            # output if the primary output will be written; otherwise, if the
            # primary output is auto-named (default_stem), derive from that.
            if neut_out_path is None and not getattr(args, '_suppress_primary_output', False):
                if getattr(args, 'out', None) and args.out != '-':
                    main_out = args.out
                    if main_out.endswith('.vcf.gz'):
                        root = main_out[:-7]; ext = '.vcf.gz'
                    else:
                        root, ext = os.path.splitext(main_out)
                    neut_out_path = root + '_neutral' + ext
                else:
                    # No explicit primary output; use default_stem and desired ext
                    ext_map_full = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
                    desired_ext = ext_map_full.get(fmt, '')
                    if desired_ext:
                        neut_out_path = default_stem + '_neutral' + desired_ext
            # Force neutral output format to mirror primary format even if user supplied a mismatched extension
            if neut_out_path and neut_out_path != '-':
                ext_map_full = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
                desired_ext = ext_map_full.get(fmt,'')
                if desired_ext:
                    if neut_out_path.endswith('.vcf.gz') and desired_ext == '.vcf.gz':
                        pass  # already correct
                    elif neut_out_path.endswith('.ms.gz') and desired_ext == '.ms.gz':
                        pass  # already correct
                    elif not neut_out_path.endswith(desired_ext):
                        # Strip any existing known extension then append desired
                        for kext in ('.vcf.gz','.vcf','.bcf','.ms.gz','.ms'):
                            if neut_out_path.endswith(kext):
                                neut_out_path = neut_out_path[:-len(kext)]
                                break
                        neut_out_path = neut_out_path + desired_ext
            # If user only wants SFS (no primary --output) and provided a value to --paired-neutral,
            # skip writing a neutral ms/vcf file; compute neutral SFS using the provided value as basename.
            skip_neut_file = (not getattr(args, '_user_out_provided', False)) and user_specified_neut and getattr(args, 'sfs', False)
            if skip_neut_file:
                # We'll still compute SFS below using neut_out_txt; adjust root for naming
                if neut_out_path:
                    neut_root = neut_out_path
                    if neut_root.endswith('.vcf.gz'):
                        neut_root = neut_root[:-7]
                    elif neut_root.endswith('.ms.gz'):
                        neut_root = neut_root[:-6]
                    else:
                        # remove single extension (e.g. .ms/.vcf/.bcf) if present
                        base_root, base_ext = os.path.splitext(neut_root)
                        if base_ext in ('.ms', '.vcf', '.bcf'):
                            neut_root = base_root
                else:
                    neut_root = 'neutral'
                sys.stderr.write('# INFO: skipping neutral output file (no --output provided); writing only neutral SFS.\n')
                try:
                    # If msprime produced VCF, convert to ms-like per replicate before SFS
                    ms_like_for_sfs = neut_out_txt
                    try:
                        if neut_meta.get('engine') == 'msprime' and neut_meta.get('produced_format') == 'vcf' and neut_out_txt:
                            vtxt = neut_out_txt
                            blocksC = []
                            if '##fileformat' in vtxt:
                                partsC = vtxt.split('##fileformat')
                                for p in partsC:
                                    p = p.strip()
                                    if not p:
                                        continue
                                    blocksC.append('##fileformat' + p if not p.startswith('##fileformat') else p)
                            else:
                                blocksC = [vtxt]
                            out_ms_chunks = []
                            for vb in blocksC:
                                vb_adj = vb
                                try:
                                    vb_adj = _adjust_vcf_contig_length(vb_adj, neut_meta.get('chromosome'), neut_meta.get('length'))
                                except Exception:
                                    pass
                                try:
                                    ms_blk = vcf_to_ms_like(vb_adj, neut_meta.get('length'), ploidy=neut_meta.get('ploidy', 1))
                                except Exception:
                                    ms_blk = 'segsites: 0\n'
                                if not ms_blk.endswith('\n'):
                                    ms_blk += '\n'
                                out_ms_chunks.append('//\n' + ms_blk)
                            ms_like_for_sfs = ''.join(out_ms_chunks)
                    except Exception:
                        pass
                    if args.sfs_mode == 'mean':
                        # We don't know n_hap cheaply here; let it infer
                        header_nsfs, lines_nsfs = compute_sfs_stream_mean(ms_like_for_sfs, n_hap_expected=None, normalized=args.sfs_normalized)
                    else:
                        header_nsfs, lines_nsfs = compute_sfs_fast(ms_like_for_sfs, normalized=args.sfs_normalized, mode=args.sfs_mode)
                    neut_sfs_path = neut_root + '.sfs'
                    nsfs_text = header_nsfs + '\n' + '\n'.join(lines_nsfs) + '\n'
                    with open(neut_sfs_path, 'w') as fns:
                        fns.write('# SFS output (neutral paired)\n')
                        fns.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                        fns.write(nsfs_text)
                    sys.stderr.write(f"# INFO: wrote neutral SFS to {neut_sfs_path}\n")
                except Exception as e:
                    sys.stderr.write(f"# ERROR: failed to compute/write neutral SFS: {e}\n")
                neut_out_path = None  # prevent file write path block below
            elif (neut_out_path is None) and getattr(args, 'sfs', False) and not getattr(args, '_user_out_provided', False):
                # No explicit neutral path given, SFS-only mode: write neutral SFS to default stem with _neutral suffix
                try:
                    if args.sfs_mode == 'mean':
                        header_nsfs, lines_nsfs = compute_sfs_stream_mean(neut_out_txt, n_hap_expected=None, normalized=args.sfs_normalized)
                    else:
                        header_nsfs, lines_nsfs = compute_sfs_fast(neut_out_txt, normalized=args.sfs_normalized, mode=args.sfs_mode)
                    neut_sfs_path = f"{default_stem}_neutral.sfs"
                    nsfs_text = header_nsfs + '\n' + '\n'.join(lines_nsfs) + '\n'
                    with open(neut_sfs_path, 'w') as fns:
                        fns.write('# SFS output (neutral paired)\n')
                        fns.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                        fns.write(nsfs_text)
                    sys.stderr.write(f"# INFO: wrote neutral SFS to {neut_sfs_path}\n")
                except Exception as e:
                    sys.stderr.write(f"# ERROR: failed to compute/write neutral SFS: {e}\n")
                neut_out_path = None
            if neut_out_path and neut_out_path != '-':
                try:
                    if fmt.startswith('vcf') or fmt == 'bcf':
                        # Helper to build neutral VCF meta lines once
                        def _build_neut_vcf_meta_lines():
                            info_pairs_n = []
                            info_pairs_n.append(('engine', neut_meta.get('engine', neut_engine)))
                            info_pairs_n.append(('species', neut_meta.get('species_id')))
                            info_pairs_n.append(('model', args.model))
                            info_pairs_n.append(('pop0', neut_meta.get('pop0')))
                            pop_order_n = neut_meta.get('pop_order') or []
                            if pop_order_n:
                                info_pairs_n.append(('order', '|'.join(map(str, pop_order_n))))
                            counts_disc_n = neut_meta.get('counts_disc') or []
                            if counts_disc_n:
                                pl_n = neut_meta.get('ploidy', 1) or 1
                                try:
                                    sample_list_n = [int(h // pl_n) for h in counts_disc_n]
                                except Exception:
                                    sample_list_n = counts_disc_n
                                info_pairs_n.append(('sample_counts', '|'.join(map(str, sample_list_n))))
                                info_pairs_n.append(('hap_counts', '|'.join(map(str, counts_disc_n))))
                            theta_val_n = rho_val_n = None
                            try:
                                toks_nt = raw_neut.split()
                                if '-t' in toks_nt:
                                    theta_val_n = toks_nt[toks_nt.index('-t')+1]
                                if '-r' in toks_nt:
                                    rho_val_n = toks_nt[toks_nt.index('-r')+1]
                            except Exception:
                                pass
                            if theta_val_n: info_pairs_n.append(('theta', theta_val_n))
                            if rho_val_n: info_pairs_n.append(('rho', rho_val_n))
                            info_pairs_n.append(('length', neut_meta.get('length')))
                            info_pairs_n.append(('N0', neut_meta.get('N0')))
                            info_pairs_n.append(('mu', neut_meta.get('mu')))
                            info_pairs_n.append(('r', neut_meta.get('rrate')))
                            info_pairs_n.append(('ploidy', neut_meta.get('ploidy')))
                            if neut_meta.get('chromosome') is not None:
                                info_pairs_n.append(('chromosome', neut_meta.get('chromosome')))
                            kvn = {k:v for k,v in info_pairs_n if v is not None}
                            core_keys_n = ['engine','species','model','pop0','order','sample_counts','hap_counts']
                            param_keys_n = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                            def build_line_n_vcf(group_name, keys):
                                vals = [f"{k}={kvn[k]}" for k in keys if k in kvn]
                                if not vals:
                                    return None
                                return '# ' + group_name + ': ' + ', '.join(vals)
                            meta_lines_neut_local = []
                            for grp_n, keys_n in (('core', core_keys_n), ('params', param_keys_n)):
                                ln_n = build_line_n_vcf(grp_n, keys_n)
                                if ln_n:
                                    meta_lines_neut_local.append(ln_n)
                            return meta_lines_neut_local

                        per_rep_vcfs_neut = []
                        if neut_meta.get('engine') == 'msprime' and neut_meta.get('produced_format') == 'vcf' and neut_out_txt:
                            vtxt_n = neut_out_txt
                            # Split concatenated VCFs on marker
                            blocks_n = []
                            if '##fileformat' in vtxt_n:
                                parts_n = vtxt_n.split('##fileformat')
                                for p in parts_n:
                                    p = p.strip()
                                    if not p:
                                        continue
                                    blocks_n.append('##fileformat' + p if not p.startswith('##fileformat') else p)
                            else:
                                blocks_n = [vtxt_n]
                            for vb in blocks_n:
                                # Skip header-only VCF blocks (no variant lines). A
                                # header-only block will only contain lines starting
                                # with '#'. We require at least one non-# line.
                                try:
                                    non_hash = any((ln.strip() and not ln.startswith('#')) for ln in vb.splitlines())
                                except Exception:
                                    non_hash = False
                                if not non_hash:
                                    continue
                                try:
                                    vb = _adjust_vcf_contig_length(vb, neut_meta.get('chromosome'), neut_meta.get('length'))
                                except Exception:
                                    pass
                                vcf_lines_neut = vb.splitlines()
                                meta_lines_neut = _build_neut_vcf_meta_lines()
                                try:
                                    header_idx_neut = next(i for i,l in enumerate(vcf_lines_neut) if l.startswith('#CHROM'))
                                except StopIteration:
                                    header_idx_neut = len(vcf_lines_neut)
                                vcf_meta_lines_neut = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines_neut]
                                vcf_lines_neut = vcf_lines_neut[:header_idx_neut] + vcf_meta_lines_neut + vcf_lines_neut[header_idx_neut:]
                                per_rep_vcfs_neut.append('\n'.join(vcf_lines_neut) + ('\n' if vcf_lines_neut and not vcf_lines_neut[-1].endswith('\n') else ''))
                        else:
                            # Convert ms-like neutral stdout to VCF per replicate
                            text_n = neut_out_txt or ''
                            blocks_n = []
                            cur_n = []
                            for ln in (text_n.splitlines(True)):
                                if ln.lstrip().startswith('//'):
                                    if cur_n:
                                        blocks_n.append(''.join(cur_n))
                                        cur_n = []
                                    continue
                                if ln.lstrip().startswith('ms ') or ln.lstrip().startswith('msms ') or ln.lstrip().startswith('discoal ') or ln.lstrip().startswith('scrm '):
                                    continue
                                cur_n.append(ln)
                            if cur_n:
                                blocks_n.append(''.join(cur_n))
                            seg_hdr_pat_n = re.compile(r'^\s*segsites\s*:', re.IGNORECASE | re.MULTILINE)
                            for blk in blocks_n:
                                # Skip blocks that don't contain any segsites header
                                if not seg_hdr_pat_n.search(blk or ''):
                                    continue
                                vtxt = ms_like_to_vcf(blk, neut_meta['length'], chrom=neut_meta.get('chromosome') or 'chr1', ploidy=neut_meta.get('ploidy',1))
                                vcf_lines = vtxt.splitlines()
                                meta_lines_neut = _build_neut_vcf_meta_lines()
                                try:
                                    header_idx_neut = next(i for i,l in enumerate(vcf_lines) if l.startswith('#CHROM'))
                                except StopIteration:
                                    header_idx_neut = len(vcf_lines)
                                vcf_meta_lines_neut = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines_neut]
                                vcf_lines_final = vcf_lines[:header_idx_neut] + vcf_meta_lines_neut + vcf_lines[header_idx_neut:]
                                per_rep_vcfs_neut.append('\n'.join(vcf_lines_final) + ('\n' if vcf_lines_final and not vcf_lines_final[-1].endswith('\n') else ''))

                        # Write outputs: per-replicate into a folder when reps > 1, else single file
                        # Optional debug: dump per-rep vcf summaries
                        try:
                            if os.getenv('RAISD_DEBUG_NEUT'):
                                dbg_vpath = os.path.join(os.getcwd(), 'neutral_vcf_blocks.log')
                                with open(dbg_vpath, 'w') as dv:
                                    dv.write(f"per_rep_vcfs_neut count: {len(per_rep_vcfs_neut)}\n")
                                    for ii, vv in enumerate(per_rep_vcfs_neut):
                                        dv.write(f"-- block {ii} len={len(vv)}\n")
                                        dv.write((vv[:500].replace('\r','') + ('\n...[truncated]\n' if len(vv) > 500 else '\n')))
                        except Exception:
                            pass
                        if args.reps and args.reps > 1:
                            parent_n = os.path.dirname(neut_out_path)
                            base_n = os.path.basename(neut_out_path)
                            # derive folder name by stripping known extensions
                            folder_base_n = base_n
                            if fmt == 'vcf.gz' and base_n.endswith('.vcf.gz'):
                                folder_base_n = base_n[:-7]
                            elif fmt == 'vcf' and base_n.endswith('.vcf'):
                                folder_base_n = base_n[:-4]
                            elif fmt == 'bcf' and base_n.endswith('.bcf'):
                                folder_base_n = base_n[:-4]
                            # ensure folder base includes _neutral suffix only when
                            # the neutral name was auto-derived (not user-specified).
                            if not user_specified_neut:
                                if not folder_base_n.endswith('_neutral'):
                                    folder_base_n = folder_base_n + '_neutral'
                            out_dir_n = os.path.join(parent_n, folder_base_n) if parent_n else folder_base_n
                            try:
                                os.makedirs(out_dir_n, exist_ok=True)
                            except Exception:
                                pass
                            root_n = folder_base_n
                            for i, vtxt in enumerate(per_rep_vcfs_neut, start=1):
                                if fmt == 'vcf.gz':
                                    path_out = os.path.join(out_dir_n, f"{root_n}_{i}.vcf.gz")
                                    vtxt_adj = _ensure_vcf_has_samples(vtxt, neut_meta)
                                    _write_bgzip_text(path_out, vtxt_adj)
                                elif fmt == 'bcf':
                                    path_out = os.path.join(out_dir_n, f"{root_n}_{i}.bcf")
                                    tmp_vcf = None
                                    try:
                                        tmpdir = _ensure_local_tmpdir()
                                        # Ensure header includes sample columns before conversion
                                        vtxt_adj = _ensure_vcf_has_samples(vtxt, neut_meta)
                                        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.vcf', dir=tmpdir) as tf:
                                            tf.write(vtxt_adj)
                                            tmp_vcf = tf.name
                                        _vcf_text_to_bcf(tmp_vcf, path_out)
                                    finally:
                                        try:
                                            if tmp_vcf is not None and os.path.exists(tmp_vcf):
                                                os.unlink(tmp_vcf)
                                        except Exception:
                                            pass
                                else:
                                    path_out = os.path.join(out_dir_n, f"{root_n}_{i}.vcf")
                                    vtxt_adj = _ensure_vcf_has_samples(vtxt, neut_meta)
                                    with open(path_out, 'w') as f2:
                                        f2.write(vtxt_adj)
                        else:
                            # Single file
                            # If multiple per-rep VCF blocks were produced, prefer
                            # the one that contains variants (heuristic: the
                            # longest block). This prevents writing a header-only
                            # VCF when another block holds the variant lines.
                            if per_rep_vcfs_neut:
                                try:
                                    single_vtxt = max(per_rep_vcfs_neut, key=len)
                                except Exception:
                                    single_vtxt = per_rep_vcfs_neut[0]
                            else:
                                single_vtxt = ''
                            if fmt == 'vcf.gz':
                                single_vtxt_adj = _ensure_vcf_has_samples(single_vtxt, neut_meta)
                                _write_bgzip_text(neut_out_path, single_vtxt_adj)  # type: ignore[arg-type]
                            elif fmt == 'bcf':
                                tmp_vcf = None
                                try:
                                    tmpdir = _ensure_local_tmpdir()
                                    # Ensure header includes sample columns before conversion
                                    single_vtxt_adj = _ensure_vcf_has_samples(single_vtxt, neut_meta)
                                    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.vcf', dir=tmpdir) as tf:
                                        tf.write(single_vtxt_adj)
                                        tmp_vcf = tf.name
                                    _vcf_text_to_bcf(tmp_vcf, neut_out_path)
                                finally:
                                    try:
                                        if tmp_vcf is not None and os.path.exists(tmp_vcf):
                                            os.unlink(tmp_vcf)
                                    except Exception:
                                        pass
                            else:
                                single_vtxt_adj = _ensure_vcf_has_samples(single_vtxt, neut_meta)
                                with open(neut_out_path, 'w') as f2:
                                    f2.write(single_vtxt_adj)
                    else:
                        # Raw neutral ms-like output (optionally gzipped)
                        f_open_neut = gzip.open if fmt == 'ms.gz' else open
                        mode_neut = 'wt' if fmt == 'ms.gz' else 'w'
                        with f_open_neut(neut_out_path, mode_neut) as f2:
                            f2.write(raw_neut + '\n')  # type: ignore[arg-type]
                            info_pairs_n = []
                            info_pairs_n.append(('engine', neut_meta.get('engine', neut_engine)))
                            info_pairs_n.append(('species', neut_meta.get('species_id')))
                            info_pairs_n.append(('model', args.model))
                            info_pairs_n.append(('pop0', neut_meta.get('pop0')))
                            pop_order_n = neut_meta.get('pop_order') or []
                            if pop_order_n:
                                info_pairs_n.append(('order', '|'.join(map(str, pop_order_n))))
                            counts_disc_n = neut_meta.get('counts_disc') or []
                            if counts_disc_n:
                                pl_n = neut_meta.get('ploidy', 1) or 1
                                try:
                                    sample_list_n = [int(h // pl_n) for h in counts_disc_n]
                                except Exception:
                                    sample_list_n = counts_disc_n
                                info_pairs_n.append(('sample_counts', '|'.join(map(str, sample_list_n))))
                                info_pairs_n.append(('hap_counts', '|'.join(map(str, counts_disc_n))))
                            theta_val_n = rho_val_n = None
                            try:
                                toks_nt = raw_neut.split()
                                if '-t' in toks_nt:
                                    theta_val_n = toks_nt[toks_nt.index('-t')+1]
                                if '-r' in toks_nt:
                                    rho_val_n = toks_nt[toks_nt.index('-r')+1]
                            except Exception:
                                pass
                            if theta_val_n: info_pairs_n.append(('theta', theta_val_n))
                            if rho_val_n: info_pairs_n.append(('rho', rho_val_n))
                            info_pairs_n.append(('length', neut_meta.get('length')))
                            info_pairs_n.append(('N0', neut_meta.get('N0')))
                            info_pairs_n.append(('mu', neut_meta.get('mu')))
                            info_pairs_n.append(('r', neut_meta.get('rrate')))
                            info_pairs_n.append(('ploidy', neut_meta.get('ploidy')))
                            if neut_meta.get('chromosome') is not None:
                                info_pairs_n.append(('chromosome', neut_meta.get('chromosome')))
                            kvn = {k:v for k,v in info_pairs_n if v is not None}
                            core_keys_n = ['engine','species','model','pop0','order','sample_counts','hap_counts']
                            param_keys_n = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                            def build_line_n(group_name, keys):
                                vals = [f"{k}={kvn[k]}" for k in keys if k in kvn]
                                if not vals: return None
                                return '# ' + group_name + ': ' + ', '.join(vals)
                            for grp_n, keys_n in (('core', core_keys_n), ('params', param_keys_n)):
                                ln_n = build_line_n(grp_n, keys_n)
                                if ln_n:
                                    f2.write(ln_n + '\n')  # type: ignore[arg-type]
                            # Body: if msprime produced VCF but requested ms/ms.gz, convert VCF->ms-like per replicate
                            if neut_meta.get('engine') == 'msprime' and neut_meta.get('produced_format') == 'vcf' and neut_out_txt:
                                vcf_text_all_n = neut_out_txt
                                # Build per-rep seed list from seeds_for_chunks_n and neut_chunks
                                seed_lines_n = []
                                try:
                                    if isinstance(neut_chunks, list) and isinstance(seeds_for_chunks_n, list):
                                        for ci, base_seed in enumerate(seeds_for_chunks_n):
                                            if base_seed is None:
                                                continue
                                            reps_in_chunk = 1 if ci >= len(neut_chunks) else int(neut_chunks[ci])
                                            for r in range(reps_in_chunk):
                                                seed_lines_n.append(str(int(base_seed) + int(r)))
                                except Exception:
                                    seed_lines_n = []
                                # Split concatenated VCFs on '##fileformat'
                                blocks_n = []
                                if '##fileformat' in vcf_text_all_n:
                                    parts_n = vcf_text_all_n.split('##fileformat')
                                    for p in parts_n:
                                        p = p.strip()
                                        if not p:
                                            continue
                                        blocks_n.append('##fileformat' + p if not p.startswith('##fileformat') else p)
                                else:
                                    blocks_n = [vcf_text_all_n]
                                for ib, vb in enumerate(blocks_n):
                                    vb_adj = vb
                                    try:
                                        vb_adj = _adjust_vcf_contig_length(vb_adj, neut_meta.get('chromosome'), neut_meta.get('length'))
                                    except Exception:
                                        pass
                                    # per-rep seed line
                                    if ib < len(seed_lines_n):
                                        f2.write(seed_lines_n[ib] + '\n')  # type: ignore[arg-type]
                                    f2.write('//\n')  # replicate separator
                                    try:
                                        ms_blk_n = vcf_to_ms_like(vb_adj, neut_meta.get('length'), ploidy=neut_meta.get('ploidy', 1))
                                    except Exception:
                                        ms_blk_n = 'segsites: 0\n'
                                    if not ms_blk_n.endswith('\n'):
                                        ms_blk_n += '\n'
                                    f2.write(ms_blk_n)  # type: ignore[arg-type]
                            else:
                                for ln2 in neut_out_txt.splitlines(True):
                                    if ln2.startswith(neut_engine+' '):
                                        continue
                                    if neut_engine == 'msms' and ln2.startswith('ms '):
                                        continue
                                    f2.write(ln2)  # type: ignore[arg-type]
                            sys.stderr.write(f"# INFO: wrote paired neutral output to {neut_out_path}\n")
                    if getattr(args, 'sfs', False):
                        try:
                            # If msprime produced VCF, convert to ms-like per replicate before SFS
                            ms_like_for_sfs = neut_out_txt
                            try:
                                if neut_meta.get('engine') == 'msprime' and neut_meta.get('produced_format') == 'vcf' and neut_out_txt:
                                    vtxt = neut_out_txt
                                    blocksC = []
                                    if '##fileformat' in vtxt:
                                        partsC = vtxt.split('##fileformat')
                                        for p in partsC:
                                            p = p.strip()
                                            if not p:
                                                continue
                                            blocksC.append('##fileformat' + p if not p.startswith('##fileformat') else p)
                                    else:
                                        blocksC = [vtxt]
                                    out_ms_chunks = []
                                    for vb in blocksC:
                                        vb_adj = vb
                                        try:
                                            vb_adj = _adjust_vcf_contig_length(vb_adj, neut_meta.get('chromosome'), neut_meta.get('length'))
                                        except Exception:
                                            pass
                                        try:
                                            ms_blk = vcf_to_ms_like(vb_adj, neut_meta.get('length'), ploidy=neut_meta.get('ploidy', 1))
                                        except Exception:
                                            ms_blk = 'segsites: 0\n'
                                        if not ms_blk.endswith('\n'):
                                            ms_blk += '\n'
                                        out_ms_chunks.append('//\n' + ms_blk)
                                    ms_like_for_sfs = ''.join(out_ms_chunks)
                            except Exception:
                                pass
                            header_nsfs, lines_nsfs = compute_sfs_fast(ms_like_for_sfs, normalized=args.sfs_normalized, mode=args.sfs_mode)
                            neut_root = neut_out_path
                            if neut_root.endswith('.vcf.gz'):
                                neut_root = neut_root[:-7]
                            elif neut_root.endswith('.ms.gz'):
                                neut_root = neut_root[:-6]
                            else:
                                neut_root = os.path.splitext(neut_root)[0]
                            neut_sfs_path = neut_root + '.sfs'
                            nsfs_text = header_nsfs + '\n' + '\n'.join(lines_nsfs) + '\n'
                            with open(neut_sfs_path, 'w') as fns:
                                fns.write('# SFS output (neutral paired)\n')
                                fns.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                                fns.write(nsfs_text)
                            sys.stderr.write(f"# INFO: wrote neutral SFS to {neut_sfs_path}\n")
                        except Exception as e:
                            sys.stderr.write(f"# ERROR: failed to compute/write neutral SFS: {e}\n")
                except Exception as e:
                    sys.stderr.write(f"# ERROR: failed to write neutral output {neut_out_path}: {e}\n")
            if neut_err_txt:
                sys.stderr.write(neut_err_txt)
        # SFS already written earlier under new naming; do not duplicate
        # If we deferred SFS writing until the end, perform the actual write now.
        try:
            if 'sfs_pending' in locals() and sfs_pending:
                sfs_path_final, sfs_text_final = sfs_pending
                try:
                    with open(sfs_path_final, 'w') as fsp:
                        fsp.write('# SFS output\n')
                        fsp.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                        fsp.write(sfs_text_final)
                        try:
                            fsp.flush()
                            os.fsync(fsp.fileno())
                        except Exception:
                            pass
                    if not sfs_written:
                        sys.stderr.write(f"# INFO: wrote SFS to {sfs_path_final}\n")
                        sfs_written = True
                except Exception as e:
                    sys.stderr.write(f"# ERROR: failed to write deferred SFS file {sfs_path_final}: {e}\n")
                    # Do not print the full SFS text to stdout per user preference.
                    # Previously the code would fall back to writing the SFS to stdout
                    # when no output file was specified; that behavior is suppressed.
                    if not args.out or args.out == '-':
                        sys.stderr.write("# INFO: deferred SFS could not be written to file; not printing SFS to stdout by user preference\n")
        except Exception:
            pass

def main():
    args = parse_args()
    return run_simulation_from_args(args)


if __name__ == "__main__":
    try:
        # Ensure main() return value is used as exit code when appropriate
        rc = main()
        if rc is None:
            sys.exit(0)
        else:
            # If main returns an int-like value, use it as exit code
            try:
                sys.exit(int(rc))
            except Exception:
                sys.exit(0)
    except ExceedsMaxRam as e:
        # Friendly single-line message for RAM aborts (no traceback)
        try:
            msg = str(e)
        except Exception:
            msg = "Aborted due to exceeding configured maximum RAM"
        sys.stderr.write(msg.rstrip() + "\n")
        # exit with a distinct non-zero code
        sys.exit(2)
    except MemoryError as e:
        # Out-of-memory inside Python: provide actionable suggestions without traceback
        try:
            msg = ("Aborted due to MemoryError (process ran out of memory). "
                   "Suggested actions: reduce --parallel; lower --sims-per-work; "
                   "reduce --replicates or --sample-individuals; or increase --max-ram-percent / add swap.")
        except Exception:
            msg = "Aborted due to MemoryError (out of memory)."
        sys.stderr.write(msg.rstrip() + "\n")
        sys.exit(2)
    except OSError as e:
        # Some low-level allocations raise OSError(errno.ENOMEM) when the OS cannot satisfy memory requests
        try:
            if getattr(e, 'errno', None) == errno.ENOMEM:
                msg = (f"Aborted due to OS OOM (ENOMEM): {e}. "
                       "Suggested actions: reduce parallelism or per-worker memory; "
                       "increase system memory or add swap; or raise --max-ram-percent.")
                sys.stderr.write(msg.rstrip() + "\n")
                sys.exit(2)
        except Exception:
            pass
        # Re-raise unexpected OSErrors
        raise