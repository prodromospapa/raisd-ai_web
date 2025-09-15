import os
import multiprocessing as mp
import argparse
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import stdpopsim
import subprocess
import time
import sys
import logging
import signal
import psutil

# CLI: species (stdpopsim id or full name) plus optional resource flags
parser = argparse.ArgumentParser()
parser.add_argument("--species", type=str, required=True,
                    help="Species id (stdpopsim short id) or full name, e.g. 'HomSap' or 'Homo sapiens'.")
parser.add_argument("--max-ram-percent", type=float, default=None,
                    help="If set, stop a simulator run when it exceeds this percent of total RAM; the run will be skipped (no lower-parallel retry). If omitted, wrapper-level monitoring is disabled (the simulator may still enforce its own threshold).")
parser.add_argument("--max-parallel", type=int, default=None,
                    help="Maximum --parallel value to start with. If omitted, defaults to min(half of CPU threads, replicates).")
parser.add_argument("--model-id", type=str, default=None,
                    help="(Optional) If set, only run this demographic model id.")
parser.add_argument("--population", type=str, default=None,
                    help="(Optional) If set with --model-id, only run this population within the model.")
parser.add_argument("--samples", type=int, default=10000,
                    help="Number of haploid samples per population to request when building the SFS (default: 10000).")
parser.add_argument("--replicates", type=int, default=5,
                    help="Number of replicates to request from the simulator for each model/population when building SFS (default: 5).")
parser.add_argument("--engine", type=str, default="msprime",
                    help="Default simulator engine to request (default: msprime).")
parser.add_argument("--chromosome", type=str, default=None,
                    help="Chromosome id to simulate; if omitted the wrapper will pick the numerically largest chromosome for the species.")
parser.add_argument("--sims-per-work", type=int, default=1,
                    help="Number of sims per worker (default: 1).")
# sfs filename and normalization are fixed for SFS generation (not CLI options)
cli_args, _ = parser.parse_known_args()
species = cli_args.species
MAX_RAM_PERCENT = cli_args.max_ram_percent
# default sample count (per-population copies to simulate when building SFS)
samples = cli_args.samples
REPLICATES = cli_args.replicates
engine_default = cli_args.engine
# chromosome will be resolved after the species object is available; if omitted, use the numerically largest chromosome
# (we set `chromosome` below once `species_std` is known).

current_sims_per_work = cli_args.sims_per_work
# sfs filename and normalization are fixed for SFS generation
# Request normalized SFS by default and use the standard filename.
sfs_filename = "sfs.sfs"
sfs_normalized_flag = True
# We'll compute a sensible default for MAX_PARALLEL: at most half the available CPUs, but no more than replicates.
cpu_count = mp.cpu_count() or 1
default_parallel = max(1, cpu_count // 2)
if cli_args.max_parallel is None:
    MAX_PARALLEL = min(default_parallel, REPLICATES)
else:
    MAX_PARALLEL = cli_args.max_parallel

# basic logger (not used for file logging here; we use append_log instead)
logger = logging.getLogger("raisd")
logger.setLevel(logging.INFO)

def biggest_chrom(species_std):
    # return the ID (e.g., "1") of the numerically largest chromosome by length (no deprecation)
    chrom = [c for c in species_std.genome.chromosomes]
    return max(chrom, key=lambda c: int(c.length)).id

def read_sfs_file(fn):
    with open(fn) as f:
        for i, l in enumerate(f):
            if i == 3:
                return np.fromstring(' '.join(l.split()[1:]), sep=' ')

species_dict = {sp.name: sp.id for sp in stdpopsim.all_species()}
if species in species_dict.values():
    species_full_name = [sp.name for sp in stdpopsim.all_species() if sp.id == species][0]
else:
    species_full_name = species
    species = species_dict[species]

species_std = stdpopsim.get_species(species)
biggest_chromosome = biggest_chrom(species_std)

# If user did not pass --chromosome, choose the numerically largest chromosome for this species
chromosome = cli_args.chromosome if cli_args.chromosome is not None else biggest_chromosome

# Build model->pop list where sampling_time == 0
demographic_models = {
    m.id: [p.name for p in species_std.get_demographic_model(m.id).populations]
    for m in species_std.demographic_models
}

# If user requested a single model/population for testing, filter accordingly
if cli_args.model_id is not None:
    if cli_args.model_id not in demographic_models:
        raise SystemExit(f"Requested model-id {cli_args.model_id} not found for species {species_full_name}")
    pops = demographic_models[cli_args.model_id]
    if cli_args.population is not None:
        if cli_args.population not in pops:
            raise SystemExit(f"Requested population {cli_args.population} not found in model {cli_args.model_id}")
        demographic_models = {cli_args.model_id: [cli_args.population]}
    else:
        demographic_models = {cli_args.model_id: pops}

species_folder_name = "".join(c if c.isalnum() or c in "-_ " else "_" for c in str(species_full_name)).strip().replace(" ", "_")

ploidy = species_std.ploidy


os.system(f'mkdir -p "data/{species_folder_name}"')

# shared record of failed parts (JSON-lines)
skipped_file = f"data/{species_folder_name}/failed_parts.jsonl"

import json
import atexit
import tempfile as _tempfile
import shutil as _shutil
import signal as _signal

def load_skipped():
    """Load skipped keys from the JSONL file.

    This function tolerates JSON objects concatenated without newlines by
    using json.JSONDecoder.raw_decode to extract objects sequentially. It
    keeps only the first-seen entry per key and atomically rewrites the
    file if duplicates or malformed concatenations were present.
    """
    skipped = set()
    try:
        if os.path.exists(skipped_file):
            text = ''
            try:
                with open(skipped_file, 'r', encoding='utf-8') as rf:
                    text = rf.read()
            except Exception:
                text = ''

            decoder = json.JSONDecoder()
            pos = 0
            first_by_key = {}
            L = len(text)
            while pos < L:
                try:
                    obj, offset = decoder.raw_decode(text, pos)
                    pos = pos + offset
                    key = obj.get('key') if isinstance(obj, dict) else None
                    if key:
                        k = str(key).strip()
                        if k not in first_by_key:
                            # store the canonical JSON string for rewrite
                            try:
                                first_by_key[k] = json.dumps(obj, ensure_ascii=False)
                            except Exception:
                                first_by_key[k] = str(obj)
                except ValueError:
                    # skip one character and continue (handles stray separators/newlines)
                    pos += 1
                    continue

            for k in first_by_key.keys():
                skipped.add(k)

            # If number of extracted objects differs from file length by lines,
            # always rewrite to keep only the first entry per key.
            try:
                # rewrite atomically
                d = os.path.dirname(skipped_file) or '.'
                fd, tmp = _tempfile.mkstemp(prefix='.failed_parts_tmp_', suffix='.jsonl', dir=d)
                try:
                    with os.fdopen(fd, 'w', encoding='utf-8') as wf:
                        for line in first_by_key.values():
                            wf.write(line + '\n')
                        try:
                            wf.flush()
                            os.fsync(wf.fileno())
                        except Exception:
                            pass
                    os.replace(tmp, skipped_file)
                except Exception:
                    try:
                        if os.path.exists(tmp):
                            os.remove(tmp)
                    except Exception:
                        pass
            except Exception:
                # best-effort only
                pass
    except Exception:
        pass
    return skipped


def safe_write_sfs(path=f"data/{species_folder_name}/sfs.csv"):
    """Atomically write the current `sfs` DataFrame to `path`.

    Writes to a temporary file in the same directory then os.replace() to
    ensure the write is atomic on POSIX filesystems. This avoids corrupting
    or truncating the existing CSV if the process is interrupted.
    """
    try:
        d = os.path.dirname(path)
        if d and not os.path.exists(d):
            os.makedirs(d, exist_ok=True)
        fd, tmp = _tempfile.mkstemp(prefix='.sfs_tmp_', suffix='.csv', dir=d or '.')
        os.close(fd)
        # Record the temp path so the signal/exit handler can clean it up if needed
        try:
            _sfs_temp_paths.add(tmp)
        except Exception:
            pass
        # Use pandas to write; follow by atomic replace
        try:
            sfs.to_csv(tmp)
            os.replace(tmp, path)
            # successful replace; remove from the tracked set
            try:
                _sfs_temp_paths.discard(tmp)
            except Exception:
                pass
        except Exception:
            # best-effort cleanup
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass
            try:
                _sfs_temp_paths.discard(tmp)
            except Exception:
                pass
    except Exception:
        try:
            logger.error(f"Failed to atomically write sfs to {path}", exc_info=True)
        except Exception:
            pass


def _on_exit_or_signal(signum=None, frame=None):
    # Attempt to persist SFS before exiting; swallow any errors to avoid
    # masking original exceptions.
    try:
        safe_write_sfs()
    except Exception:
        pass
    # Clean up any leftover temporary sfs files created earlier
    try:
        for tmp in list(_sfs_temp_paths):
            try:
                if tmp and os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass
            try:
                _sfs_temp_paths.discard(tmp)
            except Exception:
                pass
    except Exception:
        pass
    # If this was a signal handler, re-raise KeyboardInterrupt for SIGINT
    # so default behavior (exit) occurs after we've flushed.
    try:
        if signum == _signal.SIGINT:
            raise KeyboardInterrupt()
    except Exception:
        pass


# Register handlers to ensure SFS is flushed on normal exit and on common signals
try:
    atexit.register(_on_exit_or_signal)
except Exception:
    pass
try:
    for sig in (_signal.SIGINT, _signal.SIGTERM, _signal.SIGHUP):
        try:
            _signal.signal(sig, _on_exit_or_signal)
        except Exception:
            pass
except Exception:
    pass

def append_skipped(key, reason='skipped', source='3.sfs'):
    # Normalize key early
    try:
        key_norm = str(key).strip()
    except Exception:
        key_norm = key

    # If we already recorded this key (loaded at startup or appended earlier),
    # do not write another entry. This enforces one reason per demographic.
    try:
        if key_norm in existing_skipped:
            return key_norm
    except Exception:
        pass

    entry = {
        'key': key_norm,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'reason': reason,
        'source': source,
    }

    # Attempt to append to the canonical JSONL file atomically (append + fsync).
    try:
        with open(skipped_file, 'a', encoding='utf-8') as wf:
            wf.write(json.dumps(entry, ensure_ascii=False) + "\n")
            try:
                wf.flush()
                os.fsync(wf.fileno())
            except Exception:
                pass
        # record in-memory to avoid further writes for this key
        try:
            existing_skipped.add(key_norm)
        except Exception:
            pass
        return key_norm
    except Exception:
        # fallback to writing a bare line in the old file location for maximum compatibility
        try:
            fallback = skipped_file.replace('.jsonl', '.txt')
            with open(fallback, 'a') as lf:
                lf.write(str(key_norm) + "\n")
                try:
                    lf.flush()
                    os.fsync(lf.fileno())
                except Exception:
                    pass
            try:
                existing_skipped.add(key_norm)
            except Exception:
                pass
            return key_norm
        except Exception:
            # Give up on writing fallback; report to stdout instead
            try:
                print(f"Failed to write failed parts to {skipped_file}", flush=True)
            except Exception:
                pass
            return None

existing_skipped = load_skipped()

# track temporary sfs files created during atomic writes so we can clean them on exit
try:
    _sfs_temp_paths = set()
except Exception:
    _sfs_temp_paths = set()

# helper to append logs both to species folder and current working dir
def append_log(msg: str):
    # Do not write any log files. Print to stdout and also log via the logger.
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    try:
        # Do not print warnings about skips to the terminal per user request.
        # Keep a minimal info-level log via the logger only.
        logger.info(msg)
    except Exception:
        # As a last resort, fall back to logger which may be configured by caller
        try:
            logger.info(msg)
        except Exception:
            pass

if os.path.exists(f"data/{species_folder_name}/sfs.csv"):
    sfs = pd.read_csv(f"data/{species_folder_name}/sfs.csv",index_col=0)
    # Normalize index to string values (strip whitespace) so membership checks
    # like `if key in sfs.index` work reliably even if the CSV had odd types.
    try:
        sfs.index = sfs.index.astype(str)
        sfs.index = sfs.index.map(lambda x: x.strip() if isinstance(x, str) else x)
    except Exception:
        try:
            sfs.index = sfs.index.map(lambda x: str(x).strip())
        except Exception:
            pass
    samples = (sfs.shape[1] // ploidy)+1
    print(f"Resuming with {samples} samples per population")
    # Maintain a normalized set of keys for robust membership checks
    try:
        sfs_index_set = set(str(x).strip() for x in sfs.index.tolist() if x is not None)
    except Exception:
        sfs_index_set = set()
else:
    sfs = pd.DataFrame(columns= range(1,ploidy*(samples)))
    sfs_index_set = set()

# previously removed warning log; no longer used
first = True
total = sum(len(pops) for pops in demographic_models.values())
with tqdm(total=total, desc="Simulations", unit="run") as pbar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            key = f"{model_id}={population}"
            if key in existing_skipped:
                # If this demographic/population is in the skipped list, make sure
                # we don't write skipped demographics to the SFS CSV.
                # Ensure any existing row for this key is removed so skipped
                # demographics never appear in `sfs.csv`.
                try:
                    sfs.drop(index=key, inplace=True, errors='ignore')
                except Exception:
                    # best-effort: if dropping fails, ignore and continue
                    try:
                        logger.debug("Failed to drop skipped key from sfs in memory", exc_info=True)
                    except Exception:
                        pass
                try:
                    safe_write_sfs(f"data/{species_folder_name}/sfs.csv")
                except Exception:
                    # best-effort: if writing fails, continue but keep consistency in memory
                    logger.debug("Failed to persist sfs.csv for skipped key", exc_info=True)
                pbar.update(1)
                pbar.refresh()
                first = False
                append_log(f"Skipping {key} because it is listed in skipped_demographics.jsonl")
                try:
                    sfs_index_set.discard(key)
                except Exception:
                    pass
                continue

            if key in sfs_index_set or key in sfs.index:
                pbar.update(1)
                pbar.refresh()
                first = False
                continue
            if first:
                pass
            first = False

            # prepare base args; we'll modify --parallel dynamically on retries
            base_args = [
                 "simulator.py",
                "--engine", str(engine_default),
                "--species-id", str(species),
                "--model-id", str(model_id),
                "--pop-order", str(population),
                "--sample-individuals", str(samples),
                "--chromosome", str(chromosome),
                "--replicates", str(REPLICATES),
                "--parallel",
            ]
            # optionally include SFS arguments
            if sfs_filename:
                base_args += ["--sfs", str(sfs_filename)]
            if sfs_normalized_flag:
                base_args += ["--sfs-normalized"]
            # include sims-per-work
            base_args += ["--sims-per-work", str(current_sims_per_work)]
            # whether to include --sfs sfs.sfs in simulator args for this model/pop
            include_sfs = True
            current_parallel = MAX_PARALLEL
            # current_sims_per_work is set from CLI (default 1)
            skipped = False
            attempts = 0
            max_attempts = 6

            # We'll attempt runs, lowering parallel if RAM threshold is exceeded
            while attempts < max_attempts:
                attempts += 1
                args = list(base_args)

                # insert parallel value in the args list after the --parallel flag
                try:
                    pi = args.index("--parallel")
                    args.insert(pi + 1, str(current_parallel))
                except ValueError:
                    args += ["--parallel", str(current_parallel)]

                # ensure sims-per-work is set to 1 (keep simple: user requested default 1)
                # Ensure the args reflect the current desired --sims-per-work value.
                # If the flag already exists in base_args, replace its value; otherwise insert it after --parallel.
                try:
                    if "--sims-per-work" in args:
                        spi = args.index("--sims-per-work")
                        # replace existing value (or insert if missing)
                        if spi + 1 < len(args):
                            args[spi + 1] = str(current_sims_per_work)
                        else:
                            args.insert(spi + 1, str(current_sims_per_work))
                    else:
                        # insert after --parallel flag if present, otherwise append
                        try:
                            pi = args.index("--parallel")
                            args.insert(pi + 2, str(current_sims_per_work))
                            args.insert(pi + 1, "--sims-per-work")
                        except Exception:
                            args += ["--sims-per-work", str(current_sims_per_work)]
                except Exception:
                    # best-effort: fall back to appending
                    try:
                        args += ["--sims-per-work", str(current_sims_per_work)]
                    except Exception:
                        pass

                # always pass max-ram-percent flag to simulator so it enforces the same threshold
                if MAX_RAM_PERCENT is not None and "--max-ram-percent" not in args:
                    args += ["--max-ram-percent", str(MAX_RAM_PERCENT)]

                # ensure args include sfs if requested
                if include_sfs and "--sfs" not in args:
                    try:
                        pi = args.index("--parallel")
                        args.insert(pi+2, "sfs.sfs")
                        args.insert(pi+1, "--sfs")
                    except Exception:
                        args += ["--sfs", "sfs.sfs"]

                exceeded = False
                rc = 1
                proc = None
                import tempfile
                stderr_path = None
                stderr_fh = None
                try:
                    tf = tempfile.NamedTemporaryFile(prefix='sim_err_', suffix='.log', delete=False, dir=os.getcwd())
                    stderr_path = tf.name
                    stderr_fh = open(stderr_path, 'wb')
                except Exception:
                    stderr_fh = subprocess.DEVNULL
                    stderr_path = None
                try:
                    proc = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=stderr_fh)

                    # If wrapper-level monitoring requested, watch the child process's RSS
                    if MAX_RAM_PERCENT is not None and psutil is not None:
                        try:
                            p = psutil.Process(proc.pid)
                            total = psutil.virtual_memory().total
                            while proc.poll() is None:
                                try:
                                    rss = p.memory_info().rss
                                    perc = rss / total * 100.0
                                except Exception:
                                    # fallback to psutil's percent method
                                    try:
                                        perc = p.memory_percent()
                                    except Exception:
                                        perc = 0.0
                                if perc > MAX_RAM_PERCENT:
                                    append_log(f"Memory exceeded {MAX_RAM_PERCENT}% ({perc:.1f}%) for {model_id} {population} at parallel={current_parallel} (pid={proc.pid}).")
                                    exceeded = True
                                    # try to terminate the process tree
                                    try:
                                        p.kill()
                                    except Exception:
                                        try:
                                            proc.kill()
                                        except Exception:
                                            pass
                                    break
                                time.sleep(0.5)
                        except psutil.NoSuchProcess:
                            # process finished quickly
                            pass

                    # wait for process to exit if not already
                    if proc is not None:
                        try:
                            rc = proc.wait()
                        except Exception:
                            try:
                                proc.kill()
                            except Exception:
                                pass
                            rc = getattr(proc, 'returncode', 1)
                        finally:
                            # ensure stderr file is flushed/closed for inspection
                            try:
                                if stderr_fh is not None and stderr_fh is not subprocess.DEVNULL and not isinstance(stderr_fh, int) and hasattr(stderr_fh, 'flush'):
                                    try:
                                        stderr_fh.flush()
                                    except Exception:
                                        pass
                                if stderr_fh is not None and stderr_fh is not subprocess.DEVNULL and not isinstance(stderr_fh, int) and hasattr(stderr_fh, 'close'):
                                    try:
                                        stderr_fh.close()
                                    except Exception:
                                        pass
                            except Exception:
                                pass

                except Exception as e:
                    append_log(f"Exception while running simulator for {model_id} {population}: {e}")
                    exceeded = False
                    rc = getattr(proc, 'returncode', 1) if proc is not None else 1
                finally:
                    # If we redirected stderr to a temp file, try to read it for diagnostics
                    stderr_text = ''
                    try:
                        if stderr_path and os.path.isfile(stderr_path):
                            try:
                                with open(stderr_path, 'r', errors='replace') as ef:
                                    stderr_text = ef.read()
                            except Exception:
                                stderr_text = ''
                    except Exception:
                        stderr_text = ''
                    # If simulator printed a system-level RAM exceed message, treat accordingly
                    try:
                        if stderr_text and ('total system RAM usage' in stderr_text or 'Memory exceeded' in stderr_text or 'RAM usage' in stderr_text):
                            append_log(f"Simulator stderr indicates system RAM exceed for {model_id} {population} (parallel={current_parallel}).")
                            # Mark as exceeded so the retry/reduction logic below runs
                            exceeded = True
                            # If we're already at minimum parallel and sims-per-work, persist skip immediately
                            try:
                                # current_sims_per_work may be negative marker; normalize for check
                                csw = int(current_sims_per_work) if isinstance(current_sims_per_work, int) or (isinstance(current_sims_per_work, str) and current_sims_per_work.lstrip('-').isdigit()) else None
                            except Exception:
                                csw = None
                            if (current_parallel == 1) and (csw == 1 or csw is None):
                                try:
                                    k = append_skipped(key, reason=f'exceeded_system_ram_{MAX_RAM_PERCENT}', source='3.sfs')
                                    if k:
                                        existing_skipped.add(k)
                                    else:
                                        try:
                                            logger.error(f"Failed to append skip for {key} after detecting simulator stderr ram exceed")
                                        except Exception:
                                            pass
                                except Exception:
                                    try:
                                        logger.error(f"Failed to append skip for {key} after detecting simulator stderr ram exceed")
                                    except Exception:
                                        pass
                    except Exception:
                        # best-effort only
                        pass
                    # Clean up stderr temp file
                    try:
                        if stderr_path and os.path.isfile(stderr_path):
                            try:
                                os.remove(stderr_path)
                            except Exception:
                                pass
                    except Exception:
                        pass

                # Decide next action based on whether memory was exceeded or process exit code
                if exceeded:
                    # Progressive retry policy:
                    # 1) Try to reduce --sims-per-work (towards 1) while keeping parallel fixed.
                    #    We reduce by halving (faster) but always make progress.
                    # 2) When sims-per-work == 1, start lowering --parallel until 1.
                    # 3) When both are 1 and it still fails, persist skip.
                    try:
                        # normalize current_sims_per_work to integer if possible
                        try:
                            csw = int(current_sims_per_work) if (isinstance(current_sims_per_work, int) or (isinstance(current_sims_per_work, str) and current_sims_per_work.lstrip('-').isdigit())) else None
                        except Exception:
                            csw = None

                        # If we have a numeric sims-per-work > 1, reduce it (halve, but at least by 1)
                        if csw is not None and csw > 1:
                            new_csw = max(1, csw // 2)
                            if new_csw == csw:
                                new_csw = max(1, csw - 1)
                            current_sims_per_work = new_csw
                            append_log(f"Memory exceeded at parallel={current_parallel}; retrying with --sims-per-work={current_sims_per_work} (keeping parallel={current_parallel}).")
                            continue

                        # If sims-per-work is already 1 (or unknown/None), try lowering parallel
                        if current_parallel > 1:
                            current_parallel -= 1
                            current_sims_per_work = 1
                            append_log(f"Retrying with lower parallel: {current_parallel} (keeping --sims-per-work=1)")
                            continue

                        # If we've reached parallel==1 and sims-per-work==1, persist skip
                        append_log(f"Skipped {model_id} {population}: exceeded RAM {MAX_RAM_PERCENT}% even at parallel=1 and --sims-per-work=1.")
                        try:
                            k = append_skipped(key, reason=f'exceeded_ram_{MAX_RAM_PERCENT}', source='3.sfs')
                            if k:
                                existing_skipped.add(k)
                            else:
                                try:
                                    logger.error(f"Failed to append skip for {key} immediately")
                                except Exception:
                                    pass
                        except Exception:
                            try:
                                logger.error(f"Failed to append skip for {key} immediately")
                            except Exception:
                                pass
                        skipped = True
                        break
                    except Exception:
                        # best-effort fallback: if anything goes wrong, fall back to previous behavior
                        try:
                            neg_target = -int(current_parallel)
                        except Exception:
                            neg_target = -current_parallel
                        if str(current_sims_per_work) != str(neg_target):
                            current_sims_per_work = neg_target
                            append_log(f"Memory exceeded at parallel={current_parallel}; retrying with --sims-per-work={current_sims_per_work} (keeping parallel={current_parallel}).")
                            continue
                        if current_parallel > 1:
                            current_parallel -= 1
                            current_sims_per_work = 1
                            append_log(f"Retrying with lower parallel: {current_parallel} (keeping --sims-per-work=1)")
                            continue
                        skipped = True
                        try:
                            k = append_skipped(key, reason=f'exceeded_ram_{MAX_RAM_PERCENT}', source='3.sfs')
                            if k:
                                existing_skipped.add(k)
                        except Exception:
                            pass

                # process finished without exceeding memory
                if rc == 0:
                    # success
                    break
                else:
                    # non-zero exit: try swapping engine to scrm once (if not already), otherwise skip
                    try:
                        i = args.index("--engine")
                        if args[i+1] != "scrm":
                            args[i+1] = "scrm"
                            append_log(f"Warning: simulator failed (rc={rc}) for {model_id} {population} at parallel={current_parallel}. Retrying with engine=scrm.")
                            try:
                                subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                                rc = 0
                                break
                            except subprocess.CalledProcessError:
                                rc = 1
                        else:
                            append_log(f"Warning: simulator failed (rc={rc}) for {model_id} {population} at parallel={current_parallel}.")
                    except Exception:
                        append_log(f"Warning: simulator failed (rc={rc}) for {model_id} {population} at parallel={current_parallel}.")

                    # if we reach here, treat as failure and skip
                    # mark skipped and persist immediately so other runs won't retry
                    try:
                        k = append_skipped(key, reason='simulator_failed_nonzero_exit', source='3.sfs')
                        if k:
                            existing_skipped.add(k)
                        else:
                            try:
                                logger.error(f"Failed to append skip for {key} after non-zero exit")
                            except Exception:
                                pass
                    except Exception:
                        try:
                            logger.error(f"Failed to append skip for {key} after non-zero exit")
                        except Exception:
                            pass
                    skipped = True
                    break

            else:
                # exceeded max attempts
                append_log(f"Giving up after {attempts} attempts for {model_id} {population}; marking skipped")
                skipped = True
            # Safety: if we've decided to skip but no persistent record exists yet,
            # ensure we append it to the canonical skipped file. This guards against
            # code paths where an earlier append attempt may have failed.
            try:
                    if skipped and key not in existing_skipped:
                        try:
                            k = append_skipped(key, reason='failed_to_produce_sfs', source='3.sfs')
                            if k:
                                existing_skipped.add(k)
                            else:
                                try:
                                    logger.error(f"Failed to persist skip for {key} in safety net")
                                except Exception:
                                    pass
                        except Exception:
                            try:
                                logger.error(f"Failed to persist skip for {key} in safety net")
                            except Exception:
                                pass
            except Exception:
                pass

            if skipped:
                # Record the skip in the canonical JSONL skipped file and remove any existing data
                try:
                    k = append_skipped(key, reason='failed_to_produce_sfs', source='3.sfs')
                    if k:
                        existing_skipped.add(k)
                    else:
                        logger.error(f"Failed to append {key} to {skipped_file}")
                except Exception:
                    try:
                        logger.error(f"Failed to append {key} to {skipped_file}")
                    except Exception:
                        pass
                # remove any previously created data for this demographic/population (best-effort)
                try:
                    import shutil
                    target_dir = os.path.join("data", species_folder_name, str(model_id), str(population))
                    if os.path.exists(target_dir):
                        shutil.rmtree(target_dir)
                except Exception:
                    # keep going even if removal fails
                    try:
                        logger.debug("Failed to remove target_dir", exc_info=True)
                    except Exception:
                        logger.debug("Failed to remove target_dir", exc_info=True)
                # Ensure skipped demographics are not present in the SFS CSV in any form
                try:
                    sfs.drop(index=key, inplace=True, errors='ignore')
                except Exception:
                    try:
                        logger.debug("Failed to drop skipped key from sfs in memory after skip", exc_info=True)
                    except Exception:
                        pass
                try:
                    sfs_index_set.discard(key)
                except Exception:
                    pass
            else:
                # If no sfs.sfs was produced and we've already reduced to the
                # lowest resource settings (parallel=1 and sims-per-work=1),
                # treat this as a RAM/exhaustion failure and persist a skip.
                try:
                    csw = None
                    try:
                        if isinstance(current_sims_per_work, int):
                            csw = int(current_sims_per_work)
                        elif isinstance(current_sims_per_work, str) and current_sims_per_work.lstrip('-').isdigit():
                            csw = int(current_sims_per_work)
                    except Exception:
                        csw = None
                    if (not os.path.exists("sfs.sfs")) and (current_parallel == 1) and (csw == 1 or csw is None):
                        append_log(f"No sfs produced at parallel=1 and sims-per-work=1 for {model_id} {population}; marking skipped (assumed RAM/exhaustion).")
                        try:
                            k = append_skipped(key, reason=f'exceeded_ram_no_output_{MAX_RAM_PERCENT}', source='3.sfs')
                            if k:
                                existing_skipped.add(k)
                            else:
                                try:
                                    logger.error(f"Failed to append skip for {key} when no sfs was produced at minimal settings")
                                except Exception:
                                    pass
                        except Exception:
                            try:
                                logger.error(f"Failed to append skip for {key} when no sfs was produced at minimal settings")
                            except Exception:
                                pass
                        # remove any previously created data for this demographic/population (best-effort)
                        try:
                            import shutil
                            target_dir = os.path.join("data", species_folder_name, str(model_id), str(population))
                            if os.path.exists(target_dir):
                                shutil.rmtree(target_dir)
                        except Exception:
                            try:
                                logger.debug("Failed to remove target_dir after no-sfs skip", exc_info=True)
                            except Exception:
                                pass
                        # Ensure skipped demographics are not present in the SFS CSV in any form
                        try:
                            sfs.drop(index=key, inplace=True, errors='ignore')
                        except Exception:
                            try:
                                logger.debug("Failed to drop skipped key from sfs in memory after no-sfs skip", exc_info=True)
                            except Exception:
                                pass
                        try:
                            sfs_index_set.discard(key)
                        except Exception:
                            pass
                        include_sfs = False
                except Exception:
                    # best-effort only
                    pass

                # only read sfs.sfs if it exists; otherwise do not pass it in future runs
                if os.path.exists("sfs.sfs"):
                    try:
                        vals = read_sfs_file("sfs.sfs")
                        # if all zeros or empty, treat as skipped and do not add to CSV
                        if vals is None or np.allclose(vals, 0):
                                try:
                                    k = append_skipped(key, reason='all_zero_or_empty_sfs', source='3.sfs')
                                    if k:
                                        existing_skipped.add(k)
                                    else:
                                        try:
                                            logger.error(f"Failed to append skip for {key} when SFS was all zero or empty")
                                        except Exception:
                                            pass
                                except Exception:
                                    try:
                                        logger.error(f"Failed to append skip for {key} when SFS was all zero or empty")
                                    except Exception:
                                        pass
                                # remove any previously created data for this demographic/population
                                try:
                                    import shutil
                                    target_dir = os.path.join("data", species_folder_name, str(model_id), str(population))
                                    if os.path.exists(target_dir):
                                        shutil.rmtree(target_dir)
                                except Exception:
                                    try:
                                        logger.debug("Failed to remove target_dir", exc_info=True)
                                    except Exception:
                                        logger.debug("Failed to remove target_dir", exc_info=True)
                                # ensure skipped demographics are not written to sfs.csv
                                try:
                                    sfs.drop(index=key, inplace=True, errors='ignore')
                                except Exception:
                                    try:
                                        logger.debug("Failed to drop all-zero key from sfs in memory", exc_info=True)
                                    except Exception:
                                        pass
                                try:
                                    sfs_index_set.discard(key)
                                except Exception:
                                    pass
                        else:
                            sfs.loc[f"{model_id}={population}"] = vals
                            try:
                                sfs_index_set.add(f"{model_id}={population}")
                            except Exception:
                                pass
                    except Exception:
                        append_log(f"Failed to read sfs.sfs for {model_id} {population}; treating as no SFS produced.")
                        # ensure we do not add a NaN row for failed reads
                        try:
                            sfs.drop(index=key, inplace=True, errors='ignore')
                        except Exception:
                            try:
                                logger.debug("Failed to drop key from sfs after read failure", exc_info=True)
                            except Exception:
                                pass
                        try:
                            sfs_index_set.discard(key)
                        except Exception:
                            pass
                        include_sfs = False
                    try:
                        os.remove("sfs.sfs")
                    except Exception:
                        pass
                else:
                    append_log(f"No sfs.sfs produced for {model_id} {population}; will not pass --sfs on retries.")
                    # do not add a NaN row for demographics that produced no SFS
                    try:
                        sfs.drop(index=key, inplace=True, errors='ignore')
                    except Exception:
                        try:
                            logger.debug("Failed to drop key from sfs when no sfs produced", exc_info=True)
                        except Exception:
                            pass
                    try:
                        sfs_index_set.discard(key)
                    except Exception:
                        pass
                    include_sfs = False
            safe_write_sfs(f"data/{species_folder_name}/sfs.csv")
            pbar.update(1)
            pbar.refresh()
