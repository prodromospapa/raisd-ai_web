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
parser.add_argument("--samples", type=int, default=10000,
                    help="Number of haploid samples per population to request when building the SFS (default: 10000).")
parser.add_argument("--replicates", type=int, default=5,
                    help="Number of replicates to request from the simulator for each model/population when building SFS (default: 5).")
parser.add_argument("--engine", type=str, default="msprime",
                    help="Default simulator engine to request (default: msprime).")
parser.add_argument("--chromosome", type=str, default="21",
                    help="Chromosome id to simulate; defaults to the largest chromosome if not set (default: 21).")
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
chromosome = cli_args.chromosome
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

# Build model->pop list where sampling_time == 0
demographic_models = {
    m.id: [p.name for p in species_std.get_demographic_model(m.id).populations]
    for m in species_std.demographic_models
}

species_folder_name = "".join(c if c.isalnum() or c in "-_ " else "_" for c in str(species_full_name)).strip().replace(" ", "_")

ploidy = species_std.ploidy


os.system(f'mkdir -p "data/{species_folder_name}"')

# shared record of failed parts (JSON-lines)
skipped_file = f"data/{species_folder_name}/failed_parts.jsonl"

import json

def load_skipped():
    skipped = set()
    try:
        if os.path.exists(skipped_file):
            with open(skipped_file, 'r', encoding='utf-8') as rf:
                for line in rf:
                    try:
                        obj = json.loads(line)
                        key = obj.get('key')
                        if key:
                            skipped.add(key)
                    except Exception:
                        # ignore malformed lines
                        continue
    except Exception:
        pass
    return skipped

def append_skipped(key, reason='skipped', source='3.sfs'):
    entry = {
        'key': key,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'reason': reason,
        'source': source,
    }
    try:
        with open(skipped_file, 'a', encoding='utf-8') as wf:
            wf.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except Exception:
        # fallback to writing a bare line in the old file location for maximum compatibility
        try:
            fallback = skipped_file.replace('.jsonl', '.txt')
            with open(fallback, 'a') as lf:
                lf.write(key + "\n")
        except Exception:
            # Give up on writing fallback; report to stdout instead
            try:
                print(f"Failed to write failed parts to {skipped_file}", flush=True)
            except Exception:
                pass

existing_skipped = load_skipped()

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
    samples = (sfs.shape[1] // ploidy)+1
    print(f"Resuming with {samples} samples per population")
else:
    sfs = pd.DataFrame(columns= range(1,ploidy*(samples)))

# previously removed warning log; no longer used
first = True
total = sum(len(pops) for pops in demographic_models.values())
with tqdm(total=total, desc="Simulations", unit="run") as pbar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            key = f"{model_id}={population}"
            if key in existing_skipped:
                pbar.update(1)
                pbar.refresh()
                first = False
                append_log(f"Skipping {key} because it is listed in skipped_demographics.jsonl")
                continue

            if key in sfs.index:
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
                if "--sims-per-work" not in args:
                    try:
                        pi = args.index("--parallel")
                        args.insert(pi + 2, str(current_sims_per_work))
                        args.insert(pi + 1, "--sims-per-work")
                    except Exception:
                        args += ["--sims-per-work", str(current_sims_per_work)]

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
                try:
                    proc = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

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

                except Exception as e:
                    append_log(f"Exception while running simulator for {model_id} {population}: {e}")
                    exceeded = False
                    rc = getattr(proc, 'returncode', 1) if proc is not None else 1

                # Decide next action based on whether memory was exceeded or process exit code
                if exceeded:
                    # reduce parallel and retry (keep sims-per-work=1)
                    if current_parallel > 1:
                        current_parallel -= 1
                        append_log(f"Retrying with lower parallel: {current_parallel} (keeping --sims-per-work=1)")
                        continue
                    else:
                        append_log(f"Skipped {model_id} {population}: exceeded RAM {MAX_RAM_PERCENT}% even at parallel=1 and --sims-per-work=1.")
                        skipped = True
                        break

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
                    skipped = True
                    break

            else:
                # exceeded max attempts
                append_log(f"Giving up after {attempts} attempts for {model_id} {population}; marking skipped")
                skipped = True
            if skipped:
                # Record the skip in the canonical JSONL skipped file and remove any existing data
                try:
                    append_skipped(key, reason='failed_to_produce_sfs', source='3.sfs')
                    existing_skipped.add(key)
                except Exception:
                    # if we can't write the canonical skip, log an error for troubleshooting via logger
                    logger.error(f"Failed to append {key} to {skipped_file}")
                # remove any previously created data for this demographic/population (best-effort)
                try:
                    import shutil
                    target_dir = os.path.join("data", species_folder_name, str(model_id), str(population))
                    if os.path.exists(target_dir):
                        shutil.rmtree(target_dir)
                except Exception:
                    # keep going even if removal fails
                    try:
                        logger.debug(f"Failed to remove target_dir {target_dir}", exc_info=True)
                    except Exception:
                        logger.debug("Failed to remove target_dir (target_dir undefined)", exc_info=True)
            else:
                # only read sfs.sfs if it exists; otherwise do not pass it in future runs
                if os.path.exists("sfs.sfs"):
                    try:
                        vals = read_sfs_file("sfs.sfs")
                        # if all zeros or empty, treat as skipped and do not add to CSV
                        if vals is None or np.allclose(vals, 0):
                                append_skipped(key, reason='all_zero_or_empty_sfs', source='3.sfs')
                                existing_skipped.add(key)
                                # remove any previously created data for this demographic/population
                                try:
                                    import shutil
                                    target_dir = os.path.join("data", species_folder_name, str(model_id), str(population))
                                    if os.path.exists(target_dir):
                                        shutil.rmtree(target_dir)
                                except Exception:
                                    try:
                                        logger.debug(f"Failed to remove target_dir {target_dir}", exc_info=True)
                                    except Exception:
                                        logger.debug("Failed to remove target_dir (target_dir undefined)", exc_info=True)
                                sfs.loc[f"{model_id}={population}"] = np.nan
                        else:
                            sfs.loc[f"{model_id}={population}"] = vals
                    except Exception:
                        append_log(f"Failed to read sfs.sfs for {model_id} {population}; treating as no SFS produced.")
                        sfs.loc[f"{model_id}={population}"] = np.nan
                        include_sfs = False
                    try:
                        os.remove("sfs.sfs")
                    except Exception:
                        pass
                else:
                    append_log(f"No sfs.sfs produced for {model_id} {population}; will not pass --sfs on retries.")
                    sfs.loc[f"{model_id}={population}"] = np.nan
                    include_sfs = False
            sfs.to_csv(f"data/{species_folder_name}/sfs.csv")
            pbar.update(1)
            pbar.refresh()
