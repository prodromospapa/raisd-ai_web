import os
import multiprocessing as mp
import argparse
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import stdpopsim
import subprocess
import time
import logging
import signal
import psutil

species = "Homo sapiens"

samples = 10_000

# CLI: allow optional max RAM percent enforcement
parser = argparse.ArgumentParser()
parser.add_argument("--max-ram-percent", type=float,default=0.5,
                    help="If set, stop a simulator run when it exceeds this percent of total RAM; the run will be skipped (no lower-parallel retry).")
parser.add_argument("--max-parallel", type=int, default=1,
                    help="Maximum --parallel value to start with.")
cli_args, _ = parser.parse_known_args()
MAX_RAM_PERCENT = cli_args.max_ram_percent
MAX_PARALLEL = cli_args.max_parallel

# basic logger
logger = logging.getLogger("check_sfs_models")
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

if os.path.exists(f"data/{species_folder_name}/sfs.csv"):
    sfs = pd.read_csv(f"data/{species_folder_name}/sfs.csv",index_col=0)
    samples = (sfs.shape[1] // ploidy)+1
    print(f"Resuming with {samples} samples per population")
else:
    sfs = pd.DataFrame(columns= range(1,ploidy*(samples)))

#os.remove(f"data/{species_folder_name}/warning.log") if os.path.exists(f"data/{species_folder_name}/warning.log") else None
first = True
total = sum(len(pops) for pops in demographic_models.values())
with tqdm(total=total, desc="Simulations", unit="run") as pbar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            if f"{model_id}={population}" in sfs.index:
                pbar.update(1)
                pbar.refresh()
                first = False
                continue
            if first:
                os.remove(f"data/{species_folder_name}/check_sfs_models.log") if os.path.exists(f"data/{species_folder_name}/check_sfs_models.log") else None
            first = False

            # helper to append logs both to species folder and current working dir
            def append_log(msg: str):
                ts = time.strftime("%Y-%m-%d %H:%M:%S")
                line = f"[{ts}] {msg}\n"
                species_log = f"data/{species_folder_name}/check_sfs_models.log"
                local_log = os.path.join(os.getcwd(), "check_sfs_models.log")
                try:
                    with open(species_log, "a") as wf:
                        wf.write(line)
                except Exception:
                    pass
                try:
                    with open(local_log, "a") as wf:
                        wf.write(line)
                except Exception:
                    pass

            # prepare base args; we'll modify --parallel dynamically on retries
            base_args = [
                "simulator.py",
                "--engine", "scrm",
                "--species-id", str(species),
                "--model-id", str(model_id),
                "--pop-order", str(population),
                "--sample-individuals", str(samples),
                "--chromosome", "21",#str(biggest_chromosome),
                "--replicates", "5",
                "--parallel",
                "--sfs-normalized"
            ]
            # whether to include --sfs sfs.sfs in simulator args for this model/pop
            include_sfs = True
            current_parallel = MAX_PARALLEL
            skipped = False
            # We'll attempt runs, lowering parallel if RAM threshold is exceeded
            while True:
                args = list(base_args)
                # insert parallel value in the args list after the --parallel flag
                try:
                    pi = args.index("--parallel")
                    args.insert(pi + 1, str(current_parallel))
                except ValueError:
                    args += ["--parallel", str(current_parallel)]

                # start the subprocess
                # build args list; include --sfs only if requested
                if include_sfs:
                    # insert sfs args right after the --parallel value later
                    pass

                if MAX_RAM_PERCENT is None or psutil is None:
                    # no monitoring requested or psutil not available: run normally
                    try:
                        # ensure args include sfs if requested
                        if include_sfs and "--sfs" not in args:
                            try:
                                pi = args.index("--parallel")
                                args.insert(pi+2, "sfs.sfs")
                                args.insert(pi+1, "--sfs")
                            except Exception:
                                args += ["--sfs", "sfs.sfs"]
                        subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        break
                    except subprocess.CalledProcessError:
                        # fallback behavior: swap engine to scrm and retry once
                        i = args.index("--engine")
                        args[i+1] = "scrm"
                        append_log(f"Warning: msprime failed for {model_id} {population}.")
                        subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        break
                else:
                    # monitored run using psutil
                    proc = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    try:
                        p = psutil.Process(proc.pid)
                    except Exception:
                        p = None

                    exceeded = False
                    try:
                        # monitor until process ends
                        while proc.poll() is None:
                            if p is not None:
                                try:
                                    # check memory percent of process tree
                                    mem = p.memory_percent()
                                    # include children
                                    for ch in p.children(recursive=True):
                                        try:
                                            mem += ch.memory_percent()
                                        except Exception:
                                            pass
                                except Exception:
                                    mem = 0.0

                                if MAX_RAM_PERCENT is not None and mem > MAX_RAM_PERCENT:
                                    # exceed threshold: kill process tree and retry with lower parallel
                                    exceeded = True
                                    # terminate children first
                                    for ch in p.children(recursive=True):
                                        try:
                                            ch.terminate()
                                        except Exception:
                                            pass
                                    try:
                                        p.terminate()
                                    except Exception:
                                        pass
                                    # give them a moment, then kill if needed
                                    time.sleep(0.5)
                                    for ch in p.children(recursive=True):
                                        if ch.is_running():
                                            try:
                                                ch.kill()
                                            except Exception:
                                                pass
                                    if p.is_running():
                                        try:
                                            p.kill()
                                        except Exception:
                                            pass
                                    break
                            time.sleep(0.2)
                    finally:
                        # ensure subprocess is reaped
                        try:
                            proc.wait(timeout=1)
                        except Exception:
                            pass
                    # after monitoring loop, decide what happened
                    if exceeded:
                        # Do not retry with a lower parallel: skip this model/population
                        append_log(f"Skipped {model_id} {population}: exceeded RAM {MAX_RAM_PERCENT}% at parallel={current_parallel}; not retrying with lower parallel.")
                        skipped = True
                        break
                    else:
                        # process finished normally; check return code
                        if proc.returncode == 0:
                            break
                        else:
                            # non-zero exit: write warning and attempt with engine scrm once
                            try:
                                i = args.index("--engine")
                                args[i+1] = "scrm"
                            except Exception:
                                pass
                            append_log(f"Warning: simulator failed (rc={proc.returncode}) for {model_id} {population} at parallel={current_parallel}. Trying scrm.")
                            # ensure sfs flag only present if requested
                            if include_sfs and "--sfs" not in args:
                                try:
                                    pi = args.index("--parallel")
                                    args.insert(pi+2, "sfs.sfs")
                                    args.insert(pi+1, "--sfs")
                                except Exception:
                                    args += ["--sfs", "sfs.sfs"]
                            subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                            break
            if skipped:
                # mark skipped row with NaNs to preserve index and continue
                sfs.loc[f"{model_id}={population}"] = np.nan
            else:
                # only read sfs.sfs if it exists; otherwise do not pass it in future runs
                if os.path.exists("sfs.sfs"):
                    try:
                        sfs.loc[f"{model_id}={population}"] = read_sfs_file("sfs.sfs")
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
