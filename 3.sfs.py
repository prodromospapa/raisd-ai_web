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
                    help="If set, stop a simulator run when it exceeds this percent of total RAM, decrement --max-parallel and retry until 1.")
parser.add_argument("--max-parallel", type=int, default=1,
                    help="Maximum --parallel value to start with; will be decremented on RAM exceed and retried until 1.")
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
                "--sfs", "sfs.sfs",
                "--sfs-normalized"
            ]

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
                if MAX_RAM_PERCENT is None or psutil is None:
                    # no monitoring requested or psutil not available: run normally
                    try:
                        subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        break
                    except subprocess.CalledProcessError:
                        # fallback behavior: swap engine to scrm and retry once
                        i = args.index("--engine")
                        args[i+1] = "scrm"
                        with open(f"data/{species_folder_name}/check_sfs_models.log", "a") as wf:
                            wf.write(f"Warning: msprime failed for {model_id} {population}.\n")
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

                    if exceeded:
                        # reduce parallel and retry unless already at 1
                        if current_parallel > 1:
                            with open(f"data/{species_folder_name}/check_sfs_models.log", "a") as wf:
                                wf.write(f"Exceeded RAM {MAX_RAM_PERCENT}% for {model_id} {population} at parallel={current_parallel}; retrying with {current_parallel-1}\n")
                            current_parallel -= 1
                            continue
                        else:
                            # failed even at parallel==1; skip this model/population
                            with open(f"data/{species_folder_name}/check_sfs_models.log", "a") as wf:
                                wf.write(f"Skipped {model_id} {population}: exceeded RAM {MAX_RAM_PERCENT}% even at parallel=1.\n")
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
                            with open(f"data/{species_folder_name}/check_sfs_models.log", "a") as wf:
                                wf.write(f"Warning: simulator failed (rc={proc.returncode}) for {model_id} {population} at parallel={current_parallel}. Trying scrm.\n")
                            subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                            break
            if skipped:
                # mark skipped row with NaNs to preserve index and continue
                sfs.loc[f"{model_id}={population}"] = np.nan
            else:
                sfs.loc[f"{model_id}={population}"] = read_sfs_file("sfs.sfs")
                os.remove("sfs.sfs")
            sfs.to_csv(f"data/{species_folder_name}/sfs.csv")
            pbar.update(1)
            pbar.refresh()
