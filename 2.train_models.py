import stdpopsim
import pandas as pd
import warnings
from tqdm.auto import tqdm
import os
import gzip
import re
import subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
import shutil
import os, sys, contextlib, logging, warnings
import multiprocessing as mp
import traceback
import argparse

########################################################################
# --------------------
# config
# --------------------
parser = argparse.ArgumentParser(description="Simulate allele frequency spectra.")
parser.add_argument(
    "--species",
    type=str,
    required=True,
    help="Species name (e.g., 'MusMus')."
)
args = parser.parse_args()
species = args.species
pop_size = 1_000
n_simulations = 10
decrease_factor = 0.9
max_time = 600 # in seconds
slim_scaling_factor_og = 1_000_000
MAX_WORKERS = int(os.environ.get("MAX_WORKERS", os.cpu_count() or 2))

# RAiSD-AI parameters
window = 500
epochs = 100
########################################################################


# mute warnings and errors
@contextlib.contextmanager
def hard_silence():
    # Save original fds
    orig_stdout_fd = os.dup(1)
    orig_stderr_fd = os.dup(2)
    try:
        with open(os.devnull, 'w') as devnull:
            os.dup2(devnull.fileno(), 1)           # redirect fd 1 (stdout)
            os.dup2(devnull.fileno(), 2)           # redirect fd 2 (stderr)
            # Keep Python objects in sync with fds
            sys.stdout = open(os.devnull, 'w')
            sys.stderr = open(os.devnull, 'w')

            # also mute logging + warnings
            logging.disable(logging.CRITICAL)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                yield
    finally:
        # Restore fds
        os.dup2(orig_stdout_fd, 1)
        os.dup2(orig_stderr_fd, 2)
        os.close(orig_stdout_fd)
        os.close(orig_stderr_fd)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        logging.disable(logging.NOTSET)

warnings.simplefilter("ignore")
# mute warnings and errors
def run(cmd, cwd=None, quiet=True):
    stdout = subprocess.DEVNULL if quiet else None
    stderr = subprocess.DEVNULL if quiet else None
    subprocess.run(cmd, cwd=cwd, check=True, stdout=stdout, stderr=stderr)

def geo_backoff(start, factor, min_value=1):
    current = start
    while current > min_value:
        yield int(current)
        current *= factor
    yield 0

def sweep_def(contig, sweep_pos, population, event_time):
    ss_id = f"hard_sweep_{population}"
    contig.add_single_site(id=ss_id, coordinate=sweep_pos)
    return stdpopsim.selective_sweep(
        single_site_id=ss_id,
        population=population,
        selection_coeff=1,
        mutation_generation_ago=event_time,
        min_freq_at_end=1
    )

def simulate_worker(
    species, demographic_model, population, chromosome,
    type_, workdir, vcf_i=None, sweep_pos=None, window=None, mean_snp=None,
    slim_scaling_factor=None, event_time=None
):
    """Run a single replicate; returns written VCF basename or None for test runs."""
    species_std = stdpopsim.get_species(species)
    model = species_std.get_demographic_model(demographic_model)
    os.makedirs(workdir, exist_ok=True)
    samples_per_pop = {population: pop_size}

    if type_ == "test":
        contig = species_std.get_contig(chromosome, mutation_rate=model.mutation_rate, right=1_000)
        events = sweep_def(contig, 500, population, event_time)
        engine = stdpopsim.get_engine("slim")
        engine.simulate(
            model,
            contig,
            samples_per_pop,
            slim_scaling_factor=slim_scaling_factor,
            slim_burn_in=1,
            extended_events=events
        )
        return None

    if window is None or sweep_pos is None:
        raise ValueError("window and sweep_pos are required for non-test simulations")

    left_snps = 0
    right_snps = 0
    multiplier = 1.0
    m_snp = max(float(mean_snp or 1), 1.0)
    ts = None

    while left_snps <= (window // 2) or right_snps <= (window // 2):
        left = max(0, int(sweep_pos - (multiplier * (window * m_snp) // 2)))
        right = int(sweep_pos + (multiplier * (window * m_snp) // 2))
        contig = species_std.get_contig(chromosome, mutation_rate=model.mutation_rate, left=left, right=right)
        multiplier += 0.1

        engine = stdpopsim.get_engine("slim")
        if type_ == 'neutral':
            ts = engine.simulate(
                model,
                contig,
                samples_per_pop,
                slim_scaling_factor=slim_scaling_factor,
                slim_burn_in=1,
                extended_events=None
            )
        elif type_ == 'sweep':
            events = sweep_def(contig, sweep_pos, population, event_time)
            ts = engine.simulate(
                model,
                contig,
                samples_per_pop,
                slim_scaling_factor=slim_scaling_factor,
                slim_burn_in=1,
                extended_events=events
            )
        else:
            raise ValueError(f"Unknown type_: {type_}")

        positions = [v.position for v in ts.variants()]
        if not positions:
            continue
        closest = min(positions, key=lambda x: abs(x - sweep_pos))
        left_snps = positions.index(closest)
        right_snps = len(positions) - left_snps - 1

    vcf_path = os.path.join(workdir, f"{type_}_{vcf_i}.vcf.gz")
    with gzip.open(vcf_path, 'wt') as f:
        if ts is None:
            raise RuntimeError("Simulation did not produce a tree sequence (ts is None)")
        ts.write_vcf(f, contig_id=f"{chromosome}_{vcf_i}")
    return os.path.basename(vcf_path)

def vcf2bin(workdir, type_, window, sweep_pos):
    run([
        "RAiSD-AI-ZLIB",
        "-n", "bin",
        "-I", f"{type_}.vcf.gz",
        "-w", str(window),
        "-its", str(int(sweep_pos)),
        "-op", "IMG-GEN",
        "-icl", f"{type_}TR",
        "-bin", "-typ", "1"
    ], cwd=workdir)
    # clean chatter
    for f in os.listdir(workdir):
        if f.startswith("RAiSD_Info."):
            try:
                os.remove(os.path.join(workdir, f))
            except OSError:
                pass
    # clean vcf file
    try:
        os.remove(os.path.join(workdir, f"{type_}.vcf.gz"))
    except OSError:
        pass

def concat_vcfs(out_path, in_paths, workdir):
    cmd = ["bcftools", "concat", "-o", out_path] + list(in_paths)
    run(cmd, quiet=True, cwd=workdir)
    # clean up
    for f in in_paths:
        try:
            os.remove(os.path.join(workdir, f))
        except OSError:
            pass

def bin2model(workdir, epochs):
    run([
        "RAiSD-AI",
        "-n", "model",
        "-I", "RAiSD_Images.bin",
        "-f", "-op", "MDL-GEN",
        "-e", str(int(epochs)),
        "-arc", "FASTER-NN"
    ], cwd=workdir)

# tune parameters
def _r(q,f,a,k):
    try: q.put(("ok", f(*a, **(k or {}))))
    except Exception as e: q.put(("err", (repr(e), traceback.format_exc())))
    
def _run(f,a=(),k=None,t=None):
    q=mp.Queue(); p=mp.Process(target=_r, args=(q,f,a,k)); p.start(); p.join(t)
    if p.is_alive():
        p.terminate(); p.join(2)
        if p.is_alive() and hasattr(p,"kill"): p.kill(); p.join(1)
        raise TimeoutError
    if not q.empty():
        s,pay=q.get_nowait()
        if s=="ok": return pay
        raise RuntimeError(f"Worker raised: {pay[0]}\nTraceback:\n{pay[1]}")

def tune_parameters(model, population, chromosome):
    s=slim_scaling_factor_og; et_out=None
    for et in geo_backoff(int(1e9), factor=decrease_factor):
        try:
            _run(simulate_worker, (species, model.id, population, chromosome, "test", "."), {"slim_scaling_factor": s, "event_time": et}, max_time)
            et_out=et; break
        except TimeoutError:
            et_out,s=None,None; break
        except Exception as e:
            x=str(e)
            m=re.search(r'value\s+(-?\d+)\s+for a tick index', x)
            if m: et=int(et + (int(m.group(1))*s) - s); continue
            if re.search(r'zero at tick', x): s=max(1, s//2); continue
            m=re.search(r'only\s+(\d+)\s+individuals', x)
            if m: s=max(1, int(int(m.group(1))*s/pop_size)); continue
            # SLiM error when trying to target a population that does not yet exist at the chosen event_time.
            # This manifests as: 'undefined identifier p2' (or p3, etc.). In this case just try a smaller event_time.
            if re.search(r'undefined identifier p\d+', x):
                # Simply continue to the next (smaller) event_time candidate
                continue
            raise
    return et_out, s
# tune parameters

# --------------------
# main
# --------------------
if __name__ == "__main__":
    species_std = stdpopsim.get_species(species)

    # Which (model, pop) have sampling_time == 0?
    demographic_models = {
        m.id: [p.name for p in species_std.get_demographic_model(m.id).populations]
        for m in species_std.demographic_models
    }
    demographic_models_pass = {}
    for demographic_model, populations in demographic_models.items():
        model = species_std.get_demographic_model(demographic_model)
        for population in populations:
            pop_obj = next(p for p in model.populations if p.name == population)
            if pop_obj.extra_metadata.get('sampling_time', None) == 0:
                demographic_models_pass.setdefault(demographic_model, []).append(population)
    
    mean_snp_df = pd.read_csv(f"data/{species_std.name.replace(' ','_')}/mean.csv", index_col=0)

    # create chromosome.txt
    if not os.path.exists(f"data/{species_std.name.replace(' ','_')}/chromosomes.txt"):
        chromosomes = [c.id for c in species_std.genome.chromosomes]
        remove_chr = input(f"Which of {chromosomes} to remove (comma-separated): ")
        remove_chr = [c.strip() for c in remove_chr.split(",") if c.strip()]
        chromosomes = [c for c in chromosomes if c not in remove_chr]
        with open(f"data/{species_std.name.replace(' ','_')}/chromosomes.txt", 'w') as f:
            for chromosome in chromosomes:
                f.write(chromosome+"\n")
    else:
        with open(f"data/{species_std.name.replace(' ','_')}/chromosomes.txt", 'r') as f:
            chromosomes = [line.strip() for line in f]

    # Pre-compute total jobs (replicates to run), skipping already-trained chromosomes
    total_jobs = 0
    sfs = pd.read_csv(f"data/{species_std.name.replace(' ','_')}/sfs.csv", index_col=0)
    for demographic_model, populations in demographic_models_pass.items():
        for population in populations:
            for chromosome in chromosomes:
                workdir = f"data/{species_std.name.replace(' ','_')}/{demographic_model}/{population}/{chromosome}"
                if os.path.exists(f"{workdir}/RAiSD_Model.model/model.pt"):
                    continue
                if f"{demographic_model}_{population}" not in sfs.index:
                    continue
                total_jobs += n_simulations * 2  # neutral + sweep

    pbar = tqdm(total=total_jobs, desc="Simulations", unit="rep")

    for demographic_model, populations in demographic_models_pass.items():
        model = species_std.get_demographic_model(demographic_model)
        for population in populations:
            if all(os.path.exists(f"data/{species_std.name.replace(' ','_')}/{demographic_model}/{population}/{chromosome}/RAiSD_Model.model/model.pt") for chromosome in chromosomes):
                continue
            if f"{demographic_model}_{population}" not in sfs.index:
                continue
            
            with hard_silence():
                event_time, slim_scaling_factor = tune_parameters(model, population, chromosomes[-1])

            if event_time is None or slim_scaling_factor is None:
                pbar.update(len(chromosomes) * n_simulations * species_std.ploidy)
                sfs = pd.read_csv(f"data/{species_std.name.replace(' ','_')}/sfs.csv", index_col=0)
                sfs.drop(index=f"{demographic_model}_{population}", inplace=True)
                sfs.to_csv(f"data/{species_std.name.replace(' ','_')}/sfs.csv")
                continue

            # Extract scalar mean_snp robustly (row could be Series/DataFrame/scalar)
            mr = mean_snp_df.loc[f"{demographic_model}_{population}"]
            if isinstance(mr, pd.DataFrame):
                mean_snp = float(mr.values.ravel()[0])
            elif isinstance(mr, pd.Series):
                mean_snp = float(mr.iloc[0])
            else:
                mean_snp = float(mr)
            for chromosome in chromosomes:
                workdir = f"data/{species_std.name.replace(' ','_')}/{demographic_model}/{population}/{chromosome}"
                # Skip if already trained; otherwise wipe and start fresh
                if os.path.exists(f"{workdir}/RAiSD_Model.model/model.pt"):
                    continue
                shutil.rmtree(workdir, ignore_errors=True)

                contig_L = species_std.get_contig(chromosome)
                L = int(contig_L.length)
                sweep_pos = L // 2
                os.makedirs(workdir, exist_ok=True)
                
                # ---------- SWEEP (parallel) ----------
                sweep_attempts = 0
                while True:
                    try:
                        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
                            futures = [
                                ex.submit(
                                    simulate_worker,
                                    species, demographic_model, population, chromosome,
                                    "sweep", workdir, rep, sweep_pos, window, mean_snp,
                                    slim_scaling_factor, event_time
                                )
                                for rep in range(n_simulations)
                            ]
                            for _ in as_completed(futures):
                                pbar.update(1)
                        sweep_files = np.array([f"sweep_{rep}.vcf.gz" for rep in range(n_simulations)], dtype=object)
                        concat_vcfs("sweep.vcf.gz", sweep_files, workdir)
                        vcf2bin(workdir, "sweep", window, sweep_pos)
                        break
                    except Exception:
                        sweep_attempts += 1
                        if sweep_attempts > 5:
                            raise
                        continue


                # ---------- NEUTRAL (parallel) ----------
                neutral_attempts = 0
                while True:
                    try:
                        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
                            futures = [
                                ex.submit(
                                    simulate_worker,
                                    species, demographic_model, population, chromosome,
                                    "neutral", workdir, rep, sweep_pos, window, mean_snp,
                                    slim_scaling_factor, None
                                )
                                for rep in range(n_simulations)
                            ]
                            for _ in as_completed(futures):
                                pbar.update(1)
                        neutral_files = np.array([f"neutral_{rep}.vcf.gz" for rep in range(n_simulations)], dtype=object)
                        concat_vcfs("neutral.vcf.gz", neutral_files, workdir)
                        vcf2bin(workdir, "neutral", window, sweep_pos)
                        break
                    except Exception:
                        neutral_attempts += 1
                        if neutral_attempts > 5:
                            raise
                        continue
            
                # Train once per chromosome (sequential)
                bin2model(workdir, epochs)

    pbar.close()
    species_dir_name = f"data/{species_std.name.replace(' ','_')}"
    os.system(f'rm "{species_dir_name}/mean.csv"')

    if sfs.empty:
        print("I can't train any models with these demographics")
        os.system(f'rm -rf "{species_dir_name}"')
        sys.exit(1)
