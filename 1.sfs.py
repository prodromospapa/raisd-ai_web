import os
import warnings
import math
import gc
import multiprocessing as mp

import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import stdpopsim

import argparse

# Avoid thread oversubscription issues in child processes
os.environ.setdefault("OMP_NUM_THREADS", "1")

warnings.simplefilter("ignore")

# ---------------- user-configurable (kept same names) ----------------
parser = argparse.ArgumentParser(description="Simulate allele frequency spectra.")
parser.add_argument(
    "--species",
    type=str,
    required=True,
    help="Species name (e.g., 'MusMus')."
)
args = parser.parse_args()
species = args.species
pop_size = 10_000

# SLiM scaling factor (mirrors usage in 2.train_models.py for neutral sims)
SLIM_SCALING_FACTOR = int(os.getenv("SLIM_SCALING_FACTOR", "1000000"))
TUNE_MAX_ATTEMPTS = int(os.getenv("TUNE_MAX_ATTEMPTS", "12"))
TUNE_CONTIG_BP = int(os.getenv("TUNE_CONTIG_BP", "10000"))  # small test region
TUNE_TIMEOUT = int(os.getenv("TUNE_TIMEOUT", "120"))  # seconds (not enforced strictly, just guidance)

# control RAM / parallelism via env vars
WINDOW_BP = int(os.getenv("WINDOW_BP", 5_000_000))   # slice size on the largest chrom
MAX_PROCS = int(os.getenv("MAX_PROCS", 0))           # 0 => all cores
# --------------------------------------------------------------------


def biggest_chrom(species_std):
    # return the ID (e.g., "1") of the numerically largest chromosome by length (no deprecation)
    chrom = [c for c in species_std.genome.chromosomes]
    return max(chrom, key=lambda c: int(c.length)).id


def average_dist(ts):
    snp_positions = ts.tables.sites.position
    return (snp_positions[-1] - snp_positions[0]) / (len(snp_positions) - 1)


# kept for structure; main run uses the parallel window worker below
def simulate(samples_per_pop, model, contig):
    """(Unused helper) Run a single neutral simulation (SLiM)."""
    engine = stdpopsim.get_engine("slim")
    ts = engine.simulate(
        model,
        contig,
        samples_per_pop,
        slim_scaling_factor=SLIM_SCALING_FACTOR,
        slim_burn_in=1,
        extended_events=None,
    )
    return ts.allele_frequency_spectrum(polarised=True, span_normalise=False), average_dist(ts)


def _try_simulation(species_name, demog_id, population, chrom_id, scaling_factor, pop_size_local, length):
    """Run a very small neutral SLiM simulation to test a scaling_factor.
    Returns True on success, False on recoverable error, raises on fatal errors."""
    species_std_local = stdpopsim.get_species(species_name)
    model = species_std_local.get_demographic_model(demog_id)
    engine = stdpopsim.get_engine("slim")
    # Use a tiny slice from start of chromosome
    contig = species_std_local.get_contig(
        chrom_id,
        left=0,
        right=length,
        mutation_rate=model.mutation_rate,
    )
    try:
        engine.simulate(
            model,
            contig,
            {population: pop_size_local},
            slim_scaling_factor=scaling_factor,
            slim_burn_in=1,
            extended_events=None,
        )
        return True
    except Exception as e:  # noqa: BLE001
        msg = str(e)
        # Patterns indicating we should reduce scaling_factor
        if (
            "zero at tick" in msg
            or "only" in msg and "individuals" in msg
            or "not enough individuals" in msg
            or "cannot find individual" in msg
        ):
            return False
        # Unknown error -> re-raise
        raise


def tune_scaling_factor(species_name, demog_id, population, chrom_id, start_scaling, pop_size_local):
    """Tune SLiM scaling factor by progressively shrinking until a test sim succeeds.

    Strategy: Try start_scaling, then halves until success or 1. Returns working factor.
    """
    scaling = max(1, int(start_scaling))
    attempts = 0
    while attempts < TUNE_MAX_ATTEMPTS and scaling >= 1:
        success = _try_simulation(
            species_name, demog_id, population, chrom_id, scaling, pop_size_local, TUNE_CONTIG_BP
        )
        if success:
            return scaling
        # reduce (geometric) if failed; ensure at least 1
        new_scaling = scaling // 2
        if new_scaling == scaling:
            new_scaling = scaling - 1
        scaling = max(1, new_scaling)
        attempts += 1
    return scaling  # may be 1


# -------------------- per-window worker --------------------
def _window_job(args):
    """
    Simulate one (model, population, window) chunk.
    Returns raw SFS (not per-bp), and stats to compute mean SNP spacing.
    """
    species_name, demog_id, population, pop_size_local, chrom_id, left, right, slim_scaling_factor = args

    species_std_local = stdpopsim.get_species(species_name)
    # Use SLiM engine for neutral simulations (mirrors neutral path in 2.train_models.py)
    engine_local = stdpopsim.get_engine("slim")
    model = species_std_local.get_demographic_model(demog_id)

    contig = species_std_local.get_contig(
        chrom_id,
        left=int(left),
        right=int(right),
        mutation_rate=model.mutation_rate,  # model-specific mutation rate
    )

    ts = engine_local.simulate(
        model,
        contig,
        {population: pop_size_local},
        slim_scaling_factor=slim_scaling_factor,
        slim_burn_in=1,
        extended_events=None,
    )

    sfs_raw = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)

    pos = ts.tables.sites.position
    m = len(pos)
    span = float(pos[m - 1]) - float(pos[0]) if m > 1 else 0.0
    gaps = (m - 1) if m > 1 else 0

    # free asap
    del ts
    gc.collect()

    key = f"{demog_id}_{population}"
    return key, sfs_raw, span, gaps
# -----------------------------------------------------------


if __name__ == "__main__":
    # Use spawn to avoid fork-related hangs with native libs
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        pass  # already set

    species_std = stdpopsim.get_species(species)

    biggest_chromosome = biggest_chrom(species_std)
    full_contig = species_std.get_contig(biggest_chromosome)
    L = int(full_contig.length)  # full chromosome length (for per-bp normalization)

    # Build model->pop list where sampling_time == 0
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
                demographic_models_pass[demographic_model] = (
                    demographic_models_pass.get(demographic_model, []) + [population]
                )
    if demographic_models_pass == {}:
        print("No valid demographic models or populations found.")
        exit()
    # Pre-make window coordinates along the named chromosome
    if WINDOW_BP <= 0 or WINDOW_BP >= L:
        WINDOW_BP = L
    n_windows = max(1, math.ceil(L / WINDOW_BP))
    windows = [(w * WINDOW_BP, min(L, (w + 1) * WINDOW_BP)) for w in range(n_windows)]

    # Tune scaling factor per (model, population)
    scaling_factors = {}
    for demographic_model, populations in demographic_models_pass.items():
        for population in populations:
            scaling_factors[(demographic_model, population)] = tune_scaling_factor(
                species, demographic_model, population, biggest_chromosome, SLIM_SCALING_FACTOR, pop_size
            )

    # Prepare window-level jobs (balances load; avoids last-2-straggler hang)
    jobs = []
    for demographic_model, populations in demographic_models_pass.items():
        for population in populations:
            for left, right in windows:
                jobs.append((
                    species,
                    demographic_model,
                    population,
                    pop_size,
                    biggest_chromosome,
                    left,
                    right,
                    scaling_factors[(demographic_model, population)],
                ))

    # Accumulators per key
    sfs_len = species_std.ploidy * pop_size + 1
    sfs_acc = {}   # key -> np.ndarray
    span_acc = {}  # key -> float
    gaps_acc = {}  # key -> int
    keys = []      # preserve insertion order for final DataFrame

    max_workers = (os.cpu_count() or 1) if MAX_PROCS <= 0 else min(MAX_PROCS, os.cpu_count() or 1)

    # Process in parallel with spawn context; recycle workers to avoid leaks
    with tqdm(total=len(jobs), desc="Simulations", unit="win") as pbar:
        with mp.get_context("spawn").Pool(processes=max_workers) as pool:
            results = pool.imap_unordered(_window_job, jobs)
            for result in results:
                try:
                    key, sfs_raw, span, gaps = result
                except Exception as e:
                    pbar.write(f"[error] {e!r}")
                    pbar.update(1)
                    continue

                if key not in sfs_acc:
                    sfs_acc[key] = np.zeros(sfs_len, dtype=np.float64)
                    span_acc[key] = 0.0
                    gaps_acc[key] = 0
                    keys.append(key)

                sfs_arr = np.asarray(sfs_raw, dtype=np.float64)
                if sfs_arr.shape[0] != sfs_len:
                    tmp = np.zeros(sfs_len, dtype=np.float64)
                    tmp[:min(sfs_len, sfs_arr.shape[0])] = sfs_arr[:min(sfs_len, sfs_arr.shape[0])]
                    sfs_arr = tmp

                sfs_acc[key] += sfs_arr
                span_acc[key] += float(span)
                gaps_acc[key] += int(gaps)

                pbar.update(1)

    # Materialize final DataFrames (same names/structure/filenames)
    sfs_df = pd.DataFrame(columns=range(0, sfs_len))
    mean_df = pd.DataFrame(columns=["mean"])
    for key in keys:
        # per-snp SFS (normalize by number of SNPs)
        sfs_df.loc[key] = sfs_acc[key] / sfs_acc[key].sum()

        # pooled mean SNP spacing across windows
        if gaps_acc[key] > 0:
            mean_df.loc[key] = span_acc[key] / gaps_acc[key]
        else:
            mean_df.loc[key] = np.nan
    os.system(f'mkdir -p "data/{species_std.name.replace(" ","_")}"')
    sfs_df.to_csv(f"data/{species_std.name.replace(' ','_')}/sfs.csv")
    mean_df.to_csv(f"data/{species_std.name.replace(' ','_')}/mean.csv")
