import os
import multiprocessing as mp
import argparse
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import stdpopsim
import subprocess

# Hard-code your inputs here
parser = argparse.ArgumentParser(description="Simulate allele frequency spectra.")
parser.add_argument(
    "--species",
    type=str,
    required=True,
    help="Species name (e.g., 'MusMus')."
)
args = parser.parse_args()
species = args.species

parallel = 10
samples = 10_000

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
sfs = pd.DataFrame(columns= range(1,ploidy*(samples)))

os.system(f'mkdir -p "data/{species_folder_name}"')

if os.path.exists(f"data/{species_folder_name}/sfs.csv"):
    sfs = pd.read_csv(f"data/{species_folder_name}/sfs.csv",index_col=0)

total = sum(len(pops) for pops in demographic_models.values())
with tqdm(total=total, desc="Simulations", unit="run") as pbar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            if f"{model_id}={population}" in sfs.index:
                pbar.update(1)
                pbar.refresh()
                continue

            args = [
                "simulator.py",
                "--engine", "scrm",
                "--species-id", str(species),
                "--model-id", str(model_id),
                "--pop-order", str(population),
                "--sample-individuals", str(samples),
                "--chromosome", "21",#str(biggest_chromosome),
                "--replicates", "20",
                "--parallel", str(parallel),
                "--sfs", "sfs.sfs",
                "--sfs-normalized"
            ]
            try:
                subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                # switch engine to msprime and retry once
                if "--engine" in args:
                    i = args.index("--engine")
                    args[i+1] = "msprime"
                else:
                    args.extend(["--engine", "msprime"])
                    # halve parallel (integer division), ensure at least 1, and update args if present
                    if "--parallel" in args:
                        p_idx = args.index("--parallel")
                        new_parallel = max(1, int(int(args[p_idx+1]) // 2))
                        args[p_idx+1] = str(new_parallel)
                        parallel = new_parallel
                    else:
                        parallel = max(1, parallel // 2)
                        args.extend(["--parallel", str(parallel)])
                subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            sfs.loc[f"{model_id}={population}"] = read_sfs_file("sfs.sfs")
            os.remove("sfs.sfs")
            sfs.to_csv(f"data/{species_folder_name}/sfs.csv")
            pbar.update(1)
            pbar.refresh()