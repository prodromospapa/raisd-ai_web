import os
import multiprocessing as mp
import argparse
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import stdpopsim
import subprocess

# # Hard-code your inputs here
# parser = argparse.ArgumentParser(description="Simulate allele frequency spectra.")
# parser.add_argument(
#     "--species",
#     type=str,
#     required=True,
#     help="Species name (e.g., 'MusMus')."
# )
# args = parser.parse_args()
# species = args.species

species = "Homo sapiens"

parallel = 5
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

            args = [
                "simulator.py",
                "--engine", "msprime",
                "--species-id", str(species),
                "--model-id", str(model_id),
                "--pop-order", str(population),
                "--sample-individuals", str(samples),
                "--chromosome", "21",#str(biggest_chromosome),
                "--replicates", "5",
                "--parallel", str(parallel),
                "--sfs", "sfs.sfs",
                "--sfs-normalized",
                "--length", "10000000"
                ]
            try:
                subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                i = args.index("--engine")
                args[i+1] = "scrm"
                with open(f"data/{species_folder_name}/check_sfs_models.log", "a") as wf:
                    wf.write(f"Warning: msprime failed for {model_id} {population}.\n")
                subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            sfs.loc[f"{model_id}={population}"] = read_sfs_file("sfs.sfs")
            os.remove("sfs.sfs")
            sfs.to_csv(f"data/{species_folder_name}/sfs.csv")
            pbar.update(1)
            pbar.refresh()