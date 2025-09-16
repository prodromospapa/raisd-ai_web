
import stdpopsim
import pandas as pd
import os

species = "HomSap"  # example input, replace with actual input handling
species_dict = {sp.name: sp.id for sp in stdpopsim.all_species()}
if species in species_dict.values():
    species_full_name = [sp.name for sp in stdpopsim.all_species() if sp.id == species][0]
else:
    species_full_name = species

sfs_csv = pd.read_csv(f"data/{species_full_name.replace(' ', '_')}/sfs.csv", index_col=0)

not_trained = []
for demographic_model,population in sfs_csv.index.str.split('=').tolist():
    if not os.path.exists(f"data/{species_full_name.replace(' ', '_')}/{demographic_model}/{population}"):
        sfs_csv.drop(index=f"{demographic_model}={population}", inplace=True)
        not_trained.append(f"{demographic_model}={population}")

sfs_csv.to_csv(f"data/{species_full_name.replace(' ', '_')}/sfs.csv")

demographics_train = os.listdir(f"data/{species_full_name.replace(' ', '_')}")
not_sfs = []
for demographic_model in demographics_train:
    populations_train = os.listdir(f"data/{species_full_name.replace(' ', '_')}/{demographic_model}")
    for population in populations_train:
        if f"{demographic_model}={population}" not in sfs_csv.index:
            os.system(f"rm -r data/{species_full_name.replace(' ', '_')}/{demographic_model}/{population}")
            not_sfs.append(f"{demographic_model}={population}")


output_table = pd.DataFrame({
    "Not trained (no directory)": not_trained,
    "No SFS (no entry in sfs.csv)": not_sfs
})
output_table.to_csv(f"data/{species_full_name.replace(' ', '_')}/final_check.csv", index=False)# This script checks for consistency between demographic models/populations in sfs.csv
# and the directories containing trained models. It removes entries from sfs.csv without
# corresponding directories and deletes directories without corresponding sfs.csv entries.
# Finally, it outputs a summary of discrepancies to final_check.csv.
