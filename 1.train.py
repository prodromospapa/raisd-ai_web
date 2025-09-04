import os
import stdpopsim as sps
import subprocess

species = "Homo sapiens"
train_sample_individuals = 100
train_replicates = 1_000
sfs_sample_individuals = 10_000
sfs_replicates = 1_000
sel_2Ns = 500

window = 1000
ips = 11 # odd
iws = 10
epochs = 10

engine = "msms"


def run_simulation(engine,species_id,model_id,pop_order,train_sample_individuals,window,ips,iws,chromosome,train_replicates,sel_2Ns):
    # build target run directory: current_path/data/species_folder_name/model_id/population/chromosome
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(pop_order), str(chromosome))
    os.makedirs(out_dir, exist_ok=True)
    args = [
        "simulator.py",
        "--engine", str(engine),
        "--species-id", str(species_id),
        "--model-id", str(model_id),
        "--pop-order", str(pop_order),
        "--sample-individuals", str(train_sample_individuals),
        "--target-snps", str(int(window + ((ips/2) * iws)*1.1)),
        "--chromosome", str(chromosome),
        "--replicates", str(train_replicates),
        "--output", "sweep.ms",
        "--sweep-pos", "0.5",
        "--sel-2Ns", str(sel_2Ns),
        "--threads", str(os.cpu_count()),
        "--target-snps-tol", "0.1",
        "--paired-neutral",
        "--neutral-output", "neutral.ms",
        "--progress"
        ]
    
    # run inside the target directory
    try:
        subprocess.run(args, cwd=out_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Simulation failed in {out_dir}: {e}")

def get_sfs(engine,species_id,model_id,pop_order,sfs_sample_individuals,chromosome,sfs_replicates,sel_2Ns=None):
    # build target run directory: current_path/data/species_folder_name/model_id/population/chromosome
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(pop_order), str(chromosome))
    os.makedirs(out_dir, exist_ok=True)
    args = [
        "simulator.py",
        "--engine", str(engine),
        "--species-id", str(species_id),
        "--model-id", str(model_id),
        "--pop-order", str(pop_order),
        "--sample-individuals", str(sfs_sample_individuals),
        "--chromosome", str(chromosome),
        "--replicates", str(sfs_replicates),
        "--sfs",
        "--threads", str(os.cpu_count()),
        ]
    
    # run inside the target directory
    try:
        subprocess.run(args, cwd=out_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Simulation failed in {out_dir}: {e}")
    

def get_info(model_id,pop_order,chromosome):
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(pop_order), str(chromosome))
    with open(f"{out_dir}/sweep.ms", 'r') as f:
        for i in range(4):
            line = f.readline()
            if i == 2:
                params = {k: __import__('ast').literal_eval(v.strip().capitalize() if v.strip().lower() in ('true','false') else v.strip()) for k, v in (item.split('=', 1) for item in line.replace("# params: ", "").strip().split(", ") if "=" in item)}
            if i == 3:
                selection = {k: __import__('ast').literal_eval(v.strip().capitalize() if v.strip().lower() in ('true','false') else v.strip()) for k, v in (item.split('=', 1) for item in line.replace("# selection: ", "").strip().split(", ") if "=" in item)}
    return params, selection

def ms2bin(model_id,pop_order,window,ips,iws,chromosome,type_,sweep_pos,length):
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(pop_order), str(chromosome))
    args = [
        "RAiSD-AI",
        "-n", "bin",
        "-I", f"{type_}.ms",
        "-w", str(window),
        "-ips", str(ips),
        "-iws", str(iws),
        "-its", str(sweep_pos),
        "-op", "IMG-GEN",
        "-icl", f"{type_}TR",
        "-bin",
        "-typ", "1",
        "-L", str(length),
        "-f"
    ]

    try:
        subprocess.run(args, cwd=out_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Simulation failed in {out_dir}: {e}")
    #os.remove(f"{out_dir}/{type_}.ms")
    #os.remove(f"{out_dir}/RAiSD_Info.bin.{type_}TR")


def train_model(model_id, pop_order, chromosome, epochs):
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(pop_order), str(chromosome))
    args = [
        "RAiSD-AI",
        "-n", "model",
        "-I", "RAiSD_Images.bin",
        "-f", "-op", "MDL-GEN",
        "-e", str(epochs),
        "-arc", "FASTER-NN"
    ]
    try:
        subprocess.run(args, cwd=out_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Simulation failed in {out_dir}: {e}")

species_dict = {sp.name: sp.id for sp in sps.all_species()}
if species in species_dict.values():
    species_full_name = [sp.name for sp in sps.all_species() if sp.id == species][0]
else:
    species_full_name = species
species_folder_name = "".join(c if c.isalnum() or c in "-_ " else "_" for c in str(species_full_name)).strip().replace(" ", "_")
os.makedirs(os.path.join("data", species_folder_name), exist_ok=True)

if species in species_dict.keys():
    species_id = species_dict[species]
elif species in species_dict.values():
    species_id = species
else:
    raise ValueError(f"Unknown species: {species}")

species_std = sps.get_species(species_id)

if not os.path.exists(f"data/{species_folder_name}/chromosomes.txt"):
        chromosomes = [c.id for c in species_std.genome.chromosomes]
        remove_chr = input(f"Which of {chromosomes} to remove (comma-separated): ")
        remove_chr = [c.strip() for c in remove_chr.split(",") if c.strip()]
        chromosomes = [c for c in chromosomes if c not in remove_chr]
        with open(f"data/{species_folder_name}/chromosomes.txt", 'w') as f:
            for chromosome in chromosomes:
                f.write(chromosome+"\n")
else:
    with open(f"data/{species_folder_name}/chromosomes.txt", 'r') as f:
        chromosomes = [line.strip() for line in f]

demographic_models = {
        m.id: [p.name for p in species_std.get_demographic_model(m.id).populations]
        for m in species_std.demographic_models}

for model_id, populations in demographic_models.items():
    for population in populations:
        for chromosome in chromosomes:
            #for engine in ["msms", "dicoal"]:
            run_simulation(engine, species_id, model_id, population, train_sample_individuals, window, ips, iws, chromosome, train_replicates, sel_2Ns)
            params, selection = get_info(model_id, population, chromosome)
            for type_ in ["sweep","neutral"]:
                ms2bin(model_id, population, window, ips, iws, chromosome, type_, selection["sweep_bp"],params["length"])
            train_model(model_id, population, chromosome, epochs)
            exit()

# RAiSD-AI -n test -mdl data/Homo_sapiens/OutOfAfricaExtendedNeandertalAdmixturePulse_3I21/YRI/1/RAiSD_Model.model -op SWP-SCN -I test.ms -L 10000000 -frm -G 300 -pci 1 1 -R
# RAiSD-AI -n ms -f -op RSD-DEF -I test.ms -frm -G 300 -R -L 1000000