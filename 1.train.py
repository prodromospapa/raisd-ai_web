import os
import stdpopsim as sps
import subprocess
from tqdm import tqdm
import shutil

species = "Homo sapiens"
train_sample_individuals = 100
train_replicates = 1_000
sfs_sample_individuals = 10_000
sfs_replicates = 1_000
sel_s = 0.1

window = 500
ips = 11 # odd
iws = 10
epochs = 3
gpu = False

# default thread count (safe when os.cpu_count() returns None)
_cpu = os.cpu_count() or 2
parallel = max(1, _cpu // 2)

engine = "msms"


def run_subprocess(args, cwd, desc):
    """Run a subprocess, capture output, exit script on failure with detailed message."""
    try:
        return subprocess.run(
            args,
            cwd=cwd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        # print(
        #     f"\n===== SUBPROCESS ERROR: {desc} =====\n"
        #     f"Command: {' '.join(args)}\n"
        #     f"Return code: {e.returncode}\n"
        #     f"STDOUT:\n{e.stdout or '(empty)'}\n"
        #     f"STDERR:\n{e.stderr or '(empty)'}\n"
        #     f"===================================\n"
        # )
        raise SystemExit(1)


def run_simulation(engine,species_id,model_id,population,train_sample_individuals,window,ips,iws,chromosome,train_replicates,sel_s):
    # build target run directory: current_path/data/species_folder_name/model_id/population/chromosome
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(population), str(chromosome))
    os.makedirs(out_dir, exist_ok=True)
    args = [
        "simulator.py",
        "--engine", str(engine),
        "--species-id", str(species_id),
        "--model-id", str(model_id),
        "--pop-order", str(population),
        "--sample-individuals", str(train_sample_individuals),
        "--target-snps", str(int(window + ((ips/2) * iws)*1.1)),
        "--chromosome", str(chromosome),
        "--replicates", str(train_replicates),
        "--output", "sweep.ms",
        "--sweep-pos", "0.5",
        "--sel-s", str(sel_s),
        "--parallel", str(parallel),
        "--target-snps-tol", "0.1",
    "--paired-neutral", "neutral.ms",
        ]
    if gpu:
        args.append("--gpu")

    # run inside the target directory; if it fails, retry once using discoal and add --sweep-time 0.01
    try:
        run_subprocess(args, out_dir, f"simulation {model_id}/{population}/{chromosome}")
    except SystemExit:
        # prepare alternative args: switch engine to discoal and add sweep-time arg (only for this run)
        alt_args = args.copy()
        if "--engine" in alt_args:
            idx = alt_args.index("--engine")
            if idx + 1 < len(alt_args):
                alt_args[idx + 1] = "discoal"
        else:
            alt_args.extend(["--engine", "discoal"])

        if "--sweep-time" not in alt_args:
            alt_args.extend(["--sweep-time", "0.01"])

        #print(f"Engine '{engine}' failed for {model_id}/{population}/{chromosome}; retrying once with 'discoal' and --sweep-time 0.01")
        run_subprocess(alt_args, out_dir, f"simulation {model_id}/{population}/{chromosome} (discoal retry)")
        msg = f"Model '{model_id}' population '{population}'\n"
        species_dir = os.path.join(os.getcwd(), "data", species_folder_name)
        os.makedirs(species_dir, exist_ok=True)
        log_path = os.path.join(species_dir, "simulation.log")
        if os.path.exists(log_path):
            with open(log_path, "r", encoding="utf-8") as log_f:
                if msg not in log_f.read():
                    with open(log_path, "a", encoding="utf-8") as log_f:
                        log_f.write(msg)
        else:
            with open(log_path, "w", encoding="utf-8") as log_f:
                log_f.write(msg)

def get_info(model_id,population,chromosome):
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(population), str(chromosome))
    params = {}
    selection = {}
    try:
        with open(f"{out_dir}/sweep.ms", 'r') as f:
            for i in range(4):
                line = f.readline()
                if not line:
                    break
                if i == 2:
                    params = {k: __import__('ast').literal_eval(v.strip().capitalize() if v.strip().lower() in ('true','false') else v.strip()) for k, v in (item.split('=', 1) for item in line.replace("# params: ", "").strip().split(", ") if "=" in item)}
                if i == 3:
                    selection = {k: __import__('ast').literal_eval(v.strip().capitalize() if v.strip().lower() in ('true','false') else v.strip()) for k, v in (item.split('=', 1) for item in line.replace("# selection: ", "").strip().split(", ") if "=" in item)}
    except FileNotFoundError:
        print(f"Info file not found: {out_dir}/sweep.ms")

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

    run_subprocess(args, out_dir, f"RAiSD IMG-GEN {model_id}/{pop_order}/{chromosome}/{type_}")
    # remove temporary files only after success
    tmp_ms = f"{out_dir}/{type_}.ms"
    if os.path.exists(tmp_ms):
        os.remove(tmp_ms)
    info_file = f"{out_dir}/RAiSD_Info.bin.{type_}TR"
    if os.path.exists(info_file):
        os.remove(info_file)


def train_model(model_id, pop_order, chromosome, epochs):
    out_dir = os.path.join(os.getcwd(), "data", species_folder_name, str(model_id), str(pop_order), str(chromosome))
    args = [
        "RAiSD-AI",
        "-n", "model",
        "-I", "RAiSD_Images.bin",
        "-f", 
        "-op", "MDL-GEN",
        "-e", str(epochs),
        "-arc", "FASTER-NN"
    ]
    run_subprocess(args, out_dir, f"model training {model_id}/{pop_order}/{chromosome}")

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

# Progress accounting (models with varying population counts) -----------------
total_tasks = sum(len(pops) * len(chromosomes) for pops in demographic_models.values())

# Count already completed tasks (model file exists)
def _model_done(model_id, population, chromosome):
    return os.path.exists(f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/RAiSD_Model.model/model.pt")

tasks_done = sum(
    1
    for model_id, pops in demographic_models.items()
    for population in pops
    for chromosome in chromosomes
    if _model_done(model_id, population, chromosome)
)

with tqdm(total=total_tasks, initial=tasks_done, desc="total", unit="task") as total_bar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            for chromosome in chromosomes:
                if _model_done(model_id, population, chromosome):
                    total_bar.update(1)
                    total_bar.refresh()
                    # already counted in initial; just continue
                    continue
                try: 
                    sweep_path = f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/sweep.ms"
                    neutral_path = f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/neutral.ms"
                    if not os.path.exists(sweep_path) or not os.path.exists(neutral_path):
                        run_simulation(engine, species_id, model_id, population, train_sample_individuals, window, ips, iws, chromosome, train_replicates, sel_s)
                        
                    for type_ in ["sweep", "neutral"]:
                        if type_ == "sweep":
                            params, selection = get_info(model_id, population, chromosome)
                            sweep_bp = selection.get("sweep_bp")
                            length = params.get("length")
                            with open(f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/sweep_info.txt", 'w') as f:
                                f.write(f"sweep_bp: {sweep_bp}\nlength: {length}\n")
                        else:
                            with open(f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/sweep_info.txt", 'r') as f:
                                lines = f.readlines()
                                sweep_bp = int(lines[0].strip().split(": ")[1])
                                length = int(lines[1].strip().split(": ")[1])

                        img_info = f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/RAiSD_Images.bin/{type_}TR"
                        if not os.path.exists(img_info):
                            ms2bin(model_id, population, window, ips, iws, chromosome, type_, sweep_bp, length)
                        if type_ == "neutral":
                            os.remove(f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/sweep_info.txt")

                    train_model(model_id, population, chromosome, epochs)
                    img_bin = f"data/{species_folder_name}/{model_id}/{population}/{chromosome}/RAiSD_Images.bin"
                    if os.path.exists(img_bin):
                        shutil.rmtree(img_bin)
                finally:
                    # Only advance when task just completed (not pre-counted)
                    total_bar.update(1)
                    total_bar.refresh()

# RAiSD-AI -n test -mdl data/Homo_sapiens/OutOfAfricaExtendedNeandertalAdmixturePulse_3I21/YRI/1/RAiSD_Model.model -op SWP-SCN -I test.ms -L 10000000 -frm -G 300 -pci 1 1 -R
# RAiSD-AI -n ms -f -op RSD-DEF -I test.ms -frm -G 300 -R -L 1000000