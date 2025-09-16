import os
import stdpopsim as sps
import subprocess
from tqdm import tqdm
import shutil
import argparse
import time
import json

# default thread count (safe when os.cpu_count() returns None)
# CLI: accept species as stdpopsim id or full name (e.g. 'HomSap' or 'Homo sapiens')
parser = argparse.ArgumentParser(description="Train RAiSD-AI models for a species")
parser.add_argument("--species", type=str, required=True,
                    help="Species id (stdpopsim short id) or full name, e.g. 'HomSap' or 'Homo sapiens'.")
parser.add_argument("--train-sample-individuals", type=int, default=100,
                    help="Individuals per simulation (default 100)")
parser.add_argument("--train-replicates", type=int, default=1000,
                    help="Number of replicates per model (default 1000)")
parser.add_argument("--sel-s", type=float, default=0.1,
                    help="Selection coefficient used for sweep sims (default 0.1)")
parser.add_argument("--window", type=int, default=500,
                    help="Window size used to choose target SNPs (default 500)")
parser.add_argument("--ips", type=int, default=11,
                    help="ips parameter for RAiSD image generation (default 11)")
parser.add_argument("--iws", type=int, default=10,
                    help="iws parameter for RAiSD image generation (default 10)")
parser.add_argument("--epochs", type=int, default=3,
                    help="Training epochs for RAiSD-AI (default 3)")
parser.add_argument("--gpu", action='store_true', default=False,
                    help="Enable GPU mode for simulator/training when supported")
parser.add_argument("--parallel", type=int, default=max(1, os.cpu_count() // 2),
                    help=f"Parallel worker count for simulator runs (default {max(1, os.cpu_count() // 2)})")
parser.add_argument("--engine", type=str, default="msms",
                    help="Simulation engine to use (default: msms)")
args = parser.parse_args()

species = args.species
train_sample_individuals = args.train_sample_individuals
train_replicates = args.train_replicates
sel_s = args.sel_s

window = args.window
ips = args.ips
iws = args.iws
epochs = args.epochs
gpu = args.gpu
parallel = args.parallel
engine = args.engine


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
        print(
            f"\n===== SUBPROCESS ERROR: {desc} =====\n"
            f"Command: {' '.join(args)}\n"
            f"Return code: {e.returncode}\n"
            f"STDOUT:\n{e.stdout or '(empty)'}\n"
            f"STDERR:\n{e.stderr or '(empty)'}\n"
            f"===================================\n"
        )
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
        "--sweep-pos", "50",
        "--sel-s", str(sel_s),
        "--parallel", str(parallel),
        "--target-snps-tol", "10",
    "--paired-neutral", "neutral.ms",
        ]
    if gpu:
        args.append("--gpu")

    # run inside the target directory; if it fails, retry once using discoal and add --sweep-time (in generations)
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
            # Previous scripts passed 0.01 in 4N units; simulator.py now expects generations ago.
            # Assuming a typical N0 ~ 10,000, 0.01 (4N units) ≈ 0.01 * 4 * 10000 = 400 generations.
            # Use 400 generations as a conservative fallback for the retry.
            alt_args.extend(["--sweep-time", "400"])

        #print(f"Engine '{engine}' failed for {model_id}/{population}/{chromosome}; retrying once with 'discoal' and --sweep-time 0.01")
        run_subprocess(alt_args, out_dir, f"simulation {model_id}/{population}/{chromosome} (discoal retry)")
        # Report retry success via stdout. Do not write log files.
        print(f"Engine '{engine}' failed for {model_id}/{population}/{chromosome}; retry succeeded with 'discoal' and --sweep-time.")

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

# shared failed parts file (JSON-lines, used by 3.sfs.py)
skipped_file = os.path.join("data", species_folder_name, "failed_parts.jsonl")

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

def remove_from_sfs_csv_if_present(species_folder_name, key):
    """If data/<species>/sfs.csv exists, remove row matching `key` (index like 'model=pop')."""
    try:
        import pandas as _pd
    except Exception:
        return
    csv_path = os.path.join("data", species_folder_name, "sfs.csv")
    try:
        if os.path.exists(csv_path):
            df = _pd.read_csv(csv_path, index_col=0)
            if key in df.index:
                df = df.drop(index=key)
                df.to_csv(csv_path)
    except Exception:
        # best-effort: ignore any problems here
        pass

existing_skipped = load_skipped()

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

first = True
with tqdm(total=total_tasks, initial=tasks_done, desc="total", unit="task") as total_bar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            for chromosome in chromosomes:
                key = f"{model_id}={population}"
                # refresh skipped set on each iteration to catch skips written by other processes
                try:
                    latest_skipped = load_skipped()
                    # any newly recorded skipped keys should be removed from sfs.csv to keep bookkeeping consistent
                    new_skips = latest_skipped - existing_skipped
                    for nk in new_skips:
                        try:
                            remove_from_sfs_csv_if_present(species_folder_name, nk)
                        except Exception:
                            pass
                    # update the in-memory set
                    existing_skipped.update(latest_skipped)
                except Exception:
                    # if we fail to reload skipped file, fall back to previously loaded set
                    pass

                if key in existing_skipped:
                    total_bar.update(1)
                    total_bar.refresh()
                    # canonical skip list contains this key; no additional skip logging required
                    continue

                if _model_done(model_id, population, chromosome):
                    total_bar.update(1)
                    total_bar.refresh()
                    first = False
                    # already counted in initial; just continue
                    continue
                try: 
                    if first:
                        pass
                    first = False
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