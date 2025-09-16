import stdpopsim
import pandas as pd
import os
import subprocess
import json
from datetime import datetime
from tqdm import tqdm

species = "Homo sapiens"
engine ="msprime"
chromosome = "1"
samples = 10_000
replicates = 10
main_parallel = 5
max_ram_percent = 80
max_sims_per_work = 1 # can be None



species_dict = {sp.name: sp.id for sp in stdpopsim.all_species()}
if species in species_dict.values():
    species_full_name = [sp.name for sp in stdpopsim.all_species() if sp.id == species][0]
else:
    species_full_name = species
    species = species_dict[species]

species_folder_name = "".join(c if c.isalnum() or c in "-_ " else "_" for c in str(species_full_name)).strip().replace(" ", "_")
species_std = stdpopsim.get_species(species)
ploidy = species_std.ploidy

demographic_models = {
    m.id: [p.name for p in species_std.get_demographic_model(m.id).populations]
    for m in species_std.demographic_models
}

if os.path.exists(f"data/{species_folder_name}/sfs.csv"):
    sfs_csv = pd.read_csv(f"data/{species_folder_name}/sfs.csv",index_col=0)
    samples = (sfs_csv.shape[1] // ploidy)+1
else:
    os.makedirs(f"data/{species_folder_name}", exist_ok=True)
    sfs_csv = pd.DataFrame(columns= range(1,ploidy*(samples)))

# Load existing failed keys so we can skip them
failed_path = f"data/{species_folder_name}/failed_parts.jsonl"
failed_keys = set()
if os.path.exists(failed_path):
    try:
        with open(failed_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                    if "key" in obj:
                        failed_keys.add(obj["key"])
                except json.JSONDecodeError:
                    # skip malformed lines
                    continue
    except OSError:
        # If we can't read the file for some reason, proceed without skipping
        failed_keys = set()

# Count total runs (all model_id × population combinations)
total_runs = sum(len(pops) for pops in demographic_models.values())

with tqdm(total=total_runs, desc="Simulations", unit="run") as pbar:
    for model_id, populations in demographic_models.items():
        for population in populations:
            reason = "unknown"
            key = f"{model_id}={population}"
            if key in sfs_csv.index:
                # already done
                pbar.update(1)
                pbar.refresh()
                continue

            # Skip if this run previously failed and is recorded in failed_parts.jsonl
            if key in failed_keys:
                pbar.update(1)
                pbar.refresh()
                continue

            base_args = [
                "simulator.py",
                "--engine", engine,
                "--species-id", str(species),
                "--model-id", str(model_id),
                "--pop-order", str(population),
                "--sample-individuals", str(samples),
                "--chromosome", str(chromosome),
                "--replicates", str(replicates),
                "--parallel", str(main_parallel),
                "--sfs", "part.sfs",
                "--sfs-normalized",
                "--max-ram-percent", str(max_ram_percent),
            ]
            if max_sims_per_work is not None:
                base_args += ["--sims-per-work", str(max_sims_per_work)]
            first = True
            done = False
            parallel = main_parallel
            while parallel != 0:
                proc = None
                stderr = ""
                try:
                    # Start the simulator in a new session so it has its own process group.
                    # Use Popen so we can terminate the entire group if the user cancels.
                    proc = subprocess.Popen(
                        base_args,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.PIPE,
                        start_new_session=True,
                        text=True
                    )

                    # Wait for completion and capture stderr
                    _, stderr = proc.communicate()
                    rc = proc.returncode

                    # If child exit code indicates it was killed by a signal,
                    # treat it as an out-of-memory / termination and retry with
                    # reduced resources (matching previous behavior).
                    if rc != 0:
                        # Normalize to signal number when possible.
                        sig = None
                        try:
                            if rc < 0:
                                sig = -rc
                            elif rc > 128:
                                sig = rc - 128
                        except Exception:
                            sig = None

                        # If terminated by SIGKILL or SIGTERM, treat as RAM problem
                        if sig in (9, 15):
                            if first:
                                base_args += ["--sims-per-work", "1"]
                                first = False
                                continue
                            else:
                                parallel -= 1
                                base_args[base_args.index("--parallel") + 1] = str(parallel)
                                continue
                        if rc == 2:
                            if first:
                                base_args += ["--sims-per-work", "1"]
                                first = False
                                continue
                            else:
                                if parallel >= 5:
                                    parallel -= 1
                                else:
                                    parallel = parallel // 2
                                base_args[base_args.index("--parallel") + 1] = str(parallel)
                                continue
                        # Non-recoverable non-zero exit: raise to be handled below
                        raise subprocess.CalledProcessError(rc, base_args, stderr=stderr)

                    sfs_part = pd.read_csv(
                        "part.sfs",
                        index_col=0,
                        skiprows=2,
                        sep="\t",
                        header=0
                    )
                    os.remove("part.sfs")
                    values = sfs_part.loc["mean"].values
                    if sum(values) == 0:
                        reason = "all_zero_sfs"
                        break
                    sfs_csv.loc[f"{model_id}={population}"] = values
                    sfs_csv.to_csv(f"data/{species_folder_name}/sfs.csv")
                    done = True
                    break

                except KeyboardInterrupt:
                    # If the user cancels, ensure we kill the child process group
                    if proc is not None:
                        try:
                            # Kill the whole process group started by Popen (negative pid)
                            os.killpg(proc.pid, 9)
                        except Exception:
                            try:
                                proc.kill()
                            except Exception:
                                pass
                    # Re-raise to allow higher-level handlers (if any) or to stop execution
                    raise

                except subprocess.CalledProcessError as e:
                    rc = getattr(e, 'returncode', None)
                    # In our raise above we didn't include returncode, so infer from proc if possible
                    if rc is None and proc is not None:
                        rc = proc.returncode

                    if rc == 2:
                        if first:
                            base_args += ["--sims-per-work", "1"]
                            first = False
                        else:
                            parallel -= 1
                            base_args[base_args.index("--parallel") + 1] = str(parallel)
                        continue
                    # For other return codes, record and break
                    reason = f'subprocess_error "{rc}"'
                    try:
                        if e.stderr:
                            reason += f': {e.stderr.strip()[:200]}'
                    except Exception:
                        # fallback to stderr captured from proc.communicate
                        try:
                            if proc is not None and proc.stderr is not None:
                                # proc.stderr may be a pipe; we captured it above
                                if stderr:
                                    reason += f': {stderr.strip()[:200]}'
                        except Exception:
                            pass
                    break

                except ValueError as e:
                    if "non-sampling population" in str(e):
                        reason = "non_sampling_population"
                        break

                finally:
                    # Ensure no orphaned child process remains if we decided to stop
                    if proc is not None and proc.poll() is None and (not done):
                        try:
                            os.killpg(proc.pid, 9)
                        except Exception:
                            try:
                                proc.kill()
                            except Exception:
                                pass
            if parallel == 0:
                reason = "run_out_of_ram"

            if not done:
                # Record failure for this model/population
                failed_dir = f"data/{species_folder_name}"
                os.makedirs(failed_dir, exist_ok=True)
                failed_path = os.path.join(failed_dir, "failed_parts.jsonl")
                record = {
                    "key": f"{model_id}={population}",
                    "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "reason": reason,
                    "source": "3.sfs"
                }
                # Append as a JSONL line
                with open(failed_path, "a", encoding="utf-8") as fh:
                    fh.write(json.dumps(record, ensure_ascii=False) + "\n")
            pbar.update(1)
            pbar.refresh()