<div align="center">

# RAiSD‑AI Web Scanner (with Pre‑Trained Demographic Models)

Primary goal: provide a ready‑to‑use web interface to scan user VCF/BCF files for selective sweeps using **pre‑trained RAiSD‑AI models** across published demographic models & populations (e.g. those under `data/Homo_sapiens/*`).

The repository also contains the internal tooling (simulator + training scripts) that produced those models. These scripts are **optional** for end users—they are only needed if you want to regenerate / extend the pretrained collection.

</div>

---

## 1. What You Get

| Layer | Purpose | Typical User Action |
|-------|---------|--------------------|
| Web app (`website_launcher.py`) | Upload VCF/BCF → auto demographic/pop match → RAiSD‑AI sweep scan → plots & gene context | ONLY part most users run |
| Pretrained models (`data/<Species>/<Model>/<Pop>/<Chrom>/RAiSD_Model.model/`) | Ready inference weights | Consume silently |
| Expected SFS (`data/<Species>/sfs.csv`) | Basis for model/pop matching via Jensen–Shannon distance | Already provided |
| Gene annotations (`data/<Species>/annotation/*.tsv`) | Gene overlap queries in UI | Already provided |
| (Optional) Training scripts (`3.sfs.py`, `2.annotation.py`, `1.train.py`) | Rebuild / extend model library | For advanced contributors |
| (Optional) `simulator.py` | Low‑level multi‑engine demographic simulator + SFS | For method development |

> If you cloned a release that already includes the populated `data/` tree, you can ignore every “Optional / Advanced” section below.

---

## 2. Web Scanner Feature Overview

| Area | Key Features |
|------|--------------|
| Simulation | Multi‑engine (`msprime`, `ms`, `msms`, `discoal`, `scrm`), sweep + neutral pairing, target SNP length estimation, growth discretization, optional RAM cap, SFS export, parallel chunking |
| SFS Builder (`3.sfs.py`) | Iterates all (model,pop), adaptive parallel reduction on OOM, skip tracking, normalized SFS matrix persisted to `sfs.csv` |
| Annotation (`2.annotation.py`) | Multi‑host BioMart fallback, parallel per‑chromosome TSV output, minimal columns |
| Training (`1.train.py`) | Sweep + neutral simulation, IMG generation (`RAiSD-AI -op IMG-GEN`), model training (`FASTER-NN`), transparent retry with `discoal` if primary engine fails |
| Web UI | Drag & drop upload, staged file auto‑indexing, SFS projection & model selection (JSD), RAiSD run + live metric/plot endpoints, gene lookup API |
| Robustness | Safe subprocess wrappers, skip registry, cleanup threads, resumable partial runs, memory throttling, automatic neutral pairing |
| Ops | Single script install, environment reproducibility (`raisd-ai.yml`), auto stale run deletion, lightweight logging |

---

## 3. Repository Layout (Perspective: End User vs. Advanced)

```
data/                     # Generated artefacts (SFS table, models, annotations, per-run outputs)
1.train.py                # RAiSD‑AI model training pipeline
2.annotation.py           # BioMart chromosome-wise annotation fetcher
3.sfs.py                  # Expected SFS builder across demographic models
simulator.py              # Core multi-engine simulator + SFS & ms/vcf output
website_launcher.py       # Flask web server (upload → model match → scan)
static/, templates/       # Web assets
install_dependencies.sh   # One-stop installer (creates/updates env + tools)
raisd-ai.yml              # Conda environment specification
```

### File Purpose Summary
| File | Role (Primary Scope) |
|------|----------------------|
| `website_launcher.py` | Flask web server; loads expected SFS + pretrained models, matches demographic/population, runs RAiSD‑AI, serves plots & gene/metric endpoints. |
| `simulator.py` | Low-level multi-engine demographic simulator (msprime/ms/msms/discoal/scrm); produces ms-like, VCF/BCF, and SFS; internal utility for training & SFS generation. |
| `3.sfs.py` | Iterates all stdpopsim demographic models/populations; runs simulator to build normalized expected SFS matrix (`sfs.csv`). Optional for end users (already provided). |
| `2.annotation.py` | Fetches per chromosome gene annotations via BioMart (multi-host fallback) producing TSVs for web gene overlap queries. Optional if already present. |
| `1.train.py` | Automates sweep + neutral simulation, RAiSD image generation, and RAiSD‑AI model training (`RAiSD_Model.model/`). Optional; only for regenerating pretrained models. |
| `install_dependencies.sh` | One-shot installer: creates/upgrades Conda env, builds RAiSD‑AI, msms, discoal, installs simulator wrapper. |
| `raisd-ai.yml` | Conda environment spec (Python + genomics + ML deps). |
| `static/`, `templates/` | Front-end JS/CSS + Jinja templates for the web UI. |
| `data/<Species>/sfs.csv` | Expected SFS table (rows `<model>=<population>`). |
| `data/<Species>/<Model>/<Pop>/<Chrom>/RAiSD_Model.model/` | Trained RAiSD‑AI weight directory. |

Two phases in practice:
1. **Preparation / Training (offline)** – run once per species (or when hyperparameters change).
2. **Scanning (online)** – launch the web server and analyze user data with already-trained models.

---

## 4. Installation

Prerequisites: A working Conda/Mamba (Linux x86_64 tested) and internet access.

```bash
bash install_dependencies.sh
conda activate raisd-ai
```

What this does:
1. Creates / updates env `raisd-ai` from `raisd-ai.yml`.
2. Downloads & wraps `msms` (Java), builds `discoal`.
3. Clones & compiles `RAiSD-AI` (both normal & ZLIB binaries).
4. Installs a convenience `simulator` wrapper into the env `bin/`.
5. (Hotfix) Applies a small `stdpopsim` patch when needed.

Verification:
```bash
which RAiSD-AI RAiSD-AI-ZLIB simulator msms discoal
python -c "import stdpopsim, msprime; print('stdpopsim OK')"
```

---

## 5. Quick Start (Using Pre‑Trained Models)

```bash
conda activate raisd-ai

# (Skip directly to launching if data/ is already populated)
python website_launcher.py
```

Navigate to http://127.0.0.1:5000/ and upload `.vcf`, `.vcf.gz`, or `.bcf`.

The app:
1. Infers total haploid sample size & ploidy.
2. Projects every expected SFS row to that sample size.
3. Computes Jensen–Shannon distance; selects the closest (model,pop).
4. Runs RAiSD‑AI with the corresponding pretrained model.
5. Streams plots & enables gene overlap queries.

---

## 6. OPTIONAL: Regenerating / Extending the Model Library
Only needed if you (a) add new species, (b) want different hyperparameters, or (c) update stdpopsim models.

Order:
1. `3.sfs.py` – build/update `sfs.csv`.
2. `2.annotation.py` – fetch per‑chromosome gene annotation TSVs.
3. `1.train.py` – simulate sweep + neutral, generate RAiSD images, train RAiSD‑AI model weights.

All steps rely on `simulator.py` internally.

### 6.1 Expected SFS (`3.sfs.py`)
...existing description...

### 6.2 Annotation (`2.annotation.py`)
...existing description...

### 6.3 Training (`1.train.py`)
...existing description...

### 6.4 Unified Simulator (`simulator.py`)
...existing description...

---

## 7. Web Application (`website_launcher.py`)
python website_launcher.py
```

Open: http://127.0.0.1:5000/ and upload a `.vcf`, `.vcf.gz`, or `.bcf`.

The UI will: infer ploidy → project SFS of every (model,pop) to the sample size → choose best via Jensen–Shannon distance → run RAiSD with the matched model → stream plots & gene hits.

---

## 6. Preparation / Training Pipeline Details

### 6.1 Step 1 – Expected SFS (`3.sfs.py`)
Generates a normalized site frequency spectrum (SFS) per (model,population). Adaptive retries lower `--parallel` if a memory threshold (`--max-ram-percent`) is exceeded. Skips are recorded in `data/<Species>/failed_parts.jsonl`.

Key flags:
* `--species <ID|Full Name>` (required)
* `--samples <haploid>` (default 10000)
* `--replicates <n>` (default 5)
* `--engine <msprime|ms|...>` (default `msprime` here)
* `--max-parallel <n>` (default: half CPU or ≤ replicates)
* `--max-ram-percent <pct>` (optional watchdog)
* `--chromosome <id>` (auto-chooses largest if omitted)

Output: `data/<Species>/sfs.csv` (rows labelled `<model>=<population>`; columns allele count bins excluding monomorphic).

### 6.2 Step 2 – Annotation (`2.annotation.py`)
Fetches chromosome-wise gene records via BioMart multi-host strategy.
* `--species <ID|Full Name>`
Produces: `data/<Species>/annotation/<chrom>.tsv` with columns `gene_id, transcript_id, gene_name, chromosome, start, end, biotype`.

### 6.3 Step 3 – Training (`1.train.py`)
For each demographic model, population, and chromosome (from `chromosomes.txt`):
1. Simulate sweep + neutral ms-like output via `simulator.py` (retry with `discoal` + `--sweep-time` if initial engine fails).
2. Transform `.ms` files into RAiSD image binaries (`RAiSD-AI -op IMG-GEN`).
3. Train RAiSD‑AI model (`FASTER-NN`) for configured `--epochs`.
4. Remove temporary image binaries; keep `RAiSD_Model.model/model.pt`.

Flags (excerpt):
* `--species` (required)
* `--train-sample-individuals` (default 100)
* `--train-replicates` (per class; default 1000) – reduce for quick smoke tests
* `--sel-s` selection coefficient (default 0.1)
* `--window --ips --iws` image shaping (window size & RAiSD image parameters)
* `--epochs` (default 3)
* `--engine` starting simulator (default `msms`)
* `--parallel` simulation worker processes (default half of available cores)
* `--gpu` attempt GPU‑accelerated mode (if RAiSD‑AI binary supports it)

First run prompts you to optionally remove chromosomes; it writes `data/<Species>/chromosomes.txt` for future runs.

### 6.4 Unified Simulator (`simulator.py`)
Core capabilities:
* Engines: `msprime`, `ms`, `msms`, `discoal`, `scrm`
* Sweep vs neutral paired runs (`--paired-neutral <neutral.ms>`)
* Automatic length estimation from `--target-snps` (with tolerance `--target-snps-tol`)
* SFS output (`--sfs <file>`, `--sfs-normalized`)
* Parallel chunking: `--parallel <workers>` × `--sims-per-work <n>`
* RAM guard: `--max-ram-percent <pct>`
* GPU flag pass‑through (`--gpu`) where supported

Representative call (like training does internally):
```bash
simulator.py \
  --engine msms \
  --species-id HomSap \
  --model-id OutOfAfrica_3G09 \
  --pop-order YRI \
  --sample-individuals 100 \
  --target-snps 800 \
  --chromosome 21 \
  --replicates 50 \
  --parallel 4 \
  --paired-neutral neutral.ms \
  --sfs sfs.sfs --sfs-normalized
```

---

## 7. Web Application (`website_launcher.py`)

Launch:
```bash
python website_launcher.py
```
Navigate to http://127.0.0.1:5000/

Flow:
1. Upload & stage file (progress bar; pre-index if needed with `bcftools`).
2. Infer sample size & ploidy (examining early GTs).
3. Load and project expected SFS for each (model,pop) → sample size.
4. Compute Jensen–Shannon distance; select best (and keep ranking).
5. Run RAiSD with matching model (naming report accordingly).
6. Stream metrics (PNG + JSON) and support later rescans.
7. Provide gene queries via start/end coordinates.

Public endpoints (selection):
`/` (upload), `/upload`, `/analyze`, `/runs/<id>/plot`, `/runs/<id>/metric_data`, `/runs/<id>/report_csv`, `/genes`, `/chromosomes`, `/healthz`.

Metrics: exposes composite `μ` and additional sweep tracking metrics (e.g. `sweepTR`).

Auto‑cleanup removes inactive run directories after `RESULTS_RETENTION_SECONDS`.

---

## 8. Data Layout

```
data/<Species>/
  chromosomes.txt                # Selected chromosomes (user curated first run)
  sfs.csv                        # Expected SFS (rows: <model>=<population>)
  failed_parts.jsonl             # Skip / failure registry (SFS + training)
  annotation/<chrom>.tsv         # Gene annotation files
  <Model>/<Population>/<Chrom>/
      sweep.ms / neutral.ms      # Source simulations (deleted later if desired)
      RAiSD_Model.model/model.pt # Trained RAiSD‑AI weights
runs/<uuid>/                     # Per upload: VCF/BCF, RAiSD reports, logs
```

---

## 9. Model Selection (SFS‑Driven)
Observed allele count distribution (after excluding monomorphic bins) is compared to every expected SFS (projected to the same sample size). The model/pop with minimal Jensen–Shannon distance is chosen; top alternatives are kept for transparency.

---

## 10. Performance & Resource Tips

| Issue | Mitigation |
|-------|------------|
| Slow training | Lower `--train-replicates`, `--epochs`, restrict chromosomes list |
| High RAM in SFS | Set `--max-ram-percent` & reduce `--max-parallel` |
| Upload sluggish | Pre‑compress & index (bgzip + bcftools index) |
| Sparse variants | Confirm chromosome choice & simulation window size; inspect RAiSD logs |
| Suspicious model match | Rebuild SFS with higher `--samples` or more `--replicates` |

---

## 11. Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| Web error: missing SFS | `sfs.csv` absent | Run `3.sfs.py` first |
| No genes returned | Annotation TSV missing | Run `2.annotation.py` |
| RAiSD report missing | Binary/path/empty input | Verify RAiSD install & variant density |
| Many skips in logs | Memory or zero SFS | Inspect `failed_parts.jsonl`; adjust flags |
| GPU ignored | No GPU or unsupported build | Check `nvidia-smi`; rebuild RAiSD‑AI |

---

## 12. Extending / Roadmap Ideas
* Add pytest smoke tests (mini chromosome, few replicates).
* Dockerfile / container recipe.
* REST endpoint to launch training jobs remotely.
* Model registry export (JSON) for integration with other pipelines.
* SFS caching keyed by (sample size, model,pop) projections.

Contributions welcome – open focused PRs. (Add a LICENSE first if you plan broad distribution.)

---

## 13. Security & Data Handling
* Uploaded files reside in `runs/<uuid>` and are auto‑removed after inactivity (configurable via env vars: `RESULTS_RETENTION_SECONDS`).
* No authentication layer – deploy behind VPN / reverse proxy for multi‑user or sensitive datasets.
* If using human data, ensure compliance with relevant privacy / consent policies; this toolkit does not anonymize content.

---

## 14. Limitations
* Single-host parallelism only (no distributed queue).
* Limited RAiSD‑AI hyperparameter surface (architecture fixed to `FASTER-NN`).
* No incremental retraining; always simulate anew (cache your raw ms if desired).
* Missing automated tests / CI.

---

## 15. License
No license file yet – add one (MIT / BSD‑3 / Apache‑2.0) to clarify permitted use.

---

## 16. Attribution / Citations
Please cite:
* RAiSD / RAiSD‑AI original publications.
* stdpopsim, msprime, demes for demographic simulations.
* BioMart for gene annotation services.

---

## 17. Quick Cheat Sheet
```bash
# SFS (fast test)
python 3.sfs.py --species HomSap --replicates 2 --samples 2000 --max-parallel 2

# Annotation
python 2.annotation.py --species HomSap

# Training (demo)
python 1.train.py --species HomSap --train-replicates 50 --epochs 1 --parallel 2

# Web app
python website_launcher.py
```

---

Happy scanning! 🧬

