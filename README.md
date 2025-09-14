<div align="center">

# RAiSD‑AI Web Toolkit

End‑to‑end workflow for demographic simulation, RAiSD‑AI model training, and an interactive web interface to scan user VCF/BCF genomes for selective sweeps.

**Simulate → Build SFS → Annotate Genes → Train RAiSD‑AI → Upload VCF → Automatic Model Match → RAiSD Scan → Interactive Plots & Gene Context**

</div>

---

## 1. Highlights

| Area | What You Get |
|------|--------------|
| Demographic simulation | Unified `simulator.py` wrapper over ms/msms/discoal/scrm/msprime with adaptive length, resource guards, paired neutral runs, SFS extraction |
| Expected SFS table | Automated per model/population SFS generation (`3.sfs.py`) with skip‑tracking + resume |
| Gene annotations | Parallel BioMart fetcher (`2.annotation.py`) writing per‑chromosome TSVs |
| Training | Automated sweep + neutral simulation, image generation, model training (`1.train.py`) |
| Model store layout | `data/<Species>/<ModelID>/<Population>/<Chromosome>/RAiSD_Model.model/` |
| Web app | Flask UI: upload VCF/BCF, infer sample ploidy, project expected SFS, pick best `(model,pop)` via Jensen–Shannon distance, run RAiSD, visualize metrics, browse genes |
| Metrics & plots | Dynamic PNG endpoint, JSON metric streaming, CSV export, RAiSD metric selection (`μ`, `sweepTR`) |
| House‑keeping | Auto cleanup of stale runs, resumable jobs, rescan endpoint |
| Robustness | Memory capping, adaptive parallelism, retry / fallback engines, neutral pairing, detailed logs |

---

## 2. Repository Overview

```
data/                         # Generated artefacts (SFS, models, annotations)
1.train.py                    # Train RAiSD‑AI models (simulation + image + NN training)
2.annotation.py               # Fetch & store per-chromosome gene annotations
3.sfs.py                      # Build expected SFS table across demographic models
simulator.py                  # Unified demographic simulator & SFS producer
website_launcher.py           # Flask application (upload → model match → RAiSD scan)
static/, templates/           # Web assets & Jinja templates
install_dependencies.sh       # Helper bootstrap script
raisd-ai.yml                  # Conda environment definition
```

Two major phases:

1. **Training / Preparation Pipeline**: (steps 1–4 below) produce `sfs.csv`, annotation TSVs, and trained RAiSD‑AI models.
2. **Web Scanning**: run Flask server, upload user data, automatic demographic/pop match, sweep scan + visualization.

---

## 3. Installation

### 3.1 Prerequisites

Required:
* Python 3.9+ (3.11+ ok) with `conda` or `mamba` recommended
* External tools on PATH:
  * `bcftools` (indexing & allele counts)
  * One of `RAiSD-AI-ZLIB`, `RAiSD-AI`, or fallback `RAiSD` executables
* Internet access (only for first-time annotation & BioMart queries)

Core Python packages used (non‑exhaustive): `stdpopsim`, `msprime`, `demes`, `numpy`, `pandas`, `flask`, `pysam`, `tqdm`, `psutil`, `biomart`, `matplotlib`.

### 3.2 Create environment

```bash
conda env create -f raisd-ai.yml -n raisd-ai
conda activate raisd-ai
```

Optional helper (installs extra/binaries if scripted):

```bash
bash install_dependencies.sh
```

### 3.3 Verify toolchain

```bash
which bcftools
which RAiSD-AI || which RAiSD-AI-ZLIB || which RAiSD
python -c "import stdpopsim, msprime; print('stdpopsim ok')"
```

---

## 4. Quick Start (Happy Path)

```bash
# 1. Build expected SFS table
python 3.sfs.py --species "Homo sapiens"

# 2. Fetch gene annotations
python 2.annotation.py --species "Homo sapiens"

# 3. Train RAiSD‑AI models (sweeps + neutral)
python 1.train.py --species "Homo sapiens"

# 4. Launch web interface
python website_launcher.py

# Open the UI
xdg-open http://127.0.0.1:5000/
```

Upload a `.vcf`, `.vcf.gz`, or `.bcf` — the app auto:
1. Infers sample count + ploidy
2. Projects each precomputed model/pop SFS to that sample size
3. Computes Jensen–Shannon distance; selects best `(model,pop)`
4. Runs RAiSD with the corresponding trained model
5. Streams metrics & plots; allows rescan / CSV / gene queries.

---

## 5. Data & Directory Layout

Generated structure (per species):

```
data/<Species>/
  chromosomes.txt                   # User‑confirmed chromosome list (first run of training)
  sfs.csv                           # Expected (normalized or raw) SFS rows: <model>=<population>
  skipped_demographics.jsonl        # Models/pops skipped (memory / zero SFS / failure)
  <ModelID>/<Population>/<Chromosome>/
      sweep.ms / neutral.ms         # Raw simulations retained until conversion
      RAiSD_Model.model/            # Trained model (contains model.pt)
```

Annotation TSVs: `data/<Species>/annotation/<chrom>.tsv`

Web run artefacts: `runs/<uuid>/` containing uploaded file, `RAiSD_Report.*`, logs, and final assets (auto‑cleaned).

---

## 6. Core Scripts & CLI Reference

### 6.1 `3.sfs.py` – Build Expected SFS Table

Generates mean SFS per (`demographic_model`, `population`). Adaptive retries lowering parallelism if memory is exceeded; records skipped items.

Key arguments:
* `--species <ID|Full Name>` (required)
* `--samples <haploid-samples>` default 10000 (used to size SFS)
* `--replicates <n>` default 5 per model/pop
* `--engine <scrm|ms|...>` default `scrm` for SFS speed
* `--max-ram-percent <pct>` optional watchdog
* `--max-parallel <n>` upper bound parallel workers
* `--chromosome <id>` choose contig (default 21 if present)

Output: `data/<Species>/sfs.csv`

### 6.2 `2.annotation.py` – Gene Annotation Fetch

Fetches per‑chromosome records via BioMart (multi-host fallback). Writes chromosome TSVs used by the web gene window lookups.

* `--species <ID|Full Name>` (required)

Outputs: `data/<Species>/annotation/<chrom>.tsv`

### 6.3 `1.train.py` – RAiSD‑AI Model Training

Workflow per (model, population, chromosome):
1. Simulate sweep & paired neutral (`simulator.py`); fallback engine retry (e.g., to `discoal`).
2. Convert ms output → RAiSD image binaries (`RAiSD-AI -op IMG-GEN`).
3. Train neural architecture (`FASTER-NN`, epochs configurable) producing `RAiSD_Model.model/model.pt`.
4. Cleanup transient image binaries.

Key arguments:
* `--species <ID|Full Name>` required
* `--train-sample-individuals` per simulation (default 100)
* `--train-replicates` replicates per class (default 1000)
* `--sel-s <s>` selection coefficient for sweep class (default 0.1)
* `--window, --ips, --iws` RAiSD image params
* `--epochs <n>` (default 3) – increase for quality
* `--engine <msms|...>` default `msms`
* `--parallel <n>` CPU workers for simulation
* `--gpu` attempt GPU where supported

### 6.4 `simulator.py` – Unified Simulation Engine

Supports engines: `discoal`, `ms`, `msms`, `scrm`, `msprime` + SFS computation & VCF/ms output.

Features:
* Adaptive length estimation from `--target-snps`
* Sweep modeling (origin / fixation) via engine‑specific flags
* Growth discretization for exponential epochs (optional)
* Paired neutral runs (`--paired-neutral`)
* SFS calculation (mean or per-replicate) with normalization
* Memory guard (`--max-ram-percent`) & per-worker thread throttling
* Chunked parallel replication with `--parallel` & `--sims-per-work`

Invoke indirectly from training / SFS scripts, or directly for custom experiments.

---

## 7. Web Application

Start with:

```bash
python website_launcher.py
```

Open: http://127.0.0.1:5000/

### 7.1 Flow
1. Upload / select staged VCF/BCF
2. Index file (bcftools) if needed
3. Infer ploidy + sample count
4. Load expected SFS (`sfs.csv`) → normalize → project to sample size
5. Compute JSD vs observed SFS; choose best match
6. Locate trained model path; run RAiSD
7. Serve interactive plots + gene queries

### 7.2 Key Endpoints (public)

| Method | Route | Purpose |
|--------|-------|---------|
| GET | `/` | Upload form & species selector |
| POST | `/analyze` | Begin or resume analysis for an upload |
| GET | `/runs/<id>/plot` | Dynamic PNG plot of a metric (query: metric, xmin, xmax) |
| GET | `/runs/<id>/metric_data` | JSON metric series (downsampling) |
| GET | `/runs/<id>/report_csv` | Download parsed RAiSD metrics CSV |
| POST | `/runs/<id>/rescan` | Re-run RAiSD with different grid / filters |
| GET | `/runs/<id>/rescan_status` | Poll rescan job state |
| POST | `/runs/<id>/cleanup` | Schedule early deletion of run directory |
| GET | `/runs/<id>/expected_sfs` | JSON: selected model/pop expected vs observed SFS |
| GET | `/runs/<id>/final` | Finalization status (ready signal) |
| POST | `/runs/<id>/heartbeat` | Keep run alive (avoid auto cleanup) |
| GET | `/chromosomes` | List available chromosomes (species context) |
| GET | `/genes` | Query overlapping genes (params: species, chr, start, end[, biotype]) |
| POST | `/upload` | Low-level upload (AJAX chunked forms) |
| GET | `/healthz` | Liveness probe |

### 7.3 Gene Query Output

Returns overlapping genes sorted by decreasing overlap bp with fields: `label,start,end,biotype,overlap_bp`.

### 7.4 Metrics

Displayed metrics derived from RAiSD report. Primary exposed: `μ` (composite) and `sweepTR` (probability/tracking metric). JSON & PNG endpoints accept sanitized metric names (underscores instead of LaTeX).

---

## 8. Model Selection Logic

For each uploaded dataset:
1. Compute observed allele count distribution via `bcftools +fill-tags`.
2. Normalize counts to polymorphic SFS (exclude bins 0 & n).
3. For every row `<model>=<population>` in projected expected SFS, compute Jensen–Shannon distance.
4. Select minimal JSD; record top matches (first 6) plus full ranking.

Persisted meta allows rescans to name RAiSD reports as `RAiSD_Report.<Model>__<Population>` for clarity.

---

## 9. Performance & Resource Tips

| Scenario | Tip |
|----------|-----|
| Training too slow | Lower `--train-replicates`, `--epochs`; start with subset chromosomes |
| Memory errors in SFS build | Use `--max-ram-percent` and reduce `--max-parallel` |
| Large VCF upload stalls | Compress & index (`bgzip` + `bcftools index`) before upload |
| Few variants in scan | Ensure chromosome choice & model length align with data; check RAiSD logs |
| Wrong model chosen | Inspect top matches JSD values; consider regenerating SFS with larger `--samples` |

---

## 10. Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| `Missing SFS table` in web | Did not run `3.sfs.py` | Build SFS first |
| Annotation lookup fails | TSVs absent or species folder mismatch | Re-run `2.annotation.py`; verify folder name matches species selector |
| RAiSD report not produced | RAiSD binary mismatch / flags / empty region | Check `raisd_stderr.log`; ensure correct RAiSD-AI build |
| Skipped demographics in training/SFS | Memory or zero SFS | Adjust parallelism; inspect `skipped_demographics.jsonl` |
| Very slow msprime runs | Hyperthread over-subscription | Reduce `--parallel`; ensure BLAS thread env vars capped |
| GPU flag ignored | No GPU detected | Confirm `nvidia-smi` output; rebuild RAiSD with GPU support |

---

## 11. Extending & Contributing

Suggestions:
* Add CI (GitHub Actions) to build SFS & run a smoke scan.
* Provide Dockerfile for reproducible web deployment.
* Implement caching for SFS projection for repeated sample sizes.
* Add REST/JSON endpoint to trigger training jobs remotely.

To contribute: open PRs with focused commits; include brief usage notes. **Add a LICENSE** before external sharing.

---

## 12. Security & Data Handling

* Uploaded VCF/BCF files are stored under `runs/<uuid>/` and auto‑deleted after inactivity (`RESULTS_RETENTION_SECONDS`).
* No authentication layer is included; deploy behind a protected network or add auth (reverse proxy / Flask login) for multi‑user environments.
* Ensure PII stripping if using human data (the tool processes genotypes; privacy policies may apply).

---

## 13. Limitations & Future Ideas

* No built‑in distributed job queue (single host thread pool only)
* Single species processed per run instance (multi‑species fine, but one server process)
* Training hyperparameters minimal (epochs, architecture fixed to FASTER-NN)
* Lack of formal tests; consider adding pytest harness (simulate tiny model / run web in test mode)

---

## 14. License

No license file currently. Add e.g. MIT, BSD‑3, or Apache‑2.0 to clarify permitted use.

---

## 15. Attribution

* RAiSD / RAiSD‑AI: upstream algorithms and binaries.
* stdpopsim / msprime / demes: demographic simulation foundation.
* BioMart: gene annotation service.

If you publish results produced with this stack, cite the original RAiSD and stdpopsim works accordingly.

---

## 16. Quick Reference Cheat Sheet

```bash
# SFS table (fast test run)
python 3.sfs.py --species HomSap --replicates 2 --samples 2000 --max-parallel 2

# Annotation
python 2.annotation.py --species HomSap

# Training (lightweight demo)
python 1.train.py --species HomSap --train-replicates 50 --epochs 1 --parallel 2

# Launch web
python website_launcher.py
```

---

Happy scanning! 🧬

