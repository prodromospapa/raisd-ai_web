RAiSD-AI Web Toolkit
=====================

Overview
--------

This repository contains two related but distinct parts:

- Training & simulation (data preparation)
   - scripts that run population-genetic simulations, compute expected site-frequency spectra (SFS) for many stdpopsim demographic models, fetch gene annotations, and train RAiSD-AI models from simulated sweep/neutral examples.

- Website & scanning (user-facing analysis)
   - a Flask-based web application that accepts user VCF/BCF uploads, matches the uploaded data to the best precomputed demographic model using projected SFS and JSD, runs RAiSD scans with the selected model, and presents interactive plots and metadata.

This README documents repository structure, installation pointers, and common workflows based only on the files present in this repo.

Table of contents
-----------------

- Overview
- Requirements
- Installation
- Repository layout (grouped by the two parts)
- Common workflows
   - Build SFS table (training)
   - Fetch gene annotations (training)
   - Train RAiSD-AI models (training)
   - Launch the web UI and analyze an uploaded VCF (website)
- API endpoints (summary)
- Troubleshooting & tips
- Contributing
- License

Requirements
------------

The repository's scripts call a mix of Python libraries and external executables. Based on imports and subprocess calls in the repository, you will likely need:

- Python 3.9+ (recommended)
- Python packages (examples): pandas, numpy, matplotlib, flask, werkzeug, pysam, stdpopsim, msprime, demes, tqdm, psutil, biomart
- External executables:
   - `bcftools` (indexing, conversion, SFS extraction)
   - `RAiSD-AI` or `RAiSD-AI-ZLIB` / `RAiSD` (used for RAiSD scans and RAiSD-AI training)

Note: There is no pinned `requirements.txt` in this repo; consider creating one or using `install_dependencies.sh` and a reproducible environment (conda/venv).

Installation
------------

This project was developed to run inside a preconfigured conda environment named `raisd-ai` in this repository. There is an environment spec file named `raisd-ai.yml` and an install script `install_dependencies.sh` that the original developer used to install packages into the `raisd-ai` environment.

Suggested (developer) steps:

1. Create / activate the provided conda environment (example):

   conda env create -f raisd-ai.yml -n raisd-ai
   conda activate raisd-ai

   Or run the project's install script inside an existing conda base shell:

   bash install_dependencies.sh

2. Ensure external executables are available in that environment's PATH:
   - `bcftools`
   - `RAiSD-AI` or `RAiSD-AI-ZLIB` / `RAiSD`

3. The developer's environment places `simulator.py` on PATH (so top-level scripts invoke `simulator.py` directly). Ensure the working environment where you run `1.train.py` and `3.sfs.py` has `simulator.py` discoverable on PATH (for example by running from the repository root with the `raisd-ai` environment active).

4. Prepare data directories (the scripts will create missing directories as needed):
   - `data/` should contain species directories (e.g., `Homo_sapiens`) with `chromosomes.txt` and `sfs.csv` where applicable.
   - `runs/` will be created automatically and stores per-analysis run data.

Repository layout (grouped by functionality)
-------------------------------------------

Training & simulation (data preparation)
- `simulator.py` — simulation wrapper integrating several engines (ms/msms/discoal/scrm/msprime/stdpopsim) and producing `.ms` outputs and optional `sfs.sfs` files.
- `3.sfs.py` — driver that enumerates stdpopsim demographic models and runs `simulator.py` to build `data/<species>/sfs.csv` (expected SFS per demographic model/population).
- `2.annotation.py` — fetches gene/transcript annotations via Biomart and writes per-chromosome TSVs (note: script currently writes to `data_/...` while the web UI reads `data/...`; see Troubleshooting).
- `1.train.py` — orchestrates producing sweep/neutral `.ms` files, converting them to RAiSD images, and invoking `RAiSD-AI` to train models (produces `RAiSD_Model.model`).

Website & scanning (user-facing analysis)
- `website_launcher.py` — Flask app that handles uploads/staging, SFS-based best-match selection, RAiSD runs using a selected precomputed model, dynamic plotting, rescans, gene lookups, and run lifecycle management. Uses `runs/` for per-analysis storage.
- `templates/` — HTML templates used by Flask (`index.html`, `processing.html`, `result.html`).
- `static/` — JavaScript and CSS for the UI.

Common workflows
----------------

Training / data preparation (first part)

1. Build the expected SFS table for a species

   Run `3.sfs.py` to generate SFS for all stdpopsim demographic models for the chosen species. This produces `data/<Species>/sfs.csv` (expected SFS per demographic model/population).

   python3 3.sfs.py

2. Fetch gene annotations

   Use `2.annotation.py` to fetch gene/transcript annotations via Biomart and write per-chromosome TSVs. By default this script writes into `data_/.../annotation/` — if you want the web UI to load these annotations without moving files, update the script to write into `data/` instead.

   python3 2.annotation.py --species "Homo sapiens"

3. Train RAiSD-AI models

   Use `1.train.py` to run simulations that produce sweep and neutral examples, convert them for RAiSD, and run `RAiSD-AI` to train models. Trained models are expected at `data/<Species>/<dm>/<pop>/<chr>/RAiSD_Model.model`.

   python3 1.train.py

Website / scanning (second part)

1. Start the Flask app:

   python3 website_launcher.py

2. In the web UI (default http://0.0.0.0:5000/): upload a VCF/BCF and request analysis. The server will:
  - Stage the upload in `runs/<run_id>/` and optionally convert BCF -> VCF.GZ.
  - Infer sample size and ploidy from the VCF and compute the observed SFS (via `bcftools`).
  - Load `data/<species>/sfs.csv`, project expected SFS to the observed sample size, compute JSD across models/pops, and pick the best match.
  - Find the RAiSD model directory `data/<species>/<dm>/<pop>/<chr>/RAiSD_Model.model` and run RAiSD on the uploaded VCF with that model.
  - Save the RAiSD report and logs under `runs/<run_id>/` and serve dynamic metric plots and CSV exports.

API endpoints (summary)
-----------------------

The Flask app exposes endpoints (summary only). See `website_launcher.py` for full details and parameter descriptions.

- GET `/` — main index UI
- POST `/upload` — stage an uploaded VCF/BCF and return `run_id` + contigs
- POST `/analyze` — start analysis (upload or staged) and return run id + polling URLs
- GET `/runs/<run_id>/tail` — tail logs for the run
- GET `/runs/<run_id>/plot?metric=...` — dynamic PNG plot for a metric
- GET `/runs/<run_id>/metric_data` — JSON series for a metric
- GET `/runs/<run_id>/report_csv` — parsed RAiSD report as CSV
- POST `/runs/<run_id>/rescan` — rescan using a specified demographic model/population
- GET `/runs/<run_id>/rescan_status?job=<id>` — check rescan job status
- GET `/genes?species=...&chromosome=...&start=...&end=...` — return genes overlapping a window (requires annotation TSVs)

Troubleshooting & tips
----------------------

- If RAiSD fails, make sure `RAiSD-AI` (or `RAiSD-AI-ZLIB` / `RAiSD`) is installed and on PATH, or set `RAISD_AI`/`RAISD_AI_ZLIB`.
- Ensure `bcftools` is installed and available (used for indexing, conversion, and SFS extraction).
- Annotation location: `2.annotation.py` currently writes into `data_/.../annotation/` while the web UI reads annotation from `data/.../annotation/`. Move annotation TSVs or edit `2.annotation.py` to write into `data/` for the website to find them automatically.
- If indexing fails for large files, provide pre-indexed files (.tbi/.csi) or ensure the server has enough memory/disk.
- For heavy workloads (training or large RAiSD scans) use a machine with sufficient CPU and memory and run training scripts separately from the web server.

Contributing
------------

Contributions are welcome. Suggested small improvements:

- Add a pinned `requirements.txt` or `environment.yml` for reproducible environments.
- Add tests and a CI workflow.
- Reconcile `data/` vs `data_/` usage for annotations to avoid confusion (recommended).

License
-------

This repository does not contain an explicit LICENSE file. When publishing to GitHub, add an appropriate license file (for example, MIT or BSD) if you want to allow reuse.

Acknowledgements
----------------

This project integrates stdpopsim, RAiSD, and RAiSD-AI concepts for simulation-driven scan workflows. See in-repo scripts for the detailed logic and helpers used to glue the pipeline together.