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
RAiSD-AI Web Toolkit — Quick Start
=================================

What this is
------------
A compact toolkit to:
- generate expected site-frequency spectra (SFS) for stdpopsim demographic models,
- prepare training data and train RAiSD-AI models,
- run a Flask web UI to analyze user VCF/BCF files with precomputed models.

Prerequisites
-------------
- Linux/macOS with Python 3.9+ and conda recommended.
- A conda environment is provided: `raisd-ai.yml` (the repo expects a `raisd-ai` env in examples).
- External tools required on PATH for full use:
  - bcftools
  - RAiSD-AI (see: https://github.com/alachins/RAiSD-AI)

Quick install (recommended)
---------------------------
1. Create and activate the environment:

```bash
conda env create -f raisd-ai.yml -n raisd-ai
conda activate raisd-ai
```

2. Optionally run the helper script (may attempt to install binaries):

```bash
bash install_dependencies.sh
```

3. Ensure `simulator.py` and RAiSD executables are on PATH in the environment.

Basic workflows
---------------
1) Build SFS for a species (defaults to Homo sapiens):

```bash
python3 3.sfs.py --species "Homo sapiens"
```

This writes `data/<Species>/sfs.csv` containing expected SFS per demographic model/population.

2) Fetch gene annotations (Biomart):

```bash
python3 2.annotation.py --species "Homo sapiens"
```

Note: annotation TSVs should be under `data/<Species>/annotation/` for the web UI to use them.

3) Train models (simulation -> RAiSD-AI):

```bash
python3 1.train.py --species "Homo sapiens"
```

Trained models are saved under `data/<Species>/<model>/<pop>/<chr>/RAiSD_Model.model`.

4) Run the web UI (scan uploaded VCFs):

```bash
python3 website_launcher.py
```

Open http://0.0.0.0:5000/ and upload a VCF/BCF. The server will pick the best precomputed model (via projected SFS and JSD) and run RAiSD on the upload.

Notes & tips
------------
- RAiSD-AI: this project uses the RAiSD-AI training tool; see https://github.com/alachins/RAiSD-AI for the training binary and examples.
- `simulator.py` is the local simulation wrapper used by `1.train.py` and `3.sfs.py`. The scripts call `simulator.py` directly, so running inside the `raisd-ai` environment (which should expose `simulator.py` on PATH) is recommended.
- For quick tests, reduce `--replicates` and `--samples` when calling `3.sfs.py` or `1.train.py`.
- If the web UI cannot find annotations, verify `data/<Species>/annotation/*.tsv` exist (or edit `2.annotation.py` to change the output path).

Contributing
------------
- Add a LICENSE file if you plan to publish.
- Add a pinned requirements file or CI for reproducible installs.

License
-------
No license file included. Add one if you intend to re-publish.
