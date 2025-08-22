species=HomSap
python3 1.sfs.py --species $species
python3 2.train_models.py --species $species
python3 3.annotation.py --species $species