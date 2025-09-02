export species=HomSap
species_full=$(python3 - <<'PY'
import os, sys
import stdpopsim
name = os.environ.get('species')
if not name:
    sys.exit(1)
try:
    sp = stdpopsim.get_species(name)
except Exception:
    print(name.replace(" ", "_"))
    sys.exit(0)
for attr in ('full_name','display_name','latin_name','scientific_name','name','common_name'):
    if hasattr(sp, attr):
        val = getattr(sp, attr)
        if val:
            print(val.replace(" ", "_"))
            break
else:
    print(str(sp).replace(" ", "_"))
PY
)
export species_full
if [ ! -f "data_/$species_full/sfs.csv" ]; then
    python3 1.sfs.py --species $species
fi
python3 2.train_models.py --species $species
python3 3.annotation.py --species $species