import sys, io
import pandas as pd
import biomart
import stdpopsim
import concurrent.futures
from pathlib import Path
from tqdm import tqdm
import argparse

# Hard-code your inputs here
parser = argparse.ArgumentParser(description="Simulate allele frequency spectra.")
parser.add_argument(
    "--species",
    type=str,
    required=True,
    help="Species name (e.g., 'MusMus')."
)
args = parser.parse_args()
species = args.species

species_dict = {sp.name: sp.id for sp in stdpopsim.all_species()}
if species in species_dict.values():
    species = [sp.name for sp in stdpopsim.all_species() if sp.id == species][0]

HOSTS = [
    "http://www.ensembl.org/biomart",   # vertebrates
    "http://plants.ensembl.org/biomart",
    "http://fungi.ensembl.org/biomart",
    "http://metazoa.ensembl.org/biomart",
    "http://protists.ensembl.org/biomart",
    "http://bacteria.ensembl.org/biomart",
]

def base_name(sp):
    parts = sp.strip().replace("_"," ").split()
    return sp.lower() if len(parts)<2 else parts[0][0].lower()+''.join(c for c in parts[1].lower() if c.isalpha())

def find_mart(species):
    b = base_name(species)
    candidates = [f"{b}_gene_ensembl", f"{b}_eg_gene"]
    last_exc = None
    for host in HOSTS:
        try:
            s = biomart.BiomartServer(host)
            datasets = getattr(s, 'datasets', None)
            if not datasets:
                # no datasets listed for this host, try next
                continue
            for ds in candidates:
                if ds in datasets:
                    return s, ds
        except Exception as e:
            # record and continue to try other hosts (handles 404/connection errors)
            last_exc = e
            continue
    if last_exc:
        raise ValueError(f"No dataset found for '{species}'. Last error: {last_exc}")
    raise ValueError(f"No dataset found for '{species}'.")

def fetch(species, chrom):
    server, dataset = find_mart(species)
    mart = server.datasets[dataset]
    attrs = ["ensembl_gene_id","ensembl_transcript_id","external_gene_name",
             "chromosome_name","start_position","end_position","strand","gene_biotype"]
    r = mart.search({"filters":{"chromosome_name":str(chrom)}, "attributes":attrs})
    b = b"\n".join(list(r.iter_lines()))
    if not b.strip():
        raise ValueError("Empty response")
    txt = b.decode("utf-8","ignore")
    df = pd.read_csv(io.StringIO(txt), sep="\t", header=None, dtype=str, na_filter=False)
    df.columns = ["gene_id","transcript_id","gene_name","chromosome","start","end","strand","biotype"]
    for c in ("start","end","strand"):
        try: df[c]=pd.to_numeric(df[c])
        except: pass
    return df

if __name__=="__main__":
    species_std = stdpopsim.get_species(species)
    species_name = species_std.name
    annot_dir = f"data_/{species_name.replace(' ', '_')}/annotation"
    with open(f"data_/{species_name.replace(' ', '_')}/chromosomes.txt", 'r') as f:
            chromosomes = [line.strip() for line in f]
    # ensure annotation directory exists
    Path(annot_dir).mkdir(parents=True, exist_ok=True)


    def fetch_and_write(chrom_id: str):
        try:
            df = fetch(species_name, chrom_id)
            df = df.drop(columns=["strand"])
            df['biotype'] = df['biotype'].replace('_', ' ')
            out = Path(annot_dir) / f"{chrom_id}.tsv"
            df.to_csv(out, sep="\t", index=False)
            return chrom_id, None
        except Exception as e:
            return chrom_id, e

    max_workers = min(8, max(1, len(chromosomes)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as exe:
        futures = {exe.submit(fetch_and_write, cid): cid for cid in chromosomes}
        for fut in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Fetching annotations"):
            cid = futures.get(fut)
            try:
                chrom_id, err = fut.result()
                if err:
                    print(f"[ERROR] {chrom_id}: {err}", file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] {cid}: {e}", file=sys.stderr)
