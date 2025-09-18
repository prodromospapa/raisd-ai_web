import sys, io
import pandas as pd
import biomart
import stdpopsim
import concurrent.futures
from pathlib import Path
from tqdm import tqdm
import argparse
import os

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
    species_full_name = [sp.name for sp in stdpopsim.all_species() if sp.id == species][0]
else:
    species_full_name = species
    species = species_dict[species]

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
    # Write annotation files into the main `data/` tree so the website finds them.
    base_annot_dir = f"data/{species_name.replace(' ', '_')}/annotation"
    with open(f"data/{species_name.replace(' ', '_')}/chromosomes.txt", 'r') as f:
        chromosomes = [line.strip() for line in f]
    # ensure base annotation directory exists
    Path(base_annot_dir).mkdir(parents=True, exist_ok=True)
    # Auto-discover Ensembl releases from the public FTP and process each
    import urllib.request, urllib.error, urllib.parse, re, gzip, shutil, tempfile, json

    def discover_releases(base_urls=("https://ftp.ensembl.org/pub/", "http://ftp.ensembl.org/pub/")):
        found = set()
        rx = re.compile(r"release[-_/]?(\d+)|release[\-]?(\\d+)|release_(\d+)")
        # Try base URLs and parse HTML for release patterns
        for base in base_urls:
            try:
                with urllib.request.urlopen(base, timeout=20) as r:
                    txt = r.read().decode('utf-8', 'ignore')
                for m in re.finditer(r'release[-_/]?(\d+)', txt, flags=re.IGNORECASE):
                    found.add(m.group(1))
            except Exception:
                continue
        # Return sorted numeric releases (as strings)
        return sorted(found, key=lambda x: int(x))

    def process_release(release: str):
        """Download GTF for this species for a given release and write per-chromosome TSVs."""
        ver_dir = Path(base_annot_dir) / f"ensembl_{release}"
        # If directory already exists and has TSVs, skip (assume complete)
        existing = list(ver_dir.glob("*.tsv")) if ver_dir.exists() else []
        if existing:
            print(f"Skipping release {release}: annotation already present at {ver_dir}")
            return True

        species_dir = species_name.lower().replace(' ', '_')
        # Try to list gtf dir for this release
        tried_urls = []
        gtf_dir_paths = [f"https://ftp.ensembl.org/pub/release-{release}/gtf/{species_dir}/",
                         f"http://ftp.ensembl.org/pub/release-{release}/gtf/{species_dir}/",
                         f"https://ftp.ensembl.org/pub/release_{release}/gtf/{species_dir}/",
                         f"http://ftp.ensembl.org/pub/release_{release}/gtf/{species_dir}/"]
        gtf_file_url = None
        for url in gtf_dir_paths:
            tried_urls.append(url)
            try:
                with urllib.request.urlopen(url, timeout=20) as r:
                    html = r.read().decode('utf-8', 'ignore')
                # find first .gtf.gz file link
                m = re.search(r'href=["\']([^"\']+\.gtf\.gz)["\']', html, flags=re.IGNORECASE)
                if m:
                    filepart = m.group(1)
                    if filepart.startswith('http'):
                        gtf_file_url = filepart
                    else:
                        gtf_file_url = urllib.parse.urljoin(url, filepart)
                    break
            except Exception:
                continue
        if not gtf_file_url:
            # No GTF for this release; silently ignore per user request
            return False

        # Print a compact download header for the release (no long URL)
        print(f"Download ensembl_{release}")

        # stream download to temp file
        try:
            with urllib.request.urlopen(gtf_file_url, timeout=120) as r:
                tmpf = tempfile.NamedTemporaryFile(delete=False)
                shutil.copyfileobj(r, tmpf)
                tmpf.close()
        except Exception as e:
            print(f"Failed to download {gtf_file_url}: {e}", file=sys.stderr)
            return False

        # Now that we have a GTF, create the version dir and parse
        try:
            ver_dir.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass

        # parse GTF.gz and write per-chromosome TSVs.
        # Some Ensembl GTFs may not contain explicit 'gene' features. In that
        # case, aggregate coordinates from other features (exon/transcript/CDS)
        # keyed by gene_id so we still produce one entry per gene with min/max
        # coordinates.
        gene_ranges = {}  # {(chrom, gene_id): {name, biotype, start, end}}
        handlers = {}
        try:
            with gzip.open(tmpf.name, 'rt', encoding='utf-8', errors='ignore') as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue
                    chrom, src, feature, start_s, end_s, score, strand, frame, attrs = parts[:9]
                    # parse attributes into dict
                    attrd = {}
                    for a in attrs.split(';'):
                        a = a.strip()
                        if not a:
                            continue
                        if ' ' in a:
                            k, v = a.split(' ', 1)
                            v = v.strip().strip('"')
                            attrd[k] = v
                    gene_id = attrd.get('gene_id') or attrd.get('geneID') or attrd.get('gene') or ''
                    gene_name = attrd.get('gene_name') or attrd.get('gene_symbol') or attrd.get('Name') or ''
                    biotype = attrd.get('gene_biotype') or attrd.get('gene_type') or attrd.get('gene_biotype') or ''
                    if not gene_id:
                        # skip entries we can't attribute to a gene
                        continue
                    try:
                        start = int(start_s)
                        end = int(end_s)
                    except Exception:
                        continue
                    key = (chrom, gene_id)
                    if key not in gene_ranges:
                        gene_ranges[key] = {
                            'gene_name': gene_name,
                            'biotype': biotype,
                            'start': start,
                            'end': end
                        }
                    else:
                        if start < gene_ranges[key]['start']:
                            gene_ranges[key]['start'] = start
                        if end > gene_ranges[key]['end']:
                            gene_ranges[key]['end'] = end

            # now write per-chromosome TSVs, show tqdm progress while writing
            # Only create TSVs for chromosomes listed in the species' chromosomes.txt
            try:
                # `chromosomes` is loaded earlier from data/<Species>/chromosomes.txt
                allowed = set([c.strip() for c in chromosomes if c.strip()])
                allowed_l = set([c.lower() for c in allowed] + [f"chr{c.lower()}" for c in allowed])
            except Exception:
                allowed = None
                allowed_l = None

            # filter items to only those where chrom matches allowed set (case-insensitive, with/without 'chr')
            items = []
            for (chrom, gid), info in gene_ranges.items():
                chrom_l = (chrom or '').lower()
                if allowed_l is not None:
                    if chrom_l in allowed_l or chrom_l.lstrip('chr') in allowed_l:
                        items.append(((chrom, gid), info))
                else:
                    items.append(((chrom, gid), info))

            if not items:
                # No allowed chromosomes found in this GTF -> do not keep an empty version dir
                try:
                    if ver_dir.exists():
                        # attempt to remove empty dir
                        try:
                            next(ver_dir.iterdir())
                            # dir not empty, keep it
                        except StopIteration:
                            ver_dir.rmdir()
                except Exception:
                    pass
                return False

            from tqdm import tqdm as _tqdm
            for (chrom, gid), info in _tqdm(items, desc=f"ensembl_{release} tsv", unit="genes"):
                out_path = ver_dir / f"{chrom}.tsv"
                if chrom not in handlers:
                    f = open(out_path, 'w', encoding='utf-8')
                    f.write('\t'.join(["gene_id","transcript_id","gene_name","chromosome","start","end","biotype"]) + '\n')
                    handlers[chrom] = f
                # transcript_id left blank because GTF gene aggregation may merge multiple transcripts
                line_out = '\t'.join([gid, '', info.get('gene_name',''), chrom, str(info['start']), str(info['end']), info.get('biotype','')])
                handlers[chrom].write(line_out + '\n')
        except Exception as e:
            print(f"Failed to parse GTF {tmpf.name}: {e}", file=sys.stderr)
            try:
                for h in handlers.values():
                    h.close()
            except Exception:
                pass
            try:
                os.unlink(tmpf.name)
            except Exception:
                pass
            return False
        finally:
            try:
                for h in handlers.values():
                    h.close()
            except Exception:
                pass
            try:
                if os.path.exists(tmpf.name):
                    os.unlink(tmpf.name)
            except Exception:
                pass

        # concise completion
        print(f"Downloaded ensembl_{release}")
        return True

    # Load or create a small local cache of releases we've already checked so
    # we don't repeatedly probe the FTP for the same release.
    checked_path = Path(base_annot_dir) / ".checked_releases.json"
    try:
        if checked_path.exists():
            with open(checked_path, 'r', encoding='utf-8') as cf:
                cj = json.load(cf)
                checked = set(str(x) for x in (cj.get('checked') or []))
        else:
            checked = set()
    except Exception:
        checked = set()

    releases = discover_releases()
    if not releases:
        print("No Ensembl releases discovered; aborting." , file=sys.stderr)
    else:
        # Process releases sequentially (safe for bandwidth). Could be parallelized later.
        processed = []
        for r in releases:
            if str(r) in checked:
                # skip releases we've previously checked
                continue
            try:
                ok = process_release(r)
                if ok:
                    processed.append(r)
            except Exception as e:
                print(f"Error processing release {r}: {e}", file=sys.stderr)
            finally:
                # mark this release as checked regardless of outcome to avoid
                # repeated probes in subsequent runs
                try:
                    checked.add(str(r))
                    tmp = checked_path.with_suffix('.tmp')
                    with open(tmp, 'w', encoding='utf-8') as cf:
                        json.dump({'checked': sorted(list(checked))}, cf, indent=2)
                    os.replace(str(tmp), str(checked_path))
                except Exception:
                    # non-fatal; continue
                    pass

        # Write manifest
        try:
            manifest_path = Path(base_annot_dir) / "versions.json"
            manifest_path.write_text(json.dumps({"ensembl_versions": processed}, indent=2), encoding='utf-8')
            print(f"Wrote versions manifest: {manifest_path}")
        except Exception as e:
            print(f"[WARN] Failed to write versions manifest: {e}", file=sys.stderr)
