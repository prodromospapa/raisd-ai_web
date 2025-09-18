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
    import urllib.request, urllib.error, urllib.parse, re, gzip, shutil, tempfile, json, subprocess, shlex, threading

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
                m = re.search(fr'href=["\']([^"\']+\.{re.escape(release)}\.gtf\.gz)["\']', html, flags=re.IGNORECASE)
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

    # Download progress will be displayed by the tqdm progress bar below

        # Try a fast pipeline: curl/wget -> pigz -dc -> parse stream. Fall back to
        # pure-Python streaming (urllib + gzip.GzipFile) if external tools are
        # unavailable.
        p_proc = None

        def parse_fh(fh):
            # fh: text-mode file-like iterator yielding lines
            gene_ranges = {}
            handlers = {}
            try:
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
                try:
                    allowed = set([c.strip() for c in chromosomes if c.strip()])
                    allowed_l = set([c.lower() for c in allowed] + [f"chr{c.lower()}" for c in allowed])
                except Exception:
                    allowed = None
                    allowed_l = None

                items = []
                for (chrom, gid), info in gene_ranges.items():
                    chrom_l = (chrom or '').lower()
                    if allowed_l is not None:
                        if chrom_l in allowed_l or chrom_l.lstrip('chr') in allowed_l:
                            items.append(((chrom, gid), info))
                    else:
                        items.append(((chrom, gid), info))

                if not items:
                    try:
                        if ver_dir.exists():
                            try:
                                next(ver_dir.iterdir())
                            except StopIteration:
                                ver_dir.rmdir()
                    except Exception:
                        pass
                    return False

                from tqdm import tqdm as _tqdm
                for (chrom, gid), info in _tqdm(items, desc=f"Splitting ensembl_{release} to tsv", unit="genes", position=1, leave=True, ascii=True, file=sys.stderr):
                    out_path = ver_dir / f"{chrom}.tsv"
                    if chrom not in handlers:
                        f = open(out_path, 'w', encoding='utf-8')
                        f.write('\t'.join(["gene_id","transcript_id","gene_name","chromosome","start","end","biotype"]) + '\n')
                        handlers[chrom] = f
                    line_out = '\t'.join([gid, '', info.get('gene_name',''), chrom, str(info['start']), str(info['end']), info.get('biotype','')])
                    handlers[chrom].write(line_out + '\n')
                return True
            except Exception as e:
                print(f"Failed to parse GTF stream: {e}", file=sys.stderr)
                try:
                    for h in handlers.values():
                        h.close()
                except Exception:
                    pass
                return False
            finally:
                try:
                    for h in handlers.values():
                        h.close()
                except Exception:
                    pass

        # create version dir before parsing outputs
        try:
            ver_dir.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass

        # attempt pigz pipeline
        pigz = shutil.which('pigz')
        curl = shutil.which('curl')
        wget = shutil.which('wget')
        if pigz and (curl or wget):
            # Build downloader command (no shell piping). We'll stream
            # downloader.stdout -> pigz.stdin while updating a tqdm bar, and
            # parse pigz.stdout in a background thread.
            if curl:
                dl_cmd = [curl, '-L', '-s', gtf_file_url]
            else:
                dl_cmd = [wget, '-q', '-O', '-', gtf_file_url]
            pigz_cmd = [pigz, '-dc']

            # try to get Content-Length for a progress bar
            total = None
            try:
                head_req = urllib.request.Request(gtf_file_url, headers={'User-Agent': 'raisd-ai/1.0 (python urllib)'}, method='HEAD')
                with urllib.request.urlopen(head_req, timeout=20) as hr:
                    cl = hr.getheader('Content-Length')
                    try:
                        total = int(cl) if cl is not None else None
                    except Exception:
                        total = None
            except Exception:
                total = None

            try:
                # Open the HTTP response directly so we control chunk reads and
                # can reliably update tqdm. Then feed bytes into pigz stdin.
                req = urllib.request.Request(gtf_file_url, headers={'User-Agent': 'raisd-ai/1.0 (python urllib)'})
                r = urllib.request.urlopen(req, timeout=120)
                p2 = subprocess.Popen(pigz_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                if p2.stdin is None or p2.stdout is None:
                    raise RuntimeError('Failed to create pigz pipes')

                # start parser thread reading pigz stdout
                parser_result = {'ok': False, 'error': None}

                def _parser():
                    try:
                        if p2.stdout is None:
                            parser_result['error'] = 'pigz produced no stdout'
                            return
                        fh = io.TextIOWrapper(p2.stdout, encoding='utf-8', errors='ignore')
                        parser_result['ok'] = parse_fh(fh)
                        try:
                            fh.close()
                        except Exception:
                            pass
                    except Exception as e:
                        parser_result['error'] = str(e)

                th = threading.Thread(target=_parser, daemon=True)
                th.start()

                # stream bytes from downloader -> pigz.stdin while updating tqdm
                chunk_size = 64 * 1024
                desc = f"Download ensembl_{release}"
                if total:
                    pbar = tqdm(total=total, unit='B', unit_scale=True, desc=desc, position=0, leave=True, ascii=True, file=sys.stderr)
                else:
                    pbar = tqdm(unit='B', unit_scale=True, desc=desc, position=0, leave=True, ascii=True, file=sys.stderr)

                try:
                    while True:
                        chunk = r.read(chunk_size)
                        if not chunk:
                            break
                        p2.stdin.write(chunk)
                        p2.stdin.flush()
                        pbar.update(len(chunk))
                finally:
                    pbar.close()
                    try:
                        p2.stdin.close()
                    except Exception:
                        pass

                # wait for parser thread and pigz
                th.join()
                p2_stderr = b''
                try:
                    p2_stderr = p2.stderr.read() if p2.stderr is not None else b''
                except Exception:
                    pass
                p2.wait()
                try:
                    r.close()
                except Exception:
                    pass

                if not parser_result.get('ok'):
                    if parser_result.get('error'):
                        print(f"pigz parser error: {parser_result['error']}", file=sys.stderr)
                    if p2_stderr:
                        try:
                            print(f"[pigz stderr] {p2_stderr.decode('utf-8','ignore')}", file=sys.stderr)
                        except Exception:
                            pass
                    # no downloader stderr because we read directly via urllib
                    
                    # fall back to Python streaming
                else:
                    # success
                    try:
                        if p2_stderr:
                            errtxt = p2_stderr.decode('utf-8','ignore')
                            if errtxt.strip():
                                print(f"[pigz stderr] {errtxt}", file=sys.stderr)
                    except Exception:
                        pass
                    # parser succeeded
                    return True
            except Exception as e:
                print(f"pigz pipeline failed: {e}", file=sys.stderr)
                # fall through to Python streaming below

        # fallback: pure-Python streaming via urllib + gzip
        try:
            req = urllib.request.Request(gtf_file_url, headers={'User-Agent': 'raisd-ai/1.0 (python urllib)'})
            with urllib.request.urlopen(req, timeout=120) as r:
                gz = gzip.GzipFile(fileobj=r)
                fh = io.TextIOWrapper(gz, encoding='utf-8', errors='ignore')
                ok = parse_fh(fh)
                try:
                    fh.close()
                except Exception:
                    pass
                if not ok:
                    return False
        except Exception as e:
            print(f"Failed to download/stream {gtf_file_url}: {e}", file=sys.stderr)
            return False

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
                    # Persist this release as checked immediately so that if the
                    # user cancels the run later, completed releases are not lost.
                    try:
                        checked.add(str(r))
                        tmp = checked_path.with_suffix('.tmp')
                        with open(tmp, 'w', encoding='utf-8') as cf:
                            json.dump({'checked': sorted(list(checked))}, cf, indent=2)
                        os.replace(str(tmp), str(checked_path))
                    except Exception:
                        # non-fatal; we'll still attempt to persist in the loop's finally
                        pass
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
        # NOTE: per user request, do NOT create or update versions.json files.
        # The original behaviour wrote a `versions.json` manifest listing
        # processed Ensembl releases. To avoid creating this file at all,
        # skip writing it here. If in future you need to enable writing,
        # reintroduce the write behind an explicit flag.
        try:
            # Intentionally skip writing versions.json to avoid creating the file.
            pass
        except Exception:
            # keep original behaviour of not failing the whole run if manifest
            # writing would fail — but since we don't write, nothing to do.
            pass
