
import os, re, json, time, shutil, signal, base64, tempfile, threading, subprocess, queue, sys
from decimal import Decimal
from flask import Flask, render_template, request, jsonify, send_from_directory, stream_with_context

try:
    import psutil
except Exception:
    psutil = None

try:
    import stdpopsim as sps
except Exception:
    sps = None

app = Flask(__name__)
PORT = 5001
SIM_SCRIPT = "simulator.py"
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ARTIFACTS_DIR = os.path.join(BASE_DIR, "artifacts")

STATE = {}
LOCK = threading.Lock()
RUNNERS = {}

# Runner structure per-sid:
# RUNNERS[sid] = {
#   'proc': Popen,
#   'thread': Thread,
#   'queue': Queue,
#   'run_dir': str,
#   'cancel_path': str,
# }

def get_sid():
    sid = request.cookies.get("sid")
    if not sid:
        sid = f"{int(time.time()*1000)}-{os.getpid()}"
    return sid

def get_state(sid):
    with LOCK:
        if sid not in STATE:
            STATE[sid] = {
                "species_id": None, "model_id": None, "chromosome": None,
                "populations": [],
                "engine": "msms", "replicates": 1,
                "seq_length": 0, "target_snps": 0, "chrom_length_mode": False,
                "target_snps_tol_auto": True, "target_snps_tol": 10.0,
                "output_format": "ms", "output_path": "", "temp_loc": "local",
                "parallel": 1, "sims_per_work_auto": True, "sims_per_work": 1,
                "max_ram_percent": 80, "max_ram_cap": "system",
                "growth_max_fold": "1.05",
                "sfs_on": False, "sfs_output": "", "sfs_normalized": False, "sfs_mode": "mean",
                "sweep_enable": False, "sweep_pop": "", "sweep_pos_raw": 50.0,
                "sel_s": 0.1, "time_mode": "", "sweep_time": "", "fixation_time": "", "time_units": "gens",
                "progress_flag": False, "run_sfs_only": False,
                "paired_neutral": False, "paired_neutral_name": "", "neutral_engine": "",
                "_contig_len": None, "_last_run_dir": None,
            }
        return STATE[sid]

def list_models(sp):
    try:
        if hasattr(sp, "demographic_models") and sp.demographic_models:
            return sorted({m.id for m in sp.demographic_models})
        if hasattr(sp, "models"):
            return sorted({m.id for m in sp.models})
    except Exception:
        pass
    return []

def list_contigs(sp):
    try:
        genome = getattr(sp, "genome", None)
        contigs = []
        if genome is not None:
            contigs = getattr(genome, "contigs", []) or getattr(genome, "chromosomes", []) or []
        names = []
        for c in contigs:
            if isinstance(c, str):
                nm = c
            else:
                nm = getattr(c, "id", None) or getattr(c, "name", None)
            if nm: names.append(nm)
        return names
    except Exception:
        return []

def normalize_contig_names(contigs):
    def norm(nm):
        if not isinstance(nm, str):
            return str(nm)
        low = nm.lower()
        if low.startswith("chr"):
            rest = nm[3:]
            return rest.upper() if rest.upper() in {"X","Y","MT","M"} else rest
        return nm
    seen=set(); out=[]
    for nm in contigs or []:
        lab = norm(nm)
        if lab not in seen:
            seen.add(lab); out.append(lab)
    def key(s):
        su=s.upper()
        if s.isdigit(): return (0,int(s))
        if su=="X": return (1,0)
        if su=="Y": return (1,1)
        if su in {"MT","M"}: return (1,2)
        return (2,s.lower())
    return sorted(out, key=key)

def list_pop_names(sp_obj, model_id):
    try:
        model = sp_obj.get_demographic_model(model_id)
    except Exception:
        return []
    try:
        g = getattr(model, 'model', None)
        if g and hasattr(g, 'demes'):
            return [d.name for d in g.demes]
    except Exception:
        pass
    out = []
    try:
        for p in getattr(model, 'populations', []):
            nm = getattr(p, 'name', None) or getattr(p, 'id', None)
            if nm: out.append(nm)
    except Exception:
        pass
    return out

def derive_state(stt):
    species_items = []
    if sps:
        try:
            for sp in sps.all_species():
                mids = list_models(sp)
                if mids:
                    label_parts = [sp.id]
                    sci = getattr(sp, "name", None)
                    if sci: label_parts.append(f"— {sci}")
                    common = getattr(sp, "common_name", None)
                    if common: label_parts.append(f"({common})")
                    species_items.append({"id": sp.id, "label": " ".join(label_parts)})
        except Exception:
            pass
    species_items.sort(key=lambda x: x["label"].lower())

    model_options=[]; contigs=[]; contig_len=None
    if sps and stt.get("species_id"):
        try:
            sp_obj = sps.get_species(stt["species_id"])
            model_options = list_models(sp_obj)
            raw_contigs = list_contigs(sp_obj)
            contigs = normalize_contig_names(raw_contigs)
            if stt.get("chromosome"):
                genome = getattr(sp_obj, "genome", None)
                raw = getattr(genome, "contigs", []) or getattr(genome, "chromosomes", []) or []
                def norm(c):
                    nm = c if isinstance(c,str) else (getattr(c,"id",None) or getattr(c,"name",None))
                    if not isinstance(nm,str): nm=str(nm)
                    low = nm.lower()
                    if low.startswith('chr'):
                        rest=nm[3:]
                        return rest.upper() if rest.upper() in {'X','Y','MT','M'} else rest
                    return nm
                for c in raw:
                    if norm(c)==stt["chromosome"]:
                        contig_len = getattr(c,"length",None) or getattr(c,"sequence_length",None) or getattr(c,"size",None) or getattr(c,"bp_length",None)
                        break
        except Exception:
            pass
    stt["_contig_len"]=contig_len

    # populations sync
    if sps and stt.get("species_id") and stt.get("model_id") and stt.get("chromosome"):
        try:
            sp_obj = sps.get_species(stt["species_id"])
            pop_names = list_pop_names(sp_obj, stt["model_id"])
        except Exception:
            pop_names = []
        if pop_names:
            existing = {p["name"]: p for p in stt.get("populations", [])}
            stt["populations"] = [existing.get(nm, {"name": nm, "selected": False, "n": 0}) for nm in pop_names]

    chosen_idx = [i for i,p in enumerate(stt.get("populations", [])) if p.get("selected")]
    ordered_ready=False
    if chosen_idx:
        ordered_ready = ((max(chosen_idx)-min(chosen_idx)+1)==len(chosen_idx)) and all(int(stt["populations"][i]["n"])>0 for i in chosen_idx)

    try:
        logical = psutil.cpu_count(logical=True) if psutil else os.cpu_count()
    except Exception:
        logical = os.cpu_count()
    try:
        physical = psutil.cpu_count(logical=False) if psutil else None
    except Exception:
        physical = None
    engine = stt.get("engine","msms")
    engine_max = int((physical or logical or 1) if engine=="msprime" else (logical or 1))
    reps = int(stt.get("replicates",1) or 1)
    max_workers = max(1, min(engine_max, reps))

    # Collect informational/error messages mirroring the Streamlit UI
    messages = []
    try:
        # If stdpopsim available but no species with models found
        if sps and not species_items:
            messages.append("No species with demographic models found in stdpopsim. Install or update stdpopsim data.")
        # If a species selected but no models present
        if stt.get("species_id") and not model_options:
            messages.append("No models available for the selected species.")
        # If species/model selected but no contigs
        if stt.get("species_id") and stt.get("model_id") and not contigs:
            messages.append("No chromosomes/contigs available for this species in stdpopsim. The `--chromosome` argument will be omitted.")
        # Require steps 1-4 prior to preparing the engine command
        if not ordered_ready:
            messages.append("Complete steps 1–4 (species, model, ordered populations with counts) before preparing the engine command.")
        # Missing simulator script
        if not os.path.exists(SIM_SCRIPT):
            messages.append("Could not find simulator.py in the working folder.")
        # Growth max fold numeric hints
        try:
            gm = stt.get('growth_max_fold')
            if gm is not None and str(gm).strip() != '':
                try:
                    gfv = float(str(gm))
                    if not (gfv > 1.0):
                        messages.append("Enter a decimal number strictly greater than 1.0")
                except Exception:
                    messages.append("Enter a valid decimal number (e.g. 1.05).")
        except Exception:
            pass
        # Check for an existing run (approximate Streamlit 'run in progress' message)
        try:
            if RUNNERS:
                messages.append("A run is already in progress for this model (another process is active).")
        except Exception:
            pass
        # SFS / output related hints
        try:
            sfs_on = bool(stt.get('sfs_on', False))
            run_sfs_only = bool(stt.get('run_sfs_only', False))
            outp = (stt.get('output_path') or '').strip()
            sfs_out = (stt.get('sfs_output') or '').strip()
            if sfs_on and run_sfs_only and not sfs_out:
                messages.append("Run SFS only is checked: provide an SFS output path.")
            if not (run_sfs_only and sfs_out) and not outp and not run_sfs_only:
                messages.append("Primary Output path is required unless 'Run SFS only' is checked and an SFS output path is provided.")
        except Exception:
            pass
        # Target SNPs tolerance
        try:
            if not bool(stt.get('target_snps_tol_auto', True)):
                raw_tol = float(stt.get('target_snps_tol', 0.0) or 0.0)
                tol_percent = raw_tol * 100.0 if (0.0 < raw_tol <= 1.0) else raw_tol
                if tol_percent <= 0.0:
                    messages.append("Target SNPs tolerance must be a positive value greater than 0 when Custom is selected.")
        except Exception:
            pass
        # Paired neutral name
        try:
            if bool(stt.get('paired_neutral', False)):
                if not (stt.get('paired_neutral_name') and str(stt.get('paired_neutral_name')).strip()):
                    messages.append("Paired neutral is checked: provide a Paired Neutral output base name before preparing/executing the command.")
        except Exception:
            pass
        # Length/target/chromosome requirement
        try:
            chrom_mode = bool(stt.get('chrom_length_mode', False))
            seq_val = int(stt.get('seq_length', 0) or 0)
            snp_val = int(stt.get('target_snps', 0) or 0)
            if not chrom_mode and seq_val <= 0 and snp_val <= 0:
                messages.append("Set either Length (bp), Target SNPs (>0), or enable Chromosome Length. At least one is required.")
        except Exception:
            pass
        # Sweep-related hints
        try:
            sweep_enabled = bool(stt.get('sweep_enable', False))
            if sweep_enabled:
                if not stt.get('sweep_pop'):
                    messages.append("Sweep is enabled: choose a Sweep population before preparing/executing the command.")
                tm = stt.get('time_mode') or ''
                sw = stt.get('sweep_time', '')
                fx = stt.get('fixation_time', '')
                if not tm:
                    messages.append("Pick either Sweep Time or Fixation Time before running the sweep.")
                else:
                    if tm == 'Sweep Time':
                        try:
                            if isinstance(sw, str):
                                if sw.strip() == '': raise ValueError()
                            float(sw)
                        except Exception:
                            messages.append("You selected Sweep Time but did not enter a valid numeric value. Please enter a numeric Sweep Time to proceed.")
                    elif tm == 'Fixation Time':
                        try:
                            if isinstance(fx, str):
                                if fx.strip() == '': raise ValueError()
                            float(fx)
                        except Exception:
                            messages.append("You selected Fixation Time but did not enter a valid numeric value. Please enter a numeric Fixation Time to proceed.")
        except Exception:
            pass
    except Exception:
        pass

    return {
        "species_items": species_items,
        "model_options": model_options,
        "contigs": contigs,
        "contig_len": contig_len,
        "ordered_ready": ordered_ready,
        "max_workers": max_workers,
        "reps": reps,
        "engine": engine,
        "messages": messages,
    }


@app.route("/")
def index():
    sid = get_sid()
    stt = get_state(sid)
    derived = derive_state(stt)
    resp = render_template("index.html", stt=stt, derived=derived)
    from flask import make_response
    r = make_response(resp)
    if not request.cookies.get("sid"):
        r.set_cookie("sid", sid, samesite="Lax")
    return r

@app.route("/api/set", methods=["POST"])
def api_set():
    sid = get_sid()
    stt = get_state(sid)
    data = request.get_json(silent=True) or {}
    # simple update path; accepts direct field->value and pop_selected_i/pop_n_i
    def to_bool(v):
        if isinstance(v, bool): return v
        return str(v).lower() in ("1","true","on","yes")
    for k,v in data.items():
        if k.startswith("pop_selected_"):
            i = int(k.split("_")[-1])
            if 0 <= i < len(stt["populations"]):
                stt["populations"][i]["selected"] = to_bool(v)
        elif k.startswith("pop_n_"):
            i = int(k.split("_")[-1])
            if 0 <= i < len(stt["populations"]):
                try: stt["populations"][i]["n"] = int(v)
                except: stt["populations"][i]["n"] = 0
        elif k in {"replicates","seq_length","target_snps","parallel","sims_per_work","max_ram_percent"}:
            try: stt[k] = int(v)
            except: stt[k] = 0
        elif k in {"target_snps_tol","sel_s","sweep_pos_raw"}:
            try: stt[k] = float(v)
            except: stt[k] = 0.0
        elif k in {"chrom_length_mode","target_snps_tol_auto","sfs_on","sfs_normalized","sweep_enable","progress_flag","run_sfs_only","paired_neutral","sims_per_work_auto"}:
            stt[k] = to_bool(v)
        else:
            stt[k] = v

    # Apply automations similar to Streamlit UI: normalize output extension,
    # auto-name SFS and paired-neutral when appropriate.
    try:
        _normalize_and_autoname(stt)
    except Exception:
        pass

    derived = derive_state(stt)
    return jsonify({"ok": True, "stt": stt, "derived": derived})


def _normalize_and_autoname(stt: dict):
    """Mutate stt in-place to apply output/path automations:
    - normalize primary output extension to match output_format
    - ensure sfs_output ends with .sfs (auto-derive from output when empty)
    - when Run SFS only toggled off and sfs exists and output empty, derive output from sfs and chosen output_format
    - when paired_neutral checked and paired_neutral_name empty, derive default from output or sfs
    This follows the Streamlit `ui_simulator.py` behavior.
    """
    try:
        fmt_ext_map = {
            "ms": ".ms",
            "ms.gz": ".ms.gz",
            "vcf": ".vcf",
            "vcf.gz": ".vcf.gz",
            "bcf": ".bcf",
        }
        out_key = 'output_path'
        fmt_key = 'output_format'
        sfs_key = 'sfs_output'
        run_sfs_key = 'run_sfs_only'

        out_val = (stt.get(out_key) or '').strip()
        fmt = str(stt.get(fmt_key, '') or '').strip()
        sfs_val = (stt.get(sfs_key) or '').strip()
        run_sfs_only = bool(stt.get(run_sfs_key, False))

        # Normalize primary output extension to output_format
        if out_val and fmt:
            desired_ext = fmt_ext_map.get(fmt, '')
            # strip any known ext first
            for e in ('.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs'):
                if out_val.endswith(e):
                    out_val = out_val[:-len(e)]
                    break
            if desired_ext:
                out_val = out_val + desired_ext
            stt[out_key] = out_val

        # If SFS enabled and sfs_output empty but primary output exists, derive sfs
        if bool(stt.get('sfs_on', False)):
            if (not sfs_val) and out_val:
                base = os.path.basename(out_val)
                # strip known extensions
                for e in ('.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf'):
                    if base.endswith(e):
                        base = base[:-len(e)]; break
                else:
                    base = os.path.splitext(base)[0]
                sfs_val = base + '.sfs'
                stt[sfs_key] = sfs_val
            # ensure sfs has .sfs extension when user typed a base name
            if sfs_val and not sfs_val.lower().endswith('.sfs'):
                stt[sfs_key] = sfs_val + '.sfs'

        # When Run SFS only is False and sfs exists but primary output empty, derive primary output
        if not run_sfs_only and (not out_val) and sfs_val:
            base = os.path.basename(sfs_val)
            if base.lower().endswith('.sfs'):
                base = base[:-4]
            for e in ('.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf'):
                if base.endswith(e):
                    base = base[:-len(e)]; break
            desired_ext = fmt_ext_map.get(fmt, '')
            if desired_ext:
                stt[out_key] = base + desired_ext
            else:
                stt[out_key] = base

        # Paired neutral default name: set when paired_neutral is True and name empty
        try:
            if bool(stt.get('paired_neutral', False)):
                    pn_key = 'paired_neutral_name'
                    pn_val = (stt.get(pn_key) or '').strip()
                    if not pn_val:
                        # Prefer deriving paired-neutral base from the primary output
                        # path if the user provided one (even when Run SFS only is
                        # checked and we suppress passing --output). Otherwise, if
                        # running SFS-only, derive from the SFS basename. If neither
                        # is available, do not auto-generate a timestamped name so
                        # we won't pass anything to the simulator.
                        if stt.get(out_key):
                            b = os.path.basename(str(stt.get(out_key)))
                            for e in ('.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs'):
                                if b.endswith(e):
                                    b = b[:-len(e)]; break
                            else:
                                b = os.path.splitext(b)[0]
                            pn_base = b + '_neutral'
                            stt[pn_key] = pn_base
                        elif run_sfs_only and stt.get(sfs_key):
                            b = os.path.basename(str(stt.get(sfs_key)))
                            if b.lower().endswith('.sfs'):
                                b = b[:-4]
                            pn_base = os.path.splitext(b)[0] + '_neutral'
                            stt[pn_key] = pn_base
                        else:
                            # Leave paired_neutral_name empty if no basis is available
                            # so we won't pass paired-neutral to the simulator.
                            pass
        except Exception:
            pass
    except Exception:
        pass

# ---- Minimal "Show engine command" passthrough ----
def strip_ext(name):
    if not name: return name
    b = os.path.basename(name)
    for e in ('.ms.gz','.vcf.gz','.ms','.vcf','.bcf','.sfs','.zip'):
        if b.endswith(e):
            return b[:-len(e)]
    return os.path.splitext(b)[0]


def _ensure_run_dir(sid):
    ts = int(time.time())
    run_name = f"run_{ts}_{os.getpid()}"
    # Ensure base artifacts directory exists alongside app.py
    os.makedirs(ARTIFACTS_DIR, exist_ok=True)
    run_dir = os.path.join(ARTIFACTS_DIR, run_name)
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def _queue_put_safe(q, item):
    try:
        q.put(item, timeout=0.5)
    except Exception:
        pass

def _collect_validation_errors(stt: dict):
    """Collect validation errors mirroring Streamlit Build & Run gating."""
    errors = []
    try:
        # SFS/run_sfs_only vs outputs
        sfs_on = bool(stt.get('sfs_on', False))
        run_sfs_only = bool(stt.get('run_sfs_only', False))
        outp = (stt.get('output_path') or '').strip()
        sfs_out = (stt.get('sfs_output') or '').strip()
        if sfs_on and run_sfs_only and not sfs_out:
            errors.append("Run SFS only is checked: provide an SFS output path.")
        if not (run_sfs_only and sfs_out):
            if not outp and not run_sfs_only:
                errors.append("Primary Output path is required unless 'Run SFS only' is checked and an SFS output path is provided.")
    except Exception:
        pass

    # tolerance custom check
    try:
        tol_auto = bool(stt.get('target_snps_tol_auto', True))
        if not tol_auto:
            raw_tol = float(stt.get('target_snps_tol', 0.0) or 0.0)
            tol_percent = raw_tol * 100.0 if (0.0 < raw_tol <= 1.0) else raw_tol
            if tol_percent <= 0.0:
                errors.append("Target SNPs tolerance must be a positive value greater than 0 when Custom is selected.")
    except Exception:
        pass

    # paired neutral name required when PN checked
    try:
        if bool(stt.get('paired_neutral', False)):
            if not (stt.get('paired_neutral_name') and str(stt.get('paired_neutral_name')).strip()):
                errors.append("Paired neutral is checked: provide a Paired Neutral output base name before preparing/executing the command.")
    except Exception:
        pass

    # length/target/chromosome requirement (at least one)
    try:
        chrom_mode = bool(stt.get('chrom_length_mode', False))
        seq_val = int(stt.get('seq_length', 0) or 0)
        snp_val = int(stt.get('target_snps', 0) or 0)
        if not chrom_mode and seq_val <= 0 and snp_val <= 0:
            errors.append("Set either Length (bp), Target SNPs (>0), or enable Chromosome Length. At least one is required.")
    except Exception:
        pass

    # Selection-related checks
    try:
        if bool(stt.get('sweep_enable', False)):
            # require sweep population
            if not (stt.get('sweep_pop') and str(stt.get('sweep_pop')).strip()):
                errors.append("Sweep is enabled: choose a Sweep population before preparing/executing the command.")
            # time mode/value
            tm = (stt.get('time_mode') or '').strip()
            sv = stt.get('sweep_time', '')
            fv = stt.get('fixation_time', '')
            if not tm:
                errors.append("Sweep is enabled: pick either Sweep Time or Fixation Time before preparing/executing the command.")
            else:
                if tm == 'Sweep Time':
                    try:
                        if isinstance(sv, str):
                            if sv.strip() == '':
                                raise ValueError()
                        float(sv)
                    except Exception:
                        errors.append("You selected Sweep Time but did not enter a valid numeric value. Please enter a numeric Sweep Time to proceed.")
                elif tm == 'Fixation Time':
                    try:
                        if isinstance(fv, str):
                            if fv.strip() == '':
                                raise ValueError()
                        float(fv)
                    except Exception:
                        errors.append("You selected Fixation Time but did not enter a valid numeric value. Please enter a numeric Fixation Time to proceed.")
            # discoal-specific: either sweep_time>0 or fixation_time>=0
            try:
                if stt.get('engine') == 'discoal':
                    sw_ok = False; fx_ok = False
                    try:
                        if isinstance(sv, str):
                            if sv.strip() != '':
                                sw_ok = float(sv) > 0.0
                        elif sv is not None:
                            sw_ok = float(sv) > 0.0
                    except Exception:
                        sw_ok = False
                    try:
                        if isinstance(fv, str):
                            if fv.strip() != '':
                                fx_ok = float(fv) >= 0.0
                        elif fv is not None:
                            fx_ok = float(fv) >= 0.0
                    except Exception:
                        fx_ok = False
                    if not (sw_ok or fx_ok):
                        errors.append("For discoal sweeps provide a positive Sweep Time or a valid Fixation Time.")
            except Exception:
                pass
    except Exception:
        pass

    return errors

def build_cmd(stt, include_progress=False, run_dir=None):
    # Mirrors original builder but supports run_dir override and returns expected artifact paths.
    pop_order = []
    counts = []
    ch = [i for i, p in enumerate(stt.get("populations", [])) if p.get("selected")]
    if not ch:
        return None, ["Select at least one population."], {}
    if (max(ch) - min(ch) + 1) != len(ch):
        return None, ["Selected populations must form a contiguous block."], {}
    for i in ch:
        p = stt["populations"][i]
        if int(p.get("n", 0)) <= 0:
            return None, ["Selected populations must have positive sample sizes."], {}
        pop_order.append(p["name"])
        counts.append(int(p["n"]))

    engine = stt.get("engine", "msms")
    cmd = [sys.executable, SIM_SCRIPT, "--engine", engine,
        "--pop-order", ",".join(pop_order),
        "--sample-individuals", ",".join(map(str, counts))]
    for flag in ("species_id", "model_id", "chromosome"):
        if stt.get(flag):
            cmd += [f"--{flag.replace('_', '-')}", str(stt[flag])]

    # Length vs target snps
    if not stt.get("chrom_length_mode"):
        if int(stt.get("seq_length", 0)) > 0:
            cmd += ["--length", str(int(stt["seq_length"]))]
        if int(stt.get("target_snps", 0)) > 0:
            cmd += ["--target-snps", str(int(stt["target_snps"]))]
            if not stt.get("target_snps_tol_auto", True):
                raw = float(stt.get("target_snps_tol", 0.0) or 0.0)
                tol_percent = raw * 100.0 if (0.0 < raw <= 1.0) else raw
                if tol_percent > 0.0:
                    cmd += ["--target-snps-tol", str(tol_percent)]

    # Replicates
    if int(stt.get("replicates", 1)) != 1:
        cmd += ["--replicates", str(int(stt["replicates"]))]

    # Parallel workers (only when >1)
    try:
        par = int(stt.get("parallel", 1) or 1)
        if par > 1:
            cmd += ["--parallel", str(par)]
    except Exception:
        pass

    outp = (stt.get("output_path", "") or "").strip()
    fmt = stt.get("output_format", "ms")
    run_sfs_only = bool(stt.get("run_sfs_only", False))
    # If running SFS-only, ignore any primary output_path the user may have
    # filled in so we don't accidentally pass --output to the simulator.
    outp_cmd = outp if not run_sfs_only else ''
    extmap = {"ms": ".ms", "ms.gz": ".ms.gz", "vcf": ".vcf", "vcf.gz": ".vcf.gz", "bcf": ".bcf"}

    artifacts = {}
    # Output/variants
    # Validation parity: if Compute SFS only is not set, require primary output
    if outp_cmd and not run_sfs_only:
        ext = extmap.get(fmt, "")
        base = strip_ext(os.path.basename(outp))
        out_name = base + (ext or "")
        if run_dir:
            out_final = os.path.join(run_dir, out_name)
        else:
            out_final = outp if outp.endswith(ext) else (outp + ext)
        cmd += ["--output", out_final, "--output-format", fmt]
        artifacts['variants'] = out_final
    elif not outp_cmd and not run_sfs_only:
        return None, ["Primary Output path is required unless 'Run SFS only' is checked and an SFS output path is provided."], {}

    # SFS output auto-name when primary output provided
    if bool(stt.get("sfs_on", False)):
        # prefer explicit sfs_output if provided
        sfs_out = (stt.get("sfs_output", "") or "").strip()
        # If no explicit sfs_output, try to derive from the primary variants path
        # (artifacts['variants']) or fall back to the original user-specified
        # output_path (outp) even when run_sfs_only suppressed the command's
        # primary output. This ensures we still generate and pass an --sfs
        # argument when the user requested SFS.
        if not sfs_out:
            if artifacts.get('variants'):
                base = strip_ext(os.path.basename(artifacts['variants']))
                sfs_out = base + ".sfs"
            elif outp:
                base = strip_ext(os.path.basename(outp))
                sfs_out = base + ".sfs"
        if sfs_out:
            sfs_path = os.path.join(run_dir, os.path.basename(sfs_out)) if run_dir else sfs_out
            artifacts['sfs'] = sfs_path
            # Always pass --sfs to the simulator when SFS calculation is requested.
            # simulator.py accepts --sfs [path] (optional value); provide the
            # per-run sfs path so the simulator writes SFS to the run artifacts.
            cmd += ["--sfs", sfs_path]
            # SFS normalization and mode flags
            if bool(stt.get('sfs_normalized', False)):
                cmd += ["--sfs-normalized"]
            sm = stt.get('sfs_mode') or ''
            if sm:
                cmd += ["--sfs-mode", str(sm)]
        else:
            # Pass bare --sfs when enabled but no path is specified/derived
            cmd += ["--sfs"]

    # Sims-per-work: only when engine is not msprime and user selected Custom
    try:
        if engine != "msprime" and not bool(stt.get("sims_per_work_auto", True)):
            spw = int(stt.get("sims_per_work", 0) or 0)
            if spw > 0:
                cmd += ["--sims-per-work", str(spw)]
    except Exception:
        pass

    # Max RAM percent and cap
    try:
        mr = float(stt.get("max_ram_percent", 80) or 80)
        cmd += ["--max-ram-percent", str(mr)]
        mrcap = (stt.get("max_ram_cap") or "system").strip() or "system"
        cmd += ["--max-ram-cap", mrcap]
    except Exception:
        pass

    # Temp location: pass only when non-local
    try:
        tmp = (stt.get("temp_loc") or "local").strip().lower()
        if tmp and tmp not in {"local", ""}:
            cmd += ["--temp", tmp]
    except Exception:
        pass

    # paired neutral (artifacts + flags)
    if bool(stt.get("paired_neutral", False)):
        pn_name = (stt.get("paired_neutral_name") or "").strip()
        # If name not provided, prefer the user's original output_path (outp)
        # even when run_sfs_only suppressed passing --output. Next prefer the
        # generated variants artifact, then (for SFS-only runs) the SFS basename.
        # If none are present, do not auto-generate a name (we won't add a
        # paired_neutral artifact in that case).
        if not pn_name:
            if outp:
                pn_name = strip_ext(os.path.basename(outp))
            elif artifacts.get('variants'):
                pn_name = strip_ext(os.path.basename(artifacts['variants']))
            elif run_sfs_only and artifacts.get('sfs'):
                pn_name = strip_ext(os.path.basename(artifacts['sfs']))
            else:
                pn_name = ''

        if pn_name:
            pn_base = strip_ext(pn_name)
            # choose neutral engine extension if available
            neutral_engine = stt.get('neutral_engine') or stt.get('engine') or 'ms'
            pn_ext = extmap.get(neutral_engine, '.ms')
            # For SFS-only runs where we derived from the SFS basename, the UI
            # requested behavior is to NOT append an engine-specific extension.
            if run_sfs_only and artifacts.get('sfs') and not artifacts.get('variants'):
                pn_file = pn_base
            else:
                pn_file = pn_base + pn_ext
            artifacts['paired_neutral'] = os.path.join(run_dir, pn_file) if run_dir else pn_file
            # Pass paired-neutral base (without extension) to simulator
            cmd += ["--paired-neutral", pn_base]
            if (stt.get('neutral_engine') or '').strip():
                cmd += ["--neutral-engine", (stt.get('neutral_engine') or '').strip()]
        else:
            # Flag without a value when user wants paired neutral but no explicit name
            cmd += ["--paired-neutral"]

    if include_progress:
        cmd += ["--progress"]

    # Growth max fold (>1) for non-msprime engines
    try:
        gmf = str(stt.get('growth_max_fold') or '').strip()
        if gmf:
            g = float(gmf)
            if g > 1.0 and gmf != '1.05':
                cmd += ["--growth-max-fold", gmf]
    except Exception:
        pass

    # Selection / sweep parameters
    try:
        if bool(stt.get('sweep_enable', False)):
            spop = (stt.get('sweep_pop') or '').strip()
            if spop:
                cmd += ["--sweep-pop", spop]
            try:
                s = float(stt.get('sel_s', 0.0) or 0.0)
                if s != 0.0:
                    cmd += ["--sel-s", str(s)]
            except Exception:
                pass
            # Compute sweep position as percent 0..100 from raw (percent or bp)
            raw = stt.get('sweep_pos_raw', 50.0)
            try:
                rawf = float(raw)
            except Exception:
                rawf = 50.0
            # Determine if raw is bp: prefer explicit seq_length, else contig when chrom_length_mode True
            length_for_pos = None
            try:
                sl = int(stt.get('seq_length', 0) or 0)
                if sl > 0:
                    length_for_pos = sl
                elif bool(stt.get('chrom_length_mode', False)) and stt.get('_contig_len'):
                    length_for_pos = int(stt.get('_contig_len'))
            except Exception:
                length_for_pos = None
            try:
                if length_for_pos and length_for_pos > 0:
                    pct = (rawf / float(length_for_pos)) * 100.0
                else:
                    pct = rawf  # already percent
            except Exception:
                pct = 50.0
            # clamp
            try:
                if pct < 0.0: pct = 0.0
                if pct > 100.0: pct = 100.0
            except Exception:
                pct = max(0.0, min(100.0, pct))
            cmd += ["--sweep-pos", str(float(pct))]
            # time-units
            tu = (stt.get('time_units') or 'gens')
            if tu:
                cmd += ["--time-units", tu]
            # Only one of sweep-time or fixation-time depending on time_mode
            tm = (stt.get('time_mode') or '').strip()
            if engine == 'msms':
                # msms handles timing differently; skip explicit time flags
                pass
            elif tm == 'Sweep Time':
                try:
                    sv = stt.get('sweep_time', '')
                    if isinstance(sv, str):
                        if sv.strip() != '':
                            cmd += ["--sweep-time", str(float(sv))]
                    elif sv is not None:
                        cmd += ["--sweep-time", str(float(sv))]
                except Exception:
                    pass
            elif tm == 'Fixation Time':
                try:
                    fv = stt.get('fixation_time', '')
                    if isinstance(fv, str):
                        if fv.strip() != '':
                            cmd += ["--fixation-time", str(float(fv))]
                    elif fv is not None:
                        cmd += ["--fixation-time", str(float(fv))]
                except Exception:
                    pass
    except Exception:
        pass

    # return command, no errors, and artifacts mapping
    return cmd, [], artifacts

@app.route("/api/show_engine", methods=["POST"])
def api_show_engine():
    sid=get_sid(); stt=get_state(sid)
    # Include progress flag in the copyable CLI when the user requested it
    cmd, errs, _arts = build_cmd(stt, include_progress=bool(stt.get('progress_flag', False)))
    if errs:
        return jsonify({"ok": False, "errors": errs}), 400
    try:
        proc = subprocess.run(cmd+["--show-command"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode!=0:
            return jsonify({"ok": False, "errors": [proc.stderr.strip() or "Simulator error"]}), 500
        return jsonify({"ok": True, "engine": proc.stdout})
    except FileNotFoundError:
        return jsonify({"ok": False, "errors": ["simulator.py not found"]}), 404


def _parse_progress_line(line):
    # return dict with optional keys: core_pct, pn_pct, a, b
    out = {}
    if not line: return out
    # percent like ' 23%'
    m = re.search(r"(\d{1,3})%", line)
    if m:
        try:
            out['core_pct'] = int(m.group(1))
        except:
            pass
    # a/b pattern
    m2 = re.search(r"(\d+)\s*/\s*(\d+)", line)
    if m2:
        try:
            a = int(m2.group(1)); b = int(m2.group(2))
            out['a'] = a; out['b'] = b
            if b>0:
                out['pn_pct'] = int((a*100)/b)
        except:
            pass
    return out


def _runner_thread(sid, proc, q, run_dir, cancel_path):
    # read stdout lines and push parsed progress events
    try:
        for raw in iter(proc.stdout.readline, ''):
            if raw is None:
                break
            line = raw.rstrip('\n')
            _queue_put_safe(q, {'type': 'line', 'line': line})
            # check cancel sentinel
            if os.path.exists(cancel_path):
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                except Exception:
                    try:
                        proc.terminate()
                    except Exception:
                        pass
            # small sleep to yield
        proc.wait()
    except Exception as e:
        _queue_put_safe(q, {'type': 'error', 'msg': str(e)})
    finally:
        rc = None
        try:
            rc = proc.returncode
        except Exception:
            rc = None
    _queue_put_safe(q, {'type': 'done', 'returncode': rc})


@app.route('/api/execute', methods=['POST'])
def api_execute():
    sid = get_sid(); stt = get_state(sid)
    data = request.get_json(silent=True) or {}
    # allow client to pre-send simple updates
    for k,v in data.items():
        stt[k]=v

    # prepare run dir
    run_dir = _ensure_run_dir(sid)
    cancel_path = os.path.join(run_dir, 'CANCEL')

    # Server-side validation parity with Streamlit UI gating
    errors = _collect_validation_errors(stt)

    if errors:
        return jsonify({'ok': False, 'errors': errors}), 400

    cmd, errs, artifacts = build_cmd(stt, include_progress=True, run_dir=run_dir)
    if errs:
        return jsonify({'ok': False, 'errors': errs}), 400

    # ensure artifacts dict updated and stored in state
    stt['_last_run_dir'] = run_dir
    stt['_last_artifacts'] = artifacts

    # launch subprocess
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, preexec_fn=os.setsid)
    except FileNotFoundError:
        return jsonify({'ok': False, 'errors': ['simulator.py not found or failed to start']}), 500
    except Exception as e:
        return jsonify({'ok': False, 'errors': [str(e)]}), 500

    q = queue.Queue()
    th = threading.Thread(target=_runner_thread, args=(sid, proc, q, run_dir, cancel_path), daemon=True)
    th.start()

    # store runner
    with LOCK:
        # if existing runner, mark it
        RUNNERS[sid] = {'proc': proc, 'thread': th, 'queue': q, 'run_dir': run_dir, 'cancel_path': cancel_path}

    # respond with run info and expected artifact download URLs under /artifacts/<run_dir>/<file>
    download_urls = {}
    for k, p in artifacts.items():
        # if path is under artifacts/, present as artifacts/<run_dir>/<basename>
        fname = os.path.basename(p)
        download_urls[k] = f"/artifacts/{os.path.basename(run_dir)}/{fname}"

    return jsonify({'ok': True, 'run_dir': os.path.basename(run_dir), 'artifacts': download_urls})


@app.route('/api/cancel', methods=['POST'])
def api_cancel():
    sid = get_sid();
    runner = RUNNERS.get(sid)
    if not runner:
        return jsonify({'ok': False, 'error': 'No active run'}), 400
    # create cancel sentinel
    try:
        open(runner['cancel_path'], 'w').close()
    except Exception:
        pass
    # try terminate
    proc = runner.get('proc')
    try:
        if proc and proc.poll() is None:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            except Exception:
                proc.terminate()
            # give short time then kill children
            time.sleep(0.5)
            if proc.poll() is None:
                try:
                    if psutil:
                        pg = psutil.Process(proc.pid)
                        for ch in pg.children(recursive=True):
                            ch.kill()
                        pg.kill()
                    else:
                        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except Exception:
                    pass
    except Exception:
        pass
    return jsonify({'ok': True})


@app.route('/api/progress')
def api_progress():
    sid = get_sid();
    runner = RUNNERS.get(sid)
    if not runner:
        return jsonify({'ok': False, 'error': 'No active run'}), 400

    def _list_run_artifacts(run_dir):
        out = {}
        try:
            if os.path.isdir(run_dir):
                base = os.path.basename(run_dir)
                for nm in os.listdir(run_dir):
                    fp = os.path.join(run_dir, nm)
                    if os.path.isfile(fp):
                        out[nm] = f"/artifacts/{base}/{nm}"
                    elif os.path.isdir(fp):
                        # Recurse one level into subdirectories
                        try:
                            for nm2 in os.listdir(fp):
                                fp2 = os.path.join(fp, nm2)
                                if os.path.isfile(fp2):
                                    rel = f"{nm}/{nm2}"
                                    out[rel] = f"/artifacts/{base}/{rel}"
                        except Exception:
                            pass
        except Exception:
            pass
        return out

    def gen():
        q = runner['queue']
        last_core = 0; last_pn = 0; start = time.time(); last_line_time = start
        run_dir = runner.get('run_dir')
        # Track previously sent artifact names to avoid duplicates
        sent_artifacts = set()
        # Immediately send initial artifact snapshot if any
        try:
            arts = _list_run_artifacts(run_dir)
            if arts:
                pay = { 'artifacts': arts }
                sent_artifacts.update(arts.keys())
                yield f"event: artifact\ndata: {json.dumps(pay)}\n\n"
        except Exception:
            pass
        while True:
            try:
                item = q.get(timeout=0.5)
            except queue.Empty:
                # heartbeat
                # on heartbeat also attempt a lightweight artifact scan
                try:
                    arts = _list_run_artifacts(run_dir)
                    # send only new files
                    new_items = {k:v for k,v in arts.items() if k not in sent_artifacts}
                    if new_items:
                        sent_artifacts.update(new_items.keys())
                        yield f"event: artifact\ndata: {json.dumps({'artifacts': new_items})}\n\n"
                except Exception:
                    pass
                yield 'event: heartbeat\ndata: {}\n\n'
                continue
            if item.get('type') == 'line':
                parsed = _parse_progress_line(item.get('line'))
                now = time.time()
                last_line_time = now
                core = parsed.get('core_pct')
                pn = parsed.get('pn_pct')
                eta = None
                if core is not None and core>0:
                    elapsed = now - start
                    try:
                        eta = int((elapsed*(100-core)/core))
                    except Exception:
                        eta = None
                    last_core = core
                if pn is not None:
                    last_pn = pn
                payload = {'core_pct': last_core, 'pn_pct': last_pn, 'eta': eta, 'line': item.get('line')}
                yield f"data: {json.dumps(payload)}\n\n"
            elif item.get('type') == 'error':
                yield f"event: error\ndata: {json.dumps({'msg': item.get('msg')})}\n\n"
            elif item.get('type') == 'done':
                # final artifact snapshot
                try:
                    arts = _list_run_artifacts(run_dir)
                    new_items = {k:v for k,v in arts.items() if k not in sent_artifacts}
                    if new_items:
                        yield f"event: artifact\ndata: {json.dumps({'artifacts': new_items})}\n\n"
                except Exception:
                    pass
                payload = {'done': True, 'returncode': item.get('returncode', None)}
                yield f"event: done\ndata: {json.dumps(payload)}\n\n"
                break

    return app.response_class(stream_with_context(gen()), mimetype='text/event-stream')


@app.route('/api/validate', methods=['POST'])
def api_validate():
    """Lightweight validation endpoint: runs build_cmd to surface any builder
    errors without invoking the simulator. Returns JSON {ok: bool, errors: []}.
    """
    sid = get_sid(); stt = get_state(sid)
    # Use builder errors and additional UI gating errors
    cmd, errs, _arts = build_cmd(stt, include_progress=False)
    ui_errs = _collect_validation_errors(stt)
    all_errs = (errs or []) + (ui_errs or [])
    if all_errs:
        return jsonify({'ok': False, 'errors': all_errs})
    return jsonify({'ok': True, 'errors': []})

@app.route("/static/<path:p>")
def static_files(p):
    return send_from_directory("static", p)

@app.route("/artifacts/<path:p>")
def artifacts(p):
    # Serve artifacts from the absolute directory next to app.py
    return send_from_directory(ARTIFACTS_DIR, p, as_attachment=True)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=PORT, debug=False)
