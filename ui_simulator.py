# streamlit run ui_simulator.py 
import sys
import subprocess
import shlex
import re
import os
import shutil
import tempfile
import streamlit as st
import threading
import time
import atexit
import signal
import math
try:
    import psutil
except Exception:
    psutil = None

# Optional species/model catalog via stdpopsim
try:
    import stdpopsim as sps
except Exception:
    sps = None

st.set_page_config(page_title="Genomic Simulator UI", page_icon="🧬", layout="wide")

# ---- Minimal style polish ----
st.markdown(
    '''
    <style>
    .block-container {padding-top: 2rem; padding-bottom: 2rem;}
    .stTitle {font-weight: 800;}
    .card {
        padding: 1rem 1.25rem;
        /* make card background transparent and remove any border/shadow so it doesn't render as a thin line */
        box-shadow: none !important;
        border: none !important;
        background: transparent !important;
    }
    .muted {color: #6b7280;}
    .pill {display:inline-block;padding:2px 8px;border-radius:999px;background:#f3f4f6;font-size:12px;margin-left:6px;}
    </style>
    ''',
    unsafe_allow_html=True,
)

st.title("🧬 Genomic Simulator")

# simulator script is assumed to be the default filename in CWD
SIM_SCRIPT = "simulator.py"

# How long generated temporary artifacts should live (seconds).
# We keep this conservative (10 minutes) and run a background cleaner every
# CLEANUP_INTERVAL seconds to remove expired temp files. Only artifacts that
# we created in the system temp dir (zipped archives) are removed by the
# automatic cleaner. User-specified output files in the working directory are
# left alone unless the user explicitly deletes them via the UI.
CLEANUP_TTL = 60 * 10
CLEANUP_INTERVAL = 30

# Toggle to enable developer diagnostics and console prints. Set to False for
# normal user-facing UI to keep the interface clean.
SHOW_DEV_DIAG = False

def ui_diag(msg: str):
    """Emit developer diagnostics to stdout only when SHOW_DEV_DIAG is True."""
    try:
        if SHOW_DEV_DIAG:
            print(msg, flush=True)
    except Exception:
        pass


def _start_cleanup_thread():
    """Start a background thread that periodically deletes expired temp artifacts.

    The registry stored in globals()['_generated_files_registry'] maps absolute
    file paths -> expiry_timestamp (float). The worker removes files whose
    expiry is in the past. The thread is idempotent (only one worker is started).
    """
    try:
        if globals().get('_cleanup_thread_started'):
            return
        globals()['_cleanup_thread_started'] = True
    except Exception:
        pass

    def _worker():
        while True:
            try:
                now = time.time()
                reg = globals().get('_generated_files_registry', {}) or {}
                for p, expiry in list(reg.items()):
                    try:
                        if expiry is not None and now >= expiry:
                            # Only remove paths that still exist
                            if os.path.isfile(p):
                                try:
                                    os.remove(p)
                                except Exception:
                                    pass
                                # remove from registry
                                reg.pop(p, None)
                            elif os.path.isdir(p):
                                try:
                                    shutil.rmtree(p)
                                except Exception:
                                    pass
                                reg.pop(p, None)
                    except Exception:
                        # on any error with a given path, drop it from registry
                        try:
                            reg.pop(p, None)
                        except Exception:
                            pass
                # write back registry (in case it was mutated)
                globals()['_generated_files_registry'] = reg
            except Exception:
                pass
            time.sleep(CLEANUP_INTERVAL)

    t = threading.Thread(target=_worker, daemon=True)
    t.start()


# Ensure the background cleaner is running (idempotent)
try:
    _start_cleanup_thread()
except Exception:
    pass

# ---- helpers ----
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
        # Prefer contigs (newer stdpopsim), but fall back to chromosomes (older/stable API)
        contigs = []
        if genome is not None:
            contigs = getattr(genome, "contigs", []) or []
            # fallback to genome.chromosomes which provide .id on each entry
            if not contigs:
                contigs = getattr(genome, "chromosomes", []) or []

        names = []
        for c in contigs:
            # entries may be objects with .id or .name, or strings
            if isinstance(c, str):
                nm = c
            else:
                nm = getattr(c, "id", None) or getattr(c, "name", None)
            if nm:
                names.append(nm)
        return names
    except Exception:
        return []

def normalize_contig_names(contigs):
    """Return a list of contig labels with 'chr' prefix removed where applicable.
    Keep order and drop duplicates after normalization."""
    # Normalize (remove chr prefix), uppercase special contigs, dedupe and natural-sort
    def norm(nm: str) -> str:
        if not isinstance(nm, str):
            return str(nm)
        low = nm.lower()
        label = nm
        if low.startswith('chr'):
            rest = nm[3:]
            # Uppercase X/Y/MT/M; leave numeric as-is
            if rest.upper() in {'X', 'Y', 'MT', 'M'}:
                label = rest.upper()
            else:
                label = rest
        return label

    seen = set()
    normalized = []
    for nm in contigs:
        label = norm(nm)
        if label not in seen:
            seen.add(label)
            normalized.append(label)

    # natural sort: numeric contigs first (by integer), then X, Y, MT, then others alphabetically
    def sort_key(s: str):
        s_up = s.upper()
        if s.isdigit():
            return (0, int(s))
        if s_up in {'X'}:
            return (1, 0)
        if s_up in {'Y'}:
            return (1, 1)
        if s_up in {'MT', 'M'}:
            return (1, 2)
        # fallback: put after numbered and special, sort by lowercase
        return (2, s.lower())

    return sorted(normalized, key=sort_key)

def list_pop_names(sp, model_id):
    try:
        model = sp.get_demographic_model(model_id)
    except Exception:
        return []
    # demes graph first
    try:
        graph = getattr(model, "model", None)
        if graph and hasattr(graph, "demes"):
            return [d.name for d in graph.demes]
    except Exception:
        pass
    # fallback
    try:
        pops = getattr(model, "populations", [])
        names = []
        for p in pops:
            nm = getattr(p, "name", None) or getattr(p, "id", None)
            if nm:
                names.append(nm)
        if names:
            return names
    except Exception:
        pass
    return []

# species dropdown showing only species with >=1 model
st.subheader("1) Species")
if not sps:
    st.error("`stdpopsim` is required for species/model discovery. Install it to enable the dropdowns.")
    st.stop()

species_items = []
for sp in sps.all_species():
    mids = list_models(sp)
    if mids:
        label_parts = [sp.id]
        sci = getattr(sp, "name", None)
        if sci:
            label_parts.append(f"— {sci}")
        common = getattr(sp, "common_name", None)
        if common:
            label_parts.append(f"({common})")
        species_items.append((" ".join(label_parts), sp.id))
species_items.sort(key=lambda x: x[0].lower())
labels = [lbl for (lbl, sid) in species_items]
# Insert a placeholder so no species is selected by default
placeholder_species = "-- select species --"
display_species = [placeholder_species] + labels
if len(labels) == 0:
    st.error("No species with demographic models found in stdpopsim. Install or update stdpopsim data.")
    st.stop()

selected_label = st.selectbox("Species (hidden label)", display_species, label_visibility='hidden')
if selected_label == placeholder_species:
    species_id = None
    sp_obj = None
    species_selected = False
else:
    species_id = dict(species_items)[selected_label]
    species_selected = True
    try:
        sp_obj = sps.get_species(species_id)
    except Exception as e:
        st.error(f"Could not load species {species_id}: {e}")
        st.stop()

# model
st.subheader("2) Demographic Model")
# Only list models when a species has been chosen
model_options = list_models(sp_obj) if sp_obj is not None else []
if species_selected and not model_options:
    st.error("No models available for the selected species.")
    st.stop()

# Add placeholder so no model is selected by default and disable until species chosen
placeholder_model = "-- select model --"
display_models = [placeholder_model] + model_options
model_select_disabled = not species_selected


# Render the model selectbox and an "Info" link side-by-side. The Info link
# is a best-effort URL into the stdpopsim catalog that opens in a new tab.
def _make_stdpop_url(species_id: str, model_id: str) -> str:
    """Construct a best-effort stdpopsim catalog URL for a given species/model.

    The docs use anchors like:
      #sec_catalog_<something>_models_<modelname>
    We cannot perfectly reconstruct the <something> for every species, so
    build a reasonable anchor using the species_id and model_id lowercased
    and with unsafe chars removed. If that fails, fall back to the catalog
    root page.
    """
    base = "https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html"
    try:
        s = (species_id or "").lower().replace("/", "_").replace(" ", "")
        # model anchor: keep alphanum and underscores
        import re

        m = (model_id or "").lower()
        m = re.sub(r"[^a-z0-9_]+", "", m)
        if s and m:
            return f"{base}#sec_catalog_{s}_models_{m}"
    except Exception:
        pass
    return base

col_m, col_info = st.columns([7, 2])
with col_m:
    selected_model_label = st.selectbox("Demographic model (hidden label)", display_models, label_visibility='hidden', disabled=model_select_disabled)
    if selected_model_label == placeholder_model:
        model_id = None
        model_selected = False
    else:
        model_id = selected_model_label
        model_selected = True
with col_info:
    # show a more prominent "Demographic Model Info" button when a model is selected
    try:
        if model_id:
            # species_id may be None if placeholder was selected; guard for safety
            url = _make_stdpop_url(species_id or "", model_id)
            # Render a nicer button using HTML/CSS. Keep it accessible and open in new tab.
            btn_html = f"""
            <a href="{url}" target="_blank" rel="noopener noreferrer"
               style="display:inline-block;text-decoration:none;padding:8px 10px;border-radius:8px;
                      background:linear-gradient(180deg,#2563eb,#1e40af);color:#fff;font-weight:600;
                      font-size:13px;border:1px solid rgba(255,255,255,0.12);box-shadow:0 2px 6px rgba(16,24,40,0.12);"
               title="Open stdpopsim docs for {model_id}">Demographic Model Info</a>
            """
            st.markdown(btn_html, unsafe_allow_html=True)
        else:
            # keep vertical alignment tidy when nothing selected
            st.markdown('<div style="height:34px"></div>', unsafe_allow_html=True)
    except Exception:
        # best-effort UI enhancement; silently ignore on failure
        pass

# chromosome – always from available contigs (normalized, no 'chr')
st.subheader("3) Chromosome")
# Compute contigs only if we have a species object; however the chromosome
# picker should remain disabled until a model is selected (sequential flow).
raw_contigs = list_contigs(sp_obj) if sp_obj is not None else []
contigs = normalize_contig_names(raw_contigs) if raw_contigs else []
placeholder_chrom = "-- select chromosome --"
# Always show the placeholder in the selectbox so the UI doesn't look empty.
display_contigs = [placeholder_chrom] + contigs
# Disabled when the model hasn't been selected; also disable when there are no real contigs
chrom_disabled = (not model_selected) or (len(contigs) == 0)
selected_chrom = st.selectbox("Chromosome (hidden label)", display_contigs, label_visibility='hidden', disabled=chrom_disabled)
if selected_chrom == placeholder_chrom:
    chromosome = None
    chromosome_selected = False
else:
    chromosome = selected_chrom
    chromosome_selected = True
if len(contigs) == 0 and species_selected:
    st.info("No chromosomes/contigs available for this species in stdpopsim. The `--chromosome` argument will be omitted.")

# Populations & Sample Sizes
#
# UI change: population selection is represented as a toggle button per
# population (shows the population name; displays a checkmark when selected).
# Selection state is stored in `st.session_state` under keys named
# `pop_tick_{model_key}_{index}` (boolean). Sample sizes are stored under
# `pop_n_{model_key}_{index}` (integer). We preserve these keys for
# backward compatibility with the previous checkbox-based UI.
st.subheader("4) Populations & Sample Sizes")
# Enforce sequential step enabling: populations are only editable once
# Species -> Model -> Chromosome have been selected.
if not (species_selected and model_selected and chromosome_selected):
    st.info("Select previous steps to configure populations.")
    names = []
else:
    pop_names = list_pop_names(sp_obj, model_id)
    if not pop_names:
        st.warning("Could not detect populations from the selected model. You can type them manually below.")
        manual_pops = st.text_input("Manual population list (comma-separated)", value="")
        names = [p.strip() for p in manual_pops.split(",") if p.strip()] if manual_pops.strip() else []
    else:
        names = pop_names

# Wrap populations in a styled card for a cleaner visual and add a short hint
st.markdown('<div class="card">', unsafe_allow_html=True)
# Use sanitized model_id so keys remain stable across species/model switches
model_key = (model_id or "").replace("/", "_")
selected = []

# --- Helper callbacks to safely mutate session_state from button callbacks ---
def _toggle_pop(key: str):
    st.session_state[key] = not st.session_state.get(key, False)

def _inc_pop(key: str):
    try:
        st.session_state[key] = int(st.session_state.get(key, 0) or 0) + 1
    except Exception:
        st.session_state[key] = 1

def _dec_pop(key: str):
    try:
        st.session_state[key] = max(0, int(st.session_state.get(key, 0) or 0) - 1)
    except Exception:
        st.session_state[key] = 0


def _reset_number_input(key: str):
    try:
        st.session_state[key] = 0
    except Exception:
        pass


def _set_paired_default(paired_name_key: str, paired_key: str, model_key: str):
    """When the 'Run paired neutral' checkbox is turned on, set a sensible
    default name for the paired neutral artifact based on the output path
    or SFS output. Do not override a user-provided name.
    """
    try:
        # Only set when checkbox is True
        if not st.session_state.get(paired_key, False):
            return

        cur = st.session_state.get(paired_name_key, "") or ""
        if cur.strip():
            # user already provided a name -> do not override
            return

        out_path_key = f"output_path_{model_key}"
        out_val = st.session_state.get(out_path_key, "") or ""
        run_sfs_only_key = f"run_sfs_only_{model_key}"
        run_sfs_only = bool(st.session_state.get(run_sfs_only_key, False))
        sfs_key = f"sfs_output_{model_key}"
        sfs_val = st.session_state.get(sfs_key, "") or ""

        default_neutral_name = ""
        if out_val.strip():
            base = out_val
            # detect double extensions like .ms.gz or .vcf.gz
            ext = ""
            for e in [".ms.gz", ".vcf.gz"]:
                if base.endswith(e):
                    base_noext = base[:-len(e)]
                    ext = e
                    break
            else:
                base_noext, ext = os.path.splitext(base)[0], os.path.splitext(base)[1]

            if run_sfs_only:
                default_neutral_name = base_noext + "_neutral"
            else:
                default_neutral_name = base_noext + "_neutral" + (ext or "")
        elif run_sfs_only and sfs_val.strip():
            base_noext = os.path.splitext(sfs_val)[0]
            default_neutral_name = base_noext + "_neutral"

        if default_neutral_name:
            st.session_state[paired_name_key] = default_neutral_name
    except Exception:
        # best effort; don't crash the UI
        pass


# Render populations in a compact grid: a toggle button (population name) on the left and
# a small Samples input with -/+ buttons on the right. We maintain the original
# session_state keys so switching from checkbox to button is stateful and backward-compatible.
num_cols = max(1, min(3, len(names))) if names else 1
cols = st.columns(num_cols)
for idx, pname in enumerate(names):
    col = cols[idx % num_cols]
    with col:
        left, right = st.columns([3, 1])
        tick_key = f"pop_tick_{model_key}_{idx}"
        num_key = f"pop_n_{model_key}_{idx}"

        # initialize keys if missing (back-compat with previous checkbox-driven UI)
        if tick_key not in st.session_state:
            st.session_state[tick_key] = False
        if num_key not in st.session_state:
            st.session_state[num_key] = 0

        # Left: population toggle button. Use markdown-styled button to show state.
        # We store selection as a boolean in session_state under tick_key.
        # Use three columns inside this population cell: toggle | number | +/-
        col_left, col_mid, col_right = st.columns([3, 1, 1])
        with col_left:
            is_selected = st.session_state.get(tick_key, False)
            # prettier dot indicator: green when selected, hollow when not
            dot = "🟢" if is_selected else "⚪"
            btn_label = f"{dot} {pname}"
            # Use on_click callback to toggle selection safely
            st.button(btn_label, key=f"pop_btn_{model_key}_{idx}", on_click=_toggle_pop, args=(tick_key,))

        # current numeric value (ensure integer)
        cur_val = int(st.session_state.get(num_key, 0) or 0)
        with col_mid:
            # compact numeric input: bind to the same session key
            # Avoid passing `value=` when the key already exists in session_state
            if num_key in st.session_state:
                num = st.number_input("Samples (hidden)", min_value=0, step=1, key=num_key, label_visibility='hidden')
            else:
                num = st.number_input("Samples (hidden)", min_value=0, value=cur_val, step=1, key=num_key, label_visibility='hidden')
        with col_right:
            dec_key = f"pop_dec_{model_key}_{idx}"
            inc_key = f"pop_inc_{model_key}_{idx}"
            # Increase on top, decrease below for a natural vertical flow
            st.button("▲", key=inc_key, on_click=_inc_pop, args=(num_key,))
            st.button("▼", key=dec_key, on_click=_dec_pop, args=(num_key,))

        # append tuple for downstream validation
        cur_tick = st.session_state.get(tick_key, False)
        cur_num = int(st.session_state.get(num_key, 0) or 0)
        selected.append((pname, cur_tick, cur_num))

        # Inline validation: if selected but sample size is zero, show a subtle inline message
        if cur_tick and cur_num <= 0:
            # place an inline caption under the row (muted style) to avoid blocking the UI
            st.markdown(f"<div style='color:#b91c1c;font-size:13px;margin-top:4px;'>Population '{pname}' is selected but sample size is 0 — enter a positive number.</div>", unsafe_allow_html=True)

ordered_ready = False
pop_order = []
counts = []
if selected:
    # New logic: allow any contiguous block of populations to be selected (no gaps),
    # but the block may start at any index. This replaces the previous prefix-only
    # requirement which forced the first population to be selected.
    #
    # Rationale: users may want to simulate a subset of populations that doesn't
    # include the first population (for example, selecting populations 2+3 only).
    # Validate that the chosen populations form a contiguous segment and have
    # positive sample sizes.
    valid = True
    tick_indices = [i for i, (_pname, tick, _num) in enumerate(selected) if tick]
    any_tick = len(tick_indices) > 0

    # Selected indices must form a contiguous range
    if any_tick:
        if (max(tick_indices) - min(tick_indices) + 1) != len(tick_indices):
            valid = False

    # All selected populations must have positive sample counts
    if any(num <= 0 for (_pname, tick, num) in selected if tick):
        valid = False

    if valid and any_tick:
        ordered_ready = True
        for pname, tick, num in selected:
            if tick:
                pop_order.append(pname)
                counts.append(num)
    else:
        # Warning suppressed by user request: do not show the 'Set populations in order...' message
        pass

# Core & Engine (render above tabs so we can conditionally include Paired neutral tab)
    st.markdown('</div>', unsafe_allow_html=True)

# Core & Engine (render above tabs so we can conditionally include Paired neutral tab)
st.markdown("### Core & Engine")
with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        engine = st.selectbox("Engine", ["discoal","ms","msms","scrm","msprime"], index=2, disabled=not ordered_ready)
        # make replicates session-backed so changes can trigger parallel default reset
        reps_key = f"replicates_{model_key}"
        if reps_key not in st.session_state:
            # Create the widget with an explicit default value; do NOT pre-set
            # the session key via API before creating the widget to avoid
            # the Streamlit warning about defaults vs session state.
            replicates = st.number_input("Replicates", min_value=1, value=1, step=1, disabled=not ordered_ready, key=reps_key)
        else:
            # when the session key exists, do not pass `value=` (Streamlit warns)
            replicates = st.number_input("Replicates", min_value=1, step=1, disabled=not ordered_ready, key=reps_key)
    with c2:
        # Use stable keys so values persist per model and we can reference each other's state
        seq_key = f"seq_length_{model_key}"
        snp_key = f"target_snps_{model_key}"
        chrom_length_key = f"chrom_length_mode_{model_key}"

        # Chromosome Length checkbox: when True we will NOT pass --length or --target-snps
        chrom_length_default = st.session_state.get(chrom_length_key, False)
        # Attempt to resolve the picked chromosome's native length (various stdpopsim versions)
        contig_len = None
        try:
            genome = getattr(sp_obj, "genome", None)
            contigs = []
            if genome is not None:
                contigs = getattr(genome, "contigs", []) or getattr(genome, "chromosomes", []) or []

            def norm_label(c):
                if isinstance(c, str):
                    nm = c
                else:
                    nm = getattr(c, "id", None) or getattr(c, "name", None)
                if not isinstance(nm, str):
                    return str(nm)
                low = nm.lower()
                if low.startswith('chr'):
                    rest = nm[3:]
                    if rest.upper() in {'X', 'Y', 'MT', 'M'}:
                        return rest.upper()
                    else:
                        return rest
                return nm

            if chromosome and contigs:
                for c in contigs:
                    if norm_label(c) == chromosome:
                        contig_len = getattr(c, "length", None) or getattr(c, "sequence_length", None) or getattr(c, "size", None) or getattr(c, "bp_length", None)
                        break
        except Exception:
            contig_len = None

        # Build label: include picked chromosome and length in bp when available
        if chromosome:
            if contig_len:
                try:
                    contig_len_display = f"{int(contig_len):,}"
                except Exception:
                    contig_len_display = str(int(contig_len)) if contig_len is not None else ""
                chrom_label = f"Chromosome {chromosome} Length ({contig_len_display} bp)"
            else:
                chrom_label = f"Chromosome {chromosome} Length"
        else:
            # When no chromosome is selected keep the label short and clean
            chrom_label = "Chromosome Length"

        chrom_length_mode = st.checkbox(
            chrom_label,
            value=chrom_length_default,
            key=chrom_length_key,
            disabled=not ordered_ready,
            help="When enabled the simulator will use the chromosome/contig length and the UI will omit --length and --target-snps from the command."
        )

        # Determine disabled state: length cannot be set if target_snps > 0; target_snps cannot be set if length > 0
        snp_cur = st.session_state.get(snp_key, 0)
        seq_cur = st.session_state.get(seq_key, 0)
        disabled_seq = (not ordered_ready) or (int(snp_cur) > 0) or chrom_length_mode
        disabled_target = (not ordered_ready) or (int(seq_cur) > 0) or chrom_length_mode

        # use session-backed keys for the inputs to avoid Streamlit warnings
        if seq_key in st.session_state:
            seq_length = st.number_input("Length (bp)", min_value=0, step=1, disabled=disabled_seq, key=seq_key)
        else:
            seq_length = st.number_input("Length (bp)", min_value=0, value=0, step=1, disabled=disabled_seq, key=seq_key)

        # show a Reset button for Length when non-zero and not using chromosome native length
        if not chrom_length_mode and int(st.session_state.get(seq_key, 0) or 0) != 0:
            st.button("Reset Length", on_click=_reset_number_input, args=(seq_key,))

        if snp_key in st.session_state:
            target_snps = st.number_input("Target SNPs", min_value=0, step=1, disabled=disabled_target, key=snp_key)
        else:
            target_snps = st.number_input("Target SNPs", min_value=0, value=0, step=1, disabled=disabled_target, key=snp_key)

        # show a Reset button for Target SNPs when non-zero and not using chromosome native length
        if not chrom_length_mode and int(st.session_state.get(snp_key, 0) or 0) != 0:
            st.button("Reset Target SNPs", on_click=_reset_number_input, args=(snp_key,))
    with c3:
        # Only enable tolerance when a positive Target SNPs value is set (and when not using Chromosome Length)
        tol_disabled = (not ordered_ready) or chrom_length_mode or (not target_snps or int(target_snps) <= 0)

        # Use a radio control for Auto vs Custom to be more explicit and prettier
        tol_auto_key = f"target_snps_tol_auto_{model_key}"
        tol_auto_default = st.session_state.get(tol_auto_key, True)
        col_tol_a, col_tol_b = st.columns([1, 3])

        with col_tol_a:
            # Radio with two options: Auto (do not include flag) or Custom (show input)
            tol_choice = st.radio(
                "Tolerance",
                options=["Auto", "Custom"],
                index=0 if tol_auto_default else 1,
                key=f"target_snps_tol_choice_{model_key}",
                horizontal=True,
                disabled=tol_disabled,
            )
            # Update session_state boolean for backward compatibility
            target_snps_tol_auto = (tol_choice == "Auto")
            st.session_state[tol_auto_key] = target_snps_tol_auto

        with col_tol_b:
            # Only render the Custom tolerance input when the user explicitly
            # selected 'Custom' AND tolerance editing is enabled for this run.
            tol_key_ss = f"target_snps_tol_{model_key}"
            target_snps_tol = None
            if not (tol_disabled or target_snps_tol_auto):
                # When the user first switches to Custom, initialize a sensible
                # default of 10% (as requested). Create the widget with an
                # explicit default only when the session key does not exist.
                if tol_key_ss not in st.session_state:
                    target_snps_tol = st.number_input(
                        "Target SNPs tolerance (%)",
                        min_value=1e-9,
                        max_value=100.0,
                        value=10.0,
                        step=0.1,
                        key=tol_key_ss,
                        help="Tolerance for target SNPs overshoot. Enter percent (>0-100). For backward compatibility the simulator also accepts fractional values <=1."
                    )
                else:
                    # When the session key exists, do NOT pass `value=` —
                    # Streamlit will use the session state value and avoids
                    # emitting a warning about setting defaults via Session API.
                    target_snps_tol = st.number_input(
                        "Target SNPs tolerance (%)",
                        min_value=1e-9,
                        max_value=100.0,
                        step=0.1,
                        key=tol_key_ss,
                        help="Tolerance for target SNPs overshoot. Enter percent (>0-100). For backward compatibility the simulator also accepts fractional values <=1."
                    )
            else:
                # When Auto or disabled, do not show the input. Keep a local
                # placeholder variable for downstream code that may read it.
                target_snps_tol = float(st.session_state.get(tol_key_ss, 0.0) or 0.0)

        if tol_disabled:
            st.caption("Tolerance applies only when 'Target SNPs' is set to a positive value.")
        elif st.session_state.get(tol_auto_key, True):
            st.caption("Tolerance: Auto (simulator chooses). Uncheck Auto to set a custom percent.)")

        # Output format and path: use session-backed keys so we can normalize
        out_format_key = f"output_format_{model_key}"
        out_path_key = f"output_path_{model_key}"
        if out_format_key not in st.session_state:
            st.session_state[out_format_key] = "ms"
        if out_path_key not in st.session_state:
            st.session_state[out_path_key] = ""

        def _normalize_output_path():
            fmt = str(st.session_state.get(out_format_key, "")).strip()
            path = st.session_state.get(out_path_key, "") or ""
            if not path or not fmt:
                return
            fmt_ext_map = {
                "ms": ".ms",
                "ms.gz": ".ms.gz",
                "vcf": ".vcf",
                "vcf.gz": ".vcf.gz",
                "bcf": ".bcf",
            }
            desired_ext = fmt_ext_map.get(fmt, "")
            if not desired_ext:
                return
            # If already has desired ext, do nothing
            if path.endswith(desired_ext):
                return
            # Strip any known extension first
            for e in fmt_ext_map.values():
                if path.endswith(e):
                    path = path[:-len(e)]
                    break
            new_path = path + desired_ext
            if new_path != st.session_state.get(out_path_key, ""):
                st.session_state[out_path_key] = new_path

        # Disable output format selection when user chose 'Run SFS only'
        run_sfs_only_key = f"run_sfs_only_{model_key}"
        run_sfs_only_local = bool(st.session_state.get(run_sfs_only_key, False))
        out_format = st.selectbox(
            "Output format",
            ["ms","ms.gz","vcf","vcf.gz","bcf"],
            disabled=(not ordered_ready) or run_sfs_only_local,
            key=out_format_key,
            on_change=_normalize_output_path,
        )
    with c4:
        # Avoid passing `value=` when the session key exists to prevent Streamlit
        # warnings about widgets being created with a default and later set via
        # session state. Use the session-backed key directly.
        # Disable output path input when Run SFS only is enabled
        if out_path_key in st.session_state:
            out_path = st.text_input(
                "Output path",
                disabled=(not ordered_ready) or run_sfs_only_local,
                key=out_path_key,
                on_change=_normalize_output_path,
            )
        else:
            out_path = st.text_input(
                "Output path",
                value=st.session_state.get(out_path_key, ""),
                disabled=(not ordered_ready) or run_sfs_only_local,
                key=out_path_key,
                on_change=_normalize_output_path,
            )
        # Show system temp choice only if system tempdir exists and is writable
        try:
            sys_tmp = tempfile.gettempdir()
            sys_tmp_ok = os.path.isdir(sys_tmp) and os.access(sys_tmp, os.W_OK)
        except Exception:
            sys_tmp = None
            sys_tmp_ok = False
        if sys_tmp_ok:
            temp_loc = st.selectbox("Temp location", ["local","system"], disabled=not ordered_ready)
        else:
            temp_loc = 'local'
    

    st.markdown('<div class="card" style="margin-top:0.75rem;">', unsafe_allow_html=True)
    d1, d2, d3, d4 = st.columns(4)
    with d1:
        # Show a slider from 1 up to either physical cores (for msprime) or
        # logical threads (for other engines). Prefer psutil when available
        # for accurate physical core count; fall back to os.cpu_count().
        try:
            logical_count = psutil.cpu_count(logical=True) if psutil else os.cpu_count()
        except Exception:
            logical_count = os.cpu_count()
        try:
            physical_count = psutil.cpu_count(logical=False) if psutil else None
        except Exception:
            physical_count = None

        # Determine max workers based on engine
        if engine == "msprime":
            engine_max = int(physical_count or logical_count or 1)
        else:
            engine_max = int(logical_count or 1)

        # Cap parallel workers by number of replicates (can't have more workers than work units)
        try:
            reps = int(replicates) if replicates else 1
        except Exception:
            reps = 1

        max_workers = max(1, min(engine_max, reps))

        # session-backed parallel so we can reset when replicates change
        par_key = f"parallel_{model_key}"
        stored_reps_key = f"_stored_reps_{model_key}"
        # default: more than half the available workers (floor(engine_max/2)+1)
        try:
            half_plus = max(1, (engine_max // 2) + 1)
        except Exception:
            half_plus = 1
        # desired default must not exceed engine_max or replicates
        desired_default = max(1, min(half_plus, engine_max, reps))
        try:
            prev_reps = int(st.session_state.get(stored_reps_key, 0) or 0)
        except Exception:
            prev_reps = 0
        # if replicates changed, reset parallel to desired default
        if prev_reps != reps:
            st.session_state[stored_reps_key] = reps
            st.session_state[par_key] = desired_default

        if par_key not in st.session_state:
            st.session_state[par_key] = desired_default

        # If only one worker is possible, show a read-only caption instead of a slider
        if max_workers <= 1:
            parallel = 1
            st.session_state[par_key] = 1
            st.caption("Parallel workers: 1 (computed)")
        else:
            # clamp stored parallel to valid range
            try:
                stored_par = int(st.session_state.get(par_key, 1) or 1)
            except Exception:
                stored_par = 1
            if stored_par < 1:
                stored_par = 1
            if stored_par > max_workers:
                stored_par = max_workers
            # render slider using the session-backed key so its value persists
            # Avoid passing `value=` when the session key exists to prevent
            # Streamlit warnings about widgets created with a default then set
            # via session state. Use the session-backed key directly when present.
            if par_key in st.session_state:
                parallel = st.slider(
                    "Parallel workers",
                    min_value=1,
                    max_value=max_workers,
                    step=1,
                    disabled=not ordered_ready,
                    key=par_key,
                    help=("Number of parallel workers. Uses physical cores for msprime, "
                          "and logical threads for other engines. Capped at the number of replicates.")
                )
            else:
                parallel = st.slider(
                    "Parallel workers",
                    min_value=1,
                    max_value=max_workers,
                    value=stored_par,
                    step=1,
                    disabled=not ordered_ready,
                    key=par_key,
                    help=("Number of parallel workers. Uses physical cores for msprime, "
                          "and logical threads for other engines. Capped at the number of replicates.")
                )
    with d2:
        # Simulation per Work: prettier Auto vs Custom radio and renamed label
        sims_auto_key = f"sims_per_work_auto_{model_key}"
        sims_auto_default = st.session_state.get(sims_auto_key, True)

        # Hide the Simulation-per-Work controls for msprime (not applicable).
        # Ensure session state and local variables are set so later code can
        # safely reference `sims_per_work` and the session key.
        if engine == "msprime":
            # Force Auto mode and define sims_per_work as 0 (meaning Auto)
            st.session_state[sims_auto_key] = True
            sims_per_work_auto = True
            sims_per_work = 0
            # Inform the user briefly
            st.caption("Simulation per Work: not applicable for msprime; using automatic scheduling.")
        else:
            col_sims_a, col_sims_b = st.columns([1, 3])
            with col_sims_a:
                sims_choice = st.radio(
                    "Simulation per Work",
                    options=["Auto", "Custom"],
                    index=0 if sims_auto_default else 1,
                    key=f"sims_per_work_choice_{model_key}",
                    horizontal=True,
                    disabled=not ordered_ready,
                )
                sims_per_work_auto = (sims_choice == "Auto")
                st.session_state[sims_auto_key] = sims_per_work_auto
            with col_sims_b:
                # If Auto is selected we do not set the argument; if Custom is
                # selected, show a slider whose max is ceil(replicates / parallel)
                if sims_per_work_auto:
                    sims_per_work = 0
                else:
                    try:
                        par = int(parallel) if parallel else 1
                    except Exception:
                        par = 1
                    try:
                        reps = int(replicates) if replicates else 1
                    except Exception:
                        reps = 1
                    import math
                    max_val = max(1, math.ceil(reps / max(1, par)))
                    if max_val <= 1:
                        # Slider requires min < max; when only 1 is possible show a caption/read-only value
                        st.caption("Simulation per Work: 1 (computed)")
                        sims_per_work = 1
                    else:
                        sims_per_work = st.slider(
                            "Simulation per Work",
                            min_value=1,
                            max_value=max_val,
                            value=1,
                            step=1,
                            disabled=not ordered_ready,
                            help=f"Workers will process this many simulations per work (1..{max_val}).",
                            key=f"sims_per_work_slider_{model_key}",
                        )
    with d3:
        # Render slider and a small radio selector for cap mode (system/run) side-by-side
        col_ram_slider, col_ram_mode = st.columns([3, 1])
        with col_ram_slider:
            max_ram_percent = st.slider(
                "Max RAM %",
                min_value=1,
                max_value=100,
                value=80,
                step=1,
                disabled=not ordered_ready,
                help="Numeric threshold (1-100). Interpretation depends on Max RAM cap mode."
            )
        with col_ram_mode:
            # session-backed key so selection persists per model
            max_ram_cap_key = f"max_ram_cap_{model_key}"
            if max_ram_cap_key not in st.session_state:
                st.session_state[max_ram_cap_key] = 'system'
            max_ram_cap = st.radio(
                "Cap",
                options=['system', 'run'],
                index=0 if st.session_state.get(max_ram_cap_key, 'system') == 'system' else 1,
                key=max_ram_cap_key,
                horizontal=False,
                disabled=not ordered_ready,
                help="Interpretation: 'system' compares to total system RAM usage; 'run' compares the run's RSS as % of total RAM."
            )
    with d4:
        # Progress flag moved to Build & Run tab per UX change; placeholder here to preserve layout
        st.write(" ")
    

    # Build dynamic tabs list: include Paired neutral only when sweep params enabled and engine supports selection engines
    sweep_key = f"sweep_enable_{model_key}"
    # If engine is not sweep-capable, ensure sweep flag is disabled and hide sweep controls
    if engine not in {"msms", "discoal"}:
        # Force to False to avoid showing Paired neutral when switching engines
        try:
            st.session_state[sweep_key] = False
        except Exception:
            pass
    sweep_enabled_ss = bool(st.session_state.get(sweep_key, False))
    show_paired_tab = bool(sweep_enabled_ss and engine in {"msms", "discoal"})

    # Ensure paired-neutral variables exist with safe defaults even when the
    # "Paired neutral" tab is not shown. build_cmd() references these names
    # later, and Python will raise NameError if they are missing.
    paired_neutral = False
    paired_neutral_name = ""
    neutral_engine = ""

    tabs_labels = ["SFS"]
    if show_paired_tab:
        tabs_labels.append("Paired neutral")
    tabs_labels += ["Demography & Selection", "Build & Run"]
    tabs = st.tabs(tabs_labels)
    tab_map = {name: tab for name, tab in zip(tabs_labels, tabs)}

    # Site Frequency Spectrum (SFS)
    with tab_map["SFS"]:
        st.markdown("### Site Frequency Spectrum (SFS)")
        st.markdown('<div class="card">', unsafe_allow_html=True)
        sfs_key = f"sfs_enable_{model_key}"
        if sfs_key not in st.session_state: 
            st.session_state[sfs_key] = False
        sfs_on = st.checkbox("Compute SFS", value=st.session_state.get(sfs_key, False), disabled=not ordered_ready, key=sfs_key)
        # session-backed SFS output so we can auto-fill it from the Output path
        sfs_output_key = f"sfs_output_{model_key}"
        if sfs_output_key not in st.session_state:
            st.session_state[sfs_output_key] = ""

        # If sfs is empty but output path exists, derive a sensible default
        try:
            out_path_key = f"output_path_{model_key}"
            cur_sfs = st.session_state.get(sfs_output_key, "") 
            cur_out = st.session_state.get(out_path_key, "")
            if (not cur_sfs or not cur_sfs.strip()) and cur_out and cur_out.strip():
                base = cur_out
                for e in [".ms.gz", ".vcf.gz"]:
                    if base.endswith(e):
                        base = base[:-len(e)]
                        break
                else:
                    base = os.path.splitext(base)[0]
                st.session_state[sfs_output_key] = base + ".sfs"
        except Exception:
            pass

        # Ensure SFS output ends with .sfs when the user types a name without extension.
        def _ensure_sfs_ext():
            cur = st.session_state.get(sfs_output_key, "") or ""
            if cur and not cur.lower().endswith('.sfs'): 
                st.session_state[sfs_output_key] = cur + '.sfs'

        # Avoid passing `value=` when the session key already exists (prevents Streamlit warning)
        if sfs_output_key in st.session_state:
            sfs_output = st.text_input(
                "SFS output",
                disabled=not ordered_ready,
                key=sfs_output_key, 
                on_change=_ensure_sfs_ext,
                help="SFS output path. If left empty the UI derives a filename from the Output path (replaces extension with .sfs)."
            )
        else:
            sfs_output = st.text_input(
                "SFS output",
                value=st.session_state.get(sfs_output_key, ""),
                disabled=not ordered_ready,
                key=sfs_output_key,
                on_change=_ensure_sfs_ext,
                help="SFS output path. If left empty the UI derives a filename from the Output path (replaces extension with .sfs)."
            )
        sfs_normalized = st.checkbox("Normalize SFS", value=False, disabled=not ordered_ready)
        sfs_mode = st.selectbox("SFS mode", ["mean","per-rep"], disabled=not ordered_ready)
        

    # Paired neutral tab (conditionally present in tab_map)
    if "Paired neutral" in tab_map:
        with tab_map["Paired neutral"]:
            st.markdown("### Paired Neutral output")
            st.markdown('<div class="card">', unsafe_allow_html=True)
            paired_key = f"paired_neutral_{model_key}"
            paired_neutral = st.checkbox( 
                "Run paired neutral",
                value=False,
                disabled=not ordered_ready,
                key=paired_key,
                on_change=_set_paired_default,
                args=(f"paired_neutral_name_{model_key}", paired_key, model_key),
            )
            # Suggest a default name based on the primary output path or SFS output.
            # If primary output exists, insert _neutral before its extension. If SFS-only,
            # suggest the sfs base with _neutral and do not append an extension (per request).
            paired_name_key = f"paired_neutral_name_{model_key}"
            default_neutral_name = "" 
            try:
                out_path_key = f"output_path_{model_key}"
                out_val = st.session_state.get(out_path_key, "") or ""
                run_sfs_only_key = f"run_sfs_only_{model_key}"
                run_sfs_only = bool(st.session_state.get(run_sfs_only_key, False))
                sfs_key = f"sfs_output_{model_key}"
                sfs_val = st.session_state.get(sfs_key, "") or ""

                if out_val.strip():
                    base = out_val
                    # detect double extensions first
                    ext = "" 
                    for e in [".ms.gz", ".vcf.gz"]:
                        if base.endswith(e):
                            base_noext = base[:-len(e)]
                            ext = e
                            break
                    else:
                        base_noext, ext = os.path.splitext(base)[0], os.path.splitext(base)[1]

                    # Always suggest a paired-neutral base name without file extensions.
                    # The simulator --paired-neutral argument expects a base name; do not
                    # append engine-specific extensions here.
                    default_neutral_name = base_noext + "_neutral"
                elif run_sfs_only and sfs_val.strip():
                    base_noext = os.path.splitext(sfs_val)[0]
                    default_neutral_name = base_noext + "_neutral"
            except Exception:
                default_neutral_name = ""

            # Track whether the user explicitly edited the paired-neutral name.
            user_set_key = f"paired_neutral_name_user_set_{model_key}"
            if user_set_key not in st.session_state:
                st.session_state[user_set_key] = False

            def _mark_pn_user_set(): 
                # When the user edits the text input, mark the name as user-set so
                # automatic updates no longer overwrite it.
                try:
                    st.session_state[user_set_key] = True
                except Exception:
                    pass

            # If the user hasn't set a custom name, keep the session value in sync
            # with the derived default (changes when output/sfs/run_sfs_only change).
            # Only overwrite when the user-set flag is False.
            if not st.session_state.get(user_set_key, False):
                # Initialize or update the session with the current default
                if paired_name_key not in st.session_state:
                    st.session_state[paired_name_key] = default_neutral_name
                else:
                    # update only when default differs from the stored value
                    try:
                        if (st.session_state.get(paired_name_key, "") or "") != (default_neutral_name or ""):
                            st.session_state[paired_name_key] = default_neutral_name
                    except Exception:
                        pass

            # Text input: when the user edits it, _mark_pn_user_set will mark it as user-set
            if paired_name_key in st.session_state:
                paired_neutral_name = st.text_input("Paired Neutral output", disabled=not (ordered_ready and paired_neutral), key=paired_name_key, on_change=_mark_pn_user_set)
            else:
                paired_neutral_name = st.text_input("Paired Neutral output", value=default_neutral_name, disabled=not (ordered_ready and paired_neutral), key=paired_name_key, on_change=_mark_pn_user_set)
            # Neutral engine is optional (empty => use same engine as sweep). Only enable
            # when user ticked sweep parameters AND the selected engine supports sweeps
            # (msms or discoal) and paired_neutral is checked.
            neutral_engine_options = ["", "discoal", "ms", "msms", "scrm", "msprime"]
            # Determine if sweep is enabled via session_state (sweep_enable widget has its own key)
            sweep_key = f"sweep_enable_{model_key}"
            sweep_enabled_ss = bool(st.session_state.get(sweep_key, False))
            neutral_engine_enabled = bool(ordered_ready and st.session_state.get(paired_key, False) and sweep_enabled_ss and engine in {"msms", "discoal"})
            # Default the neutral engine to the current sweep/core engine when
            # available. Persist choice in session_state using a model-scoped key.
            neutral_engine_key = f"neutral_engine_{model_key}"
            try:
                default_index = neutral_engine_options.index(engine) if engine in neutral_engine_options else 0
            except Exception:
                default_index = 0
            neutral_engine = st.selectbox(
                "Neutral engine",
                neutral_engine_options,
                index=default_index,
                key=neutral_engine_key,
                disabled=not neutral_engine_enabled,
                help="Leave empty to use the same engine as the sweep. Enabled only when sweep params are on and engine is msms/discoal."
            )
            if not neutral_engine_enabled:
                if not sweep_enabled_ss:
                    st.caption("Enable sweep parameters to choose a neutral engine (defaults to sweep engine).")
                elif engine not in {"msms", "discoal"}:
                    st.caption("Paired neutral engine selection is meaningful mainly for msms or discoal engines.")
            

    with tab_map["Demography & Selection"]:
        st.markdown("### Demography & Selection")
        st.markdown('<div class="card">', unsafe_allow_html=True)
        e1, e2 = st.columns(2)
        with e1:
            no_en_ladder = st.checkbox("Disable growth discretization", value=False, disabled=not ordered_ready)
            growth_max_fold = st.number_input("Growth max fold", min_value=1.0, value=1.05, step=0.01, disabled=not ordered_ready)
        with e2:
            # Allow user to enable sweep/selection parameters. When enabled, the UI
            # reveals the sweep population (dropdown of model populations) and other params.
            # Use a stable key so other tabs can read the value from session_state.
            sweep_key = f"sweep_enable_{model_key}"
            # Only show sweep enable when engine supports sweep simulations
            if engine in {"msms", "discoal"}:
                sweep_enable = st.checkbox("Enable sweep/selection parameters", value=False, disabled=not ordered_ready, key=sweep_key)
            else:
                # Ensure session flag off and present explanatory caption
                try:
                    st.session_state[sweep_key] = False
                except Exception:
                    pass
                sweep_enable = False
                st.caption("Sweep/selection parameters are available only for msms or discoal engines.")
            if not ordered_ready:
                st.caption("Set populations in order and sample sizes (step 4) to enable sweep/selection parameters.")
            if sweep_enable:
                if names:
                    # Do not pre-select a population. Insert a placeholder
                    # option and require the user to pick one explicitly.
                    sweep_choices = [""] + names
                    sweep_label_map = {"": "-- select population --"}
                    # Render display labels so the blank option is user-friendly
                    display_labels = [sweep_label_map.get(x, x) for x in sweep_choices]
                    sel_index = 0
                    sweep_sel = st.selectbox("Sweep population", display_labels, index=sel_index, disabled=not ordered_ready, key=f"sweep_pop_{model_key}_sel")
                    # Map back to actual population name (empty string means none chosen)
                    sweep_pop = "" if sweep_sel == sweep_label_map.get("") else sweep_sel
                else:
                    sweep_pop = st.text_input("Sweep population", value="", disabled=not ordered_ready, key=f"sweep_pop_manual_{model_key}")
                # Sweep position: show a bp slider when chromosome native length
                # or explicit length is available, otherwise show percentage (0–100).
                # Internally we keep `sweep_pos_pct` as percent (0..100) so build_cmd
                # can convert to fractional 0..1 as before.
                length_for_pos = None
                try:
                    # prefer explicit sequence length when set
                    if 'seq_length' in globals() or 'seq_length' in locals():
                        # seq_length may be int or 0
                        sl = int(st.session_state.get(f"seq_length_{model_key}", 0) or 0)
                        if sl and sl > 0:
                            length_for_pos = sl
                except Exception:
                    length_for_pos = None
                # fallback to contig/native length if chromosome length mode enabled
                try:
                    if not length_for_pos and chrom_length_mode:
                        # contig_len is computed earlier when resolving chromosome; if present use it
                        if 'contig_len' in globals() and contig_len:
                            length_for_pos = int(contig_len)
                except Exception:
                    pass

                if length_for_pos and length_for_pos > 0:
                    default_bp = int(min(max(0, length_for_pos // 2), length_for_pos))
                    sweep_pos_bp = st.slider(
                        "Sweep position (bp)",
                        min_value=0,
                        max_value=length_for_pos,
                        value=default_bp,
                        step=max(1, int(length_for_pos // 1000) or 1),
                        disabled=not ordered_ready,
                        key=f"sweep_pos_{model_key}",
                    )
                    # convert bp to percent for internal consistency
                    try:
                        sweep_pos_pct = float(sweep_pos_bp) / float(length_for_pos) * 100.0
                    except Exception:
                        sweep_pos_pct = 50.0
                else:
                    sweep_pos_pct = st.slider(
                        "Sweep position (%)",
                        min_value=0.0,
                        max_value=100.0,
                        value=50.0,
                        step=0.1,
                        disabled=not ordered_ready,
                        key=f"sweep_pos_{model_key}",
                    )
                sel_s = st.number_input("Selection Coefficient (s)", value=0.1, step=0.0001, format="%g", disabled=not ordered_ready, key=f"sel_s_{model_key}")
                # Engine-specific sweep times and time-units: render time units
                # inline next to the relevant time field so the UI mirrors the
                # placement of "Cap" next to Max RAM.
                sweep_time = None
                fixation_time = None
                if ordered_ready and engine in {"discoal", "msms"}:
                    # use two columns: time input (wider) and units selector (narrower)
                    col_time, col_units = st.columns([3, 2])
                    if engine == "discoal":
                        with col_time:
                            sweep_time = st.number_input("Sweep time", min_value=0.0, value=0.0, step=1.0, disabled=not ordered_ready, key=f"sweep_time_{model_key}")
                        with col_units:
                            # Compact vertical radio similar to "Cap" control
                            tu_key = f"time_units_{model_key}"
                            if tu_key not in st.session_state:
                                st.session_state[tu_key] = 'gens'
                            time_units = st.radio(
                                "Units",
                                ["gens", "4N"],
                                index=0 if st.session_state.get(tu_key, 'gens') == 'gens' else 1,
                                key=tu_key,
                                disabled=not ordered_ready,
                            )
                    elif engine == "msms":
                        with col_time:
                            # default fixation time as 0.0 (float) to keep numeric types consistent
                            fixation_time = st.number_input("Fixation time", min_value=0.0, value=0.0, step=1.0, disabled=not ordered_ready, key=f"fixation_time_{model_key}")
                        with col_units:
                            # Compact vertical radio similar to "Cap" control
                            tu_key = f"time_units_{model_key}"
                            if tu_key not in st.session_state:
                                st.session_state[tu_key] = 'gens'
                            time_units = st.radio(
                                "Units",
                                ["gens", "4N"],
                                index=0 if st.session_state.get(tu_key, 'gens') == 'gens' else 1,
                                key=tu_key,
                                disabled=not ordered_ready,
                            )
                else:
                    # When sweep params are disabled or engine not sweep-capable,
                    # default the units to 'gens' so build_cmd has a sensible value.
                    time_units = "gens"
            else:
                # default placeholders when sweep params disabled
                sweep_pop = ""
                # default in percent (50%) -> converted later
                sweep_pos_pct = 50.0
                sel_s = 0.0
                time_units = "gens"
        


    with tab_map["Build & Run"]:
        st.markdown("### Build & Run")
        st.markdown('<div class="card">', unsafe_allow_html=True)
        # Include progress flag in the Build & Run tab so users can choose whether
        # the CLI includes --progress (UI shows progress regardless).
        progress_flag = st.checkbox("Include progress flag in command", value=False, disabled=not ordered_ready, help="Include --progress in the copyable CLI (UI always shows a progress bar during execution).")
        # Allow running only the SFS branch: enabled only when SFS is requested
        run_sfs_only_key = f"run_sfs_only_{model_key}"
        if run_sfs_only_key not in st.session_state:
            st.session_state[run_sfs_only_key] = False
        run_sfs_only = st.checkbox(
            "Run SFS only (omit primary output)",
            value=st.session_state.get(run_sfs_only_key, False),
            key=run_sfs_only_key,
            disabled=not ordered_ready or not sfs_on,
            help="When enabled the simulator will produce only the SFS artifact (no primary ms/vcf output). Requires Compute SFS to be checked."
        )
    # Button area for preparing or executing the engine command. Disabled
    # when required UI steps aren't complete. msprime is supported by the
    # simulator backend (it runs in-process), so do not disable Execute
    # when msprime is selected. The UI still hides Simulation-per-Work
    # controls for msprime above (they're not applicable).
    btn_disabled = not ordered_ready
    def build_cmd():
        if not ordered_ready:
            return None
        # Prefer the session-backed parallel value (slider uses a session key)
        # so that the built CLI always reflects the current widget state.
        try:
            parallel_effective = st.session_state.get(par_key, parallel)
        except Exception:
            parallel_effective = parallel
        # Helper: read from session_state if present, else fallback to module globals
        def _ss_glob_read(ss_key: str, glob_name: str, cast, default):
            try:
                if ss_key in st.session_state and st.session_state.get(ss_key) is not None:
                    return cast(st.session_state.get(ss_key))
            except Exception:
                pass
            try:
                v = globals().get(glob_name, None)
                if v is None:
                    return default
                return cast(v)
            except Exception:
                return default

        seq_length = _ss_glob_read(f"seq_length_{model_key}", 'seq_length', int, 0)
        target_snps = _ss_glob_read(f"target_snps_{model_key}", 'target_snps', int, 0)
        target_snps_tol = _ss_glob_read(f"target_snps_tol_{model_key}", 'target_snps_tol', float, 0.0)
        replicates = _ss_glob_read(f"replicates_{model_key}", 'replicates', int, 1)
        out_path = _ss_glob_read(f"output_path_{model_key}", 'out_path', str, None)
        out_format = _ss_glob_read(f"output_format_{model_key}", 'out_format', str, None)
        temp_loc = _ss_glob_read('temp_loc', 'temp_loc', str, 'local')
        sims_per_work = _ss_glob_read(f"sims_per_work_slider_{model_key}", 'sims_per_work', int, 0)
        max_ram_percent = _ss_glob_read(f"max_ram_percent_{model_key}", 'max_ram_percent', int, 80)
        sfs_output = _ss_glob_read(f"sfs_output_{model_key}", 'sfs_output', str, None)
        sfs_normalized = _ss_glob_read(f"sfs_normalized_{model_key}", 'sfs_normalized', bool, False)
        sfs_mode = _ss_glob_read(f"sfs_mode_{model_key}", 'sfs_mode', str, None)
        paired_neutral_name = _ss_glob_read(f"paired_neutral_name_{model_key}", 'paired_neutral_name', str, None)
        paired_neutral = _ss_glob_read(f"paired_neutral_{model_key}", 'paired_neutral', bool, False)
        neutral_engine = _ss_glob_read(f"neutral_engine_{model_key}", 'neutral_engine', str, '')
        no_en_ladder = _ss_glob_read('no_en_ladder', 'no_en_ladder', bool, False)
        growth_max_fold = _ss_glob_read('growth_max_fold', 'growth_max_fold', float, 1.05)
        sweep_enabled_ss_local = bool(st.session_state.get(f"sweep_enable_{model_key}", False))
        # sweep_pop: prefer selectbox value, then manual entry, then globals fallback
        try:
            sel_val = st.session_state.get(f"sweep_pop_{model_key}_sel", None)
            manual_val = st.session_state.get(f"sweep_pop_manual_{model_key}", None)
            if sel_val and sel_val != "-- select population --":
                sweep_pop = str(sel_val)
            elif manual_val:
                sweep_pop = str(manual_val)
            else:
                sweep_pop = globals().get('sweep_pop', '')
        except Exception:
            sweep_pop = globals().get('sweep_pop', '')
        # Determine sweep position percentage (0..100).
        # The UI stores the slider under key `sweep_pos_{model_key}`. When a
        # bp-based slider was shown this stores a base-pair integer; when a
        # percent slider was shown it stores a percent float. Recompute the
        # percentage here so the command builder always gets a percent.
        sweep_pos_pct = 50.0
        try:
            # First, detect if an explicit sequence length was provided via UI
            seq_len_ss = st.session_state.get(f"seq_length_{model_key}", None)
            seq_len_val = int(seq_len_ss) if seq_len_ss not in (None, "") else 0
        except Exception:
            seq_len_val = 0
        try:
            chrom_length_mode_ss = bool(st.session_state.get(f"chrom_length_mode_{model_key}", False))
        except Exception:
            chrom_length_mode_ss = False
        # Determine contig/native length fallback from globals (same logic as UI)
        length_for_pos = None
        try:
            if seq_len_val and seq_len_val > 0:
                length_for_pos = seq_len_val
            elif not seq_len_val and chrom_length_mode_ss:
                if 'contig_len' in globals() and contig_len:
                    length_for_pos = int(contig_len)
        except Exception:
            length_for_pos = None

        # Read raw slider value from session if present; otherwise fall back to global
        raw_val = None
        try:
            if f"sweep_pos_{model_key}" in st.session_state:
                raw_val = st.session_state.get(f"sweep_pos_{model_key}")
            else:
                raw_val = globals().get('sweep_pos_pct', None)
        except Exception:
            raw_val = globals().get('sweep_pos_pct', None)

        try:
            if length_for_pos and raw_val is not None:
                # raw_val is bp position -> convert to percent
                sweep_pos_pct = float(raw_val) / float(length_for_pos) * 100.0
            else:
                # raw_val is percent already (or None)
                sweep_pos_pct = float(raw_val if raw_val is not None else 50.0)
        except Exception:
            sweep_pos_pct = 50.0
        sel_s = _ss_glob_read(f"sel_s_{model_key}", 'sel_s', float, 0.0)
        time_units = _ss_glob_read(f"time_units_{model_key}", 'time_units', str, 'gens')
        sweep_time = _ss_glob_read(f"sweep_time_{model_key}", 'sweep_time', float, 0.0)
        fixation_time = _ss_glob_read(f"fixation_time_{model_key}", 'fixation_time', float, 0.0)
        cmd = [sys.executable, SIM_SCRIPT,
               "--engine", engine,
               "--pop-order", ",".join(pop_order),
               "--sample-individuals", ",".join(str(n) for n in counts)]
        # Only include species/model arguments when the user explicitly selected them
        if species_selected and species_id:
            cmd += ["--species-id", species_id]
        if model_selected and model_id:
            cmd += ["--model-id", model_id]
        if chromosome:
            cmd += ["--chromosome", str(chromosome)]
        # If chromosome/native length mode is enabled, ignore explicit length/target-snps
        chrom_length_mode_ss = bool(st.session_state.get(f"chrom_length_mode_{model_key}", False))
        if not chrom_length_mode_ss:
            if seq_length and seq_length > 0:
                cmd += ["--length", str(int(seq_length))]
            if target_snps and target_snps > 0:
                cmd += ["--target-snps", str(int(target_snps))]
                # Only include tolerance flag when user explicitly set a custom value (Auto unchecked)
                tol_auto_ss = bool(st.session_state.get(f"target_snps_tol_auto_{model_key}", True))
                if not tol_auto_ss:
                    # Read the session-backed tolerance value (percent). Accept
                    # fractional values <=1 for backward compatibility (convert
                    # to percent by *100). Only include flag when value > 0.
                    try:
                        raw_tol = float(st.session_state.get(f"target_snps_tol_{model_key}", 0.0) or 0.0)
                    except Exception:
                        raw_tol = float(target_snps_tol or 0.0)
                    # Normalize fractional inputs (<=1) into percent
                    tol_percent = raw_tol
                    try:
                        if 0.0 < raw_tol <= 1.0:
                            tol_percent = raw_tol * 100.0
                    except Exception:
                        pass
                    # Only include when strictly positive
                    try:
                        if float(tol_percent) > 0.0:
                            cmd += ["--target-snps-tol", str(float(tol_percent))]
                    except Exception:
                        pass
        if replicates and replicates != 1:
            cmd += ["--replicates", str(int(replicates))]
        # Only include primary output path when not running SFS-only
        if out_path and not st.session_state.get(f"run_sfs_only_{model_key}", False):
            # Ensure output path has the proper extension for chosen output format
            out_path_final = out_path
            try:
                fmt = str(out_format or "").strip()
            except Exception:
                fmt = ""
            # Determine canonical extensions for known formats
            fmt_ext_map = {
                "ms": ".ms",
                "ms.gz": ".ms.gz",
                "vcf": ".vcf",
                "vcf.gz": ".vcf.gz",
                "bcf": ".bcf",
            }
            ext = fmt_ext_map.get(fmt, "")
            if ext and not out_path_final.endswith(ext):
                # If user supplied an output without the chosen extension, append it
                out_path_final = out_path_final + ext
            cmd += ["--output", out_path_final]
        # Only include output-format when a primary output is being produced
        # (skip when user selected "Run SFS only").
        if out_format and not st.session_state.get(f"run_sfs_only_{model_key}", False):
            cmd += ["--output-format", out_format]
        # Only include --temp when the user explicitly chose a non-local temp
        # location. The simulator's default for local temp should be used when
        # temp_loc is 'local' or empty.
        try:
            if temp_loc and str(temp_loc).strip().lower() not in {"local", ""}:
                cmd += ["--temp", str(temp_loc)]
        except Exception:
            pass
        if progress_flag:
            cmd += ["--progress"]  # only included if user ticked it
        # Use the effective parallel value (session-backed when available).
        try:
            par_val = int(parallel_effective) if parallel_effective is not None else None
        except Exception:
            par_val = None
        if par_val and par_val != 1:
            cmd += ["--parallel", str(int(par_val))]
        # Only include sims-per-work when user selected a custom value (Auto unchecked).
        # Do NOT pass --sims-per-work to the simulator when engine is msprime
        # because msprime performs its own in-process batching.
        sims_auto_ss = bool(st.session_state.get(f"sims_per_work_auto_{model_key}", True))
        try:
            if engine == "msprime":
                pass
            else:
                if not sims_auto_ss and sims_per_work and sims_per_work > 0:
                    cmd += ["--sims-per-work", str(int(sims_per_work))]
        except Exception:
            # Defensive: if engine or session state access fails, fall back to previous behaviour
            if not sims_auto_ss and sims_per_work and sims_per_work > 0:
                cmd += ["--sims-per-work", str(int(sims_per_work))]
        # Always pass --max-ram-percent so the simulator receives an explicit
        # value even when the UI default (80) is selected.
        try:
            if max_ram_percent is None:
                # fallback to default slider value if missing
                mr = 80.0
            else:
                mr = float(max_ram_percent)
            cmd += ["--max-ram-percent", str(mr)]
            # also include the cap interpretation selected in the UI (system|run)
            try:
                max_ram_cap_val = _ss_glob_read(f"max_ram_cap_{model_key}", 'max_ram_cap', str, 'system')
                if max_ram_cap_val:
                    cmd += ["--max-ram-cap", str(max_ram_cap_val)]
            except Exception:
                # ignore failure and continue
                pass
        except Exception:
            # best effort: don't crash the UI if conversion fails
            try:
                cmd += ["--max-ram-percent", "80"]
            except Exception:
                pass
        if sfs_on:
            # If user provided an explicit SFS output, use it; otherwise derive
            # from the output filename by replacing its extension with .sfs when possible.
            sfs_arg = None
            try:
                if sfs_output and sfs_output.strip():
                    sfs_arg = sfs_output.strip()
                else:
                    # derive from out_path_final if available
                    if out_path:
                        # strip known extensions then append .sfs
                        base = out_path
                        # remove double extensions like .vcf.gz/.ms.gz first
                        for e in [".ms.gz", ".vcf.gz"]:
                            if base.endswith(e):
                                base = base[:-len(e)]
                                break
                        else:
                            # remove single extension
                            base = os.path.splitext(base)[0]
                        sfs_arg = base + ".sfs"
            except Exception:
                sfs_arg = None

            if sfs_arg:
                cmd += ["--sfs", sfs_arg]
            else:
                # leave simulator to use its default behavior
                cmd += ["--sfs"]
            if sfs_normalized:
                cmd += ["--sfs-normalized"]
            if sfs_mode:
                cmd += ["--sfs-mode", sfs_mode]
        if paired_neutral:
            # The paired_neutral_name stored in session_state is expected to be
            # a base name (without extensions). For robustness strip common
            # extensions if the user provided them accidentally.
            pn_val = (paired_neutral_name or "").strip()
            if pn_val:
                for e in [".ms.gz", ".vcf.gz", ".ms", ".vcf", ".bcf", ".sfs"]:
                    if pn_val.endswith(e):
                        pn_val = pn_val[: -len(e)]
                        break
                cmd += ["--paired-neutral", pn_val]
            else:
                cmd += ["--paired-neutral"]
            if neutral_engine.strip():
                cmd += ["--neutral-engine", neutral_engine.strip()]
        if no_en_ladder:
            cmd += ["--disable-growth-discretization"]
        if growth_max_fold and growth_max_fold != 1.05:
            cmd += ["--growth-max-fold", str(float(growth_max_fold))]
        # selection (only include when sweep parameters are enabled)
        sweep_enabled_ss_local = bool(st.session_state.get(f"sweep_enable_{model_key}", False))
        if sweep_enabled_ss_local:
            if sweep_pop and str(sweep_pop).strip():
                cmd += ["--sweep-pop", str(sweep_pop).strip()]
            if sel_s and float(sel_s) != 0.0:
                cmd += ["--sel-s", str(float(sel_s))]
            # sweep_pos_pct holds percentage (0..100); pass it directly to the
            # simulator as a percentage value (0..100) per updated engine API.
            # Ensure sweep_pos_pct is a float percentage (0..100). When sweep
            # parameters are enabled we always pass --sweep-pos (clamped to the
            # valid range) so downstream engines receive an explicit position.
            try:
                sp_pct = float(sweep_pos_pct)
            except Exception:
                sp_pct = 50.0
            # clamp into valid range
            try:
                if sp_pct < 0.0:
                    sp_pct = 0.0
                elif sp_pct > 100.0:
                    sp_pct = 100.0
            except Exception:
                sp_pct = max(0.0, min(100.0, float(sp_pct or 50.0)))
            # pass percentage (0..100) as the engine expects
            cmd += ["--sweep-pos", str(float(sp_pct))]
            if time_units:
                cmd += ["--time-units", time_units]
            # For msms engine, do not pass --sweep-time (msms expects selection
            # strength via --sel-s and handles timing differently). For other
            # engines include --sweep-time when provided.
            try:
                current_engine = engine
            except Exception:
                current_engine = None
            if sweep_time is not None and (current_engine is None or current_engine != "msms"):
                cmd += ["--sweep-time", str(float(sweep_time))]
            # For discoal, do not pass --fixation-time (discoal handles
            # fixation timing differently). For other engines, include it when set.
            if fixation_time is not None and (current_engine is None or current_engine != "discoal"):
                cmd += ["--fixation-time", str(float(fixation_time))]
        # Optional developer toggle to include --show-command
        # include_show_flag variable is set lower when building for Execute;
        # here we default to False when invoked directly.
        return cmd

    # helper to run a command and show progress + stdout/stderr
    def run_and_report(cmd_list):
        """Run the simulator subprocess and stream its stdout into the UI.
        Parse tqdm-like progress (percent or a/b) and update a Streamlit
        progress bar in real time. Capture full output and show it at the end.
        """
    # Do not display the full running command or raw terminal output in the UI.
    # We only present progress bars (percent/frac -> progress) and final
    # success/error state. The copyable CLI shown elsewhere remains available.
        import time, re
        try:
            ui_diag(f"[ui_simulator] run_and_report invoked for model_key={model_key} cmd={' '.join(cmd_list) if isinstance(cmd_list, (list,tuple)) else str(cmd_list)}")
        except Exception:
            pass

        # Start the subprocess, merging stderr into stdout so we only need one stream
        # Use a simple file lock in artifacts/ to avoid concurrent runs for the same model_key
        try:
            workspace_dir = os.path.abspath(os.getcwd())
            artifacts_dir = os.path.join(workspace_dir, 'artifacts')
            os.makedirs(artifacts_dir, exist_ok=True)
            lock_path = os.path.join(artifacts_dir, f".runlock_{model_key}")
            cancel_path = os.path.join(artifacts_dir, f".runlock_{model_key}.cancel")
        except Exception:
            lock_path = None
            cancel_path = None

        # track whether we actually created the lock so we can remove it on exit
        lock_created = False

        try:
            # if a lock exists, try to detect whether it's stale by reading PID
            if lock_path and os.path.exists(lock_path):
                try:
                    with open(lock_path, 'r') as fh:
                        content = fh.read().strip()
                    # expected format: <pid>\n<timestamp>
                    parts = [p for p in content.split() if p.strip()]
                    pid_existing = int(parts[0]) if parts else None
                except Exception:
                    pid_existing = None

                alive = False
                if pid_existing:
                    try:
                        # use psutil when available for a reliable check
                        if psutil:
                            alive = psutil.pid_exists(pid_existing)
                        else:
                            # os.kill(0) raises if not exists
                            os.kill(pid_existing, 0)
                            alive = True
                    except Exception:
                        alive = False

                if alive:
                    st.error("A run is already in progress for this model (another process is active).", icon="❌")
                    return
                else:
                    # stale lock: remove it and continue
                    try:
                        os.remove(lock_path)
                    except Exception:
                        pass

            # Start the subprocess in its own process group so we can signal it and its children
            try:
                # Prefer start_new_session for portability (creates separate process group)
                proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True)
            except TypeError:
                # Older Python/platforms might not accept start_new_session; fallback to preexec_fn
                try:
                    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, preexec_fn=os.setsid)
                except Exception:
                    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

            # create the lock containing PID and timestamp so other attempts see an active run
            if lock_path:
                try:
                    with open(lock_path, 'w') as fh:
                        fh.write(f"{proc.pid}\n{time.time()}")
                    lock_created = True
                except Exception:
                    lock_created = False
                # store running pid in session_state so the UI Cancel handler
                # can locate and signal the process group immediately
                try:
                    pid_key = f"_running_pid_{model_key}"
                    st.session_state[pid_key] = int(proc.pid)
                except Exception:
                    pass
                # artifact diagnostic: record run start (write to artifacts/ui_simulator_diag.log)
                try:
                    diag_path = os.path.join(artifacts_dir, 'ui_simulator_diag.log')
                    with open(diag_path, 'a') as dh:
                        dh.write(f"START {model_key} pid={proc.pid} time={time.time()}\n")
                except Exception:
                    pass
                    try:
                        pg = None
                        try:
                            pg = os.getpgid(proc.pid)
                        except Exception:
                            pg = None
                        # check psutil existence if available
                        alive = None
                        try:
                            if psutil:
                                alive = psutil.pid_exists(proc.pid)
                        except Exception:
                            alive = None
                        ui_diag(f"[ui_simulator] Started subprocess pid={proc.pid} pgid={pg} lock_created={lock_created} lock_path={lock_path} pid_alive={alive}")
                    except Exception:
                        pass
        except FileNotFoundError:
            # ensure we remove the lock if we created it
            try:
                if lock_created and lock_path and os.path.exists(lock_path):
                    os.remove(lock_path)
            except Exception:
                pass
            st.error("Could not find simulator.py in the working folder.", icon="❌")
            return
        except Exception as e:
            try:
                if lock_created and lock_path and os.path.exists(lock_path):
                    os.remove(lock_path)
            except Exception:
                pass
            st.error(f"Failed to start process: {e}", icon="❌")
            return
        captured_lines = []

        # Setup cancellation token in session_state so UI button can request stop
        cancel_key = f"_cancel_run_{model_key}"
        try:
            st.session_state[cancel_key] = False
        except Exception:
            # if session_state unavailable for some reason, continue without cancel
            pass

        # Primary progress bar (overall) and secondary for neutral if present.
        # Use labeled containers so the bars are clearly identified.
        status = st.empty()
        # Create distinct placeholders for labels, bars and ETAs so updating one
        # does not overwrite the others.
        # Helper to render a custom HTML progress bar with centered percentage text.
        def _render_progress_bar(ph, pct, bar_id: str = ""):
            try:
                # accept float or int; coerce to int percent for display
                pi = max(0, min(100, int(float(pct))))
            except Exception:
                pi = 0
            # text color depends on fill for legibility
            text_color = '#ffffff' if pi >= 50 else '#111827'
            # gradient fill for a nicer look
            bar_html = f"""
<div style="width:100%;background:#f1f5f9;border-radius:8px;overflow:hidden;border:1px solid #e2e8f0;padding:3px 3px;">
  <div style="position:relative;height:20px;border-radius:6px;overflow:hidden;background:#e6eef8;">
    <div style="width:{pi}%;height:100%;background:linear-gradient(90deg,#6366f1,#10b981);transition:width 0.3s ease;"></div>
    <div style="position:absolute;left:0;top:0;right:0;bottom:0;display:flex;align-items:center;justify-content:center;font-weight:600;color:{text_color};font-size:12px;">{pi}%</div>
  </div>
</div>
"""
            try:
                ph.markdown(bar_html, unsafe_allow_html=True)
            except Exception:
                # fallback: use built-in progress if HTML rendering fails
                try:
                    ph.progress(pi)
                except Exception:
                    pass

        sweep_label_ph = st.empty()
        # Use a concise, accurate label — it's not always a "sweep" phase.
        sweep_label_ph.markdown("**Core progress**")
        sweep_bar_ph = st.empty()
        _render_progress_bar(sweep_bar_ph, 0, bar_id=f"sweep_{model_key}")
        sweep_eta_ph = st.empty()

        # Do not create the paired-neutral bar until we see the first neutral
        # progress update. Creating it on-demand prevents showing an empty
        # disabled bar before the neutral phase actually starts.
        pn_label_ph = None
        pn_bar_ph = None
        pn_eta_ph = None

        percent_re = re.compile(r"(\d{1,3})\s*%")
        frac_re = re.compile(r"(\d+)\s*/\s*(\d+)")

        # Determine paired-neutral identifiers from the expected PN name
        pn_expected = st.session_state.get(f"_expected_pn_{model_key}") or ""
        def _strip_known_ext_local(name):
            if not name:
                return name
            b = os.path.basename(name)
            for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs']:
                if b.endswith(e):
                    return b[:-len(e)]
            return os.path.splitext(b)[0]

        pn_base = _strip_known_ext_local(pn_expected).lower() if pn_expected else ""
        # friendly local alias for mid-run detection (may be empty)
        try:
            pn_base_local = pn_base
        except Exception:
            pn_base_local = ""

        def _is_neutral_line(text: str) -> bool:
            if not text:
                return False
            low = text.lower()
            # explicit match on pn_base (preferred)
            if pn_base and pn_base in low:
                return True
            # fallback keywords: 'paired' or 'neutral' often appear in PN logs
            if 'paired' in low or 'neutral' in low:
                return True
            return False

        # Track separate ETA/progress for sweep and paired-neutral
        sweep_start = None
        sweep_last_progress = 0.0
        pn_start = None
        pn_last_progress = 0.0

        # Helper to compute ETA and store progress in session_state
        def _update_phase(pct_float, start_time, phase_prefix):
            # pct_float: 0..100
            try:
                now = time.time()
                if start_time is None:
                    start_time = now
                progress_frac = max(1e-9, pct_float / 100.0)
                elapsed = now - start_time
                est_total = elapsed / progress_frac
                remaining = max(0.0, est_total - elapsed)
                # store numeric values in session_state for external use
                try:
                    st.session_state[f"_progress_pct_{model_key}{phase_prefix}"] = float(pct_float)
                    st.session_state[f"_progress_eta_{model_key}{phase_prefix}"] = float(remaining)
                except Exception:
                    pass
                return start_time, int(remaining), int(est_total)
            except Exception:
                return start_time, None, None

        # Helper to format seconds into a human-friendly duration string.
        def _format_duration(sec):
            try:
                s = int(max(0, int(sec)))
            except Exception:
                return "0s"
            parts = []
            days, rem = divmod(s, 86400)
            if days:
                parts.append(f"{days}d")
            hours, rem = divmod(rem, 3600)
            if hours:
                parts.append(f"{hours}h")
            minutes, seconds = divmod(rem, 60)
            if minutes:
                parts.append(f"{minutes}m")
            if seconds or not parts:
                parts.append(f"{seconds}s")
            return " ".join(parts)

        # Stream the process output line-by-line and parse progress updates.
        # This loop checks session_state cancel flag and attempts to terminate
        # the subprocess gracefully if requested.
        try:
            if proc.stdout is None:
                # No stream available; wait for completion and capture output
                out_text = proc.communicate()[0] or ""
                captured_lines.append(out_text)
            else:
                while True:
                    # check cancel request (session flag or cancel sentinel file)
                    cancel_flag = False
                    try:
                        if st.session_state.get(cancel_key, False):
                            cancel_flag = True
                    except Exception:
                        cancel_flag = False
                    try:
                        if not cancel_flag and cancel_path and os.path.exists(cancel_path):
                            cancel_flag = True
                    except Exception:
                        pass
                    if cancel_flag:
                        # Robust cancellation: try graceful termination of the
                        # whole process group, fallback to terminating children
                        # (via psutil) and finally force-kill the group.
                        try:
                            status.text("Cancelling…")
                        except Exception:
                            pass

                        try:
                            pgid = None
                            # Prefer signaling the process group so children are included
                            try:
                                pgid = os.getpgid(proc.pid)
                                ui_diag(f"[ui_simulator] Cancelling run pid={proc.pid} pgid={pgid}")
                                try:
                                    os.killpg(pgid, signal.SIGTERM)
                                except Exception:
                                    # fallback to terminating the process directly
                                    proc.terminate()
                            except Exception:
                                # getpgid may fail; fallback
                                try:
                                    proc.terminate()
                                except Exception:
                                    pass

                            # Allow a short grace period for processes to exit
                            try:
                                proc.wait(timeout=5)
                                ui_diag(f"[ui_simulator] Process {proc.pid} terminated after SIGTERM.")
                            except Exception:
                                # If psutil is available, try terminating children first
                                try:
                                    if psutil:
                                        try:
                                            p = psutil.Process(proc.pid)
                                            children = p.children(recursive=True)
                                            for c in children:
                                                try:
                                                    ui_diag(f"[ui_simulator] terminating child pid={c.pid}")
                                                    c.terminate()
                                                except Exception:
                                                    pass
                                            gone, alive = psutil.wait_procs(children, timeout=3)
                                            if alive:
                                                for c in alive:
                                                    try:
                                                        ui_diag(f"[ui_simulator] killing child pid={c.pid}")
                                                        c.kill()
                                                    except Exception:
                                                        pass
                                        except Exception:
                                            pass
                                except Exception:
                                    pass

                                # Finally, force kill the entire process group
                                try:
                                    if pgid is not None:
                                        os.killpg(pgid, signal.SIGKILL)
                                    else:
                                        proc.kill()
                                    ui_diag(f"[ui_simulator] Sent SIGKILL to pid/group {proc.pid}")
                                except Exception:
                                    try:
                                        proc.kill()
                                    except Exception:
                                        pass

                        except Exception as e:
                            ui_diag(f"[ui_simulator] Error during cancel: {e}")

                        # Ensure lock cleaned up immediately so future runs are allowed
                        try:
                            if lock_created and lock_path and os.path.exists(lock_path):
                                try:
                                    os.remove(lock_path)
                                    ui_diag(f"[ui_simulator] Removed run lock {lock_path} after cancel")
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        # remove cancel sentinel if present
                        try:
                            if cancel_path and os.path.exists(cancel_path):
                                try:
                                    os.remove(cancel_path)
                                except Exception:
                                    pass
                        except Exception:
                            pass

                        # artifact diagnostic: record cancel completion
                        try:
                            diag_path = os.path.join(artifacts_dir, 'ui_simulator_diag.log')
                            with open(diag_path, 'a') as dh:
                                dh.write(f"CANCELLED {model_key} pid={proc.pid} time={time.time()}\n")
                        except Exception:
                            pass

                        break

                    raw = proc.stdout.readline()
                    if raw == '' and proc.poll() is not None:
                        break
                    if not raw:
                        time.sleep(0.01)
                        continue
                    captured_lines.append(raw)
                    # tqdm often uses carriage returns to update the same line; split and take last token
                    last = raw.rstrip('\n').split('\r')[-1].strip()
                    if not last:
                        continue

                    # percent like ' 42%'
                    m = percent_re.search(last)
                    pct = None
                    if m:
                        try:
                            pct = float(m.group(1))
                            pct = max(0.0, min(100.0, pct))
                        except Exception:
                            pct = None

                    # fraction like '3/10'
                    if pct is None:
                        m2 = frac_re.search(last)
                        if m2:
                            try:
                                a = int(m2.group(1))
                                b = int(m2.group(2))
                                if b > 0:
                                    pct = (a / float(b)) * 100.0
                                    pct = max(0.0, min(100.0, pct))
                                else:
                                    pct = None
                            except Exception:
                                pct = None

                    if pct is not None:
                        try:
                            # Decide whether this update is for the neutral phase
                            if _is_neutral_line(last):
                                # create the paired-neutral bar on first sighting
                                if pn_bar_ph is None:
                                    try:
                                        pn_label_ph = st.empty()
                                        pn_label_ph.markdown("**Paired neutral progress**")
                                        pn_bar_ph = st.empty()
                                        _render_progress_bar(pn_bar_ph, 0, bar_id=f"pn_{model_key}")
                                        pn_eta_ph = st.empty()
                                    except Exception:
                                        pn_bar_ph = None
                                if pn_bar_ph is not None:
                                    # Only update the numeric progress bar and ETA;
                                    # do not echo the raw terminal text into the UI.
                                    _render_progress_bar(pn_bar_ph, pct, bar_id=f"pn_{model_key}")
                                    if pn_start is None:
                                        pn_start = time.time()
                                    pn_start, remaining, est_total = _update_phase(pct, pn_start, "_pn")
                                    try:
                                        # Explicit remaining-time label per request
                                            if remaining is not None and pn_eta_ph is not None:
                                                pn_eta_ph.text(f"Remaining Time: {_format_duration(remaining)}")
                                    except Exception:
                                        pass
                            else:
                                # core progress
                                # Only update numeric core progress and ETA; do not echo raw text.
                                _render_progress_bar(sweep_bar_ph, pct, bar_id=f"sweep_{model_key}")
                                if sweep_start is None:
                                    sweep_start = time.time()
                                sweep_start, remaining, est_total = _update_phase(pct, sweep_start, "")
                                try:
                                    # Explicit remaining-time label per request
                                    if remaining is not None:
                                        sweep_eta_ph.text(f"Remaining Time: {_format_duration(remaining)}")
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        continue

                    # Mid-run artifact discovery: if a sweep file (heuristic) appears, register it so user can download early
                    try:
                        mpath = re.search(r"\b(\S+\.(?:vcf|vcf\.gz|ms|ms\.gz|sfs|bcf))\b", last)
                        if mpath:
                            candidate = mpath.group(1)
                            # If path exists on disk and looks like a sweep artifact, register it
                            if os.path.exists(candidate):
                                reg_key_mid = f"generated_files_{model_key}"
                                cur = st.session_state.get(reg_key_mid, [])
                                absp = os.path.abspath(candidate)
                                already = any(os.path.abspath(x[1]) == absp for x in cur)
                                if not already:
                                    # compute label for candidate
                                    try:
                                        # normalize comparison by stripping known extensions
                                        def _base_noext(p):
                                            b = os.path.basename(p)
                                            for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs']:
                                                if b.endswith(e):
                                                    return b[:-len(e)]
                                            return os.path.splitext(b)[0]

                                        cand_base = _base_noext(candidate)
                                        out_base = _base_noext(out_path or "") if out_path else ""
                                        pn_base_local = _base_noext((paired_neutral_name or "").strip()) if paired_neutral_name else ""

                                        # First, check if this is a paired-neutral artifact
                                        # (match by base without extension). If so, label
                                        # accordingly and skip the generic SFS check.
                                        if pn_base_local and cand_base == pn_base_local:
                                            if candidate.lower().endswith('.sfs'):
                                                lab = "Paired Neutral SFS Output"
                                            else:
                                                lab = "Paired Neutral Output"
                                        # Next, check for main SFS by extension or exact
                                        # filename match to the configured sfs_output.
                                        elif candidate.lower().endswith('.sfs') or (sfs_output and os.path.basename(candidate) == os.path.basename(sfs_output)):
                                            lab = "SFS Output"
                                        elif out_base and cand_base == out_base:
                                            lab = "Variant Output"
                                        else:
                                            lab = os.path.basename(candidate) or candidate
                                    except Exception:
                                        lab = os.path.basename(candidate) or candidate

                                    # schedule zip only if directory; otherwise add directly
                                    if os.path.isdir(candidate):
                                        try:
                                            # derive a clean base name (strip known extensions)
                                            base = os.path.basename(candidate.rstrip(os.path.sep))
                                            for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs']:
                                                if base.endswith(e):
                                                    base = base[:-len(e)]
                                                    break
                                            # If this candidate is classified as Paired Neutral,
                                            # prefer the paired-neutral base name for the archive
                                            zip_base = pn_base_local or base
                                            # create the zip inside a workspace 'artifacts' folder
                                            workspace_dir = os.path.abspath(os.getcwd())
                                            artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                                            try:
                                                os.makedirs(artifacts_dir, exist_ok=True)
                                            except Exception:
                                                pass
                                            archive_name = f"{zip_base}.zip"
                                            archive_base = os.path.join(artifacts_dir, zip_base)
                                            zip_path = os.path.join(artifacts_dir, archive_name)
                                            # avoid overwriting an existing archive by appending a timestamp if needed
                                            if os.path.exists(zip_path):
                                                archive_base = os.path.join(artifacts_dir, f"{zip_base}_{int(time.time())}")
                                            else:
                                                archive_base = os.path.join(artifacts_dir, zip_base)
                                            zip_path = shutil.make_archive(archive_base, 'zip', root_dir=candidate)
                                            # after creating the archive, remove the original folder
                                            try:
                                                shutil.rmtree(candidate)
                                            except Exception:
                                                pass
                                            cur.append((lab, candidate, zip_path, True))
                                        except Exception:
                                            pass
                                    else:
                                        cur.append((lab, candidate, candidate, False))
                                    st.session_state[reg_key_mid] = cur
                    except Exception:
                        pass
                    # Fallback: do not display raw terminal lines in the UI. We keep
                    # progress bars and ETA updates only.
                    pass
        except Exception:
            # If streaming fails, fall back to waiting for the process
            try:
                proc.wait()
            except Exception:
                pass

    # Finalize
        try:
            proc.poll()
        except Exception:
            pass
        try:
            # final render at 100%
            _render_progress_bar(sweep_bar_ph, 100, bar_id=f"sweep_{model_key}")
        except Exception:
            pass
        # Minimal final text: indicate finalizing (kept brief) and then show
        # a success/error banner below. Do not display the full process stdout.
        try:
            sweep_eta_ph.text("Finalizing…")
        except Exception:
            pass
        try:
            # Do not display a success popup when the run finishes; keep
            # the final progress bar and finalizing text visible for the
            # user. Only show an error banner when the process exits
            # with a non-zero return code.
            if proc.returncode != 0:
                st.error(f"Process exited with code {proc.returncode}.", icon="❌")
        except Exception:
            pass

        # Ensure UI running/cancel flags are cleared so the Execute button
        # reappears after the run (or cancellation) completes.
        try:
            running_key = f"_running_{model_key}"
            cancel_key = f"_cancel_run_{model_key}"
            st.session_state[running_key] = False
            # reset cancel flag too so subsequent runs start clean
            st.session_state[cancel_key] = False
        except Exception:
            pass
        # clear start_request/processed so Execute becomes available again
        try:
            start_request_key = f"_start_request_{model_key}"
            processed_key = f"_exec_processed_{model_key}"
            st.session_state[start_request_key] = False
            st.session_state[processed_key] = False
        except Exception:
            pass
        try:
            immediate_key = f"_exec_start_immediate_{model_key}"
            st.session_state[immediate_key] = False
        except Exception:
            pass
        try:
            ui_lock_key = f"_ui_lock_{model_key}"
            st.session_state[ui_lock_key] = False
        except Exception:
            pass
        # signal to the UI render that a run finished so it can perform a
        # safe rerun from the top-level script (some Streamlit contexts do
        # not allow experimental_rerun() directly from background threads).
        try:
            finished_key = f"_exec_finished_{model_key}"
            st.session_state[finished_key] = True
        except Exception:
            pass
        # Opportunistically clean up any stale lock files for ANY model key so
        # the UI never remains disabled due to crashes or aborted runs.
        try:
            workspace_dir = os.path.abspath(os.getcwd())
            artifacts_dir = os.path.join(workspace_dir, 'artifacts')
            for name in os.listdir(artifacts_dir):
                if not name.startswith('.runlock_'):  # only inspect runlock files
                    continue
                path = os.path.join(artifacts_dir, name)
                try:
                    with open(path, 'r') as fh:
                        content = fh.read().strip()
                    parts = [p for p in content.split() if p.strip()]
                    pid_existing = int(parts[0]) if parts else None
                except Exception:
                    pid_existing = None
                alive = False
                if pid_existing:
                    try:
                        if psutil:
                            alive = psutil.pid_exists(pid_existing)
                        else:
                            os.kill(pid_existing, 0)
                            alive = True
                    except Exception:
                        alive = False
                if not alive:
                    try:
                        os.remove(path)
                    except Exception:
                        pass
        except Exception:
            pass
        # ---------------------
        # Discover produced artifacts and register them for download/cleanup
        # ---------------------
        try:
            gen_list = []
            # Use the expected output names (stashed before execution) as the
            # authoritative list of artifacts we expect to produce. This avoids
            # brittle heuristics and ensures the UI always lists SFS/primary/PN
            # artifacts immediately after Execute (even if the files are not
            # yet present on disk). We will still scan the filesystem for any
            # additional artifacts produced.
            exp_out = st.session_state.get(f"_expected_out_{model_key}")
            exp_sfs = st.session_state.get(f"_expected_sfs_{model_key}")
            exp_pn = st.session_state.get(f"_expected_pn_{model_key}")

            # primary output (may be None when running SFS-only)
            if exp_out:
                # store the expected final path (may not exist yet)
                # When output format is VCF/BCF and replicates >= 2 the
                # simulator produces a directory of per-replicate files
                # named after the output base (without extension). In that
                # case, prefer the directory base so the UI will zip it for
                # download.
                try:
                    fmt = str(st.session_state.get(f"output_format_{model_key}", "") or out_format or "").strip()
                except Exception:
                    fmt = str(out_format) if 'out_format' in locals() else ""
                try:
                    reps = int(replicates or 1)
                except Exception:
                    reps = 1

                if fmt in {"vcf", "vcf.gz", "bcf"} and reps >= 2:
                    base = exp_out
                    for e in [".ms.gz", ".vcf.gz"]:
                        if base.endswith(e):
                            base = base[:-len(e)]
                            break
                    else:
                        base = os.path.splitext(base)[0]
                    gen_list.append(base)
                else:
                    gen_list.append(exp_out)

            # sfs output: only add expected SFS candidates when Compute SFS is enabled
            if sfs_on:
                if exp_sfs:
                    sfs_expected = exp_sfs
                    if not sfs_expected.lower().endswith('.sfs'):
                        sfs_expected = sfs_expected + '.sfs'
                    gen_list.append(sfs_expected)
                else:
                    # derive from exp_out when available
                    if exp_out:
                        base = exp_out
                        for e in [".ms.gz", ".vcf.gz"]:
                            if base.endswith(e):
                                base = base[:-len(e)]
                                break
                        else:
                            base = os.path.splitext(base)[0]
                        gen_list.append(base + '.sfs')

            # paired neutral: user-provided base; include only the most
            # likely artifact to avoid noisy duplicate entries in the UI.
            if exp_pn:
                pn_val = exp_pn
                try:
                    fmt = str(st.session_state.get(f"output_format_{model_key}", "") or out_format or "").strip()
                except Exception:
                    fmt = str(out_format) if 'out_format' in locals() else ""
                try:
                    reps = int(replicates or 1)
                except Exception:
                    reps = 1

                # If main output is a multi-replicate directory, prefer the directory base
                if fmt in {"vcf", "vcf.gz", "bcf"} and reps >= 2:
                    base = pn_val
                    for e in [".ms.gz", ".vcf.gz"]:
                        if base.endswith(e):
                            base = base[:-len(e)]
                            break
                    else:
                        base = os.path.splitext(base)[0]
                    gen_list.append(base)
                else:
                    # prefer the explicit PN value (may include extension)
                    gen_list.append(pn_val)

            # Also include any files discovered on disk that match common artifact
            # extensions (fallback/extra artifacts). This captures cases where
            # the simulator created files we didn't predict. Prefer files under
            # ./artifacts/ when present.
            try:
                # scan captured stdout for explicit paths first (existing logic)
                for ln in captured_lines:
                    mpath = re.search(r"\b(\S+\.(?:vcf|vcf\.gz|ms|ms\.gz|sfs|bcf))\b", ln)
                    if mpath:
                        candidate = mpath.group(1)
                        try:
                            # prefer artifact under ./artifacts/
                            workspace_dir = os.path.abspath(os.getcwd())
                            artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                            art_candidate = os.path.join(artifacts_dir, os.path.basename(candidate))
                        except Exception:
                            art_candidate = None
                        if art_candidate and os.path.exists(art_candidate) and art_candidate not in gen_list:
                            gen_list.append(art_candidate)
                        elif os.path.exists(candidate) and candidate not in gen_list:
                            gen_list.append(candidate)
            except Exception:
                pass

            # Build canonical download items from expected outputs first so
            # labels are deterministic and not subject to heuristics.
            # We store tuples as (label, original_path, download_path, zipped_flag)
            download_items = []

            def _maybe_add(path, label):
                if not path:
                    return
                # avoid duplicates (by absolute path)
                try:
                    ab = os.path.abspath(path)
                except Exception:
                    ab = path
                for _, _, dpath, _ in download_items:
                    try:
                        if os.path.abspath(dpath) == ab:
                            return
                    except Exception:
                        pass
                # if it's a directory, zip it into temp and mark zipped
                if os.path.isdir(path):
                    try:
                        base = os.path.basename(path.rstrip(os.sep))
                        for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs']:
                            if base.endswith(e):
                                base = base[:-len(e)]
                                break
                        # create zips inside workspace 'artifacts' folder
                        workspace_dir = os.path.abspath(os.getcwd())
                        artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                        try:
                            os.makedirs(artifacts_dir, exist_ok=True)
                        except Exception:
                            pass
                        archive_name = f"{base}.zip"
                        archive_base = os.path.join(artifacts_dir, base)
                        zip_path = os.path.join(artifacts_dir, archive_name)
                        # avoid overwriting by appending timestamp if needed
                        if os.path.exists(zip_path):
                            archive_base = os.path.join(artifacts_dir, f"{base}_{int(time.time())}")
                        else:
                            archive_base = os.path.join(artifacts_dir, base)
                        zip_path = shutil.make_archive(archive_base, 'zip', root_dir=path)
                        # remove the original folder after zipping so artifacts/ contains only archives
                        try:
                            shutil.rmtree(path)
                        except Exception:
                            pass
                        download_items.append((label, path, zip_path, True))
                        return
                    except Exception:
                        pass
                # otherwise add the file path (may not exist yet)
                download_items.append((label, path, path, False))

            # Variant output (only when not running SFS-only)
            run_sfs_only_ss = bool(st.session_state.get(f"run_sfs_only_{model_key}", False))
            if exp_out and not run_sfs_only_ss:
                # If main output is a multi-replicate directory (vcf/bcf with reps>=2)
                # use the directory base so _maybe_add will zip it for download.
                try:
                    fmt = str(st.session_state.get(f"output_format_{model_key}", "") or out_format or "").strip()
                except Exception:
                    fmt = str(out_format) if 'out_format' in locals() else ""
                try:
                    reps = int(replicates or 1)
                except Exception:
                    reps = 1
                if fmt in {"vcf", "vcf.gz", "bcf"} and reps >= 2:
                    base = exp_out
                    for e in [".ms.gz", ".vcf.gz"]:
                        if base.endswith(e):
                            base = base[:-len(e)]
                            break
                    else:
                        base = os.path.splitext(base)[0]
                    variant_add_path = base
                else:
                    variant_add_path = exp_out
                _maybe_add(variant_add_path, "Variant Output")

            # SFS output
            if sfs_on:
                try:
                    workspace_dir = os.path.abspath(os.getcwd())
                    artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                except Exception:
                    artifacts_dir = None

                if exp_sfs:
                    sfs_path = exp_sfs
                    if not sfs_path.lower().endswith('.sfs'):
                        sfs_path = sfs_path + '.sfs'
                    # prefer artifacts/<basename>.sfs when present
                    try:
                        if artifacts_dir:
                            art_sfs = os.path.join(artifacts_dir, os.path.basename(sfs_path))
                            if os.path.exists(art_sfs):
                                _maybe_add(art_sfs, "SFS Output")
                            else:
                                _maybe_add(sfs_path, "SFS Output")
                        else:
                            _maybe_add(sfs_path, "SFS Output")
                    except Exception:
                        _maybe_add(sfs_path, "SFS Output")
                else:
                    # derive from exp_out if present
                    if exp_out:
                        base = exp_out
                        for e in [".ms.gz", ".vcf.gz"]:
                            if base.endswith(e):
                                base = base[:-len(e)]
                                break
                        else:
                            base = os.path.splitext(base)[0]
                        sfs_candidate = base + '.sfs'
                        try:
                            if artifacts_dir:
                                art_sfs = os.path.join(artifacts_dir, os.path.basename(sfs_candidate))
                                if os.path.exists(art_sfs):
                                    _maybe_add(art_sfs, "SFS Output")
                                else:
                                    _maybe_add(sfs_candidate, "SFS Output")
                            else:
                                _maybe_add(sfs_candidate, "SFS Output")
                        except Exception:
                            _maybe_add(sfs_candidate, "SFS Output")

            # Paired neutral expected artifacts
            if exp_pn:
                # the user-provided PN base may be a bare base or a path; add
                # the base itself and common engine extensions we might produce
                pn_dir = os.path.dirname(exp_pn) or '.'
                # strip any extension the user may have included so we always
                # build candidates from the plain base name
                def _strip_known_ext(name):
                    b = name
                    for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs']:
                        if b.endswith(e):
                            return b[:-len(e)]
                    return os.path.splitext(b)[0]
                pn_base = _strip_known_ext(os.path.basename(exp_pn))
                # prefer the same extension as the main variant output when available
                variant_ext = None
                try:
                    if exp_out:
                        # detect double extensions first
                        for ex in ['.ms.gz', '.vcf.gz']:
                            if exp_out.endswith(ex):
                                variant_ext = ex
                                break
                        if variant_ext is None:
                            variant_ext = os.path.splitext(exp_out)[1]
                            if variant_ext == '':
                                variant_ext = None
                    # If no explicit expected out path, fall back to the selected
                    # output format (session-backed). This produces deterministic
                    # PN filenames like a_neutral.ms when user selected output-format ms.
                    if variant_ext is None:
                        try:
                            fmt = str(st.session_state.get(f"output_format_{model_key}", "") or "").strip()
                            fmt_map = {"ms": ".ms", "ms.gz": ".ms.gz", "vcf": ".vcf", "vcf.gz": ".vcf.gz", "bcf": ".bcf"}
                            variant_ext = fmt_map.get(fmt)
                        except Exception:
                            variant_ext = None
                except Exception:
                    variant_ext = None

                # Prefer adding PN with an explicit extension. Use the PN base
                # (stripped of any extension the user might have typed) and append
                # the chosen variant extension so we produce a file like
                # '<pn_base>.ms' or '<pn_base>.vcf.gz'. Avoid adding the bare base
                # without extension unless no candidates exist.
                added_pn = False
                run_sfs_only_ss = bool(st.session_state.get(f"run_sfs_only_{model_key}", False))
                # Do not add paired-neutral variant outputs when running SFS-only
                if not run_sfs_only_ss:
                    if variant_ext:
                        _maybe_add(os.path.join(pn_dir, pn_base + variant_ext), "Paired Neutral Variant Output")
                        added_pn = True
                    else:
                        # fallback: consider common extensions and add the first that exists or all candidates
                        for e in ['.vcf', '.vcf.gz', '.ms', '.ms.gz', '.bcf']:
                            _maybe_add(os.path.join(pn_dir, pn_base + e), "Paired Neutral Variant Output")
                            added_pn = True
                else:
                    # when run_sfs_only is True we still want PN SFS if requested
                    added_pn = False

                # If none of the extension candidates were usable (rare), then
                # finally add the bare PN base so the UI still shows an entry.
                if not added_pn and not run_sfs_only_ss:
                    _maybe_add(os.path.join(pn_dir, pn_base), "Paired Neutral Variant Output")

                # also add paired neutral SFS candidate only when SFS computation is enabled
                if sfs_on:
                    try:
                        workspace_dir = os.path.abspath(os.getcwd())
                        artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                    except Exception:
                        artifacts_dir = None
                    pn_sfs = os.path.join(pn_dir, pn_base + '.sfs')
                    try:
                        if artifacts_dir:
                            # prefer artifacts/<pn_base>.sfs explicitly
                            art_pn_sfs = os.path.join(artifacts_dir, f"{pn_base}.sfs")
                            if os.path.exists(art_pn_sfs):
                                _maybe_add(art_pn_sfs, "Paired Neutral SFS Output")
                            else:
                                art_alt = os.path.join(artifacts_dir, os.path.basename(pn_sfs))
                                if os.path.exists(art_alt):
                                    _maybe_add(art_alt, "Paired Neutral SFS Output")
                                else:
                                    _maybe_add(pn_sfs, "Paired Neutral SFS Output")
                        else:
                            _maybe_add(pn_sfs, "Paired Neutral SFS Output")
                    except Exception:
                        _maybe_add(pn_sfs, "Paired Neutral SFS Output")

            # Append any additional discovered files from captured_lines that are
            # present on disk and not already included.
            for ln in captured_lines:
                mpath = re.search(r"\b(\S+\.(?:vcf|vcf\.gz|ms|ms\.gz|sfs|bcf))\b", ln)
                if mpath:
                    candidate = mpath.group(1)
                    try:
                        if os.path.exists(candidate):
                            # attempt to classify: skip if already included
                            _maybe_add(candidate, os.path.basename(candidate) or candidate)
                    except Exception:
                        pass

            # Finally, include any files that happen to exist in gen_list but
            # were not part of the expected names (fallback/extra artifacts).
            # Prefer artifacts/<basename> when present.
            try:
                workspace_dir = os.path.abspath(os.getcwd())
                artifacts_dir = os.path.join(workspace_dir, 'artifacts')
            except Exception:
                artifacts_dir = None

            for p in gen_list:
                try:
                    if not p:
                        continue
                    chosen = p
                    if artifacts_dir:
                        art_alt = os.path.join(artifacts_dir, os.path.basename(p))
                        if os.path.exists(art_alt):
                            chosen = art_alt
                    if chosen and os.path.exists(chosen):
                        _maybe_add(chosen, os.path.basename(chosen) or chosen)
                    else:
                        # include anyway (may be not-yet-existing path)
                        _maybe_add(p, os.path.basename(p) or p)
                except Exception:
                    try:
                        _maybe_add(p, os.path.basename(p) or p)
                    except Exception:
                        pass

            # Filter out noisy duplicates: if we have expected exp_out or exp_pn
            # prefer the labeled entries and drop other paths that share the
            # same canonical base (e.g., 'test', 'test_neutral.ms').
            def _canon_base(path):
                try:
                    b = os.path.basename(path or "")
                except Exception:
                    b = str(path or "")
                for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs', '.zip']:
                    if b.endswith(e):
                        return b[:-len(e)]
                return os.path.splitext(b)[0]

            try:
                exp_out_canon = _canon_base(exp_out) if exp_out else None
            except Exception:
                exp_out_canon = None
            try:
                exp_pn_canon = _canon_base(exp_pn) if exp_pn else None
            except Exception:
                exp_pn_canon = None

            filtered_items = []
            for item in download_items:
                lab, orig, dpath, zflag = item
                try:
                    b = _canon_base(orig)
                except Exception:
                    b = None
                if exp_pn_canon and b == exp_pn_canon and lab not in ("Paired Neutral Variant Output", "Paired Neutral SFS Output"):
                    # skip noisy non-PN labeled entries
                    continue
                # Keep SFS Output entries even when their canonical base matches
                # the primary base. Only skip duplicates that are neither the
                # primary Variant Output nor SFS Output.
                if exp_out_canon and b == exp_out_canon and lab not in ("Variant Output", "SFS Output"):
                    # skip noisy non-primary labeled entries that match the primary base
                    continue
                filtered_items.append(item)
            download_items = filtered_items

            # Build canonical ordering: keep first occurrence per canonical label
            order = ["Variant Output", "SFS Output", "Paired Neutral Variant Output", "Paired Neutral SFS Output"]
            by_label = {}
            others = []
            for item in download_items:
                lab = item[0] or ""
                if lab in order:
                    if lab not in by_label:
                        by_label[lab] = item
                else:
                    others.append(item)

            ordered = [by_label[l] for l in order if l in by_label]
            ordered.extend(others)

            # store in session for UI and attempt to cleanup at process exit
            reg_key = f"generated_files_{model_key}"
            st.session_state[reg_key] = ordered
            # global registry for atexit cleanup
            try:
                if '_generated_files_registry' not in globals():
                    # registry maps path -> expiry_timestamp (float)
                    globals()['_generated_files_registry'] = {}
                # set expiry only for files created in the system tempdir
                reg = globals()['_generated_files_registry']
                for _, dpath, _, _ in download_items:
                    try:
                        dpath_abs = os.path.abspath(dpath)
                        # default: do not set expiry for user files (in repo cwd)
                        expiry = None
                        try:
                            tmpdir = os.path.abspath(tempfile.gettempdir())
                            if dpath_abs.startswith(tmpdir + os.sep) or dpath_abs == tmpdir:
                                expiry = time.time() + CLEANUP_TTL
                        except Exception:
                            expiry = time.time() + CLEANUP_TTL
                        reg[dpath_abs] = expiry
                    except Exception:
                        pass
                globals()['_generated_files_registry'] = reg
            except Exception:
                pass

            def _cleanup_at_exit():
                try:
                    for p in list(globals().get('_generated_files_registry', [])):
                        try:
                            if os.path.isfile(p):
                                os.remove(p)
                            elif os.path.isdir(p):
                                shutil.rmtree(p)
                        except Exception:
                            pass
                except Exception:
                    pass

            try:
                atexit.register(_cleanup_at_exit)
            except Exception:
                pass
        except Exception:
            pass

        # remove lock file if we created it so future runs are allowed
        try:
            if lock_created and lock_path and os.path.exists(lock_path):
                os.remove(lock_path)
        except Exception:
            pass
        # Robust rerun attempt so the Execute button state refreshes immediately.
        # Prefer stable st.rerun when available; fallback to experimental_rerun.
        reran = False
        try:
            if hasattr(st, 'rerun'):
                st.rerun()
                reran = True
        except Exception:
            reran = False
        if not reran:
            try:
                if hasattr(st, 'experimental_rerun'):
                    st.experimental_rerun()
                    reran = True
            except Exception:
                reran = False
        # If we still couldn't rerun (older Streamlit), set a refresh flag the
        # top-level UI will notice on the next natural interaction.
        if not reran:
            try:
                st.session_state[f"_needs_refresh_{model_key}"] = True
            except Exception:
                pass

    with tab_map["Build & Run"]:
        # Button that prepares (but does not execute) the engine command and
        # makes it available as a copyable string.
        # Validate SFS output requirement: if Compute SFS is enabled and no primary
        # Output path was provided, require a specific SFS output path. When the
        # requirement is not met, disable the Prepare/Execute buttons.
        sfs_missing_req = False
        try:
            out_val = (st.session_state.get(f"output_path_{model_key}", "") or "").strip()
            sfs_val = (st.session_state.get(f"sfs_output_{model_key}", "") or "").strip()
        except Exception:
            out_val = ""
            sfs_val = ""
        # Determine sfs_on and run_sfs_only from session_state (stable keys)
        sfs_on_local = bool(st.session_state.get(f"sfs_enable_{model_key}", False))
        run_sfs_only_local = bool(st.session_state.get(f"run_sfs_only_{model_key}", False))

        # Validation rules:
        # - If Compute SFS is enabled and no primary output is provided, require an SFS output path.
        # - If NOT running SFS-only, require a primary output path.
        # - If running SFS-only but SFS output is not provided, require a primary output path.
        sfs_missing_req = False
        primary_missing_req = False

        if sfs_on_local and not out_val and not sfs_val:
            sfs_missing_req = True
            st.error("Compute SFS is enabled: provide an SFS output path when no primary Output path is set.")

        # require primary output unless (run_sfs_only_local is True AND sfs_val is set)
        if not (run_sfs_only_local and sfs_val):
            # primary output is required
            if not out_val:
                primary_missing_req = True
                st.error("Primary Output path is required unless 'Run SFS only' is checked and an SFS output path is provided.")

        # If engine is discoal and sweep parameters are enabled, ensure sweep_time > 0
        try:
            sweep_enabled_ss_local = bool(st.session_state.get(f"sweep_enable_{model_key}", False))
            if engine == 'discoal' and sweep_enabled_ss_local:
                sweep_time_val = float(st.session_state.get(f"sweep_time_{model_key}", 0.0) or 0.0)
                if sweep_time_val <= 0.0:
                    st.error("discoal requires a positive Sweep time. Set 'Sweep time' > 0.")
                    primary_missing_req = True
        except Exception:
            pass

        # If user selected a custom target-snps tolerance, require it to be > 0
        tol_invalid_missing = False
        try:
            tol_auto_ss = bool(st.session_state.get(f"target_snps_tol_auto_{model_key}", True))
            if not tol_auto_ss:
                # read the session-backed tolerance value (percent)
                tol_val = float(st.session_state.get(f"target_snps_tol_{model_key}", 0.0) or 0.0)
                # disallow zero or negative values for custom tolerance
                if tol_val <= 0.0:
                    tol_invalid_missing = True
                    st.error("Target SNPs tolerance must be a positive value greater than 0 when Custom is selected.")
        except Exception:
            tol_invalid_missing = False

        # Require at least one of: explicit Length >0, Target SNPs >0, or Chromosome Length mode
        try:
            chrom_mode = bool(st.session_state.get(f"chrom_length_mode_{model_key}", False))
            seq_val = int(st.session_state.get(f"seq_length_{model_key}", 0) or 0)
            snp_val = int(st.session_state.get(f"target_snps_{model_key}", 0) or 0)
            if not chrom_mode and seq_val <= 0 and snp_val <= 0:
                primary_missing_req = True
                st.error("Set either Length (bp), Target SNPs (>0), or enable Chromosome Length. At least one is required.")
        except Exception:
            pass

        engine_cmd_key = f"engine_cmd_{model_key}"
        # include all validation failures in the effective disabled flag so
        # Prepare/Execute (and Show engine command) are consistently disabled
        # when required inputs or engine-specific checks fail.
        # If sweep parameters are enabled, require a sweep population to be selected.
        sweep_pop_required_missing = False
        try:
            sweep_enabled_ss_local = bool(st.session_state.get(f"sweep_enable_{model_key}", False))
            if sweep_enabled_ss_local:
                # check both selectbox and manual keys
                sp_key_sel = f"sweep_pop_{model_key}_sel"
                sp_key_manual = f"sweep_pop_manual_{model_key}"
                sp_val = ""
                if sp_key_sel in st.session_state:
                    try:
                        sel_val = st.session_state.get(sp_key_sel, "") or ""
                        # map the placeholder back to empty
                        if sel_val == "-- select population --":
                            sel_val = ""
                        sp_val = sel_val
                    except Exception:
                        sp_val = ""
                elif sp_key_manual in st.session_state:
                    sp_val = st.session_state.get(sp_key_manual, "") or ""

                if not str(sp_val).strip():
                    sweep_pop_required_missing = True
                    st.error("Sweep is enabled: choose a Sweep population before preparing/executing the command.")
        except Exception:
            sweep_pop_required_missing = False

        # btn_disabled_effective collects validation failures that disable
        # Prepare/Execute and the copyable CLI. Do NOT include engine==msprime
        # here — msprime should only disable the *Show engine command* button
        # because it runs in-process and doesn't emit an external engine command.
        btn_disabled_effective = (
            btn_disabled
            or sfs_missing_req
            or primary_missing_req
            or sweep_pop_required_missing
            or tol_invalid_missing
        )

        # Separate flag used solely for the Show engine command button so
        # we can disable that button for msprime while leaving execution and
        # the copyable CLI available when other inputs are valid.
        show_engine_disabled = btn_disabled_effective or (engine == "msprime")

        if st.button("Show engine command", key=f"show_engine_cmd_btn_{model_key}", disabled=show_engine_disabled):
            cmd_show = build_cmd()
            if cmd_show is None:
                st.error("Complete steps 1–4 (species, model, ordered populations with counts) before showing the engine command.")
            else:
                # Ask simulator.py to print the engine command (use --show-command)
                cmd_proc = cmd_show + ["--show-command"]
                try:
                    proc = subprocess.run(cmd_proc, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    if proc.returncode != 0:
                        st.error(f"Simulator returned non-zero exit code {proc.returncode} while preparing engine command.")
                        st.subheader("Stderr")
                        st.text(proc.stderr)
                    else:
                        # Trim whitespace and store the engine command as returned by the simulator
                        engine_text = proc.stdout.strip()
                        st.session_state[engine_cmd_key] = engine_text
                except FileNotFoundError:
                    st.error("Could not find simulator.py in the working folder.")
                except Exception as e:
                    st.error(f"Failed to run simulator to retrieve engine command: {e}")

        # NOTE: removed 'Include --show-command' checkbox per user request.

        cmd = build_cmd()
        if cmd is None:
            st.info("Complete steps 1–4 (species, model, ordered populations with counts) to build the command.", icon="ℹ️")
        else:
            # Only display the copyable CLI when execution is available. When
            # execution is disabled (btn_disabled_effective True) we avoid
            # showing the CLI to prevent confusion.
            if not btn_disabled_effective:
                st.markdown("**Copyable CLI:**")
                # When showing the copyable CLI, drop the absolute python executable path
                try:
                    cmd_display = cmd.copy()
                    py_exec = sys.executable
                    if cmd_display and cmd_display[0] == py_exec:
                        cmd_display = cmd_display[1:]
                    st.code(" ".join(shlex.quote(c) for c in cmd_display), language="bash")
                except Exception:
                    st.code(" ".join(shlex.quote(c) for c in cmd), language="bash")

                # If the engine command was prepared via the Show engine command
                # button, display it in the same style as the copyable CLI.
                if st.session_state.get(engine_cmd_key):
                    st.markdown("**Engine command:**")
                    engine_text = st.session_state.get(engine_cmd_key) or ""
                    # Strip leading python executable if present
                    try:
                        py_exec = sys.executable
                        if engine_text and engine_text.startswith(py_exec + " "):
                            engine_text = engine_text[len(py_exec) + 1 :]
                    except Exception:
                        pass

                    # Show the engine command as code (same presentation as copyable CLI)
                    st.code(engine_text, language="bash")

                    # (Custom HTML copy button removed — Streamlit's code block provides copy UI.)
            else:
                # Execution isn't possible right now; avoid showing the copyable CLI.
                st.info("Execution is disabled; copyable CLI is not available.", icon="ℹ️")

            # Execute immediately when the Execute button is pressed
            cancel_key_ui = f"_cancel_run_{model_key}"
            # ensure cancel key exists in session_state for the UI to toggle
            if cancel_key_ui not in st.session_state:
                st.session_state[cancel_key_ui] = False

            running_key = f"_running_{model_key}"
            # Ensure running key exists
            if running_key not in st.session_state:
                try:
                    st.session_state[running_key] = False
                except Exception:
                    pass

            # Render Execute and Cancel side-by-side so both widgets always exist
            # Cancel is enabled only while running; Execute is disabled while running
            # Also disable Execute when an artifacts lock file exists for this model
            try:
                workspace_dir = os.path.abspath(os.getcwd())
                artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                lock_path_check = os.path.join(artifacts_dir, f".runlock_{model_key}")
                lock_exists = os.path.exists(lock_path_check)
            except Exception:
                lock_exists = False

            # Ensure start_request and processed keys exist with safe defaults
            # so stale/missing keys don't permanently disable Execute.
            start_request_key = f"_start_request_{model_key}"
            processed_key = f"_exec_processed_{model_key}"
            ui_lock_key = f"_ui_lock_{model_key}"
            if ui_lock_key not in st.session_state:
                try:
                    st.session_state[ui_lock_key] = False
                except Exception:
                    pass
            if start_request_key not in st.session_state:
                st.session_state[start_request_key] = False
            if processed_key not in st.session_state:
                st.session_state[processed_key] = False

            

            # Compute Execute disabled state. Only rely on running flag and lock
            # presence so stale start_request flags cannot permanently disable Execute.
            # If a run just finished in a different execution context, perform
            # a safe UI-side rerun here: clear the finished flag and request
            # an experimental rerun so widget states are refreshed.
            try:
                finished_key = f"_exec_finished_{model_key}"
                if st.session_state.get(finished_key, False):
                    # clear the finished marker and ensure running/cancel/ui_lock are false
                    st.session_state[finished_key] = False
                    st.session_state[running_key] = False
                    st.session_state[cancel_key_ui] = False
                    st.session_state[ui_lock_key] = False
                    # If a lock file exists for this model, check whether the
                    # PID recorded in it is still alive. If it's stale, remove
                    # the lock so Execute becomes available.
                    try:
                        try:
                            workspace_dir = os.path.abspath(os.getcwd())
                            artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                            lock_path_check = os.path.join(artifacts_dir, f".runlock_{model_key}")
                        except Exception:
                            lock_path_check = None
                        if lock_path_check and os.path.exists(lock_path_check):
                            try:
                                with open(lock_path_check, 'r') as _lf:
                                    content = _lf.read().strip()
                                parts = [p for p in content.split() if p.strip()]
                                pid_existing = int(parts[0]) if parts else None
                            except Exception:
                                pid_existing = None
                            alive = False
                            if pid_existing:
                                try:
                                    if psutil:
                                        alive = psutil.pid_exists(pid_existing)
                                    else:
                                        os.kill(pid_existing, 0)
                                        alive = True
                                except Exception:
                                    alive = False
                            if not alive:
                                try:
                                    os.remove(lock_path_check)
                                except Exception:
                                    pass
                    except Exception:
                        pass
                    try:
                        if hasattr(st, 'experimental_rerun'):
                            st.experimental_rerun()
                    except Exception:
                        pass
            except Exception:
                pass
            # Consume any refresh flag left by run completion when rerun APIs
            # were not available. This ensures we clear UI locks on next render.
            try:
                refresh_key = f"_needs_refresh_{model_key}"
                if st.session_state.get(refresh_key, False):
                    st.session_state[refresh_key] = False
                    st.session_state[running_key] = False
                    st.session_state[cancel_key_ui] = False
                    st.session_state[ui_lock_key] = False
                    # remove stale lock file for this model if PID dead
                    try:
                        workspace_dir = os.path.abspath(os.getcwd())
                        artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                        lock_path_check = os.path.join(artifacts_dir, f".runlock_{model_key}")
                        if os.path.exists(lock_path_check):
                            try:
                                with open(lock_path_check, 'r') as _lf:
                                    content = _lf.read().strip()
                                parts = [p for p in content.split() if p.strip()]
                                pid_existing = int(parts[0]) if parts else None
                            except Exception:
                                pid_existing = None
                            alive = False
                            if pid_existing:
                                try:
                                    if psutil:
                                        alive = psutil.pid_exists(pid_existing)
                                    else:
                                        os.kill(pid_existing, 0)
                                        alive = True
                                except Exception:
                                    alive = False
                            if not alive:
                                try:
                                    os.remove(lock_path_check)
                                except Exception:
                                    pass
                    except Exception:
                        pass
            except Exception:
                pass
            exec_disabled = (
                bool(btn_disabled_effective)
                or bool(st.session_state.get(running_key, False))
                or bool(st.session_state.get(ui_lock_key, False))
                or bool(lock_exists)
            )
            cancel_disabled = not bool(st.session_state.get(running_key, False))
            c_exec, c_cancel = st.columns([1, 1])
            def _set_cancel_local():
                """UI Cancel handler: create a cancel sentinel, try to signal the process/group,
                fallback to psutil cleanup, update session state and remove lock files.
                """
                workspace_dir = os.path.abspath(os.getcwd())
                artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                lock_path = os.path.join(artifacts_dir, f".runlock_{model_key}")
                cancel_path = os.path.join(artifacts_dir, f".runlock_{model_key}.cancel")

                # ensure artifacts dir exists and write a cancel sentinel so the runner can see it
                try:
                    os.makedirs(artifacts_dir, exist_ok=True)
                    with open(cancel_path, 'w') as _cf:
                        _cf.write(f"{os.getpid()}\n{time.time()}\n")
                    try:
                        ui_diag(f"[ui_simulator] UI Cancel wrote cancel sentinel: {cancel_path}")
                    except Exception:
                        pass
                except Exception as e:
                    try:
                        ui_diag(f"[ui_simulator] UI Cancel failed to write sentinel: {e}")
                    except Exception:
                        pass

                # attempt to find pid from session or lock
                pid_val = None
                try:
                    pid_key = f"_running_pid_{model_key}"
                    pid_val = st.session_state.get(pid_key)
                except Exception:
                    pid_val = None
                if not pid_val:
                    try:
                        if os.path.exists(lock_path):
                            with open(lock_path, 'r') as fh:
                                content = fh.read().strip()
                            parts = [p for p in content.split() if p.strip()]
                            pid_val = int(parts[0]) if parts else None
                    except Exception:
                        pid_val = None

                if pid_val:
                    try:
                        try:
                            pg = os.getpgid(int(pid_val))
                        except Exception:
                            pg = None

                        if pg is not None:
                            try:
                                ui_diag(f"[ui_simulator] UI Cancel sending SIGTERM to pgid={pg}")
                                os.killpg(pg, signal.SIGTERM)
                            except Exception as e:
                                ui_diag(f"[ui_simulator] UI Cancel os.killpg failed: {e}")

                        try:
                            os.kill(int(pid_val), signal.SIGTERM)
                        except Exception as e:
                            ui_diag(f"[ui_simulator] UI Cancel os.kill(pid) failed: {e}")

                        time.sleep(0.5)
                        still = False
                        try:
                            if psutil:
                                still = psutil.pid_exists(int(pid_val))
                            else:
                                os.kill(int(pid_val), 0)
                                still = True
                        except Exception:
                            still = False

                        if still and psutil:
                            try:
                                p = psutil.Process(int(pid_val))
                                children = p.children(recursive=True)
                                for c in children:
                                    try:
                                        ui_diag(f"[ui_simulator] UI Cancel terminating child pid={c.pid}")
                                        c.terminate()
                                    except Exception:
                                        pass
                                gone, alive = psutil.wait_procs(children, timeout=3)
                                for c in alive:
                                    try:
                                        ui_diag(f"[ui_simulator] UI Cancel killing child pid={c.pid}")
                                        c.kill()
                                    except Exception:
                                        pass
                                try:
                                    p.kill()
                                except Exception:
                                    pass
                            except Exception as e:
                                ui_diag(f"[ui_simulator] UI Cancel psutil cleanup failed: {e}")
                    except Exception as e:
                        ui_diag(f"[ui_simulator] UI Cancel outer error: {e}")

                # set session flags so UI reflects cancellation
                try:
                    st.session_state[cancel_key_ui] = True
                except Exception:
                    pass
                try:
                    st.session_state[running_key] = False
                except Exception:
                    pass

                # clear start_request/processed flags so Execute becomes available again
                try:
                    start_request_key = f"_start_request_{model_key}"
                    processed_key = f"_exec_processed_{model_key}"
                    st.session_state[start_request_key] = False
                    st.session_state[processed_key] = False
                    # clear transient immediate flag
                    try:
                        immediate_key = f"_exec_start_immediate_{model_key}"
                        st.session_state[immediate_key] = False
                    except Exception:
                        pass
                except Exception:
                    pass

                # remove lock and cancel sentinel if present (cleanup)
                try:
                    if os.path.exists(lock_path):
                        os.remove(lock_path)
                except Exception:
                    pass
                try:
                    if os.path.exists(cancel_path):
                        os.remove(cancel_path)
                except Exception:
                    pass

                # clear start_request/processed/ui_lock so Execute can be pressed again
                try:
                    start_request_key = f"_start_request_{model_key}"
                    processed_key = f"_exec_processed_{model_key}"
                    st.session_state[start_request_key] = False
                    st.session_state[processed_key] = False
                except Exception:
                    pass
                try:
                    ui_lock_key = f"_ui_lock_{model_key}"
                    st.session_state[ui_lock_key] = False
                except Exception:
                    pass

                # Force an immediate rerun so the UI updates state (enables Execute)
                try:
                    if hasattr(st, 'experimental_rerun'):
                        st.experimental_rerun()
                except Exception:
                    pass

            with c_exec:
                # Execute triggers a start-request so the app reruns before
                # starting the long-running job. We set running immediately in
                # the on_click handler so the Execute button is disabled right away.
                start_request_key = f"_start_request_{model_key}"
                if start_request_key not in st.session_state:
                    st.session_state[start_request_key] = False

                def _request_exec_local():
                    try:
                        # clear any stale processed flag so the start block will execute
                        try:
                            processed_key = f"_exec_processed_{model_key}"
                            st.session_state[processed_key] = False
                        except Exception:
                            pass
                        # mark immediate start so the run starter will kick off
                        try:
                            immediate_key = f"_exec_start_immediate_{model_key}"
                            st.session_state[immediate_key] = True
                        except Exception:
                            pass
                        # set a transient UI lock so repeated clicks before file-lock
                        # creation do not spawn multiple runs
                        try:
                            st.session_state[ui_lock_key] = True
                        except Exception:
                            pass
                        # Set the start_request and mark running immediately so
                        # the Cancel button is enabled right away in the UI.
                        try:
                            running_key_local = f"_running_{model_key}"
                            st.session_state[running_key_local] = True
                        except Exception:
                            pass
                        st.session_state[start_request_key] = True
                    except Exception:
                        pass
                    try:
                        # clear any previous cancel flag so the new run starts clean
                        st.session_state[cancel_key_ui] = False
                    except Exception:
                        pass
                    # Do NOT set running here; the actual run will set running when
                    # it begins. Setting running immediately prevents the start
                    # condition below from firing. Log the request for debugging.
                    try:
                        ui_diag(f"[ui_simulator] Execute requested for model_key={model_key} (start_request set)")
                        # Try to force an immediate rerun so the start_request is
                        # processed in the same click (avoids double-click symptom).
                        try:
                            if hasattr(st, 'experimental_rerun'):
                                st.experimental_rerun()
                        except Exception:
                            pass
                    except Exception:
                        pass

                st.button("Execute", disabled=exec_disabled, key=f"exec_btn_{model_key}", on_click=_request_exec_local)
                # Diagnostics: intentionally hide any user-facing message when Execute is disabled
                try:
                    # Keep diagnostics internal-only via ui_diag (controlled by SHOW_DEV_DIAG)
                    pass
                except Exception:
                    pass

                # Developer diagnostics expander: show session_state keys and lock contents
                try:
                    # ensure lock_path_check exists for diagnostics
                    try:
                        workspace_dir = os.path.abspath(os.getcwd())
                        artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                        lock_path_check = os.path.join(artifacts_dir, f".runlock_{model_key}")
                    except Exception:
                        lock_path_check = None

                    if SHOW_DEV_DIAG:
                        with st.expander("Diagnostics (dev)", expanded=False):
                            diag_state = {
                                'btn_disabled_effective': bool(btn_disabled_effective),
                                'exec_disabled': bool(exec_disabled),
                                'running': bool(st.session_state.get(running_key, False)),
                                'start_request': bool(st.session_state.get(start_request_key, False)),
                                'exec_processed': bool(st.session_state.get(processed_key, False)),
                                'running_pid': st.session_state.get(f"_running_pid_{model_key}", None),
                            }

                            # lock file contents
                            lock_info = {'exists': False, 'content': None, 'content_error': None}
                            try:
                                if lock_path_check and os.path.exists(lock_path_check):
                                    lock_info['exists'] = True
                                    try:
                                        with open(lock_path_check, 'r') as _lf:
                                            lock_info['content'] = _lf.read().strip()
                                    except Exception as e:
                                        lock_info['content_error'] = str(e)
                            except Exception:
                                lock_info['exists'] = False

                            st.write(diag_state)
                            st.write(lock_info)

                            try:
                                # allow clearing a stale lock from UI during debugging
                                if lock_info.get('exists'):
                                    if st.button("Clear lock (dev)", key=f"clear_lock_{model_key}"):
                                        try:
                                            if lock_path_check:
                                                os.remove(lock_path_check)
                                        except Exception:
                                            pass
                                        # reset session flags that could have been left stale
                                        try:
                                            st.session_state[processed_key] = False
                                            st.session_state[start_request_key] = False
                                            st.session_state[running_key] = False
                                            st.session_state[cancel_key_ui] = False
                                            # clear transient UI lock too
                                            try:
                                                st.session_state[ui_lock_key] = False
                                            except Exception:
                                                pass
                                        except Exception:
                                            pass
                                        try:
                                            ui_diag(f"[ui_simulator] Dev cleared lock for model_key={model_key}")
                                        except Exception:
                                            pass
                            except Exception:
                                pass
                except Exception:
                    pass
            with c_cancel:
                # Cancel button: only enabled while running. Use on_click to reliably
                # set the cancel session flag (avoids issues with immediate reruns).
                st.button("Cancel", disabled=cancel_disabled, key=f"cancel_btn_{model_key}", on_click=_set_cancel_local)

            # Start the run after rendering the Execute/Cancel buttons so Cancel
            # is visible immediately while the run is starting.
            start_request_key = f"_start_request_{model_key}"
            processed_key = f"_exec_processed_{model_key}"
            immediate_key = f"_exec_start_immediate_{model_key}"
            immediate_flag = bool(st.session_state.get(immediate_key, False))
            if (st.session_state.get(start_request_key, False) and not st.session_state.get(processed_key, False) and not st.session_state.get(running_key, False)) or immediate_flag:
                try:
                    st.session_state[processed_key] = True
                except Exception:
                    pass
                try:
                    # clear the transient immediate flag so we don't re-enter
                    try:
                        if immediate_flag:
                            st.session_state[immediate_key] = False
                    except Exception:
                        pass
                    st.session_state[start_request_key] = False
                except Exception:
                    pass
                try:
                    st.session_state[cancel_key_ui] = False
                except Exception:
                    pass
                try:
                    st.session_state[running_key] = True
                except Exception:
                    pass

                # stash expected outputs for later artifact discovery
                try:
                    out_path_final = None
                    run_sfs_only_ss = bool(st.session_state.get(f"run_sfs_only_{model_key}", False))
                    if out_path and not run_sfs_only_ss:
                        try:
                            fmt = str(out_format or "").strip()
                        except Exception:
                            fmt = ""
                        fmt_ext_map = {"ms": ".ms", "ms.gz": ".ms.gz", "vcf": ".vcf", "vcf.gz": ".vcf.gz", "bcf": ".bcf"}
                        ext = fmt_ext_map.get(fmt, "")
                        out_path_final = out_path
                        if ext and not out_path_final.endswith(ext):
                            out_path_final = out_path_final + ext
                    sfs_arg = None
                    try:
                        if sfs_output and sfs_output.strip():
                            sfs_arg = sfs_output.strip()
                        else:
                            if out_path:
                                base = out_path
                                for e in [".ms.gz", ".vcf.gz"]:
                                    if base.endswith(e):
                                        base = base[:-len(e)]
                                        break
                                else:
                                    base = os.path.splitext(base)[0]
                                sfs_arg = base + ".sfs"
                    except Exception:
                        sfs_arg = None
                    paired_key_local = f"paired_neutral_{model_key}"
                    try:
                        pn_arg = (paired_neutral_name or "").strip() if bool(st.session_state.get(paired_key_local, False)) else ""
                    except Exception:
                        pn_arg = ""
                    st.session_state[f"_expected_out_{model_key}"] = out_path_final
                    st.session_state[f"_expected_sfs_{model_key}"] = sfs_arg
                    st.session_state[f"_expected_pn_{model_key}"] = pn_arg
                except Exception:
                    pass

                # Build runtime command and run (wrap paths into artifacts dir)
                try:
                    cmd_exec = list(cmd) if isinstance(cmd, (list, tuple)) else cmd
                    try:
                        workspace_dir = os.path.abspath(os.getcwd())
                        artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                        os.makedirs(artifacts_dir, exist_ok=True)
                    except Exception:
                        artifacts_dir = os.path.join(os.path.abspath(os.getcwd()), 'artifacts')
                    runtime_cmd = []
                    skip_next = False
                    path_flags = {'--output', '--sfs', '--paired-neutral'}
                    for i, token in enumerate(cmd_exec):
                        if skip_next:
                            skip_next = False
                            try:
                                t = str(token)
                                if os.path.isabs(t) or os.path.commonpath([os.path.abspath(t), artifacts_dir]) == artifacts_dir:
                                    runtime_cmd.append(t)
                                else:
                                    runtime_cmd.append(os.path.join(artifacts_dir, t))
                            except Exception:
                                runtime_cmd.append(token)
                            continue
                        if token in path_flags:
                            runtime_cmd.append(token)
                            skip_next = True
                            continue
                        appended = False
                        for pf in path_flags:
                            if isinstance(token, str) and token.startswith(pf + "="):
                                try:
                                    key, val = token.split('=', 1)
                                    if os.path.isabs(val) or os.path.commonpath([os.path.abspath(val), artifacts_dir]) == artifacts_dir:
                                        runtime_cmd.append(token)
                                    else:
                                        runtime_cmd.append(f"{key}={os.path.join(artifacts_dir, val)}")
                                    appended = True
                                except Exception:
                                    runtime_cmd.append(token)
                                    appended = True
                                break
                        if appended:
                            continue
                        runtime_cmd.append(token)
                    cmd_exec = runtime_cmd
                    if '--progress' not in cmd_exec:
                        cmd_exec = cmd_exec + ['--progress'] if isinstance(cmd_exec, list) else cmd_exec + ['--progress']
                except Exception:
                    cmd_exec = cmd
                # start the run (blocks until finished inside run_and_report)
                run_and_report(cmd_exec)
                try:
                    st.session_state[processed_key] = False
                except Exception:
                    pass
            if st.session_state.get(running_key, False):
                try:
                    st.caption("Running…")
                except Exception:
                    st.markdown("**Running…**")

            # Download panel: show any generated artifacts for this model.
            # Prefer files placed under ./artifacts/ in the workspace for downloads
            # (do not change the copyable CLI paths used by build_cmd()).
            reg_key = f"generated_files_{model_key}"
            gen = st.session_state.get(reg_key, [])
            if gen:
                st.markdown("---")
                st.markdown("**Generated artifacts**")

                # ensure artifacts folder exists in the workspace so users can
                # manually inspect or remove files if needed
                try:
                    workspace_dir = os.path.abspath(os.getcwd())
                    artifacts_dir = os.path.join(workspace_dir, 'artifacts')
                    os.makedirs(artifacts_dir, exist_ok=True)
                except Exception:
                    artifacts_dir = os.path.join(os.path.abspath(os.getcwd()), 'artifacts')

                # helper to compute canonical base name
                def _canon_base_local(path):
                    try:
                        b = os.path.basename(path or "")
                    except Exception:
                        b = str(path or "")
                    for e in ['.ms.gz', '.vcf.gz', '.ms', '.vcf', '.bcf', '.sfs', '.zip']:
                        if b.endswith(e):
                            return b[:-len(e)]
                    return os.path.splitext(b)[0]

                for label, original_path, download_path, zipped in gen:
                    try:
                        # Prefer a matching file under ./artifacts/ when available.
                        chosen_path = download_path
                        chosen_zipped = zipped
                        base_name = os.path.basename(download_path) or os.path.basename(original_path) or label

                        # If a zipped artifact with the canonical base exists in artifacts/, prefer it.
                        try:
                            canon = _canon_base_local(original_path or download_path or base_name)
                        except Exception:
                            canon = None

                        # Candidate zip in artifacts folder (e.g. artifacts/<base>.zip)
                        zip_candidate = None
                        if canon:
                            zip_candidate = os.path.join(artifacts_dir, f"{canon}.zip")
                        # Special handling: if this item is a Paired Neutral candidate
                        # and the run produced multi-replicate VCF/BCF outputs, prefer
                        # artifacts/<paired_neutral_base>.zip when present.
                        try:
                            is_pn = str(label).lower().startswith('paired neutral') or 'paired neutral' in str(label).lower()
                        except Exception:
                            is_pn = False

                        # If PN zip exists, prefer it.
                        if is_pn and zip_candidate and os.path.exists(zip_candidate):
                            chosen_path = zip_candidate
                            chosen_zipped = True
                        else:
                            # Otherwise, if a general zip exists for the canonical base, prefer it.
                            if zip_candidate and os.path.exists(zip_candidate):
                                chosen_path = zip_candidate
                                chosen_zipped = True
                            else:
                                # If the download_path itself is not present, check if an
                                # identically-named file exists inside artifacts/ and use it.
                                try:
                                    art_alt = os.path.join(artifacts_dir, os.path.basename(download_path))
                                    if os.path.exists(art_alt):
                                        chosen_path = art_alt
                                        # if the alternative is a .zip file, mark zipped
                                        chosen_zipped = art_alt.lower().endswith('.zip')
                                except Exception:
                                    pass

                        dl_name = os.path.basename(chosen_path) or os.path.basename(original_path) or label
                        exists = os.path.exists(chosen_path)

                        with st.container():
                            col_a, col_b = st.columns([4, 1])
                            with col_a:
                                if not exists:
                                    st.write(f"{label}  — (not available yet)")
                                else:
                                    st.write(label)
                            with col_b:
                                if exists:
                                    try:
                                        # small files: read into memory for download_button
                                        with open(chosen_path, 'rb') as fh:
                                            data = fh.read()
                                        st.download_button(f"Download {dl_name}", data=data, file_name=dl_name)
                                    except Exception:
                                        st.write("(file not available)")
                                else:
                                    st.button(f"Download {dl_name}", disabled=True)
                    except Exception:
                        # ignore issues enumerating this item
                        pass

                # Manual delete control removed per user request.

        

# (Removed tip about progress bar per user request)
