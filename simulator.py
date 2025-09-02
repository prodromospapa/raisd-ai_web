#!/usr/bin/env python3
import argparse, math, sys, subprocess, shlex, concurrent.futures, gzip, re
import os
import stdpopsim as sps
import demes
from demes.ms import to_ms
from tqdm import tqdm

# (tolerance constant removed; configured via --target-snps-tol argument only)

def present_size_of(deme: demes.Deme):
    for ep in deme.epochs:
        if getattr(ep, 'end_time', None) == 0:
            size = getattr(ep, 'end_size', None) or getattr(ep, 'start_size', None)
            if size is None:
                raise SystemExit(f"ERROR: cannot determine present size for deme {deme.name}")
            return float(size)
    ep = deme.epochs[-1]
    size = getattr(ep, 'end_size', None) or getattr(ep, 'start_size', None)
    if size is None:
        raise SystemExit(f"ERROR: cannot determine present size for deme {deme.name}")
    return float(size)

def reorder_put_first(order_list, first_name):
    if first_name is None:
        return order_list[:]
    if first_name not in order_list:
        raise SystemExit(f"ERROR: requested name '{first_name}' not in list {order_list}")
    return [first_name] + [d for d in order_list if d != first_name]

def harmonic_number(n):
    return sum(1.0/i for i in range(1, n+1)) if n > 0 else 0.0

# ---------- ms/discoal demography helpers ----------

def strip_I_block(ms_cmd):
    toks = ms_cmd.split()
    if '-I' not in toks:
        return ms_cmd
    i = toks.index('-I')
    if i+1 >= len(toks):
        return ms_cmd
    try:
        npop = int(toks[i+1])
    except ValueError:
        return ms_cmd
    j = i + 2 + npop
    return ' '.join(toks[:i] + toks[j:])

def strip_all_growth(ms_cmd):
    toks = ms_cmd.split()
    out = []
    i = 0
    while i < len(toks):
        if toks[i] in ('-G','-g'):
            i += 2
            continue
        out.append(toks[i]); i += 1
    return ' '.join(out)

def remap_indices_ms_1based(ms_cmd, mapping):
    toks = ms_cmd.split(); out = []; i = 0
    while i < len(toks):
        t = toks[i]
        if t == '-n' and i+2 < len(toks):
            try:
                idx = int(float(toks[i+1]))
                out += ['-n', str(mapping.get(idx, idx)), toks[i+2]]
                i += 3; continue
            except Exception: pass
        if t == '-en' and i+3 < len(toks):
            try:
                idx = int(float(toks[i+2]))
                out += ['-en', toks[i+1], str(mapping.get(idx, idx)), toks[i+3]]
                i += 4; continue
            except Exception: pass
        if t == '-m' and i+3 < len(toks):
            try:
                ii = int(float(toks[i+1])); jj = int(float(toks[i+2]))
                out += ['-m', str(mapping.get(ii, ii)), str(mapping.get(jj, jj)), toks[i+3]]
                i += 4; continue
            except Exception: pass
        if t == '-em' and i+4 < len(toks):
            try:
                ii = int(float(toks[i+2])); jj = int(float(toks[i+3]))
                out += ['-em', toks[i+1], str(mapping.get(ii, ii)), str(mapping.get(jj, jj)), toks[i+4]]
                i += 5; continue
            except Exception: pass
        if t == '-ej' and i+3 < len(toks):
            try:
                ii = int(float(toks[i+2])); jj = int(float(toks[i+3]))
                out += ['-ej', toks[i+1], str(mapping.get(ii, ii)), str(mapping.get(jj, jj))]
                i += 4; continue
            except Exception: pass
        out.append(t); i += 1
    return ' '.join(out)

def map_ej_to_ed(ms_cmd):
    toks = ms_cmd.split(); out = []; i = 0
    while i < len(toks):
        if toks[i] == '-ej' and i + 3 < len(toks):
            out += ['-ed', toks[i+1], toks[i+2], toks[i+3]]
            i += 4; continue
        out.append(toks[i]); i += 1
    return ' '.join(out)

def shift_to_0_based(ms_cmd):
    toks = ms_cmd.split(); out = []; i = 0
    while i < len(toks):
        t = toks[i]
        if t == '-en' and i + 3 < len(toks):
            try:
                out += ['-en', toks[i+1], str(int(float(toks[i+2]))-1), toks[i+3]]
                i += 4; continue
            except ValueError: pass
        if t == '-em' and i + 4 < len(toks):
            try:
                out += ['-em', toks[i+1], str(int(float(toks[i+2]))-1), str(int(float(toks[i+3]))-1), toks[i+4]]
                i += 5; continue
            except ValueError: pass
        if t == '-ed' and i + 3 < len(toks):
            try:
                out += ['-ed', toks[i+1], str(int(float(toks[i+2]))-1), str(int(float(toks[i+3]))-1)]
                i += 4; continue
            except ValueError: pass
        if t == '-n' and i + 2 < len(toks):
            try:
                out += ['-n', str(int(float(toks[i+1]))-1), toks[i+2]]
                i += 3; continue
            except ValueError: pass
        out.append(t); i += 1
    return ' '.join(out)

def stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step=1.05):
    en = []
    graph_order = [d.name for d in graph.demes]
    name_to_idx = {name: i+1 for i, name in enumerate(graph_order)}
    for deme in graph.demes:
        i1 = name_to_idx[deme.name]
        for ep in deme.epochs:
            t0, t1 = ep.start_time, ep.end_time
            n0, n1 = ep.start_size, ep.end_size
            if ep.size_function == 'exponential' and t0 > t1:
                fold = n0 / n1 if n1 > 0 else float('inf')
                fold = max(fold, 1.0 + 1e-9)
                steps = max(1, math.ceil(math.log(fold) / math.log(max_fold_per_step)))
                for s in range(1, steps+1):
                    frac = s/steps
                    n_t = n0 * (n1/n0)**frac
                    t_gen = t0 - frac*(t0 - t1)
                    t_ms = t_gen / (4.0 * N0)
                    en.append((t_ms, i1, n_t / N0))
    en.sort(key=lambda x: -x[0])
    return en

def sanitize_ms_demography(ms_demog_str):
    toks = ms_demog_str.split(); out = []; i = 0
    arg_counts = {'-n':2,'-en':3,'-m':3,'-em':4,'-ej':3}
    while i < len(toks):
        t = toks[i]
        if t == '-eg':
            i += 4 if i+3 < len(toks) else 1
            continue
        if t in arg_counts:
            need = arg_counts[t]
            if i + need < len(toks)+1:
                out.extend(toks[i:i+1+need]); i += 1 + need; continue
            else:
                break
        if t.startswith('-'): out.append(t)
        i += 1
    return ' '.join(out)

# ---------- output converters ----------

def ms_like_to_vcf(ms_text, length, chrom='chr1', ploidy=1):
    """Convert single-replicate ms/discoal/msms output to a minimal VCF string.

    Groups consecutive haplotypes into individuals of size=ploidy to form genotypes.
    Assumes haplotypes are emitted in consistent order across loci.
    """
    # Ensure only one replicate (one segsites block)
    seg_blocks = [ln for ln in ms_text.splitlines() if ln.startswith('segsites:')]
    if len(seg_blocks) == 0:
        raise SystemExit('ERROR: Cannot locate segsites line in simulator output for VCF conversion.')
    if len(seg_blocks) > 1:
        raise SystemExit('ERROR: VCF export currently supports only a single replicate; use --replicates 1.')
    lines = ms_text.splitlines()
    # Find segsites line index
    seg_idx = None
    for i,ln in enumerate(lines):
        if ln.startswith('segsites:'):
            seg_idx = i; break
    if seg_idx is None:
        raise SystemExit('ERROR: segsites line not found.')
    try:
        segsites = int(lines[seg_idx].split()[1])
    except Exception:
        raise SystemExit('ERROR: Failed to parse segsites count.')
    if segsites == 0:
        header = [
            '##fileformat=VCFv4.2',
            '##source=stdpopsim2ms.py',
            f'##contig=<ID={chrom},length={length}>',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        ]
        return '\n'.join(header) + '\n'
    if seg_idx + 1 >= len(lines) or not lines[seg_idx+1].startswith('positions:'):
        # ms-like outputs sometimes have '//' before positions; search ahead
        pos_line = None; pos_line_idx = None
        for j in range(seg_idx+1, min(seg_idx+5, len(lines))):
            if lines[j].startswith('positions:'):
                pos_line = lines[j]; pos_line_idx = j; break
        if pos_line is None:
            raise SystemExit('ERROR: positions line not found after segsites line.')
    else:
        pos_line = lines[seg_idx+1]; pos_line_idx = seg_idx+1
    parts = pos_line.split()
    pos_floats = [float(p) for p in parts[1:]]
    if len(pos_floats) != segsites:
        # some simulators label with 'positions:' then floats
        # tolerate mismatch by using min length
        pos_floats = pos_floats[:segsites]
    # Haplotypes start after positions line until blank or next '//' or end
    hap_lines = []
    if pos_line_idx is None:
        raise SystemExit('ERROR: Internal: positions line index unresolved.')
    for ln in lines[pos_line_idx+1:]:
        if not ln.strip():
            break
        if ln.startswith('//'):
            break
        if ln.startswith('segsites:'):
            break
        if any(kw in ln for kw in ('segsites:', 'positions:')):
            break
        # haplotype line: sequence of 0/1
        seq = ln.strip()
        if set(seq) <= {'0','1'}:
            hap_lines.append(seq)
        else:
            # Ignore non 0/1 lines (e.g., comments)
            continue
    if not hap_lines:
        raise SystemExit('ERROR: No haplotype lines parsed for VCF conversion.')
    n_hap = len(hap_lines)
    # Ensure all hap lines same length
    Lset = {len(h) for h in hap_lines}
    if len(Lset) != 1:
        raise SystemExit('ERROR: Inconsistent haplotype lengths in simulator output.')
    seq_len = Lset.pop()
    if seq_len < segsites:
        raise SystemExit('ERROR: Fewer columns in haplotypes than segsites.')
    # Map fractional positions to integer coordinates (1..length)
    used = set()
    pos_ints = []
    for f in pos_floats:
        p = int(round(f * length))
        if p < 1: p = 1
        if p > length: p = length
        # ensure uniqueness by shifting forward
        while p in used and p < length:
            p += 1
        used.add(p)
        pos_ints.append(p)
    # Build VCF
    header = [
        '##fileformat=VCFv4.2',
        '##source=stdpopsim2ms.py',
        f'##contig=<ID={chrom},length={length}>',
        '##INFO=<ID=.,Number=0,Type=Flag,Description="Placeholder">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    ]
    if ploidy < 1:
        ploidy = 1
    if n_hap % ploidy != 0:
        raise SystemExit(f'ERROR: Total haplotypes {n_hap} not divisible by ploidy {ploidy} for VCF conversion.')
    n_ind = n_hap // ploidy
    sample_names = [f'Ind{i+1}' for i in range(n_ind)]
    header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sample_names))
    body_lines = []
    for vidx,(p) in enumerate(pos_ints, start=1):
        ref = 'A'; alt = 'T'
        geno_chars = [hap_lines[h][vidx-1] for h in range(n_hap)]
        # Group into ploidy-sized genotypes (phased)
        genos = []
        for i in range(0, n_hap, ploidy):
            alleles = ['1' if c == '1' else '0' for c in geno_chars[i:i+ploidy]]
            genos.append('|'.join(alleles))
        body_lines.append(f'{chrom}\t{p}\tsnp{vidx}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t' + '\t'.join(genos))
    return '\n'.join(header + body_lines) + '\n'

# ---------- SFS utilities ----------

def _parse_ms_like_replicates(ms_text):
    """Yield haplotype matrix (list of strings) per replicate from ms-like output."""
    lines = ms_text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith('//'):
            i += 1
            continue
        if ln.startswith('segsites:'):
            # parse segsites
            try:
                segsites = int(ln.split()[1])
            except Exception:
                # skip malformed replicate
                i += 1; continue
            pos_line_idx = None
            # look ahead for positions line
            j = i + 1
            while j < len(lines) and j <= i + 10:
                if lines[j].startswith('positions:'):
                    pos_line_idx = j; break
                if lines[j].startswith('segsites:'):
                    break
                j += 1
            if pos_line_idx is None:
                i += 1; continue
            hap_start = pos_line_idx + 1
            hap_lines = []
            k = hap_start
            while k < len(lines):
                hl = lines[k].strip()
                if not hl or hl.startswith('//') or hl.startswith('segsites:'):
                    break
                if hl.startswith('positions:'):
                    break
                if set(hl) <= {'0','1'}:
                    hap_lines.append(hl)
                k += 1
            if hap_lines:
                yield hap_lines
            i = k
        else:
            i += 1

def compute_sfs(ms_text, normalized=False, mode='mean'):
    """Compute unfolded SFS from ms-like output.

    Returns (header_line, data_lines) where data_lines is list of strings.
    mode: 'mean' averages across replicates; 'per-rep' outputs each replicate separately.
    """
    replicates = list(_parse_ms_like_replicates(ms_text))
    if not replicates:
        return ('#SFS', ['# No replicates parsed'])
    # assume all replicates have same haplotype count
    n_hap = len(replicates[0])
    for haps in replicates:
        if len(haps) != n_hap:
            raise SystemExit('ERROR: Inconsistent haplotype counts across replicates for SFS.')
    bins = n_hap - 1  # k=1..n_hap-1
    if mode == 'per-rep':
        lines = []
        for ridx, haps in enumerate(replicates, start=1):
            if not haps:
                continue
            L = len(haps[0])
            # verify sequences length
            for hp in haps:
                if len(hp) != L:
                    raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
            sfs = [0]*(n_hap+1)
            segsites = L
            for col in range(L):
                k = sum(1 for hp in haps if hp[col] == '1')
                if 0 < k < n_hap:
                    sfs[k] += 1
            if segsites == 0:
                norm = [0]*(n_hap-1)
            else:
                if normalized:
                    norm = [sfs[k]/segsites for k in range(1,n_hap)]
                else:
                    norm = [sfs[k] for k in range(1,n_hap)]
            lines.append('rep'+str(ridx)+'\t'+'\t'.join(str(x) for x in norm))
        header = '#SFS\t' + '\t'.join(str(k) for k in range(1,n_hap))
        return header, lines
    else:  # mean
        accum = [0.0]*(n_hap-1)
        count_used = 0
        for haps in replicates:
            if not haps:
                continue
            L = len(haps[0])
            for hp in haps:
                if len(hp) != L:
                    raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
            sfs = [0]*(n_hap+1)
            segsites = L
            if segsites == 0:
                continue
            for col in range(L):
                k = sum(1 for hp in haps if hp[col] == '1')
                if 0 < k < n_hap:
                    sfs[k] += 1
            if normalized:
                row = [sfs[k]/segsites for k in range(1,n_hap)]
            else:
                row = [sfs[k] for k in range(1,n_hap)]
            accum = [a + r for a,r in zip(accum,row)]
            count_used += 1
        if count_used == 0:
            mean_vals = [0]*(n_hap-1)
        else:
            mean_vals = [a / count_used for a in accum]
        header = '#SFS\t' + '\t'.join(str(k) for k in range(1,n_hap))
        line = 'mean\t' + '\t'.join(str(x) for x in mean_vals)
        return header, [line]

def compute_sfs_fast(ms_text, normalized=False, mode='mean'):
    """Accelerated SFS computation using NumPy where available.

    Falls back to compute_sfs if NumPy is not installed. Interface-compatible.
    """
    try:
        import numpy as _np  # local import to avoid hard dependency
    except Exception:
        return compute_sfs(ms_text, normalized=normalized, mode=mode)

    replicates = list(_parse_ms_like_replicates(ms_text))
    if not replicates:
        return ('#SFS', ['# No replicates parsed'])
    n_hap = len(replicates[0])
    for haps in replicates:
        if len(haps) != n_hap:
            raise SystemExit('ERROR: Inconsistent haplotype counts across replicates for SFS.')
    header = '#SFS\t' + '\t'.join(str(k) for k in range(1, n_hap))
    if mode == 'per-rep':
        out_lines = []
        for ridx, haps in enumerate(replicates, start=1):
            if not haps:
                continue
            L = len(haps[0])
            for hp in haps:
                if len(hp) != L:
                    raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
            if L == 0:
                vals = [0]*(n_hap-1)
            else:
                # Build contiguous byte buffer for all haplotypes and reshape
                buf = ''.join(haps).encode('ascii')
                mat = _np.frombuffer(buf, dtype='S1').reshape(n_hap, L)
                counts = (mat == b'1').sum(axis=0)  # per-site derived allele count
                binc = _np.bincount(counts, minlength=n_hap+1)
                raw = binc[1:n_hap]  # exclude monomorphic counts
                if normalized:
                    denom = raw.sum()
                    if denom > 0:
                        vals = (raw / denom).tolist()
                    else:
                        vals = [0.0]*(n_hap-1)
                else:
                    vals = raw.tolist()
            out_lines.append(f"rep{ridx}\t" + '\t'.join(str(v) for v in vals))
        return header, out_lines
    else:  # mean mode
        accum = _np.zeros(n_hap-1, dtype=float)
        used = 0
        for haps in replicates:
            if not haps:
                continue
            L = len(haps[0])
            for hp in haps:
                if len(hp) != L:
                    raise SystemExit('ERROR: Inconsistent haplotype sequence lengths.')
            if L == 0:
                continue
            buf = ''.join(haps).encode('ascii')
            mat = _np.frombuffer(buf, dtype='S1').reshape(n_hap, L)
            counts = (mat == b'1').sum(axis=0)
            binc = _np.bincount(counts, minlength=n_hap+1)
            raw = binc[1:n_hap]
            if normalized:
                denom = raw.sum()
                row = raw/denom if denom > 0 else _np.zeros_like(raw)
            else:
                row = raw
            accum += row
            used += 1
        if used == 0:
            mean_vals = [0]*(n_hap-1)
        else:
            mean_vals = (accum / used).tolist()
        return header, ['mean\t' + '\t'.join(str(v) for v in mean_vals)]

# ---------- discoal builder ----------

def build_discoal_command(*, species, model_id, user_order, individual_counts, discoal_pop0,
                          reps, length, max_fold_per_step, sweep_time, x, a, min_snps, chr_name=None,
                          disable_en_ladder=False, debug=False, seed=None):
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try:
        graph = raw_demog.to_demes()
    except AttributeError:
        graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')

    graph_names = [d.name for d in graph.demes]
    ploidy = sp.ploidy
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if discoal_pop0 is None:
        discoal_pop0 = user_order[0]
    desired_disc_order = reorder_put_first(user_order, discoal_pop0)
    N0 = deme_present_sizes[discoal_pop0]

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try:
            chosen_contig = sp.get_contig(chr_name)
        except Exception:
            chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig, 'mutation_rate'):
                    mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig, 'recombination_map') and hasattr(chosen_contig.recombination_map, 'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    counts_disc = [name_to_hap[n] for n in desired_disc_order]

    total_hap = sum(counts_disc)
    H = harmonic_number(total_hap - 1)
    if length is None:
        if min_snps is not None:
            # discoal theta later set as theta = 2 * N0 * mu * ploidy * L
            # Expected segregating sites S ≈ theta * H (same harmonic expectation as ms) =>
            # S ≈ (2 * N0 * mu * ploidy * L) * H
            # Solve L = min_snps / (2 * N0 * mu * ploidy * H)
            theta_factor = 2.0 * N0 * mu * ploidy
            if theta_factor * H <= 0:
                raise SystemExit('ERROR: Non-positive theta factor * H, cannot derive length.')
            safety = 1.15  # bias upward so realized S more likely >= target
            length = int(math.ceil(min_snps / (theta_factor * H) * safety))
        else:
            raise SystemExit('ERROR: Provide --length or --min-snps.')

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n, 0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_no_growth = strip_all_growth(ms_no_I)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_disc_idx = {n: i+1 for i,n in enumerate(desired_disc_order)}
    mapping = {name_to_graph_idx[n]: name_to_disc_idx[n] for n in desired_disc_order}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_remap = map_ej_to_ed(demog_remap)
    demog_0b = shift_to_0_based(demog_remap)

    if chr_name:
        theta = 2.0 * N0 * mu * ploidy * length
        rho   = 2.0 * N0 * rrate * ploidy * length
    else:
        theta = 2.0 * ploidy * N0 * mu * length
        rho   = 2.0 * ploidy * N0 * rrate * length

    sel_args = []
    # discoal -ws t x a : t interpreted here as ORIGIN time (when allele became beneficial) per user spec.
    if (a is not None) and (sweep_time is not None):
        if x is None:
            x = 0.5
        sel_args = ['-ws', str(sweep_time), '-x', str(x), '-a', str(a)]

    npop_disc = len(desired_disc_order)
    total_hap_disc = sum(counts_disc)
    seed_part = f" -seed {seed}" if seed is not None else ''
    discoal_cmd = (
        f"# pop0 (sweep-capable) deme: {discoal_pop0}; indices 0-based; order: {desired_disc_order}\n"
        f"discoal {total_hap_disc} {reps} {length} -t {theta} -r {rho} -p {npop_disc} {' '.join(map(str,counts_disc))}{seed_part} "
        f"{' '.join(sel_args)} {demog_0b}".strip()
    )
    meta = {
        'N0': N0, 'mu': mu, 'rrate': rrate, 'counts_disc': counts_disc, 'npop': npop_disc,
        'demog': demog_0b, 'pop_order': desired_disc_order, 'pop0': discoal_pop0,
        'ploidy': ploidy, 'sel_args': sel_args, 'length': length,
        'chromosome': chr_name if chr_name else None, 'species_id': species,
        'sweep_pos': x if sel_args else None,
        'sel_2Ns': a if sel_args else None,
        'sweep_time': sweep_time if sel_args else None,
        'fixation_time': None,
        'base_comments': [f"# pop0 (sweep-capable) deme: {discoal_pop0}; indices 0-based; order: {desired_disc_order}",
                          f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    if sel_args and x is not None:
        try:
            sweep_bp = int(round(float(x) * length)) if length else None
        except Exception:
            sweep_bp = None
        loc_line = f"# sweep-pos={x}"
        if sweep_bp is not None:
            loc_line += f" sweep-bp={sweep_bp}"
        if sweep_time is not None:
            loc_line += f" sweep-time={sweep_time}"
        meta['base_comments'].append(loc_line)
    if debug:
        sys.stderr.write('# DEBUG discoal_demog_ms(after_remap_0b) '+demog_0b+'\n')
    if seed is not None:
        meta['seed'] = seed
    return discoal_cmd, meta

# ---------- ms (neutral) builder ----------

def build_ms_command(*, species, model_id, user_order, individual_counts, pop0,
                     reps, length, max_fold_per_step, chr_name=None,
                     min_snps=None, disable_en_ladder=False, debug=False, seed=None):
    engine = 'ms'
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try: graph = raw_demog.to_demes()
    except AttributeError: graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')
    ploidy = sp.ploidy
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if pop0 is None: pop0 = user_order[0]
    if pop0 not in deme_present_sizes:
        raise SystemExit(f"ERROR: pop0 '{pop0}' not in model demes {list(deme_present_sizes)}")
    N0 = deme_present_sizes[pop0]

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try: chosen_contig = sp.get_contig(chr_name)
        except Exception: chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig,'mutation_rate'): mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig,'recombination_map') and hasattr(chosen_contig.recombination_map,'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    total_hap_placeholder = sum(ploidy * c for c in individual_counts)
    if length is None:
        if min_snps is not None:
            H = harmonic_number(total_hap_placeholder - 1)
            theta_per_bp = 4.0 * N0 * mu  # theta = 4 N0 mu L
            if theta_per_bp * H <= 0:
                raise SystemExit('ERROR: Cannot derive length for ms (theta_per_bp*H<=0).')
            safety = 1.10
            length = int(math.ceil(min_snps / (theta_per_bp * H) * safety))
        elif chr_name and chosen_contig is not None:
            length = int(getattr(chosen_contig,'length'))
        else:
            raise SystemExit('ERROR: Provide --length or --min-snps for engine=ms/msms.')

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    desired_order = user_order[:]
    counts_order = [name_to_hap[n] for n in desired_order]
    total_hap = sum(counts_order)

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n,0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    ms_no_growth = strip_all_growth(ms_no_I)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_user_idx = {n: i+1 for i,n in enumerate(desired_order)}
    mapping = {name_to_graph_idx[n]: name_to_user_idx[n] for n in desired_order if n in name_to_graph_idx}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_final = sanitize_ms_demography(demog_remap)

    theta = 4.0 * N0 * mu * length
    rho = 4.0 * N0 * rrate * length
    npop = len(desired_order)

    base_parts = [engine, str(total_hap), str(reps),
                  '-t', str(theta), '-r', str(rho), str(length),
                  '-I', str(npop)] + list(map(str, counts_order))
    # Intentionally NOT adding -seeds to avoid external seed file creation per user request.
    base_cmd_line = ' '.join(base_parts) + ' ' + demog_final
    full_cmd = (f"# engine={engine}\n# pop0={pop0}\n# order={desired_order}\n# N0={N0} mu={mu} r={rrate} length={length}\n{base_cmd_line}")
    meta = {
        'engine': engine, 'N0': N0, 'mu': mu, 'rrate': rrate,
        'counts_disc': counts_order, 'npop': npop, 'demog': demog_final,
        'pop_order': desired_order, 'pop0': pop0, 'ploidy': ploidy,
        'sel_args': [], 'length': length,
        'chromosome': chr_name if chr_name else None,
        'species_id': species,
        'sweep_pos': None,
        'sel_2Ns': None,
        'sweep_time': None,
        'fixation_time': None,
    'base_comments': ['# engine='+engine, f"# pop0={pop0}", f"# order={desired_order}",
              f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    if debug:
        sys.stderr.write('# DEBUG ms demog(after_remap) '+demog_final+'\n')
    if seed is not None:
        meta['seed'] = seed
    return full_cmd, meta

# ---------- msms (selection) builder with SFC + Smark ----------

def build_msms_command(*, species, model_id, user_order, individual_counts, pop0,
                       reps, length, max_fold_per_step, chr_name=None,
                       min_snps=None, disable_en_ladder=False, a=None, x=None,
                       fixation_time=0.0, debug=False, seed=None):
    sp = sps.get_species(species)
    model = sp.get_demographic_model(model_id)
    raw_demog = model.model
    try: graph = raw_demog.to_demes()
    except AttributeError: graph = raw_demog
    if not hasattr(graph, 'demes'):
        raise SystemExit('ERROR: Could not obtain demes.Graph from demographic model.')

    ploidy = sp.ploidy
    if pop0 is None:
        pop0 = user_order[0]
    deme_present_sizes = {d.name: present_size_of(d) for d in graph.demes}
    if pop0 not in deme_present_sizes:
        raise SystemExit(f"ERROR: pop0 '{pop0}' not in model demes {list(deme_present_sizes)}")
    N0 = float(deme_present_sizes[pop0])

    genome = sp.genome
    mu = None; rrate = None; chosen_contig = None
    if chr_name:
        try: chosen_contig = sp.get_contig(chr_name)
        except Exception: chosen_contig = None
        if chosen_contig is not None:
            try:
                if hasattr(chosen_contig,'mutation_rate'): mu = float(chosen_contig.mutation_rate)
            except Exception: pass
            try:
                if hasattr(chosen_contig,'recombination_map') and hasattr(chosen_contig.recombination_map,'mean_rate'):
                    rrate = float(chosen_contig.recombination_map.mean_rate)
            except Exception: pass
    if mu is None:
        for attr in ('per_site_mutation_rate','mutation_rate','mu'):
            if hasattr(genome, attr):
                try: mu = float(getattr(genome, attr)); break
                except Exception: pass
    if rrate is None:
        for attr in ('per_site_recombination_rate','recombination_rate','r'):
            if hasattr(genome, attr):
                try: rrate = float(getattr(genome, attr)); break
                except Exception: pass
    if mu is None: mu = 1e-8
    if rrate is None: rrate = 1e-8

    total_hap_placeholder = sum(ploidy * c for c in individual_counts)
    if length is None:
        if min_snps is not None:
            H = harmonic_number(total_hap_placeholder - 1)
            theta_per_bp = 4.0 * N0 * mu  # neutral expectation
            if theta_per_bp * H <= 0:
                raise SystemExit('ERROR: Cannot derive length for msms (theta_per_bp*H<=0).')
            # Under selective sweeps segregating sites are reduced; inflate length more aggressively
            base_len = min_snps / (theta_per_bp * H)
            safety = 1.25 if (a is not None and a > 0) else 1.10
            length = int(math.ceil(base_len * safety))
        elif chr_name and chosen_contig is not None:
            length = int(getattr(chosen_contig,'length'))
        else:
            raise SystemExit('ERROR: Provide --length or --min-snps for engine=msms.')

    hap_counts_input_order = [ploidy * c for c in individual_counts]
    name_to_hap = {n: hap_counts_input_order[i] for i,n in enumerate(user_order)}
    desired_order = user_order[:]
    counts_order = [name_to_hap[n] for n in desired_order]
    total_hap = sum(counts_order)

    graph_order = [d.name for d in graph.demes]
    name_to_hap_all = {n: (ploidy * individual_counts[user_order.index(n)]) for n in user_order}
    samples_graph_order = [name_to_hap_all.get(n,0) for n in graph_order]
    ms_cmd = to_ms(graph, N0=N0, samples=samples_graph_order)
    ms_no_I = strip_I_block(ms_cmd)
    ms_no_growth = strip_all_growth(ms_no_I)
    en_ladder = [] if disable_en_ladder else stepwise_from_exponential_epochs_graph_order(graph, N0, max_fold_per_step)
    en_ladder_str = ' '.join(f"-en {t} {i} {n}" for t,i,n in en_ladder)
    ms_aug = (ms_no_growth + (' ' + en_ladder_str if en_ladder_str else '')).strip()
    name_to_graph_idx = {n: i+1 for i,n in enumerate(graph_order)}
    name_to_user_idx = {n: i+1 for i,n in enumerate(desired_order)}
    mapping = {name_to_graph_idx[n]: name_to_user_idx[n] for n in desired_order if n in name_to_graph_idx}
    demog_remap = remap_indices_ms_1based(ms_aug, mapping)
    demog_final = sanitize_ms_demography(demog_remap)

    # time-invariant detection (for completeness; we will prefer SFC anyway)
    time_invariant = not any(flag in (' ' + demog_final + ' ')
                             for flag in (' -en ', ' -em ', ' -ej ', ' -eg ', ' -eG ', ' -es ', ' -ema '))

    theta = 4.0 * N0 * mu * length
    rho   = 4.0 * N0 * rrate * length
    npop  = len(desired_order)

    parts = ['msms', str(total_hap), str(reps),
             '-t', str(theta), '-r', str(rho), str(length),
             '-N', str(int(round(N0))),
             '-I', str(npop), *map(str, counts_order)]

    sel_args = []
    if (a is not None) or (x is not None):
        if a is None:
            raise SystemExit('ERROR: engine=msms sweep requested but --sel-2Ns missing.')
        # genic (h≈0.5)
        sel_args += ['-SaA', str(a), '-SAA', str(2.0 * a)]
        # selected site position + mark it in output
        if x is None:
            x = 0.5
        sel_args += ['-Sp', str(x), '-Smark']
        # Derive origin time (backwards) given fixation at present (fixation_time must be 0 per spec).
        # Approx fixation duration in coalescent units (duration from origin to fixation):
        # T_duration ≈ ln(a)/a for a>1 else ln(1+a)/a. (a = 2Ns)
        # If fixation_time > 0 (fixed in the past), origin_time = fixation_time + T_duration.
        if a <= 0:
            raise SystemExit('ERROR: --sel-2Ns must be > 0 for msms selection.')
        if a > 1.0:
            t_duration = math.log(a) / a
        else:
            t_duration = math.log(1.0 + a) / a
        tstart = t_duration + (fixation_time or 0.0)
        f0 = max(1.0 / (2.0 * N0), 1e-6)
        init_freqs = ['0'] * npop
        sel_deme_index = desired_order.index(pop0)  # 0-based index into our order list
        init_freqs[sel_deme_index] = str(f0)
        sel_args += ['-SI', str(tstart), str(npop), *init_freqs, '-SFC']

    if sel_args:
        parts.extend(sel_args)
    parts.extend(demog_final.split())

    if seed is not None:
        parts.extend(['-seed', str(seed)])
    full_cmd = (f"# engine=msms\n# pop0={pop0}\n# order={desired_order}\n# N0={N0} mu={mu} r={rrate} length={length}\n{' '.join(parts)}")
    meta = {
        'engine': 'msms', 'N0': N0, 'mu': mu, 'rrate': rrate,
        'counts_disc': counts_order, 'npop': npop, 'demog': demog_final,
        'pop_order': desired_order, 'pop0': pop0, 'ploidy': ploidy,
        'sel_args': sel_args, 'length': length,
        'chromosome': chr_name if chr_name else None,
        'species_id': species,
        'sweep_pos': x if sel_args else None,
        'sel_2Ns': a if sel_args else None,
        'sweep_time': None,
        'fixation_time': fixation_time if sel_args else None,
    'base_comments': ['# engine=msms', f"# pop0={pop0}", f"# order={desired_order}",
              f"# theta={theta} rho={rho} length={length} N0={N0} mu={mu} r={rrate} ploidy={ploidy}"],
    }
    if sel_args and x is not None:
        try:
            sweep_bp = int(round(float(x) * length)) if length else None
        except Exception:
            sweep_bp = None
        loc_line = f"# sweep-pos={x}"
        if sweep_bp is not None:
            loc_line += f" sweep-bp={sweep_bp}"
        if fixation_time:
            loc_line += f" fixation-time={fixation_time}"
        meta['base_comments'].append(loc_line)
    if debug:
        sys.stderr.write('# DEBUG msms demog(after_remap) '+demog_final+'\n')
    if seed is not None:
        meta['seed'] = seed
    return full_cmd, meta

# ---------- CLI ----------

class _RequiredAnnotatingFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """HelpFormatter that appends '(required)' to help text for required options."""
    def add_argument(self, action):
        if action.required:
            if action.help:
                if '(required)' not in action.help:
                    action.help += ' (required)'
            else:
                action.help = '(required)'
        super().add_argument(action)

def parse_args():
    """Unified parser: always show selection arguments; validate post parse based on --engine."""
    ap = argparse.ArgumentParser(
        description=("Export stdpopsim demography to ms-format using: discoal (selection origin), ms (neutral), msms (selection fixation present)."),
        formatter_class=_RequiredAnnotatingFormatter,
    )
    g = ap.add_argument_group('Core')
    g.add_argument('--engine', required=True, choices=['discoal','ms','msms'],
                   help="Simulator engine: discoal (selection origin), ms (neutral), msms (selection fixation).")
    g.add_argument('--species-id', dest='species', required=True,
                   help="Species identifier or name as accepted by stdpopsim (e.g. 'Homo sapiens').")
    g.add_argument('--model-id', dest='model', required=True,
                   help="Demographic model id for the chosen species (see stdpopsim models).")
    g.add_argument('--pop-order', dest='demes', required=True,
                   help="Comma list of deme/population names specifying the sample order (e.g. 'popA,popB').")
    g.add_argument('--sample-individuals', dest='samples_config', required=True,
                   help="Comma integers: individuals sampled per deme in same order as --pop-order.")
    g.add_argument('--length', type=int, default=None,
                   help="Sequence length (bp). If omitted, may be derived from --target-snps or contig info.")
    g.add_argument('--target-snps', dest='min_snps', type=int, default=None,
                   help="Target/minimum segregating sites; derives --length when not provided.")
    g.add_argument('--target-snps-tol', dest='min_snps_tol', type=float, default=0.0,
                   help=("Fractional overshoot tolerance for --target-snps. After all replicates reach >= target, "
                         "if max(segsites) > target*(1+tol) the sequence length is iteratively shrunk while keeping every replicate >= target. "
                         "Set 0 to disable shrink (default 0). Example: 0.02 allows ~2% overshoot before shrink attempts."))
    g.add_argument('--chromosome', dest='chromosome', default=None,
                   help="Chromosome/contig name for mutation/recombination rates/length (optional).")
    g.add_argument('--replicates', dest='reps', type=int, default=1,
                   help="Number of simulation replicates to run.")
    g.add_argument('--output', dest='out', default=None,
                   help=("Output file path. If omitted, modelid_pop0[_chrom][_sweeppop].EXT is written (EXT per --output-format). "
                         "Use '-' for stdout."))
    g.add_argument('--output-format', dest='format', default='ms', choices=['ms','ms.gz','vcf','vcf.gz','bcf'],
                   help="Output format. 'ms'/'ms.gz' write raw ms-like (optionally gzipped); 'vcf'/'vcf.gz' convert single replicate to VCF.")
    g.add_argument('--threads', dest='threads', type=int, default=1,
                   help="Parallel threads for replicate chunks.")
    g.add_argument('--debug', action='store_true', help="Enable debug logging to stderr.")
    g.add_argument('--progress', action='store_true', help="Show progress bar (requires tqdm).")
    g.add_argument('--sfs', nargs='?', const=True, default=False,
                   help="Compute unfolded SFS. Optional value sets output basename; default name modelid_pop0[_chrom][_sweeppop].sfs.")
    g.add_argument('--sfs-normalized', action='store_true', help='Normalize SFS counts to frequencies.')
    g.add_argument('--sfs-mode', choices=['mean','per-rep'], default='mean', help="SFS aggregation mode.")

    pn = ap.add_argument_group('Paired neutral baseline')
    pn.add_argument('--paired-neutral', action='store_true', help='Also run a neutral simulation with identical parameters (selection removed).')
    pn.add_argument('--neutral-engine', choices=['discoal','ms','msms'], help='Engine to use for neutral run (default: same as --engine).')
    pn.add_argument('--neutral-output', dest='neutral_out', help='Neutral output path. Default: main output name with _neutral before extension.')
    pn.add_argument('--neutral-same-seed', action='store_true', help='Reuse primary run seed for neutral run when engines match; otherwise a different random seed is used.')

    gd = ap.add_argument_group('Demography')
    gd.add_argument('--disable-growth-discretization', dest='no_en_ladder', action='store_true', help='Disable exponential growth discretization.')
    gd.add_argument('--growth-max-fold', dest='max_fold_per_step', type=float, default=1.05, help='Max fold change per discretized growth step.')

    sg = ap.add_argument_group('Selection (engine-specific), required for discoal/msms')
    sg.add_argument('--sweep-pop', dest='discoal_pop0', help='Sweep population (discoal/msms). Ignored for ms.')
    sg.add_argument('--sweep-pos', dest='x', type=float, help='Selected site position in (0,1) (required for discoal/msms). Ignored for ms.')
    sg.add_argument('--sel-2Ns', dest='a', type=float, help='Selection strength 2Ns (required for discoal/msms). Ignored for ms.')
    sg.add_argument('--sweep-time', dest='sweep_time', type=float, help='Origin time (4N gens) allele becomes beneficial (discoal only).')
    sg.add_argument('--fixation-time', dest='fix_time', type=float, default=0.0, help='Fixation time back in 4N units (msms only; default 0=present).')

    args = ap.parse_args()
    engine = args.engine
    argv_full = sys.argv[1:]
    # Misuse checks
    if engine == 'msms' and any(a.startswith('--sweep-time') for a in argv_full):
        raise SystemExit('ERROR: engine=msms uses --fixation-time; do not use --sweep-time.')
    if engine == 'discoal' and any(a.startswith('--fixation-time') for a in argv_full):
        raise SystemExit('ERROR: engine=discoal uses --sweep-time; do not use --fixation-time.')
    # Required selection params
    if engine == 'discoal':
        if args.a is None or args.x is None or args.sweep_time is None:
            raise SystemExit('ERROR: discoal requires --sel-2Ns, --sweep-pos, and --sweep-time.')
    elif engine == 'msms':
        if args.a is None or args.x is None:
            raise SystemExit('ERROR: msms requires --sel-2Ns and --sweep-pos.')
        if args.sweep_time is not None:
            raise SystemExit('ERROR: msms uses --fixation-time not --sweep-time.')
    else:  # ms
        sel_flags = ['--sweep-pos', '--sel-2Ns', '--sweep-time', '--fixation-time', '--sweep-pop']
        provided = [fl for fl in sel_flags if any(a == fl or a.startswith(fl + '=') for a in argv_full)]
        if provided:
            sys.stderr.write('# INFO: Ignoring selection flags ' + ', '.join(provided) + ' because engine=ms is neutral only. Use --engine discoal or --engine msms for selection.\n')
    # target-snps tolerance dependency (no hard error: just ignore if target absent)
    tol_flag_provided = any(a == '--target-snps-tol' or a.startswith('--target-snps-tol=') for a in argv_full)
    if args.min_snps is None and tol_flag_provided:
        # Inform user and neutralize tolerance value
        sys.stderr.write('# INFO: Ignoring --target-snps-tol (no --target-snps provided). It only applies when deriving/ enforcing target SNP count.\n')
        setattr(args, 'min_snps_tol', 0.0)
    # SFS dependency (convert prior hard error into informational ignore)
    have_sfs = bool(getattr(args, 'sfs', False))
    sfs_norm_flag = any(a == '--sfs-normalized' or a.startswith('--sfs-normalized=') for a in argv_full)
    sfs_mode_flag = any(a == '--sfs-mode' or a.startswith('--sfs-mode=') for a in argv_full)
    if not have_sfs and (sfs_norm_flag or sfs_mode_flag):
        sys.stderr.write('# INFO: Ignoring --sfs-normalized/--sfs-mode (require --sfs). Provide --sfs to enable SFS calculation.\n')
        # Leave defaults; ensure normalization flag off
        if hasattr(args, 'sfs_normalized'):
            setattr(args, 'sfs_normalized', False)
    return args

def main():
    args = parse_args()
    user_order = [d.strip() for d in args.demes.split(',') if d.strip()]

    # Automatic primary seed (used for reproducibility and optionally neutral run)
    # Primary seed retained only for optional reuse; individual parallel chunks get unique seeds later.
    try:
        import random as _rnd
        primary_seed = _rnd.randint(1, 2**31 - 1)
    except Exception:
        primary_seed = None

    # Resolve species name -> id if needed
    try:
        _species_map = {sp.name: sp.id for sp in sps.all_species()}
    except Exception:
        _species_map = {}
    if _species_map:
        parts = args.species.strip().split()
        if len(parts) == 2:
            args.species = parts[0].title() + ' ' + parts[1].lower()
        if args.species in _species_map.values():
            pass
        elif args.species in _species_map:
            args.species = _species_map[args.species]
        else:
            sys.exit(f"ERROR: Unknown species '{args.species}'.")

    # threads
    try: max_cpu = os.cpu_count() or 1
    except Exception: max_cpu = 1
    args.threads = min(max(1, args.threads), max_cpu)

    # validate model & complete demes list
    try:
        species_obj = sps.get_species(args.species)
    except Exception as e:
        sys.exit(f"ERROR: Unknown species '{args.species}': {e}")
    valid_models = [m.id for m in species_obj.demographic_models]
    if args.model not in valid_models:
        sys.exit(f"ERROR: Unknown model '{args.model}'. Available models: {','.join(valid_models)}")
    try:
        model_obj = species_obj.get_demographic_model(args.model)
        raw_demog_m = model_obj.model
        try: graph_m = raw_demog_m.to_demes()
        except AttributeError: graph_m = raw_demog_m
        model_demes = [d.name for d in getattr(graph_m, 'demes', [])]
    except Exception:
        model_demes = []
    missing = []
    if model_demes:
        unknown_user = [d for d in user_order if d not in model_demes]
        if unknown_user:
            sys.exit(f"ERROR: Unknown deme(s) {unknown_user}. Available: {','.join(model_demes)}")
        missing = [d for d in model_demes if d not in user_order]
        if missing:
            user_order = user_order + missing

    # parse samples
    raw_sample_tokens = [t.strip() for t in args.samples_config.split(',') if t.strip()!='']
    try:
        individual_counts = [int(x) for x in raw_sample_tokens]
    except ValueError:
        sys.exit('ERROR: --sample-individuals must be a comma list of integers.')
    if len(individual_counts) > len(user_order):
        sys.exit('ERROR: More sample counts provided than demes.')
    padded = False
    if len(individual_counts) < len(user_order):
        individual_counts += [0] * (len(user_order) - len(individual_counts))
        padded = True

    engine = args.engine

    # sweep flag validation (engine specific)
    if engine == 'discoal':
        if any(getattr(args, n) is not None for n in ('sweep_time','a','x')):
            if args.sweep_time is None or args.a is None:
                sys.exit('ERROR: discoal sweep requires --sweep-time and --sel-2Ns.')
            if args.x is None:
                args.x = 0.5
            if args.sweep_time < 0:
                sys.exit('ERROR: --sweep-time must be >= 0.')
    elif engine == 'msms':
        if any(getattr(args, n) is not None for n in ('fix_time','a','x')):
            if args.a is None:
                sys.exit('ERROR: msms sweep requires --sel-2Ns (and optionally --sweep-pos).')
            if args.x is None:
                args.x = 0.5
            if args.fix_time is not None and args.fix_time < 0:
                sys.exit('ERROR: --fixation-time must be >= 0 (time back in 4N units).')

    # Adaptive target-S segregating sites refinement if user gave --target-snps without explicit --length
    def _pilot_build(engine_name, length_val, reps_override=1, min_snps_forward=None):
        if engine_name == 'discoal':
            cmd, meta_local = build_discoal_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                discoal_pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, sweep_time=getattr(args,'sweep_time',None),
                x=getattr(args,'x',None), a=getattr(args,'a',None), min_snps=min_snps_forward, chr_name=args.chromosome,
                disable_en_ladder=args.no_en_ladder, debug=False, seed=None,
            )
        elif engine_name == 'msms':
            cmd, meta_local = build_msms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=min_snps_forward, disable_en_ladder=args.no_en_ladder,
                a=getattr(args,'a',None), x=getattr(args,'x',None), fixation_time=getattr(args,'fix_time',0.0), debug=False, seed=None,
            )
        else:  # ms
            cmd, meta_local = build_ms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=reps_override, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=min_snps_forward, disable_en_ladder=args.no_en_ladder, debug=False, seed=None,
            )
        return cmd, meta_local

    if args.length is None and args.min_snps is not None:
        # Initial builder-derived length (using existing formula logic in builders)
        init_cmd, init_meta = _pilot_build(engine, None, reps_override=1, min_snps_forward=args.min_snps)
        candidate_len = init_meta['length']
        target_snps = args.min_snps
        max_iter = 6
        last_obs = None
        for it in range(max_iter):
            pilot_cmd, _m = _pilot_build(engine, candidate_len, reps_override=1, min_snps_forward=None)
            raw_line = None
            for l in pilot_cmd.splitlines():
                if l.startswith(engine+' '):
                    raw_line = l; break
            if raw_line is None:
                break
            toks = shlex.split(raw_line)
            try:
                r = subprocess.run(toks, capture_output=True, text=True, check=True)
                out_txt = r.stdout
            except Exception:
                break
            m = re.search(r'^segsites:\s*(\d+)', out_txt, re.MULTILINE)
            segs = int(m.group(1)) if m else 0
            last_obs = segs
            if segs >= target_snps:
                break
            if segs == 0:
                candidate_len = int(candidate_len * 2)
            else:
                factor = target_snps / float(segs)
                inflation = 1.10 if engine == 'ms' else (1.35 if engine == 'msms' else 1.50)
                candidate_len = int(math.ceil(candidate_len * factor * inflation))
        # Optional shrink if overshoot beyond internal tolerance
        tol = getattr(args,'min_snps_tol', 0.0)
        if tol > 0 and last_obs is not None and last_obs > target_snps * (1.0 + tol):
            shrink_iters = 3
            for _ in range(shrink_iters):
                if last_obs <= target_snps * (1.0 + tol):
                    break
                ratio = target_snps / float(last_obs)
                candidate_len = max(100, int(candidate_len * ratio * 0.98))
                pilot_cmd, _m = _pilot_build(engine, candidate_len, reps_override=1, min_snps_forward=None)
                raw_line = None
                for l in pilot_cmd.splitlines():
                    if l.startswith(engine+' '):
                        raw_line = l; break
                if raw_line is None:
                    break
                toks = shlex.split(raw_line)
                try:
                    r = subprocess.run(toks, capture_output=True, text=True, check=True)
                    out_txt = r.stdout
                except Exception:
                    break
                m = re.search(r'^segsites:\s*(\d+)', out_txt, re.MULTILINE)
                last_obs = int(m.group(1)) if m else 0
                if last_obs < target_snps:
                    candidate_len = int(math.ceil(candidate_len * (target_snps / max(1,last_obs)) * 1.05))
                    break
        refined_len = candidate_len
        if refined_len != init_meta['length'] and args.debug:
            sys.stderr.write(f"# DEBUG: auto-adjust length {init_meta['length']} -> {refined_len} (target_snps={target_snps}, last_observed={last_obs})\n")
        args.length = refined_len

    # build command + run (final, using possibly refined args.length)
    if engine == 'discoal':
        sim_cmd, meta = build_discoal_command(
            species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
            discoal_pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
            max_fold_per_step=args.max_fold_per_step, sweep_time=getattr(args,'sweep_time',None),
            x=getattr(args,'x',None), a=getattr(args,'a',None), min_snps=args.min_snps, chr_name=args.chromosome,
            disable_en_ladder=args.no_en_ladder, debug=args.debug, seed=None,
        )
    elif engine == 'msms':
        sim_cmd, meta = build_msms_command(
            species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
            pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
            max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
            min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder,
            a=getattr(args,'a',None), x=getattr(args,'x',None), fixation_time=getattr(args,'fix_time',0.0), debug=args.debug, seed=None,
        )
    else:
        # ms neutral only: selection flags were already warned & ignored in parse_args.
        sim_cmd, meta = build_ms_command(
            species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
            # For engine=ms there is no --sweep-pop/--sweep-pos group; pop0 falls back to first deme.
            pop0=getattr(args, 'discoal_pop0', None), reps=args.reps, length=args.length,
            max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
            min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, debug=args.debug, seed=None,
        )

    # ---------------- Simulation runner + enforcement helpers ----------------
    def _simulate_current_length(length_val, show_progress=True):
        """Build and execute simulation at given length. Returns meta, raw_cmd_line, stdout, stderr."""
        # rebuild command for current length (meta updated)
        if engine == 'discoal':
            sim_cmd_loc, meta_loc = build_discoal_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                discoal_pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, sweep_time=getattr(args,'sweep_time',None),
                x=getattr(args,'x',None), a=getattr(args,'a',None), min_snps=args.min_snps, chr_name=args.chromosome,
                disable_en_ladder=args.no_en_ladder, debug=args.debug, seed=None,
            )
        elif engine == 'msms':
            sim_cmd_loc, meta_loc = build_msms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder,
                a=getattr(args,'a',None), x=getattr(args,'x',None), fixation_time=getattr(args,'fix_time',0.0), debug=args.debug, seed=None,
            )
        else:
            sim_cmd_loc, meta_loc = build_ms_command(
                species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=length_val,
                max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                min_snps=args.min_snps, disable_en_ladder=args.no_en_ladder, debug=args.debug, seed=None,
            )
        raw_cmd_line_loc = None
        for ln in sim_cmd_loc.splitlines():
            if ln.startswith(engine+' '):
                raw_cmd_line_loc = ln.strip(); break
        if raw_cmd_line_loc is None:
            sys.exit('ERROR: Failed to locate simulator command line (builder).')
        cmd_tokens_loc = shlex.split(raw_cmd_line_loc)
        total_reps_loc = int(cmd_tokens_loc[2])
        # progress handling (optionally suppressed for enforcement iterations)
        threads_loc = max(1, min(int(getattr(args,'threads',1)), total_reps_loc))
        per_rep_mode_loc = bool(args.progress and show_progress and total_reps_loc > 1)
        if per_rep_mode_loc:
            chunks_loc = [1]*total_reps_loc
        else:
            if threads_loc == 1 or total_reps_loc == 1:
                chunks_loc = [total_reps_loc]
            else:
                base = total_reps_loc // threads_loc; rem = total_reps_loc % threads_loc
                chunks_loc = [base + (1 if i<rem else 0) for i in range(threads_loc) if base + (1 if i<rem else 0) > 0]
        # Prepare unique seeds per chunk for engines supporting explicit seeding
        chunk_process_count = None
        per_rep_engine = engine in ('discoal','msms') and per_rep_mode_loc
        # Consolidate ms into single process to avoid duplicate seeds while honoring user request not to use -seeds
        if engine == 'ms' and threads_loc > 1 and not per_rep_mode_loc and len(chunks_loc) > 1:
            chunks_loc = [total_reps_loc]
        if engine in ('discoal','msms'):
            if per_rep_mode_loc:
                chunk_process_count = len(chunks_loc)  # each chunk = 1 replicate
            else:
                chunk_process_count = len(chunks_loc)
            try:
                unique_seeds = set()
                seeds_for_chunks = []
                import random as _rnd2
                while len(seeds_for_chunks) < chunk_process_count:
                    s = _rnd2.randint(1, 2**31 - 1)
                    if s in unique_seeds: continue
                    unique_seeds.add(s); seeds_for_chunks.append(s)
            except Exception:
                seeds_for_chunks = [None]*chunk_process_count
        else:
            seeds_for_chunks = [None]*len(chunks_loc)

        def run_chunk_loc(rc, replicate_index=None, chunk_idx=None):
            toks = cmd_tokens_loc[:]; toks[2] = str(rc)
            if engine in ('discoal','msms'):
                seed_val = seeds_for_chunks[chunk_idx] if chunk_idx is not None and chunk_idx < len(seeds_for_chunks) else None
                if seed_val is not None:
                    # append -seed unless already present
                    if '-seed' not in toks:
                        toks.extend(['-seed', str(seed_val)])
            try:
                r = subprocess.run(toks, capture_output=True, text=True, check=True)
                return (rc, r.stdout, r.stderr, None, replicate_index)
            except subprocess.CalledProcessError as e:
                return (rc, e.stdout, e.stderr, e.returncode, replicate_index)
        results_loc = []
        bar_loc = None
        progress_done_loc = 0
        if args.progress and show_progress:
            if tqdm is not None:
                bar_loc = tqdm(total=total_reps_loc, desc='replicates', unit='rep')
            else:
                if total_reps_loc == 1:
                    sys.stderr.write('# progress: 0/1\n')
                else:
                    sys.stderr.write('# WARNING: --progress requested but tqdm not installed; using textual updates.\n')
        if per_rep_mode_loc:
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads_loc) as ex:
                futs = []
                rep_counter = 0
                for idx,c in enumerate(chunks_loc):
                    futs.append(ex.submit(run_chunk_loc, c, rep_counter, idx))
                    rep_counter += 1
                for fut in concurrent.futures.as_completed(futs):
                    rc,out,err,code,rep_idx = fut.result()
                    results_loc.append((rc,out,err,code,rep_idx))
                    if bar_loc is not None:
                        bar_loc.update(1)
                    elif args.progress and show_progress:
                        progress_done_loc += 1
                        sys.stderr.write(f'# progress: {progress_done_loc}/{total_reps_loc}\n')
            results_loc.sort(key=lambda x: (x[4] if len(x)>4 else 0))
        else:
            if len(chunks_loc) == 1:
                rc,out,err,code,rep_idx = run_chunk_loc(chunks_loc[0], 0, 0)
                results_loc.append((rc,out,err,code,rep_idx))
                if bar_loc is not None:
                    bar_loc.update(total_reps_loc)
                elif args.progress and show_progress and total_reps_loc == 1:
                    sys.stderr.write('# progress: 1/1\n')
            else:
                with concurrent.futures.ThreadPoolExecutor(max_workers=len(chunks_loc)) as ex:
                    futs = [ex.submit(run_chunk_loc, c, i, i) for i,c in enumerate(chunks_loc)]
                    for fut in concurrent.futures.as_completed(futs):
                        rc,out,err,code,rep_idx = fut.result()
                        results_loc.append((rc,out,err,code,rep_idx))
                        if bar_loc is not None:
                            try:
                                bar_loc.update(int(rc))
                            except Exception:
                                bar_loc.update(1)
                        elif args.progress and show_progress:
                            progress_done_loc += rc
                            sys.stderr.write(f'# progress: {progress_done_loc}/{total_reps_loc}\n')
        if bar_loc is not None:
            bar_loc.close()
        for tpl in results_loc:
            if len(tpl) == 5:
                rc,out,err,code,_r = tpl
            else:
                rc,out,err,code = tpl
            if code is not None:
                sys.stderr.write(f'# ERROR: {engine} exited with code {code}.\n')
                if out: sys.stderr.write(out)
                if err: sys.stderr.write(err)
                sys.exit(code)
        concatenated_stdout_loc = ''.join(t[1] for t in results_loc)
        concatenated_stderr_loc = ''.join(t[2] for t in results_loc)
        return meta_loc, raw_cmd_line_loc, concatenated_stdout_loc, concatenated_stderr_loc

    # Initial simulation at refined/derived length
    meta, raw_cmd_line, concatenated_stdout, concatenated_stderr = _simulate_current_length(args.length, show_progress=True)

    # Enforcement: ensure ALL replicates meet/exceed target and optionally shrink overshoot
    if args.min_snps is not None:
        target_snps = args.min_snps
        tol_internal = getattr(args,'min_snps_tol', 0.0)
        max_enforce_iters = 5
        enforce_iter = 0  # count of inflation iterations
        length_current = args.length
        pattern_segs = re.compile(r'^segsites:\s*(\d+)', re.MULTILINE)
        while enforce_iter < max_enforce_iters:
            seg_list = [int(x) for x in pattern_segs.findall(concatenated_stdout)]
            if not seg_list:
                break
            min_segs = min(seg_list)
            if min_segs < target_snps:
                factor = target_snps / max(1, min_segs)
                inflation = 1.05 if engine == 'ms' else (1.15 if engine == 'msms' else 1.20)
                new_length = int(math.ceil(length_current * factor * inflation))
                if new_length <= length_current:
                    new_length = length_current + 1
                if args.debug:
                    sys.stderr.write(f"# DEBUG: enforce iter{enforce_iter} length {length_current} -> {new_length} (min_segs={min_segs} target={target_snps})\n")
                length_current = new_length
                meta, raw_cmd_line, concatenated_stdout, concatenated_stderr = _simulate_current_length(length_current, show_progress=False)
                enforce_iter += 1
                continue
            break
        args.length = length_current
        if args.debug:
            seg_list_final = [int(x) for x in pattern_segs.findall(concatenated_stdout)]
            if seg_list_final:
                sys.stderr.write(f"# DEBUG: enforcement complete length={args.length} segsites_min={min(seg_list_final)} segsites_max={max(seg_list_final)} target={target_snps} tol={tol_internal}\n")

        overshoot_status = None
        shrink_bracket_attempts = 0
        shrink_extra_reductions = 0
        shrink_bs_iters = 0
        if tol_internal > 0:
            seg_list_curr = [int(x) for x in pattern_segs.findall(concatenated_stdout)]
            if seg_list_curr:
                hi_len = args.length
                hi_stdout = concatenated_stdout
                hi_meta = meta
                hi_seg_list = seg_list_curr
                lo_len = None
                lo_seg_list = None
                attempts = 0
                while attempts < 6:
                    max_segs_curr = max(hi_seg_list)
                    target_cap = target_snps * (1.0 + tol_internal)
                    scale_cap = target_cap / float(max_segs_curr)
                    trial_len = int(math.ceil(hi_len * max(0.60, min(0.95, scale_cap))))
                    if trial_len >= hi_len or trial_len < 100:
                        break
                    meta_trial, raw_cmd_line_trial, stdout_trial, stderr_trial = _simulate_current_length(trial_len, show_progress=False)
                    seg_list_trial = [int(x) for x in pattern_segs.findall(stdout_trial)]
                    if not seg_list_trial:
                        break
                    if min(seg_list_trial) >= target_snps:
                        hi_len, hi_stdout, hi_meta, hi_seg_list = trial_len, stdout_trial, meta_trial, seg_list_trial
                        attempts += 1
                        shrink_bracket_attempts = attempts
                        continue
                    else:
                        lo_len = trial_len
                        lo_seg_list = seg_list_trial
                        break
                if lo_len is None:
                    probe_len = int(hi_len * 0.85)
                    extra = 0
                    while extra < 6 and probe_len >= 100:
                        meta_trial, raw_cmd_line_trial, stdout_trial, stderr_trial = _simulate_current_length(probe_len, show_progress=False)
                        seg_list_trial = [int(x) for x in pattern_segs.findall(stdout_trial)]
                        if not seg_list_trial:
                            break
                        if min(seg_list_trial) >= target_snps:
                            hi_len, hi_stdout, hi_meta, hi_seg_list = probe_len, stdout_trial, meta_trial, seg_list_trial
                            probe_len = int(probe_len * 0.85)
                            shrink_extra_reductions = extra + 1
                        else:
                            lo_len = probe_len
                            lo_seg_list = seg_list_trial
                            break
                        extra += 1
                if lo_len is not None:
                    bs_iter = 0
                    while bs_iter < 12 and lo_len + 1 < hi_len:
                        mid = (lo_len + hi_len) // 2
                        if mid == hi_len:
                            break
                        meta_mid, raw_cmd_line_mid, stdout_mid, stderr_mid = _simulate_current_length(mid, show_progress=False)
                        seg_list_mid = [int(x) for x in pattern_segs.findall(stdout_mid)]
                        if not seg_list_mid:
                            break
                        if min(seg_list_mid) >= target_snps:
                            hi_len, hi_stdout, hi_meta, hi_seg_list = mid, stdout_mid, meta_mid, seg_list_mid
                        else:
                            lo_len = mid
                        bs_iter += 1
                        shrink_bs_iters = bs_iter
                args.length = hi_len
                concatenated_stdout = hi_stdout
                meta = hi_meta
                max_final = max(hi_seg_list)
                if max_final <= target_snps * (1.0 + tol_internal):
                    overshoot_status = 'met'
                else:
                    overshoot_status = 'unmet_minimal_feasible'
                    if args.debug:
                        sys.stderr.write(f"# DEBUG: tolerance unmet at minimal feasible length (len={hi_len}) max_segs={max_final} threshold={target_snps*(1.0+tol_internal):.2f}\n")
                if args.debug:
                    min_final = min(hi_seg_list)
                    sys.stderr.write(f"# DEBUG: shrink-binary finalized length={hi_len} segsites_min={min_final} segsites_max={max_final} status={overshoot_status} tol={tol_internal}\n")
            else:
                overshoot_status = 'no_data'
        final_seg_list = [int(x) for x in pattern_segs.findall(concatenated_stdout)] if 'pattern_segs' in locals() else []
        if final_seg_list:
            meta.setdefault('enforcement_stats', {})
            meta['enforcement_stats'].update({
                'target_snps': target_snps,
                'tol': tol_internal,
                'segsites_min': min(final_seg_list),
                'segsites_max': max(final_seg_list),
                'overshoot': overshoot_status or ('n/a' if tol_internal == 0 else 'unknown'),
                'enforce_iters': enforce_iter,
                'shrink_bracket_attempts': shrink_bracket_attempts,
                'shrink_extra_reductions': shrink_extra_reductions,
                'shrink_bs_iters': shrink_bs_iters,
            })
    # Decide default naming pieces
    model_id = getattr(args, 'model', 'model')
    pop0 = None
    if engine == 'discoal':
        pop0 = getattr(args, 'discoal_pop0', None)
    elif engine == 'msms':
        pop0 = getattr(args, 'discoal_pop0', None)
    else:  # ms
        pop0 = getattr(args, 'discoal_pop0', None)
    if not pop0:
        # fallback to first deme in user_order
        pop0 = user_order[0] if user_order else 'pop0'
    chrom = getattr(args, 'chromosome', None)
    sweeppop = None
    if engine in ('discoal','msms'):
        sweeppop = getattr(args, 'discoal_pop0', None)
    # Base stem
    parts_name = [model_id, pop0]
    if chrom:
        parts_name.append(str(chrom))
    if sweeppop and sweeppop != pop0:
        parts_name.append(sweeppop)
    default_stem = '_'.join(parts_name)

    if args.sfs:
        # Use accelerated implementation if NumPy is available; falls back automatically.
        header_sfs, lines_sfs = compute_sfs_fast(concatenated_stdout, normalized=args.sfs_normalized, mode=args.sfs_mode)
        sfs_text = header_sfs + '\n' + '\n'.join(lines_sfs) + '\n'
        # Determine SFS file path
        if isinstance(args.sfs, str) and args.sfs is not True and args.sfs != 'True':
            # User provided custom base/path; honor it (append .sfs if missing)
            if args.sfs.endswith('.sfs'):
                sfs_path = args.sfs
            else:
                sfs_path = args.sfs + '.sfs'
        else:
            # No custom SFS name: if --output supplied (and not '-') derive from it; else use default stem
            if getattr(args, 'out', None) and args.out not in (None, '-'):
                # Strip known output extensions then add .sfs
                base = args.out
                for ext_candidate in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                    if base.endswith(ext_candidate):
                        base = base[:-len(ext_candidate)]
                        break
                sfs_path = base + '.sfs'
            else:
                sfs_path = default_stem + '.sfs'
        try:
            with open(sfs_path, 'w') as f:
                f.write('# SFS output\n')
                f.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                f.write(sfs_text)
            sys.stderr.write(f"# INFO: wrote SFS to {sfs_path}\n")
        except Exception as e:
            sys.stderr.write(f"# ERROR: failed to write SFS file {sfs_path}: {e}\n")
            if not args.out or args.out == '-':
                sys.stdout.write(sfs_text)
    # Handle primary output: if args.out is None -> write file with default stem; if '-' -> stdout
    if args.out is None:
        # build file name according to format
        fmt_ext = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
        ext = fmt_ext.get(args.format, '.ms')
        args.out = default_stem + ext
    if args.out == '-':
        if concatenated_stdout:
            sys.stdout.write(concatenated_stdout)
    elif not args.sfs and not concatenated_stdout and args.format in ('ms','ms.gz'):
        # nothing to write
        pass
    if concatenated_stderr:
        sys.stderr.write(concatenated_stderr)
    if args.out and args.out != '-':
        # Normalize output filename extension based on requested format if user omitted one.
        ext_map = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
        desired_ext = ext_map.get(args.format, '')
        if desired_ext and not args.out.endswith(desired_ext):
            # Strip any recognized extension (handle multi-part .vcf.gz first) then append desired
            for kext in ('.vcf.gz', '.vcf', '.bcf', '.ms.gz', '.ms'):
                if args.out.endswith(kext):
                    args.out = args.out[:-len(kext)]
                    break
            args.out = args.out + desired_ext
        fmt = args.format
        # Write primary output first
        if fmt.startswith('vcf'):
            vcf_text = ms_like_to_vcf(concatenated_stdout, meta['length'], chrom=meta.get('chromosome') or 'chr1', ploidy=meta.get('ploidy',1))
            # Build metadata lines identical to ms output style
            info_pairs = []
            info_pairs.append(('engine', meta.get('engine', engine)))
            info_pairs.append(('species', meta.get('species_id')))
            info_pairs.append(('model', model_id))
            info_pairs.append(('pop0', meta.get('pop0')))
            if meta.get('pop_order'):
                info_pairs.append(('order', '|'.join(map(str, meta.get('pop_order')))))
            if meta.get('counts_disc'):
                info_pairs.append(('sample_counts', '|'.join(map(str, meta.get('counts_disc')))))
            theta_val = rho_val = None
            try:
                toks_tmp = raw_cmd_line.split()
                if '-t' in toks_tmp:
                    theta_val = toks_tmp[toks_tmp.index('-t')+1]
                if '-r' in toks_tmp:
                    rho_val = toks_tmp[toks_tmp.index('-r')+1]
            except Exception:
                pass
            if theta_val: info_pairs.append(('theta', theta_val))
            if rho_val: info_pairs.append(('rho', rho_val))
            info_pairs.append(('length', meta.get('length')))
            info_pairs.append(('N0', meta.get('N0')))
            info_pairs.append(('mu', meta.get('mu')))
            info_pairs.append(('r', meta.get('rrate')))
            info_pairs.append(('ploidy', meta.get('ploidy')))
            if meta.get('chromosome') is not None:
                info_pairs.append(('chromosome', meta.get('chromosome')))
            if meta.get('sweep_pos') is not None:
                sp = meta.get('sweep_pos')
                info_pairs.append(('sweep_pos', sp))
                try:
                    if sp is not None and meta.get('length'):
                        sweep_bp = int(round(float(sp) * int(meta.get('length'))))
                        info_pairs.append(('sweep_bp', sweep_bp))
                except Exception:
                    pass
            if meta.get('sel_2Ns') is not None:
                info_pairs.append(('sel_2Ns', meta.get('sel_2Ns')))
            if meta.get('sweep_time') is not None:
                info_pairs.append(('sweep_time', meta.get('sweep_time')))
            if meta.get('fixation_time') is not None:
                info_pairs.append(('fixation_time', meta.get('fixation_time')))
            kv_dict = {k:v for k,v in info_pairs if v is not None}
            enf_stats = meta.get('enforcement_stats') or {}
            for k,v in enf_stats.items():
                kv_dict[k] = v
            core_keys = ['engine','species','model','pop0','order','sample_counts']
            param_keys = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
            sel_keys = ['sweep_pos','sweep_bp','sel_2Ns','sweep_time','fixation_time']
            def build_group_line(prefix, keys):
                vals = [f"{k}={kv_dict[k]}" for k in keys if k in kv_dict]
                if not vals:
                    return None
                return '# ' + prefix + ': ' + ', '.join(vals)
            meta_comment_lines = []
            for grp, keys in (('core', core_keys), ('params', param_keys), ('selection', sel_keys)):
                ln = build_group_line(grp, keys)
                if ln:
                    meta_comment_lines.append(ln)
            # Inject these lines into VCF header (convert '# ' to '##' for VCF compliance while preserving human readability)
            vcf_lines = vcf_text.splitlines()
            try:
                header_idx = next(i for i,l in enumerate(vcf_lines) if l.startswith('#CHROM'))
            except StopIteration:
                header_idx = len(vcf_lines)
            vcf_meta_lines = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_comment_lines]
            vcf_lines = vcf_lines[:header_idx] + vcf_meta_lines + vcf_lines[header_idx:]
            vcf_text_final = '\n'.join(vcf_lines) + ('\n' if not vcf_lines[-1].endswith('\n') else '')
            if fmt == 'vcf.gz':
                with gzip.open(args.out, 'wt') as f:
                    f.write(vcf_text_final)
            else:
                with open(args.out, 'w') as f:
                    f.write(vcf_text_final)
        else:
            # Raw ms-like output (optionally gzipped)
            f_open = gzip.open if fmt == 'ms.gz' else open
            mode = 'wt' if fmt == 'ms.gz' else 'w'
            with f_open(args.out, mode) as f:
                f.write(raw_cmd_line + '\n')
                # Build unified comma-separated info line
                info_pairs = []
                info_pairs.append(('engine', meta.get('engine', engine)))
                info_pairs.append(('species', meta.get('species_id')))
                info_pairs.append(('model', model_id))
                info_pairs.append(('pop0', meta.get('pop0')))
                if meta.get('pop_order'):
                    info_pairs.append(('order', '|'.join(map(str, meta.get('pop_order')))))
                if meta.get('counts_disc'):
                    info_pairs.append(('sample_counts', '|'.join(map(str, meta.get('counts_disc')))))
                # parse theta/rho from command
                theta_val = rho_val = None
                try:
                    toks_tmp = raw_cmd_line.split()
                    if '-t' in toks_tmp:
                        theta_val = toks_tmp[toks_tmp.index('-t')+1]
                    if '-r' in toks_tmp:
                        rho_val = toks_tmp[toks_tmp.index('-r')+1]
                except Exception:
                    pass
                if theta_val: info_pairs.append(('theta', theta_val))
                if rho_val: info_pairs.append(('rho', rho_val))
                info_pairs.append(('length', meta.get('length')))
                info_pairs.append(('N0', meta.get('N0')))
                info_pairs.append(('mu', meta.get('mu')))
                info_pairs.append(('r', meta.get('rrate')))
                info_pairs.append(('ploidy', meta.get('ploidy')))
                if meta.get('chromosome') is not None:
                    info_pairs.append(('chromosome', meta.get('chromosome')))
                if meta.get('sweep_pos') is not None:
                    sp = meta.get('sweep_pos')
                    info_pairs.append(('sweep_pos', sp))
                    try:
                        if sp is not None and meta.get('length'):
                            sweep_bp = int(round(float(sp) * int(meta.get('length'))))
                            info_pairs.append(('sweep_bp', sweep_bp))
                    except Exception: pass
                if meta.get('sel_2Ns') is not None:
                    info_pairs.append(('sel_2Ns', meta.get('sel_2Ns')))
                if meta.get('sweep_time') is not None:
                    info_pairs.append(('sweep_time', meta.get('sweep_time')))
                if meta.get('fixation_time') is not None:
                    info_pairs.append(('fixation_time', meta.get('fixation_time')))
                # demog intentionally omitted from metadata output per user request
                # Group related keys: core, params, selection (enforcement line removed per user request)
                core_keys = ['engine','species','model','pop0','order','sample_counts']
                param_keys = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                sel_keys = ['sweep_pos','sweep_bp','sel_2Ns','sweep_time','fixation_time']
                # Enforcement diagnostics retained internally but not emitted.
                kv_dict = {k:v for k,v in info_pairs if v is not None}
                # Inject enforcement stats if present
                enf = meta.get('enforcement_stats') or {}
                for k,v in enf.items():
                    kv_dict[k] = v
                def build_line(group_name, keys):
                    vals = [f"{k}={kv_dict[k]}" for k in keys if k in kv_dict]
                    if not vals:
                        return None
                    return '# ' + group_name + ': ' + ', '.join(vals)
                lines_info = []
                for grp_name, keys in (
                    ('core', core_keys),
                    ('params', param_keys),
                    ('selection', sel_keys),
                ):
                    ln = build_line(grp_name, keys)
                    if ln:
                        lines_info.append(ln)
                for ln in lines_info:
                    f.write(ln + '\n')
                # format line removed per user request
                lines = concatenated_stdout.splitlines(True)
                for ln in lines:
                    # Avoid duplicating echoed command lines. For ms and discoal we skip the engine line.
                    # For msms, some versions may echo both an 'msms ...' line and a compatibility 'ms ...' line;
                    # keep only the 'msms' line we already wrote (raw_cmd_line) and drop any stray 'ms ' line.
                    if ln.startswith(engine+' '):
                        continue
                    if engine == 'msms' and ln.startswith('ms '):
                        continue
                    f.write(ln)
    # Optional paired neutral run
    if getattr(args, 'paired_neutral', False):
        fmt = args.format  # ensure available for neutral output generation
        neut_engine = args.neutral_engine or engine
        if getattr(args, 'neutral_same_seed', False) and (neut_engine == engine):
            neut_seed = primary_seed
        else:
            try:
                import random as _rnd
                neut_seed = _rnd.randint(1, 2**31 - 1)
                if primary_seed is not None and neut_seed == primary_seed:
                    neut_seed = (neut_seed + 1) % (2**31 - 1) or 1
            except Exception:
                neut_seed = None
        # If user wants exact neutral twin and engine matches, derive neutral command by stripping selection flags
        neut_cmd = None; neut_meta = None
        if neut_engine == engine and engine in ('discoal','msms'):
            # Start from original raw command line (post-enforcement) to guarantee identical theta/rho/length/demog/sample order
            base_line = raw_cmd_line
            toks = shlex.split(base_line)
            def _strip_selection_tokens_discoal(tokens):
                out = []
                skip_next = 0
                sel_flags = {'-ws':1,'-x':1,'-a':1}
                i=0
                while i < len(tokens):
                    t = tokens[i]
                    if t in sel_flags:
                        i += 1 + sel_flags[t]  # skip flag + its args
                        continue
                    out.append(t); i += 1
                return out
            def _strip_selection_tokens_msms(tokens):
                out = []
                i = 0
                while i < len(tokens):
                    t = tokens[i]
                    if t in ('-SaA','-SAA','-Sp') and i+1 < len(tokens):
                        i += 2; continue
                    if t == '-Smark':
                        i += 1; continue
                    if t == '-SI':
                        # pattern: -SI <tstart> <npop> <freq1> ... <freqN>
                        if i+2 < len(tokens):
                            try:
                                npop_local = int(tokens[i+2])
                                i += 3 + npop_local
                                continue
                            except Exception:
                                # fallback: skip just flag
                                i += 1
                                continue
                        else:
                            i += 1
                            continue
                    if t == '-SFC':
                        i += 1; continue
                    out.append(t); i += 1
                return out
            if engine == 'discoal':
                stripped = _strip_selection_tokens_discoal(toks)
            else:  # msms
                stripped = _strip_selection_tokens_msms(toks)
            # Append seed if requested and not already there
            if neut_seed is not None and '-seed' not in stripped:
                stripped += ['-seed', str(neut_seed)]
            neut_line = ' '.join(stripped)
            # Reconstruct full command block (reuse comments but remove sweep-specific comment lines)
            base_comments = [c for c in meta.get('base_comments', []) if 'sweep-pos=' not in c and 'sweep-time=' not in c and 'fixation-time=' not in c]
            neut_cmd = '\n'.join(base_comments + [neut_line])
            # Clone meta and blank selection fields
            neut_meta = dict(meta)
            for k in ('sweep_pos','sel_2Ns','sweep_time','fixation_time'):
                neut_meta[k] = None
            neut_meta['base_comments'] = base_comments
            neut_meta['engine'] = neut_engine
        else:
            # Fallback: build via neutral engine-specific builder (may differ if user intentionally changed engine)
            if neut_engine == 'discoal':
                neut_cmd, neut_meta = build_discoal_command(
                    species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    discoal_pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
                    max_fold_per_step=args.max_fold_per_step, sweep_time=None, x=None, a=None, min_snps=None, chr_name=args.chromosome,
                    disable_en_ladder=args.no_en_ladder, debug=args.debug, seed=neut_seed,
                )
            elif neut_engine == 'msms':
                neut_cmd, neut_meta = build_msms_command(
                    species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
                    max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                    min_snps=None, disable_en_ladder=args.no_en_ladder,
                    a=None, x=None, fixation_time=0.0, debug=args.debug, seed=neut_seed,
                )
            else:
                neut_cmd, neut_meta = build_ms_command(
                    species=args.species, model_id=args.model, user_order=user_order, individual_counts=individual_counts,
                    pop0=getattr(args,'discoal_pop0',None), reps=args.reps, length=args.length,
                    max_fold_per_step=args.max_fold_per_step, chr_name=args.chromosome,
                    min_snps=None, disable_en_ladder=args.no_en_ladder, debug=args.debug, seed=neut_seed,
                )
        raw_neut = None
        neut_err_txt = ''
        for ln in neut_cmd.splitlines():
            if ln.startswith(neut_engine+' '):
                raw_neut = ln.strip(); break
        if raw_neut is None:
            sys.stderr.write('# ERROR: failed to extract neutral command line; skipping neutral run.\n')
        else:
            try:
                rN = subprocess.run(shlex.split(raw_neut), capture_output=True, text=True, check=True)
                neut_out_txt = rN.stdout; neut_err_txt = rN.stderr
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"# ERROR: neutral {neut_engine} exited code {e.returncode}.\n")
                if e.stdout: sys.stderr.write(e.stdout)
                if e.stderr: sys.stderr.write(e.stderr)
                neut_out_txt = ''; neut_err_txt = ''
            neut_out_path = getattr(args,'neutral_out', None)
            if neut_out_path is None and args.out and args.out != '-':
                main_out = args.out
                if main_out.endswith('.vcf.gz'):
                    root = main_out[:-7]; ext = '.vcf.gz'
                else:
                    root, ext = os.path.splitext(main_out)
                neut_out_path = root + '_neutral' + ext
            # Force neutral output format to mirror primary format even if user supplied a mismatched extension
            if neut_out_path and neut_out_path != '-':
                ext_map_full = {'ms':'.ms','ms.gz':'.ms.gz','vcf':'.vcf','vcf.gz':'.vcf.gz','bcf':'.bcf'}
                desired_ext = ext_map_full.get(fmt,'')
                if desired_ext:
                    if neut_out_path.endswith('.vcf.gz') and desired_ext == '.vcf.gz':
                        pass  # already correct
                    elif neut_out_path.endswith('.ms.gz') and desired_ext == '.ms.gz':
                        pass  # already correct
                    elif not neut_out_path.endswith(desired_ext):
                        # Strip any existing known extension then append desired
                        for kext in ('.vcf.gz','.vcf','.bcf','.ms.gz','.ms'):
                            if neut_out_path.endswith(kext):
                                neut_out_path = neut_out_path[:-len(kext)]
                                break
                        neut_out_path = neut_out_path + desired_ext
            if neut_out_path and neut_out_path != '-':
                try:
                    if fmt.startswith('vcf'):
                        vcf_text_neut = ms_like_to_vcf(neut_out_txt, neut_meta['length'], chrom=neut_meta.get('chromosome') or 'chr1', ploidy=neut_meta.get('ploidy',1))
                        # Build neutral metadata lines analogous to ms output
                        info_pairs_n = []
                        info_pairs_n.append(('engine', neut_meta.get('engine', neut_engine)))
                        info_pairs_n.append(('species', neut_meta.get('species_id')))
                        info_pairs_n.append(('model', args.model))
                        info_pairs_n.append(('pop0', neut_meta.get('pop0')))
                        pop_order_n = neut_meta.get('pop_order') or []
                        if pop_order_n:
                            info_pairs_n.append(('order', '|'.join(map(str, pop_order_n))))
                        counts_disc_n = neut_meta.get('counts_disc') or []
                        if counts_disc_n:
                            info_pairs_n.append(('sample_counts', '|'.join(map(str, counts_disc_n))))
                        theta_val_n = rho_val_n = None
                        try:
                            toks_nt = raw_neut.split()
                            if '-t' in toks_nt:
                                theta_val_n = toks_nt[toks_nt.index('-t')+1]
                            if '-r' in toks_nt:
                                rho_val_n = toks_nt[toks_nt.index('-r')+1]
                        except Exception:
                            pass
                        if theta_val_n: info_pairs_n.append(('theta', theta_val_n))
                        if rho_val_n: info_pairs_n.append(('rho', rho_val_n))
                        info_pairs_n.append(('length', neut_meta.get('length')))
                        info_pairs_n.append(('N0', neut_meta.get('N0')))
                        info_pairs_n.append(('mu', neut_meta.get('mu')))
                        info_pairs_n.append(('r', neut_meta.get('rrate')))
                        info_pairs_n.append(('ploidy', neut_meta.get('ploidy')))
                        if neut_meta.get('chromosome') is not None:
                            info_pairs_n.append(('chromosome', neut_meta.get('chromosome')))
                        kvn = {k:v for k,v in info_pairs_n if v is not None}
                        core_keys_n = ['engine','species','model','pop0','order','sample_counts']
                        param_keys_n = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                        def build_line_n_vcf(group_name, keys):
                            vals = [f"{k}={kvn[k]}" for k in keys if k in kvn]
                            if not vals:
                                return None
                            return '# ' + group_name + ': ' + ', '.join(vals)
                        meta_lines_neut = []
                        for grp_n, keys_n in (('core', core_keys_n), ('params', param_keys_n)):
                            ln_n = build_line_n_vcf(grp_n, keys_n)
                            if ln_n:
                                meta_lines_neut.append(ln_n)
                        # Inject before #CHROM
                        vcf_lines_neut = vcf_text_neut.splitlines()
                        try:
                            header_idx_neut = next(i for i,l in enumerate(vcf_lines_neut) if l.startswith('#CHROM'))
                        except StopIteration:
                            header_idx_neut = len(vcf_lines_neut)
                        vcf_meta_lines_neut = [('##' + l[2:]) if l.startswith('# ') else ('##'+l.lstrip('# ')) for l in meta_lines_neut]
                        vcf_lines_neut = vcf_lines_neut[:header_idx_neut] + vcf_meta_lines_neut + vcf_lines_neut[header_idx_neut:]
                        vcf_text_neut_final = '\n'.join(vcf_lines_neut) + ('\n' if not vcf_lines_neut[-1].endswith('\n') else '')
                        if fmt == 'vcf.gz':
                            with gzip.open(neut_out_path, 'wt') as f2:
                                f2.write(vcf_text_neut_final)
                        else:
                            with open(neut_out_path, 'w') as f2:
                                f2.write(vcf_text_neut_final)
                    else:
                        # Raw neutral ms-like output (optionally gzipped)
                        f_open_neut = gzip.open if fmt == 'ms.gz' else open
                        mode_neut = 'wt' if fmt == 'ms.gz' else 'w'
                        with f_open_neut(neut_out_path, mode_neut) as f2:
                            f2.write(raw_neut + '\n')
                            info_pairs_n = []
                            info_pairs_n.append(('engine', neut_meta.get('engine', neut_engine)))
                            info_pairs_n.append(('species', neut_meta.get('species_id')))
                            info_pairs_n.append(('model', args.model))
                            info_pairs_n.append(('pop0', neut_meta.get('pop0')))
                            pop_order_n = neut_meta.get('pop_order') or []
                            if pop_order_n:
                                info_pairs_n.append(('order', '|'.join(map(str, pop_order_n))))
                            counts_disc_n = neut_meta.get('counts_disc') or []
                            if counts_disc_n:
                                info_pairs_n.append(('sample_counts', '|'.join(map(str, counts_disc_n))))
                            theta_val_n = rho_val_n = None
                            try:
                                toks_nt = raw_neut.split()
                                if '-t' in toks_nt:
                                    theta_val_n = toks_nt[toks_nt.index('-t')+1]
                                if '-r' in toks_nt:
                                    rho_val_n = toks_nt[toks_nt.index('-r')+1]
                            except Exception:
                                pass
                            if theta_val_n: info_pairs_n.append(('theta', theta_val_n))
                            if rho_val_n: info_pairs_n.append(('rho', rho_val_n))
                            info_pairs_n.append(('length', neut_meta.get('length')))
                            info_pairs_n.append(('N0', neut_meta.get('N0')))
                            info_pairs_n.append(('mu', neut_meta.get('mu')))
                            info_pairs_n.append(('r', neut_meta.get('rrate')))
                            info_pairs_n.append(('ploidy', neut_meta.get('ploidy')))
                            if neut_meta.get('chromosome') is not None:
                                info_pairs_n.append(('chromosome', neut_meta.get('chromosome')))
                            kvn = {k:v for k,v in info_pairs_n if v is not None}
                            core_keys_n = ['engine','species','model','pop0','order','sample_counts']
                            param_keys_n = ['theta','rho','length','N0','mu','r','ploidy','chromosome']
                            def build_line_n(group_name, keys):
                                vals = [f"{k}={kvn[k]}" for k in keys if k in kvn]
                                if not vals: return None
                                return '# ' + group_name + ': ' + ', '.join(vals)
                            for grp_n, keys_n in (('core', core_keys_n), ('params', param_keys_n)):
                                ln_n = build_line_n(grp_n, keys_n)
                                if ln_n:
                                    f2.write(ln_n + '\n')
                            for ln2 in neut_out_txt.splitlines(True):
                                if ln2.startswith(neut_engine+' '):
                                    continue
                                if neut_engine == 'msms' and ln2.startswith('ms '):
                                    continue
                                f2.write(ln2)
                    sys.stderr.write(f"# INFO: wrote paired neutral output to {neut_out_path}\n")
                    if getattr(args, 'sfs', False):
                        try:
                            header_nsfs, lines_nsfs = compute_sfs_fast(neut_out_txt, normalized=args.sfs_normalized, mode=args.sfs_mode)
                            neut_root = neut_out_path
                            if neut_root.endswith('.vcf.gz'):
                                neut_root = neut_root[:-7]
                            elif neut_root.endswith('.ms.gz'):
                                neut_root = neut_root[:-6]
                            else:
                                neut_root = os.path.splitext(neut_root)[0]
                            neut_sfs_path = neut_root + '.sfs'
                            nsfs_text = header_nsfs + '\n' + '\n'.join(lines_nsfs) + '\n'
                            with open(neut_sfs_path, 'w') as fns:
                                fns.write('# SFS output (neutral paired)\n')
                                fns.write('# normalized=' + str(args.sfs_normalized) + ' mode=' + args.sfs_mode + '\n')
                                fns.write(nsfs_text)
                            sys.stderr.write(f"# INFO: wrote neutral SFS to {neut_sfs_path}\n")
                        except Exception as e:
                            sys.stderr.write(f"# ERROR: failed to compute/write neutral SFS: {e}\n")
                except Exception as e:
                    sys.stderr.write(f"# ERROR: failed to write neutral output {neut_out_path}: {e}\n")
            if neut_err_txt:
                sys.stderr.write(neut_err_txt)
        # Write separate SFS file if requested
    # SFS already written earlier under new naming; do not duplicate

if __name__ == "__main__":
    main()