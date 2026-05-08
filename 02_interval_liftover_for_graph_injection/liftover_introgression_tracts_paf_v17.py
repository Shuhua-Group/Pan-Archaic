#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PAF-only true-subset liftover (JOIN by default, inject-safe, with detailed rejections)

- Determine the target subinterval inside the GRCh38 query window, then map it
  back to sample-path query coordinates.
- By default, emit one joined interval per path as [min_start, max_end], allowing
  internal graph bubbles or SVs.
- Endpoint snapping falls back from expand to nearest to raw true-subset window.
- Defaults are permissive but non-expanding: min-mapq=5 and 30 bp / 1% GRCh38
  boundary margin.

Example:
  python3 liftover_introgression_tracts_paf_v17.py -p graph.paf -i intro.tsv -o out.bed
"""

import argparse, sys, multiprocessing, math, bisect, json, hashlib
from collections import defaultdict
import pandas as pd
from intervaltree import IntervalTree
from tqdm import tqdm

# ---------- helpers ----------

def paf_id_to_og_path(paf_id: str) -> str:
    if paf_id.startswith('id='): paf_id = paf_id[3:]
    if '.' in paf_id and '|' in paf_id:
        try:
            sample_part, rest = paf_id.split('.', 1)
            hap_num, contig = rest.split('|', 1)
            return f"{sample_part}#{hap_num}#{contig}#0"
        except ValueError:
            pass
    if '|' in paf_id:
        try:
            ref_name, chrom_name = paf_id.split('|', 1)
            return f"{ref_name}#0#{chrom_name}"
        except ValueError:
            pass
    return f"{paf_id}#0"

def parse_minigraph_id(t_name: str) -> str:
    if t_name.startswith('id='): t_name = t_name[3:]
    try:
        return t_name.split('|', 1)[1]
    except IndexError:
        return t_name

def md5_short(s: str, k=6) -> str:
    return hashlib.md5(s.encode()).hexdigest()[:k]

def path_matches_hap(path: str, hap: str) -> bool:
    # Accept exact haplotype prefixes such as SAMPLE#1, SAMPLE#2, or SAMPLE.
    return path == hap or path.startswith(hap) or path.startswith(hap + '#')

def merge_intervals(iv_list):
    if not iv_list: return []
    iv_list = sorted((int(s), int(e)) for s,e in iv_list if e > s)
    out = [list(iv_list[0])]
    for s,e in iv_list[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s,e])
    return [(s,e) for s,e in out]

# q (GRCh38) -> t (node)
def q2t(q, g_qs, g_qe, t_s, t_e, strand):
    if g_qe == g_qs: return None
    Lq = g_qe - g_qs; Lt = t_e - t_s
    if strand == '+':
        return t_s + (q - g_qs) * Lt / Lq
    else:
        return t_e - (q - g_qs) * Lt / Lq

# t (node) -> q (sample path)
def t2q(t, s_qs, s_qe, t_s, t_e, strand):
    Lt = t_e - t_s
    if Lt == 0: return None
    Lq = s_qe - s_qs
    if strand == '+':
        return s_qs + (t - t_s) * Lq / Lt
    else:
        return s_qs + (t_e - t) * Lq / Lt

# ---------- globals ----------

GT   = None   # chr -> IntervalTree(begin=q_start,end=q_end,data=node_id) for GRCh38
GREC = None   # node_id -> list of dict( GRCh38 paf rows on this node )
NM   = None   # node_id -> list of dict( sample paf rows on this node )
PLEN = None   # path -> max q_len
BND  = None   # path -> sorted list of PAF boundaries (qs/qe)

MIN_MAPQ         = 5
SNAP_MODE        = "expand"   # expand|shrink|nearest
G38_MARGIN_BP    = 30
G38_MARGIN_FRAC  = 0.01
MAX_RATIO        = 1e9        # no ratio limit by default in JOIN mode
CLIP             = True
MIN_OUTPUT_BP    = 1
EMIT_MODE        = "join"     # join | fragments
FORCE_SNAP_ONLY  = False      # True disables the raw fallback
REJECT_VERBOSE   = True       # include detailed rejection records

def init_worker(gt, grec, nm, plen, bnd,
                min_mapq, snap_mode,
                g38_bp, g38_frac,
                max_ratio, clip, min_out_bp,
                emit_mode, force_snap_only, reject_verbose):
    global GT, GREC, NM, PLEN, BND
    global MIN_MAPQ, SNAP_MODE, G38_MARGIN_BP, G38_MARGIN_FRAC
    global MAX_RATIO, CLIP, MIN_OUTPUT_BP, EMIT_MODE, FORCE_SNAP_ONLY, REJECT_VERBOSE
    GT, GREC, NM, PLEN, BND = gt, grec, nm, plen, bnd
    MIN_MAPQ = min_mapq
    SNAP_MODE = snap_mode
    G38_MARGIN_BP = g38_bp
    G38_MARGIN_FRAC = g38_frac
    MAX_RATIO = max_ratio
    CLIP = clip
    MIN_OUTPUT_BP = min_out_bp
    EMIT_MODE = emit_mode
    FORCE_SNAP_ONLY = force_snap_only
    REJECT_VERBOSE = reject_verbose

# ---------- PAF parse ----------

def build_from_paf(paf_file: str, min_mapq: int = 5):
    gt = defaultdict(IntervalTree)
    grec = defaultdict(list)
    nm = defaultdict(list)
    plen = defaultdict(int)
    paf_bnd = defaultdict(set)

    cols = ['q_name','q_len','q_start','q_end','strand',
            't_name','t_len','t_start','t_end','matches','aln_len','map_q']
    try:
        total = sum(1 for _ in open(paf_file, 'r'))
    except FileNotFoundError:
        print(f"[ERROR] PAF not found: {paf_file}", file=sys.stderr)
        sys.exit(1)

    reader = pd.read_csv(paf_file, sep='\t', header=None,
                         chunksize=200_000, usecols=range(12), names=cols)

    with tqdm(total=total, unit=" lines", desc="PAF parse") as bar:
        for chunk in reader:
            chunk = chunk[chunk['map_q'] >= min_mapq]
            for _, r in chunk.iterrows():
                node   = parse_minigraph_id(str(r['t_name']))
                qid    = str(r['q_name'])
                qs     = int(r['q_start']);  qe   = int(r['q_end'])
                ql     = int(r['q_len'])
                ts     = int(r['t_start']);  te   = int(r['t_end'])
                tl     = int(r['t_len'])
                strand = '+' if str(r['strand']) == '+' else '-'

                if qid.startswith('id=GRCh38|'):
                    chrom = qid.split('|', 1)[1]
                    gt[chrom].addi(qs, qe, node)
                    grec[node].append({
                        "chrom": chrom, "qs": qs, "qe": qe,
                        "ts": ts, "te": te, "t_len": tl, "strand": strand
                    })
                else:
                    path = paf_id_to_og_path(qid)
                    nm[node].append({
                        "path": path, "qs": qs, "qe": qe, "q_len": ql,
                        "ts": ts, "te": te, "t_len": tl, "strand": strand
                    })
                    plen[path] = max(plen[path], ql)
                    paf_bnd[path].add(qs); paf_bnd[path].add(qe)
            bar.update(len(chunk))

    bnd = {p: sorted(v) for p, v in paf_bnd.items()}
    print(f"[INFO] Built maps: GRCh38_chroms={len(gt)} nodes={len(nm)} paths={len(plen)}", file=sys.stderr)
    return gt, grec, nm, dict(plen), bnd

# ---------- snap to PAF boundaries ----------

def _snap_endpoint(boundaries, x, mode, is_left):
    if not boundaries: return int(x)
    i = bisect.bisect_left(boundaries, x)
    if mode == "nearest":
        if i == 0: return boundaries[0]
        if i == len(boundaries): return boundaries[-1]
        L = boundaries[i-1]; R = boundaries[i]
        return L if (x - L) <= (R - x) else R
    if mode == "expand":
        if is_left:
            if i == 0: return boundaries[0]
            if i < len(boundaries) and boundaries[i] == x: return boundaries[i]
            return boundaries[i-1]
        else:
            return boundaries[i] if i < len(boundaries) else boundaries[-1]
    if mode == "shrink":
        if is_left:
            return boundaries[i] if i < len(boundaries) else boundaries[-1]
        else:
            if i < len(boundaries) and boundaries[i] == x: return boundaries[i]
            return boundaries[i-1] if i > 0 else boundaries[0]
    return int(x)

def snap_interval_to_paf(path, s, e, mode, lo=None, hi=None):
    b = BND.get(path)
    info = {"snap_mode": mode, "had_bounds": bool(b)}
    if not b:
        return int(s), int(e), "NO_BOUNDARIES", info
    if lo is not None or hi is not None:
        lo = -10**18 if lo is None else int(lo)
        hi =  10**18 if hi is None else int(hi)
        L = bisect.bisect_left(b, lo)
        R = bisect.bisect_right(b, hi)
        info["window"] = [lo, hi]
        info["n_bounds_in_window"] = R - L
        b = b[L:R]
        if not b:
            return int(s), int(e), "NO_BOUNDS_IN_WINDOW", info
    s2 = _snap_endpoint(b, int(math.floor(s)), mode, is_left=True)
    e2 = _snap_endpoint(b, int(math.ceil(e)),  mode, is_left=False)
    return int(s2), int(e2), "SNAPPED", info

# ---------- core ----------

def lift_one(job):
    chrom, g_s, g_e, itype, hap = job
    g_s = int(g_s); g_e = int(g_e)
    if g_e <= g_s:
        return [], [(chrom, g_s, g_e, hap, "ZERO_OR_NEG", "{}")], [], 0

    orig_len = g_e - g_s
    margin   = max(G38_MARGIN_BP, int(G38_MARGIN_FRAC * max(orig_len, 1)))
    G_L = g_s - margin
    G_R = g_e + margin

    if chrom not in GT:
        return [], [(chrom, g_s, g_e, hap, "NO_GRCH38", "{}")], [], 0

    overlaps = list(GT[chrom][g_s:g_e])
    if not overlaps:
        return [], [(chrom, g_s, g_e, hap, "NO_OVERLAP", "{}")], [], 0

    nodes = [iv.data for iv in overlaps]

    # node -> allowed target segments from GRCh38 true-subset
    per_path_spans = defaultdict(list)  # path -> [(qs,qe)]
    ns_grch_nodes  = 0
    ns_hap_hits    = 0

    for node in nodes:
        grec_list = [r for r in GREC.get(node, []) if r["chrom"] == chrom]
        if not grec_list:
            continue

        target_allowed = []
        for g in grec_list:
            loq = max(G_L, g["qs"]); hiq = min(G_R, g["qe"])
            if hiq <= loq: continue
            t0 = q2t(loq, g["qs"], g["qe"], g["ts"], g["te"], g["strand"])
            t1 = q2t(hiq, g["qs"], g["qe"], g["ts"], g["te"], g["strand"])
            if t0 is None or t1 is None: continue
            ta, tb = min(t0, t1), max(t0, t1)
            t_min = min(g["ts"], g["te"]); t_max = max(g["ts"], g["te"])
            ta = max(t_min, min(ta, t_max)); tb = max(t_min, min(tb, t_max))
            if tb > ta:
                target_allowed.append((ta, tb))
        target_allowed = merge_intervals(target_allowed)
        if not target_allowed:
            continue
        ns_grch_nodes += 1

        hit_this_node = False
        for m in NM.get(node, []):
            p = m["path"]
            if not path_matches_hap(p, hap):
                continue
            st_min = min(m["ts"], m["te"]); st_max = max(m["ts"], m["te"])
            for ta, tb in target_allowed:
                ia = max(ta, st_min); ib = min(tb, st_max)
                if ib <= ia: continue
                qa = t2q(ia, m["qs"], m["qe"], m["ts"], m["te"], m["strand"])
                qb = t2q(ib, m["qs"], m["qe"], m["ts"], m["te"], m["strand"])
                if qa is None or qb is None: continue
                s = int(math.floor(min(qa, qb))); e = int(math.ceil(max(qa, qb)))
                if e > s:
                    per_path_spans[p].append((s, e))
                    hit_this_node = True
        if hit_this_node:
            ns_hap_hits += 1

    if not per_path_spans:
        detail = {"ns_grch_nodes": ns_grch_nodes, "ns_hap_hit_nodes": ns_hap_hits}
        return [], [(chrom, g_s, g_e, hap, "NO_SAMPLE_HITS_ON_NODES", json.dumps(detail))], [], 0

    bed_out, rejects, stats_rows = [], [], []
    covered_this_job = 0

    for p, spans in per_path_spans.items():
        merged = merge_intervals(spans)
        if not merged:
            continue
        ql = PLEN.get(p, None)

        if EMIT_MODE == "fragments":
            for s0, e0 in merged:
                if e0 - s0 < MIN_OUTPUT_BP: continue
                s, e, snap_status, info = snap_interval_to_paf(p, s0, e0, SNAP_MODE, lo=s0, hi=e0)
                used = SNAP_MODE
                if e <= s:
                    s, e, snap_status2, info2 = snap_interval_to_paf(p, s0, e0, "nearest", lo=s0, hi=e0)
                    used = "nearest"
                    if e <= s and not FORCE_SNAP_ONLY:
                        s, e = s0, e0
                        used = "raw"
                        snap_status2 = "RAW_FALLBACK"
                if ql is not None and CLIP:
                    s = max(0, min(s, ql)); e = max(0, min(e, ql))
                if e - s < MIN_OUTPUT_BP:
                    if REJECT_VERBOSE:
                        detail = {"path": p, "mode": used, "snap_status": snap_status}
                        rejects.append((chrom, g_s, g_e, hap, "FRAG_DEGENERATE_AFTER_SNAP", json.dumps(detail)))
                    continue
                if orig_len > 0 and (e - s) > MAX_RATIO * orig_len:
                    if REJECT_VERBOSE:
                        detail = {"path": p, "span_len": e-s, "orig_len": orig_len}
                        rejects.append((chrom, g_s, g_e, hap, "FRAG_EXCEED_MAX_RATIO", json.dumps(detail)))
                    continue
                try:
                    contig = p.split('#')[2]
                except Exception:
                    contig = p.replace('#','_')
                tag = "FRAG_RAWFALLBACK" if used == "raw" else "FRAG"
                name = f"{hap}_{itype}_{chrom}_{g_s}_{g_e}_{contig}_{tag}_{md5_short(f'{p}_{s}_{e}')}"
                bed_out.append(f"{p}\t{s}\t{e}\t{name}\n")
            covered_this_job = 1 if merged else covered_this_job
        else:
            # JOIN mode emits one interval per path as [min_start, max_end].
            s0 = merged[0][0]; e0 = merged[-1][1]
            s, e, snap_status, info = snap_interval_to_paf(p, s0, e0, SNAP_MODE, lo=s0, hi=e0)
            used = SNAP_MODE
            if e <= s:
                s, e, snap_status2, info2 = snap_interval_to_paf(p, s0, e0, "nearest", lo=s0, hi=e0)
                used = "nearest"
                if e <= s and not FORCE_SNAP_ONLY:
                    s, e = s0, e0
                    used = "raw"
                    snap_status2 = "RAW_FALLBACK"
            if ql is not None and CLIP:
                s = max(0, min(s, ql)); e = max(0, min(e, ql))
            if e - s >= MIN_OUTPUT_BP:
                try:
                    contig = p.split('#')[2]
                except Exception:
                    contig = p.replace('#','_')
                tag = "JOIN_RAWFALLBACK" if used == "raw" else "JOIN"
                name = f"{hap}_{itype}_{chrom}_{g_s}_{g_e}_{contig}_{tag}_{md5_short(f'{p}_{s}_{e}')}"
                bed_out.append(f"{p}\t{s}\t{e}\t{name}\n")
                covered_this_job = 1
            else:
                if REJECT_VERBOSE:
                    detail = {"path": p, "mode": used, "snap_status": snap_status}
                    rejects.append((chrom, g_s, g_e, hap, "JOIN_DEGENERATE_AFTER_SNAP", json.dumps(detail)))

    if covered_this_job == 0 and not bed_out:
        return [], [(chrom, g_s, g_e, hap, "ALL_SPANS_FILTERED", "{}")], [], 0

    # Lightweight per-job statistics.
    stats_rows.append({
        "chr": chrom, "g_start": g_s, "g_end": g_e, "hap": hap, "type": itype,
        "n_paths": len(per_path_spans), "emit_mode": EMIT_MODE
    })

    return bed_out, rejects, stats_rows, covered_this_job

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="PAF-only true-subset liftover (JOIN default, inject-safe)")
    ap.add_argument('-p', '--paf', required=True)
    ap.add_argument('-i', '--introgression', required=True, help="chr start end type hap (whitespace-delimited)")
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count())

    # Permissive defaults without expanding the final output interval.
    ap.add_argument('--min-mapq', type=int, default=5)
    ap.add_argument('--g38-margin-bp', type=int, default=30)
    ap.add_argument('--g38-margin-frac', type=float, default=0.01)

    ap.add_argument('--max-ratio', type=float, default=1e9)
    ap.add_argument('--snap-mode', choices=['expand','shrink','nearest'], default='expand')
    ap.add_argument('--no-clip', action='store_true')
    ap.add_argument('--min-output-bp', type=int, default=1)

    ap.add_argument('--emit', choices=['join','fragments'], default='join',
                    help="join=one interval per path (default); fragments=emit all merged fragments")
    ap.add_argument('--force-snap-only', action='store_true',
                    help="disable the raw fallback; endpoints must snap to PAF boundaries")
    ap.add_argument('--no-reject-detail', action='store_true',
                    help="omit detailed rejection metadata to reduce output size")

    ap.add_argument('--rejected', default=None)
    ap.add_argument('--stats', default=None)

    args = ap.parse_args()

    rejected_path = args.rejected if args.rejected else f"{args.output}.rejected.tsv"
    stats_path    = args.stats if args.stats else f"{args.output}.stats.tsv"

    gt, grec, nm, plen, bnd = build_from_paf(args.paf, args.min_mapq)

    try:
        intro = pd.read_csv(args.introgression, sep=r'\s+', header=None,
                            names=['chr','start','end','type','hap'], engine='python')
    except FileNotFoundError:
        print(f"[ERROR] Introgression file not found: {args.introgression}", file=sys.stderr)
        sys.exit(1)
    jobs = [tuple(r) for r in intro.itertuples(index=False)]

    bed_lines, rej_rows, stats_rows = [], [], []
    covered_jobs = 0
    with multiprocessing.Pool(
        processes=args.threads,
        initializer=init_worker,
        initargs=(gt, grec, nm, plen, bnd,
                  args.min_mapq, args.snap_mode,
                  args.g38_margin_bp, args.g38_margin_frac,
                  args.max_ratio, not args.no_clip, args.min_output_bp,
                  args.emit, args.force_snap_only, not args.no_reject_detail)
    ) as pool:
        for bed, rej, st, covered in tqdm(pool.imap_unordered(lift_one, jobs),
                                          total=len(jobs), desc="Lifting"):
            if bed: bed_lines.extend(bed)
            if rej: rej_rows.extend(rej)
            if st:  stats_rows.extend(st)
            covered_jobs += covered

    with open(args.output, 'w') as f:
        f.writelines(bed_lines)

    # rejected with detail
    if rej_rows:
        rej_df = pd.DataFrame(rej_rows, columns=["chr","g_start","g_end","hap","reason","detail"])
    else:
        rej_df = pd.DataFrame(columns=["chr","g_start","g_end","hap","reason","detail"])
    rej_df.to_csv(rejected_path, sep='\t', index=False)

    # stats
    if stats_rows:
        pd.DataFrame(stats_rows).to_csv(stats_path, sep='\t', index=False)
    else:
        pd.DataFrame(columns=["chr","g_start","g_end","hap","type","n_paths","emit_mode"]).to_csv(stats_path, sep='\t', index=False)

    total_jobs = len(jobs)
    records = len(bed_lines)
    fail = len(rej_rows)
    job_success_rate = (covered_jobs / total_jobs * 100.0) if total_jobs else 0.0
    record_ratio = (records / total_jobs * 100.0) if total_jobs else 0.0

    print(f"[DONE] BED: {args.output} (records={records})", file=sys.stderr)
    print(f"[DONE] Rejected: {rejected_path} (rows={fail})", file=sys.stderr)
    print(f"[DONE] Stats: {stats_path} (rows={len(stats_rows)})", file=sys.stderr)
    print(f"[INFO] Jobs covered: {covered_jobs}/{total_jobs} ({job_success_rate:.1f}%)", file=sys.stderr)
    print(f"[INFO] Records/job: {record_ratio:.1f}% (can exceed 100% if multiple paths per job)", file=sys.stderr)

if __name__ == "__main__":
    main()
