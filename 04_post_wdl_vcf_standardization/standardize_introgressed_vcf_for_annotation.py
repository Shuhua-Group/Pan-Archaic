#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import gzip
import argparse
import hashlib
import re
from typing import List, Tuple, Optional, Dict

DNA_RE = re.compile(r'^[ACGTNacgtn]+$')
_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Retain only standard SV classes; same-length non-inversion replacements are dropped.
STANDARD_SVTYPES = {"INS", "DEL", "DUP", "INV", "BND"}

def open_text(path: str):
    """Open plain or gzipped text file. path == '-' means stdin."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")

def is_dna(seq: str) -> bool:
    return bool(DNA_RE.match(seq))

def revcomp(seq: str) -> str:
    return seq.translate(_TRANS)[::-1]

def sanitize_chrom(chrom: str) -> str:
    return re.sub(r"\s+", "_", chrom)

def parse_info_keys(info: str) -> Dict[str, str]:
    """Parse INFO into dict for key existence checks."""
    d = {}
    if info in {".", ""}:
        return d
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = ""
    return d

def add_info_field(info: str, key: str, value: str) -> str:
    """Append key=value to INFO safely, respecting '.'."""
    if info in {".", ""}:
        return f"{key}={value}"
    return f"{info};{key}={value}"

def make_stable_id(chrom: str, pos: str, main_svtype: str, ref: str, alt: str, orig_id: str) -> str:
    """
    Stable, near-zero collision ID.
    chr_pos_mainSVTYPE_hash (hash over chrom/pos/ref/alt/origid).
    """
    chrom_s = sanitize_chrom(chrom)
    payload = f"{chrom}\t{pos}\t{ref}\t{alt}\t{orig_id}".encode("utf-8", errors="ignore")
    h = hashlib.sha1(payload).hexdigest()[:10]
    main_svtype = re.sub(r"[^A-Za-z0-9_]+", "_", main_svtype if main_svtype else "SV")
    return f"{chrom_s}_{pos}_{main_svtype}_{h}"

def choose_main_svtype(svtypes: List[Optional[str]], svlens: List[Optional[int]]) -> str:
    """
    Choose main SVTYPE for record-level SVTYPE/ID.
    """
    priority = {"INV": 5, "DUP": 4, "DEL": 3, "INS": 2, "BND": 0}
    best = None
    best_score = (-1, -1)  # (size, priority)
    for t, l in zip(svtypes, svlens):
        if t is None:
            continue
        if l is None:
            size = 0
        else:
            size = abs(l) if t in {"INS", "DEL", "DUP"} else int(l)
        score = (size, priority.get(t, 0))
        if score > best_score:
            best_score = score
            best = t
    return best if best is not None else "SV"

def infer_svtype_for_allele(
    ref: str,
    alt: str,
    min_inv_len: int,
    min_dup_event_len: int
) -> Tuple[Optional[str], Optional[int], bool]:
    """
    Infer SVTYPE and SVLEN for one ALT allele.

    Key change vs previous version:
    - diff==0 (same-length) will be kept ONLY if it is a clear INV (revcomp and length>=min_inv_len).
    - otherwise it is treated as SNV/MNV-like replacement and DROPPED (no DELINS output).

    Returns (svtype, svlen, keep_allele)
      - INS/DUP: svlen = len(alt) - len(ref) (positive)
      - DEL:     svlen = len(alt) - len(ref) (negative)
      - INV:     svlen = event size (positive)
      - BND:     svlen = None
    """
    if alt in {"*", "<*>", "."}:
        return (None, None, False)

    # Breakend
    if ("[" in alt) or ("]" in alt):
        return ("BND", None, True)

    # Symbolic ALT
    if alt.startswith("<") and alt.endswith(">") and len(alt) > 2:
        inner = alt[1:-1].split(":")[0].upper()
        if inner in STANDARD_SVTYPES:
            return (inner, None, True)
        # unknown symbolic: keep but don't guess
        return (None, None, True)

    # Non-DNA allele: keep but don't guess
    if not (is_dna(ref) and is_dna(alt)):
        return (None, None, True)

    diff = len(alt) - len(ref)

    # Length-changing: INS/DEL; DUP only when long enough & strong tandem signature
    if diff != 0:
        if diff < 0:
            return ("DEL", diff, True)

        # diff > 0: default INS; conservative DUP only when inserted segment is long enough
        svtype = "INS"
        inserted_len = diff

        if inserted_len >= min_dup_event_len and alt.startswith(ref):
            inserted = alt[len(ref):]

            # Conservative DUP heuristics:
            # H1: inserted equals REF without anchor (requires long enough due to threshold)
            if len(ref) > 1 and len(inserted) == (len(ref) - 1) and inserted == ref[1:]:
                svtype = "DUP"
            # H2: inserted equals full REF (rare; also requires long enough)
            elif len(inserted) == len(ref) and inserted == ref:
                svtype = "DUP"
            else:
                # motif repeat and motif matches suffix of REF
                max_motif = min(50, len(inserted) // 2)
                for m in range(1, max_motif + 1):
                    if len(inserted) % m != 0:
                        continue
                    motif = inserted[:m]
                    if motif * (len(inserted) // m) != inserted:
                        continue
                    if len(ref) >= m and ref[-m:] == motif:
                        svtype = "DUP"
                        break

        return (svtype, diff, True)

    # diff == 0: keep ONLY if strong INV
    ref_len = len(ref)
    if ref_len >= min_inv_len and alt.upper() == revcomp(ref.upper()):
        size = ref_len - 1 if ref_len > 1 else ref_len
        return ("INV", size, True)

    # Otherwise drop same-length replacements (no DELINS)
    return (None, None, False)

def main():
    ap = argparse.ArgumentParser(
        description="Preprocess PanGenie/pangenome VCF to be AnnotSV-friendly. Preserves original INFO/FORMAT/samples and multi-allelic structure. Filters out same-length non-INV replacements (no DELINS)."
    )
    ap.add_argument("-i", "--input", required=True, help="Input VCF file (.vcf or .vcf.gz). Use '-' for stdin.")
    ap.add_argument("--min-inv-len", type=int, default=50,
                    help="Minimum same-length allele size (bp) to consider INV by revcomp (default: 50).")
    ap.add_argument("--min-dup-event-len", type=int, default=20,
                    help="Minimum inserted length (bp) to consider calling DUP (default: 20). Prevents 1bp A->AA being called DUP.")
    ap.add_argument("--rewrite-id", action="store_true", help="Rewrite ID to chr_pos_mainSVTYPE_hash and store original in ORIGID (default: on).")
    ap.add_argument("--no-rewrite-id", dest="rewrite_id", action="store_false", help="Do not rewrite ID.")
    ap.set_defaults(rewrite_id=True)
    ap.add_argument("--keep-all-nonstandard", action="store_true",
                    help="If set, keep records even when no ALT allele is kept (default: off).")
    args = ap.parse_args()

    existing_info_ids = set()
    injected = False

    kept_records = 0
    dropped_records = 0
    warnings = 0

    ins_called = 0
    del_called = 0
    dup_called = 0
    inv_called = 0
    bnd_called = 0

    out = sys.stdout
    err = sys.stderr

    with open_text(args.input) as f:
        for raw in f:
            if raw.startswith("#"):
                line = raw.rstrip("\n")
                if line.startswith("##INFO=<ID="):
                    m = re.match(r"##INFO=<ID=([^,>]+)", line)
                    if m:
                        existing_info_ids.add(m.group(1))

                if line.startswith("#CHROM") and not injected:
                    injections = []
                    if "SVTYPE" not in existing_info_ids:
                        injections.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (inferred if absent)">')
                    if "SVLEN" not in existing_info_ids:
                        injections.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length (inferred); DEL negative; INS/DUP positive; INV event size">')
                    if "END" not in existing_info_ids:
                        injections.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant on the reference allele (POS + len(REF) - 1 if inferred)">')
                    if "SVTYPEA" not in existing_info_ids:
                        injections.append('##INFO=<ID=SVTYPEA,Number=A,Type=String,Description="Per-ALT inferred SVTYPE (INS/DEL/DUP/INV/BND or .)">')
                    if "SVLENA" not in existing_info_ids:
                        injections.append('##INFO=<ID=SVLENA,Number=A,Type=Integer,Description="Per-ALT inferred SVLEN (DEL negative; INS/DUP positive; INV event size; . if unknown)">')
                    if "ORIGID" not in existing_info_ids:
                        injections.append('##INFO=<ID=ORIGID,Number=1,Type=String,Description="Original variant ID before preprocessing (if ID is rewritten)">')

                    for h in injections:
                        out.write(h + "\n")
                    injected = True

                out.write(line + "\n")
                continue

            line = raw.rstrip("\n")
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                # Malformed, pass through unchanged (safest)
                out.write(line + "\n")
                warnings += 1
                continue

            chrom, pos, vid, ref, alt, qual, flt, info = cols[:8]
            rest = cols[8:]  # FORMAT + samples etc (preserve verbatim)

            alt_alleles = alt.split(",")

            svtype_list: List[Optional[str]] = []
            svlen_list: List[Optional[int]] = []
            keep_any = False

            for a in alt_alleles:
                t, l, keep = infer_svtype_for_allele(
                    ref=ref,
                    alt=a,
                    min_inv_len=args.min_inv_len,
                    min_dup_event_len=args.min_dup_event_len
                )
                svtype_list.append(t)
                svlen_list.append(l)
                if keep:
                    keep_any = True
                    if t == "INS": ins_called += 1
                    elif t == "DEL": del_called += 1
                    elif t == "DUP": dup_called += 1
                    elif t == "INV": inv_called += 1
                    elif t == "BND": bnd_called += 1

            if not keep_any and not args.keep_all_nonstandard:
                dropped_records += 1
                continue

            info_dict = parse_info_keys(info)

            # END: infer from reference allele if missing
            if "END" not in info_dict:
                try:
                    end_val = str(int(pos) + len(ref) - 1)
                    info = add_info_field(info, "END", end_val)
                    info_dict["END"] = end_val
                except Exception:
                    pass

            # per-allele annotations (do not overwrite)
            if "SVTYPEA" not in info_dict:
                svtypea = ",".join([t if t is not None else "." for t in svtype_list])
                info = add_info_field(info, "SVTYPEA", svtypea)
                info_dict["SVTYPEA"] = svtypea

            if "SVLENA" not in info_dict:
                svlena = ",".join([str(l) if l is not None else "." for l in svlen_list])
                info = add_info_field(info, "SVLENA", svlena)
                info_dict["SVLENA"] = svlena

            # record-level SVTYPE/SVLEN (do not overwrite existing)
            main_svtype = choose_main_svtype(svtype_list, svlen_list)

            if "SVTYPE" not in info_dict and main_svtype in STANDARD_SVTYPES:
                info = add_info_field(info, "SVTYPE", main_svtype)
                info_dict["SVTYPE"] = main_svtype

            if "SVLEN" not in info_dict:
                best_l = None
                best_size = -1
                for t, l in zip(svtype_list, svlen_list):
                    if t is None or l is None:
                        continue
                    size = abs(l) if t in {"INS", "DEL", "DUP"} else int(l)
                    if size > best_size:
                        best_size = size
                        best_l = l
                if best_l is not None:
                    info = add_info_field(info, "SVLEN", str(int(best_l)))
                    info_dict["SVLEN"] = str(int(best_l))

            # rewrite ID if requested
            if args.rewrite_id:
                if vid != "." and "ORIGID" not in info_dict:
                    info = add_info_field(info, "ORIGID", vid)
                svtype_for_id = info_dict.get("SVTYPE", main_svtype)
                vid = make_stable_id(chrom, pos, svtype_for_id, ref, alt, vid)

            out.write("\t".join([chrom, pos, vid, ref, alt, qual, flt, info] + rest) + "\n")
            kept_records += 1

    err.write(
        f"[standardize_introgressed_vcf_for_annotation] done\n"
        f"  kept_records: {kept_records}\n"
        f"  dropped_records: {dropped_records}\n"
        f"  warnings_passed_through: {warnings}\n"
        f"  allele_type_counts (kept alleles): INS={ins_called}, DEL={del_called}, DUP={dup_called}, INV={inv_called}, BND={bnd_called}\n"
        f"  params: min_inv_len={args.min_inv_len}, min_dup_event_len={args.min_dup_event_len}\n"
    )

if __name__ == "__main__":
    main()
