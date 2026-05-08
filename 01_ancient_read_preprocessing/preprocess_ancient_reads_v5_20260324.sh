#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# aDNA AutoPrep for PanGenie input
#
# Main flow:
#   1. AdapterRemoval for adapter trimming / quality trimming / fixed trimming
#   2. Collapse overlapping PE reads with AdapterRemoval
#   3. Pool retained reads into one final FASTQ for PanGenie
#   4. Optional fastp dedup-only step on the pooled FASTQ
#
# Default trimming policy:
#   - fixed trim: --trim5p 2 --trim3p 2
#   - minimum final read length: 35 bp AFTER trimming
#
# Note on the reference paper:
#   - the paper used fixed trimming of 5 bp from both ends plus fastp -l 31
#   - fastp -l / --length_required is a post-processing length threshold
#   - here we use trim2 + minlength 35 as a more practical compromise
#
# Note on quality trimming:
#   - trim2 is fixed trimming and is separate from quality trimming
#   - --trimqualities in AdapterRemoval uses the default --minquality=2
#   - we intentionally keep that default and do not raise it here
#
# Note on adapter handling:
#   - AdapterRemoval --identify-adapters is documented for PE overlap-based
#     adapter inference from fully overlapping paired-end reads
#   - therefore we use it only for PE
#   - for SE we do NOT use --identify-adapters
#   - instead, SE input relies on AdapterRemoval's official default adapter1
#
# Optional fastp step:
#   - fastp is used only for FASTQ-level dedup if enabled
#   - dedup-only means: --dedup with -A -Q -L to disable adapter trimming,
#     quality filtering, and length filtering
#   - this is not BAM coordinate deduplication
# ============================================================

READS_ROOT=""
OUTDIR="work_processed_final"
THREADS=8
MINLEN=35
FIXED_TRIM=2
AR_BIN="AdapterRemoval"
FASTP_BIN="fastp"
INCLUDE_UNMERGED_IN_FINAL="yes"
RUN_FASTP_DEDUP="no"

FINAL_GZ_FILENAME="Altai_Final_Input_for_PanGenie_and_GraphSE.fq.gz"
FINAL_FQ_FILENAME="Altai_Final_Input_for_PanGenie_and_GraphSE.fq"
FINAL_DEDUP_FQ_FILENAME="Altai_Final_Input_for_PanGenie_and_GraphSE.dedup.fq"
FINAL_DEDUP_GZ_FILENAME="Altai_Final_Input_for_PanGenie_and_GraphSE.dedup.fq.gz"
FASTP_JSON_FILENAME="Altai_Final_Input_for_PanGenie_and_GraphSE.fastp_dedup.json"
FASTP_HTML_FILENAME="Altai_Final_Input_for_PanGenie_and_GraphSE.fastp_dedup.html"

die() { echo "[FATAL] $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

usage() {
  cat <<EOF
Usage:
  bash $(basename "$0") -R <reads_dir> -o <outdir> [options]

Required:
  -R, --reads-root <dir>          Directory containing raw FASTQ(.gz) files

Optional:
  -o, --outdir <dir>              Output directory (default: ${OUTDIR})
  -t, --threads <int>             Threads for AdapterRemoval / fastp (default: ${THREADS})
  --minlen <int>                  Minimum final read length AFTER trimming (default: ${MINLEN})
  --fixed-trim <int>              Fixed trim on both ends (default: ${FIXED_TRIM})
  --ar-bin <path>                 AdapterRemoval binary (default: ${AR_BIN})
  --fastp-bin <path>              fastp binary (default: ${FASTP_BIN})
  --include-unmerged-in-final yes|no
                                  Include unmerged pair1/pair2 in final pooled FASTQ
                                  (default: ${INCLUDE_UNMERGED_IN_FINAL})
  --run-fastp-dedup yes|no
                                  Run fastp dedup-only after final pooling
                                  (default: ${RUN_FASTP_DEDUP})
  --help                          Show this help

Outputs:
  - clean/reads/*.merged.fq.gz
  - clean/reads/*.singletons.fq.gz
  - clean/reads/*.pair1.fq.gz
  - clean/reads/*.pair2.fq.gz
  - clean/reads/*.clean.fq.gz
  - summary/unmerged_PE_summary.tsv
  - ${FINAL_GZ_FILENAME}
  - ${FINAL_FQ_FILENAME}
  - optional:
      ${FINAL_DEDUP_FQ_FILENAME}
      ${FINAL_DEDUP_GZ_FILENAME}
      ${FASTP_JSON_FILENAME}
      ${FASTP_HTML_FILENAME}
EOF
}

count_fq_reads() {
  local fq="$1"
  [[ -s "$fq" ]] || { echo 0; return 0; }
  gzip -dc "$fq" 2>/dev/null | awk 'END{print NR/4}'
}

move_first_existing() {
  local dst="$1"
  shift
  local f
  for f in "$@"; do
    if [[ -f "$f" ]]; then
      mv "$f" "$dst"
      return 0
    fi
  done
  return 1
}

make_tag() {
  local name="$1"
  name="${name// /_}"
  name="${name//[^A-Za-z0-9._-]/_}"
  printf "%s" "$name"
}

# Global vars returned:
#   DETECTED_A1
#   DETECTED_A2
identify_adapters_pe() {
  local f1="$1"
  local f2="$2"
  local tag="$3"
  local ident_log="${OUTDIR}/temp/${tag}.ident.log"

  log "  > Identifying adapters for PE library ${tag} ..."

  "$AR_BIN" \
    --file1 "$f1" \
    --file2 "$f2" \
    --identify-adapters \
    --threads "$THREADS" \
    > "$ident_log" 2>&1 || true

  DETECTED_A1=$(grep -m 1 -- "--adapter1:" "$ident_log" | awk '{print $2}' || true)
  DETECTED_A2=$(grep -m 1 -- "--adapter2:" "$ident_log" | awk '{print $2}' || true)

  if [[ -z "$DETECTED_A1" ]]; then
    log "  [WARN] PE adapter identification failed for ${tag}. Using AdapterRemoval defaults."
    DETECTED_A1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
    DETECTED_A2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
  elif [[ -z "$DETECTED_A2" ]]; then
    DETECTED_A2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
  fi

  rm -f "$ident_log"
}

ARGS=$(getopt -o R:o:t: -l reads-root:,outdir:,threads:,minlen:,fixed-trim:,ar-bin:,fastp-bin:,include-unmerged-in-final:,run-fastp-dedup:,help -- "$@") || exit 1
eval set -- "$ARGS"
while true; do
  case "$1" in
    -R|--reads-root) READS_ROOT="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --minlen) MINLEN="$2"; shift 2;;
    --fixed-trim) FIXED_TRIM="$2"; shift 2;;
    --ar-bin) AR_BIN="$2"; shift 2;;
    --fastp-bin) FASTP_BIN="$2"; shift 2;;
    --include-unmerged-in-final) INCLUDE_UNMERGED_IN_FINAL="$2"; shift 2;;
    --run-fastp-dedup) RUN_FASTP_DEDUP="$2"; shift 2;;
    --help) usage; exit 0;;
    --) shift; break;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -n "$READS_ROOT" ]] || die "Need -R/--reads-root"
[[ -d "$READS_ROOT" ]] || die "READS_ROOT does not exist: $READS_ROOT"
command -v "$AR_BIN" >/dev/null 2>&1 || die "AdapterRemoval not found: $AR_BIN"
[[ "$INCLUDE_UNMERGED_IN_FINAL" == "yes" || "$INCLUDE_UNMERGED_IN_FINAL" == "no" ]] || \
  die "--include-unmerged-in-final must be yes or no"
[[ "$RUN_FASTP_DEDUP" == "yes" || "$RUN_FASTP_DEDUP" == "no" ]] || \
  die "--run-fastp-dedup must be yes or no"
if [[ "$RUN_FASTP_DEDUP" == "yes" ]]; then
  command -v "$FASTP_BIN" >/dev/null 2>&1 || die "fastp not found: $FASTP_BIN"
fi

mkdir -p \
  "$OUTDIR"/logs \
  "$OUTDIR"/temp \
  "$OUTDIR"/clean/reads \
  "$OUTDIR"/summary

SUMMARY_TSV="${OUTDIR}/summary/unmerged_PE_summary.tsv"
printf "tag\traw_pairs\tmerged_pairs\tsingleton_reads\tunmerged_pairs\tunmerged_pair_fraction_of_raw_pairs\tunmerged_read_fraction_in_all_retained_reads\n" > "$SUMMARY_TSV"

run_pe_smart() {
  local R1="$1"
  local R2="$2"
  local tag="$3"

  local base="${OUTDIR}/temp/${tag}"
  local logf="${OUTDIR}/logs/${tag}.PE.log"

  log "[PE] Processing ${tag} ..."

  identify_adapters_pe "$R1" "$R2" "$tag"
  log "     A1: $DETECTED_A1"
  log "     A2: $DETECTED_A2"

  local raw_pairs
  raw_pairs=$(count_fq_reads "$R1")

  local AR_ARGS=(
    --file1 "$R1"
    --file2 "$R2"
    --basename "$base"
    --threads "$THREADS"
    --gzip
    --trimns
    --trimqualities
    --minlength "$MINLEN"
    --qualitybase 33
    --qualitymax 93
    --minalignmentlength 11
    --mm 3
    --collapse
    --trim5p "$FIXED_TRIM"
    --trim3p "$FIXED_TRIM"
    --adapter1 "$DETECTED_A1"
    --adapter2 "$DETECTED_A2"
  )

  "$AR_BIN" "${AR_ARGS[@]}" > "$logf" 2>&1

  move_first_existing \
    "${OUTDIR}/clean/reads/${tag}.merged.fq.gz" \
    "${base}.collapsed.truncated.gz" \
    "${base}.collapsed.gz" || true

  move_first_existing \
    "${OUTDIR}/clean/reads/${tag}.singletons.fq.gz" \
    "${base}.singleton.truncated.gz" \
    "${base}.singleton.gz" || true

  move_first_existing \
    "${OUTDIR}/clean/reads/${tag}.pair1.fq.gz" \
    "${base}.pair1.truncated.gz" \
    "${base}.pair1.gz" || true

  move_first_existing \
    "${OUTDIR}/clean/reads/${tag}.pair2.fq.gz" \
    "${base}.pair2.truncated.gz" \
    "${base}.pair2.gz" || true

  local merged_pairs singleton_reads pair1_reads pair2_reads unmerged_pairs
  merged_pairs=$(count_fq_reads "${OUTDIR}/clean/reads/${tag}.merged.fq.gz")
  singleton_reads=$(count_fq_reads "${OUTDIR}/clean/reads/${tag}.singletons.fq.gz")
  pair1_reads=$(count_fq_reads "${OUTDIR}/clean/reads/${tag}.pair1.fq.gz")
  pair2_reads=$(count_fq_reads "${OUTDIR}/clean/reads/${tag}.pair2.fq.gz")

  if [[ "$pair1_reads" -ne "$pair2_reads" ]]; then
    log "[WARN] ${tag}: pair1_reads (${pair1_reads}) != pair2_reads (${pair2_reads})"
  fi

  if [[ "$pair1_reads" -le "$pair2_reads" ]]; then
    unmerged_pairs="$pair1_reads"
  else
    unmerged_pairs="$pair2_reads"
  fi

  awk \
    -v tag="$tag" \
    -v raw="$raw_pairs" \
    -v merged="$merged_pairs" \
    -v singles="$singleton_reads" \
    -v unmerged="$unmerged_pairs" \
    'BEGIN{
      pair_frac = (raw > 0 ? unmerged / raw : 0);
      read_frac = ((merged + singles + 2*unmerged) > 0 ? (2*unmerged) / (merged + singles + 2*unmerged) : 0);
      printf "%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\n", tag, raw, merged, singles, unmerged, pair_frac, read_frac;
    }' >> "$SUMMARY_TSV"

  rm -f "${base}"*
}

run_se_smart() {
  local SE="$1"
  local tag="$2"

  local base="${OUTDIR}/temp/${tag}"
  local logf="${OUTDIR}/logs/${tag}.SE.log"

  log "[SE] Processing ${tag} ..."
  log "  > Using AdapterRemoval default adapter1 for SE input ${tag}."

  local AR_ARGS=(
    --file1 "$SE"
    --basename "$base"
    --threads "$THREADS"
    --gzip
    --trimns
    --trimqualities
    --minlength "$MINLEN"
    --qualitybase 33
    --qualitymax 93
    --trim5p "$FIXED_TRIM"
    --trim3p "$FIXED_TRIM"
  )

  "$AR_BIN" "${AR_ARGS[@]}" > "$logf" 2>&1

  move_first_existing \
    "${OUTDIR}/clean/reads/${tag}.clean.fq.gz" \
    "${base}.truncated.gz" \
    "${base}.gz" || true

  rm -f "${base}"*
}

shopt -s nullglob

log ">>> Starting PE Processing..."
for R1 in "${READS_ROOT}"/*_1.fastq.gz "${READS_ROOT}"/*_1.fq.gz; do
  [[ -e "$R1" ]] || continue
  R2="${R1/_1.fastq/_2.fastq}"
  R2="${R2/_1.fq/_2.fq}"
  if [[ -f "$R2" ]]; then
    bn=$(basename "$R1")
    tag="${bn%_1.fastq.gz}"
    tag="${tag%_1.fq.gz}"
    tag=$(make_tag "$tag")
    run_pe_smart "$R1" "$R2" "$tag"
  fi
done

log ">>> Starting SE Processing..."
for SE in "${READS_ROOT}"/*.fastq.gz "${READS_ROOT}"/*.fq.gz; do
  [[ -e "$SE" ]] || continue
  bn=$(basename "$SE")
  if [[ "$bn" != *_1.fastq.gz && "$bn" != *_2.fastq.gz && "$bn" != *_1.fq.gz && "$bn" != *_2.fq.gz ]]; then
    tag="${bn%.fastq.gz}"
    tag="${tag%.fq.gz}"
    tag=$(make_tag "$tag")
    run_se_smart "$SE" "$tag"
  fi
done

log ">>> Performing FINAL GRAND MERGE..."

FINAL_GZ="${OUTDIR}/${FINAL_GZ_FILENAME}"
FINAL_FQ="${OUTDIR}/${FINAL_FQ_FILENAME}"
FINAL_DEDUP_FQ="${OUTDIR}/${FINAL_DEDUP_FQ_FILENAME}"
FINAL_DEDUP_GZ="${OUTDIR}/${FINAL_DEDUP_GZ_FILENAME}"
FASTP_JSON="${OUTDIR}/${FASTP_JSON_FILENAME}"
FASTP_HTML="${OUTDIR}/${FASTP_HTML_FILENAME}"

files=()

for f in "${OUTDIR}/clean/reads/"*.merged.fq.gz; do
  [[ -e "$f" ]] && files+=("$f")
done
for f in "${OUTDIR}/clean/reads/"*.singletons.fq.gz; do
  [[ -e "$f" ]] && files+=("$f")
done
for f in "${OUTDIR}/clean/reads/"*.clean.fq.gz; do
  [[ -e "$f" ]] && files+=("$f")
done

if [[ "$INCLUDE_UNMERGED_IN_FINAL" == "yes" ]]; then
  for f in "${OUTDIR}/clean/reads/"*.pair1.fq.gz; do
    [[ -e "$f" ]] && files+=("$f")
  done
  for f in "${OUTDIR}/clean/reads/"*.pair2.fq.gz; do
    [[ -e "$f" ]] && files+=("$f")
  done
fi

if (( ${#files[@]} == 0 )); then
  die "No cleaned reads found for final merge."
fi

cat "${files[@]}" > "$FINAL_GZ"
gzip -dc "$FINAL_GZ" > "$FINAL_FQ"

log ">>> SUCCESS! Final pooled outputs generated:"
ls -lh "$FINAL_GZ" "$FINAL_FQ"

if [[ "$RUN_FASTP_DEDUP" == "yes" ]]; then
  log ">>> Running fastp dedup-only on pooled FASTQ..."
  "$FASTP_BIN" \
    --in1 "$FINAL_FQ" \
    --out1 "$FINAL_DEDUP_FQ" \
    --thread "$THREADS" \
    --dedup \
    -A \
    -Q \
    -L \
    --json "$FASTP_JSON" \
    --html "$FASTP_HTML"

  gzip -c "$FINAL_DEDUP_FQ" > "$FINAL_DEDUP_GZ"
  log ">>> fastp dedup-only outputs generated:"
  ls -lh "$FINAL_DEDUP_FQ" "$FINAL_DEDUP_GZ" "$FASTP_JSON" "$FASTP_HTML"
fi

cat <<EOF

============================================================
Done.

Key outputs:
  1) Pooled FASTQ.GZ      : ${FINAL_GZ}
  2) Pooled FASTQ         : ${FINAL_FQ}
  3) Unmerged summary     : ${SUMMARY_TSV}
EOF

if [[ "$RUN_FASTP_DEDUP" == "yes" ]]; then
  cat <<EOF
  4) Dedup FASTQ          : ${FINAL_DEDUP_FQ}
  5) Dedup FASTQ.GZ       : ${FINAL_DEDUP_GZ}
  6) fastp JSON           : ${FASTP_JSON}
  7) fastp HTML           : ${FASTP_HTML}
EOF
fi

cat <<EOF

Current setting:
  fixed_trim                 = ${FIXED_TRIM}
  minlen_after_trimming      = ${MINLEN}
  include_unmerged_in_final  = ${INCLUDE_UNMERGED_IN_FINAL}
  run_fastp_dedup            = ${RUN_FASTP_DEDUP}

Interpretation:
  trim2 is fixed trimming.
  --trimqualities is enabled with AdapterRemoval default minquality=2.
  minlength is enforced after trimming.

============================================================
EOF

