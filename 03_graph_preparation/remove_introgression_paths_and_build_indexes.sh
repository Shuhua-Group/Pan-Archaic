#!/usr/bin/env bash
# Remove injected introgression paths and build graph indexes for the WDL.
#
# This preserves the node partition created by odgi inject, but prevents the
# injected paths from being retained as haplotype/reference-like paths in the
# genotyping graph indexes.

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input.injected.og> <output_prefix> <output_dir> [threads]" >&2
    echo "Optional environment: ODGI_SIF, VG_SIF, SINGULARITY_BIND_ARGS, DROP_PATH_REGEX, REF_PATH_REGEX" >&2
    exit 1
fi

INJECTED_OG="$1"
OUT_PREFIX="$2"
OUT_DIR="$3"
THREADS="${4:-8}"
DROP_PATH_REGEX="${DROP_PATH_REGEX:-Neanderthal|Denisovan|Mosaic}"
REF_PATH_REGEX="${REF_PATH_REGEX:-CHM13#0#}"

if [ ! -f "$INJECTED_OG" ]; then
    echo "Error: injected graph not found: $INJECTED_OG" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
TMP_DIR=$(mktemp -d "${OUT_DIR}/tmp.${OUT_PREFIX}.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

CLEAN_OG="${OUT_DIR}/${OUT_PREFIX}.remove_introgression_paths.og"
GFA="${OUT_DIR}/${OUT_PREFIX}.gfa"
XG="${OUT_DIR}/${OUT_PREFIX}.xg"
GBZ="${OUT_DIR}/${OUT_PREFIX}.gbz"
DIST="${OUT_DIR}/${OUT_PREFIX}.dist"
RI="${OUT_DIR}/${OUT_PREFIX}.ri"
HAPL="${OUT_DIR}/${OUT_PREFIX}.hapl"
SNARLS="${OUT_DIR}/${OUT_PREFIX}.snarls"
SURJECT_PATHS="${OUT_DIR}/${OUT_PREFIX}.surject.paths.txt"
DROP_LIST="${TMP_DIR}/paths_to_drop.txt"

run_odgi() {
    if [ -n "${ODGI_SIF:-}" ]; then
        singularity exec ${SINGULARITY_BIND_ARGS:-} "$ODGI_SIF" odgi "$@"
    else
        odgi "$@"
    fi
}

run_vg() {
    if [ -n "${VG_SIF:-}" ]; then
        singularity exec ${SINGULARITY_BIND_ARGS:-} "$VG_SIF" vg "$@"
    else
        vg "$@"
    fi
}

echo "[INFO] Preparing boundary-aware graph indexes"
echo "  injected graph: $INJECTED_OG"
echo "  output prefix : $OUT_PREFIX"
echo "  output dir    : $OUT_DIR"
echo "  drop regex    : $DROP_PATH_REGEX"

run_odgi paths -i "$INJECTED_OG" -L | grep -E "$DROP_PATH_REGEX" > "$DROP_LIST" || true

if [ -s "$DROP_LIST" ]; then
    DROP_COUNT=$(wc -l < "$DROP_LIST")
    echo "[INFO] Removing $DROP_COUNT injected introgression paths"
    run_odgi paths -i "$INJECTED_OG" -X "$DROP_LIST" -o "$CLEAN_OG"
else
    echo "[INFO] No matching injected paths found; copying input graph"
    cp "$INJECTED_OG" "$CLEAN_OG"
fi

echo "[INFO] Converting cleaned graph to GFA"
run_odgi view -i "$CLEAN_OG" -g > "$GFA"

echo "[INFO] Building XG"
run_vg convert -t "$THREADS" -x -g "$GFA" > "$XG"

echo "[INFO] Building GBZ"
run_vg gbwt -p --gbz-format --graph-name "$GBZ" --gfa-input "$GFA"

echo "[INFO] Building distance index"
run_vg index -t "$THREADS" -j "$DIST" --no-nested-distance "$GBZ"

echo "[INFO] Building R-index"
run_vg gbwt -p --num-threads "$THREADS" -r "$RI" -Z "$GBZ"

echo "[INFO] Building haplotype index"
run_vg haplotypes -v 2 -t "$THREADS" -d "$DIST" -r "$RI" -H "$HAPL" "$GBZ"

echo "[INFO] Calling snarls"
run_vg snarls -t "$THREADS" "$XG" > "$SNARLS"

echo "[INFO] Writing optional surjection path list"
run_vg paths -L -x "$XG" | grep -E "$REF_PATH_REGEX" > "$SURJECT_PATHS" || true

echo "[OK] Boundary-aware graph preparation complete"
echo "  cleaned graph : $CLEAN_OG"
echo "  WDL GBZ input : $GBZ"
echo "  WDL HAPL input: $HAPL"
echo "  WDL RI input  : $RI"
echo "  WDL snarls    : $SNARLS"
