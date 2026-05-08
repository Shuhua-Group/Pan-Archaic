#!/usr/bin/env bash
# Inject lifted introgression intervals into an ODGI graph.
#
# The injected graph is used to align graph node boundaries to introgression
# tract boundaries and to retain path-level topology for later tracing.

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input.og> <introgression_intervals.bed> <output.injected.og> [threads]" >&2
    echo "Optional environment: ODGI_SIF, SINGULARITY_BIND_ARGS" >&2
    exit 1
fi

INPUT_OG="$1"
INTERVAL_BED="$2"
OUTPUT_OG="$3"
THREADS="${4:-8}"

if [ ! -f "$INPUT_OG" ]; then
    echo "Error: input graph not found: $INPUT_OG" >&2
    exit 1
fi
if [ ! -f "$INTERVAL_BED" ]; then
    echo "Error: interval BED not found: $INTERVAL_BED" >&2
    exit 1
fi

run_odgi() {
    if [ -n "${ODGI_SIF:-}" ]; then
        singularity exec ${SINGULARITY_BIND_ARGS:-} "$ODGI_SIF" odgi "$@"
    else
        odgi "$@"
    fi
}

echo "[INFO] Injecting introgression intervals into graph"
echo "  input graph : $INPUT_OG"
echo "  intervals   : $INTERVAL_BED"
echo "  output graph: $OUTPUT_OG"

run_odgi inject \
    -i "$INPUT_OG" \
    -b "$INTERVAL_BED" \
    -o "$OUTPUT_OG" \
    -t "$THREADS" \
    -P

echo "[OK] Injected graph written: $OUTPUT_OG"
