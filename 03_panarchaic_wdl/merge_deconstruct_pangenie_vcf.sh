#!/usr/bin/env bash
# Merge deconstruct-derived variant records with PanGenie genotype fields.
#
# The script uses a text-normalized merge path to avoid gzip/bgzip compatibility
# issues across environments. Input paths are resolved before changing into the
# output directory.

set -euo pipefail

get_abs_path() {
  if [ -f "$1" ]; then
    echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
  else
    echo "$1"
  fi
}

if [ $# -lt 3 ]; then
    echo "Usage: $0 <deconstruct_vcf> <pangenie_vcf> <sample> [out_dir]"
    exit 1
fi

DECON_VCF=$(get_abs_path "$1")
PANG_VCF=$(get_abs_path "$2")
SAMPLE="$3"
OUTPUT_DIR="${4:-merge_deconstruct_pangenie_vcf_out}"

echo "[INFO] Starting VCF merge"
echo "  Source (deconstruct): $DECON_VCF"
echo "  Target (PanGenie):    $PANG_VCF"

if [ ! -f "$DECON_VCF" ]; then
    echo "Error: deconstruct VCF not found: $DECON_VCF" >&2
    exit 1
fi
if [ ! -f "$PANG_VCF" ]; then
    echo "Error: PanGenie VCF not found: $PANG_VCF" >&2
    exit 1
fi

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "[INFO] Step 1: decompressing inputs"
if [[ "$DECON_VCF" == *.gz ]]; then
    zcat "$DECON_VCF" > decon.vcf
else
    cat "$DECON_VCF" > decon.vcf
fi

if [[ "$PANG_VCF" == *.gz ]]; then
    zcat "$PANG_VCF" > pang.vcf
else
    cat "$PANG_VCF" > pang.vcf
fi

echo "[INFO] Step 2: loading PanGenie genotype fields"
grep -v "^#" pang.vcf | awk -v OFS="\t" '{key=$2":"$4":"$5; print key, $10}' > pang_data.txt

echo "[INFO] Step 3: merging records"
OUTPUT_VCF="${SAMPLE}.merged.vcf"

awk -v OFS="\t" '
    BEGIN {
        while ((getline < "pang_data.txt") > 0) {
            split($0, a, "\t")
            pang_map[a[1]] = a[2]
        }
        close("pang_data.txt")
    }
    /^##/ {
        print $0
        next
    }
    /^#CHROM/ {
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
        print "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihoods\">"
        print "##FORMAT=<ID=KC,Number=1,Type=Float,Description=\"K-mer counts\">"
        print $0
        next
    }
    {
        key = $2":"$4":"$5
        if (key in pang_map) {
            $9 = "GT:GQ:GL:KC"
            $10 = pang_map[key]
            print $0
        } else {
            print $0
        }
    }
' decon.vcf > "$OUTPUT_VCF"

echo "[INFO] Step 4: compressing and indexing"
bgzip -f "$OUTPUT_VCF"
tabix -f -p vcf "${OUTPUT_VCF}.gz"
rm decon.vcf pang.vcf pang_data.txt

echo "[SUCCESS] Merge completed: $(pwd)/${OUTPUT_VCF}.gz"

echo "[CHECK] Verifying header and first data row"
if zgrep -q "ID=KC," "${OUTPUT_VCF}.gz"; then
    echo "  [OK] KC FORMAT definition present"
else
    echo "  [WARN] KC FORMAT definition missing"
fi

FIRST_DATA=$(zgrep -v "^#" "${OUTPUT_VCF}.gz" | head -n 1 || true)
if [ -n "$FIRST_DATA" ]; then
    echo "  First sample field: $(echo "$FIRST_DATA" | awk '{print $9, $10}')"
else
    echo "  First sample field: no data lines found"
fi

exit 0
