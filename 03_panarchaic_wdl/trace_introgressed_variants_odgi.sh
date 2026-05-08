#!/usr/bin/env bash
# Trace genotype-supported variant alleles through introgressed ODGI graph paths.
#
# Multi-allelic records are evaluated allele-specifically: the script parses the
# sample GT field, validates only ALT alleles carried by the sample, and checks
# whether candidate introgressed paths match the allele-specific graph topology.

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 ODGI_SIF OG_FILE FULL_VCF_GZ INTROGRESSED_PATHS_FILE SAMPLE [THREADS]"
    exit 1
fi

ODGI_SIF="$1"
OG_FILE="$2"
FULL_VCF="$3"
INTRO_PATHS="$4"
SAMPLE="$5"
THREADS="${6:-8}"

MIN_SV_LEN=2

GFA_FILE="${OG_FILE%.og}.gfa"
NODE2PATHS_FILE="node_to_paths.txt"

BASE_VCF="${SAMPLE}.introgressed_svs.vcf"
BASE_VCF_GZ="${BASE_VCF}.gz"
SV_VCF="${SAMPLE}.introgressed_svs.SV.vcf"
SV_VCF_GZ="${SV_VCF}.gz"
SV_STRICT_VCF="${SAMPLE}.introgressed_svs.SVstrict.vcf"
SV_STRICT_VCF_GZ="${SV_STRICT_VCF}.gz"

echo "=========================================="
echo "  ODGI introgressed variant tracer"
echo "=========================================="
echo "  Sample   : $SAMPLE"
echo "  Min len  : ${MIN_SV_LEN} bp"
echo "  Output   : ${BASE_VCF_GZ}"
echo "  Logic    : genotype-specific topology check"
echo ""

echo "[1/3] Building node-to-introgressed-path map"
if [ -s "$GFA_FILE" ]; then
    echo "  [OK] Reusing existing GFA: $GFA_FILE"
else
    echo "  [INFO] Generating GFA from ODGI graph"
    if [ -n "${SINGULARITY_CONTAINER:-}" ] || command -v odgi &> /dev/null; then
        odgi view -i "$OG_FILE" -g > "$GFA_FILE"
    else
        SINGULARITY_BIND_ARGS="${SINGULARITY_BIND_ARGS:-}"
        singularity exec ${SINGULARITY_BIND_ARGS} "$ODGI_SIF" odgi view -i "$OG_FILE" -g > "$GFA_FILE"
    fi
fi

if [ -s "$NODE2PATHS_FILE" ]; then
    echo "  [OK] Reusing existing node map: $NODE2PATHS_FILE"
else
    echo "  [INFO] Parsing GFA path records"
    awk -v path_list="$INTRO_PATHS" '
      BEGIN {
        while ((getline p < path_list) > 0) {
            gsub(/\r$/, "", p); if(p!="") want[p]=1
        }
      }
      /^P/ {
        path_name = $2
        if (!(path_name in want)) next
        n = split($3, nodes, ",")
        for(i=1; i<=n; i++) {
            node = nodes[i]
            gsub(/[+-]/, "", node)
            print node "\t" path_name
        }
      }
    ' "$GFA_FILE" | sort -u > "$NODE2PATHS_FILE"
fi

echo "[2/3] Running genotype-aware topology validation"
cat <<'EOF' > analyze_sv_topology_v6.py
import sys
import gzip
from collections import defaultdict

node_map_file = sys.argv[1]
vcf_file = sys.argv[2]
min_sv_len = int(sys.argv[3])

node_db = defaultdict(set)
with open(node_map_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            node_db[parts[0]].add(parts[1])

def parse_at_path(at_str):
    nodes = set()
    clean_str = at_str.replace('<', '>').replace(',', '')
    parts = clean_str.split('>')
    for p in parts:
        if p.isdigit():
            nodes.add(p)
    return nodes

def get_gt_alleles(gt_str):
    """Return ALT allele indices carried by a genotype string."""
    alleles = []
    if '.' in gt_str:
        return []

    gt_clean = gt_str.replace('|', '/')
    parts = gt_clean.split('/')

    for p in parts:
        if p.isdigit():
            val = int(p)
            if val > 0:
                alleles.append(val)
    return list(set(alleles))

opener = gzip.open if vcf_file.endswith('.gz') else open

with opener(vcf_file, 'rt') as fin:
    for line in fin:
        if line.startswith('#'):
            if line.startswith('##INFO=<ID=AT'):
                 print('##INFO=<ID=INTROG_CONFIRMED,Number=1,Type=String,Description="TRUE if a candidate introgressed path matches the topology of a GT-carried ALT allele">')
                 print('##INFO=<ID=SUPPORTING_PATHS,Number=.,Type=String,Description="Candidate introgressed paths supporting the variant">')
            print(line.strip())
            continue

        parts = line.strip().split('\t')

        fmt_fields = parts[8].split(':')
        try:
            gt_idx = fmt_fields.index('GT')
            sample_data = parts[9].split(':')
            gt_str = sample_data[gt_idx]
            alt_indices = get_gt_alleles(gt_str)
            if not alt_indices:
                continue
        except ValueError:
            continue

        ref_seq = parts[3]
        alt_seqs_raw = parts[4].split(',')
        all_seqs = [ref_seq] + alt_seqs_raw

        info = parts[7]
        at_val = None
        info_parts = info.split(';')
        for field in info_parts:
            if field.startswith('AT='):
                at_val = field[3:]
                break

        if not at_val:
            continue
        at_paths_all = at_val.split(',')

        valid_site = False
        all_confirmed_paths = set()

        for alt_idx in alt_indices:
            if alt_idx >= len(all_seqs) or alt_idx >= len(at_paths_all):
                continue

            target_alt_seq = all_seqs[alt_idx]
            len_ref = len(ref_seq)
            len_alt = len(target_alt_seq)
            diff = abs(len_ref - len_alt)
            max_len = max(len_ref, len_alt)
            if diff < min_sv_len and max_len < min_sv_len:
                continue

            ref_path_nodes = parse_at_path(at_paths_all[0])
            target_alt_path_nodes = parse_at_path(at_paths_all[alt_idx])

            unique_ref = ref_path_nodes - target_alt_path_nodes
            unique_alt = target_alt_path_nodes - ref_path_nodes
            anchors = ref_path_nodes & target_alt_path_nodes

            if not anchors:
                continue

            potential_paths = set()
            for node in anchors:
                if node in node_db:
                    potential_paths.update(node_db[node])

            for path in potential_paths:
                is_valid = True

                if unique_ref:
                    for n in unique_ref:
                        if path in node_db.get(n, set()):
                            is_valid = False
                            break
                if not is_valid:
                    continue

                if unique_alt:
                    for n in unique_alt:
                        if path not in node_db.get(n, set()):
                            is_valid = False
                            break
                if not is_valid:
                    continue

                valid_site = True
                all_confirmed_paths.add(path)

        if valid_site and all_confirmed_paths:
            path_str = ",".join(sorted(list(all_confirmed_paths)))
            new_info = f"{info};INTROG_CONFIRMED=TRUE;SUPPORTING_PATHS={path_str}"
            parts[7] = new_info
            print("\t".join(parts))
EOF

python3 analyze_sv_topology_v6.py "$NODE2PATHS_FILE" "$FULL_VCF" "$MIN_SV_LEN" > "$BASE_VCF"

echo "[3/3] Writing final VCF outputs"
if [ ! -s "$BASE_VCF" ]; then
    echo "[WARN] No qualifying variants were found; writing header-only VCF"
    zcat "$FULL_VCF" | grep "^#" > "$BASE_VCF"
fi

bgzip -@ "$THREADS" -c "$BASE_VCF" > "$BASE_VCF_GZ"
tabix -p vcf "$BASE_VCF_GZ"

cp "$BASE_VCF_GZ" "$SV_VCF_GZ"
cp "$BASE_VCF_GZ.tbi" "$SV_VCF_GZ.tbi"
cp "$BASE_VCF_GZ" "$SV_STRICT_VCF_GZ"
cp "$BASE_VCF_GZ.tbi" "$SV_STRICT_VCF_GZ.tbi"

echo "=========================================="
count=$(grep -v "^#" "$BASE_VCF" | wc -l || true)
echo "Completed genotype-aware introgressed variant tracing"
echo "Retained variant count: $count"
echo "Output VCF: $BASE_VCF_GZ"
echo "=========================================="
