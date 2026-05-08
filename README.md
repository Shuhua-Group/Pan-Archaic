# PanArchaic indel/SV discovery scripts

This repository contains the core scripts used for pangenome-resolved archaic indel and structural-variant discovery. It is intended as a compact methods and reproducibility package for reviewers. It does not include sequencing reads, graph files, container images, or scheduler-specific launch scripts.

## Workflow

1. `01_ancient_read_preprocessing/preprocess_ancient_reads_v5_20260324.sh`
   - Preprocess ancient sequencing reads for PanGenie and graph-based genotyping.

2. `02_interval_liftover_for_graph_injection/format_introgression_tracts.awk`
   and `02_interval_liftover_for_graph_injection/liftover_introgression_tracts_paf_v17.py`
   - Prepare and lift introgression-tract intervals onto graph/sample-path coordinates.

3. `03_graph_preparation/inject_introgression_paths_odgi.sh`
   and `03_graph_preparation/remove_introgression_paths_and_build_indexes.sh`
   - Inject lifted introgression intervals into the graph, which aligns graph node boundaries to introgression-tract boundaries.
   - Remove the injected introgression paths before building graph indexes, preserving the refined node granularity without retaining the injected paths as genotyping/haplotype evidence.

4. `03_panarchaic_wdl/panarchaic_introgression_indel_sv_v4.2.2.wdl`
   - Run the main tract-first, graph-aware workflow.
   - The WDL performs read preparation, k-mer counting, haplotype sampling, graph deconstruction, PanGenie genotyping, deconstruct/PanGenie VCF merging, introgressed-path extraction, and ODGI-based tracing of retained events.

5. `03_panarchaic_wdl/merge_deconstruct_pangenie_vcf.sh`
   and `03_panarchaic_wdl/trace_introgressed_variants_odgi.sh`
   - Helper scripts used by the WDL to merge VCF annotations and trace events through introgressed graph paths.

6. `04_post_wdl_vcf_standardization/standardize_introgressed_vcf_for_annotation.py`
   - Standardize final VCF records for downstream annotation and summary. This script assigns stable IDs, adds record-level/per-allele SV annotations, and removes same-length non-inversion replacements from the final indel/SV denominator.

See `03_graph_preparation/README.md` for the short rationale behind the inject-then-remove graph preparation step.

## Required inputs

The WDL expects these inputs to be supplied through a Cromwell-compatible JSON file:

- Ancient-sample FASTQ input, optionally with an already prepared uncompressed FASTQ for PanGenie.
- Boundary-aware graph/index inputs generated after injecting introgression intervals and then removing the injected paths: `original_gbz_file`, `original_hapl_file`, `original_ri_file`, and `original_snarls_file`.
- The introgression-injected ODGI graph before path removal: `injected_og_file`.
- Reference FASTA and index files.
- Container images for `vg`, PanGenie, ODGI, and a basic Unix utility environment.
- Helper script paths pointing to `merge_deconstruct_pangenie_vcf.sh` and `trace_introgressed_variants_odgi.sh`.

If a container runtime needs explicit bind mounts, set the optional WDL input `singularity_bind_args`. Otherwise it can be left as an empty string.

## Minimal execution pattern

```bash
java -jar cromwell.jar run \
  03_panarchaic_wdl/panarchaic_introgression_indel_sv_v4.2.2.wdl \
  -i inputs.json
```

The strict event VCF produced by the WDL is:

```text
<sample>.introgressed_svs.SVstrict.vcf.gz
```

## Post-WDL filtering and standardization

The analysis used genotype-quality and k-mer support filters before final VCF standardization:

```bash
bcftools view -i 'FORMAT/GQ>30 && FORMAT/KC>5' \
  input.SVstrict.vcf.gz \
  -Oz -o output.GQ30.KC5.vcf.gz

tabix -p vcf output.GQ30.KC5.vcf.gz

python 04_post_wdl_vcf_standardization/standardize_introgressed_vcf_for_annotation.py \
  -i output.GQ30.KC5.vcf.gz \
  > output.GQ30.KC5.preprocessed.vcf

bgzip -c output.GQ30.KC5.preprocessed.vcf \
  > output.GQ30.KC5.preprocessed.vcf.gz

tabix -p vcf output.GQ30.KC5.preprocessed.vcf.gz
```

## Repository contents

`SCRIPT_MANIFEST.tsv` summarizes each script, its role, expected inputs, and expected outputs. The package intentionally excludes scheduler wrappers, absolute-path input JSON files, plotting code, and downstream study-specific analyses.
