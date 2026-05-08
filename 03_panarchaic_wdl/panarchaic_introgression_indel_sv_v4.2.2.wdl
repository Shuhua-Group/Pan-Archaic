version 1.0

# ============================================================================
# PanArchaic v4.2.2: Fast Mode - FIXED (Removed BuildMinimizerIndex) + LV=0 Fix
# ============================================================================
# Based on v4.2, streamlined for rapid SV detection without alignment visualization
#
# v4.2.2 Changes (vs v4.2):
# - **Removed BuildMinimizerIndex task** - not used by any downstream tasks!
# - **Removed minimizer_k/minimizer_w parameters** - no longer needed
# - **Removed original_dist_file input** - no longer needed
# - **Result**: Saves 5-15 minutes per task, eliminates segfault issue
#
# Fix for PanGenie (Current):
# - Added LV=0 filtering for PanGenie input to prevent overlapping nested variants
#
# Retained workflow:
# 1. PrepareReadsForPanGenie (decompress for PanGenie)
# 2. KmerCountingKMC (k-mer counting, uses .gz)
# 3. HaplotypeSampling (diploid graph construction)
# 4. Deconstruct (VCF from sampled graph)
# 5. PanGenieIndex (build genotyping index)
# 6. PanGenie (k-mer based genotyping)
# 7. MergeVCF (combine Deconstruct + PanGenie)
# 8. ExtractIntrogressedPaths (get archaic paths)
# 9. TraceSV_odgi (trace SVs to introgressed segments)
# ============================================================================

workflow AncientIntrogressionSVDetection_FastMode {
    input {
        # ===== Input Reads (Single-End, can be .fq.gz) =====
        File input_fastq
        String sample = "AncientSample"


        # Optional: if provided, skip PrepareReadsForPanGenie and use this as PanGenie reads (must be uncompressed .fq)
        File? prepared_fastq_for_pangenie

        # ===== Boundary-aware graph/index files =====
        # These are built after introgression-path injection and subsequent
        # removal of the injected paths. The node partition remains aligned to
        # introgression boundaries, but the injected paths are not used as
        # genotyping/haplotype evidence.
        File original_gbz_file
        File original_hapl_file
        File original_ri_file
        File original_snarls_file
        # Note: original_dist_file removed - not needed in Fast Mode

        # ===== Injected Graph Files =====
        File injected_og_file

        # ===== Introgressed Path Patterns =====
        Array[String] introgressed_path_patterns = ["Neanderthal", "Denisovan", "Mosaic"]

        # ===== Reference Files =====
        File in_ref_file
        File in_ref_index_file

        # ===== Parameters (aDNA Optimized) =====
        String pangenie_prefix = "pangenie_out"
        Int threads = 16
        String chromosome = "chr22"
        Int kmer_length = 29
        # Note: minimizer_k and minimizer_w removed - not needed

        # ===== Singularity Images =====
        File vg_sif
        File ubuntu_sif
        File pangenie_sif
        File odgi_sif
        String singularity_bind_args = ""

        # ===== Scripts =====
        File merge_vcf_script
        File trace_sv_script
    }

    # ========================================================================
    # PHASE 1: Genotyping Pipeline (No Alignment, No Minimizer Index)
    # ========================================================================

    if (!defined(prepared_fastq_for_pangenie)) {
        call PrepareReadsForPanGenie {
            input:
                input_fastq = input_fastq,
                sample = sample,
                ubuntu_sif = ubuntu_sif
        }
    }

    # PanGenie requires an uncompressed FASTQ; reuse a pre-prepared FASTQ if provided
    File pangenie_reads = select_first([prepared_fastq_for_pangenie, PrepareReadsForPanGenie.prepared_fastq])

    call KmerCountingKMC {
        input:
            input_fastq = input_fastq,
            kmer_length = kmer_length,
            sample = sample
    }

    call HaplotypeSampling {
        input:
            kff_file = KmerCountingKMC.kff_file,
            gbz_file = original_gbz_file,
            hapl_file = original_hapl_file,
            sample = sample,
            vg_sif = vg_sif,
    }

    call Deconstruct {
        input:
            sampled_gbz_file = HaplotypeSampling.sampled_gbz,
            sample = sample,
            chromosome = chromosome,
            threads = threads,
            vg_sif = vg_sif,
            singularity_bind_args = singularity_bind_args
    }

    call PanGenieIndex {
        input:
            in_ref_file = in_ref_file,
            in_ref_index_file = in_ref_index_file,
            variants_vcf = Deconstruct.pangenie_input_vcf,
            pangenie_prefix = pangenie_prefix,
            threads = threads,
            pangenie_sif = pangenie_sif
    }

    call PanGenie {
        input:
            reads = pangenie_reads,
            index_prefix = pangenie_prefix,
            graph_cereal_files = PanGenieIndex.graph_cereal_files,
            kmers_tsv_gz_files = PanGenieIndex.kmers_tsv_gz_files,
            unique_kmers_map = PanGenieIndex.unique_kmers_map,
            path_segments_fasta = PanGenieIndex.path_segments_fasta,
            pangenie_prefix = pangenie_prefix,
            sample = sample,
            threads = threads,
            pangenie_sif = pangenie_sif
    }

    call MergeVCF {
        input:
            deconstruct_vcf = Deconstruct.full_vcf_gz,
            pangenie_vcf = PanGenie.genotype_results,
            sample = sample,
            merge_script = merge_vcf_script
    }

    # ========================================================================
    # PHASE 2: Archaic Introgression Tracing
    # ========================================================================

    call ExtractIntrogressedPaths {
        input:
            injected_og_file = injected_og_file,
            path_patterns = introgressed_path_patterns,
            sample = sample,
            odgi_sif = odgi_sif
    }

    call TraceSV_odgi {
        input:
            injected_og_file = injected_og_file,
            merged_vcf = MergeVCF.merged_vcf_gz,
            introgressed_paths_file = ExtractIntrogressedPaths.paths_file,
            sample = sample,
            threads = threads,
            odgi_sif = odgi_sif,
            trace_script = trace_sv_script,
            singularity_bind_args = singularity_bind_args
    }

    output {
        # Core outputs (VCF only, no BAM, no minimizer index)
        File sampled_gbz = HaplotypeSampling.sampled_gbz
        File deconstruct_vcf = Deconstruct.full_vcf_gz
        File pangenie_vcf = PanGenie.genotype_results
        File merged_vcf = MergeVCF.merged_vcf_gz

        # Introgression-specific outputs
        File introgressed_paths = ExtractIntrogressedPaths.paths_file

        # Final SV outputs (3-tier)
        File introgressed_svs_base_vcf = TraceSV_odgi.base_vcf_gz
        File introgressed_svs_sv_vcf = TraceSV_odgi.sv_vcf_gz
        File introgressed_svs_strict_vcf = TraceSV_odgi.strict_vcf_gz

        # Intermediate file
        File node_to_paths = TraceSV_odgi.node_to_paths
    }
}

# ============================================================================
# Task Definitions
# ============================================================================

task PrepareReadsForPanGenie {
    input {
        File input_fastq
        String sample
        File ubuntu_sif
    }

    command <<<
        set -euo pipefail

        echo "[INFO] Preparing reads for PanGenie (requires uncompressed)"

        if [[ "~{input_fastq}" == *.gz ]]; then
            echo "[INFO] Decompressing .gz file..."
            zcat ~{input_fastq} > ~{sample}.fq
        else
            echo "[INFO] Copying uncompressed file..."
            cat ~{input_fastq} > ~{sample}.fq
        fi

        echo "[INFO] PanGenie input prepared"
        ls -lh ~{sample}.fq
    >>>

    output {
        File prepared_fastq = "~{sample}.fq"
    }

    runtime {
        req_cpu: 2
        req_memory: "8Gi"
        singularity_image: ubuntu_sif
    }
}

task KmerCountingKMC {
    input {
        File input_fastq
        Int kmer_length
        String sample
    }

    command <<<
        set -euo pipefail

        echo "[INFO] Running KMC with k=~{kmer_length}"
        echo ~{input_fastq} > file_list.txt
        kmc -k~{kmer_length} -m64 -okff -t16 @file_list.txt ~{sample} .
        rm file_list.txt

        echo "[INFO] KMC complete"
        ls -lh ~{sample}.kff
    >>>

    output {
        File kff_file = "~{sample}.kff"
    }

    runtime {
        req_cpu: 16
        req_memory: "64Gi"
    }
}

task HaplotypeSampling {
    input {
        File kff_file
        File gbz_file
        File hapl_file
        String sample
        File vg_sif
    }

    command <<<
        set -euo pipefail

        BASE=$(basename ~{gbz_file} .gbz)
        ln -s ~{gbz_file} ./$BASE.gbz
        ln -s ~{hapl_file} ./$BASE.hapl

        echo "[INFO] Running haplotype sampling"

        vg haplotypes \
            -v 2 \
            -i $BASE.hapl \
            -k ~{kff_file} \
            --include-reference \
            --set-reference "GRCh38" \
            --diploid-sampling \
            -g ~{sample}.gbz \
            $BASE.gbz

        echo "[INFO] Sampled graph created"
        ls -lh ~{sample}.gbz

        echo "[INFO] Verifying reference paths:"
        if vg paths -L -x ~{sample}.gbz | grep -q "GRCh38#0#"; then
            echo "[OK] Reference paths preserved"
        else
            echo "[WARN] No GRCh38#0# paths found"
            vg paths -L -x ~{sample}.gbz
        fi
    >>>

    output {
        File sampled_gbz = "~{sample}.gbz"
    }

    runtime {
        req_cpu: 4
        req_memory: "40Gi"
        singularity_image: vg_sif
    }
}

task Deconstruct {
    input {
        File sampled_gbz_file
        String sample
        String chromosome
        Int threads
        File vg_sif
        String singularity_bind_args
    }

    command <<<
        set -euo pipefail

        echo "[INFO] Running vg deconstruct"
        SINGULARITY_BIND_ARGS="~{singularity_bind_args}"
        singularity exec ${SINGULARITY_BIND_ARGS} ~{vg_sif} vg deconstruct \
            -L 0.75 \
            -a \
            -t ~{threads} \
            -v ~{sampled_gbz_file} \
            -P "GRCh38#0#" \
            > ~{sample}.gbz.raw.vcf

        sed 's/GRCh38#0#//g' ~{sample}.gbz.raw.vcf \
            | bgzip --threads ~{threads} \
            > ~{sample}.gbz.vcf.gz

        tabix -p vcf ~{sample}.gbz.vcf.gz

        # Keep only LV=0 records to avoid nested-overlap errors in PanGenie.
        bcftools view --threads ~{threads} \
            -r "~{chromosome}" \
            -s ^CHM13 \
            --force-samples \
            ~{sample}.gbz.vcf.gz \
            | bcftools filter -i 'INFO/LV=0' \
            | bcftools view --threads ~{threads} --min-ac 1 \
            | bcftools view --threads ~{threads} --trim-alt-alleles \
            > ~{sample}.pangenie_input.vcf

        echo "[INFO] Deconstruct complete"
        echo "Full VCF variants:"
        bcftools view -H ~{sample}.gbz.vcf.gz | wc -l
        echo "PanGenie input variants (LV=0 only):"
        grep -v "^#" ~{sample}.pangenie_input.vcf | wc -l || echo "0"
    >>>

    output {
        File full_vcf_gz = "~{sample}.gbz.vcf.gz"
        File full_vcf_tbi = "~{sample}.gbz.vcf.gz.tbi"
        File pangenie_input_vcf = "~{sample}.pangenie_input.vcf"
    }

    runtime {
        req_cpu: threads
        req_memory: "64Gi"
    }
}

task PanGenieIndex {
    input {
        File in_ref_file
        File in_ref_index_file
        File variants_vcf
        String pangenie_prefix
        Int threads
        File pangenie_sif
    }

    command <<<
        set -euo pipefail

        ln -s ~{in_ref_file} ref.fa
        ln -s ~{in_ref_index_file} ref.fa.fai

        PanGenie-index \
            -r ref.fa \
            -v ~{variants_vcf} \
            -o ~{pangenie_prefix}.pangenie_input \
            -t ~{threads}

        echo "[INFO] PanGenie index built"
        ls -lh ~{pangenie_prefix}.pangenie_input*
    >>>

    output {
        Array[File] graph_cereal_files = glob("~{pangenie_prefix}.pangenie_input_chr*_Graph.cereal")
        Array[File] kmers_tsv_gz_files = glob("~{pangenie_prefix}.pangenie_input_chr*_kmers.tsv.gz")
        File unique_kmers_map = "~{pangenie_prefix}.pangenie_input_UniqueKmersMap.cereal"
        File path_segments_fasta = "~{pangenie_prefix}.pangenie_input_path_segments.fasta"
    }

    runtime {
        req_cpu: threads
        req_memory: "40Gi"
        singularity_image: pangenie_sif
    }
}

task PanGenie {
    input {
        File reads
        String index_prefix
        Array[File] graph_cereal_files
        Array[File] kmers_tsv_gz_files
        File unique_kmers_map
        File path_segments_fasta
        String pangenie_prefix
        String sample
        Int threads
        File pangenie_sif
    }

    command <<<
        set -euo pipefail

        ln -sf ~{sep=" " graph_cereal_files} ./
        ln -sf ~{sep=" " kmers_tsv_gz_files} ./
        ln -sf ~{unique_kmers_map} ./
        ln -sf ~{path_segments_fasta} ./
        ln -sf ~{reads} ./

        echo "[INFO] Running PanGenie"

        PanGenie \
            -f ~{pangenie_prefix}.pangenie_input \
            -i ~{reads} \
            -s ~{sample} \
            -o ~{pangenie_prefix} \
            -t ~{threads} \
            -j ~{threads} \
            -u

        echo "[INFO] PanGenie complete"
        ls -lh ~{pangenie_prefix}_genotyping.vcf
    >>>

    output {
        File genotype_results = "~{pangenie_prefix}_genotyping.vcf"
    }

    runtime {
        req_cpu: threads
        req_memory: "40Gi"
        singularity_image: pangenie_sif
        continueOnReturnCode: [0, 1]
    }
}

task MergeVCF {
    input {
        File deconstruct_vcf
        File pangenie_vcf
        String sample
        File merge_script
    }

    command <<<
        set -euo pipefail

        echo "[INFO] Merging VCFs"

        cp ~{merge_script} ./merge_deconstruct_pangenie_vcf.sh
        chmod +x ./merge_deconstruct_pangenie_vcf.sh

        MERGE_OUT_DIR="merge_deconstruct_pangenie_vcf_out"
        bash ./merge_deconstruct_pangenie_vcf.sh \
            ~{deconstruct_vcf} \
            ~{pangenie_vcf} \
            ~{sample} \
            "${MERGE_OUT_DIR}"

        mv "${MERGE_OUT_DIR}"/~{sample}.merged.vcf.gz ./
        mv "${MERGE_OUT_DIR}"/~{sample}.merged.vcf.gz.tbi ./

        echo "[INFO] Merge complete"
        ls -lh ~{sample}.merged.vcf.gz*
    >>>

    output {
        File merged_vcf_gz = "~{sample}.merged.vcf.gz"
        File merged_vcf_tbi = "~{sample}.merged.vcf.gz.tbi"
    }

    runtime {
        req_cpu: 4
        req_memory: "16Gi"
    }
}

task ExtractIntrogressedPaths {
    input {
        File injected_og_file
        Array[String] path_patterns
        String sample
        File odgi_sif
    }

    command <<<
        set -euo pipefail

        echo "[INFO] Extracting introgressed paths"

        odgi paths \
            -i ~{injected_og_file} \
            -L \
            > all_paths.txt

        echo "[INFO] Total paths: $(wc -l < all_paths.txt)"

        PATTERN="~{sep='|' path_patterns}"
        grep -E "$PATTERN" all_paths.txt > ~{sample}.introgressed_paths.txt || true

        COUNT=$(wc -l < ~{sample}.introgressed_paths.txt)
        echo "[INFO] Introgressed paths found: $COUNT"

        if [ "$COUNT" -eq 0 ]; then
            echo "[ERROR] No introgressed paths found"
            head -20 all_paths.txt
            exit 1
        fi

        echo "First 10 introgressed paths:"
        head -10 ~{sample}.introgressed_paths.txt
    >>>

    output {
        File paths_file = "~{sample}.introgressed_paths.txt"
        File all_paths_file = "all_paths.txt"
    }

    runtime {
        req_cpu: 2
        req_memory: "8Gi"
        singularity_image: odgi_sif
    }
}

task TraceSV_odgi {
    input {
        File injected_og_file
        File merged_vcf
        File introgressed_paths_file
        String sample
        Int threads
        File odgi_sif
        File trace_script
        String singularity_bind_args
    }

    command <<<
        set -euo pipefail

        echo "[INFO] Tracing introgressed SVs"

        cp ~{trace_script} ./trace_introgressed_variants_odgi.sh
        chmod +x ./trace_introgressed_variants_odgi.sh

        export SINGULARITY_BIND_ARGS="~{singularity_bind_args}"
        bash -x ./trace_introgressed_variants_odgi.sh \
            ~{odgi_sif} \
            ~{injected_og_file} \
            ~{merged_vcf} \
            ~{introgressed_paths_file} \
            ~{sample} \
            ~{threads} || {
            echo "[ERROR] Script failed"
            ls -lha
            exit 1
        }

        echo "[INFO] SV tracing complete"
        ls -lh ~{sample}.introgressed_svs*.vcf.gz
        ls -lh node_to_paths.txt
    >>>

    output {
        File base_vcf_gz = "~{sample}.introgressed_svs.vcf.gz"
        File base_vcf_tbi = "~{sample}.introgressed_svs.vcf.gz.tbi"
        File sv_vcf_gz = "~{sample}.introgressed_svs.SV.vcf.gz"
        File sv_vcf_tbi = "~{sample}.introgressed_svs.SV.vcf.gz.tbi"
        File strict_vcf_gz = "~{sample}.introgressed_svs.SVstrict.vcf.gz"
        File strict_vcf_tbi = "~{sample}.introgressed_svs.SVstrict.vcf.gz.tbi"
        File node_to_paths = "node_to_paths.txt"
    }

    runtime {
        req_cpu: threads
        req_memory: "32Gi"
    }
}
