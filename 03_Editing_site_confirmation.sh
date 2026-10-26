#!/usr/bin/env bash
# RNA editing site selection and visualization from snpEff-annotated VCFs.

# =============================================================================
# 00. Paths, sample metadata and software
# =============================================================================
set -euo pipefail
export LC_ALL=C
# Locate bundled scripts relative to this entrypoint; input paths remain configurable.
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
SCRIPT_DIR="${PIPELINE_DIR}/code/03"
export INPUT_ROOT="${INPUT_ROOT:-${PIPELINE_DIR}/../04_reproduce_original_20260913/inputs}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-${PIPELINE_DIR}/results}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${SCRIPT_DIR}/snpeff_samples.tsv}"
[[ -s "${SAMPLE_TABLE}" ]] || {
    echo "Missing sample table: ${SAMPLE_TABLE}" >&2
    exit 1
}
G_ALIGN="${INPUT_ROOT}/new_snpeff/g_align"
CHAIN_DIR="${CHAIN_DIR:-/path/to/Proj1_RNAediting/02_Bcftools_SnpEff/g_align}"
LIFTOVER="${LIFTOVER:-${CHAIN_DIR}/liftOver}"
RSCRIPT="${RSCRIPT:-/path/to/.local/share/mamba/envs/r_env/bin/Rscript}"
PYTHON="${PYTHON:-python3}"
PLOT_RSCRIPT="${PLOT_RSCRIPT:-Rscript}"
PLOT_R_LIBS="${PLOT_R_LIBS:-${PIPELINE_DIR}/../04_reproduce_original_20260913/Rlib}"
OUT_DIR="${OUT_DIR:-${PIPELINE_DIR}/results/03_editing_sites_test2_union}"
WORKERS="${WORKERS:-8}"
N_UPSET="${N_UPSET:-40}"
ISO1_MODE=all
SPECIES_MODE=test2
FIGURE_DIR="${OUT_DIR}/figure1"
FIGURE_DATA_DIR="${FIGURE_DIR}/data"
RECLASSIFY_SEC_UTR="${RECLASSIFY_SEC_UTR:-1}"

# Input: code/03/snpeff_samples.tsv (tab-separated, with a header).
# First rows:
# library	group	annotation_relative_path	table_prefix	chain	chromosome_substitution
# Adar5G1_R1	Adar5G1	Dmel_Adar5G1_rep1_mRNA/05_variants/Dmel_Adar5G1_rep1_mRNA.annotated.vcf	old_head/Dmel_Adar5G1_rep1_mRNA.uniq1.vcf.Filter.ann	-	-
# Adar5G1_R2	Adar5G1	Dmel_Adar5G1_rep2_mRNA/05_variants/Dmel_Adar5G1_rep2_mRNA.annotated.vcf	old_head/Dmel_Adar5G1_rep2_mRNA.uniq1.vcf.Filter.ann	-	-
# WT_F1	WT	WT-F1-LJD7606/05_variants/WT-F1-LJD7606.annotated.vcf	new_snpeff/WT-F1-LJD7606.uniq1.vcf.Filter.ann	-	-
# Read annotation_relative_path under ANNOTATION_ROOT; write table_prefix under INPUT_ROOT.

# =============================================================================
# 01. Parse snpEff ANN records into 11-column test2 tables
# =============================================================================
while IFS=$'\t' read -r library group annotation_relative_path table_prefix chain chromosome_substitution; do
    source_vcf="${ANNOTATION_ROOT}/${annotation_relative_path}"
    ann="${INPUT_ROOT}/${table_prefix}"
    [[ -f "${source_vcf}" ]] || {
        echo "Missing annotation: ${source_vcf}" >&2
        exit 1
    }
    mkdir -p "$(dirname "${ann}")"
    echo "[01 parse] ${library}: ${source_vcf}"
    grep 'ANN=' "${source_vcf}" | perl "${SCRIPT_DIR}/VCFAnnParse_step1.pl" > "${ann}.test"
    perl "${SCRIPT_DIR}/VCFAnnParse_step2.pl" "${ann}.test" |
        awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' \
            > "${ann}.test2"
done < <(tail -n +2 "${SAMPLE_TABLE}")

# =============================================================================
# 02. Filter ISO1 sites and extract transcript-oriented A>G records
# =============================================================================
# ISO1 uses at2-filtered tables; others use unfiltered test2 records.
# sim/sec/ere use test2_A2G without additional DP/MQ filtering.
while IFS=$'\t' read -r library group annotation_relative_path table_prefix chain chromosome_substitution; do
    ann="${INPUT_ROOT}/${table_prefix}"
    echo "[02 arrange] ${library}"
    "${RSCRIPT}" "${SCRIPT_DIR}/arrange_step3.R" "${ann}.test2"
    awk -F '\t' '$10 ~ /A>G/' "${ann}.test2" > "${ann}.test2_A2G"
done < <(tail -n +2 "${SAMPLE_TABLE}")

# =============================================================================
# 03. Lift sim/sec/ere coordinates to dm6 before candidate selection
# =============================================================================
# Process the 22 libraries with a chain entry; Dmel coordinates are unchanged.
# Convert 1-based POS to BED [POS-1, POS), then liftOver and export .sites.
# Payload columns 12/16 are empty for test2; the selector recalculates read depth.
mkdir -p "${G_ALIGN}"
while IFS=$'\t' read -r library group annotation_relative_path table_prefix chain chromosome_substitution; do
    input="${INPUT_ROOT}/${table_prefix}.test2_A2G"
    output_name="$(basename "${table_prefix}").test2_A2G"
    bed="${G_ALIGN}/${output_name}.bed"
    mapped="${G_ALIGN}/${output_name}.${chain}"
    echo "[03 liftover] ${library}: ${chain}"
    awk -F '\t' '{print $1"\t"$2-1"\t"$2"\t"$1","$2","$3","$4","$5","$6","$7","$8","$9","$11","$12","$16}' "${input}" |
        grep -v '^chr' | sed "${chromosome_substitution}" > "${bed}"
    "${LIFTOVER}" "${bed}" "${CHAIN_DIR}/${chain}.over.chain" "${mapped}" "${mapped}.unMapped"
    sed 's/chr//1' "${mapped}" | sed 's/,/\t/g' > "${mapped}.sites"
done < <(awk -F '\t' 'NR > 1 && $5 != "-"' "${SAMPLE_TABLE}")

# =============================================================================
# 04. Select candidate sites and export editing levels and UpSet plots
# =============================================================================
# Candidate = ISO1_union intersect (WBG_min2 union Published_union), after
# Adar background subtraction applied uniformly to all sites.
# sim/sec/ere use library unions after one-to-one mapping checks; WT require
# support from at least two libraries. The three-source UpSet covers candidates.
mkdir -p "${OUT_DIR}"
R_LIBS_USER="${PLOT_R_LIBS}" PLOT_RSCRIPT="${PLOT_RSCRIPT}" \
    "${PYTHON}" "${SCRIPT_DIR}/rebuild_candidate_corrected_WBG_change2_species_union.py" \
    "${OUT_DIR}" "${WORKERS}" "${N_UPSET}" "${ISO1_MODE}" "${SPECIES_MODE}" \
    > "${OUT_DIR}/candidate.run.log" 2>&1

# =============================================================================
# 05. Prepare plotting inputs from final result tables
# =============================================================================
mkdir -p "${FIGURE_DIR}" "${FIGURE_DATA_DIR}"
R_LIBS_USER="${PLOT_R_LIBS}" "${PLOT_RSCRIPT}" "${SCRIPT_DIR}/prepare_figure1_inputs.R" \
    "${OUT_DIR}" "${FIGURE_DATA_DIR}" > "${FIGURE_DIR}/prepare_inputs.log" 2>&1

# =============================================================================
# 06. Generate Figure 1 panels and All site types
# =============================================================================
# RECLASSIFY_SEC_UTR=1 assigns Dsec Others to homologous ISO1 UTR classes for display.
# Heatmaps display missing values and the Adar row as zero. MYA values are schematic.
mkdir -p "${FIGURE_DIR}/panels"
R_LIBS_USER="${PLOT_R_LIBS}" "${PLOT_RSCRIPT}" "${SCRIPT_DIR}/plot_figure1_panels.R" \
    "${FIGURE_DATA_DIR}" "${FIGURE_DIR}/panels" "${RECLASSIFY_SEC_UTR}" \
    > "${FIGURE_DIR}/figure1.run.log" 2>&1

# =============================================================================
# 07. Optional: redraw UpSet plots from existing candidate tables
# =============================================================================
# R_LIBS_USER="${PLOT_R_LIBS}" PLOT_RSCRIPT="${PLOT_RSCRIPT}" \
#     "${PYTHON}" "${SCRIPT_DIR}/plot_candidate_upsets.py" "${OUT_DIR}" "${N_UPSET}" \
#     > "${OUT_DIR}/upset.run.log" 2>&1

# =============================================================================
# 08. Figure S1: sequence motifs around editing sites
# =============================================================================
# Configure the candidate, one-to-one species, quantitative-table and genome inputs.
# The original analysis exports both kept and min2 branches; S1 uses current_1v1_min2.
export FS1_CANDIDATE_DIR="${FS1_CANDIDATE_DIR:-<path/to/candidate_results>}"
export FS1_SPECIES_DIR="${FS1_SPECIES_DIR:-<path/to/species_1v1_results>}"
export FS1_PLOT_DATA_DIR="${FS1_PLOT_DATA_DIR:-<path/to/editing_level_table_directory>}"
export REF_DIR="${REF_DIR:-<path/to/species_genome_fasta_directory>}"
FS1_OUTPUT_DIR="${FS1_OUTPUT_DIR:-${PIPELINE_DIR}/results/03/FS1}"

"${PLOT_RSCRIPT:-Rscript}" "${SCRIPT_DIR}/FS1/plot_ADAR_motif_1v1.R" \
    "${FS1_OUTPUT_DIR}" "${MIN_DEPTH:-10}" "${MIN_EDITING:-0.1}" \
    "${PLOT_MODE:-both}" "${FLANK:-1}"
