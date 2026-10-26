#!/usr/bin/env bash
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

# Mass-spectrometry analysis workflow for detecting and quantifying
# RNA-editing-derived peptide products.
#
# This file is organized as independent modules. Run one section at a time
# after replacing values enclosed in <...>. It is not intended to be launched
# from top to bottom as a single batch job.
#
# Dependency order:
#   00-01  Configure paths and validate site/reference inputs.
#   02-03  Build the two edited/unedited sequence databases (parallel strategies).
#   04     Optional DDA/MaxQuant search using either database.
#   05-06  Predict DIA libraries, then search the DIA raw files.
#   07-09  Quantify peptide states; summarize multi-site states and site-level Lp.
#   10-12  Calculate Lm/Lr, join with Lp, and plot the four-species results.
#   13     Figure S19 hierarchical upstream/downstream peptide analysis.
#
# Sections 02 and 03 are alternatives to compare, not inputs to each other.
# Section 04 is a DDA branch; its output is not used by the DIA-NN commands below.
# Section 08 describes a separate multi-site summary and is not required for Lp.
# Section 08B reproduces the RNA/RPF Figure S9 independently of MS.
# Section 10 can run independently once the RNA/Ribo input tables are ready.
# Section 13 reproduces Figure S19 from retained full-proteome peptide x sample tables.

# =============================================================================
# 00. Project paths and species-specific inputs
# =============================================================================

PROJECT_DIR="${PROJECT_DIR:-<path/to/project>}"
RESULT_DIR="${PROJECT_DIR}/03_MS"
ANALYSIS_DIR="${RESULT_DIR}/DIA_NN_analysis"
SCRIPT_DIR="${PIPELINE_DIR}/code/06"
REF_DIR="${REF_DIR:-<path/to/reference_genomes>}"

# Select one species for each run: dmel, dsim, dsec, or dere.
SPECIES="${SPECIES:-<species>}"
GROUP_NAME="${GROUP_NAME:-<ISO1|sim|sec|ere>}"
RELEASE="${RELEASE:-<r6.12|r2.02|r1.3|r1.05>}"

FBTR_FBPP_MAP="${FBTR_FBPP_MAP:-<path/to/${SPECIES}_${RELEASE}_fbtr_fbpp.tsv>}"
CDS_FASTA="${CDS_FASTA:-${REF_DIR}/<species-CDS.fasta>}"
TRYPTIC_PEPTIDES="${TRYPTIC_PEPTIDES:-<path/to/${SPECIES}_${RELEASE}_tryptic_peptides.tsv>}"
MISSENSE_POS="${MISSENSE_POS:-<path/to/${GROUP_NAME}.missense.pos>}"
SOURCE_ANNOTATIONS="${SOURCE_ANNOTATIONS:-<path/to/${GROUP_NAME}.candidate.annotations.long.tsv>}"

mkdir -p "${RESULT_DIR}" "${ANALYSIS_DIR}"

# Reference configurations used in this study:
#   dmel  ISO1  dmel-all-CDS-r6.12.fasta
#   dsim  sim   dsim-all-CDS-r2.02.fasta
#   dsec  sec   dsec-all-CDS-r1.3.fasta
#   dere  ere   dere-all-CDS-r1.05.fasta

# =============================================================================
# 01. Prepare editing sites and verify reference coordinates
# =============================================================================

# No conversion step is included in this workflow. Prepare MISSENSE_POS as a
# tab-delimited file with the following header and one editing site per row:
#
# site_id       gene     FBgn         FBtr         refAA  altAA  trans_pos  trans_len  cds_pos  cds_len  aa_pos  aa_len
# 2L_11159274   Ca-beta  FBgn0259822  FBtr0310279  Tyr    Cys    3618       4514       2531     2535     844     844
# 2L_11159893   Ca-beta  FBgn0259822  FBtr0310279  Met    Val    3149       4514       2062     2535     688     844
#
# Coordinates in trans_pos, cds_pos, and aa_pos are 1-based. site_id identifies
# the genomic editing site; refAA and altAA denote the unedited and edited amino
# acids, respectively. The peptide-database script requires site_id, gene, FBgn,
# FBtr, refAA, altAA, and aa_pos. The 100-aa workflow additionally uses cds_pos.

# Before database construction, ensure FBTR_FBPP_MAP, CDS_FASTA,
# TRYPTIC_PEPTIDES and the amino-acid coordinates use the same reference release
# and the same transcript/protein pairing. aa_pos must refer to the protein
# used to build the search database.
#
# A genome-derived CDS can differ from the curated protein because of a
# reference-sequence error, a strain indel or a translation exception. Such
# records need a validated CDS-to-protein codon map before amino-acid editing
# consequences are used in either database strategy. --strict-ref-check checks
# the edited site's reference residue; it does not validate the entire frame.

# =============================================================================
# 02. Build edited and unedited tryptic peptide databases
# =============================================================================

TRYPTIC_OUT="${RESULT_DIR}/${SPECIES}_${RELEASE}_edit_peptides"

# Generate reference and edited tryptic peptides for sites that produce
# missense substitutions. --strict-ref-check requires the reported reference
# amino acid to agree with the reference peptide sequence.
python3 "${SCRIPT_DIR}/make_edit_peptide_fasta_v2.py" \
    --map "${FBTR_FBPP_MAP}" \
    --edits "${MISSENSE_POS}" \
    --peptides "${TRYPTIC_PEPTIDES}" \
    --out-prefix "${TRYPTIC_OUT}" \
    --strict-ref-check

# Generate peptide sequences containing multiple editable codons when more
# than one candidate editing site occurs in the same tryptic peptide.
python3 "${SCRIPT_DIR}/make_multi_codon_edit_peptide_fasta.py" \
    --map "${FBTR_FBPP_MAP}" \
    --source-annotations "${SOURCE_ANNOTATIONS}" \
    --peptides "${TRYPTIC_PEPTIDES}" \
    --out-prefix "${TRYPTIC_OUT}.multi_codon"

# Combine reference, single-edit, and multi-edit peptide sequences into the
# FASTA database used for the tryptic-peptide search.
cat \
    "${TRYPTIC_OUT}.wt.fasta" \
    "${TRYPTIC_OUT}.edit.fasta" \
    "${TRYPTIC_OUT}.multi_codon.edit.fasta" \
    > "${TRYPTIC_OUT}.edit_wt_multi.fasta"

# =============================================================================
# 03. Build edited and unedited 100-aa sequence databases
# =============================================================================

AA100_PREFIX="${RESULT_DIR}/${SPECIES}_candidate_missense"

# Introduce all candidate missense edits into their source CDS sequences and
# translate the edited CDS. For D. sechellia, use 1.edit_more_site_Dsec.py.
python3 "${SCRIPT_DIR}/1.edit_more_site.py" \
    -i "${MISSENSE_POS}" \
    -f "${CDS_FASTA}" \
    -o "${AA100_PREFIX}.edited.CDS.fasta" \
    -p "${AA100_PREFIX}.edited.protein.fasta"

# Extract codon-aligned windows with 150 nt of flanking sequence on each side
# of the edited-site cluster (the nominal "100-aa" database design).
# Windows longer than 600 nt are skipped, not truncated to 600 nt.
python3 "${SCRIPT_DIR}/2.get_100peptide.py" \
    -i "${AA100_PREFIX}.edited.CDS.fasta" \
    -o "${AA100_PREFIX}.edited.CDS.100AA.fasta" \
    --flank 150 \
    --max_len 600

# Reconstruct the matching unedited sequence from the reference CDS.
# For D. sechellia, use 3.bp_unedit_Dsec.py.
python3 "${SCRIPT_DIR}/3.bp_unedit.py" \
    -w "${AA100_PREFIX}.edited.CDS.100AA.fasta" \
    -r "${CDS_FASTA}" \
    -o "${AA100_PREFIX}.unedited.CDS.100AA.fasta"

cat \
    "${AA100_PREFIX}.edited.CDS.100AA.fasta" \
    "${AA100_PREFIX}.unedited.CDS.100AA.fasta" \
    > "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.fasta"

# Pad the sequence end with N when its length is not divisible by three,
# translate, then standardize FASTA identifiers for downstream DIA-NN matching.
# Padding only changes the terminal length; it does not repair an internal
# frameshift or establish the correct protein reading frame.
python3 "${SCRIPT_DIR}/4.check_3.py" \
    "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.fasta" \
    "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.pad.fasta"

python3 "${SCRIPT_DIR}/4.cds_translate.py" \
    "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.pad.fasta" \
    "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.protein.fasta"

python3 "${SCRIPT_DIR}/rename_fasta_header_old_style_patched.py" \
    -i "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.protein.fasta" \
    -o "${AA100_PREFIX}.edit_and_unedited.CDS.100AA.search.fasta"

# =============================================================================
# 04. Optional DDA branch: search with MaxQuant
# =============================================================================

# This optional DDA branch starts after database construction. Configure the
# FASTA and raw-file paths in MAXQUANT_XML; no MaxQuant output is consumed by
# the DIA-NN branch below.
# Run each MaxQuant parameter file in a separate working directory because
# MaxQuant writes a fixed combined/ output directory.
MAXQUANT_DLL="${MAXQUANT_DLL:-<path/to/MaxQuantCmd.dll>}"
# Example: code/06/config/mqpar.example.xml (Dsim 100-aa, MaxQuant 2.7.5.0).
# Copy and edit all input/output paths before use; see code/06/config/README.txt.
MAXQUANT_XML="${MAXQUANT_XML:-<path/to/edited_mqpar.xml>}"
MAXQUANT_RUN_DIR="${MAXQUANT_RUN_DIR:-${RESULT_DIR}/MaxQuant/<analysis_name>}"

mkdir -p "${MAXQUANT_RUN_DIR}"
(
    cd "${MAXQUANT_RUN_DIR}"
    dotnet "${MAXQUANT_DLL}" "${MAXQUANT_XML}" \
        > "${MAXQUANT_RUN_DIR}/maxquant.log" 2>&1
)

# The tryptic and 100-aa databases were searched separately. The MaxQuant
# version and search parameters must be reported from the corresponding XML.

# =============================================================================
# 05. Build predicted spectral libraries for DIA-NN
# =============================================================================

DIA_RUN_ROOT="${RESULT_DIR}/DIA_NN"
TRYPTIC_FASTA="${TRYPTIC_OUT}.edit_wt_multi.fasta"
AA100_FASTA="${AA100_PREFIX}.edit_and_unedited.CDS.100AA.search.fasta"

mkdir -p "${DIA_RUN_ROOT}/trypsin" "${DIA_RUN_ROOT}/100AA"

# Pass the exact FASTAs generated above, rather than relying on a host-specific glob.
case "$SPECIES" in
    dmel) DIA_SPECIES_DIR=Dmel ;;
    dsim) DIA_SPECIES_DIR=Dsim ;;
    dsec) DIA_SPECIES_DIR=Dsec ;;
    dere) DIA_SPECIES_DIR=Dere ;;
    *)
        echo "SPECIES must be dmel, dsim, dsec or dere" >&2
        exit 2
        ;;
esac
FASTA_FILE="$TRYPTIC_FASTA" OUT_ROOT="$DIA_RUN_ROOT/trypsin" \
    bash "${SCRIPT_DIR}/01_build_diann_predicted_libraries_multi_species_trypsin.sh" "$SPECIES"
FASTA_FILE="$AA100_FASTA" OUT_ROOT="$DIA_RUN_ROOT/100AA" \
    bash "${SCRIPT_DIR}/01_build_diann_predicted_libraries_multi_species_100AA.sh" "$SPECIES"

# =============================================================================
# 06. Process DIA raw data with DIA-NN
# =============================================================================

ANALYSIS_NAME="${ANALYSIS_NAME:-<analysis_name>}"
RAW_PREFIX="${RAW_PREFIX:-<raw-file-prefix>}"
TRYPTIC_LIBRARY="${TRYPTIC_LIBRARY:-${DIA_RUN_ROOT}/trypsin/${DIA_SPECIES_DIR}/$(basename "${TRYPTIC_FASTA}" .fasta).predicted.speclib}"
AA100_LIBRARY="${AA100_LIBRARY:-${DIA_RUN_ROOT}/100AA/${DIA_SPECIES_DIR}/$(basename "${AA100_FASTA}" .fasta).predicted.speclib}"

# Search against the edited/unedited tryptic-peptide database.
OUT_ROOT_BASE="${DIA_RUN_ROOT}/trypsin" \
    bash "${SCRIPT_DIR}/02_run_diann_astral_multi_species_trypsin.sh" \
    "${SPECIES}" \
    "${ANALYSIS_NAME}" \
    "${TRYPTIC_LIBRARY}" \
    "${TRYPTIC_FASTA}" \
    "${RAW_PREFIX}" \
    > "${DIA_RUN_ROOT}/trypsin/${ANALYSIS_NAME}.log" 2>&1

# Search against the edited/unedited 100-aa sequence database.
OUT_ROOT_BASE="${DIA_RUN_ROOT}/100AA" \
    bash "${SCRIPT_DIR}/02_run_diann_astral_multi_species_100AA.sh" \
    "${SPECIES}" \
    "${ANALYSIS_NAME}" \
    "${AA100_LIBRARY}" \
    "${AA100_FASTA}" \
    "${RAW_PREFIX}" \
    > "${DIA_RUN_ROOT}/100AA/${ANALYSIS_NAME}.log" 2>&1

# =============================================================================
# 07. Quantify edited and unedited peptide states at each site
# =============================================================================

mkdir -p "${ANALYSIS_DIR}"

# Use the exact FASTA that produced each corresponding DIA-NN report.
# Repeat for each species before the cross-species plots in Section 12.
# Quantification using the 100-aa search results.
python3 "${SCRIPT_DIR}/DIA_NN_proteoform_state_quant_v2.py" \
    --sites "${MISSENSE_POS}" \
    --matrix "<path/to/100AA/report.pr_matrix.tsv>" \
    --fasta "${AA100_FASTA}" \
    --out "${ANALYSIS_DIR}/${SPECIES}.site_intensity.100AA.tsv" \
    --detail-out "${ANALYSIS_DIR}/${SPECIES}.site_intensity.100AA.detail.tsv"

# Quantification using the in-silico tryptic-peptide search results.
python3 "${SCRIPT_DIR}/DIA_NN_proteoform_state_quant_FBpp_v4.py" \
    --sites "${MISSENSE_POS}" \
    --matrix "<path/to/trypsin/report.pr_matrix.tsv>" \
    --fasta "${TRYPTIC_FASTA}" \
    --id-map "${FBTR_FBPP_MAP}" \
    --out "${ANALYSIS_DIR}/${SPECIES}.site_intensity.trypsin.tsv" \
    --detail-out "${ANALYSIS_DIR}/${SPECIES}.site_intensity.trypsin.detail.tsv"

# =============================================================================
# 08. Optional multi-site peptide-state summary
# =============================================================================

# This optional branch uses peptide/allele-level intensities from the 100-aa
# search and quantification in Sections 06-07. It is separate from the
# single-site Lp summary in Section 09.
#
# Run this section independently after configuring the Table S10 inputs.
# See code/06/table_s10_multiedit/README.txt for evidence provenance and definitions.
# MS detected means any IA/AI/II peptide state, not II alone. Updated Dere
# evidence includes AA and IA; an unlisted-alt whole-site filter omits it.
#
# For peptides containing multiple editable codons, allele states (e.g.
# AA, IA, AI, and II) are assigned within each sample. Intensities are first
# summed within each site x genotype x sex x replicate x allele combination.
# Replicate-level proportions are then calculated as:
#
#   allele proportion = allele intensity / sum of all allele intensities
#
# Mean output values are obtained across biological replicates within each
# site x genotype x sex x allele group. This step is performed with
# summarize_allele_intensity_proportion_100AA.py using the desired DIA-NN
# intensity column (for example, samplesum_intensity).

bash "${SCRIPT_DIR}/run_table_s10.sh"

# =============================================================================
# 08B. Figure S9: equal-replicate pooled RNA/RPF editing-state fractions
# =============================================================================

# RNA/RPF total >=5 per paired sample; at the codon level, IA+AI+II must
# be observed in each assay across qualifying samples (not necessarily the
# same sample). Rescale each assay in each sample to total=1, pool rescaled
# counts, calculate state fractions, then log2(pRPF/pRNA). No pseudocounts;
# preserve -infinity, +infinity, and NA. The publication layout has 84 codons.
# FS9_COUNTS_DIR supplies double_edit_codons.tsv and
# double_edit_allele_counts.by_sample.tsv from the upstream candidate counter.
# Run this section independently after configuring FS9_COUNTS_DIR.
bash "${SCRIPT_DIR}/run_figure_s9.sh"

# =============================================================================
# 08C. Figure 4B / S14: final DIA-NN fractions, Dmel WT = WBg only
# =============================================================================
bash "${PIPELINE_DIR}/code/06/figure4_s14_diann/run.sh"

# =============================================================================
# 09. Compare database strategies and calculate mean Lp
# =============================================================================

# These site-intensity tables come directly from Section 07. The optional
# multi-site allele summary in Section 08 is not an input to this calculation.
# First assess the overlap between sites detected by the 100-aa and tryptic
# searches.
python3 "${SCRIPT_DIR}/compare_100AA_trypsin_site_overlap.py" \
    --aa100 "${ANALYSIS_DIR}/${SPECIES}.site_intensity.100AA.tsv" \
    --trypsin "${ANALYSIS_DIR}/${SPECIES}.site_intensity.trypsin.tsv" \
    --out-prefix "${ANALYSIS_DIR}/${SPECIES}.100AA_vs_trypsin"

# Combine site-level measurements hierarchically across genotype, sex, and
# biological replicates. Change --intensity-col to median_intensity or
# sum_intensity only for the corresponding sensitivity analysis.
python3 "${SCRIPT_DIR}/calculate_Lp_hierarchical_genotype_sex_mean_v4.py" \
    --aa100 "${ANALYSIS_DIR}/${SPECIES}.site_intensity.100AA.tsv" \
    --trypsin "${ANALYSIS_DIR}/${SPECIES}.site_intensity.trypsin.tsv" \
    --out-prefix "${ANALYSIS_DIR}/${SPECIES}.Lp_union"

# =============================================================================
# 10. Calculate RNA/Ribo mean Lm and Lr
# =============================================================================

# This all-species RNA/Ribo summary is independent of MS database generation.
# Its output is first consumed by the Lp/Lm/Lr join in Section 11.
# Recalculate the site-level mRNA (Lm) and ribosome-profiling (Lr) measurements.
# A site/sample measurement is retained when either mRNA coverage or RPF A-site
# coverage is at least three reads.
python3 "${SCRIPT_DIR}/recalculate_Lm_Lr_min3_or_coverage_v2.py" \
    --input-dir "<path/to/all-editing-site-results>" \
    --out-prefix "${ANALYSIS_DIR}/all_species_LmLr_min3" \
    --min-coverage 3

# =============================================================================
# 11. Match peptide-level Lp with RNA/Ribo Lm and Lr
# =============================================================================

# Inputs: species-specific Lp from Section 09 and the all-species Lm/Lr
# summary from Section 10.
python3 "${SCRIPT_DIR}/match_Lp_with_min3_LmLr_v5.py" \
    --lp-wide "${ANALYSIS_DIR}/${SPECIES}.Lp_union.combined_sex_Lp_wide.tsv" \
    --lmlr "${ANALYSIS_DIR}/all_species_LmLr_min3.sex_site_mean_LmLr.tsv" \
    --cdspos-dir "<path/to/species-CDS-position-tables>" \
    --out "${ANALYSIS_DIR}/${SPECIES}.Lp_union.min3_LmLr.tsv"

# =============================================================================
# 12. Plot cross-layer editing measurements
# =============================================================================

# Plot Lp versus Lm and Lp versus Lr for the four species. Supply the four
# species-specific tables generated in Section 11.
Rscript "${SCRIPT_DIR}/plot_male_Lp_vs_Lm_Lr_longestTx.R" \
    "<dmel.Lp_LmLr.tsv>" \
    "<dsim.Lp_LmLr.tsv>" \
    "<dsec.Lp_LmLr.tsv>" \
    "<dere.Lp_LmLr.tsv>" \
    "${ANALYSIS_DIR}/Lp_Lm_Lr_male_figures"

# Compare the distributions of protein- and mRNA-level editing measurements.
Rscript "${SCRIPT_DIR}/plot_Ipos_Lp_Lm_3sex_2methods_A4_v4.R" \
    "${ANALYSIS_DIR}/Ipos_Lp_Lm_3sex" \
    "<path/to/species-CDS-position-tables>" \
    "<dmel.Lp_LmLr.tsv>" \
    "<dsim.Lp_LmLr.tsv>" \
    "<dsec.Lp_LmLr.tsv>" \
    "<dere.Lp_LmLr.tsv>"

# =============================================================================
# 12B. Final counter-rule Lp/Lm/Lr figures (independent consolidated workflow)
# =============================================================================
# Configure SOURCE_PROJECT, MS_ROOT and CDSPOS_DIR for this tested workflow.
# This is an alternative consolidated entry to Sections 09-12, not an extra filter.
export OUTPUT_ROOT="${OUTPUT_ROOT:-${PIPELINE_DIR}/results/06/Lp_Lm_Lr}"
bash "${PIPELINE_DIR}/code/06/run_Lp_Lm_Lr_counter_rule_figures.sh"

# =============================================================================
# 13. Figure S19: peptide abundance upstream/downstream of editing sites
# =============================================================================

# Independent branch: full-protein library -> full-protein DIA search ->
# editing-site versus random-site upstream/downstream peptide abundance.
# It does not consume the site-level Lp/Lm/Lr tables from Sections 09-12.
# Run this module on the MS host, where /data/dongxn/MS is available.
# If needed, deactivate conda in the interactive shell before running 13A/13B:
#   conda deactivate
# Select a Python environment containing pandas, numpy, matplotlib and scipy.

UPDOWN_MS_DIR="${UPDOWN_MS_DIR:-/data/dongxn/MS}"
# Local bundled code; raw files, FASTAs and outputs remain on the configured MS host.
UPDOWN_CODE_DIR="${PIPELINE_DIR}/code/06"
# The full-protein library builder uses longest-CDS protein FASTAs.
FULLPROTEIN_LIBRARY_SCRIPT="${FULLPROTEIN_LIBRARY_SCRIPT:-${UPDOWN_CODE_DIR}/01_build_diann_predicted_libraries_multi_species_old.sh}"
# Keep new full-protein runs separate: some historical DIA_NN_RUN reports were reused for 100-aa searches.
UPDOWN_RUN_ROOT="${UPDOWN_RUN_ROOT:-${UPDOWN_MS_DIR}/DIA_NN_RUN_full_protein}"
UPDOWN_RAW_DIR="${UPDOWN_RAW_DIR:-${UPDOWN_MS_DIR}/Rawdata_Astral}"
UPDOWN_REF_DIR="${UPDOWN_REF_DIR:-/path/to/RNA_editing_DATAprocessing_new/ref_genome}"
UPDOWN_PYTHON="${UPDOWN_PYTHON:-python}"

# Wait for every job in a stage and propagate any failure before the next stage.
# A bare 'wait' would not reliably report failures from every background job.
ms_wait_for_jobs() {
    local ms_pid ms_failed=0
    for ms_pid in "$@"; do
        if wait "${ms_pid}"; then
            :
        else
            echo "MS job failed (PID ${ms_pid}); check its log." >&2
            ms_failed=1
        fi
    done
    return "${ms_failed}"
}

# -----------------------------------------------------------------------------
# 13A. Optional reconstruction of the historical full-protein DIA-NN inputs
#
# The published Figure S19 uses retained deduplicated peptide x sample tables
# under DIA_NN_RUN. Set FS19_REBUILD_DIANN=1 only to reconstruct upstream
# search products for provenance checks; the default Figure S19 run skips it.
#
if [[ "${FS19_REBUILD_DIANN:-0}" == "1" ]]; then
    mkdir -p "${UPDOWN_RUN_ROOT}/logs"

    # The species-specific longest-CDS-per-gene protein FASTAs used in 13B are
    # pre-existing inputs; their generation is not included in this section.
    # These libraries contain the longest-CDS protein per gene and provide the
    # reference-protein searches used for the upstream/downstream comparison.
    # All four species run concurrently; finish this stage before starting 13B.
    # -----------------------------------------------------------------------------

    ms_library_pids=()
    for ms_species in dmel dsim dsec dere; do
        FASTA_DIR="$UPDOWN_MS_DIR" OUT_ROOT="$UPDOWN_RUN_ROOT" \
            bash "${FULLPROTEIN_LIBRARY_SCRIPT}" \
            "${ms_species}" \
            > "${UPDOWN_RUN_ROOT}/logs/build_diann_${ms_species}.log" 2>&1 &
        ms_library_pids+=("$!")
    done
    ms_wait_for_jobs "${ms_library_pids[@]}"

    # -----------------------------------------------------------------------------
    # 13B. Process Astral RAW files against the full-protein libraries
    #
    # Use the stable .predicted.speclib alias supplied by the builder.
    # The actual DIA-NN 2.6.1 file retains its .predicted.predicted.speclib name.
    # Dmel combines WBg, ISO1, AS and AG RAW prefixes in the same analysis.
    # Finish all four searches before reading their report.pr_matrix.tsv in 13C.
    # -----------------------------------------------------------------------------

    ms_search_pids=()

    RAW_DIR="$UPDOWN_RAW_DIR" OUT_ROOT_BASE="$UPDOWN_RUN_ROOT" \
        bash "${UPDOWN_CODE_DIR}/02_run_diann_astral_multi_species.sh" \
        dmel \
        WBg_ag_as_ISO1 \
        "${UPDOWN_RUN_ROOT}/Dmel/dmel-all-CDS-r6.12.longest_CDS_per_gene.protein.predicted.speclib" \
        "${UPDOWN_MS_DIR}/dmel-all-CDS-r6.12.longest_CDS_per_gene.protein.fasta" \
        XQH_20260320_WBg_ \
        DXN_20260401_iso_ \
        XQH_20260320_AS_ \
        XQH_20260320_AG_ \
        > "${UPDOWN_RUN_ROOT}/logs/run_diann_WBg_ISO1.log" 2>&1 &
    ms_search_pids+=("$!")

    RAW_DIR="$UPDOWN_RAW_DIR" OUT_ROOT_BASE="$UPDOWN_RUN_ROOT" \
        bash "${UPDOWN_CODE_DIR}/02_run_diann_astral_multi_species.sh" \
        dsim \
        Dsim_diann \
        "${UPDOWN_RUN_ROOT}/Dsim/dsim-all-CDS-r2.02.longest_CDS_per_gene.protein.predicted.speclib" \
        "${UPDOWN_MS_DIR}/dsim-all-CDS-r2.02.longest_CDS_per_gene.protein.fasta" \
        DXN_20260401_sim_ \
        > "${UPDOWN_RUN_ROOT}/logs/run_diann_Dsim.log" 2>&1 &
    ms_search_pids+=("$!")

    RAW_DIR="$UPDOWN_RAW_DIR" OUT_ROOT_BASE="$UPDOWN_RUN_ROOT" \
        bash "${UPDOWN_CODE_DIR}/02_run_diann_astral_multi_species.sh" \
        dsec \
        Dsec_diann \
        "${UPDOWN_RUN_ROOT}/Dsec/dsec-all-CDS-r1.3.longest_CDS_per_gene.protein.predicted.speclib" \
        "${UPDOWN_MS_DIR}/dsec-all-CDS-r1.3.longest_CDS_per_gene.protein.fasta" \
        DXN_20260401_sec_ \
        > "${UPDOWN_RUN_ROOT}/logs/run_diann_Dsec.log" 2>&1 &
    ms_search_pids+=("$!")

    RAW_DIR="$UPDOWN_RAW_DIR" OUT_ROOT_BASE="$UPDOWN_RUN_ROOT" \
        bash "${UPDOWN_CODE_DIR}/02_run_diann_astral_multi_species.sh" \
        dere \
        Dere_diann \
        "${UPDOWN_RUN_ROOT}/Dere/dere-all-CDS-r1.05.longest_CDS_per_gene.protein.predicted.speclib" \
        "${UPDOWN_MS_DIR}/dere-all-CDS-r1.05.longest_CDS_per_gene.protein.fasta" \
        DXN_20260401_ere_ \
        > "${UPDOWN_RUN_ROOT}/logs/run_diann_Dere.log" 2>&1 &
    ms_search_pids+=("$!")

    ms_wait_for_jobs "${ms_search_pids[@]}"

fi

# -----------------------------------------------------------------------------
# 13C. Reproduce the published Figure S19 analysis
#
# The retained peptide x sample tables are reclassified against the updated
# counter-rule representative sites. For every gene, sample-level
# log2(downstream/upstream) values are averaged within sex and the female and
# male means are then averaged with equal weight. A gene observed in only one
# sex retains that sex-specific mean; missing gene-sample rows are not imputed.
# Dmel ISO1 and WBg are separate units, and Dmel AG/AS samples are excluded.
# Each final gene contributes once to the Edited-versus-Random two-sample KS
# test. The script exports sample-, sex- and final gene-level tables.
# -----------------------------------------------------------------------------

FS19_MS_RUN_ROOT="${FS19_MS_RUN_ROOT:-${UPDOWN_MS_DIR}/DIA_NN_RUN}"
FS19_CDSPOS_DIR="${FS19_CDSPOS_DIR:-/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed}"
FS19_REF_DIR="${FS19_REF_DIR:-${UPDOWN_REF_DIR}}"
FS19_OUT_DIR="${FS19_OUT_DIR:-${ANALYSIS_DIR}/Figure_S19_truncated_peptide_updown}"
FS19_MPLCONFIGDIR="${FS19_MPLCONFIGDIR:-/tmp/mpl_figure_s19}"
mkdir -p "${FS19_OUT_DIR}" "${FS19_MPLCONFIGDIR}"

MPLCONFIGDIR="${FS19_MPLCONFIGDIR}" \
    "${UPDOWN_PYTHON}" "${UPDOWN_CODE_DIR}/figure_s19_truncated_peptide_updown.py" \
    --ms-run-root "${FS19_MS_RUN_ROOT}" \
    --cdspos-dir "${FS19_CDSPOS_DIR}" \
    --ref-dir "${FS19_REF_DIR}" \
    --out-dir "${FS19_OUT_DIR}" \
    --random-multiplier 10 \
    > "${FS19_OUT_DIR}/figure_s19.log" 2>&1
