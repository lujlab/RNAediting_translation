#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
shopt -s nullglob

# Figures 2/3 and S6-S8/S10-S11: stalling, PausePred, gmRC and mixed models.
# Configure paths in section 00 and execute the required sections in order.
# Editing-site inputs must include transcript and CDS coordinate projections.

# ============================================================================
# 00. Input files, output directories, and analysis parameters
# ============================================================================
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
# Optional explicit handoff from the completed 05 project/edited-cds/handoff stages.
# Example: SITE_HANDOFF=/path/to/05_output/sites_for_06.env bash this_script.sh
if [[ -n "${SITE_HANDOFF:-}" ]]; then
    [[ -s "$SITE_HANDOFF" ]] || {
        echo "Missing SITE_HANDOFF: $SITE_HANDOFF" >&2
        exit 1
    }
    source "$SITE_HANDOFF"
fi
PROJECT_DIR="$(cd -- "${PIPELINE_DIR}/.." && pwd)"
CODE_DIR="${PIPELINE_DIR}/code/05"
PEAK_CODE="${CODE_DIR}/figure2_3_s6_edger"
PAUSE_CODE="${CODE_DIR}/figure_s6_s8_pausepred"
ORIGINAL_STALLING="${ORIGINAL_STALLING:-/path/to/Proj1_editing/codex/02_stalling}"
PEAK_INPUT="${PEAK_INPUT:-${ORIGINAL_STALLING}/3flank/peak_find}"
RESCUE_ROOT="${RESCUE_ROOT:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean}"
EDIT_ROOT="${EDIT_ROOT:-${RESCUE_ROOT}/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed}"
OUT_DIR="${OUT_DIR:-${PIPELINE_DIR}/results/05_stalling}"
REF_DIR="${REF_DIR:-${PROJECT_DIR}/ref_genome}"
SAMPLE_CONFIG_DIR="${SAMPLE_CONFIG_DIR:-${CODE_DIR}/config}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
PROJECT_R_LIBRARY="${PROJECT_R_LIBRARY:-${PIPELINE_DIR}/.runtime/05/R_library}"
if [[ -d "${PROJECT_R_LIBRARY}" ]]; then
    export R_LIBS_USER="${PROJECT_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"
fi
JOBS="${JOBS:-2}"
R_THREADS="${R_THREADS:-4}"
OVERWRITE="${OVERWRITE:-0}"
FLANK_DIR="${FLANK_DIR:-${OUT_DIR}/flank}"
MASTER_DIR="${MASTER_DIR:-${OUT_DIR}/master}"
EDGER_DIR="${EDGER_DIR:-${OUT_DIR}/edgeR_peak_only_20261002}"
RUN_LABEL="${RUN_LABEL:-peak_only_CDS150_m6A_20261002}"
TRANSCRIPT_MEDIAN_RULE="${TRANSCRIPT_MEDIAN_RULE:-any_both}"

RNALFOLD_CODE="${CODE_DIR}/figure3_s8_s11_rna_structure/run_RNALfold_T25Z3.sh"
RNALFOLD_INPUT_DIR="${RNALFOLD_INPUT_DIR:-${ORIGINAL_STALLING}/../RNALfold/37C_rnalfold/input}"
RNALFOLD_BIN="${RNALFOLD_BIN:-RNALfold}"
RNALFOLD_JOBS="${RNALFOLD_JOBS:-${JOBS}}"
HAIRPIN_DIR="${HAIRPIN_DIR:-${OUT_DIR}/RNALfold}"

PSITE_BAM_DIR="${PSITE_BAM_DIR:-${PROJECT_DIR}/05.psite/psite_2times}"
MERGED_BAM_DIR="${MERGED_BAM_DIR:-${OUT_DIR}/pausepred/merged_bam}"
MASTER_SHEET="${MASTER_SHEET:-${MERGED_BAM_DIR}/pausepred_sample_sheet.offset_plus3.tsv}"
# Validated mergefix batch; score and downstream analyses use 11 samples.
PAUSEPRED_EXCLUDE_SAMPLES="${PAUSEPRED_EXCLUDE_SAMPLES-AdarG_F,AdarG_M,AdarS_F,AdarS_M}"
PAUSEPRED_BATCH_ROOT="${PAUSEPRED_BATCH_ROOT:-${PIPELINE_DIR}/results/05_pausepred_mergefix_20260914}"
SHEET_DIR="${SHEET_DIR:-${PAUSEPRED_BATCH_ROOT}/sample_sheets}"
PAUSEPRED_DIR="${PAUSEPRED_DIR:-${PAUSEPRED_BATCH_ROOT}/pausepred_results}"
PAUSEPRED_QC_DIR="${PAUSEPRED_QC_DIR:-$(dirname "$PAUSEPRED_DIR")/qc}"
PAUSEPRED_SCRIPT="${PAUSEPRED_SCRIPT:-${PAUSE_CODE}/offline_pausepred_no_flag_filter_mergefix.pl}"
PAUSEPRED_PYTHON="${PAUSEPRED_PYTHON:-${PYTHON_BIN}}"
WINDOW_SIZE="${WINDOW_SIZE:-500}"
PAUSE_CUTOFF="${PAUSE_CUTOFF:-0}"
WINDOW_COVERAGE="${WINDOW_COVERAGE:-10}"
UPSTREAM_SEQ="${UPSTREAM_SEQ:-90}"
DOWNSTREAM_SEQ="${DOWNSTREAM_SEQ:-90}"
PAUSE_RUNTIME="${PAUSE_RUNTIME:-${PIPELINE_DIR}/.runtime/06}"
PYTHON_BIN="$(command -v "${PYTHON_BIN}")"
export PATH="$(dirname -- "${PYTHON_BIN}"):${PATH}"
export PATH="${PAUSE_RUNTIME}/bin:${PATH}"
export PERL5LIB="${PAUSE_CODE}/lib:${PAUSE_RUNTIME}/lib${PERL5LIB:+:${PERL5LIB}}"

TXINFO_DIR="${TXINFO_DIR:-${PROJECT_DIR}/05.psite}"
ANNOTATION_DIR="${ANNOTATION_DIR:-${FLANK_DIR}/annotation}"
MASTER_MANIFEST="${MASTER_MANIFEST:-${MASTER_DIR}/edgeR_master_manifest_counter_rule.tsv}"
ALL_SUMMARY="${EDGER_DIR}/edited_unedited_summary_codonA"
HAIRPIN_SUMMARY="${EDGER_DIR}/edited_unedited_hairpin_summary"
RUN_TAG="${RUN_TAG:-seq90_up90_down90_counterRule_newalledit}"

# Peak-only revision: CDS >=150 excluding stop, focal first15/last5 and m6A exclusion.
# Edited codons remain focal candidates and are excluded from flank backgrounds.
EDGER_CODE="${PEAK_CODE}/peak_only"
FOCAL_REFERENCE_MANIFEST="${FOCAL_REFERENCE_MANIFEST:-${EDGER_DIR}/inputs/focal_reference_manifest.tsv}"
RNA_SCALE_MANIFEST="${RNA_SCALE_MANIFEST:-${PIPELINE_DIR}/config/rna_count_scaling.tsv}"

# ============================================================================
# 00b. Predict stable local hairpins with RNALfold and merge CDS intervals
# ============================================================================
# Run section 00 then 00b independently to generate only RNALfold results.
# Use the original unedited longest-CDS FASTAs; output coordinates are 1-based
# and closed within those CDS sequences. Do not substitute edited/transcript FASTAs.
# Preserve --noLP --temp=25 --zscore=-3, one-hairpin-loop filtering and overlap merging.
# Species run sequentially; RNALFOLD_JOBS controls parallel chunks within each species.
# Sections 04 and 08 read HAIRPIN_DIR; section 10 reads its Dmel table by default.
command -v "$RNALFOLD_BIN" > /dev/null || {
    echo "RNALfold not found: set RNALFOLD_BIN to its executable path." >&2
    exit 1
}
for species in Dmel Dsim Dsec Dere; do
    [[ -s "$RNALFOLD_INPUT_DIR/$species.longest_CDS.fasta" ]] || {
        echo "Missing RNALfold FASTA: $RNALFOLD_INPUT_DIR/$species.longest_CDS.fasta" >&2
        exit 1
    }
done
mkdir -p "$HAIRPIN_DIR/logs"
for species in Dmel Dsim Dsec Dere; do
    env RNALFOLD_SPECIES="$species" \
        RNALFOLD_INPUT_FASTA="$RNALFOLD_INPUT_DIR/$species.longest_CDS.fasta" \
        RNALFOLD_BIN="$RNALFOLD_BIN" RNALFOLD_OUT_DIR="$HAIRPIN_DIR" \
        RNALFOLD_JOBS="$RNALFOLD_JOBS" RECORDS_PER_CHUNK="${RECORDS_PER_CHUNK:-100}" \
        bash "$RNALFOLD_CODE" > "$HAIRPIN_DIR/logs/$species.T25Z3.log" 2>&1
    [[ -s "$HAIRPIN_DIR/$species.fa.T25Z3.2.merged.tsv" ]]
    echo "RNALfold completed: $species"
done

# ============================================================================
# 01. Select representative transcripts and calculate focal/flanking counts
# ============================================================================
# Select the longest CDS, then the longest transcript, then the smallest ID.
# Flanking coverage uses 30 codons on each side; terminal exclusions are
# 15 codons at the CDS start and 5 codons at the CDS end.
mkdir -p "${ANNOTATION_DIR}" "${FLANK_DIR}/samples"
overwrite_args=()
[[ "$OVERWRITE" == 0 ]] || overwrite_args=(--overwrite)
for sp in dmel dere dsec dsim; do
    case "$sp" in
        dmel)
            tag=dmelr6.12
            txinfo=dmel_r6.12
            m6a=m6A_list.handle
            edit=ISO1
            sample_table=samples.dmel.all.tsv
            ;;
        dere)
            tag=derer1.05
            txinfo=dere_all_r1.05
            m6a=m6A_list.droEre2ToDm6.handle
            edit=ere
            sample_table=samples.dere.tsv
            ;;
        dsec)
            tag=dsecr1.3
            txinfo=dsec_all_r1.3
            m6a=m6A_list.droSec1ToDm6.handle
            edit=sec
            sample_table=samples.dsec.tsv
            ;;
        dsim)
            tag=dsimr2.02
            txinfo=dsim_all_r2.02
            m6a=m6A_list.droSim1.handle
            edit=sim
            sample_table=samples.dsim.tsv
            ;;
    esac
    "$PYTHON_BIN" "${PEAK_CODE}/select_longest_cds_transcripts_counter_rule.py" \
        --txinfo "${TXINFO_DIR}/${txinfo}_txinfo.biotype.tsv" \
        --ids-out "${ANNOTATION_DIR}/${sp}.representative_transcripts.txt" \
        --table-out "${ANNOTATION_DIR}/${sp}.representative_transcripts.tsv"
    "$PYTHON_BIN" "${PEAK_CODE}/codon_flank_pipeline.py" prepare \
        --annotation "${PEAK_INPUT}/${tag}.gtf.cds.bed.ordered.1codon" \
        --longest-transcripts "${ANNOTATION_DIR}/${sp}.representative_transcripts.txt" \
        --m6a "${PEAK_INPUT}/${m6a}" \
        --output-annotation "${ANNOTATION_DIR}/${tag}.gtf.cds.bed.ordered.1codon.longest.tsv" \
        --output-m6a-map "${ANNOTATION_DIR}/${sp}.m6a.tsv" "${overwrite_args[@]}"
    sample_input="${SAMPLE_CONFIG_DIR}/${sample_table}"
    if [[ "$sp" == dmel ]]; then
        awk -F '\t' 'NR == 1 || $1 ~ /^(WBg_F[123]|WBg_M[123]|AdarG_[FM][123]|AdarS_[FM][123]|Dmel_ISO1_(female|male)|Dmel_Adar5G1_rep[12])$/' \
            "$sample_input" > "${FLANK_DIR}/samples/samples.dmel.tsv"
        sample_input="${FLANK_DIR}/samples/samples.dmel.tsv"
    fi
    "$PYTHON_BIN" "${PEAK_CODE}/codon_flank_pipeline.py" count \
        --annotation-longest "${ANNOTATION_DIR}/${tag}.gtf.cds.bed.ordered.1codon.longest.tsv" \
        --m6a-map "${ANNOTATION_DIR}/${sp}.m6a.tsv" \
        --editing-sites "${EDIT_ROOT}/${edit}.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv" \
        --sample-sheet "$sample_input" --output-dir "${FLANK_DIR}/results_${sp}" \
        --flanks 30 --exclude-start 15 --exclude-end 5 --jobs "$JOBS" "${overwrite_args[@]}"
done

# ============================================================================
# 02. Build group-specific master tables and the combined ISO1 female/male table
# ============================================================================
# The combined ISO1 model uses the female and male libraries as separate columns.
mkdir -p "${MASTER_DIR}"
awk -F '\t' -v OFS='\t' -v root="${FLANK_DIR}" '
  NR == 1 {print; next}
  $1 != "dmel_ISO1_F" && $1 != "dmel_ISO1_M" {
    $8 = root "/results_" $2; print
  }' "${SAMPLE_CONFIG_DIR}/stalling_groups.tsv" > "${MASTER_DIR}/stalling_groups.tsv"
master_args=(--groups "${MASTER_DIR}/stalling_groups.tsv" --output-dir "$MASTER_DIR" --flanks 30 --jobs "$JOBS")
[[ "$OVERWRITE" == 0 ]] || master_args+=(--overwrite)
"$PYTHON_BIN" "${PEAK_CODE}/build_stalling_master_tables.py" "${master_args[@]}"
"$PYTHON_BIN" "${PEAK_CODE}/combine_iso1_FM_master_counter_rule.py" \
    --female "${FLANK_DIR}/results_dmel/Dmel_ISO1_female.flank30.tsv" \
    --male "${FLANK_DIR}/results_dmel/Dmel_ISO1_male.flank30.tsv" \
    --output "${MASTER_DIR}/dmel_ISO1_FM.flank30.master.tsv" \
    --manifest "${MASTER_DIR}/dmel_ISO1_FM.master_manifest.tsv"
{
    cat "${MASTER_DIR}/master_table_manifest.tsv"
    awk 'NR > 1' "${MASTER_DIR}/dmel_ISO1_FM.master_manifest.tsv"
} > "$MASTER_MANIFEST"

# ============================================================================
# 03. Correct paired-end RNA counts and fit edgeR quasi-likelihood models
# ============================================================================
# Test the ribosome-specific focal effect with library-specific intercepts.
# The design is ~ library + region + ribo_focal.
# Call peaks at positive logFC and one-sided BH FDR <= 0.05; retain zero focal RPF.
# Transcript coverage is filtered using the any_both rule.
# Apply sample-specific RNA factors to focal/flank counts and library sizes
# before filtering and normalization. Supply unscaled master tables.
# Replicates remain separate; fractional corrected counts are retained.
# Reuse unscaled masters by setting MASTER_MANIFEST in section 00 and running
# this block independently. OVERWRITE=0 preserves completed model outputs.
"$PYTHON_BIN" "$EDGER_CODE/prepare_references.py" --annotation-dir "$ANNOTATION_DIR" \
    --master-manifest "$MASTER_MANIFEST" --output-dir "$EDGER_DIR/inputs"
env RSCRIPT_BIN="$RSCRIPT_BIN" MASTER_MANIFEST="$MASTER_MANIFEST" \
    FOCAL_REFERENCE_MANIFEST="$FOCAL_REFERENCE_MANIFEST" \
    RNA_SCALE_MANIFEST="$RNA_SCALE_MANIFEST" \
    EDGER_V2_DIR="$EDGER_DIR" EDGER_JOBS="${EDGER_JOBS:-4}" R_THREADS="$R_THREADS" \
    PROJECT_R_LIBRARY="$PROJECT_R_LIBRARY" OVERWRITE_EDGER="$OVERWRITE" \
    GROUP_IDS='dmel_Adar5G1_M,dmel_ISO1_FM,dmel_WBg_F,dmel_WBg_M,dmel_AdarG_F,dmel_AdarG_M,dmel_AdarS_F,dmel_AdarS_M,dere_Dere_F,dere_Dere_M,dsec_Dsec_F,dsec_Dsec_M,dsim_w501_F,dsim_w501_M' \
    SAMPLE_FRACTION=1 SAMPLE_SEED=20260805 NORMALIZATION_METHOD=TMM \
    MIN_STALLING_LOG2FC=0 \
    MIN_TRANSCRIPT_MEDIAN_RNA=0 MIN_TRANSCRIPT_MEDIAN_RPF=0 EXCLUDE_HEAD_REPLICATES=0 \
    TRANSCRIPT_MEDIAN_RULE="$TRANSCRIPT_MEDIAN_RULE" RUN_LABEL="$RUN_LABEL" \
    bash "${EDGER_CODE}/run_peak_only.sh"

# ============================================================================
# 04. Summarize edited/unedited event rates and hairpin-stratified enrichment
# ============================================================================
# Restrict unedited controls to reference codons containing adenine.
inputs=("${EDGER_DIR}"/*.flank30."${RUN_LABEL}".edgeR.tsv.gz)
((${#inputs[@]})) || {
    echo "No edgeR outputs for RUN_LABEL=${RUN_LABEL}" >&2
    exit 1
}
"$RSCRIPT_BIN" "${EDGER_CODE}/summarize_edited_unedited_peaks.R" \
    --result-dir "$EDGER_DIR" --manifest "$MASTER_MANIFEST" --run-label "$RUN_LABEL" \
    --codon-dir "$REF_DIR" --editing-dir "$EDIT_ROOT" --output-dir "$ALL_SUMMARY"
"$RSCRIPT_BIN" "${EDGER_CODE}/summarize_edited_unedited_hairpin_peaks.R" \
    --result-dir "$EDGER_DIR" --manifest "$MASTER_MANIFEST" --run-label "$RUN_LABEL" \
    --annotation-dir "$ANNOTATION_DIR" --hairpin-dir "$HAIRPIN_DIR" \
    --codon-dir "$REF_DIR" --editing-dir "$EDIT_ROOT" --output-dir "$HAIRPIN_SUMMARY"
# ============================================================================
# 04b. Peak-only enrichment figures and editable manuscript panels/tables
# ============================================================================
EDGER_PLOT_PYTHON="${EDGER_PLOT_PYTHON:-${PIPELINE_DIR}/.runtime/06/venv/bin/python}"
if [[ ! -x "$EDGER_PLOT_PYTHON" ]]; then EDGER_PLOT_PYTHON="$PAUSEPRED_PYTHON"; fi
EDGER_PLOT_LIB="${PIPELINE_DIR}/.runtime/06/python_library${PYTHONPATH:+:$PYTHONPATH}"
env PYTHONPATH="$EDGER_PLOT_LIB" "$EDGER_PLOT_PYTHON" "${EDGER_CODE}/plot_peak_enrichment.py" \
    --edgeR-dir "$EDGER_DIR"
env PYTHONPATH="$EDGER_PLOT_LIB" "$EDGER_PLOT_PYTHON" "${EDGER_CODE}/plot_peak_enrichment.py" \
    --edgeR-dir "$EDGER_DIR" --include-adarg-adars
env PYTHONPATH="$EDGER_PLOT_LIB" "$EDGER_PLOT_PYTHON" "${EDGER_CODE}/publication_figures.py" \
    --edgeR-dir "$EDGER_DIR" --out-dir "$EDGER_DIR/publication"
env PYTHONPATH="$EDGER_PLOT_LIB" "$EDGER_PLOT_PYTHON" "${EDGER_CODE}/publication_tables.py" \
    --edgeR-dir "$EDGER_DIR" --out-dir "$EDGER_DIR/publication"

# ============================================================================
# 05. Merge BAMs and derive read-length offsets from modal P-site offsets +3 nt
# ============================================================================
"$PYTHON_BIN" "${PAUSE_CODE}/merge_bams_and_build_pausepred.py" \
    --source-dir "$PSITE_BAM_DIR" --output-dir "$MERGED_BAM_DIR" \
    --jobs "$JOBS" --threads "$R_THREADS" --stage all

# ============================================================================
# 06. Validate BAM/reference inputs and prepare species-specific sample sheets
# ============================================================================
env MASTER_SHEET="$MASTER_SHEET" SHEET_DIR="$SHEET_DIR" \
    SCRIPT="$PAUSEPRED_SCRIPT" REF_DIR="$REF_DIR" \
    OUTDIR="$PAUSEPRED_DIR" WINDOW_SIZE="$WINDOW_SIZE" PAUSE_CUTOFF="$PAUSE_CUTOFF" \
    WINDOW_COVERAGE="$WINDOW_COVERAGE" UPSTREAM_SEQ="$UPSTREAM_SEQ" DOWNSTREAM_SEQ="$DOWNSTREAM_SEQ" \
    SPECIES_JOBS="${SPECIES_JOBS:-2}" MAX_JOBS="${MAX_JOBS:-1}" SKIP_EXISTING="${SKIP_EXISTING:-1}" \
    bash "${PAUSE_CODE}/run_pausepred_4species_safe_parallel.sh" --review

# ============================================================================
# 07. Run PausePred with 500-nt windows and 90-nt flanking sequences
# ============================================================================
env MASTER_SHEET="$MASTER_SHEET" SHEET_DIR="$SHEET_DIR" \
    SCRIPT="$PAUSEPRED_SCRIPT" REF_DIR="$REF_DIR" \
    OUTDIR="$PAUSEPRED_DIR" WINDOW_SIZE="$WINDOW_SIZE" PAUSE_CUTOFF="$PAUSE_CUTOFF" \
    WINDOW_COVERAGE="$WINDOW_COVERAGE" UPSTREAM_SEQ="$UPSTREAM_SEQ" DOWNSTREAM_SEQ="$DOWNSTREAM_SEQ" \
    SPECIES_JOBS="${SPECIES_JOBS:-2}" MAX_JOBS="${MAX_JOBS:-1}" SKIP_EXISTING="${SKIP_EXISTING:-1}" \
    bash "${PAUSE_CODE}/run_pausepred_4species_safe_parallel.sh" --run

# ============================================================================
# 07b. Pause-score tail distribution: reads >=2, coverage >=20%, score >=0
# ============================================================================
# Execute 00 then 07b–08b to reuse the latest validated CSVs without BAM recounting.
# Requires completed RNALfold structures for section 08. To reuse the structures
# from the reference figures, set HAIRPIN_DIR to the original
# ${ORIGINAL_STALLING}/3flank/hairpin directory before executing section 00.
# To reproduce historical inputs, override PAUSEPRED_DIR and SHEET_DIR in 00.
# The original plotting script retains score zero and uses all eligible positions
# as its denominator; it does not filter at the later enrichment cutoff of 5.
"$PAUSEPRED_PYTHON" -c 'import numpy, pandas, matplotlib, scipy'
"$PAUSEPRED_PYTHON" "${PAUSE_CODE}/02_03_01_plot_pause_score_distribution.py" \
    --indir "$PAUSEPRED_DIR" --pattern '*.pausepred.csv' --exclude-samples "$PAUSEPRED_EXCLUDE_SAMPLES" \
    --outprefix "$PAUSEPRED_DIR/pause_score_distribution.reads2_cov20_pause0" \
    --xmin 0 --xmax 30 --ymin 0 --ymax 0.7 --step 0.5 \
    --min-reads 2 --min-coverage 20

# ============================================================================
# 08. Calculate editing enrichment and generate hairpin-stratified figures
# ============================================================================
env PATH="$(dirname "$(command -v "$PAUSEPRED_PYTHON")"):$PATH" \
    PAUSEPRED_EXCLUDE_SAMPLES="$PAUSEPRED_EXCLUDE_SAMPLES" \
    PAUSEPRED_DIR="$PAUSEPRED_DIR" SAMPLE_SHEETS_DIR="$SHEET_DIR" PAUSEPRED_QC_DIR="$PAUSEPRED_QC_DIR" \
    RESCUE_ROOT="$RESCUE_ROOT" EDIT_ROOT="$EDIT_ROOT" IPOS_DIR="${CODE_DIR}/figure_s6_s8_pausepred/enrichment" \
    HAIRPIN_DIR="$HAIRPIN_DIR" \
    RUN_TAG="$RUN_TAG" PAUSE_CUTOFFS=5,8,10 ENRICH_WINDOWS=0,2 \
    ENRICH_JOBS="$JOBS" HAIRPIN_JOBS="$JOBS" \
    bash "${PAUSE_CODE}/run_pausepred_edit_unedit_hairpin_workflow.sh" --no-wait --windows 2 --offset 5

# 08b. Export Table S8 for the 2-nt window in the specified library order.
infile="$PAUSEPRED_DIR/enrichment_${RUN_TAG}_unedit_A_base/plots/all_samples.pausepred_enrichment.merged.tsv"
outfile="$PAUSEPRED_DIR/enrichment_${RUN_TAG}_unedit_A_base/plots/tableS8_pausepred_all_samples.pausepred_enrichment.merged_2nt_ordered.tsv"
[[ -s "$infile" ]] || {
    echo "Missing PausePred merged table: $infile" >&2
    exit 1
}

awk -F '\t' -v OFS='\t' -v outfile="$outfile" '
BEGIN {
  split("Dmel_Adar5G1_M Dmel_ISO1_head_female_ribo Dmel_ISO1_head_male_ribo WBg_F WBg_M Dsim_F Dsim_M Dsec_F Dsec_M Dere_F Dere_M", want, " ")
}
NR == 1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "sample") s = i
    if ($i == "window_nt") w = i
    if ($i == "nonediting_not_in_peak") non_out = i
    if ($i == "nonediting_in_peak") non_in = i
    if ($i == "nonediting_peak_fraction") non_frac = i
    if ($i == "editing_not_in_peak") edit_out = i
    if ($i == "editing_in_peak") edit_in = i
    if ($i == "editing_peak_fraction") edit_frac = i
    if ($i == "p_value") pv = i
  }
  print "Library", "Unedited_outside_peak", "Unedited_inside_peak", "Unedited_peak_fraction", "Edited_outside_peak", "Edited_inside_peak", "Edited_peak_fraction", "P_value" > outfile
  next
}
$w == 2 {
  rec[$s] = $s OFS $(non_out) OFS $(non_in) OFS $(non_frac) OFS $(edit_out) OFS $(edit_in) OFS $(edit_frac) OFS $(pv)
}
END {
  for (i = 1; i <= length(want); i++) if (want[i] in rec) print rec[want[i]] >> outfile
}
' "$infile"

# ============================================================================
# 08c. Redraw the two A-base bar charts in the supplied Figure 10 style
# ============================================================================
"$PAUSEPRED_PYTHON" "${PAUSE_CODE}/plot_pausepred_reference_style.py" \
    --pausepred-dir "$PAUSEPRED_DIR" --run-tag "$RUN_TAG" --cutoff 5 --window 2

# ============================================================================
# 08d. Edited peak overlap: PausePred and RNA /2 edgeR, 11 selected libraries
# ============================================================================
# Main table counts edited nucleotide sites. Also exports CDS-only and codon
# tables. The ISO1 F/M rows share the existing jointly fitted edgeR ISO1 model.
"$PAUSEPRED_PYTHON" "${EDGER_CODE}/compare_pausepred_edgeR_edited_peaks.py" \
    --pausepred-dir "$PAUSEPRED_DIR" --run-tag "$RUN_TAG" \
    --edgeR-dir "$EDGER_DIR" --edgeR-run-label "$RUN_LABEL" \
    --editing-dir "$EDIT_ROOT" \
    --out-dir "${EDGER_DIR}/pausepred_edited_peak_overlap" \
    --cutoff 5 --window 2 --jobs "$JOBS"

# =============================================================================
# 09. Main gmRC: input paths and plotting environment
# =============================================================================
GMRC_ROOT="$PIPELINE_DIR/code/05/figure2_3_s7_s10_s11_gmrc"
WORKSPACE_ROOT="$(cd "$GMRC_ROOT/../../.." && pwd)"
PROJECT_ROOT="$(dirname -- "$WORKSPACE_ROOT")"
FIGURE_ROOT="${FIGURE_ROOT:-$GMRC_ROOT/figure_update}"
PUBLICATION_MEAN_DIR="${PUBLICATION_MEAN_DIR:-$GMRC_ROOT/group_sex_mean}"
RESCUE_ROOT="${RESCUE_ROOT:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean}"
EDITED_REFERENCE_DIR="${EDITED_REFERENCE_DIR:-$RESCUE_ROOT/09_Ipos/rebuilt_edited_references_3R/split}"
M6A_REFERENCE_CONFIG="${M6A_REFERENCE_CONFIG:-$GMRC_ROOT/../config/m6A_references.tsv}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
RENDER_PYTHON="${RENDER_PYTHON:-/usr/bin/python3}"
RENDER_PYTHONPATH="${RENDER_PYTHONPATH:-$GMRC_ROOT/runtime/python36}"
JOBS="${JOBS:-4}"
OVERWRITE="${OVERWRITE:-0}"

RSCU_FILE="${RSCU_FILE:-/path/to/RNAi/240719drosophila_check/RSCU_flybase.txt}"
CODON_TABLE_FILE="${CODON_TABLE_FILE:-/path/to/RNAi/240719drosophila_check/Codon_table.new.txt}"
GLM_BIGWIG="${GLM_BIGWIG:-/path/to/RNA_editing_DATAprocessing/ref_genome/dm6_phyloP27way/dm6.phyloP27way.bw}"
GLM_TRANS_FILE="${GLM_TRANS_FILE:-/path/to/Proj2_allele_specific_translation/02_codon/dmel-r6.12_CDS_position}"
GLM_CODON_FILE="${GLM_CODON_FILE:-$PROJECT_ROOT/ref_genome/dmel-all-CDS-r6.12.fasta.edited.tab.split}"
GLM_TX_GENE_FILE="${GLM_TX_GENE_FILE:-$PROJECT_ROOT/05.psite/dmel_r6.12_txinfo.biotype.tsv}"
GLM_TE_FILE="${GLM_TE_FILE:-$PROJECT_ROOT/10.stalling/04_GLM/log2TE_matrix.tsv}"
GLM_STRUCTURE_FILE="${GLM_STRUCTURE_FILE:-$PROJECT_ROOT/10.stalling/04_GLM/RNALfold_result/Dmel.fa.T25Z3.2.merged.tsv}"
GLM_EDIT_WBG="${GLM_EDIT_WBG:-/path/to/Proj1_editing/code/01_identify_sites/sites_final/WBG.candidate.AtoG.annotations.long.tsv.new}"

export R_LIBS_USER="$FIGURE_ROOT/R_library${R_LIBS_USER:+:$R_LIBS_USER}"
export M6A_REFERENCE_CONFIG FIGURE_ROOT

require_file() {
    [[ -s "$1" ]] || {
        echo "Missing required file: $1" >&2
        exit 1
    }
}

render_python() {
    env PYTHONPATH="${RENDER_PYTHONPATH}${PYTHONPATH:+:${PYTHONPATH}}" "$RENDER_PYTHON" "$@"
}

# =============================================================================
# 10. Calculate replicate gmRC and group means
# =============================================================================
echo "[compute] corrected single-library mRC/gmRC and group means"
"$PYTHON_BIN" "$GMRC_ROOT/code/recompute.py" prepare
"$PYTHON_BIN" "$GMRC_ROOT/code/recompute.py" run --jobs "$JOBS"
"$PYTHON_BIN" "$GMRC_ROOT/code/prepare_figure_manifest.py"
"$PYTHON_BIN" "$GMRC_ROOT/code/prepare_glm_inputs.py"

# =============================================================================
# 11. Metagene profiles and Figure 3C/D
# =============================================================================
echo "[profiles] metagene, F3C and F3D; FS10 uses the effective-replicate median threshold below"
require_file "$PUBLICATION_MEAN_DIR/mean_gmRC_table_manifest.tsv"
require_file "$M6A_REFERENCE_CONFIG"

env METAGENE_SOURCE_SCRIPT="$FIGURE_ROOT/code/All_metagene_mean_existingMean_any1libLm02_mean_sd.R" \
    METAGENE_MEAN_DIR="$PUBLICATION_MEAN_DIR" \
    METAGENE_OUTDIR="$FIGURE_ROOT/metagene" \
    METAGENE_ALLEDIT_DIR="$RESCUE_ROOT/04_alledit_new/results" \
    METAGENE_CDSPOS_DIR="$RESCUE_ROOT/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed" \
    METAGENE_REFERENCE_CONFIG="$GMRC_ROOT/reference/metagene_references.tsv" \
    METAGENE_CONTROL_CENTER_A_ONLY=true \
    M6A_REFERENCE_CONFIG="$M6A_REFERENCE_CONFIG" \
    "$RSCRIPT_BIN" --vanilla "$FIGURE_ROOT/code/run_metagene_counter_rule_streamed.R"

M6A_REFERENCE_CONFIG="$M6A_REFERENCE_CONFIG" "$RSCRIPT_BIN" --vanilla \
    "$FIGURE_ROOT/code/F3C_precompute.R" "$PUBLICATION_MEAN_DIR" "$RESCUE_ROOT" "$FIGURE_ROOT/F3C"
F3C_TABLE_DIR="$FIGURE_ROOT/F3C" "$RSCRIPT_BIN" --vanilla \
    "$FIGURE_ROOT/code/F3C_comparison.R"

AA_UNEDIT_A_ONLY=false M6A_REFERENCE_CONFIG="$M6A_REFERENCE_CONFIG" \
    "$RSCRIPT_BIN" --vanilla "$FIGURE_ROOT/code/F3D_classify.R" \
    "$PUBLICATION_MEAN_DIR" "$FIGURE_ROOT/AA" "$EDITED_REFERENCE_DIR" \
    "$CODON_TABLE_FILE" "$RSCU_FILE"
"$RSCRIPT_BIN" --vanilla "$FIGURE_ROOT/code/F3D_pool_summary.R" "$FIGURE_ROOT/AA"
"$RSCRIPT_BIN" --vanilla "$FIGURE_ROOT/code/F3D_models.R" "$FIGURE_ROOT/AA"

# =============================================================================
# 12. Figure S10: main gmRC, 3N threshold and sex model
# =============================================================================
RSCRIPT="$RSCRIPT_BIN" bash "$WORKSPACE_ROOT/code/05/figure_s10_codon_preference/run.sh"

# =============================================================================
# 13. Figure 3E and Figure S11 statistical models
# =============================================================================
echo "[models] F3E and the final 1:20 FS11 with the editing-level term"
require_file "$GLM_BIGWIG"
command -v bigWigAverageOverBed > /dev/null || {
    echo "bigWigAverageOverBed not found" >&2
    exit 1
}
mkdir -p "$FIGURE_ROOT/phyloP_annotation"
for sample in WBg_F{1..3} WBg_M{1..3}; do
    input="$GMRC_ROOT/glm_inputs/$sample.codonTE_AEP.AA"
    output="$FIGURE_ROOT/phyloP_annotation/$sample.codonTE_AEP.AA.phyloP.txt"
    require_file "$input"
    if [[ "$OVERWRITE" == 1 || ! -s "$output" ]]; then
        bash "$FIGURE_ROOT/code/annotate_codon_phyloP.sh" "$input" "$GLM_BIGWIG" "$output"
    fi
done

env GMRC_OUTPUT_DIR="$FIGURE_ROOT/F3E" \
    GLM_PHYLOP_DIR="$FIGURE_ROOT/phyloP_annotation" \
    GLM_TRANS_FILE="$GLM_TRANS_FILE" GLM_CODON_FILE="$GLM_CODON_FILE" \
    GLM_TX_GENE_FILE="$GLM_TX_GENE_FILE" GLM_TE_FILE="$GLM_TE_FILE" \
    GLM_RSCU_FILE="$RSCU_FILE" GLM_STRUCTURE_FILE="$GLM_STRUCTURE_FILE" \
    GLM_EDIT_WBG="$GLM_EDIT_WBG" M6A_REFERENCE_CONFIG="$M6A_REFERENCE_CONFIG" \
    GLM_FAST_AGGREGATION_HELPER="$FIGURE_ROOT/code/collapse_replicates_fast.R" \
    "$RSCRIPT_BIN" --vanilla "$FIGURE_ROOT/code/F3E_group_mean_model_WT_m6A.R"
GMRC_OUTPUT_DIR="$FIGURE_ROOT/F3E" "$RSCRIPT_BIN" --vanilla \
    "$FIGURE_ROOT/code/plot_WT_fixed_effects.R"
FS11_SOURCE_MODEL_DIR="$FIGURE_ROOT/F3E" \
    FS11_EDITLEVEL_OUTPUT_DIR="$FIGURE_ROOT/FS11_with_editing_level" \
    "$RSCRIPT_BIN" --vanilla "$FIGURE_ROOT/code/FS11_matched_with_editing_level.R"

# =============================================================================
# 14. Export individual native-vector panels
# =============================================================================
echo "[render] native-vector publication panels"
export GLM_RESULT_DIR="$FIGURE_ROOT/F3E"
for part in metagene f3c f3d glm; do
    render_python "$FIGURE_ROOT/code/render_figures.py" "$part"
done
export PUBLICATION_GROUP_MEAN_DIR="$PUBLICATION_MEAN_DIR"
render_python "$FIGURE_ROOT/code/render_F2D_FS7BC_A.py"
render_python "$FIGURE_ROOT/code/render_panels.py"
render_python "$FIGURE_ROOT/code/render_FS11_with_editing_level.py"
cp "$FIGURE_ROOT/FS11_with_editing_level/FS11_with_editing_level_fixed_effects.tsv" \
    "$FIGURE_ROOT/tables/FS11_fixed_effects.tsv"

# Full-page composition is not required; individual panels are the deliverables.
