# Input paths and output locations for the Figure 5 analysis.
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)}"
PROJECT_DIR="${PROJECT_DIR:-$(dirname -- "$PIPELINE_DIR")}"
CODE_DIR="${CODE_DIR:-$PIPELINE_DIR/code/08}"
RESCUE_ROOT="${RESCUE_ROOT:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean}"
OUT_DIR="${OUT_DIR:-$PIPELINE_DIR/results/08_editing_expression}"
FEATURECOUNTS_DIR="${FEATURECOUNTS_DIR:-$OUT_DIR/featurecounts}"
RNA_ALIGNMENT_ROOT="${RNA_ALIGNMENT_ROOT:-$PIPELINE_DIR/results}"
RNA_SAMPLE_TABLE="${RNA_SAMPLE_TABLE:-$CODE_DIR/setup/rna_samples.tsv}"
RNA_GTF="${RNA_GTF:-$PROJECT_DIR/ref_genome/dmel-all-r6.12.gtf}"
FEATURECOUNTS_BIN="${FEATURECOUNTS_BIN:-featureCounts}"
FEATURECOUNTS_THREADS="${FEATURECOUNTS_THREADS:-16}"
EDIT_INPUT_DIR="${EDIT_INPUT_DIR:-$RESCUE_ROOT/13_AdarSG_Lm_cluster}"
SITE_ROOT="${SITE_ROOT:-$RESCUE_ROOT/01_rebuild/result/sites_candidate_rescue_3R_18907675_all_test2}"
ALL_SITE_FILE="${ALL_SITE_FILE:-$SITE_ROOT/raw_rows_by_group/ISO1.candidate.AtoG.annotations.long.tsv}"
SITE_ANNOTATION_FILE="${SITE_ANNOTATION_FILE:-$SITE_ROOT/candidate_lists_by_group/ISO1.candidate.detected.annotated.tsv}"
CDSPOS_FILE="${CDSPOS_FILE:-$RESCUE_ROOT/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed/ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv}"
# The publication clustering panels use degPatterns with six clusters.
CLUSTER_METHOD="${CLUSTER_METHOD:-degPatterns}"
[[ "$CLUSTER_METHOD" == degPatterns ]] || {
    echo "Publication workflow requires CLUSTER_METHOD=degPatterns" >&2
    return 2
}
RNA_OUT_DIR="${RNA_OUT_DIR:-$OUT_DIR/RNA_LRT_$CLUSTER_METHOD}"
EDIT_OUT_DIR="${EDIT_OUT_DIR:-$OUT_DIR/editing_within_sex}"
INTERACTION_OUT_DIR="${INTERACTION_OUT_DIR:-$OUT_DIR/editing_interaction}"
# Preserve the original expression-background definition for the GO panel.
EDIT_GO_EXPRESSION_FILE="${EDIT_GO_EXPRESSION_FILE:-$PROJECT_DIR/06_Xtail_TE/xtail_input/RNA_count_matrix_WBg_Adar_merged.txt}"
PROTEIN_INPUT="${PROTEIN_INPUT:-$RESCUE_ROOT/04_alledit_new/03_Lp_Lm_Lr_MS_scatter/00_inputs/dmel_Lp_union.counter_rule_LmLr.tsv}"
PROTEIN_OUT_DIR="${PROTEIN_OUT_DIR:-$OUT_DIR/protein_output_bias}"
FIGURE_DIR="${FIGURE_DIR:-$OUT_DIR/figures}"
ANALYSIS_SEED="${ANALYSIS_SEED:-20260819}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
R07_LIBRARY="${R07_LIBRARY:-$PIPELINE_DIR/.runtime/07/R_library}"
export R_LIBS_USER="$R07_LIBRARY${R_LIBS_EXTRA:+:$R_LIBS_EXTRA}${R_LIBS_USER:+:$R_LIBS_USER}"
export PIPELINE_DIR PROJECT_DIR CODE_DIR RESCUE_ROOT OUT_DIR FEATURECOUNTS_DIR EDIT_INPUT_DIR
export SITE_ROOT ALL_SITE_FILE SITE_ANNOTATION_FILE CDSPOS_FILE CLUSTER_METHOD RNA_OUT_DIR
export EDIT_OUT_DIR INTERACTION_OUT_DIR EDIT_GO_EXPRESSION_FILE PROTEIN_INPUT PROTEIN_OUT_DIR
export FIGURE_DIR ANALYSIS_SEED R07_LIBRARY

ADAR_ALL_EDIT_DIR="${ADAR_ALL_EDIT_DIR:-$PROJECT_DIR/05.psite/psite_2times}"
ADAR_OUT_DIR="${ADAR_OUT_DIR:-$OUT_DIR/adar_autosite}"
export ADAR_ALL_EDIT_DIR ADAR_OUT_DIR

export PYTHONPATH="$PIPELINE_DIR/.runtime/07/python_library${PYTHONPATH:+:$PYTHONPATH}"
