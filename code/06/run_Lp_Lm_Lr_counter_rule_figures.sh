#!/usr/bin/env bash
set -euo pipefail

# Rebuild Lm/Lr-linked DIA-NN / MaxQuant figures after the counter-rule
# alledit recalculation.  All outputs are written beneath 04_alledit_new.
#
# Step 1: derive per-replicate and sex-mean Lm/Lr from new A-site alledit.
#   Lm = alt RNA / (ref RNA + alt RNA)
#   Lr = alt Ribo / (ref Ribo + alt Ribo)
#   Each metric independently requires coverage >= 3.
#
# Step 2: map each DIA-NN Lp site to new CDS coordinates and append Lm/Lr.
# Step 3: create (a) ISO1/sim DIA-NN/MaxQuant/combined Lp-vs-Lm,Lr plots,
#         (b) I-position plots, (c) Dmel genotype violin plots, and
#         (d) 7-panel male/female Lp-vs-Lm,Lr scatter plots.

PROJECT=${SOURCE_PROJECT:?Set SOURCE_PROJECT}
MS_ROOT=${MS_ROOT:?Set MS_ROOT}
OUTPUT_ROOT=${OUTPUT_ROOT:-$PROJECT/04_alledit_new/03_Lp_Lm_Lr_MS_scatter}
CODE_DIR=${CODE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}
ALLEDIT_DIR=${ALLEDIT_DIR:-$PROJECT/04_alledit_new/results}
CDSPOS_DIR=${CDSPOS_DIR:-$PROJECT/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed}
INPUT_DIR=$OUTPUT_ROOT/00_inputs
# These defaults must be visible to the R subprocesses as well as Python.
export MS_ROOT CDSPOS_DIR

mkdir -p "$INPUT_DIR" "$OUTPUT_ROOT/01_ISO1_sim_union" \
    "$OUTPUT_ROOT/02_Ipos_Lp_Lm_3sex" "$OUTPUT_ROOT/03_Dmel_AG_AS_WBg_violin" \
    "$OUTPUT_ROOT/04_Lp_Lm_Lr_scatter"

for code in recalculate_Lm_Lr_min3_or_coverage_v2.py match_Lp_with_min3_LmLr_v5.py \
    plot_ISO1_sim_combined_MaxQuant_DIANN_Lp_vs_LmLr.R \
    plot_Ipos_Lp_Lm_3sex_2methods_A4_v4.R \
    plot_Lp_ratio_difference_violin_paired_by_sex_SX31style_current_wide_AG_WBg_AS_A4quarter.R \
    plot_male_Lp_vs_Lm_Lr_longestTx.R; do
    [[ -s "$CODE_DIR/$code" ]] || {
        echo "Missing code: $CODE_DIR/$code" >&2
        exit 1
    }
done

# 1. New Lm/Lr table from all 42 eligible A-site alledit samples.
"${PYTHON:-python3}" "$CODE_DIR/recalculate_Lm_Lr_min3_or_coverage_v2.py" \
    --input-dir "$ALLEDIT_DIR" --out-prefix "$INPUT_DIR/all_species_LmLr_min3" \
    --min-coverage 3
LMLR="$INPUT_DIR/all_species_LmLr_min3.sex_site_mean_LmLr.tsv"

# 2. Attach the new Lm/Lr sex means to each DIA-NN Lp wide table.
match_lp() {
    local input=$1 output=$2
    "${PYTHON:-python3}" "$CODE_DIR/match_Lp_with_min3_LmLr_v5.py" \
        --lp-wide "$input" --lmlr "$LMLR" --cdspos-dir "$CDSPOS_DIR" --out "$output"
}
match_lp "$MS_ROOT/ISO1_sim_combined_MaxQuant_DIANN_Lp_scatter/intermediate/dmel_iso_Lp_wide.tsv" "$INPUT_DIR/dmel_iso_Lp_wide.counter_rule_LmLr.tsv"
match_lp "$MS_ROOT/DIA_NN_analysis/dmel_Lp_union.combined_sex_Lp_wide.tsv" "$INPUT_DIR/dmel_Lp_union.counter_rule_LmLr.tsv"
match_lp "$MS_ROOT/DIA_NN_analysis/dsim_Lp_union.combined_sex_Lp_wide.tsv" "$INPUT_DIR/dsim_Lp_union.counter_rule_LmLr.tsv"
match_lp "$MS_ROOT/DIA_NN_analysis/dsec_Lp_union.combined_sex_Lp_wide.tsv" "$INPUT_DIR/dsec_Lp_union.counter_rule_LmLr.tsv"
match_lp "$MS_ROOT/DIA_NN_analysis/dere_Lp_union.combined_sex_Lp_wide.tsv" "$INPUT_DIR/dere_Lp_union.counter_rule_LmLr.tsv"

# 3a. ISO1/sim: DIA-NN, MaxQuant, and equal-weight combined Lp plots.
LMLR_TABLE="$LMLR" \
    DIA_ISO="$INPUT_DIR/dmel_iso_Lp_wide.counter_rule_LmLr.tsv" \
    DIA_SIM="$INPUT_DIR/dsim_Lp_union.counter_rule_LmLr.tsv" \
    "${RSCRIPT:-Rscript}" "$CODE_DIR/plot_ISO1_sim_combined_MaxQuant_DIANN_Lp_vs_LmLr.R" "$OUTPUT_ROOT/01_ISO1_sim_union"

# 3b--d. Ipos, Dmel genotype violin, and sex-stratified Lm/Lr scatter figures.
"${RSCRIPT:-Rscript}" "$CODE_DIR/plot_Ipos_Lp_Lm_3sex_2methods_A4_v4.R" \
    "$OUTPUT_ROOT/02_Ipos_Lp_Lm_3sex" "$CDSPOS_DIR" \
    "$INPUT_DIR/dmel_Lp_union.counter_rule_LmLr.tsv" "$INPUT_DIR/dsim_Lp_union.counter_rule_LmLr.tsv" \
    "$INPUT_DIR/dsec_Lp_union.counter_rule_LmLr.tsv" "$INPUT_DIR/dere_Lp_union.counter_rule_LmLr.tsv"
"${RSCRIPT:-Rscript}" "$CODE_DIR/plot_Lp_ratio_difference_violin_paired_by_sex_SX31style_current_wide_AG_WBg_AS_A4quarter.R" \
    "$OUTPUT_ROOT/03_Dmel_AG_AS_WBg_violin" "$INPUT_DIR/dmel_Lp_union.counter_rule_LmLr.tsv"
"${RSCRIPT:-Rscript}" "$CODE_DIR/plot_male_Lp_vs_Lm_Lr_longestTx.R" \
    "$INPUT_DIR/dmel_Lp_union.counter_rule_LmLr.tsv" "$INPUT_DIR/dsim_Lp_union.counter_rule_LmLr.tsv" \
    "$INPUT_DIR/dsec_Lp_union.counter_rule_LmLr.tsv" "$INPUT_DIR/dere_Lp_union.counter_rule_LmLr.tsv" \
    "$OUTPUT_ROOT/04_Lp_Lm_Lr_scatter"
