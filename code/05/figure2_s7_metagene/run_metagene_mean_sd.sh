#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
METAGENE_CODE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$METAGENE_CODE/../../.." && pwd)}"
RESCUE_ROOT="${RESCUE_ROOT:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean}"
OUT_DIR="${OUT_DIR:-${PIPELINE_DIR}/results/05_stalling}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
METAGENE_MEAN_DIR="${METAGENE_MEAN_DIR:-${RESCUE_ROOT}/09_Ipos/species_sex_mean_gmRC_tables}"
METAGENE_OUTDIR="${METAGENE_OUTDIR:-${OUT_DIR}/metagene/profiles_existingMean_mean_sd}"
METAGENE_ALLEDIT_DIR="${METAGENE_ALLEDIT_DIR:-${RESCUE_ROOT}/04_alledit_new/results}"
METAGENE_CDSPOS_DIR="${METAGENE_CDSPOS_DIR:-${EDIT_ROOT:-${RESCUE_ROOT}/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed}}"
METAGENE_REFERENCE_CONFIG="${METAGENE_REFERENCE_CONFIG:-${SAMPLE_CONFIG_DIR:-${PIPELINE_DIR}/code/05/config}/metagene_references.tsv}"
# Same existing means and mean/SD algorithm as sections 13-14 of the main Bash.
exec env METAGENE_SOURCE_SCRIPT="$METAGENE_CODE/All_metagene_mean_existingMean_any1libLm02_mean_sd.R" \
    METAGENE_MEAN_DIR="$METAGENE_MEAN_DIR" METAGENE_OUTDIR="$METAGENE_OUTDIR" \
    METAGENE_ALLEDIT_DIR="$METAGENE_ALLEDIT_DIR" METAGENE_CDSPOS_DIR="$METAGENE_CDSPOS_DIR" \
    METAGENE_REFERENCE_CONFIG="$METAGENE_REFERENCE_CONFIG" METAGENE_CONTROL_CENTER_A_ONLY=true \
    "$RSCRIPT_BIN" --vanilla "$METAGENE_CODE/run_metagene_counter_rule_streamed.R"
