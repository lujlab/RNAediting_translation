#!/usr/bin/env bash
# Final publication workflow for corrected mRC/gmRC and every downstream gmRC panel.
# This supersedes the historical gmRC/plotting sections in 05_Stalling_analysis.sh.
set -euo pipefail
export LC_ALL=C

STAGE="${1:-all}"
case "$STAGE" in
    all | compute | profiles | fs10 | models | render) ;;
    *)
        echo "Usage: $0 [all|compute|profiles|fs10|models|render]" >&2
        exit 2
        ;;
esac

GMRC_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
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

run_stage() {
    [[ "$STAGE" == all || "$STAGE" == "$1" ]]
}

require_file() {
    [[ -s "$1" ]] || {
        echo "Missing required file: $1" >&2
        exit 1
    }
}

render_python() {
    env PYTHONPATH="${RENDER_PYTHONPATH}${PYTHONPATH:+:${PYTHONPATH}}" "$RENDER_PYTHON" "$@"
}

if run_stage compute; then
    echo "[compute] corrected single-library mRC/gmRC and group means"
    "$PYTHON_BIN" "$GMRC_ROOT/code/recompute.py" prepare
    "$PYTHON_BIN" "$GMRC_ROOT/code/recompute.py" run --jobs "$JOBS"
    "$PYTHON_BIN" "$GMRC_ROOT/code/prepare_figure_manifest.py"
    "$PYTHON_BIN" "$GMRC_ROOT/code/prepare_glm_inputs.py"
fi

if run_stage profiles; then
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

fi

# FS10: MAIN gmRC (RPF50 among RNA>=5 internal codons), sum of A-site medians >=3*n_effective.
# Keep sexes separate; model gmRC ~ AA + sex + preference. No additional editing/m6A mask.
# Report one-sided codon-clustered CR1 P values without multiple-testing correction.
if run_stage fs10 || [[ "$STAGE" == profiles ]]; then
    RSCRIPT="$RSCRIPT_BIN" bash "$WORKSPACE_ROOT/code/05/figure_s10_codon_preference/run.sh"
fi

if run_stage models; then
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
fi

if run_stage render; then
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
fi



echo "Publication gmRC workflow completed: stage=$STAGE"
