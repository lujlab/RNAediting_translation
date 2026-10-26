#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# RNA-editing levels, genotype-dependent gene expression, and Figure 5 panels.
# Configure inputs in section 00, then execute the numbered sections in order.
# Biological replicates enter the models separately; means are used for plots.
# To reuse existing counts, configure FEATURECOUNTS_DIR and skip section 00a.
# Code layout: setup/ for configuration, analysis/ for models, figures/ for plots.

# ============================================================================
# 00. Paths and analysis settings
# ============================================================================
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
source "$PIPELINE_DIR/code/08/setup/paths.sh"
mkdir -p "$OUT_DIR/inputs" "$OUT_DIR/logs" "$FIGURE_DIR"
export OPENBLAS_NUM_THREADS="${ANALYSIS_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${ANALYSIS_BLAS_THREADS:-1}"

# ============================================================================
# 00a. Count gene-level RNA-seq fragments from the BAMs produced by workflow 02
# ============================================================================
# Reverse-stranded paired-end libraries; three female and three male replicates
# per genotype. Preserve raw fragment counts for expression analysis.
[[ -f "$RNA_GTF" ]] || {
    echo "Missing GTF: $RNA_GTF" >&2
    exit 1
}
[[ -s "$RNA_SAMPLE_TABLE" ]] || {
    echo "Missing sample table: $RNA_SAMPLE_TABLE" >&2
    exit 1
}
mkdir -p "$FEATURECOUNTS_DIR"
for genotype in WBG AG AS; do
    mapfile -t bams < <(awk -F '\t' -v g="$genotype" -v root="$RNA_ALIGNMENT_ROOT" \
        'NR > 1 && $1 == g {print root "/" $3}' "$RNA_SAMPLE_TABLE")
    [[ ${#bams[@]} -eq 6 ]] || {
        echo "Expected six BAMs for $genotype" >&2
        exit 1
    }
    for bam in "${bams[@]}"; do
        [[ -f "$bam" ]] || {
            echo "Missing BAM: $bam" >&2
            exit 1
        }
    done
    "$FEATURECOUNTS_BIN" \
        -T "$FEATURECOUNTS_THREADS" \
        -p --countReadPairs \
        -s 2 \
        -t exon \
        -g gene_id \
        -Q 10 \
        -a "$RNA_GTF" \
        -o "$FEATURECOUNTS_DIR/${genotype}_STAR_featureCounts_new" \
        "${bams[@]}" > "$OUT_DIR/logs/00a_${genotype}.featureCounts.log" 2>&1
done

# ============================================================================
# 01. Validate dependencies and input files
# ============================================================================
# Install missing packages once, if required:
# "$RSCRIPT_BIN" "$CODE_DIR/setup/install_packages.R"
"$RSCRIPT_BIN" "$CODE_DIR/setup/preflight.R" > "$OUT_DIR/logs/01_preflight.log" 2>&1

# ============================================================================
# 02. Merge raw featureCounts fragment counts and export sample metadata
# ============================================================================
# These files were counted with --countReadPairs. Do not divide them by two.
"$RSCRIPT_BIN" "$CODE_DIR/analysis/expression/01_merge_featurecounts.R" > "$OUT_DIR/logs/02_counts.log" 2>&1

# ============================================================================
# 03. Test genotype effects on editing within each sex; cluster and plot
# ============================================================================
# NB and beta-binomial LRTs; >=2 replicates per genotype, depth >=10.
# Apply model-specific BH FDR <0.05, mean editing >0.10 and pairwise difference >0.05.
"$RSCRIPT_BIN" "$CODE_DIR/analysis/editing/01_within_sex_NB_BB.R" \
    "$EDIT_INPUT_DIR" "$EDIT_OUT_DIR" > "$OUT_DIR/logs/03_editing.log" 2>&1

# ============================================================================
# 03b. Enrich GO terms for genes with altered editing levels
# ============================================================================
"$RSCRIPT_BIN" "$CODE_DIR/analysis/editing/03_editing_GO.R" "$EDIT_OUT_DIR" \
    > "$OUT_DIR/logs/03b_editing_GO.log" 2>&1

# ============================================================================
# 04. Test editing genotype-by-sex interactions
# ============================================================================
# Compare additive and interaction models using the complete 3-by-2 design.
"$RSCRIPT_BIN" "$CODE_DIR/analysis/editing/02_genotype_sex_interaction.R" \
    "$EDIT_INPUT_DIR" "$INTERACTION_OUT_DIR" > "$OUT_DIR/logs/04_interaction.log" 2>&1

# ============================================================================
# 05. Test gene expression by DESeq2 LRT; generate PCA, clusters, and cluster GO
# ============================================================================
# Within each sex: full model ~genotype; reduced model ~1; BH FDR <0.01.
# degPatterns produces the six-cluster publication panels.
"$RSCRIPT_BIN" "$CODE_DIR/analysis/expression/02_lrt_${CLUSTER_METHOD}.R" > "$OUT_DIR/logs/05_rna_lrt.log" 2>&1

# ============================================================================
# 05b. Enrich GO terms for the expression clusters
# ============================================================================
"$RSCRIPT_BIN" "$CODE_DIR/analysis/expression/06_cluster_GO_${CLUSTER_METHOD}.R" > "$OUT_DIR/logs/05b_cluster_GO.log" 2>&1

# ============================================================================
# 06. Run DESeq2 Wald contrasts, UpSet plots, effect-size plots, and directional GO
# ============================================================================
# AdarG vs WBg, AdarS vs WBg, and AdarS vs AdarG are tested separately for F and M.
"$RSCRIPT_BIN" "$CODE_DIR/analysis/expression/03_pairwise_DE_and_GO.R" > "$OUT_DIR/logs/06_pairwise.log" 2>&1

# ============================================================================
# 07. Test enrichment of ISO1-defined edited genes among DE genes and clusters
# ============================================================================
"$RSCRIPT_BIN" "$CODE_DIR/analysis/expression/04_edited_gene_enrichment_${CLUSTER_METHOD}.R" \
    > "$OUT_DIR/logs/07_edited_gene_enrichment.log" 2>&1

# ============================================================================
# 08. Plot VST expression for all genes entered into each sex-specific LRT
# ============================================================================
"$RSCRIPT_BIN" "$CODE_DIR/analysis/expression/05_lrt_heatmaps.R" > "$OUT_DIR/logs/08_heatmaps.log" 2>&1

# ============================================================================
# 09. Plot editing-associated protein-output bias (Lp minus Lm)
# ============================================================================
# Genotype comparisons are paired by site, with BH correction within each sex.
"$RSCRIPT_BIN" "$CODE_DIR/figures/plot_protein_output_bias.R" \
    "$PROTEIN_OUT_DIR" "$PROTEIN_INPUT" difference > "$OUT_DIR/logs/09_protein_bias.log" 2>&1

# ============================================================================
# 10. Plot ADAR auto-editing from A-site alledit counts
# ============================================================================
"$RSCRIPT_BIN" "$CODE_DIR/figures/plot_adar_autosite.R" \
    "$ADAR_ALL_EDIT_DIR" "$ADAR_OUT_DIR" > "$OUT_DIR/logs/10_adar_autosite.log" 2>&1

# ============================================================================
# 11. Figure S26: CDS editing groups versus DESeq2 model log2 fold changes
# ============================================================================
RSCRIPT="$RSCRIPT_BIN" S26_BASE="$OUT_DIR" \
    S26_CANDIDATES="$EDIT_OUT_DIR/AG_AS_WBG_sample_sites.long.standardized.tsv" \
    S26_TARGETS="$EDIT_OUT_DIR/NB_BetaBinomial_overlap_all_clusters_Lm_sites_for_AdarG_AdarS_delta.tsv" \
    S26_CDS_MAP="$RNA_OUT_DIR/RNA_editing_gene_DE_cluster_enrichment/alledit_CDS.site_to_gene.tsv" \
    bash "$CODE_DIR/figures/figure_s26/run.sh" > "$OUT_DIR/logs/11_S26.log" 2>&1
