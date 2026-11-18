#!/usr/bin/env bash
# Configure section 00, then run the required numbered sections.
set -euo pipefail

# =============================================================================
# 00. Paths and software
# =============================================================================
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
RSCRIPT="${RSCRIPT:-Rscript}"
PHENOTYPE_DATA="${PHENOTYPE_DATA:-$PIPELINE_DIR/data/10/phenotype_data.tsv}"
PHENOTYPE_OUT="${PHENOTYPE_OUT:-$PIPELINE_DIR/results/10}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

# =============================================================================
# 01. Figure 6: phenotype plots and statistical tests
# =============================================================================
# Sex-specific Kruskal-Wallis and Dunn tests with BH correction.
"$RSCRIPT" --vanilla "$PIPELINE_DIR/code/10/figure6.R" "$PHENOTYPE_DATA" "$PHENOTYPE_OUT"
