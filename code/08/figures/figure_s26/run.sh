#!/usr/bin/env bash
# Reproduce only the CDS editing-group analysis and four panels of Figure S26.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$HERE/../../../.." && pwd)"
BASE="${S26_BASE:-$ROOT/results/08_editing_expression}"
OUT="${S26_OUT:-$ROOT/results/08_publication/S26}"
RSCRIPT="${RSCRIPT:-Rscript}"
export OPENBLAS_NUM_THREADS="${ANALYSIS_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${ANALYSIS_BLAS_THREADS:-1}"
COUNTS="${S26_COUNTS:-$BASE/inputs/WBg_AdarG_AdarS_featureCounts_raw_merged.tsv}"
CANDIDATES="${S26_CANDIDATES:-$BASE/editing_within_sex/AG_AS_WBG_sample_sites.long.standardized.tsv}"
TARGETS="${S26_TARGETS:-$BASE/editing_within_sex/NB_BetaBinomial_overlap_all_clusters_Lm_sites_for_AdarG_AdarS_delta.tsv}"
CDS_MAP="${S26_CDS_MAP:-$BASE/RNA_LRT_degPatterns/RNA_editing_gene_DE_cluster_enrichment/alledit_CDS.site_to_gene.tsv}"
STAGE="${1:-all}"
case "$STAGE" in all | fit | groups | plot) ;; *)
    echo 'Usage: bash run.sh [all|fit|groups|plot]' >&2
    exit 2
    ;;
esac
mkdir -p "$OUT/logs"
if [[ "$STAGE" == all || "$STAGE" == fit ]]; then
    "$RSCRIPT" --vanilla "$HERE/01_deseq2.R" "$COUNTS" "$OUT/deseq2" > "$OUT/logs/01_deseq2.log" 2>&1
fi
if [[ "$STAGE" == all || "$STAGE" == groups ]]; then
    "$RSCRIPT" --vanilla "$HERE/02_groups.R" "$CANDIDATES" "$TARGETS" "$CDS_MAP" "$OUT/groups" > "$OUT/logs/02_groups.log" 2>&1
fi
if [[ "$STAGE" == all || "$STAGE" == plot ]]; then
    "$RSCRIPT" --vanilla "$HERE/03_plot.R" "$OUT/deseq2" "$OUT/groups" "$OUT" > "$OUT/logs/03_plot.log" 2>&1
fi
printf 'Figure S26 output: %s\n' "$OUT"
