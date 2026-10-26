#!/usr/bin/env bash
# Reproduce Supplementary Figure 18 using embryo CDS sites (>=1 site/gene).
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR="$SCRIPT_DIR"
PIPELINE_DIR="$(cd -- "$ANALYSIS_DIR/../../.." && pwd)"
cd "$ANALYSIS_DIR"

RSCRIPT=${RSCRIPT:-Rscript}
# Pin this figure to the embryo CDS input and >=1 threshold, even if an old
# head-site configuration is exported in the calling shell.
export HALF_LIFE_EDITING_FILE="$ANALYSIS_DIR/inputs/published/Duan2017_Graveley2011_embryo_CDS_editing_sites_long.tsv"
export HALF_LIFE_MIN_EDITING_SITES=1
export HALF_LIFE_OUT_DIR="${HALF_LIFE_OUT_DIR:-$PIPELINE_DIR/output/duan2017_embryo_cds_min1}"
HALF_LIFE_PDF=${HALF_LIFE_PDF:-$PIPELINE_DIR/output/Supplementary Figure 18 Half-life CDS embroy sites min1.pdf}

command -v "$RSCRIPT" > /dev/null || {
    echo "Rscript not found: $RSCRIPT (activate your R environment or set RSCRIPT)." >&2
    exit 1
}
command -v gs > /dev/null || {
    echo "Ghostscript is required." >&2
    exit 1
}
[[ -s "$HALF_LIFE_EDITING_FILE" ]] || {
    echo "Missing embryo CDS input: $HALF_LIFE_EDITING_FILE" >&2
    exit 1
}
mkdir -p "$HALF_LIFE_OUT_DIR" "$(dirname -- "$HALF_LIFE_PDF")"

"$RSCRIPT" --vanilla run_all.R
"$RSCRIPT" --vanilla scripts/plot_suppfig18_duan2017_embryo_cds.R \
    "$HALF_LIFE_OUT_DIR" "$HALF_LIFE_PDF" "$HALF_LIFE_MIN_EDITING_SITES"
printf '[DONE] %s\n' "$HALF_LIFE_PDF"
