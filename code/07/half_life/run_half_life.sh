#!/usr/bin/env bash
set -euo pipefail
figure_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
figure_rscript="${RSCRIPT:-Rscript}"
command -v "$figure_rscript" >/dev/null || {
    echo "Rscript not found; activate the R environment or set RSCRIPT." >&2
    exit 1
}
export HALF_LIFE_EDITING_FILE="$figure_script_dir/inputs/published/Duan2017_Graveley2011_embryo_CDS_editing_sites_long.tsv"
export HALF_LIFE_MIN_EDITING_SITES=1
"$figure_rscript" --vanilla "$figure_script_dir/run_all.R"
