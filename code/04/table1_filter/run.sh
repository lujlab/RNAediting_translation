#!/usr/bin/env bash
# Generate Table1 filter scatter data and individual panels.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$HERE/../../.." && pwd)"
OUT="${TABLE1_FILTER_OUT:-$ROOT/results/04_RNA_Ribo_alledit_Table1/04_Table1_filter}"
mkdir -p "$OUT"
"${RSCRIPT:-Rscript}" --vanilla "$HERE/run_analysis.R" "$ROOT" "$OUT" > "$OUT/analysis.log" 2>&1
printf 'Table1 filter outputs: %s\n' "$OUT"
