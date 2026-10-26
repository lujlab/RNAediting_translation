#!/usr/bin/env bash
# Reproduce the selected S16 from the bundled 5/30/60-minute Ct observations.
# Usage: bash code/07/run_S16.sh [OUTPUT_DIRECTORY]
set -euo pipefail
s16_code_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
s16_root="$(cd -- "$s16_code_dir/../.." && pwd)"
s16_output_dir="${1:-$s16_root/results/S16}"
s16_input_dir="${S16_INPUT_DIR:-$s16_root/data/07}"
s16_python="${PYTHON_BIN:-python3}"
s16_rscript="${RSCRIPT:-Rscript}"
mkdir -p "$s16_output_dir"
s16_work_dir="$(mktemp -d "${TMPDIR:-/tmp}/s16.XXXXXXXX")"
trap 'rm -rf -- "$s16_work_dir"' EXIT
"$s16_python" "$s16_code_dir/extract_s16_qpcr.py" \
    --input-dir "$s16_input_dir" --output-dir "$s16_work_dir"
"$s16_rscript" --vanilla "$s16_code_dir/figure_S16_qPCR.R" \
    "$s16_work_dir/raw_Ct_5_30_60min.tsv" "$s16_output_dir"
