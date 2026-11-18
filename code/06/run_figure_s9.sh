#!/usr/bin/env bash
set -euo pipefail
FS9_CODE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
FS9_PIPELINE_DIR="$(cd -- "${FS9_CODE_DIR}/../.." && pwd)"
FS9_COUNTS_DIR="${FS9_COUNTS_DIR:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean/04_alledit_new/02_double_edit_and_MS/01_counts}"
FS9_OUT_DIR="${FS9_OUT_DIR:-${FS9_PIPELINE_DIR}/results/06/FigureS9}"
exec "${FS9_PYTHON:-python3}" "${FS9_CODE_DIR}/figure_s9_equalrep.py" \
    --catalog "${FS9_COUNTS_DIR}/double_edit_codons.tsv" \
    --counts "${FS9_COUNTS_DIR}/double_edit_allele_counts.by_sample.tsv" \
    --outdir "${FS9_OUT_DIR}" "$@"
