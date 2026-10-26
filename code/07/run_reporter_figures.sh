#!/usr/bin/env bash
# Figures 4C-F and S15D-G from the bundled experimental observations.
set -euo pipefail
reporter_code="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
reporter_root="$(cd -- "$reporter_code/../.." && pwd)"
cd "$reporter_root"
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
export REPORTER_OUT_DIR="${1:-$reporter_root/results/07}"
mkdir -p "$REPORTER_OUT_DIR"
REPORTER_OUT_DIR="$(cd -- "$REPORTER_OUT_DIR" && pwd)"
export REPORTER_WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/reporter.XXXXXXXX")"
trap 'rm -rf -- "$REPORTER_WORK_DIR"' EXIT
LUC_PLOT_OUTPUT_DIR="$REPORTER_OUT_DIR" "$RSCRIPT" --vanilla "$reporter_code/figure4_CEF.R"
"$RSCRIPT" --vanilla "$reporter_code/standalone_panel_data.R"
"$PYTHON_BIN" "$reporter_code/draw_standalone_panels.py"
