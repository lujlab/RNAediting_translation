#!/usr/bin/env bash
# Selected 1016 table: separate any-edited MS support and target-codon II support.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${TABLE_S10_INPUTS:?Set TABLE_S10_INPUTS to the upstream codon/evidence tables}"
: "${DIA_NN_ANALYSIS_DIR:?Set DIA_NN_ANALYSIS_DIR to the four full 100-aa DIA-NN intensity tables}"
export TABLE_S10_INPUTS DIA_NN_ANALYSIS_DIR
"${PYTHON_BIN:-python3}" "$HERE/table_s10_multiedit/build_table_s10.py"
"${PYTHON_BIN:-python3}" "$HERE/table_s10_multiedit/render_table.py"
