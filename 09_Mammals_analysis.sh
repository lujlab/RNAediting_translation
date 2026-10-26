#!/usr/bin/env bash
# Configure section 00, then run the required numbered sections.
set -euo pipefail
# Human/mouse only; use existing NH-weighted codon counts.

# =============================================================================
# 00. Paths and software
# =============================================================================
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CODE_DIR="$PIPELINE_DIR/code/09"
cd "$CODE_DIR"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PLOT_PYTHON="${PLOT_PYTHON:-$PYTHON_BIN}"
CXX="${CXX:-g++}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$CODE_DIR/.cache/matplotlib}"
mkdir -p "$MPLCONFIGDIR"

# =============================================================================
# 01. Prepare reference and candidate inputs
# =============================================================================
"$PYTHON_BIN" prepare.py

# =============================================================================
# 02. Select the longest CDS per gene
# =============================================================================
# Section 01 has just regenerated the full reference/candidate inputs.
# For a selection-only rerun on unchanged inputs, omit --refresh-reference.
"$PYTHON_BIN" select_longest.py --refresh-reference

# =============================================================================
# 03. Compile codon-counting programs
# =============================================================================
"$CXX" -O3 -std=c++17 count.cpp -o count
"$CXX" -O3 -std=c++17 count_all_CDS.cpp -o count_all_CDS

# =============================================================================
# 04. Count candidate and all-CDS positions
# =============================================================================
"$PYTHON_BIN" run_counts.py
"$PYTHON_BIN" recount_all_CDS.py

# =============================================================================
# 05. Assemble Table S12 from existing counts
# =============================================================================
# If counts already exist, start here after section 00.
"$PYTHON_BIN" assemble_alt2.py

# =============================================================================
# 06. Draw Figure S28
# =============================================================================
"$PLOT_PYTHON" plot_scatter.py

# =============================================================================
# 07. Export Table S13
# =============================================================================
"$PLOT_PYTHON" TableS13/build_tableS13.py

