#!/usr/bin/env bash
# Configure the external S17/S18/S27 inputs before running all sections.
set -euo pipefail

# 00. Paths and software
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PIPELINE_DIR"
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OUT="${PUBLICATION_07_OUT:-${PUBLICATION_09_OUT:-$PIPELINE_DIR/results/07}}"
mkdir -p "$OUT"
OUT="$(cd -- "$OUT" && pwd)"

# 01. Figures 4C-F and S15D-G: bundled experimental inputs
RSCRIPT="$RSCRIPT" PYTHON_BIN="$PYTHON_BIN" bash code/07/run_reporter_figures.sh "$OUT"

# 02. Figure S16: bundled 5/30/60-minute qPCR measurements
RSCRIPT="$RSCRIPT" PYTHON_BIN="$PYTHON_BIN" bash code/07/run_S16.sh "$OUT"

# Intermediate tables are consumed by subsequent plotting steps, then removed.
rna_work_dir="$(mktemp -d "${TMPDIR:-/tmp}/rna_figures.XXXXXXXX")"
trap 'rm -rf -- "$rna_work_dir"' EXIT

# 03. Figure S17: nascent and mature RNA
S17_OUT_DIR="$rna_work_dir/S17" "$RSCRIPT" --vanilla code/07/figure_s17_nascent_mature.R
cp "$rna_work_dir/S17/FigureS17_Nascent_Mature_rawP.pdf" "$OUT/FigureS17.pdf"

# 04. Figure S18: RNA half-life
RSCRIPT="$RSCRIPT" HALF_LIFE_OUT_DIR="$rna_work_dir/S18" HALF_LIFE_PDF="$OUT/FigureS18.pdf" \
    bash code/07/half_life/run_half_life.sh

# 05. Figure S27: protein/RNA ratio panel
"$RSCRIPT" --vanilla code/07/figure_s27_protein_rna.R --pseudocount=0 "--out-dir=$rna_work_dir/S27"
cp "$rna_work_dir/S27/FigureS27_CDSge2_protein_RNA_ratio_only.pdf" "$OUT/FigureS27.pdf"
