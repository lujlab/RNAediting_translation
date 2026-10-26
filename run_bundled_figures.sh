#!/usr/bin/env bash
# Reproduce the figures whose experimental inputs are included in this package.
set -euo pipefail
bundle_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$bundle_root"
export RSCRIPT="${RSCRIPT:-Rscript}"
export PYTHON_BIN="${PYTHON_BIN:-python3}"
bundle_out="${1:-$bundle_root/results}"
mkdir -p "$bundle_out"
bundle_out="$(cd -- "$bundle_out" && pwd)"
bash code/07/run_reporter_figures.sh "$bundle_out"
bash code/07/run_S16.sh "$bundle_out"
PHENOTYPE_OUT="$bundle_out" bash 10_Phenotype_analysis.sh
