#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
"${RSCRIPT:-${RSCRIPT_BIN:-Rscript}}" filter.R "$PWD"
"${RSCRIPT:-${RSCRIPT_BIN:-Rscript}}" model.R "$PWD"
"${PYTHON_BIN:-python3}" render.py
