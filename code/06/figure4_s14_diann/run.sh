#!/usr/bin/env bash
set -euo pipefail
AUDIT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
"${PYTHON}" "${AUDIT_DIR}/build_data.py"
"${PYTHON}" "${AUDIT_DIR}/draw_figures.py"
