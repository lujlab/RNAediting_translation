#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${R_PREFIX:-}" ]]; then export PATH="$R_PREFIX/bin:$PATH"; fi
pipeline_dir="$(cd -- "$here/../../.." && pwd)"
library_dir="${R05_LIBRARY:-$pipeline_dir/.runtime/05/R_library}"
export R_LIBS_USER="$library_dir${R_LIBS_USER:+:$R_LIBS_USER}"
mkdir -p "$library_dir"
: "${R_PACKAGE_SOURCE_DIR:?Provide source archives for statmod 1.5.2, locfit, and edgeR 4.8.2}"
for src in statmod_1.5.2.tar.gz locfit.tar.gz edgeR_4.8.2.tar.gz; do
    R CMD INSTALL -l "$library_dir" "$R_PACKAGE_SOURCE_DIR/$src"
done
Rscript --vanilla -e 'library(edgeR); library(data.table); library(ggplot2); sessionInfo()'
