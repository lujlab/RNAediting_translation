#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript.")
base_dir <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
source(file.path(base_dir, "scripts", "figure_s18", "run_s18.R"))
data_dir <- if (length(args) >= 1L) args[[1L]] else Sys.getenv("HALF_LIFE_OUT_DIR", "output/duan2017_embryo_cds_min1")
output_pdf <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("HALF_LIFE_PDF",
    file.path(base_dir, "..", "..", "..", "output", "Supplementary_Figure_18_r6_12_transcripts_editable.pdf"))
if (length(args) >= 3L && args[[3L]] != "1") stop("Current S18 uses the >=1 CDS-site definition.")
if (!file.exists(file.path(data_dir, "plot_data.json"))) stop("Run run_half_life.sh first; current plot_data.json is missing.")
s18_render(base_dir, data_dir, output_pdf)
