#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args) >= 1L) args[[1L]] else "../results"
manifest_path <- if (length(args) >= 2L) args[[2L]] else {
    "../../stalling_master_tables/master_table_manifest.tsv"
}
summary_dir <- file.path(result_dir, "run_summaries")
summary_files <- list.files(summary_dir, pattern = "\\.tsv$", full.names = TRUE)
if (!length(summary_files)) stop("No per-run summaries found in ", summary_dir, call. = FALSE)

summaries <- rbindlist(lapply(summary_files, fread), fill = TRUE)
if (file.exists(manifest_path)) {
    manifest <- fread(manifest_path, integer64 = "double")
    manifest <- unique(manifest[, .(group_id, species, genotype, sex, flank_size)])
    summaries <- merge(summaries, manifest, by = c("group_id", "flank_size"), all.x = TRUE)
}
setcolorder(summaries, c("group_id", "species", "genotype", "sex", "flank_size", setdiff(names(summaries),
    c("group_id", "species", "genotype", "sex", "flank_size"))))
setorder(summaries, species, genotype, sex, flank_size)
fwrite(summaries, file.path(result_dir, "edgeR_run_summary.tsv"), sep = "\t")

completed <- summaries[status %in% c("completed", "existing") & !is.na(n_tested) & n_tested >
    0]
if (nrow(completed)) {
    completed[, `:=`(peak_ratio = n_significant/n_tested)]
}
fwrite(completed, file.path(result_dir, "overall_peak_by_library.tsv"), sep = "\t")

totals <- completed[, .(libraries = .N, tested_codons_sum = sum(n_tested), peak_sum = sum(n_significant),
    peak_ratio_over_all_library_sites = sum(n_significant)/sum(n_tested)), by = flank_size]
fwrite(totals, file.path(result_dir, "overall_peak_totals.tsv"), sep = "\t")

print(summaries)
print(totals)
