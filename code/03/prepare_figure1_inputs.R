#!/usr/bin/env Rscript
# Only arrange the current final tables into the filenames expected by original R.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_figure1_inputs.R FINAL_DIR DATA_DIR")
final <- normalizePath(args[1], mustWork = TRUE)
dest <- args[2]
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
master <- fread(file.path(final, "candidate.final.tsv"))
stopifnot(!anyDuplicated(master$handle))
groups <- c("ISO1", "WBG", "sim", "sec", "ere")
sets <- list()
manifest <- list()
for (g in groups) {
    src <- file.path(final, "candidate_lists_by_group", paste0(g, ".candidate.detected.annotated.tsv"))
    dt <- fread(src)
    stopifnot(all(c("handle", "type") %in% names(dt)), !anyDuplicated(dt$handle), all(dt$handle %in%
        master$handle))
    if (g == "sec")
        stopifnot("ISO1_type" %in% names(dt))
    sets[[g]] <- dt$handle
    out <- file.path(dest, paste0(g, ".annotation.tsv"))
    stopifnot(file.copy(src, out, overwrite = TRUE))
    manifest[[g]] <- data.table(group = g, input = src, output = out, rows = nrow(dt))
}
stopifnot(setequal(union(sets$ISO1, sets$WBG), master$handle))
src <- file.path(final, "heatmap_inputs", "candidate.editing_level.observed.long.tsv")
q <- fread(src)
stopifnot(all(c("group", "library", "handle", "editing_level") %in% names(q)))
stopifnot(!anyDuplicated(q[, .(group, library, handle)]), all(is.finite(q$editing_level)),
    all(q$editing_level >= 0 & q$editing_level <= 1))
out <- file.path(dest, basename(src))
stopifnot(file.copy(src, out, overwrite = TRUE))
manifest[["quant"]] <- data.table(group = "quant", input = src, output = out, rows = nrow(q))
fwrite(rbindlist(manifest), file.path(dest, "input_manifest.tsv"), sep = "\t")
shared <- Reduce(intersect, list(union(sets$ISO1, sets$WBG), sets$sim, sets$sec, sets$ere))
accepted_shared <- fread(file.path(final, "site_lists", "candidate.shared_Dmel_Dsim_Dsec_Dere.handles.txt"),
    header = FALSE)[[1]]
stopifnot(setequal(shared, accepted_shared))
fwrite(data.table(set = c("Dmel", "Dsim", "Dsec", "Dere", "shared"), n = c(nrow(master),
    length(sets$sim), length(sets$sec), length(sets$ere), length(shared))), file.path(dest,
    "accepted_set_counts.tsv"), sep = "\t")
cat("Prepared current Figure1 inputs:", nrow(master), "Dmel candidates;", length(shared),
    "shared sites\n")
