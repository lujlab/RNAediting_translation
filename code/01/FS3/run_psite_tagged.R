#!/usr/bin/env Rscript
# Figure S3: one species per invocation; use existing transcriptome BAM PS/XF tags.
# Usage: Rscript run_psite_tagged.R samples.tsv annotation.tsv output_dir
# samples.tsv: sample, group, bam (absolute path); one biological replicate per row.
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(riboWaltz)
})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 3L) stop("Usage: run_psite_tagged.R samples.tsv annotation.tsv output_dir")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(script))
source(file.path(here, "ribowaltz_function_code.R"))
source(file.path(here, "ribowaltz_psite_bam_patch.R"))

manifest <- fread(argv[1])
stopifnot(all(c("sample", "group", "bam") %in% names(manifest)))
if (!nrow(manifest) || anyNA(manifest[, .(sample, group, bam)]) || anyDuplicated(manifest$sample) ||
    anyDuplicated(manifest$bam) || any(!nzchar(manifest$sample)) || any(!nzchar(manifest$group)) ||
    any(!file.exists(manifest$bam))) stop("Invalid sample manifest or missing BAM.")
annotation <- fread(argv[2])
if (anyDuplicated(annotation$transcript)) stop("Duplicate transcript IDs in annotation.")
out <- argv[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
groups <- split(manifest$sample, factor(manifest$group, levels = unique(manifest$group)))

# Preserve patch coordinate convention: SAM POS + PS; no new offset estimation.
# Exclude unmapped and reverse-strand reads, as in the forward-only main analysis.
reads <- read_psite_bams_as_ribowaltz(bams = manifest$bam, annotation = annotation, sample_names = manifest$sample,
    samtools = Sys.getenv("SAMTOOLS", "samtools"), exclude_flag = 20L, unique_only = FALSE,
    psite_shift = 0L, check_frame = TRUE)
for (nm in names(reads)) {
    if (!nrow(reads[[nm]]))
        stop("No usable reads: ", nm)
    cds <- reads[[nm]][psite_region == "cds"]
    if (!nrow(cds) || anyNA(cds$frame_bam) || any(!cds$frame_bam %in% 0:2))
        stop("Missing/invalid CDS XF tags: ", nm)
}

save_result <- function(result, prefix) {
    for (nm in names(result)) {
        x <- result[[nm]]
        safe <- gsub("[^A-Za-z0-9_.-]", "_", nm)
        path <- file.path(out, paste0(prefix, "_", safe))
        if (inherits(x, "ggplot")) {
            ggsave(paste0(path, ".pdf"), plot = x, width = 7, height = 4, device = "pdf",
                useDingbats = FALSE)
        } else if (is.data.frame(x)) {
            fwrite(x, paste0(path, ".tsv"), sep = "\t")
        }
    }
}

# A: package regional counting, using the PS-derived coordinates, not new offsets.
save_result(riboWaltz::region_psite(reads, annotation = annotation, sample = groups, multisamples = "average",
    plot_style = "dodge", cl = 85, colour = c("#333f50", "gray70", "#39827c")), "S3A")
# B: explicitly read XF; do not silently substitute calculated frame.
save_result(frame_psite_bam(reads, annotation = annotation, sample = groups, multisamples = "average",
    plot_style = "facet", region = "cds", cl = 100, frame_source = "XF", colour = c("#333f50",
        "#39827c")), "S3B")
# C: the supplied custom metaprofile function, preserving its normalization.
save_result(metaprofile_psite_bam(reads, annotation = annotation, sample = groups, multisamples = "average",
    utr5l = 20, cdsl = 40, utr3l = 20, cl = 100, colour = "#333f50"), "S3C")

fwrite(manifest, file.path(out, "samples_used.tsv"), sep = "\t")
fwrite(manifest[, .(n = .N), by = group], file.path(out, "biological_replicate_counts.tsv"),
    sep = "\t")
invisible(NULL)
