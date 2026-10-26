#!/usr/bin/env Rscript
# Fit the sex-specific RNA models used by Figure S26.
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 01_deseq2.R merged_featurecounts.tsv output_dir")
input <- args[1]
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- fread(input)
samples <- unlist(lapply(c("WBg", "AdarG", "AdarS"), function(g) paste0(g, "_", rep(c("F",
    "M"), each = 3), rep(1:3, 2))))
stopifnot(all(c("Geneid", samples) %in% names(x)), !anyDuplicated(x$Geneid))
m <- as.matrix(x[, ..samples])
storage.mode(m) <- "double"
stopifnot(all(is.finite(m)), all(m >= 0))
rownames(m) <- x$Geneid
metadata <- data.frame(sample = samples, genotype = sub("_.*", "", samples), sex = sub(".*_([FM])[0-9]+$",
    "\\1", samples), row.names = samples)
fwrite(as.data.table(metadata), file.path(out, "samples.tsv"), sep = "\t")
for (sx in c("F", "M")) {
    meta <- metadata[metadata$sex == sx, , drop = FALSE]
    mat <- round(m[, meta$sample, drop = FALSE])
    storage.mode(mat) <- "integer"
    mat <- mat[rowSums(mat) >= 100, , drop = FALSE]
    meta$genotype <- factor(meta$genotype, levels = c("WBg", "AdarG", "AdarS"))
    dds <- DESeqDataSetFromMatrix(mat, meta, design = ~genotype)
    dds <- DESeq(dds, test = "Wald", fitType = "local", quiet = TRUE)
    for (mut in c("AdarG", "AdarS")) {
        r <- as.data.table(as.data.frame(results(dds, contrast = c("genotype", mut, "WBg"),
            alpha = 0.01)), keep.rownames = "gene_id")
        r[, `:=`(sex = sx, mutant = mut, denominator = "WBg")]
        fwrite(r, file.path(out, paste0("RNA_", sx, "_", mut, "_vs_WBg.tsv")), sep = "\t")
    }
}
invisible(NULL)
