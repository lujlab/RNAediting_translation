source(file.path(Sys.getenv("CODE_DIR"), "analysis", "common.R"))
dir.create(dirname(input_rds), recursive = TRUE, showWarnings = FALSE)
# merge counts files
suppressPackageStartupMessages({
    library(data.table)
})

# ============================================================
# 1. Input files
# ============================================================


file_WBg <- file.path(base_dir, "WBG_STAR_featureCounts_new")
file_AdarG <- file.path(base_dir, "AG_STAR_featureCounts_new")
file_AdarS <- file.path(base_dir, "AS_STAR_featureCounts_new")

# Output
out_file <- file.path(
    dirname(input_rds),
    "WBg_AdarG_AdarS_featureCounts_raw_merged.tsv"
)

# ============================================================
# 2. Desired sample order
# ============================================================

sampleNames <- c(
    "WBg_F1", "WBg_F2", "WBg_F3",
    "WBg_M1", "WBg_M2", "WBg_M3",
    "AdarG_F1", "AdarG_F2", "AdarG_F3",
    "AdarG_M1", "AdarG_M2", "AdarG_M3",
    "AdarS_F1", "AdarS_F2", "AdarS_F3",
    "AdarS_M1", "AdarS_M2", "AdarS_M3"
)

sample_WBg <- sampleNames[1:6]
sample_AdarG <- sampleNames[7:12]
sample_AdarS <- sampleNames[13:18]

# featureCounts annotation columns
info_cols <- c(
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
)

# ============================================================
# 3. Function: read one featureCounts file
# ============================================================

read_featurecounts <- function(file, sample_names) {
    cat("[INFO] Reading:", file, "\n")

    if (!file.exists(file)) {
        stop("Input file does not exist: ", file)
    }

    # skip featureCounts comment/header command line
    dt <- fread(
        file,
        skip = "Geneid",
        header = TRUE,
        sep = "\t",
        check.names = FALSE
    )

    cat(
        "[INFO] Rows:", nrow(dt),
        " Columns:", ncol(dt), "\n"
    )

    # Check basic structure
    missing_info <- setdiff(info_cols, names(dt))

    if (length(missing_info) > 0) {
        stop(
            "Missing annotation columns in ", file, ": ",
            paste(missing_info, collapse = ", ")
        )
    }

    count_cols <- setdiff(names(dt), info_cols)

    if (length(count_cols) != length(sample_names)) {
        stop(
            "Expected ", length(sample_names),
            " sample columns in ", file,
            ", but found ", length(count_cols),
            "\nDetected count columns:\n",
            paste(count_cols, collapse = "\n")
        )
    }

    # Resolve sample identities from BAM names; column order is not an identity.
    bam_names <- basename(count_cols)
    group_code <- sub("[-_].*$", "", bam_names)
    group_map <- c(WBG = "WBg", AG = "AdarG", AS = "AdarS", WBg = "WBg", AdarG = "AdarG", AdarS = "AdarS")
    suffix <- sub("^.*?[-_]([FM][123])(?:[-_.].*|Aligned.*|$)", "\\1", bam_names, perl = TRUE)
    actual_names <- paste(unname(group_map[group_code]), suffix, sep = "_")
    if (anyNA(group_map[group_code]) || anyDuplicated(actual_names) ||
        !setequal(actual_names, sample_names)) {
        stop("Cannot map featureCounts BAM columns to expected samples: ", paste(bam_names, collapse = ", "))
    }
    setnames(dt, old = count_cols, new = actual_names)
    setcolorder(dt, c(info_cols, sample_names))

    # Check duplicate genes
    if (anyDuplicated(dt$Geneid)) {
        dup <- unique(dt$Geneid[duplicated(dt$Geneid)])
        stop(
            "Duplicated Geneid detected in ", file,
            ". Example: ", paste(head(dup, 10), collapse = ", ")
        )
    }

    return(dt)
}

# ============================================================
# 4. Read three genotypes
# ============================================================

WBg <- read_featurecounts(
    file_WBg,
    sample_WBg
)

AdarG <- read_featurecounts(
    file_AdarG,
    sample_AdarG
)

AdarS <- read_featurecounts(
    file_AdarS,
    sample_AdarS
)

# ============================================================
# 5. Check gene sets
# ============================================================

cat("\n[INFO] Gene number:\n")
cat("  WBg   :", nrow(WBg), "\n")
cat("  AdarG :", nrow(AdarG), "\n")
cat("  AdarS :", nrow(AdarS), "\n")

genes_all <- Reduce(
    union,
    list(WBg$Geneid, AdarG$Geneid, AdarS$Geneid)
)

genes_common <- Reduce(
    intersect,
    list(WBg$Geneid, AdarG$Geneid, AdarS$Geneid)
)

cat("  Union :", length(genes_all), "\n")
cat("  Common:", length(genes_common), "\n")

if (
    length(genes_all) != nrow(WBg) ||
        length(genes_all) != nrow(AdarG) ||
        length(genes_all) != nrow(AdarS)
) {
    stop(
        "Gene sets are not identical among WBg / AdarG / AdarS.\n",
        "Please check whether all three files were generated using the same GTF."
    )
}

# ============================================================
# 6. Use WBg gene annotation as reference
# ============================================================

gene_info <- WBg[, ..info_cols]

# Match other genotype tables to WBg Geneid order
idx_AdarG <- match(gene_info$Geneid, AdarG$Geneid)
idx_AdarS <- match(gene_info$Geneid, AdarS$Geneid)

if (anyNA(idx_AdarG) || anyNA(idx_AdarS)) {
    stop("Failed to match Geneid among the three featureCounts files.")
}

AdarG <- AdarG[idx_AdarG]
AdarS <- AdarS[idx_AdarS]

# ============================================================
# 7. Check annotation consistency
# ============================================================

check_annotation <- function(ref, query, name) {
    for (cc in info_cols) {
        x <- as.character(ref[[cc]])
        y <- as.character(query[[cc]])

        bad <- which(
            !(is.na(x) & is.na(y)) &
                (is.na(x) | is.na(y) | x != y)
        )

        if (length(bad) > 0) {
            warning(
                "Annotation column '", cc,
                "' differs between WBg and ", name,
                " for ", length(bad), " genes."
            )
        }
    }
}

check_annotation(WBg, AdarG, "AdarG")
check_annotation(WBg, AdarS, "AdarS")

# ============================================================
# 8. Combine counts
# ============================================================

count_WBg <- WBg[, ..sample_WBg]
count_AdarG <- AdarG[, ..sample_AdarG]
count_AdarS <- AdarS[, ..sample_AdarS]

merged <- cbind(
    gene_info,
    count_WBg,
    count_AdarG,
    count_AdarS
)

# Force exact final column order
final_cols <- c(info_cols, sampleNames)

setcolorder(
    merged,
    final_cols
)

# ============================================================
# 9. Sanity checks
# ============================================================

stopifnot(
    identical(names(merged)[1:6], info_cols)
)

stopifnot(
    identical(names(merged)[7:24], sampleNames)
)

cat("\n[INFO] Final table:\n")
cat("  genes   :", nrow(merged), "\n")
cat("  columns :", ncol(merged), "\n")

cat("\n[INFO] Sample order:\n")
cat(
    paste(
        sprintf("%2d. %s", seq_along(sampleNames), sampleNames),
        collapse = "\n"
    ),
    "\n"
)

# ============================================================
# 10. Write output
# ============================================================

fwrite(
    merged,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

cat("\n[DONE] Output written to:\n")
cat(out_file, "\n")

rna <- as.matrix(merged[, ..sampleNames])
rownames(rna) <- merged$Geneid
if (any(!is.finite(rna)) || any(rna < 0) || any(rna != round(rna)))
    stop("DESeq2 requires finite, nonnegative integer raw counts.")
sample_info <- data.frame(
    sample = colnames(rna),
    genotype = sub("_.*$", "", colnames(rna)),
    sex = substr(sub("^[^_]+_", "", colnames(rna)), 1, 1),
    rep = sub("^[^_]+_[FM]", "", colnames(rna))
)
rownames(sample_info) <- sample_info$sample
sample_info$sex <- factor(sample_info$sex, levels = c("M", "F"))
sample_info$genotype <- factor(sample_info$genotype, levels = c("WBg", "AdarG", "AdarS"))
sample_info$group <- paste(sample_info$genotype, sample_info$sex, sep = "_")
if (anyNA(sample_info) || any(table(sample_info$genotype, sample_info$sex) != 3L))
    stop("Expected three biological replicates in every genotype-by-sex cell.")
rna_f <- rna[rowSums(rna) >= 100, , drop = FALSE]
saveRDS(list(rna_f = rna_f, sample_info = sample_info), input_rds)
fwrite(data.table(gene_id = rownames(rna), rna), file.path(dirname(input_rds), "RNA_gene_counts.tsv"), sep = "\t")
fwrite(as.data.table(sample_info), file.path(dirname(input_rds), "RNA_sample_metadata.tsv"), sep = "\t")
