source(file.path(Sys.getenv("CODE_DIR"), "analysis", "common.R"))
list2env(readRDS(analysis_rds), envir = environment())
# ============================================================
# Heatmap of ALL genes entered into DESeq2 LRT
#
# IMPORTANT:
#   genes = ALL genes actually present in dds and tested by LRT
#   NOT only LRT-significant genes
#   NOT only degPatterns genes
#
# Female and Male are plotted separately.
# Rows    = biological samples
# Columns = genes
# Values  = VST-normalized expression
# ============================================================

suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
    library(gridExtra)
    library(grid)
    library(data.table)
})


# ============================================================
# 1. Output directory
# ============================================================

heatmap_out_dir <- file.path(
    count_out_dir,
    "LRT_all_gene_heatmap"
)

dir.create(
    heatmap_out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 2. Get EXACT genes entered into LRT
#
# rownames(dds) = genes actually used in DESeq2 LRT
# ============================================================

lrt_genes_F <- rownames(
    rna_F_res$dds
)

lrt_genes_M <- rownames(
    rna_M_res$dds
)


cat(
    "\nFemale genes entered into LRT:",
    length(lrt_genes_F),
    "\n"
)

cat(
    "Male genes entered into LRT:",
    length(lrt_genes_M),
    "\n"
)


# ============================================================
# 3. QC
#
# vst_mat should contain exactly the genes in dds
# ============================================================

stopifnot(
    setequal(
        lrt_genes_F,
        rownames(rna_F_res$vst_mat)
    )
)

stopifnot(
    setequal(
        lrt_genes_M,
        rownames(rna_M_res$vst_mat)
    )
)

# LRT result table should also have one row per tested gene
stopifnot(
    length(lrt_genes_F) ==
        nrow(rna_F_res$res_tbl)
)

stopifnot(
    length(lrt_genes_M) ==
        nrow(rna_M_res$res_tbl)
)


# ============================================================
# 4. Extract VST expression matrices
#
# Original:
#
# Heatmap:
# ============================================================

mat_F_gene_sample <- rna_F_res$vst_mat[
    lrt_genes_F, ,
    drop = FALSE
]

mat_M_gene_sample <- rna_M_res$vst_mat[
    lrt_genes_M, ,
    drop = FALSE
]


# transpose
mat_F <- t(mat_F_gene_sample)
mat_M <- t(mat_M_gene_sample)


# Do NOT remove genes here.
# Every LRT-tested gene must remain in the heatmap.
if (any(!is.finite(mat_F))) {
    stop(
        "Non-finite VST values detected in Female LRT matrix."
    )
}

if (any(!is.finite(mat_M))) {
    stop(
        "Non-finite VST values detected in Male LRT matrix."
    )
}


cat(
    "Female heatmap dimensions:",
    nrow(mat_F),
    "samples <U+00D7>",
    ncol(mat_F),
    "genes\n"
)

cat(
    "Male heatmap dimensions:",
    nrow(mat_M),
    "samples <U+00D7>",
    ncol(mat_M),
    "genes\n"
)


# ============================================================
# 5. Save exact LRT gene lists
# ============================================================

fwrite(
    data.table(
        gene_id = lrt_genes_F
    ),
    file.path(
        heatmap_out_dir,
        "RNA_F_genes_entered_into_LRT.tsv"
    ),
    sep = "\t"
)

fwrite(
    data.table(
        gene_id = lrt_genes_M
    ),
    file.path(
        heatmap_out_dir,
        "RNA_M_genes_entered_into_LRT.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 6. Pretty sample names
#
# WBg_F1   -> WT_rep1
# AdarG_F1 -> AdarG_rep1
# AdarS_F1 -> AdarS_rep1
# ============================================================

pretty_sample_name <- function(x) {
    genotype <- sub(
        "_.*$",
        "",
        x
    )

    rep_num <- sub(
        "^.*_[FM]([0-9]+)$",
        "\\1",
        x
    )

    genotype[
        genotype == "WBg"
    ] <- "WT"

    paste0(
        genotype,
        "_rep",
        rep_num
    )
}


labels_F <- pretty_sample_name(
    rownames(mat_F)
)

labels_M <- pretty_sample_name(
    rownames(mat_M)
)


# ============================================================
# 7. Sample annotation
# ============================================================

make_sample_annotation <- function(mat) {
    sample_names <- rownames(mat)

    genotype <- sub(
        "_.*$",
        "",
        sample_names
    )

    annotation <- data.frame(
        Genotype = factor(
            genotype,
            levels = c(
                "AdarG",
                "WBg",
                "AdarS"
            )
        )
    )

    rownames(annotation) <- sample_names

    annotation
}


annotation_F <- make_sample_annotation(
    mat_F
)

annotation_M <- make_sample_annotation(
    mat_M
)


# Similar to the reference figure:
# AdarG = blue
# WT    = yellow
# AdarS = pink

annotation_colors <- list(
    Genotype = c(
        "AdarG" = "#86D4EE",
        "WBg"   = "#F2EA00",
        "AdarS" = "#F3B5C8"
    )
)


# ============================================================
# 8. Shared Female/Male expression color scale
#
# Keep VST values.
# DO NOT perform row-wise or gene-wise Z-score.
# ============================================================

all_values <- c(
    as.numeric(mat_F),
    as.numeric(mat_M)
)

all_values <- all_values[
    is.finite(all_values)
]


# Use the same scale in Female and Male.
# 1%-99% clipping is ONLY for display to prevent a few
# extreme genes from compressing the color scale.
heat_low <- as.numeric(
    quantile(
        all_values,
        0.01,
        na.rm = TRUE
    )
)

heat_high <- as.numeric(
    quantile(
        all_values,
        0.99,
        na.rm = TRUE
    )
)


cat(
    "Shared VST display range:",
    round(heat_low, 3),
    "-",
    round(heat_high, 3),
    "\n"
)


# Clip only values used for visualization.
# Raw VST values are retained separately below.

mat_F_plot <- pmax(
    pmin(
        mat_F,
        heat_high
    ),
    heat_low
)

mat_M_plot <- pmax(
    pmin(
        mat_M,
        heat_high
    ),
    heat_low
)


heat_colors <- colorRampPalette(
    c(
        "#2166AC",
        "#67A9CF",
        "#D1E5F0",
        "#FFFFBF",
        "#FDAE61",
        "#D73027"
    )
)(100)


heat_breaks <- seq(
    heat_low,
    heat_high,
    length.out =
        length(heat_colors) + 1
)


# ============================================================
# 9. Female heatmap
#
# NOTE:
# Euclidean distance is used for genes so genes with
# nearly constant expression are NOT removed.
#
# Therefore ALL LRT-tested genes remain displayed.
# ============================================================

p_F <- pheatmap(
    mat_F_plot,
    color = heat_colors,
    breaks = heat_breaks,
    scale = "none",

    # samples
    cluster_rows = TRUE,

    # all genes
    cluster_cols = TRUE,
    clustering_distance_rows =
        "euclidean",
    clustering_distance_cols =
        "euclidean",
    clustering_method =
        "complete",
    annotation_row =
        annotation_F,
    annotation_colors =
        annotation_colors,
    labels_row =
        labels_F,
    show_rownames = TRUE,
    show_colnames = FALSE,
    border_color = NA,
    treeheight_row = 40,
    treeheight_col = 45,
    fontsize = 10,
    fontsize_row = 10,
    annotation_names_row = FALSE,
    main = paste0(
        "Female gene expression (",
        ncol(mat_F_plot),
        " LRT-tested genes)"
    ),
    silent = TRUE
)


# ============================================================
# 10. Male heatmap
# ============================================================

p_M <- pheatmap(
    mat_M_plot,
    color = heat_colors,
    breaks = heat_breaks,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows =
        "euclidean",
    clustering_distance_cols =
        "euclidean",
    clustering_method =
        "complete",
    annotation_row =
        annotation_M,
    annotation_colors =
        annotation_colors,
    labels_row =
        labels_M,
    show_rownames = TRUE,
    show_colnames = FALSE,
    border_color = NA,
    treeheight_row = 40,
    treeheight_col = 45,
    fontsize = 10,
    fontsize_row = 10,
    annotation_names_row = FALSE,
    main = paste0(
        "Male gene expression (",
        ncol(mat_M_plot),
        " LRT-tested genes)"
    ),
    silent = TRUE
)


# ============================================================
# 11. Save Female separately
# ============================================================

pdf(
    file.path(
        heatmap_out_dir,
        "RNA_F_ALL_LRT_tested_genes_heatmap.pdf"
    ),
    width = 10,
    height = 6.5,
    useDingbats = FALSE
)

grid.newpage()
grid.draw(
    p_F$gtable
)

dev.off()


# ============================================================
# 12. Save Male separately
# ============================================================

pdf(
    file.path(
        heatmap_out_dir,
        "RNA_M_ALL_LRT_tested_genes_heatmap.pdf"
    ),
    width = 10,
    height = 6.5,
    useDingbats = FALSE
)

grid.newpage()
grid.draw(
    p_M$gtable
)

dev.off()


# ============================================================
# 13. Female + Male combined PDF
# ============================================================

combined_pdf <- file.path(
    heatmap_out_dir,
    "RNA_FM_ALL_LRT_tested_genes_heatmap.pdf"
)


pdf(
    combined_pdf,
    width = 10,
    height = 12,
    useDingbats = FALSE
)


gridExtra::grid.arrange(
    gridExtra::arrangeGrob(
        grid::textGrob(
            "A",
            x = 0.01,
            hjust = 0,
            gp = grid::gpar(
                fontface = "bold",
                fontsize = 14
            )
        ),
        p_F$gtable,
        ncol = 1,
        heights = c(
            0.04,
            0.96
        )
    ),
    gridExtra::arrangeGrob(
        grid::textGrob(
            "B",
            x = 0.01,
            hjust = 0,
            gp = grid::gpar(
                fontface = "bold",
                fontsize = 14
            )
        ),
        p_M$gtable,
        ncol = 1,
        heights = c(
            0.04,
            0.96
        )
    ),
    ncol = 1
)


dev.off()


# ============================================================
# 14. Save RAW VST matrices used for the heatmap
#
# These contain ALL genes entered into the LRT.
# ============================================================

write.table(
    mat_F,
    file = file.path(
        heatmap_out_dir,
        "RNA_F_ALL_LRT_tested_genes_VST_matrix.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


write.table(
    mat_M,
    file = file.path(
        heatmap_out_dir,
        "RNA_M_ALL_LRT_tested_genes_VST_matrix.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


cat("\n")
cat("============================================================\n")
cat("ALL-LRT-gene expression heatmap finished\n")
cat("============================================================\n")

cat(
    "Female LRT-tested genes:",
    length(lrt_genes_F),
    "\n"
)

cat(
    "Male LRT-tested genes:",
    length(lrt_genes_M),
    "\n"
)

cat(
    "Combined PDF:",
    combined_pdf,
    "\n"
)
