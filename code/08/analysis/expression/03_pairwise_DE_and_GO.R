source(file.path(Sys.getenv("CODE_DIR"), "analysis", "common.R"))
list2env(readRDS(input_rds), envir = environment())
suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(UpSetR)
    library(clusterProfiler)
    library(GOSemSim)
    library(GO.db)
    library(org.Dm.eg.db)
    library(AnnotationDbi)
    library(gridExtra)
})

# ============================================================
# ============================================================

required_objects <- c(
    "rna_f",
    "sample_info",
    "count_out_dir"
)

missing_objects <- required_objects[
    !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
    stop(
        "Run this section after the current RNA LRT/Mfuzz analysis. ",
        "Missing objects: ",
        paste(missing_objects, collapse = ", ")
    )
}

# ============================================================
# 2. Analysis settings
# ============================================================

pairwise_out_dir <- file.path(
    count_out_dir,
    "RNA_pairwise_DE_current"
)

dir.create(
    pairwise_out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

min_sex_row_sum <- 100

# Keep pairwise calls consistent with the preceding LRT analysis.
pairwise_padj_cutoff <- 0.01

pairwise_fitType <- "local"

# Direction is determined only by sign:
#   log2FC > 0 -> Up
#   log2FC < 0 -> Down
lfc_cutoff <- 0

# Display-only cap for the violin plot. The paired Wilcoxon tests below are
# calculated before this cap is applied and therefore use the raw values.
violin_cap <- 0.75

# GO settings
go_ont <- "BP"
go_pvalue_cutoff <- 0.05
go_qvalue_cutoff <- 0.2
go_top_n <- 10
go_similarity_cutoff <- 0.60
# For display only: terms deeper than this are first mapped to a broad
# biological-process ancestor before semantic clustering.
go_parent_max_depth <- 4L
go_overlap_fraction_cutoff <- 0.50
go_overlap_jaccard_cutoff <- 0.25

if (!exists("analysis_seed", inherits = TRUE)) {
    analysis_seed <- as.integer(format(Sys.Date(), "%Y%m%d"))
}

set.seed(analysis_seed)

# ============================================================
# 3. Helpers
# ============================================================

format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 2.2e-16) return("< 2.2e-16")
    if (p < 0.001) return(formatC(p, format = "e", digits = 2))
    formatC(p, format = "f", digits = 3)
}

safe_wilcox_paired <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]

    if (length(x) < 2) {
        return(list(p.value = NA_real_, n = length(x)))
    }

    wt <- suppressWarnings(
        wilcox.test(
            x,
            y,
            paired = TRUE,
            exact = FALSE
        )
    )

    list(
        p.value = wt$p.value,
        n = length(x)
    )
}

# ============================================================
# 4. Core function:
#    run Wald pairwise DE for one sex
# ============================================================

run_pairwise_DE_by_sex <- function(
  count_mat,
  sample_info,
  sex_use,
  layer_name = "RNA",
  genotype_ref = "WBg",
  padj_cutoff = 0.05,
  min_row_sum = 100,
  fitType = "local",
  out_dir = "."
) {
    message("")
    message("============================================================")
    message("PAIRWISE DE: ", layer_name, " ", sex_use)
    message("============================================================")

    # ----------------------------------------------------------
    # ----------------------------------------------------------

    sample_info_sex <- sample_info[
        sample_info$sex == sex_use, ,
        drop = FALSE
    ]

    sample_info_sex <- sample_info_sex[
        colnames(count_mat)[
            colnames(count_mat) %in% sample_info_sex$sample
        ], ,
        drop = FALSE
    ]

    mat_sex <- count_mat[
        ,
        sample_info_sex$sample,
        drop = FALSE
    ]

    mat_sex <- round(as.matrix(mat_sex))
    mode(mat_sex) <- "integer"

    stopifnot(
        all(colnames(mat_sex) == rownames(sample_info_sex))
    )

    # ----------------------------------------------------------
    #     rowSums(mat_sex) >= 100
    # ----------------------------------------------------------

    keep <- rowSums(mat_sex) >= min_row_sum
    mat_sex <- mat_sex[
        keep, ,
        drop = FALSE
    ]

    message(
        layer_name, " ", sex_use,
        ": genes after count filter = ",
        nrow(mat_sex)
    )

    # ----------------------------------------------------------
    # 4.3 DESeq2 Wald model
    #
    # Do NOT reuse the LRT p-values for pairwise testing.
    # A fresh Wald fit is used for true pairwise tests.
    # ----------------------------------------------------------

    coldata <- sample_info_sex

    expected_genotypes <- c(
        "WBg",
        "AdarG",
        "AdarS"
    )

    missing_genotypes <- setdiff(
        expected_genotypes,
        unique(as.character(coldata$genotype))
    )

    if (length(missing_genotypes) > 0) {
        stop(
            sex_use,
            ": missing genotype(s): ",
            paste(missing_genotypes, collapse = ", ")
        )
    }

    coldata$genotype <- factor(
        coldata$genotype,
        levels = expected_genotypes
    )

    coldata$genotype <- relevel(
        coldata$genotype,
        ref = genotype_ref
    )

    dds <- DESeqDataSetFromMatrix(
        countData = mat_sex,
        colData = coldata,
        design = ~genotype
    )

    dds <- DESeq(
        dds,
        test = "Wald",
        fitType = fitType
    )

    normalized_counts <- counts(
        dds,
        normalized = TRUE
    )

    normalized_file <- file.path(
        out_dir,
        paste0(
            layer_name,
            "_",
            sex_use,
            "_pairwise_DESeq2_normalized_counts.tsv"
        )
    )

    normalized_dt <- as.data.table(
        normalized_counts,
        keep.rownames = "gene_id"
    )

    fwrite(
        normalized_dt,
        normalized_file,
        sep = "\t"
    )

    # ----------------------------------------------------------
    # 4.4 All three pairwise comparisons
    #
    # Main downstream figures will use:
    #   AdarG_vs_WBg
    #   AdarS_vs_WBg
    #
    # AdarG_vs_AdarS is also output for completeness.
    # ----------------------------------------------------------

    comparisons <- list(
        AdarG_vs_WBg = c("AdarG", "WBg"),
        AdarS_vs_WBg = c("AdarS", "WBg"),
        AdarG_vs_AdarS = c("AdarG", "AdarS")
    )

    result_list <- list()
    sig_list <- list()
    summary_list <- list()
    tested_genes_list <- list()

    for (comp_name in names(comparisons)) {
        numerator <- comparisons[[comp_name]][1]
        denominator <- comparisons[[comp_name]][2]

        res <- results(
            dds,
            contrast = c(
                "genotype",
                numerator,
                denominator
            ),
            alpha = padj_cutoff
        )

        res_dt <- as.data.table(
            as.data.frame(res),
            keep.rownames = "gene_id"
        )

        res_dt[
            ,
            `:=`(
                layer = layer_name,
                sex = sex_use,
                comparison = comp_name,
                numerator = numerator,
                denominator = denominator
            )
        ]

        # DEG status:
        # padj NA remains NA / NS, rather than being forced to 1.
        res_dt[
            ,
            direction := fifelse(
                !is.na(padj) &
                    is.finite(log2FoldChange) &
                    padj < padj_cutoff &
                    log2FoldChange > lfc_cutoff,
                "Up",
                fifelse(
                    !is.na(padj) &
                        is.finite(log2FoldChange) &
                        padj < padj_cutoff &
                        log2FoldChange < -lfc_cutoff,
                    "Down",
                    "NS"
                )
            )
        ]

        # Useful ordering
        setcolorder(
            res_dt,
            c(
                "gene_id",
                "layer",
                "sex",
                "comparison",
                "numerator",
                "denominator",
                "baseMean",
                "log2FoldChange",
                "lfcSE",
                "stat",
                "pvalue",
                "padj",
                "direction"
            )
        )

        res_dt[
            ,
            padj_na_tmp := is.na(padj)
        ]

        setorder(
            res_dt,
            padj_na_tmp,
            padj
        )

        res_dt[
            ,
            padj_na_tmp := NULL
        ]

        sig_dt <- res_dt[
            direction %in% c("Up", "Down")
        ]

        all_file <- file.path(
            out_dir,
            paste0(
                layer_name,
                "_",
                sex_use,
                "_",
                comp_name,
                "_DESeq2_Wald_all.tsv"
            )
        )

        sig_file <- file.path(
            out_dir,
            paste0(
                layer_name,
                "_",
                sex_use,
                "_",
                comp_name,
                "_DESeq2_Wald_sig_padj",
                padj_cutoff,
                ".tsv"
            )
        )

        fwrite(
            res_dt,
            all_file,
            sep = "\t"
        )

        fwrite(
            sig_dt,
            sig_file,
            sep = "\t"
        )

        tested_n <- res_dt[
            !is.na(padj),
            .N
        ]

        tested_genes_list[[comp_name]] <- res_dt[
            !is.na(padj),
            gene_id
        ]

        up_n <- sig_dt[
            direction == "Up",
            .N
        ]

        down_n <- sig_dt[
            direction == "Down",
            .N
        ]

        summary_list[[comp_name]] <- data.table(
            layer = layer_name,
            sex = sex_use,
            comparison = comp_name,
            numerator = numerator,
            denominator = denominator,
            genes_after_count_filter = nrow(mat_sex),
            genes_with_valid_padj = tested_n,
            DEG_total = up_n + down_n,
            DEG_up = up_n,
            DEG_down = down_n,
            padj_cutoff = padj_cutoff,
            lfc_cutoff = lfc_cutoff
        )

        result_list[[comp_name]] <- res_dt
        sig_list[[comp_name]] <- sig_dt

        message(
            sex_use, " ", comp_name,
            ": up = ", up_n,
            ", down = ", down_n,
            ", total = ", up_n + down_n
        )
    }

    summary_dt <- rbindlist(
        summary_list,
        use.names = TRUE,
        fill = TRUE
    )

    fwrite(
        summary_dt,
        file.path(
            out_dir,
            paste0(
                layer_name,
                "_",
                sex_use,
                "_pairwise_DE_summary.tsv"
            )
        ),
        sep = "\t"
    )

    return(
        list(
            layer = layer_name,
            sex = sex_use,
            dds = dds,
            normalized_counts = normalized_counts,
            tested_genes = unique(unlist(tested_genes_list)),
            tested_genes_by_comparison = tested_genes_list,
            results = result_list,
            sig = sig_list,
            summary = summary_dt
        )
    )
}

# ============================================================
# 5. Run Female and Male
# ============================================================

pairwise_F <- run_pairwise_DE_by_sex(
    count_mat = rna_f,
    sample_info = sample_info,
    sex_use = "F",
    layer_name = "RNA",
    genotype_ref = "WBg",
    padj_cutoff = pairwise_padj_cutoff,
    min_row_sum = min_sex_row_sum,
    fitType = pairwise_fitType,
    out_dir = pairwise_out_dir
)

pairwise_M <- run_pairwise_DE_by_sex(
    count_mat = rna_f,
    sample_info = sample_info,
    sex_use = "M",
    layer_name = "RNA",
    genotype_ref = "WBg",
    padj_cutoff = pairwise_padj_cutoff,
    min_row_sum = min_sex_row_sum,
    fitType = pairwise_fitType,
    out_dir = pairwise_out_dir
)

# ============================================================
# 6. Combined summary
# ============================================================

pairwise_summary_all <- rbindlist(
    list(
        pairwise_F$summary,
        pairwise_M$summary
    ),
    use.names = TRUE,
    fill = TRUE
)

fwrite(
    pairwise_summary_all,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_DE_summary.ALL.tsv"
    ),
    sep = "\t"
)

message("")
message("Pairwise DE summary:")
print(pairwise_summary_all)

# ============================================================
# 7. Main DEG sets:
#    mutant vs WBg only
#
#   AdarG_F_up
#   AdarS_F_up
#   AdarG_F_down
#   AdarS_F_down
#   AdarG_M_up
#   AdarS_M_up
#   AdarG_M_down
#   AdarS_M_down
# ============================================================

# data.table NSE is clearer here without helper shadowing
AdarG_F_up <- unique(
    pairwise_F$sig[["AdarG_vs_WBg"]][
        direction == "Up",
        gene_id
    ]
)

AdarG_F_down <- unique(
    pairwise_F$sig[["AdarG_vs_WBg"]][
        direction == "Down",
        gene_id
    ]
)

AdarS_F_up <- unique(
    pairwise_F$sig[["AdarS_vs_WBg"]][
        direction == "Up",
        gene_id
    ]
)

AdarS_F_down <- unique(
    pairwise_F$sig[["AdarS_vs_WBg"]][
        direction == "Down",
        gene_id
    ]
)

AdarG_M_up <- unique(
    pairwise_M$sig[["AdarG_vs_WBg"]][
        direction == "Up",
        gene_id
    ]
)

AdarG_M_down <- unique(
    pairwise_M$sig[["AdarG_vs_WBg"]][
        direction == "Down",
        gene_id
    ]
)

AdarS_M_up <- unique(
    pairwise_M$sig[["AdarS_vs_WBg"]][
        direction == "Up",
        gene_id
    ]
)

AdarS_M_down <- unique(
    pairwise_M$sig[["AdarS_vs_WBg"]][
        direction == "Down",
        gene_id
    ]
)

sets_direction <- list(
    AdarG_F_up = AdarG_F_up,
    AdarS_F_up = AdarS_F_up,
    AdarG_F_down = AdarG_F_down,
    AdarS_F_down = AdarS_F_down,
    AdarG_M_up = AdarG_M_up,
    AdarS_M_up = AdarS_M_up,
    AdarG_M_down = AdarG_M_down,
    AdarS_M_down = AdarS_M_down
)

# Use the original UpSetR layout so that the intersection bars, set-size bars,
# and membership matrix share the package's native coordinate system.
# `newpage = FALSE` is important for direct PDF output: UpSetR already creates
# the page layout internally, so allowing print.upset() to open another page
# can leave a blank first page in non-interactive runs.
save_upset_classic_pdf <- function(
  upset_data,
  nsets,
  pdf_file,
  width,
  height,
  text_scale,
  max_intersections = 40L
) {
    upset_data <- as.data.frame(upset_data)
    set_names <- colnames(upset_data)

    if (length(set_names) < 2L) {
        stop("At least two non-empty sets are required for an UpSet plot.")
    }

    if (length(set_names) != nsets) {
        warning(
            "nsets does not match the number of columns; using all set columns."
        )
    }

    upset_data[set_names] <- lapply(
        upset_data[set_names],
        function(x) as.integer(as.numeric(x) > 0)
    )

    upset_data <- upset_data[
        rowSums(as.matrix(upset_data[set_names])) > 0, ,
        drop = FALSE
    ]

    if (nrow(upset_data) == 0L) {
        stop("The UpSet input contains no rows with set membership.")
    }

    grDevices::pdf(
        file = pdf_file,
        width = width,
        height = height,
        onefile = TRUE,
        useDingbats = FALSE,
        compress = TRUE
    )
    on.exit(grDevices::dev.off(), add = TRUE)

    upset_plot <- UpSetR::upset(
        upset_data,
        sets = set_names,
        nsets = nsets,
        nintersects = max_intersections,
        order.by = "freq",
        keep.order = FALSE,
        mb.ratio = c(0.6, 0.4),
        text.scale = text_scale,
        number.angles = 0,
        point.size = 5,
        line.size = 2,
        mainbar.y.label = "Genes Intersections",
        sets.x.label = "Genes Per Library"
    )

    print(upset_plot, newpage = FALSE)
    invisible(upset_plot)
}

# ============================================================
# 8. UpSet 1:
#    8 direction-specific sets
# ============================================================

nonempty_direction <- sets_direction[
    lengths(sets_direction) > 0
]

if (length(nonempty_direction) >= 2) {
    upset_data_direction <- UpSetR::fromList(
        nonempty_direction
    )

    upset_direction_pdf <- file.path(
        pairwise_out_dir,
        "RNA_pairwise_upset_8_direction_sets.pdf"
    )

    save_upset_classic_pdf(
        upset_data = upset_data_direction,
        nsets = length(nonempty_direction),
        pdf_file = upset_direction_pdf,
        width = 12,
        height = 8,
        text_scale = 1.3
    )

    message(
        "Saved vector PDF: ", upset_direction_pdf
    )
}

# ============================================================
# 9. UpSet 2:
#    merge Up + Down for each mutant-vs-WBg comparison
# ============================================================

sets_any_DE <- list(
    AdarG_F = unique(
        c(
            AdarG_F_up,
            AdarG_F_down
        )
    ),
    AdarS_F = unique(
        c(
            AdarS_F_up,
            AdarS_F_down
        )
    ),
    AdarG_M = unique(
        c(
            AdarG_M_up,
            AdarG_M_down
        )
    ),
    AdarS_M = unique(
        c(
            AdarS_M_up,
            AdarS_M_down
        )
    )
)

nonempty_any <- sets_any_DE[
    lengths(sets_any_DE) > 0
]

if (length(nonempty_any) >= 2) {
    upset_data_any <- UpSetR::fromList(
        nonempty_any
    )

    upset_any_pdf <- file.path(
        pairwise_out_dir,
        "RNA_pairwise_upset_any_DE.pdf"
    )

    save_upset_classic_pdf(
        upset_data = upset_data_any,
        nsets = length(nonempty_any),
        pdf_file = upset_any_pdf,
        width = 9,
        height = 7,
        text_scale = 1.4
    )

    message(
        "Saved vector PDF: ", upset_any_pdf
    )
}

# ============================================================
# 10. |log2FoldChange| violin plot using EXACTLY the same
#     8 significant DEG sets used in the direction-specific UpSet
#
# IMPORTANT:
#   * The plotted genes are exactly:
#       AdarG_F_up, AdarS_F_up, AdarG_F_down, AdarS_F_down,
#       AdarG_M_up, AdarS_M_up, AdarG_M_down, AdarS_M_down
#   * These are the same objects stored in sets_direction and used by UpSet 1.
#   * y = |log2FC vs WBg|, so Up and Down are both shown as perturbation magnitude.
#   * Because AdarG and AdarS DEG sets are not necessarily the same genes,
#     The optional comparison below therefore uses an unpaired Wilcoxon rank-sum test.
# ============================================================

# ------------------------------------------------------------
# 10.1 Metadata for the exact 8 UpSet sets
# ------------------------------------------------------------

upset8_meta <- data.table(
    set_name = c(
        "AdarG_F_up",
        "AdarS_F_up",
        "AdarG_F_down",
        "AdarS_F_down",
        "AdarG_M_up",
        "AdarS_M_up",
        "AdarG_M_down",
        "AdarS_M_down"
    ),
    genotype = c(
        "AdarG", "AdarS", "AdarG", "AdarS",
        "AdarG", "AdarS", "AdarG", "AdarS"
    ),
    sex_code = c(
        "F", "F", "F", "F",
        "M", "M", "M", "M"
    ),
    sex = c(
        "Female", "Female", "Female", "Female",
        "Male", "Male", "Male", "Male"
    ),
    direction = c(
        "Up", "Up", "Down", "Down",
        "Up", "Up", "Down", "Down"
    ),
    comparison = c(
        "AdarG_vs_WBg", "AdarS_vs_WBg",
        "AdarG_vs_WBg", "AdarS_vs_WBg",
        "AdarG_vs_WBg", "AdarS_vs_WBg",
        "AdarG_vs_WBg", "AdarS_vs_WBg"
    )
)

# Sanity check: these names must be exactly the 8 sets used by UpSet.
stopifnot(
    setequal(
        upset8_meta$set_name,
        names(sets_direction)
    )
)


# ------------------------------------------------------------
# 10.2 Extract |log2FC| for one exact UpSet DEG set
# ------------------------------------------------------------

extract_upset8_abs_lfc <- function(meta_row) {
    set_name_i <- as.character(meta_row$set_name)
    sex_code_i <- as.character(meta_row$sex_code)
    comparison_i <- as.character(meta_row$comparison)

    gene_set_i <- unique(
        sets_direction[[set_name_i]]
    )

    res_obj_i <- if (sex_code_i == "F") {
        pairwise_F
    } else {
        pairwise_M
    }

    res_i <- copy(
        res_obj_i$results[[comparison_i]]
    )

    out_i <- res_i[
        gene_id %in% gene_set_i &
            is.finite(log2FoldChange),
        .(
            gene_id,
            log2FC = log2FoldChange,
            abs_log2FC = abs(log2FoldChange),
            padj
        )
    ]

    out_i[
        ,
        `:=`(
            set_name = set_name_i,
            genotype = as.character(meta_row$genotype),
            sex_code = sex_code_i,
            sex = as.character(meta_row$sex),
            direction = as.character(meta_row$direction),
            comparison = comparison_i
        )
    ]

    # The plotting set must match the UpSet set exactly.
    # All genes in sets_direction came from the corresponding significant table,
    # so a mismatch indicates an unexpected upstream inconsistency.
    if (uniqueN(out_i$gene_id) != length(gene_set_i)) {
        missing_i <- setdiff(gene_set_i, out_i$gene_id)
        stop(
            "UpSet/violin gene-set mismatch for ", set_name_i,
            ": UpSet n = ", length(gene_set_i),
            ", violin n = ", uniqueN(out_i$gene_id),
            ". Missing examples: ",
            paste(head(missing_i, 10), collapse = ", ")
        )
    }

    out_i
}


# ------------------------------------------------------------
# 10.3 Build plotting data from the SAME 8 UpSet sets
# ------------------------------------------------------------

upset8_violin_list <- lapply(
    seq_len(nrow(upset8_meta)),
    function(i) {
        extract_upset8_abs_lfc(
            upset8_meta[i]
        )
    }
)

lfc_plot <- rbindlist(
    upset8_violin_list,
    use.names = TRUE,
    fill = TRUE
)

lfc_plot[
    ,
    genotype := factor(
        genotype,
        levels = c(
            "AdarS",
            "AdarG"
        )
    )
]

lfc_plot[
    ,
    sex := factor(
        sex,
        levels = c(
            "Female",
            "Male"
        )
    )
]

lfc_plot[
    ,
    direction := factor(
        direction,
        levels = c(
            "Up",
            "Down"
        )
    )
]


# ------------------------------------------------------------
# 10.4 QC: violin n must equal the UpSet set size
# ------------------------------------------------------------

upset8_expected_n <- data.table(
    set_name = names(sets_direction),
    upset_n = lengths(sets_direction)
)

upset8_observed_n <- lfc_plot[
    ,
    .(
        violin_n = uniqueN(gene_id)
    ),
    by = set_name
]

upset8_qc <- merge(
    upset8_expected_n,
    upset8_observed_n,
    by = "set_name",
    all = TRUE
)

upset8_qc[
    ,
    identical_n := upset_n == violin_n
]

print(upset8_qc)

if (any(!upset8_qc$identical_n)) {
    stop(
        "The 8 violin gene sets are not identical to the 8 UpSet gene sets."
    )
}

fwrite(
    upset8_qc,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_UPSET8_violin_gene_set_QC.tsv"
    ),
    sep = "\t"
)

fwrite(
    lfc_plot,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_UPSET8_absolute_log2FC_violin_data.tsv"
    ),
    sep = "\t"
)


# ------------------------------------------------------------
# 10.5 Summary of the exact 8 sets
# Use RAW |log2FC| for summary; display cap does not alter statistics.
# ------------------------------------------------------------

violin_summary <- lfc_plot[
    is.finite(abs_log2FC),
    .(
        n = uniqueN(gene_id),
        mean_abs_log2FC = mean(abs_log2FC),
        median_abs_log2FC = median(abs_log2FC),
        sd_abs_log2FC = sd(abs_log2FC)
    ),
    by = .(
        set_name,
        sex,
        direction,
        genotype
    )
]

violin_summary[
    ,
    label := sprintf(
        "italic(n) == %d ~ ',' ~ %.3f %%+-%% %.3f",
        n,
        mean_abs_log2FC,
        sd_abs_log2FC
    )
]

fwrite(
    violin_summary,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_UPSET8_absolute_log2FC_violin_summary.tsv"
    ),
    sep = "\t"
)

print(violin_summary)


# ------------------------------------------------------------
# 10.6 AdarS vs AdarG distribution comparison within each
#
# NOTE:
# The two DEG sets are not guaranteed to contain the same genes,
# ------------------------------------------------------------

safe_wilcox_unpaired <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    if (length(x) < 2 || length(y) < 2) {
        return(NA_real_)
    }

    suppressWarnings(
        wilcox.test(
            x,
            y,
            paired = FALSE,
            exact = FALSE
        )$p.value
    )
}

wilcox_summary <- lfc_plot[
    ,
    {
        genes_S <- unique(gene_id[genotype == "AdarS"])
        genes_G <- unique(gene_id[genotype == "AdarG"])

        .(
            n_AdarS = length(genes_S),
            n_AdarG = length(genes_G),
            n_overlap_genes = length(intersect(genes_S, genes_G)),
            p_value = safe_wilcox_unpaired(
                abs_log2FC[genotype == "AdarS"],
                abs_log2FC[genotype == "AdarG"]
            )
        )
    },
    by = .(
        sex,
        direction
    )
]

wilcox_summary[
    ,
    padj_BH := p.adjust(
        p_value,
        method = "BH"
    )
]

fwrite(
    wilcox_summary,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_UPSET8_absolute_log2FC_AdarS_vs_AdarG_wilcox_unpaired.tsv"
    ),
    sep = "\t"
)

print(wilcox_summary)


# ------------------------------------------------------------
# 10.7 Display cap ONLY for plotting
# ------------------------------------------------------------

if (is.null(violin_cap)) {
    lfc_plot[
        ,
        value_plot := abs_log2FC
    ]

    y_max <- max(
        lfc_plot$value_plot,
        na.rm = TRUE
    ) * 1.08

    y_lab <- expression("|log"[2] * "FC vs WT|")
} else {
    lfc_plot[
        ,
        value_plot := pmin(
            abs_log2FC,
            violin_cap
        )
    ]

    y_max <- violin_cap

    y_lab <- bquote(
        "|log"[2] * "FC vs WT|" ~ "(display capped at" ~ .(violin_cap) * ")"
    )
}


# ------------------------------------------------------------
# 10.8 P-value and n labels
# ------------------------------------------------------------

plot_summary <- copy(violin_summary)
plot_summary[
    ,
    label_y := y_max * 0.97
]

p_annotation_df <- copy(wilcox_summary)
p_annotation_df[
    ,
    `:=`(
        x = 1.5,
        y = y_max * 0.80,
        label = paste0(
            "italic(P) == ",
            vapply(
                p_value,
                format_p,
                character(1)
            )
        )
    )
]


# ------------------------------------------------------------
# ------------------------------------------------------------

p_lfc <- ggplot(
    lfc_plot,
    aes(
        x = genotype,
        y = value_plot,
        fill = genotype
    )
) +
    geom_violin(
        trim = TRUE,
        color = NA,
        alpha = 0.75
    ) +
    geom_boxplot(
        width = 0.18,
        outlier.shape = NA,
        alpha = 0.8
    ) +
    facet_grid(
        direction ~ sex
    ) +
    geom_text(
        data = plot_summary,
        aes(
            x = genotype,
            y = label_y,
            label = label
        ),
        inherit.aes = FALSE,
        parse = TRUE,
        vjust = 1,
        size = 3.3
    ) +
    geom_text(
        data = p_annotation_df,
        aes(
            x = x,
            y = y,
            label = label
        ),
        inherit.aes = FALSE,
        parse = TRUE,
        vjust = 0,
        size = 4
    ) +
    scale_y_continuous(
        limits = c(0, y_max),
        expand = expansion(
            mult = c(0, 0.02)
        )
    ) +
    theme_classic(
        base_size = 14
    ) +
    labs(
        x = NULL,
        y = y_lab
    ) +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(
            size = 14
        )
    )


# ------------------------------------------------------------
# 10.10 Save
# ------------------------------------------------------------

lfc_pdf <- file.path(
    pairwise_out_dir,
    "RNA_pairwise_UPSET8_absolute_log2FC_violin.pdf"
)

lfc_png <- file.path(
    pairwise_out_dir,
    "RNA_pairwise_UPSET8_absolute_log2FC_violin.png"
)

ggsave(
    lfc_pdf,
    p_lfc,
    width = 9,
    height = 7
)

ggsave(
    lfc_png,
    p_lfc,
    width = 9,
    height = 7,
    dpi = 300
)

message(
    "Saved UpSet-consistent |log2FC| violin: ",
    lfc_pdf
)

# ============================================================
# 11. GO enrichment:
#     8 mutant-vs-WBg direction-specific sets
#
# IMPORTANT:
#   Use sex-specific tested genes as the background universe,
# ============================================================

go_group_info <- list(
    AdarG_F_up = list(
        genes = AdarG_F_up,
        sex = "Female",
        sex_code = "F",
        genotype = "AdarG",
        direction = "Up",
        comparison = "AdarG_vs_WBg"
    ),
    AdarS_F_up = list(
        genes = AdarS_F_up,
        sex = "Female",
        sex_code = "F",
        genotype = "AdarS",
        direction = "Up",
        comparison = "AdarS_vs_WBg"
    ),
    AdarG_F_down = list(
        genes = AdarG_F_down,
        sex = "Female",
        sex_code = "F",
        genotype = "AdarG",
        direction = "Down",
        comparison = "AdarG_vs_WBg"
    ),
    AdarS_F_down = list(
        genes = AdarS_F_down,
        sex = "Female",
        sex_code = "F",
        genotype = "AdarS",
        direction = "Down",
        comparison = "AdarS_vs_WBg"
    ),
    AdarG_M_up = list(
        genes = AdarG_M_up,
        sex = "Male",
        sex_code = "M",
        genotype = "AdarG",
        direction = "Up",
        comparison = "AdarG_vs_WBg"
    ),
    AdarS_M_up = list(
        genes = AdarS_M_up,
        sex = "Male",
        sex_code = "M",
        genotype = "AdarS",
        direction = "Up",
        comparison = "AdarS_vs_WBg"
    ),
    AdarG_M_down = list(
        genes = AdarG_M_down,
        sex = "Male",
        sex_code = "M",
        genotype = "AdarG",
        direction = "Down",
        comparison = "AdarG_vs_WBg"
    ),
    AdarS_M_down = list(
        genes = AdarS_M_down,
        sex = "Male",
        sex_code = "M",
        genotype = "AdarS",
        direction = "Down",
        comparison = "AdarS_vs_WBg"
    )
)

background_by_comparison <- list(
    F_AdarG_vs_WBg = unique(
        pairwise_F$tested_genes_by_comparison[["AdarG_vs_WBg"]]
    ),
    F_AdarS_vs_WBg = unique(
        pairwise_F$tested_genes_by_comparison[["AdarS_vs_WBg"]]
    ),
    M_AdarG_vs_WBg = unique(
        pairwise_M$tested_genes_by_comparison[["AdarG_vs_WBg"]]
    ),
    M_AdarS_vs_WBg = unique(
        pairwise_M$tested_genes_by_comparison[["AdarS_vs_WBg"]]
    )
)

ego_list <- list()

for (group_name in names(go_group_info)) {
    info <- go_group_info[[group_name]]

    genes_use <- unique(
        info$genes
    )

    background_key <- paste(
        info$sex_code,
        info$comparison,
        sep = "_"
    )

    background_use <- unique(
        background_by_comparison[[background_key]]
    )

    genes_use <- intersect(
        genes_use,
        background_use
    )

    message(
        "GO ", group_name,
        ": genes = ", length(genes_use),
        ", background = ", length(background_use)
    )

    if (length(genes_use) < 3) {
        message(
            "Skip ", group_name,
            ": fewer than 3 genes."
        )
        next
    }

    ego <- enrichGO(
        gene = genes_use,
        universe = background_use,
        OrgDb = org.Dm.eg.db,
        keyType = "FLYBASE",
        ont = go_ont,
        pAdjustMethod = "BH",
        minGSSize = 3,
        pvalueCutoff = go_pvalue_cutoff,
        qvalueCutoff = go_qvalue_cutoff,
        readable = TRUE
    )

    ego_df <- as.data.frame(
        ego
    )

    if (nrow(ego_df) == 0) {
        message(
            "No significant GO terms: ",
            group_name
        )
        next
    }

    ego <- clusterProfiler::simplify(
        ego,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min,
        measure = "Wang"
    )

    ego_df <- as.data.frame(
        ego
    )

    if (nrow(ego_df) == 0) {
        next
    }

    ego_df$group <- group_name
    ego_df$direction <- info$direction
    ego_df$sex <- info$sex
    ego_df$genotype <- info$genotype
    ego_df$comparison <- info$comparison
    ego_df$background_key <- background_key
    ego_df$input_gene_n <- length(genes_use)
    ego_df$background_gene_n <- length(background_use)

    ego_list[[group_name]] <- ego_df
}

if (length(ego_list) > 0) {
    go_all <- bind_rows(
        ego_list
    )

    go_all <- go_all %>%
        mutate(
            sex = factor(
                sex,
                levels = c(
                    "Female",
                    "Male"
                )
            ),
            genotype = factor(
                genotype,
                levels = c(
                    "AdarG",
                    "AdarS"
                )
            ),
            direction = factor(
                direction,
                levels = c(
                    "Up",
                    "Down"
                )
            ),
            minus_log10_FDR = -log10(
                p.adjust
            )
        )

    fwrite(
        as.data.table(go_all),
        file.path(
            pairwise_out_dir,
            "RNA_pairwise_GO_all.tsv"
        ),
        sep = "\t"
    )

    # ----------------------------------------------------------
    # Select terms that are shared across the most groups first,
    # then ranked by the best adjusted P value.
    # ----------------------------------------------------------

    select_shared_GO <- function(go_data,
                                 direction_use,
                                 top_n = 10) {
        term_rank <- go_data %>%
            filter(
                direction == direction_use
            ) %>%
            group_by(
                ID,
                Description
            ) %>%
            summarise(
                group_n = n_distinct(
                    group
                ),
                best_FDR = min(
                    p.adjust,
                    na.rm = TRUE
                ),
                .groups = "drop"
            ) %>%
            arrange(
                desc(group_n),
                best_FDR
            ) %>%
            slice_head(
                n = top_n
            )

        plot_data <- go_data %>%
            filter(
                direction == direction_use
            ) %>%
            semi_join(
                term_rank,
                by = c(
                    "ID",
                    "Description"
                )
            )

        plot_data$ID <- factor(
            plot_data$ID,
            levels = rev(
                term_rank$ID
            )
        )

        list(
            ranking = term_rank,
            plot_data = plot_data
        )
    }

    # ----------------------------------------------------------
    # Up-regulated GO
    # ----------------------------------------------------------

    go_up <- select_shared_GO(
        go_data = go_all,
        direction_use = "Up",
        top_n = go_top_n
    )

    if (nrow(go_up$plot_data) > 0) {
        up_label <- setNames(
            go_up$ranking$Description,
            go_up$ranking$ID
        )

        p_up <- ggplot(
            go_up$plot_data,
            aes(
                x = genotype,
                y = ID
            )
        ) +
            geom_point(
                aes(
                    size = minus_log10_FDR,
                    color = Count
                )
            ) +
            facet_grid(
                . ~ sex,
                scales = "free_x",
                space = "free_x"
            ) +
            scale_y_discrete(
                labels = up_label
            ) +
            scale_size_continuous(
                name = "Gene number",
                range = c(2, 10)
            ) +
            scale_color_gradientn(
                name = expression(-log[10](FDR)),
                colours = c(
                    "#3F7FB5",
                    "#8EB5D3",
                    "#F7F7F7",
                    "#D7898F",
                    "#C43D4B"
                )
            ) +
            theme_classic(
                base_size = 14
            ) +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(
                    size = 15
                ),
                axis.text.x = element_text(
                    size = 14,
                    face = "italic"
                ),
                axis.text.y = element_text(
                    size = 11
                ),
                axis.title = element_blank(),
                plot.title = element_text(
                    hjust = 0.5,
                    size = 17
                )
            ) +
            ggtitle(
                "BP GO: Up-regulated genes"
            )

        ggsave(
            file.path(
                pairwise_out_dir,
                "GO_up_source_terms_top10.pdf"
            ),
            p_up,
            width = 9,
            height = 6
        )
    }

    # ----------------------------------------------------------
    # Down-regulated GO
    # ----------------------------------------------------------

    go_down <- select_shared_GO(
        go_data = go_all,
        direction_use = "Down",
        top_n = go_top_n
    )

    if (nrow(go_down$plot_data) > 0) {
        down_label <- setNames(
            go_down$ranking$Description,
            go_down$ranking$ID
        )

        p_down <- ggplot(
            go_down$plot_data,
            aes(
                x = genotype,
                y = ID
            )
        ) +
            geom_point(
                aes(
                    size = minus_log10_FDR,
                    color = Count
                )
            ) +
            facet_grid(
                . ~ sex,
                scales = "free_x",
                space = "free_x"
            ) +
            scale_y_discrete(
                labels = down_label
            ) +
            scale_size_continuous(
                name = "Gene number",
                range = c(2, 10)
            ) +
            scale_color_gradientn(
                name = expression(-log[10](FDR)),
                colours = c(
                    "#3F7FB5",
                    "#8EB5D3",
                    "#F7F7F7",
                    "#D7898F",
                    "#C43D4B"
                )
            ) +
            theme_classic(
                base_size = 14
            ) +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(
                    size = 15
                ),
                axis.text.x = element_text(
                    size = 14,
                    face = "italic"
                ),
                axis.text.y = element_text(
                    size = 11
                ),
                axis.title = element_blank(),
                plot.title = element_text(
                    hjust = 0.5,
                    size = 17
                )
            ) +
            ggtitle(
                "BP GO: Down-regulated genes"
            )

        ggsave(
            file.path(
                pairwise_out_dir,
                "GO_down_source_terms_top10.pdf"
            ),
            p_down,
            width = 11,
            height = 6
        )
    }
} else {
    warning(
        "No significant GO results were found for any pairwise DEG set."
    )
}

# ============================================================
# 12. Cross-group GO semantic collapse for the final GO bubble plots
#
# This section collapses similar terms across all eight DEG sets,
# prefers the broadest available representative term, and then
# keeps the original term-by-group GO bubble-plot format.
# ============================================================

if (exists("go_all", inherits = TRUE) && nrow(go_all) > 0) {
    get_go_isa_parents <- function(go_id) {
        parents <- mget(
            go_id,
            GO.db::GOBPPARENTS,
            ifnotfound = NA
        )[[1]]

        if (length(parents) == 0 || all(is.na(parents))) {
            return(character(0))
        }

        parents <- unname(
            parents[
                names(parents) %in% c("isa", "is_a")
            ]
        )

        parents <- parents[
            !is.na(parents) &
                parents != "all"
        ]

        unique(parents)
    }

    get_go_ancestors <- function(go_id) {
        ancestors <- mget(
            go_id,
            GO.db::GOBPANCESTOR,
            ifnotfound = NA
        )[[1]]

        if (length(ancestors) == 0 || all(is.na(ancestors))) {
            return(character(0))
        }

        ancestors <- unname(ancestors)
        unique(
            ancestors[
                !is.na(ancestors) &
                    ancestors != "all"
            ]
        )
    }

    go_depth_cache <- new.env(parent = emptyenv())

    get_go_depth <- function(go_id) {
        if (exists(go_id, envir = go_depth_cache, inherits = FALSE)) {
            return(get(go_id, envir = go_depth_cache, inherits = FALSE))
        }

        if (go_id == "GO:0008150") {
            assign(go_id, 0L, envir = go_depth_cache)
            return(0L)
        }

        parents <- get_go_isa_parents(go_id)

        if (length(parents) == 0) {
            assign(go_id, NA_integer_, envir = go_depth_cache)
            return(NA_integer_)
        }

        parent_depths <- vapply(
            parents,
            get_go_depth,
            integer(1)
        )

        parent_depths <- parent_depths[
            is.finite(parent_depths)
        ]

        depth <- if (length(parent_depths) == 0) {
            NA_integer_
        } else {
            as.integer(min(parent_depths) + 1L)
        }

        assign(go_id, depth, envir = go_depth_cache)
        depth
    }

    get_go_descendant_n <- function(go_id) {
        descendants <- mget(
            go_id,
            GO.db::GOBPOFFSPRING,
            ifnotfound = NA
        )[[1]]

        if (length(descendants) == 0 || all(is.na(descendants))) {
            return(0L)
        }

        length(unique(unname(descendants)))
    }

    get_go_term_description <- function(go_id) {
        term_obj <- mget(
            go_id,
            GO.db::GOTERM,
            ifnotfound = NA
        )[[1]]

        if (
            length(term_obj) == 0 ||
                (!isS4(term_obj) && all(is.na(term_obj)))
        ) {
            return(NA_character_)
        }

        as.character(AnnotationDbi::Term(term_obj))
    }

    choose_go_parent <- function(go_id,
                                 max_depth = go_parent_max_depth) {
        candidates <- unique(
            c(
                go_id,
                get_go_ancestors(go_id)
            )
        )

        candidate_depth <- vapply(
            candidates,
            get_go_depth,
            integer(1)
        )

        candidates <- candidates[
            is.finite(candidate_depth) &
                candidate_depth <= max_depth
        ]
        candidate_depth <- candidate_depth[
            is.finite(candidate_depth) &
                candidate_depth <= max_depth
        ]

        if (length(candidates) == 0) {
            return(
                data.frame(
                    parent_ID = go_id,
                    parent_depth = NA_integer_,
                    stringsAsFactors = FALSE
                )
            )
        }

        # Keep the deepest ancestor still considered broad. If multiple
        # branches have the same depth, prefer the narrower term so that,
        # for example, amino-acid terms collapse to "amino acid metabolic
        # process" rather than the generic "metabolic process".
        deepest <- max(candidate_depth)
        candidates <- candidates[candidate_depth == deepest]
        descendant_n <- vapply(
            candidates,
            get_go_descendant_n,
            integer(1)
        )
        parent_id <- sort(
            candidates[
                descendant_n == min(descendant_n)
            ]
        )[1]

        data.frame(
            parent_ID = parent_id,
            parent_depth = as.integer(deepest),
            stringsAsFactors = FALSE
        )
    }

    safe_min <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) NA_real_ else min(x)
    }

    collapse_gene_ids <- function(x) {
        x <- x[!is.na(x)]

        if (length(x) == 0) {
            return("")
        }

        gene_ids <- unique(
            unlist(
                strsplit(x, "/", fixed = TRUE),
                use.names = FALSE
            )
        )

        gene_ids <- gene_ids[nzchar(gene_ids)]

        paste(sort(gene_ids), collapse = "/")
    }

    split_gene_ids <- function(x) {
        x <- x[!is.na(x)]

        if (length(x) == 0) {
            return(character(0))
        }

        gene_ids <- unique(
            unlist(
                strsplit(x, "/", fixed = TRUE),
                use.names = FALSE
            )
        )

        gene_ids[nzchar(gene_ids)]
    }

    count_gene_ids <- function(x) {
        x <- collapse_gene_ids(x)
        if (!nzchar(x)) 0L else length(strsplit(x, "/", fixed = TRUE)[[1]])
    }

    collapse_go_terms <- function(go_data,
                                  similarity_cutoff = 0.70) {
        go_data <- as.data.frame(go_data)
        term_ids <- unique(as.character(go_data$ID))

        if (length(term_ids) == 0) {
            return(go_data)
        }

        parent_map <- do.call(
            rbind,
            lapply(
                term_ids,
                function(id) {
                    parent_row <- choose_go_parent(id)
                    data.frame(
                        ID = id,
                        parent_ID = parent_row$parent_ID,
                        parent_depth = parent_row$parent_depth,
                        stringsAsFactors = FALSE
                    )
                }
            )
        )
        rownames(parent_map) <- NULL
        parent_map$parent_Description <- vapply(
            parent_map$parent_ID,
            get_go_term_description,
            character(1)
        )

        # Prefer the description already returned by enrichGO when the
        # parent term is one of the input terms; otherwise use GO.db.
        raw_description <- go_data[, c("ID", "Description"), drop = FALSE]
        raw_description <- raw_description[
            !duplicated(raw_description$ID), ,
            drop = FALSE
        ]
        description_by_id <- setNames(
            as.character(raw_description$Description),
            as.character(raw_description$ID)
        )
        parent_description_from_input <- unname(
            description_by_id[parent_map$parent_ID]
        )
        use_input_description <- !is.na(parent_description_from_input) &
            nzchar(parent_description_from_input)
        parent_map$parent_Description[use_input_description] <-
            parent_description_from_input[use_input_description]

        # GO.db does not always encode modified/L/non-proteinogenic amino-acid
        # terms under one common parent. For display, consolidate those explicit
        # amino-acid metabolic/biosynthetic/catabolic labels into the standard
        # amino acid metabolic process term, while retaining all source IDs.
        source_description <- unname(
            description_by_id[parent_map$ID]
        )
        amino_acid_family <- !is.na(source_description) &
            grepl("amino acid", source_description, ignore.case = TRUE) &
            grepl(
                "metabolic|biosynthetic|catabolic",
                source_description,
                ignore.case = TRUE
            )
        parent_map$display_rule <- "GO_parent_depth"
        parent_map$display_rule[amino_acid_family] <-
            "amino_acid_family"
        parent_map$parent_ID[amino_acid_family] <- "GO:0006520"
        parent_map$parent_depth[amino_acid_family] <- get_go_depth(
            "GO:0006520"
        )
        parent_map$parent_Description[amino_acid_family] <-
            get_go_term_description("GO:0006520")

        # Some broad GO terms remain visually redundant with the amino-acid
        # family even after GO-parent mapping. For the final bubble plot only,
        # hide a broad term when at least half of its enriched genes overlap the
        # explicit amino-acid family and the two gene sets have sufficient
        # Jaccard overlap. The complete enrichment table is retained, and the
        # reason is recorded below, so this is a display filter rather than a
        # change to the enrichment statistics.
        go_data$plot_suppressed <- FALSE
        go_data$plot_suppress_reason <- ""
        broad_metabolism <- grepl(
            "organic acid metabolic|amide metabolic|small molecule (catabolic|biosynthetic) process",
            as.character(go_data$Description),
            ignore.case = TRUE
        )
        group_direction <- interaction(
            go_data$group,
            go_data$direction,
            drop = TRUE,
            sep = "\r"
        )

        for (group_direction_i in unique(group_direction)) {
            row_i <- which(group_direction == group_direction_i)
            amino_i <- row_i[
                grepl(
                    "amino acid",
                    as.character(go_data$Description[row_i]),
                    ignore.case = TRUE
                ) &
                    grepl(
                        "metabolic|biosynthetic|catabolic",
                        as.character(go_data$Description[row_i]),
                        ignore.case = TRUE
                    )
            ]

            if (length(amino_i) == 0) {
                next
            }

            amino_genes <- split_gene_ids(go_data$geneID[amino_i])
            candidate_i <- row_i[broad_metabolism[row_i]]

            for (term_i in candidate_i) {
                term_genes <- split_gene_ids(go_data$geneID[term_i])

                if (length(term_genes) == 0 || length(amino_genes) == 0) {
                    next
                }

                overlap_n <- length(intersect(term_genes, amino_genes))
                overlap_fraction <- overlap_n / length(term_genes)
                union_n <- length(union(term_genes, amino_genes))
                jaccard <- if (union_n == 0) 0 else overlap_n / union_n

                if (
                    overlap_fraction >= go_overlap_fraction_cutoff &&
                        jaccard >= go_overlap_jaccard_cutoff
                ) {
                    go_data$plot_suppressed[term_i] <- TRUE
                    go_data$plot_suppress_reason[term_i] <- paste0(
                        "broad_term_overlaps_amino_acid_family",
                        ":overlap_fraction=",
                        sprintf("%.2f", overlap_fraction),
                        ";jaccard=",
                        sprintf("%.2f", jaccard)
                    )
                }
            }
        }

        go_data <- go_data %>%
            left_join(parent_map, by = "ID")

        parent_term_ids <- unique(as.character(go_data$parent_ID))

        if (length(parent_term_ids) == 1) {
            semantic_cluster <- "GO_cluster_1"
        } else {
            sem_data <- suppressWarnings(
                GOSemSim::godata(
                    OrgDb = org.Dm.eg.db,
                    ont = go_ont,
                    computeIC = TRUE
                )
            )

            sim <- GOSemSim::mgoSim(
                parent_term_ids,
                parent_term_ids,
                semData = sem_data,
                measure = "Wang",
                combine = NULL
            )

            sim[is.na(sim)] <- 0
            diag(sim) <- 1

            hc <- hclust(
                as.dist(1 - sim),
                method = "average"
            )

            cluster_id <- cutree(
                hc,
                h = 1 - similarity_cutoff
            )

            semantic_cluster <- paste0(
                "GO_cluster_",
                cluster_id[match(as.character(go_data$parent_ID), parent_term_ids)]
            )
        }

        go_data$semantic_cluster <- semantic_cluster

        term_meta <- data.frame(
            parent_ID = parent_term_ids,
            parent_depth = vapply(
                parent_term_ids,
                get_go_depth,
                integer(1)
            ),
            descendant_n = vapply(
                parent_term_ids,
                get_go_descendant_n,
                integer(1)
            ),
            stringsAsFactors = FALSE
        )

        representative_rank <- go_data %>%
            group_by(
                semantic_cluster,
                parent_ID,
                parent_Description
            ) %>%
            summarise(
                group_n = n_distinct(group),
                best_FDR = safe_min(p.adjust),
                .groups = "drop"
            ) %>%
            left_join(
                term_meta,
                by = "parent_ID"
            ) %>%
            group_by(semantic_cluster) %>%
            arrange(
                is.na(parent_depth),
                parent_depth,
                desc(descendant_n),
                desc(group_n),
                best_FDR,
                .by_group = TRUE
            ) %>%
            slice_head(n = 1) %>%
            ungroup()

        representative_map <- representative_rank %>%
            transmute(
                semantic_cluster,
                representative_ID = parent_ID,
                representative_Description = parent_Description,
                representative_depth = parent_depth
            )

        go_data <- go_data %>%
            left_join(
                representative_map,
                by = "semantic_cluster"
            )

        go_display <- go_data %>%
            group_by(
                representative_ID,
                representative_Description,
                representative_depth,
                group,
                sex,
                genotype,
                direction,
                comparison,
                background_key
            ) %>%
            summarise(
                geneID = collapse_gene_ids(geneID),
                source_term_n = n_distinct(ID),
                source_term_IDs = paste(sort(unique(ID)), collapse = ";"),
                source_parent_n = n_distinct(parent_ID),
                display_rule = paste(
                    sort(unique(display_rule)),
                    collapse = ";"
                ),
                plot_suppressed = any(plot_suppressed),
                plot_suppress_reason = paste(
                    sort(unique(plot_suppress_reason[nzchar(plot_suppress_reason)])),
                    collapse = ";"
                ),
                input_gene_n = max(input_gene_n),
                background_gene_n = max(background_gene_n),
                best_pvalue = safe_min(pvalue),
                p.adjust = safe_min(p.adjust),
                .groups = "drop"
            ) %>%
            mutate(
                Count = vapply(geneID, count_gene_ids, integer(1)),
                GeneRatio_num = Count / input_gene_n,
                minus_log10_FDR = -log10(p.adjust),
                fdr_basis = "minimum source-term FDR; display collapse only"
            )

        list(
            display = go_display,
            representatives = representative_rank
        )
    }

    go_collapsed <- collapse_go_terms(
        go_all,
        similarity_cutoff = go_similarity_cutoff
    )

    go_display <- go_collapsed$display

    fwrite(
        as.data.table(go_display),
        file.path(
            pairwise_out_dir,
            "RNA_pairwise_GO_parent_collapsed.tsv"
        ),
        sep = "\t"
    )

    fwrite(
        as.data.table(
            go_display %>%
                filter(!plot_suppressed)
        ),
        file.path(
            pairwise_out_dir,
            "RNA_pairwise_GO_parent_plot_input.tsv"
        ),
        sep = "\t"
    )

    fwrite(
        as.data.table(go_collapsed$representatives),
        file.path(
            pairwise_out_dir,
            "RNA_pairwise_GO_parent_representatives.tsv"
        ),
        sep = "\t"
    )

    select_parent_GO <- function(go_data,
                                 direction_use,
                                 top_n = 10) {
        # Broad terms that are gene-set redundant with the amino-acid family are
        # retained in RNA_pairwise_GO_parent_collapsed.tsv but omitted from the
        # final display to prevent an amino-acid-heavy panel.
        go_data <- go_data %>%
            filter(!plot_suppressed)

        term_rank <- go_data %>%
            filter(direction == direction_use) %>%
            group_by(
                representative_ID,
                representative_Description,
                representative_depth
            ) %>%
            summarise(
                group_n = n_distinct(group),
                best_FDR = safe_min(p.adjust),
                source_term_n = max(source_term_n),
                .groups = "drop"
            ) %>%
            arrange(
                is.na(representative_depth),
                representative_depth,
                desc(group_n),
                best_FDR
            )

        # Do not show both a GO term and its ancestor in the final heatmap.
        # The table above remains complete; this only reduces visual redundancy.
        selected_index <- integer(0)
        selected_ids <- character(0)
        if (nrow(term_rank) > 0) {
            for (i in seq_len(nrow(term_rank))) {
                candidate_id <- as.character(term_rank$representative_ID[i])
                redundant <- vapply(
                    selected_ids,
                    function(selected_id) {
                        candidate_id == selected_id ||
                            candidate_id %in% get_go_ancestors(selected_id) ||
                            selected_id %in% get_go_ancestors(candidate_id)
                    },
                    logical(1)
                )
                if (!any(redundant)) {
                    selected_index <- c(selected_index, i)
                    selected_ids <- c(selected_ids, candidate_id)
                }
                if (length(selected_index) >= top_n) {
                    break
                }
            }
        }
        term_rank <- term_rank[selected_index, , drop = FALSE]

        plot_data <- go_data %>%
            semi_join(
                term_rank,
                by = c(
                    "representative_ID",
                    "representative_Description",
                    "representative_depth"
                )
            )

        plot_data$representative_ID <- factor(
            plot_data$representative_ID,
            levels = rev(term_rank$representative_ID)
        )

        list(
            ranking = term_rank,
            plot_data = plot_data
        )
    }

    plot_parent_GO_bubble <- function(selected,
                                      direction_label,
                                      out_file,
                                      width = 9,
                                      height = 6) {
        if (nrow(selected$plot_data) == 0) {
            return(invisible(NULL))
        }

        term_labels <- setNames(
            selected$ranking$representative_Description,
            selected$ranking$representative_ID
        )

        p <- ggplot(
            selected$plot_data,
            aes(
                x = genotype,
                y = representative_ID
            )
        ) +
            geom_point(
                aes(
                    size = Count,
                    color = minus_log10_FDR
                )
            ) +
            facet_grid(
                . ~ sex,
                scales = "free_x",
                space = "free_x"
            ) +
            scale_y_discrete(
                labels = term_labels
            ) +
            scale_size_continuous(
                name = "Gene number",
                range = c(2, 10)
            ) +
            scale_color_gradientn(
                name = "-log10(min. source-term FDR)",
                colours = c(
                    "#3F7FB5",
                    "#8EB5D3",
                    "#F7F7F7",
                    "#D7898F",
                    "#C43D4B"
                )
            ) +
            theme_classic(
                base_size = 14
            ) +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(size = 15),
                axis.text.x = element_text(size = 12, face = "italic", angle = 45, hjust = 1),
                axis.text.y = element_text(size = 11),
                axis.title = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 17)
            ) +
            ggtitle(
                paste0(
                    "BP GO: ",
                    direction_label
                )
            )

        ggsave(
            out_file,
            p,
            width = width,
            height = height
        )

        invisible(p)
    }

    go_parent_up <- select_parent_GO(
        go_display,
        direction_use = "Up",
        top_n = go_top_n
    )

    go_parent_down <- select_parent_GO(
        go_display,
        direction_use = "Down",
        top_n = go_top_n
    )

    plot_parent_GO_bubble(
        go_parent_up,
        direction_label = "Up-regulated genes",
        out_file = file.path(
            pairwise_out_dir,
            "GO_up_shared_top10.pdf"
        ),
        width = 9,
        height = 6
    )

    plot_parent_GO_bubble(
        go_parent_down,
        direction_label = "Down-regulated genes",
        out_file = file.path(
            pairwise_out_dir,
            "GO_down_shared_top10.pdf"
        ),
        width = 11,
        height = 6
    )
}


# ============================================================
# 13. Final summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PAIRWISE DE DONE\n")
cat("============================================================\n")

cat("\nOutput directory:\n")
cat(
    normalizePath(
        pairwise_out_dir,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n"
)

cat("\nPairwise summary:\n")
print(
    pairwise_summary_all
)

cat("\nPaired |log2FC| Wilcoxon:\n")
print(
    wilcox_summary
)


# ============================================================
# Additional violin:
# Up + Down combined
#
# DIRECTLY use the same 4 gene sets as:
# RNA_pairwise_upset_any_DE.pdf
#
#   AdarG_F = AdarG_F_up + AdarG_F_down
#   AdarS_F = AdarS_F_up + AdarS_F_down
#   AdarG_M = AdarG_M_up + AdarG_M_down
#   AdarS_M = AdarS_M_up + AdarS_M_down
#
# Y axis = |log2FoldChange vs WBg|
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})


# ============================================================
# 1. Check required objects
# ============================================================

required_objects_allDEG <- c(
    "pairwise_F",
    "pairwise_M",
    "sets_any_DE",
    "pairwise_out_dir"
)

missing_objects_allDEG <- required_objects_allDEG[
    !vapply(
        required_objects_allDEG,
        exists,
        logical(1),
        inherits = TRUE
    )
]

if (length(missing_objects_allDEG) > 0) {
    stop(
        "Please run the complete pairwise analysis first. Missing objects: ",
        paste(missing_objects_allDEG, collapse = ", ")
    )
}


# ============================================================
# 2. Extract one exact UpSet gene set
# ============================================================

extract_one_anyDE_group <- function(
  gene_set,
  result_table,
  genotype_label,
  sex_label
) {
    result_table <- as.data.table(
        copy(result_table)
    )

    gene_set <- unique(
        as.character(gene_set)
    )

    gene_set <- gene_set[
        !is.na(gene_set) &
            gene_set != ""
    ]


    # Keep EXACTLY genes from the corresponding UpSet set
    out <- result_table[
        gene_id %in% gene_set &
            is.finite(log2FoldChange),
        .(
            gene_id,
            log2FC = log2FoldChange,
            abs_log2FC = abs(log2FoldChange),
            padj
        )
    ]


    # QC:
    # every UpSet gene should be recovered from the Wald result
    missing_gene <- setdiff(
        gene_set,
        out$gene_id
    )

    if (length(missing_gene) > 0) {
        warning(
            genotype_label, " ", sex_label,
            ": ", length(missing_gene),
            " UpSet genes did not have finite log2FC and were not plotted."
        )
    }


    out[
        ,
        `:=`(
            genotype = genotype_label,
            sex = sex_label
        )
    ]

    out
}


# ============================================================
# 3. Build the four groups
#
# These gene sets are EXACTLY the ones used in:
# RNA_pairwise_upset_any_DE.pdf
# ============================================================

AdarG_F_violin <- extract_one_anyDE_group(
    gene_set = sets_any_DE[["AdarG_F"]],
    result_table = pairwise_F$results[["AdarG_vs_WBg"]],
    genotype_label = "AdarG",
    sex_label = "Female"
)


AdarS_F_violin <- extract_one_anyDE_group(
    gene_set = sets_any_DE[["AdarS_F"]],
    result_table = pairwise_F$results[["AdarS_vs_WBg"]],
    genotype_label = "AdarS",
    sex_label = "Female"
)


AdarG_M_violin <- extract_one_anyDE_group(
    gene_set = sets_any_DE[["AdarG_M"]],
    result_table = pairwise_M$results[["AdarG_vs_WBg"]],
    genotype_label = "AdarG",
    sex_label = "Male"
)


AdarS_M_violin <- extract_one_anyDE_group(
    gene_set = sets_any_DE[["AdarS_M"]],
    result_table = pairwise_M$results[["AdarS_vs_WBg"]],
    genotype_label = "AdarS",
    sex_label = "Male"
)


# ============================================================
# 4. Combine
# ============================================================

allDEG_violin_data <- rbindlist(
    list(
        AdarS_F_violin,
        AdarG_F_violin,
        AdarS_M_violin,
        AdarG_M_violin
    ),
    use.names = TRUE,
    fill = TRUE
)


allDEG_violin_data[
    ,
    genotype := factor(
        genotype,
        levels = c(
            "AdarS",
            "AdarG"
        )
    )
]


allDEG_violin_data[
    ,
    sex := factor(
        sex,
        levels = c(
            "Female",
            "Male"
        )
    )
]


# ============================================================
# 5. Strong QC:
# violin n should equal corresponding UpSet-set size
# ============================================================

upset_violin_QC <- data.table(
    sex = c(
        "Female",
        "Female",
        "Male",
        "Male"
    ),
    genotype = c(
        "AdarS",
        "AdarG",
        "AdarS",
        "AdarG"
    ),
    upset_n = c(
        length(unique(sets_any_DE[["AdarS_F"]])),
        length(unique(sets_any_DE[["AdarG_F"]])),
        length(unique(sets_any_DE[["AdarS_M"]])),
        length(unique(sets_any_DE[["AdarG_M"]]))
    )
)


violin_n <- allDEG_violin_data[
    ,
    .(
        violin_n = uniqueN(gene_id)
    ),
    by = .(
        sex,
        genotype
    )
]


upset_violin_QC <- merge(
    upset_violin_QC,
    violin_n,
    by = c(
        "sex",
        "genotype"
    ),
    all.x = TRUE
)


upset_violin_QC[
    ,
    same_gene_number :=
        upset_n == violin_n
]


cat("\n")
cat("============================================================\n")
cat("UpSet vs violin gene-set QC\n")
cat("============================================================\n")

print(
    upset_violin_QC
)


if (
    any(
        upset_violin_QC$same_gene_number == FALSE,
        na.rm = TRUE
    )
) {
    warning(
        "Some violin gene numbers differ from UpSet because some genes have non-finite log2FC."
    )
}


fwrite(
    upset_violin_QC,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_anyDE_UPSET_vs_violin_QC.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 6. Raw summary
# ============================================================

allDEG_summary <- allDEG_violin_data[
    ,
    .(
        n = uniqueN(gene_id),
        mean_abs_log2FC = mean(
            abs_log2FC,
            na.rm = TRUE
        ),
        median_abs_log2FC = median(
            abs_log2FC,
            na.rm = TRUE
        ),
        sd_abs_log2FC = sd(
            abs_log2FC,
            na.rm = TRUE
        )
    ),
    by = .(
        sex,
        genotype
    )
]


cat("\nCombined Up + Down DEG summary:\n")

print(
    allDEG_summary
)


fwrite(
    allDEG_summary,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_anyDE_absolute_log2FC_summary.tsv"
    ),
    sep = "\t"
)


fwrite(
    allDEG_violin_data,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_anyDE_absolute_log2FC_violin_data.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 7. Wilcoxon rank-sum:
# AdarS DEG distribution vs AdarG DEG distribution
#
# Not paired:
# AdarS and AdarG have different DEG sets.
# ============================================================

safe_wilcox_anyDE <- function(
  x,
  y
) {
    x <- x[
        is.finite(x)
    ]

    y <- y[
        is.finite(y)
    ]

    if (
        length(x) < 2 ||
            length(y) < 2
    ) {
        return(
            NA_real_
        )
    }

    suppressWarnings(
        wilcox.test(
            x,
            y,
            paired = FALSE,
            exact = FALSE
        )$p.value
    )
}


allDEG_test <- allDEG_violin_data[
    ,
    .(
        n_AdarS = sum(
            genotype == "AdarS"
        ),
        n_AdarG = sum(
            genotype == "AdarG"
        ),
        p_value = safe_wilcox_anyDE(
            abs_log2FC[
                genotype == "AdarS"
            ],
            abs_log2FC[
                genotype == "AdarG"
            ]
        )
    ),
    by = sex
]


allDEG_test[
    ,
    p_adjust_BH := p.adjust(
        p_value,
        method = "BH"
    )
]


cat("\nAdarS vs AdarG Wilcoxon:\n")

print(
    allDEG_test
)


fwrite(
    allDEG_test,
    file.path(
        pairwise_out_dir,
        "RNA_pairwise_anyDE_AdarS_vs_AdarG_absolute_log2FC_wilcox.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 8. Format P values
# ============================================================

format_p_anyDE <- function(p) {
    if (!is.finite(p)) {
        return("NA")
    }

    if (p < 2.2e-16) {
        return("<2.2e-16")
    }

    if (p < 0.001) {
        return(
            formatC(
                p,
                format = "e",
                digits = 2
            )
        )
    }

    formatC(
        p,
        format = "f",
        digits = 3
    )
}


# ============================================================
# 9. DISPLAY cap
#
# Statistical tests and mean/median above use RAW values.
# Only the displayed violin is capped.
# ============================================================

VIOLIN_CAP_ANYDE <- 0.75


allDEG_violin_data[
    ,
    value_plot := pmin(
        abs_log2FC,
        VIOLIN_CAP_ANYDE
    )
]


# ============================================================
# ============================================================

plot_summary_anyDE <- allDEG_violin_data[
    ,
    .(
        n = .N,
        mean_abs_log2FC = mean(
            abs_log2FC,
            na.rm = TRUE
        ),
        sd_abs_log2FC = sd(
            abs_log2FC,
            na.rm = TRUE
        )
    ),
    by = .(
        sex,
        genotype
    )
]


plot_summary_anyDE[
    ,
    label := sprintf(
        "italic(n) == %d ~ ',' ~ %.3f %%+-%% %.3f",
        n,
        mean_abs_log2FC,
        sd_abs_log2FC
    )
]


# ============================================================
# 11. P annotation
# ============================================================

p_annotation_anyDE <- copy(
    allDEG_test
)


p_annotation_anyDE[
    ,
    label := paste0(
        "italic(P) == ",
        vapply(
            p_value,
            format_p_anyDE,
            character(1)
        )
    )
]


p_annotation_anyDE[
    ,
    `:=`(
        x = 1.5,
        y = 0.60
    )
]


# ============================================================
# 12. Plot
# ============================================================

cols_anyDE <- c(
    "AdarS" = "#F8766D",
    "AdarG" = "#00BFC4"
)


p_anyDE <- ggplot(
    allDEG_violin_data,
    aes(
        x = genotype,
        y = value_plot,
        fill = genotype
    )
) +
    geom_violin(
        trim = TRUE,
        color = NA,
        alpha = 0.75
    ) +
    geom_boxplot(
        width = 0.16,
        outlier.shape = NA,
        color = "black",
        linewidth = 0.30,
        alpha = 0.85
    ) +
    facet_wrap(
        ~sex,
        nrow = 1
    ) +
    geom_text(
        data = plot_summary_anyDE,
        aes(
            x = genotype,
            y = 0.73,
            label = label
        ),
        inherit.aes = FALSE,
        parse = TRUE,
        size = 3.2,
        vjust = 1
    ) +
    geom_text(
        data = p_annotation_anyDE,
        aes(
            x = x,
            y = y,
            label = label
        ),
        inherit.aes = FALSE,
        parse = TRUE,
        size = 4
    ) +
    scale_fill_manual(
        values = cols_anyDE
    ) +
    scale_y_continuous(
        limits = c(
            0,
            VIOLIN_CAP_ANYDE
        ),
        breaks = seq(
            0,
            VIOLIN_CAP_ANYDE,
            by = 0.25
        ),
        expand = expansion(
            mult = c(
                0,
                0.02
            )
        )
    ) +
    labs(
        x = NULL,
        y = expression(
            "|" * log[2] * "FC vs WT" * "|"
        )
    ) +
    theme_classic(
        base_size = 14
    ) +
    theme(
        legend.position = "none",
        strip.background =
            element_blank(),
        strip.text =
            element_text(
                size = 14
            ),
        axis.text.x =
            element_text(
                size = 12
            )
    )


# ============================================================
# 13. Save PDF + PNG
# ============================================================

anyDE_pdf <- file.path(
    pairwise_out_dir,
    "RNA_pairwise_UPSET_anyDE_absolute_log2FC_violin.pdf"
)


anyDE_png <- file.path(
    pairwise_out_dir,
    "RNA_pairwise_UPSET_anyDE_absolute_log2FC_violin.png"
)


ggsave(
    filename = anyDE_pdf,
    plot = p_anyDE,
    width = 7.5,
    height = 5
)


ggsave(
    filename = anyDE_png,
    plot = p_anyDE,
    width = 7.5,
    height = 5,
    dpi = 300
)


message("")
message("Combined Up + Down violin finished.")
message("PDF: ", anyDE_pdf)
message("PNG: ", anyDE_png)
