#!/usr/bin/env Rscript

# Recreate the CDS-site >=2 versus unedited DESeq2 violin plot, replacing the
# mature WT samples with the three male samples from WBG_STAR_featureCounts_new.
# The mature mutant remains the two historical Adar5G1 samples, as in the
# original corrected_5G1_vs_Iso1_w1118 analysis.

suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(ggplot2)
})
Sys.setenv(TZ = "Asia/Shanghai")

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
script_path <- if (length(file_arg)) file_arg[1] else
    "make_DESeq2_CDS_sites_ge2_vs_unedited_WBG_new_male.R"
script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
input_dir <- Sys.getenv("S17_INPUT_DIR", unset = file.path(script_dir, "inputs"))
out_dir <- Sys.getenv("S17_OUT_DIR", unset = file.path(dirname(dirname(script_dir)), "results/07_publication/S17"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fc_dir <- input_dir
adar_fc_file <- file.path(fc_dir, "oldDmel_STAR_featureCounts_new")
wbg_fc_file <- file.path(fc_dir, "WBG_STAR_featureCounts_new")
nascent_fc_file <- file.path(input_dir, "nascent.featureCounts.txt")
nascent_meta_file <- file.path(input_dir, "nascent_sample_metadata.tsv")
cds_site_file <- file.path(input_dir, "ISO1.candidate.AtoG.annotations.long.tsv.new")

id_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
group_levels <- c("unedited", "edited")
group_labels <- c(unedited = "Unedit", edited = "Edited")
group_colors <- c(unedited = "grey70", edited = "#B2182B")

inputs <- c(adar_fc_file, wbg_fc_file, nascent_fc_file, nascent_meta_file, cds_site_file)
invisible(lapply(inputs, function(x) if (!file.exists(x)) stop("Missing input: ", x)))

read_fc <- function(path) {
    x <- fread(path, skip = "Geneid")
    if (!all(c("Geneid", "Length") %in% names(x))) stop("Invalid featureCounts file: ", path)
    x[, `:=`(Geneid = as.character(Geneid), Length = as.numeric(Length))]
    if (anyNA(x$Geneid) || any(x$Geneid == "") || anyDuplicated(x$Geneid)) {
        stop("Missing or duplicated Geneid values: ", path)
    }
    if (any(!is.finite(x$Length) | x$Length <= 0)) stop("Invalid gene lengths: ", path)
    x
}

pick_columns <- function(fc, wanted) {
    candidates <- setdiff(names(fc), id_cols)
    found <- candidates[basename(candidates) %in% wanted]
    found <- found[match(wanted, basename(found))]
    if (length(found) != length(wanted) || anyNA(found)) {
        stop("Missing expected columns: ", paste(wanted, collapse = ", "))
    }
    found
}

make_count_matrix <- function(fc, cols, sample_names) {
    x <- as.matrix(fc[, ..cols])
    storage.mode(x) <- "double"
    if (any(!is.finite(x)) || any(x < 0) || any(abs(x - round(x)) > 1e-8)) {
        stop("Counts must be finite, non-negative integers")
    }
    storage.mode(x) <- "integer"
    rownames(x) <- fc$Geneid
    colnames(x) <- sample_names
    x
}

make_deseq_result <- function(count_mat, case_idx, wt_idx, case_label, wt_label) {
    keep <- rowSums(count_mat >= 10L) >= 3L
    condition <- factor(c(rep(case_label, length(case_idx)), rep(wt_label, length(wt_idx))))
    dds <- DESeqDataSetFromMatrix(
        countData = count_mat[keep, , drop = FALSE],
        colData = data.frame(row.names = colnames(count_mat), condition = condition),
        design = ~condition
    )
    dds <- estimateSizeFactors(dds)
    normalized <- counts(dds, normalized = TRUE)
    out <- data.table(
        gene_id = rownames(normalized),
        mean_case = rowMeans(normalized[, case_idx, drop = FALSE]),
        mean_wt = rowMeans(normalized[, wt_idx, drop = FALSE])
    )
    out[, log2FC := fifelse(
        mean_case > 0 & mean_wt > 0,
        log2(mean_case / mean_wt),
        NA_real_
    )]
    list(
        result = out,
        normalized_counts = normalized,
        size_factors = data.table(
            sample = colnames(normalized),
            condition = as.character(condition),
            size_factor = as.numeric(sizeFactors(dds)),
            all_gene_count_sum = colSums(count_mat),
            retained_gene_count_sum = colSums(count_mat[keep, , drop = FALSE])
        ),
        n_genes_before = nrow(count_mat),
        n_genes_after = sum(keep)
    )
}

format_p_plotmath <- function(p) {
    if (!is.finite(p)) return("italic(P) == NA")
    if (p == 0) return("italic(P) < 2.2 %*% 10^{-16}")
    if (p < 0.001) {
        exponent <- floor(log10(p))
        coefficient <- p / 10^exponent
        return(sprintf("italic(P) == %.2f %%*%% 10^{%d}", coefficient, exponent))
    }
    sprintf("italic(P) == %.3f", p)
}

make_clipped_violin_polygons <- function(d, y_min = -3.25, y_max = 3.05) {
    # Compute each kernel density from every finite observation, then clip only
    # the polygon coordinates. This preserves the full-data density estimate
    # while preventing extreme tails from extending outside the reference panel.
    density_plot <- ggplot(
        d,
        aes(x = as.numeric(edit_group), y = log2FC, group = edit_group)
    ) +
        geom_violin(trim = FALSE, scale = "area", width = 0.9)
    built <- as.data.table(ggplot_build(density_plot)$data[[1]])
    built <- built[y >= y_min & y <= y_max]

    rbindlist(lapply(split(built, built$group), function(z) {
        setorder(z, y)
        half_width <- z$violinwidth * z$width / 2
        data.table(
            polygon_x = c(z$x - half_width, rev(z$x + half_width)),
            polygon_y = c(z$y, rev(z$y)),
            polygon_group = unique(z$group),
            edit_group = factor(group_levels[unique(z$group)], levels = group_levels)
        )
    }))
}

make_violin <- function(d, assay_title, y_axis_label, p_value, panel_tag, tag_x) {
    d <- copy(d[is.finite(log2FC) & !is.na(edit_group)])
    d[, edit_group := factor(edit_group, levels = group_levels)]
    d[, plot_x := as.numeric(edit_group)]
    violin_polygons <- make_clipped_violin_polygons(d)
    ann <- d[, .(
        n = .N,
        mean_log2FC = mean(log2FC),
        sd_log2FC = sd(log2FC)
    ), by = edit_group]
    ann[, `:=`(
        x = as.numeric(edit_group),
        n_y = 3.89,
        stat_y = 3.55,
        n_label = paste0("italic(n) == ", n),
        stat_label = sprintf("%.2f \u00B1 %.2f", mean_log2FC, sd_log2FC)
    )]
    bracket <- data.table(
        x1 = 1, x2 = 2, y = 4.62,
        label = format_p_plotmath(p_value)
    )

    ggplot(d, aes(plot_x, log2FC)) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
        geom_polygon(
            data = violin_polygons,
            aes(x = polygon_x, y = polygon_y, group = polygon_group, fill = edit_group),
            inherit.aes = FALSE, alpha = 0.85, colour = "grey20", linewidth = 0.35
        ) +
        geom_boxplot(aes(group = edit_group),
            width = 0.12, outlier.shape = NA,
            fill = "white", alpha = 0.9, linewidth = 0.35
        ) +
        geom_segment(
            data = bracket, aes(x = x1, xend = x2, y = y, yend = y),
            inherit.aes = FALSE, linewidth = 0.35
        ) +
        geom_segment(
            data = bracket, aes(x = x1, xend = x1, y = y - 0.16, yend = y),
            inherit.aes = FALSE, linewidth = 0.35
        ) +
        geom_segment(
            data = bracket, aes(x = x2, xend = x2, y = y - 0.16, yend = y),
            inherit.aes = FALSE, linewidth = 0.35
        ) +
        geom_text(
            data = bracket,
            aes(x = (x1 + x2) / 2 + 0.037, y = y + 0.18, label = label),
            inherit.aes = FALSE, size = 2.82, family = "Arial", parse = TRUE
        ) +
        geom_text(
            data = ann, aes(x = x, y = n_y, label = n_label),
            inherit.aes = FALSE, size = 2.82, family = "Arial", parse = TRUE
        ) +
        geom_text(
            data = ann, aes(x = x, y = stat_y, label = stat_label),
            inherit.aes = FALSE, size = 2.82, family = "Arial"
        ) +
        scale_x_continuous(
            breaks = 1:2,
            labels = unname(group_labels),
            limits = c(0.38, 2.62),
            expand = expansion(mult = 0)
        ) +
        scale_fill_manual(values = group_colors) +
        scale_y_continuous(breaks = c(-2, 0, 2)) +
        coord_cartesian(ylim = c(-3.25, 3.05), clip = "off") +
        labs(
            x = NULL,
            y = y_axis_label,
            title = assay_title,
            tag = panel_tag
        ) +
        theme_classic(base_size = 8, base_family = "Arial") +
        theme(
            legend.position = "none",
            plot.title = element_text(
                hjust = 0.512, size = 10,
                margin = margin(b = 38)
            ),
            plot.tag = element_text(face = "bold", size = 12),
            plot.tag.position = c(tag_x, 1),
            axis.title.y = element_text(size = 11, margin = margin(r = 1)),
            axis.text = element_text(colour = "black", size = 8),
            axis.text.x = element_text(size = 8.8, margin = margin(t = 1)),
            axis.ticks.x = element_blank(),
            plot.margin = margin(t = 18, r = 12, b = 4, l = 12)
        )
}

# Mature RNA: historical Adar5G1 mutant versus current WBG male WT.
adar_fc <- read_fc(adar_fc_file)
wbg_fc <- read_fc(wbg_fc_file)
if (!identical(adar_fc$Geneid, wbg_fc$Geneid) || !identical(adar_fc$Length, wbg_fc$Length)) {
    stop("Adar5G1 and WBG featureCounts gene rows or lengths differ")
}
adar_cols <- pick_columns(adar_fc, c(
    "Dmel_Adar5G1_rep1_mRNAAligned.sortedByCoord.out.bam",
    "Dmel_Adar5G1_rep2_mRNAAligned.sortedByCoord.out.bam"
))
wbg_cols <- pick_columns(wbg_fc, c(
    "WBG-M1-LJD7603Aligned.sortedByCoord.out.bam",
    "WBG-M2-LJD7604Aligned.sortedByCoord.out.bam",
    "WBG-M3-LJD7605Aligned.sortedByCoord.out.bam"
))
mature_counts <- cbind(
    make_count_matrix(adar_fc, adar_cols, c("Adar5G1_rep1", "Adar5G1_rep2")),
    make_count_matrix(wbg_fc, wbg_cols, c("WBG_M1", "WBG_M2", "WBG_M3"))
)
mature_fit <- make_deseq_result(mature_counts, 1:2, 3:5, "Adar5G1", "WBG_new_male")
mature_base <- mature_fit$result
mature_base[, `:=`(assay = "Mature mRNA", comparison = "Adar5G1 / WBG new male")]

# Mature-sample QC uses the same retained genes and normalized counts as the
# plotted comparison. Correlation and PCA help expose replicate or batch issues.
mature_log_norm <- log2(mature_fit$normalized_counts + 1)
mature_cor <- as.data.table(
    cor(mature_log_norm, method = "spearman"),
    keep.rownames = "sample"
)
mature_pca <- prcomp(t(mature_log_norm), center = TRUE, scale. = FALSE)
mature_pca_scores <- data.table(
    sample = rownames(mature_pca$x),
    PC1 = mature_pca$x[, 1],
    PC2 = mature_pca$x[, 2]
)
mature_pca_variance <- data.table(
    component = colnames(mature_pca$x),
    percent_variance = 100 * mature_pca$sdev^2 / sum(mature_pca$sdev^2)
)

# Nascent RNA is unchanged from the source analysis.
nascent_fc <- read_fc(nascent_fc_file)
meta <- fread(nascent_meta_file)[assay == "nascent", .(sample_id, group, replicate)]
nas_cols <- setdiff(names(nascent_fc), id_cols)
extract_srr <- function(x) {
    m <- regexpr("SRR[0-9]+", basename(x), perl = TRUE)
    if (m[1] < 0) NA_character_ else regmatches(basename(x), m)
}
nas_map <- data.table(
    feature_col = nas_cols,
    sample_id = vapply(nas_cols, extract_srr, character(1))
)
nas_map <- merge(nas_map, meta, by = "sample_id", all.x = TRUE, sort = FALSE)[
    group %in% c("ADAR0", "yw")
]
if (anyNA(nas_map$group) || !any(nas_map$group == "ADAR0") || !any(nas_map$group == "yw")) {
    stop("Could not select ADAR0/yw nascent samples")
}
nas_map[, group_order := match(group, c("ADAR0", "yw"))]
setorder(nas_map, group_order, replicate, sample_id)
nascent_counts <- make_count_matrix(
    nascent_fc,
    nas_map$feature_col,
    make.unique(paste(nas_map$group, nas_map$sample_id, sep = "_"))
)
nas_case_idx <- which(nas_map$group == "ADAR0")
nas_wt_idx <- which(nas_map$group == "yw")
nascent_fit <- make_deseq_result(nascent_counts, nas_case_idx, nas_wt_idx, "ADAR0", "yw")
nascent_base <- nascent_fit$result
nascent_base[, `:=`(assay = "Nascent RNA", comparison = "ADAR0 / yw")]

# Define unedited and CDS-edited groups exactly as in the requested source PDF.
site_rows <- fread(cds_site_file, select = c("handle", "dmel_chr", "dmel_pos", "gene"))
site_rows <- site_rows[
    !is.na(gene) & gene != "" & !is.na(handle) & handle != "",
    .(gene_id = as.character(gene), handle, dmel_chr, dmel_pos)
]
unique_gene_sites <- unique(site_rows, by = c("gene_id", "handle"))
site_counts <- unique_gene_sites[, .(n_sites = .N), by = gene_id]
setorder(site_counts, -n_sites, gene_id)

expression_genes <- union(rownames(mature_counts), rownames(nascent_counts))
gene_groups <- merge(data.table(gene_id = expression_genes), site_counts,
    by = "gene_id", all.x = TRUE
)
gene_groups[is.na(n_sites), n_sites := 0L]
gene_groups[, edit_group := fcase(
    n_sites == 0L, "unedited",
    n_sites >= 2L, "edited",
    default = NA_character_
)]
gene_groups[, edit_group := factor(edit_group, levels = group_levels)]
setorder(gene_groups, gene_id)

nascent_result <- merge(nascent_base, gene_groups, by = "gene_id", all.x = TRUE)
mature_result <- merge(mature_base, gene_groups, by = "gene_id", all.x = TRUE)
plot_data <- rbindlist(list(nascent_result, mature_result), use.names = TRUE)
plot_data <- plot_data[is.finite(log2FC)]

group_summary <- plot_data[!is.na(edit_group), .(
    n = .N,
    mean_log2FC = mean(log2FC),
    median_log2FC = median(log2FC),
    sd_log2FC = sd(log2FC)
), by = .(assay, comparison, edit_group)]

wilcox_results <- plot_data[!is.na(edit_group),
    {
        x <- log2FC[edit_group == "edited"]
        y <- log2FC[edit_group == "unedited"]
        wt <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
        .(
            n_edited = length(x),
            n_unedited = length(y),
            median_edited = median(x),
            median_unedited = median(y),
            median_difference_edited_minus_unedited = median(x) - median(y),
            W = unname(wt$statistic),
            p_value = wt$p.value
        )
    },
    by = .(assay, comparison)
]

filter_summary <- data.table(
    assay = c("Mature mRNA", "Nascent RNA"),
    genes_before_count_filter = c(mature_fit$n_genes_before, nascent_fit$n_genes_before),
    genes_after_count_filter = c(mature_fit$n_genes_after, nascent_fit$n_genes_after),
    genes_with_finite_log2FC = c(
        sum(is.finite(mature_base$log2FC)),
        sum(is.finite(nascent_base$log2FC))
    )
)

get_p <- function(assay_name) wilcox_results[assay == assay_name, p_value][1]
panels <- list(
    make_violin(
        nascent_result,
        "Nascent",
        expression(log[2] * "FC (ADAR0 / yw)"),
        get_p("Nascent RNA"),
        "A",
        0.017
    ),
    make_violin(
        mature_result,
        "Mature",
        expression(log[2] * "FC (Adar5G1 / WBG new male)"),
        get_p("Mature mRNA"),
        "B",
        -0.011
    )
)

output_pdf <- file.path(
    out_dir,
    "FigureS17_Nascent_Mature_rawP.pdf"
)
ggsave(
    output_pdf,
    patchwork::wrap_plots(
        patchwork::wrap_plots(panels, ncol = 2, byrow = TRUE),
        patchwork::plot_spacer(),
        ncol = 1,
        heights = c(0.206, 0.794)
    ),
    width = 595.276 / 72,
    height = 841.890 / 72,
    units = "in",
    device = cairo_pdf
)

message("Wrote: ", output_pdf)
