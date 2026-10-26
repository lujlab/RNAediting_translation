#!/usr/bin/env Rscript

# Plot DESeq2 model-estimated log2 fold changes (AdarG or AdarS versus WBg)
# for RNA-editing-change gene groups, separately for CDS sites and all sites.
#
# Group rule used in every panel:
#   AdarG: exclude editing_up
#   AdarS: exclude editing_down
#   retain all editing_nochange groups
#
# Usage:
#   Rscript RNA_DESeq2_log2FC_by_editing_change_groups_CDS_ALL.R [output_root]

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("The patchwork package is required to assemble the four-panel figures.")
}

args <- commandArgs(trailingOnly = TRUE)

default_output_root <- file.path(getwd(), "RNA_DESeq2_log2FC_by_editing_change_groups_results")
output_root <- getOption("deseq2_lfc_output_root", if (length(args) >= 1L) args[[1L]] else default_output_root)

prepared_group_files <- c(CDS = paste0("/path/to/outputs/", "RNA_expression_by_editing_change_groups_CDS_direction_filtered_results/",
    "RNA_expression_and_editing_group_all_genes.tsv"), ALL = paste0("/path/to/outputs/",
    "RNA_expression_by_editing_change_groups_ALL_direction_filtered_results/", "RNA_expression_and_editing_group_all_genes.tsv"))

deseq2_result_dir <- paste0("/path/to/F5_AdarSG/RNA_LRT_degPatterns/", "RNA_pairwise_DE_current")

group_levels <- c("unedit", "editing_up", "editing_down", "editing_nochange")

group_colours <- c(unedit = "#BDBDBD", editing_up = "#F8766D", editing_down = "#00BFC4",
    editing_nochange = "#E6AB02")

group_axis_labels_all <- c(unedit = "Unedited", editing_up = "Editing\nup", editing_down = "Editing\ndown",
    editing_nochange = "Editing\nno change")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

assert_files_exist <- function(paths) {
    missing <- paths[!file.exists(paths)]
    if (length(missing)) {
        stop("Missing required input file(s):\n", paste(missing, collapse = "\n"))
    }
}

assert_columns <- function(x, required, label) {
    missing <- setdiff(required, names(x))
    if (length(missing)) {
        stop(label, " is missing required column(s): ", paste(missing, collapse = ", "))
    }
}

format_p <- function(p) {
    if (!length(p) || !is.finite(p))
        return("NA")
    if (p < 0.001)
        return(formatC(p, format = "e", digits = 2))
    formatC(p, format = "f", digits = 3)
}

safe_kruskal <- function(value, group) {
    keep <- is.finite(value) & !is.na(group)
    value <- value[keep]
    group <- droplevels(factor(group[keep]))
    if (length(value) < 2L || nlevels(group) < 2L)
        return(NA_real_)
    tryCatch(suppressWarnings(kruskal.test(value ~ group)$p.value), error = function(e) NA_real_)
}

cliffs_delta <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    nx <- length(x)
    ny <- length(y)
    if (nx < 1L || ny < 1L)
        return(NA_real_)
    pooled_ranks <- rank(c(x, y), ties.method = "average")
    u_x <- sum(pooled_ranks[seq_len(nx)]) - nx * (nx + 1)/2
    2 * u_x/(nx * ny) - 1
}

cliffs_magnitude <- function(delta) {
    out <- rep(NA_character_, length(delta))
    a <- abs(delta)
    out[is.finite(a) & a < 0.147] <- "negligible"
    out[is.finite(a) & a >= 0.147 & a < 0.33] <- "small"
    out[is.finite(a) & a >= 0.33 & a < 0.474] <- "medium"
    out[is.finite(a) & a >= 0.474] <- "large"
    out
}

run_dunn_test <- function(d, ordered_groups) {
    d <- copy(d[is.finite(log2FoldChange) & !is.na(plot_group)])
    ordered_groups <- ordered_groups[ordered_groups %chin% unique(as.character(d$plot_group))]
    if (length(ordered_groups) < 2L)
        return(data.table())

    d[, dunn_rank := rank(log2FoldChange, ties.method = "average")]
    n_total <- nrow(d)
    tie_sizes <- as.numeric(table(d$log2FoldChange))
    rank_variance <- n_total * (n_total + 1)/12 - sum(tie_sizes^3 - tie_sizes)/(12 * (n_total -
        1))

    rank_summary <- d[, .(n = .N, mean_rank = mean(dunn_rank)), by = plot_group]
    rank_summary[, plot_group := as.character(plot_group)]

    result <- rbindlist(lapply(combn(ordered_groups, 2L, simplify = FALSE), function(pair) {
        x <- d[as.character(plot_group) == pair[[1L]], log2FoldChange]
        y <- d[as.character(plot_group) == pair[[2L]], log2FoldChange]
        rank_1 <- rank_summary[plot_group == pair[[1L]]]
        rank_2 <- rank_summary[plot_group == pair[[2L]]]
        se <- sqrt(rank_variance * (1/rank_1$n + 1/rank_2$n))
        z <- (rank_1$mean_rank - rank_2$mean_rank)/se
        p <- 2 * pnorm(abs(z), lower.tail = FALSE)
        delta <- cliffs_delta(x, y)

        data.table(group_1 = pair[[1L]], group_2 = pair[[2L]], n_group_1 = length(x), n_group_2 = length(y),
            mean_rank_group_1 = rank_1$mean_rank, mean_rank_group_2 = rank_2$mean_rank, dunn_z = z,
            dunn_p = p, median_log2FC_group_1 = median(x), median_log2FC_group_2 = median(y),
            median_log2FC_difference_group1_minus_group2 = median(x) - median(y), fold_change_ratio_of_group_medians = 2^(median(x) -
                median(y)), cliffs_delta_group1_vs_group2 = delta, cliffs_delta_magnitude = cliffs_magnitude(delta))
    }))

    result[, dunn_BH_FDR_within_plot := p.adjust(dunn_p, method = "BH")]
    result[]
}

is_excluded_group <- function(mutant, gene_set_status) {
    (mutant == "AdarG" & gene_set_status == "editing_up") | (mutant == "AdarS" & gene_set_status ==
        "editing_down")
}

# Read the DESeq2 Wald coefficients once. These are the direct model estimates
# from results(dds, contrast = c('genotype', mutant, 'WBg')).
deseq2_files <- unlist(lapply(c("F", "M"), function(sex_use) {
    vapply(c("AdarG", "AdarS"), function(mutant_use) {
        file.path(deseq2_result_dir, paste0("RNA_", sex_use, "_", mutant_use, "_vs_WBg_DESeq2_Wald_all.tsv"))
    }, character(1L))
}))

assert_files_exist(c(prepared_group_files, deseq2_files))

deseq2_lfc <- rbindlist(lapply(c("F", "M"), function(sex_use) {
    rbindlist(lapply(c("AdarG", "AdarS"), function(mutant_use) {
        file_use <- file.path(deseq2_result_dir, paste0("RNA_", sex_use, "_", mutant_use,
            "_vs_WBg_DESeq2_Wald_all.tsv"))
        x <- fread(file_use, check.names = FALSE)
        assert_columns(x, c("gene_id", "sex", "numerator", "denominator", "baseMean", "log2FoldChange",
            "lfcSE", "stat", "pvalue", "padj", "direction"), basename(file_use))
        x <- x[, .(gene_id = as.character(gene_id), sex = as.character(sex), mutant = as.character(numerator),
            denominator = as.character(denominator), baseMean = as.numeric(baseMean), log2FoldChange = as.numeric(log2FoldChange),
            lfcSE = as.numeric(lfcSE), Wald_stat = as.numeric(stat), Wald_pvalue = as.numeric(pvalue),
            Wald_padj = as.numeric(padj), DE_direction = as.character(direction), DESeq2_result_file = basename(file_use))]
        if (any(x$sex != sex_use) || any(x$mutant != mutant_use) || any(x$denominator !=
            "WBg")) {
            stop("Unexpected contrast metadata in ", basename(file_use))
        }
        x
    }))
}))

if (deseq2_lfc[, anyDuplicated(paste(gene_id, sex, mutant))] != 0L) {
    stop("Duplicate gene/sex/mutant keys in the DESeq2 Wald tables.")
}

all_plot_data <- list()
all_group_summaries <- list()
all_global_tests <- list()
all_dunn_tests <- list()
all_membership_summaries <- list()
all_validations <- list()

for (scope_use in c("CDS", "ALL")) {
    scope_out <- file.path(output_root, scope_use)
    dir.create(scope_out, recursive = TRUE, showWarnings = FALSE)

    audit <- fread(prepared_group_files[[scope_use]], check.names = FALSE)
    assert_columns(audit, c("gene_id", "sex", "mutant", "normalization", "gene_set_status",
        "editing_group", "n_target_sites", "mean_site_delta_vs_WBG"), paste(scope_use, "prepared editing-group table"))
    audit <- audit[normalization == "DESeq2"]
    audit <- audit[, .(gene_id = as.character(gene_id), sex = as.character(sex), mutant = as.character(mutant),
        gene_set_status = as.character(gene_set_status), editing_group = as.character(editing_group),
        n_target_sites = as.integer(n_target_sites), mean_site_delta_vs_WBG = as.numeric(mean_site_delta_vs_WBG))]
    audit[, scope := scope_use]

    if (audit[, anyDuplicated(paste(gene_id, sex, mutant))] != 0L) {
        stop("Duplicate gene/sex/mutant keys in the ", scope_use, " audit table.")
    }

    audit[, excluded_by_Adar_direction_rule := is_excluded_group(mutant, gene_set_status)]
    audit[, included_in_plot_and_tests := gene_set_status %chin% group_levels & !excluded_by_Adar_direction_rule]

    membership <- audit[, .(n_genes_in_prepared_DESeq2_table = uniqueN(gene_id)), by = .(scope,
        sex, mutant, gene_set_status, excluded_by_Adar_direction_rule, included_in_plot_and_tests)]

    eligible <- audit[included_in_plot_and_tests == TRUE]
    joined <- merge(eligible, deseq2_lfc, by = c("gene_id", "sex", "mutant"), all.x = TRUE)
    joined[, scope := scope_use]
    joined[, fold_change_Adar_vs_WBg := 2^log2FoldChange]

    missing_or_nonfinite <- joined[!is.finite(log2FoldChange)]
    plot_data <- joined[is.finite(log2FoldChange)]
    plot_data[, plot_group := factor(gene_set_status, levels = group_levels)]
    setorder(plot_data, sex, mutant, plot_group, gene_id)

    membership_after_lfc <- plot_data[, .(n_genes_with_finite_DESeq2_log2FC = uniqueN(gene_id)),
        by = .(scope, sex, mutant, gene_set_status)]
    membership <- merge(membership, membership_after_lfc, by = c("scope", "sex", "mutant",
        "gene_set_status"), all.x = TRUE)
    membership[is.na(n_genes_with_finite_DESeq2_log2FC), n_genes_with_finite_DESeq2_log2FC :=
        0L]

    group_summary <- plot_data[, .(n_genes = uniqueN(gene_id), mean_log2FoldChange = mean(log2FoldChange),
        sd_log2FoldChange = sd(log2FoldChange), median_log2FoldChange = median(log2FoldChange),
        q1_log2FoldChange = quantile(log2FoldChange, 0.25), q3_log2FoldChange = quantile(log2FoldChange,
            0.75), median_fold_change_Adar_vs_WBg = 2^median(log2FoldChange), proportion_log2FC_above_zero = mean(log2FoldChange >
            0), n_DESeq2_padj_below_0_01 = sum(Wald_padj < 0.01, na.rm = TRUE)), by = .(scope,
        sex, mutant, plot_group)]

    global_tests <- plot_data[, .(n_genes = uniqueN(gene_id), kruskal_wallis_p = safe_kruskal(log2FoldChange,
        plot_group)), by = .(scope, sex, mutant)]

    dunn_tests <- rbindlist(lapply(c("F", "M"), function(sex_use) {
        rbindlist(lapply(c("AdarG", "AdarS"), function(mutant_use) {
            d <- plot_data[sex == sex_use & mutant == mutant_use]
            retained_groups <- group_levels[group_levels %chin% unique(as.character(d$plot_group))]
            z <- run_dunn_test(d, retained_groups)
            z[, `:=`(scope = scope_use, sex = sex_use, mutant = mutant_use)]
            setcolorder(z, c("scope", "sex", "mutant", setdiff(names(z), c("scope", "sex",
                "mutant"))))
            z
        }))
    }))

    fwrite(plot_data, file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_plot_data.tsv")),
        sep = "\t")
    fwrite(membership, file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_membership_summary.tsv")),
        sep = "\t")
    fwrite(membership[excluded_by_Adar_direction_rule == TRUE], file.path(scope_out, paste0("RNA_DESeq2_log2FC_",
        scope_use, "_excluded_groups.tsv")), sep = "\t")
    fwrite(group_summary, file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_group_summary.tsv")),
        sep = "\t")
    fwrite(global_tests, file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_Kruskal_Wallis.tsv")),
        sep = "\t")
    fwrite(dunn_tests, file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_Dunn_BH_effect_sizes.tsv")),
        sep = "\t")

    data_min <- min(plot_data$log2FoldChange, na.rm = TRUE)
    data_max <- max(plot_data$log2FoldChange, na.rm = TRUE)
    data_span <- data_max - data_min
    y_unit <- max(data_span * 0.03, 0.2)
    label_y <- data_max + y_unit
    bracket_base <- data_max + 2.35 * y_unit
    bracket_step <- 1.35 * y_unit
    p_y <- bracket_base + 2 * bracket_step + 1.3 * y_unit
    panel_ymax <- p_y + 0.9 * y_unit
    panel_ymin <- data_min - 0.45 * y_unit

    axis_labels <- group_axis_labels_all
    if (scope_use == "CDS") {
        axis_labels[["unedit"]] <- "No CDS\nediting"
    }
    scope_label <- if (scope_use == "CDS")
        "CDS sites" else "all sites"

    make_plot <- function(sex_use, mutant_use) {
        d <- copy(plot_data[sex == sex_use & mutant == mutant_use])
        retained_groups <- group_levels[group_levels %chin% unique(as.character(d$plot_group))]
        d[, plot_group := factor(as.character(plot_group), levels = retained_groups)]

        s <- copy(group_summary[sex == sex_use & mutant == mutant_use])
        s[, plot_group := factor(as.character(plot_group), levels = retained_groups)]
        s[, label := sprintf("italic(n) == %d ~ ',' ~ %.2f %%+-%% %.2f", n_genes, mean_log2FoldChange,
            sd_log2FoldChange)]
        s[, label_y := label_y]

        g <- global_tests[sex == sex_use & mutant == mutant_use]
        a <- copy(dunn_tests[sex == sex_use & mutant == mutant_use])
        a[, `:=`(x_start = match(group_1, retained_groups), x_end = match(group_2, retained_groups))]
        a[, span := x_end - x_start]
        setorder(a, span, x_start)
        a[, annotation_row := seq_len(.N)]
        a[, `:=`(x_mid = (x_start + x_end)/2, bracket_y = bracket_base + (annotation_row -
            1L) * bracket_step)]
        a[, `:=`(bracket_tip_y = bracket_y - 0.2 * y_unit, bracket_label_y = bracket_y +
            0.1 * y_unit, bracket_label = paste0("Dunn BH P = ", vapply(dunn_BH_FDR_within_plot,
            format_p, character(1L))))]

        sex_label <- if (sex_use == "F")
            "Female" else "Male"

        ggplot(d, aes(x = plot_group, y = log2FoldChange, fill = plot_group)) + geom_hline(yintercept = 0,
            linetype = "dashed", linewidth = 0.55, colour = "#555555") + geom_violin(trim = TRUE,
            scale = "width", colour = NA, alpha = 0.75, width = 0.92) + geom_boxplot(width = 0.17,
            outlier.shape = NA, colour = "#333333", linewidth = 0.45, alpha = 0.85) + geom_text(data = s,
            aes(x = plot_group, y = label_y, label = label), inherit.aes = FALSE, parse = TRUE,
            size = 3.15, vjust = 0.5) + geom_segment(data = a, aes(x = x_start, xend = x_end,
            y = bracket_y, yend = bracket_y), inherit.aes = FALSE, linewidth = 0.45, colour = "black") +
            geom_segment(data = a, aes(x = x_start, xend = x_start, y = bracket_y, yend = bracket_tip_y),
                inherit.aes = FALSE, linewidth = 0.45, colour = "black") + geom_segment(data = a,
            aes(x = x_end, xend = x_end, y = bracket_y, yend = bracket_tip_y), inherit.aes = FALSE,
            linewidth = 0.45, colour = "black") + geom_text(data = a, aes(x = x_mid, y = bracket_label_y,
            label = bracket_label), inherit.aes = FALSE, size = 2.95, vjust = 0) + annotate("text",
            x = (length(retained_groups) + 1)/2, y = p_y, label = paste0("Kruskal-Wallis P = ",
                format_p(g$kruskal_wallis_p)), size = 4) + scale_fill_manual(values = group_colours[retained_groups],
            drop = TRUE) + scale_x_discrete(labels = axis_labels[retained_groups], drop = TRUE) +
            scale_y_continuous(breaks = pretty(c(data_min, data_max), n = 7), expand = expansion(mult = c(0,
                0))) + coord_cartesian(ylim = c(panel_ymin, panel_ymax), clip = "off") +
            labs(title = sex_label, subtitle = paste(mutant_use, "vs WBg | DESeq2 |", scope_label),
                x = NULL, y = expression(log[2] ~ fold ~ change ~ ("Adar"/"WBg"))) + theme_classic(base_size = 14) +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16,
                face = "plain"), plot.subtitle = element_text(hjust = 0.5, size = 11), axis.text.x = element_text(size = 11,
                colour = "black"), axis.text.y = element_text(size = 11, colour = "black"),
                axis.title.y = element_text(size = 13), axis.line = element_line(linewidth = 0.75,
                  colour = "black"), axis.ticks = element_line(linewidth = 0.65, colour = "black"),
                plot.margin = margin(8, 12, 8, 8))
    }

    plot_list <- list()
    for (sex_use in c("F", "M")) {
        for (mutant_use in c("AdarG", "AdarS")) {
            key <- paste(sex_use, mutant_use, sep = "_")
            p <- make_plot(sex_use, mutant_use)
            plot_list[[key]] <- p
            prefix <- paste0("RNA_DESeq2_log2FC_", scope_use, "_", ifelse(sex_use == "F",
                "Female", "Male"), "_", mutant_use, "_by_editing_change_group")
            ggsave(file.path(scope_out, paste0(prefix, ".pdf")), p, width = 8.6, height = 6.2,
                device = cairo_pdf)
            ggsave(file.path(scope_out, paste0(prefix, ".png")), p, width = 8.6, height = 6.2,
                dpi = 320, bg = "white")
        }
    }

    combined <- patchwork::wrap_plots(plot_list[["F_AdarG"]], plot_list[["F_AdarS"]], plot_list[["M_AdarG"]],
        plot_list[["M_AdarS"]], ncol = 2)
    ggsave(file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_all_four_panels.pdf")),
        combined, width = 17.2, height = 12.4, device = cairo_pdf)
    ggsave(file.path(scope_out, paste0("RNA_DESeq2_log2FC_", scope_use, "_all_four_panels.png")),
        combined, width = 17.2, height = 12.4, dpi = 260, bg = "white")

    png_files <- list.files(scope_out, pattern = "\\.png$", full.names = TRUE)
    pdf_files <- list.files(scope_out, pattern = "\\.pdf$", full.names = TRUE)
    forbidden_rows <- plot_data[is_excluded_group(mutant, gene_set_status), .N]
    nochange_expected <- audit[gene_set_status == "editing_nochange", uniqueN(paste(sex,
        mutant))]
    nochange_observed <- plot_data[gene_set_status == "editing_nochange", uniqueN(paste(sex,
        mutant))]

    validation <- data.table(scope = scope_use, check = c("no duplicate prepared gene/sex/mutant keys",
        "all eligible genes matched a finite DESeq2 log2FoldChange", "no forbidden Adar-direction groups enter plots",
        "all editing_nochange panel groups are retained", "all expected PNG files are nonempty",
        "all expected PDF files are nonempty"), observed = c(as.character(audit[, anyDuplicated(paste(gene_id,
        sex, mutant))]), paste0(nrow(plot_data), "/", nrow(eligible), " genes"), paste0(forbidden_rows,
        " rows"), paste0(nochange_observed, "/", nochange_expected, " panel groups"), paste0(length(png_files),
        " files"), paste0(length(pdf_files), " files")), pass = c(audit[, anyDuplicated(paste(gene_id,
        sex, mutant))] == 0L, nrow(missing_or_nonfinite) == 0L && nrow(plot_data) == nrow(eligible),
        forbidden_rows == 0L, nochange_observed == nochange_expected, length(png_files) ==
            5L && all(file.info(png_files)$size > 0), length(pdf_files) == 5L && all(file.info(pdf_files)$size >
            0)))
    invisible(NULL)
    if (!all(validation$pass)) {
        stop("Validation failed for ", scope_use, ". See its validation table.")
    }

    readme <- c(paste0("DESeq2 log2 fold change by editing-change group [", scope_use, "]"),
        "===============================================================", "", "Y axis: the DESeq2 Wald model log2FoldChange for AdarG or AdarS versus WBg.",
        "The coefficient comes directly from the supplied DESeq2 results table; it is", "not a ratio calculated from normalized-count means and is not LFC-shrunk.",
        "A value of 0 means no modeled expression change, +1 means 2-fold expression,", "and -1 means one-half expression in the named Adar genotype relative to WBg.",
        "", paste0("Editing-site scope: ", ifelse(scope_use == "CDS", "CDS only", "all sites"),
            "."), "AdarG editing_up and AdarS editing_down are excluded from plots and tests.",
        "Every editing_nochange group is retained.", "", "Each panel uses a Kruskal-Wallis global test followed by pairwise Dunn tests.",
        "Dunn P values are BH-adjusted within the panel. Cliff's delta and median", "log2FC differences are reported in the pairwise table.")
    invisible(NULL)

    all_plot_data[[scope_use]] <- plot_data
    all_group_summaries[[scope_use]] <- group_summary
    all_global_tests[[scope_use]] <- global_tests
    all_dunn_tests[[scope_use]] <- dunn_tests
    all_membership_summaries[[scope_use]] <- membership
    all_validations[[scope_use]] <- validation
}

fwrite(rbindlist(all_group_summaries), file.path(output_root, "RNA_DESeq2_log2FC_CDS_ALL_group_summary.tsv"),
    sep = "\t")
fwrite(rbindlist(all_global_tests), file.path(output_root, "RNA_DESeq2_log2FC_CDS_ALL_Kruskal_Wallis.tsv"),
    sep = "\t")
fwrite(rbindlist(all_dunn_tests), file.path(output_root, "RNA_DESeq2_log2FC_CDS_ALL_Dunn_BH_effect_sizes.tsv"),
    sep = "\t")
fwrite(rbindlist(all_membership_summaries, fill = TRUE), file.path(output_root, "RNA_DESeq2_log2FC_CDS_ALL_membership_summary.tsv"),
    sep = "\t")
invisible(NULL)

# ============================================================
# FINAL A4 FIGURE: match the supplied 2-column x 4-row example
# ============================================================

compact_letters_from_dunn <- function(groups, pairwise, alpha = 0.05) {
    groups <- as.character(groups)
    n_groups <- length(groups)
    nonsignificant <- diag(TRUE, n_groups)
    dimnames(nonsignificant) <- list(groups, groups)

    if (nrow(pairwise)) {
        for (i in seq_len(nrow(pairwise))) {
            g1 <- pairwise$group_1[[i]]
            g2 <- pairwise$group_2[[i]]
            if (g1 %chin% groups && g2 %chin% groups) {
                is_ns <- is.finite(pairwise$dunn_BH_FDR_within_plot[[i]]) && pairwise$dunn_BH_FDR_within_plot[[i]] >=
                  alpha
                nonsignificant[g1, g2] <- is_ns
                nonsignificant[g2, g1] <- is_ns
            }
        }
    }

    candidate_cliques <- lapply(seq_len(2^n_groups - 1L), function(mask) {
        members <- which(as.logical(intToBits(mask)[seq_len(n_groups)]))
        if (all(nonsignificant[members, members, drop = FALSE]))
            members else NULL
    })
    candidate_cliques <- Filter(Negate(is.null), candidate_cliques)

    maximal <- vapply(seq_along(candidate_cliques), function(i) {
        !any(vapply(seq_along(candidate_cliques), function(j) {
            i != j && length(candidate_cliques[[j]]) > length(candidate_cliques[[i]]) &&
                all(candidate_cliques[[i]] %in% candidate_cliques[[j]])
        }, logical(1L)))
    }, logical(1L))
    cliques <- candidate_cliques[maximal]
    clique_order <- order(vapply(cliques, min, integer(1L)), -vapply(cliques, length, integer(1L)))
    cliques <- cliques[clique_order]

    letter_pool <- c(letters, paste0("a", letters), paste0("b", letters))
    group_letters <- rep("", n_groups)
    for (i in seq_along(cliques)) {
        group_letters[cliques[[i]]] <- paste0(group_letters[cliques[[i]]], letter_pool[[i]])
    }

    data.table(plot_group = groups, significance_letter = group_letters)
}

make_compact_A4_panel <- function(scope_use, sex_use, mutant_use) {
    d <- copy(all_plot_data[[scope_use]][sex == sex_use & mutant == mutant_use])
    s <- copy(all_group_summaries[[scope_use]][sex == sex_use & mutant == mutant_use])
    pw <- copy(all_dunn_tests[[scope_use]][sex == sex_use & mutant == mutant_use])

    retained_groups <- group_levels[group_levels %chin% unique(as.character(d$plot_group))]
    d[, plot_group := factor(as.character(plot_group), levels = retained_groups)]
    d[, plotted_log2FoldChange := pmax(-0.75, pmin(0.75, log2FoldChange))]
    s[, plot_group := as.character(plot_group)]

    s[, original_group_order := match(plot_group, group_levels)]
    groups_high_to_low <- s[order(-median_log2FoldChange, original_group_order), plot_group]
    letters_one <- compact_letters_from_dunn(groups_high_to_low, pw, alpha = 0.05)
    s <- merge(s, letters_one, by = "plot_group", all.x = TRUE)
    s[, plot_group := factor(plot_group, levels = retained_groups)]
    s[, display_mean_log2FC := fifelse(abs(mean_log2FoldChange) < 0.005, 0, mean_log2FoldChange)]
    s[, top_label := sprintf("atop(italic(n) == %d, %.2f %%+-%% %.2f)", n_genes, display_mean_log2FC,
        sd_log2FoldChange)]

    # Display-only winsorization requested by the user. Tests, compact letters,
    # effect sizes, and the mean +/- SD labels continue to use raw log2FoldChange.
    y_low_core <- -0.75
    y_high_core <- 0.75
    letter_y <- 0.82
    n_label_y <- 0.99

    s[, `:=`(letter_y = letter_y, n_label_y = n_label_y)]

    axis_labels <- group_axis_labels_all
    if (scope_use == "CDS") {
        axis_labels[["unedit"]] <- "No CDS\nediting"
    }
    scope_label <- if (scope_use == "CDS")
        "CDS sites" else "all sites"
    sex_label <- if (sex_use == "F")
        "Female" else "Male"
    adar_superscript <- if (mutant_use == "AdarG")
        "G" else "S"
    panel_title <- bquote(italic(Adar)^.(adar_superscript) ~ .(sex_label) ~ "|" ~ .(scope_label))

    ggplot(d, aes(x = plot_group, y = plotted_log2FoldChange, fill = plot_group)) + geom_hline(yintercept = 0,
        linetype = "dashed", linewidth = 0.28, colour = "#666666") + geom_violin(trim = TRUE,
        scale = "width", colour = NA, alpha = 0.78, width = 0.92) + geom_boxplot(width = 0.17,
        outlier.shape = NA, colour = "#333333", linewidth = 0.3, alpha = 0.62) + geom_text(data = s,
        aes(x = plot_group, y = n_label_y, label = top_label), inherit.aes = FALSE, parse = TRUE,
        size = 2.3, vjust = 0.5, lineheight = 0.92) + geom_text(data = s, aes(x = plot_group,
        y = letter_y, label = significance_letter), inherit.aes = FALSE, size = 3, fontface = "bold",
        vjust = 0.5) + scale_fill_manual(values = group_colours[retained_groups], drop = TRUE) +
        scale_x_discrete(labels = axis_labels[retained_groups], drop = TRUE) + scale_y_continuous(breaks = seq(-0.75,
        0.75, by = 0.25), expand = expansion(mult = c(0, 0))) + coord_cartesian(ylim = c(y_low_core,
        y_high_core), clip = "off") + labs(title = panel_title, x = NULL, y = expression(DESeq2 ~
        log[2] ~ fold ~ change ~ ("Adar"/"WBg"))) + theme_classic(base_size = 9) + theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 9.5, face = "plain", margin = margin(b = 27,
            unit = "pt")), axis.title.y = element_text(size = 8.5, margin = margin(r = 2,
            unit = "pt")), axis.text.x = element_text(size = 8.2, colour = "black", lineheight = 0.88),
        axis.text.y = element_text(size = 8, colour = "black"), axis.line = element_line(linewidth = 0.45,
            colour = "black"), axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks.length = unit(2.2, "pt"), plot.margin = margin(t = 3, r = 3, b = 3, l = 3,
            unit = "pt"))
}

a4_panels <- list(make_compact_A4_panel("CDS", "F", "AdarG"), make_compact_A4_panel("CDS",
    "F", "AdarS"), make_compact_A4_panel("CDS", "M", "AdarG"), make_compact_A4_panel("CDS",
    "M", "AdarS"), make_compact_A4_panel("ALL", "F", "AdarG"), make_compact_A4_panel("ALL",
    "F", "AdarS"), make_compact_A4_panel("ALL", "M", "AdarG"), make_compact_A4_panel("ALL",
    "M", "AdarS"))

a4_grid <- patchwork::wrap_plots(a4_panels, ncol = 2, byrow = TRUE)

# The example uses two compact columns occupying roughly two thirds of an A4
# portrait page, with a reserved blank column on the right.
a4_page <- (a4_grid | patchwork::plot_spacer()) + patchwork::plot_layout(widths = c(2.05,
    0.95)) & theme(plot.margin = margin(2, 2, 2, 2, unit = "mm"))

a4_pdf <- file.path(output_root, "RNA_DESeq2_log2FC_CDS_all_sites_A4_portrait_2col_4row_9pt.pdf")
a4_png <- file.path(output_root, "RNA_DESeq2_log2FC_CDS_all_sites_A4_portrait_2col_4row_9pt.png")

ggsave(a4_pdf, a4_page, width = 210, height = 297, units = "mm", device = cairo_pdf, bg = "white")
ggsave(a4_png, a4_page, width = 210, height = 297, units = "mm", dpi = 220, bg = "white")

a4_letters <- rbindlist(lapply(c("CDS", "ALL"), function(scope_use) {
    rbindlist(lapply(c("F", "M"), function(sex_use) {
        rbindlist(lapply(c("AdarG", "AdarS"), function(mutant_use) {
            d <- all_plot_data[[scope_use]][sex == sex_use & mutant == mutant_use]
            pw <- all_dunn_tests[[scope_use]][sex == sex_use & mutant == mutant_use]
            s <- all_group_summaries[[scope_use]][sex == sex_use & mutant == mutant_use]
            retained_groups <- group_levels[group_levels %chin% unique(as.character(d$plot_group))]
            s[, plot_group := as.character(plot_group)]
            s[, original_group_order := match(plot_group, group_levels)]
            groups_high_to_low <- s[order(-median_log2FoldChange, original_group_order),
                plot_group]
            x <- compact_letters_from_dunn(groups_high_to_low, pw, alpha = 0.05)
            x[, `:=`(scope = scope_use, sex = sex_use, mutant = mutant_use)]
            x <- merge(x, s[, .(plot_group, median_log2FoldChange)], by = "plot_group", all.x = TRUE)
            x[, letter_order_basis := "median raw DESeq2 log2FoldChange, descending"]
            setcolorder(x, c("scope", "sex", "mutant", "plot_group", "median_log2FoldChange",
                "significance_letter", "letter_order_basis"))
            x
        }))
    }))
}))
fwrite(a4_letters, file.path(output_root, "RNA_DESeq2_log2FC_CDS_ALL_compact_letters.tsv"),
    sep = "\t")

a4_axis_cap_counts <- rbindlist(lapply(c("CDS", "ALL"), function(scope_use) {
    all_plot_data[[scope_use]][, .(n_total = .N, n_capped_at_minus_0_75 = sum(log2FoldChange <
        -0.75), n_capped_at_plus_0_75 = sum(log2FoldChange > 0.75)), by = .(scope, sex, mutant,
        plot_group)]
}))
fwrite(a4_axis_cap_counts, file.path(output_root, "RNA_DESeq2_log2FC_CDS_ALL_axis_cap_counts.tsv"),
    sep = "\t")

if (!file.exists(a4_pdf) || file.info(a4_pdf)$size <= 0) {
    stop("The final A4 PDF was not created correctly.")
}

cat("Completed. Output root:\n", normalizePath(output_root), "\n", sep = "")
