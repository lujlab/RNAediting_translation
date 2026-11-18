#!/usr/bin/env Rscript
# Join CDS editing groups to model-estimated log2FC and generate Figure S26.
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 03_plot.R deseq2_dir groups_dir output_dir")
de_dir <- args[1]
group_dir <- args[2]
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "plot_helpers.R"))
groups <- fread(file.path(group_dir, "05_gene_editing_change_groups.tsv"))
background <- fread(file.path(group_dir, "candidate_CDS_genes.tsv"))
stopifnot(!anyDuplicated(groups[, .(gene_id, sex, mutant)]))
lfc <- rbindlist(lapply(c("F", "M"), function(sx) rbindlist(lapply(c("AdarG", "AdarS"), function(mut) {
    d <- fread(file.path(de_dir, paste0("RNA_", sx, "_", mut, "_vs_WBg.tsv")))
    stopifnot(all(d$sex == sx), all(d$mutant == mut), all(d$denominator == "WBg"))
    d
}))))
stopifnot(!anyDuplicated(lfc[, .(gene_id, sex, mutant)]))
membership <- merge(lfc, groups, by = c("gene_id", "sex", "mutant"), all.x = TRUE)
background[, has_CDS_candidate := TRUE]
membership <- merge(membership, background, by = c("gene_id", "sex"), all.x = TRUE)
membership[, gene_set_status := fcase(!is.na(editing_group), editing_group, is.na(has_CDS_candidate),
    "unedit", default = "edited_not_in_target_sites")]
group_levels <- c("unedit", "editing_up", "editing_down", "editing_nochange")
membership[, included := gene_set_status %chin% group_levels & !(mutant == "AdarG" & gene_set_status ==
    "editing_up") & !(mutant == "AdarS" & gene_set_status == "editing_down") & is.finite(log2FoldChange)]
fwrite(membership, file.path(out, "gene_membership.tsv"), sep = "\t")
plot_data <- membership[included == TRUE]
plot_data[, plot_group := factor(gene_set_status, levels = group_levels)]
setorder(plot_data, sex, mutant, plot_group, gene_id)
fwrite(plot_data, file.path(out, "plot_data.tsv"), sep = "\t")
group_summary <- plot_data[, .(n_genes = .N, mean_log2FoldChange = mean(log2FoldChange),
    sd_log2FoldChange = sd(log2FoldChange), median_log2FoldChange = median(log2FoldChange)),
    by = .(sex, mutant, plot_group)]
kw <- plot_data[, .(kruskal_wallis_p = kruskal.test(log2FoldChange ~ plot_group)$p.value),
    by = .(sex, mutant)]

# Dunn's two-sided normal approximation includes the pooled-rank tie correction.
run_dunn <- function(d) {
    d <- copy(d)
    d[, pooled_rank := rank(log2FoldChange, ties.method = "average")]
    n <- nrow(d)
    ties <- as.numeric(table(d$log2FoldChange))
    variance <- n * (n + 1)/12 - sum(ties^3 - ties)/(12 * (n - 1))
    rs <- d[, .(n = .N, mean_rank = mean(pooled_rank)), by = plot_group]
    rs[, plot_group := as.character(plot_group)]
    retained <- group_levels[group_levels %chin% rs$plot_group]
    stopifnot(length(retained) >= 2, is.finite(variance), variance > 0)
    ans <- rbindlist(lapply(combn(retained, 2, simplify = FALSE), function(pair) {
        a <- rs[plot_group == pair[1]]
        b <- rs[plot_group == pair[2]]
        z <- (a$mean_rank - b$mean_rank)/sqrt(variance * (1/a$n + 1/b$n))
        data.table(group_1 = pair[1], group_2 = pair[2], dunn_z = z, dunn_p = 2 * pnorm(abs(z),
            lower.tail = FALSE))
    }))
    ans[, dunn_BH_FDR_within_plot := p.adjust(dunn_p, method = "BH")]
    ans
}
dunn <- plot_data[, run_dunn(.SD), by = .(sex, mutant)]
fwrite(group_summary, file.path(out, "group_summary.tsv"), sep = "\t")
fwrite(kw, file.path(out, "Kruskal_Wallis.tsv"), sep = "\t")
fwrite(dunn, file.path(out, "Dunn_BH.tsv"), sep = "\t")
letters_table <- rbindlist(lapply(c("F", "M"), function(sx) rbindlist(lapply(c("AdarG", "AdarS"),
    function(mut) {
        s <- copy(group_summary[sex == sx & mutant == mut])
        s[, original_order := match(plot_group, group_levels)]
        order_groups <- s[order(-median_log2FoldChange, original_order), as.character(plot_group)]
        ans <- compact_letters_from_dunn(order_groups, dunn[sex == sx & mutant == mut])
        ans[, `:=`(sex = sx, mutant = mut)]
        ans
    }))))
fwrite(letters_table, file.path(out, "significance_letters.tsv"), sep = "\t")
group_colours <- c(unedit = "#BDBDBD", editing_up = "#F8766D", editing_down = "#00BFC4",
    editing_nochange = "#E6AB02")
group_axis_labels <- c(unedit = "No CDS\nediting", editing_up = "Editing\nup", editing_down = "Editing\ndown",
    editing_nochange = "Editing\nno change")
panels <- list(make_panel("F", "AdarG"), make_panel("F", "AdarS"), make_panel("M", "AdarG"),
    make_panel("M", "AdarS"))
p <- wrap_plots(panels, ncol = 2, byrow = TRUE)
ggsave(file.path(out, "FigureS26.pdf"), p, width = 143.5, height = 148.5, units = "mm", device = cairo_pdf,
    bg = "white")
fwrite(plot_data[, .(n = .N, n_below_minus075 = sum(log2FoldChange < -0.75), n_above_plus075 = sum(log2FoldChange >
    0.75)), by = .(sex, mutant, plot_group)], file.path(out, "display_caps.tsv"), sep = "\t")
invisible(NULL)
