#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
})

theme_prism_manual <- function(base_size = 11) {
    theme_classic(base_size = base_size) + theme(panel.border = element_rect(colour = "black",
        fill = NA, linewidth = 0.8), panel.grid = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"), axis.ticks.length = unit(0.18, "cm"))
}

input_dir <- Sys.getenv("F3C_TABLE_DIR")
output_dir <- Sys.getenv("F3C_TABLE_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list(Dmel = "Dmel_df_means_AA_density_Asite_mean_dp10_mRNA5", Dsim = "Dsim_df_means_AA_density_Asite_mean_dp10_mRNA5",
    Dsec = "Dsec_df_means_AA_density_Asite_mean_dp10_mRNA5", Dere = "Dere_df_means_AA_density_Asite_mean_dp10_mRNA5")

aa <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
    "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

read_one <- function(sp) {

    x <- as.data.table(readRDS(file.path(input_dir, files[[sp]])))

    if (sp == "Dmel") {
        x <- x[table_name %chin% c("Dmel_ISO1_WBg_F_mean", "Dmel_ISO1_WBg_M_mean", "Dmel_Adar5G1_mean")]
    } else {
        x <- x[grepl(paste0("^", sp, "_[FM]_mean$"), table_name)]
    }

    x[, species := sp]
    x[, plot_group := ifelse(sp == "Dmel" & table_name == "Dmel_Adar5G1_mean", "adar5G1",
        sp)]
    x[, sex := ifelse(grepl("_F_", table_name), "Female", "Male")]
    x
}

x <- rbindlist(lapply(names(files), read_one), fill = TRUE)

available_aa <- intersect(c("edited", aa), names(x))

long <- melt(x, id.vars = c("table_name", "species", "plot_group", "sex"), measure.vars = available_aa,
    variable.name = "AA", value.name = "gmRC")

long[, AA := as.character(AA)]
long[, gmRC := as.numeric(gmRC)]
long <- long[is.finite(gmRC) & gmRC >= 0]  # zero category means retained


# AA mean background
AA_mean <- long[AA %in% aa, .(gmRC = mean(gmRC, na.rm = TRUE)), by = .(table_name, species,
    plot_group, sex)]

AA_mean[, AA := "AA_mean"]


plot_data <- rbindlist(list(long[, .(table_name, species, plot_group, sex, AA, gmRC)], AA_mean[,
    .(table_name, species, plot_group, sex, AA, gmRC)]), fill = TRUE)

plot_data[, AA := as.character(AA)]
fwrite(plot_data, file.path(output_dir, "F3C_plot_data.tsv"), sep = "\t")


# 8 vs 8 paired Wilcoxon
test_data <- rbindlist(list(long[AA == "edited" & plot_group != "adar5G1", .(table_name,
    species, plot_group, sex, AA, gmRC)], AA_mean[plot_group != "adar5G1", .(table_name,
    species, plot_group, sex, AA, gmRC)]))

test_wide <- dcast(test_data, table_name + species + plot_group + sex ~ AA, value.var = "gmRC")

stat <- test_wide[, {
    ok <- is.finite(edited) & is.finite(AA_mean)

    p <- if (sum(ok) >= 2) {
        wilcox.test(edited[ok], AA_mean[ok], paired = TRUE, exact = FALSE)$p.value
    } else {
        NA_real_
    }

    .(n = sum(ok), p = p)
}]

stat[, p_adj := p.adjust(p, method = "BH")]

format_p <- function(p) {
    if (is.na(p)) {
        return("italic(P)==NA")
    }

    if (p == 0) {
        return("italic(P)<2.2%*%10^-308")
    }

    e <- floor(log10(p))
    m <- signif(p/(10^e), 3)

    sprintf("italic(P)==%.3f %%*%% 10^%d", m, e)
}

stat[, label := vapply(p_adj, format_p, character(1))]

fwrite(stat, file.path(output_dir, "AA_gmRC_v10_edited_vs_AAmean_pairedWilcox.tsv"), sep = "\t")


# Median order of canonical AA
median_order <- plot_data[AA %in% aa & plot_group != "adar5G1", .(median_gmRC = median(gmRC,
    na.rm = TRUE)), by = .(AA = as.character(AA))]

median_order <- median_order[order(-median_gmRC)]

plot_order <- c("edited", "AA_mean", as.character(median_order$AA))

plot_data[, AA := factor(AA, levels = plot_order)]


p <- ggplot(plot_data, aes(x = AA, y = log2(gmRC))) + geom_boxplot(data = plot_data[plot_group !=
    "adar5G1"], outlier.shape = NA, fill = "grey85", colour = "black", width = 0.65) + geom_point(aes(colour = plot_group,
    shape = sex), position = position_jitter(width = 0.12), size = 2) + scale_colour_manual(values = c(Dmel = "#F2B37B",
    adar5G1 = "#000000", Dsim = "#D6456C", Dsec = "#7bbcd5", Dere = "#6e7cb8")) + scale_shape_manual(values = c(Female = 16,
    Male = 17)) + labs(x = "Amino acid", y = expression(log[2] * "(mean gmRC)")) + theme_prism_manual()


ymax <- max(log2(plot_data$gmRC), na.rm = TRUE)

p <- p + geom_segment(aes(x = 1, xend = 2, y = ymax + 0.2, yend = ymax + 0.2), inherit.aes = FALSE) +
    geom_text(data = stat, aes(x = 1.5, y = ymax + 0.35, label = label), inherit.aes = FALSE,
        parse = TRUE, size = 4)


pdf(file.path(output_dir, "AA_gmRC_boxplot_v10_DmelWT_adar5G1_Dsim_Dsec_Dere_AAmean_ordered.pdf"),
    width = 11, height = 5)

print(p)
dev.off()

cat("[DONE]\n")
print(plot_order)
