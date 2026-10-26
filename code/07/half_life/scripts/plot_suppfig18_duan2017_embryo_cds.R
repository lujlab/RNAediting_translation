#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1L) args[[1L]] else "output/duan2017_embryo_cds_min1"
output_pdf <- if (length(args) >= 2L) args[[2L]] else
    "output/Supplementary Figure 18 Half-life CDS sites min1.pdf"

min_sites <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1L
if (is.na(min_sites) || min_sites < 1L) stop("min_sites must be a positive integer")

font_family <- "Arial"
red <- "#E45756"
blue <- "#4C78A8"
green <- "#54A24B"
light_grid <- "#DEDEDE"

meta <- fread(file.path(analysis_dir, "10_tables", "all_genes_merged_metadata.tsv"))
fisher <- fread(file.path(
    analysis_dir, "03_2010_stable_decay", "2010_stable_decay_fisher.tsv"
))[1]
reg <- fread(file.path(
    analysis_dir, "05_2010_regression", "2010_half_life_regression_all_models.tsv"
))

meta <- meta[editing_binary_status %chin% c("Edited", "No_detected_editing")]
meta[, plot_group := factor(
    fifelse(editing_binary_status == "Edited", "Edited", "Unedited"),
    levels = c("Edited", "Unedited")
)]
decay <- meta[
    decay_class %chin% c("II", "III", "IV", "V") &
        is.finite(half_life_min) & is.finite(transcript_length)
]

# Panel A includes all class I-V genes with a positive finite transcript length.
# Half-life availability is not an eligibility condition for this panel.
length_genes <- meta[decay_class %chin% c("I", "II", "III", "IV", "V") &
    is.finite(transcript_length) & transcript_length > 0]
n_length_edited <- length_genes[plot_group == "Edited", .N]
n_length_unedited <- length_genes[plot_group == "Unedited", .N]
n_length_edited_label <- comma(n_length_edited, accuracy = 1)
n_length_unedited_label <- comma(n_length_unedited, accuracy = 1)
fwrite(length_genes, file.path(analysis_dir, "10_tables", "FigureS18_panel_A_genes.tsv"), sep = "\t")
fwrite(
    data.table(
        panel = c("A", "A", "B", "B"), group = rep(c("Edited", "Unedited"), 2),
        n = c(n_length_edited, n_length_unedited, fisher$N_edited, fisher$N_control)
    ),
    file.path(analysis_dir, "10_tables", "FigureS18_panel_AB_counts.tsv"),
    sep = "\t"
)

n_edited <- decay[plot_group == "Edited", .N]
n_unedited <- decay[plot_group == "Unedited", .N]
n_edited_label <- comma(n_edited, accuracy = 1)
n_unedited_label <- comma(n_unedited, accuracy = 1)

theme_ac <- function() {
    theme_classic(base_family = font_family, base_size = 8) +
        theme(
            text = element_text(family = font_family, colour = "black"),
            axis.title = element_text(size = 7.2),
            axis.text = element_text(size = 6.4, colour = "black"),
            axis.line = element_line(linewidth = 0.45, colour = "black"),
            axis.ticks = element_line(linewidth = 0.35, colour = "black"),
            axis.ticks.length = unit(1.5, "mm"),
            plot.margin = margin(4, 3, 1, 2, unit = "pt"),
            legend.position = "none"
        )
}

add_italic_n_axis <- function(p) {
    # Grid text is used for the second line so the sample-size n is explicitly
    # italic while group names and numbers remain roman, as in the reference PDF.
    p + scale_x_discrete(labels = c("Edited", "Unedited")) +
        theme(axis.text.x = element_text(margin = margin(t = 2)))
}

length_stats <- length_genes[, .(
    mean = mean(transcript_length), SD = sd(transcript_length),
    lower = mean(transcript_length) - sd(transcript_length),
    upper = mean(transcript_length) + sd(transcript_length)
), by = plot_group]
length_stats[, hjust := fifelse(plot_group == "Edited", 0.50, 0.48)]
length_stats[, annotation := sprintf(
    "%s \u00B1 %s nt", comma(round(mean), accuracy = 1), comma(round(SD), accuracy = 1)
)]

pA <- ggplot(length_genes, aes(plot_group, transcript_length, fill = plot_group)) +
    geom_hline(yintercept = 3000, colour = light_grid, linewidth = 0.35) +
    geom_violin(trim = FALSE, linewidth = 0.45, colour = "black") +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.4, alpha = 0.35) +
    geom_errorbar(
        data = length_stats, aes(x = plot_group, ymin = lower, ymax = upper),
        inherit.aes = FALSE, width = 0.09, linewidth = 0.55
    ) +
    geom_point(
        data = length_stats, aes(x = plot_group, y = mean), inherit.aes = FALSE,
        shape = 23, size = 3.0, stroke = 0.55, fill = "white", colour = "black"
    ) +
    geom_text(
        data = length_stats, aes(x = plot_group, y = 60000, label = annotation, hjust = hjust),
        inherit.aes = FALSE, family = font_family, size = 7.0 / ggplot2::.pt
    ) +
    scale_fill_manual(values = c(Edited = red, Unedited = blue)) +
    scale_y_log10(
        limits = c(200, 100000), breaks = c(300, 1000, 3000, 10000, 30000),
        labels = label_number(big.mark = ",")
    ) +
    labs(x = NULL, y = "Transcript length (nt, log scale)") +
    theme_ac()
pA <- add_italic_n_axis(pA)

stack <- data.table(
    plot_group = factor(c("Edited", "Edited", "Unedited", "Unedited"),
        levels = c("Edited", "Unedited")
    ),
    status = factor(c("Stable", "Decay", "Stable", "Decay"),
        levels = c("Decay", "Stable")
    ),
    n = c(
        fisher$Edited_stable, fisher$Edited_decay,
        fisher$Control_stable, fisher$Control_decay
    )
)
stack[, fraction := n / sum(n), by = plot_group]
stack[, midpoint := cumsum(fraction) - fraction / 2, by = plot_group]

p_label <- sprintf("%.3g", fisher$one_sided_greater_P)
pB <- ggplot(stack, aes(plot_group, fraction, fill = status)) +
    geom_col(width = 0.72) +
    geom_text(
        aes(y = midpoint + 0.045, label = sprintf("italic(n) == \"%s\"", comma(n, accuracy = 1))),
        parse = TRUE, family = font_family, size = 6.5 / ggplot2::.pt
    ) +
    geom_text(
        aes(y = midpoint - 0.055, label = sprintf("%.1f%%", 100 * fraction)),
        family = font_family, size = 6.5 / ggplot2::.pt
    ) +
    annotate("segment",
        x = 0.88, xend = 2.12, y = 1.065, yend = 1.065,
        linewidth = 0.35
    ) +
    annotate("segment",
        x = 0.88, xend = 0.88, y = 1.04, yend = 1.065,
        linewidth = 0.35
    ) +
    annotate("segment",
        x = 2.12, xend = 2.12, y = 1.04, yend = 1.065,
        linewidth = 0.35
    ) +
    annotate(
        "text",
        x = 1.5, y = 1.12,
        label = sprintf("'One-sided Fisher exact'~italic(P) == %s", p_label),
        parse = TRUE, family = font_family, size = 6.0 / ggplot2::.pt
    ) +
    scale_fill_manual(values = c(Decay = red, Stable = blue), breaks = c("Decay", "Stable")) +
    scale_y_continuous(
        limits = c(0, 1.17), breaks = seq(0, 1, 0.25), labels = label_percent(accuracy = 1),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = "Fraction", fill = NULL) +
    theme_ac() +
    theme(
        legend.position = "top", legend.direction = "horizontal",
        legend.justification = "center", legend.margin = margin(0, 0, -2, 0),
        legend.key.size = unit(2.8, "mm"), legend.spacing.x = unit(1.2, "mm"),
        legend.text = element_text(size = 6.0),
        plot.margin = margin(0, 3, 1, 1, unit = "pt")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

half_stats <- decay[, .(
    mean = mean(half_life_min), SD = sd(half_life_min),
    lower = mean(half_life_min) - sd(half_life_min),
    upper = mean(half_life_min) + sd(half_life_min)
), by = plot_group]
half_stats[, hjust := fifelse(plot_group == "Edited", 0.58, 0.42)]
half_stats[, annotation := sprintf("%.1f \u00B1 %.1f min", mean, SD)]

pC <- ggplot(decay, aes(plot_group, half_life_min, fill = plot_group)) +
    geom_violin(trim = FALSE, linewidth = 0.45, colour = "black") +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.4, alpha = 0.35) +
    geom_errorbar(
        data = half_stats,
        aes(x = plot_group, ymin = pmax(lower, 5), ymax = upper),
        inherit.aes = FALSE, width = 0.09, linewidth = 0.55
    ) +
    geom_point(
        data = half_stats, aes(x = plot_group, y = mean), inherit.aes = FALSE,
        shape = 23, size = 3.0, stroke = 0.55, fill = "white", colour = "black"
    ) +
    geom_text(
        data = half_stats, aes(x = plot_group, y = 26000, label = annotation, hjust = hjust),
        inherit.aes = FALSE, family = font_family, size = 7.0 / ggplot2::.pt
    ) +
    scale_fill_manual(values = c(Edited = red, Unedited = blue)) +
    scale_y_log10(
        limits = c(5, 50000), breaks = c(10, 100, 1000, 10000),
        labels = label_number(big.mark = ",")
    ) +
    labs(x = NULL, y = "Half-life (min, log scale)") +
    theme_ac()
pC <- add_italic_n_axis(pC)

label_map <- c(
    edited = "Editing: Edited vs No detected",
    z_log2_transcript_length = "Transcript length (per 1 SD)",
    z_log2_published_expression = "Published expression (per 1 SD)",
    z_log2_CDS_length = "CDS length (per 1 SD)",
    z_log2_UTR5_length = "5'UTR length (per 1 SD)",
    z_log2_UTR3_length = "3'UTR length (per 1 SD)",
    z_GC_content = "GC content (per 1 SD)",
    z_A_content = "A content (per 1 SD)",
    decay_classIII = "Decay class III vs II",
    decay_classIV = "Decay class IV vs II",
    decay_classV = "Decay class V vs II"
)
top_to_bottom <- c(
    "Decay class IV vs II", "Editing: Edited vs No detected",
    "A content (per 1 SD)", "GC content (per 1 SD)",
    "5'UTR length (per 1 SD)", "3'UTR length (per 1 SD)",
    "CDS length (per 1 SD)", "Transcript length (per 1 SD)",
    "Published expression (per 1 SD)", "Decay class V vs II",
    "Decay class III vs II"
)
forest <- reg[
    model %chin% c("H3", "H4_component_length_model", "H5_extended_model") &
        variable != "(Intercept)" & is.finite(estimate)
]
forest[, label := unname(label_map[variable])]
forest[, label := factor(label, levels = rev(top_to_bottom))]
forest[, model_label := factor(
    model,
    levels = c("H3", "H4_component_length_model", "H5_extended_model"),
    labels = c("H3: total-length", "H4: component-length", "H5: sequence-extended")
)]
model_n <- unique(forest[, .(N_total, N_edited, N_unedited)])[1]

pD <- ggplot(
    forest,
    aes(estimate, label, xmin = CI_low95, xmax = CI_high95, colour = model_label)
) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
    geom_pointrange(
        position = position_dodge(width = 0.42), linewidth = 0.38,
        fatten = 1.6
    ) +
    scale_colour_manual(values = c(red, green, blue)) +
    scale_x_continuous(
        limits = c(-1.22, 1.30), breaks = seq(-1, 1, 0.5),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Regression coefficient beta", y = NULL, colour = "Model") +
    theme_minimal(base_family = font_family, base_size = 6.4) +
    theme(
        text = element_text(family = font_family, colour = "black"),
        panel.grid.major = element_line(colour = light_grid, linewidth = 0.32),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", linewidth = 0.45),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.35),
        axis.ticks.length = unit(1.5, "mm"),
        axis.text = element_text(size = 6.8, colour = "black"),
        axis.title.x = element_text(size = 7.2, margin = margin(t = 4)),
        legend.position = "none",
        plot.margin = margin(14, 2, 1, 0, unit = "pt")
    ) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE))

burden <- decay[
    plot_group == "Edited" & is.finite(gene_editing_burden) & gene_editing_burden > 0 & is.finite(log2_half_life)
]
burden[, log2_gene_editing_burden := log2(gene_editing_burden)]
burden_test <- cor.test(
    burden$gene_editing_burden, burden$log2_half_life,
    method = "spearman", exact = FALSE
)
burden_rho <- unname(burden_test$estimate)
burden_p <- burden_test$p.value
x_limits <- c(-7.5, -1.0)
x_breaks <- -7:-1
burden_n <- nrow(burden)
burden_rho_label <- sprintf("%.3f", burden_rho)
burden_p_label <- sprintf("%.3g", burden_p)
burden_title_prefix <- sprintf("Edited genes (>= %d CDS site%s),", min_sites, if (min_sites == 1L) "" else "s")
burden_title <- bquote(.(burden_title_prefix) ~ italic(n) == .(burden_n))
burden_subtitle <- bquote("Spearman's rho =" ~ .(burden_rho_label) * "," ~ italic(P) == .(burden_p_label))

pE <- ggplot(burden, aes(log2_gene_editing_burden, log2_half_life)) +
    geom_smooth(
        method = "lm", formula = y ~ x, se = TRUE,
        colour = "#A63D40", fill = "#BDBDBD", linewidth = 0.35, alpha = 0.45
    ) +
    geom_point(colour = red, size = 0.9) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = 4:9) +
    coord_cartesian(xlim = x_limits, ylim = c(3.8, 9.2)) +
    labs(
        x = expression(log[2] ~ "pooled editing allele fraction"), y = "log2 half-life",
        title = burden_title,
        subtitle = burden_subtitle
    ) +
    theme_ac() +
    theme(
        plot.title = element_text(size = 7.0, face = "plain", hjust = 0, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 7.0, face = "plain", hjust = 0, margin = margin(b = -1)),
        plot.margin = margin(0, 2, 3, 1, unit = "pt")
    )

# Match the Illustrator reference page exactly: 165 x 200 mm.  The figure
# occupies the upper half of the page; the lower half remains intentionally
# blank because this is how the supplied manuscript PDF is laid out.
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
raw_pdf <- paste0(output_pdf, ".raw.pdf")
cairo_pdf(
    raw_pdf,
    width = 165 / 25.4, height = 200 / 25.4,
    family = font_family, onefile = FALSE
)
grid.newpage()

draw_plot <- function(plot, x, y, width, height) {
    print(plot, vp = viewport(x = x, y = y, width = width, height = height))
}

draw_plot(pA, x = 0.180, y = 0.864, width = 0.350, height = 0.203)
draw_plot(pB, x = 0.503, y = 0.868, width = 0.292, height = 0.211)
draw_plot(pC, x = 0.802, y = 0.864, width = 0.350, height = 0.203)
draw_plot(pD, x = 0.320, y = 0.6085, width = 0.580, height = 0.241)
draw_plot(pE, x = 0.822, y = 0.614, width = 0.310, height = 0.248)

# Panel letters use the same size and weight as the supplied PDF.
for (item in list(
    list("A", 0.028, 0.979), list("B", 0.368, 0.979),
    list("C", 0.654, 0.979), list("D", 0.055, 0.740),
    list("E", 0.663, 0.740)
)) {
    grid.text(
        item[[1]],
        x = item[[2]], y = item[[3]], just = c("left", "top"),
        gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 10)
    )
}

# Add the italic n sample-size lines below panels A and C.  They are separate
# grobs so only n is italic, exactly matching the reference typography.
grid.text(
    bquote(italic(n) == .(n_length_edited_label)),
    x = 0.160, y = 0.759,
    gp = gpar(fontfamily = font_family, fontsize = 6.8)
)
grid.text(
    bquote(italic(n) == .(n_length_unedited_label)),
    x = 0.278, y = 0.759,
    gp = gpar(fontfamily = font_family, fontsize = 6.8)
)
grid.text(
    bquote(italic(n) == .(n_edited_label)),
    x = 0.783, y = 0.759,
    gp = gpar(fontfamily = font_family, fontsize = 6.8)
)
grid.text(
    bquote(italic(n) == .(n_unedited_label)),
    x = 0.901, y = 0.759,
    gp = gpar(fontfamily = font_family, fontsize = 6.8)
)

grid.text(
    "Model",
    x = 0.080, y = 0.731, just = "left",
    gp = gpar(fontfamily = font_family, fontsize = 6.8)
)
legend_items <- list(
    list(red, 0.137, 0.148, "H3: total-length"),
    list(green, 0.270, 0.281, "H4: component-length"),
    list(blue, 0.450, 0.461, "H5: sequence-extended")
)
for (item in legend_items) {
    grid.segments(
        x0 = item[[2]] - 0.006, x1 = item[[2]] + 0.006,
        y0 = 0.731, y1 = 0.731,
        gp = gpar(col = item[[1]], lwd = 1.5)
    )
    grid.points(
        x = item[[2]], y = 0.731, pch = 16, size = unit(1.8, "mm"),
        gp = gpar(col = item[[1]])
    )
    grid.text(
        item[[4]],
        x = item[[3]], y = 0.731, just = "left",
        gp = gpar(fontfamily = font_family, fontsize = 6.8)
    )
}

# Model N line is placed independently so the D legend and the plot panel keep
# the same vertical alignment as the supplied figure.
grid.text(
    sprintf(
        "N = %s (Edited = %s; control = %s)",
        comma(model_n$N_total, accuracy = 1), comma(model_n$N_edited, accuracy = 1),
        comma(model_n$N_unedited, accuracy = 1)
    ),
    x = 0.345, y = 0.715,
    gp = gpar(fontfamily = font_family, fontsize = 7.2)
)

dev.off()

gs_status <- system2(
    "gs",
    args = c(
        "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
        "-dFIXEDMEDIA", "-dPDFFitPage",
        "-dDEVICEWIDTHPOINTS=467.717", "-dDEVICEHEIGHTPOINTS=566.929",
        paste0("-sOutputFile=", shQuote(output_pdf)), shQuote(raw_pdf)
    )
)
if (!identical(gs_status, 0L)) stop("Ghostscript failed while setting the exact page size")
unlink(raw_pdf)

cat("Wrote:", normalizePath(output_pdf, mustWork = FALSE), "\n")
cat(sprintf(
    "Panel samples: decay Edited=%d, Unedited=%d; all Edited=%d, control=%d\n",
    n_edited, n_unedited, fisher$N_edited, fisher$N_control
))
cat(sprintf("Panel E: n=%d, Spearman rho=%.6f, P=%.6g\n", nrow(burden), burden_rho, burden_p))

cat(sprintf("Panel A: Edited=%d, Unedited=%d; classes I-V, length available, irrespective of half-life\n", n_length_edited, n_length_unedited))
