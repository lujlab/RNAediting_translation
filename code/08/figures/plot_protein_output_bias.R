#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(2, 3)) {
    stop("Usage:\n", "  Rscript plot_protein_output_bias.R outdir input.tsv [ratio|difference]\n\n",
        "Example:\n", "  Rscript plot_protein_output_bias.R ./Lp_Lm_violin_by_sex \\\n",
        "    dmel_Lp_union.combined_sex_Lp_wide.min3_LmLr.tsv\n")
}

outdir <- args[1]
input_file <- args[2]
metric <- if (length(args) == 3) tolower(args[3]) else "ratio"
if (!metric %in% c("ratio", "difference")) {
    stop("Metric must be 'ratio' or 'difference'.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

metric_tag <- if (metric == "ratio") "log2_Lp_over_Lm" else "Lp_minus_Lm"
prefix <- file.path(outdir, if (metric == "ratio") "DIANN_Lp_vs_Lm" else "DIANN_Lp_minus_Lm")

# ============================================================
# Parameters
# ============================================================

group_levels <- c("AG", "WBG", "AS")
sex_levels <- c("M", "F")
sex_labels <- c(M = "Male", F = "Female")
min_altcount <- 1
min_dp <- 1
log_ratio_pseudocount <- 1e-04

fill_cols <- c(AG = "#E777A7", WBG = "#6AACDD", AS = "#F5A966")

# ============================================================
# Helper functions
# ============================================================

first_nonempty <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0)
        NA_character_ else x[1]
}

mean_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0)
        NA_real_ else mean(x)
}

format_p_plotmath <- function(p, digits = 2) {
    if (!is.finite(p))
        return("\"NA\"")
    if (p <= 0)
        return("0")

    exponent <- floor(log10(p))
    mantissa <- p/(10^exponent)
    mantissa_txt <- formatC(mantissa, format = "f", digits = digits)
    mantissa_txt <- sub("0+$", "", mantissa_txt)
    mantissa_txt <- sub("\\.$", "", mantissa_txt)

    paste0(mantissa_txt, "%*%10^{", exponent, "}")
}

parse_sample <- function(x) {
    x <- as.character(x)
    x <- basename(gsub("\\\\", "/", x))
    x <- sub("\\.[^.]+$", "", x)

    m <- regexec("^(AG|AS|WBG|WBg)_([FM][0-9]+)$", x, perl = TRUE)
    r <- regmatches(x, m)[[1]]

    if (length(r) < 3) {
        return(data.table(sample = x, group = NA_character_, rep = NA_character_, sex = NA_character_))
    }

    group <- toupper(r[2])
    rep <- r[3]

    data.table(sample = x, group = group, rep = rep, sex = substr(rep, 1, 1))
}

pairwise_paired_wilcox <- function(dt, sex_value, g1, g2) {
    d1 <- dt[sex == sex_value & as.character(group) == g1, .(site_id, value1 = plot_value)]
    d2 <- dt[sex == sex_value & as.character(group) == g2, .(site_id, value2 = plot_value)]

    paired <- merge(d1, d2, by = "site_id", all = FALSE)
    paired <- paired[is.finite(value1) & is.finite(value2)]
    n_shared <- nrow(paired)

    if (n_shared < 3) {
        p_value <- NA_real_
    } else if (all(paired$value1 == paired$value2)) {
        p_value <- 1
    } else {
        wt <- tryCatch(suppressWarnings(wilcox.test(paired$value1, paired$value2, paired = TRUE,
            alternative = "two.sided", exact = FALSE, correct = TRUE)), error = function(e) NULL)
        p_value <- if (is.null(wt))
            NA_real_ else as.numeric(wt$p.value)
    }

    data.table(sex = sex_value, group1 = g1, group2 = g2, n_shared_sites = n_shared, median_group1 = if (n_shared >
        0)
        median(paired$value1) else NA_real_, median_group2 = if (n_shared > 0)
        median(paired$value2) else NA_real_, median_paired_difference_group1_minus_group2 = if (n_shared > 0) {
        median(paired$value1 - paired$value2)
    } else {
        NA_real_
    }, p_value = p_value)
}

make_violin_plot_by_sex <- function(plot_dt, pairwise_dt, summary_dt, sex_value, common_range,
    output_stem) {
    sex_name <- unname(sex_labels[sex_value])
    d <- copy(plot_dt[sex == sex_value])
    pw <- copy(pairwise_dt[sex == sex_value])
    sm <- copy(summary_dt[sex == sex_value])

    if (nrow(d) == 0) {
        return(NULL)
    }

    plot_lower <- min(common_range[1], 0)
    plot_upper_base <- max(common_range[2], 0)
    y_span <- plot_upper_base - plot_lower
    if (!is.finite(y_span) || y_span == 0)
        y_span <- 1

    if (nrow(pw) > 0) {
        pw[, y := plot_upper_base + (0.1 + 0.12 * comparison_order) * y_span]
        y_upper <- max(pw$y, na.rm = TRUE) + 0.13 * y_span
    } else {
        y_upper <- plot_upper_base + 0.2 * y_span
    }

    count_dt <- sm[, .(group, x = as.numeric(group), y = plot_lower, label = paste0("italic(n)==",
        n_sites))]

    p <- ggplot(d, aes(x = group_num, y = plot_value, group = group, fill = group)) + geom_hline(yintercept = 0,
        linetype = 2, linewidth = 0.45, color = "grey40") + geom_violin(width = 0.82, scale = "width",
        trim = TRUE, alpha = 0.95, color = NA, linewidth = 0) + geom_boxplot(width = 0.16,
        outlier.shape = NA, fill = "#F7F7F7", alpha = 0.7, color = "grey35", linewidth = 0.45)

    if (nrow(pw) > 0) {
        p <- p + geom_segment(data = pw, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE,
            linewidth = 0.45, color = "black") + geom_segment(data = pw, aes(x = x1, xend = x1,
            y = y - 0.02 * y_span, yend = y), inherit.aes = FALSE, linewidth = 0.45, color = "black") +
            geom_segment(data = pw, aes(x = x2, xend = x2, y = y - 0.02 * y_span, yend = y),
                inherit.aes = FALSE, linewidth = 0.45, color = "black") + geom_text(data = pw,
            aes(x = xmid, y = y + 0.015 * y_span, label = p_label), inherit.aes = FALSE,
            hjust = 0.5, vjust = 0, size = 3, color = "black", lineheight = 0.95, parse = TRUE)
    }

    p <- p + geom_text(data = count_dt, aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.55, size = 3.1, color = "black", parse = TRUE) + scale_fill_manual(values = fill_cols,
        drop = FALSE) + scale_x_continuous(breaks = seq_along(group_levels), labels = expression(italic(Adar)^G,
        plain(WBg), italic(Adar)^S), limits = c(0.55, 3.45)) + scale_y_continuous(expand = expansion(mult = c(0.01,
        0.01))) + coord_cartesian(ylim = c(plot_lower, y_upper), clip = "off") + labs(x = NULL,
        y = if (metric == "ratio") {
            expression(log[2]((L[p] + 10^{
                -4
            })/(L[m] + 10^{
                -4
            })))
        } else {
            expression(L[p] - L[m])
        }, title = sex_name) + theme_classic(base_size = 10) + theme(axis.title.y = element_text(size = 10,
        margin = margin(r = 7)), axis.text.x = element_text(size = 9, color = "black", margin = margin(t = 4)),
        axis.text.y = element_text(size = 9, color = "black"), axis.line = element_line(color = "black",
            linewidth = 0.55), axis.ticks = element_line(color = "black", linewidth = 0.45),
        axis.ticks.length = unit(0.14, "cm"), legend.position = "none", plot.title = element_text(size = 10,
            face = "bold", hjust = 0.5), plot.margin = margin(12, 12, 20, 12))

    ggsave(paste0(output_stem, ".", sex_name, ".", metric_tag, "_violin_paired_wilcox.SX31style.pdf"),
        p, width = 4, height = 5.7)
    ggsave(paste0(output_stem, ".", sex_name, ".", metric_tag, "_violin_paired_wilcox.SX31style.png"),
        p, width = 4, height = 5.7, dpi = 300)
    p
}

save_combined_plot <- function(p_male, p_female, output_stem) {
    if (is.null(p_male) || is.null(p_female))
        return(invisible(NULL))

    # A4 portrait canvas.
    #
    # Male/Female are still arranged left-to-right in ONE row, but the
    # complete plotting row occupies only ~1/4 of the A4 page height.
    # This avoids vertically stretching the two panels to fill the page.
    pdf(paste0(output_stem, ".Male_Female.combined.A4_portrait.SX31style.pdf"), width = 8.27,
        height = 11.69)
    grid.newpage()

    pushViewport(viewport(x = 0.5, y = 0.5, width = unit(0.98, "npc"), height = unit(0.25,
        "npc"), just = c("center", "center"), layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1,
        1), "null")), clip = "off"))

    print(p_male, vp = viewport(layout.pos.row = 1, layout.pos.col = 1, clip = "off"))

    print(p_female, vp = viewport(layout.pos.row = 1, layout.pos.col = 2, clip = "off"))

    popViewport()
    dev.off()

    png(paste0(output_stem, ".Male_Female.combined.A4_portrait.SX31style.png"), width = 2481,
        height = 3507, res = 300)
    grid.newpage()

    pushViewport(viewport(x = 0.5, y = 0.5, width = unit(0.98, "npc"), height = unit(0.25,
        "npc"), just = c("center", "center"), layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1,
        1), "null")), clip = "off"))

    print(p_male, vp = viewport(layout.pos.row = 1, layout.pos.col = 1, clip = "off"))

    print(p_female, vp = viewport(layout.pos.row = 1, layout.pos.col = 2, clip = "off"))

    popViewport()
    dev.off()
}

# ============================================================
# Read and standardize input
# ============================================================
#
# Current input format:
#   site
#   genotype
#   Lp_mean_F
#   Lp_mean_M
#   mean_Lm_F
#   mean_Lm_M
#
# Only AG / WBg / AS are retained.
# The current table is already summarized at genotype x sex level, so no
# sample parsing, altcount/dp filtering, or extra replicate averaging is done.
# ============================================================

if (!file.exists(input_file)) {
    stop("Missing input file: ", input_file)
}

dt <- fread(input_file, sep = "\t", header = TRUE, quote = "", fill = TRUE, check.names = FALSE,
    na.strings = c("NA", "NaN", ""))

if (!"site_id" %in% names(dt) && "site" %in% names(dt)) {
    setnames(dt, "site", "site_id")
}

if (!"genotype" %in% names(dt) && "group" %in% names(dt)) {
    dt[, genotype := as.character(group)]
}

required <- c("site_id", "genotype", "Lp_mean_F", "Lp_mean_M", "mean_Lm_F", "mean_Lm_M")

missing_required <- setdiff(required, names(dt))
if (length(missing_required) > 0) {
    stop("Input file is missing required columns: ", paste(missing_required, collapse = ", "))
}

optional_character <- c("gene", "FBgn", "FBtr", "refAA", "known_alt_AAs")

for (nm in optional_character) {
    if (!nm %in% names(dt)) {
        dt[, (nm) := NA_character_]
    }
}

for (nm in c("Lp_mean_F", "Lp_mean_M", "mean_Lm_F", "mean_Lm_M")) {
    dt[, (nm) := suppressWarnings(as.numeric(get(nm)))]
}

dt[, site_id := trimws(as.character(site_id))]
dt[, genotype := trimws(as.character(genotype))]
dt[, group := toupper(genotype)]

# Only AG / WBg / AS.
dt <- dt[group %in% group_levels]

if (nrow(dt) == 0) {
    stop("No AG/WBg/AS rows remain in the current input table.")
}

male_dt <- dt[, .(site_id, group, sex = "M", gene, FBgn, FBtr, refAA, altAA = known_alt_AAs,
    Lp_use = Lp_mean_M, Lm_use = mean_Lm_M)]

female_dt <- dt[, .(site_id, group, sex = "F", gene, FBgn, FBtr, refAA, altAA = known_alt_AAs,
    Lp_use = Lp_mean_F, Lm_use = mean_Lm_F)]

long_dt <- rbindlist(list(male_dt, female_dt), use.names = TRUE, fill = TRUE)

# Keep true zero values because the plotted ratio contains the 1e-4 pseudocount.
long_dt[, passes_filter := !is.na(group) & group %in% group_levels & !is.na(sex) & sex %in%
    c("F", "M") & !is.na(site_id) & site_id != "" & site_id != "NA" & is.finite(Lp_use) &
    Lp_use >= 0 & Lp_use <= 1 & is.finite(Lm_use) & Lm_use >= 0 & Lm_use <= 1]

filter_qc <- long_dt[, .(n_input_rows = .N, n_unique_sites_input = uniqueN(site_id), n_rows_Lp_valid_0_to_1 = sum(is.finite(Lp_use) &
    Lp_use >= 0 & Lp_use <= 1, na.rm = TRUE), n_rows_Lp_zero = sum(is.finite(Lp_use) & Lp_use ==
    0, na.rm = TRUE), n_unique_sites_with_Lp_zero = uniqueN(site_id[is.finite(Lp_use) & Lp_use ==
    0]), n_rows_Lm_valid_0_to_1 = sum(is.finite(Lm_use) & Lm_use >= 0 & Lm_use <= 1, na.rm = TRUE),
    n_rows_Lm_zero = sum(is.finite(Lm_use) & Lm_use == 0, na.rm = TRUE), n_rows_pass = sum(passes_filter,
        na.rm = TRUE), n_unique_sites_pass = uniqueN(site_id[passes_filter])), by = .(group,
    sex)][order(group, sex)]

fwrite(filter_qc, paste0(prefix, ".filter_QC_by_sample.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

fwrite(long_dt, paste0(prefix, ".all_rows_with_filter_status.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

work <- long_dt[passes_filter == TRUE]

if (nrow(work) == 0) {
    stop("No rows remain after requiring finite 0 <= Lp <= 1 ", "and finite 0 <= Lm <= 1.")
}

plot_dt <- work[, .(gene = first_nonempty(gene), FBgn = first_nonempty(FBgn), FBtr = first_nonempty(FBtr),
    refAA = first_nonempty(refAA), altAA = first_nonempty(altAA), Lm_mean = mean_finite(Lm_use),
    Lp_mean = mean_finite(Lp_use), n_input_rows = .N), by = .(sex, group, site_id)]

plot_dt[, log2_Lp_over_Lm := log2((Lp_mean + log_ratio_pseudocount)/(Lm_mean + log_ratio_pseudocount))]

plot_dt[, Lp_minus_Lm := Lp_mean - Lm_mean]
plot_dt[, plot_value := if (metric == "ratio") log2_Lp_over_Lm else Lp_minus_Lm]

plot_dt[, ratio_of_group_means := fifelse(is.finite(Lm_mean) & Lm_mean > 0, Lp_mean/Lm_mean,
    NA_real_)]

plot_dt[, `:=`(n_samples = n_input_rows, samples_used = "pre_averaged_current_wide_table")]

plot_dt <- plot_dt[is.finite(plot_value)]
plot_dt[, sex := factor(sex, levels = sex_levels)]
plot_dt[, group := factor(group, levels = group_levels)]
plot_dt[, group_num := as.numeric(group)]

if (nrow(plot_dt) == 0) {
    stop("No finite site-level values remain for metric: ", metric)
}

# Keep the original auxiliary output filenames.
fwrite(work, paste0(prefix, ".sample_level_", metric_tag, "_by_sex.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

fwrite(plot_dt, paste0(prefix, ".site_sex_group_mean_", metric_tag, ".tsv"), sep = "\t",
    quote = FALSE, na = "NA")

comparison_pairs <- list(c("AG", "WBG"), c("WBG", "AS"), c("AG", "AS"))
sex_present <- sex_levels[sex_levels %in% as.character(unique(plot_dt$sex))]

pairwise_dt <- rbindlist(lapply(sex_present, function(sx) {
    rbindlist(lapply(comparison_pairs, function(z) pairwise_paired_wilcox(plot_dt, sx, z[1],
        z[2])), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

pairwise_dt[, p_adjust_BH := p.adjust(p_value, method = "BH"), by = sex]
pairwise_dt[, comparison := paste(group1, group2, sep = " vs ")]
pairwise_dt[, p_label := paste0("italic(P)==", vapply(p_adjust_BH, format_p_plotmath,
    character(1)), "*','~~italic(n)==", n_shared_sites)]
pairwise_dt[, x1 := match(group1, group_levels)]
pairwise_dt[, x2 := match(group2, group_levels)]
pairwise_dt[, xmid := (x1 + x2)/2]
pairwise_dt[, comparison_order := seq_len(.N), by = sex]

fwrite(pairwise_dt, paste0(prefix, ".", metric_tag, ".pairwise_paired_wilcox_BH_by_sex.tsv"),
    sep = "\t", quote = FALSE, na = "NA")

summary_dt <- plot_dt[, .(n_sites = uniqueN(site_id), mean_value = mean(plot_value), median_value = median(plot_value),
    q25_value = quantile(plot_value, 0.25), q75_value = quantile(plot_value, 0.75), n_below_reference = sum(plot_value <
        0), n_equal_reference = sum(plot_value == 0), n_above_reference = sum(plot_value >
        0)), by = .(sex, group)][order(sex, group)]

fwrite(summary_dt, paste0(prefix, ".", metric_tag, ".group_summary_by_sex.tsv"), sep = "\t",
    quote = FALSE, na = "NA")

common_range <- range(plot_dt$plot_value, na.rm = TRUE)
p_male <- if ("M" %in% sex_present) make_violin_plot_by_sex(plot_dt, pairwise_dt, summary_dt,
    "M", common_range, prefix) else NULL
p_female <- if ("F" %in% sex_present) make_violin_plot_by_sex(plot_dt, pairwise_dt, summary_dt,
    "F", common_range, prefix) else NULL
save_combined_plot(p_male, p_female, prefix)

cat("\nDone.\n")
cat("Input:\n", input_file, "\n\n", sep = "")
cat("Input format:\n")
cat("  current genotype/sex wide table\n")
cat("  only AG, WBg and AS are retained\n")
cat("  Female uses Lp_mean_F and mean_Lm_F\n")
cat("  Male uses Lp_mean_M and mean_Lm_M\n\n")
cat("Filter:\n")
cat("  finite 0 <= Lp <= 1; Lp = 0 is retained exactly\n")
cat("  finite 0 <= Lm <= 1; Lm = 0 is retained exactly\n\n")
cat("Metric:\n")
if (metric == "ratio") {
    cat("  log2((Lp+1e-4)/(Lm+1e-4)): uses log2((Lp + ", log_ratio_pseudocount, ") / (Lm + ",
        log_ratio_pseudocount, "))\n", sep = "")
} else {
    cat("  Lp-Lm: direct difference of the current sex-level means; no pseudocount\n")
}
cat("  Lp and Lm are read directly from the current sex-level mean columns\n\n")
cat("Pairwise tests:\n")
cat("  paired Wilcoxon signed-rank tests on shared site IDs within the same sex\n")
cat("  BH correction is applied separately within each sex (3 comparisons)\n\n")
cat("Main outputs:\n")
cat(paste0(prefix, ".site_sex_group_mean_", metric_tag, ".tsv\n"))
cat(paste0(prefix, ".", metric_tag, ".pairwise_paired_wilcox_BH_by_sex.tsv\n"))
cat(paste0(prefix, ".Male.", metric_tag, "_violin_paired_wilcox.SX31style.pdf\n"))
cat(paste0(prefix, ".Female.", metric_tag, "_violin_paired_wilcox.SX31style.pdf\n"))
cat(paste0(prefix, ".Male_Female.combined.A4_portrait.SX31style.pdf\n"))
