#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
summary_dir <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else getwd()

pool_summary_file <- file.path(summary_dir, "AA_three_group_four_species_pool_summary.tsv")
codon_file <- file.path(summary_dir, "AA_four_group_species_sex_codon_long.tsv.gz")
unit_mean_file <- file.path(summary_dir, "AA_three_group_main8_species_sex_AA_group_means.tsv")

for (input_file in c(pool_summary_file, codon_file)) {
    if (!file.exists(input_file))
        stop("Missing input: ", input_file)
}

pdf_file <- file.path(summary_dir, "AA_three_group_mean_heatmaps_two_way_and_species_sex_adjusted_A4.pdf")
legacy_pdf_file <- file.path(summary_dir, "AA_three_group_pool_mean_median_heatmaps_only_blue0_white1.5_red3_wilcox_A4.pdf")
alias_pdf_file <- file.path(summary_dir, "AA_three_group_pool_mean_median_heatmaps_blue0_white1_red2_ANOVA_A4_portrait.pdf")

group_order <- c("Unedited_preferred", "Unedited_unpreferred", "Edited")
group_labels <- c(Unedited_preferred = "Unedited\npreferred", Unedited_unpreferred = "Unedited\nunpreferred",
    Edited = "Edited")
amino_acid_order <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
font_family <- "Arial"

significance_code <- function(p) {
    fcase(is.na(p), "NA", p < 1e-04, "****", p < 0.001, "***", p < 0.01, "**", p < 0.05,
        "*", default = "ns")
}

format_axis_p <- function(p) {
    if (!is.finite(p))
        return("NA")
    if (p < 1e-300)
        return("<1e-300")
    if (p < 0.001)
        return(formatC(p, format = "e", digits = 1))
    formatC(p, format = "f", digits = 3L)
}

displayed_one_decimal <- function(x) {
    suppressWarnings(as.numeric(sprintf("%.1f", x)))
}

quadratic_f_test <- function(delta, covariance, df_denominator) {
    rank_c <- qr(covariance)$rank
    inverse <- tryCatch(solve(covariance), error = function(e) NULL)
    if (is.null(inverse) || rank_c < 1L || !is.finite(df_denominator)) {
        return(c(F_statistic = NA_real_, df1 = NA_real_, p_value = NA_real_))
    }
    f_value <- drop(t(delta) %*% inverse %*% delta)/rank_c
    c(F_statistic = f_value, df1 = rank_c, p_value = pf(f_value, rank_c, df_denominator,
        lower.tail = FALSE))
}

# -----------------------------------------------------------------------------
# Panel 1: codon-pooled two-way model, gmRC ~ AA * group.
# The model is reconstructed exactly from cell n/mean/SD. The full interaction
# model has one fitted mean per AA x group cell, so its residual MSE is the
# pooled within-cell variance from the original codon observations.
# -----------------------------------------------------------------------------
pool <- fread(pool_summary_file)
pool[, `:=`(AA = gsub("^\"|\"$", "", as.character(AA)), group = gsub("^\"|\"$", "", as.character(group)),
    n_codons = as.numeric(n_codons), mean_gmRC = as.numeric(mean_gmRC), sd_gmRC = as.numeric(sd_gmRC))]
pool <- pool[AA %chin% amino_acid_order & group %chin% group_order & n_codons >= 2 & is.finite(mean_gmRC) &
    is.finite(sd_gmRC)]
present_aa <- amino_acid_order[amino_acid_order %chin% unique(pool$AA)]
pool <- pool[AA %chin% present_aa]
if (pool[, uniqueN(paste(AA, group))] != length(present_aa) * length(group_order)) {
    stop("Pooled summary does not contain a complete AA x group table")
}

pool[, `:=`(AA = factor(AA, levels = present_aa), group = factor(group, levels = group_order))]

ss_within_pool <- pool[, sum((n_codons - 1) * sd_gmRC^2)]
df_within_pool <- pool[, sum(n_codons) - .N]
mse_pool <- ss_within_pool/df_within_pool

weighted_between_rss <- function(formula_object) {
    fit <- lm(formula_object, data = pool, weights = n_codons)
    sum(pool$n_codons * residuals(fit)^2)
}

rss_additive <- weighted_between_rss(mean_gmRC ~ AA + group)
rss_group_only <- weighted_between_rss(mean_gmRC ~ group)
rss_aa_only <- weighted_between_rss(mean_gmRC ~ AA)
ss_aa <- max(0, rss_group_only - rss_additive)
ss_group <- max(0, rss_aa_only - rss_additive)
ss_interaction <- max(0, rss_additive)

pool_global <- data.table(term = c("AA (adjusted for group)", "group (adjusted for AA)",
    "AA:group"), df1 = c(length(present_aa) - 1L, length(group_order) - 1L, (length(present_aa) -
    1L) * (length(group_order) - 1L)), sum_squares = c(ss_aa, ss_group, ss_interaction),
    residual_df = df_within_pool, residual_MSE = mse_pool, method = "Two-way ANOVA reconstructed from codon-level cell n/mean/SD")
pool_global[, F_statistic := (sum_squares/df1)/residual_MSE]
pool_global[, p_value := pf(F_statistic, df1, residual_df, lower.tail = FALSE)]
fwrite(pool_global, file.path(summary_dir, "AA_three_group_codon_pooled_two_way_ANOVA.tsv"),
    sep = "\t", na = "NA")

run_pool_simple_effect <- function(aa_name) {
    d <- copy(pool[as.character(AA) == aa_name])
    d[, group_order_value := match(as.character(group), group_order)]
    setorder(d, group_order_value)
    d[, group_order_value := NULL]
    mu <- d$mean_gmRC
    n <- d$n_codons
    contrast_matrix <- rbind(`Edited - Unedited_preferred` = c(-1, 0, 1), `Edited - Unedited_unpreferred` = c(0,
        -1, 1))
    delta <- drop(contrast_matrix %*% mu)
    covariance_mu <- diag(mse_pool/n, nrow = 3L)
    covariance_delta <- contrast_matrix %*% covariance_mu %*% t(contrast_matrix)
    test <- quadratic_f_test(delta, covariance_delta, df_within_pool)
    data.table(AA = aa_name, test = "Simple group effect within AA from two-way AA*group model",
        df1 = unname(test[["df1"]]), df2 = df_within_pool, F_statistic = unname(test[["F_statistic"]]),
        p_value = unname(test[["p_value"]]))
}

pool_simple <- rbindlist(lapply(present_aa, run_pool_simple_effect))
pool_simple[, p_adj_BH := p.adjust(p_value, method = "BH")]
pool_simple[, significance := significance_code(p_adj_BH)]

run_pool_posthoc <- function(aa_name) {
    d <- copy(pool[as.character(AA) == aa_name])
    edited <- d[as.character(group) == "Edited"]
    rbindlist(lapply(group_order[1:2], function(control_name) {
        control <- d[as.character(group) == control_name]
        delta <- edited$mean_gmRC - control$mean_gmRC
        se_tukey <- sqrt(mse_pool/2 * (1/edited$n_codons + 1/control$n_codons))
        q_value <- abs(delta)/se_tukey
        p_tukey <- ptukey(q_value, nmeans = length(group_order), df = df_within_pool, lower.tail = FALSE)
        edited_display <- displayed_one_decimal(edited$mean_gmRC)
        control_display <- displayed_one_decimal(control$mean_gmRC)
        data.table(AA = aa_name, contrast = paste("Edited", control_name, sep = " - "), edited_mean = edited$mean_gmRC,
            control_mean = control$mean_gmRC, edited_heatmap_value = edited_display, control_heatmap_value = control_display,
            estimate_difference = delta, standard_error_TukeyKramer = se_tukey, q_statistic = q_value,
            p_adj_TukeyKramer = p_tukey, edited_higher_on_heatmap = edited_display > control_display,
            edited_significantly_higher = edited_display > control_display & p_tukey < 0.05)
    }))
}

pool_posthoc <- rbindlist(lapply(present_aa, run_pool_posthoc))
pool_direction <- pool_posthoc[, .(edited_higher_than_any_unedited = any(edited_significantly_higher)),
    by = AA]
pool_simple <- merge(pool_simple, pool_direction, by = "AA", all.x = TRUE)
pool_simple[, aa_order := match(AA, present_aa)]
setorder(pool_simple, aa_order)
pool_simple[, aa_order := NULL]
pool_posthoc[, aa_order := match(AA, present_aa)]
setorder(pool_posthoc, aa_order, contrast)
pool_posthoc[, aa_order := NULL]
fwrite(pool_simple, file.path(summary_dir, "AA_three_group_codon_pooled_simple_group_effect_by_AA.tsv"),
    sep = "\t", na = "NA")
fwrite(pool_posthoc, file.path(summary_dir, "AA_three_group_codon_pooled_Tukey_within_AA.tsv"),
    sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# Panel 2: each species-sex unit contributes one AA x group mean. Model:
# unit_mean_gmRC ~ AA * group + species + sex.
# -----------------------------------------------------------------------------
if (!file.exists(unit_mean_file)) {
    message("Creating species-sex AA x group means from the codon-level table...")
    extract_command <- sprintf("gzip -dc %s | cut -f1,2,6,8,9,30,31,44", shQuote(codon_file))
    unit_raw <- fread(cmd = extract_command, showProgress = TRUE)
    unit_raw <- unit_raw[analysis_role == "main_8" & valid_codon == TRUE & AA %chin% present_aa &
        group %chin% group_order & is.finite(gmRC) & gmRC >= 0]
    unit_means <- unit_raw[, .(unit_mean_gmRC = mean(gmRC), n_codons = .N, unit_sd_gmRC = sd(gmRC)),
        by = .(unit_id, species, sex, AA, group)]
    rm(unit_raw)
    invisible(gc())
    fwrite(unit_means, unit_mean_file, sep = "\t", na = "NA")
} else {
    unit_means <- fread(unit_mean_file)
}

unit_means[, `:=`(unit_id = as.character(unit_id), species = as.character(species), sex = as.character(sex),
    AA = as.character(AA), group = as.character(group), unit_mean_gmRC = as.numeric(unit_mean_gmRC))]
unit_means <- unit_means[AA %chin% present_aa & group %chin% group_order & is.finite(unit_mean_gmRC)]
unit_means[, `:=`(AA = factor(AA, levels = present_aa), group = factor(group, levels = group_order),
    species = factor(species), sex = factor(sex))]

unit_fit <- lm(unit_mean_gmRC ~ AA * group + species + sex, data = unit_means)
if (anyNA(coef(unit_fit))) {
    stop("Species-sex model is rank-deficient; adjusted means are not unique")
}

unit_anova <- as.data.table(anova(unit_fit), keep.rownames = "term")
setnames(unit_anova, old = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"), new = c("df",
    "sum_squares", "mean_squares", "F_statistic", "p_value"), skip_absent = TRUE)
unit_anova[, method := "lm: unit_mean_gmRC ~ AA*group + species + sex"]
fwrite(unit_anova, file.path(summary_dir, "AA_three_group_species_sex_unit_mean_model_ANOVA.tsv"),
    sep = "\t", na = "NA")

unit_grid_base <- CJ(species = levels(unit_means$species), sex = levels(unit_means$sex),
    unique = TRUE)
model_terms <- delete.response(terms(unit_fit))
unit_beta <- coef(unit_fit)
unit_vcov <- vcov(unit_fit)

get_unit_l_vector <- function(aa_name, group_name) {
    newdata <- copy(unit_grid_base)
    newdata[, `:=`(AA = factor(aa_name, levels = present_aa), group = factor(group_name,
        levels = group_order), species = factor(species, levels = levels(unit_means$species)),
        sex = factor(sex, levels = levels(unit_means$sex)))]
    x_matrix <- model.matrix(model_terms, data = newdata)
    colMeans(x_matrix)
}

unit_l_vectors <- list()
unit_emmeans <- rbindlist(lapply(present_aa, function(aa_name) {
    rbindlist(lapply(group_order, function(group_name) {
        key <- paste(aa_name, group_name, sep = "||")
        l_vector <- get_unit_l_vector(aa_name, group_name)
        unit_l_vectors[[key]] <<- l_vector
        estimate <- drop(l_vector %*% unit_beta)
        standard_error <- sqrt(drop(l_vector %*% unit_vcov %*% l_vector))
        data.table(AA = aa_name, group = group_name, adjusted_mean_gmRC = estimate, standard_error = standard_error,
            df = df.residual(unit_fit), n_species_sex_cells = unit_means[as.character(AA) ==
                aa_name & as.character(group) == group_name, .N])
    }))
}))

run_unit_posthoc <- function(aa_name) {
    edited_key <- paste(aa_name, "Edited", sep = "||")
    edited_l <- unit_l_vectors[[edited_key]]
    edited_mean <- unit_emmeans[AA == aa_name & group == "Edited", adjusted_mean_gmRC]
    ans <- rbindlist(lapply(group_order[1:2], function(control_name) {
        control_key <- paste(aa_name, control_name, sep = "||")
        control_l <- unit_l_vectors[[control_key]]
        contrast_l <- edited_l - control_l
        delta <- drop(contrast_l %*% unit_beta)
        se <- sqrt(drop(contrast_l %*% unit_vcov %*% contrast_l))
        t_value <- delta/se
        p_value <- 2 * pt(abs(t_value), df.residual(unit_fit), lower.tail = FALSE)
        control_mean <- unit_emmeans[AA == aa_name & group == control_name, adjusted_mean_gmRC]
        edited_display <- displayed_one_decimal(edited_mean)
        control_display <- displayed_one_decimal(control_mean)
        data.table(AA = aa_name, contrast = paste("Edited", control_name, sep = " - "), edited_adjusted_mean = edited_mean,
            control_adjusted_mean = control_mean, edited_heatmap_value = edited_display,
            control_heatmap_value = control_display, estimate_difference = delta, standard_error = se,
            t_statistic = t_value, df = df.residual(unit_fit), p_value = p_value, edited_higher_on_heatmap = edited_display >
                control_display)
    }))
    ans[, p_adj_BH_within_AA := p.adjust(p_value, method = "BH")]
    ans[, edited_significantly_higher := edited_higher_on_heatmap & p_adj_BH_within_AA <
        0.05]
    ans
}

unit_posthoc <- rbindlist(lapply(present_aa, run_unit_posthoc))

run_unit_simple_effect <- function(aa_name) {
    l_edited <- unit_l_vectors[[paste(aa_name, "Edited", sep = "||")]]
    l_pref <- unit_l_vectors[[paste(aa_name, "Unedited_preferred", sep = "||")]]
    l_unpref <- unit_l_vectors[[paste(aa_name, "Unedited_unpreferred", sep = "||")]]
    contrast_matrix <- rbind(l_edited - l_pref, l_edited - l_unpref)
    delta <- drop(contrast_matrix %*% unit_beta)
    covariance <- contrast_matrix %*% unit_vcov %*% t(contrast_matrix)
    test <- quadratic_f_test(delta, covariance, df.residual(unit_fit))
    data.table(AA = aa_name, test = "Simple group effect within AA from AA*group+species+sex model",
        df1 = unname(test[["df1"]]), df2 = df.residual(unit_fit), F_statistic = unname(test[["F_statistic"]]),
        p_value = unname(test[["p_value"]]))
}

unit_simple <- rbindlist(lapply(present_aa, run_unit_simple_effect))
unit_simple[, p_adj_BH := p.adjust(p_value, method = "BH")]
unit_simple[, significance := significance_code(p_adj_BH)]
unit_direction <- unit_posthoc[, .(edited_higher_than_any_unedited = any(edited_significantly_higher)),
    by = AA]
unit_simple <- merge(unit_simple, unit_direction, by = "AA", all.x = TRUE)
for (object_name in c("unit_emmeans", "unit_posthoc", "unit_simple")) {
    object_value <- get(object_name)
    object_value[, aa_order := match(AA, present_aa)]
    setorder(object_value, aa_order)
    object_value[, aa_order := NULL]
    assign(object_name, object_value)
}
fwrite(unit_emmeans, file.path(summary_dir, "AA_three_group_species_sex_adjusted_marginal_means.tsv"),
    sep = "\t", na = "NA")
fwrite(unit_simple, file.path(summary_dir, "AA_three_group_species_sex_simple_group_effect_by_AA.tsv"),
    sep = "\t", na = "NA")
fwrite(unit_posthoc, file.path(summary_dir, "AA_three_group_species_sex_Edited_posthoc_by_AA.tsv"),
    sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# Heatmaps and axis annotations.
# -----------------------------------------------------------------------------
make_axis_labels <- function(stats_table) {
    setNames(stats_table[, paste0(AA, "  P_BH", fifelse(p_adj_BH < 1e-300, "", "="), vapply(p_adj_BH,
        format_axis_p, character(1)), "  ", significance)], stats_table$AA)
}

pool_axis_labels <- make_axis_labels(pool_simple)
unit_axis_labels <- make_axis_labels(unit_simple)

pool_plot_dt <- pool[, .(AA = as.character(AA), group = as.character(group), value = mean_gmRC)]
unit_plot_dt <- unit_emmeans[, .(AA, group, value = adjusted_mean_gmRC)]

make_heatmap <- function(plot_dt, axis_labels, title_text, subtitle_text) {
    d <- copy(plot_dt)
    d[, `:=`(AA = factor(AA, levels = present_aa), group = factor(group, levels = group_order),
        color_value = pmin(pmax(value, 0), 2), value_label = fifelse(is.finite(value), sprintf("%.1f",
            value), "NA"))]
    ggplot(d, aes(x = AA, y = group, fill = color_value)) + geom_tile(color = "white", linewidth = 0.55) +
        geom_text(aes(label = value_label), family = font_family, size = 8/2.845276, color = "black") +
        scale_x_discrete(limits = present_aa, labels = axis_labels, drop = FALSE) + scale_y_discrete(limits = rev(group_order),
        labels = group_labels, drop = FALSE) + scale_fill_gradient2(low = "#2166AC", mid = "white",
        high = "#B2182B", midpoint = 1, limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2),
        labels = c("0", "0.5", "1.0", "1.5", "2.0"), name = "gmRC") + labs(title = title_text,
        subtitle = subtitle_text, x = NULL, y = NULL) + theme_minimal(base_size = 8, base_family = font_family) +
        theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 2)),
            plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 5)),
            axis.text.x = element_text(size = 8, angle = 270, hjust = 0, vjust = 0.5, color = "black"),
            axis.text.y = element_text(size = 8, color = "black"), panel.grid = element_blank(),
            legend.position = "right", legend.title = element_text(size = 8), legend.text = element_text(size = 8),
            legend.key.height = grid::unit(0.27, "in"), legend.key.width = grid::unit(0.11,
                "in"), plot.margin = margin(t = 4, r = 6, b = 3, l = 6))
}

style_direction_stars <- function(plot_object, stats_table) {
    g <- ggplotGrob(plot_object)
    axis_index <- which(g$layout$name == "axis-b")
    axis_grob <- g$grobs[[axis_index]]
    axis_table <- axis_grob$children[["axis"]]
    title_index <- which(vapply(axis_table$grobs, inherits, logical(1), what = "titleGrob"))
    title_grob <- axis_table$grobs[[title_index]]
    original_text <- title_grob$children[[1L]]
    base_labels <- stats_table[, paste0(AA, "  P_BH", fifelse(p_adj_BH < 1e-300, "", "="),
        vapply(p_adj_BH, format_axis_p, character(1)))]
    star_labels <- fifelse(stats_table$edited_higher_than_any_unedited == TRUE, stats_table$significance,
        "")
    star_gp <- original_text$gp
    star_gp$col <- "#8B0000"
    base_widths <- do.call(grid::unit.c, lapply(base_labels, function(label_i) {
        grid::grobWidth(grid::textGrob(label_i, gp = original_text$gp))
    }))
    star_text <- grid::textGrob(name = "direction_star_overlay", label = star_labels, x = original_text$x,
        y = original_text$y - base_widths - grid::unit(2, "pt"), rot = 270, hjust = 0, vjust = 0.5,
        gp = star_gp)
    title_grob <- grid::addGrob(title_grob, star_text)
    axis_table$grobs[[title_index]] <- title_grob
    axis_grob$children[["axis"]] <- axis_table
    g$grobs[[axis_index]] <- axis_grob
    g
}

p_pool <- make_heatmap(pool_plot_dt, pool_axis_labels, "Codon-pooled mean gmRC", "AA <U+00D7> group two-way model; BH-adjusted within-AA P")
p_unit <- make_heatmap(unit_plot_dt, unit_axis_labels, "Species<U+2013>sex adjusted marginal mean gmRC",
    "AA <U+00D7> group + species + sex; BH-adjusted within-AA P")
g_pool <- style_direction_stars(p_pool, pool_simple)
g_unit <- style_direction_stars(p_unit, unit_simple)

grDevices::cairo_pdf(pdf_file, width = 8.27, height = 11.69, family = font_family, pointsize = 8)
grid::grid.newpage()
grid::pushViewport(grid::viewport(x = 0.5, y = 0.625, width = 0.5, height = 0.25, just = c("center",
    "center")))
grid::grid.draw(g_pool)
grid::upViewport()
grid::pushViewport(grid::viewport(x = 0.5, y = 0.375, width = 0.5, height = 0.25, just = c("center",
    "center")))
grid::grid.draw(g_unit)
grid::upViewport()
grDevices::dev.off()

file.copy(pdf_file, legacy_pdf_file, overwrite = TRUE)
file.copy(pdf_file, alias_pdf_file, overwrite = TRUE)

message("PDF: ", normalizePath(pdf_file, mustWork = FALSE))
message("Legacy PDF updated: ", normalizePath(legacy_pdf_file, mustWork = FALSE))
message("Unit means: ", normalizePath(unit_mean_file, mustWork = FALSE))
