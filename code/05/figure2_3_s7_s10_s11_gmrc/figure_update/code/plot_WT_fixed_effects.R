#!/usr/bin/env Rscript

# WT male/female fixed-effect bar plots with BH correction across all non-intercept fixed effects within each model and sex.
#
# This block can be pasted directly at the end of
# gmRC_groupMean_WBg_Adar5G1_pseudoEdit_model.R.  It uses the `all_fixed`
# and `output_dir` objects created by that script.  It can also be sourced
# independently after the model has finished; in that case it reads the
# exported fixed-effect table.

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required.")
}

if (!exists("output_dir", inherits = TRUE)) {
    output_dir <- Sys.getenv("GMRC_OUTPUT_DIR", unset = paste0("/path/to/RNA_editing_DATAprocessing_new/",
        "10.stalling/04_GLM/", "gmRC_groupMean_WBg_Adar5G1_merged_tsv_sites_final"))
}

if (!exists("all_fixed", inherits = TRUE)) {
    fixed_file <- file.path(output_dir, "all_group_mean_fixed_effects.tsv")
    if (!file.exists(fixed_file)) {
        stop("Missing fixed-effect table: ", fixed_file)
    }
    all_fixed <- data.table::fread(fixed_file)
}

plot_output_dir <- file.path(output_dir, "WT_sex_BH_plots", "latest_all_variables_BH")
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c("group_id", "genotype", "sex", "model_type", "term", "estimate", "p_value")
missing_columns <- setdiff(required_columns, names(all_fixed))
if (length(missing_columns) > 0L) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

term_label_map_bh <- c(edit_statusedited = "Editing status", edit_status_matchededited = "Editing status",
    category_beforeoptimal = "Focal codon preference", phyloP_z = "Codon phyloP", secondary_structureinside = "Stable hairpin",
    logTE_z = "Translation efficiency", transcript_length_z = "Transcript length", upstream_categoryoptimal = "Upstream codon preference",
    downstream_categoryoptimal = "Downstream codon preference")

variable_order_bh <- c("Editing status", "Focal codon preference", "Codon phyloP", "Stable hairpin",
    "Translation efficiency", "Transcript length", "Upstream codon preference", "Downstream codon preference")

all_wt_fixed_bh <- data.table::copy(all_fixed[group_id %in% c("WBg_M_mean", "WBg_F_mean") &
    model_type %in% c("GroupMean_Full", "GroupMean_Matched") & term != "(Intercept)"])

all_wt_fixed_bh[, sex_label := data.table::fcase(group_id == "WBg_M_mean", "Male", group_id ==
    "WBg_F_mean", "Female", default = NA_character_)]
all_wt_fixed_bh[, sex_label := factor(sex_label, levels = c("Male", "Female"))]
all_wt_fixed_bh[, p_adjust_bh_all := p.adjust(p_value, method = "BH"), by = .(model_type,
    sex_label)]
all_wt_fixed_bh[, significance_bh_all := data.table::fifelse(is.finite(p_adjust_bh_all) &
    p_adjust_bh_all < 0.05, "*", "")]

data.table::fwrite(all_wt_fixed_bh, file.path(plot_output_dir, "all_WT_fixed_effects_BH_all_variables.tsv"),
    sep = "\t")

plot_dat_bh <- data.table::copy(all_wt_fixed_bh[term %in% names(term_label_map_bh)])

plot_dat_bh[, variable := unname(term_label_map_bh[term])]
plot_dat_bh[, variable := factor(variable, levels = rev(variable_order_bh))]

# Place the star just beyond the end of each bar.
plot_dat_bh[, star_position := estimate + data.table::fifelse(estimate >= 0, 0.004, -0.004)]
plot_dat_bh[, star_hjust := data.table::fifelse(estimate >= 0, 0, 1)]

expected_counts_bh <- plot_dat_bh[, .N, by = .(model_type, sex_label)]
if (nrow(expected_counts_bh) != 4L || any(expected_counts_bh$N != length(variable_order_bh))) {
    stop("Expected exactly eight displayed terms for each model-by-sex family. ", "Observed: ",
        paste(paste0(expected_counts_bh$model_type, "/", expected_counts_bh$sex_label, "=",
            expected_counts_bh$N), collapse = "; "))
}

data.table::setorder(plot_dat_bh, model_type, sex_label, variable)
data.table::fwrite(plot_dat_bh, file.path(plot_output_dir, "plot_group_mean_common_effects_WT_sex_bar_BH_all_variables.tsv"),
    sep = "\t")

sex_colors_bh <- c(Male = "#4C78A8", Female = "#E78AC3")

make_fixed_effect_bar_bh <- function(model_type_i, title_i) {
    dat_i <- plot_dat_bh[model_type == model_type_i]
    dodge_i <- ggplot2::position_dodge(width = 0.7)

    ggplot2::ggplot(dat_i, ggplot2::aes(x = variable, y = estimate, fill = sex_label, group = sex_label)) +
        ggplot2::geom_col(width = 0.6, position = dodge_i) + ggplot2::geom_text(ggplot2::aes(y = star_position,
        label = significance_bh_all, hjust = star_hjust), position = dodge_i, size = 5, fontface = "bold",
        color = "#231815", show.legend = FALSE) + ggplot2::coord_flip(ylim = c(-0.05, 0.2),
        clip = "off") + ggplot2::scale_fill_manual(values = sex_colors_bh, breaks = c("Male",
        "Female"), drop = FALSE) + ggplot2::scale_y_continuous(breaks = seq(-0.05, 0.2, by = 0.05),
        labels = function(x) sprintf("%.2f", x), expand = ggplot2::expansion(mult = c(0,
            0))) + ggplot2::labs(x = NULL, y = "Fixed-effect coefficient for ln(1 + gmRC)",
        fill = "Sex", title = title_i, caption = NULL) + ggplot2::theme_classic(base_size = 14,
        base_family = "Arial") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        face = "bold", size = 16), axis.title.x = ggplot2::element_text(size = 15), axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 14, color = "#231815"), axis.line = ggplot2::element_line(color = "#231815",
            linewidth = 0.6), axis.ticks = ggplot2::element_line(color = "#231815", linewidth = 0.6),
        legend.position = c(0.82, 0.53), legend.title = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 13), plot.caption = ggplot2::element_text(hjust = 1,
            size = 10, color = "#231815"), plot.margin = ggplot2::margin(12, 24, 8, 12))
}

p_full_bh <- make_fixed_effect_bar_bh("GroupMean_Full", "Full dataset")
p_matched_bh <- make_fixed_effect_bar_bh("GroupMean_Matched", "1:20 edited-to-unedited sampling ratio")

ggplot2::ggsave(filename = file.path(plot_output_dir, "group_mean_common_effects_WT_sex_full_bar_BH_all_variables.pdf"),
    plot = p_full_bh, device = grDevices::cairo_pdf, width = 7.2, height = 4.7, units = "in")
ggplot2::ggsave(filename = file.path(plot_output_dir, "group_mean_common_effects_WT_sex_1_20_bar_BH_all_variables.pdf"),
    plot = p_matched_bh, device = grDevices::cairo_pdf, width = 7.2, height = 4.7, units = "in")

message("BH-adjusted WT sex plots written to: ", normalizePath(plot_output_dir))


# ============================================================
# Amino-acid effects: same BH correction, amino acids on x-axis
# ============================================================

aa_order_bh <- c("Ala", "Pro", "Arg", "His", "Asp", "Leu", "Thr", "Ile", "Lys", "Asn", "Gly",
    "Met", "Tyr", "Ser", "Glu", "Gln", "Val")

aa_dat_bh <- all_wt_fixed_bh[grepl("^AA", term), .(group_id, genotype, sex, sex_label, model_type,
    term, estimate, std_error, df, statistic, p_value, p_adjust_bh_all, significance_bh_all)]
aa_dat_bh[, amino_acid := sub("^AA", "", term)]

# Ala is the reference: coefficient 0, with no coefficient-level P value.
aa_reference_bh <- unique(aa_dat_bh[, .(group_id, genotype, sex, sex_label, model_type)])
aa_reference_bh[, `:=`(term = "AAAla_reference", estimate = 0, std_error = NA_real_, df = NA_real_,
    statistic = NA_real_, p_value = NA_real_, p_adjust_bh_all = NA_real_, significance_bh_all = "",
    amino_acid = "Ala")]

aa_dat_bh <- data.table::rbindlist(list(aa_dat_bh, aa_reference_bh), use.names = TRUE, fill = TRUE)
aa_dat_bh[, amino_acid := factor(amino_acid, levels = aa_order_bh)]
aa_dat_bh[, star_position := estimate + data.table::fifelse(estimate > 0, 0.008, data.table::fifelse(estimate <
    0, -0.008, 0))]
aa_dat_bh[, star_vjust := data.table::fifelse(estimate >= 0, 0, 1)]

aa_counts_bh <- aa_dat_bh[, .N, by = .(model_type, sex_label)]
if (nrow(aa_counts_bh) != 4L || any(aa_counts_bh$N != length(aa_order_bh))) {
    stop("Expected 17 amino-acid levels, including Ala, per family.")
}

data.table::fwrite(aa_dat_bh, file.path(plot_output_dir, "plot_group_mean_amino_acid_effects_WT_sex_bar_BH_all_variables.tsv"),
    sep = "\t")

make_aa_bar_bh <- function(model_type_i, title_i) {
    dat_i <- aa_dat_bh[model_type == model_type_i]
    dodge_i <- ggplot2::position_dodge(width = 0.78)

    ggplot2::ggplot(dat_i, ggplot2::aes(x = amino_acid, y = estimate, fill = sex_label, group = sex_label)) +
        ggplot2::geom_hline(yintercept = 0, color = "#231815", linewidth = 0.45) + ggplot2::geom_col(width = 0.68,
        position = dodge_i) + ggplot2::geom_text(ggplot2::aes(y = star_position, label = significance_bh_all,
        vjust = star_vjust), position = dodge_i, size = 4.6, fontface = "bold", color = "#231815",
        show.legend = FALSE) + ggplot2::scale_fill_manual(values = sex_colors_bh, breaks = c("Male",
        "Female"), drop = FALSE) + ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_continuous(breaks = seq(-0.15,
        0.15, by = 0.05), labels = function(x) sprintf("%.2f", x), limits = c(-0.18, 0.18),
        expand = ggplot2::expansion(mult = c(0, 0))) + ggplot2::labs(x = "Amino acid", y = paste0("Fixed-effect coefficient for ln(1 + gmRC)\n",
        "relative to Ala"), fill = "Sex", title = title_i, caption = NULL) + ggplot2::theme_classic(base_size = 14,
        base_family = "Arial") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        face = "bold", size = 16), axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 11), axis.line = ggplot2::element_line(color = "#231815",
            linewidth = 0.6), axis.ticks = ggplot2::element_line(color = "#231815", linewidth = 0.6),
        legend.position = "top", legend.title = ggplot2::element_text(size = 14), legend.text = ggplot2::element_text(size = 13),
        plot.caption = ggplot2::element_text(hjust = 1, size = 9.2, color = "#231815"), plot.margin = ggplot2::margin(12,
            24, 8, 12))
}

p_aa_full_bh <- make_aa_bar_bh("GroupMean_Full", "Full dataset")
p_aa_matched_bh <- make_aa_bar_bh("GroupMean_Matched", "1:20 edited-to-unedited sampling ratio")

ggplot2::ggsave(filename = file.path(plot_output_dir, "group_mean_amino_acid_effects_WT_sex_full_bar_BH_all_variables.pdf"),
    plot = p_aa_full_bh, device = grDevices::cairo_pdf, width = 10.5, height = 5.2, units = "in")
ggplot2::ggsave(filename = file.path(plot_output_dir, "group_mean_amino_acid_effects_WT_sex_1_20_bar_BH_all_variables.pdf"),
    plot = p_aa_matched_bh, device = grDevices::cairo_pdf, width = 10.5, height = 5.2, units = "in")

message("All-variable BH amino-acid plots written to: ", normalizePath(plot_output_dir))


# ============================================================
# Additional Full plot including editing level
# ============================================================

full_level_label_map_bh <- c(edit_statusedited = "Editing status", editing_level_within_editSD = "Editing level (+1 SD; edited)",
    category_beforeoptimal = "Focal codon preference", phyloP_z = "Codon phyloP", secondary_structureinside = "Stable hairpin",
    logTE_z = "Translation efficiency", transcript_length_z = "Transcript length", upstream_categoryoptimal = "Upstream codon preference",
    downstream_categoryoptimal = "Downstream codon preference")

full_level_order_bh <- c("Editing status", "Editing level (+1 SD; edited)", "Focal codon preference",
    "Codon phyloP", "Stable hairpin", "Translation efficiency", "Transcript length", "Upstream codon preference",
    "Downstream codon preference")

full_level_dat_bh <- data.table::copy(all_wt_fixed_bh[model_type == "GroupMean_Full" & term %in%
    names(full_level_label_map_bh)])
full_level_dat_bh[, variable := unname(full_level_label_map_bh[term])]
full_level_dat_bh[, variable := factor(variable, levels = rev(full_level_order_bh))]
full_level_dat_bh[, star_position := estimate + data.table::fifelse(estimate >= 0, 0.004,
    -0.004)]
full_level_dat_bh[, star_hjust := data.table::fifelse(estimate >= 0, 0, 1)]

full_level_counts_bh <- full_level_dat_bh[, .N, by = sex_label]
if (nrow(full_level_counts_bh) != 2L || any(full_level_counts_bh$N != length(full_level_order_bh))) {
    stop("Expected nine Full-model terms for each sex.")
}

data.table::fwrite(full_level_dat_bh, file.path(plot_output_dir, "plot_group_mean_common_effects_WT_sex_full_with_editing_level_BH_all_variables.tsv"),
    sep = "\t")

dodge_full_level_bh <- ggplot2::position_dodge(width = 0.7)
p_full_with_level_bh <- ggplot2::ggplot(full_level_dat_bh, ggplot2::aes(x = variable, y = estimate,
    fill = sex_label, group = sex_label)) + ggplot2::geom_col(width = 0.6, position = dodge_full_level_bh) +
    ggplot2::geom_text(ggplot2::aes(y = star_position, label = significance_bh_all, hjust = star_hjust),
        position = dodge_full_level_bh, size = 5, fontface = "bold", color = "#231815", show.legend = FALSE) +
    ggplot2::coord_flip(ylim = c(-0.05, 0.2), clip = "off") + ggplot2::scale_fill_manual(values = sex_colors_bh,
    breaks = c("Male", "Female"), drop = FALSE) + ggplot2::scale_y_continuous(breaks = seq(-0.05,
    0.2, by = 0.05), labels = function(x) sprintf("%.2f", x), expand = ggplot2::expansion(mult = c(0,
    0))) + ggplot2::labs(x = NULL, y = "Fixed-effect coefficient for ln(1 + gmRC)", fill = "Sex",
    title = "Full dataset", caption = NULL) + ggplot2::theme_classic(base_size = 14, base_family = "Arial") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = ggplot2::element_text(size = 15), axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 13, color = "#231815"), axis.line = ggplot2::element_line(color = "#231815",
            linewidth = 0.6), axis.ticks = ggplot2::element_line(color = "#231815", linewidth = 0.6),
        legend.position = c(0.82, 0.53), legend.title = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 13), plot.caption = ggplot2::element_text(hjust = 1,
            size = 9.2, color = "#231815"), plot.margin = ggplot2::margin(12, 24, 8, 12))

ggplot2::ggsave(filename = file.path(plot_output_dir, "group_mean_common_effects_WT_sex_full_with_editing_level_BH_all_variables.pdf"),
    plot = p_full_with_level_bh, device = grDevices::cairo_pdf, width = 7.6, height = 5.2,
    units = "in")

message("Full plot with editing level written to: ", normalizePath(plot_output_dir))
