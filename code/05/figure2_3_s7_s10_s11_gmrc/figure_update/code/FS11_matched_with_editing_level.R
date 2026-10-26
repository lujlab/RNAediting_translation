#!/usr/bin/env Rscript

# Figure S11 sensitivity model with the Figure 3E editing-level covariate.
#
# The model input and the 1:20 transcript-matched pairs are reused verbatim
# from the validated Figure 3E/old Figure S11 run.  The only model change from
# the old Figure S11 is addition of editing_level_within_editSD, encoded in
# exactly the same way as Figure 3E: z-scored among edited codons separately
# for each sex and set to zero for unedited codons.

.libPaths(c("/path/to/RNA_editing_DATAprocessing_new/code_upload/new_mRC_gmRC/figure_update/R_library",
    .libPaths()))

suppressPackageStartupMessages({
    library(data.table)
    library(lme4)
    library(lmerTest)
})

workspace <- "/path/to/RNA_editing_DATAprocessing_new/code_upload"
figure_root <- file.path(workspace, "new_mRC_gmRC", "figure_update")
source_model_dir <- Sys.getenv("FS11_SOURCE_MODEL_DIR", unset = file.path(figure_root, "F3E"))
output_dir <- Sys.getenv("FS11_EDITLEVEL_OUTPUT_DIR", unset = file.path(figure_root, "FS11_with_editing_level"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "per_group"), recursive = TRUE, showWarnings = FALSE)

groups <- data.table(group_id = c("WBg_M_mean", "WBg_F_mean"), genotype = "WBg", sex = c("M",
    "F"), sex_label = c("Male", "Female"), n_replicates = 3L)

model_formula <- log_gmRC ~ edit_status_matched + editing_level_within_editSD + category_before +
    AA + phyloP_z + logTE_z + transcript_length_z + upstream_category + downstream_category +
    secondary_structure + (1 | FBtr) + (1 | match_id)

extract_fixed_effects <- function(model, meta, editing_definition) {
    coef_matrix <- coef(summary(model))
    result <- data.table(group_id = meta$group_id, genotype = meta$genotype, sex = meta$sex,
        n_replicates = meta$n_replicates, editing_definition = editing_definition, model_type = "GroupMean_Matched_with_editing_level",
        term = rownames(coef_matrix), estimate = coef_matrix[, "Estimate"], std_error = coef_matrix[,
            "Std. Error"], df = coef_matrix[, "df"], statistic = coef_matrix[, "t value"],
        p_value = coef_matrix[, "Pr(>|t|)"])
    result[, critical_value := qt(0.975, df = df)]
    result[, `:=`(conf_low = estimate - critical_value * std_error, conf_high = estimate +
        critical_value * std_error)]
    result[, `:=`(effect_percent = 100 * (exp(estimate) - 1), lower_percent = 100 * (exp(conf_low) -
        1), upper_percent = 100 * (exp(conf_high) - 1))]
    result[, critical_value := NULL]
    result
}

max_gradient <- function(model) {
    gradient <- model@optinfo$derivs$gradient
    if (is.null(gradient))
        return(NA_real_)
    max(abs(gradient))
}

convergence_message <- function(model) {
    value <- model@optinfo$conv$lme4$messages
    if (is.null(value))
        "OK" else paste(value, collapse = "; ")
}

all_fixed <- list()
all_qc <- list()

for (ii in seq_len(nrow(groups))) {
    meta <- groups[ii]
    group_id <- meta$group_id
    group_dir <- file.path(source_model_dir, "per_group", group_id)
    input_file <- file.path(group_dir, "group_mean_model_input.tsv.gz")
    pair_file <- file.path(group_dir, "matched_pairs.tsv.gz")

    if (!file.exists(input_file))
        stop("Missing model input: ", input_file)
    if (!file.exists(pair_file))
        stop("Missing matched pairs: ", pair_file)

    model_input <- fread(cmd = paste("gzip -dc", shQuote(input_file)))
    model_input[, row_id := .I]
    pairs <- fread(cmd = paste("gzip -dc", shQuote(pair_file)))

    case_map <- unique(pairs[, .(match_id, row_id = case_row_id, matched_role = "edited")])
    control_map <- pairs[, .(match_id, row_id = control_row_id, matched_role = "unedited")]
    matched_map <- rbindlist(list(case_map, control_map), use.names = TRUE)
    matched_data <- merge(matched_map, model_input, by = "row_id", all.x = TRUE, sort = FALSE)

    if (anyNA(matched_data$FBtr)) {
        stop(group_id, ": at least one saved match row did not map to model input")
    }

    matched_data[, edit_status_matched := factor(matched_role, levels = c("unedited",
        "edited"))]
    matched_data[, category_before := factor(category_before, levels = c("nonoptimal",
        "optimal"))]
    matched_data[, upstream_category := factor(upstream_category, levels = c("nonoptimal",
        "optimal"))]
    matched_data[, downstream_category := factor(downstream_category, levels = c("nonoptimal",
        "optimal"))]
    matched_data[, secondary_structure := factor(secondary_structure, levels = c("outside",
        "inside"))]
    matched_data[, AA := relevel(factor(AA), ref = "Ala")]
    matched_data[, FBtr := factor(FBtr)]
    matched_data[, match_id := factor(match_id)]
    setorder(matched_data, match_id, edit_status_matched)

    edited_levels <- matched_data[edit_status_matched == "edited", editing_level_within_editSD]
    input_edited_levels <- model_input[edit_status == "edited", editing_level_within_editSD]
    if (!isTRUE(all.equal(mean(input_edited_levels), 0, tolerance = 1e-10))) {
        stop(group_id, ": saved Figure 3E editing-level variable is not centered")
    }
    if (!isTRUE(all.equal(sd(input_edited_levels), 1, tolerance = 1e-10))) {
        stop(group_id, ": saved Figure 3E editing-level variable is not scaled")
    }
    if (any(matched_data[edit_status_matched == "unedited", editing_level_within_editSD] !=
        0)) {
        stop(group_id, ": unedited rows do not have editing-level value zero")
    }

    message("Fitting ", group_id, " with ", nrow(matched_data), " matched rows")
    model <- lmerTest::lmer(model_formula, data = matched_data, REML = FALSE, na.action = na.omit,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)))

    fixed <- extract_fixed_effects(model, meta, unique(model_input$editing_definition))
    all_fixed[[group_id]] <- fixed

    per_group_out <- file.path(output_dir, "per_group", group_id)
    dir.create(per_group_out, recursive = TRUE, showWarnings = FALSE)
    fwrite(fixed, file.path(per_group_out, "fixed_effects.tsv"), sep = "\t")
    fwrite(pairs, file.path(per_group_out, "matched_pairs.tsv.gz"), sep = "\t")
    saveRDS(model, file.path(per_group_out, "model.rds"))
    capture.output(summary(model), file = file.path(per_group_out, "summary.txt"))

    counts_per_case <- pairs[, .N, by = match_id]
    qc <- data.table(group_id = group_id, n_model_input = nrow(model_input), n_matched_rows = nrow(matched_data),
        n_matched_edited = sum(matched_data$edit_status_matched == "edited"), n_matched_unedited = sum(matched_data$edit_status_matched ==
            "unedited"), n_match_sets = uniqueN(matched_data$match_id), n_sets_with_20_controls = sum(counts_per_case$N ==
            20L), min_controls_per_set = min(counts_per_case$N), median_controls_per_set = median(counts_per_case$N),
        max_controls_per_set = max(counts_per_case$N), n_transcripts = uniqueN(matched_data$FBtr),
        mean_editing_level_SD_matched_edited = mean(edited_levels), sd_editing_level_SD_matched_edited = sd(edited_levels),
        n_fixed_effect_coefficients = nrow(fixed) - 1L, singular_fit = lme4::isSingular(model,
            tol = 1e-04), max_gradient = max_gradient(model), convergence = convergence_message(model),
        logLik = as.numeric(logLik(model)), AIC = AIC(model), BIC = BIC(model))
    all_qc[[group_id]] <- qc
}

fixed_all <- rbindlist(all_fixed, use.names = TRUE)
qc_all <- rbindlist(all_qc, use.names = TRUE)

fixed_all[, sex_label := fifelse(sex == "M", "Male", "Female")]
fixed_all[, p_adjust_bh_all := NA_real_]
fixed_all[term != "(Intercept)", p_adjust_bh_all := p.adjust(p_value, method = "BH"),
    by = sex_label]
fixed_all[, significance_bh_all := fifelse(is.finite(p_adjust_bh_all) & p_adjust_bh_all <
    0.05, "*", "")]

term_labels <- c(edit_status_matchededited = "Editing status", editing_level_within_editSD = "Editing level (+1 SD; edited)",
    category_beforeoptimal = "Focal codon preference", phyloP_z = "Codon phyloP", secondary_structureinside = "Stable hairpin",
    logTE_z = "Translation efficiency", transcript_length_z = "Transcript length", upstream_categoryoptimal = "Upstream codon preference",
    downstream_categoryoptimal = "Downstream codon preference")
display_order <- unname(term_labels)

displayed <- copy(fixed_all[term %in% names(term_labels)])
displayed[, display_variable := unname(term_labels[term])]
displayed[, display_index := match(display_variable, display_order)]
displayed[, sex_index := match(sex_label, c("Male", "Female"))]
setorder(displayed, display_index, sex_index)
displayed[, c("display_index", "sex_index") := NULL]

expected <- displayed[, .N, by = sex_label]
if (nrow(expected) != 2L || any(expected$N != 9L)) {
    stop("Expected nine displayed terms for each sex")
}

fwrite(fixed_all, file.path(output_dir, "all_fixed_effects_with_BH.tsv"), sep = "\t")
fwrite(displayed, file.path(output_dir, "FS11_with_editing_level_fixed_effects.tsv"), sep = "\t")
fwrite(qc_all, file.path(output_dir, "model_QC.tsv"), sep = "\t")

metadata <- data.table(parameter = c("analysis", "source_model_directory", "response", "fixed_effect_formula",
    "random_effects", "editing_level_definition", "length_covariate", "matched_sampling",
    "matched_pairs_reused", "BH_correction", "timestamp"), value = c("FS11_1_to_20_matched_with_F3E_editing_level",
    normalizePath(source_model_dir), "ln(1 + group-mean gmRC)", paste(deparse(nobars(model_formula)),
        collapse = " "), "random intercepts for FBtr and match_id", paste0("same saved variable as F3E: z-score among all eligible edited codons ",
        "within each sex; unedited codons set to 0"), "same F3E z-score of ln(full mature transcript length in nt), including both UTRs",
    paste0("each edited codon matched to up to 20 eligible unedited codons from ", "the same transcript; without replacement"),
    "yes; exact matched_pairs.tsv.gz files from the validated F3E/old FS11 run", paste0("BH across every non-intercept fixed-effect coefficient separately ",
        "within male and female models"), as.character(Sys.time())))
fwrite(metadata, file.path(output_dir, "run_metadata.tsv"), sep = "\t")
invisible(NULL)

readme <- c("Figure S11 with the Figure 3E editing-level term", "", "Updated 2026-10-05: the length covariate includes both UTRs.",
    "It reuses the exact saved 1:20 transcript-matched edited/unedited pairs from", "the validated F3E/old-FS11 analysis. The response, eligibility filters,",
    "reference levels and gmRC inputs are unchanged in the length-only update.", "Length uses the saved F3E transcript_length_z: standardized ln(tx_len in nt).",
    "", "The added editing-level covariate is identical to Figure 3E:", "editing_level_within_editSD is standardized among all eligible edited codons",
    "separately for WT male and WT female; unedited codons are assigned zero.", "The matched model includes random intercepts for transcript (FBtr) and match set.",
    "Raw coefficient P values and BH-adjusted P values are both in the result table;", "stars in the PDF use BH-adjusted P < 0.05 across all non-intercept fixed effects",
    "within each sex-specific model.")
invisible(NULL)

message("Completed: ", normalizePath(output_dir))
