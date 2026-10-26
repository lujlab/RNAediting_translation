source("config.R")
source("R/utils.R")
source("R/read_xlsx_base.R")
make_output_dirs(out_dir)
prep_file <- file.path(out_dir, "01_gene_mapping", "prepared_annotation_and_editing.rds")
if (!file.exists(prep_file)) stop("Run 01_build_features.R first.")
obj <- readRDS(prep_file)
gene_features <- as.data.table(obj$gene_features)
transcript_features <- as.data.table(obj$transcript_features)
editing_gene <- as.data.table(obj$editing_gene)
callable <- as.data.table(obj$callable)
cg_mapping <- as.data.table(obj$cg_mapping)

zscore <- function(x) {
    s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else as.numeric((x - mean(x, na.rm = TRUE)) / s)
}
map_to_fbgn <- function(id1, id2 = NULL) {
    lookup <- unique(gene_features[, .(FBgn, gene_symbol)])
    out <- clean_id(id1); direct <- grepl("^FBgn", out)
    if (any(!direct)) {
        old <- out[!direct]
        mapped <- if (nrow(cg_mapping)) cg_mapping$FBgn[match(old, cg_mapping$annotation_id)] else rep(NA_character_, length(old))
        missing <- is.na(mapped); mapped[missing] <- lookup$FBgn[match(old[missing], lookup$gene_symbol)]
        out[!direct] <- mapped
    }
    if (!is.null(id2)) {
        missing <- is.na(out) | out == ""; alt <- clean_id(id2); out[missing] <- lookup$FBgn[match(alt[missing], lookup$gene_symbol)]
    }
    out
}
add_editing <- function(dt) {
    x <- merge(dt, editing_gene, by = "FBgn", all.x = TRUE)
    x[is.na(editing_site_count), editing_site_count := 0L]
    x[, edited := fifelse(
        editing_site_count >= min_editing_sites_per_gene, 1L,
        fifelse(editing_site_count == 0L, 0L, NA_integer_)
    )]
    x[, editing_binary_status := fifelse(
        editing_site_count >= min_editing_sites_per_gene, "Edited",
        fifelse(
            editing_site_count == 0L, "No_detected_editing",
            sprintf("Excluded_below_%d_site_threshold", min_editing_sites_per_gene)
        )
    )]
    x
}
read_background_B2 <- function() {
    if (is.null(editing_expression_file) || !file.exists(editing_expression_file)) return(data.table(FBgn = character()))
    x <- standardize_gene_column(fread(editing_expression_file), "editing-expression B2")
    ec <- find_col(x, c("mean_normalized_count", "normalized_count", "mean_count", "TPM", "FPKM", "expression"))
    x[, editing_experiment_expression := suppressWarnings(as.numeric(get(ec)))]
    unique(x[is.finite(editing_experiment_expression) & editing_experiment_expression >= editing_expression_threshold, .(FBgn)])
}
b2_genes <- read_background_B2()
b3_genes <- if (nrow(callable)) unique(callable[callable == TRUE, .(FBgn)]) else data.table(FBgn = character())
background_names <- c("B1_All_detected_background", if (nrow(b2_genes)) "B2_Expression_detectable_background", if (nrow(b3_genes)) "B3_A_site_coverage_callable_background")
apply_background <- function(dt, bg) {
    x <- copy(dt)
    # Detected genes below the strict binary threshold are not controls.
    x <- x[!is.na(edited)]
    if (startsWith(bg, "B2")) x <- x[FBgn %in% b2_genes$FBgn]
    if (startsWith(bg, "B3")) x <- x[FBgn %in% b3_genes$FBgn]
    control <- if (startsWith(bg, "B1")) "No_detected_editing" else if (startsWith(bg, "B2")) "Expression_detectable_unedited" else "Callable_unedited"
    x[, editing_group := factor(ifelse(edited == 1, "Edited", control), levels = c(control, "Edited"))]
    x[, background_definition := bg]
    x
}
standardize_model_vars <- function(x) {
    x[, `:=`(z_log2_transcript_length = zscore(safe_log2(transcript_length)), z_log2_published_expression = zscore(safe_log2(published_expression, 1)), z_log2_CDS_length = zscore(safe_log2(CDS_length)), z_log2_UTR5_length = zscore(safe_log2(UTR5_length, 1)), z_log2_UTR3_length = zscore(safe_log2(UTR3_length, 1)), z_GC_content = zscore(GC_content), z_A_content = zscore(A_content))]
    x
}

# ---- 2010 processed microarray, decay classes, half-lives and expression ----
parse_2010 <- function(file) {
    if (!file.exists(file)) return(NULL)
    raw <- fread(file)
    class_col <- find_col(raw, c("class.I.to.V", "decay_class", "class")); cg_col <- find_col(raw, c("CG.gene.ID", "gene", "gene_id")); sym_col <- find_col(raw, c("Gene.symbol", "symbol"), FALSE)
    raw[, decay_class := toupper(trimws(as.character(get(class_col))))]
    raw[!decay_class %chin% c("I", "II", "III", "IV", "V"), decay_class := NA_character_]
    u2 <- grep("^U2[._]", names(raw), value = TRUE); u3 <- grep("^U3[._]", names(raw), value = TRUE); e2 <- grep("^E2[._]", names(raw), value = TRUE); e3 <- grep("^E3[._]", names(raw), value = TRUE)
    if (!all(lengths(list(u2, u3, e2, e3)) > 0)) stop("2010 processed table lacks U2/U3/E2/E3 replicate columns.")
    for (v in c(u2, u3, e2, e3)) set(raw, j = v, value = suppressWarnings(as.numeric(raw[[v]])))
    raw[, `:=`(U2_log2 = rowMeans(.SD, na.rm = TRUE)), .SDcols = u2]
    raw[, U3_log2 := rowMeans(.SD, na.rm = TRUE), .SDcols = u3]
    raw[, E2_log2 := rowMeans(.SD, na.rm = TRUE), .SDcols = e2]
    raw[, E3_log2 := rowMeans(.SD, na.rm = TRUE), .SDcols = e3]
    raw[, half_life_min := fifelse(decay_class %chin% c("II", "III"), 60 / (U2_log2 - U3_log2), fifelse(decay_class %chin% c("IV", "V"), 60 / (E2_log2 - E3_log2), NA_real_))]
    raw[!is.finite(half_life_min) | half_life_min <= 0, half_life_min := NA_real_]
    raw[, log2_baseline_expression := fifelse(decay_class %chin% c("II", "III"), U2_log2, E2_log2)]
    raw[, published_expression := 2^log2_baseline_expression]
    id2 <- if (!is.na(sym_col)) raw[[sym_col]] else NULL
    raw[, FBgn := map_to_fbgn(get(cg_col), id2)]
    raw[, source_gene := as.character(get(cg_col))]
    conflicts <- raw[!is.na(FBgn) & !is.na(decay_class), .(n_probe_sets = .N, n_decay_classes = uniqueN(decay_class), classes = paste(sort(unique(decay_class)), collapse = ";")), by = FBgn][n_decay_classes > 1]
    write_tsv(conflicts, file.path(out_dir, "00_input_QC", "2010_conflicting_probe_classifications.tsv"))
    valid <- raw[!is.na(FBgn) & !FBgn %in% conflicts$FBgn & !is.na(decay_class)]
    gene <- valid[, .(gene = source_gene[1], decay_class = decay_class[1], half_life_min = median(half_life_min, na.rm = TRUE), published_expression = median(published_expression, na.rm = TRUE), log2_baseline_expression = median(log2_baseline_expression, na.rm = TRUE), expression_source = ifelse(decay_class[1] %chin% c("II", "III"), "E-MEXP-2580 processed U2 mean (class-specific decay-series baseline)", "E-MEXP-2580 processed E2 mean (embryo-series baseline)"), n_probe_sets = .N), by = FBgn]
    for (v in c("half_life_min", "published_expression", "log2_baseline_expression")) gene[!is.finite(get(v)), (v) := NA_real_]
    gene <- merge(gene, gene_features, by = "FBgn", all.x = TRUE)
    gene <- add_editing(gene); gene[, decay_status := as.integer(decay_class != "I")]; gene[, log2_half_life := safe_log2(half_life_min)]; gene[, decay_class := factor(decay_class, levels = c("II", "III", "IV", "V", "I"))]
    gene <- standardize_model_vars(gene)
    write_tsv(
        gene[is.na(edited), .(FBgn, gene_symbol, editing_site_count, editing_binary_status, decay_class, half_life_min)],
        file.path(out_dir, "00_input_QC", "2010_editing_threshold_excluded_genes.tsv")
    )
    write_tsv(gene[, .(gene = FBgn, baseline_expression = published_expression, log2_baseline_expression, expression_source)], file.path(out_dir, "00_input_QC", "2010_published_expression.tsv"))
    write_tsv(raw[, .(source_gene, FBgn, decay_class, half_life_min, log2_baseline_expression)], file.path(out_dir, "00_input_QC", "2010_probe_level_processed.tsv"))
    probe_medians <- raw[decay_class %chin% c("II", "III", "IV", "V") & is.finite(half_life_min), .(N_probe_sets = .N, median_half_life_min = median(half_life_min)), by = decay_class]
    probe_medians[, expected_approximate_median_from_paper := c(64, 31, 133, 38)[match(decay_class, c("II", "III", "IV", "V"))]]
    write_tsv(probe_medians, file.path(out_dir, "00_input_QC", "2010_published_median_reproduction.tsv"))
    write_tsv(data.table(metric = c("processed_probe_rows", "mapped_probe_rows", "gene_rows", "conflicting_genes_excluded", "finite_half_lives"), value = c(nrow(raw), sum(!is.na(raw$FBgn)), nrow(gene), nrow(conflicts), sum(is.finite(gene$half_life_min)))), file.path(out_dir, "00_input_QC", "2010_input_QC.tsv"))
    gene
}

# Add group-specific N to model tables and retain exact complete-case membership.
fit_lm_record <- function(dt, formula, name, dataset, cc_store, diag_store, vif_store) {
    vars <- all.vars(formula); report <- complete_case_report(dt, vars, dataset, name); cc_store[[length(cc_store) + 1L]] <- report
    d <- dt[report$included]
    if (nrow(d) < max(8, length(vars) + 3) || uniqueN(d$edited) < 2) return(list(result = data.table(dataset = dataset, model = name, variable = "MODEL_SKIPPED", estimate = NA_real_, SE = NA_real_, ordinary_P = NA_real_, robust_SE_HC3 = NA_real_, robust_P_HC3 = NA_real_, CI_low95 = NA_real_, CI_high95 = NA_real_, N_total = nrow(d), adjusted_R2 = NA_real_, AIC = NA_real_, N_edited = sum(d$edited == 1), N_unedited = sum(d$edited == 0), note = "Insufficient complete cases or only one editing group"), model = NULL))
    m <- lm(formula, data = d); tab <- tidy_lm(m, name, dataset); tab[, `:=`(N_edited = sum(d$edited == 1), N_unedited = sum(d$edited == 0), note = "")]
    vif_store[[length(vif_store) + 1L]] <- vif_manual(m, name, dataset)
    diag_store[[length(diag_store) + 1L]] <- diagnose_lm(m, name, dataset, file.path(out_dir, "05_2010_regression", paste0(gsub("[^A-Za-z0-9]+", "_", paste(dataset, name)), "_diagnostics.pdf")))
    list(result = tab, model = m)
}
fit_glm_record <- function(dt, formula, name, dataset, cc_store) {
    vars <- all.vars(formula); report <- complete_case_report(dt, vars, dataset, name); cc_store[[length(cc_store) + 1L]] <- report; d <- dt[report$included]
    if (nrow(d) < max(10, length(vars) + 4) || uniqueN(d$edited) < 2 || uniqueN(d$decay_status) < 2) return(data.table(dataset = dataset, model = name, variable = "MODEL_SKIPPED", estimate = NA_real_, SE = NA_real_, OR = NA_real_, OR_low95 = NA_real_, OR_high95 = NA_real_, P = NA_real_, N_total = nrow(d), AIC = NA_real_, N_edited = sum(d$edited == 1), N_unedited = sum(d$edited == 0), note = "Insufficient complete cases or outcome variation"))
    m <- glm(formula, family = binomial(), data = d); z <- tidy_glm(m, name, dataset); z[, `:=`(N_edited = sum(d$edited == 1), N_unedited = sum(d$edited == 0), note = "")]; z
}

run_2010_background <- function(all2010, bg, full_outputs = FALSE) {
    x <- apply_background(all2010, bg); decay <- x[decay_class %chin% c("II", "III", "IV", "V") & is.finite(half_life_min)]
    cc <- list(); diags <- list(); vifs <- list(); lm_results <- list(); lm_models <- list(); glm_results <- list()
    # Stable/decay contingency and logistic models.
    tab <- table(factor(x$edited, levels = c(1, 0)), factor(x$decay_status, levels = c(1, 0))); ft <- fisher.test(tab); ft_one <- fisher.test(tab, alternative = "greater")
    fisher <- data.table(background_definition = bg, N_total = nrow(x), N_edited = sum(x$edited == 1), N_control = sum(x$edited == 0), Edited_decay = tab[1, 1], Edited_stable = tab[1, 2], Control_decay = tab[2, 1], Control_stable = tab[2, 2], odds_ratio = unname(ft$estimate), OR_low95 = ft$conf.int[1], OR_high95 = ft$conf.int[2], two_sided_P = ft$p.value, one_sided_greater_P = ft_one$p.value, one_sided_alternative = "Edited odds of Decay > control odds of Decay", BH_adjusted_two_sided_P = p.adjust(ft$p.value, "BH"))
    glm_results[["L0"]] <- fit_glm_record(x, decay_status ~ edited, "L0", bg, cc)
    glm_results[["L1"]] <- fit_glm_record(x, decay_status ~ edited + z_log2_transcript_length, "L1", bg, cc)
    if (sum(is.finite(x$published_expression)) > 0) glm_results[["L2"]] <- fit_glm_record(x, decay_status ~ edited + z_log2_transcript_length + z_log2_published_expression, "L2", bg, cc) else glm_results[["L2"]] <- data.table(dataset = bg, model = "L2", variable = "MODEL_SKIPPED", note = "Published expression unavailable; no substitute used")
    # Half-life simple comparison.
    ed <- decay[edited == 1, half_life_min]; co <- decay[edited == 0, half_life_min]
    wt <- if (length(ed) && length(co)) wilcox.test(ed, co, exact = FALSE) else NULL
    simple <- decay[, .(n = .N, median = median(half_life_min), mean = mean(half_life_min), SD = sd(half_life_min), IQR = IQR(half_life_min), min = min(half_life_min), max = max(half_life_min)), by = editing_group]
    simple[, `:=`(background_definition = bg, Wilcoxon_P = if (is.null(wt)) NA_real_ else wt$p.value, rank_biserial_edited_vs_control = if (length(ed) && length(co)) (sum(outer(ed, co, ">")) + 0.5 * sum(outer(ed, co, "=="))) / (length(ed) * length(co)) * 2 - 1 else NA_real_)]
    forms <- list(H0 = log2_half_life ~ edited, H1 = log2_half_life ~ edited + z_log2_transcript_length)
    if (sum(is.finite(decay$published_expression)) > 0) {
        forms$H2 <- log2_half_life~edited + z_log2_transcript_length + z_log2_published_expression
        forms$H3 <- log2_half_life~edited + z_log2_transcript_length + z_log2_published_expression + decay_class
        forms$H4_component_length_model <- log2_half_life~edited + z_log2_CDS_length + z_log2_UTR5_length + z_log2_UTR3_length + z_log2_published_expression + decay_class
        forms$H5_extended_model <- log2_half_life~edited + z_log2_transcript_length + z_log2_published_expression + z_GC_content + z_A_content + decay_class
    }
    for (nm in names(forms)) {
        z <- fit_lm_record(decay, forms[[nm]], nm, bg, cc, diags, vifs); lm_results[[nm]] <- z$result; lm_models[[nm]] <- z$model
    }
    # Interaction is exploratory and kept separate from primary main-effect models.
    int_form <- if (sum(is.finite(decay$published_expression)) > 0) log2_half_life ~ edited * decay_class + z_log2_transcript_length + z_log2_published_expression else log2_half_life ~ edited * decay_class + z_log2_transcript_length
    z <- fit_lm_record(decay, int_form, "H_interaction_exploratory", bg, cc, diags, vifs); lm_results[["interaction"]] <- z$result; lm_models[["interaction"]] <- z$model
    lm_tab <- rbindlist(lm_results, fill = TRUE); lm_tab[, BH_adjusted_P := p.adjust(ordinary_P, "BH")]
    glm_tab <- rbindlist(glm_results, fill = TRUE); if ("P" %in% names(glm_tab)) glm_tab[, BH_adjusted_P := p.adjust(P, "BH")]
    if (full_outputs) {
        cc_out <- rbindlist(c(lapply(names(forms), function(nm) complete_case_report(decay, all.vars(forms[[nm]]), bg, nm)), list(complete_case_report(decay, all.vars(int_form), bg, "H_interaction_exploratory"))), fill = TRUE)
        vif_out <- rbindlist(lapply(names(lm_models), function(nm) if (!is.null(lm_models[[nm]])) vif_manual(lm_models[[nm]], nm, bg) else NULL), fill = TRUE)
        diag_out <- rbindlist(lapply(names(lm_models), function(nm) if (!is.null(lm_models[[nm]])) diagnose_lm(lm_models[[nm]], nm, bg, file.path(out_dir, "05_2010_regression", paste0(nm, "_diagnostics.pdf"))) else NULL), fill = TRUE)
        if (!is.null(lm_models$H3) && !is.null(lm_models$interaction)) {
            ia <- anova(lm_models$H3, lm_models$interaction)
            write_tsv(data.table(background_definition = bg, N_total = nobs(lm_models$interaction), N_edited = sum(model.frame(lm_models$interaction)$edited == 1), N_control = sum(model.frame(lm_models$interaction)$edited == 0), reduced_model = "H3", interaction_model = "H_interaction_exploratory", numerator_df = ia$Df[2], F = ia$F[2], interaction_global_P = ia$`Pr(>F)`[2]), file.path(out_dir, "05_2010_regression", "2010_editing_decay_class_interaction_test.tsv"))
        }
        write_tsv(fisher, file.path(out_dir, "03_2010_stable_decay", "2010_stable_decay_fisher.tsv")); write_tsv(glm_tab, file.path(out_dir, "03_2010_stable_decay", "2010_decay_logistic_regression.tsv")); write_tsv(simple, file.path(out_dir, "04_2010_half_life", "2010_half_life_edited_vs_unedited.tsv")); write_tsv(lm_tab, file.path(out_dir, "05_2010_regression", "2010_half_life_regression_all_models.tsv")); write_tsv(vif_out, file.path(out_dir, "05_2010_regression", "2010_regression_VIF.tsv")); write_tsv(diag_out, file.path(out_dir, "05_2010_regression", "2010_regression_diagnostics.tsv")); write_tsv(cc_out, file.path(out_dir, "05_2010_regression", "model_complete_cases.tsv"))
    }
    primary_name <- if ("H3" %in% names(lm_models) && !is.null(lm_models$H3)) "H3" else "H1"
    pe <- lm_tab[model == primary_name & variable == "edited"]
    sens <- data.table(background_definition = bg, N_edited = length(ed), N_control = length(co), median_half_life_edited = if (length(ed)) median(ed) else NA_real_, median_half_life_control = if (length(co)) median(co) else NA_real_, difference = if (length(ed) && length(co)) median(ed) - median(co) else NA_real_, Wilcoxon_P = if (is.null(wt)) NA_real_ else wt$p.value, regression_model = primary_name, editing_beta = pe$estimate[1], CI_low = pe$CI_low95[1], CI_high = pe$CI_high95[1], regression_P = pe$ordinary_P[1], robust_P = pe$robust_P_HC3[1])
    list(data = x, decay = decay, fisher = fisher, logistic = glm_tab, simple = simple, regression = lm_tab, models = lm_models, background = sens)
}

data2010 <- parse_2010(half_life_2010_file)
res2010 <- list()
if (!is.null(data2010)) {
    for (bg in background_names) res2010[[bg]] <- run_2010_background(data2010, bg, full_outputs = startsWith(bg, "B1"))
    write_tsv(res2010[[1]]$data, file.path(out_dir, "10_tables", "all_genes_merged_metadata.tsv"))
    bg_sensitivity <- rbindlist(lapply(res2010, `[[`, "background"), fill = TRUE)
    write_tsv(bg_sensitivity, file.path(out_dir, "08_cross_dataset", "background_sensitivity.tsv"))
} else {
    warning("2010 processed data unavailable: 2010 analysis skipped."); bg_sensitivity <- data.table()
}

# ---- 2010 burden, sensitivity, matching and descriptive figures (B1) ----
burden_results <- data.table(); sensitivity_results <- data.table(); matched_pairs <- data.table()
if (length(res2010)) {
    primary <- res2010[["B1_All_detected_background"]]$decay
    base_covars <- if (sum(is.finite(primary$published_expression)) > 0) "+z_log2_transcript_length+z_log2_published_expression+decay_class" else "+z_log2_transcript_length+decay_class"
    burden_specs <- list(pooled_allele_fraction = list(data = primary[edited == 1], term = "z_gene_editing_burden", source = "gene_editing_burden"), site_number = list(data = primary, term = "z_log1p_editing_site_count", source = "editing_site_count"), mean_site_editing = list(data = primary[edited == 1], term = "z_mean_editing", source = "mean_editing"), median_site_editing = list(data = primary[edited == 1], term = "z_median_editing", source = "median_editing"), max_site_editing = list(data = primary[edited == 1], term = "z_max_editing", source = "max_editing"))
    for (nm in names(burden_specs)) {
        d <- copy(burden_specs[[nm]]$data); src <- burden_specs[[nm]]$source; term <- burden_specs[[nm]]$term
        d[, (term) := zscore(if (src == "editing_site_count") log1p(get(src)) else get(src))]
        f <- as.formula(paste("log2_half_life~", term, base_covars)); vars <- all.vars(f); q <- d[complete.cases(d[, ..vars])]
        if (nrow(q) >= length(vars) + 4 && sd(q[[term]]) > 0) {
            m <- lm(f, data = q); t <- tidy_lm(m, paste0("burden_", nm), "2010_B1"); t <- t[variable == term]; t[, `:=`(burden_variable = src, N_edited = sum(q$edited == 1), N_unedited = sum(q$edited == 0))]; burden_results <- rbind(burden_results, t, fill = TRUE)
        }
    }
    write_tsv(burden_results, file.path(out_dir, "05_2010_regression", "2010_editing_burden_regression.tsv"))

    # Published valid, top-1%-trimmed, and winsorized responses; no value is silently discarded.
    primary_formula <- if (sum(is.finite(primary$published_expression)) > 0) log2_half_life ~ edited + z_log2_transcript_length + z_log2_published_expression + decay_class else log2_half_life ~ edited + z_log2_transcript_length + decay_class
    cut <- quantile(primary$half_life_min, 1 - top_outlier_fraction, na.rm = TRUE)
    sens_sets <- list(A_all_published_valid = copy(primary), B_remove_top_1pct = primary[half_life_min <= cut], C_winsorized = copy(primary))
    sens_sets$C_winsorized[, half_life_min := pmin(half_life_min, cut)]; sens_sets$C_winsorized[, log2_half_life := log2(half_life_min)]
    for (nm in names(sens_sets)) {
        d <- sens_sets[[nm]]; vars <- all.vars(primary_formula); q <- d[complete.cases(d[, ..vars])]; if (nrow(q) > 10 && uniqueN(q$edited) == 2) {
            m <- lm(primary_formula, data = q);t <- tidy_lm(m, nm, "2010_B1")[variable == "edited"];t[, `:=`(sensitivity_analysis = nm, outlier_cutoff_min = cut, N_edited = sum(q$edited == 1), N_unedited = sum(q$edited == 0))];sensitivity_results <- rbind(sensitivity_results, t, fill = TRUE)
        }
    }
    write_tsv(sensitivity_results, file.path(out_dir, "05_2010_regression", "2010_half_life_outlier_sensitivity.tsv"))

    # Within-decay-class nearest-neighbour matching on length and published expression.
    match_vars <- c("z_log2_transcript_length", if (sum(is.finite(primary$published_expression)) > 0) "z_log2_published_expression")
    edrows <- primary[edited == 1 & complete.cases(primary[, ..match_vars])]; controls <- primary[edited == 0 & complete.cases(primary[, ..match_vars])]
    pair_list <- lapply(seq_len(nrow(edrows)), function(i) {
        e <- edrows[i];pool <- controls[decay_class == e$decay_class];if (!nrow(pool)) return(NULL);dist <- sqrt(rowSums((as.matrix(pool[, ..match_vars]) - matrix(as.numeric(e[, ..match_vars]), nrow(pool), length(match_vars), byrow = TRUE))^2));k <- min(max_matching_controls, length(dist));o <- order(dist)[seq_len(k)];data.table(edited_FBgn = e$FBgn, control_FBgn = pool$FBgn[o], decay_class = as.character(e$decay_class), match_rank = seq_len(k), distance = dist[o], edited_half_life = e$half_life_min, control_half_life = pool$half_life_min[o])
    })
    matched_pairs <- rbindlist(pair_list, fill = TRUE)
    write_tsv(matched_pairs, file.path(out_dir, "06_matched_analysis", "matched_pairs.tsv"))
    if (nrow(matched_pairs)) {
        paired <- matched_pairs[, .(edited_half_life = first(edited_half_life), matched_control_mean = mean(control_half_life)), by = edited_FBgn];mw <- wilcox.test(paired$edited_half_life, paired$matched_control_mean, paired = TRUE, exact = FALSE);write_tsv(data.table(N_edited = nrow(paired), controls_per_edited_min = min(matched_pairs[, .N, by = edited_FBgn]$N), controls_per_edited_max = max(matched_pairs[, .N, by = edited_FBgn]$N), median_edited = median(paired$edited_half_life), median_matched_control_mean = median(paired$matched_control_mean), paired_Wilcoxon_P = mw$p.value), file.path(out_dir, "06_matched_analysis", "matched_half_life_comparison.tsv"))
    }

    # Figures A/B/D and descriptive predictor plots.
    xall <- res2010[["B1_All_detected_background"]]$data
    counts <- xall[, .N, by = .(editing_group, decay_label = ifelse(decay_status == 1, "Decay", "Stable"))]
    counts[, fraction := N / sum(N), by = editing_group]
    group_totals <- xall[, .(group_N = .N), by = editing_group]; counts <- merge(counts, group_totals, by = "editing_group");counts[, group_label := paste0(as.character(editing_group), "\nN = ", group_N)]
    fisher_one_p <- res2010[["B1_All_detected_background"]]$fisher$one_sided_greater_P[1]
    pA <- ggplot(counts, aes(group_label, fraction, fill = decay_label)) +
        geom_col() +
        geom_text(aes(label = sprintf("n = %d\n%.1f%%", N, 100 * fraction)), position = position_stack(vjust = .5), size = 3.5) +
        annotate("segment", x = 1, xend = 2, y = 1.055, yend = 1.055, linewidth = .5) +
        annotate("segment", x = 1, xend = 1, y = 1.03, yend = 1.055, linewidth = .5) +
        annotate("segment", x = 2, xend = 2, y = 1.03, yend = 1.055, linewidth = .5) +
        annotate("text", x = 1.5, y = 1.095, label = sprintf("One-sided Fisher exact P = %.3g", fisher_one_p), size = 3.8) +
        scale_y_continuous(labels = function(z) paste0(round(100 * z), "%"), expand = expansion(mult = c(0, .15))) +
        labs(x = NULL, y = "Fraction", fill = NULL, title = "2010 stable versus decay", subtitle = sprintf("B1 exploratory background; Edited: >=%d sites; one-site genes excluded", min_editing_sites_per_gene)) +
        theme_bw()
    save_plot(pA, file.path(out_dir, "09_figures", "Figure_A_2010_stable_decay.pdf"))
    bstats <- primary[, .(n = .N, mean = mean(half_life_min), SD = sd(half_life_min), min_value = min(half_life_min)), by = editing_group];bstats[, `:=`(lower = pmax(mean - SD, min_value), upper = mean + SD, group_label = paste0(as.character(editing_group), "\nn = ", n))];primary <- merge(primary, bstats[, .(editing_group, group_label)], by = "editing_group", all.x = TRUE)
    pB <- ggplot(primary, aes(group_label, half_life_min, fill = editing_group)) +
        geom_violin(trim = FALSE, alpha = .45) +
        geom_boxplot(width = .18, outlier.shape = NA) +
        geom_errorbar(data = bstats, aes(x = group_label, ymin = lower, ymax = upper), inherit.aes = FALSE, width = .10, linewidth = .8) +
        geom_point(data = bstats, aes(x = group_label, y = mean), inherit.aes = FALSE, shape = 23, size = 3.5, fill = "white", color = "black") +
        geom_text(data = bstats, aes(x = group_label, y = upper, label = sprintf("mean <U+00B1> SD\n%.1f <U+00B1> %.1f min", mean, SD)), inherit.aes = FALSE, vjust = -.35, size = 3.2) +
        scale_y_log10(expand = expansion(mult = c(.05, .22))) +
        labs(x = NULL, y = "Published half-life (min, log scale)", title = "2010 decay-gene half-lives", subtitle = "Diamond and vertical line: arithmetic mean <U+00B1> SD; box: median and IQR", caption = NULL) +
        theme_bw() +
        theme(legend.position = "none")
    save_plot(pB, file.path(out_dir, "09_figures", "Figure_B_2010_half_life_distribution.pdf"))

    # Separate violin figure for the representative transcript-length distribution.
    # This does not alter Figure C or add transcript length a second time to any model.
    length_plot_data <- primary[
        is.finite(transcript_length) & transcript_length > 0,
        .(FBgn, editing_group, transcript_length)
    ]
    length_wilcox <- wilcox.test(transcript_length ~ editing_group,
        data = length_plot_data, exact = FALSE
    )
    length_stats <- length_plot_data[, .(
        n = .N,
        mean = mean(transcript_length),
        SD = sd(transcript_length),
        median = median(transcript_length),
        IQR = IQR(transcript_length),
        min_value = min(transcript_length),
        max_value = max(transcript_length)
    ), by = editing_group]
    length_stats[, `:=`(
        lower = pmax(mean - SD, min_value),
        upper = mean + SD,
        length_group_label = paste0(as.character(editing_group), "\nn = ", n),
        Wilcoxon_P = length_wilcox$p.value
    )]
    length_stats[, length_group_label := factor(
        length_group_label,
        # Keep all Edited-vs-control violins consistent: Edited on the left.
        levels = length_group_label[order(as.integer(editing_group), decreasing = TRUE)]
    )]
    length_plot_data <- merge(
        length_plot_data,
        length_stats[, .(editing_group, length_group_label)],
        by = "editing_group", all.x = TRUE
    )
    length_plot_data[, length_group_label := factor(
        length_group_label,
        levels = levels(length_stats$length_group_label)
    )]
    write_tsv(length_stats[, .(
        editing_group, n,
        mean_transcript_length_nt = mean,
        SD_transcript_length_nt = SD, median_transcript_length_nt = median,
        IQR_transcript_length_nt = IQR, min_transcript_length_nt = min_value,
        max_transcript_length_nt = max_value, Wilcoxon_P
    )], file.path(out_dir, "02_gene_features", "2010_transcript_length_edited_vs_unedited.tsv"))

    p_length <- ggplot(
        length_plot_data,
        aes(length_group_label, transcript_length, fill = editing_group)
    ) +
        geom_violin(trim = FALSE, alpha = .45) +
        geom_boxplot(width = .18, outlier.shape = NA) +
        geom_errorbar(
            data = length_stats,
            aes(x = length_group_label, ymin = lower, ymax = upper),
            inherit.aes = FALSE, width = .10, linewidth = .8
        ) +
        geom_point(
            data = length_stats, aes(x = length_group_label, y = mean),
            inherit.aes = FALSE, shape = 23, size = 3.5,
            fill = "white", color = "black"
        ) +
        geom_text(
            data = length_stats,
            aes(
                x = length_group_label, y = upper,
                label = sprintf("mean <U+00B1> SD\n%.0f <U+00B1> %.0f nt", mean, SD)
            ),
            inherit.aes = FALSE, vjust = -.35, size = 3.2
        ) +
        scale_y_log10(expand = expansion(mult = c(.05, .23))) +
        labs(
            x = NULL, y = "Representative transcript length (nt, log scale)",
            title = "Transcript-length distribution of 2010 decay genes",
            subtitle = sprintf("Edited versus No_detected_editing; Wilcoxon P = %.3g", length_wilcox$p.value),
            caption = NULL
        ) +
        theme_bw() +
        theme(legend.position = "none", plot.caption = element_text(hjust = 0, size = 8.5))
    save_plot(p_length,
        file.path(out_dir, "09_figures", "Figure_H_2010_transcript_length_distribution.pdf"),
        width = 7.5, height = 5.8
    )

    burden_plot_data <- primary[
        edited == 1 & is.finite(gene_editing_burden) & is.finite(log2_half_life)
    ]
    if (nrow(burden_plot_data)) {
        burden_cor <- cor.test(
            burden_plot_data$gene_editing_burden,
            burden_plot_data$log2_half_life,
            method = "spearman", exact = FALSE
        )
        burden_rho <- unname(burden_cor$estimate)
        burden_cor_p <- burden_cor$p.value
        burden_cor_label <- sprintf(
            "Spearman rho = %.3f; P = %.3g", burden_rho, burden_cor_p
        )
        pD <- ggplot(burden_plot_data, aes(gene_editing_burden, log2_half_life)) +
            geom_point(alpha = .65, color = "#2C7FB8") +
            geom_smooth(method = "lm", se = TRUE, color = "#D95F0E") +
            scale_x_continuous(labels = function(z) sprintf("%.0f%%", 100 * z)) +
            labs(
                x = "Pooled gene editing burden (sum altcount / sum depth)",
                y = "log2 half-life",
                title = "Pooled editing allele fraction and published half-life",
                subtitle = sprintf(
                    "Edited genes (>=%d detected sites): n = %d; %s",
                    min_editing_sites_per_gene, nrow(burden_plot_data), burden_cor_label
                ),
                caption = NULL
            ) +
            theme_bw() +
            theme(plot.caption = element_text(hjust = 0, size = 8.5, margin = margin(t = 8)))
        save_plot(
            pD, file.path(out_dir, "09_figures", "Figure_D_2010_editing_burden.pdf"),
            8, 5.5
        )
    }
    pred_map <- c(transcript_length = "Transcript length", CDS_length = "CDS length", UTR3_length = "3'UTR length", published_expression = "Published baseline expression", editing_site_count = "Editing site count", mean_editing = "Mean editing level")
    descriptive_results <- data.table()
    for (v in names(pred_map)) {
        d <- primary[is.finite(get(v)) & is.finite(log2_half_life)];if (nrow(d) >= 3) {
            ct <- cor.test(d[[v]], d$log2_half_life, method = "spearman", exact = FALSE);descriptive_results <- rbind(descriptive_results, data.table(predictor = v, N = nrow(d), Spearman_rho = unname(ct$estimate), P = ct$p.value));p <- ggplot(d, aes(x = .data[[v]], y = log2_half_life)) +
                geom_point(alpha = .5) +
                geom_smooth(method = "lm", se = TRUE) +
                labs(x = pred_map[[v]], y = "log2 half-life", subtitle = sprintf("Spearman rho = %.3f; P = %.3g; n = %d", unname(ct$estimate), ct$p.value, nrow(d))) +
                theme_bw();save_plot(p, file.path(out_dir, "09_figures", paste0("descriptive_", v, ".pdf")))
        }
    }
    write_tsv(descriptive_results, file.path(out_dir, "04_2010_half_life", "2010_descriptive_correlations.tsv"))
}

# ---- 2023 PLOS independent validation and matched expression from S1 ----
parse_2023 <- function(s2_file, s1_file) {
    if (!file.exists(s2_file)) return(NULL)
    s2 <- read_xlsx_sheet_base(s2_file, 1)
    gene_col <- find_col(s2, c("gene_id", "FBgn", "gene")); tx_col <- find_col(s2, c("tr_id", "transcript_id", "target_id")); hl_col <- find_col(s2, c("Half_life (mins)", "half_life_min", "half-life")); cat_col <- find_col(s2, c("category", "class"), FALSE)
    half <- data.table(FBgn = clean_id(s2[[gene_col]]), transcript_id = clean_id(s2[[tx_col]]), half_life_min = suppressWarnings(as.numeric(s2[[hl_col]])), half_life_category = if (!is.na(cat_col)) as.character(s2[[cat_col]]) else NA_character_)
    half <- half[is.finite(half_life_min) & half_life_min > 0]
    expr_cache <- file.path(out_dir, "00_input_QC", "2023_published_expression.tsv")
    expr_out <- data.table(transcript_id = character(), initial_TPM = numeric(), max_TPM = numeric(), mean_TPM = numeric(), median_TPM = numeric())
    if (file.exists(expr_cache)) {
        expr_out <- fread(expr_cache)
    } else if (file.exists(s1_file)) {
        message("Reading 2023 S1 transcript TPM worksheet; first run may take about two minutes.")
        e <- read_xlsx_sheet_base(s1_file, 2); idc <- find_col(e, c("target_id", "transcript_id", "tr_id")); valcols <- setdiff(names(e), idc)
        for (v in valcols) set(e, j = v, value = suppressWarnings(as.numeric(e[[v]])))
        initial_cols <- grep("^95[_-]", valcols, value = TRUE)
        m <- as.matrix(e[, ..valcols]); initm <- if (length(initial_cols)) as.matrix(e[, ..initial_cols]) else m[, 1:min(3, ncol(m)), drop = FALSE]
        expr_out <- data.table(transcript_id = clean_id(e[[idc]]), initial_TPM = rowMeans(initm, na.rm = TRUE), max_TPM = apply(m, 1, max, na.rm = TRUE), mean_TPM = rowMeans(m, na.rm = TRUE), median_TPM = apply(m, 1, median, na.rm = TRUE))
        for (v in c("initial_TPM", "max_TPM", "mean_TPM", "median_TPM")) expr_out[!is.finite(get(v)), (v) := NA_real_]
    } else warning("2023 S1 Table unavailable: expression models will be skipped and no other-tissue substitute will be used.")
    write_tsv(expr_out, file.path(out_dir, "00_input_QC", "2023_published_expression.tsv"))
    tx <- merge(half, expr_out, by = "transcript_id", all.x = TRUE); tx <- merge(tx, transcript_features, by = c("transcript_id", "FBgn"), all.x = TRUE)
    write_tsv(tx, file.path(out_dir, "07_2023_validation", "2023_transcript_half_life_merged.tsv"))
    gene <- tx[, .(gene_symbol = gene_symbol[which(!is.na(gene_symbol))[1]], n_transcripts = .N, min_half_life = as.numeric(min(half_life_min)), half_life_min = as.numeric(median(half_life_min)), max_half_life = as.numeric(max(half_life_min)), transcript_length = as.numeric(median(transcript_length, na.rm = TRUE)), CDS_length = as.numeric(median(CDS_length, na.rm = TRUE)), UTR5_length = as.numeric(median(UTR5_length, na.rm = TRUE)), UTR3_length = as.numeric(median(UTR3_length, na.rm = TRUE)), exon_number = as.numeric(median(exon_number, na.rm = TRUE)), GC_content = as.numeric(median(GC_content, na.rm = TRUE)), A_content = as.numeric(median(A_content, na.rm = TRUE)), initial_TPM = as.numeric(median(initial_TPM, na.rm = TRUE)), max_TPM = as.numeric(median(max_TPM, na.rm = TRUE)), mean_TPM = as.numeric(median(mean_TPM, na.rm = TRUE)), median_TPM = as.numeric(median(median_TPM, na.rm = TRUE))), by = FBgn]
    for (v in c("transcript_length", "CDS_length", "UTR5_length", "UTR3_length", "exon_number", "GC_content", "A_content", "initial_TPM", "max_TPM", "mean_TPM", "median_TPM")) gene[!is.finite(get(v)), (v) := NA_real_]
    gene[, published_expression := mean_TPM];gene[, log2_half_life := log2(half_life_min)];gene <- add_editing(gene);gene <- standardize_model_vars(gene)
    write_tsv(
        gene[is.na(edited), .(FBgn, gene_symbol, editing_site_count, editing_binary_status, half_life_min)],
        file.path(out_dir, "00_input_QC", "2023_editing_threshold_excluded_genes.tsv")
    )
    gene
}

run_2023_background <- function(all2023, bg, full_outputs = FALSE) {
    x <- apply_background(all2023, bg);ed <- x[edited == 1, half_life_min];co <- x[edited == 0, half_life_min];wt <- if (length(ed) && length(co)) wilcox.test(ed, co, exact = FALSE) else NULL
    simple <- x[, .(n = .N, median = median(half_life_min), mean = mean(half_life_min), SD = sd(half_life_min), IQR = IQR(half_life_min), min = min(half_life_min), max = max(half_life_min)), by = editing_group];simple[, `:=`(background_definition = bg, Wilcoxon_P = if (is.null(wt)) NA_real_ else wt$p.value)]
    specs <- list(V1_primary_length = log2_half_life ~ edited + z_log2_transcript_length)
    if (sum(is.finite(x$mean_TPM)) > 0) specs$V2_secondary_mean_TPM <- log2_half_life~edited + z_log2_transcript_length + z_log2_published_expression
    # Sensitivity definitions use the same transcript/gene time series, never editing-experiment expression.
    for (v in c("initial_TPM", "max_TPM", "median_TPM")) {
        if (sum(is.finite(x[[v]])) > 0) {
            zv <- paste0("z_log2_", v);x[, (zv) := zscore(safe_log2(get(v), 1))];specs[[paste0("V_expression_sensitivity_", v)]] <- as.formula(paste("log2_half_life~edited+z_log2_transcript_length+", zv))
        }
    }
    results <- list();models <- list();complete <- list()
    for (nm in names(specs)) {
        f <- specs[[nm]]; vars <- all.vars(f)
        complete[[nm]] <- complete_case_report(x, vars, bg, nm)
        d <- x[complete.cases(x[, ..vars])]
        if (nrow(d) >= max(8, length(vars) + 3) && uniqueN(d$edited) == 2) {
            m <- lm(f, data = d); t <- tidy_lm(m, nm, bg)
            expr_def <- if (grepl("initial", nm)) "initial_TPM" else if (grepl("max", nm)) "max_TPM" else if (grepl("median", nm)) "median_TPM" else if (grepl("mean", nm)) "mean_TPM" else "none"
            t[, `:=`(N_edited = sum(d$edited == 1), N_unedited = sum(d$edited == 0), expression_definition = expr_def)]
            results[[nm]] <- t; models[[nm]] <- m
        } else results[[nm]] <- data.table(dataset = bg, model = nm, variable = "MODEL_SKIPPED", N_total = nrow(d), N_edited = sum(d$edited == 1), N_unedited = sum(d$edited == 0), note = "Insufficient complete cases or only one editing group")
    }
    tab <- rbindlist(results, fill = TRUE);tab[, BH_adjusted_P := p.adjust(ordinary_P, "BH")]
    if (full_outputs) {
        write_tsv(simple, file.path(out_dir, "07_2023_validation", "2023_half_life_analysis.tsv"));write_tsv(tab, file.path(out_dir, "07_2023_validation", "2023_regression.tsv"));write_tsv(rbindlist(complete, fill = TRUE), file.path(out_dir, "07_2023_validation", "2023_model_complete_cases.tsv"));diag <- rbindlist(lapply(names(models), function(nm) diagnose_lm(models[[nm]], nm, bg, file.path(out_dir, "07_2023_validation", paste0(nm, "_diagnostics.pdf")))), fill = TRUE);write_tsv(diag, file.path(out_dir, "07_2023_validation", "2023_regression_diagnostics.tsv"))
    }
    preferred <- if ("V2_secondary_mean_TPM" %in% names(models)) "V2_secondary_mean_TPM" else "V1_primary_length";pe <- tab[model == preferred & variable == "edited"]
    sens <- data.table(background_definition = bg, N_edited = length(ed), N_control = length(co), median_half_life_edited = if (length(ed)) median(ed) else NA_real_, median_half_life_control = if (length(co)) median(co) else NA_real_, difference = if (length(ed) && length(co)) median(ed) - median(co) else NA_real_, Wilcoxon_P = if (is.null(wt)) NA_real_ else wt$p.value, regression_model = preferred, editing_beta = pe$estimate[1], CI_low = pe$CI_low95[1], CI_high = pe$CI_high95[1], regression_P = pe$ordinary_P[1], robust_P = pe$robust_P_HC3[1])
    list(data = x, simple = simple, regression = tab, models = models, background = sens)
}

data2023 <- parse_2023(half_life_2023_file, expression_2023_file);res2023 <- list()
if (!is.null(data2023)) {
    for (bg in background_names) res2023[[bg]] <- run_2023_background(data2023, bg, startsWith(bg, "B1"))
    write_tsv(data2023, file.path(out_dir, "07_2023_validation", "2023_gene_half_life_merged.tsv"))
    bg23 <- rbindlist(lapply(res2023, `[[`, "background"), fill = TRUE);bg23[, dataset := "2023"]
    if (nrow(bg_sensitivity)) {
        bg_sensitivity[, dataset := "2010"];allbg <- rbind(bg_sensitivity, bg23, fill = TRUE)
    } else allbg <- bg23
    write_tsv(allbg, file.path(out_dir, "08_cross_dataset", "background_sensitivity.tsv"))
    eplot <- copy(res2023[[1]]$data);estats <- eplot[, .(n = .N, mean = mean(half_life_min), SD = sd(half_life_min), min_value = min(half_life_min)), by = editing_group];estats[, `:=`(lower = pmax(mean - SD, min_value), upper = mean + SD, group_label = paste0(as.character(editing_group), "\nn = ", n))];eplot <- merge(eplot, estats[, .(editing_group, group_label)], by = "editing_group", all.x = TRUE)
    pE <- ggplot(eplot, aes(group_label, half_life_min, fill = editing_group)) +
        geom_violin(trim = FALSE, alpha = .45) +
        geom_boxplot(width = .18, outlier.shape = NA) +
        geom_jitter(width = .12, size = .8, alpha = .5) +
        geom_errorbar(data = estats, aes(x = group_label, ymin = lower, ymax = upper), inherit.aes = FALSE, width = .10, linewidth = .8) +
        geom_point(data = estats, aes(x = group_label, y = mean), inherit.aes = FALSE, shape = 23, size = 3.5, fill = "white", color = "black") +
        geom_text(data = estats, aes(x = group_label, y = upper, label = sprintf("mean <U+00B1> SD\n%.1f <U+00B1> %.1f min", mean, SD)), inherit.aes = FALSE, vjust = -.35, size = 3.2) +
        scale_y_log10(expand = expansion(mult = c(.05, .22))) +
        labs(x = NULL, y = "Published half-life (min, log scale)", title = "2023 independent validation", subtitle = "Diamond and vertical line: arithmetic mean <U+00B1> SD; box: median and IQR", caption = NULL) +
        theme_bw() +
        theme(legend.position = "none");save_plot(pE, file.path(out_dir, "09_figures", "Figure_E_2023_validation.pdf"))
} else {
    warning("2023 S2 Table unavailable: independent validation skipped.");allbg <- bg_sensitivity
}

# ---- Cross-dataset half-life distributions and overlapping genes ----
# The two studies are displayed side-by-side but are not numerically pooled.
# Black points identify genes with a finite gene-level half-life in both studies.
if (!is.null(data2010) && !is.null(data2023)) {
    hl2010 <- unique(data2010[
        decay_class %chin% c("II", "III", "IV", "V") & is.finite(half_life_min),
        .(FBgn, gene_symbol, decay_class = as.character(decay_class), half_life_2010_min = half_life_min)
    ])
    hl2023 <- unique(data2023[
        is.finite(half_life_min),
        .(FBgn, gene_symbol_2023 = gene_symbol, half_life_2023_min = half_life_min)
    ])
    overlap_hl <- merge(hl2010, hl2023, by = "FBgn", all = FALSE)
    overlap_hl[, gene_symbol := fifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_symbol_2023)]
    overlap_hl[, `:=`(
        log2_half_life_2010 = log2(half_life_2010_min),
        log2_half_life_2023 = log2(half_life_2023_min),
        ratio_2010_to_2023 = half_life_2010_min / half_life_2023_min,
        within_2fold = pmax(
            half_life_2010_min / half_life_2023_min,
            half_life_2023_min / half_life_2010_min
        ) <= 2
    )]
    setcolorder(overlap_hl, c(
        "FBgn", "gene_symbol", "decay_class", "half_life_2010_min",
        "half_life_2023_min", "log2_half_life_2010", "log2_half_life_2023",
        "ratio_2010_to_2023", "within_2fold"
    ))
    write_tsv(
        overlap_hl[, !"gene_symbol_2023"],
        file.path(out_dir, "08_cross_dataset", "2010_2023_overlapping_gene_half_lives.tsv")
    )

    spearman_test <- if (nrow(overlap_hl) >= 3) cor.test(overlap_hl$half_life_2010_min,
        overlap_hl$half_life_2023_min,
        method = "spearman", exact = FALSE
    ) else NULL
    pearson_test <- if (nrow(overlap_hl) >= 3) cor.test(overlap_hl$log2_half_life_2010,
        overlap_hl$log2_half_life_2023,
        method = "pearson"
    ) else NULL
    overlap_summary <- data.table(
        N_overlap = nrow(overlap_hl),
        median_2010_min = median(overlap_hl$half_life_2010_min),
        median_2023_min = median(overlap_hl$half_life_2023_min),
        Spearman_rho = if (is.null(spearman_test)) NA_real_ else unname(spearman_test$estimate),
        Spearman_P = if (is.null(spearman_test)) NA_real_ else spearman_test$p.value,
        Pearson_log2_r = if (is.null(pearson_test)) NA_real_ else unname(pearson_test$estimate),
        Pearson_log2_P = if (is.null(pearson_test)) NA_real_ else pearson_test$p.value,
        N_within_2fold = sum(overlap_hl$within_2fold),
        proportion_within_2fold = mean(overlap_hl$within_2fold)
    )
    write_tsv(
        overlap_summary,
        file.path(out_dir, "08_cross_dataset", "2010_2023_half_life_consistency.tsv")
    )

    violin_data <- rbind(
        hl2010[, .(FBgn, dataset = "2010 Genome Biology", half_life_min = half_life_2010_min)],
        hl2023[, .(FBgn, dataset = "2023 PLOS Biology", half_life_min = half_life_2023_min)]
    )
    violin_data[, overlap_gene := FBgn %chin% overlap_hl$FBgn]
    violin_stats <- violin_data[, .(
        n = .N,
        overlap_n = sum(overlap_gene),
        mean = mean(half_life_min),
        SD = sd(half_life_min),
        min_value = min(half_life_min)
    ), by = dataset]
    violin_stats[, `:=`(
        lower = pmax(mean - SD, min_value),
        upper = mean + SD,
        group_label = sprintf("%s\nTotal n = %d; overlap n = %d", dataset, n, overlap_n)
    )]
    violin_data <- merge(violin_data, violin_stats[, .(dataset, group_label)], by = "dataset", all.x = TRUE)
    violin_data[, group_label := factor(group_label, levels = violin_stats$group_label)]
    violin_stats[, group_label := factor(group_label, levels = violin_stats$group_label)]

    p_overlap_violin <- ggplot(violin_data, aes(group_label, half_life_min, fill = dataset)) +
        geom_violin(trim = FALSE, alpha = 0.48, color = "grey35") +
        geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.70) +
        geom_point(
            data = violin_data[overlap_gene == TRUE],
            position = position_jitter(width = 0.075, height = 0, seed = 20260827),
            shape = 21, size = 2.2, stroke = 0.25, fill = "black", color = "black"
        ) +
        geom_errorbar(
            data = violin_stats,
            aes(x = group_label, ymin = lower, ymax = upper),
            inherit.aes = FALSE, width = 0.10, linewidth = 0.8
        ) +
        geom_point(
            data = violin_stats, aes(x = group_label, y = mean),
            inherit.aes = FALSE, shape = 23, size = 3.5, fill = "white", color = "black"
        ) +
        geom_text(
            data = violin_stats,
            aes(
                x = group_label, y = upper,
                label = sprintf("mean +/- SD\n%.1f +/- %.1f min", mean, SD)
            ),
            inherit.aes = FALSE, vjust = -0.35, size = 3.2
        ) +
        scale_y_log10(expand = expansion(mult = c(0.05, 0.23))) +
        scale_fill_manual(values = c("2010 Genome Biology" = "#5B8DB8", "2023 PLOS Biology" = "#E69F5B")) +
        labs(
            x = NULL, y = "Published half-life (min, log scale)",
            title = "Published half-life distributions in the 2010 and 2023 datasets",
            subtitle = sprintf("Black points: the same %d genes measured in both datasets; diamond and line: arithmetic mean +/- SD", nrow(overlap_hl)),
            caption = NULL
        ) +
        theme_bw() +
        theme(legend.position = "none", plot.caption = element_text(hjust = 0, size = 8.5))
    save_plot(p_overlap_violin,
        file.path(out_dir, "09_figures", "Figure_G_2010_2023_overlap_half_life_violin.pdf"),
        width = 8.5, height = 6
    )
}

# ---- Cross-model summaries and forest plots ----
model_summary <- data.table();cross_effect <- data.table()
if (length(res2010)) {
    r <- res2010[[1]]$regression;model_summary <- rbind(model_summary, r[variable == "edited", .(dataset = "2010", model, N = N_total, N_edited, N_unedited, editing_beta = estimate, CI_low = CI_low95, CI_high = CI_high95, P = ordinary_P, robust_P = robust_P_HC3, adjusted_R2, AIC)], fill = TRUE)
    pick <- if (any(r$model == "H3")) "H3" else "H1";cross_effect <- rbind(cross_effect, r[model == pick & variable == "edited", .(dataset = "2010 Genome Biology", model, effect = estimate, CI_low = CI_low95, CI_high = CI_high95, P = ordinary_P, robust_P = robust_P_HC3, N = N_total, N_edited, N_control = N_unedited)], fill = TRUE)
    forest <- r[model %chin% c("H3", "H4_component_length_model", "H5_extended_model") & variable != "(Intercept)" & is.finite(estimate)]
    if (nrow(forest)) {
        label_map <- c(edited = "Editing: Edited vs No detected", z_log2_transcript_length = "Transcript length (per 1 SD)", z_log2_published_expression = "Published expression (per 1 SD)", z_log2_CDS_length = "CDS length (per 1 SD)", z_log2_UTR5_length = "5'UTR length (per 1 SD)", z_log2_UTR3_length = "3'UTR length (per 1 SD)", z_GC_content = "GC content (per 1 SD)", z_A_content = "A content (per 1 SD)", decay_classIII = "Decay class III vs II", decay_classIV = "Decay class IV vs II", decay_classV = "Decay class V vs II")
        forest[, label := unname(label_map[variable])];forest[is.na(label), label := variable]
        forest[, model_label := factor(model, levels = c("H3", "H4_component_length_model", "H5_extended_model"), labels = c("H3: total-length", "H4: component-length", "H5: sequence-extended"))]
        model_n <- unique(forest[, .(N_total, N_edited, N_unedited)])[1]
        pC <- ggplot(forest, aes(estimate, reorder(label, estimate), xmin = CI_low95, xmax = CI_high95, color = model_label)) +
            geom_vline(xintercept = 0, linetype = 2) +
            geom_pointrange(position = position_dodge(width = .5)) +
            labs(x = "Regression coefficient beta (continuous covariates standardized)", y = NULL, color = "Model", title = "2010 multivariable coefficients", subtitle = sprintf("Each model: N = %d (Edited n = %d; No_detected_editing n = %d)", model_n$N_total, model_n$N_edited, model_n$N_unedited), caption = NULL) +
            theme_bw();save_plot(pC, file.path(out_dir, "09_figures", "Figure_C_multivariable_coefficients.pdf"), 9, 7)
    }
}
if (length(res2023)) {
    r <- res2023[[1]]$regression;model_summary <- rbind(model_summary, r[variable == "edited", .(dataset = "2023", model, N = N_total, N_edited, N_unedited, editing_beta = estimate, CI_low = CI_low95, CI_high = CI_high95, P = ordinary_P, robust_P = robust_P_HC3, adjusted_R2, AIC)], fill = TRUE)
    pick <- if (any(r$model == "V2_secondary_mean_TPM")) "V2_secondary_mean_TPM" else "V1_primary_length";cross_effect <- rbind(cross_effect, r[model == pick & variable == "edited", .(dataset = "2023 PLOS Biology", model, effect = estimate, CI_low = CI_low95, CI_high = CI_high95, P = ordinary_P, robust_P = robust_P_HC3, N = N_total, N_edited, N_control = N_unedited)], fill = TRUE)
}
write_tsv(model_summary, file.path(out_dir, "10_tables", "model_summary.tsv"));write_tsv(cross_effect, file.path(out_dir, "08_cross_dataset", "cross_dataset_effect_summary.tsv"))
if (nrow(cross_effect)) {
    cross_effect[, plot_label := sprintf("%s\nN=%d; Edited n=%d; control n=%d", dataset, N, N_edited, N_control)];pF <- ggplot(cross_effect, aes(effect, plot_label, xmin = CI_low, xmax = CI_high)) +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_pointrange() +
        labs(x = "Editing coefficient for log2 half-life", y = NULL, title = "Dataset-level editing associations", caption = NULL) +
        theme_bw();save_plot(pF, file.path(out_dir, "09_figures", "Figure_F_cross_dataset_forest.pdf"), 8, 4.5)
}

# ---- One-page, logically ordered A4 composite ----
# Logical order: raw outcomes (A/B), major confounder (C), adjusted model (D),
# editing burden (E), cross-study comparability (F), validation (G), and the
# dataset-level adjusted-effect summary (H). Patchwork aligns plotting panels
# within the shared design grid; D spans two columns because it has many terms.
if (requireNamespace("patchwork", quietly = TRUE) &&
    all(vapply(
        c("pA", "pB", "p_length", "pC", "pD", "p_overlap_violin", "pE", "pF"),
        function(nm) exists(nm, envir = .GlobalEnv, inherits = FALSE), logical(1)
    ))) {
    composite_theme <- theme_minimal(base_size = 8) + theme(
        text = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, face = "bold", colour = "black", margin = margin(b = 2)),
        plot.subtitle = element_text(size = 8, colour = "black", margin = margin(b = 2)),
        axis.title = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.key.height = grid::unit(3.0, "mm"),
        legend.key.width = grid::unit(4.0, "mm"),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(colour = "black", linewidth = 0.35),
        axis.line.y.left = element_line(colour = "black", linewidth = 0.35),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        plot.margin = margin(4, 5, 4, 5)
    )
    # geom_text sizes are specified independently of theme text. Force every
    # annotation layer to the same 8-point size (8 pt / 2.845 pt/mm).
    force_annotation_text_8pt <- function(p) {
        for (i in seq_along(p$layers)) {
            if (inherits(p$layers[[i]]$geom, "GeomText") ||
                inherits(p$layers[[i]]$geom, "GeomLabel")) {
                p$layers[[i]]$aes_params$size <- 8 / 2.845276
            }
        }
        p
    }
    force_opaque_vector_layers <- function(p) {
        for (i in seq_along(p$layers)) {
            # Explicit alpha=1 prevents PDF transparency ExtGState objects. All
            # shapes remain ordinary editable vector paths rather than soft masks.
            p$layers[[i]]$aes_params$alpha <- 1
        }
        p
    }
    compact_plot <- function(p, title, subtitle = NULL, legend = "none") {
        q <- p + labs(title = title, subtitle = subtitle, caption = NULL) +
            composite_theme + theme(legend.position = legend) +
            guides(
                colour = guide_legend(nrow = 1, byrow = TRUE),
                fill = guide_legend(nrow = 1, byrow = TRUE)
            )
        force_opaque_vector_layers(force_annotation_text_8pt(q))
    }
    editing_palette <- c("No_detected_editing" = "#4C78A8", "Edited" = "#E45756")
    dataset_palette <- c("2010 Genome Biology" = "#4C78A8", "2023 PLOS Biology" = "#F2A65A")
    comp_A <- compact_plot(pA, "A  Stable versus decay", legend = "bottom") +
        scale_fill_manual(values = c("Decay" = "#E45756", "Stable" = "#4C78A8"))
    comp_B <- compact_plot(pB, "B  2010 half-life distribution") +
        scale_fill_manual(values = editing_palette) +
        labs(y = "Half-life (min, log scale)")
    comp_C <- compact_plot(p_length, "C  Representative transcript length") +
        scale_fill_manual(values = editing_palette) +
        labs(y = "Transcript length (nt, log scale)")
    comp_D <- compact_plot(
        pC, "D  Multivariable half-life models",
        sprintf("N=3,566 (Edited=63; control=3,503); Edited >=%d sites", min_editing_sites_per_gene),
        legend = "bottom"
    ) + scale_color_manual(values = c("#E45756", "#54A24B", "#4C78A8")) +
        labs(x = "Regression coefficient beta")
    comp_E <- compact_plot(
        pD, "E  Editing burden",
        sprintf(
            "n=%d; Spearman rho=%.3f\nP=%.3g",
            nrow(burden_plot_data), burden_rho, burden_cor_p
        )
    ) + labs(x = "Pooled editing allele fraction", y = "log2 half-life")
    comp_E$layers[[1]]$aes_params$colour <- "#E45756"
    comp_E$layers[[2]]$aes_params$colour <- "#A63D40"
    # With an opaque HC ribbon, draw the smooth/ribbon first and points second so
    # no observations are covered in the editable PDF.
    comp_E$layers <- comp_E$layers[c(2, 1)]
    comp_F <- compact_plot(
        p_overlap_violin, "F  2010 versus 2023 half-life",
        sprintf("Black points mark the same %d genes in both datasets", nrow(overlap_hl))
    ) + scale_fill_manual(values = dataset_palette) +
        scale_x_discrete(labels = function(z) ifelse(
            grepl("^2010", z), "2010\nn=3,637; overlap=26", "2023\nn=186; overlap=26"
        )) + labs(y = "Half-life (min, log scale)")
    comp_G <- compact_plot(pE, "G  2023 independent validation") +
        scale_fill_manual(values = editing_palette) +
        labs(y = "Half-life (min, log scale)")
    comp_H <- compact_plot(pF, "H  Adjusted editing effects") +
        labs(x = "Editing coefficient (log2 half-life)")

    a4_design <- "
ABC
DDE
FGH
"
    a4_composite <- patchwork::wrap_plots(
        A = comp_A, B = comp_B, C = comp_C, D = comp_D,
        E = comp_E, F = comp_F, G = comp_G, H = comp_H,
        design = a4_design, guides = "keep"
    ) + patchwork::plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
    ggsave(
        file.path(out_dir, "09_figures", "Figure_A-H_A4_composite.pdf"),
        a4_composite,
        width = 11.69, height = 8.27, units = "in",
        device = "pdf", limitsize = FALSE
    )
    write_tsv(data.table(
        panel = LETTERS[1:8],
        source_figure = c(
            "Figure_A_2010_stable_decay.pdf",
            "Figure_B_2010_half_life_distribution.pdf",
            "Figure_H_2010_transcript_length_distribution.pdf",
            "Figure_C_multivariable_coefficients.pdf",
            "Figure_D_2010_editing_burden.pdf",
            "Figure_G_2010_2023_overlap_half_life_violin.pdf",
            "Figure_E_2023_validation.pdf",
            "Figure_F_cross_dataset_forest.pdf"
        ),
        panel_topic = c(
            "Stable versus decay", "2010 half-life distribution",
            "Representative transcript length", "Multivariable half-life models",
            "Editing burden", "2010 versus 2023 half-life distributions",
            "2023 independent validation", "Adjusted editing effects across datasets"
        )
    ), file.path(out_dir, "09_figures", "Figure_A-H_A4_panel_order.tsv"))
} else {
    warning("A4 composite skipped: patchwork or one or more component plots are unavailable.")
}

if (nrow(allbg)) {
    allbg[, plot_label := sprintf("%s | %s\nEdited n=%d; control n=%d", dataset, background_definition, N_edited, N_control)];pbg <- ggplot(allbg[is.finite(editing_beta)], aes(editing_beta, plot_label, xmin = CI_low, xmax = CI_high, color = dataset)) +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_pointrange() +
        labs(x = "Editing coefficient", y = NULL, title = "Background-definition sensitivity", caption = NULL) +
        theme_bw();save_plot(pbg, file.path(out_dir, "09_figures", "background_sensitivity_forest.pdf"), 10, 5)
}

# Plain-language, non-causal summary. Numerical conclusions are only populated
# from actually fitted models and every statement preserves condition mismatch.
fmt_effect <- function(x) {
    if (!nrow(x) || !is.finite(x$estimate[1])) "not estimable" else sprintf("beta = %.3f (95%% CI %.3f to %.3f; ordinary P = %.3g; HC3 P = %.3g; N = %d)", x$estimate[1], x$CI_low95[1], x$CI_high95[1], x$ordinary_P[1], x$robust_P_HC3[1], x$N_total[1])
}
sum_lines <- c("Scope and interpretation:", "Genes with detected A-to-I RNA editing in our dataset were tested for association with published embryonic mRNA half-lives.", "The editing and half-life measurements may differ in developmental stage, tissue, genotype, and experimental condition. The analysis does not establish that these transcripts were edited in the embryos used for half-life measurement and cannot establish that editing caused decay.", "", "Background:", paste0("B1 is exploratory: Edited requires >=", min_editing_sites_per_gene, " detected sites; No_detected_editing requires zero sites; genes below the threshold are explicitly excluded. B2 available: ", nrow(b2_genes) > 0, ". B3 available: ", nrow(b3_genes) > 0, "."))
if (length(res2010)) {
    r <- res2010[[1]];f <- r$fisher[1];h2 <- r$regression[model == "H2" & variable == "edited"];h3 <- r$regression[model == "H3" & variable == "edited"];h4 <- r$regression[model == "H4_component_length_model" & variable == "edited"]
    l2 <- r$logistic[model == "L2" & variable == "edited"]
    sm_ed <- r$simple[as.character(editing_group) == "Edited"];sm_co <- r$simple[as.character(editing_group) != "Edited"]
    allele_burden <- burden_results[burden_variable == "gene_editing_burden"];site_burden <- burden_results[burden_variable == "editing_site_count"]
    sum_lines <- c(sum_lines, "", "2010 results:", sprintf("Q1 stable/decay Fisher: OR = %.3f, 95%% CI %.3f to %.3f, two-sided P = %.3g, one-sided P (Edited has higher Decay odds) = %.3g. Published-expression-adjusted L2 editing OR = %.3f (95%% CI %.3f to %.3f; P = %.3g).", f$odds_ratio, f$OR_low95, f$OR_high95, f$two_sided_P, f$one_sided_greater_P, l2$OR[1], l2$OR_low95[1], l2$OR_high95[1], l2$P[1]), sprintf("Q2 decay-gene medians: Edited %.2f min versus No_detected_editing %.2f min; Wilcoxon P = %.3g.", sm_ed$median[1], sm_co$median[1], sm_ed$Wilcoxon_P[1]), paste0("Q3 H2 (published expression): ", fmt_effect(h2)), paste0("Q4/Q6 H3 (published expression and decay class): ", fmt_effect(h3)), paste0("Q5 H4 component-length model editing term: ", fmt_effect(h4)), paste0("Q7 pooled allele-fraction burden: ", fmt_effect(allele_burden), "; site-count burden (secondary): ", fmt_effect(site_burden), ". Burden coefficients are modeled separately."))
}
if (length(res2023)) {
    r <- res2023[[1]]$regression;summary_pick <- if (any(r$model == "V2_secondary_mean_TPM")) "V2_secondary_mean_TPM" else "V1_primary_length";v <- r[model == summary_pick & variable == "edited"];sum_lines <- c(sum_lines, "", paste0("Q8 2023 independent validation: ", fmt_effect(v)))
}
if (!file.exists(expression_2023_file)) sum_lines <- c(sum_lines, "WARNING: 2023 published expression was unavailable; expression-adjusted validation models were skipped.")
invisible(NULL)
message("Analysis complete. Outputs written under: ", normalizePath(out_dir, mustWork = FALSE))
