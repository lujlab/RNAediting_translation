suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

dir_create <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)
make_output_dirs <- function(out_dir) {
    dirs <- file.path(out_dir, c("00_input_QC", "01_gene_mapping", "02_gene_features", "03_2010_stable_decay",
        "04_2010_half_life", "05_2010_regression", "06_matched_analysis", "07_2023_validation",
        "08_cross_dataset", "09_figures", "10_tables"))
    vapply(dirs, dir_create, logical(1))
    dirs
}
write_tsv <- function(x, file) {
    dir_create(dirname(file))
    fwrite(as.data.table(x), file, sep = "\t", na = "")
}
clean_id <- function(x) sub("\\.[0-9]+$", "", trimws(as.character(x)))
or_default <- function(x, default) if (is.null(x) || !length(x) || is.na(x)) default else x
find_col <- function(x, choices, required = TRUE) {
    n <- names(x)
    hit <- which(tolower(n) %in% tolower(choices))
    if (!length(hit))
        hit <- which(vapply(n, function(z) any(vapply(choices, function(q) grepl(q, z, ignore.case = TRUE,
            fixed = TRUE), logical(1))), logical(1)))
    if (!length(hit) && required)
        stop("Could not find a column matching: ", paste(choices, collapse = ", "), ". Available: ",
            paste(n, collapse = ", "))
    if (length(hit))
        n[hit[1]] else NA_character_
}
standardize_gene_column <- function(dt, dataset) {
    col <- find_col(dt, c("FBgn", "gene_id", "gene", "gene_symbol", "symbol"))
    dt[, FBgn := clean_id(get(col))]
    if (all(!grepl("^FBgn", dt$FBgn)))
        message(dataset, ": non-FBgn IDs will be mapped via annotation gene_symbol when possible.")
    dt
}
safe_log2 <- function(x, offset = 0) ifelse(is.finite(x) & x + offset > 0, log2(x + offset),
    NA_real_)
rank_biserial <- function(x, group) {
    ok <- is.finite(x) & !is.na(group)
    x <- x[ok]
    group <- as.character(group[ok])
    lev <- unique(group)
    if (length(lev) != 2 || !all(table(group) > 0))
        return(NA_real_)
    a <- x[group == lev[1]]
    b <- x[group == lev[2]]
    (sum(outer(a, b, ">")) + 0.5 * sum(outer(a, b, "==")))/(length(a) * length(b)) * 2 -
        1
}
hc3_vcov <- function(model) {
    X <- model.matrix(model)
    e <- residuals(model)
    h <- hatvalues(model)
    bread <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(bread))
        return(matrix(NA_real_, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X))))
    bread %*% crossprod(X, X * (e/pmax(1 - h, .Machine$double.eps))^2) %*% bread
}
tidy_lm <- function(model, model_name, dataset) {
    sm <- summary(model)$coefficients
    vc <- hc3_vcov(model)
    robust_se <- sqrt(diag(vc))
    robust_z <- sm[, "Estimate"]/robust_se
    data.table(dataset = dataset, model = model_name, variable = rownames(sm), estimate = sm[,
        "Estimate"], SE = sm[, "Std. Error"], ordinary_P = sm[, ncol(sm)], robust_SE_HC3 = robust_se,
        robust_P_HC3 = 2 * pnorm(abs(robust_z), lower.tail = FALSE), CI_low95 = sm[, "Estimate"] -
            qnorm(0.975) * sm[, "Std. Error"], CI_high95 = sm[, "Estimate"] + qnorm(0.975) *
            sm[, "Std. Error"], N_total = nobs(model), adjusted_R2 = summary(model)$adj.r.squared,
        AIC = AIC(model))
}
tidy_glm <- function(model, model_name, dataset) {
    sm <- summary(model)$coefficients
    data.table(dataset = dataset, model = model_name, variable = rownames(sm), estimate = sm[,
        "Estimate"], SE = sm[, "Std. Error"], OR = exp(sm[, "Estimate"]), OR_low95 = exp(sm[,
        "Estimate"] - qnorm(0.975) * sm[, "Std. Error"]), OR_high95 = exp(sm[, "Estimate"] +
        qnorm(0.975) * sm[, "Std. Error"]), P = sm[, ncol(sm)], N_total = nobs(model), AIC = AIC(model))
}
vif_manual <- function(model, model_name, dataset) {
    X <- model.matrix(model)[, -1, drop = FALSE]
    if (!ncol(X))
        return(data.table(dataset = dataset, model = model_name, variable = character(),
            VIF = numeric(), warning = character()))
    rbindlist(lapply(seq_len(ncol(X)), function(i) {
        z <- X[, -i, drop = FALSE]
        r2 <- if (!ncol(z))
            0 else summary(lm(X[, i] ~ z))$r.squared
        v <- 1/max(1 - r2, .Machine$double.eps)
        data.table(dataset = dataset, model = model_name, variable = colnames(X)[i], VIF = v,
            warning = ifelse(v >= 5, "VIF >= 5; inspect collinearity", ""))
    }))
}
diagnose_lm <- function(model, model_name, dataset, path) {
    d <- data.table(dataset = dataset, model = model_name, N = nobs(model), residual_shapiro_P = if (nobs(model) >=
        3 && nobs(model) <= 5000)
        shapiro.test(residuals(model))$p.value else NA_real_, max_cooks_distance = max(cooks.distance(model), na.rm = TRUE), n_cooks_gt_4_over_n = sum(cooks.distance(model) >
        4/nobs(model), na.rm = TRUE))
    pdf(path, width = 8, height = 8)
    par(mfrow = c(2, 2))
    plot(model)
    dev.off()
    d
}
save_plot <- function(p, filename, width = 7, height = 5) {
    ggsave(filename, p, width = width, height = height, device = "pdf")
    if (isTRUE(write_png))
        ggsave(sub("\\.pdf$", ".png", filename), p, width = width, height = height, dpi = 300)
}
complete_case_report <- function(dt, vars, dataset, model) {
    ok <- complete.cases(dt[, ..vars])
    data.table(dataset = dataset, model = model, FBgn = dt$FBgn, included = ok, missing_variables = vapply(seq_len(nrow(dt)),
        function(i) paste(vars[is.na(unlist(dt[i, ..vars]))], collapse = ";"), character(1)))
}
