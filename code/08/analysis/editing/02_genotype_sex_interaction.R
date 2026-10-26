#!/usr/bin/env Rscript

# ============================================================================
# Two-factor genotype-by-sex RNA-editing interaction analysis
# Updated to follow the filtering / QC logic used in
# ANOVA_editinglevel_NB_BB_overlap_pairdiff005_full.R
#
# Main biological question:
#   Does the genotype effect on RNA-editing level depend on sex?
#
# Statistical models are fitted BEFORE biological effect-size filtering.
# Final high-confidence interaction sites are defined AFTER model fitting as:
#   1) valid in both NB and beta-binomial models;
#   2) NB interaction BH FDR < 0.05;
#   3) beta-binomial interaction BH FDR < 0.05;
#   4) at least one of the six genotype x sex mean editing levels > 0.10;
#   5) maximum absolute difference-of-differences > 0.05.
#
# Usage:
#   Rscript two_factor_editing_NB_BB_genotype_sex_interaction_UPDATED.R INPUT_DIR OUTPUT_DIR
#
# Expected files in INPUT_DIR:
#   AG.candidate.AtoG.annotations.long.tsv
#   WBG.candidate.AtoG.annotations.long.tsv
#   AS.candidate.AtoG.annotations.long.tsv
# ============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(MASS)
})

if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop(
        "VGAM is required for beta-binomial regression. Install it first, e.g. install.packages('VGAM')."
    )
}

# ============================================================================
# 0. Analysis settings
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1) {
    args[1]
} else {
    Sys.getenv("EDIT_INPUT_DIR")
}


outdir <- if (length(args) >= 2L) {
    args[[2L]]
} else {
    file.path(input_dir, "two_factor_NB_BB_genotype_sex_interaction_UPDATED")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

GENOTYPES <- c("AG", "WBG", "AS")
SEXES <- c("F", "M")

MIN_REPS_PER_CELL <- 2L

# Statistical threshold.
FDR_CUTOFF <- 0.05

# Biological effect-size thresholds, applied only AFTER interaction testing.
LEVEL_CUTOFF <- 0.10
INTERACTION_DIFF_CUTOFF <- 0.05

# Basic coverage filter.
MIN_DP <- 10L

# ============================================================================
# 1. Locate input files
# ============================================================================

find_group_file <- function(input_dir, group_name) {
    expected <- file.path(
        input_dir,
        paste0(group_name, ".candidate.AtoG.annotations.long.tsv")
    )

    if (!file.exists(expected)) {
        stop("Missing input file: ", expected)
    }

    expected
}

input_files <- setNames(
    vapply(
        GENOTYPES,
        function(g) find_group_file(input_dir, g),
        character(1)
    ),
    GENOTYPES
)

cat("Input files:\n")
print(input_files)

# ============================================================================
# 2. Read and standardize one annotations.long file
# ============================================================================

read_group <- function(path, group) {
    x <- fread(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        fill = TRUE,
        check.names = FALSE,
        na.strings = c("NA", "NaN", "")
    )

    required <- c(
        "handle",
        "chr",
        "pos",
        "refcount",
        "altcount",
        "depth",
        "lib"
    )

    missing <- setdiff(required, names(x))

    if (length(missing) > 0L) {
        stop(
            "Missing columns in ", path, ": ",
            paste(missing, collapse = ", ")
        )
    }

    if (!"gene" %in% names(x)) {
        x[, gene := NA_character_]
    }

    if (!"annotation_source" %in% names(x)) {
        x[, annotation_source := NA_character_]
    }

    x[, `:=`(
        coord_handle = trimws(as.character(handle)),
        chr = as.character(chr),
        pos = suppressWarnings(as.numeric(pos)),
        refcount = suppressWarnings(as.numeric(refcount)),
        altcount = suppressWarnings(as.numeric(altcount)),
        dp = suppressWarnings(as.numeric(depth)),
        sample_raw = trimws(as.character(lib)),
        annotation_source_raw = trimws(as.character(annotation_source)),
        genotype = group,
        gene = trimws(as.character(gene))
    )]

    # Keep supplied handle; fall back to chr_pos only when handle is missing.
    x[
        (is.na(coord_handle) | coord_handle == "") &
            !is.na(chr) & chr != "" & is.finite(pos),
        coord_handle := paste0(chr, "_", as.integer(pos))
    ]

    # AG_F1 / WBG_M3 / AS_F2 -> F1 / M3 / F2.
    x[, sample := sub(
        "^.*_([FM][0-9]+)$",
        "\\1",
        sample_raw,
        perl = TRUE
    )]

    x[!grepl("^[FM][0-9]+$", sample), sample := NA_character_]

    x[, sex := fifelse(
        !is.na(sample),
        substr(sample, 1L, 1L),
        NA_character_
    )]

    # Fallback only if annotation_source itself contains F/M.
    x[
        (is.na(sex) | sex == "") &
            toupper(annotation_source_raw) %in% SEXES,
        sex := toupper(annotation_source_raw)
    ]

    x[, rep := fifelse(
        !is.na(sample),
        sub("^[FM]", "", sample),
        NA_character_
    )]

    # Use depth as dp; recover from ref+alt only when depth is missing.
    x[
        !is.finite(dp) & is.finite(refcount) & is.finite(altcount),
        dp := refcount + altcount
    ]

    # QC only: does supplied depth equal ref+alt?
    x[, dp_count_sum := fifelse(
        is.finite(refcount) & is.finite(altcount),
        refcount + altcount,
        NA_real_
    )]

    x[, dp_mismatch := (
        is.finite(dp) &
            is.finite(dp_count_sum) &
            dp != dp_count_sum
    )]

    x[, source_file := basename(path)]

    x
}

raw_all <- rbindlist(
    lapply(
        GENOTYPES,
        function(g) read_group(input_files[[g]], g)
    ),
    use.names = TRUE,
    fill = TRUE
)

raw_all[, genotype := factor(as.character(genotype), levels = GENOTYPES)]

fwrite(
    raw_all,
    file.path(outdir, "AG_WBG_AS.long.standardized.before_filter.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

cat("\nDepth mismatch summary before filtering:\n")
print(
    raw_all[, .(
        n_rows = .N,
        n_dp_mismatch = sum(dp_mismatch, na.rm = TRUE),
        fraction_dp_mismatch = mean(dp_mismatch, na.rm = TRUE)
    ), by = genotype]
)

# ============================================================================
# 3. Basic filtering
#    IMPORTANT: altcount == 0 is retained as a valid observation.
# ============================================================================

raw <- raw_all[
    !is.na(coord_handle) & coord_handle != "" &
        !is.na(genotype) &
        !is.na(sample) &
        !is.na(rep) &
        sex %in% SEXES &
        is.finite(refcount) &
        is.finite(altcount) &
        is.finite(dp) &
        dp >= MIN_DP &
        refcount >= 0 &
        altcount >= 0 &
        altcount <= dp &
        (refcount + altcount) > 1
]

raw[, `:=`(
    refcount = as.integer(round(refcount)),
    altcount = as.integer(round(altcount)),
    dp = as.integer(round(dp))
)]

raw[, editing_level := altcount / dp]

# Also save the exact proportion modeled by beta-binomial.
raw[, bb_editing_fraction := altcount / (altcount + refcount)]

cat("\nRows after basic filtering: ", nrow(raw), "\n", sep = "")
cat("Unique sites after basic filtering: ", uniqueN(raw$coord_handle), "\n", sep = "")

cat("\nSample summary after basic filtering:\n")
print(
    raw[, .N, by = .(genotype, sex, sample, rep)][order(sex, genotype, rep)]
)

# ============================================================================
# 4. Collapse annotation duplicates to one quantitative record per site/library
# ============================================================================

site_rep <- raw[, .(
    n_annotation_rows = .N,
    n_unique_measurements = uniqueN(
        paste(refcount, altcount, dp, sep = "|")
    ),
    refcount = first(refcount),
    altcount = first(altcount),
    dp = first(dp),
    editing_level = first(editing_level),
    bb_editing_fraction = first(bb_editing_fraction)
), by = .(
    coord_handle,
    genotype,
    sex,
    sample,
    rep
)]

if (any(site_rep$n_unique_measurements > 1L)) {
    bad <- site_rep[n_unique_measurements > 1L]
    fwrite(
        bad,
        file.path(outdir, "ERROR.conflicting_site_library_measurements.tsv"),
        sep = "\t"
    )
    stop(
        "Conflicting refcount/altcount/dp records were found within at least one site/library."
    )
}

fwrite(
    site_rep,
    file.path(outdir, "site_replicate_measurement_summary.tsv"),
    sep = "\t"
)

# Optional site -> gene annotation table for downstream interpretation.
site_gene_map <- unique(
    raw_all[
        !is.na(coord_handle) & coord_handle != "" &
            !is.na(gene) & gene != "" & gene != "NA",
        .(coord_handle, gene)
    ]
)

if (nrow(site_gene_map) > 0L) {
    fwrite(
        site_gene_map,
        file.path(outdir, "site_gene_map.tsv"),
        sep = "\t"
    )
}

# ============================================================================
# 5. Eligibility: require >=2 independent replicates in EACH of the 6 cells
#       AG_F, AG_M, WBG_F, WBG_M, AS_F, AS_M
# ============================================================================

rep_counts <- site_rep[, .(
    n_reps = uniqueN(rep)
), by = .(
    coord_handle,
    genotype,
    sex
)]

expected_cells <- CJ(
    coord_handle = unique(site_rep$coord_handle),
    genotype = GENOTYPES,
    sex = SEXES,
    unique = TRUE
)

rep_counts <- merge(
    expected_cells,
    rep_counts,
    by = c("coord_handle", "genotype", "sex"),
    all.x = TRUE,
    sort = FALSE
)

rep_counts[is.na(n_reps), n_reps := 0L]

eligibility <- rep_counts[, .(
    min_reps_across_6_cells = min(n_reps),
    complete_3x2 = all(n_reps >= MIN_REPS_PER_CELL)
), by = coord_handle]

fwrite(
    rep_counts,
    file.path(outdir, "replicate_count_by_site_genotype_sex.tsv"),
    sep = "\t"
)

fwrite(
    eligibility,
    file.path(outdir, "site_eligibility_3x2.tsv"),
    sep = "\t"
)

eligible <- eligibility[complete_3x2 == TRUE, coord_handle]

dat <- site_rep[coord_handle %in% eligible]

dat[, `:=`(
    genotype = factor(as.character(genotype), levels = GENOTYPES),
    sex = factor(as.character(sex), levels = SEXES)
)]

cat(
    "\nEligible sites with >=", MIN_REPS_PER_CELL,
    " replicates in every genotype x sex cell: ",
    length(eligible), "\n",
    sep = ""
)

if (length(eligible) == 0L) {
    stop("No sites passed the complete 3x2 replicate filter.")
}

# ============================================================================
# 6. Robust model helpers
# ============================================================================

fit_nb_one <- function(formula_use, x) {
    # First try MASS default initialization.
    fit <- tryCatch(
        suppressWarnings(
            glm.nb(
                formula_use,
                data = x,
                control = glm.control(maxit = 500)
            )
        ),
        error = function(e) NULL
    )

    if (!is.null(fit)) {
        return(fit)
    }

    # Retry reasonable finite theta starts for sparse / near-Poisson sites.
    theta_candidates <- c(0.1, 1, 10, 100, 1000, 10000, 1e6)

    for (theta_init in theta_candidates) {
        fit <- tryCatch(
            suppressWarnings(
                glm.nb(
                    formula_use,
                    data = x,
                    init.theta = theta_init,
                    control = glm.control(maxit = 500)
                )
            ),
            error = function(e) NULL
        )

        if (!is.null(fit)) {
            return(fit)
        }
    }

    NULL
}

fit_bb_one <- function(formula_use, x, irho_init = NULL) {
    fam <- if (is.null(irho_init)) {
        VGAM::betabinomial(zero = 2)
    } else {
        VGAM::betabinomial(zero = 2, irho = irho_init)
    }

    tryCatch(
        suppressWarnings(
            VGAM::vglm(
                formula_use,
                family = fam,
                data = x,
                control = VGAM::vglm.control(
                    maxit = 500,
                    trace = FALSE
                )
            )
        ),
        error = function(e) NULL
    )
}

# Fit two nested beta-binomial models using the same initialization attempt.
fit_bb_pair <- function(null_formula, full_formula, x) {
    irho_candidates <- list(
        NULL,
        1e-4,
        1e-3,
        0.01,
        0.05,
        0.10,
        0.25,
        0.50
    )

    for (irho_try in irho_candidates) {
        null_fit <- fit_bb_one(null_formula, x, irho_try)
        full_fit <- fit_bb_one(full_formula, x, irho_try)

        if (!is.null(null_fit) && !is.null(full_fit)) {
            ll0 <- tryCatch(as.numeric(logLik(null_fit)), error = function(e) NA_real_)
            ll1 <- tryCatch(as.numeric(logLik(full_fit)), error = function(e) NA_real_)

            if (is.finite(ll0) && is.finite(ll1)) {
                return(list(
                    null = null_fit,
                    full = full_fit,
                    irho_init = irho_try
                ))
            }
        }
    }

    NULL
}

lrt_info <- function(null_fit, full_fit) {
    empty <- list(
        p = NA_real_,
        statistic = NA_real_,
        df = NA_real_
    )

    if (is.null(null_fit) || is.null(full_fit)) {
        return(empty)
    }

    ll0_obj <- tryCatch(logLik(null_fit), error = function(e) NULL)
    ll1_obj <- tryCatch(logLik(full_fit), error = function(e) NULL)

    if (is.null(ll0_obj) || is.null(ll1_obj)) {
        return(empty)
    }

    ll0 <- suppressWarnings(as.numeric(ll0_obj))
    ll1 <- suppressWarnings(as.numeric(ll1_obj))

    if (!is.finite(ll0) || !is.finite(ll1)) {
        return(empty)
    }

    df0 <- suppressWarnings(as.numeric(attr(ll0_obj, "df")))
    df1 <- suppressWarnings(as.numeric(attr(ll1_obj, "df")))

    if (length(df0) != 1L || !is.finite(df0)) {
        df0 <- length(stats::coef(null_fit))
    }

    if (length(df1) != 1L || !is.finite(df1)) {
        df1 <- length(stats::coef(full_fit))
    }

    df <- df1 - df0
    statistic <- 2 * (ll1 - ll0)

    if (
        length(df) != 1L ||
            length(statistic) != 1L ||
            !is.finite(df) ||
            !is.finite(statistic) ||
            df <= 0 ||
            statistic < -1e-6
    ) {
        return(empty)
    }

    statistic <- max(statistic, 0)

    list(
        p = pchisq(statistic, df = df, lower.tail = FALSE),
        statistic = statistic,
        df = df
    )
}

# ============================================================================
# 7. Fit one site
#
# Main effects are adjusted for the other factor:
#   genotype:  sex-only  vs genotype + sex
#   sex:       genotype-only vs genotype + sex
# Interaction:
#   genotype + sex vs genotype * sex
# ============================================================================

fit_one_site <- function(x, key) {
    # Safety checks inside each site.
    cell_n <- x[, .(n_reps = uniqueN(rep)), by = .(genotype, sex)]

    if (nrow(cell_n) != length(GENOTYPES) * length(SEXES) ||
        any(cell_n$n_reps < MIN_REPS_PER_CELL)) {
        return(data.table(
            coord_handle = key,
            fit_status = "failed_complete_3x2_filter"
        ))
    }

    if (any(!is.finite(x$altcount)) ||
        any(!is.finite(x$refcount)) ||
        any(!is.finite(x$dp)) ||
        any(x$altcount < 0) ||
        any(x$refcount < 0) ||
        any(x$dp <= 0) ||
        any(x$altcount > x$dp)) {
        return(data.table(
            coord_handle = key,
            fit_status = "invalid_counts"
        ))
    }

    if (sum(x$altcount, na.rm = TRUE) == 0) {
        return(data.table(
            coord_handle = key,
            fit_status = "all_altcount_zero"
        ))
    }

    # ---------------- NB ----------------
    nb_sex <- fit_nb_one(
        altcount ~ sex + offset(log(dp)),
        x
    )

    nb_geno <- fit_nb_one(
        altcount ~ genotype + offset(log(dp)),
        x
    )

    nb_add <- fit_nb_one(
        altcount ~ genotype + sex + offset(log(dp)),
        x
    )

    nb_int <- fit_nb_one(
        altcount ~ genotype * sex + offset(log(dp)),
        x
    )

    nb_genotype_lrt <- lrt_info(nb_sex, nb_add)
    nb_sex_lrt <- lrt_info(nb_geno, nb_add)
    nb_interaction_lrt <- lrt_info(nb_add, nb_int)

    # -------------- Beta-binomial --------------
    bb_genotype_pair <- fit_bb_pair(
        cbind(altcount, refcount) ~ sex,
        cbind(altcount, refcount) ~ genotype + sex,
        x
    )

    bb_sex_pair <- fit_bb_pair(
        cbind(altcount, refcount) ~ genotype,
        cbind(altcount, refcount) ~ genotype + sex,
        x
    )

    bb_interaction_pair <- fit_bb_pair(
        cbind(altcount, refcount) ~ genotype + sex,
        cbind(altcount, refcount) ~ genotype * sex,
        x
    )

    bb_genotype_lrt <- if (is.null(bb_genotype_pair)) {
        list(p = NA_real_, statistic = NA_real_, df = NA_real_)
    } else {
        lrt_info(bb_genotype_pair$null, bb_genotype_pair$full)
    }

    bb_sex_lrt <- if (is.null(bb_sex_pair)) {
        list(p = NA_real_, statistic = NA_real_, df = NA_real_)
    } else {
        lrt_info(bb_sex_pair$null, bb_sex_pair$full)
    }

    bb_interaction_lrt <- if (is.null(bb_interaction_pair)) {
        list(p = NA_real_, statistic = NA_real_, df = NA_real_)
    } else {
        lrt_info(bb_interaction_pair$null, bb_interaction_pair$full)
    }

    data.table(
        coord_handle = key,
        fit_status = "ok",
        NB_p_genotype = nb_genotype_lrt$p,
        NB_LRT_genotype = nb_genotype_lrt$statistic,
        NB_df_genotype = nb_genotype_lrt$df,
        NB_p_sex = nb_sex_lrt$p,
        NB_LRT_sex = nb_sex_lrt$statistic,
        NB_df_sex = nb_sex_lrt$df,
        NB_p_interaction = nb_interaction_lrt$p,
        NB_LRT_interaction = nb_interaction_lrt$statistic,
        NB_df_interaction = nb_interaction_lrt$df,
        BB_p_genotype = bb_genotype_lrt$p,
        BB_LRT_genotype = bb_genotype_lrt$statistic,
        BB_df_genotype = bb_genotype_lrt$df,
        BB_p_sex = bb_sex_lrt$p,
        BB_LRT_sex = bb_sex_lrt$statistic,
        BB_df_sex = bb_sex_lrt$df,
        BB_p_interaction = bb_interaction_lrt$p,
        BB_LRT_interaction = bb_interaction_lrt$statistic,
        BB_df_interaction = bb_interaction_lrt$df,
        NB_fit_additive = !is.null(nb_add),
        NB_fit_interaction = !is.null(nb_int),
        BB_fit_additive = !is.null(bb_interaction_pair),
        BB_fit_interaction = !is.null(bb_interaction_pair)
    )
}

# ============================================================================
# 8. Fit all eligible sites with checkpointing
# ============================================================================

checkpoint_file <- file.path(
    outdir,
    "two_factor_interaction_fit_checkpoint_UPDATED.tsv"
)

# Reuse model fits only when inputs, script, parameters, and package versions match.
checkpoint_signature_file <- paste0(checkpoint_file, ".signature.rds")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
checkpoint_signature <- list(
    inputs = tools::md5sum(unname(input_files)), script = tools::md5sum(script_path),
    min_dp = MIN_DP, min_reps = MIN_REPS_PER_CELL,
    MASS = as.character(packageVersion("MASS")), VGAM = as.character(packageVersion("VGAM"))
)
if (file.exists(checkpoint_file) &&
    (!file.exists(checkpoint_signature_file) ||
        !identical(readRDS(checkpoint_signature_file), checkpoint_signature)))
    stop("Checkpoint provenance differs. Use a new output directory for this analysis.")
saveRDS(checkpoint_signature, checkpoint_signature_file)

previous <- if (file.exists(checkpoint_file)) {
    fread(checkpoint_file)
} else {
    data.table()
}

if (nrow(previous) > 0L) {
    previous <- previous[coord_handle %in% eligible]
    previous <- unique(previous, by = "coord_handle")
    setorder(previous, coord_handle)
}

remaining <- setdiff(eligible, previous$coord_handle)

cat(
    "\nResuming with ", nrow(previous),
    " completed and ", length(remaining),
    " remaining sites.\n",
    sep = ""
)

new_results <- vector("list", length(remaining))

for (i in seq_along(remaining)) {
    key <- remaining[[i]]

    new_results[[i]] <- fit_one_site(
        dat[coord_handle == key],
        key
    )

    if (i %% 25L == 0L || i == length(remaining)) {
        completed <- rbindlist(
            c(
                list(previous),
                new_results[seq_len(i)]
            ),
            use.names = TRUE,
            fill = TRUE
        )

        fwrite(
            completed,
            checkpoint_file,
            sep = "\t"
        )

        cat(
            "Fitted ", nrow(previous) + i,
            "/", length(eligible),
            " sites; checkpoint saved.\n",
            sep = ""
        )
    }
}

results <- rbindlist(
    c(list(previous), new_results),
    use.names = TRUE,
    fill = TRUE
)

results <- unique(results, by = "coord_handle")

cat("\nFit-status summary:\n")
print(table(results$fit_status, useNA = "ifany"))

# ============================================================================
# 9. Common-valid testing universe + BH correction
#
# For each effect, NB and beta-binomial P values are BH-corrected over the same
# pipeline and makes the two model-specific FDR families directly comparable.
# ============================================================================

effects <- c("genotype", "sex", "interaction")

for (effect in effects) {
    nb_p <- paste0("NB_p_", effect)
    bb_p <- paste0("BB_p_", effect)

    valid_col <- paste0("common_valid_", effect)
    nb_fdr <- paste0("NB_FDR_", effect)
    bb_fdr <- paste0("BB_FDR_", effect)

    nb_sig <- paste0("NB_significant_", effect)
    bb_sig <- paste0("BB_significant_", effect)
    both_sig <- paste0("both_models_FDR005_", effect)
    nb_only <- paste0("NB_only_FDR005_", effect)
    bb_only <- paste0("BB_only_FDR005_", effect)
    either_sig <- paste0("either_model_FDR005_", effect)

    results[, (valid_col) := (
        is.finite(get(nb_p)) &
            is.finite(get(bb_p))
    )]

    results[, (nb_fdr) := NA_real_]
    results[, (bb_fdr) := NA_real_]

    common_idx <- which(results[[valid_col]] == TRUE)

    if (length(common_idx) > 0L) {
        results[common_idx, (nb_fdr) := p.adjust(
            get(nb_p),
            method = "BH"
        )]

        results[common_idx, (bb_fdr) := p.adjust(
            get(bb_p),
            method = "BH"
        )]
    }

    results[, (nb_sig) := (
        get(valid_col) &
            is.finite(get(nb_fdr)) &
            get(nb_fdr) < FDR_CUTOFF
    )]

    results[, (bb_sig) := (
        get(valid_col) &
            is.finite(get(bb_fdr)) &
            get(bb_fdr) < FDR_CUTOFF
    )]

    results[, (both_sig) := get(nb_sig) & get(bb_sig)]
    results[, (nb_only) := get(nb_sig) & !get(bb_sig)]
    results[, (bb_only) := !get(nb_sig) & get(bb_sig)]
    results[, (either_sig) := get(nb_sig) | get(bb_sig)]
}

# ============================================================================
# 10. Mean editing levels and interaction effect sizes
#
# Effect-size filtering is applied AFTER statistical testing.
# ============================================================================

means <- dat[, .(
    mean_editing_level = mean(editing_level, na.rm = TRUE),
    median_editing_level = median(editing_level, na.rm = TRUE),
    mean_bb_fraction = mean(bb_editing_fraction, na.rm = TRUE),
    mean_altcount = mean(altcount, na.rm = TRUE),
    mean_refcount = mean(refcount, na.rm = TRUE),
    mean_dp = mean(dp, na.rm = TRUE),
    n_reps = uniqueN(rep)
), by = .(
    coord_handle,
    genotype,
    sex
)]

# Mean Lm matrix: AG_F, AG_M, WBG_F, WBG_M, AS_F, AS_M.
means_wide <- dcast(
    means,
    coord_handle ~ genotype + sex,
    value.var = "mean_editing_level"
)

required_mean_cols <- as.vector(
    outer(
        GENOTYPES,
        SEXES,
        paste,
        sep = "_"
    )
)

for (cc in required_mean_cols) {
    if (!cc %in% names(means_wide)) {
        means_wide[, (cc) := NA_real_]
    }
}

# Replicate counts by cell.
reps_wide <- dcast(
    means,
    coord_handle ~ genotype + sex,
    value.var = "n_reps"
)

rep_cols <- setdiff(names(reps_wide), "coord_handle")
setnames(
    reps_wide,
    rep_cols,
    paste0("nrep_", rep_cols)
)

results <- merge(
    results,
    means_wide,
    by = "coord_handle",
    all.x = TRUE,
    sort = FALSE
)

results <- merge(
    results,
    reps_wide,
    by = "coord_handle",
    all.x = TRUE,
    sort = FALSE
)

results <- merge(
    results,
    eligibility,
    by = "coord_handle",
    all.x = TRUE,
    sort = FALSE
)

# At least one of the six cells has mean editing level > 0.10.
results[, pass_level_cutoff := rowSums(
    .SD > LEVEL_CUTOFF,
    na.rm = TRUE
) > 0,
.SDcols = required_mean_cols
]

# Within-sex genotype contrasts, useful for interpreting an interaction.
results[, `:=`(
    delta_AG_WBG_F = AG_F - WBG_F,
    delta_AG_WBG_M = AG_M - WBG_M,
    delta_AS_WBG_F = AS_F - WBG_F,
    delta_AS_WBG_M = AS_M - WBG_M,
    delta_AG_AS_F = AG_F - AS_F,
    delta_AG_AS_M = AG_M - AS_M,
    sex_delta_AG = AG_F - AG_M,
    sex_delta_WBG = WBG_F - WBG_M,
    sex_delta_AS = AS_F - AS_M
)]

# Difference-of-differences (DoD): the genotype contrast in Female minus the
# same genotype contrast in Male. This is the biologically relevant magnitude
# of genotype x sex interaction.
results[, `:=`(
    interaction_DID_AG_WBG = delta_AG_WBG_F - delta_AG_WBG_M,
    interaction_DID_AS_WBG = delta_AS_WBG_F - delta_AS_WBG_M,
    interaction_DID_AG_AS = delta_AG_AS_F - delta_AG_AS_M
)]

results[, max_abs_interaction_DID := pmax(
    abs(interaction_DID_AG_WBG),
    abs(interaction_DID_AS_WBG),
    abs(interaction_DID_AG_AS),
    na.rm = FALSE
)]

results[, pass_interaction_diff_cutoff := (
    is.finite(max_abs_interaction_DID) &
        max_abs_interaction_DID > INTERACTION_DIFF_CUTOFF
)]

# Label the contrast contributing the largest observed interaction magnitude.
results[, dominant_interaction_contrast := fifelse(
    !is.finite(max_abs_interaction_DID),
    NA_character_,
    fifelse(
        abs(interaction_DID_AG_WBG) == max_abs_interaction_DID,
        "AG_vs_WBG",
        fifelse(
            abs(interaction_DID_AS_WBG) == max_abs_interaction_DID,
            "AS_vs_WBG",
            "AG_vs_AS"
        )
    )
)]

# Final high-confidence interaction definition.
results[, high_confidence_interaction := (
    common_valid_interaction == TRUE &
        NB_significant_interaction == TRUE &
        BB_significant_interaction == TRUE &
        pass_level_cutoff == TRUE &
        pass_interaction_diff_cutoff == TRUE
)]

# ============================================================================
# 11. Outputs
# ============================================================================

all_results_file <- file.path(
    outdir,
    "two_factor_NB_BetaBinomial_all_eligible_sites_UPDATED.tsv"
)

fwrite(results, all_results_file, sep = "\t")

# Common-valid sites for interaction.
fwrite(
    results[common_valid_interaction == TRUE],
    file.path(
        outdir,
        "two_factor_interaction_NB_BetaBinomial_common_valid.tsv"
    ),
    sep = "\t"
)

# Statistical interaction significance before biological effect-size filters.
fwrite(
    results[both_models_FDR005_interaction == TRUE],
    file.path(
        outdir,
        "two_factor_interaction_both_models_FDR005.before_effect_filter.tsv"
    ),
    sep = "\t"
)

# Final high-confidence interaction sites.
final_interaction_file <- file.path(
    outdir,
    paste0(
        "two_factor_interaction_both_models_FDR005.level_gt",
        LEVEL_CUTOFF,
        ".max_abs_DID_gt",
        INTERACTION_DIFF_CUTOFF,
        ".tsv"
    )
)

fwrite(
    results[high_confidence_interaction == TRUE],
    final_interaction_file,
    sep = "\t"
)

# Also retain main-effect significant sites without applying interaction-specific
# biological effect filters.
for (effect in c("genotype", "sex")) {
    keep <- paste0("both_models_FDR005_", effect)
    fwrite(
        results[get(keep) == TRUE],
        file.path(
            outdir,
            paste0("two_factor_", effect, "_both_models_FDR005.tsv")
        ),
        sep = "\t"
    )
}

# ============================================================================
# 12. Summary tables
# ============================================================================

summary_effects <- rbindlist(
    lapply(effects, function(effect) {
        valid_col <- paste0("common_valid_", effect)
        nb_sig <- paste0("NB_significant_", effect)
        bb_sig <- paste0("BB_significant_", effect)
        both_sig <- paste0("both_models_FDR005_", effect)
        nb_only <- paste0("NB_only_FDR005_", effect)
        bb_only <- paste0("BB_only_FDR005_", effect)
        either_sig <- paste0("either_model_FDR005_", effect)

        data.table(
            effect = effect,
            n_eligible_sites = nrow(results),
            n_common_valid = sum(results[[valid_col]], na.rm = TRUE),
            n_NB_FDR005 = sum(results[[nb_sig]], na.rm = TRUE),
            n_BetaBinomial_FDR005 = sum(results[[bb_sig]], na.rm = TRUE),
            n_both_models_FDR005 = sum(results[[both_sig]], na.rm = TRUE),
            n_NB_only = sum(results[[nb_only]], na.rm = TRUE),
            n_BB_only = sum(results[[bb_only]], na.rm = TRUE),
            n_either_model = sum(results[[either_sig]], na.rm = TRUE)
        )
    })
)

fwrite(
    summary_effects,
    file.path(outdir, "two_factor_effect_significance_summary_UPDATED.tsv"),
    sep = "\t"
)

interaction_filter_summary <- data.table(
    n_eligible_sites = nrow(results),
    n_common_valid_interaction = sum(
        results$common_valid_interaction,
        na.rm = TRUE
    ),
    n_both_models_interaction_FDR005 = sum(
        results$both_models_FDR005_interaction,
        na.rm = TRUE
    ),
    n_both_sig_pass_level = sum(
        results$both_models_FDR005_interaction &
            results$pass_level_cutoff,
        na.rm = TRUE
    ),
    n_both_sig_pass_DID = sum(
        results$both_models_FDR005_interaction &
            results$pass_interaction_diff_cutoff,
        na.rm = TRUE
    ),
    n_final_high_confidence_interaction = sum(
        results$high_confidence_interaction,
        na.rm = TRUE
    )
)

fwrite(
    interaction_filter_summary,
    file.path(outdir, "interaction_effect_filter_summary_UPDATED.tsv"),
    sep = "\t"
)

cat("\nEffect significance summary:\n")
print(summary_effects)

cat("\nInteraction effect-filter summary:\n")
print(interaction_filter_summary)

# ============================================================================
# 13. README + session info
# ============================================================================

writeLines(
    c(
        "Two-factor genotype x sex RNA-editing interaction analysis.",
        "",
        "BASIC FILTERING:",
        paste0("- depth >= ", MIN_DP),
        "- refcount >= 0 and altcount >= 0",
        "- altcount == 0 is retained; zero editing is a valid observation",
        "- altcount <= depth",
        "- refcount + altcount > 1",
        "- annotation duplicate rows are collapsed only when counts agree",
        paste0("- each of the 6 genotype x sex cells must have >= ", MIN_REPS_PER_CELL, " independent replicates"),
        "",
        "MODELS:",
        "- NB uses altcount with offset(log(depth))",
        "- beta-binomial uses cbind(altcount, refcount)",
        "- genotype main effect: sex-only vs genotype + sex",
        "- sex main effect: genotype-only vs genotype + sex",
        "- interaction: genotype + sex vs genotype * sex",
        "",
        "MULTIPLE TESTING:",
        "- for each effect, NB and beta-binomial are BH-corrected over the same common-valid site universe",
        paste0("- statistical significance cutoff: FDR < ", FDR_CUTOFF),
        "",
        "FINAL HIGH-CONFIDENCE INTERACTION FILTER (applied AFTER model fitting):",
        "- interaction significant in BOTH NB and beta-binomial models",
        paste0("- at least one of AG_F/AG_M/WBG_F/WBG_M/AS_F/AS_M has mean Lm > ", LEVEL_CUTOFF),
        paste0("- max absolute genotype-contrast difference-of-differences > ", INTERACTION_DIFF_CUTOFF),
        "",
        "INTERACTION EFFECT SIZE:",
        "- AG_vs_WBG DoD = (AG_F - WBG_F) - (AG_M - WBG_M)",
        "- AS_vs_WBG DoD = (AS_F - WBG_F) - (AS_M - WBG_M)",
        "- AG_vs_AS  DoD = (AG_F - AS_F)  - (AG_M - AS_M)",
        "- max_abs_interaction_DID is the maximum absolute value among these three",
        "",
        "IMPORTANT:",
        "- No within-sex genotype-significance filter is applied before interaction testing.",
        "- The 0.10 level and 0.05 interaction-effect cutoffs are downstream biological-effect filters, not model eligibility criteria.",
        "- editing_level/Lm summaries use altcount/depth to remain consistent with the previous pipeline; beta-binomial itself models altcount/(altcount+refcount)."
    ),
    file.path(outdir, "README_UPDATED.txt")
)

invisible(NULL)

cat("\nDone.\n")
cat("All-site result: ", all_results_file, "\n", sep = "")
cat("Final interaction result: ", final_interaction_file, "\n", sep = "")
