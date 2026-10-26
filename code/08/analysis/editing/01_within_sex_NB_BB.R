suppressPackageStartupMessages({
    library(data.table)
    library(MASS)
    library(ggplot2)
    library(gridExtra)
    library(UpSetR)
    library(ggprism)
})
if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop(
        "The beta-binomial analysis requires the R package 'VGAM'. ",
        "Install it before running this script, for example with: ",
        "install.packages('VGAM')"
    )
}

# ============================================================
# Input and output paths
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1) {
    args[1]
} else {
    Sys.getenv("EDIT_INPUT_DIR")
}

outdir <- if (length(args) >= 2) {
    args[2]
} else {
    Sys.getenv("EDIT_INPUT_DIR")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

MIN_REP_PER_GENOTYPE <- 2
GENOTYPE_LEVELS <- c("AG", "WBG", "AS")
TOP_N_CLUSTERS <- 6L

# Statistical and biological-effect thresholds.
# F1/F2/F3 and M1/M2/M3 are independent biological replicates, so replicate
# number is not included as a cross-genotype block term in either model.
FDR_CUTOFF <- 0.05
LEVEL_CUTOFF <- 0.10
PAIRWISE_DIFF_CUTOFF <- 0.05

# ============================================================
# 1. Locate AG, AS and WBG files
# ============================================================

find_group_file <- function(input_dir, group_name) {
    matched_files <- list.files(
        input_dir,
        pattern = paste0(
            "^",
            group_name,
            "\\.candidate\\.AtoG\\.annotations\\.long\\.tsv$"
        ),
        full.names = TRUE,
        ignore.case = TRUE
    )

    if (length(matched_files) == 0) {
        stop(
            "Cannot find input file for ",
            group_name,
            " in:\n",
            input_dir,
            "\nExpected file: ",
            group_name,
            ".candidate.AtoG.annotations.long.tsv"
        )
    }

    if (length(matched_files) > 1) {
        stop(
            "More than one file was found for ",
            group_name,
            ":\n",
            paste(matched_files, collapse = "\n")
        )
    }

    matched_files[1]
}

input_files <- c(
    AG = find_group_file(input_dir, "AG"),
    WBG = find_group_file(input_dir, "WBG"),
    AS = find_group_file(input_dir, "AS")
)

cat("Input files:\n")
print(input_files)

# ============================================================
# 2. Read and standardize one annotations.long file
# ============================================================

read_new_candidate_file <- function(file, genotype_name) {
    cat(
        "Reading ",
        genotype_name,
        ": ",
        file,
        "\n",
        sep = ""
    )

    x <- fread(
        file,
        sep = "\t",
        header = TRUE,
        quote = "",
        fill = TRUE,
        data.table = TRUE,
        check.names = FALSE,
        na.strings = c("NA", "NaN", "")
    )

    required_new_cols <- c(
        "handle",
        "chr",
        "pos",
        "Ref",
        "Alt",
        "refcount",
        "altcount",
        "gene",
        "lib",
        "annotation_source",
        "depth"
    )

    missing_new_cols <- setdiff(
        required_new_cols,
        colnames(x)
    )

    if (length(missing_new_cols) > 0) {
        stop(
            "Missing required columns in:\n",
            file,
            "\nMissing columns: ",
            paste(missing_new_cols, collapse = ", "),
            "\nAvailable columns: ",
            paste(colnames(x), collapse = ", ")
        )
    }

    # e.g. WBG_M3, and `annotation_source` contains the library/genotype label.
    # Sex and replicate are therefore recovered from the lib suffix below.
    x[, sample_raw := trimws(as.character(lib))]
    x[, annotation_source_raw := trimws(as.character(annotation_source))]

    x[, `:=`(
        genotype = genotype_name,
        chr = as.character(chr),
        pos = suppressWarnings(as.numeric(pos)),
        refcount = suppressWarnings(as.numeric(refcount)),
        altcount = suppressWarnings(as.numeric(altcount)),
        depth = suppressWarnings(as.numeric(depth)),
        sex = NA_character_
    )]

    # Use the supplied site handle.  Fall back to chr + pos only if handle is
    x[, coord_handle := trimws(as.character(handle))]
    x[
        (is.na(coord_handle) | coord_handle == "") &
            !is.na(chr) &
            chr != "" &
            is.finite(pos),
        coord_handle := paste0(
            chr,
            "_",
            as.integer(pos)
        )
    ]

    # WBG_M3 / AG_F1 / AS_M2 -> M3 / F1 / M2
    x[
        ,
        sample := sub(
            "^.*_([FM][0-9]+)$",
            "\\1",
            sample_raw,
            perl = TRUE
        )
    ]

    x[
        !grepl(
            "^[FM][0-9]+$",
            sample
        ),
        sample := NA_character_
    ]

    # Recover sex from lib.  annotation_source is retained as provenance and
    # is only used as a fallback if it itself contains F/M.
    x[
        !is.na(sample),
        sex := substr(sample, 1, 1)
    ]
    x[
        is.na(sex) &
            toupper(annotation_source_raw) %in% c("F", "M"),
        sex := toupper(annotation_source_raw)
    ]

    x[
        is.na(sex) |
            sex == "" |
            !sex %in% c("F", "M"),
        sex := fifelse(
            !is.na(sample),
            substr(sample, 1, 1),
            NA_character_
        )
    ]

    # Use depth as dp.
    x[, dp := depth]

    # Only recover dp from counts when depth is missing.
    x[
        !is.finite(dp) &
            is.finite(refcount) &
            is.finite(altcount),
        dp := refcount + altcount
    ]

    # Check whether depth agrees with refcount + altcount.
    x[
        ,
        dp_count_sum := fifelse(
            is.finite(refcount) &
                is.finite(altcount),
            refcount + altcount,
            NA_real_
        )
    ]

    x[
        ,
        dp_mismatch := (
            is.finite(dp) &
                is.finite(dp_count_sum) &
                dp != dp_count_sum
        )
    ]

    x[, source_file := basename(file)]

    x
}

# ============================================================
# 3. Combine AG, WBG and AS
# ============================================================

dt <- rbindlist(
    lapply(
        names(input_files),
        function(g) {
            read_new_candidate_file(
                file = input_files[[g]],
                genotype_name = g
            )
        }
    ),
    use.names = TRUE,
    fill = TRUE
)

cat("Combined input rows:\n")
print(nrow(dt))

cat("Combined input columns:\n")
print(colnames(dt))

cat("Depth mismatch summary:\n")
print(
    dt[
        ,
        .(
            n_rows = .N,
            n_dp_mismatch = sum(
                dp_mismatch,
                na.rm = TRUE
            )
        ),
        by = genotype
    ]
)

# Save the standardized input for checking.
standardized_input_file <- file.path(
    outdir,
    "AG_AS_WBG_sample_sites.long.standardized.tsv"
)

fwrite(
    dt,
    standardized_input_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

cat(
    "Written standardized input table:\n",
    standardized_input_file,
    "\n"
)

# ============================================================
# 3.1 Arrange standardized columns
# ============================================================

required_cols <- c(
    "genotype",
    "coord_handle",
    "refcount",
    "altcount",
    "dp",
    "sample",
    "sex"
)

missing_cols <- setdiff(
    required_cols,
    colnames(dt)
)

if (length(missing_cols) > 0) {
    stop(
        "Missing required standardized columns: ",
        paste(missing_cols, collapse = ", ")
    )
}

dt[, genotype := factor(
    as.character(genotype),
    levels = GENOTYPE_LEVELS
)]

dt[, sex := as.character(sex)]
dt[, sample := as.character(sample)]
dt[, coord_handle := as.character(coord_handle)]

dt[, refcount := as.numeric(refcount)]
dt[, altcount := as.numeric(altcount)]
dt[, dp := as.numeric(dp)]

dt[
    is.na(dp) &
        !is.na(refcount) &
        !is.na(altcount),
    dp := refcount + altcount
]

# F1 -> 1; M3 -> 3
dt[
    ,
    rep := sub(
        "^[FM]",
        "",
        sample
    )
]

dt[
    is.na(sex) |
        sex == "",
    sex := substr(
        sample,
        1,
        1
    )
]

cat("Genotype table:\n")
print(table(dt$genotype, useNA = "ifany"))

cat("Sex table:\n")
print(table(dt$sex, useNA = "ifany"))

cat("Sample table:\n")
print(table(dt$sample, useNA = "ifany"))

cat("Rep table:\n")
print(table(dt$rep, useNA = "ifany"))

# ============================================================
# 4. Filter
# ============================================================

dt <- dt[!is.na(coord_handle)]
dt <- dt[!is.na(genotype)]
dt <- dt[!is.na(sex)]
dt <- dt[!is.na(sample)]
dt <- dt[!is.na(rep)]

dt <- dt[!is.na(refcount) & !is.na(altcount) & !is.na(dp)]
dt <- dt[dp >= 10]
dt <- dt[refcount >= 0]
dt <- dt[altcount >= 0]
dt <- dt[altcount <= dp]
dt <- dt[(refcount + altcount) > 1]

dt[, refcount := as.integer(round(refcount))]
dt[, altcount := as.integer(round(altcount))]
dt[, dp := as.integer(round(dp))]

# Keep the original definition of Lm/editing level for all summaries and plots.
dt[, editing_level := altcount / dp]

cat("Rows after basic filtering:\n")
print(nrow(dt))

cat("Sample summary after filtering:\n")
print(dt[, .N, by = .(genotype, sex, sample, rep)][order(sex, genotype, rep)])

# The new annotations.long files contain one row per annotation, so a single
# site/sample can occur multiple times.  Counts are repeated across those
# annotation rows; collapse to one quantitative record per site and library
# before replicate counting, NB/beta-binomial modeling, level summaries, and plots.
site_rep_measurements <- dt[
    ,
    .(
        n_annotation_rows = .N,
        n_unique_measurements = uniqueN(
            paste(refcount, altcount, dp, sep = "|")
        ),
        refcount = first(refcount),
        altcount = first(altcount),
        dp = first(dp),
        editing_level = first(editing_level)
    ),
    by = .(coord_handle, sex, genotype, sample, rep)
]

if (any(site_rep_measurements$n_unique_measurements > 1L)) {
    stop(
        "Conflicting refcount/altcount/depth values occur within at least one ",
        "coord_handle + sex + genotype + sample record."
    )
}

fwrite(
    site_rep_measurements,
    file.path(outdir, "site_replicate_measurement_summary.tsv"),
    sep = "\t"
)

dt_site_rep <- site_rep_measurements

cat("Quantitative site-replicate rows after annotation collapse:\n")
print(nrow(dt_site_rep))
cat("Maximum annotation rows per site-replicate:\n")
print(max(dt_site_rep$n_annotation_rows))

group_rep_n <- dt_site_rep[, .(
    n_records = .N,
    n_reps = uniqueN(rep),
    mean_level = mean(editing_level, na.rm = TRUE)
), by = .(coord_handle, sex, genotype)]

rep_wide <- dcast(
    group_rep_n,
    coord_handle + sex ~ genotype,
    value.var = "n_reps",
    fill = 0
)

for (g in GENOTYPE_LEVELS) {
    if (!g %in% colnames(rep_wide)) {
        rep_wide[, (g) := 0]
    }
}

rep_wide[, complete_3genotypes := AG >= MIN_REP_PER_GENOTYPE &
    WBG >= MIN_REP_PER_GENOTYPE &
    AS >= MIN_REP_PER_GENOTYPE]

eligible_keys <- rep_wide[complete_3genotypes == TRUE, .(coord_handle, sex)]
ineligible_keys <- rep_wide[complete_3genotypes == FALSE]

cat("Eligibility summary:\n")
cat("Total coord_handle + sex combinations: ", nrow(rep_wide), "\n")
cat("Eligible complete combinations: ", nrow(eligible_keys), "\n")
cat("Ineligible combinations: ", nrow(ineligible_keys), "\n")

cat("Ineligible reason summary by available reps:\n")
print(rep_wide[, .N, by = .(AG, WBG, AS, complete_3genotypes)][order(complete_3genotypes, AG, WBG, AS)])

fwrite(
    rep_wide,
    file.path(outdir, "eligibility_replicate_count_by_site_sex.tsv"),
    sep = "\t"
)

dt_eligible <- merge(
    dt_site_rep,
    eligible_keys,
    by = c("coord_handle", "sex"),
    all = FALSE
)

cat("Rows used in NB and beta-binomial models after complete-3-genotype filtering:\n")
print(nrow(dt_eligible))


# -------------------------------
# -------------------------------
run_nb_genotype_test_complete <- function(dt_input) {
    required_cols <- c(
        "coord_handle",
        "sex",
        "genotype",
        "rep",
        "altcount",
        "dp"
    )

    missing_cols <- setdiff(
        required_cols,
        colnames(dt_input)
    )

    if (length(missing_cols) > 0) {
        stop(
            "Missing required columns: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    if (nrow(dt_input) == 0) {
        return(
            data.table(
                coord_handle = character(),
                sex = character(),
                pvalue = numeric(),
                status = character(),
                warning_message = character()
            )
        )
    }

    keys <- unique(
        dt_input[
            ,
            .(
                coord_handle,
                sex
            )
        ]
    )

    res_list <- vector(
        mode = "list",
        length = nrow(keys)
    )

    for (i in seq_len(nrow(keys))) {
        h <- keys$coord_handle[i]
        sx <- keys$sex[i]

        x <- copy(
            dt_input[
                coord_handle == h &
                    sex == sx
            ]
        )

        x[
            ,
            genotype := factor(
                as.character(genotype),
                levels = GENOTYPE_LEVELS
            )
        ]

        # ========================================================
        # Check replicate number in each genotype
        # ========================================================

        check_n <- x[
            ,
            .(
                n_reps = uniqueN(
                    rep[
                        !is.na(rep)
                    ]
                )
            ),
            by = .(
                genotype = as.character(genotype)
            )
        ]

        check_n <- merge(
            data.table(
                genotype = GENOTYPE_LEVELS
            ),
            check_n,
            by = "genotype",
            all.x = TRUE,
            sort = FALSE
        )

        check_n[
            is.na(n_reps),
            n_reps := 0L
        ]

        if (
            any(
                check_n$n_reps <
                    MIN_REP_PER_GENOTYPE
            )
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                status = "failed_complete_filter_inside_model",
                warning_message = NA_character_
            )

            next
        }

        # ========================================================
        # Check altcount
        # ========================================================

        if (
            any(!is.finite(x$altcount)) ||
                any(x$altcount < 0)
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                status = "invalid_altcount",
                warning_message = NA_character_
            )

            next
        }

        if (
            sum(
                x$altcount,
                na.rm = TRUE
            ) == 0
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                status = "all_altcount_zero",
                warning_message = NA_character_
            )

            next
        }

        # ========================================================
        # Check sequencing depth
        # ========================================================

        if (
            any(!is.finite(x$dp)) ||
                any(x$dp <= 0)
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                status = "invalid_dp",
                warning_message = NA_character_
            )

            next
        }

        if (
            any(
                x$altcount >
                    x$dp
            )
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                status = "altcount_greater_than_dp",
                warning_message = NA_character_
            )

            next
        }

        # ========================================================
        # Fit negative-binomial models
        # ========================================================

        warning_messages <- character()

        fit_result <- tryCatch(
            withCallingHandlers(
                {
                    null_model <- glm.nb(
                        altcount ~
                            1 +
                            offset(log(dp)),
                        data = x,
                        control = glm.control(
                            maxit = 500
                        )
                    )

                    # Fit the full model with MASS's default dispersion initialization
                    # first.  For sparse or near-Poisson sites, theta.ml() can fail
                    # numerically even though the counts are valid.  In that case retry
                    # the same NB model with a small set of finite starting values; do
                    # not silently switch to a different test family.
                    fit_full_model <- function(init_theta = NULL) {
                        if (is.null(init_theta)) {
                            glm.nb(
                                altcount ~
                                    genotype +
                                    offset(log(dp)),
                                data = x,
                                control = glm.control(
                                    maxit = 500
                                )
                            )
                        } else {
                            glm.nb(
                                altcount ~
                                    genotype +
                                    offset(log(dp)),
                                data = x,
                                init.theta = init_theta,
                                control = glm.control(
                                    maxit = 500
                                )
                            )
                        }
                    }

                    full_model <- tryCatch(
                        fit_full_model(),
                        error = function(e) {
                            NULL
                        }
                    )

                    if (is.null(full_model)) {
                        theta_candidates <- suppressWarnings(
                            as.numeric(
                                c(
                                    null_model$theta,
                                    0.1,
                                    1,
                                    10,
                                    100,
                                    1000,
                                    10000,
                                    1e6
                                )
                            )
                        )
                        theta_candidates <- unique(
                            theta_candidates[
                                is.finite(theta_candidates) &
                                    theta_candidates > 0
                            ]
                        )

                        for (theta_init in theta_candidates) {
                            candidate_model <- tryCatch(
                                fit_full_model(theta_init),
                                error = function(e) {
                                    NULL
                                }
                            )

                            if (!is.null(candidate_model)) {
                                full_model <- candidate_model
                                warning_messages <<- c(
                                    warning_messages,
                                    "full_model_init_theta_fallback"
                                )
                                break
                            }
                        }
                    }

                    if (is.null(full_model)) {
                        stop(
                            "NB full model failed for all dispersion starting values"
                        )
                    }

                    lrt <- anova(
                        null_model,
                        full_model,
                        test = "Chisq"
                    )

                    # Find the P-value column robustly.
                    p_column <- grep(
                        "^Pr\\(",
                        colnames(lrt),
                        value = TRUE
                    )

                    if (length(p_column) == 0) {
                        p_value <- NA_real_
                    } else {
                        p_value <- suppressWarnings(
                            as.numeric(
                                lrt[
                                    2,
                                    p_column[1]
                                ]
                            )
                        )
                    }

                    data.table(
                        coord_handle = h,
                        sex = sx,
                        pvalue = p_value,
                        status = ifelse(
                            is.finite(p_value),
                            "ok",
                            "invalid_lrt_pvalue"
                        ),
                        warning_message = if (
                            length(warning_messages) == 0
                        ) {
                            NA_character_
                        } else {
                            paste(
                                unique(warning_messages),
                                collapse = " | "
                            )
                        }
                    )
                },
                warning = function(w) {
                    warning_messages <<- c(
                        warning_messages,
                        conditionMessage(w)
                    )

                    invokeRestart(
                        "muffleWarning"
                    )
                }
            ),
            error = function(e) {
                data.table(
                    coord_handle = h,
                    sex = sx,
                    pvalue = NA_real_,
                    status = paste0(
                        "error: ",
                        conditionMessage(e)
                    ),
                    warning_message = if (
                        length(warning_messages) == 0
                    ) {
                        NA_character_
                    } else {
                        paste(
                            unique(warning_messages),
                            collapse = " | "
                        )
                    }
                )
            }
        )

        # Directly assign one complete result per site.
        res_list[[i]] <- fit_result
    }

    # ==========================================================
    # Combine results
    # ==========================================================

    res <- rbindlist(
        res_list,
        use.names = TRUE,
        fill = TRUE
    )

    # This should never be needed, but ensures a stable schema.
    if (!"pvalue" %in% colnames(res)) {
        res[
            ,
            pvalue := NA_real_
        ]
    }

    res
}

nb_result <- run_nb_genotype_test_complete(dt_eligible)

cat("NB model status summary:\n")
print(table(nb_result$status, useNA = "ifany"))
cat("NB finite P values: ", sum(is.finite(nb_result$pvalue)), "\n", sep = "")

# Retain a model-specific BH result for reference.  The final overlap filter
# below uses FDR values recalculated over the common valid testing universe.
nb_result[
    ,
    FDR_all_NB_valid := NA_real_
]

nb_result_valid <- nb_result[
    status == "ok" &
        is.finite(pvalue)
]

nb_result_valid[
    ,
    FDR_all_NB_valid := p.adjust(
        pvalue,
        method = "BH"
    ),
    by = sex
]

nb_result[
    nb_result_valid,
    on = .(
        coord_handle,
        sex
    ),
    FDR_all_NB_valid := i.FDR_all_NB_valid
]

fwrite(
    nb_result,
    file.path(
        outdir,
        "NB_genotype_difference_complete3groups.all_sites.with_sexwise_FDR.tsv"
    ),
    sep = "\t"
)


# -------------------------------
# 5.1 Beta-binomial model:
#     cbind(edited reads, unedited reads) ~ genotype
#
# Biological replicates are independent observations.  Replicate labels
# F1/F2/F3 or M1/M2/M3 are therefore not treated as matched blocks across
# genotypes.  The beta-binomial dispersion parameter models extra-binomial
# variation among biological replicates.
# -------------------------------

run_betabinomial_genotype_test_complete <- function(dt_input) {
    required_cols <- c(
        "coord_handle",
        "sex",
        "genotype",
        "rep",
        "altcount",
        "refcount"
    )

    missing_cols <- setdiff(
        required_cols,
        colnames(dt_input)
    )

    if (length(missing_cols) > 0) {
        stop(
            "Missing required columns for beta-binomial model: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    if (nrow(dt_input) == 0) {
        return(
            data.table(
                coord_handle = character(),
                sex = character(),
                pvalue = numeric(),
                LRT_statistic = numeric(),
                LRT_df = numeric(),
                status = character(),
                warning_message = character()
            )
        )
    }

    keys <- unique(
        dt_input[
            ,
            .(
                coord_handle,
                sex
            )
        ]
    )

    res_list <- vector(
        mode = "list",
        length = nrow(keys)
    )

    for (i in seq_len(nrow(keys))) {
        h <- keys$coord_handle[i]
        sx <- keys$sex[i]

        x <- copy(
            dt_input[
                coord_handle == h &
                    sex == sx
            ]
        )

        x[
            ,
            genotype := factor(
                as.character(genotype),
                levels = GENOTYPE_LEVELS
            )
        ]

        # Check the number of independent biological replicates per genotype.
        check_n <- x[
            ,
            .(
                n_reps = uniqueN(
                    rep[
                        !is.na(rep)
                    ]
                )
            ),
            by = .(
                genotype = as.character(genotype)
            )
        ]

        check_n <- merge(
            data.table(
                genotype = GENOTYPE_LEVELS
            ),
            check_n,
            by = "genotype",
            all.x = TRUE,
            sort = FALSE
        )

        check_n[
            is.na(n_reps),
            n_reps := 0L
        ]

        if (
            any(
                check_n$n_reps <
                    MIN_REP_PER_GENOTYPE
            )
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                LRT_statistic = NA_real_,
                LRT_df = NA_real_,
                status = "failed_complete_filter_inside_model",
                warning_message = NA_character_
            )

            next
        }

        # Validate edited and unedited read counts.
        if (
            any(!is.finite(x$altcount)) ||
                any(!is.finite(x$refcount)) ||
                any(x$altcount < 0) ||
                any(x$refcount < 0)
        ) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                LRT_statistic = NA_real_,
                LRT_df = NA_real_,
                status = "invalid_alt_or_ref_count",
                warning_message = NA_character_
            )

            next
        }

        x[
            ,
            `:=`(
                altcount = as.integer(round(altcount)),
                refcount = as.integer(round(refcount))
            )
        ]

        if (any((x$altcount + x$refcount) <= 1)) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                LRT_statistic = NA_real_,
                LRT_df = NA_real_,
                status = "alt_plus_ref_total_le_1",
                warning_message = NA_character_
            )

            next
        }

        if (sum(x$altcount, na.rm = TRUE) == 0) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                LRT_statistic = NA_real_,
                LRT_df = NA_real_,
                status = "all_altcount_zero",
                warning_message = NA_character_
            )

            next
        }

        if (sum(x$refcount, na.rm = TRUE) == 0) {
            res_list[[i]] <- data.table(
                coord_handle = h,
                sex = sx,
                pvalue = NA_real_,
                LRT_statistic = NA_real_,
                LRT_df = NA_real_,
                status = "all_refcount_zero",
                warning_message = NA_character_
            )

            next
        }

        warning_messages <- character()

        fit_result <- tryCatch(
            withCallingHandlers(
                {
                    # The mean/rho beta-binomial parameterization is used.  zero = 2
                    # keeps rho as an intercept-only overdispersion parameter, while
                    # genotype tests the mean editing proportion.
                    fit_one_bb_model <- function(formula_use,
                                                 irho_init = NULL) {
                        family_use <- if (is.null(irho_init)) {
                            VGAM::betabinomial(
                                zero = 2
                            )
                        } else {
                            VGAM::betabinomial(
                                zero = 2,
                                irho = irho_init
                            )
                        }

                        VGAM::vglm(
                            formula_use,
                            family = family_use,
                            data = x,
                            control = VGAM::vglm.control(
                                maxit = 500,
                                trace = FALSE
                            )
                        )
                    }

                    # VGAM recommends supplying an initial rho when the automatic
                    # initialization fails.  Both nested models are fitted using the
                    # same starting rho in each attempt.
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

                    null_model <- NULL
                    full_model <- NULL
                    used_irho <- NULL
                    last_fit_error <- NA_character_

                    for (irho_try in irho_candidates) {
                        candidate_pair <- tryCatch(
                            {
                                candidate_null <- fit_one_bb_model(
                                    cbind(altcount, refcount) ~ 1,
                                    irho_init = irho_try
                                )

                                candidate_full <- fit_one_bb_model(
                                    cbind(altcount, refcount) ~ genotype,
                                    irho_init = irho_try
                                )

                                candidate_ll_null <- logLik(
                                    candidate_null
                                )
                                candidate_ll_full <- logLik(
                                    candidate_full
                                )

                                if (
                                    !is.finite(
                                        suppressWarnings(
                                            as.numeric(candidate_ll_null)
                                        )
                                    ) ||
                                        !is.finite(
                                            suppressWarnings(
                                                as.numeric(candidate_ll_full)
                                            )
                                        )
                                ) {
                                    stop(
                                        "Non-finite log-likelihood"
                                    )
                                }

                                list(
                                    null = candidate_null,
                                    full = candidate_full
                                )
                            },
                            error = function(e) {
                                last_fit_error <<- conditionMessage(e)
                                NULL
                            }
                        )

                        if (!is.null(candidate_pair)) {
                            null_model <- candidate_pair$null
                            full_model <- candidate_pair$full
                            used_irho <- irho_try
                            break
                        }
                    }

                    if (
                        is.null(null_model) ||
                            is.null(full_model)
                    ) {
                        stop(
                            "Beta-binomial null/full model fitting failed for all rho initializations: ",
                            last_fit_error
                        )
                    }

                    if (!is.null(used_irho)) {
                        warning_messages <<- c(
                            warning_messages,
                            paste0(
                                "beta_binomial_irho_fallback=",
                                format(
                                    used_irho,
                                    scientific = TRUE
                                )
                            )
                        )
                    }

                    ll_null <- logLik(null_model)
                    ll_full <- logLik(full_model)

                    loglik_null <- suppressWarnings(
                        as.numeric(ll_null)
                    )
                    loglik_full <- suppressWarnings(
                        as.numeric(ll_full)
                    )

                    df_null <- suppressWarnings(
                        as.numeric(
                            attr(
                                ll_null,
                                "df"
                            )
                        )
                    )
                    df_full <- suppressWarnings(
                        as.numeric(
                            attr(
                                ll_full,
                                "df"
                            )
                        )
                    )

                    if (
                        length(df_null) != 1L ||
                            !is.finite(df_null)
                    ) {
                        df_null <- length(
                            stats::coef(null_model)
                        )
                    }

                    if (
                        length(df_full) != 1L ||
                            !is.finite(df_full)
                    ) {
                        df_full <- length(
                            stats::coef(full_model)
                        )
                    }

                    lrt_df <- df_full - df_null
                    lrt_statistic <- 2 * (
                        loglik_full -
                            loglik_null
                    )

                    if (
                        !is.finite(loglik_null) ||
                            !is.finite(loglik_full) ||
                            !is.finite(lrt_statistic) ||
                            !is.finite(lrt_df) ||
                            lrt_df <= 0
                    ) {
                        stop(
                            "Invalid beta-binomial likelihood-ratio statistic or degrees of freedom"
                        )
                    }

                    # Tiny negative values can occur through numerical rounding.
                    if (lrt_statistic < -1e-7) {
                        stop(
                            "Beta-binomial full model has a lower log-likelihood than the null model"
                        )
                    }

                    lrt_statistic <- max(
                        lrt_statistic,
                        0
                    )

                    p_value <- stats::pchisq(
                        lrt_statistic,
                        df = lrt_df,
                        lower.tail = FALSE
                    )

                    data.table(
                        coord_handle = h,
                        sex = sx,
                        pvalue = p_value,
                        LRT_statistic = lrt_statistic,
                        LRT_df = lrt_df,
                        status = ifelse(
                            is.finite(p_value),
                            "ok",
                            "invalid_lrt_pvalue"
                        ),
                        warning_message = if (
                            length(warning_messages) == 0
                        ) {
                            NA_character_
                        } else {
                            paste(
                                unique(warning_messages),
                                collapse = " | "
                            )
                        }
                    )
                },
                warning = function(w) {
                    warning_messages <<- c(
                        warning_messages,
                        conditionMessage(w)
                    )

                    invokeRestart(
                        "muffleWarning"
                    )
                }
            ),
            error = function(e) {
                data.table(
                    coord_handle = h,
                    sex = sx,
                    pvalue = NA_real_,
                    LRT_statistic = NA_real_,
                    LRT_df = NA_real_,
                    status = paste0(
                        "error: ",
                        conditionMessage(e)
                    ),
                    warning_message = if (
                        length(warning_messages) == 0
                    ) {
                        NA_character_
                    } else {
                        paste(
                            unique(warning_messages),
                            collapse = " | "
                        )
                    }
                )
            }
        )

        res_list[[i]] <- fit_result
    }

    res <- rbindlist(
        res_list,
        use.names = TRUE,
        fill = TRUE
    )

    if (!"pvalue" %in% colnames(res)) {
        res[
            ,
            pvalue := NA_real_
        ]
    }

    res
}


bb_result <- run_betabinomial_genotype_test_complete(
    dt_eligible
)

cat("Beta-binomial model status summary:\n")
print(table(bb_result$status, useNA = "ifany"))
cat(
    "Beta-binomial finite P values: ",
    sum(is.finite(bb_result$pvalue)),
    "\n",
    sep = ""
)

bb_result[
    ,
    FDR_all_BB_valid := NA_real_
]

bb_result_valid <- bb_result[
    status == "ok" &
        is.finite(pvalue)
]

bb_result_valid[
    ,
    FDR_all_BB_valid := p.adjust(
        pvalue,
        method = "BH"
    ),
    by = sex
]

bb_result[
    bb_result_valid,
    on = .(
        coord_handle,
        sex
    ),
    FDR_all_BB_valid := i.FDR_all_BB_valid
]

fwrite(
    bb_result,
    file.path(
        outdir,
        "BetaBinomial_genotype_difference_complete3groups.all_sites.with_sexwise_FDR.tsv"
    ),
    sep = "\t"
)


# -------------------------------
# 5.2 Put both models on the same valid testing universe
#     and perform BH correction separately within each sex.
# -------------------------------

nb_result_named <- copy(nb_result)
setnames(
    nb_result_named,
    old = c(
        "pvalue",
        "status",
        "warning_message",
        "FDR_all_NB_valid"
    ),
    new = c(
        "NB_pvalue",
        "NB_status",
        "NB_warning_message",
        "NB_FDR_all_valid"
    )
)

bb_result_named <- copy(bb_result)
setnames(
    bb_result_named,
    old = c(
        "pvalue",
        "LRT_statistic",
        "LRT_df",
        "status",
        "warning_message",
        "FDR_all_BB_valid"
    ),
    new = c(
        "BB_pvalue",
        "BB_LRT_statistic",
        "BB_LRT_df",
        "BB_status",
        "BB_warning_message",
        "BB_FDR_all_valid"
    )
)

model_result_all <- merge(
    nb_result_named,
    bb_result_named,
    by = c(
        "coord_handle",
        "sex"
    ),
    all = TRUE
)

model_result_all[
    ,
    common_valid := (
        NB_status == "ok" &
            BB_status == "ok" &
            is.finite(NB_pvalue) &
            is.finite(BB_pvalue)
    )
]

model_result_common <- model_result_all[
    common_valid == TRUE
]

# Both P-value families are corrected over the identical common-valid set.
model_result_common[
    ,
    `:=`(
        NB_FDR = p.adjust(
            NB_pvalue,
            method = "BH"
        ),
        BB_FDR = p.adjust(
            BB_pvalue,
            method = "BH"
        )
    ),
    by = sex
]

model_result_all[
    ,
    `:=`(
        NB_FDR = NA_real_,
        BB_FDR = NA_real_
    )
]

model_result_all[
    model_result_common,
    on = .(
        coord_handle,
        sex
    ),
    `:=`(
        NB_FDR = i.NB_FDR,
        BB_FDR = i.BB_FDR
    )
]

model_result_all[
    ,
    `:=`(
        NB_significant = (
            common_valid == TRUE &
                is.finite(NB_FDR) &
                NB_FDR < FDR_CUTOFF
        ),
        BB_significant = (
            common_valid == TRUE &
                is.finite(BB_FDR) &
                BB_FDR < FDR_CUTOFF
        )
    )
]

model_result_all[
    ,
    `:=`(
        both_significant = (
            NB_significant == TRUE &
                BB_significant == TRUE
        ),
        NB_only_significant = (
            NB_significant == TRUE &
                BB_significant == FALSE
        ),
        BB_only_significant = (
            NB_significant == FALSE &
                BB_significant == TRUE
        ),
        either_significant = (
            NB_significant == TRUE |
                BB_significant == TRUE
        )
    )
]

fwrite(
    model_result_all,
    file.path(
        outdir,
        "NB_BetaBinomial_model_comparison.all_sites.tsv"
    ),
    sep = "\t"
)

fwrite(
    model_result_all[
        common_valid == TRUE
    ],
    file.path(
        outdir,
        "NB_BetaBinomial_common_valid_tests.with_sexwise_FDR.tsv"
    ),
    sep = "\t"
)

fwrite(
    model_result_all[
        both_significant == TRUE
    ],
    file.path(
        outdir,
        "NB_BetaBinomial_both_significant.before_effect_filter.tsv"
    ),
    sep = "\t"
)

summarize_model_overlap <- function(x) {
    x[
        ,
        .(
            n_common_valid = .N,
            n_NB_significant = sum(
                NB_significant,
                na.rm = TRUE
            ),
            n_BB_significant = sum(
                BB_significant,
                na.rm = TRUE
            ),
            n_both_significant = sum(
                both_significant,
                na.rm = TRUE
            ),
            n_NB_only = sum(
                NB_only_significant,
                na.rm = TRUE
            ),
            n_BB_only = sum(
                BB_only_significant,
                na.rm = TRUE
            ),
            n_significant_union = sum(
                either_significant,
                na.rm = TRUE
            )
        )
    ]
}

overlap_summary_by_sex <- model_result_all[
    common_valid == TRUE,
    summarize_model_overlap(.SD),
    by = sex
]

overlap_summary_total <- summarize_model_overlap(
    model_result_all[
        common_valid == TRUE
    ]
)
overlap_summary_total[
    ,
    sex := "All"
]
setcolorder(
    overlap_summary_total,
    c(
        "sex",
        setdiff(
            colnames(overlap_summary_total),
            "sex"
        )
    )
)

model_overlap_summary <- rbindlist(
    list(
        overlap_summary_by_sex,
        overlap_summary_total
    ),
    use.names = TRUE,
    fill = TRUE
)

cat("NB versus beta-binomial overlap summary:\n")
print(model_overlap_summary)

fwrite(
    model_overlap_summary,
    file.path(
        outdir,
        "NB_BetaBinomial_significance_overlap_summary.tsv"
    ),
    sep = "\t"
)


# -------------------------------
# 6. Calculate AG/WBG/AS mean editing levels and retain:
#    1) significant by both NB and beta-binomial;
#    2) at least one group with mean editing level > 0.10;
#    3) at least one pairwise absolute mean difference > 0.05.
# -------------------------------

level_summary <- dt_eligible[
    ,
    .(
        mean_level = mean(
            editing_level,
            na.rm = TRUE
        ),
        median_level = median(
            editing_level,
            na.rm = TRUE
        ),
        mean_altcount = mean(
            altcount,
            na.rm = TRUE
        ),
        mean_refcount = mean(
            refcount,
            na.rm = TRUE
        ),
        mean_dp = mean(
            dp,
            na.rm = TRUE
        ),
        n_records = .N,
        n_reps = uniqueN(rep)
    ),
    by = .(
        coord_handle,
        sex,
        genotype
    )
]

level_wide <- dcast(
    level_summary,
    coord_handle + sex ~ genotype,
    value.var = "mean_level"
)

for (g in GENOTYPE_LEVELS) {
    if (!g %in% colnames(level_wide)) {
        level_wide[
            ,
            (g) := NA_real_
        ]
    }
}

level_wide <- merge(
    level_wide,
    model_result_all[
        ,
        .(
            coord_handle,
            sex,
            NB_pvalue,
            NB_FDR_all_valid,
            NB_FDR,
            NB_status,
            NB_warning_message,
            BB_pvalue,
            BB_FDR_all_valid,
            BB_FDR,
            BB_LRT_statistic,
            BB_LRT_df,
            BB_status,
            BB_warning_message,
            common_valid,
            NB_significant,
            BB_significant,
            both_significant,
            NB_only_significant,
            BB_only_significant,
            either_significant
        )
    ],
    by = c(
        "coord_handle",
        "sex"
    ),
    all.x = TRUE
)

rep_info <- copy(rep_wide)
setnames(
    rep_info,
    c(
        "AG",
        "WBG",
        "AS"
    ),
    c(
        "nrep_AG",
        "nrep_WBG",
        "nrep_AS"
    )
)

level_wide <- merge(
    level_wide,
    rep_info[
        ,
        .(
            coord_handle,
            sex,
            nrep_AG,
            nrep_WBG,
            nrep_AS,
            complete_3genotypes
        )
    ],
    by = c(
        "coord_handle",
        "sex"
    ),
    all.x = TRUE
)

level_wide[
    ,
    pass_level_cutoff := rowSums(
        .SD > LEVEL_CUTOFF,
        na.rm = TRUE
    ) > 0,
    .SDcols = c(
        "AG",
        "WBG",
        "AS"
    )
]

level_wide[
    ,
    `:=`(
        abs_diff_AG_WBG = abs(
            AG -
                WBG
        ),
        abs_diff_AS_WBG = abs(
            AS -
                WBG
        ),
        abs_diff_AG_AS = abs(
            AG -
                AS
        )
    )
]

level_wide[
    ,
    max_pairwise_diff := pmax(
        abs_diff_AG_WBG,
        abs_diff_AS_WBG,
        abs_diff_AG_AS,
        na.rm = FALSE
    )
]

level_wide[
    ,
    `:=`(
        pass_AG_WBG_diff = (
            is.finite(abs_diff_AG_WBG) &
                abs_diff_AG_WBG > PAIRWISE_DIFF_CUTOFF
        ),
        pass_AS_WBG_diff = (
            is.finite(abs_diff_AS_WBG) &
                abs_diff_AS_WBG > PAIRWISE_DIFF_CUTOFF
        ),
        pass_AG_AS_diff = (
            is.finite(abs_diff_AG_AS) &
                abs_diff_AG_AS > PAIRWISE_DIFF_CUTOFF
        ),
        pass_pairwise_diff_cutoff = (
            is.finite(max_pairwise_diff) &
                max_pairwise_diff > PAIRWISE_DIFF_CUTOFF
        )
    )
]

fwrite(
    level_wide,
    file.path(
        outdir,
        "NB_BetaBinomial_complete3groups.with_mean_levels_and_effects.all_sites.tsv"
    ),
    sep = "\t"
)

sig_level_wide <- level_wide[
    is.finite(AG) &
        is.finite(WBG) &
        is.finite(AS) &
        common_valid == TRUE &
        NB_significant == TRUE &
        BB_significant == TRUE &
        pass_level_cutoff == TRUE &
        pass_pairwise_diff_cutoff == TRUE
]

cat(
    "Final sites significant in both models and passing level/effect cutoffs:\n"
)
print(
    sig_level_wide[
        ,
        .(
            n_sites = uniqueN(
                coord_handle
            )
        ),
        by = sex
    ]
)

effect_filter_summary <- level_wide[
    common_valid == TRUE,
    .(
        n_common_valid = .N,
        n_both_significant = sum(
            both_significant,
            na.rm = TRUE
        ),
        n_both_sig_pass_level = sum(
            both_significant &
                pass_level_cutoff,
            na.rm = TRUE
        ),
        n_both_sig_pass_pairwise_diff = sum(
            both_significant &
                pass_pairwise_diff_cutoff,
            na.rm = TRUE
        ),
        n_final = sum(
            both_significant &
                pass_level_cutoff &
                pass_pairwise_diff_cutoff,
            na.rm = TRUE
        )
    ),
    by = sex
]

cat("Effect-size filtering summary:\n")
print(effect_filter_summary)

fwrite(
    effect_filter_summary,
    file.path(
        outdir,
        "NB_BetaBinomial_effect_filter_summary.tsv"
    ),
    sep = "\t"
)

fwrite(
    sig_level_wide,
    file.path(
        outdir,
        paste0(
            "NB_BetaBinomial_both_FDR005.level_gt",
            LEVEL_CUTOFF,
            ".any_pairwise_diff_gt",
            PAIRWISE_DIFF_CUTOFF,
            ".with_mean_levels.tsv"
        )
    ),
    sep = "\t"
)

if (nrow(sig_level_wide) == 0L) {
    stop(
        "No sites passed the joint NB + beta-binomial significance and effect-size filters."
    )
}


# -------------------------------
# 7.PLOT VIOLINE
# -------------------------------

plot_top_clusters_violin_by_sex <- function(
  sig_level_wide,
  sex_use = "M",
  top_n_cluster = 4,
  outdir = ".",
  width = 10,
  height = 3.2
) {
    dt_plot <- as.data.table(copy(sig_level_wide))

    dt_plot <- dt_plot[sex == sex_use]

    dt_plot <- dt_plot[!is.na(AG) & !is.na(WBG) & !is.na(AS)]

    if (nrow(dt_plot) == 0L) {
        warning(
            "No final jointly significant sites are available for sex = ",
            sex_use,
            "."
        )

        if (!"shape_trend" %in% colnames(dt_plot)) {
            dt_plot[
                ,
                shape_trend := character()
            ]
        }

        empty_top <- data.table(
            sex = character(),
            shape_trend = character(),
            n_sites = integer(),
            cluster_rank = integer(),
            cluster_id = character()
        )

        fwrite(
            empty_top,
            file.path(
                outdir,
                paste0(
                    "NB_BB_overlap_complete3groups_top",
                    top_n_cluster,
                    "_clusters_",
                    sex_use,
                    ".tsv"
                )
            ),
            sep = "\t"
        )

        fwrite(
            dt_plot,
            file.path(
                outdir,
                paste0(
                    "dt_plot_NB_BB_overlap_complete3groups_top",
                    top_n_cluster,
                    "_clusters_plot_data_",
                    sex_use,
                    ".tsv"
                )
            ),
            sep = "\t"
        )

        return(
            list(
                plot = NULL,
                top_clusters = empty_top,
                plot_data = data.table(),
                dt_plot = dt_plot
            )
        )
    }

    if (!"shape_trend" %in% colnames(dt_plot)) {
        dt_plot[, shape_trend := fifelse(
            AG < WBG & WBG < AS,
            "AG < WBG < AS",
            fifelse(
                AG > WBG & WBG > AS,
                "AG > WBG > AS",
                fifelse(
                    AG > WBG & WBG < AS,
                    "AG > WBG < AS",
                    fifelse(
                        AG < WBG & WBG >= AS,
                        "AG < WBG >= AS",
                        fifelse(
                            AG == WBG & WBG < AS,
                            "AG = WBG < AS",
                            fifelse(
                                AG == WBG & WBG > AS,
                                "AG = WBG > AS",
                                fifelse(
                                    AG < WBG & WBG == AS,
                                    "AG < WBG = AS",
                                    fifelse(
                                        AG > WBG & WBG == AS,
                                        "AG > WBG = AS",
                                        "other"
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )]
    }

    dt_plot <- dt_plot[!is.na(shape_trend)]

    cluster_count <- dt_plot[, .(
        n_sites = uniqueN(coord_handle)
    ), by = .(sex, shape_trend)]

    setorder(cluster_count, -n_sites)
    cluster_count[, cluster_rank := seq_len(.N)]
    top_clusters <- cluster_count[cluster_rank <= top_n_cluster]
    top_clusters[, cluster_id := paste0("C", cluster_rank)]

    print(top_clusters)

    dt_top <- merge(
        dt_plot,
        top_clusters[, .(sex, shape_trend, n_sites, cluster_rank, cluster_id)],
        by = c("sex", "shape_trend"),
        all = FALSE
    )

    dt_top[, facet_label := paste0(cluster_id, "_", sex, " (n = ", n_sites, ")")]
    dt_top[, facet_label := factor(
        facet_label,
        levels = top_clusters[order(cluster_rank), paste0(cluster_id, "_", sex, " (n = ", n_sites, ")")]
    )]

    long_plot <- melt(
        dt_top,
        measure.vars = c(
            "AG",
            "WBG",
            "AS"
        ),
        variable.name = "genotype",
        value.name = "Lm"
    )

    long_plot[, genotype := factor(genotype, levels = c("AG", "WBG", "AS"))]

    mean_df <- long_plot[, .(
        mean_Lm = mean(Lm, na.rm = TRUE)
    ), by = .(sex, shape_trend, cluster_id, cluster_rank, n_sites, facet_label, genotype)]

    cols_fill <- c(
        "AS"  = "#D95F02",
        "WBG" = "#E6A23C",
        "AG"  = "#D8B44A"
    )

    p <- ggplot(long_plot, aes(x = genotype, y = Lm)) +
        geom_line(
            aes(group = coord_handle),
            color = "grey75",
            alpha = 0.55,
            linewidth = 0.25
        ) +
        geom_point(
            aes(group = coord_handle),
            color = "grey75",
            alpha = 0.55,
            size = 0.45
        ) +
        geom_violin(
            aes(fill = genotype),
            width = 0.75,
            trim = FALSE,
            alpha = 0.9,
            color = NA
        ) +
        geom_boxplot(
            aes(fill = genotype),
            width = 0.11,
            outlier.shape = NA,
            color = "black",
            linewidth = 0.3,
            alpha = 0.95
        ) +
        geom_line(
            data = mean_df,
            aes(x = genotype, y = mean_Lm, group = 1),
            inherit.aes = FALSE,
            color = "grey35",
            linewidth = 0.6
        ) +
        geom_point(
            data = mean_df,
            aes(x = genotype, y = mean_Lm, fill = genotype),
            inherit.aes = FALSE,
            shape = 24,
            size = 2.3,
            color = "grey35"
        ) +
        facet_wrap(~facet_label, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = cols_fill) +
        scale_x_discrete(
            labels = c(
                "AS" = expression(Adar^S),
                "WBG" = "WT",
                "AG" = expression(Adar^G)
            )
        ) +
        labs(
            x = "Genotype",
            y = expression(L[m])
        ) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = "white", color = NA),
            strip.text = element_text(size = 10, face = "bold"),
            panel.grid = element_blank(),
            axis.title = element_text(size = 11),
            axis.text = element_text(size = 9),
            legend.position = "none",
            plot.margin = margin(5, 5, 5, 5)
        ) +
        theme_prism(
            palette = "candy_bright",
            base_fontface = "plain",
            base_size = 14,
            base_line_size = 0,
            axis_text_angle = 0
        ) +
        theme(legend.position = "none", axis.text.x = element_text(size = 12))

    print(p)

    ggsave(
        filename = file.path(outdir, paste0("NB_BB_overlap_complete3groups_top", top_n_cluster, "_clusters_violin_lineplot_", sex_use, ".pdf")),
        plot = p,
        width = width,
        height = height
    )

    fwrite(
        top_clusters,
        file.path(outdir, paste0("NB_BB_overlap_complete3groups_top", top_n_cluster, "_clusters_", sex_use, ".tsv")),
        sep = "\t"
    )

    fwrite(
        long_plot,
        file.path(outdir, paste0("NB_BB_overlap_complete3groups_top", top_n_cluster, "_clusters_plot_data_", sex_use, ".tsv")),
        sep = "\t"
    )

    fwrite(
        dt_plot,
        file.path(outdir, paste0("dt_plot_NB_BB_overlap_complete3groups_top", top_n_cluster, "_clusters_plot_data_", sex_use, ".tsv")),
        sep = "\t"
    )

    return(list(
        plot = p,
        top_clusters = top_clusters,
        plot_data = long_plot,
        dt_plot = dt_plot
    ))
}

res_M <- plot_top_clusters_violin_by_sex(
    sig_level_wide = sig_level_wide,
    sex_use = "M",
    top_n_cluster = TOP_N_CLUSTERS,
    outdir = outdir,
    width = 12,
    height = 3.2
)

res_F <- plot_top_clusters_violin_by_sex(
    sig_level_wide = sig_level_wide,
    sex_use = "F",
    top_n_cluster = TOP_N_CLUSTERS,
    outdir = outdir,
    width = 12,
    height = 3.2
)


# -------------------------------
# 8. UPSET PLOT
# -------------------------------
dt_plot_f <- fread(
    file.path(
        outdir,
        paste0(
            "dt_plot_NB_BB_overlap_complete3groups_top",
            TOP_N_CLUSTERS,
            "_clusters_plot_data_F.tsv"
        )
    )
)
dt_plot_m <- fread(
    file.path(
        outdir,
        paste0(
            "dt_plot_NB_BB_overlap_complete3groups_top",
            TOP_N_CLUSTERS,
            "_clusters_plot_data_M.tsv"
        )
    )
)
dt_plot <- rbind(dt_plot_f, dt_plot_m)
sig_for_downstream <- dt_plot[
    ,
    .(
        coord_handle,
        sex,
        AG,
        WBG,
        AS,
        NB_pvalue,
        NB_FDR,
        BB_pvalue,
        BB_FDR,
        abs_diff_AG_WBG,
        abs_diff_AS_WBG,
        abs_diff_AG_AS,
        max_pairwise_diff,
        pass_pairwise_diff_cutoff,
        shape_trend,
        nrep_AG,
        nrep_WBG,
        nrep_AS
    )
]

fwrite(
    sig_for_downstream,
    file.path(outdir, "NB_BetaBinomial_both_significant_effect_filtered.downstream_table.tsv"),
    sep = "\t"
)

combined <- sig_for_downstream[, .(coord_handle, sex, shape_trend)]
combined[, trend_label := paste(sex, shape_trend, sep = "_")]

upset_matrix <- dcast(
    combined,
    coord_handle ~ trend_label,
    fun.aggregate = length,
    fill = 0
)

upset_matrix[, (2:ncol(upset_matrix)) := lapply(.SD, function(x) ifelse(x > 0, 1, 0)),
    .SDcols = 2:ncol(upset_matrix)
]

pdf(
    file = file.path(outdir, "Upset_NB_BetaBinomial_overlap_clusters.pdf"),
    width = 15,
    height = 9
)

upset_plot <- upset(
    upset_matrix,
    sets = names(upset_matrix)[-1],
    nsets = ncol(upset_matrix) - 1,
    nintersects = 20,
    order.by = "freq",
    keep.order = FALSE,
    mb.ratio = c(0.6, 0.4),
    text.scale = 3.5,
    number.angles = 0,
    point.size = 5,
    line.size = 2,
    mainbar.y.label = "Site Intersections",
    sets.x.label = "Sites Per Cluster"
)
print(upset_plot, newpage = FALSE)

dev.off()


# -------------------------------
# 10. AdarG versus AdarS:
#     absolute Lm change relative to WBG
# -------------------------------

plot_lm_delta_violin <- function() {
    # Use all sites in sig_level_wide.  This is the complete set of sites with
    # significant Lm differences across all trend clusters; do not use res_M$dt_plot
    # or res_F$dt_plot here because those objects contain only the top clusters.
    if (!exists("sig_level_wide")) {
        stop("sig_level_wide is required for the AdarG/AdarS delta violin plot.")
    }

    if (!exists("anno_gene")) {
        stop("anno_gene is required to restrict the analysis to annotated genes.")
    }

    lm_sig_all_sites <- as.data.table(copy(sig_level_wide))
    lm_sig_all_sites <- lm_sig_all_sites[
        is.finite(AG) &
            is.finite(WBG) &
            is.finite(AS) &
            is.finite(NB_FDR) &
            NB_FDR < FDR_CUTOFF &
            is.finite(BB_FDR) &
            BB_FDR < FDR_CUTOFF &
            pass_pairwise_diff_cutoff == TRUE
    ]

    site_gene_map <- unique(
        anno_gene[
            !is.na(coord_handle) &
                !is.na(gene) &
                gene != "",
            .(coord_handle, gene)
        ]
    )

    annotated_sig_coord <- unique(
        site_gene_map[
            coord_handle %chin% lm_sig_all_sites$coord_handle,
            coord_handle
        ]
    )

    lm_sig_all_sites <- lm_sig_all_sites[
        coord_handle %chin% annotated_sig_coord
    ]

    if (nrow(lm_sig_all_sites) == 0) {
        stop(
            "No annotated significant Lm sites are available for the delta violin plot."
        )
    }

    lm_sig_site_gene_map <- site_gene_map[
        coord_handle %chin% unique(lm_sig_all_sites$coord_handle)
    ]

    cat(
        "Annotated significant Lm sites used for AdarG/AdarS delta violin:\n"
    )
    print(lm_sig_all_sites[, .(n_sites = uniqueN(coord_handle)), by = sex])
    cat(
        "Annotated genes represented in these sites: ",
        uniqueN(lm_sig_site_gene_map$gene),
        "\n",
        sep = ""
    )

    fwrite(
        lm_sig_all_sites,
        file.path(
            outdir,
            "NB_BetaBinomial_overlap_all_clusters_Lm_sites_for_AdarG_AdarS_delta.tsv"
        ),
        sep = "\t"
    )

    fwrite(
        lm_sig_site_gene_map,
        file.path(
            outdir,
            "NB_BetaBinomial_overlap_all_clusters_Lm_site_gene_map.tsv"
        ),
        sep = "\t"
    )

    lm_delta_wide <- lm_sig_all_sites[
        ,
        .(
            coord_handle,
            sex,
            AdarG_abs_delta = abs(AG - WBG),
            AdarS_abs_delta = abs(AS - WBG)
        )
    ]

    lm_delta_wide <- lm_delta_wide[
        is.finite(AdarG_abs_delta) &
            is.finite(AdarS_abs_delta)
    ]

    # Paired Wilcoxon tests use the untruncated values.  Each pair is the same
    # coord_handle within one sex, so the comparison is site-paired.
    safe_wilcox_paired_lm_delta <- function(x, y) {
        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]
        y <- y[ok]

        if (length(x) < 2) {
            return(NA_real_)
        }

        tryCatch(
            suppressWarnings(
                wilcox.test(
                    x,
                    y,
                    paired = TRUE,
                    exact = FALSE
                )$p.value
            ),
            error = function(e) NA_real_
        )
    }

    lm_delta_test <- lm_delta_wide[
        ,
        .(
            n_pairs = .N,
            p_value = safe_wilcox_paired_lm_delta(
                AdarG_abs_delta,
                AdarS_abs_delta
            )
        ),
        by = sex
    ]

    lm_delta_test[, p_adjust_BH := p.adjust(p_value, method = "BH")]

    fwrite(
        lm_delta_test,
        file.path(
            outdir,
            "NB_BetaBinomial_overlap_AdarG_AdarS_absolute_Lm_change_vs_WBG_paired_wilcoxon.tsv"
        ),
        sep = "\t"
    )

    lm_delta_plot <- melt(
        lm_delta_wide,
        id.vars = c("coord_handle", "sex"),
        measure.vars = c("AdarS_abs_delta", "AdarG_abs_delta"),
        variable.name = "genotype",
        value.name = "delta_raw"
    )

    lm_delta_plot[
        ,
        genotype := factor(
            genotype,
            levels = c("AdarS_abs_delta", "AdarG_abs_delta"),
            labels = c("AdarS", "AdarG")
        )
    ]

    lm_delta_plot[
        ,
        sex := factor(
            sex,
            levels = c("F", "M"),
            labels = c("Female", "Male")
        )
    ]

    LM_DELTA_PLOT_CAP <- 0.4

    # Only the displayed values are capped.  Raw values remain in delta_raw for
    # the Wilcoxon test and in the summary table below.
    lm_delta_plot[
        ,
        value_plot := pmin(delta_raw, LM_DELTA_PLOT_CAP)
    ]

    lm_delta_summary <- lm_delta_plot[
        is.finite(delta_raw),
        .(
            n = .N,
            mean_abs_delta = mean(delta_raw),
            sd_abs_delta = sd(delta_raw)
        ),
        by = .(sex, genotype)
    ]

    lm_delta_summary[
        ,
        `:=`(
            label = sprintf(
                "italic(n) == %d ~ ',' ~ %.3f %%+-%% %.3f",
                n,
                mean_abs_delta,
                sd_abs_delta
            ),
            label_x = fifelse(
                as.character(genotype) == "AdarS",
                0.75,
                2.25
            ),
            label_hjust = fifelse(
                as.character(genotype) == "AdarS",
                0,
                1
            )
        )
    ]

    format_p_lm_delta <- function(p) {
        if (!is.finite(p)) {
            return("NA")
        }

        if (p < 0.001) {
            return(formatC(p, format = "e", digits = 2))
        }

        formatC(p, format = "f", digits = 3)
    }

    lm_delta_annotation <- data.table(
        sex = factor(
            c("Female", "Male"),
            levels = c("Female", "Male")
        ),
        x = 1.5,
        y = 0.30,
        label = c(
            paste0(
                "italic(P) == ",
                format_p_lm_delta(
                    lm_delta_test[sex == "F", p_value][1]
                )
            ),
            paste0(
                "italic(P) == ",
                format_p_lm_delta(
                    lm_delta_test[sex == "M", p_value][1]
                )
            )
        )
    )

    cols_lm_delta <- c(
        "AdarS" = "#F8766D",
        "AdarG" = "#00BFC4"
    )

    p_lm_delta <- ggplot(
        lm_delta_plot,
        aes(
            x = genotype,
            y = value_plot,
            fill = genotype
        )
    ) +
        geom_violin(
            trim = TRUE,
            color = NA,
            alpha = 0.75
        ) +
        geom_boxplot(
            width = 0.18,
            outlier.shape = NA,
            alpha = 0.8
        ) +
        facet_wrap(
            ~sex,
            nrow = 1
        ) +
        geom_text(
            data = lm_delta_summary,
            aes(
                x = label_x,
                y = 0.385,
                label = label,
                hjust = label_hjust
            ),
            inherit.aes = FALSE,
            parse = TRUE,
            vjust = 1,
            size = 3.3
        ) +
        geom_text(
            data = lm_delta_annotation,
            aes(
                x = x,
                y = y,
                label = label
            ),
            inherit.aes = FALSE,
            parse = TRUE,
            vjust = 0,
            size = 4
        ) +
        scale_fill_manual(values = cols_lm_delta) +
        scale_y_continuous(
            limits = c(0, LM_DELTA_PLOT_CAP),
            breaks = seq(0, LM_DELTA_PLOT_CAP, by = 0.1),
            expand = expansion(mult = c(0, 0.02))
        ) +
        theme_classic(
            base_size = 14
        ) +
        labs(
            x = NULL,
            y = expression("|" * Delta * L[m] * " vs WBG" * "|")
        ) +
        theme(
            legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_text(size = 14)
        )

    lm_delta_pdf <- file.path(
        outdir,
        "NB_BetaBinomial_overlap_AdarG_AdarS_absolute_Lm_change_vs_WBG_violin.pdf"
    )

    grDevices::pdf(
        file = lm_delta_pdf,
        width = 9,
        height = 5.5,
        onefile = TRUE,
        useDingbats = FALSE,
        compress = TRUE
    )
    print(p_lm_delta)
    grDevices::dev.off()

    fwrite(
        lm_delta_plot,
        file.path(
            outdir,
            "NB_BetaBinomial_overlap_AdarG_AdarS_absolute_Lm_change_vs_WBG_violin_data.tsv"
        ),
        sep = "\t"
    )

    cat(
        "Saved AdarG/AdarS absolute Lm-change violin plot:\n",
        lm_delta_pdf,
        "\n"
    )
}


fwrite(
    upset_matrix,
    file.path(outdir, "NB_BetaBinomial_both_significant_effect_filtered.upset_matrix.tsv"),
    sep = "\t"
)


# -------------------------------
# 9. GO PLOT
# -------------------------------

# Export the model context for the separate GO analysis.
save(list = ls(envir = .GlobalEnv), file = file.path(outdir, "editing_models.RData"), envir = .GlobalEnv)
anno_gene <- unique(dt[!is.na(gene) & grepl("^FBgn", gene), .(coord_handle, gene)])
plot_lm_delta_violin()
