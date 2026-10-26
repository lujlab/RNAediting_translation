# ============================================================
# Combined plotting script
#   1) RNA digest inosine quantification
#   2) Panel D: dual-luciferase assay
#   3) Panel E: WB quantification normalized to 0 mM
#   4) Panel F: S2-cell Renilla activity
#
# Statistics:
#   - Digest: two-sided Pearson correlation + linear regression
#             on dose means (linear concentration-response question).
#   - D/E/F: one-sided Spearman correlation (alternative = 'less'),
#          using the exact null distribution from pspearman and the original ITP concentration.
# ============================================================

required_packages <- c("readxl", "dplyr", "stringr", "tibble", "readr", "ggplot2", "scales",
    "data.table", "tidyr", "gridExtra", "pspearman")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE,
    FUN.VALUE = logical(1))]
if (length(missing_packages) > 0L) {
    stop("Please install required packages first: install.packages(c(", paste(sprintf("\"%s\"",
        missing_packages), collapse = ", "), "))")
}

suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(stringr)
    library(tibble)
    library(readr)
    library(ggplot2)
    library(scales)
    library(data.table)
    library(tidyr)
    library(gridExtra)
})


# ============================================================
# Exact Spearman helper and P-value formatter
# pspearman provides the precomputed exact null distribution for n <= 22.
# This avoids the AS89/Edgeworth approximation used by base R cor.test()
# for Spearman tests with n >= 10, which can return a numerical P value of 0
# for extremely strong correlations.
# ============================================================
run_exact_spearman <- function(x, y, alternative = "less") {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    n <- length(x)

    if (n < 3L)
        stop("At least 3 complete dose means are required for Spearman correlation.")
    if (n > 22L) {
        stop("Exact Spearman helper is configured for n <= 22; found n = ", n, ".")
    }

    pspearman::spearman.test(x = x, y = y, alternative = alternative, approximation = "exact")
}

# Always display P values in scientific notation.
# If a value is below the display floor (or numerically returned as 0),
# report it as 'P < ...' rather than 'P = 0'.
P_DISPLAY_FLOOR <- 1e-16

format_p_value <- function(p, digits = 2L, display_floor = P_DISPLAY_FLOOR) {
    if (!is.finite(p))
        return("P = NA")

    floor_text <- formatC(display_floor, format = "e", digits = digits, flag = "")

    if (p <= 0 || p < display_floor) {
        return(paste0("P < ", floor_text))
    }

    paste0("P = ", formatC(p, format = "e", digits = digits, flag = ""))
}

format_spearman_text <- function(test) {
    paste0("Spearman's rho = ", sprintf("%.3f", unname(test$estimate)), "\n", format_p_value(test$p.value))
}

# ============================================================
# User paths / output settings
# ============================================================
DIGEST_FILE <- Sys.getenv("INOSINE_INPUT_FILE", unset = "data/07/20260322data-ITP.xlsx")
DUAL_FILE <- Sys.getenv("DUAL_INPUT_FILE", "data/07/dual_plot.txt")
WB_FILE <- Sys.getenv("WB_INPUT_FILE", "data/07/WB.txt")

OUT_DIR <- Sys.getenv("LUC_PLOT_OUTPUT_DIR", "results/07")
SAVE_PLOTS <- TRUE

if (SAVE_PLOTS && !dir.exists(OUT_DIR)) {
    dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# Common theme for panels D/E
# ============================================================
my_theme <- theme_classic(base_size = 14) + theme(axis.line = element_line(color = "black",
    linewidth = 1), axis.text = element_text(color = "black", size = 12), axis.text.x = element_text(angle = 60,
    hjust = 1), axis.title = element_text(color = "black", size = 14, face = "bold"), plot.title = element_text(hjust = 0.5,
    face = "bold"))

# ============================================================
# PART 1. RNA digest inosine quantification
# ============================================================

sheet_name <- "Sheet1"
I_RESP_COL <- 25L  # Excel Y: I Resp.
I_FINAL_CONC_COL <- 26L  # Excel Z: I Final Conc.

FORCE_STANDARD_THROUGH_ZERO <- FALSE
CLIP_NEGATIVE_TO_ZERO <- TRUE
APPLY_OUTLIER_FILTER <- TRUE
OUTLIER_RATIO_CUTOFF <- 10
ZERO_TOL <- 1e-12
DROP_IF_ONLY_ONE_NONZERO_DEV <- FALSE

INOSINE_MW_G_MOL <- 268.23
REPORT_UNITS_PER_PG_UL <- 1000
ADDITIONAL_DILUTION_FACTOR <- 1

CONC_MAP <- tibble(rrl_conc_id = 1:14, target_f = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
    0.4, 0.45, 0.5, 0.55, 0.6, 0.65), concentration = c(0, 3.95, 8.33, 13.24, 18.75, 25,
    32.14, 40.38, 50, 61.36, 75, 91.67, 112.5, 139.29))

BASE_SIZE_DIGEST <- 12
SHOW_DIGEST_STATS <- TRUE
DIGEST_X_LABEL <- "ITP working-stock concentration (mM)"
DIGEST_Y_LABEL <- "Inosine concentration in RNA digests (mM)"

if (!file.exists(DIGEST_FILE)) stop("Digest input file not found: ", DIGEST_FILE)

for (value in c(INOSINE_MW_G_MOL, REPORT_UNITS_PER_PG_UL, ADDITIONAL_DILUTION_FACTOR, OUTLIER_RATIO_CUTOFF)) {
    if (!is.finite(value) || value <= 0) {
        stop("Molecular weight, conversion factors and cutoff must be positive.")
    }
}

if (anyDuplicated(CONC_MAP$rrl_conc_id) || any(!is.finite(CONC_MAP$concentration) | CONC_MAP$concentration <
    0)) {
    stop("CONC_MAP IDs must be unique and concentrations must be finite and non-negative.")
}

to_num <- function(x) {
    suppressWarnings(readr::parse_number(as.character(x)))
}

pg_uL_to_mM <- function(x) x/(INOSINE_MW_G_MOL * 1000)

mM_axis_labels <- function(x) {
    format(x, scientific = FALSE, trim = TRUE, digits = 5)
}

stopifnot(isTRUE(all.equal(pg_uL_to_mM(268.23), 0.001)))

flag_triplicate_outlier <- function(.x, .y) {
    .x <- .x %>%
        mutate(outlier_median = NA_real_, outlier_abs_dev = NA_real_, outlier_ratio = NA_real_,
            remove_as_outlier = FALSE)

    if (!isTRUE(.y$sample_type == "RRL"))
        return(.x)

    idx <- which(is.finite(.x$final_conc_recalc))
    if (length(idx) < 3)
        return(.x)

    med <- median(.x$final_conc_recalc[idx], na.rm = TRUE)
    abs_dev <- abs(.x$final_conc_recalc[idx] - med)

    .x$outlier_median[idx] <- med
    .x$outlier_abs_dev[idx] <- abs_dev

    nonzero_dev <- abs_dev[abs_dev > ZERO_TOL]

    if (length(nonzero_dev) >= 2) {
        min_nonzero <- min(nonzero_dev, na.rm = TRUE)
        max_dev <- max(nonzero_dev, na.rm = TRUE)
        ratio <- max_dev/min_nonzero
        .x$outlier_ratio[idx] <- ratio

        if (is.finite(ratio) && ratio >= OUTLIER_RATIO_CUTOFF) {
            max_local <- which(abs_dev == max(abs_dev, na.rm = TRUE))
            if (length(max_local) == 1) {
                .x$remove_as_outlier[idx[max_local]] <- TRUE
            }
        }
    }

    if (length(nonzero_dev) == 1 && DROP_IF_ONLY_ONE_NONZERO_DEV) {
        max_local <- which(abs_dev == max(abs_dev, na.rm = TRUE))
        .x$outlier_ratio[idx] <- Inf
        if (length(max_local) == 1 && max(abs_dev, na.rm = TRUE) > ZERO_TOL) {
            .x$remove_as_outlier[idx[max_local]] <- TRUE
        }
    }

    .x
}

raw_digest <- readxl::read_excel(DIGEST_FILE, sheet = sheet_name, col_names = FALSE, .name_repair = "minimal")

if (ncol(raw_digest) < I_FINAL_CONC_COL || nrow(raw_digest) < 3L) {
    stop("Digest table dimensions are not as expected; check sheet name and column indices.")
}

dat_digest <- tibble(source_excel_row = seq_len(nrow(raw_digest)), sample = as.character(raw_digest[[3]]),
    data_file = as.character(raw_digest[[4]]), resp = to_num(raw_digest[[I_RESP_COL]]), final_conc_reported = to_num(raw_digest[[I_FINAL_CONC_COL]]),
    accuracy_reported = if (ncol(raw_digest) >= I_FINAL_CONC_COL + 1L) {
        to_num(raw_digest[[I_FINAL_CONC_COL + 1L]])
    } else {
        rep(NA_real_, nrow(raw_digest))
    }) %>%
    slice(-(1:2)) %>%
    mutate(sample = str_trim(sample)) %>%
    filter(!is.na(sample), sample != "")

std_df_I <- dat_digest %>%
    filter(str_detect(sample, regex("^(P0|Std\\s*\\d+)$", ignore_case = TRUE))) %>%
    filter(is.finite(resp), is.finite(final_conc_reported))

if (nrow(std_df_I) < 3L || n_distinct(std_df_I$resp) < 2L) {
    stop("Insufficient valid P0/Std records for inosine concentration recalculation.")
}

cal_fit_I <- if (FORCE_STANDARD_THROUGH_ZERO) {
    lm(final_conc_reported ~ 0 + resp, data = std_df_I)
} else {
    lm(final_conc_reported ~ resp, data = std_df_I)
}

if (!is.finite(coef(cal_fit_I)[["resp"]]) || coef(cal_fit_I)[["resp"]] <= 0) {
    stop("Invalid inosine response-to-concentration slope; check the input data.")
}

calibration_summary_I <- tibble(n_standard_points = nrow(std_df_I), intercept = if (FORCE_STANDARD_THROUGH_ZERO) 0 else unname(coef(cal_fit_I)[1]),
    slope = unname(coef(cal_fit_I)[["resp"]]), r_squared = summary(cal_fit_I)$r.squared,
    status = "Legacy software-report mapping; calibration not revalidated")

rrl_df_I_raw <- dat_digest %>%
    filter(str_detect(sample, regex("^RRL\\d+", ignore_case = TRUE))) %>%
    mutate(base = "I", sample_type = "RRL", rrl_conc_id = as.integer(str_match(sample, regex("^RRL(\\d+)",
        ignore_case = TRUE))[, 2]), replicate = if_else(str_detect(sample, "_\\d+$"), as.integer(str_match(sample,
        "_(\\d+)$")[, 2]), 1L)) %>%
    left_join(CONC_MAP, by = "rrl_conc_id")

if (nrow(rrl_df_I_raw) == 0L) stop("No RRL samples found.")
if (anyNA(rrl_df_I_raw$concentration)) {
    stop("CONC_MAP is missing RRL IDs: ", paste(unique(rrl_df_I_raw$rrl_conc_id[is.na(rrl_df_I_raw$concentration)]),
        collapse = ", "))
}

rrl_df_I_raw$final_conc_recalc_raw <- as.numeric(predict(cal_fit_I, newdata = rrl_df_I_raw))

rrl_df_I_raw <- rrl_df_I_raw %>%
    mutate(final_conc_recalc = if (CLIP_NEGATIVE_TO_ZERO) {
        pmax(0, final_conc_recalc_raw)
    } else {
        final_conc_recalc_raw
    }, final_conc_pg_uL_assayed = final_conc_recalc/REPORT_UNITS_PER_PG_UL, final_conc_mM_assayed = pg_uL_to_mM(final_conc_pg_uL_assayed),
        final_conc_mM = final_conc_mM_assayed * ADDITIONAL_DILUTION_FACTOR)

rrl_df_I_checked <- rrl_df_I_raw %>%
    group_by(base, sample_type, rrl_conc_id) %>%
    group_modify(flag_triplicate_outlier) %>%
    ungroup() %>%
    mutate(excluded_from_plot = APPLY_OUTLIER_FILTER & remove_as_outlier)

removed_outliers_I <- rrl_df_I_checked %>%
    filter(excluded_from_plot)

rrl_df_I <- rrl_df_I_checked %>%
    filter(!excluded_from_plot, is.finite(concentration), is.finite(final_conc_mM))

if (nrow(rrl_df_I) == 0L) stop("No RRL measurements remain after filtering.")

message("[Digest] measurements used: ", nrow(rrl_df_I), "; removed by configured outlier rule: ",
    nrow(removed_outliers_I))

rrl_summary_I <- rrl_df_I %>%
    group_by(rrl_conc_id, target_f, concentration) %>%
    summarise(n = n(), mean_conc_mM = mean(final_conc_mM), sd_conc_mM = sd(final_conc_mM),
        se_conc_mM = sd_conc_mM/sqrt(n), .groups = "drop") %>%
    arrange(concentration)

can_fit_I <- nrow(rrl_summary_I) >= 3L && n_distinct(rrl_summary_I$concentration) >= 2L
can_cor_I <- can_fit_I && n_distinct(rrl_summary_I$mean_conc_mM) >= 2L

# Pearson is retained here intentionally because this panel asks whether the
# measured inosine concentration follows a linear concentration-response.
# The 'exact' option discussed for Spearman is not applicable to Pearson.
pearson_test_I <- if (can_cor_I) {
    cor.test(x = rrl_summary_I$concentration, y = rrl_summary_I$mean_conc_mM, method = "pearson",
        alternative = "two.sided")
} else {
    NULL
}

lm_fit_I <- if (can_fit_I) {
    lm(mean_conc_mM ~ concentration, data = rrl_summary_I)
} else {
    NULL
}

mean_stats_I <- tibble(analysis = "Dose means after configured outlier handling", n_doses = nrow(rrl_summary_I),
    pearson_r = if (is.null(pearson_test_I)) NA_real_ else unname(pearson_test_I$estimate),
    pearson_p = if (is.null(pearson_test_I)) NA_real_ else pearson_test_I$p.value, lm_r_squared = if (is.null(lm_fit_I)) NA_real_ else summary(lm_fit_I)$r.squared)

label_I <- NULL
if (SHOW_DIGEST_STATS && !is.null(pearson_test_I)) {
    label_I <- paste0("Pearson r = ", sprintf("%.3f", unname(pearson_test_I$estimate)), "\n",
        format_p_value(pearson_test_I$p.value))
}

x_breaks_I <- sort(unique(rrl_summary_I$concentration))
min_gap_I <- if (length(x_breaks_I) >= 2L) min(diff(x_breaks_I)) else 1

p_I_mM <- ggplot() + geom_point(data = rrl_df_I, aes(x = concentration, y = final_conc_mM),
    position = position_jitter(width = min_gap_I * 0.05, height = 0, seed = 1), alpha = 0.65,
    size = 2) + geom_errorbar(data = rrl_summary_I, aes(x = concentration, ymin = mean_conc_mM -
    sd_conc_mM, ymax = mean_conc_mM + sd_conc_mM), width = min_gap_I * 0.12, na.rm = TRUE) +
    geom_point(data = rrl_summary_I, aes(x = concentration, y = mean_conc_mM), size = 3)

if (can_fit_I) {
    p_I_mM <- p_I_mM + geom_smooth(data = rrl_summary_I, aes(x = concentration, y = mean_conc_mM),
        method = "lm", formula = y ~ x, se = TRUE, level = 0.95)
}

p_I_mM <- p_I_mM + scale_x_continuous(breaks = x_breaks_I, labels = scales::label_number(accuracy = 0.01,
    trim = TRUE), expand = expansion(mult = c(0.02, 0.05))) + scale_y_continuous(labels = mM_axis_labels,
    expand = expansion(mult = c(0.03, 0.12))) + theme_classic(base_size = BASE_SIZE_DIGEST) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(8, 12,
        8, 8)) + labs(title = NULL, x = DIGEST_X_LABEL, y = DIGEST_Y_LABEL)

if (!is.null(label_I)) {
    p_I_mM <- p_I_mM + annotate("text", x = Inf, y = Inf, label = label_I, hjust = 1.05,
        vjust = 1.15, size = 4.5)
}

# ============================================================
# PART 2. Panel D: dual-luciferase, 90 min
# ============================================================

if (!file.exists(DUAL_FILE)) stop("Dual-luciferase input file not found: ", DUAL_FILE)

data_D <- data.table::fread(DUAL_FILE)

# Three biological replicates, each represented by three technical wells.
well_D <- data_D %>%
    mutate(Concentration = as.numeric(gsub(" mM", "", Inosine_concentration))) %>%
    pivot_longer(starts_with("well"), names_to = "Well", values_to = "Value") %>%
    mutate(Well_number = as.integer(sub("well_", "", Well)), Replicate = case_when(Well_number %in%
        1:3 ~ "rep1", Well_number %in% 4:6 ~ "rep2", Well_number %in% 7:9 ~ "rep3", TRUE ~
        NA_character_), Included_in_analysis = !is.na(Replicate))

long_D <- well_D %>%
    filter(Included_in_analysis) %>%
    group_by(Concentration, Replicate) %>%
    summarise(Technical_SD = sd(Value), N_technical_wells = n(), Value = mean(Value), Source_wells = paste(Well,
        collapse = "; "), .groups = "drop")
stopifnot(all(long_D$N_technical_wells == 3L), all(is.finite(long_D$Value)), all(table(long_D$Concentration) ==
    3L))

processed_D <- long_D %>%
    group_by(Concentration) %>%
    summarise(mean = mean(Value), sd = sd(Value), .groups = "drop") %>%
    arrange(Concentration)

cor_D <- run_exact_spearman(x = processed_D$Concentration, y = processed_D$mean, alternative = "less")

label_D <- format_spearman_text(cor_D)

# Three independent-replicate means per dose; each averages three technical wells.
p_D <- ggplot(processed_D, aes(x = Concentration, y = mean)) + geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.12, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "brown") + scale_x_continuous(breaks = processed_D$Concentration,
    labels = processed_D$Concentration) + labs(x = "ITP input (mM)", y = "Renilla / firefly activity") +
    annotate("text", x = Inf, y = Inf, label = label_D, size = 5, hjust = 1.05, vjust = 1.15) +
    my_theme + theme(plot.margin = margin(10, 14, 10, 10))

# ============================================================
# PART 3. Panel E: WB quantification, normalize 0 mM to 1
# ============================================================

if (!file.exists(WB_FILE)) stop("WB input file not found: ", WB_FILE)

data_E <- data.table::fread(WB_FILE)

long_E <- data_E %>%
    mutate(Concentration = as.numeric(gsub(" mM", "", Inosine_concentration))) %>%
    pivot_longer(cols = starts_with("well"), names_to = "Well", values_to = "Value") %>%
    filter(is.finite(Value), Value > 0)

baseline_0 <- mean(long_E$Value[long_E$Concentration == 0], na.rm = TRUE)

if (!is.finite(baseline_0) || baseline_0 <= 0) {
    stop("Cannot normalize panel E: 0 mM WB mean is missing or non-positive.")
}

processed_E <- long_E %>%
    mutate(Value_norm = Value/baseline_0) %>%
    group_by(Concentration) %>%
    summarise(mean = mean(Value_norm, na.rm = TRUE), sd = sd(Value_norm, na.rm = TRUE), .groups = "drop") %>%
    arrange(Concentration)

cor_E <- run_exact_spearman(x = processed_E$Concentration, y = processed_E$mean, alternative = "less")

label_E <- format_spearman_text(cor_E)

# Dynamic upper limit so that error bars are not clipped; explicitly include 1
# among the y-axis tick labels to show the 0 mM normalization baseline.
y_max_E <- max(processed_E$mean + processed_E$sd, long_E$Value/baseline_0, na.rm = TRUE)
y_upper_E <- max(1.2, ceiling(y_max_E * 10)/10)
y_breaks_E <- sort(unique(c(pretty(c(0, y_upper_E)), 1)))
y_breaks_E <- y_breaks_E[y_breaks_E >= 0 & y_breaks_E <= y_upper_E]

p_E <- ggplot(processed_E, aes(x = Concentration, y = mean)) + geom_hline(yintercept = 1,
    linetype = "dashed", linewidth = 0.6, color = "grey60") + geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.12, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "brown") + scale_x_continuous(breaks = processed_E$Concentration,
    labels = processed_E$Concentration) + scale_y_continuous(breaks = y_breaks_E, limits = c(0,
    y_upper_E), expand = expansion(mult = c(0, 0.08))) + labs(x = "ITP input (mM)", y = "Relative Renilla / firefly\nprotein abundance") +
    annotate("text", x = Inf, y = Inf, label = label_E, size = 5, hjust = 1.05, vjust = 1.15) +
    my_theme + theme(plot.margin = margin(10, 14, 10, 10))

# ============================================================
# Statistical results for main panels C-E
# ============================================================

cat("\n=== Digest: Pearson correlation ===\n")
print(pearson_test_I)

cat("\n=== Panel D: one-sided exact Spearman ===\n")
print(cor_D)

cat("\n=== Panel E: one-sided exact Spearman ===\n")
print(cor_E)

# ============================================================
# PART 4. Panel F: S2 cells, 48 h after transfection
# NOTE:
#   - Spearman test uses the ORIGINAL ITP concentration.
#   - log2(concentration + 1) is used ONLY for x-axis spacing.
#   - Raw RLU is divided by 1e6 only for y-axis display.
#   - Error bars are mean +/- SD across the three biological replicates.
# ============================================================

data_F <- data.frame(Inosine_concentration = c("0 mM", "0.1 mM", "1 mM", "2.5 mM", "5 mM",
    "10 mM"), well_1 = c(1466821, 1562016, 980864, 858581, 849538, 443171), well_2 = c(1845707,
    1417781, 1168401, 1043189, 938257, 492345), well_3 = c(1623857, 1504468, 1501635, 1170317,
    752628, 407052))

long_F <- data_F %>%
    mutate(Concentration = as.numeric(gsub(" mM", "", Inosine_concentration))) %>%
    pivot_longer(starts_with("well"), names_to = "Well", values_to = "Value") %>%
    mutate(log2_conc = log2(Concentration + 1), Value_plot = Value/1e+06)

processed_F <- data_F %>%
    mutate(Concentration = as.numeric(gsub(" mM", "", Inosine_concentration))) %>%
    pivot_longer(cols = starts_with("well"), names_to = "Well", values_to = "Value") %>%
    group_by(Concentration) %>%
    summarise(mean_RLU = mean(Value, na.rm = TRUE), sd_RLU = sd(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(log2_conc = log2(Concentration + 1), mean = mean_RLU/1e+06, sd = sd_RLU/1e+06) %>%
    arrange(Concentration)

# Statistics use raw concentration and raw RLU means.
# Scaling RLU by 1e6 would not alter Spearman rho or P value.
cor_F <- run_exact_spearman(x = processed_F$Concentration, y = processed_F$mean_RLU, alternative = "less")

label_F <- format_spearman_text(cor_F)

p_F <- ggplot(processed_F, aes(x = log2_conc, y = mean)) + geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.08, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "brown") + scale_x_continuous(breaks = processed_F$log2_conc,
    labels = processed_F$Concentration) + labs(title = "S2 cells, 48 h after transfection",
    x = "ITP input (mM)", y = expression(paste("Renilla activity (", 10^6, " RLU)"))) + annotate("text",
    x = Inf, y = Inf, label = label_F, size = 5, hjust = 1.05, vjust = 1.15) + my_theme +
    theme(plot.margin = margin(10, 14, 10, 10))

cat("\n=== Panel F: one-sided exact Spearman ===\n")
print(cor_F)

# ============================================================
# PART 5. Supplementary Figure S15, panels D-G
# Source: Westernblot_and_luciferase.rds
#
# Panel order in the final supplementary figure:
#   D = luciferase activity, mRNA with 5' UTR       (lfs_with_utr)
#   E = WB quantification, mRNA with 5' UTR         (wb_with_utr)
#   F = luciferase activity, mRNA without 5' UTR    (lfs_no_utr)
#   G = WB quantification, mRNA without 5' UTR      (wb_no_utr)
#
# Statistics for all four panels:
#   - Spearman correlation uses ORIGINAL ITP concentration.
#   - alternative = 'less' tests the pre-specified negative association.
#   - exact = TRUE is used for the six dose means.
#   - log2(1 + concentration) is used for x-axis spacing only.
#   - The lm line is shown only as a visual guide; inference is from Spearman.
# ============================================================

SUPP_RDS_FILE <- Sys.getenv("LUC_WB_RDS_FILE", unset = "data/07/Westernblot_and_luciferase.rds")

if (!file.exists(SUPP_RDS_FILE)) {
    stop("Supplementary RDS file not found: ", SUPP_RDS_FILE)
}

wb_supp <- readRDS(SUPP_RDS_FILE)

# Original object order from new_Fig3E-FS13_lfs_and_WB.R
supp_table_names <- c("wb_no_utr", "wb_with_utr", "lfs_no_utr", "lfs_with_utr")
concentration_order_supp <- c("c_0mM", "c_0.1mM", "c_1mM", "c_2.5mM", "c_5mM", "c_10mM")

prepare_supp_table <- function(index, n_rep) {
    d <- data.frame(wb_supp[index])

    expected_ncol <- n_rep + 2L  # replicate columns + original mean + se
    if (ncol(d) != expected_ncol) {
        stop("Unexpected number of columns for ", supp_table_names[index], ": expected ",
            expected_ncol, ", observed ", ncol(d))
    }

    # Keep the replicate-level measurements and ignore the stored mean/SE for plotting.
    # Mean and SD are recalculated directly from the replicate values so that all panels
    # use mean +/- SD consistently.
    colnames(d) <- c(paste0("rep", seq_len(n_rep)), "mean_original", "se_original")

    missing_rows <- setdiff(concentration_order_supp, rownames(d))
    if (length(missing_rows) > 0L) {
        stop("Missing concentration rows in ", supp_table_names[index], ": ", paste(missing_rows,
            collapse = ", "))
    }

    d <- d[concentration_order_supp, , drop = FALSE]
    d$Concentration <- as.numeric(gsub("c_|mM", "", rownames(d)))

    rep_cols <- paste0("rep", seq_len(n_rep))
    rep_mat <- as.matrix(d[, rep_cols, drop = FALSE])
    storage.mode(rep_mat) <- "numeric"

    d$mean <- rowMeans(rep_mat, na.rm = TRUE)
    d$sd <- apply(rep_mat, 1, stats::sd, na.rm = TRUE)

    d
}

format_spearman_label <- function(test) format_spearman_text(test)

make_supp_panel <- function(data, y_label, panel_title = NULL) {
    test <- run_exact_spearman(x = data$Concentration, y = data$mean, alternative = "less")

    label <- format_spearman_label(test)
    y_top <- max(data$mean + data$sd, na.rm = TRUE)
    y_min <- min(data$mean - data$sd, na.rm = TRUE)
    y_span <- y_top - y_min
    if (!is.finite(y_span) || y_span <= 0)
        y_span <- max(abs(y_top), 1)

    points <- data %>%
        pivot_longer(matches("^rep[0-9]+$"), names_to = "Replicate", values_to = "Value")
    # Match the manuscript S15 caption; transform spacing only, not the test.
    data$plot_x <- log2(data$Concentration + 1)
    points$plot_x <- log2(points$Concentration + 1)
    # Symmetric horizontal offsets expose equal-valued repeats (e.g. both WB controls).
    # These display offsets never enter the dose summaries or statistical tests.
    rep_ids <- sort(unique(points$Replicate))
    offsets <- seq(-0.05, 0.05, length.out = length(rep_ids))
    points$point_x <- points$plot_x + offsets[match(points$Replicate, rep_ids)]
    p <- ggplot(data, aes(x = plot_x, y = mean)) + geom_point(data = points, aes(x = point_x,
        y = Value), inherit.aes = FALSE, size = 1.8, alpha = 0.8) + geom_point(size = 3,
        color = "blue") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15,
        color = "black") + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "brown") +
        scale_x_continuous(breaks = data$plot_x, labels = data$Concentration) + labs(title = panel_title,
        x = "ITP input (mM)", y = y_label) + annotate("text", x = Inf, y = Inf, label = label,
        size = 4.6, hjust = 1.05, vjust = 1.15) + my_theme + theme(plot.margin = margin(10,
        14, 10, 10))

    list(plot = p, test = test, data = data, points = points)
}

# The RDS contains two WB replicates and three luciferase-assay replicates.
wb_no_utr_supp <- prepare_supp_table(index = 1, n_rep = 2)
wb_with_utr_supp <- prepare_supp_table(index = 2, n_rep = 2)
lfs_no_utr_supp <- prepare_supp_table(index = 3, n_rep = 3)
lfs_with_utr_supp <- prepare_supp_table(index = 4, n_rep = 3)

# Use the final supplementary-figure panel order.
res_S15D <- make_supp_panel(lfs_with_utr_supp, y_label = "Relative luciferase activity",
    panel_title = "mRNA with 5' UTR")

res_S15E <- make_supp_panel(wb_with_utr_supp, y_label = "Normalized luciferase expression",
    panel_title = NULL)

res_S15F <- make_supp_panel(lfs_no_utr_supp, y_label = "Relative luciferase activity", panel_title = "mRNA without 5' UTR")

res_S15G <- make_supp_panel(wb_no_utr_supp, y_label = "Normalized luciferase expression",
    panel_title = NULL)

p_S15D <- res_S15D$plot
p_S15E <- res_S15E$plot
p_S15F <- res_S15F$plot
p_S15G <- res_S15G$plot

cat("\n=== Supplementary Fig. S15D: with 5' UTR, luciferase activity ===\n")
print(res_S15D$test)
cat("\n=== Supplementary Fig. S15E: with 5' UTR, WB quantification ===\n")
print(res_S15E$test)
cat("\n=== Supplementary Fig. S15F: without 5' UTR, luciferase activity ===\n")
print(res_S15F$test)
cat("\n=== Supplementary Fig. S15G: without 5' UTR, WB quantification ===\n")
print(res_S15G$test)

# ============================================================

# Export only the C/E/F panels here; D and S15 use the native-vector renderer.
for (panel in c("C", "E", "F")) {
    p <- switch(panel, C = p_I_mM, E = p_E, F = p_F)
    ggplot2::ggsave(file.path(OUT_DIR, paste0("Figure4", panel, ".pdf")), p, width = 7, height = 5)
}
