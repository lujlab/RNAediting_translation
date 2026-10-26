#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(grid)
})

# ============================================================
# 0. Command-line arguments
# ============================================================
#
# Usage:
#
# Rscript plot_male_Lp_vs_Lm_Lr_longestTx.R \
#   Dmel.longestTx_LmLr.tsv \
#   Dsim.longestTx_LmLr.tsv \
#   Dsec.longestTx_LmLr.tsv \
#   Dere.longestTx_LmLr.tsv \
#   (Dmel input may contain Dmel_ISO1/Dmel_WBg/AG/AS rows)
#   output_directory
#
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        paste0(
            "\nUsage:\n",
            "Rscript plot_male_Lp_vs_Lm_Lr_longestTx.R \\\n",
            "  <Dmel_Lp_wide_with_LmLr.tsv> \\\n",
            "  <Dsim_Lp_wide_with_LmLr.tsv> \\\n",
            "  <Dsec_Lp_wide_with_LmLr.tsv> \\\n",
            "  <Dere_Lp_wide_with_LmLr.tsv> \\\n",
            "  <output_dir>\n"
        )
    )
}

DMEL_FILE <- args[1]
DSIM_FILE <- args[2]
DSEC_FILE <- args[3]
DERE_FILE <- args[4]
OUT_DIR <- args[5]

INPUT_FILES <- c(
    Dmel = DMEL_FILE,
    Dsim = DSIM_FILE,
    Dsec = DSEC_FILE,
    Dere = DERE_FILE
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Plot settings
# ============================================================

# A4 portrait
A4_WIDTH_IN <- 8.27
A4_HEIGHT_IN <- 11.69

# Font.
# If Arial is unavailable on the server, change this to "sans".
BASE_FAMILY <- "Arial"

# Plot colors.
POINT_COLOR <- "#A9D3EA" # light blue
ADAR_COLOR <- "#9E2F2F" # dark red
DIAG_COLOR <- "#B33A3A" # y=x dashed line

# Panel order.
PANEL_ORDER <- c(
    "Dmel_ISO1",
    "Dmel_WBg",
    "Dmel_AG",
    "Dmel_AS",
    "Dsim",
    "Dsec",
    "Dere"
)
PANEL_NCOL <- 5

# Adar sites in the LONGEST-TRANSCRIPT coordinate system.
# These are matched directly against `new_handle`.
adar_sites_by_species <- c(
    Dmel = "FBtr0339272_1288",
    Dsim = "FBtr0358708_1399",
    Dsec = "FBtr0202151_1165",
    Dere = "FBtr0413873_1519"
)

# ============================================================
# 2. Utility functions
# ============================================================

check_file <- function(x) {
    if (!file.exists(x)) {
        stop("Input file does not exist: ", x)
    }
}

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(x))
}

canonical <- function(x) {
    tolower(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

# ------------------------------------------------------------
# Assign each row to the desired panels.
#
# For Dmel:
#   iso / iso1 -> Dmel_ISO1
#   WBg        -> Dmel_WBg
#   AG         -> Dmel_AG
#   AS         -> Dmel_AS
#
# For the other species, dataset identity determines the panel.
# ------------------------------------------------------------

assign_panel <- function(dataset_species, genotype) {
    if (dataset_species == "Dmel") {
        z <- canonical(genotype)

        if (z %in% c("iso", "iso1", "dmeliso1")) {
            return("Dmel_ISO1")
        }

        if (z %in% c("wbg", "dmelwbg")) {
            return("Dmel_WBg")
        }

        if (z %in% c("ag", "adarg", "dmelag")) {
            return("Dmel_AG")
        }

        if (z %in% c("as", "adars", "dmelas")) {
            return("Dmel_AS")
        }

        return(NA_character_)
    }

    if (dataset_species == "Dsim") return("Dsim")
    if (dataset_species == "Dsec") return("Dsec")
    if (dataset_species == "Dere") return("Dere")

    NA_character_
}

panel_species <- function(panel_id) {
    if (panel_id %in% c("Dmel_ISO1", "Dmel_WBg", "Dmel_AG", "Dmel_AS")) return("Dmel")
    if (panel_id == "Dsim") return("Dsim")
    if (panel_id == "Dsec") return("Dsec")
    if (panel_id == "Dere") return("Dere")
    NA_character_
}

# ------------------------------------------------------------
# new_handle can contain one or several semicolon-separated
# longest-transcript handles, e.g.
#
# FBtr0343368_3089;FBtr0343368_3088
#
# Adar is considered matched when its exact target handle occurs.
# ------------------------------------------------------------

split_handles <- function(x) {
    if (is.na(x) || !nzchar(trimws(as.character(x)))) {
        return(character())
    }

    z <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
    z <- trimws(z)
    z[nzchar(z)]
}

is_adar_site <- function(new_handle, species_key) {
    target <- adar_sites_by_species[[species_key]]

    if (is.null(target) || is.na(target) || !nzchar(target)) {
        return(FALSE)
    }

    target %in% split_handles(new_handle)
}

# ------------------------------------------------------------
# Paired Wilcoxon:
#
# Each point is one RNA-editing site.
# x and y for the same site are therefore paired observations.
#
# NA is removed.
# 0 and 1 are retained.
# ------------------------------------------------------------

paired_wilcox <- function(x, y) {
    keep <- is.finite(x) & is.finite(y)

    x <- x[keep]
    y <- y[keep]

    n <- length(x)

    if (n < 2) {
        return(list(n = n, p = NA_real_))
    }

    p <- tryCatch(
        suppressWarnings(
            wilcox.test(
                y,
                x,
                paired = TRUE,
                alternative = "two.sided",
                exact = FALSE
            )$p.value
        ),
        error = function(e) NA_real_
    )

    list(n = n, p = p)
}

# ------------------------------------------------------------
# Plotmath expression:
#
# italic(n) = ...
# italic(P) = 1.23 <U+00D7> 10^-4
# ------------------------------------------------------------

make_stat_expression <- function(n, p) {
    if (is.na(p)) {
        return(
            bquote(
                italic(n) == .(n) * "," ~ ~
                    italic(P) == "NA"
            )
        )
    }

    if (p == 0) {
        return(
            bquote(
                italic(n) == .(n) * "," ~ ~
                    italic(P) < 2.2 %*% 10^-308
            )
        )
    }

    sci <- formatC(p, format = "e", digits = 2)
    ss <- strsplit(sci, "e", fixed = TRUE)[[1]]

    mantissa <- as.numeric(ss[1])
    exponent <- as.integer(ss[2])

    bquote(
        italic(n) == .(n) * "," ~ ~
            italic(P) == .(mantissa) %*% 10^.(exponent)
    )
}

make_title_expression <- function(panel_id) {
    base <- switch(panel_id,
        Dmel_ISO1 =
            expression(italic(mel) ~ "(iso-1)"),
        Dmel_WBg =
            expression(italic(mel) ~ "(WBg)"),
        Dmel_AG =
            expression(italic(mel) ~ "(AG)"),
        Dmel_AS =
            expression(italic(mel) ~ "(AS)"),
        Dsim =
            expression(italic(sim)),
        Dsec =
            expression(italic(sec)),
        Dere =
            expression(italic(ere)),
        panel_id
    )

    base
}

# ============================================================
# 3. Read data
# ============================================================

for (ff in INPUT_FILES) {
    check_file(ff)
}

dt_list <- list()

for (sp in names(INPUT_FILES)) {
    f <- INPUT_FILES[[sp]]

    x <- fread(
        f,
        sep = "\t",
        header = TRUE,
        quote = "\"",
        na.strings = c("", "NA", "NaN", "nan")
    )

    x[, dataset_species := sp]
    x[, source_file := f]

    dt_list[[sp]] <- x
}

dt <- rbindlist(
    dt_list,
    use.names = TRUE,
    fill = TRUE
)

# ============================================================
# 4. Validate columns
# ============================================================

required_columns <- c(
    "site",
    "genotype",
    "new_handle",
    "Lp_mean_M",
    "Lp_mean_F",
    "mean_Lm_M",
    "mean_Lm_F",
    "mean_Lr_M",
    "mean_Lr_F"
)

missing_columns <- setdiff(required_columns, names(dt))

if (length(missing_columns) > 0) {
    stop(
        paste0(
            "Missing required columns:\n  ",
            paste(missing_columns, collapse = "\n  "),
            "\n\nDetected columns:\n  ",
            paste(names(dt), collapse = "\n  ")
        )
    )
}

# ============================================================
# 5. Prepare panels (shared by sex)
# ============================================================

dt[, panel_id := mapply(
    assign_panel,
    dataset_species,
    genotype
)]

dt <- dt[
    !is.na(panel_id) &
        panel_id %in% PANEL_ORDER
]

dt[, panel_id := factor(
    panel_id,
    levels = PANEL_ORDER
)]

# Explicit numeric conversion.
for (cc in c(
    "Lp_mean_M",
    "Lp_mean_F",
    "mean_Lm_M",
    "mean_Lm_F",
    "mean_Lr_M",
    "mean_Lr_F"
)) {
    set(
        dt,
        j = cc,
        value = safe_numeric(dt[[cc]])
    )
}

# Identify species for Adar lookup.
dt[, species_key := vapply(
    as.character(panel_id),
    panel_species,
    character(1)
)]

# Adar flag by LONGEST-TRANSCRIPT `new_handle`.
dt[, is_Adar := mapply(
    is_adar_site,
    new_handle,
    species_key
)]

# ------------------------------------------------------------
# Defensive de-duplication.
#
# Keep one row per panel <U+00D7> site.
# If duplicated rows exist, prefer the row with more complete
# male/female values.
# ------------------------------------------------------------

dt[, completeness := rowSums(
    !is.na(.SD)
),
.SDcols = c(
    "Lp_mean_M",
    "Lp_mean_F",
    "mean_Lm_M",
    "mean_Lm_F",
    "mean_Lr_M",
    "mean_Lr_F"
)
]

setorder(
    dt,
    panel_id,
    site,
    -completeness
)

dt <- unique(
    dt,
    by = c("panel_id", "site")
)

cat("\nRows retained before x/y complete-case filtering:\n")
print(
    dt[, .N, by = panel_id][
        match(PANEL_ORDER, panel_id)
    ]
)

cat("\nAdar matches in retained data:\n")
print(
    dt[
        is_Adar == TRUE,
        .(site, genotype, panel_id, new_handle)
    ]
)

# ============================================================
# 6. Draw one panel
# ============================================================

make_one_panel <- function(
  full_dt,
  panel_name,
  x_col,
  y_col,
  show_y_axis = FALSE
) {
    d <- full_dt[
        as.character(panel_id) == panel_name &
            is.finite(get(x_col)) &
            is.finite(get(y_col))
    ]

    xvals <- d[[x_col]]
    yvals <- d[[y_col]]

    st <- paired_wilcox(
        xvals,
        yvals
    )

    title_expr <- make_title_expression(panel_name)
    stat_expr <- make_stat_expression(st$n, st$p)

    # Ordinary points first.
    p <- ggplot(
        d,
        aes(
            x = .data[[x_col]],
            y = .data[[y_col]]
        )
    ) +
        geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed",
            linewidth = 0.50,
            color = DIAG_COLOR
        ) +
        geom_point(
            color = POINT_COLOR,
            size = 1.70,
            alpha = 0.90,
            shape = 16
        ) +
        coord_fixed(
            ratio = 1,
            xlim = c(-0.02, 1.02),
            ylim = c(-0.02, 1.02),
            expand = FALSE,
            clip = "off"
        ) +
        scale_x_continuous(
            breaks = c(
                0,
                0.25,
                0.50,
                0.75,
                1.00
            ),
            labels = sprintf(
                "%.2f",
                c(
                    0,
                    0.25,
                    0.50,
                    0.75,
                    1.00
                )
            )
        ) +
        scale_y_continuous(
            breaks = c(
                0,
                0.25,
                0.50,
                0.75,
                1.00
            ),
            labels = if (show_y_axis) {
                sprintf(
                    "%.2f",
                    c(
                        0,
                        0.25,
                        0.50,
                        0.75,
                        1.00
                    )
                )
            } else {
                rep("", 5)
            }
        ) +
        labs(
            title = title_expr,
            subtitle = stat_expr,
            x = NULL,
            y = NULL
        ) +
        theme_classic(
            base_size = 8.5,
            base_family = BASE_FAMILY
        ) +
        theme(
            plot.title = element_text(
                size = 10.0,
                hjust = 0.5,
                margin = margin(
                    b = 1.5
                )
            ),
            plot.subtitle = element_text(
                size = 7.8,
                hjust = 0.5,
                margin = margin(
                    b = 4.5
                )
            ),
            axis.text = element_text(
                size = 7.3,
                color = "black"
            ),
            axis.ticks = element_line(
                linewidth = 0.35,
                color = "black"
            ),
            axis.ticks.length = unit(
                1.3,
                "mm"
            ),
            axis.line = element_line(
                linewidth = 0.45,
                color = "black"
            ),
            plot.margin = margin(
                t = 4,
                r = 2,
                b = 5,
                l = 2,
                unit = "pt"
            ),
            legend.position = "none"
        )

    # ----------------------------------------------------------
    # Highlight Adar.
    # ----------------------------------------------------------

    d_adar <- d[
        is_Adar == TRUE
    ]

    if (nrow(d_adar) > 0) {
        # In principle there should be only one.
        d_adar <- d_adar[1]

        x_adar <- d_adar[[x_col]]
        y_adar <- d_adar[[y_col]]

        # Put label near the point and keep it within plot region.
        label_x <- min(
            x_adar + 0.035,
            0.82
        )

        label_y <- min(
            y_adar + 0.030,
            0.95
        )

        p <- p +

            geom_point(
                data = d_adar,
                aes(
                    x = .data[[x_col]],
                    y = .data[[y_col]]
                ),
                inherit.aes = FALSE,
                color = ADAR_COLOR,
                size = 2.8,
                shape = 16
            ) +

            annotate(
                "text",
                x = label_x,
                y = label_y,
                label = "Adar",
                family = BASE_FAMILY,
                fontface = "italic",
                color = ADAR_COLOR,
                size = 3.0,
                hjust = 0
            )
    }

    stats_dt <- data.table(
        panel_id = panel_name,
        x_variable = x_col,
        n = st$n,
        p_value = st$p,
        adar_found = nrow(d_adar) > 0
    )

    list(
        plot = p,
        stats = stats_dt
    )
}

# ============================================================
# 7. Save one A4 portrait PDF
# ============================================================

save_a4_figure <- function(
  combined_plot,
  out_pdf,
  x_label_expression
) {
    cairo_pdf(
        filename = out_pdf,
        width = A4_WIDTH_IN,
        height = A4_HEIGHT_IN,
        family = BASE_FAMILY
    )

    grid.newpage()

    # Panels.
    pushViewport(
        viewport(
            x = 0.5,
            y = 0.5,
            width = 0.965,
            height = 0.86,
            just = c(
                "center",
                "center"
            )
        )
    )

    grid.draw(
        patchworkGrob(
            combined_plot
        )
    )

    popViewport()

    # Shared x title.
    grid.text(
        x_label_expression,
        x = 0.5,
        y = 0.058,
        gp = gpar(
            fontsize = 13,
            fontfamily = BASE_FAMILY
        )
    )

    # Shared y title.
    grid.text(
        expression(L[p]),
        x = 0.022,
        y = 0.50,
        rot = 90,
        gp = gpar(
            fontsize = 13,
            fontfamily = BASE_FAMILY
        )
    )

    dev.off()
}

# ============================================================
# 8. Create one multi-panel figure
# ============================================================

make_figure <- function(
  x_col,
  y_col,
  x_label_expression,
  pdf_name,
  stats_name
) {
    results <- vector(
        "list",
        length(PANEL_ORDER)
    )

    for (i in seq_along(PANEL_ORDER)) {
        results[[i]] <- make_one_panel(
            full_dt = dt,
            panel_name = PANEL_ORDER[i],
            x_col = x_col,
            y_col = y_col,
            show_y_axis = (i == 1)
        )
    }

    plots <- lapply(
        results,
        `[[`,
        "plot"
    )

    stats <- rbindlist(
        lapply(
            results,
            `[[`,
            "stats"
        )
    )

    combined <- wrap_plots(
        plots,
        ncol = PANEL_NCOL
    )

    out_pdf <- file.path(
        OUT_DIR,
        pdf_name
    )

    out_stats <- file.path(
        OUT_DIR,
        stats_name
    )

    save_a4_figure(
        combined_plot = combined,
        out_pdf = out_pdf,
        x_label_expression = x_label_expression
    )

    fwrite(
        stats,
        out_stats,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )

    cat("\nSaved:\n")
    cat("  ", out_pdf, "\n", sep = "")
    cat("  ", out_stats, "\n", sep = "")
}

# ============================================================
# 9. Lm vs Lp: Male + Female
# ============================================================

make_figure(
    x_col = "mean_Lm_M",
    y_col = "Lp_mean_M",
    x_label_expression = expression(L[m]),
    pdf_name = "Male_Lm_vs_Lp_7panels_A4_portrait.pdf",
    stats_name = "Male_Lm_vs_Lp_7panels_stats.tsv"
)

make_figure(
    x_col = "mean_Lm_F",
    y_col = "Lp_mean_F",
    x_label_expression = expression(L[m]),
    pdf_name = "Female_Lm_vs_Lp_7panels_A4_portrait.pdf",
    stats_name = "Female_Lm_vs_Lp_7panels_stats.tsv"
)

# ============================================================
# 10. Lr vs Lp: Male + Female
# ============================================================

make_figure(
    x_col = "mean_Lr_M",
    y_col = "Lp_mean_M",
    x_label_expression = expression(L[r]),
    pdf_name = "Male_Lr_vs_Lp_7panels_A4_portrait.pdf",
    stats_name = "Male_Lr_vs_Lp_7panels_stats.tsv"
)

make_figure(
    x_col = "mean_Lr_F",
    y_col = "Lp_mean_F",
    x_label_expression = expression(L[r]),
    pdf_name = "Female_Lr_vs_Lp_7panels_A4_portrait.pdf",
    stats_name = "Female_Lr_vs_Lp_7panels_stats.tsv"
)

cat("\nAll figures finished.\n")
