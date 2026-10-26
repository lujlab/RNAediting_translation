# Plotting and statistical helpers for Table1 filter scatter panels.
format_stat <- function(x, digits = 3L) {
    x <- as.numeric(x)
    out <- rep("NA", length(x))
    keep <- is.finite(x)
    out[keep] <- formatC(x[keep], format = "f", digits = digits)
    out
}

format_p_plotmath <- function(p) {
    if (!is.finite(p)) return("italic(P)==NA")
    if (p < 0.001) return("italic(P)<0.001")
    paste0("italic(P)==", formatC(p, format = "f", digits = 3L))
}

adar_sites_by_species <- c(
    Dmel = "FBtr0339272_1288",
    Dsim = "FBtr0358708_1399",
    Dsec = "FBtr0202151_1165",
    Dere = "FBtr0413873_1519"
)
make_paired_scatter <- function(
  d, x_col, y_col, x_axis_label, y_axis_label,
  metric_name, fixed_limits = NULL
) {
    d <- copy(d[is.finite(get(x_col)) & is.finite(get(y_col))])
    x <- d[[x_col]]
    y <- d[[y_col]]
    n_points <- nrow(d)
    n_edit_sites <- sum(d$n_total_edit_sites)
    delta_median <- if (n_points > 0L) median(y - x) else NA_real_
    paired_test <- if (n_points >= 2L) {
        tryCatch(
            wilcox.test(y, x, paired = TRUE, exact = FALSE),
            error = function(e) NULL
        )
    } else {
        NULL
    }
    p_value <- if (is.null(paired_test)) NA_real_ else paired_test$p.value

    upper <- if (is.null(fixed_limits)) {
        observed_max <- max(c(x, y), na.rm = TRUE)
        if (!is.finite(observed_max) || observed_max <= 0) 1 else observed_max * 1.06
    } else {
        fixed_limits[[2L]]
    }
    lower <- if (is.null(fixed_limits)) 0 else fixed_limits[[1L]]
    axis_span <- upper - lower

    species_name <- unique(as.character(d$species))
    cohort_name <- unique(as.character(d$cohort))
    panel_name <- unique(as.character(d$panel_id))
    sex_name <- unique(as.character(d$sex))
    if (length(species_name) != 1L || length(cohort_name) != 1L ||
        length(panel_name) != 1L || length(sex_name) != 1L) {
        stop("Each scatter panel must contain one cohort-sex unit")
    }
    title_text <- paste0(
        ifelse(panel_name == "Dmel_ISO1", "Dmel (ISO-1)",
            ifelse(panel_name == "Dmel_WBg", "Dmel (WBg)", panel_name)
        ), " ", sex_name
    )

    adar_site <- adar_sites_by_species[[species_name]]
    adar_pattern <- paste0("(^|,)", adar_site, "(,|$)")
    d[, is_adar := grepl(adar_pattern, cdspos_ids)]
    d[, `:=`(
        in_plot_window =
            get(x_col) >= lower & get(x_col) <= upper &
                get(y_col) >= lower & get(y_col) <= upper,
        # Preserve raw values for the statistical test, but display ratios beyond
        # the fixed 0-1 window at the boundary (e.g. Ar/Am = 1.4 is drawn at 1).
        x_plot = pmin(upper, pmax(lower, get(x_col))),
        y_plot = pmin(upper, pmax(lower, get(y_col))),
        x_highlight = pmin(upper, pmax(lower, get(x_col))),
        y_highlight = pmin(upper, pmax(lower, get(y_col)))
    )]
    adar_dt <- d[is_adar == TRUE]
    adar_dt[, `:=`(
        highlight_shape = fifelse(in_plot_window, 16L, 17L),
        highlight_label = "ADAR",
        label_hjust = fifelse(x_highlight >= upper - 0.08 * axis_span, 1.05, -0.12),
        label_vjust = fifelse(y_highlight >= upper - 0.08 * axis_span, 1.35, 1.15)
    )]

    n_label <- paste0("italic(n)==", n_points)
    delta_label <- paste0(
        "Delta*median==", format_stat(delta_median, digits = 4L)
    )
    p_label <- format_p_plotmath(p_value)
    # Keep the summary block away from the highlighted ADAR site.  Several
    # species have the ADAR point in the upper-left part of the ratio panel.
    stat_x <- upper - 0.04 * axis_span

    p <- ggplot(d, aes(x = x_plot, y = y_plot)) +
        geom_abline(
            intercept = 0, slope = 1,
            linetype = 2, linewidth = 1.05, color = "#8B1A1A"
        ) +
        geom_point(
            data = d[is_adar == FALSE],
            shape = 16, size = 1.8, alpha = 0.78, color = "#8FCBE6"
        ) +
        geom_point(
            data = adar_dt,
            aes(x = x_highlight, y = y_highlight, shape = highlight_shape),
            inherit.aes = FALSE,
            size = 3.4, alpha = 1, color = "#8B0000"
        ) +
        geom_text(
            data = adar_dt,
            aes(
                x = x_highlight, y = y_highlight, label = highlight_label,
                hjust = label_hjust, vjust = label_vjust
            ),
            inherit.aes = FALSE,
            color = "#8B0000", fontface = "bold",
            size = 2.8
        ) +
        annotate(
            "text",
            x = stat_x, y = lower + 0.28 * axis_span,
            label = n_label, parse = TRUE, hjust = 1, vjust = 1, size = 2.8
        ) +
        annotate(
            "text",
            x = stat_x, y = lower + 0.19 * axis_span,
            label = delta_label, parse = TRUE, hjust = 1, vjust = 1, size = 2.8
        ) +
        annotate(
            "text",
            x = stat_x, y = lower + 0.10 * axis_span,
            label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 2.8
        ) +
        scale_shape_identity() +
        scale_x_continuous(
            breaks = seq(lower, upper, length.out = 5L),
            expand = expansion(mult = 0)
        ) +
        scale_y_continuous(
            breaks = seq(lower, upper, length.out = 5L),
            expand = expansion(mult = 0)
        ) +
        coord_equal(xlim = c(lower, upper), ylim = c(lower, upper), clip = "on") +
        labs(
            title = title_text,
            x = x_axis_label,
            y = y_axis_label
        ) +
        theme_bw(base_size = 8) +
        theme(
            axis.line = element_line(color = "black", linewidth = 0.6),
            axis.title = element_text(size = 8, face = "plain", color = "black"),
            axis.text = element_text(size = 8, color = "black"),
            axis.ticks = element_line(linewidth = 0.6),
            axis.ticks.length = grid::unit(0.08, "cm"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(size = 8, hjust = 0.5),
            plot.margin = margin(3, 8, 3, 3),
            legend.position = "none"
        )

    list(
        plot = p,
        summary = data.table(
            species = species_name,
            cohort = cohort_name,
            panel_id = panel_name,
            sex = sex_name,
            metric = metric_name,
            n_points = n_points,
            n_edit_sites = n_edit_sites,
            n_visible_in_axis_window = sum(d$in_plot_window),
            n_outside_axis_window = sum(!d$in_plot_window),
            n_x_capped_at_upper = sum(d[[x_col]] > upper),
            n_y_capped_at_upper = sum(d[[y_col]] > upper),
            median_x = median(x),
            median_y = median(y),
            delta_median_y_minus_x = delta_median,
            wilcoxon_V = if (is.null(paired_test)) NA_real_ else unname(paired_test$statistic),
            p_value = p_value,
            p_adjustment = "none",
            adar_site = adar_site,
            adar_found = nrow(adar_dt) > 0L,
            adar_outside_axis_window = any(!adar_dt$in_plot_window),
            x_min = lower,
            x_max = upper,
            y_min = lower,
            y_max = upper
        )
    )
}
