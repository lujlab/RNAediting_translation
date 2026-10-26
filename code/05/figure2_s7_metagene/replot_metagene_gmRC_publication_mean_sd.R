#!/usr/bin/env Rscript

# Replot the completed metagene summary without rereading the large source
# mean-gmRC tables. The full analysis script contains the same plotting code.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
    normalizePath(dirname(sub("^--file=", "", script_arg[[1L]])),
        mustWork = TRUE
    )
} else {
    getwd()
}
outdir_name <- Sys.getenv("METAGENE_OUTDIR", unset = "metagene_gmRC")
outdir <- if (grepl("^/", outdir_name)) {
    outdir_name
} else {
    file.path(script_dir, outdir_name)
}
summary_file <- file.path(outdir, "metagene_gmRC_summary.tsv")
unit_file <- file.path(outdir, "metagene_units.tsv")
if (!file.exists(summary_file) || !file.exists(unit_file)) {
    stop("Completed metagene summary is missing in ", outdir)
}

meta_dt <- fread(summary_file, na.strings = c("", "NA", "NaN"))
unit_info <- fread(unit_file, na.strings = c("", "NA", "NaN"))
required_mean_sd <- c(
    "edited_sd", "control_sd",
    "edited_n_after_trim", "control_n_after_trim"
)
if (!all(required_mean_sd %in% names(meta_dt))) {
    stop(
        "Summary lacks mean/SD columns; recompute with ",
        "All_metagene_mean_existingMean_any1libLm02_mean_sd.R"
    )
}
desired_order <- c(
    "Dmel_Adar5G1__male",
    "Dmel_ISO1__female", "Dmel_ISO1__male",
    "Dmel_WBg__female", "Dmel_WBg__male",
    "Dmel_ISO1_WBg__female", "Dmel_ISO1_WBg__male",
    "Dsim__female", "Dsim__male",
    "Dsec__female", "Dsec__male",
    "Dere__female", "Dere__male"
)
unit_info[, unit_order := match(unit_id, desired_order)]
setorder(unit_info, unit_order)
meta_dt[, unit_order := match(unit_id, desired_order)]
setorder(meta_dt, unit_order, Position)

plot_data <- melt(
    meta_dt,
    id.vars = c(
        "unit_id", "species", "cohort", "sex", "label", "Position",
        "edited_sd", "control_sd",
        "edited_n_before_trim", "control_n_before_trim",
        "edited_n_after_trim", "control_n_after_trim", "n_edit_centers"
    ),
    measure.vars = c("edited", "control"),
    variable.name = "profile",
    value.name = "gmRC"
)
plot_data[, profile := factor(
    profile,
    levels = c("edited", "control"),
    labels = c("Edited sites", "Control sites")
)]

format_count <- function(x) {
    if (length(x) == 0L || !is.finite(x)) return("NA")
    formatC(as.numeric(x), format = "d", big.mark = ",")
}

format_gmrc <- function(x) {
    if (length(x) == 0L || !is.finite(x)) return("NA")
    formatC(as.numeric(x), format = "f", digits = 3L)
}

make_unit_title <- function(info) {
    species_name <- as.character(info$species[[1L]])
    species_plotmath <- if (species_name %in% c("Dsim", "Dsec", "Dere")) {
        sub("^D", "", species_name)
    } else {
        species_name
    }
    cohort_name <- as.character(info$cohort[[1L]])
    cohort_label <- if (cohort_name == "Dmel_ISO1") {
        "ISO1"
    } else if (cohort_name == "Dmel_WBg") {
        "WBg"
    } else if (cohort_name == "Dmel_ISO1_WBg") {
        "ISO1 + WBg"
    } else if (cohort_name == "Dmel_Adar5G1") {
        "Adar5G1"
    } else {
        ""
    }
    sex_label <- if (as.character(info$sex[[1L]]) == "female") "F" else "M"
    title_tail <- paste(
        c(cohort_label, sex_label)[nzchar(c(cohort_label, sex_label))],
        collapse = " "
    )
    parse(text = paste0(
        "italic(", species_plotmath, ")~'", title_tail, "'"
    ))[[1L]]
}

plot_objects <- lapply(seq_len(nrow(unit_info)), function(i) {
    info <- unit_info[i]
    current_unit_id <- info$unit_id[[1L]]
    d <- plot_data[unit_id == current_unit_id]
    center <- meta_dt[unit_id == current_unit_id & Position == 0][1L]
    center_label <- paste0(
        "Edited: ", format_gmrc(center$edited), " \u00b1 ",
        format_gmrc(center$edited_sd),
        " (n=", format_count(center$edited_n_before_trim), ")\n",
        "bg: ", format_gmrc(center$control), " \u00b1 ",
        format_gmrc(center$control_sd)
    )
    ggplot(d, aes(x = Position, y = gmRC, color = profile)) +
        geom_vline(
            xintercept = 0, linetype = 2, linewidth = 0.45,
            color = "grey45"
        ) +
        geom_line(linewidth = 0.45, na.rm = TRUE) +
        scale_color_manual(values = c(
            "Edited sites" = "#C0392B", "Control sites" = "#2C7FB8"
        ), drop = FALSE) +
        scale_x_continuous(
            breaks = c(-100, -50, 0, 50, 100),
            labels = c("-100", "-50", "0", "50", "100"),
            expand = expansion(mult = c(0.01, 0.01))
        ) +
        scale_y_continuous(
            breaks = seq(0.2, 2.2, by = 0.4),
            expand = expansion(mult = c(0, 0))
        ) +
        annotate(
            "text",
            x = 0, y = 2.5, label = center_label,
            size = 2.64, lineheight = 0.9, hjust = 0.5, vjust = 0.5
        ) +
        coord_cartesian(ylim = c(0.2, 2.2), clip = "off") +
        labs(
            title = make_unit_title(info),
            x = "Codon position relative to edit",
            y = "gmRC",
            color = NULL
        ) +
        theme_classic(base_size = 7.5) +
        theme(
            plot.title = element_text(
                size = 8, face = "bold", hjust = 0.5,
                margin = margin(t = 0, r = 0, b = 20, l = 0)
            ),
            axis.title = element_text(size = 7.5),
            axis.text = element_text(size = 7.5, color = "black"),
            axis.text.x = element_text(size = 7.5),
            axis.ticks = element_line(linewidth = 0.3, color = "black"),
            axis.line = element_line(linewidth = 0.35, color = "black"),
            legend.position = "bottom",
            legend.text = element_text(size = 7.5),
            legend.key.width = grid::unit(0.65, "cm"),
            legend.key.height = grid::unit(0.3, "cm"),
            plot.margin = margin(6, 3, 4, 3)
        )
})

names(plot_objects) <- unit_info$unit_id
metagene_pdf <- file.path(outdir, "metagene_gmRC_13_units.pdf")
pdf(
    metagene_pdf,
    width = 8.27, height = 11.69, onefile = TRUE, useDingbats = FALSE
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(
    y = grid::unit(0.5, "npc"),
    height = grid::unit(0.94, "npc"),
    layout = grid::grid.layout(
        6L, 3L,
        heights = grid::unit(c(0.18, 0.92, 0.92, 0.92, 0.92, 0.92), "null")
    ),
    gp = grid::gpar(fontfamily = "sans")
))
for (i in seq_along(plot_objects)) {
    row <- (i - 1L) %/% 3L + 2L
    col <- (i - 1L) %% 3L + 1L
    print(
        plot_objects[[i]],
        vp = grid::viewport(layout.pos.row = row, layout.pos.col = col)
    )
}
grid::popViewport()
dev.off()
file.copy(
    metagene_pdf,
    file.path(outdir, "metagene_gmRC_11_units.pdf"),
    overwrite = TRUE
)
message(
    "Replotted publication-style metagene PDF: ",
    metagene_pdf
)
