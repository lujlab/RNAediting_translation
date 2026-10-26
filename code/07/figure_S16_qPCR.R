#!/usr/bin/env Rscript
# Selected publication S16: 5/30/60 min, original unadjusted P values.
# Run via code/07/run_S16.sh; models and layout match the selected final PDF.
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript code/07/figure_S16_qPCR.R CT_TSV OUTPUT_DIRECTORY")
out <- normalizePath(args[2], mustWork = TRUE)
options(digits = 16)
raw <- read.delim(args[1], stringsAsFactors = FALSE)
raw$Group <- factor(raw$Group, levels = c("R0", "R10", "R30", "R50"))
raw$Time <- as.numeric(raw$Time)
stopifnot(nrow(raw) == 144L, all(is.finite(raw$Ct)), all(table(
    raw$Region, raw$Gene, raw$Group, raw$Time) == 3L))

# The released source measurements contain only the selected time points.
included_times <- c(5, 30, 60)
stopifnot(identical(sort(unique(raw$Time)), included_times))
dat <- raw
stopifnot(nrow(dat) == 144L, identical(sort(unique(dat$Time)), included_times))
normalize_ct <- function(d) {
    baseline <- aggregate(Ct ~ Region + Gene + Group, d[d$Time == 5, ], mean)
    names(baseline)[4] <- "baseline_Ct_5min"
    d <- merge(d, baseline, by = c("Region", "Gene", "Group"), sort = FALSE)
    d$DeltaCt <- d$Ct - d$baseline_Ct_5min
    d$ITP_mM <- unname(c(R0 = 0, R10 = 0.8, R30 = 3.2, R50 = 7.5)[as.character(d$Group)])
    d[order(d$Region, d$Gene, d$Group, d$Time, d$Technical_well), ]
}
dat <- normalize_ct(dat)

keys <- unique(dat[c("Region", "Gene", "Group", "ITP_mM", "Time")])
summary_rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    x <- dat[dat$Region == k$Region & dat$Gene == k$Gene &
             dat$Group == k$Group & dat$Time == k$Time, ]
    cbind(k, n_technical_wells = nrow(x), mean_Ct = mean(x$Ct),
          baseline_Ct_5min = x$baseline_Ct_5min[1],
          mean_DeltaCt = mean(x$DeltaCt), SD_DeltaCt = sd(x$DeltaCt),
          SEM_DeltaCt = sd(x$DeltaCt) / sqrt(nrow(x)))
})
summ <- do.call(rbind, summary_rows)
stopifnot(nrow(summ) == 48L, all(summ$n_technical_wells == 3L),
          max(abs(summ$mean_DeltaCt[summ$Time == 5])) < 1e-12)

# Preserve the original models: categorical Group; continuous numeric Time;
# Type-I sequential ANOVA; unadjusted P values; individual technical wells.
# The interaction is an omnibus slope-difference test, not a directional test
# of increasing RNA loss with increasing inosine incorporation.
fit_panel <- function(x, time_window) {
    main <- lm(DeltaCt ~ Group + Time, data = x)
    interaction <- lm(DeltaCt ~ Group * Time, data = x)
    a <- anova(main)
    b <- anova(interaction)
    data.frame(Region = x$Region[1], Gene = x$Gene[1], times_min = time_window,
               n_technical_wells = nrow(x), p_concentration = a["Group", "Pr(>F)"],
               p_time = a["Time", "Pr(>F)"],
               p_concentration_time = b["Group:Time", "Pr(>F)"],
               F_concentration = a["Group", "F value"],
               df_concentration = a["Group", "Df"],
               residual_df_main = df.residual(main),
               F_concentration_time = b["Group:Time", "F value"],
               df_concentration_time = b["Group:Time", "Df"],
               residual_df_interaction = df.residual(interaction))
}
panels <- expand.grid(Gene = c("firefly", "renilla"),
                      Region = c("3prime", "5prime"), stringsAsFactors = FALSE)
stats <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
    fit_panel(dat[dat$Region == panels$Region[i] & dat$Gene == panels$Gene[i], ], "5,30,60")
}))
stopifnot(all(stats$n_technical_wells == 36L), all(stats$residual_df_main == 31L),
          all(stats$residual_df_interaction == 28L))

fmt_p <- function(p) {
    ifelse(p < 0.0001, formatC(p, format = "e", digits = 2),
           formatC(p, format = "fg", digits = 3, flag = "#"))
}
stats$label <- paste0("Conc P = ", fmt_p(stats$p_concentration),
                     "\nConc \u00d7 time P = ", fmt_p(stats$p_concentration_time))
summ$Gene <- factor(summ$Gene, levels = c("firefly", "renilla"),
                    labels = c("Firefly", "Renilla"))
stats$Gene <- factor(stats$Gene, levels = c("firefly", "renilla"),
                     labels = c("Firefly", "Renilla"))
palette <- c(R0 = "#E41A1C", R10 = "#377EB8", R30 = "#4DAF4A", R50 = "#984EA3")
make_region_plot <- function(region, show_y_title = TRUE) {
    ggplot(summ[summ$Region == region, ], aes(Time, mean_DeltaCt, color = Group, group = Group)) +
        geom_hline(yintercept = 0, color = "grey65", linetype = "dashed", linewidth = 0.3) +
        geom_line(linewidth = 0.55) +
        geom_errorbar(aes(ymin = mean_DeltaCt - SEM_DeltaCt,
                          ymax = mean_DeltaCt + SEM_DeltaCt), width = 2.2, linewidth = 0.4) +
        geom_point(size = 1.5) +
        geom_text(data = stats[stats$Region == region, ],
                  aes(x = 5, y = 10.0, label = label), inherit.aes = FALSE,
                  hjust = 0, vjust = 1, size = 2.65, lineheight = 1.15, color = "black") +
        facet_wrap(~Gene, nrow = 1) +
        scale_color_manual(values = palette, labels = c("0", "0.8", "3.2", "7.5"),
                           name = "ITP used for Renilla transcription (mM)") +
        scale_x_continuous(breaks = included_times, limits = c(2, 63), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), expand = c(0, 0)) +
        coord_cartesian(ylim = c(-0.6, 10.5)) +
        labs(title = if (region == "3prime") "mRNA 3\u2032 end" else "mRNA 5\u2032 end",
             x = "Translation time (min)",
             y = if (show_y_title) expression(Delta*"Ct relative to 5 min") else NULL) +
        theme_classic(base_size = 9, base_family = "sans") +
        theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
              strip.background = element_blank(),
              strip.text = element_text(size = 10, face = "bold"),
              axis.text = element_text(color = "black", size = 8.5),
              axis.title = element_text(size = 9),
              axis.line = element_line(linewidth = 0.35),
              axis.ticks = element_line(linewidth = 0.35),
              panel.spacing = grid::unit(3.5, "mm"),
              legend.position = "bottom", legend.title = element_text(size = 8.5),
              legend.text = element_text(size = 8.5),
              legend.key.width = grid::unit(7, "mm"),
              legend.margin = margin(0, 0, 0, 0),
              plot.margin = margin(4, 4, 4, 4))
}
fig <- (make_region_plot("3prime") + make_region_plot("5prime", FALSE) +
        plot_layout(nrow = 1, guides = "collect")) & theme(legend.position = "bottom")
ggsave(file.path(out, "Figure_S16_no90min.pdf"), fig, device = cairo_pdf,
       width = 210, height = 88, units = "mm", bg = "white")
ggsave(file.path(out, "Figure_S16_no90min.png"), fig, device = "png", type = "cairo",
       width = 210, height = 88, units = "mm", dpi = 300, bg = "white")
cat("Saved Figure S16 PDF/PNG.\n")
