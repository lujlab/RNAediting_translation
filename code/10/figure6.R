#!/usr/bin/env Rscript
# Reproduce Figure 6 using sex-specific Kruskal-Wallis and Dunn tests.
suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
})
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(script))
input <- if (length(args)) args[1] else file.path(here, "phenotype_data.tsv")
out <- if (length(args) > 1) args[2] else here
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.delim(input, stringsAsFactors = FALSE)
stopifnot(all(is.finite(d$value)), all(d$Genotype %in% c("AdarG", "WT", "AdarS")))
groups <- c("AdarG", "WT", "AdarS")
colours <- c(AdarG = "#00BFC4", WT = "#FDAE61", AdarS = "#E777AE")

# Dunn's two-sided normal approximation with pooled ranks and tie correction.
dunn_bh <- function(x) {
    n <- nrow(x)
    ranks <- rank(x$value, ties.method = "average")
    sizes <- table(factor(x$Genotype, levels = groups))
    means <- tapply(ranks, factor(x$Genotype, levels = groups), mean)
    ties <- as.numeric(table(x$value))
    variance <- n * (n + 1)/12 - sum(ties^3 - ties)/(12 * (n - 1))
    pairs <- combn(groups, 2)
    ans <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
        a <- pairs[1, i]
        b <- pairs[2, i]
        z <- if (variance == 0)
            0 else (means[a] - means[b])/sqrt(variance * (1/sizes[a] + 1/sizes[b]))
        data.frame(group_1 = a, group_2 = b, z = unname(z), P = 2 * pnorm(-abs(z)))
    }))
    ans$P_BH <- p.adjust(ans$P, method = "BH")
    ans
}

# Maximal nonsignificant cliques provide an exact compact-letter display
# for three groups. Letters are ordered by descending mean phenotype value.
compact_letters <- function(pw, ordered_groups) {
    ns <- matrix(TRUE, 3, 3, dimnames = list(ordered_groups, ordered_groups))
    for (i in seq_len(nrow(pw))) ns[pw$group_1[i], pw$group_2[i]] <- ns[pw$group_2[i], pw$group_1[i]] <- pw$P_BH[i] >=
        0.05
    sets <- lapply(1:7, function(mask) which(as.logical(intToBits(mask)[1:3])))
    sets <- Filter(function(s) all(ns[s, s, drop = FALSE]), sets)
    maximal <- Filter(function(s) !any(vapply(sets, function(t) length(t) > length(s) &&
        all(s %in% t), logical(1))), sets)
    maximal <- maximal[order(vapply(maximal, min, integer(1)), -lengths(maximal))]
    lab <- setNames(rep("", 3), ordered_groups)
    for (i in seq_along(maximal)) lab[maximal[[i]]] <- paste0(lab[maximal[[i]]], letters[i])
    for (i in seq_len(nrow(pw))) {
        shared <- length(intersect(strsplit(lab[pw$group_1[i]], "")[[1]], strsplit(lab[pw$group_2[i]],
            "")[[1]])) > 0
        stopifnot(shared == (pw$P_BH[i] >= 0.05))
    }
    lab
}
kw_all <- pw_all <- summary_all <- list()
for (panel in LETTERS[1:8]) for (sex in unique(d$sex[d$panel == panel])) {
    x <- d[d$panel == panel & d$sex == sex, ]
    stopifnot(all(table(factor(x$Genotype, levels = groups)) > 1))
    kw <- kruskal.test(value ~ Genotype, data = x)
    pw <- dunn_bh(x)
    s <- do.call(rbind, lapply(groups, function(g) {
        v <- x$value[x$Genotype == g]
        data.frame(Genotype = g, n = length(v), mean = mean(v), SD = sd(v), maximum = max(v))
    }))
    lab <- compact_letters(pw, s$Genotype[order(-s$mean, match(s$Genotype, groups))])
    s$letter <- unname(lab[s$Genotype])
    s$panel <- panel
    s$sex <- sex
    pw$panel <- panel
    pw$sex <- sex
    key <- paste(panel, sex)
    kw_all[[key]] <- data.frame(panel = panel, sex = sex, n = nrow(x), statistic = unname(kw$statistic),
        df = unname(kw$parameter), P = kw$p.value)
    pw_all[[key]] <- pw
    summary_all[[key]] <- s
}
kw_all <- do.call(rbind, kw_all)
pw_all <- do.call(rbind, pw_all)
summary_all <- do.call(rbind, summary_all)
titles <- c(A = "Sector-boundary crossings / 2 min", B = "Climbing success (%)", C = "Survival after cold exposure (%)",
    D = "Recovery time after hypoxia (min)", E = "Adult offspring number", F = "Embryo hatching (%)",
    G = "Larval pupation (%)", H = "Pupal eclosion (%)")
tops <- c(A_Female = 80, A_Male = 75, B_Female = 80, B_Male = 90, C_Female = 112, C_Male = 112,
    D_Female = 22, D_Male = 26, E_Female = 250, F_Pooled = 115, G_Pooled = 115, H_Pooled = 115)
make_axis <- function(panel, sex) {
    x <- d[d$panel == panel & d$sex == sex, ]
    s <- summary_all[summary_all$panel == panel & summary_all$sex == sex, ]
    x$Genotype <- factor(x$Genotype, levels = groups)
    s$Genotype <- factor(s$Genotype, levels = groups)
    top <- unname(tops[paste(panel, sex, sep = "_")])
    s$label_y <- pmax(s$maximum, s$mean + s$SD) + 0.045 * top
    top <- max(top, max(s$label_y) * 1.055)
    genotype_text <- c("italic(Adar)^G", "WT", "italic(Adar)^S")
    xlabs <- parse(text = paste0("atop(", genotype_text, ",italic(n)==", s$n, ")"))
    ggplot(s, aes(Genotype, mean, fill = Genotype)) + geom_col(width = 0.63) + geom_errorbar(aes(ymin = mean -
        SD, ymax = mean + SD), width = 0.13, linewidth = 0.28) + geom_point(data = x, aes(y = value),
        position = position_jitter(width = 0.11, height = 0, seed = 1011), size = 0.75, alpha = 0.7,
        colour = "black") + geom_text(aes(y = label_y, label = letter), colour = "#3535CA",
        size = 2.6) + scale_fill_manual(values = colours) + scale_x_discrete(labels = xlabs,
        drop = FALSE) + scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = pretty(c(0,
        top), n = 4)) + coord_cartesian(ylim = c(0, top), clip = "off") + labs(x = NULL,
        y = NULL, title = if (panel %in% LETTERS[1:4])
            sex else NULL) + theme_classic(base_size = 8, base_family = "sans") + theme(legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5), axis.text = element_text(colour = "black",
            size = 7.5), axis.text.x = element_text(size = 8), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), axis.ticks.length = grid::unit(1, "mm"),
        plot.margin = margin(3, 6, 2, 5))
}
plots <- list()
for (panel in LETTERS[1:8]) {
    sexes <- if (panel %in% LETTERS[1:4])
        c("Female", "Male") else if (panel == "E")
        "Female" else "Pooled"
    axes <- lapply(sexes, function(s) make_axis(panel, s))
    body <- plot_grid(plotlist = axes, nrow = 1, align = "h", axis = "tb")
    plots[[panel]] <- ggdraw() + draw_plot(body, y = 0, height = 0.92) + draw_label(titles[panel],
        x = 0.55, y = 0.97, size = 8.5) + draw_label(panel, x = 0.012, y = 0.97, hjust = 0,
        fontface = "bold", size = 11)
}
for (panel in names(plots)) {
    width <- if (panel %in% LETTERS[1:4]) 3.6 else 1.8
    stem <- paste0("Figure6", panel, "_sexwise_KW_DunnBH")
    ggsave(file.path(out, paste0(stem, ".pdf")), plots[[panel]], width = width,
        height = 2.6, device = cairo_pdf)
    ggsave(file.path(out, paste0(stem, ".png")), plots[[panel]], width = width,
        height = 2.6, dpi = 240, bg = "white")
}
invisible(NULL)
message("Figure 6 panels A-H written to ", out)
