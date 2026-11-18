# Selected Figure S26 panel geometry and compact-letter display.
compact_letters_from_dunn <- function(groups, pairwise, alpha = 0.05) {
    groups <- as.character(groups)
    n_groups <- length(groups)
    nonsignificant <- diag(TRUE, n_groups)
    dimnames(nonsignificant) <- list(groups, groups)

    if (nrow(pairwise)) {
        for (i in seq_len(nrow(pairwise))) {
            g1 <- pairwise$group_1[[i]]
            g2 <- pairwise$group_2[[i]]
            if (g1 %chin% groups && g2 %chin% groups) {
                is_ns <- is.finite(pairwise$dunn_BH_FDR_within_plot[[i]]) && pairwise$dunn_BH_FDR_within_plot[[i]] >=
                  alpha
                nonsignificant[g1, g2] <- is_ns
                nonsignificant[g2, g1] <- is_ns
            }
        }
    }

    candidate_cliques <- lapply(seq_len(2^n_groups - 1L), function(mask) {
        members <- which(as.logical(intToBits(mask)[seq_len(n_groups)]))
        if (all(nonsignificant[members, members, drop = FALSE]))
            members else NULL
    })
    candidate_cliques <- Filter(Negate(is.null), candidate_cliques)

    maximal <- vapply(seq_along(candidate_cliques), function(i) {
        !any(vapply(seq_along(candidate_cliques), function(j) {
            i != j && length(candidate_cliques[[j]]) > length(candidate_cliques[[i]]) &&
                all(candidate_cliques[[i]] %in% candidate_cliques[[j]])
        }, logical(1L)))
    }, logical(1L))
    cliques <- candidate_cliques[maximal]
    clique_order <- order(vapply(cliques, min, integer(1L)), -vapply(cliques, length, integer(1L)))
    cliques <- cliques[clique_order]

    letter_pool <- c(letters, paste0("a", letters), paste0("b", letters))
    group_letters <- rep("", n_groups)
    for (i in seq_along(cliques)) {
        group_letters[cliques[[i]]] <- paste0(group_letters[cliques[[i]]], letter_pool[[i]])
    }

    data.table(plot_group = groups, significance_letter = group_letters)
}

make_panel <- function(sex_use, mutant_use) {
    d <- copy(plot_data[sex == sex_use & mutant == mutant_use])
    s <- copy(group_summary[sex == sex_use & mutant == mutant_use])
    pw <- copy(dunn[sex == sex_use & mutant == mutant_use])

    retained_groups <- group_levels[group_levels %chin% unique(as.character(d$plot_group))]
    d[, plot_group := factor(as.character(plot_group), levels = retained_groups)]
    d[, plotted_log2FoldChange := pmax(-0.75, pmin(0.75, log2FoldChange))]
    s[, plot_group := as.character(plot_group)]

    s[, original_group_order := match(plot_group, group_levels)]
    groups_high_to_low <- s[order(-median_log2FoldChange, original_group_order), plot_group]
    letters_one <- compact_letters_from_dunn(groups_high_to_low, pw, alpha = 0.05)
    s <- merge(s, letters_one, by = "plot_group", all.x = TRUE)
    s[, plot_group := factor(plot_group, levels = retained_groups)]
    s[, display_mean_log2FC := fifelse(abs(mean_log2FoldChange) < 0.005, 0, mean_log2FoldChange)]
    s[, top_label := sprintf("atop(italic(n) == %d, %.2f %%+-%% %.2f)", n_genes, display_mean_log2FC,
        sd_log2FoldChange)]

    # Cap displayed values only; tests and labels use the original log2FC.
    y_low_core <- -0.75
    y_high_core <- 0.75
    letter_y <- 0.82
    n_label_y <- 0.99

    s[, `:=`(letter_y = letter_y, n_label_y = n_label_y)]

    axis_labels <- group_axis_labels
    axis_labels[["unedit"]] <- "No CDS\nediting"
    scope_label <- "CDS sites"
    sex_label <- if (sex_use == "F")
        "Female" else "Male"
    adar_superscript <- if (mutant_use == "AdarG")
        "G" else "S"
    panel_title <- bquote(italic(Adar)^.(adar_superscript) ~ .(sex_label) ~ "|" ~ .(scope_label))

    ggplot(d, aes(x = plot_group, y = plotted_log2FoldChange, fill = plot_group)) + geom_hline(yintercept = 0,
        linetype = "dashed", linewidth = 0.28, colour = "#666666") + geom_violin(trim = TRUE,
        scale = "width", colour = NA, alpha = 0.78, width = 0.92) + geom_boxplot(width = 0.17,
        outlier.shape = NA, colour = "#333333", linewidth = 0.3, alpha = 0.62) + geom_text(data = s,
        aes(x = plot_group, y = n_label_y, label = top_label), inherit.aes = FALSE, parse = TRUE,
        size = 2.3, vjust = 0.5, lineheight = 0.92) + geom_text(data = s, aes(x = plot_group,
        y = letter_y, label = significance_letter), inherit.aes = FALSE, size = 3, fontface = "bold",
        vjust = 0.5) + scale_fill_manual(values = group_colours[retained_groups], drop = TRUE) +
        scale_x_discrete(labels = axis_labels[retained_groups], drop = TRUE) + scale_y_continuous(breaks = seq(-0.75,
        0.75, by = 0.25), expand = expansion(mult = c(0, 0))) + coord_cartesian(ylim = c(y_low_core,
        y_high_core), clip = "off") + labs(title = panel_title, x = NULL, y = expression(DESeq2 ~
        log[2] ~ fold ~ change ~ ("Adar"/"WBg"))) + theme_classic(base_size = 9) + theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 9.5, face = "plain", margin = margin(b = 27,
            unit = "pt")), axis.title.y = element_text(size = 8.5, margin = margin(r = 2,
            unit = "pt")), axis.text.x = element_text(size = 8.2, colour = "black", lineheight = 0.88),
        axis.text.y = element_text(size = 8, colour = "black"), axis.line = element_line(linewidth = 0.45,
            colour = "black"), axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks.length = unit(2.2, "pt"), plot.margin = margin(t = 3, r = 3, b = 3, l = 3,
            unit = "pt"))
}
