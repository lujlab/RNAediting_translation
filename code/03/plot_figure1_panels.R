#!/usr/bin/env Rscript
# Draw Figure 1 panels from the input tables.
# Usage: plot_figure1_panels.R DATA_DIR FIGURE_DIR/panels [reclassify=1]
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2L) stop("Usage: plot_figure1_panels.R DATA_DIR FIGURE_DIR/panels [1|0]")
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
helper_dir <- dirname(normalizePath(script_path))
reclassify <- if (length(a) >= 3L) a[3] == "1" else TRUE
set.seed(20260913)
# The original UpSet code draws on the active device before capturing its grob.
scratch_device <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(scratch_device, width = 9, height = 5.6, family = "Arial")
options(warn = 1)
source(file.path(helper_dir, if (reclassify) "plot_figure1_from_raw_reclassified.R" else "plot_figure1_from_raw.R"))
figure_root <- dirname(output_dir)
# Export the exact matrix used by the original heatmap, including display zeros.
fwrite(as.data.table(heat_wide), file.path(table_dir, "heatmap.mean_matrix.display_zero.tsv"), sep = "\t")
fwrite(replicate_counts, file.path(table_dir, "replicate_counts.tsv"), sep = "\t")
# Record the missing cells independently, before the old display convention fills 0.
observed_keys <- unique(heat_mean[, .(sample_id, handle)])
missing_cells <- fsetdiff(CJ(sample_id = sample_order, handle = shared_handles), observed_keys)
fwrite(missing_cells, file.path(table_dir, "heatmap.missing_cells_filled_zero.tsv"), sep = "\t")

# N/S/Others stacking uses the ISO1 representative annotation, as in the old
# Python display code. No site memberships or annotation tables are changed.
iso_type <- fread(file.path(data_dir, "ISO1.annotation.tsv"))[, .(handle, type)]
iso_type[, display_class := fifelse(grepl("missense_variant", type), "N", fifelse(grepl("synonymous_variant", type), "S", "Others"))]
stacked <- merge(membership, iso_type[, .(handle, display_class)], by = "handle", all.x = TRUE)
stacked[is.na(display_class), display_class := "Others"]
intersections <- stacked[, .(total = .N, N = sum(display_class == "N"), S = sum(display_class == "S"), Others = sum(display_class == "Others")), by = .(pattern, Dmel, Dsim, Dsec, Dere)]
setorder(intersections, -total, pattern)
fwrite(intersections, file.path(table_dir, "upset.intersections.with_NS.tsv"), sep = "\t")

line <- function(x0, y0, x1, y1, col = "black", lwd = .8) grid.lines(c(x0, x1), c(y0, y1), gp = gpar(col = col, lwd = lwd))
label <- function(txt, x, y, size = 9, just = "centre", rot = 0, col = "black", fontface = 1) grid.text(txt, x, y, just = just, rot = rot, gp = gpar(fontfamily = font_family, fontsize = size, col = col, fontface = fontface))
box <- function(x, y, w, h, fill) grid.rect(x, y, w, h, gp = gpar(fill = fill, col = NA))

# Fixed phylogeny coordinates are retained from draw_phylogeny() in the source
# old plot_figure1_python.py. They are a schematic, not estimated divergence ages.
draw_A <- function() {
    k <- nrow(intersections)
    x0 <- .43
    x1 <- .985
    yb <- .40
    yt <- .91
    xx <- x0 + ((seq_len(k) - .5) / k) * (x1 - x0)
    bw <- (x1 - x0) / k * .60
    ymax <- ceiling(max(intersections$total) * 1.13 / 250) * 250
    yy <- function(v) yb + v / ymax * (yt - yb)
    for (i in seq_len(k)) {
        base <- 0
        for (cl in c("N", "S", "Others")) {
            n <- intersections[[cl]][i]
            box(xx[i], (yy(base) + yy(base + n)) / 2, bw, yy(base + n) - yy(base), c(N = "#E97274", S = "#58B6D9", Others = "#3B3B3B")[[cl]])
            base <- base + n
        }
        label(format_number(base), xx[i], yy(base) + .017, 9)
    }
    line(x0, yb, x1, yb)
    line(x0, yb, x0, yt)
    for (v in seq(0, ymax, 250)) {
        line(x0 - .005, yy(v), x0, yy(v))
        label(v, x0 - .009, yy(v), 8, just = "right")
    }
    label("Editing Sites Intersections", x0 - .055, (yb + yt) / 2, 10, rot = 90)
    for (j in 1:3) {
        lx <- .725 + (j - 1) * .078
        box(lx, .958, .018, .016, c("#58B6D9", "#E97274", "#3B3B3B")[j])
        label(c("S", "N", "Others")[j], lx + .013, .958, 8, just = "left")
    }
    sy <- c(.305, .235, .165, .095)
    for (j in 1:4) {
        if (j %% 2 == 0) box((x0 + x1) / 2, sy[j], x1 - x0, .068, "#F1F1F1")
        grid.points(xx, rep(sy[j], k), pch = 16, size = unit(1.8, "mm"), gp = gpar(col = "#E5E5E5"))
    }
    flags <- as.matrix(intersections[, .(Dmel, Dsim, Dsec, Dere)])
    for (i in seq_len(k)) {
        sel <- which(flags[i, ] == 1L)
        if (length(sel) > 1L) line(xx[i], min(sy[sel]), xx[i], max(sy[sel]), "#3B3B3B", 1.0)
        grid.points(rep(xx[i], length(sel)), sy[sel], pch = 16, size = unit(1.9, "mm"), gp = gpar(col = "#3B3B3B"))
    }
    sizes <- c(length(upset_handles$Dmel), length(upset_handles$Dsim), length(upset_handles$Dsec), length(upset_handles$Dere))
    end <- .208
    span <- .12
    limit <- ceiling(max(sizes) / 1000) * 1000
    for (j in 1:4) {
        w <- sizes[j] / limit * span
        box(end - w / 2, sy[j], w, .042, "#BDBDBD")
        label(format_number(sizes[j]), end - .004, sy[j], 8.5, just = "right")
        label(c("D. melanogaster", "D. simulans", "D. sechellia", "D. erecta")[j], .22, sy[j], 8.5, just = "left", fontface = 3)
    }
    line(end - span, .052, end, .052)
    for (v in c(0, 2500)) {
        x <- end - v / limit * span
        line(x, .052, x, .044)
        label(v, x, .031, 7.5)
    }
    label("Editing sites", end - span / 2, .006, 8)
    tx <- function(v) .006 + v / 10 * .067
    ty <- function(v) sy[4] + v * (sy[1] - sy[4]) / 3
    seg <- function(a, b, c, d) line(tx(a), ty(b), tx(c), ty(d))
    seg(1.3, 0, 9.1, 0)
    seg(1.3, 0, 1.3, 2.25)
    seg(1.3, 2.25, 4.3, 2.25)
    seg(4.3, 1.5, 4.3, 3)
    seg(4.3, 3, 9.1, 3)
    seg(4.3, 1.5, 7.1, 1.5)
    seg(7.1, 1, 7.1, 2)
    seg(7.1, 2, 9.1, 2)
    seg(7.1, 1, 9.1, 1)
    line(tx(1), ty(3.6), tx(9), ty(3.6), "#1F5AA6")
    for (j in 1:3) {
        x <- tx(c(1, 5, 9)[j])
        line(x, ty(3.6), x, ty(3.76), "#1F5AA6")
        label(c("10", "5", "0")[j], x, ty(3.95), 7, col = "#1F5AA6")
    }
    label("MYA", tx(5), ty(4.22), 7.5, col = "#1F5AA6")
}

# The hclust object and ordered matrix are the original R results, unchanged.
draw_B <- function() {
    n <- length(ordered_handles)
    left <- .025
    right <- .755
    bottom <- .18
    top <- .73
    xp <- function(pos) left + (pos - .5) / n * (right - left)
    maxh <- max(hc$height)
    yh <- function(h) top + .012 + h / maxh * .205
    xpos <- numeric(n - 1)
    height <- numeric(n - 1)
    leafpos <- match(seq_len(n), hc$order)
    for (i in seq_len(n - 1)) {
        a <- hc$merge[i, 1]
        b <- hc$merge[i, 2]
        xa <- if (a < 0) leafpos[-a] else xpos[a]
        xb <- if (b < 0) leafpos[-b] else xpos[b]
        ha <- if (a < 0) 0 else height[a]
        hb <- if (b < 0) 0 else height[b]
        xpos[i] <- (xa + xb) / 2
        height[i] <- hc$height[i]
        line(xp(xa), yh(ha), xp(xa), yh(height[i]), lwd = .35)
        line(xp(xb), yh(hb), xp(xb), yh(height[i]), lwd = .35)
        line(xp(xa), yh(height[i]), xp(xb), yh(height[i]), lwd = .35)
    }
    pal <- colorRampPalette(c("#48489F", "#F7F7F7", "#FEC17D"))(256)
    vals <- t(ordered_matrix)
    cols <- matrix(pal[pmax(1, pmin(256, round(vals * 255) + 1))], nrow = nrow(vals), ncol = ncol(vals))
    grid.raster(as.raster(cols), x = (left + right) / 2, y = (top + bottom) / 2, width = right - left, height = top - bottom, interpolate = FALSE)
    ys <- top - (seq_along(heatmap_sample_order) - .5) / length(heatmap_sample_order) * (top - bottom)
    reps <- heat_metadata$n_replicates[match(heatmap_sample_order, heat_metadata$sample_id)]
    titles <- c("Adar5G1 M", "iso-1 M", "iso-1 F", "WT M", "WT F", "Dsim M", "Dsim F", "Dsec M", "Dsec F", "Dere M", "Dere F")
    for (i in seq_along(ys)) label(paste0(titles[i], " (n = ", reps[i], ")"), right + .014, ys[i], 8.5, just = "left")
    for (i in c(1, 5, 7, 9)) {
        y <- top - i / 11 * (top - bottom)
        line(left, y, right, y, "white", .75)
    }
    for (j in seq_along(group_counts)) {
        inds <- which(ordered_group == paste0("Group ", j))
        l <- left + (min(inds) - 1) / n * (right - left)
        r <- left + max(inds) / n * (right - left)
        box((l + r) / 2, bottom - .023, r - l, .023, c("#F38EA5", "#8BC63F")[j])
        label(paste0("Cluster ", j, "\n(n = ", format_number(length(inds)), ")"), if (j == 1) l else (l + r) / 2, bottom - .075, 9, just = if (j == 1) "left" else "centre")
    }
    for (j in 1:4) {
        ii <- list(2:5, 6:7, 8:9, 10:11)[[j]]
        hi <- max(ys[ii]) + (top - bottom) / 22
        lo <- min(ys[ii]) - (top - bottom) / 22
        line(.966, lo, .966, hi, lwd = .55)
        line(.958, lo, .966, lo, lwd = .55)
        line(.958, hi, .966, hi, lwd = .55)
        label(c("Dmel", "Dsim", "Dsec", "Dere")[j], .978, (hi + lo) / 2, 8, rot = 90)
    }
    grid.raster(as.raster(matrix(rev(pal), nrow = 1)), x = .838, y = .84, width = .12, height = .028, interpolate = FALSE)
    for (j in 1:3) label(c("1.0", "0.5", "0.0")[j], .778 + (j - 1) * .06, .875, 8)
    label(expression(italic(L)[m]), .925, .84, 10)
}

# Panel C uses exactly the original site means / quartiles / polygons; formatting
# is matched to the supplied PDF (two decimals, WT label, ticks through 1.0).
pretty_sample <- function(id, n) {
    parts <- strsplit(id, "_", fixed = TRUE)[[1]]
    nm <- c(ISO1 = "iso-1", WBG = "WT", sim = "Dsim", sec = "Dsec", ere = "Dere")[[parts[1]]]
    parse(text = sprintf("atop(plain('%s')~plain('%s'),group('(',italic(n)==%d,')'))", nm, parts[2], n))[[1]]
}
cs <- copy(stats)
cs[, summary_label := paste0("(", sprintf("%.2f", mean), " \u00b1 ", sprintf("%.2f", sd), ", ", format_number(sites), " sites)")]
ys_label <- as.expression(Map(pretty_sample, sample_meta$sample_id, sample_meta$n_replicates))
p_C <- ggplot() +
    geom_polygon(data = violin_polygons, aes(x = x, y = y, group = interaction(sample_id, ns_class), fill = ns_class), colour = NA, alpha = .92) +
    geom_segment(data = cs, aes(x = whisker_low, xend = whisker_high, y = box_y, yend = box_y), linewidth = .32) +
    geom_rect(data = cs, aes(xmin = q1, xmax = q3, ymin = box_y - .060, ymax = box_y + .060), fill = "white", colour = "black", linewidth = .30) +
    geom_segment(data = cs, aes(x = median, xend = median, y = box_y - .060, yend = box_y + .060), linewidth = .34) +
    geom_point(data = cs, aes(x = mean, y = box_y), size = 1.4) +
    geom_text(data = cs, aes(x = 1.03, y = box_y, label = summary_label), hjust = 0, size = 2.65, family = font_family) +
    scale_fill_manual(values = palette_ns, breaks = c("S", "N"), name = NULL) +
    scale_x_continuous(limits = c(0, 1.88), breaks = seq(0, 1, .25), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(.45, 10.6), breaks = sample_meta$base_y, labels = ys_label, expand = expansion(mult = 0)) +
    labs(x = expression("Editing level (" * italic(L)[m] * ")"), y = NULL) +
    annotate("segment", x = 0, xend = 1, y = .45, yend = .45, linewidth = .4) +
    theme_classic(base_size = 10, base_family = font_family) +
    theme(
        axis.line.x = element_blank(), axis.text.y = element_text(size = 9, colour = "black", margin = margin(r = 6)), axis.text.x = element_text(size = 9), axis.title.x = element_text(size = 10),
        legend.position = "top", legend.direction = "horizontal", legend.background = element_rect(fill = "white", colour = "grey60", linetype = "dashed"),
        legend.key.size = unit(.28, "cm"), plot.margin = margin(5, 7, 6, 5, "mm")
    )

save_grid_panel <- function(fun, path, w, h, tag = NULL) {
    cairo_pdf(path, width = w, height = h, family = font_family)
    grid.newpage()
    fun()
    if (!is.null(tag)) label(tag, .003, .993, 16, just = c("left", "top"), fontface = 2)
    dev.off()
}
save_grid_panel(draw_A, file.path(output_dir, "panel_A_upset_1v1_R.pdf"), 10, 6, "A")
save_grid_panel(draw_B, file.path(output_dir, "panel_B_heatmap_1v1_R.pdf"), 11, 6, "B")
save_plot(p_C, file.path(output_dir, "panel_C_NS_half_violin_1v1_R.pdf"), file.path(output_dir, "panel_C_NS_half_violin_1v1_R.png"), 8.9, 11.1)
all_prefix <- if (reclassify) "Figure1_all_site_types_violin_1v1_R_reclassified" else "Figure1_all_site_types_violin_1v1_R"
for (ext in c("pdf", "png")) {
    stopifnot(file.copy(file.path(output_dir, paste0(all_prefix, ".", ext)), file.path(figure_root, paste0("All_site_types.", ext)), overwrite = TRUE))
    # These redundant files were just generated by this run, and are now renamed.
    unlink(file.path(output_dir, paste0(all_prefix, ".", ext)))
}
if (dev.cur() > 1) dev.off()
unlink(scratch_device)
cat("[DONE] Figure 1 panels: ", figure_root, "\n", sep = "")
