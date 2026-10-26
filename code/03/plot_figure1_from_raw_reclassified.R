#!/usr/bin/env Rscript

# Rebuild the 1v1 Figure 1 panels directly from the raw plotting tables.
#
# Required data files in DATA_DIR:
#   candidate.editing_level.observed.long.tsv
#   ISO1.annotation.tsv, WBG.annotation.tsv
#   sim.annotation.tsv, sec.annotation.tsv, ere.annotation.tsv
#
# The primary annotation remains species specific, with two display
# display classification: only N, S, 5'UTR, 3'UTR, Intron, and ncRNA are kept
# as named classes; every other annotation is pooled into 'Others'. A
# D. sechellia site in Others is provisionally displayed as 3'UTR/5'UTR when
# its homologous D. melanogaster annotation in ISO1_type is 3'UTR/5'UTR. Every
# such reassignment is written to an audit TSV.
#
# Usage:
#   Rscript plot_figure1_from_raw.R DATA_DIR OUTPUT_DIR

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(UpSetR)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
    getwd()
}

data_dir <- if (length(args) >= 1L) {
    args[[1L]]
} else {
    file.path(script_dir, "../data")
}
output_dir <- if (length(args) >= 2L) {
    args[[2L]]
} else {
    file.path(script_dir, "../result/R_run_current")
}
data_dir <- normalizePath(data_dir, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
table_dir <- file.path(output_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
group_order <- c("ISO1", "WBG", "sim", "sec", "ere")
sample_order <- c("ISO1_M", "ISO1_F", "WBG_M", "WBG_F", "sim_M", "sim_F", "sec_M", "sec_F",
    "ere_M", "ere_F")
heatmap_sample_order <- c("Adar5G1_M", sample_order)
variant_class_order <- c("N", "S", "5'UTR", "3'UTR", "Intron", "ncRNA", "Others")
variant_class_colors <- c(N = "#CC3A63", S = "#D2848D", `5'UTR` = "#E99C81", `3'UTR` = "#F5DB98",
    Intron = "#7BBCD5", ncRNA = "#6E7CB8", Others = "#A48CBE")
palette_ns <- c(N = "#E97274", S = "#58B6D9")
group_colors <- c(`Group 1` = "#E69F00", `Group 2` = "#2CA25F")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

stop_missing <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) {
        stop("Missing or empty file: ", path)
    }
}

require_columns <- function(dt, required, source_file) {
    missing <- setdiff(required, names(dt))
    if (length(missing) > 0L) {
        stop("Missing column(s) in ", source_file, ": ", paste(missing, collapse = ", "))
    }
}

format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L)
        return(NA_real_)
    sd(x)
}

format_sd <- function(x) {
    ifelse(is.finite(x), sprintf("%.3f", x), "NA")
}

classify_type <- function(type) {
    type <- fcoalesce(as.character(type), "")
    fcase(grepl("missense_variant", type, ignore.case = TRUE), "N", grepl("synonymous_variant",
        type, ignore.case = TRUE), "S", grepl("3_prime_UTR_variant", type, ignore.case = TRUE),
        "3'UTR", grepl("5_prime_UTR_variant", type, ignore.case = TRUE), "5'UTR", grepl("intron_variant",
            type, ignore.case = TRUE), "Intron", grepl("non_coding", type, ignore.case = TRUE),
        "ncRNA", default = "Others")
}

derive_sex <- function(library) {
    library <- as.character(library)
    fifelse(grepl("(^|[_-])F([0-9]+)?($|[_-])", library, ignore.case = TRUE), "F", fifelse(grepl("(^|[_-])M([0-9]+)?($|[_-])",
        library, ignore.case = TRUE), "M", NA_character_))
}

sample_plotmath <- function(sample_id) {
    prefix <- sub("_.*$", "", sample_id)
    suffix <- sub("^.*_", "", sample_id)
    if (prefix == "ISO1")
        return(sprintf("plain('iso-1')~plain('%s')", suffix))
    if (prefix == "WBG")
        return(sprintf("plain('WBG')~plain('%s')", suffix))
    if (prefix %chin% c("sim", "sec", "ere")) {
        return(sprintf("italic(%s)~plain('%s')", prefix, suffix))
    }
    sprintf("plain('%s')", sample_id)
}

facet_label_code <- function(sample_id, n_replicates, total_sites) {
    sprintf(paste0("atop(%s, paste('(', italic(n), ", "' = %s, sites = %s', ')'))"), sample_plotmath(sample_id),
        n_replicates, format_number(total_sites))
}

axis_label_expression <- function(sample_id, n_replicates) {
    parse(text = sprintf("atop(%s, group('(', italic(n)==%s, ')'))", sample_plotmath(sample_id),
        n_replicates))[[1L]]
}

save_plot <- function(plot_object, pdf_file, png_file, width, height) {
    if (!capabilities("cairo"))
        stop("Cairo support is required for PDF output.")
    cairo_pdf(pdf_file, width = width, height = height, family = font_family)
    print(plot_object)
    invisible(dev.off())
    grDevices::png(png_file, width = width, height = height, units = "in", res = 300, type = "cairo")
    print(plot_object)
    invisible(dev.off())
}

save_grob <- function(grob_object, pdf_file, width, height) {
    cairo_pdf(pdf_file, width = width, height = height, family = font_family)
    grid.newpage()
    grid.draw(grob_object)
    invisible(dev.off())
}

make_half_violin <- function(values, baseline, direction, max_width = 0.27) {
    values <- values[is.finite(values)]
    if (length(values) == 0L)
        return(data.table(x = numeric(), y = numeric()))
    if (length(values) < 2L || uniqueN(values) < 2L) {
        x0 <- values[[1L]]
        return(data.table(x = c(x0 - 0.002, x0 + 0.002, x0), y = c(baseline, baseline, baseline +
            direction * 0.03)))
    }
    lower <- max(0, min(values))
    upper <- min(1, max(values))
    if (upper <= lower)
        upper <- min(1, lower + 0.002)
    dens <- density(values, from = lower, to = upper, n = 384, cut = 0)
    scaled_width <- dens$y/max(dens$y) * max_width
    outer_y <- baseline + direction * scaled_width
    data.table(x = c(dens$x, rev(dens$x)), y = c(rep(baseline, length(dens$x)), rev(outer_y)))
}


# -----------------------------------------------------------------------------
# Read raw plotting data and species-specific annotations
# -----------------------------------------------------------------------------

quant_file <- file.path(data_dir, "candidate.editing_level.observed.long.tsv")
stop_missing(quant_file)
quant <- fread(quant_file, showProgress = FALSE)
require_columns(quant, c("group", "library", "handle", "editing_level"), quant_file)
quant[, group := as.character(group)]
quant[, library := as.character(library)]
quant[, handle := as.character(handle)]
quant[, editing_level := suppressWarnings(as.numeric(editing_level))]
if (!"sex" %in% names(quant)) quant[, sex := NA_character_]
quant[, sex := toupper(as.character(sex))]
quant[!sex %chin% c("F", "M"), sex := derive_sex(library)]

annotation_paths <- file.path(data_dir, paste0(group_order, ".annotation.tsv"))
names(annotation_paths) <- group_order
invisible(lapply(annotation_paths, stop_missing))

annotation_list <- lapply(group_order, function(group_name) {
    path <- annotation_paths[[group_name]]
    dt <- fread(path, showProgress = FALSE)
    require_columns(dt, c("handle", "type"), path)
    if (group_name == "sec") {
        require_columns(dt, "ISO1_type", path)
    }
    dt[, handle := as.character(handle)]
    dt[, type := as.character(type)]
    dt <- dt[nzchar(handle) & nzchar(type)]
    dt[, group := group_name]
    if (!"ISO1_type" %in% names(dt))
        dt[, ISO1_type := NA_character_]
    dt[, ISO1_type := as.character(ISO1_type)]
    dt[, original_type_class := classify_type(type)]
    dt[, dmel_homolog_type_class := NA_character_]
    if (group_name == "sec") {
        dt[, dmel_homolog_type_class := classify_type(ISO1_type)]
    }
    dt[, type_class := original_type_class]
    dt[group == "sec" & original_type_class == "Others" & dmel_homolog_type_class %chin%
        c("5'UTR", "3'UTR"), type_class := dmel_homolog_type_class]
    unique(dt[, .(group, handle, type, ISO1_type, original_type_class, dmel_homolog_type_class,
        type_class)], by = c("group", "handle"))
})
annotations <- rbindlist(annotation_list, use.names = TRUE)

sec_utr_reclassification <- unique(annotations[group == "sec" & original_type_class == "Others" &
    type_class %chin% c("5'UTR", "3'UTR"), .(handle, sec_type = type, dmel_ISO1_type = ISO1_type,
    original_type_class, dmel_homolog_type_class, final_type_class = type_class)])
fwrite(sec_utr_reclassification[order(final_type_class, handle)], file.path(table_dir, "sec_other_to_dmel_UTR.reclassification.tsv"),
    sep = "\t")

species_handles <- lapply(group_order, function(group_name) {
    unique(annotations[group == group_name, handle])
})
names(species_handles) <- group_order
dmel_handles <- sort(unique(c(species_handles[["ISO1"]], species_handles[["WBG"]])))
upset_handles <- list(Dmel = dmel_handles, Dsim = species_handles[["sim"]], Dsec = species_handles[["sec"]],
    Dere = species_handles[["ere"]])
union_handles <- sort(unique(unlist(upset_handles, use.names = FALSE)))
shared_handles <- sort(Reduce(intersect, upset_handles))

quant <- quant[group %chin% group_order & handle %chin% union_handles & is.finite(editing_level) &
    sex %chin% c("F", "M")]
quant <- quant[group != "ISO1" | library %chin% c("ISO1_head_F", "ISO1_head_M")]
quant[, sample_id := paste(group, sex, sep = "_")]

replicate_counts <- quant[, .(n_replicates = uniqueN(library)), by = .(sample_id, group,
    sex)]
replicate_counts[, sample_order__ := match(sample_id, sample_order)]
setorder(replicate_counts, sample_order__)
replicate_counts[, sample_order__ := NULL]
if (any(!sample_order %in% replicate_counts$sample_id)) {
    stop("One or more M/F sample groups are missing from the raw table.")
}

site_means <- quant[, .(site_mean = mean(editing_level, na.rm = TRUE)), by = .(group, sample_id,
    handle)]
site_means <- merge(site_means, annotations, by = c("group", "handle"), all.x = TRUE, sort = FALSE)
missing_annotation <- site_means[is.na(type_class), .(group, sample_id, handle)]
fwrite(missing_annotation, file.path(table_dir, "species_annotation_missing.site_mean.tsv"),
    sep = "\t")
site_means <- site_means[!is.na(type_class) & is.finite(site_mean)]
if (nrow(site_means) == 0L) stop("No annotated site means remain after the species merge.")

if (site_means[group == "sec" & type_class %chin% c("5'UTR", "3'UTR") & dmel_homolog_type_class !=
    type_class, .N] > 0L) {
    stop("A sec UTR display class is not supported by its Dmel homolog annotation.")
}

fwrite(site_means[order(group, sample_id, type_class, handle)], file.path(table_dir, "site_mean.by_species_type.tsv"),
    sep = "\t")
fwrite(unique(annotations)[order(group, handle)], file.path(table_dir, "species_specific_type_annotations.tsv"),
    sep = "\t")


# -----------------------------------------------------------------------------
# UpSet data and panel
# -----------------------------------------------------------------------------

membership <- data.table(handle = union_handles)
membership[, Dmel := as.integer(handle %chin% upset_handles[["Dmel"]])]
membership[, Dsim := as.integer(handle %chin% upset_handles[["Dsim"]])]
membership[, Dsec := as.integer(handle %chin% upset_handles[["Dsec"]])]
membership[, Dere := as.integer(handle %chin% upset_handles[["Dere"]])]
membership[, pattern := paste0(Dmel, Dsim, Dsec, Dere)]
fwrite(membership, file.path(table_dir, "upset.input.tsv"), sep = "\t")
fwrite(data.table(set = c("Dmel", "Dsim", "Dsec", "Dere"), n = c(sum(membership$Dmel), sum(membership$Dsim),
    sum(membership$Dsec), sum(membership$Dere))), file.path(table_dir, "upset.set_count.tsv"),
    sep = "\t")

grid.newpage()
UpSetR::upset(as.data.frame(membership[, .(Dmel, Dsim, Dsec, Dere)]), sets = c("Dmel", "Dsim",
    "Dsec", "Dere"), nsets = 4, nintersects = 8, order.by = "freq", decreasing = TRUE, keep.order = TRUE,
    main.bar.color = "#3B3B3B", sets.bar.color = "#3B3B3B", matrix.color = "#222222", point.size = 3,
    line.size = 0.7, text.scale = c(1.15, 1.05, 0.9, 0.9, 1, 0.9))
upset_grob <- grid.grab(wrap.grobs = TRUE)
save_grob(upset_grob, file.path(output_dir, "panel_A_upset_1v1_R.pdf"), width = 9, height = 5.6)


# -----------------------------------------------------------------------------
# Heatmap and in-strip group labels
# -----------------------------------------------------------------------------

heat_quant <- quant[handle %chin% shared_handles]
heat_mean <- heat_quant[, .(editing_level = mean(editing_level, na.rm = TRUE)), by = .(sample_id,
    handle)]
heat_grid <- CJ(sample_id = sample_order, handle = shared_handles, unique = TRUE)
heat_grid <- merge(heat_grid, heat_mean, by = c("sample_id", "handle"), all.x = TRUE, sort = FALSE)
heat_grid[is.na(editing_level), editing_level := 0]
heat_wide <- dcast(heat_grid, handle ~ sample_id, value.var = "editing_level")
heat_wide <- heat_wide[match(shared_handles, handle)]
heat_wide[, Adar5G1_M := 0]
heat_wide <- heat_wide[, c("handle", heatmap_sample_order), with = FALSE]
heat_matrix <- as.matrix(heat_wide[, -1])
storage.mode(heat_matrix) <- "numeric"
rownames(heat_matrix) <- heat_wide$handle

hc <- hclust(dist(heat_matrix), method = "complete")
leaf_order <- hc$order
ordered_handles <- rownames(heat_matrix)[leaf_order]
raw_group <- cutree(hc, k = 2L)
encountered <- unique(raw_group[leaf_order])
group_map <- setNames(seq_along(encountered), encountered)
ordered_group <- paste0("Group ", unname(group_map[as.character(raw_group[leaf_order])]))
group_counts <- table(factor(ordered_group, levels = c("Group 1", "Group 2")))

fwrite(data.table(handle = ordered_handles, ordered_position = seq_along(ordered_handles),
    group = ordered_group), file.path(table_dir, "heatmap.site_groups.tsv"), sep = "\t")
fwrite(data.table(handle = ordered_handles), file.path(table_dir, "shared_four_species_heatmap_1v1.handles.txt"),
    sep = "\t", col.names = FALSE)

ordered_matrix <- heat_matrix[leaf_order, , drop = FALSE]
heat_plot_dt <- data.table(handle = ordered_handles)
for (sample_id in heatmap_sample_order) {
    heat_plot_dt[, (sample_id) := ordered_matrix[, sample_id]]
}
heat_long <- melt(heat_plot_dt, id.vars = "handle", variable.name = "sample_id", value.name = "editing_level")
heat_long[, handle := factor(handle, levels = ordered_handles)]
heat_long[, sample_id := factor(sample_id, levels = rev(heatmap_sample_order))]

heat_metadata <- copy(replicate_counts[, .(sample_id, n_replicates)])
heat_metadata <- rbind(data.table(sample_id = "Adar5G1_M", n_replicates = 2L), heat_metadata,
    fill = TRUE)
heat_metadata <- heat_metadata[match(heatmap_sample_order, sample_id)]
heat_labels <- paste0(heat_metadata$n_replicates, " replicates")
names(heat_labels) <- heat_metadata$sample_id

p_heatmap <- ggplot(heat_long, aes(x = handle, y = sample_id, fill = editing_level)) + geom_raster() +
    scale_fill_gradientn(colours = c("#41469B", "#F7F7F7", "#E97274", "#8F1D2C"), limits = c(0,
        1), name = "Editing level") + scale_y_discrete(labels = function(x) {
    paste0(x, " (n = ", heat_labels[x], ")")
}) + labs(title = paste0("Mean editing level across four-species shared 1v1 sites (n = ",
    format_number(length(shared_handles)), ")"), x = "Clustered 1v1 site order", y = NULL) +
    theme_classic(base_size = 10, base_family = font_family) + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 8.5),
    axis.line = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
    legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 8),
    plot.title = element_text(size = 11, hjust = 0.5), plot.margin = margin(2, 5, 5, 5, unit = "mm"))

group_positions <- data.table(group = c("Group 1", "Group 2"), xmin = c(1, sum(group_counts[["Group 1"]]) +
    1), xmax = c(sum(group_counts[["Group 1"]]), length(ordered_handles)))
group_positions[, center := (xmin + xmax)/2]
group_positions[, label := paste0(group, " (n = ", format_number(as.integer(group_counts[group])),
    ")")]

p_group <- ggplot() + geom_rect(data = group_positions, aes(xmin = xmin - 0.5, xmax = xmax +
    0.5, ymin = 0, ymax = 1, fill = group), colour = "white", linewidth = 0.25) + geom_text(data = group_positions,
    aes(x = center, y = 0.5, label = label), colour = c("black", "white"), fontface = "bold",
    size = 3) + scale_fill_manual(values = group_colors, guide = "none") + scale_x_continuous(limits = c(0.5,
    length(ordered_handles) + 0.5), expand = c(0, 0)) + scale_y_continuous(limits = c(0,
    1), expand = c(0, 0)) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0,
    5, 0, 40, unit = "pt"))

heatmap_grob <- arrangeGrob(ggplotGrob(p_group), ggplotGrob(p_heatmap), ncol = 1, heights = c(0.12,
    1))
save_grob(heatmap_grob, file.path(output_dir, "panel_B_heatmap_1v1_R.pdf"), width = 11, height = 5.8)


# -----------------------------------------------------------------------------
# All type violin: one facet per sex/sample, all species in one figure
# -----------------------------------------------------------------------------

total_sites <- site_means[, .(total_sites = uniqueN(handle)), by = sample_id]
total_sites <- merge(total_sites, replicate_counts[, .(sample_id, n_replicates)], by = "sample_id",
    all.x = TRUE)
total_sites[, sample_order__ := match(sample_id, sample_order)]
setorder(total_sites, sample_order__)
total_sites[, sample_order__ := NULL]
total_sites[, facet_label := mapply(facet_label_code, sample_id, n_replicates, total_sites,
    USE.NAMES = FALSE)]

facet_levels <- total_sites$facet_label
site_means[, facet_label := factor(total_sites$facet_label[match(sample_id, total_sites$sample_id)],
    levels = facet_levels)]
site_means[, variant_class := factor(type_class, levels = variant_class_order)]

class_summary <- site_means[, .(sites = uniqueN(handle), mean = mean(site_mean), sd = safe_sd(site_mean)),
    by = .(sample_id, facet_label, variant_class)]
class_summary[, summary_label := paste0("n=", format_number(sites), "; mean=", sprintf("%.3f",
    mean), "; sd=", format_sd(sd))]

p_all_types <- ggplot(site_means, aes(x = variant_class, y = site_mean, fill = variant_class)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.28, colour = "grey25") +
    geom_boxplot(width = 0.13, fill = "white", alpha = 0.94, outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.055, height = 0, colour = "black", size = 0.27, alpha = 0.12) +
    stat_summary(fun = mean, geom = "point", shape = 16, size = 1.35, colour = "black") +
    geom_text(data = class_summary, aes(x = variant_class, y = 1.23, label = summary_label),
        inherit.aes = FALSE, family = font_family, size = 1.8, angle = 90, hjust = 0.5, vjust = 0.5) +
    facet_wrap(~facet_label, ncol = 2, scales = "free_x", drop = FALSE, labeller = label_parsed) +
    scale_fill_manual(values = variant_class_colors, limits = variant_class_order, drop = TRUE) +
    scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1, by = 0.25), labels = sprintf("%.2f",
        seq(0, 1, by = 0.25)), expand = expansion(mult = c(0, 0))) + labs(x = "Variant type",
    y = expression("Editing level (" * italic(L)[italic(m)] * ")")) + theme_classic(base_size = 10,
    base_family = font_family) + theme(strip.background = element_blank(), strip.text = element_text(size = 10.2,
    lineheight = 1.02, margin = margin(3.5, 0, 5, 0)), axis.text.x = element_text(angle = 32,
    hjust = 1, vjust = 1, size = 7.8), axis.text.y = element_text(size = 8.5), axis.title = element_text(size = 10.5),
    legend.position = "none", panel.grid = element_blank(), panel.spacing = unit(0.78, "lines"),
    plot.margin = margin(7, 7, 7, 7, unit = "mm"))

save_plot(p_all_types, file.path(output_dir, "Figure1_all_site_types_violin_1v1_R_reclassified.pdf"),
    file.path(output_dir, "Figure1_all_site_types_violin_1v1_R_reclassified.png"), width = 11.4,
    height = 15)
fwrite(class_summary[, .(sample_id = as.character(sample_id), variant_class = as.character(variant_class),
    mean, sd, sites)], file.path(table_dir, "p1_multitype_summary.tsv"), sep = "\t")


# -----------------------------------------------------------------------------
# Figure C: N/S half-violins
# -----------------------------------------------------------------------------

ns_values <- site_means[as.character(type_class) %chin% c("N", "S"), .(sample_id, handle,
    site_mean, ns_class = as.character(type_class))]
sample_meta <- copy(replicate_counts[, .(sample_id, n_replicates)])
sample_meta[, sample_order__ := match(sample_id, sample_order)]
setorder(sample_meta, sample_order__)
sample_meta[, sample_order__ := NULL]
sample_meta[, base_y := rev(seq_along(sample_order))]
sample_meta[, axis_label := Map(axis_label_expression, sample_id, n_replicates)]

ns_values <- merge(ns_values, sample_meta[, .(sample_id, base_y)], by = "sample_id", all.x = TRUE)
ns_values[, direction := fifelse(ns_class == "N", 1, -1)]
ns_values[, violin_baseline := base_y + direction * 0.1]
ns_values[, box_y := base_y + direction * 0.22]

stats <- ns_values[, {
    x <- site_mean[is.finite(site_mean)]
    q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
    iqr <- q[[3L]] - q[[1L]]
    lower <- x[x >= q[[1L]] - 1.5 * iqr]
    upper <- x[x <= q[[3L]] + 1.5 * iqr]
    if (!length(lower))
        lower <- x
    if (!length(upper))
        upper <- x
    list(sites = uniqueN(handle), mean = mean(x), sd = safe_sd(x), q1 = q[[1L]], median = q[[2L]],
        q3 = q[[3L]], whisker_low = min(lower), whisker_high = max(upper), box_y = unique(box_y))
}, by = .(sample_id, ns_class, base_y, direction, violin_baseline)]
stats[, summary_label := paste0("(", sprintf("%.3f", mean), " +/- ", format_sd(sd), ", n=",
    format_number(sites), " sites)")]
violin_polygons <- ns_values[, make_half_violin(site_mean, unique(violin_baseline), unique(direction)),
    by = .(sample_id, ns_class)]

axis_labels <- as.expression(sample_meta$axis_label)
x_ticks <- seq(0, 1.25, by = 0.25)
axis_y <- 0.52
p_panel_c <- ggplot() + geom_polygon(data = violin_polygons, aes(x = x, y = y, group = interaction(sample_id,
    ns_class), fill = ns_class), colour = NA, alpha = 0.92) + geom_segment(data = stats,
    aes(x = whisker_low, xend = whisker_high, y = box_y, yend = box_y), linewidth = 0.38) +
    geom_rect(data = stats, aes(xmin = q1, xmax = q3, ymin = box_y - 0.06, ymax = box_y +
        0.06), fill = "white", colour = "black", linewidth = 0.34) + geom_segment(data = stats,
    aes(x = median, xend = median, y = box_y - 0.06, yend = box_y + 0.06), linewidth = 0.42) +
    geom_point(data = stats, aes(x = mean, y = box_y), shape = 16, size = 1.65) + geom_text(data = stats,
    aes(x = 1.04, y = box_y, label = summary_label), hjust = 0, size = 3, family = font_family) +
    geom_segment(aes(x = 0, xend = 0, y = axis_y, yend = max(sample_meta$base_y) + 0.52),
        linewidth = 0.42) + geom_segment(aes(x = 0, xend = max(x_ticks), y = axis_y, yend = axis_y),
    linewidth = 0.42) + geom_segment(data = data.table(x = x_ticks), aes(x = x, xend = x,
    y = axis_y, yend = axis_y - 0.055), linewidth = 0.36) + geom_segment(data = sample_meta,
    aes(x = -0.018, xend = 0, y = base_y, yend = base_y), linewidth = 0.36) + scale_fill_manual(values = palette_ns,
    breaks = c("S", "N"), labels = c("S", "N"), name = NULL) + scale_x_continuous(limits = c(-0.02,
    1.72), breaks = x_ticks, labels = sprintf("%.2f", x_ticks), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0.42, max(sample_meta$base_y) + 0.64), breaks = sample_meta$base_y,
        labels = axis_labels, expand = expansion(mult = 0)) + labs(x = expression("Editing level (" *
    italic(L)[italic(m)] * ")"), y = NULL) + coord_cartesian(clip = "off") + theme_classic(base_size = 10,
    base_family = font_family) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 9.5, margin = margin(t = 4)), axis.text.y = element_text(size = 10.2,
        lineheight = 0.95, margin = margin(r = 8)), axis.title.x = element_text(size = 10.5,
        margin = margin(t = 7)), legend.position = "top", legend.direction = "horizontal",
    legend.justification = "center", legend.key = element_rect(fill = "white", colour = NA),
    legend.key.width = unit(0.42, "cm"), legend.key.height = unit(0.3, "cm"), legend.spacing.x = unit(0.08,
        "cm"), legend.text = element_text(size = 10), legend.background = element_rect(fill = "white",
        colour = "grey70", linewidth = 0.45, linetype = "dashed"), legend.margin = margin(2,
        4, 2, 4), plot.margin = margin(5, 8, 7, 13, unit = "mm"))

save_plot(p_panel_c, file.path(output_dir, "panel_C_NS_half_violin_1v1_R.pdf"), file.path(output_dir,
    "panel_C_NS_half_violin_1v1_R.png"), width = 8.9, height = 11.1)
fwrite(stats[, .(sample_id = as.character(sample_id), class = as.character(ns_class), mean,
    sd, sites)], file.path(table_dir, "panel_C_NS_summary.tsv"), sep = "\t")


cat("[DONE] output=", output_dir, "\n", sep = "")
cat("[INFO] Dmel=", length(species_handles[["ISO1"]]), " Dsim=", length(species_handles[["sim"]]),
    " Dsec=", length(species_handles[["sec"]]), " Dere=", length(species_handles[["ere"]]),
    " union=", length(union_handles), " shared=", length(shared_handles), "\n", sep = "")
cat("[INFO] sec types=", paste(sort(unique(as.character(site_means[group == "sec", type_class]))),
    collapse = ", "), "\n", sep = "")
cat("[INFO] sec Other-to-Dmel-UTR reclassified handles=", nrow(sec_utr_reclassification),
    " (5'UTR=", sec_utr_reclassification[final_type_class == "5'UTR", .N], ", 3'UTR=", sec_utr_reclassification[final_type_class ==
        "3'UTR", .N], ")\n", sep = "")
