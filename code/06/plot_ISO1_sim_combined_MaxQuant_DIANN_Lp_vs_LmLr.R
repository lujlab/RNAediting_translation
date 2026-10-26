#!/usr/bin/env Rscript

# Combine DIA-NN and MaxQuant site-level Lp for ISO-1 and D. simulans, then
# compare the combined peptide editing proportion (Lp) with RNA Lm or Ribo Lr.
# Within each platform, 100AA and Trypsin are averaged with equal weight;
# afterwards DIA-NN and MaxQuant are averaged with equal weight.

suppressPackageStartupMessages(library(ggplot2))

root <- Sys.getenv("MS_ROOT")
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) {
    args[[1L]]
} else {
    file.path(root, "ISO1_sim_combined_MaxQuant_DIANN_Lp_scatter")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ipos_all <- "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/09_Ipos/Ipos_multispecies_results/species_sex_site_means_all.tsv"
ipos_all <- Sys.getenv("LMLR_TABLE", unset = ipos_all)
cdspos_dir <- Sys.getenv("CDSPOS_DIR")

inputs <- list(list(label = "mel (ISO-1)", genotype = "iso", dia = Sys.getenv("DIA_ISO",
    unset = file.path(root, "DIA_NN_analysis/dmel_iso_Lp_union.combined_sex_Lp_wide.longestTx_LmLr_all.tsv")),
    mq = file.path(root, "ISO1_MaxQuant_100AA_Trypsin_analysis/results/iso1_mean_maxquant_intensity.tsv"),
    coordinate = file.path(cdspos_dir, "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"),
    ipos_species = "Dmel"), list(label = "sim", genotype = "sim", dia = Sys.getenv("DIA_SIM",
    unset = file.path(root, "DIA_NN_analysis/dsim_Lp_union.combined_sex_Lp_wide.longestTx_LmLr_all.tsv")),
    mq = file.path(root, "SIM_MaxQuant_100AA_Trypsin_analysis/results/sim_mean_maxquant_intensity.tsv"),
    coordinate = file.path(cdspos_dir, "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"),
    ipos_species = "Dsim"))

as_num <- function(x) suppressWarnings(as.numeric(x))
mean_present <- function(x) {
    x <- as_num(x)
    if (all(is.na(x)))
        NA_real_ else mean(x, na.rm = TRUE)
}

split_source_ids <- function(x) {
    x <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
    trimws(x[nzchar(trimws(x))])
}

mq_lmlr_lookup <- function(cfg, handles) {
    # MaxQuant handles are source genomic site IDs, whereas Ipos values use
    # longest-transcript CDS coordinates.  Map source ID -> new coordinate and
    # accept both `handle` and alternate `cdspos_ids` in the Ipos table.
    coordinate <- read.delim(cfg$coordinate, check.names = FALSE, stringsAsFactors = FALSE)
    coordinate <- coordinate[coordinate$site_key %in% handles & !is.na(coordinate$new_transcript) &
        !is.na(coordinate$new_cds_pos_1based), c("site_key", "new_transcript", "new_cds_pos_1based"),
        drop = FALSE]
    coordinate$coord_key <- paste(coordinate$new_transcript, coordinate$new_cds_pos_1based,
        sep = "_")

    ipos <- read.delim(ipos_all, check.names = FALSE, stringsAsFactors = FALSE)
    if ("genotype" %in% colnames(ipos)) {
        ipos <- ipos[ipos$species == cfg$ipos_species & ipos$genotype == cfg$genotype & ipos$sex %in%
            c("female", "male"), c("sex", "handle", "mean_Lm", "mean_Lr"), drop = FALSE]
        ipos$cdspos_ids <- ipos$handle
    } else {
        ipos <- ipos[ipos$species == cfg$ipos_species & ipos$sex %in% c("female", "male"),
            c("sex", "handle", "cdspos_ids", "mean_Lm", "mean_Lr"), drop = FALSE]
    }
    ipos$mean_Lm <- as_num(ipos$mean_Lm)
    ipos$mean_Lr <- as_num(ipos$mean_Lr)
    # `cdspos_ids` can contain multiple positions within one codon.
    alternate <- do.call(rbind, lapply(seq_len(nrow(ipos)), function(i) {
        keys <- trimws(unlist(strsplit(as.character(ipos$cdspos_ids[i]), ",", fixed = TRUE)))
        data.frame(key = keys[nzchar(keys)], sex = ipos$sex[i], mean_Lm = ipos$mean_Lm[i],
            mean_Lr = ipos$mean_Lr[i], stringsAsFactors = FALSE)
    }))
    primary <- data.frame(key = ipos$handle, sex = ipos$sex, mean_Lm = ipos$mean_Lm, mean_Lr = ipos$mean_Lr,
        stringsAsFactors = FALSE)
    ipos_keyed <- unique(rbind(primary, alternate))

    out <- lapply(handles, function(handle) {
        coord_keys <- unique(coordinate$coord_key[coordinate$site_key == handle])
        matched <- ipos_keyed[ipos_keyed$key %in% coord_keys, , drop = FALSE]
        female <- matched[matched$sex == "female", , drop = FALSE]
        male <- matched[matched$sex == "male", , drop = FALSE]
        data.frame(handle = handle, mq_Lm = mean_present(c(mean_present(female$mean_Lm),
            mean_present(male$mean_Lm))), mq_Lr = mean_present(c(mean_present(female$mean_Lr),
            mean_present(male$mean_Lr))), mq_LmLr_coordinate_matches = length(unique(coord_keys)),
            stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

read_dia <- function(cfg) {
    d <- read.delim(cfg$dia, check.names = FALSE, stringsAsFactors = FALSE)
    d <- d[d$genotype == cfg$genotype, , drop = FALSE]
    if (!nrow(d))
        stop("No DIA-NN rows for ", cfg$label)

    rows <- lapply(seq_len(nrow(d)), function(i) {
        z <- d[i, , drop = FALSE]
        ids <- split_source_ids(z$source_site_ids)
        # Sex means are calculated within each method before equal-weight method averaging.
        dia_100aa <- mean_present(c(z$Lp_100AA_F, z$Lp_100AA_M))
        dia_trypsin <- mean_present(c(z$Lp_trypsin_F, z$Lp_trypsin_M))
        data.frame(handle = ids, dia_gene = z$gene, dia_transcript = z$FBtr, dia_aa_pos = as_num(z$aa_pos),
            dia_ref_aa = z$refAA, dia_alt_aa = z$known_alt_AAs, dia_Lp_100AA = dia_100aa,
            dia_Lp_trypsin = dia_trypsin, dia_Lp = mean_present(c(dia_100aa, dia_trypsin)),
            Lm = mean_present(c(z$mean_Lm_F, z$mean_Lm_M)), Lr = mean_present(c(z$mean_Lr_F,
                z$mean_Lr_M)), stringsAsFactors = FALSE)
    })
    x <- do.call(rbind, rows)
    # A source coordinate can map to more than one peptide/transcript row.
    # Keep one source site, averaging only duplicate mappings.
    do.call(rbind, lapply(split(x, x$handle), function(z) data.frame(handle = z$handle[1],
        dia_gene = z$dia_gene[1], dia_transcript = z$dia_transcript[1], dia_aa_pos = z$dia_aa_pos[1],
        dia_ref_aa = z$dia_ref_aa[1], dia_alt_aa = z$dia_alt_aa[1], dia_Lp_100AA = mean_present(z$dia_Lp_100AA),
        dia_Lp_trypsin = mean_present(z$dia_Lp_trypsin), dia_Lp = mean_present(z$dia_Lp),
        Lm = mean_present(z$Lm), Lr = mean_present(z$Lr), stringsAsFactors = FALSE)))
}

read_mq <- function(cfg) {
    d <- read.delim(cfg$mq, check.names = FALSE, stringsAsFactors = FALSE)
    total_trypsin <- as_num(d$total_intensity_Trypsin)
    total_100aa <- as_num(d$total_intensity_100AA)
    ratio_trypsin <- as_num(d$ratio_Trypsin)
    ratio_100aa <- as_num(d$ratio_100AA)
    # A zero ratio is valid when a reference peptide is observed, so use total
    # intensity (not Lp > 0) as the detection criterion.
    mq_trypsin <- ifelse(!is.na(total_trypsin) & total_trypsin > 0, ratio_trypsin, NA_real_)
    mq_100aa <- ifelse(!is.na(total_100aa) & total_100aa > 0, ratio_100aa, NA_real_)
    out <- data.frame(handle = d$handle, mq_gene = d$gene, mq_transcript = d$FBtr_id, mq_aa_pos = as_num(d$aa_pos),
        mq_ref_aa = d$ref_aa, mq_alt_aa = d$alt_aa, mq_Lp_100AA = mq_100aa, mq_Lp_trypsin = mq_trypsin,
        mq_Lp = vapply(seq_len(nrow(d)), function(i) mean_present(c(mq_100aa[i], mq_trypsin[i])),
            numeric(1)), stringsAsFactors = FALSE)
    merge(out, mq_lmlr_lookup(cfg, out$handle), by = "handle", all.x = TRUE)
}

combine_one <- function(cfg) {
    dia <- read_dia(cfg)
    mq <- read_mq(cfg)
    z <- merge(dia, mq, by = "handle", all = TRUE)
    z$species <- cfg$label
    z$DIA_NN_detected <- !is.na(z$dia_Lp)
    z$MaxQuant_detected <- !is.na(z$mq_Lp)
    # Use the Lm/Lr already matched to DIA-NN when available; for a genuine
    # MaxQuant-only site, use its direct source-site -> Ipos coordinate match.
    z$Lm <- ifelse(is.na(z$Lm), z$mq_Lm, z$Lm)
    z$Lr <- ifelse(is.na(z$Lr), z$mq_Lr, z$Lr)
    z$gene <- ifelse(!is.na(z$dia_gene) & nzchar(z$dia_gene), z$dia_gene, z$mq_gene)
    z$transcript <- ifelse(!is.na(z$dia_transcript) & nzchar(z$dia_transcript), z$dia_transcript,
        z$mq_transcript)
    z$aa_position <- ifelse(!is.na(z$dia_aa_pos), z$dia_aa_pos, z$mq_aa_pos)
    z$ref_aa <- ifelse(!is.na(z$dia_ref_aa) & nzchar(z$dia_ref_aa), z$dia_ref_aa, z$mq_ref_aa)
    z$alt_aa <- ifelse(!is.na(z$dia_alt_aa) & nzchar(z$dia_alt_aa), z$dia_alt_aa, z$mq_alt_aa)
    z$aa_substitution <- ifelse(!is.na(z$ref_aa) & nzchar(z$ref_aa) & !is.na(z$alt_aa) &
        nzchar(z$alt_aa), paste0(z$ref_aa, ">", z$alt_aa), NA_character_)
    z$Lp_combined <- vapply(seq_len(nrow(z)), function(i) {
        mean_present(c(z$dia_Lp[i], z$mq_Lp[i]))
    }, numeric(1))
    z$detection_group <- ifelse(z$DIA_NN_detected & z$MaxQuant_detected, "Both DIA-NN and MaxQuant",
        ifelse(z$DIA_NN_detected, "DIA-NN only", ifelse(z$MaxQuant_detected, "MaxQuant only",
            "Neither")))
    z
}

all_sites <- do.call(rbind, lapply(inputs, combine_one))
all_sites$detection_group <- factor(all_sites$detection_group, levels = c("Both DIA-NN and MaxQuant",
    "DIA-NN only", "MaxQuant only", "Neither"))
write.table(all_sites, file.path(outdir, "ISO1_sim_combined_Lp_site_table.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA")

annotated_cols <- c("species", "handle", "gene", "transcript", "aa_position", "aa_substitution",
    "Lp_combined", "dia_Lp", "mq_Lp", "detection_group", "Lm", "Lr")
write.table(all_sites[, annotated_cols], file.path(outdir, "ISO1_sim_combined_Lp_annotated_sites.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

mq_only_audit <- all_sites[all_sites$detection_group == "MaxQuant only" & is.finite(all_sites$Lp_combined),
    , drop = FALSE]
mq_only_audit$included_in_Lm_scatter <- is.finite(mq_only_audit$Lm)
mq_only_audit$included_in_Lr_scatter <- is.finite(mq_only_audit$Lr)
write.table(mq_only_audit, file.path(outdir, "ISO1_sim_MaxQuant_only_sites_LmLr_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

paired_stats <- function(d, x_col, y_col) {
    z <- d[is.finite(d[[y_col]]) & is.finite(d[[x_col]]), , drop = FALSE]
    if (nrow(z) < 2L)
        return(data.frame(n = nrow(z), p_value = NA_real_, delta_median = NA_real_))
    wt <- suppressWarnings(wilcox.test(z[[y_col]], z[[x_col]], paired = TRUE, exact = FALSE))
    data.frame(n = nrow(z), p_value = wt$p.value, delta_median = median(z[[y_col]] - z[[x_col]]))
}

make_plot <- function(y_col, y_label, x_col, x_label, file_stem) {
    d <- all_sites[is.finite(all_sites[[y_col]]) & is.finite(all_sites[[x_col]]), , drop = FALSE]
    # `Neither` is retained in the full audit table but cannot have a finite Lp
    # and therefore is not a plotted detection class.
    d$detection_group <- factor(as.character(d$detection_group), levels = c("Both DIA-NN and MaxQuant",
        "DIA-NN only", "MaxQuant only"))
    stats <- do.call(rbind, lapply(split(d, d$species), function(z) {
        s <- paired_stats(z, x_col, y_col)
        s$species <- z$species[1]
        s$DIA_NN_n <- sum(z$DIA_NN_detected)
        s$MaxQuant_n <- sum(z$MaxQuant_detected)
        s$Both_n <- sum(z$DIA_NN_detected & z$MaxQuant_detected)
        s
    }))
    stats$label <- sprintf("DIA-NN n = %d; MaxQuant n = %d\nBoth n = %d; paired Wilcoxon n = %d\nP = %s; Delta median = %.3f",
        stats$DIA_NN_n, stats$MaxQuant_n, stats$Both_n, stats$n, format.pval(stats$p_value,
            digits = 2, eps = 0.001), stats$delta_median)
    stats$panel_title <- paste(stats$species, stats$label, sep = "\n")
    d$panel_title <- stats$panel_title[match(d$species, stats$species)]
    write.table(stats, file.path(outdir, paste0(file_stem, "_paired_wilcox_stats.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    p <- ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]], colour = detection_group)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#777777") +
        geom_point(size = 1.8, alpha = 0.82) + facet_wrap(~panel_title, nrow = 1) + scale_colour_manual(values = c(`Both DIA-NN and MaxQuant` = "#1B9E77",
        `DIA-NN only` = "#4D4D4D", `MaxQuant only` = "#A8DDA0"), name = "Peptide detection",
        drop = FALSE) + coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + labs(x = x_label,
        y = y_label) + theme_classic(base_size = 8) + theme(strip.background = element_rect(fill = "grey95",
        colour = "grey70"), strip.text = element_text(size = 8, lineheight = 1.05), legend.position = "bottom",
        panel.spacing = grid::unit(1.1, "lines"))
    ggsave(file.path(outdir, paste0(file_stem, ".pdf")), p, width = 11, height = 5.5)
    ggsave(file.path(outdir, paste0(file_stem, ".png")), p, width = 11, height = 5.5, dpi = 300)
    list(data = d, stats = stats)
}

combined_lm_result <- make_plot("Lp_combined", "Combined peptide editing proportion (Lp)",
    "Lm", "RNA editing proportion (Lm)", "ISO1_sim_combined_Lp_vs_Lm_scatter")
combined_lr_result <- make_plot("Lp_combined", "Combined peptide editing proportion (Lp)",
    "Lr", "Ribo editing proportion (Lr)", "ISO1_sim_combined_Lp_vs_Lr_scatter")
dia_lm_result <- make_plot("dia_Lp", "DIA-NN peptide editing proportion (Lp)", "Lm", "RNA editing proportion (Lm)",
    "ISO1_sim_DIA_NN_Lp_vs_Lm_scatter")
dia_lr_result <- make_plot("dia_Lp", "DIA-NN peptide editing proportion (Lp)", "Lr", "Ribo editing proportion (Lr)",
    "ISO1_sim_DIA_NN_Lp_vs_Lr_scatter")
mq_lm_result <- make_plot("mq_Lp", "MaxQuant peptide editing proportion (Lp)", "Lm", "RNA editing proportion (Lm)",
    "ISO1_sim_MaxQuant_Lp_vs_Lm_scatter")
mq_lr_result <- make_plot("mq_Lp", "MaxQuant peptide editing proportion (Lp)", "Lr", "Ribo editing proportion (Lr)",
    "ISO1_sim_MaxQuant_Lp_vs_Lr_scatter")

# A compact A4 counterpart to the all-species Ar/Am-Ir/Im figure: four
# panels in one row (ISO-1 Lm, sim Lm, ISO-1 Lr, sim Lr), 8-pt typography.
make_a4_panel <- function(result, y_col, y_label, x_col, x_label, panel_title, species_label) {
    d <- result$data[result$data$species == species_label, , drop = FALSE]
    s <- result$stats[result$stats$species == species_label, , drop = FALSE]
    label <- sprintf("DIA n = %d; MQ n = %d\nBoth n = %d; paired n = %d\nP = %s; Delta med = %.3f",
        s$DIA_NN_n, s$MaxQuant_n, s$Both_n, s$n, format.pval(s$p_value, digits = 2, eps = 0.001),
        s$delta_median)
    ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]], colour = detection_group)) + geom_abline(slope = 1,
        intercept = 0, linetype = "dashed", colour = "#777777", linewidth = 0.55) + geom_point(size = 1.8,
        alpha = 0.82) + scale_colour_manual(values = c(`Both DIA-NN and MaxQuant` = "#1B9E77",
        `DIA-NN only` = "#4D4D4D", `MaxQuant only` = "#A8DDA0"), guide = "none", drop = FALSE) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + labs(title = panel_title,
        subtitle = label, x = x_label, y = y_label) + theme_classic(base_size = 8) + theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8), plot.title = element_text(size = 8, hjust = 0.5),
        plot.subtitle = element_text(size = 6.2, hjust = 0.5, lineheight = 1), plot.margin = margin(3,
            5, 3, 3))
}

write_a4_four_panels <- function(lm_result, lr_result, y_col, y_label, method_label, file_name) {
    panels <- list(make_a4_panel(lm_result, y_col, y_label, "Lm", "Lm", paste0("mel (ISO-1), ",
        method_label, ", Lm"), "mel (ISO-1)"), make_a4_panel(lm_result, y_col, y_label, "Lm",
        "Lm", paste0("sim, ", method_label, ", Lm"), "sim"), make_a4_panel(lr_result, y_col,
        y_label, "Lr", "Lr", paste0("mel (ISO-1), ", method_label, ", Lr"), "mel (ISO-1)"),
        make_a4_panel(lr_result, y_col, y_label, "Lr", "Lr", paste0("sim, ", method_label,
            ", Lr"), "sim"))
    pdf(file.path(outdir, file_name), width = 8.27, height = 11.69, onefile = TRUE)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.62, width = grid::unit(8.05, "inches"),
        height = grid::unit(2.6, "inches"), just = c("center", "center"), layout = grid::grid.layout(1,
            4)))
    for (i in seq_along(panels)) {
        print(panels[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
    }
    dev.off()
}

write_a4_four_panels(combined_lm_result, combined_lr_result, "Lp_combined", "Lp", "combined",
    "ISO1_sim_combined_Lp_vs_Lm_Lr_scatter_A4_4panels.pdf")
write_a4_four_panels(dia_lm_result, dia_lr_result, "dia_Lp", "DIA-NN Lp", "DIA-NN", "ISO1_sim_DIA_NN_Lp_vs_Lm_Lr_scatter_A4_4panels.pdf")
write_a4_four_panels(mq_lm_result, mq_lr_result, "mq_Lp", "MaxQuant Lp", "MaxQuant", "ISO1_sim_MaxQuant_Lp_vs_Lm_Lr_scatter_A4_4panels.pdf")

invisible(NULL)
