#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

indir <- ifelse(length(args) >= 1, args[1], "/path/to/RNA_editing_DATAprocessing/05.psite/psite_2times")

outdir <- ifelse(length(args) >= 2, args[2], file.path("/path/to/RNA_editing_DATAprocessing/all_code/plots",
    "alledit_repmean_scatter.Asites_only_psite2times"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

run_date <- Sys.Date()
run_date_seed <- as.integer(Sys.getenv("ANALYSIS_SEED", "20260819"))
set.seed(run_date_seed)
fwrite(data.table(run_date = as.character(run_date), random_seed = run_date_seed), file.path(outdir,
    "Adar_autosite.run_metadata.tsv"), sep = "\t")

min_m_dp <- 3
min_r_dp <- 3
min_alt_rna <- 1

highlight_sites_by_group <- c(ISO1 = "FBtr0339272_1288", w1118 = "FBtr0339272_1288", Adar5G1 = "FBtr0339272_1288",
    WBG = "FBtr0339272_1288", AG = "FBtr0339272_1288", AS = "FBtr0339272_1288", sim = "FBtr0358708_1399",
    sec = "FBtr0202151_1165", ere = "FBtr0413873_1519")

# Keep the group-to-target mapping as a table for the sample-level join below.
target_dt <- data.table(group = names(highlight_sites_by_group), target_key = unname(highlight_sites_by_group))

highlight_label <- "Adar"

# ISO1 and WBG are both D. melanogaster and are pooled before taking the
# Female/Male species-level mean. AG and AS are deliberately excluded.
species_by_group <- c(ISO1 = "Dmel", WBG = "Dmel", sim = "Dsim", sec = "Dsec", ere = "Dere")

species_levels <- c("Dmel", "Dsim", "Dsec", "Dere")
selected_groups <- names(species_by_group)

# Explicit sample plan. The old *_head_ribo files are the fourth replicate
# for Dsim (F/M) and Dsec male; Dere is restricted to F1--F3 and M1--M3.
# ISO1 female/male are pooled with WBG F1--F3/M1--M3 as Dmel.
expected_meta <- rbindlist(list(CJ(group = "ISO1", sex = c("F", "M"), rep = "4", unique = TRUE),
    CJ(group = "WBG", sex = c("F", "M"), rep = as.character(1:3), unique = TRUE), CJ(group = "sim",
        sex = c("F", "M"), rep = as.character(1:4), unique = TRUE), data.table(group = "sec",
        sex = "F", rep = as.character(1:3)), data.table(group = "sec", sex = "M", rep = as.character(1:4)),
    CJ(group = "ere", sex = c("F", "M"), rep = as.character(1:3), unique = TRUE)), fill = TRUE)

replicate_plan <- expected_meta[, .(n_expected = .N, reps_expected = paste(sort(unique(rep)),
    collapse = ",")), by = .(group, sex)]

species_colors <- c(Dmel = "#4C78A8", Dsim = "#F58518", Dsec = "#54A24B", Dere = "#B279A2")

sex_shapes <- c(Female = 16, Male = 17)

mean_no_nan <- function(x) {
    y <- mean(x, na.rm = TRUE)
    if (is.nan(y))
        NA_real_ else y
}

format_p <- function(p) {
    if (is.na(p))
        return("NA")
    if (p < 2.2e-16)
        return(format(signif(p, 3), scientific = TRUE))
    format(signif(p, 3), scientific = TRUE)
}

paired_wilcox_summary <- function(dt) {
    keep <- is.finite(dt$Lm) & is.finite(dt$Lr)
    x <- dt$Lm[keep]
    y <- dt$Lr[keep]

    if (length(x) < 2L) {
        p_value <- NA_real_
    } else if (all(x == y)) {
        p_value <- 1
    } else {
        p_value <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value)
    }

    data.table(n_pair = length(x), mean_Lm = if (length(x))
        mean(x) else NA_real_, mean_Lr = if (length(y))
        mean(y) else NA_real_, paired_wilcox_p = p_value)
}

safe_name <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
}

normalize_group <- function(x) {
    z <- tolower(x)
    z <- gsub("_mono$", "", z)
    z <- gsub("-mono$", "", z)
    z <- gsub("_head$", "", z)
    z <- gsub("-head$", "", z)
    z <- gsub("^dmel[_-]", "", z)

    if (grepl("iso[-_]?1|iso1", z))
        return("ISO1")
    if (grepl("w1118", z))
        return("w1118")
    if (grepl("adar5g1", z))
        return("Adar5G1")
    if (grepl("adarg", z))
        return("AG")
    if (grepl("adars", z))
        return("AS")
    if (grepl("wbg", z))
        return("WBG")
    if (grepl("w501", z))
        return("sim")
    if (grepl("^d?sim", z))
        return("sim")
    if (grepl("^d?sec", z))
        return("sec")
    if (grepl("^d?ere|dere", z))
        return("ere")
    x
}

parse_meta_one <- function(f, source) {
    bn <- basename(f)

    clean <- sub("\\.(Asite|codon_Asite)_df_alledit$", "", bn, ignore.case = TRUE)

    group_raw <- NA_character_
    sex <- NA_character_
    rep <- NA_character_

    m <- regexec("^(.+?)[_-]([FM])([0-9]+)$", clean, perl = TRUE)
    hit <- regmatches(clean, m)[[1]]

    if (length(hit) >= 4) {
        group_raw <- hit[2]
        sex <- hit[3]
        rep <- hit[4]

    } else {

        group_raw <- clean

        ## sex
        if (grepl("female", tolower(clean))) {
            sex <- "F"
        } else if (grepl("male", tolower(clean))) {
            sex <- "M"
        } else if (grepl("w1118|adar5g1", tolower(clean))) {
            sex <- "M"
        }

        ## rep
        rep_match <- regexpr("rep([0-9]+)", clean, ignore.case = TRUE)
        if (rep_match > 0) {
            rep <- regmatches(clean, rep_match)[1]
            rep <- sub("rep", "", rep)
        } else {
            rep <- "4"
        }

        group_raw <- sub("_?head_?(female|male).*$", "", group_raw, ignore.case = TRUE)
        group_raw <- sub("_ribo.*$", "", group_raw, ignore.case = TRUE)
        group_raw <- sub("_rep[0-9]+.*$", "", group_raw, ignore.case = TRUE)
    }

    group_raw <- sub("_mono$", "", group_raw, ignore.case = TRUE)
    group_raw <- sub("-mono$", "", group_raw, ignore.case = TRUE)

    group <- normalize_group(group_raw)

    data.table(file = f, source = source, basename = bn, group = group, sex = sex, rep = rep)
}

list_asite_alledit <- function(dir) {
    files <- list.files(dir, pattern = "alledit$", full.names = TRUE, ignore.case = TRUE)

    if (length(files) == 0)
        return(character())

    bn <- basename(files)

    keep <- (grepl("Asite_df_alledit$", bn, ignore.case = TRUE) | grepl("codon_Asite_df_alledit$",
        bn, ignore.case = TRUE))

    keep <- keep & !grepl("Esite|Psite", bn, ignore.case = TRUE)
    # Choose one explicit file convention when site and codon tables coexist.
    suffix <- Sys.getenv("ADAR_ASITE_SUFFIX", ".Asite_df_alledit")
    keep <- keep & endsWith(bn, suffix)
    files[keep]
}

read_alledit_one <- function(f, group, sex, rep, sample_id, source, target_key) {
    message("Reading: ", basename(f), " -> ", sample_id)

    x <- fread(f, sep = "\t", header = TRUE, quote = "", fill = TRUE, check.names = FALSE)

    req_cols <- c("handle", "ref_number_rna", "alt_number_rna", "ref_number_ribo", "alt_number_ribo",
        "edit_pos", "FBt_cdspos")

    miss <- setdiff(req_cols, colnames(x))
    if (length(miss) > 0) {
        stop("Missing columns in file: ", f, "\n", paste(miss, collapse = ", "))
    }

    num_cols <- c("ref_number_rna", "alt_number_rna", "ref_number_ribo", "alt_number_ribo",
        "edit_pos")

    for (cc in num_cols) {
        x[, (cc) := as.numeric(get(cc))]
    }

    x[, site_key := fifelse(!is.na(FBt_cdspos) & FBt_cdspos != "", as.character(FBt_cdspos),
        as.character(handle))]

    x <- x[site_key == target_key]

    if (nrow(x) == 0) {
        stop("Target Adar site was not found in file: ", f, "\nExpected target key: ", target_key)
    }

    x[, m_dp := ref_number_rna + alt_number_rna]
    x[, r_dp := ref_number_ribo + alt_number_ribo]

    x[, valid_pair := (!is.na(m_dp) & !is.na(r_dp) & m_dp >= min_m_dp & r_dp >= min_r_dp &
        alt_number_rna >= min_alt_rna)]

    x[, Lm := fifelse(valid_pair, alt_number_rna/m_dp, NA_real_)]
    x[, Lr := fifelse(valid_pair, alt_number_ribo/r_dp, NA_real_)]

    setorder(x, site_key, -r_dp, -alt_number_ribo, -m_dp)
    x <- x[, .SD[1], by = site_key]

    x[, `:=`(group = group, sex = sex, rep = rep, sample_id = sample_id, source = source,
        file = basename(f))]

    x[, .(site_key, target_key, handle, FBt_cdspos, edit_pos, group, sex, rep, sample_id,
        source, file, ref_number_rna, alt_number_rna, ref_number_ribo, alt_number_ribo, m_dp,
        r_dp, valid_pair, Lm, Lr)]
}

make_plot <- function(d, title, out_prefix) {
    d <- d[is.finite(Lm) & is.finite(Lr)]
    n_use <- nrow(d)

    if (n_use >= 3) {
        wt <- tryCatch(wilcox.test(d$Lm, d$Lr, paired = TRUE, exact = FALSE), error = function(e) NULL)
        p_value <- if (is.null(wt))
            NA_real_ else wt$p.value
    } else {
        p_value <- NA_real_
    }

    p_str <- format_p(p_value)
    this_group <- unique(d$group)

    if (length(this_group) != 1) {
        warning("Cannot uniquely determine group for plot: ", title)
        this_highlight_site <- NA_character_
    } else if (this_group %in% names(highlight_sites_by_group)) {
        this_highlight_site <- highlight_sites_by_group[[this_group]]
    } else {
        this_highlight_site <- NA_character_
    }

    hit_dt <- d[site_key == this_highlight_site]

    p <- ggplot(d, aes(x = Lm, y = Lr)) + geom_point(shape = 16, size = 3.8, alpha = 0.9,
        color = "#4D4D4D") + xlab("Lm") + ylab("Lr") + scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1)) + ggtitle(title) + geom_abline(intercept = 0,
        slope = 1, color = "#b81316", linewidth = 1.2, linetype = 2) + theme_bw() + theme(axis.line = element_line(color = "black"),
        axis.title = element_text(size = 20, face = "plain", color = "black"), axis.text = element_text(size = 16,
            color = "black"), axis.ticks = element_line(linewidth = 1), axis.ticks.length = unit(0.25,
            "cm"), axis.title.y = element_text(margin = margin(0, 15, 0, 0)), axis.title.x = element_text(margin = margin(15,
            0, 0, 0)), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5), axis.line.x = element_line(color = "black",
            linewidth = 1.2), axis.line.y = element_line(color = "black", linewidth = 1.2)) +
        annotate("text", x = 0.8, y = 0.12, label = paste0("italic(n) == ", n_use), size = 6,
            color = "black", parse = TRUE) + annotate("text", x = 0.7, y = 0.05, label = paste0("italic(P) == '",
        p_str, "'"), size = 6, color = "black", parse = TRUE)

    if (!is.na(this_highlight_site) && nrow(hit_dt) > 0) {
        message("Highlight found in ", title, ": ", this_highlight_site)
        p <- p + geom_point(data = hit_dt, aes(x = Lm, y = Lr), size = 6, color = "#b81316") +
            geom_text(data = hit_dt, aes(x = Lm, y = Lr, label = highlight_label), color = "#b81316",
                vjust = 1, hjust = -0.4, size = 5)
    } else {
        message("Highlight not found in ", title, ": ", this_highlight_site)
    }

    out_pdf <- paste0(out_prefix, ".pdf")
    out_png <- paste0(out_prefix, ".png")

    if (capabilities("cairo")) {
        ggsave(out_pdf, p, width = 8, height = 8, device = cairo_pdf)
    } else {
        ggsave(out_pdf, p, width = 8, height = 8)
    }

    ggsave(out_png, p, width = 8, height = 8, dpi = 300)

    list(plot = p, summary = data.table(title = title, group = ifelse(length(this_group) ==
        1, this_group, NA_character_), n_site = n_use, p_value = p_value, p_method = "paired Wilcoxon signed-rank test after replicate mean",
        highlight_site = this_highlight_site, highlight_found = !is.na(this_highlight_site) &&
            nrow(hit_dt) > 0, pdf = out_pdf, png = out_png))
}

############################
############################

files <- list_asite_alledit(indir)

cat("Selected A-site alledit files:\n")
print(basename(files))

if (length(files) == 0) {
    stop("No A-site alledit files found in: ", indir)
}

meta_all <- rbindlist(lapply(files, parse_meta_one, source = "auto"), fill = TRUE)

meta_all[, rep := as.character(rep)]
fwrite(meta_all, file.path(outdir, "sample_metadata_all_files.tsv"), sep = "\t", quote = FALSE)
fwrite(expected_meta, file.path(outdir, "Adar_autosite_expected_sample_plan.tsv"), sep = "\t",
    quote = FALSE)
fwrite(replicate_plan, file.path(outdir, "Adar_autosite_replicate_plan.tsv"), sep = "\t",
    quote = FALSE)

meta_candidates <- meta_all[group %in% selected_groups & sex %in% c("F", "M")]

meta_counts <- meta_candidates[, .N, by = .(group, sex, rep)]
missing_meta <- expected_meta[!meta_counts, on = .(group, sex, rep)]
duplicate_meta <- meta_counts[N > 1]

fwrite(missing_meta, file.path(outdir, "Adar_autosite_missing_expected_rep.tsv"), sep = "\t",
    quote = FALSE)
fwrite(duplicate_meta, file.path(outdir, "Adar_autosite_duplicate_expected_rep.tsv"), sep = "\t",
    quote = FALSE)

if (nrow(missing_meta) > 0 || nrow(duplicate_meta) > 0) {
    stop("The expected species/sex/replicate sample plan is incomplete or duplicated. ",
        "See Adar_autosite_missing_expected_rep.tsv and ", "Adar_autosite_duplicate_expected_rep.tsv.")
}

# Keep only the specified samples. This excludes, for example,
# the available Dere_male_head_ribo file because Dere male is defined as M1--M3.
meta <- merge(expected_meta, meta_candidates, by = c("group", "sex", "rep"), all.x = TRUE,
    sort = FALSE)

meta <- merge(meta, target_dt[, .(group, target_key)], by = "group", all.x = TRUE, sort = FALSE)

meta[, sample_id := paste(group, sex, paste0("R", rep), sep = "_")]

bad_meta <- meta[is.na(group) | is.na(sex) | is.na(target_key)]
if (nrow(bad_meta) > 0) {
    fwrite(bad_meta, file.path(outdir, "BAD_selected_metadata.tsv"), sep = "\t", quote = FALSE)
    stop("Some selected files could not be parsed. See: ", file.path(outdir, "BAD_selected_metadata.tsv"))
}

fwrite(meta, file.path(outdir, "sample_metadata_selected_4species.tsv"), sep = "\t", quote = FALSE)

message("Files to process:")
print(meta[, .N, by = .(group, sex)])

############################
############################

rep_dt <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
    read_alledit_one(f = meta$file[i], group = meta$group[i], sex = meta$sex[i], rep = meta$rep[i],
        sample_id = meta$sample_id[i], source = meta$source[i], target_key = meta$target_key[i])
}), fill = TRUE)

fwrite(rep_dt, file.path(outdir, "alledit_rep_level_Lm_Lr.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

rep_summary <- rep_dt[, .(n_site = .N, n_valid = sum(valid_pair, na.rm = TRUE), n_Lm_Lr = sum(is.finite(Lm) &
    is.finite(Lr), na.rm = TRUE), mean_Lm = mean_no_nan(Lm), mean_Lr = mean_no_nan(Lr)),
    by = .(group, sex, rep, sample_id, file)]

fwrite(rep_summary, file.path(outdir, "alledit_rep_level_summary.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

############################
############################

mean_dt <- rep_dt[, .(Lm = mean_no_nan(Lm), Lr = mean_no_nan(Lr), n_rep_total = uniqueN(sample_id),
    n_rep_valid = uniqueN(sample_id[is.finite(Lm) & is.finite(Lr)]), reps_used = paste(sort(unique(sample_id[is.finite(Lm) &
        is.finite(Lr)])), collapse = ",")), by = .(group, sex, site_key)]

mean_dt <- mean_dt[is.finite(Lm) & is.finite(Lr)]
setorder(mean_dt, group, sex, site_key)

fwrite(mean_dt, file.path(outdir, "alledit_group_sex_repmean_Lm_Lr.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

############################
############################

plot_rep_dt <- rep_dt[group %in% selected_groups & sex %in% c("F", "M") & is.finite(Lm) &
    is.finite(Lr)]

# Assign the species before averaging so ISO1 and WBG contribute to the same
# Dmel Female/Male point rather than producing two points per sex.
plot_rep_dt[, species := unname(species_by_group[as.character(group)])]

rep_mean_plot_dt <- plot_rep_dt[, .(Lm = mean(Lm), Lr = mean(Lr), n_rep_valid = .N, reps_used = paste(sort(unique(rep)),
    collapse = ",")), by = .(species, sex)]

rep_mean_plot_dt[, species := factor(species, levels = species_levels)]

rep_mean_plot_dt <- rep_mean_plot_dt[!is.na(species) & is.finite(Lm) & is.finite(Lr)]

expected_plot <- CJ(species = species_levels, sex = c("F", "M"), unique = TRUE)

expected_plot[, species := factor(species, levels = species_levels)]

missing_plot <- expected_plot[!rep_mean_plot_dt, on = .(species, sex)]

duplicate_plot <- rep_mean_plot_dt[, .N, by = .(species, sex)][N > 1]

fwrite(missing_plot, file.path(outdir, "Adar_autosite_missing_final_8_points.tsv"), sep = "\t",
    quote = FALSE)
fwrite(duplicate_plot, file.path(outdir, "Adar_autosite_duplicate_final_8_points.tsv"), sep = "\t",
    quote = FALSE)

if (nrow(missing_plot) > 0 || nrow(duplicate_plot) > 0) {
    stop("The final plot requires exactly 8 species-sex pairs. See ", "Adar_autosite_missing_final_8_points.tsv and ",
        "Adar_autosite_duplicate_final_8_points.tsv.")
}

setorder(rep_mean_plot_dt, species, sex)

rep_mean_plot_dt[, sex_label := factor(sex, levels = c("F", "M"), labels = c("Female",
    "Male"))]

fwrite(rep_mean_plot_dt, file.path(outdir, "Adar_autosite_4species_2sex_repmean.tsv"), sep = "\t",
    quote = FALSE, na = "NA")

plot_stats <- paired_wilcox_summary(rep_mean_plot_dt)
plot_stats[, comparison := "Eight species-sex rep means"]
plot_stats[, p_method := "Two-sided paired Wilcoxon signed-rank test"]

fwrite(plot_stats, file.path(outdir, "Adar_autosite_4species_2sex_repmean_paired_wilcox.tsv"),
    sep = "\t", quote = FALSE, na = "NA")

long_plot_dt <- melt(rep_mean_plot_dt, id.vars = c("species", "sex", "sex_label"), measure.vars = c("Lm",
    "Lr"), variable.name = "level", value.name = "editing_level")

long_plot_dt[, level := factor(level, levels = c("Lm", "Lr"))]

fwrite(long_plot_dt, file.path(outdir, "Adar_autosite_4species_2sex_repmean_plot_data_long.tsv"),
    sep = "\t", quote = FALSE, na = "NA")

plot_y_max <- max(1.05, max(long_plot_dt$editing_level, na.rm = TRUE) + 0.08)

p_label <- paste0("P = ", format_p(plot_stats$paired_wilcox_p[1]), "\nn = ", plot_stats$n_pair[1],
    " pairs")

final_plot <- ggplot(long_plot_dt, aes(x = level, y = editing_level, group = interaction(species,
    sex), color = species, shape = sex_label)) + geom_line(linewidth = 0.9, alpha = 0.82) +
    geom_point(size = 3.8, stroke = 0.45) + annotate("text", x = 1.5, y = 0.13, label = p_label,
    size = 4.2, hjust = 0.5, vjust = 0) + scale_color_manual(values = species_colors, breaks = names(species_colors),
    name = "Species") + scale_shape_manual(values = sex_shapes, name = "Sex") + scale_x_discrete(labels = c(Lm = expression(L[m]),
    Lr = expression(L[r]))) + scale_y_continuous(limits = c(0, plot_y_max), breaks = seq(0,
    1, by = 0.2), expand = expansion(mult = c(0, 0.02))) + labs(title = "Adar A-site", x = NULL,
    y = "Editing level") + theme_classic(base_size = 14) + theme(plot.title = element_text(hjust = 0.5,
    size = 18, face = "plain"), axis.title.y = element_text(size = 15), axis.text.x = element_text(size = 16,
    color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.line = element_line(linewidth = 0.9,
    color = "black"), axis.ticks = element_line(linewidth = 0.8, color = "black"), legend.title = element_text(size = 11),
    legend.text = element_text(size = 10), legend.position = "right", plot.margin = margin(12,
        15, 10, 10))

final_pdf <- file.path(outdir, "Adar_autosite_4species_2sex_repmean_paired.pdf")
final_png <- file.path(outdir, "Adar_autosite_4species_2sex_repmean_paired.png")

ggsave(final_pdf, final_plot, width = 6.8, height = 6, device = if (capabilities("cairo")) cairo_pdf else NULL)

ggsave(final_png, final_plot, width = 6.8, height = 6, dpi = 300)

message("====================================================")
message("Final species-sex pairs: ", nrow(rep_mean_plot_dt))
message("Replicate plan: ", file.path(outdir, "Adar_autosite_replicate_plan.tsv"))
message("Paired Wilcoxon P: ", format_p(plot_stats$paired_wilcox_p[1]))
message("Final PDF: ", final_pdf)
message("Final PNG: ", final_png)
message("Output directory: ", outdir)
message("Done.")
