#!/usr/bin/env Rscript

# I-position analysis for multiple species and biological replicates.
#
# Usage:
#   Rscript test_Ipos_multispecies.R ALL
#   Rscript test_Ipos_multispecies.R Dmel Dmel_test_results
#
# Processing order:
#   1. Use each species' complete cdspos_for_alledit table to classify codons
#      as pos1, pos2, pos3, or multi.
#   2. Match those positions to A-site *_df_alledit files only.
#   3. In every biological replicate, apply RNA/RPF/alt-RNA thresholds first,
#      then calculate Lm, Lr, Ar/Am, and Ir/Im.
#   4. Average those ratios across biological replicates within species/sex,
#      then calculate Ipos=log2(mean Lr/mean Lm).
#   5. Plot one retained CDS codon per point in sex-separated violin/box and
#      I-position plots; scatter analysis is in table1_filter/run.sh.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
    getwd()
}
pipeline_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
project_root <- Sys.getenv("SOURCE_PROJECT", unset = pipeline_root)
args <- commandArgs(trailingOnly = TRUE)

requested_species <- if (length(args) >= 1L && nzchar(args[[1L]])) {
    trimws(strsplit(args[[1L]], ",", fixed = TRUE)[[1L]])
} else {
    c("Dmel", "Dsim", "Dsec", "Dere")
}
if (length(requested_species) == 1L && toupper(requested_species) == "ALL") {
    requested_species <- c("Dmel", "Dsim", "Dsec", "Dere")
}

valid_species <- c("Dmel", "Dsim", "Dsec", "Dere")
species_lookup <- setNames(valid_species, tolower(valid_species))
requested_species <- unname(species_lookup[tolower(requested_species)])
if (anyNA(requested_species)) {
    stop("Species must be one or more of: ", paste(valid_species, collapse = ", "))
}
requested_species <- unique(requested_species)

outdir <- if (length(args) >= 2L && nzchar(args[[2L]])) {
    args[[2L]]
} else {
    file.path(pipeline_root, "results", "05_RNA_Ribo_alledit_Table1", "04_multispecies_plots")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

alledit_dir <- file.path(
    project_root, "04_alledit", "result",
    "alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed_newdata",
    "results"
)
alledit_dir <- Sys.getenv("ALLEDIT_DIR", unset = alledit_dir)
cdspos_dir <- file.path(
    project_root, "03_cdspos", "result",
    "cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed"
)
cdspos_dir <- Sys.getenv("CDSPOS_DIR", unset = cdspos_dir)

cdspos_files <- c(
    Dmel = "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    Dsim = "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    Dsec = "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    Dere = "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
)

# Dmel includes ISO1 and WBg only. Adar5G1, AdarS, AdarG, and w1118 are
# excluded. Other species include their F/M replicates and head-ribo inputs.
sample_patterns <- c(
    Dmel = "^(Dmel_ISO1_head_(female|male)_ribo|WBg_[FM][0-9]+)$",
    Dsim = "^Dsim_(F[0-9]+|M[0-9]+|female_head_ribo|male_head_ribo)$",
    Dsec = "^Dsec_(F[0-9]+|M[0-9]+|male_head_ribo)$",
    Dere = "^Dere_(F[0-9]+|M[0-9]+|male_head_ribo)$"
)

min_rep_rna_total <- 3
min_rep_ribo_total <- 3
min_rep_alt_rna <- 1
eps <- 1e-6
group_levels <- c("pos1", "pos2", "pos3", "multi")

mean_no_na <- function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

median_no_na <- function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

sd_no_na <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else sd(x)
}

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

sum_no_na <- function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

first_no_na <- function(x) {
    z <- x[!is.na(x)]
    if (length(z) == 0L) NA_real_ else as.numeric(z[[1L]])
}

as_numeric_na <- function(x) {
    z <- trimws(as.character(x))
    z[tolower(z) %in% c("", "na", "nan", "null", ".")] <- NA_character_
    suppressWarnings(as.numeric(z))
}

read_cds_classification <- function(species) {
    file <- file.path(cdspos_dir, cdspos_files[[species]])
    if (!file.exists(file)) stop("Missing cdspos file: ", file)

    x <- fread(file, sep = "\t", na.strings = c("NA", "na", "NaN", ""))
    needed <- c("new_transcript", "new_cds_pos_1based")
    missing_cols <- setdiff(needed, names(x))
    if (length(missing_cols) > 0L) {
        stop(
            "Missing columns in ", basename(file), ": ",
            paste(missing_cols, collapse = ", ")
        )
    }

    x <- x[, .(
        transcript = as.character(new_transcript),
        cds_pos = as_numeric_na(new_cds_pos_1based)
    )]
    x <- unique(x[!is.na(transcript) & transcript != "" & !is.na(cds_pos)])
    x[, cds_pos := as.integer(cds_pos)]
    x <- x[cds_pos >= 1L]
    x[, codon_start := cds_pos - ((cds_pos - 1L) %% 3L)]
    x[, codon_pos := ((cds_pos - 1L) %% 3L) + 1L]
    x[, `:=`(
        FBt_cdspos = paste(transcript, cds_pos, sep = "_"),
        cds_handle = paste(transcript, codon_start, sep = "_")
    )]

    codons <- x[, .(
        n_total_edit_sites = uniqueN(FBt_cdspos),
        total_edit_positions = paste(sort(unique(codon_pos)), collapse = ","),
        cdspos_ids = paste(sort(unique(FBt_cdspos)), collapse = ",")
    ), by = cds_handle]
    codons[, group := ifelse(
        n_total_edit_sites > 1L,
        "multi",
        paste0("pos", total_edit_positions)
    )]

    out <- merge(x, codons, by = "cds_handle", all.x = TRUE, sort = FALSE)
    out[, species := species]
    out[]
}

discover_samples <- function(species) {
    files <- list.files(
        alledit_dir,
        pattern = "\\.Asite_df_alledit$",
        full.names = TRUE
    )
    sample_ids <- sub("\\.Asite_df_alledit$", "", basename(files))
    keep <- grepl(sample_patterns[[species]], sample_ids, perl = TRUE)
    data.table(
        species = species,
        sample_id = sample_ids[keep],
        file = files[keep],
        cohort = if (species == "Dmel") {
            ifelse(grepl("^WBg_", sample_ids[keep]), "WBg", "ISO1")
        } else {
            species
        }
    )[order(sample_id)]
}

read_one_sample <- function(file, species, cohort, sample_id, sex, cds_sites) {
    x <- fread(
        file,
        sep = "\t",
        na.strings = c("NA", "na", "NaN", "nan", "")
    )
    needed <- c(
        "handle", "ref_number_rna", "alt_number_rna",
        "ref_number_ribo", "alt_number_ribo", "FBt_cdspos"
    )
    missing_cols <- setdiff(needed, names(x))
    if (length(missing_cols) > 0L) {
        stop(
            "Missing columns in ", basename(file), ": ",
            paste(missing_cols, collapse = ", ")
        )
    }

    x <- x[, .(
        alledit_handle = as.character(handle),
        FBt_cdspos = as.character(FBt_cdspos),
        ref_rna = as_numeric_na(ref_number_rna),
        alt_rna = as_numeric_na(alt_number_rna),
        ref_ribo = as_numeric_na(ref_number_ribo),
        alt_ribo = as_numeric_na(alt_number_ribo)
    )]

    # Matching to the complete CDS-site list both restricts the analysis to
    # candidate editing sites and supplies a sample-independent multi label.
    x <- merge(
        x,
        cds_sites[, .(
            FBt_cdspos, cds_handle, group, cdspos_ids,
            n_total_edit_sites, total_edit_positions
        )],
        by = "FBt_cdspos",
        all = FALSE,
        sort = FALSE,
        allow.cartesian = FALSE
    )
    if (nrow(x) == 0L) return(NULL)

    n_handle_mismatch <- x[
        !is.na(alledit_handle) & alledit_handle != cds_handle, .N
    ]
    if (n_handle_mismatch > 0L) {
        warning(
            sample_id, ": ", n_handle_mismatch,
            " rows have an alledit/CDS codon-handle mismatch; CDS handle used"
        )
    }

    # There is one output row per codon per biological replicate. In a multi
    # codon, reference counts are shared and alternative counts are summed.
    ans <- x[, .(
        ref_rna = first_no_na(ref_rna),
        alt_rna = sum_no_na(alt_rna),
        ref_ribo = first_no_na(ref_ribo),
        alt_ribo = sum_no_na(alt_ribo),
        n_matched_edit_sites = uniqueN(FBt_cdspos)
    ), by = .(
        handle = cds_handle, group,
        cdspos_ids, n_total_edit_sites, total_edit_positions
    )]
    ans[, `:=`(
        species = species, cohort = cohort, sample_id = sample_id, sex = sex
    )]
    ans[, `:=`(
        tot_rna = ref_rna + alt_rna,
        tot_ribo = ref_ribo + alt_ribo
    )]
    ans[, valid_pair :=
        is.finite(tot_rna) & tot_rna >= min_rep_rna_total &
            is.finite(tot_ribo) & tot_ribo >= min_rep_ribo_total &
            is.finite(alt_rna) & alt_rna >= min_rep_alt_rna]
    # Ratios are calculated inside each valid replicate. Biological replicates
    # are averaged only later, so every valid replicate has equal weight.
    ans[, Lm := fifelse(valid_pair, alt_rna / tot_rna, NA_real_)]
    ans[, Lr := fifelse(valid_pair, alt_ribo / tot_ribo, NA_real_)]
    ans[, Ar_Am := fifelse(valid_pair & ref_rna > 0, ref_ribo / ref_rna, NA_real_)]
    ans[, Ir_Im := fifelse(valid_pair & alt_rna > 0, alt_ribo / alt_rna, NA_real_)]
    setcolorder(ans, c(
        "species", "cohort", "sample_id", "sex", "handle", "group", "cdspos_ids",
        "n_total_edit_sites", "total_edit_positions", "n_matched_edit_sites",
        "ref_rna", "alt_rna", "tot_rna", "Lm",
        "ref_ribo", "alt_ribo", "tot_ribo", "Lr",
        "Ar_Am", "Ir_Im", "valid_pair"
    ))
    ans[]
}

alledit_sample_info <- rbindlist(lapply(requested_species, discover_samples), fill = TRUE)
if (nrow(alledit_sample_info) == 0L) {
    stop("No requested A-site alledit files were found in: ", alledit_dir)
}
missing_species <- setdiff(requested_species, unique(alledit_sample_info$species))
if (length(missing_species) > 0L) {
    stop("No A-site alledit files found for: ", paste(missing_species, collapse = ", "))
}
alledit_sample_info[, sex := fifelse(
    grepl("female|(^|_)F[0-9]+$", sample_id, ignore.case = TRUE, perl = TRUE),
    "female",
    fifelse(
        grepl("male|(^|_)M[0-9]+$", sample_id, ignore.case = TRUE, perl = TRUE),
        "male",
        NA_character_
    )
)]
if (anyNA(alledit_sample_info$sex)) {
    stop(
        "Could not infer sex from sample name(s): ",
        paste(alledit_sample_info[is.na(sex), sample_id], collapse = ", ")
    )
}
fwrite(
    alledit_sample_info,
    file.path(outdir, "input_Asite_alledit_files.tsv"),
    sep = "\t"
)

cds_list <- setNames(lapply(requested_species, read_cds_classification), requested_species)
cds_codons <- rbindlist(lapply(requested_species, function(sp) {
    z <- unique(cds_list[[sp]][, .(
        species, cds_handle, group, cdspos_ids,
        n_total_edit_sites, total_edit_positions
    )])
    z
}))
setorder(cds_codons, species, group, cds_handle)
fwrite(
    cds_codons,
    file.path(outdir, "cds_codon_classification.tsv"),
    sep = "\t", na = "NA"
)
cds_summary <- cds_codons[, .(n_codons = .N), by = .(species, group)]
fwrite(cds_summary, file.path(outdir, "cds_group_summary.tsv"), sep = "\t")

# Final analysis path: every alledit biological replicate is filtered first;
# Lm, Lr, Ar/Am, and Ir/Im are calculated within that replicate. Those ratios
# are then averaged across replicates within species and sex. Ipos is finally
# calculated as log2(mean Lr / mean Lm). No pooled count ratio is used.
replicate_data <- rbindlist(lapply(seq_len(nrow(alledit_sample_info)), function(i) {
    info <- alledit_sample_info[i]
    sp <- info$species[[1L]]
    message("Reading alledit ", sp, ": ", info$sample_id[[1L]])
    read_one_sample(
        info$file[[1L]], sp, info$cohort[[1L]], info$sample_id[[1L]], info$sex[[1L]],
        cds_list[[sp]]
    )
}), fill = TRUE)
if (nrow(replicate_data) == 0L) stop("No alledit rows matched the CDS coordinates")
fwrite(
    replicate_data,
    file.path(outdir, "replicate_level_Lm_Lr.tsv"),
    sep = "\t", na = "NA"
)
fwrite(
    replicate_data,
    file.path(outdir, "alledit_replicate_ratios.tsv"),
    sep = "\t", na = "NA"
)

sex_sample_n <- alledit_sample_info[, .(
    n_input_samples = uniqueN(sample_id)
), by = .(species, cohort, sex)]
sex_site_means <- replicate_data[, .(
    mean_Lm = mean_no_na(Lm),
    mean_Lr = mean_no_na(Lr),
    Ar_Am = mean_no_na(Ar_Am),
    Ir_Im = mean_no_na(Ir_Im),
    n_Ipos_replicates = sum(is.finite(Lm) & is.finite(Lr)),
    n_Lm_Lr_replicates = sum(is.finite(Lm) & is.finite(Lr)),
    n_ArAm_IrIm_replicates = sum(is.finite(Ar_Am) & is.finite(Ir_Im)),
    replicates_used = paste(
        sort(unique(sample_id[is.finite(Lm) & is.finite(Lr)])),
        collapse = ","
    )
), by = .(
    species, cohort, sex, handle, group, cdspos_ids,
    n_total_edit_sites, total_edit_positions
)]
sex_site_means[, lr_lm := log2((mean_Lr + eps) / (mean_Lm + eps))]
sex_site_means[, lr_minus_lm := mean_Lr - mean_Lm]
sex_site_means <- merge(
    sex_site_means, sex_sample_n,
    by = c("species", "cohort", "sex"), all.x = TRUE
)
sex_site_means[, pass_filter :=
    is.finite(lr_lm) & is.finite(mean_Lm) & is.finite(mean_Lr) &
        is.finite(Ar_Am) & is.finite(Ir_Im)]
sex_site_means[, group := factor(group, levels = group_levels)]
setorder(sex_site_means, species, sex, group, handle)
fwrite(
    sex_site_means,
    file.path(outdir, "species_sex_site_means_all.tsv"),
    sep = "\t", na = "NA"
)
sex_plot_dt <- sex_site_means[pass_filter == TRUE & !is.na(group)]
if (nrow(sex_plot_dt) == 0L) {
    stop("No species-sex codons passed all Ipos/scatter filters")
}
fwrite(
    sex_plot_dt,
    file.path(outdir, "species_sex_site_means_filtered.tsv"),
    sep = "\t", na = "NA"
)

# Species-only provenance summary. The plotted site distributions remain
# separated by sex and use sex_plot_dt directly.
plot_dt <- sex_plot_dt[, .(
    lr_lm = mean_no_na(lr_lm),
    n_total_edit_sites = max(n_total_edit_sites),
    cdspos_ids = first(cdspos_ids),
    total_edit_positions = first(total_edit_positions)
), by = .(species, cohort, handle, group)]
plot_dt[, group := factor(group, levels = group_levels)]
setorder(plot_dt, species, group, handle)
fwrite(
    plot_dt,
    file.path(outdir, "species_site_means_all.tsv"),
    sep = "\t", na = "NA"
)
fwrite(
    plot_dt,
    file.path(outdir, "species_site_means_filtered.tsv"),
    sep = "\t", na = "NA"
)

sample_info <- copy(alledit_sample_info)

species_order <- valid_species[valid_species %in% requested_species]
plot_dt[, species := factor(species, levels = species_order)]
sex_plot_dt[, species := factor(species, levels = species_order)]

group_cols <- c(
    pos1 = "#4E9AD3", pos2 = "#F5A623",
    pos3 = "#7ED321", multi = "#D0021B"
)
species_cols <- c(
    Dmel = "#4E79A7", Dsim = "#F28E2B",
    Dsec = "#59A14F", Dere = "#E15759"
)
adar_sites_by_species <- c(
    Dmel = "FBtr0339272_1288",
    Dsim = "FBtr0358708_1399",
    Dsec = "FBtr0202151_1165",
    Dere = "FBtr0413873_1519"
)
combined_y_limits <- c(-2, 2)
combined_y_n <- -1.75
combined_y_sig <- 1.75

# Figure 1: one retained CDS codon per point. Female and male are drawn in
# separate figures. The all-sex panel uses the union of female and male codons.
# When both sexes have valid values, mean_Lr and mean_Lm are averaged equally
# across sex before subtraction; a codon observed in only one sex retains that
# sex's value. A multi-edit point can contain multiple cdspos.

all_sex_plot_dt <- sex_plot_dt[, .(
    mean_Lm = mean_no_na(mean_Lm),
    mean_Lr = mean_no_na(mean_Lr),
    n_sexes = uniqueN(as.character(sex)),
    sexes_used = paste(sort(unique(as.character(sex))), collapse = ","),
    n_total_edit_sites = max(n_total_edit_sites),
    cdspos_ids = first(cdspos_ids),
    total_edit_positions = first(total_edit_positions)
), by = .(species, cohort, handle, group)]
all_sex_plot_dt <- all_sex_plot_dt[
    n_sexes >= 1L & is.finite(mean_Lm) & is.finite(mean_Lr)
]
all_sex_plot_dt[, lr_minus_lm := mean_Lr - mean_Lm]
setorder(all_sex_plot_dt, species, group, handle)
all_sex_plot_dt_out <- copy(all_sex_plot_dt)
all_sex_plot_dt_out[, `:=`(
    species = as.character(species),
    group = as.character(group)
)]
fwrite(
    all_sex_plot_dt_out,
    file.path(outdir, "species_site_sex_mean_Lr_Lm.tsv"),
    sep = "\t", na = "NA"
)

# Tie-corrected Dunn test for the pooled site distributions. The BH correction
# is applied across the six pairwise editing-group comparisons within a panel.
site_dunn_bh <- function(dt) {
    d <- copy(dt)
    d[, group := as.character(group)]
    d <- d[is.finite(lr_minus_lm) & group %in% group_levels]
    groups <- group_levels[group_levels %in% unique(d$group)]
    if (length(groups) < 2L) {
        return(data.table(
            group1 = character(), group2 = character(), statistic = numeric(),
            p_value = numeric(), p_adj = numeric(), method = character()
        ))
    }

    d[, rank_value := rank(lr_minus_lm, ties.method = "average")]
    n_total <- nrow(d)
    tie_sizes <- as.numeric(table(d$lr_minus_lm))
    tie_correction <- if (n_total > 1L) {
        1 - sum(tie_sizes^3 - tie_sizes) / (n_total^3 - n_total)
    } else {
        1
    }
    rank_summary <- d[, .(n = .N, mean_rank = mean(rank_value)), by = group]
    ans <- rbindlist(lapply(combn(groups, 2L, simplify = FALSE), function(pair) {
        g1 <- rank_summary[group == pair[[1L]]]
        g2 <- rank_summary[group == pair[[2L]]]
        denominator <- sqrt(
            (n_total * (n_total + 1) / 12) *
                (1 / g1$n + 1 / g2$n) * tie_correction
        )
        z_value <- (g1$mean_rank - g2$mean_rank) / denominator
        data.table(
            group1 = pair[[1L]], group2 = pair[[2L]],
            statistic = z_value,
            p_value = 2 * pnorm(-abs(z_value))
        )
    }))
    ans[, p_adj := p.adjust(p_value, method = "BH")]
    ans[, method := "Dunn test; Benjamini-Hochberg adjusted"]
    ans[]
}

# Compact-letter display: groups sharing at least one letter are not
# significantly different at the selected alpha level.
site_compact_letters <- function(groups, pairwise, alpha = 0.05) {
    groups <- as.character(groups)
    n_groups <- length(groups)
    if (n_groups == 0L) {
        return(data.table(group = character(), cld = character()))
    }
    significant <- matrix(
        FALSE,
        nrow = n_groups, ncol = n_groups,
        dimnames = list(groups, groups)
    )
    if (nrow(pairwise) > 0L) {
        for (i in seq_len(nrow(pairwise))) {
            g1 <- as.character(pairwise$group1[[i]])
            g2 <- as.character(pairwise$group2[[i]])
            if (g1 %in% groups && g2 %in% groups &&
                is.finite(pairwise$p_adj[[i]]) && pairwise$p_adj[[i]] < alpha) {
                significant[g1, g2] <- TRUE
                significant[g2, g1] <- TRUE
            }
        }
    }
    candidate_cliques <- list()
    if (n_groups >= 2L) {
        for (size in 2:n_groups) {
            candidate_cliques <- c(
                candidate_cliques,
                combn(groups, size, simplify = FALSE)
            )
        }
    }
    candidate_cliques <- Filter(function(clique) {
        idx <- match(clique, groups)
        all(!significant[idx, idx, drop = FALSE])
    }, candidate_cliques)
    if (length(candidate_cliques) > 0L) {
        is_subsumed <- vapply(seq_along(candidate_cliques), function(i) {
            any(vapply(seq_along(candidate_cliques), function(j) {
                j != i && length(candidate_cliques[[j]]) > length(candidate_cliques[[i]]) &&
                    all(candidate_cliques[[i]] %in% candidate_cliques[[j]])
            }, logical(1)))
        }, logical(1))
        candidate_cliques <- candidate_cliques[!is_subsumed]
    }
    covered <- unique(unlist(candidate_cliques, use.names = FALSE))
    uncovered <- setdiff(groups, covered)
    if (length(uncovered) > 0L) {
        candidate_cliques <- c(candidate_cliques, as.list(uncovered))
    }
    if (length(candidate_cliques) == 0L) candidate_cliques <- as.list(groups)
    label_letters <- letters[seq_along(candidate_cliques)]
    data.table(
        group = groups,
        cld = vapply(groups, function(g) {
            paste(
                label_letters[vapply(
                    candidate_cliques, function(clique) g %in% clique, logical(1)
                )],
                collapse = ""
            )
        }, character(1))
    )
}

make_site_violin <- function(d, plot_title, panel_name) {
    d <- copy(d)
    if (nrow(d) == 0L) stop("No filtered sites for panel: ", panel_name)
    d[, group := factor(as.character(group), levels = group_levels)]

    site_stats <- d[, .(
        n_points = .N,
        n_edit_sites = sum(n_total_edit_sites),
        mean_lr_minus_lm = mean(lr_minus_lm),
        sd_lr_minus_lm = sd_no_na(lr_minus_lm)
    ), by = group]
    y_min <- min(d$lr_minus_lm, na.rm = TRUE)
    y_max <- max(d$lr_minus_lm, na.rm = TRUE)
    y_range <- y_max - y_min
    if (!is.finite(y_range) || y_range <= 0) y_range <- 1

    kw_test <- tryCatch(
        kruskal.test(lr_minus_lm ~ group, data = d),
        error = function(e) NULL
    )
    kw_out <- data.table(
        panel = panel_name,
        test = "Kruskal-Wallis",
        statistic = if (is.null(kw_test)) NA_real_ else unname(as.numeric(kw_test$statistic)),
        df = if (is.null(kw_test)) NA_real_ else unname(as.numeric(kw_test$parameter)),
        p_value = if (is.null(kw_test)) NA_real_ else kw_test$p.value,
        n_points = nrow(d)
    )
    dunn_out <- site_dunn_bh(d)
    dunn_out[, panel := panel_name]
    setcolorder(dunn_out, c(
        "panel", "group1", "group2", "statistic", "p_value", "p_adj", "method"
    ))
    letters_dt <- site_compact_letters(
        group_levels[group_levels %in% as.character(site_stats$group)],
        dunn_out
    )
    # Letter names are arbitrary statistically; remap them so the group with
    # the largest mean starts at "a", using the reference figure style.
    mean_order <- site_stats[
        order(-mean_lr_minus_lm), as.character(group)
    ]
    ordered_cld <- letters_dt[
        match(mean_order, as.character(group)), cld
    ]
    old_letters <- unique(unlist(strsplit(ordered_cld, "", fixed = TRUE)))
    letter_map <- setNames(letters[seq_along(old_letters)], old_letters)
    letters_dt[, cld := vapply(
        strsplit(cld, "", fixed = TRUE),
        function(x) paste(letter_map[x], collapse = ""),
        character(1)
    )]
    letters_dt[, `:=`(
        panel = panel_name,
        group = factor(group, levels = group_levels),
        y_letter = y_max + 0.22 * y_range
    )]
    kw_label <- if (is.null(kw_test)) {
        "Kruskal-Wallis P = NA"
    } else if (kw_test$p.value < 0.001) {
        "Kruskal-Wallis P < 0.001"
    } else {
        paste0("Kruskal-Wallis P = ", formatC(kw_test$p.value, format = "f", digits = 3L))
    }
    site_stats[, `:=`(
        y_n = y_max + 0.045 * y_range,
        y_mean_sd = y_max + 0.115 * y_range,
        n_label = paste0("italic(n)==", n_points),
        mean_sd_label = paste0(
            format_stat(mean_lr_minus_lm, digits = 2L),
            " <U+00B1> ", format_stat(sd_lr_minus_lm, digits = 2L)
        )
    )]
    y_bracket <- y_min - 0.16 * y_range
    y_bracket_tick <- y_min - 0.13 * y_range
    y_class <- y_min - 0.23 * y_range

    set.seed(18907675L)
    p <- ggplot(d, aes(x = group, y = lr_minus_lm, fill = group)) +
        geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, color = "grey65") +
        geom_violin(
            scale = "width", trim = FALSE, alpha = 0.38,
            color = "grey30", linewidth = 0.35
        ) +
        geom_boxplot(
            width = 0.18, fill = "white", color = "black",
            outlier.shape = NA, linewidth = 0.35
        ) +
        geom_jitter(
            aes(color = group),
            width = 0.08, height = 0,
            size = 0.75, alpha = 0.28, show.legend = FALSE
        ) +
        geom_text(
            data = site_stats,
            aes(x = group, y = y_n, label = n_label),
            inherit.aes = FALSE, vjust = 0, size = 3, parse = TRUE
        ) +
        geom_text(
            data = site_stats,
            aes(x = group, y = y_mean_sd, label = mean_sd_label),
            inherit.aes = FALSE, vjust = 0, size = 2.55
        ) +
        geom_text(
            data = letters_dt,
            aes(x = group, y = y_letter, label = cld),
            inherit.aes = FALSE, color = "#2447A8", fontface = "bold", size = 4
        ) +
        annotate("segment", x = 0.70, xend = 3.30, y = y_bracket, yend = y_bracket) +
        annotate("segment", x = 0.70, xend = 0.70, y = y_bracket, yend = y_bracket_tick) +
        annotate("segment", x = 3.30, xend = 3.30, y = y_bracket, yend = y_bracket_tick) +
        annotate("text", x = 2, y = y_class, label = "Single-editing", vjust = 1) +
        annotate("text", x = 4, y = y_class, label = "Multi-editing", vjust = 1) +
        scale_fill_manual(values = group_cols, drop = FALSE) +
        scale_color_manual(values = group_cols, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.26, 0.25))) +
        coord_cartesian(clip = "off") +
        labs(
            title = plot_title,
            x = "Editing-position group",
            y = expression(mean ~ italic(Lr) - mean ~ italic(Lm)),
            caption = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
            legend.position = "none",
            strip.background = element_rect(fill = "grey95", color = "grey70"),
            strip.text = element_text(face = "plain"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.caption = element_text(hjust = 0, size = 8.5),
            plot.margin = margin(8, 8, 34, 8)
        )

    list(
        plot = p, stats = site_stats,
        kruskal_wallis = kw_out, dunn = dunn_out, letters = letters_dt
    )
}

site_plots <- list(
    female = make_site_violin(
        sex_plot_dt[as.character(sex) == "female"], "Female", "female"
    ),
    male = make_site_violin(
        sex_plot_dt[as.character(sex) == "male"], "Male", "male"
    ),
    all = make_site_violin(
        all_sex_plot_dt, "All (female/male union)", "all"
    )
)
site_stats_out <- rbindlist(lapply(names(site_plots), function(panel_name) {
    z <- copy(site_plots[[panel_name]]$stats)
    z[, panel := panel_name]
    z[, .(
        panel,
        group = as.character(group), n_points, n_edit_sites,
        mean_lr_minus_lm, sd_lr_minus_lm
    )]
}))
setorder(site_stats_out, panel, group)
fwrite(
    site_stats_out,
    file.path(outdir, "site_violin_group_n_mean_sd.tsv"),
    sep = "\t", na = "NA"
)

site_kw_out <- rbindlist(lapply(site_plots, `[[`, "kruskal_wallis"))
site_dunn_out <- rbindlist(lapply(site_plots, `[[`, "dunn"))
site_letters_out <- rbindlist(lapply(site_plots, function(x) {
    z <- copy(x$letters)
    z[, .(panel, group = as.character(group), cld)]
}))
fwrite(
    site_kw_out,
    file.path(outdir, "site_violin_kruskal_wallis.tsv"),
    sep = "\t", na = "NA"
)
fwrite(
    site_dunn_out,
    file.path(outdir, "site_violin_dunn_bh.tsv"),
    sep = "\t", na = "NA"
)
fwrite(
    site_letters_out,
    file.path(outdir, "site_violin_compact_letters.tsv"),
    sep = "\t", na = "NA"
)

site_width <- 7.2
site_height <- 5.5
for (panel_name in names(site_plots)) {
    ggsave(
        file.path(outdir, paste0("site_violin_box_by_species_", panel_name, ".pdf")),
        site_plots[[panel_name]]$plot,
        width = site_width, height = site_height
    )
}

# Keep the established filename as a three-page PDF: female, male, all-sex.
pdf(
    file.path(outdir, "site_violin_box_by_species.pdf"),
    width = site_width, height = site_height, onefile = TRUE
)
print(site_plots$female$plot)
print(site_plots$male$plot)
print(site_plots$all$plot)
dev.off()

# Save the previous all-sample species medians as a separate provenance table.
all_group_medians <- plot_dt[, .(
    median_lr_lm = median_no_na(lr_lm),
    n_sites = sum(!is.na(lr_lm))
), by = .(species, group)]
all_group_medians[, species := factor(species, levels = species_order)]
all_group_medians[, group := factor(group, levels = group_levels)]
setorder(all_group_medians, species, group)
all_group_medians_out <- copy(all_group_medians)
all_group_medians_out[, `:=`(
    species = as.character(species),
    group = as.character(group)
)]
fwrite(
    all_group_medians_out,
    file.path(outdir, "species_group_medians_all_samples.tsv"),
    sep = "\t", na = "NA"
)

# Figure 2: one site median per species-sex unit for each editing group.
group_medians <- sex_plot_dt[, .(
    median_lr_lm = median_no_na(lr_lm),
    n_sites = sum(!is.na(lr_lm))
), by = .(species, sex, group)]
group_medians[, species := factor(species, levels = species_order)]
group_medians[, sex := factor(sex, levels = c("female", "male"))]
group_medians[, group := factor(group, levels = group_levels)]
group_medians[, unit_id := paste(as.character(species), as.character(sex), sep = "__")]
setorder(group_medians, species, sex, group)
median_points <- copy(group_medians)
median_points[, median_lr_lm_display := median_lr_lm]
median_points[
    as.character(group) == "pos3" & median_lr_lm < combined_y_limits[[1L]],
    median_lr_lm_display := combined_y_limits[[1L]]
]
group_medians_out <- copy(group_medians)
group_medians_out[, `:=`(
    species = as.character(species),
    sex = as.character(sex),
    group = as.character(group)
)]
group_medians_out[, plot_median_lr_lm := median_lr_lm]
group_medians_out[
    group == "pos3" & median_lr_lm < combined_y_limits[[1L]],
    plot_median_lr_lm := combined_y_limits[[1L]]
]
fwrite(
    group_medians_out,
    file.path(outdir, "species_group_medians.tsv"),
    sep = "\t", na = "NA"
)

# Requested unpaired analysis: Kruskal-Wallis followed by a tie-corrected
# Dunn test. No optional post-hoc package is required for reproducibility.
dunn_bh <- function(dt) {
    d <- copy(dt)
    d[, group := as.character(group)]
    d <- d[is.finite(median_lr_lm) & !is.na(group)]
    groups <- group_levels[group_levels %in% unique(d$group)]
    if (length(groups) < 2L) {
        return(data.table(
            group1 = character(), group2 = character(), statistic = numeric(),
            p_value = numeric(), p_adj = numeric(), method = character()
        ))
    }

    d[, rank_value := rank(median_lr_lm, ties.method = "average")]
    n_total <- nrow(d)
    tie_sizes <- as.numeric(table(d$median_lr_lm))
    tie_correction <- if (n_total > 1L) {
        1 - sum(tie_sizes^3 - tie_sizes) / (n_total^3 - n_total)
    } else {
        1
    }
    rank_summary <- d[, .(n = .N, mean_rank = mean(rank_value)), by = group]
    pair_list <- combn(groups, 2L, simplify = FALSE)
    ans <- rbindlist(lapply(pair_list, function(pair) {
        g1 <- rank_summary[group == pair[[1L]]]
        g2 <- rank_summary[group == pair[[2L]]]
        denominator <- sqrt(
            (n_total * (n_total + 1) / 12) *
                (1 / g1$n + 1 / g2$n) * tie_correction
        )
        z_value <- (g1$mean_rank - g2$mean_rank) / denominator
        data.table(
            group1 = pair[[1L]], group2 = pair[[2L]],
            statistic = z_value,
            p_value = 2 * pnorm(-abs(z_value))
        )
    }))
    ans[, p_adj := p.adjust(p_value, method = "BH")]
    ans[, method := "Dunn"]
    setcolorder(ans, c("group1", "group2", "statistic", "p_value", "p_adj", "method"))
    ans[]
}

combined_for_stats <- copy(group_medians_out)
combined_for_stats[, group := factor(group, levels = group_levels)]
combined_for_stats <- combined_for_stats[
    is.finite(median_lr_lm) & !is.na(group)
]
kw_test <- tryCatch(
    kruskal.test(median_lr_lm ~ group, data = combined_for_stats),
    error = function(e) NULL
)
kw_out <- if (is.null(kw_test)) {
    data.table(
        test = "Kruskal-Wallis", statistic = NA_real_, df = NA_real_,
        p_value = NA_real_, n_points = nrow(combined_for_stats)
    )
} else {
    data.table(
        test = "Kruskal-Wallis",
        statistic = unname(as.numeric(kw_test$statistic)),
        df = unname(as.numeric(kw_test$parameter)),
        p_value = kw_test$p.value,
        n_points = nrow(combined_for_stats)
    )
}
dunn_out <- dunn_bh(combined_for_stats)
fwrite(kw_out, file.path(outdir, "combined_median_kruskal_wallis.tsv"), sep = "\t", na = "NA")
fwrite(dunn_out, file.path(outdir, "combined_median_dunn_bh.tsv"), sep = "\t", na = "NA")

# Primary analysis for the plot: each species-sex unit contributes one value
# to every group, so use a Friedman omnibus test and paired Wilcoxon tests.
paired_wilcoxon_bh <- function(dt) {
    d <- copy(dt)
    d[, group := as.character(group)]
    d[, unit_id := as.character(unit_id)]
    groups <- group_levels[group_levels %in% unique(d$group)]
    wide <- dcast(d, unit_id ~ group, value.var = "median_lr_lm")
    groups <- groups[groups %in% names(wide)]
    if (length(groups) < 2L) {
        return(data.table(
            group1 = character(), group2 = character(), statistic = numeric(),
            p_value = numeric(), p_adj = numeric(), n_pairs = integer(),
            method = character()
        ))
    }
    pair_list <- combn(groups, 2L, simplify = FALSE)
    ans <- rbindlist(lapply(pair_list, function(pair) {
        x <- wide[[pair[[1L]]]]
        y <- wide[[pair[[2L]]]]
        keep <- is.finite(x) & is.finite(y)
        x <- x[keep]
        y <- y[keep]
        test <- if (length(x) >= 2L) {
            tryCatch(
                wilcox.test(x, y, paired = TRUE, exact = FALSE),
                error = function(e) NULL
            )
        } else {
            NULL
        }
        data.table(
            group1 = pair[[1L]], group2 = pair[[2L]],
            statistic = if (is.null(test)) NA_real_ else unname(as.numeric(test$statistic)),
            p_value = if (is.null(test)) NA_real_ else test$p.value,
            n_pairs = length(x)
        )
    }))
    ans[, p_adj := p.adjust(p_value, method = "BH")]
    ans[, method := "Paired Wilcoxon signed-rank"]
    setcolorder(ans, c("group1", "group2", "statistic", "p_value", "p_adj", "n_pairs", "method"))
    ans[]
}

friedman_test <- tryCatch(
    {
        complete_units <- combined_for_stats[, .(n_groups = uniqueN(group)), by = unit_id][
            n_groups == length(group_levels), unit_id
        ]
        friedman_dt <- combined_for_stats[unit_id %in% complete_units]
        if (length(complete_units) < 2L || uniqueN(friedman_dt$group) < 2L) {
            NULL
        } else {
            friedman_dt[, group_factor := factor(group, levels = group_levels)]
            friedman_dt[, unit_factor := factor(unit_id)]
            friedman.test(
                median_lr_lm ~ group_factor | unit_factor,
                data = friedman_dt
            )
        }
    },
    error = function(e) NULL
)
friedman_out <- if (is.null(friedman_test)) {
    data.table(
        test = "Friedman", statistic = NA_real_, df = NA_real_,
        p_value = NA_real_, n_units = NA_integer_
    )
} else {
    data.table(
        test = "Friedman",
        statistic = unname(as.numeric(friedman_test$statistic)),
        df = unname(as.numeric(friedman_test$parameter)),
        p_value = friedman_test$p.value,
        n_units = uniqueN(combined_for_stats$unit_id)
    )
}
paired_out <- paired_wilcoxon_bh(combined_for_stats)
fwrite(friedman_out, file.path(outdir, "combined_median_friedman.tsv"), sep = "\t", na = "NA")
fwrite(paired_out, file.path(outdir, "combined_median_paired_wilcoxon_bh.tsv"), sep = "\t", na = "NA")

# Convert the BH-adjusted paired comparisons into compact letters for the
# figure. Groups sharing a letter are not significantly different at alpha=.05.
compact_letters <- function(groups, pairwise, alpha = 0.05) {
    groups <- as.character(groups)
    n_groups <- length(groups)
    if (n_groups == 0L) return(data.table(group = character(), cld = character()))
    significant <- matrix(FALSE,
        nrow = n_groups, ncol = n_groups,
        dimnames = list(groups, groups)
    )
    if (nrow(pairwise) > 0L) {
        for (i in seq_len(nrow(pairwise))) {
            g1 <- as.character(pairwise$group1[[i]])
            g2 <- as.character(pairwise$group2[[i]])
            if (g1 %in% groups && g2 %in% groups &&
                !is.na(pairwise$p_adj[[i]]) && pairwise$p_adj[[i]] < alpha) {
                significant[g1, g2] <- TRUE
                significant[g2, g1] <- TRUE
            }
        }
    }
    candidate_cliques <- list()
    if (n_groups >= 2L) {
        for (size in 2:n_groups) {
            candidate_cliques <- c(
                candidate_cliques,
                combn(groups, size, simplify = FALSE)
            )
        }
    }
    candidate_cliques <- Filter(function(clique) {
        idx <- match(clique, groups)
        all(!significant[idx, idx, drop = FALSE])
    }, candidate_cliques)
    if (length(candidate_cliques) > 0L) {
        is_subsumed <- vapply(seq_along(candidate_cliques), function(i) {
            any(vapply(seq_along(candidate_cliques), function(j) {
                j != i && length(candidate_cliques[[j]]) > length(candidate_cliques[[i]]) &&
                    all(candidate_cliques[[i]] %in% candidate_cliques[[j]])
            }, logical(1)))
        }, logical(1))
        candidate_cliques <- candidate_cliques[!is_subsumed]
    }
    covered <- unique(unlist(candidate_cliques, use.names = FALSE))
    uncovered <- setdiff(groups, covered)
    if (length(uncovered) > 0L) {
        candidate_cliques <- c(candidate_cliques, as.list(uncovered))
    }
    if (length(candidate_cliques) == 0L) {
        candidate_cliques <- as.list(groups)
    }
    label_letters <- letters[seq_along(candidate_cliques)]
    data.table(
        group = groups,
        cld = vapply(groups, function(g) {
            paste(label_letters[vapply(candidate_cliques, function(clique) g %in% clique, logical(1))], collapse = "")
        }, character(1))
    )
}

stats_groups <- group_levels[group_levels %in% unique(as.character(group_medians$group))]
annotation_pairwise <- paired_out
if (nrow(annotation_pairwise) == 0L || all(is.na(annotation_pairwise$p_adj))) {
    annotation_pairwise <- dunn_out
}
letters_dt <- compact_letters(stats_groups, annotation_pairwise)
letters_dt[, group := factor(group, levels = group_levels)]

median_y <- group_medians[, .(
    y_min = min(median_lr_lm, na.rm = TRUE),
    y_max = max(median_lr_lm, na.rm = TRUE)
), by = group]
median_range <- max(median_y$y_max, na.rm = TRUE) - min(median_y$y_min, na.rm = TRUE)
if (!is.finite(median_range) || median_range <= 0) median_range <- 1
median_y[, `:=`(
    y_sig = combined_y_sig,
    y_n = -1.69,
    y_mean_sd = -1.90
)]
letters_dt <- merge(letters_dt, median_y[, .(group, y_sig)], by = "group", all.x = TRUE)
n_medians <- group_medians[, .(
    n_points = sum(is.finite(median_lr_lm)),
    mean_median_lr_lm = mean(median_lr_lm[is.finite(median_lr_lm)]),
    sd_median_lr_lm = sd_no_na(median_lr_lm)
), by = group]
n_medians <- merge(
    n_medians,
    median_y[, .(group, y_n, y_mean_sd)],
    by = "group", all.x = TRUE
)
n_medians[, `:=`(
    n_label = paste0("italic(n)==", n_points),
    mean_sd_label = paste0(
        format_stat(mean_median_lr_lm, digits = 2L),
        " <U+00B1> ", format_stat(sd_median_lr_lm, digits = 2L)
    )
)]
combined_n_out <- group_medians[, .(n_points = uniqueN(unit_id)), by = group]
fwrite(
    combined_n_out,
    file.path(outdir, "combined_median_n.tsv"),
    sep = "\t", na = "NA"
)

friedman_caption <- if (is.null(friedman_test)) {
    paste0(
        "Friedman test unavailable; letters show post-hoc groups\n",
        "pos3 values < -2 are displayed at -2; tests use original values; ",
        "values below n are mean <U+00B1> SD"
    )
} else {
    paste0(
        "Letters: paired Wilcoxon signed-rank, BH-adjusted\n",
        "pos3 values < -2 are displayed at -2; tests use original values; ",
        "values below n are mean <U+00B1> SD"
    )
}
friedman_p_label <- if (is.null(friedman_test)) {
    "Friedman~italic(P)==NA"
} else if (friedman_test$p.value < 0.001) {
    "Friedman~italic(P)<0.001"
} else {
    paste0(
        "Friedman~italic(P)==",
        formatC(friedman_test$p.value, format = "f", digits = 3L)
    )
}

# Figure 2 is intentionally a boxplot rather than a violin plot. Each point
# is a species-sex median; circles and triangles distinguish female and male.
set.seed(18907675L)
p_medians <- ggplot(group_medians, aes(x = group, y = median_lr_lm)) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, color = "grey65") +
    annotate(
        "text",
        x = "pos1", y = 1.98,
        label = friedman_p_label, parse = TRUE,
        hjust = 0, vjust = 1, size = 3.4
    ) +
    geom_boxplot(
        aes(group = group),
        width = 0.48, fill = "white", color = "grey25",
        outlier.shape = NA, linewidth = 0.45
    ) +
    geom_point(
        data = median_points,
        aes(x = group, y = median_lr_lm_display, color = species, shape = sex),
        inherit.aes = FALSE,
        position = position_jitter(width = 0.07, height = 0),
        size = 3, alpha = 0.9
    ) +
    geom_text(
        data = letters_dt,
        aes(x = group, y = y_sig, label = cld),
        inherit.aes = FALSE, fontface = "bold", size = 4
    ) +
    geom_text(
        data = n_medians,
        aes(x = group, y = y_n, label = n_label),
        inherit.aes = FALSE, vjust = 1, size = 3.2, parse = TRUE
    ) +
    geom_text(
        data = n_medians,
        aes(x = group, y = y_mean_sd, label = mean_sd_label),
        inherit.aes = FALSE, vjust = 1, size = 2.8
    ) +
    scale_color_manual(
        values = species_cols, drop = FALSE,
        labels = function(x) parse(text = paste0("italic(", x, ")"))
    ) +
    scale_shape_manual(values = c(female = 16, male = 17), drop = FALSE) +
    scale_y_continuous(
        breaks = seq(combined_y_limits[[1L]], combined_y_limits[[2L]], by = 1),
        expand = expansion(mult = 0)
    ) +
    coord_cartesian(ylim = combined_y_limits, clip = "off") +
    labs(
        x = "Editing-position group",
        y = expression(Median ~ log[2](mean[rep] ~ italic(Lr) / mean[rep] ~ italic(Lm))),
        color = "Species",
        shape = "Sex",
        caption = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.05),
        plot.margin = margin(8, 8, 38, 8)
    )

ggsave(
    file.path(outdir, "group_median_boxplot_all_species.pdf"),
    p_medians,
    width = 7.5, height = 6.2
)

message("I-position analysis complete: ", normalizePath(outdir))
