#!/usr/bin/env Rscript
# Revalidate target sites and derive the CDS editing groups used in Figure S26.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 02_groups.R candidates.tsv targets.tsv CDS_map.tsv output_dir")
candidate_file <- args[1]
target_site_file <- args[2]
site_gene_map_file <- args[3]
output_dir <- args[4]
scope <- "CDS"
min_depth <- 10
min_reps_per_genotype <- 2L
editing_change_cutoff <- 0.05
genotype_levels <- c("AG", "WBG", "AS")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- function(x, required, label) {
    missing <- setdiff(required, names(x))
    if (length(missing)) {
        stop(label, " is missing columns: ", paste(missing, collapse = ", "))
    }
}

for (f in c(candidate_file, target_site_file, site_gene_map_file)) {
    if (!file.exists(f))
        stop("Input file does not exist: ", f)
}

# -----------------------------------------------------------------------------
# 2. Replicate-level QC
# -----------------------------------------------------------------------------

candidate <- fread(candidate_file, select = c("coord_handle", "genotype", "sex", "sample",
    "refcount", "altcount", "dp"), check.names = FALSE)

required_columns(candidate, c("coord_handle", "genotype", "sex", "sample", "refcount", "altcount",
    "dp"), "candidate editing table")

candidate[, `:=`(coord_handle = trimws(as.character(coord_handle)), genotype = trimws(as.character(genotype)),
    sex = trimws(as.character(sex)), sample = trimws(as.character(sample)), refcount = suppressWarnings(as.numeric(refcount)),
    altcount = suppressWarnings(as.numeric(altcount)), dp = suppressWarnings(as.numeric(dp)))]

# The standardized table uses F1/F2/F3 or M1/M2/M3. The fallback also accepts
# full library names such as AG_F1.
candidate[, sample_short := fifelse(grepl("^[FM][0-9]+$", sample), sample, sub("^.*_([FM][0-9]+)$",
    "\\1", sample, perl = TRUE))]
candidate[!grepl("^[FM][0-9]+$", sample_short), sample_short := NA_character_]
candidate[, rep := sub("^[FM]", "", sample_short)]

candidate <- candidate[!is.na(coord_handle) & coord_handle != "" & genotype %in% genotype_levels &
    sex %in% c("F", "M") & !is.na(sample_short) & !is.na(rep) & rep != "" & is.finite(refcount) &
    is.finite(altcount) & is.finite(dp) & dp >= min_depth & refcount >= 0 & altcount >= 0 &
    altcount <= dp & (refcount + altcount) > 1]

candidate[, Lm := altcount/dp]

# Multiple rows can represent different annotations of the same quantitative
# site/library measurement. Collapse them before counting biological replicates.
site_rep <- candidate[, .(n_annotation_rows = .N, n_unique_measurements = uniqueN(paste(refcount,
    altcount, dp, sep = "|")), refcount = first(refcount), altcount = first(altcount), dp = first(dp),
    Lm = first(Lm)), by = .(coord_handle, sex, genotype, sample_short, rep)]

if (site_rep[, any(n_unique_measurements > 1L)]) {
    stop("Conflicting counts exist within at least one site/library measurement")
}

fwrite(site_rep, file.path(output_dir, "01_site_replicate_measurements_after_QC.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. Require >=2 valid replicates in AG, WBG and AS separately
# -----------------------------------------------------------------------------

replicate_summary <- site_rep[, .(n_valid_replicates = uniqueN(rep), mean_Lm = mean(Lm)),
    by = .(coord_handle, sex, genotype)]

nrep_wide <- dcast(replicate_summary, coord_handle + sex ~ genotype, value.var = "n_valid_replicates",
    fill = 0)

mean_wide <- dcast(replicate_summary, coord_handle + sex ~ genotype, value.var = "mean_Lm",
    fill = NA_real_)

for (g in genotype_levels) {
    if (!g %in% names(nrep_wide))
        nrep_wide[, (g) := 0L]
    if (!g %in% names(mean_wide))
        mean_wide[, (g) := NA_real_]
}

setnames(nrep_wide, genotype_levels, paste0("nrep_", genotype_levels))
setnames(mean_wide, genotype_levels, paste0("mean_Lm_", genotype_levels))

eligibility <- merge(nrep_wide, mean_wide, by = c("coord_handle", "sex"), all = TRUE, sort = FALSE)

eligibility[, complete_3genotypes := nrep_AG >= min_reps_per_genotype & nrep_WBG >= min_reps_per_genotype &
    nrep_AS >= min_reps_per_genotype]

setorder(eligibility, sex, coord_handle)

fwrite(eligibility, file.path(output_dir, "02_replicate_eligibility_and_site_means.tsv"),
    sep = "\t")

eligible_keys <- eligibility[complete_3genotypes == TRUE, .(coord_handle, sex, recalculated_nrep_AG = nrep_AG,
    recalculated_nrep_WBG = nrep_WBG, recalculated_nrep_AS = nrep_AS, recalculated_mean_Lm_AG = mean_Lm_AG,
    recalculated_mean_Lm_WBG = mean_Lm_WBG, recalculated_mean_Lm_AS = mean_Lm_AS)]

# -----------------------------------------------------------------------------
# 4. Revalidate the final NB+beta-binomial target-site universe
# -----------------------------------------------------------------------------

target <- fread(target_site_file, check.names = FALSE)
required_columns(target, c("coord_handle", "sex", "AG", "WBG", "AS", "NB_FDR", "BB_FDR",
    "pass_level_cutoff", "pass_pairwise_diff_cutoff"), "target-site table")

target[, `:=`(coord_handle = trimws(as.character(coord_handle)), sex = trimws(as.character(sex)),
    AG = as.numeric(AG), WBG = as.numeric(WBG), AS = as.numeric(AS), NB_FDR = as.numeric(NB_FDR),
    BB_FDR = as.numeric(BB_FDR))]

target <- target[sex %in% c("F", "M") & is.finite(AG) & is.finite(WBG) & is.finite(AS) &
    is.finite(NB_FDR) & NB_FDR < 0.05 & is.finite(BB_FDR) & BB_FDR < 0.05 & pass_level_cutoff ==
    TRUE & pass_pairwise_diff_cutoff == TRUE]

if (target[, anyDuplicated(paste(coord_handle, sex))] != 0L) {
    stop("Duplicated coord_handle + sex keys in target-site table")
}

target_checked <- merge(target, eligible_keys, by = c("coord_handle", "sex"), all.x = TRUE,
    sort = FALSE)

if (target_checked[, any(is.na(recalculated_nrep_AG) | is.na(recalculated_nrep_WBG) | is.na(recalculated_nrep_AS))]) {
    stop("At least one final target site fails the recalculated replicate rule")
}

if (all(c("nrep_AG", "nrep_WBG", "nrep_AS") %in% names(target_checked))) {
    replicate_count_match <- target_checked[, all(nrep_AG == recalculated_nrep_AG & nrep_WBG ==
        recalculated_nrep_WBG & nrep_AS == recalculated_nrep_AS)]
    if (!replicate_count_match) {
        stop("Stored and recalculated biological-replicate counts do not match")
    }
}

# The final target table's AG/WBG/AS columns must equal the means recovered
# directly from the valid biological replicates.
mean_differences <- target_checked[, .(max_abs_difference = max(c(abs(AG - recalculated_mean_Lm_AG),
    abs(WBG - recalculated_mean_Lm_WBG), abs(AS - recalculated_mean_Lm_AS)), na.rm = TRUE))]

if (!is.finite(mean_differences$max_abs_difference) || mean_differences$max_abs_difference >
    1e-12) {
    stop("Target-site means do not match replicate-derived means; max difference = ", mean_differences$max_abs_difference)
}

fwrite(target_checked, file.path(output_dir, "03_final_target_sites_replicate_validated.tsv"),
    sep = "\t")

# -----------------------------------------------------------------------------
# 5. Map target sites to genes
# -----------------------------------------------------------------------------

site_gene_map <- fread(site_gene_map_file, check.names = FALSE)
required_columns(site_gene_map, c("site_key", "gene"), "CDS site-gene map")
setnames(site_gene_map, "site_key", "coord_handle")

site_gene_map[, `:=`(coord_handle = trimws(as.character(coord_handle)), gene = trimws(as.character(gene)))]

site_gene_map <- unique(site_gene_map[!is.na(coord_handle) & coord_handle != "" & !is.na(gene) &
    grepl("^FBgn[0-9]+$", gene), .(coord_handle, gene_id = gene)])

if (site_gene_map[, anyDuplicated(coord_handle)] != 0L) stop("The CDS map assigns a site to multiple genes")

# No-CDS-editing background uses all basic-QC candidates in each sex,
# including candidates outside the significant NB/BB target set.
background <- merge(unique(candidate[, .(coord_handle, sex)]), site_gene_map, by = "coord_handle",
    all = FALSE)
fwrite(unique(background[, .(sex, gene_id)]), file.path(output_dir, "candidate_CDS_genes.tsv"),
    sep = "\t")

target_gene <- merge(target_checked, site_gene_map, by = "coord_handle", all = FALSE, allow.cartesian = TRUE)

if (!nrow(target_gene)) stop("No target sites mapped to genes")

# -----------------------------------------------------------------------------
# 6. Calculate signed site changes and classify genes
# -----------------------------------------------------------------------------

site_changes <- rbindlist(list(target_gene[, .(scope, coord_handle, gene_id, sex, mutant = "AdarG",
    mutant_mean_Lm = AG, WBg_mean_Lm = WBG, site_delta_vs_WBg = AG - WBG, n_mutant_reps = recalculated_nrep_AG,
    n_WBg_reps = recalculated_nrep_WBG)], target_gene[, .(scope, coord_handle, gene_id, sex,
    mutant = "AdarS", mutant_mean_Lm = AS, WBg_mean_Lm = WBG, site_delta_vs_WBg = AS - WBG,
    n_mutant_reps = recalculated_nrep_AS, n_WBg_reps = recalculated_nrep_WBG)]), use.names = TRUE)

site_changes[, site_direction := fcase(site_delta_vs_WBg > editing_change_cutoff, "site_up",
    site_delta_vs_WBg < -editing_change_cutoff, "site_down", default = "site_nochange")]

gene_groups <- site_changes[, .(n_target_sites = uniqueN(coord_handle), mean_site_delta_vs_WBg = mean(site_delta_vs_WBg),
    median_site_delta_vs_WBg = median(site_delta_vs_WBg), sd_site_delta_vs_WBg = if (.N >
        1L) sd(site_delta_vs_WBg) else NA_real_, min_site_delta_vs_WBg = min(site_delta_vs_WBg),
    max_site_delta_vs_WBg = max(site_delta_vs_WBg), n_site_up = uniqueN(coord_handle[site_delta_vs_WBg >
        editing_change_cutoff]), n_site_down = uniqueN(coord_handle[site_delta_vs_WBg < -editing_change_cutoff]),
    n_site_nochange = uniqueN(coord_handle[abs(site_delta_vs_WBg) <= editing_change_cutoff]),
    target_site_ids = paste(sort(unique(coord_handle)), collapse = ";")), by = .(scope, gene_id,
    sex, mutant)]

# Classification is performed AFTER averaging signed changes across all target
# sites in a gene. Exactly +0.05 and -0.05 are included in editing_nochange.
gene_groups[, editing_group := fcase(mean_site_delta_vs_WBg > editing_change_cutoff, "editing_up",
    mean_site_delta_vs_WBg < -editing_change_cutoff, "editing_down", default = "editing_nochange")]

setorder(site_changes, sex, mutant, gene_id, coord_handle)
setorder(gene_groups, sex, mutant, editing_group, gene_id)

fwrite(site_changes, file.path(output_dir, "04_site_signed_changes.tsv"), sep = "\t")

fwrite(gene_groups, file.path(output_dir, "05_gene_editing_change_groups.tsv"), sep = "\t")

