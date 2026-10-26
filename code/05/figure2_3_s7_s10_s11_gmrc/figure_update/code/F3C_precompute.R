#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
args <- commandArgs(TRUE)
mean_dir <- args[1]
rescue <- args[2]
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(mean_dir, "mean_gmRC_table_manifest.tsv"))
prefix <- c(Dmel = "ISO1", Dsim = "sim", Dsec = "sec", Dere = "ere")
reference_config <- Sys.getenv("M6A_REFERENCE_CONFIG", "")
if (!nzchar(reference_config) || !file.exists(reference_config)) stop("Missing M6A_REFERENCE_CONFIG")
reference <- fread(reference_config)
stopifnot(all(c("species", "m6a_file") %in% names(reference)), !anyDuplicated(reference$species))
m6a_by_species <- setNames(lapply(reference$m6a_file, function(path) {
    if (!file.exists(path))
        stop("Missing m6A file: ", path)
    unique(scan(path, what = character(), quiet = TRUE))
}), reference$species)
all_stats <- list()
m6a_qc <- list()
for (sp in names(prefix)) {
    z <- fread(file.path(rescue, "03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed",
        paste0(prefix[[sp]], ".candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv")),
        select = c("new_transcript", "target_region", "new_cds_pos_1based"))
    z <- unique(z[toupper(target_region) == "CDS" & is.finite(new_cds_pos_1based) & new_cds_pos_1based >=
        1, .(FBtr = new_transcript, cp = (as.integer(new_cds_pos_1based) - 1L)%/%3L + 1L)])
    edited <- z[, paste(FBtr, cp, sep = "_")]
    excluded <- unique(z[, .(pos = cp + seq.int(-5L, 5L)), by = .(FBtr, cp)][pos >= 1, paste(FBtr,
        pos, sep = "_")])
    stats <- list()
    for (i in which(manifest$species == sp)) {
        info <- manifest[i]
        message("F3C ", info$unit_id)
        d <- fread(info$mean_gmRC_file, select = c("FBtr_codonpos", "AA", "gmRC", "chr",
            "pos"))
        n_before_m6a <- nrow(d)
        d[, handle := fifelse(!is.na(chr) & chr != "" & is.finite(pos), paste(chr, format(pos,
            scientific = FALSE, trim = TRUE), sep = "_"), NA_character_)]
        n_m6a_removed <- d[, sum(!is.na(handle) & handle %chin% m6a_by_species[[sp]])]
        d <- d[is.na(handle) | !(handle %chin% m6a_by_species[[sp]])]
        m6a_qc[[length(m6a_qc) + 1L]] <- data.table(species = sp, unit_id = info$unit_id,
            n_before_m6a = n_before_m6a, n_m6a_removed = n_m6a_removed, n_after_m6a = nrow(d))
        d <- d[is.finite(gmRC) & (FBtr_codonpos %chin% edited | !FBtr_codonpos %chin% excluded)]
        d[, category := fifelse(FBtr_codonpos %chin% edited, "edited", AA)]
        s <- d[, .(n_codons = .N, mean_gmRC = mean(gmRC), sd_gmRC = sd(gmRC)), by = category]
        s[, `:=`(unit_id = info$unit_id, species = sp, sex = info$sex)]
        s[, table_name := if (info$unit_id == "Dmel_Adar5G1__male")
            "Dmel_Adar5G1_mean" else paste0(if (sp == "Dmel")
            "Dmel_ISO1_WBg" else sp, "_", if (info$sex == "female")
            "F" else "M", "_mean")]
        stats[[length(stats) + 1L]] <- s
        rm(d)
        gc()
    }
    s <- rbindlist(stats)
    saveRDS(dcast(s, table_name ~ category, value.var = "mean_gmRC"), file.path(out, paste0(sp,
        "_df_means_AA_density_Asite_mean_dp10_mRNA5")))
    all_stats[[sp]] <- s
}
fwrite(rbindlist(all_stats), file.path(out, "F3C_category_summary.tsv"), sep = "\t")
fwrite(rbindlist(m6a_qc), file.path(out, "F3C_m6A_filter_QC.tsv"), sep = "\t")
