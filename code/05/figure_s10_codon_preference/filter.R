suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)
args3 <- commandArgs(TRUE)
root3 <- normalizePath(args3[1])
repo3 <- dirname(root3)
cache3 <- file.path(root3, "main_Asite_median_cache")
dir.create(cache3, showWarnings = FALSE)
out3 <- file.path(root3, "results_3N")
dir.create(out3, showWarnings = FALSE)
qc3 <- list()
read_filtered_mean <- function(info) {
    unit <- info$unit_id[[1]]
    message("3N filter: ", unit)
    x <- fread(info$mean_gmRC_file[[1]], select = c("FBtr_codonpos", "FBtr", "codon", "gmRC",
        "n_gmRC_reps_finite", "sample_ids_finite"))
    before <- nrow(x)
    samples <- strsplit(info$input_samples[[1]], ",", fixed = TRUE)[[1]]
    med <- rbindlist(lapply(samples, function(s) {
        path <- file.path(cache3, paste0(s, ".tsv"))
        if (!file.exists(path)) {
            message("Extract main gmRC median: ", s)
            z <- fread(file.path(repo3, "figure2_3_s7_s10_s11_gmrc/per_library", paste0(s, ".codonTE_AEP.AA.gz")),
                select = c("FBtr", "Asite_median"))
            z <- unique(z)
            stopifnot(!anyDuplicated(z$FBtr), all(is.finite(z$Asite_median)))
            fwrite(z, path, sep = "\t")
        }
        z <- fread(path)
        z[, sample := s]
        z
    }))
    u <- unique(x[, .(FBtr, sample_ids_finite, n_effective = n_gmRC_reps_finite)])
    members <- u[, .(sample = strsplit(sample_ids_finite, ",", fixed = TRUE)[[1]]), by = .(FBtr,
        sample_ids_finite, n_effective)]
    stopifnot(all(members$sample %chin% samples))
    joined <- merge(members, med, by = c("FBtr", "sample"), all.x = TRUE, sort = FALSE)
    sums <- joined[, .(n_verified = .N, n_missing = sum(!is.finite(Asite_median)), Asite_median_sum = if (all(is.finite(Asite_median)))
        sum(Asite_median) else NA_real_), by = .(FBtr, sample_ids_finite, n_effective)]
    stopifnot(all(sums$n_verified == sums$n_effective), all(sums$n_missing == 0))
    sums[, `:=`(threshold = 3 * n_effective, pass_3N = Asite_median_sum >= 3 * n_effective)]
    invisible(NULL)
    x <- merge(x, sums, by = c("FBtr", "sample_ids_finite"), all.x = TRUE, sort = FALSE)
    stopifnot(nrow(x) == before, all(x$n_effective == x$n_gmRC_reps_finite), !anyNA(x$pass_3N))
    fwrite(x, file.path(out3, paste0(unit, ".all_positions.tsv.gz")), sep = "\t",
        compress = "gzip")
    qc3[[unit]] <<- data.table(unit_id = unit, N_before = before, N_pass = sum(x$pass_3N),
        N_fail = sum(!x$pass_3N), N_zero_pass = sum(x$pass_3N & x$gmRC == 0), all_memberships_verified = TRUE,
        all_medians_available = TRUE)
    fwrite(rbindlist(qc3), file.path(out3, "filter_QC.tsv"), sep = "\t")
    x <- x[pass_3N == TRUE]
    # Group means remain unchanged. Filter the group, not individual replicates.
    stopifnot(all(x$Asite_median_sum >= 3 * x$n_effective))
    x[, .(FBtr_codonpos, FBtr, codon, gmRC, n_gmRC_reps_finite)]
}
{
    m3 <- fread(file.path(repo3, "figure2_3_s7_s10_s11_gmrc/group_sex_mean/mean_gmRC_table_manifest.tsv"))[analysis_role ==
        "main_8"]
    fwrite(m3, file.path(out3, "input_manifest.tsv"), sep = "\t")
    ann <- fread(file.path(repo3, "config/fs10_codon_annotation.tsv"))
    fwrite(ann, file.path(out3, "codon_annotation.tsv"), sep = "\t")
    for (ii in seq_len(nrow(m3))) {
        invisible(read_filtered_mean(m3[ii]))
        gc()
    }
    quit(save = "no", status = 0)
}
