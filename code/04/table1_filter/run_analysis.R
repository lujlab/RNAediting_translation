suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else "..")
outdir <- normalizePath(if (length(args) > 1) args[2] else ".")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "plot_helpers.R"))
base <- Sys.getenv("TABLE1_SCALED_ROOT", file.path(root, "results/04_Table1_TableS4S5_RNAdiv2/RNAdiv2"))
table_dir <- Sys.getenv("TABLE1_RESULT_DIR", file.path(base, "Table1_edgeR_QL_Fisher_minrep2_pairedcoverage"))
samples <- fread(file.path(table_dir, "included_independent_samples.tsv"))
universe <- fread(cmd = paste("gzip -dc", shQuote(file.path(table_dir, "common_test_universe.tsv.gz"))))
universe <- universe[position == "A" & analysis_level == "sex_specific"]
counts <- rbindlist(lapply(seq_len(nrow(samples)), function(i) {
    s <- samples[i]
    x <- fread(file.path(base, "alledit", paste0(s$sample_id, ".Asite_df_alledit")))
    x[, site_id := as.character(FBt_cdspos)]
    x[is.na(site_id) | !nzchar(site_id), site_id := as.character(handle)]
    stopifnot(!anyNA(x$site_id), all(nzchar(x$site_id)))
    cols <- c("ref_number_rna", "alt_number_rna", "ref_number_ribo", "alt_number_ribo")
    for (cc in cols) set(x, which(is.na(x[[cc]])), cc, 0)
    y <- x[, .(Am = sum(ref_number_rna), Im = sum(alt_number_rna), Ar = sum(ref_number_ribo),
        Ir = sum(alt_number_ribo)), by = site_id]
    y[, `:=`(sample_id = s$sample_id, analysis_group = paste0("D", s$species, "_", s$sex))]
    y
}))
mean_finite <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
values <- list()
repvalues <- list()
qc <- list()
for (g in unique(universe$analysis_group)) {
    u <- universe[analysis_group == g]
    z <- counts[analysis_group == g]
    reps <- sort(unique(z$sample_id))
    d <- merge(CJ(site_id = sort(unique(z$site_id)), sample_id = reps), z, by = c("site_id",
        "sample_id"), all.x = TRUE)
    for (cc in c("Am", "Im", "Ar", "Ir")) set(d, which(is.na(d[[cc]])), cc, 0)
    d[, `:=`(RNA_total = Am + Im, Ribo_total = Ar + Ir, analysis_group = g)]
    f <- d[, .(n_group_replicates = .N, RNA_sum = sum(RNA_total), Ribo_sum = sum(Ribo_total),
        pooled_Lm = if (sum(RNA_total) > 0)
            sum(Im)/sum(RNA_total) else NA_real_, paired_coverage_reps = sum(RNA_total >= 5 & Ribo_total >= 5)), by = site_id]
    f[, test_universe := n_group_replicates >= 2 & paired_coverage_reps >= 2 & is.finite(pooled_Lm) &
        pooled_Lm >= 0.05]
    stopifnot(setequal(f$site_id, u$site_id))
    a <- merge(f, u, by = "site_id", suffixes = c("_new", "_Table1"))
    stopifnot(all(a$test_universe_new == a$test_universe_Table1), max(abs(a$RNA_sum_new -
        a$RNA_sum_Table1)) < 1e-08, max(abs(a$Ribo_sum_new - a$Ribo_sum_Table1)) < 1e-08,
        all(a$paired_coverage_reps_new == a$paired_coverage_reps_Table1))
    qc[[g]] <- data.table(analysis_group = g, n_replicates = length(reps), n_input_sites = nrow(f),
        n_selected = sum(f$test_universe), exact_Table1_universe_match = TRUE)
    d <- d[site_id %in% f[test_universe == TRUE, site_id]]
    # Undefined ratios remain missing; no additional coverage or alternate-count gate.
    d[, `:=`(Lm = fifelse(RNA_total > 0, Im/RNA_total, NA_real_), Lr = fifelse(Ribo_total >
        0, Ir/Ribo_total, NA_real_), Ar_Am = fifelse(Am > 0, Ar/Am, NA_real_), Ir_Im = fifelse(Im >
        0, Ir/Im, NA_real_))]
    repvalues[[g]] <- d
    v <- d[, .(mean_Lm = mean_finite(Lm), mean_Lr = mean_finite(Lr), Ar_Am = mean_finite(Ar_Am),
        Ir_Im = mean_finite(Ir_Im), n_Lm = sum(is.finite(Lm)), n_Lr = sum(is.finite(Lr)),
        n_ArAm = sum(is.finite(Ar_Am)), n_IrIm = sum(is.finite(Ir_Im))), by = site_id]
    v <- merge(v, f[test_universe == TRUE], by = "site_id")
    v[, `:=`(analysis_group = g, species = sub("_.*", "", g), cohort = "Table1", panel_id = sub("_.*",
        "", g), sex = if (grepl("_F$", g))
        "female" else "male", cdspos_ids = site_id, n_total_edit_sites = 1L)]
    values[[g]] <- v
}
scatter_dt <- rbindlist(values)
invisible(NULL)
fwrite(rbindlist(repvalues), file.path(outdir, "selected_site_replicate_counts_and_ratios.tsv"),
    sep = "\t", na = "NA")
fwrite(scatter_dt, file.path(outdir, "scatter_site_values.tsv"), sep = "\t", na = "NA")
summary_rows <- list()
specs <- list(Lm_Lr = list(x = "mean_Lm", y = "mean_Lr", xl = expression(italic(Lm)), yl = expression(italic(Lr))),
    ArAm_IrIm = list(x = "Ar_Am", y = "Ir_Im", xl = expression(italic(Ar)/italic(Am)), yl = expression(italic(Ir)/italic(Im))))
for (metric in names(specs)) {
    spec <- specs[[metric]]
    for (sp in c("Dmel", "Dsim", "Dsec", "Dere")) for (sx in c("male", "female")) {
        d <- scatter_dt[species == sp & sex == sx]
        ans <- make_paired_scatter(d, spec$x, spec$y, spec$xl, spec$yl, metric, c(0, 1))
        key <- paste(sp, sx, sep = "_")
        ggsave(file.path(outdir, paste0("Table1_filter_", metric, "_", key, ".pdf")),
            ans$plot, width = 3, height = 3, device = cairo_pdf)
        ans$summary[, `:=`(n_Table1_selected = nrow(d), n_undefined_omitted = nrow(d) - n_points)]
        summary_rows[[paste(metric, key)]] <- ans$summary
    }

}
fwrite(rbindlist(summary_rows), file.path(outdir, "scatter_statistics.tsv"), sep = "\t",
    na = "NA")
invisible(NULL)
