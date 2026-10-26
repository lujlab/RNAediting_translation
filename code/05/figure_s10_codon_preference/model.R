suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)
root <- normalizePath(commandArgs(TRUE)[1])
out <- file.path(root, "results_3N_sex_model")
dir.create(out, showWarnings = FALSE)
manifest <- fread(file.path(root, "results_3N/input_manifest.tsv"))
ann <- fread(file.path(root, "results_3N/codon_annotation.tsv"))
qc <- fread(file.path(root, "results_3N/filter_QC.tsv"))
fwrite(manifest, file.path(out, "input_manifest.tsv"), sep = "\t")
stats <- densities <- models <- audit <- list()
for (sp in c("Dmel", "Dsim", "Dsec", "Dere")) {
    message("Processing ", sp)
    sm <- manifest[species == sp]
    d <- rbindlist(lapply(seq_len(nrow(sm)), function(i) {
        unit <- sm$unit_id[i]
        x <- fread(file.path(root, "results_3N", paste0(unit, ".all_positions.tsv.gz")))
        stopifnot(!anyDuplicated(x$FBtr_codonpos), all(is.finite(x$gmRC)), all(x$n_missing ==
            0), all(x$n_verified == x$n_effective), all(x$n_effective == x$n_gmRC_reps_finite),
            all(x$threshold == 3 * x$n_effective), all(x$pass_3N == (x$Asite_median_sum >=
                x$threshold)))
        x <- x[pass_3N == TRUE]
        stopifnot(nrow(x) == qc[unit_id == unit, N_pass])
        x[, sex := sm$sex[i]]
        x
    }))
    stopifnot(!anyDuplicated(d[, .(FBtr_codonpos, sex)]))
    n_before <- nrow(d)
    d <- merge(d, ann, by = "codon", all.x = TRUE, sort = FALSE)
    stopifnot(nrow(d) == n_before, all(d[is.na(RSCU), codon] %chin% c("TAA", "TAG", "TGA") |
        !grepl("^[ACGT]{3}$", d[is.na(RSCU), codon])))
    fwrite(d, file.path(out, paste0(sp, "_position_by_sex_source_data.tsv.gz")), sep = "\t",
        compress = "gzip")
    fwrite(d[is.na(RSCU)], file.path(out, paste0(sp, "_codons_not_classifiable.tsv")), sep = "\t")
    n_unclassified <- d[is.na(RSCU), .N]
    d <- d[!is.na(RSCU)]
    audit[[sp]] <- data.table(species = sp, N_pass_3N = n_before, N_unclassified = n_unclassified,
        N_observations = nrow(d), N_unique_positions = uniqueN(d$FBtr_codonpos), N_female = sum(d$sex ==
            "female"), N_male = sum(d$sex == "male"), sex_averaged = FALSE, all_3N_memberships_verified = TRUE)
    st <- d[, .(n_positions = .N, mean_gmRC = mean(gmRC), sd_gmRC = sd(gmRC), median_gmRC = median(gmRC),
        q1_gmRC = as.numeric(quantile(gmRC, 0.25)), q3_gmRC = as.numeric(quantile(gmRC, 0.75))),
        by = preference]
    top <- max(as.numeric(quantile(d$gmRC, 0.99)), 0.1)
    st[, `:=`(species = sp, y_upper = top)]
    stats[[sp]] <- st
    for (j in 1:2) {
        pref <- c("Preferred", "Unpreferred")[j]
        den <- density(d[preference == pref, gmRC], n = 512, from = 0, to = top)
        w <- den$y/max(den$y) * 0.4
        densities[[paste(sp, pref)]] <- data.table(species = sp, preference = pref, x = c(j -
            w, rev(j + w)), y = c(den$x, rev(den$x)), polygon_order = seq_len(1024), y_upper = top)
    }
    # Exact observation-level OLS using codon-by-sex sufficient statistics.
    # Cluster scores are summed across BOTH sexes within each codon identity.
    g <- d[, .(n = .N, sy = sum(gmRC)), by = .(codon, AA, sex, preference)]
    g[, `:=`(sex = factor(sex, levels = c("female", "male")), preference = factor(preference,
        levels = c("Preferred", "Unpreferred")))]
    X <- model.matrix(~AA + sex + preference, g)
    N <- nrow(d)
    p <- ncol(X)
    G <- uniqueN(g$codon)
    bread <- solve(crossprod(X, X * g$n))
    beta <- bread %*% crossprod(X, g$sy)
    scores <- rowsum(X * as.vector(g$sy - g$n * (X %*% beta)), g$codon, reorder = FALSE)
    stopifnot(nrow(scores) == G, G == 61)
    vc <- (G/(G - 1)) * ((N - 1)/(N - p)) * bread %*% crossprod(scores) %*% bread
    k <- match("preferenceUnpreferred", colnames(X))
    b <- beta[k]
    se <- sqrt(vc[k, k])
    models[[sp]] <- data.table(species = sp, model = "AA_and_sex_adjusted_primary", SE_type = "cluster_robust_codon_CR1",
        estimate = b, SE = se, df = G - 1, n_clusters = G, P_one_sided_unpreferred_greater = pt(b/se,
            G - 1, lower.tail = FALSE), CI_low95 = b - qt(0.975, G - 1) * se, CI_high95 = b +
            qt(0.975, G - 1) * se, N_positions = N)
    fit <- lm(sy/n ~ AA + sex + preference, g, weights = n)
    stopifnot(max(abs(coef(fit) - as.vector(beta))) < 1e-08)
    rm(d, g, X)
    gc()
}
fwrite(rbindlist(stats), file.path(out, "FS10_plot_stats.tsv"), sep = "\t")
fwrite(rbindlist(densities), file.path(out, "FS10_density_polygons.tsv"), sep = "\t")
fwrite(rbindlist(models), file.path(out, "RSCU_site_level_AA_adjusted_models.tsv"), sep = "\t")
invisible(NULL)
print(rbindlist(models))
