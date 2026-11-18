#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
out <- commandArgs(TRUE)[1]
d <- fread(file.path(out, "AA_four_group_species_sex_codon_long.tsv.gz"), select = c("unit_id",
    "species", "sex", "analysis_role", "AA", "group", "gmRC", "valid_codon"))
d <- d[analysis_role == "main_8" & valid_codon & is.finite(gmRC) & gmRC >= 0]
pool <- d[, .(n_codons = .N, mean_gmRC = mean(gmRC), sd_gmRC = sd(gmRC), median_gmRC = median(gmRC)),
    by = .(AA, group)]
fwrite(pool, file.path(out, "AA_three_group_four_species_pool_summary.tsv"), sep = "\t")
units <- d[, .(unit_mean_gmRC = mean(gmRC), n_codons = .N, unit_sd_gmRC = sd(gmRC)), by = .(unit_id,
    species, sex, AA, group)]
fwrite(units, file.path(out, "AA_three_group_main8_species_sex_AA_group_means.tsv"), sep = "\t")
