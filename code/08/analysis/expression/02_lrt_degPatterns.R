source(file.path(Sys.getenv("CODE_DIR"), "analysis", "common.R"))
list2env(readRDS(input_rds), envir = environment())
source(file.path(Sys.getenv("CODE_DIR"), "analysis", "expression", "functions_degPatterns.R"))
padj_cutoff_main <- 0.01

# without changing the LRT or downstream statistical framework.
n_clusters_main <- 6

# Display order only. The LRT treats genotype as categorical.
genotype_order_plot <- c("AdarG", "WBg", "AdarS")

# ------------------------------------------------------------
# RNA FEMALE
# ------------------------------------------------------------
rna_F_res <- run_count_LRT_degPatterns_by_sex(count_mat = rna_f, sample_info = sample_info,
    sex_use = "F", layer_name = "RNA", padj_cutoff = padj_cutoff_main, n_clusters = n_clusters_main,
    genotype_ref = "WBg", genotype_order = genotype_order_plot, out_dir = count_out_dir,
    fitType = "local", scale_patterns = TRUE, plot_width = 11, plot_height = 7.5)

# ------------------------------------------------------------
# RNA MALE
# ------------------------------------------------------------
rna_M_res <- run_count_LRT_degPatterns_by_sex(count_mat = rna_f, sample_info = sample_info,
    sex_use = "M", layer_name = "RNA", padj_cutoff = padj_cutoff_main, n_clusters = n_clusters_main,
    genotype_ref = "WBg", genotype_order = genotype_order_plot, out_dir = count_out_dir,
    fitType = "local", scale_patterns = TRUE, plot_width = 11, plot_height = 7.5)

# ------------------------------------------------------------
# RNA GO enrichment
# genes with a valid DESeq2 LRT adjusted P value for each sex.
# ------------------------------------------------------------

saveRDS(list(rna_F_res = rna_F_res, rna_M_res = rna_M_res), analysis_rds)
