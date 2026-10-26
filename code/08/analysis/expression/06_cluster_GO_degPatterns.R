source(file.path(Sys.getenv("CODE_DIR"), "analysis", "common.R"))
library(clusterProfiler)
list2env(readRDS(analysis_rds), envir = environment())
source(file.path(Sys.getenv("CODE_DIR"), "analysis", "expression", "functions_degPatterns.R"))
rna_background_genes <- list(F = unique(rna_F_res$res_tbl$gene_id[!is.na(rna_F_res$res_tbl$padj)]),
    M = unique(rna_M_res$res_tbl$gene_id[!is.na(rna_M_res$res_tbl$padj)]))

cluster_files_main <- c(F = file.path(count_out_dir, "RNA_F_degPatterns_cluster_results.tsv"),
    M = file.path(count_out_dir, "RNA_M_degPatterns_cluster_results.tsv"))

go_RNA <- run_count_degPatterns_cluster_GO(cluster_files = cluster_files_main, sex_names = c("F",
    "M"), background_genes = rna_background_genes, layer_name = "RNA", ont = "BP", pvalueCutoff = 0.05,
    qvalueCutoff = 0.2, top_n = 10, out_prefix = file.path(count_out_dir, "RNA_degPatterns"))


