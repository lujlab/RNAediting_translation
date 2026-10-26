#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 03_editing_GO.R EDIT_OUT_DIR")
requested_outdir <- args[[1]]
load(file.path(requested_outdir, "editing_models.RData"))
outdir <- requested_outdir
library(clusterProfiler)
library(org.Dm.eg.db)
library(GOSemSim)
library(data.table)
library(ggplot2)

anno <- dt

head(anno)
colnames(anno)
table(grepl("^FBgn", anno$gene), useNA = "ifany")
head(unique(anno$gene), 20)

anno_gene <- anno[!is.na(coord_handle) & !is.na(gene) & grepl("^FBgn", gene), .(coord_handle,
    gene)]

anno_gene <- unique(anno_gene)

head(anno_gene)

sig_sites <- as.data.table(sig_level_wide)

sig_sites <- sig_sites[!is.na(coord_handle) & !is.na(sex) & is.finite(NB_FDR) & NB_FDR <
    FDR_CUTOFF & is.finite(BB_FDR) & BB_FDR < FDR_CUTOFF & pass_pairwise_diff_cutoff == TRUE]

sig_site_gene <- merge(sig_sites, anno_gene, by = "coord_handle", all.x = TRUE)

sig_site_gene[, .(n_sites = uniqueN(coord_handle), n_gene_annotated_sites = uniqueN(coord_handle[!is.na(gene)]),
    n_genes = uniqueN(gene[!is.na(gene)])), by = sex]


rna_dt <- fread(Sys.getenv("EDIT_GO_EXPRESSION_FILE"))
count_cols <- setdiff(colnames(rna_dt), "gene_id")
MIN_COUNT <- 10
MIN_SAMPLE <- 3

rna_dt[, all_pass_n := rowSums(.SD >= MIN_COUNT, na.rm = TRUE), .SDcols = count_cols]

background_genes_expr <- unique(rna_dt[all_pass_n >= MIN_SAMPLE, gene_id])

length(background_genes_expr)

sig_site_gene_id_female <- sig_site_gene[sig_site_gene$sex == "F"]$gene
sig_site_gene_id_male <- sig_site_gene[sig_site_gene$sex == "M"]$gene

female_ego_BP <- enrichGO(gene = sig_site_gene_id_female, universe = background_genes_expr,
    OrgDb = org.Dm.eg.db, keyType = "FLYBASE", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
    qvalueCutoff = 0.2, readable = TRUE)

male_ego_BP <- enrichGO(gene = sig_site_gene_id_male, universe = background_genes_expr, OrgDb = org.Dm.eg.db,
    keyType = "FLYBASE", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
    readable = TRUE)
male_ego_BP_table <- as.data.frame(male_ego_BP)
head(male_ego_BP_table)
female_ego_BP_table <- as.data.frame(female_ego_BP)
head(female_ego_BP_table)

male_ego_BP_table2 <- male_ego_BP_table
male_ego_BP_table2$sex <- "male"

female_ego_BP_table2 <- female_ego_BP_table
female_ego_BP_table2$sex <- "female"

ego_BP_table <- rbind(male_ego_BP_table, female_ego_BP_table)
ego_BP_table2 <- rbind(male_ego_BP_table2, female_ego_BP_table2)

fwrite(ego_BP_table2, file.path(outdir, "NB_BetaBinomial_overlap_GO_BP_enrichment_table.tsv"),
    sep = "\t")

bp_sem <- godata("org.Dm.eg.db", ont = "BP")

ego_female_BP_s <- simplify(female_ego_BP, cutoff = 0.7, by = "p.adjust", select_fun = min,
    measure = "Wang", semData = bp_sem)

ego_male_BP_s <- simplify(male_ego_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang",
    semData = bp_sem)

female_tab <- as.data.table(as.data.frame(ego_female_BP_s))
male_tab <- as.data.table(as.data.frame(ego_male_BP_s))

# Export all retained terms after sex-specific GO semantic simplification.
# This table is saved before top-30 selection, female/male intersection,
# and top-10 plotting. Count and geneID belong to each retained GO term;
# they are not unions of genes from the removed redundant terms.
go_BP_simplified_table <- rbindlist(list(female = female_tab, male = male_tab), use.names = TRUE,
    fill = TRUE, idcol = "sex")

if (nrow(go_BP_simplified_table) > 0L) {
    setorderv(go_BP_simplified_table, c("sex", "p.adjust", "ID"))
}

go_BP_simplified_file <- file.path(outdir, "NB_BetaBinomial_overlap_GO_BP_simplified_enrichment_table.tsv")
fwrite(go_BP_simplified_table, go_BP_simplified_file, sep = "\t")
message("Saved all post-simplify GO BP terms: ", go_BP_simplified_file, " (female: ", nrow(female_tab),
    "; male: ", nrow(male_tab), ")")

female_top30 <- head(female_tab[order(p.adjust)], 30L)
male_top30 <- head(male_tab[order(p.adjust)], 30L)

overlap_ids <- intersect(female_top30$ID, male_top30$ID)

length(overlap_ids)
overlap_ids

female_overlap <- female_top30[ID %in% overlap_ids]
male_overlap <- male_top30[ID %in% overlap_ids]

female_overlap2 <- female_overlap[, .(ID, Description, Count_female = Count, padj_female = p.adjust)]

male_overlap2 <- male_overlap[, .(ID, Description, Count_male = Count, padj_male = p.adjust)]

overlap_summary <- merge(female_overlap2, male_overlap2, by = c("ID", "Description"), all = FALSE)

overlap_summary[, neglogFDR_female := -log10(padj_female)]
overlap_summary[, neglogFDR_male := -log10(padj_male)]

overlap_summary[, total_Count := Count_female + Count_male]
overlap_summary[, best_sig := pmax(neglogFDR_female, neglogFDR_male)]

top10_overlap <- head(overlap_summary[order(-total_Count, -best_sig)], 10L)

top10_overlap
plot_dt <- rbindlist(list(top10_overlap[, .(ID, Description, sex = "Female", Count = Count_female,
    p.adjust = padj_female, neglogFDR = neglogFDR_female)], top10_overlap[, .(ID, Description,
    sex = "Male", Count = Count_male, p.adjust = padj_male, neglogFDR = neglogFDR_male)]),
    fill = TRUE)

term_order <- rev(top10_overlap$Description)

plot_dt[, Description := factor(Description, levels = term_order)]
plot_dt[, sex := factor(sex, levels = c("Female", "Male"))]

shorten_term <- function(x, width = 42) {
    ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

top10_overlap[, Description_short := shorten_term(Description, 42)]

plot_dt <- merge(plot_dt, top10_overlap[, .(ID, Description, Description_short)], by = c("ID",
    "Description"), all.x = TRUE)

term_order_short <- rev(top10_overlap$Description_short)
plot_dt[, Description_short := factor(Description_short, levels = term_order_short)]

p_overlap <- ggplot(plot_dt, aes(x = sex, y = Description_short)) + geom_point(aes(size = Count,
    color = neglogFDR)) + scale_size_continuous(name = "Count", range = c(3, 10)) + scale_color_gradient(name = "-log10(FDR)",
    low = "#E7A07A", high = "#2B8CBE") + labs(title = expression("BP GO: Genes with altered " *
    L[m]), x = NULL, y = NULL) + theme_bw() + theme(plot.title = element_text(face = "bold",
    hjust = 0.5, size = 12), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 10),
    legend.text = element_text(size = 9))


pdf(file = file.path(outdir, "NB_BetaBinomial_overlap_GO_FM_clusters.pdf"), width = 6, height = 5)
print(p_overlap)
dev.off()


# Run this plot after GO annotation has created anno_gene.


cat("All done.\n")
cat("Output directory:\n")
cat(outdir, "\n")
