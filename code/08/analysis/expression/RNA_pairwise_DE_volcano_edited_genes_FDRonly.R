#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

# Re-analysis of the AdarG / AdarS / WBg RNA-seq featureCounts data.
# The DE model and thresholds deliberately match
# RNA_pairwise_DE_current_UPSET8_violin_fixed.R:
#   - female and male samples are analysed separately
#   - genes must have a summed raw count >= 100 within that sex
#   - DESeq2 Wald test, design = ~ genotype, fitType = 'local'
#   - DE genes are defined only by padj < 0.01; no fold-change cutoff

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else file.path(getwd(), "RNA_pairwise_DE_volcano_edited_genes_FDRonly")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_dir <- paste0("/path/to/RNA_editing_DATAprocessing_new/", "RNA_processing_JWZ/03_featurecounts")
featurecount_files <- c(WBg = file.path(input_dir, "WBG_STAR_featureCounts_new"), AdarG = file.path(input_dir,
    "AG_STAR_featureCounts_new"), AdarS = file.path(input_dir, "AS_STAR_featureCounts_new"))
all_editing_file <- paste0("/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/",
    "01_rebuild/result/sites_candidate_rescue_3R_18907675_all_test2/", "raw_rows_by_group/ISO1.candidate.AtoG.annotations.long.tsv.new")
cds_editing_file <- paste0("/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/",
    "03_cdspos/result/", "cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed/",
    "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.", "cdspos_for_alledit.tsv")
original_editing_summary_file <- paste0("/path/to/RNA_editing_DATAprocessing_new/2012MC_nascent/",
    "analysis_editing_stability/results/gene_stability_editing.summary.tsv")
gtf_file <- paste0("/path/to/RNA_editing_DATAprocessing_new/", "ref_genome/dmel-all-r6.12.gtf")

padj_cutoff <- 0.01
min_sex_row_sum <- 100L
fit_type <- "local"

sample_names <- c("WBg_F1", "WBg_F2", "WBg_F3", "WBg_M1", "WBg_M2", "WBg_M3", "AdarG_F1",
    "AdarG_F2", "AdarG_F3", "AdarG_M1", "AdarG_M2", "AdarG_M3", "AdarS_F1", "AdarS_F2", "AdarS_F3",
    "AdarS_M1", "AdarS_M2", "AdarS_M3")
sample_by_genotype <- list(WBg = sample_names[1:6], AdarG = sample_names[7:12], AdarS = sample_names[13:18])
info_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

assert_files_exist <- function(paths) {
    missing <- paths[!file.exists(paths)]
    if (length(missing))
        stop("Missing input file(s):\n", paste(missing, collapse = "\n"))
}
assert_files_exist(c(featurecount_files, original_editing_summary_file, all_editing_file,
    cds_editing_file, gtf_file))

read_featurecounts <- function(path, new_sample_names) {
    d <- fread(path, skip = "Geneid", header = TRUE, sep = "\t", check.names = FALSE)
    missing_info <- setdiff(info_cols, names(d))
    if (length(missing_info))
        stop("Missing annotation fields in ", path)
    count_cols <- setdiff(names(d), info_cols)
    if (length(count_cols) != length(new_sample_names)) {
        stop("Expected ", length(new_sample_names), " count columns in ", path, ", found ",
            length(count_cols))
    }
    setnames(d, count_cols, new_sample_names)
    if (anyDuplicated(d$Geneid))
        stop("Duplicated Geneid in ", path)
    d
}

counts_by_genotype <- Map(read_featurecounts, unname(featurecount_files), sample_by_genotype[names(featurecount_files)])
names(counts_by_genotype) <- names(featurecount_files)

gene_ids <- counts_by_genotype$WBg$Geneid
for (nm in names(counts_by_genotype)[-1L]) {
    if (!setequal(gene_ids, counts_by_genotype[[nm]]$Geneid)) {
        stop("Gene sets differ between WBg and ", nm)
    }
    setkeyv(counts_by_genotype[[nm]], "Geneid")
    counts_by_genotype[[nm]] <- counts_by_genotype[[nm]][gene_ids]
}

gene_info <- counts_by_genotype$WBg[, ..info_cols]
wbg_cols <- sample_by_genotype$WBg
adarg_cols <- sample_by_genotype$AdarG
adars_cols <- sample_by_genotype$AdarS
count_dt <- cbind(gene_info[, .(gene_id = Geneid)], counts_by_genotype$WBg[, ..wbg_cols],
    counts_by_genotype$AdarG[, ..adarg_cols], counts_by_genotype$AdarS[, ..adars_cols])
stopifnot(identical(names(count_dt)[-1L], sample_names))
fwrite(cbind(gene_info, count_dt[, -1L]), file.path(out_dir, "WBg_AdarG_AdarS_featureCounts_raw_merged.tsv"),
    sep = "\t")

sample_info <- data.frame(sample = sample_names, genotype = sub("^([^_]+)_.*$", "\\1", sample_names),
    sex = sub("^[^_]+_([FM]).*$", "\\1", sample_names), rep = sub("^[^_]+_[FM]([0-9]+)$",
        "\\1", sample_names), stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$sample
fwrite(as.data.table(sample_info), file.path(out_dir, "sample_information.tsv"), sep = "\t")

# Gene symbols from the same FlyBase release used for featureCounts.
gtf <- fread(gtf_file, sep = "\t", header = FALSE, quote = "", fill = TRUE, select = c(3L,
    9L), col.names = c("feature", "attribute"))
gtf <- gtf[feature == "gene"]
gtf[, gene_id := sub(".*gene_id \"([^\"]+)\".*", "\\1", attribute)]
gtf[, gene_symbol := sub(".*gene_symbol \"([^\"]+)\".*", "\\1", attribute)]
gtf[!grepl("gene_symbol \"", attribute), gene_symbol := gene_id]
gene_symbols <- unique(gtf[, .(gene_id, gene_symbol)], by = "gene_id")

# Primary edited-gene definition: exactly the editing summary used by
# 06_AG_AS_WBG_RPKM_violin_by_sex.R.
original_editing <- fread(original_editing_summary_file)
if (length(setdiff(c("gene_id", "n_sites"), names(original_editing)))) {
    stop("Original editing summary lacks gene_id or n_sites")
}
original_editing <- original_editing[!is.na(gene_id) & gene_id != "" & !is.na(n_sites) &
    n_sites > 0, .(n_editing_sites_06_source = max(as.integer(n_sites), na.rm = TRUE)), by = .(gene_id = as.character(gene_id))]
original_editing[, has_editing_site_06_source := TRUE]
original_editing <- merge(original_editing, gene_symbols, by = "gene_id", all.x = TRUE)
setcolorder(original_editing, c("gene_id", "gene_symbol", "n_editing_sites_06_source", "has_editing_site_06_source"))
setorder(original_editing, -n_editing_sites_06_source, gene_id)
fwrite(original_editing, file.path(out_dir, "editing_genes_from_06_script_source.tsv"), sep = "\t")

# Current editing-site list. Multiple transcript annotations of one genomic site
# are collapsed before sites are counted per gene.
editing_raw <- fread(all_editing_file)
required_edit <- c("dmel_chr", "dmel_pos", "gene", "target_region")
if (length(setdiff(required_edit, names(editing_raw)))) {
    stop("Current editing table is missing required columns")
}
editing_raw[, `:=`(site_id = paste(dmel_chr, dmel_pos, sep = "_"), gene_id = as.character(gene),
    target_region = toupper(as.character(target_region)))]
editing_site_gene <- unique(editing_raw[!is.na(gene_id) & gene_id != "", .(site_id, gene_id,
    target_region)])
editing_gene_summary <- editing_site_gene[, .(n_editing_sites = uniqueN(site_id), n_target_CDS_editing_sites = uniqueN(site_id[target_region ==
    "CDS"])), by = gene_id]

# Keep a separate flag based on the supplied CDS-position projection table.
cds_positions <- fread(cds_editing_file)
if (!"site_key" %in% names(cds_positions)) stop("site_key missing from CDS table")
cds_site_gene <- merge(unique(cds_positions[, .(site_id = as.character(site_key))]), unique(editing_site_gene[,
    .(site_id, gene_id)]), by = "site_id", all.x = TRUE, allow.cartesian = TRUE)
if (cds_site_gene[is.na(gene_id), uniqueN(site_id)] > 0L) {
    stop("Some CDS-position sites could not be mapped to genes in the current site list")
}
cds_gene_summary <- cds_site_gene[, .(n_cds_position_list_sites = uniqueN(site_id)), by = gene_id]
editing_gene_summary <- merge(editing_gene_summary, cds_gene_summary, by = "gene_id", all = TRUE)
for (cc in c("n_editing_sites", "n_target_CDS_editing_sites", "n_cds_position_list_sites")) {
    set(editing_gene_summary, which(is.na(editing_gene_summary[[cc]])), cc, 0L)
}
editing_gene_summary[, `:=`(has_editing_site = n_editing_sites > 0L, has_target_CDS_editing_site = n_target_CDS_editing_sites >
    0L, has_cds_position_list_site = n_cds_position_list_sites > 0L)]
editing_gene_summary <- merge(editing_gene_summary, gene_symbols, by = "gene_id", all.x = TRUE)
setcolorder(editing_gene_summary, c("gene_id", "gene_symbol", "n_editing_sites", "has_editing_site",
    "n_target_CDS_editing_sites", "has_target_CDS_editing_site", "n_cds_position_list_sites",
    "has_cds_position_list_site"))
setorder(editing_gene_summary, -n_editing_sites, gene_id)
fwrite(editing_gene_summary, file.path(out_dir, "current_editing_genes_site_counts.tsv"),
    sep = "\t")

comparisons <- list(AdarG_vs_WBg = c("AdarG", "WBg"), AdarS_vs_WBg = c("AdarS", "WBg"), AdarG_vs_AdarS = c("AdarG",
    "AdarS"))

fisher_direction_stats <- function(d, direction_class) {
    if (!direction_class %in% c("Up", "Down", "Up+Down")) {
        stop("Unknown Fisher direction class: ", direction_class)
    }
    event <- if (direction_class == "Up+Down") {
        d$direction != "NS"
    } else {
        d$direction == direction_class
    }
    edited <- d$has_editing_site_06_source
    edited_event <- sum(edited & event)
    edited_not_event <- sum(edited & !event)
    nonedited_event <- sum(!edited & event)
    nonedited_not_event <- sum(!edited & !event)
    contingency <- matrix(c(edited_event, edited_not_event, nonedited_event, nonedited_not_event),
        nrow = 2L, byrow = TRUE, dimnames = list(editing = c("edited", "nonedited"), DE_status = c("event",
            "not_event")))
    ft <- fisher.test(contingency, alternative = "two.sided")
    data.table(direction_class = direction_class, total_genes = nrow(d), DEG_genes = sum(event),
        total_edited_genes = sum(edited), edited_DEG_genes = edited_event, edited_not_DEG = edited_not_event,
        nonedited_DEG = nonedited_event, nonedited_not_DEG = nonedited_not_event, odds_ratio = unname(ft$estimate),
        fisher_p_two_sided = ft$p.value)
}

format_fisher_p <- function(p) {
    if (!is.finite(p))
        return("NA")
    if (p < 0.001)
        return(formatC(p, format = "e", digits = 2))
    formatC(p, format = "f", digits = 3)
}

plot_volcano <- function(d, sex_label, comparison_label) {
    pdat <- d[!is.na(padj) & is.finite(log2FoldChange)]
    positive_padj <- pdat[padj > 0, padj]
    zero_replacement <- if (length(positive_padj))
        min(positive_padj)/10 else 1e-300
    pdat[, plot_padj := fifelse(padj > 0, padj, zero_replacement)]
    # Display-only truncation. Original log2FC and padj remain unchanged in tables.
    pdat[, plot_log2FoldChange := pmax(-2.5, pmin(2.5, log2FoldChange))]
    pdat[, plot_neglog10_padj := pmin(10, -log10(plot_padj))]
    pdat[, highlight := fifelse(has_editing_site_06_source, "Gene with editing site",
        fifelse(direction != "NS", "DEG, no editing site", "Not significant"))]
    pdat[, highlight := factor(highlight, levels = c("Not significant", "DEG, no editing site",
        "Gene with editing site"))]
    fisher_stats <- rbindlist(lapply(c("Up", "Down", "Up+Down"), function(z) fisher_direction_stats(pdat,
        z)))
    line_for <- function(z) {
        q <- fisher_stats[direction_class == z]
        sprintf("%s: N=%s; DEG=%s; edited=%s; edited DEG=%s; Fisher P=%s", toupper(z), format(q$total_genes,
            big.mark = ","), format(q$DEG_genes, big.mark = ","), format(q$total_edited_genes,
            big.mark = ","), format(q$edited_DEG_genes, big.mark = ","), format_fisher_p(q$fisher_p_two_sided))
    }
    ggplot(pdat, aes(plot_log2FoldChange, plot_neglog10_padj)) + geom_hline(yintercept = -log10(padj_cutoff),
        linewidth = 0.35, linetype = "dashed", colour = "#555555") + geom_point(data = pdat[highlight ==
        "Not significant"], aes(colour = highlight), size = 1, alpha = 0.45) + geom_point(data = pdat[highlight ==
        "DEG, no editing site"], aes(colour = highlight), size = 1.15, alpha = 0.75) + geom_point(data = pdat[highlight ==
        "Gene with editing site"], aes(colour = highlight), size = 1.45, alpha = 0.9) + scale_colour_manual(values = c(`Not significant` = "#BDBDBD",
        `DEG, no editing site` = "#4C78A8", `Gene with editing site` = "#D73027"), drop = FALSE) +
        labs(title = paste0(comparison_label, " (", sex_label, ")"), subtitle = paste(line_for("Up"),
            line_for("Down"), line_for("Up+Down"), sep = "\n"), x = expression(log[2] ~ fold ~
            change), y = expression(-log[10] ~ adjusted ~ italic(P)), colour = NULL) + scale_x_continuous(limits = c(-2.5,
        2.5), breaks = c(-2.5, -1.25, 0, 1.25, 2.5), expand = expansion(mult = c(0.02, 0.02))) +
        scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), expand = expansion(mult = c(0.01,
            0.02))) + theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold",
        size = 12), plot.subtitle = element_text(size = 7.2, lineheight = 1.05), legend.position = "bottom",
        legend.text = element_text(size = 8), axis.line = element_line(linewidth = 0.4))
}

count_mat <- as.matrix(count_dt[, -1L])
rownames(count_mat) <- count_dt$gene_id
mode(count_mat) <- "integer"

summary_list <- list()
plot_list <- list()
all_result_list <- list()

for (sex_use in c("F", "M")) {
    samples_use <- sample_info$sample[sample_info$sex == sex_use]
    mat_sex <- count_mat[, samples_use, drop = FALSE]
    keep <- rowSums(mat_sex) >= min_sex_row_sum
    mat_sex <- mat_sex[keep, , drop = FALSE]

    coldata <- sample_info[samples_use, , drop = FALSE]
    coldata$genotype <- factor(coldata$genotype, levels = c("WBg", "AdarG", "AdarS"))
    dds <- DESeqDataSetFromMatrix(countData = mat_sex, colData = coldata, design = ~genotype)
    dds <- DESeq(dds, test = "Wald", fitType = fit_type, quiet = TRUE)

    norm <- as.data.table(counts(dds, normalized = TRUE), keep.rownames = "gene_id")
    fwrite(norm, file.path(out_dir, paste0("RNA_", sex_use, "_pairwise_DESeq2_normalized_counts.tsv")),
        sep = "\t")

    for (comp_name in names(comparisons)) {
        numerator <- comparisons[[comp_name]][1L]
        denominator <- comparisons[[comp_name]][2L]
        res <- results(dds, contrast = c("genotype", numerator, denominator), alpha = padj_cutoff)
        d <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
        d[, `:=`(layer = "RNA", sex = sex_use, comparison = comp_name, numerator = numerator,
            denominator = denominator, direction = fifelse(!is.na(padj) & is.finite(log2FoldChange) &
                padj < padj_cutoff & log2FoldChange > 0, "Up", fifelse(!is.na(padj) & is.finite(log2FoldChange) &
                padj < padj_cutoff & log2FoldChange < 0, "Down", "NS")))]
        d <- merge(d, gene_symbols, by = "gene_id", all.x = TRUE)
        d <- merge(d, original_editing[, setdiff(names(original_editing), "gene_symbol"),
            with = FALSE], by = "gene_id", all.x = TRUE)
        d <- merge(d, editing_gene_summary[, setdiff(names(editing_gene_summary), "gene_symbol"),
            with = FALSE], by = "gene_id", all.x = TRUE)
        for (cc in c("n_editing_sites", "n_target_CDS_editing_sites", "n_cds_position_list_sites")) {
            set(d, which(is.na(d[[cc]])), cc, 0L)
        }
        set(d, which(is.na(d$n_editing_sites_06_source)), "n_editing_sites_06_source", 0L)
        for (cc in c("has_editing_site", "has_target_CDS_editing_site", "has_cds_position_list_site")) {
            set(d, which(is.na(d[[cc]])), cc, FALSE)
        }
        set(d, which(is.na(d$has_editing_site_06_source)), "has_editing_site_06_source",
            FALSE)
        setcolorder(d, c("gene_id", "gene_symbol", "layer", "sex", "comparison", "numerator",
            "denominator", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
            "direction", "n_editing_sites_06_source", "has_editing_site_06_source", "n_editing_sites",
            "has_editing_site", "n_target_CDS_editing_sites", "has_target_CDS_editing_site",
            "n_cds_position_list_sites", "has_cds_position_list_site"))
        d[, padj_is_na := is.na(padj)]
        setorder(d, padj_is_na, padj)
        d[, padj_is_na := NULL]

        all_file <- file.path(out_dir, paste0("RNA_", sex_use, "_", comp_name, "_DESeq2_Wald_all_annotated.tsv"))
        sig_file <- file.path(out_dir, paste0("RNA_", sex_use, "_", comp_name, "_DESeq2_Wald_sig_padj0.01.tsv"))
        fwrite(d, all_file, sep = "\t")
        fwrite(d[direction != "NS"], sig_file, sep = "\t")

        tested <- d[!is.na(padj)]
        summary_list[[paste(sex_use, comp_name, sep = "_")]] <- data.table(sex = sex_use,
            comparison = comp_name, numerator = numerator, denominator = denominator, genes_after_count_filter = nrow(d),
            genes_with_nonNA_padj = nrow(tested), DEG_FDR_lt_0.01 = sum(tested$direction !=
                "NS"), Up_FDR_lt_0.01 = sum(tested$direction == "Up"), Down_FDR_lt_0.01 = sum(tested$direction ==
                "Down"), edited_genes_tested = sum(tested$has_editing_site_06_source), edited_genes_DEG = sum(tested$has_editing_site_06_source &
                tested$direction != "NS"), edited_genes_Up = sum(tested$has_editing_site_06_source &
                tested$direction == "Up"), edited_genes_Down = sum(tested$has_editing_site_06_source &
                tested$direction == "Down"))

        p <- plot_volcano(d, ifelse(sex_use == "F", "female", "male"), comp_name)
        plot_key <- paste(sex_use, comp_name, sep = "_")
        plot_list[[plot_key]] <- p
        all_result_list[[plot_key]] <- d
        ggsave(file.path(out_dir, paste0("volcano_", plot_key, ".pdf")), p, width = 6.3,
            height = 5.2, device = cairo_pdf)
        ggsave(file.path(out_dir, paste0("volcano_", plot_key, ".png")), p, width = 6.3,
            height = 5.2, dpi = 300)
    }
}

summary_dt <- rbindlist(summary_list, use.names = TRUE)
fwrite(summary_dt, file.path(out_dir, "DE_and_editing_gene_summary.tsv"), sep = "\t")
fwrite(rbindlist(all_result_list, use.names = TRUE), file.path(out_dir, "all_six_comparisons_DESeq2_Wald_annotated.tsv"),
    sep = "\t")

fisher_dt <- rbindlist(lapply(names(all_result_list), function(plot_key) {
    d <- all_result_list[[plot_key]]
    d <- d[!is.na(padj) & is.finite(log2FoldChange)]
    z <- rbindlist(lapply(c("Up", "Down", "Up+Down"), function(direction_class) fisher_direction_stats(d,
        direction_class)))
    z[, `:=`(sex = unique(d$sex), comparison = unique(d$comparison), numerator = unique(d$numerator),
        denominator = unique(d$denominator))]
    z
}), use.names = TRUE)
fisher_dt[, fisher_BH_FDR_across_18_tests := p.adjust(fisher_p_two_sided, method = "BH")]
setcolorder(fisher_dt, c("sex", "comparison", "numerator", "denominator", "direction_class",
    "total_genes", "DEG_genes", "total_edited_genes", "edited_DEG_genes", "edited_not_DEG",
    "nonedited_DEG", "nonedited_not_DEG", "odds_ratio", "fisher_p_two_sided", "fisher_BH_FDR_across_18_tests"))
fwrite(fisher_dt, file.path(out_dir, "Fisher_edited_vs_nonedited_DE_enrichment_Up_Down_All.tsv"),
    sep = "\t")

combined <- wrap_plots(plot_list, ncol = 3, guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out_dir, "volcano_all_pairwise_female_male.pdf"), combined, width = 18,
    height = 10.3, device = cairo_pdf)
ggsave(file.path(out_dir, "volcano_all_pairwise_female_male.png"), combined, width = 18,
    height = 10.3, dpi = 300)

manifest <- data.table(item = c("featureCounts_WBg", "featureCounts_AdarG", "featureCounts_AdarS",
    "editing_genes_for_red_highlight", "current_all_editing_sites", "current_CDS_position_sites",
    "gene_annotation_GTF", "sex_specific_count_filter", "DE_model", "DE_test", "fitType",
    "DE_significance", "fold_change_cutoff", "volcano_red_definition"), value = c(unname(featurecount_files),
    original_editing_summary_file, all_editing_file, cds_editing_file, gtf_file, "rowSums(raw counts within each sex) >= 100",
    "~ genotype", "DESeq2 Wald", fit_type, "BH-adjusted P value (padj) < 0.01", "none (direction uses FC sign only)",
    "gene has n_sites > 0 in the editing summary used by 06_AG_AS_WBG_RPKM_violin_by_sex.R"))
fwrite(manifest, file.path(out_dir, "analysis_parameters_and_inputs.tsv"), sep = "\t")

invisible(NULL)

message("Completed. Output: ", normalizePath(out_dir))
print(summary_dt)
