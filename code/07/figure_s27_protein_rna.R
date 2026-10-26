#!/usr/bin/env Rscript

# Complete reproducible workflow for Supplementary Figure 27:
# featureCounts raw counts -> joint DESeq2 normalization -> DIA-NN protein
# quantities -> CDS edited-gene definition -> paired AS/AG analysis -> PDF.
#
# Run without pseudocounts:
#   Rscript run_from_initial_quantification_to_SuppFig30.R --pseudocount=0
#
# Optional output directory:
#   --out-dir=/absolute/path/to/output

suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(ggplot2)
    library(grid)
})

set.seed(20200827)

# =============================================================================
# 1. Configuration
# =============================================================================

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- if (length(script_arg)) script_arg[[1L]] else
    "run_from_initial_quantification_to_SuppFig30.R"
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
user_args <- commandArgs(trailingOnly = TRUE)

get_option <- function(option_name, default) {
    prefix <- paste0("^", option_name, "=")
    hit <- grep(prefix, user_args, value = TRUE)
    if (!length(hit)) return(default)
    sub(prefix, "", hit[[length(hit)]])
}

pseudocount <- as.numeric(get_option("--pseudocount", "0"))
if (length(pseudocount) != 1L || !is.finite(pseudocount) ||
    pseudocount != 0) {
    stop("This publication version requires --pseudocount=0")
}
mode_stub <- if (pseudocount == 0) "no_pseudocount" else "pseudocount1"

out_dir <- get_option(
    "--out-dir",
    file.path(dirname(dirname(script_dir)), "results/07_publication/S27")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_dir <- Sys.getenv("S27_INPUT_DIR", unset = file.path(script_dir, "inputs"))
input <- list(
    RNA_featureCounts = c(
        AG = file.path(input_dir, "AG_STAR_featureCounts_new"),
        WBG = file.path(input_dir, "WBG_STAR_featureCounts_new"),
        AS = file.path(input_dir, "AS_STAR_featureCounts_new")
    ),
    DIA_NN_report = file.path(input_dir, "DIA_NN_Report.tsv"),
    GTF = file.path(input_dir, "dmel-all-r6.12.CDS.clean.gtf"),
    cds_editing_sites = file.path(input_dir, "ISO1.candidate.AtoG.annotations.long.tsv.new")
)

missing_inputs <- unlist(input, use.names = FALSE)
missing_inputs <- missing_inputs[!file.exists(missing_inputs)]
if (length(missing_inputs)) {
    stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))
}

genotype_levels <- c("AG", "WBG", "AS")
plot_genotypes <- c("AG", "AS")
sex_levels <- c("female", "male")

fwrite(
    data.table(
        parameter = c(
            paste0("RNA_featureCounts_", names(input$RNA_featureCounts)),
            "DIA_NN_report", "GTF", "cds_editing_sites",
            "DESeq2_filter", "DESeq2_design", "pseudocount", "output_directory"
        ),
        value = c(
            unname(input$RNA_featureCounts), input$DIA_NN_report, input$GTF,
            input$cds_editing_sites,
            "raw count >=10 in at least 3 of 18 RNA samples",
            "~ genotype + sex; only estimateSizeFactors() is required",
            as.character(pseudocount), normalizePath(out_dir, mustWork = FALSE)
        )
    ),
    file.path(out_dir, "00_analysis_configuration.tsv"),
    sep = "\t", quote = FALSE
)

# =============================================================================
# 2. Read featureCounts tables and construct the 18-sample count matrix
# =============================================================================

read_featurecounts <- function(path, genotype_name) {
    fixed_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    d <- fread(path, skip = "Geneid")
    if (!all(fixed_columns %in% names(d))) {
        stop("Unexpected featureCounts columns: ", path)
    }
    sample_columns <- setdiff(names(d), fixed_columns)
    if (!length(sample_columns)) stop("No RNA samples in: ", path)

    long <- melt(
        d,
        id.vars = "Geneid",
        measure.vars = sample_columns,
        variable.name = "sample_col",
        value.name = "raw_count"
    )
    long[, raw_count := as.numeric(raw_count)]
    long[!is.finite(raw_count), raw_count := 0]
    if (any(long$raw_count < 0)) stop("Negative RNA count in: ", path)
    long[, `:=`(
        gene_id = as.character(Geneid),
        genotype = genotype_name,
        sex = fifelse(
            grepl("-[Ff][0-9]", sample_col), "female",
            fifelse(grepl("-[Mm][0-9]", sample_col), "male", NA_character_)
        )
    )]
    if (anyNA(long$sex)) stop("Cannot infer sex from RNA sample names in: ", path)
    long[, .(gene_id, sample_col, genotype, sex, raw_count)]
}

message("[1/7] Reading featureCounts tables...")
rna_long <- rbindlist(Map(
    read_featurecounts,
    path = unname(input$RNA_featureCounts),
    genotype_name = names(input$RNA_featureCounts)
))
rna_sample_info <- unique(rna_long[, .(sample_col, genotype, sex)])
setorder(rna_sample_info, genotype, sex, sample_col)

rna_replicate_qc <- rna_sample_info[, .N, by = .(genotype, sex)]
if (nrow(rna_sample_info) != 18L || any(rna_replicate_qc$N != 3L)) {
    stop("Expected exactly 3 RNA replicates for each of 3 genotypes x 2 sexes")
}
fwrite(
    rna_sample_info,
    file.path(out_dir, "01_RNA_sample_information.tsv"),
    sep = "\t", quote = FALSE
)

# =============================================================================
# 3. Joint DESeq2 size-factor normalization and gene-level means
# =============================================================================

message("[2/7] DESeq2 joint size-factor normalization...")
rna_wide <- dcast(rna_long, gene_id ~ sample_col, value.var = "raw_count", fill = 0)
rna_sample_columns <- setdiff(names(rna_wide), "gene_id")
rna_sample_info <- rna_sample_info[match(rna_sample_columns, sample_col)]
if (anyNA(rna_sample_info$genotype) || anyNA(rna_sample_info$sex)) {
    stop("RNA sample metadata do not match the count matrix")
}

count_matrix <- as.matrix(rna_wide[, ..rna_sample_columns])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- rna_wide$gene_id
keep_gene <- rowSums(count_matrix >= 10L) >= 3L
if (!any(keep_gene)) stop("No genes pass the DESeq2 prefilter")

col_data <- data.frame(
    genotype = factor(rna_sample_info$genotype, levels = genotype_levels),
    sex = factor(rna_sample_info$sex, levels = sex_levels),
    row.names = rna_sample_columns
)
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix[keep_gene, , drop = FALSE],
    colData = col_data,
    design = ~ genotype + sex
)
dds <- estimateSizeFactors(dds)
normalized_matrix <- counts(dds, normalized = TRUE)

size_factor_table <- data.table(
    sample_col = names(sizeFactors(dds)),
    DESeq2_size_factor = as.numeric(sizeFactors(dds))
)
size_factor_table <- merge(
    size_factor_table, rna_sample_info,
    by = "sample_col", all.x = TRUE, sort = FALSE
)
fwrite(
    size_factor_table,
    file.path(out_dir, "02_RNA_DESeq2_joint_size_factors.tsv"),
    sep = "\t", quote = FALSE
)
fwrite(
    data.table(
        total_genes = nrow(count_matrix),
        retained_genes = sum(keep_gene),
        removed_genes = sum(!keep_gene),
        n_RNA_samples = ncol(count_matrix),
        filter_rule = "raw count >=10 in at least 3 samples"
    ),
    file.path(out_dir, "03_RNA_DESeq2_filter_QC.tsv"),
    sep = "\t", quote = FALSE
)

normalized_long <- melt(
    as.data.table(normalized_matrix, keep.rownames = "gene_id"),
    id.vars = "gene_id",
    variable.name = "sample_col",
    value.name = "DESeq2_normalized_count"
)
normalized_long <- merge(
    normalized_long, rna_sample_info,
    by = "sample_col", all.x = TRUE, sort = FALSE
)
rna_gene <- normalized_long[, .(
    mean_DESeq2_normalized_count = mean(DESeq2_normalized_count),
    n_RNA_samples = .N
), by = .(gene_id, genotype, sex)]
fwrite(
    rna_gene,
    file.path(out_dir, "04_RNA_gene_mean_DESeq2_normalized_count.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
)
retained_genes <- unique(rna_gene$gene_id)

# =============================================================================
# 4. Map DIA-NN gene symbols to FBgn and summarize direct PG.Quantity
# =============================================================================

read_unique_symbol_map <- function(gtf_path) {
    gtf <- fread(
        cmd = sprintf("grep -v '^#' %s", shQuote(gtf_path)),
        sep = "\t", header = FALSE, quote = "", fill = TRUE,
        showProgress = FALSE
    )
    if (ncol(gtf) < 9L) stop("GTF has fewer than 9 columns")
    attributes <- as.character(gtf[[9L]])
    gene_id <- ifelse(
        grepl('gene_id "', attributes),
        sub('.*gene_id "([^"]+)".*', '\\1', attributes), NA_character_
    )
    gene_symbol <- ifelse(
        grepl('gene_symbol "', attributes),
        sub('.*gene_symbol "([^"]+)".*', '\\1', attributes), NA_character_
    )
    gene_name <- ifelse(
        grepl('gene_name "', attributes),
        sub('.*gene_name "([^"]+)".*', '\\1', attributes), NA_character_
    )
    symbol <- fifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_name)
    mapping <- unique(data.table(gene_symbol = symbol, gene_id = gene_id))
    mapping <- mapping[
        !is.na(gene_symbol) & gene_symbol != "" & !is.na(gene_id) & gene_id != ""
    ]
    mapping[, n_gene_ids := uniqueN(gene_id), by = gene_symbol]
    unique(mapping[n_gene_ids == 1L, .(gene_symbol, gene_id)])
}

read_dia_nn <- function(report_path, symbol_map) {
    d <- fread(report_path)
    if (!"PG.Genes" %in% names(d)) stop("PG.Genes is absent from DIA-NN report")
    quantity_columns <- grep("\\.PG\\.Quantity$", names(d), value = TRUE)
    if (length(quantity_columns) != 18L) {
        stop("Expected 18 DIA-NN PG.Quantity columns; observed ", length(quantity_columns))
    }

    sample_info <- data.table(sample_col = quantity_columns)
    sample_info[, `:=`(
        genotype = fifelse(
            grepl("_AG_", sample_col), "AG",
            fifelse(
                grepl("_AS_", sample_col), "AS",
                fifelse(grepl("_WBg_", sample_col, ignore.case = TRUE), "WBG", NA_character_)
            )
        ),
        sex = fifelse(
            grepl("_[Ff][123]\\.raw", sample_col), "female",
            fifelse(grepl("_[Mm][123]\\.raw", sample_col), "male", NA_character_)
        )
    )]
    if (anyNA(sample_info$genotype) || anyNA(sample_info$sex)) {
        stop("Cannot infer genotype/sex from DIA-NN quantity columns")
    }

    d[, `:=`(
        report_row = .I,
        gene_symbol_raw = trimws(as.character(PG.Genes))
    )]
    d[, gene_symbol := gsub("^Dmel\\\\", "", gene_symbol_raw)]
    d[, single_gene_symbol :=
        !is.na(gene_symbol) & gene_symbol != "" & !grepl(";", gene_symbol, fixed = TRUE)]
    d <- merge(d, symbol_map, by = "gene_symbol", all.x = TRUE, sort = FALSE)
    mapping_qc <- d[, .N, by = .(
        single_gene_symbol,
        mapping = fifelse(is.na(gene_id), "unmapped_or_ambiguous", "uniquely_mapped")
    )]

    long <- melt(
        d,
        id.vars = c("report_row", "gene_id", "gene_symbol_raw", "single_gene_symbol"),
        measure.vars = quantity_columns,
        variable.name = "sample_col",
        value.name = "protein_quantity"
    )
    long[, protein_quantity := as.numeric(protein_quantity)]
    long[!is.finite(protein_quantity) | protein_quantity < 0, protein_quantity := 0]
    long <- long[single_gene_symbol & !is.na(gene_id)]

    # If several protein groups map to one FBgn, retain the group with the largest
    # total PG.Quantity over all 18 DIA-NN samples.
    group_score <- long[, .(
        total_quantity = sum(protein_quantity),
        n_positive_samples = sum(protein_quantity > 0)
    ), by = .(report_row, gene_id, gene_symbol_raw)]
    group_score <- group_score[
        n_positive_samples > 0
    ][order(gene_id, -total_quantity, report_row)]
    selected_groups <- group_score[, .SD[1L], by = gene_id]
    long <- merge(
        long,
        selected_groups[, .(report_row, gene_id)],
        by = c("report_row", "gene_id"), all = FALSE, sort = FALSE
    )
    long <- merge(long, sample_info, by = "sample_col", all.x = TRUE, sort = FALSE)

    list(
        long = long,
        sample_info = sample_info,
        mapping_qc = mapping_qc,
        selected_groups = selected_groups
    )
}

message("[3/7] Reading and mapping DIA-NN protein quantities...")
symbol_map <- read_unique_symbol_map(input$GTF)
dia <- read_dia_nn(input$DIA_NN_report, symbol_map)
protein_gene <- dia$long[, .(
    mean_protein_quantity = mean(protein_quantity),
    n_MS_samples = .N
), by = .(gene_id, genotype, sex)]
protein_gene <- protein_gene[
    is.finite(mean_protein_quantity) & mean_protein_quantity > 0
]

fwrite(dia$sample_info, file.path(out_dir, "05_DIA_NN_sample_information.tsv"),
    sep = "\t", quote = FALSE
)
fwrite(dia$mapping_qc, file.path(out_dir, "06_DIA_NN_gene_mapping_QC.tsv"),
    sep = "\t", quote = FALSE
)
fwrite(dia$selected_groups, file.path(out_dir, "07_DIA_NN_selected_protein_group_per_FBgn.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
)
fwrite(protein_gene, file.path(out_dir, "08_protein_gene_mean_direct_quantity.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
)

# =============================================================================
# 5. Define genes with >=2 CDS editing sites
# =============================================================================

message("[4/7] Defining CDS edited genes...")
editing_sites <- fread(
    input$cds_editing_sites,
    select = c("handle", "gene", "target_region")
)
editing_sites <- editing_sites[
    target_region == "CDS" & grepl("^FBgn", gene) & !is.na(handle) & handle != "",
    .(gene_id = as.character(gene), handle = as.character(handle))
]
fwrite(unique(editing_sites, by = c("gene_id", "handle")),
    file.path(out_dir, "09_CDS_gene_site_membership.tsv"),
    sep = "\t"
)
editing_counts <- unique(editing_sites, by = c("gene_id", "handle"))[, .(
    n_cds_editing_sites = .N
), by = gene_id]

gene_groups <- merge(
    data.table(gene_id = retained_genes),
    editing_counts,
    by = "gene_id", all.x = TRUE, sort = FALSE
)
gene_groups[is.na(n_cds_editing_sites), n_cds_editing_sites := 0L]
gene_groups[, edit_group := fifelse(
    n_cds_editing_sites == 0L, "unedit",
    fifelse(n_cds_editing_sites >= 2L, "CDS edit >=2", NA_character_)
)]
setorder(gene_groups, gene_id)
fwrite(
    gene_groups,
    file.path(out_dir, "09_gene_edit_groups.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
)

# =============================================================================
# 6. Construct the three final per-gene plotted variables
# =============================================================================

message("[5/7] Constructing RNA, protein, and protein/RNA variables...")

rna_values <- merge(
    rna_gene, gene_groups,
    by = "gene_id", all = FALSE, sort = FALSE
)
rna_values <- rna_values[
    edit_group == "CDS edit >=2" & genotype %in% plot_genotypes &
        is.finite(mean_DESeq2_normalized_count) &
        mean_DESeq2_normalized_count > 0
]
rna_values[, `:=`(
    metric = "RNA",
    log2_value = log2(mean_DESeq2_normalized_count)
)]

protein_values <- merge(
    protein_gene[gene_id %in% retained_genes], gene_groups,
    by = "gene_id", all = FALSE, sort = FALSE
)
protein_values <- protein_values[
    edit_group == "CDS edit >=2" & genotype %in% plot_genotypes &
        is.finite(mean_protein_quantity) & mean_protein_quantity > 0
]
protein_values[, `:=`(
    metric = "protein",
    log2_value = log2(mean_protein_quantity)
)]

ratio_values <- merge(
    protein_gene[gene_id %in% retained_genes],
    rna_gene,
    by = c("gene_id", "genotype", "sex"),
    all = FALSE, sort = FALSE
)
ratio_values <- merge(
    ratio_values, gene_groups,
    by = "gene_id", all = FALSE, sort = FALSE
)
ratio_values <- ratio_values[
    edit_group == "CDS edit >=2" & genotype %in% plot_genotypes &
        is.finite(mean_protein_quantity) &
        is.finite(mean_DESeq2_normalized_count) &
        mean_protein_quantity > 0 &
        mean_DESeq2_normalized_count > 0
]
ratio_values[, `:=`(
    metric = "protein_over_RNA",
    log2_value = log2(mean_protein_quantity) -
        log2(mean_DESeq2_normalized_count)
)]

plot_data <- list(RNA = rna_values, protein = protein_values, ratio = ratio_values)
for (metric_name in names(plot_data)) {
    fwrite(
        plot_data[[metric_name]],
        file.path(out_dir, paste0("10_", metric_name, "_edited_AS_AG_plot_data.tsv")),
        sep = "\t", quote = FALSE, na = "NA"
    )
}

# =============================================================================
# 7. Summary, paired Wilcoxon tests, and Supplementary Figure 27 style PDF
# =============================================================================

font_family <- "Arial"
plot_genotype_colours <- c(AG = "#2572A9", AS = "#D5231E")
pt_to_mm <- ggplot2::.pt

axis_labels <- list(
    RNA = expression(log[2] * "(mean DESeq2-normalized count)"),
    protein = expression(log[2] * "(mean direct protein quantity)"),
    ratio = expression(log[2] * "(protein / RNA)")
)

plot_config <- list(
    RNA = list(
        panel = "A", y_label = axis_labels$RNA,
        limits = c(3.7, 22.2), breaks = c(5, 10, 15, 20),
        summary_n_y = 20.65, summary_value_y = 19.90,
        bracket_y = 21.60, p_y = 22.15, zero_line = FALSE
    ),
    protein = list(
        panel = "B", y_label = axis_labels$protein,
        limits = c(5.5, 26.0), breaks = c(10, 15, 20, 25),
        summary_n_y = 23.85, summary_value_y = 23.02,
        bracket_y = 25.70, p_y = 26.20, zero_line = FALSE
    ),
    ratio = list(
        panel = "C", y_label = axis_labels$ratio,
        limits = c(-8.0, 10.2), breaks = c(-5, 0, 5, 10),
        summary_n_y = 8.35, summary_value_y = 7.58,
        bracket_y = 9.58, p_y = 10.15, zero_line = TRUE
    )
)

format_p_expression <- function(p_value) {
    exponent <- floor(log10(p_value))
    mantissa <- p_value / 10^exponent
    sprintf('italic(P)*" = "*"%.2f"%%*%%10^{%d}', mantissa, exponent)
}

calculate_plot_statistics <- function(d, metric_name, conf) {
    d[, `:=`(
        genotype = factor(as.character(genotype), levels = plot_genotypes),
        sex = factor(as.character(sex), levels = sex_levels)
    )]
    summary_dt <- d[, .(
        n = .N,
        mean = mean(log2_value),
        sd = sd(log2_value),
        median = median(log2_value)
    ), by = .(sex, genotype)]
    summary_dt[, `:=`(
        metric = metric_name,
        x = as.numeric(genotype),
        n_y = conf$summary_n_y,
        value_y = conf$summary_value_y,
        n_label = paste0("italic(n)==", n),
        value_label = sprintf("%.2f \u00B1 %.2f", mean, sd)
    )]

    paired_dt <- rbindlist(lapply(sex_levels, function(sex_name) {
        ag <- d[sex == sex_name & genotype == "AG", .(gene_id, AG = log2_value)]
        as <- d[sex == sex_name & genotype == "AS", .(gene_id, AS = log2_value)]
        paired <- merge(ag, as, by = "gene_id")
        fwrite(paired, file.path(out_dir, paste0("paired_genes_", metric_name, "_", sex_name, ".tsv")), sep = "\t")
        test <- suppressWarnings(
            wilcox.test(paired$AG, paired$AS, paired = TRUE, exact = FALSE)
        )
        data.table(
            metric = metric_name,
            sex = sex_name,
            comparison = "AG vs AS",
            n_paired_genes = nrow(paired),
            median_AG = median(paired$AG),
            median_AS = median(paired$AS),
            median_paired_difference_AG_minus_AS = median(paired$AG - paired$AS),
            statistic = unname(test$statistic),
            p_value = test$p.value
        )
    }))
    paired_dt[, `:=`(
        sex = factor(sex, levels = sex_levels),
        x_start = 1, x_end = 2,
        bracket_y = conf$bracket_y,
        tick_bottom = conf$bracket_y - 0.34,
        p_y = conf$p_y,
        p_label = vapply(p_value, format_p_expression, character(1))
    )]
    list(data = d, summary = summary_dt, paired = paired_dt)
}

metric_results <- Map(
    calculate_plot_statistics,
    d = plot_data,
    metric_name = names(plot_data),
    conf = plot_config
)

all_summary <- rbindlist(lapply(metric_results, `[[`, "summary"), fill = TRUE)
all_paired <- rbindlist(lapply(metric_results, `[[`, "paired"), fill = TRUE)
fwrite(
    all_summary[, .(metric, sex, genotype, n, mean, sd, median)],
    file.path(out_dir, "11_final_group_summary.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
    all_paired[, .(
        metric, sex, comparison, n_paired_genes, median_AG, median_AS,
        median_paired_difference_AG_minus_AS, statistic, p_value
    )],
    file.path(out_dir, "12_final_paired_Wilcoxon.tsv"),
    sep = "\t", quote = FALSE, na = "NA"
)

make_plot <- function(result, conf) {
    d <- result$data
    summary_dt <- result$summary
    paired_dt <- result$paired
    p <- ggplot(d, aes(genotype, log2_value, fill = genotype))
    if (isTRUE(conf$zero_line)) {
        p <- p + geom_hline(
            yintercept = 0, linetype = "dashed", linewidth = 0.247, colour = "grey45"
        )
    }
    p +
        geom_violin(
            trim = FALSE, width = 0.82, alpha = 0.62,
            colour = "#1A1A1A", linewidth = 0.247
        ) +
        geom_boxplot(
            width = 0.12, outlier.shape = NA, fill = "white",
            alpha = 0.85, linewidth = 0.247
        ) +
        stat_summary(
            fun.data = function(z) {
                mean_value <- mean(z)
                sd_value <- sd(z)
                data.frame(
                    y = mean_value,
                    ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value
                )
            },
            geom = "errorbar", width = 0.10, linewidth = 0.247, colour = "black"
        ) +
        stat_summary(
            fun = mean, geom = "point", shape = 21, size = 1.45,
            fill = "black", colour = "black"
        ) +
        geom_text(
            data = summary_dt,
            aes(x = x, y = n_y, label = n_label),
            inherit.aes = FALSE, parse = TRUE,
            family = font_family, size = 7.3 / pt_to_mm
        ) +
        geom_text(
            data = summary_dt,
            aes(x = x, y = value_y, label = value_label),
            inherit.aes = FALSE, family = font_family, size = 7.3 / pt_to_mm
        ) +
        geom_segment(
            data = paired_dt,
            aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
            inherit.aes = FALSE, linewidth = 0.247
        ) +
        geom_segment(
            data = paired_dt,
            aes(x = x_start, xend = x_start, y = tick_bottom, yend = bracket_y),
            inherit.aes = FALSE, linewidth = 0.247
        ) +
        geom_segment(
            data = paired_dt,
            aes(x = x_end, xend = x_end, y = tick_bottom, yend = bracket_y),
            inherit.aes = FALSE, linewidth = 0.247
        ) +
        geom_text(
            data = paired_dt,
            aes(x = 1.5, y = p_y, label = p_label),
            inherit.aes = FALSE, parse = TRUE,
            family = font_family, size = 8 / pt_to_mm
        ) +
        facet_wrap(
            ~sex,
            nrow = 1,
            labeller = as_labeller(c(female = "Female", male = "Male"))
        ) +
        scale_fill_manual(values = plot_genotype_colours) +
        scale_x_discrete(labels = c(
            AG = expression(italic(Adar)^italic(G)),
            AS = expression(italic(Adar)^italic(S))
        )) +
        scale_y_continuous(breaks = conf$breaks, expand = expansion(mult = 0)) +
        coord_cartesian(ylim = conf$limits, clip = "off") +
        labs(x = NULL, y = conf$y_label) +
        theme_classic(base_size = 8, base_family = font_family) +
        theme(
            legend.position = "none",
            panel.spacing.x = unit(5.5, "pt"),
            axis.line = element_line(linewidth = 0.247, colour = "black"),
            axis.ticks = element_line(linewidth = 0.247, colour = "black"),
            axis.ticks.length = unit(2.2, "pt"),
            axis.text = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8.7, colour = "black", margin = margin(r = 4)),
            strip.background = element_blank(),
            strip.text = element_text(
                size = 10, face = "plain", colour = "#1A1A1A",
                margin = margin(b = if (conf$panel == "C") 9 else 8)
            ),
            plot.margin = margin(
                t = if (conf$panel == "C") 10 else 7,
                r = 4, b = 1, l = 10, unit = "pt"
            )
        )
}

message("[6/7] Drawing the S27 protein/RNA ratio panel...")
plots <- Map(make_plot, result = metric_results, conf = plot_config)
ggsave(file.path(out_dir, "FigureS27_CDSge2_protein_RNA_ratio_only.pdf"),
    plots$ratio,
    width = 4.1, height = 3.0, device = cairo_pdf
)

invisible(NULL)

message("[7/7] Complete")
message("PDF: ", file.path(out_dir, "FigureS27_CDSge2_protein_RNA_ratio_only.pdf"))
