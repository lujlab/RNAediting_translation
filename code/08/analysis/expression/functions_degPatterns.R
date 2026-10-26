# ============================================================
# Figure 5 core functions
# DESeq2 LRT + VST + DEGreport::degPatterns
# ============================================================
suppressPackageStartupMessages({
    library(DESeq2)
    library(DEGreport)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tibble)
    library(org.Dm.eg.db)
    library(AnnotationDbi)
})

# ------------------------------------------------------------
# Helper: extract one hard cluster assignment per gene from
# DEGreport::degPatterns() output.
# This is deliberately robust to DEGreport versions that name
# the gene column 'genes', 'gene', 'name', or 'gene_id'.
# ------------------------------------------------------------
extract_degpatterns_assignment <- function(deg_res) {

    if (is.null(deg_res$df)) {
        stop("degPatterns result does not contain $df.")
    }

    d <- as.data.frame(deg_res$df)

    gene_candidates <- c("gene_id", "genes", "gene", "name")
    gene_col <- gene_candidates[gene_candidates %in% colnames(d)][1]

    if (length(gene_col) == 0 || is.na(gene_col)) {
        stop("Cannot identify the gene column in degPatterns$df. Columns are: ", paste(colnames(d),
            collapse = ", "))
    }

    if (!"cluster" %in% colnames(d)) {
        stop("degPatterns$df does not contain a 'cluster' column. Columns are: ", paste(colnames(d),
            collapse = ", "))
    }

    ans <- unique(d[, c(gene_col, "cluster"), drop = FALSE])

    colnames(ans)[1] <- "gene_id"
    ans$gene_id <- as.character(ans$gene_id)

    # Convert cluster labels to integers when possible.
    cluster_chr <- as.character(ans$cluster)
    cluster_num <- suppressWarnings(as.integer(cluster_chr))

    if (anyNA(cluster_num)) {
        cluster_num2 <- suppressWarnings(as.integer(gsub("[^0-9]+", "", cluster_chr)))
        cluster_num[is.na(cluster_num)] <- cluster_num2[is.na(cluster_num)]
    }

    if (anyNA(cluster_num)) {
        stop("Failed to convert some degPatterns cluster labels to integers: ", paste(unique(cluster_chr[is.na(cluster_num)]),
            collapse = ", "))
    }

    ans$cluster <- cluster_num

    # Each gene must belong to one and only one hard cluster.
    dup_check <- ans %>%
        dplyr::count(gene_id, name = "n_cluster") %>%
        dplyr::filter(n_cluster > 1)

    if (nrow(dup_check) > 0) {
        stop("Some genes were assigned to more than one degPatterns cluster. Example: ",
            paste(head(dup_check$gene_id, 10), collapse = ", "))
    }

    ans <- ans[order(ans$cluster, ans$gene_id), , drop = FALSE]
    rownames(ans) <- NULL
    ans
}


# ------------------------------------------------------------
# Core analysis for one sex
#
# Important:
#   * LRT has NO time ordering.
#   * full model    = ~ genotype
#   * reduced model = ~ 1
#   * genotype is treated as a categorical factor.
#   * degPatterns is used only after LRT to group genes with
#     similar expression patterns across the three genotypes.
# ------------------------------------------------------------
run_count_LRT_degPatterns_by_sex <- function(count_mat, sample_info, sex_use = "F", layer_name = "RNA",
    padj_cutoff = 0.01, n_clusters = 6, genotype_ref = "WBg", genotype_order = c("AdarG",
        "WBg", "AdarS"), out_dir = ".", fitType = "local", scale_patterns = TRUE, plot_width = 11,
    plot_height = 7.5) {

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # =========================
    # 1. Subset one sex
    # =========================
    sample_info_sex <- sample_info[sample_info$sex == sex_use, , drop = FALSE]

    sample_info_sex <- sample_info_sex[colnames(count_mat)[colnames(count_mat) %in% sample_info_sex$sample],
        , drop = FALSE]

    mat_sex <- count_mat[, sample_info_sex$sample, drop = FALSE]
    mat_sex <- round(as.matrix(mat_sex))
    mode(mat_sex) <- "integer"

    stopifnot(all(colnames(mat_sex) == rownames(sample_info_sex)))

    # Keep the same sex-specific low-count filter used in the
    keep <- rowSums(mat_sex) >= 100
    mat_sex <- mat_sex[keep, , drop = FALSE]

    message(layer_name, " ", sex_use, ": genes after count filter = ", nrow(mat_sex))

    # =========================
    # 2. DESeq2 LRT
    # =========================
    coldata <- sample_info_sex

    # DESeq2 model reference is WBg, but the LRT P value itself is
    # invariant to which genotype is chosen as the reference level.
    coldata$genotype <- factor(coldata$genotype)
    coldata$genotype <- relevel(coldata$genotype, ref = genotype_ref)

    dds <- DESeqDataSetFromMatrix(countData = mat_sex, colData = coldata, design = ~genotype)

    dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = fitType)

    res_LRT <- results(dds, alpha = padj_cutoff)

    res_tbl <- as.data.frame(res_LRT) %>%
        rownames_to_column(var = "gene_id") %>%
        as_tibble() %>%
        arrange(padj)

    sig_tbl <- res_tbl %>%
        filter(!is.na(padj), padj < padj_cutoff)

    sig_genes <- sig_tbl$gene_id

    message(layer_name, " ", sex_use, ": significant LRT genes padj < ", padj_cutoff, " = ",
        length(sig_genes))

    lrt_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_DESeq2_LRT_all_results.tsv"))

    sig_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_DESeq2_LRT_sig_padj",
        padj_cutoff, ".tsv"))

    fwrite(res_tbl, file = lrt_file, sep = "\t")
    fwrite(sig_tbl, file = sig_file, sep = "\t")

    # =========================
    # 3. VST for PCA and clustering
    # =========================
    vsd <- vst(dds, blind = FALSE, fitType = fitType)
    vst_mat <- assay(vsd)

    # =========================
    # 4. PCA
    # =========================
    pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)

    percentVar <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(PC1, PC2, color = genotype, label = name)) + geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 3) + theme_classic(base_size = 14) + labs(title = paste0(layer_name,
        " ", sex_use, " PCA"), x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ",
        percentVar[2], "% variance"))

    pca_pdf <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_PCA.pdf"))

    ggsave(pca_pdf, p_pca, width = 7, height = 6)

    # =========================
    # 5. Check number of LRT genes
    # =========================
    if (length(sig_genes) < 2) {
        warning(layer_name, " ", sex_use, ": too few significant genes for degPatterns.")

        return(list(layer = layer_name, sex = sex_use, dds = dds, res_LRT = res_LRT, res_tbl = res_tbl,
            sig_tbl = sig_tbl, sig_genes = sig_genes, vsd = vsd, vst_mat = vst_mat, pca_plot = p_pca,
            pca_pdf = pca_pdf, degpatterns_input = NULL, degpatterns = NULL, cluster_assignment = NULL,
            cluster_result_table = NULL, cluster_result_file = NULL, cluster_pdf = NULL))
    }

    # =========================
    # 6. LRT-significant genes for degPatterns
    # =========================
    sig_genes <- intersect(sig_genes, rownames(vst_mat))
    cluster_vst <- vst_mat[sig_genes, , drop = FALSE]

    # Remove rows with non-finite or zero-variance VST values.
    finite_keep <- apply(cluster_vst, 1, function(z) all(is.finite(z)))

    cluster_vst <- cluster_vst[finite_keep, , drop = FALSE]

    var_keep <- apply(cluster_vst, 1, sd, na.rm = TRUE) > 0
    cluster_vst <- cluster_vst[var_keep, , drop = FALSE]

    if (nrow(cluster_vst) < 2) {
        stop(layer_name, " ", sex_use, ": fewer than two clusterable LRT-significant genes remain.")
    }

    # Metadata used for degPatterns. genotype_order affects display only;
    # it does NOT introduce a time model into the LRT.
    pattern_coldata <- as.data.frame(coldata)
    pattern_coldata <- pattern_coldata[colnames(cluster_vst), , drop = FALSE]
    pattern_coldata$genotype <- factor(as.character(pattern_coldata$genotype), levels = genotype_order)

    if (anyNA(pattern_coldata$genotype)) {
        stop("Some genotype labels are absent from genotype_order: ", paste(unique(as.character(coldata$genotype)[!as.character(coldata$genotype) %in%
            genotype_order]), collapse = ", "))
    }

    # If the requested number of clusters exceeds the number of genes,
    # reduce it automatically rather than failing cryptically.
    n_clusters_use <- min(as.integer(n_clusters), nrow(cluster_vst))

    if (n_clusters_use < 2) {
        stop("n_clusters must be at least 2.")
    }

    # =========================
    # 7. degPatterns clustering
    # =========================
    # edited-gene enrichment analyses.
    deg_res <- DEGreport::degPatterns(cluster_vst, pattern_coldata, time = "genotype", col = NULL,
        nClusters = n_clusters_use, scale = scale_patterns)

    cluster_assignment <- extract_degpatterns_assignment(deg_res)

    # Only retain genes that were actually assigned by degPatterns.
    clustered_genes <- cluster_assignment$gene_id

    # =========================
    # 8. Genotype-average VST values
    # =========================
    mat_clustered <- cluster_vst[clustered_genes, , drop = FALSE]

    mat_avg <- sapply(genotype_order, function(g) {
        cols <- rownames(pattern_coldata)[pattern_coldata$genotype == g]
        rowMeans(mat_clustered[, cols, drop = FALSE], na.rm = TRUE)
    })

    mat_avg <- as.data.frame(mat_avg)
    mat_avg$gene_id <- rownames(mat_avg)

    # =========================
    # 9. Build complete cluster result table
    # =========================
    cluster_result <- merge(cluster_assignment, mat_avg, by = "gene_id", all.x = TRUE, sort = FALSE)

    cluster_result <- merge(cluster_result, as.data.frame(res_tbl), by = "gene_id", all.x = TRUE,
        sort = FALSE)

    cluster_result <- cluster_result[order(cluster_result$cluster, cluster_result$padj, cluster_result$gene_id),
        , drop = FALSE]

    cluster_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_cluster_results.tsv"))

    fwrite(cluster_result, cluster_file, sep = "\t")

    # Raw long-format data returned by degPatterns, useful for QC/replotting.
    cluster_plot_data_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_plot_data.tsv"))

    fwrite(as.data.frame(deg_res$df), cluster_plot_data_file, sep = "\t")

    # Cluster sizes.
    cluster_sizes <- cluster_assignment %>%
        dplyr::count(cluster, name = "n_genes") %>%
        dplyr::arrange(cluster)

    cluster_sizes_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_cluster_sizes.tsv"))

    fwrite(cluster_sizes, cluster_sizes_file, sep = "\t")

    # One file per cluster for downstream inspection.
    for (k in sort(unique(cluster_result$cluster))) {
        cluster_k <- cluster_result[cluster_result$cluster == k, , drop = FALSE]

        fwrite(cluster_k, file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_cluster",
            k, ".tsv")), sep = "\t")
    }

    # =========================
    # 10. Save old-style degPatterns figure
    # =========================
    cluster_pdf <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_c",
        n_clusters_use, ".pdf"))

    if (!is.null(deg_res$plot)) {
        ggsave(cluster_pdf, plot = deg_res$plot, width = plot_width, height = plot_height)
    } else {
        warning("degPatterns result does not contain $plot; no cluster PDF was saved.")
        cluster_pdf <- NA_character_
    }

    # =========================
    # 11. QC: cluster-centre similarity
    # =========================
    # This is NOT used to redefine clusters. It only flags whether the
    # requested nClusters may have split very similar expression patterns.
    centre_df <- cluster_result %>%
        dplyr::select(gene_id, cluster, dplyr::all_of(genotype_order)) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(genotype_order), names_to = "genotype",
            values_to = "value") %>%
        dplyr::group_by(cluster, genotype) %>%
        dplyr::summarise(centre = mean(value, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = genotype, values_from = centre) %>%
        dplyr::arrange(cluster)

    centre_mat <- as.matrix(centre_df[, genotype_order, drop = FALSE])
    rownames(centre_mat) <- paste0("Cluster", centre_df$cluster)

    cluster_center_cor <- cor(t(centre_mat), use = "pairwise.complete.obs")

    cluster_center_distance <- as.matrix(dist(centre_mat))

    cluster_center_cor_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_cluster_center_correlation.tsv"))

    cluster_center_distance_file <- file.path(out_dir, paste0(layer_name, "_", sex_use, "_degPatterns_cluster_center_distance.tsv"))

    fwrite(data.frame(cluster = rownames(cluster_center_cor), cluster_center_cor, check.names = FALSE),
        cluster_center_cor_file, sep = "\t")

    fwrite(data.frame(cluster = rownames(cluster_center_distance), cluster_center_distance,
        check.names = FALSE), cluster_center_distance_file, sep = "\t")

    similar_pairs <- which(upper.tri(cluster_center_cor) & cluster_center_cor >= 0.9, arr.ind = TRUE)

    if (nrow(similar_pairs) > 0) {
        similar_pair_labels <- apply(similar_pairs, 1, function(z) {
            paste0(rownames(cluster_center_cor)[z[1]], "-", colnames(cluster_center_cor)[z[2]])
        })

        warning(layer_name, " ", sex_use, ": highly similar degPatterns cluster centres (cor >= 0.90): ",
            paste(similar_pair_labels, collapse = ", "), ". Consider testing a smaller n_clusters as a sensitivity analysis.")
    }

    message("Saved: ", lrt_file)
    message("Saved: ", pca_pdf)
    message("Saved: ", cluster_file)
    message("Saved: ", cluster_pdf)

    return(list(layer = layer_name, sex = sex_use, dds = dds, res_LRT = res_LRT, res_tbl = res_tbl,
        sig_tbl = sig_tbl, sig_genes = sig_genes, vsd = vsd, vst_mat = vst_mat, pca_plot = p_pca,
        pca_pdf = pca_pdf, degpatterns_input = cluster_vst, degpatterns_coldata = pattern_coldata,
        degpatterns = deg_res, cluster_assignment = cluster_assignment, cluster_result_table = cluster_result,
        cluster_result_file = cluster_file, cluster_plot_data_file = cluster_plot_data_file,
        cluster_sizes = cluster_sizes, cluster_sizes_file = cluster_sizes_file, cluster_pdf = cluster_pdf,
        cluster_center_cor = cluster_center_cor, cluster_center_distance = cluster_center_distance,
        cluster_center_cor_file = cluster_center_cor_file, cluster_center_distance_file = cluster_center_distance_file))
}


# ------------------------------------------------------------
# Convert gene IDs to FlyBase IDs for clusterProfiler.
# ------------------------------------------------------------
convert_to_flybase <- function(gene_ids) {

    gene_ids <- unique(gene_ids)
    gene_ids <- gene_ids[!is.na(gene_ids)]
    gene_ids <- gene_ids[gene_ids != ""]

    flybase_direct <- gene_ids[grepl("^FBgn", gene_ids)]
    symbol_ids <- gene_ids[!grepl("^FBgn", gene_ids)]

    if (length(symbol_ids) > 0) {
        symbol_to_fb <- AnnotationDbi::mapIds(org.Dm.eg.db, keys = unique(symbol_ids), keytype = "SYMBOL",
            column = "FLYBASE", multiVals = "first")

        flybase_from_symbol <- unname(symbol_to_fb)
        flybase_from_symbol <- flybase_from_symbol[!is.na(flybase_from_symbol)]
    } else {
        flybase_from_symbol <- character(0)
    }

    flybase_all <- unique(c(flybase_direct, flybase_from_symbol))
    flybase_all <- flybase_all[!is.na(flybase_all)]

    flybase_all
}


# ------------------------------------------------------------
# GO enrichment for degPatterns clusters.
#   ont = BP, BH, pvalueCutoff = 0.05, qvalueCutoff = 0.2.
# ------------------------------------------------------------
run_count_degPatterns_cluster_GO <- function(cluster_files, sex_names = c("F", "M"), background_genes,
    layer_name = "RNA", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, top_n = 10,
    out_prefix = "RNA_degPatterns_GO") {

    if (length(cluster_files) != length(sex_names)) {
        stop("cluster_files and sex_names must have the same length.")
    }

    if (!is.list(background_genes) || is.null(names(background_genes)) || !all(sex_names %in%
        names(background_genes))) {
        stop("background_genes must be a named list containing one gene set per sex: ", paste(sex_names,
            collapse = ", "))
    }

    bg_flybase_by_sex <- lapply(background_genes[sex_names], convert_to_flybase)
    names(bg_flybase_by_sex) <- sex_names

    all_go_df <- list()

    for (i in seq_along(cluster_files)) {

        cluster_file <- cluster_files[i]
        sex_use <- sex_names[i]

        background_genes_use <- background_genes[[sex_use]]
        bg_flybase <- bg_flybase_by_sex[[sex_use]]

        message(layer_name, " ", sex_use, " background genes: ", length(background_genes_use))

        message(layer_name, " ", sex_use, " background FlyBase genes: ", length(bg_flybase))

        if (!file.exists(cluster_file)) {
            warning("File does not exist: ", cluster_file)
            next
        }

        cluster_df <- fread(cluster_file)
        cluster_df <- as.data.frame(cluster_df)

        required_cols <- c("gene_id", "cluster")
        missing_cols <- setdiff(required_cols, colnames(cluster_df))

        if (length(missing_cols) > 0) {
            stop("degPatterns cluster result file is missing columns: ", paste(missing_cols,
                collapse = ", "))
        }

        clusters <- sort(unique(cluster_df$cluster))

        for (cl in clusters) {

            sub_df <- cluster_df %>%
                dplyr::filter(cluster == cl)

            gene_cl <- convert_to_flybase(sub_df$gene_id)
            gene_cl <- unique(gene_cl)

            message(layer_name, " ", sex_use, " cluster ", cl, ": input genes = ", nrow(sub_df),
                ", FlyBase genes = ", length(gene_cl))

            if (length(gene_cl) < 5) {
                message("Skip ", layer_name, " ", sex_use, " cluster ", cl, ": too few genes.")
                next
            }

            ego <- enrichGO(gene = gene_cl, universe = bg_flybase, OrgDb = org.Dm.eg.db,
                keyType = "FLYBASE", ont = ont, pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff,
                qvalueCutoff = qvalueCutoff, readable = TRUE)

            ego_df <- as.data.frame(ego)

            if (nrow(ego_df) == 0) {
                message("No enriched GO terms for ", layer_name, " ", sex_use, " cluster ",
                  cl)
                next
            }

            ego_df$layer <- layer_name
            ego_df$sex <- sex_use
            ego_df$cluster <- cl
            ego_df$group <- paste0(layer_name, "_", sex_use, "_C", cl)
            ego_df$input_gene_n <- length(gene_cl)
            ego_df$input_raw_n <- nrow(sub_df)

            all_go_df[[paste0(layer_name, "_", sex_use, "_C", cl)]] <- ego_df
        }
    }

    if (length(all_go_df) == 0) {
        warning("No GO enrichment results found for ", layer_name)
        return(NULL)
    }

    go_all <- dplyr::bind_rows(all_go_df)

    go_outfile <- paste0(out_prefix, "_all_GO_", ont, ".tsv")

    write.table(go_all, file = go_outfile, sep = "\t", quote = FALSE, row.names = FALSE)

    message("GO result table saved to: ", go_outfile)

    go_plot_df <- go_all %>%
        dplyr::group_by(group) %>%
        dplyr::arrange(p.adjust, .by_group = TRUE) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::ungroup()

    go_plot_df$GeneRatio_num <- sapply(go_plot_df$GeneRatio, function(x) {
        tmp <- strsplit(x, "/")[[1]]
        as.numeric(tmp[1])/as.numeric(tmp[2])
    })

    go_plot_df$Description <- factor(go_plot_df$Description, levels = rev(unique(go_plot_df$Description)))

    go_plot_df$group <- factor(go_plot_df$group, levels = sort(unique(go_plot_df$group)))

    p <- ggplot(go_plot_df, aes(x = GeneRatio_num, y = Description, size = Count, color = -log10(p.adjust))) +
        geom_point() + facet_wrap(~group, scales = "free_y", ncol = 3) + theme_bw(base_size = 12) +
        labs(x = "GeneRatio", y = NULL, size = "Gene count", color = "-log10(FDR)", title = paste0(layer_name,
            " GO enrichment of degPatterns clusters: ", ont)) + theme(strip.text = element_text(size = 11,
        face = "bold"), axis.text.y = element_text(size = 8))

    pdf_file <- paste0(out_prefix, "_GO_", ont, ".pdf")

    ggsave(pdf_file, p, width = 20, height = 10)
    message("GO plot saved to: ", pdf_file)

    return(list(go_all = go_all, go_plot_df = go_plot_df, plot = p, background_flybase = bg_flybase_by_sex,
        go_table_file = go_outfile, go_pdf = pdf_file))
}
