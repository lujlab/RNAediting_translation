source(file.path(Sys.getenv("CODE_DIR"), "analysis", "common.R"))
# cluster enrichment
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

# ============================================================
# 1. Analysis settings
# ============================================================

# ------------------------------------------------------------
# Editing-site files
# ------------------------------------------------------------

all_site_file <- required_env("ALL_SITE_FILE")

site_annotation_file <- required_env("SITE_ANNOTATION_FILE")

# alledit CDS-position table
alledit_file <- required_env("CDSPOS_FILE")

# ------------------------------------------------------------
# ------------------------------------------------------------


# ------------------------------------------------------------
# ------------------------------------------------------------

lrt_files <- c(
    F = file.path(
        count_out_dir,
        "RNA_F_DESeq2_LRT_all_results.tsv"
    ),
    M = file.path(
        count_out_dir,
        "RNA_M_DESeq2_LRT_all_results.tsv"
    )
)

# ------------------------------------------------------------
# ------------------------------------------------------------

cluster_files <- c(
    F = file.path(
        count_out_dir,
        "RNA_F_degPatterns_cluster_results.tsv"
    ),
    M = file.path(
        count_out_dir,
        "RNA_M_degPatterns_cluster_results.tsv"
    )
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

out_dir <- file.path(
    count_out_dir,
    "RNA_editing_gene_DE_cluster_enrichment"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

padj_cutoff <- 0.01

# degPatterns is a hard-clustering method, so there is no degPatterns membership filter.


# ============================================================
# 2. Helper
# ============================================================

clean_gene_id <- function(x) {
    x <- as.character(x)

    x <- trimws(x)

    x[
        is.na(x) |
            x == "" |
            x == "NA" |
            x == "."
    ] <- NA_character_

    x
}


safe_fraction <- function(x, n) {
    if (n == 0) {
        return(NA_real_)
    }

    x / n
}


# ============================================================
# 3. Extract edited genes:
#    Definition 1 = all editing sites
# ============================================================

cat("\n")
cat("============================================================\n")
cat("[1] Extract edited genes from ALL editing sites\n")
cat("============================================================\n")

if (!file.exists(all_site_file)) {
    stop("File does not exist: ", all_site_file)
}

all_sites <- fread(
    all_site_file,
    select = c(
        "handle",
        "gene",
        "change2"
    )
)

all_sites[, gene := clean_gene_id(gene)]

all_sites <- all_sites[
    change2 == "A>G"
]

n_all_rows <- nrow(all_sites)

n_all_sites <- uniqueN(
    all_sites$handle
)

n_all_sites_with_gene <- uniqueN(
    all_sites[
        !is.na(gene),
        handle
    ]
)

n_all_sites_without_gene <-
    n_all_sites -
    n_all_sites_with_gene


handle_gene_check <- unique(
    all_sites[
        !is.na(gene),
        .(
            handle,
            gene
        )
    ]
)

handle_gene_conflict <- handle_gene_check[
    ,
    .(
        n_gene = uniqueN(gene)
    ),
    by = handle
][
    n_gene > 1
]

if (nrow(handle_gene_conflict) > 0) {
    stop(
        "Some editing sites map to multiple genes. ",
        "Please inspect before continuing."
    )
}


edited_genes_all_sites <- sort(
    unique(
        all_sites[
            !is.na(gene),
            gene
        ]
    )
)


cat(
    "Rows in all-site file          : ",
    n_all_rows,
    "\n"
)

cat(
    "Unique editing sites           : ",
    n_all_sites,
    "\n"
)

cat(
    "Sites with gene annotation     : ",
    n_all_sites_with_gene,
    "\n"
)

cat(
    "Sites without gene annotation  : ",
    n_all_sites_without_gene,
    "\n"
)

cat(
    "Unique edited genes [all_sites]: ",
    length(edited_genes_all_sites),
    "\n"
)


fwrite(
    data.table(
        gene_id = edited_genes_all_sites
    ),
    file.path(
        out_dir,
        "edited_genes.all_sites.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 4. Build handle -> gene mapping
#    for alledit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("[2] Build editing-site -> gene map\n")
cat("============================================================\n")


if (!file.exists(site_annotation_file)) {
    stop(
        "File does not exist: ",
        site_annotation_file
    )
}


site_annotation <- fread(
    site_annotation_file,
    select = c(
        "handle",
        "gene"
    )
)

site_annotation[
    ,
    gene := clean_gene_id(gene)
]


site_gene_map <- unique(
    site_annotation[
        !is.na(handle) &
            handle != "",
        .(
            handle,
            gene
        )
    ]
)


# ------------------------------------------------------------
# Ensure each handle maps to only one non-NA gene
# ------------------------------------------------------------

map_conflict <- site_gene_map[
    !is.na(gene),
    .(
        n_gene = uniqueN(gene)
    ),
    by = handle
][
    n_gene > 1
]


if (nrow(map_conflict) > 0) {
    fwrite(
        map_conflict,
        file.path(
            out_dir,
            "ERROR.handle_multiple_gene.tsv"
        ),
        sep = "\t"
    )

    stop(
        "At least one handle maps to multiple genes."
    )
}


# Keep one row per handle
site_gene_map <- site_gene_map[
    order(
        handle,
        is.na(gene)
    )
]

site_gene_map <- site_gene_map[
    !duplicated(handle)
]


cat(
    "Unique site annotation handles: ",
    nrow(site_gene_map),
    "\n"
)


# ============================================================
# 5. Extract edited genes:
#    Definition 2 = CDS sites from alledit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("[3] Extract CDS-edited genes from alledit\n")
cat("============================================================\n")


if (!file.exists(alledit_file)) {
    stop(
        "File does not exist: ",
        alledit_file
    )
}


alledit <- fread(
    alledit_file
)


required_alledit_cols <- c(
    "site_key",
    "target_region"
)

missing_cols <- setdiff(
    required_alledit_cols,
    names(alledit)
)

if (length(missing_cols) > 0) {
    stop(
        "alledit file missing columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


# ------------------------------------------------------------
# IMPORTANT:
# only CDS
#
# ------------------------------------------------------------

alledit_CDS <- unique(
    alledit[
        target_region == "CDS",
        .(
            site_key
        )
    ]
)


cat(
    "Total alledit rows             : ",
    nrow(alledit),
    "\n"
)

cat(
    "Unique CDS alledit sites       : ",
    nrow(alledit_CDS),
    "\n"
)


# ------------------------------------------------------------
# Map CDS site_key -> gene
# ------------------------------------------------------------

alledit_CDS_gene <- merge(
    alledit_CDS,
    site_gene_map,
    by.x = "site_key",
    by.y = "handle",
    all.x = TRUE,
    sort = FALSE
)


# ------------------------------------------------------------
# QC: unmatched sites
# ------------------------------------------------------------

unmatched_CDS <- alledit_CDS_gene[
    is.na(gene)
]


if (nrow(unmatched_CDS) > 0) {
    fwrite(
        unmatched_CDS,
        file.path(
            out_dir,
            "alledit_CDS.unmatched_sites.tsv"
        ),
        sep = "\t"
    )

    stop(
        "Some CDS alledit sites cannot be mapped to a gene: ",
        nrow(unmatched_CDS)
    )
}


edited_genes_CDS_alledit <- sort(
    unique(
        alledit_CDS_gene$gene
    )
)


cat(
    "Mapped CDS sites               : ",
    nrow(alledit_CDS_gene),
    "\n"
)

cat(
    "Unique CDS-edited genes        : ",
    length(edited_genes_CDS_alledit),
    "\n"
)


# ------------------------------------------------------------
# Output site -> gene table
# ------------------------------------------------------------

fwrite(
    alledit_CDS_gene,
    file.path(
        out_dir,
        "alledit_CDS.site_to_gene.tsv"
    ),
    sep = "\t"
)


# ------------------------------------------------------------
# Output CDS edited gene list
# ------------------------------------------------------------

fwrite(
    data.table(
        gene_id = edited_genes_CDS_alledit
    ),
    file.path(
        out_dir,
        "edited_genes.CDS_alledit.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 6. Edited-gene-set QC
# ============================================================

editing_set_qc <- data.table(
    editing_definition = c(
        "all_sites",
        "CDS_alledit"
    ),
    unique_sites = c(
        n_all_sites,
        nrow(alledit_CDS)
    ),
    unique_genes = c(
        length(edited_genes_all_sites),
        length(edited_genes_CDS_alledit)
    )
)


fwrite(
    editing_set_qc,
    file.path(
        out_dir,
        "edited_gene_set.QC.tsv"
    ),
    sep = "\t"
)


cat("\nEdited gene-set QC:\n")
print(editing_set_qc)


# ============================================================
# 7. Core function
# ============================================================

run_one_enrichment <- function(
  lrt_file,
  cluster_file,
  edited_genes,
  sex_use,
  editing_definition,
  padj_cutoff = 0.01,
  out_dir = "."
) {
    cat("\n")
    cat("============================================================\n")
    cat(
        "[RUN] sex = ",
        sex_use,
        " | definition = ",
        editing_definition,
        "\n",
        sep = ""
    )
    cat("============================================================\n")


    # ==========================================================
    # 7.1 Read DESeq2 LRT
    # ==========================================================

    if (!file.exists(lrt_file)) {
        stop(
            "LRT file does not exist: ",
            lrt_file
        )
    }


    lrt <- fread(lrt_file)


    required_lrt_cols <- c(
        "gene_id",
        "padj"
    )

    missing_lrt <- setdiff(
        required_lrt_cols,
        names(lrt)
    )

    if (length(missing_lrt) > 0) {
        stop(
            "LRT file missing columns: ",
            paste(
                missing_lrt,
                collapse = ", "
            )
        )
    }


    # one row per gene
    if (anyDuplicated(lrt$gene_id)) {
        stop(
            "Duplicated gene_id found in LRT file."
        )
    }


    # ==========================================================
    # IMPORTANT:
    # padj = NA is excluded
    #
    # No_DE is explicitly:
    # padj >= 0.01
    # ==========================================================

    gene_status <- lrt[
        !is.na(padj),
        .(
            gene_id,
            padj
        )
    ]


    gene_status[
        ,
        DE_status := fifelse(
            padj < padj_cutoff,
            "DE",
            "No_DE"
        )
    ]


    gene_status[
        ,
        edited := gene_id %in% edited_genes
    ]


    cat(
        "Genes with valid padj : ",
        nrow(gene_status),
        "\n"
    )

    cat(
        "DE genes              : ",
        gene_status[
            DE_status == "DE",
            .N
        ],
        "\n"
    )

    cat(
        "No-DE genes           : ",
        gene_status[
            DE_status == "No_DE",
            .N
        ],
        "\n"
    )


    # ==========================================================
    # 7.2 PART 1:
    #     DE vs No-DE
    # ==========================================================

    de_edit <- gene_status[
        DE_status == "DE",
        sum(edited)
    ]

    de_nonedit <- gene_status[
        DE_status == "DE",
        sum(!edited)
    ]

    node_edit <- gene_status[
        DE_status == "No_DE",
        sum(edited)
    ]

    node_nonedit <- gene_status[
        DE_status == "No_DE",
        sum(!edited)
    ]


    fisher_DE_matrix <- matrix(
        c(
            de_edit,
            de_nonedit,
            node_edit,
            node_nonedit
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            c(
                "DE",
                "No_DE"
            ),
            c(
                "Edited",
                "NonEdited"
            )
        )
    )


    cat("\nDE vs No-DE contingency table:\n")
    print(fisher_DE_matrix)


    # ----------------------------------------------------------
    # Primary test:
    # H1 = edited genes enriched in DE
    # ----------------------------------------------------------

    fisher_DE_enrichment <- fisher.test(
        fisher_DE_matrix,
        alternative = "greater"
    )


    # ----------------------------------------------------------
    # Also output conventional two-sided test + 95% CI
    # ----------------------------------------------------------

    fisher_DE_two <- fisher.test(
        fisher_DE_matrix,
        alternative = "two.sided"
    )


    de_total <- de_edit + de_nonedit
    node_total <- node_edit + node_nonedit


    DE_result <- data.table(
        sex = sex_use,
        editing_definition =
            editing_definition,
        comparison =
            "DE_vs_No_DE",
        DE_total =
            de_total,
        DE_edited =
            de_edit,
        DE_edited_fraction =
            safe_fraction(
                de_edit,
                de_total
            ),
        NoDE_total =
            node_total,
        NoDE_edited =
            node_edit,
        NoDE_edited_fraction =
            safe_fraction(
                node_edit,
                node_total
            ),
        odds_ratio =
            unname(
                fisher_DE_two$estimate
            ),
        OR_conf_low =
            fisher_DE_two$conf.int[1],
        OR_conf_high =
            fisher_DE_two$conf.int[2],
        fisher_p_enrichment =
            fisher_DE_enrichment$p.value,
        fisher_p_two_sided =
            fisher_DE_two$p.value
    )


    # ==========================================================
    # 7.3 Read degPatterns cluster assignments
    # ==========================================================

    if (!file.exists(cluster_file)) {
        stop(
            "degPatterns cluster file does not exist: ",
            cluster_file
        )
    }


    cluster_input <- fread(
        cluster_file
    )


    required_cluster_cols <- c(
        "gene_id",
        "cluster"
    )


    missing_cluster <- setdiff(
        required_cluster_cols,
        names(cluster_input)
    )


    if (length(missing_cluster) > 0) {
        stop(
            "degPatterns cluster file missing columns: ",
            paste(
                missing_cluster,
                collapse = ", "
            )
        )
    }


    cluster_use <- unique(
        cluster_input[
            ,
            .(
                gene_id,
                cluster
            )
        ],
        by = "gene_id"
    )


    # ==========================================================
    # 7.4 Attach edited status
    # ==========================================================

    cluster_df <- merge(
        cluster_use,
        gene_status[
            ,
            .(
                gene_id,
                DE_status,
                edited
            )
        ],
        by = "gene_id",
        all = FALSE
    )


    # degPatterns cluster file should contain only LRT-significant DE genes
    cluster_df <- cluster_df[
        DE_status == "DE"
    ]


    cluster_df[
        ,
        group := paste0(
            "Cluster",
            cluster
        )
    ]


    # ==========================================================
    # 7.5 DE genes that were not clusterable
    # ==========================================================

    all_DE_genes <- gene_status[
        DE_status == "DE",
        gene_id
    ]


    clustered_DE_genes <- unique(
        cluster_df$gene_id
    )


    DE_unclustered <- setdiff(
        all_DE_genes,
        clustered_DE_genes
    )


    cat(
        "DE genes assigned to cluster   : ",
        length(clustered_DE_genes),
        "\n"
    )

    cat(
        "DE genes not in degPatterns clusters: ",
        length(DE_unclustered),
        "\n"
    )


    if (length(DE_unclustered) > 0) {
        fwrite(
            data.table(
                gene_id = DE_unclustered
            ),
            file.path(
                out_dir,
                paste0(
                    "RNA_",
                    sex_use,
                    ".",
                    editing_definition,
                    ".DE_unclustered_genes.tsv"
                )
            ),
            sep = "\t"
        )
    }


    # ==========================================================
    # 7.6 Build No_DE + cluster table
    # ==========================================================

    noDE_df <- gene_status[
        DE_status == "No_DE",
        .(
            gene_id,
            edited
        )
    ]


    noDE_df[
        ,
        group := "No_DE"
    ]


    cluster_for_test <- cluster_df[
        ,
        .(
            gene_id,
            edited,
            group
        )
    ]


    group_df <- rbindlist(
        list(
            noDE_df,
            cluster_for_test
        ),
        use.names = TRUE
    )


    # ==========================================================
    # 7.7 Summary:
    #     edited fraction in each group
    # ==========================================================

    group_summary <- group_df[
        ,
        .(
            total_genes =
                .N,
            edited_genes =
                sum(edited),
            nonedited_genes =
                sum(!edited),
            edited_fraction =
                mean(edited)
        ),
        by = group
    ]


    group_summary[, order_tmp := 0L]

    group_summary[
        group != "No_DE",
        order_tmp := as.integer(
            sub(
                "^Cluster",
                "",
                group
            )
        )
    ]


    setorder(
        group_summary,
        order_tmp
    )


    group_summary[
        ,
        order_tmp := NULL
    ]


    group_summary[
        ,
        `:=`(
            sex =
                sex_use,
            editing_definition =
                editing_definition
        )
    ]


    # ==========================================================
    # 7.8 Global Fisher:
    #
    #            Edited  NonEdited
    # No_DE
    # Cluster1
    # Cluster2
    # ...
    # ==========================================================

    global_matrix <- as.matrix(
        group_summary[
            ,
            .(
                edited_genes,
                nonedited_genes
            )
        ]
    )


    rownames(global_matrix) <-
        group_summary$group

    colnames(global_matrix) <-
        c(
            "Edited",
            "NonEdited"
        )


    cat("\nGlobal contingency table:\n")
    print(global_matrix)


    set.seed(analysis_seed)


    global_fisher <- tryCatch(
        fisher.test(
            global_matrix,
            workspace = 2e8
        ),
        error = function(e) {
            message(
                "Exact global Fisher failed; ",
                "switch to Fisher Monte Carlo."
            )

            fisher.test(
                global_matrix,
                simulate.p.value = TRUE,
                B = 1000000
            )
        }
    )


    global_result <- data.table(
        sex =
            sex_use,
        editing_definition =
            editing_definition,
        test =
            "No_DE_vs_all_clusters_global",
        fisher_p =
            global_fisher$p.value,
        method =
            global_fisher$method
    )


    # ==========================================================
    # 7.9 Pairwise:
    #     each cluster vs No_DE
    # ==========================================================

    noDE_row <- group_summary[
        group == "No_DE"
    ]


    cluster_names <- group_summary[
        group != "No_DE",
        group
    ]


    pairwise_list <- vector(
        "list",
        length(cluster_names)
    )


    for (
        i in seq_along(
            cluster_names
        )
    ) {
        g <- cluster_names[i]


        g_row <- group_summary[
            group == g
        ]


        pair_mat <- matrix(
            c(
                g_row$edited_genes,
                g_row$nonedited_genes,
                noDE_row$edited_genes,
                noDE_row$nonedited_genes
            ),
            nrow = 2,
            byrow = TRUE,
            dimnames = list(
                c(
                    g,
                    "No_DE"
                ),
                c(
                    "Edited",
                    "NonEdited"
                )
            )
        )


        # H1:
        # edited genes enriched in this cluster
        ft_enrich <- fisher.test(
            pair_mat,
            alternative = "greater"
        )


        ft_two <- fisher.test(
            pair_mat,
            alternative = "two.sided"
        )


        pairwise_list[[i]] <- data.table(
            sex =
                sex_use,
            editing_definition =
                editing_definition,
            group =
                g,
            group_total =
                g_row$total_genes,
            group_edited =
                g_row$edited_genes,
            group_edited_fraction =
                g_row$edited_fraction,
            NoDE_total =
                noDE_row$total_genes,
            NoDE_edited =
                noDE_row$edited_genes,
            NoDE_edited_fraction =
                noDE_row$edited_fraction,
            odds_ratio =
                unname(
                    ft_two$estimate
                ),
            OR_conf_low =
                ft_two$conf.int[1],
            OR_conf_high =
                ft_two$conf.int[2],
            fisher_p_enrichment =
                ft_enrich$p.value,
            fisher_p_two_sided =
                ft_two$p.value
        )
    }


    pairwise_result <- rbindlist(
        pairwise_list
    )


    # ----------------------------------------------------------
    # BH correction:
    # ----------------------------------------------------------

    pairwise_result[
        ,
        fisher_p_enrichment_BH :=
            p.adjust(
                fisher_p_enrichment,
                method = "BH"
            )
    ]


    pairwise_result[
        ,
        fisher_p_two_sided_BH :=
            p.adjust(
                fisher_p_two_sided,
                method = "BH"
            )
    ]


    pairwise_result[
        ,
        enriched :=
            odds_ratio > 1 &
                fisher_p_enrichment_BH < 0.05
    ]


    # ==========================================================
    # ==========================================================

    prefix <- paste0(
        "RNA_",
        sex_use,
        ".",
        editing_definition
    )


    fwrite(
        gene_status,
        file.path(
            out_dir,
            paste0(
                prefix,
                ".gene_DE_editing_status.tsv"
            )
        ),
        sep = "\t"
    )


    fwrite(
        DE_result,
        file.path(
            out_dir,
            paste0(
                prefix,
                ".DE_vs_NoDE.fisher.tsv"
            )
        ),
        sep = "\t"
    )


    fwrite(
        group_summary,
        file.path(
            out_dir,
            paste0(
                prefix,
                ".NoDE_cluster.edited_fraction.tsv"
            )
        ),
        sep = "\t"
    )


    fwrite(
        global_result,
        file.path(
            out_dir,
            paste0(
                prefix,
                ".NoDE_cluster.global_fisher.tsv"
            )
        ),
        sep = "\t"
    )


    fwrite(
        pairwise_result,
        file.path(
            out_dir,
            paste0(
                prefix,
                ".cluster_vs_NoDE.fisher.tsv"
            )
        ),
        sep = "\t"
    )


    cat("\nDE vs No-DE result:\n")
    print(DE_result)

    cat("\nGroup summary:\n")
    print(group_summary)

    cat("\nCluster vs No-DE:\n")
    print(pairwise_result)


    return(
        list(
            DE =
                DE_result,
            group_summary =
                group_summary,
            global =
                global_result,
            pairwise =
                pairwise_result
        )
    )
}


# ============================================================
# 8. Two editing definitions
# ============================================================

editing_definitions <- list(
    all_sites =
        edited_genes_all_sites,
    CDS_alledit =
        edited_genes_CDS_alledit
)


# ============================================================
# ============================================================

all_results <- list()

for (
    sex_use in c(
        "F",
        "M"
    )
) {
    for (
        definition in names(
            editing_definitions
        )
    ) {
        key <- paste(
            sex_use,
            definition,
            sep = "_"
        )


        all_results[[key]] <-
            run_one_enrichment(
                lrt_file =
                    lrt_files[[sex_use]],
                cluster_file =
                    cluster_files[[sex_use]],
                edited_genes =
                    editing_definitions[[definition]],
                sex_use =
                    sex_use,
                editing_definition =
                    definition,
                padj_cutoff =
                    padj_cutoff,
                out_dir =
                    out_dir
            )
    }
}


# ============================================================
# 10. Merge all DE vs No-DE results
# ============================================================

DE_all <- rbindlist(
    lapply(
        all_results,
        `[[`,
        "DE"
    ),
    fill = TRUE
)


# 4 primary tests:
DE_all[
    ,
    fisher_p_enrichment_BH :=
        p.adjust(
            fisher_p_enrichment,
            method = "BH"
        )
]


DE_all[
    ,
    fisher_p_two_sided_BH :=
        p.adjust(
            fisher_p_two_sided,
            method = "BH"
        )
]


DE_all[
    ,
    enriched_in_DE :=
        odds_ratio > 1 &
            fisher_p_enrichment_BH < 0.05
]


fwrite(
    DE_all,
    file.path(
        out_dir,
        "ALL.DE_vs_NoDE.fisher.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 11. Merge all group summaries
# ============================================================

group_summary_all <- rbindlist(
    lapply(
        all_results,
        `[[`,
        "group_summary"
    ),
    fill = TRUE
)


fwrite(
    group_summary_all,
    file.path(
        out_dir,
        "ALL.NoDE_cluster.edited_fraction.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 12. Merge global Fisher results
# ============================================================

global_all <- rbindlist(
    lapply(
        all_results,
        `[[`,
        "global"
    ),
    fill = TRUE
)


global_all[
    ,
    fisher_p_BH :=
        p.adjust(
            fisher_p,
            method = "BH"
        )
]


fwrite(
    global_all,
    file.path(
        out_dir,
        "ALL.NoDE_cluster.global_fisher.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 13. Merge cluster vs No-DE
# ============================================================

pairwise_all <- rbindlist(
    lapply(
        all_results,
        `[[`,
        "pairwise"
    ),
    fill = TRUE
)


fwrite(
    pairwise_all,
    file.path(
        out_dir,
        "ALL.cluster_vs_NoDE.fisher.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 14. Grouped bar plots
#     Make the enrichment results directly interpretable
# ============================================================

# ------------------------------------------------------------
# Which P value should be printed on the bars?
#
# "two_sided":
#   Fisher two-sided test; asks whether the edited-gene
#   proportion differs from No-DE in either direction.
#
# "enrichment":
#   one-sided Fisher test; asks specifically whether edited
#   genes are enriched relative to No-DE.
#
# Because the figure is intended to show "difference P",
# two-sided Fisher is used by default.
# ------------------------------------------------------------

plot_p_mode <- "two_sided"

if (!plot_p_mode %in% c("two_sided", "enrichment")) {
    stop("plot_p_mode must be 'two_sided' or 'enrichment'.")
}


# ------------------------------------------------------------
# Pretty P-value labels
# ------------------------------------------------------------

format_plot_p <- function(p) {
    if (is.na(p)) {
        return("")
    }

    if (p < 2.2e-16) {
        return("BH P < 2.2e-16")
    }

    if (p < 0.001) {
        return(
            paste0(
                "BH P = ",
                formatC(
                    p,
                    format = "e",
                    digits = 1
                )
            )
        )
    }

    paste0(
        "BH P = ",
        formatC(
            p,
            format = "f",
            digits = 3
        )
    )
}


# ------------------------------------------------------------
# Common display labels
# ------------------------------------------------------------

editing_definition_levels <- c(
    "all_sites",
    "CDS_alledit"
)

editing_definition_labels <- c(
    "All editing genes",
    "CDS-edited genes"
)

sex_levels <- c(
    "F",
    "M"
)

sex_labels <- c(
    "Female",
    "Male"
)


# ============================================================
# 14.1 Figure A:
#      DE vs No-DE
# ============================================================

DE_plot_df <- rbindlist(
    list(
        DE_all[
            ,
            .(
                sex,
                editing_definition,
                group = "No_DE",
                total_genes = NoDE_total,
                edited_genes = NoDE_edited,
                edited_fraction = NoDE_edited_fraction
            )
        ],
        DE_all[
            ,
            .(
                sex,
                editing_definition,
                group = "DE",
                total_genes = DE_total,
                edited_genes = DE_edited,
                edited_fraction = DE_edited_fraction
            )
        ]
    ),
    use.names = TRUE
)


# P value is attached only to the DE bar.
if (plot_p_mode == "two_sided") {
    DE_p_df <- DE_all[
        ,
        .(
            sex,
            editing_definition,
            group = "DE",
            plot_p = fisher_p_two_sided_BH
        )
    ]
} else {
    DE_p_df <- DE_all[
        ,
        .(
            sex,
            editing_definition,
            group = "DE",
            plot_p = fisher_p_enrichment_BH
        )
    ]
}


DE_plot_df <- merge(
    DE_plot_df,
    DE_p_df,
    by = c(
        "sex",
        "editing_definition",
        "group"
    ),
    all.x = TRUE,
    sort = FALSE
)


DE_plot_df[
    ,
    edited_percent :=
        100 * edited_fraction
]


DE_plot_df[
    ,
    n_label :=
        paste0(
            edited_genes,
            "/",
            total_genes
        )
]


DE_plot_df[
    ,
    p_label :=
        vapply(
            plot_p,
            format_plot_p,
            character(1)
        )
]


DE_plot_df[
    ,
    bar_label :=
        ifelse(
            p_label == "",
            n_label,
            paste0(
                n_label,
                "\n",
                p_label
            )
        )
]


DE_plot_df[
    ,
    group_plot := factor(
        group,
        levels = c(
            "No_DE",
            "DE"
        ),
        labels = c(
            "No-DE",
            "DE"
        )
    )
]


DE_plot_df[
    ,
    editing_definition_plot := factor(
        editing_definition,
        levels = editing_definition_levels,
        labels = editing_definition_labels
    )
]


DE_plot_df[
    ,
    sex_plot := factor(
        sex,
        levels = sex_levels,
        labels = sex_labels
    )
]


fwrite(
    DE_plot_df,
    file.path(
        out_dir,
        "PLOT_DATA.DE_vs_NoDE.grouped_bar.tsv"
    ),
    sep = "\t"
)


p_DE <- ggplot(
    DE_plot_df,
    aes(
        x = group_plot,
        y = edited_percent,
        fill = editing_definition_plot
    )
) +
    geom_col(
        position = position_dodge(
            width = 0.80
        ),
        width = 0.70
    ) +
    geom_text(
        aes(
            label = bar_label
        ),
        position = position_dodge(
            width = 0.80
        ),
        vjust = -0.25,
        size = 3.2,
        lineheight = 0.90
    ) +
    facet_wrap(
        ~sex_plot,
        nrow = 1
    ) +
    scale_y_continuous(
        expand = expansion(
            mult = c(
                0,
                0.28
            )
        )
    ) +
    labs(
        x = NULL,
        y = "Edited genes (%)",
        fill = NULL,
        title = "Edited-gene enrichment in differentially expressed genes",
        subtitle = "Numbers above bars: edited genes / total genes; P: Fisher exact test vs No-DE"
    ) +
    theme_classic(
        base_size = 12
    ) +
    theme(
        legend.position = "top",
        legend.justification = "left",
        strip.background = element_blank(),
        strip.text = element_text(
            face = "bold",
            size = 12
        ),
        axis.text.x = element_text(
            size = 11
        ),
        axis.text.y = element_text(
            size = 10
        ),
        axis.title.y = element_text(
            size = 11
        ),
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9
        )
    )


ggsave(
    filename = file.path(
        out_dir,
        "RNA_editing_enrichment.DE_vs_NoDE.grouped_bar.pdf"
    ),
    plot = p_DE,
    width = 9.0,
    height = 5.5
)


# ============================================================
# 14.2 Figure B:
#      No-DE + degPatterns clusters
# ============================================================

cluster_plot_df <- copy(
    group_summary_all
)


# ------------------------------------------------------------
# Select the P values used for plotting
# ------------------------------------------------------------

if (plot_p_mode == "two_sided") {
    cluster_p_df <- pairwise_all[
        ,
        .(
            sex,
            editing_definition,
            group,
            plot_p = fisher_p_two_sided_BH
        )
    ]
} else {
    cluster_p_df <- pairwise_all[
        ,
        .(
            sex,
            editing_definition,
            group,
            plot_p = fisher_p_enrichment_BH
        )
    ]
}


# ------------------------------------------------------------
# Merge P onto the corresponding cluster bar.
# No-DE has no P because it is the reference group.
# ------------------------------------------------------------

cluster_plot_df <- merge(
    cluster_plot_df,
    cluster_p_df,
    by = c(
        "sex",
        "editing_definition",
        "group"
    ),
    all.x = TRUE,
    sort = FALSE
)


cluster_plot_df[
    ,
    edited_percent :=
        100 * edited_fraction
]


cluster_plot_df[
    ,
    n_label :=
        paste0(
            edited_genes,
            "/",
            total_genes
        )
]


cluster_plot_df[
    ,
    p_label :=
        vapply(
            plot_p,
            format_plot_p,
            character(1)
        )
]


cluster_plot_df[
    ,
    bar_label :=
        ifelse(
            group == "No_DE",
            n_label,
            paste0(
                n_label,
                "\n",
                p_label
            )
        )
]


# ------------------------------------------------------------
# Determine cluster order numerically
# ------------------------------------------------------------

cluster_names_plot <- unique(
    cluster_plot_df[
        group != "No_DE",
        group
    ]
)

cluster_numbers_plot <- as.integer(
    sub(
        "^Cluster",
        "",
        cluster_names_plot
    )
)

cluster_names_plot <- cluster_names_plot[
    order(
        cluster_numbers_plot
    )
]

group_levels_plot <- c(
    "No_DE",
    cluster_names_plot
)

group_labels_plot <- c(
    "No-DE",
    paste0(
        "Cluster ",
        as.integer(
            sub(
                "^Cluster",
                "",
                cluster_names_plot
            )
        )
    )
)


cluster_plot_df[
    ,
    group_plot := factor(
        group,
        levels = group_levels_plot,
        labels = group_labels_plot
    )
]


cluster_plot_df[
    ,
    editing_definition_plot := factor(
        editing_definition,
        levels = editing_definition_levels,
        labels = editing_definition_labels
    )
]


cluster_plot_df[
    ,
    sex_plot := factor(
        sex,
        levels = sex_levels,
        labels = sex_labels
    )
]


fwrite(
    cluster_plot_df,
    file.path(
        out_dir,
        "PLOT_DATA.NoDE_cluster.grouped_bar.tsv"
    ),
    sep = "\t"
)


p_cluster <- ggplot(
    cluster_plot_df,
    aes(
        x = group_plot,
        y = edited_percent,
        fill = editing_definition_plot
    )
) +
    geom_col(
        position = position_dodge(
            width = 0.82
        ),
        width = 0.72
    ) +
    geom_text(
        aes(
            label = bar_label
        ),
        position = position_dodge(
            width = 0.82
        ),
        vjust = -0.25,
        size = 2.9,
        lineheight = 0.90
    ) +
    facet_wrap(
        ~sex_plot,
        nrow = 1
    ) +
    scale_y_continuous(
        expand = expansion(
            mult = c(
                0,
                0.33
            )
        )
    ) +
    labs(
        x = NULL,
        y = "Edited genes (%)",
        fill = NULL,
        title = "Edited-gene enrichment across expression-response clusters",
        subtitle = "Each Cluster is compared with No-DE by Fisher exact test; labels show edited genes / total genes and BH-adjusted P"
    ) +
    theme_classic(
        base_size = 12
    ) +
    theme(
        legend.position = "top",
        legend.justification = "left",
        strip.background = element_blank(),
        strip.text = element_text(
            face = "bold",
            size = 12
        ),
        axis.text.x = element_text(
            size = 10,
            angle = 25,
            hjust = 1
        ),
        axis.text.y = element_text(
            size = 10
        ),
        axis.title.y = element_text(
            size = 11
        ),
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9
        )
    )


ggsave(
    filename = file.path(
        out_dir,
        "RNA_editing_enrichment.NoDE_clusters.grouped_bar.pdf"
    ),
    plot = p_cluster,
    width = 11.0,
    height = 5.8
)


# ============================================================
# 14.3 Two-page summary PDF
# ============================================================

summary_pdf <- file.path(
    out_dir,
    "RNA_editing_enrichment.grouped_bar_summary.pdf"
)

pdf(
    summary_pdf,
    width = 11.0,
    height = 5.8
)

print(p_DE)
print(p_cluster)

dev.off()


cat("\nGrouped bar plots saved:\n")
cat(
    "  ",
    file.path(
        out_dir,
        "RNA_editing_enrichment.DE_vs_NoDE.grouped_bar.pdf"
    ),
    "\n",
    sep = ""
)

cat(
    "  ",
    file.path(
        out_dir,
        "RNA_editing_enrichment.NoDE_clusters.grouped_bar.pdf"
    ),
    "\n",
    sep = ""
)

cat(
    "  ",
    summary_pdf,
    "\n",
    sep = ""
)


# ============================================================
# 15. Final console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")

cat("\nEdited gene sets:\n")
print(editing_set_qc)

cat("\nDE vs No-DE:\n")
print(DE_all)

cat("\nGlobal No-DE / cluster Fisher:\n")
print(global_all)

cat("\nCluster vs No-DE:\n")
print(pairwise_all)

cat("\nOutput directory:\n")
cat(
    normalizePath(
        out_dir,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n"
)
