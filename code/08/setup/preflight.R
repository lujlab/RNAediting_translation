#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
method <- Sys.getenv("CLUSTER_METHOD", "degPatterns")
if (method != "degPatterns") stop("Publication workflow requires CLUSTER_METHOD=degPatterns")
packages <- c("DESeq2", "DEGreport", "data.table", "MASS", "VGAM", "dplyr", "tidyr", "tibble",
    "ggplot2", "gridExtra", "UpSetR", "ggprism", "pheatmap", "RColorBrewer", "Biobase", "clusterProfiler",
    "GOSemSim", "GO.db", "org.Dm.eg.db", "AnnotationDbi")
versions <- data.table(package = packages, available = vapply(packages, requireNamespace,
    logical(1), quietly = TRUE))
versions[, version := vapply(package, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_,
    character(1))]
invisible(NULL)
inputs <- c(setNames(file.path(Sys.getenv("FEATURECOUNTS_DIR"), paste0(c("WBG", "AG", "AS"),
    "_STAR_featureCounts_new")), c("WBG_counts", "AG_counts", "AS_counts")), setNames(file.path(Sys.getenv("EDIT_INPUT_DIR"),
    paste0(c("WBG", "AG", "AS"), ".candidate.AtoG.annotations.long.tsv")), c("WBG_editing",
    "AG_editing", "AS_editing")), sapply(c("ALL_SITE_FILE", "SITE_ANNOTATION_FILE", "CDSPOS_FILE",
    "EDIT_GO_EXPRESSION_FILE", "PROTEIN_INPUT"), Sys.getenv))
manifest <- data.table(input = names(inputs), path = unname(inputs), exists = file.exists(inputs))
manifest[, md5 := as.character(tools::md5sum(path))]
invisible(NULL)
print(versions)
print(manifest)
if (any(!manifest$exists)) stop("Missing inputs: ", paste(manifest[exists == FALSE, input],
    collapse = ", "))
if (any(!versions$available)) stop("Missing or unloadable R packages: ", paste(versions[available ==
    FALSE, package], collapse = ", "), ". See setup/install_packages.R.")
if (!dir.exists(Sys.getenv("ADAR_ALL_EDIT_DIR"))) stop("ADAR_ALL_EDIT_DIR does not exist")
invisible(NULL)
