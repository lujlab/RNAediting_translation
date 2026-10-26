#!/usr/bin/env Rscript
# Install analysis dependencies into the library specified by R07_LIBRARY.
lib <- Sys.getenv("R07_LIBRARY")
if (!nzchar(lib)) stop("Set R07_LIBRARY to the project R library directory.")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), BioC_mirror = "https://bioconductor.posit.co",
    timeout = 900)
cran <- c("MASS", "VGAM", "data.table", "dplyr", "tidyr", "tibble", "ggplot2", "gridExtra",
    "UpSetR", "ggprism", "pheatmap", "RColorBrewer", "BiocManager")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, Ncpus = 4)
bioc <- c("DESeq2", "DEGreport", "Biobase", "clusterProfiler", "GOSemSim", "GO.db", "org.Dm.eg.db",
    "AnnotationDbi")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, lib = lib, ask = FALSE, update = FALSE,
    Ncpus = 4)
missing <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Dependencies unavailable: ", paste(missing, collapse = ", "))
invisible(NULL)
