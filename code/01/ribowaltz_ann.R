library(data.table)
library(txdbmaker)
# Annotation preparation for Figure S3; this script does not draw QC panels.
ref_dir <- Sys.getenv("REF_DIR")
out_dir <- Sys.getenv("FS3_ANNOTATION_DIR")
if (!nzchar(ref_dir) || !nzchar(out_dir)) {
    stop("Set REF_DIR and FS3_ANNOTATION_DIR before running this script.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gtf_file_dmel <- file.path(ref_dir, "dmel-all-r6.12.gtf")

gtf_file_dsec <- file.path(ref_dir, "dsec-all-r1.3.with_lifted_UTR.rebuilt_exon.sorted.gtf")

gtf_file_dere <- file.path(ref_dir, "dere-all-r1.05.gtf")

gtf_file_dsim <- file.path(ref_dir, "dsim-all-r2.02.gtf")

create_annotation <- function(gtfpath = NULL, txdb = NULL, dataSource = NA, organism = NA) {
    if (length(gtfpath) != 0 & length(txdb) != 0) {
        cat("\n")
        warning("gtfpath and txdb are both specified. Only gtfpath will be considered\n")
        txdb <- NULL
    }

    if (length(gtfpath) == 0 & length(txdb) == 0) {
        cat("\n")
        stop("neither gtfpath nor txdb is specified\n\n")
    }

    if (length(gtfpath) != 0) {
        path_to_gtf <- gtfpath
        txdbanno <- txdbmaker::makeTxDbFromGFF(file = path_to_gtf, format = "gtf", dataSource = dataSource,
            organism = organism)
    } else {
        if (txdb %in% rownames(installed.packages())) {
            library(txdb, character.only = TRUE)
        } else {
            source("https://bioconductor.org/biocLite.R")
            biocLite(txdb, suppressUpdates = TRUE)
            library(txdb, character.only = TRUE)
        }
        txdbanno <- get(txdb)
    }

    exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx", use.names = TRUE))
    utr5 <- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdbanno, use.names = TRUE))
    cds <- suppressWarnings(GenomicFeatures::cdsBy(txdbanno, by = "tx", use.names = TRUE))
    utr3 <- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdbanno, use.names = TRUE))

    exon <- data.table::as.data.table(exon[unique(names(exon))])
    utr5 <- data.table::as.data.table(utr5[unique(names(utr5))])
    cds <- data.table::as.data.table(cds[unique(names(cds))])
    utr3 <- data.table::as.data.table(utr3[unique(names(utr3))])

    # Use the data-table width column, not GenomicRanges::width.
    anno_df <- exon[, .(l_tr = sum(width)), by = .(transcript = group_name)]
    l_utr5 <- utr5[, .(l_utr5 = sum(width)), by = .(transcript = group_name)]
    l_cds <- cds[, .(l_cds = sum(width)), by = .(transcript = group_name)]
    l_utr3 <- utr3[, .(l_utr3 = sum(width)), by = .(transcript = group_name)]

    merge_allx <- function(x, y) merge(x, y, all.x = TRUE)
    anno_df <- Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
    anno_df[is.na(anno_df)] <- 0

    anno_df
}

anno_dmel <- create_annotation(gtfpath = gtf_file_dmel, dataSource = "FlyBase", organism = "Drosophila melanogaster")

fwrite(anno_dmel, file = file.path(out_dir, "anno_dmel.tsv"), sep = "\t", quote = FALSE,
    na = "NA")


anno_dsec <- create_annotation(gtfpath = gtf_file_dsec, dataSource = "FlyBase", organism = "Drosophila sechellia")
fwrite(anno_dsec, file = file.path(out_dir, "anno_dsec_utr.tsv"), sep = "\t", quote = FALSE,
    na = "NA")

anno_dere <- create_annotation(gtfpath = gtf_file_dere, dataSource = "FlyBase", organism = "Drosophila erecta")

fwrite(anno_dere, file = file.path(out_dir, "anno_dere.tsv"), sep = "\t", quote = FALSE,
    na = "NA")


anno_dsim <- create_annotation(gtfpath = gtf_file_dsim, dataSource = "FlyBase", organism = "Drosophila simulans")

fwrite(anno_dsim, file = file.path(out_dir, "anno_dsim.tsv"), sep = "\t", quote = FALSE,
    na = "NA")
