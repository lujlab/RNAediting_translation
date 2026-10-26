source("config.R")
source("R/utils.R")
make_output_dirs(out_dir)
stopifnot(file.exists(editing_file), file.exists(gtf_file))

attr_value <- function(x, key) {
    m <- regexec(paste0("(?:^|;[[:space:]]*)", key, " \\\"([^\\\"]+)\\\""), x, perl = TRUE)
    vapply(regmatches(x, m), function(z) if (length(z) >= 2)
        z[2] else NA_character_, character(1))
}
message("Reading GTF: ", gtf_file)
gtf <- fread(gtf_file, sep = "\t", header = FALSE, quote = "", comment.char = "#", fill = TRUE)
setnames(gtf, paste0("V", 1:9), c("seqname", "source", "feature", "start", "end", "score",
    "strand", "frame", "attribute"))
gtf[, `:=`(gene_id = clean_id(attr_value(attribute, "gene_id")), transcript_id = clean_id(attr_value(attribute,
    "transcript_id")), gene_symbol = attr_value(attribute, "gene_name"))]
gtf[is.na(gene_symbol), gene_symbol := attr_value(attribute, "gene_symbol")]
gtf[, width_nt := as.numeric(end) - as.numeric(start) + 1]
tx <- gtf[!is.na(transcript_id) & !is.na(gene_id) & feature %in% c("exon", "transcript",
    "mRNA")]
exon <- gtf[!is.na(transcript_id) & feature == "exon", .(exon_length = sum(width_nt), exon_number = .N),
    by = .(gene_id, transcript_id)]
if (!nrow(exon)) stop("No exon features were found in GTF; cannot calculate transcript lengths.")
cds <- gtf[!is.na(transcript_id) & toupper(feature) == "CDS", .(CDS_length = sum(width_nt)),
    by = .(gene_id, transcript_id)]
utr5 <- gtf[!is.na(transcript_id) & grepl("five_prime_utr|5utr|five_prime_UTR", feature,
    ignore.case = TRUE), .(UTR5_length = sum(width_nt)), by = .(gene_id, transcript_id)]
utr3 <- gtf[!is.na(transcript_id) & grepl("three_prime_utr|3utr|three_prime_UTR", feature,
    ignore.case = TRUE), .(UTR3_length = sum(width_nt)), by = .(gene_id, transcript_id)]
exon_span <- gtf[!is.na(transcript_id) & feature == "exon", .(span_length = max(as.numeric(end)) -
    min(as.numeric(start)) + 1), by = .(gene_id, transcript_id)]
introns <- merge(exon[, .(gene_id, transcript_id, exon_number, exon_length)], exon_span,
    by = c("gene_id", "transcript_id"))
introns[, `:=`(intron_number = pmax(exon_number - 1L, 0L), total_intron_length = pmax(span_length -
    exon_length, 0))]
introns <- introns[, .(gene_id, transcript_id, intron_number, total_intron_length)]
gene_names <- gtf[!is.na(gene_id), .(gene_symbol = gene_symbol[which(!is.na(gene_symbol))[1]]),
    by = gene_id]
features_tx <- Reduce(function(x, y) merge(x, y, by = c("gene_id", "transcript_id"), all.x = TRUE),
    list(exon, cds, utr5, utr3, introns))
features_tx <- merge(features_tx, gene_names, by = "gene_id", all.x = TRUE)
for (x in c("CDS_length", "UTR5_length", "UTR3_length", "intron_number", "total_intron_length")) set(features_tx,
    which(is.na(features_tx[[x]])), x, 0)
setnames(features_tx, "gene_id", "FBgn")
setnames(features_tx, "exon_length", "transcript_length")
setorder(features_tx, FBgn, -transcript_length, transcript_id)
features_tx[, representative_transcript := seq_len(.N) == 1L, by = FBgn]

# Optional sequence composition, calculated only for representative transcript IDs.
features_tx[, `:=`(GC_content = NA_real_, A_content = NA_real_, editable_A_count = NA_integer_)]
cg_mapping <- data.table(annotation_id = character(), FBgn = character())
if (file.exists(transcript_fasta)) {
    target_ids <- features_tx[representative_transcript == TRUE, transcript_id]
    message("Calculating GC/A content for ", length(target_ids), " representative transcript sequences.")
    lines <- readLines(transcript_fasta, warn = FALSE)
    hdr <- grep("^>", lines)
    ends <- c(hdr[-1] - 1L, length(lines))
    header_text <- lines[hdr]
    has_cg <- grepl("FlyBase_Annotation_IDs:CG[0-9]+", header_text) & grepl("parent=FBgn[0-9]+",
        header_text)
    if (any(has_cg)) {
        cg_mapping <- unique(data.table(annotation_id = sub(".*FlyBase_Annotation_IDs:(CG[0-9]+).*",
            "\\1", header_text[has_cg]), FBgn = sub(".*parent=(FBgn[0-9]+).*", "\\1", header_text[has_cg])))
        write_tsv(cg_mapping, file.path(out_dir, "01_gene_mapping", "r6.12_annotation_ID_to_FBgn.tsv"))
    }
    seqs <- vector("list", length(hdr))
    names(seqs) <- sub(" .*", "", sub("^>", "", header_text))
    for (i in seq_along(hdr)) {
        id <- names(seqs)[i]
        if (id %in% target_ids)
            seqs[[i]] <- toupper(paste0(lines[(hdr[i] + 1L):ends[i]], collapse = ""))
    }
    seqs <- seqs[!vapply(seqs, is.null, logical(1))]
    if (length(seqs)) {
        comp <- rbindlist(lapply(names(seqs), function(id) {
            s <- seqs[[id]]
            n <- nchar(s)
            data.table(transcript_id = clean_id(id), GC_content = (nchar(gsub("[^GC]", "",
                s))/n), A_content = nchar(gsub("[^A]", "", s))/n, editable_A_count = nchar(gsub("[^A]",
                "", s)))
        }))
        features_tx <- merge(features_tx, comp, by = "transcript_id", all.x = TRUE, suffixes = c("",
            ".seq"))
        for (z in c("GC_content", "A_content", "editable_A_count")) {
            features_tx[is.na(get(z)), (z) := get(paste0(z, ".seq"))]
            features_tx[, (paste0(z, ".seq")) := NULL]
        }
    }
} else message("Transcript FASTA not found: sequence covariates will remain NA.")
write_tsv(features_tx, file.path(out_dir, "02_gene_features", "transcript_structure_features.tsv"))
gene_features <- features_tx[representative_transcript == TRUE]
write_tsv(gene_features, file.path(out_dir, "02_gene_features", "gene_structure_features.tsv"))
write_tsv(unique(gene_features[, .(FBgn, gene_symbol)]), file.path(out_dir, "01_gene_mapping",
    "annotation_gene_mapping.tsv"))

message("Reading and aggregating editing sites: ", editing_file)
sites <- fread(editing_file)
gene_col <- find_col(sites, c("gene", "FBgn", "gene_id"))
ratio_col <- find_col(sites, c("ratio", "editing_ratio", "editing_level"), FALSE)
depth_col <- find_col(sites, c("depth", "total_depth"), FALSE)
alt_col <- find_col(sites, c("altcount", "alt_count", "edited_allele_count"), FALSE)
ref_col <- find_col(sites, c("refcount", "ref_count", "reference_allele_count"), FALSE)
lib_col <- find_col(sites, c("lib", "library", "sample"), FALSE)
sites[, FBgn := clean_id(get(gene_col))]
if (!is.na(ratio_col)) sites[, editing_ratio := suppressWarnings(as.numeric(get(ratio_col)))] else sites[,
    editing_ratio := NA_real_]
if (!is.na(depth_col)) sites[, sequencing_depth := suppressWarnings(as.numeric(get(depth_col)))] else sites[,
    sequencing_depth := NA_real_]
if (!is.na(alt_col)) sites[, edited_allele_count_row := suppressWarnings(as.numeric(get(alt_col)))] else sites[,
    edited_allele_count_row := NA_real_]
if (!is.na(ref_col)) sites[, reference_allele_count_row := suppressWarnings(as.numeric(get(ref_col)))] else sites[,
    reference_allele_count_row := NA_real_]
if (all(is.na(sites$sequencing_depth)) && !all(is.na(sites$edited_allele_count_row)) && !all(is.na(sites$reference_allele_count_row))) sites[,
    sequencing_depth := edited_allele_count_row + reference_allele_count_row]
sites <- sites[grepl("^FBgn", FBgn) & (is.na(sequencing_depth) | sequencing_depth >= editing_min_depth) &
    (is.na(editing_ratio) | editing_ratio >= editing_min_ratio)]
site_coord_cols <- if (all(c("dmel_chr", "dmel_pos") %in% names(sites))) c("dmel_chr", "dmel_pos") else c("chr",
    "pos")
site_library_key <- c("FBgn", site_coord_cols, if (!is.na(lib_col)) lib_col)
# Remove only duplicate annotation rows for the same gene/site/library. Counts
# from distinct libraries are retained and pooled for the gene burden.
setorderv(sites, c(site_library_key, "sequencing_depth"), c(rep(1, length(site_library_key)),
    -1), na.last = TRUE)
site_library <- unique(sites, by = site_library_key)
site_summary <- site_library[, .(site_edited_alleles = sum(edited_allele_count_row, na.rm = TRUE),
    site_total_alleles = sum(sequencing_depth, na.rm = TRUE), fallback_ratio = mean(editing_ratio,
        na.rm = TRUE)), by = c("FBgn", site_coord_cols)]
site_summary[, site_editing_fraction := ifelse(site_total_alleles > 0, site_edited_alleles/site_total_alleles,
    fallback_ratio)]
editing_gene <- site_summary[, .(editing_site_count = .N, edited_allele_count = sum(site_edited_alleles,
    na.rm = TRUE), total_allele_count = sum(site_total_alleles, na.rm = TRUE), mean_editing = mean(site_editing_fraction,
    na.rm = TRUE), median_editing = median(site_editing_fraction, na.rm = TRUE), max_editing = max(site_editing_fraction,
    na.rm = TRUE), sum_editing = sum(site_editing_fraction, na.rm = TRUE), mean_depth = mean(site_total_alleles,
    na.rm = TRUE)), by = FBgn]
editing_gene[, gene_editing_burden := ifelse(total_allele_count > 0, edited_allele_count/total_allele_count,
    NA_real_)]
editing_gene[!is.finite(mean_editing), c("mean_editing", "median_editing", "max_editing",
    "sum_editing", "mean_depth", "gene_editing_burden") := NA_real_]
editing_gene[, edited := as.integer(editing_site_count >= min_editing_sites_per_gene)]
editing_gene[, binary_analysis_status := fifelse(editing_site_count >= min_editing_sites_per_gene,
    "Edited", sprintf("Excluded_below_%d_site_threshold", min_editing_sites_per_gene))]
write_tsv(editing_gene, file.path(out_dir, "01_gene_mapping", "editing_gene_summary.tsv"))
write_tsv(editing_gene[editing_site_count < min_editing_sites_per_gene], file.path(out_dir,
    "00_input_QC", "editing_binary_threshold_excluded_genes.tsv"))
write_tsv(data.table(metric = c("input_rows", "rows_after_filters", "unique_gene_site_library_rows",
    "unique_gene_sites", "genes_with_at_least_1_detected_site", sprintf("genes_with_at_least_%d_detected_sites",
        min_editing_sites_per_gene), sprintf("genes_excluded_with_1_to_%d_sites", min_editing_sites_per_gene -
        1L), "duplicate_annotation_rows_removed"), value = c(nrow(fread(editing_file, select = gene_col)),
    nrow(sites), nrow(site_library), nrow(site_summary), uniqueN(editing_gene$FBgn), sum(editing_gene$editing_site_count >=
        min_editing_sites_per_gene), sum(editing_gene$editing_site_count > 0 & editing_gene$editing_site_count <
        min_editing_sites_per_gene), nrow(sites) - nrow(site_library))), file.path(out_dir,
    "00_input_QC", "editing_QC.tsv"))

# Callable background: a one-column FBgn list means all rows are callable; an
# explicit callable flag is respected.  Without it, the analysis is exploratory.
callable <- data.table(FBgn = character(), callable = logical())
if (!is.null(callable_gene_file) && file.exists(callable_gene_file)) {
    callable <- standardize_gene_column(fread(callable_gene_file), "callable_gene_file")
    ccol <- find_col(callable, c("callable", "is_callable", "eligible"), FALSE)
    if (is.na(ccol))
        callable[, callable := TRUE] else callable[, callable := tolower(as.character(get(ccol))) %in% c("true", "t", "1",
        "yes", "y")]
    callable <- unique(callable[, .(FBgn, callable)])
}
saveRDS(list(gene_features = gene_features, transcript_features = features_tx, cg_mapping = cg_mapping,
    editing_gene = editing_gene, callable = callable), file.path(out_dir, "01_gene_mapping",
    "prepared_annotation_and_editing.rds"))
message("Feature construction complete: ", nrow(gene_features), " representative transcripts; ",
    nrow(editing_gene), " genes with >=1 detected site; ", sum(editing_gene$editing_site_count >=
        min_editing_sites_per_gene), " Edited genes with >=", min_editing_sites_per_gene,
    " sites; callable background supplied: ", nrow(callable) > 0)
