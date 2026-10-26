.libPaths(c("/lustre/user/lulab/lulabrotation/R/x86_64-conda-linux-gnu-library/4.1", .libPaths()))

suppressPackageStartupMessages(library(data.table))
message("[fallback-alledit] script started")

# Recalculate one A-site alledit file.  Candidate sites are reported in the
# selected (new/longest-CDS) coordinate system.  If that transcript has no
# RNA coverage at the candidate codon in this sample, counts are instead
# taken from the source-transcript CDS coordinate retained by the cdspos file.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
    stop("Usage: Rscript recompute_alledit_with_original_tx_fallback.R <rna.num> <asite.num> <cdspos.tsv> [output_prefix]")
}
rna_file <- args[[1L]]
ribo_file <- args[[2L]]
edit_file <- args[[3L]]
output_prefix <- if (length(args) >= 4L) args[[4L]] else sub("\\.[^.]+$", "", ribo_file)

edit <- fread(edit_file, header = TRUE, stringsAsFactors = FALSE)
needed <- c("genomic_key", "target_region", "new_transcript", "new_cds_pos_1based", "source_transcript",
    "source_cds_pos_1based")
missing <- setdiff(needed, names(edit))
if (length(missing) > 0L) stop("cdspos file is missing: ", paste(missing, collapse = ", "))

edit <- edit[target_region == "CDS"]
edit[, `:=`(new_transcript = as.character(new_transcript), source_transcript = as.character(source_transcript),
    new_cds_pos = suppressWarnings(as.integer(new_cds_pos_1based)), source_cds_pos = suppressWarnings(as.integer(source_cds_pos_1based)))]
edit <- edit[!is.na(new_transcript) & new_transcript != "" & new_transcript != "NA" & !is.na(new_cds_pos) &
    new_cds_pos >= 1L]
edit[, `:=`(new_handle = paste(new_transcript, new_cds_pos, sep = "_"), new_codon_start = new_cds_pos -
    ((new_cds_pos - 1L)%%3L), target_handle = paste(new_transcript, new_cds_pos - ((new_cds_pos -
    1L)%%3L), sep = "_"), source_handle = fifelse(!is.na(source_transcript) & source_transcript !=
    "" & source_transcript != "NA" & !is.na(source_cds_pos) & source_cds_pos >= 1L, paste(source_transcript,
    source_cds_pos, sep = "_"), NA_character_))]

# The raw count files are large.  Filter them by the candidate transcripts
# while streaming, before they enter R memory.  Both coordinate systems are
# kept here because source coordinates may be selected per sample below.
candidate_tx_pre <- unique(c(edit$new_transcript, edit$source_transcript[!is.na(edit$source_handle)]))
candidate_tx_pre <- candidate_tx_pre[!is.na(candidate_tx_pre) & candidate_tx_pre != ""]
candidate_count_patterns <- unique(c(paste(edit$new_transcript, edit$new_codon_start, sep = "\t"),
    paste(edit$source_transcript[!is.na(edit$source_handle)], edit$source_cds_pos[!is.na(edit$source_handle)] -
        ((edit$source_cds_pos[!is.na(edit$source_handle)] - 1L)%%3L), sep = "\t")))
candidate_count_patterns <- candidate_count_patterns[!is.na(candidate_count_patterns) & !grepl("^NA\\t|\\tNA$",
    candidate_count_patterns)]
read_selected_counts <- function(file, value_col) {
    pattern_file <- tempfile(pattern = "candidate_transcripts_", fileext = ".txt")
    on.exit(unlink(pattern_file), add = TRUE)
    writeLines(paste0(candidate_count_patterns, "\t"), pattern_file)
    # Only candidate codons are needed.  This keeps a 100+ MB raw count table
    # to a few thousand rows before it enters R memory.
    cmd <- sprintf("grep -F -f %s %s", shQuote(pattern_file), shQuote(file))
    x <- fread(cmd = cmd, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) != 4L)
        stop("No usable candidate-transcript rows in ", file)
    setnames(x, c("FBtr", "CDSpos_codonstart", "codon", value_col))
    x[, CDSpos_codonstart := as.integer(CDSpos_codonstart)]
    set(x, j = value_col, value = as.numeric(x[[value_col]]))
    x[, handle := paste(FBtr, CDSpos_codonstart, sep = "_")]
    x[, handle2 := paste(handle, codon, sep = "_")]
    x
}
rna <- read_selected_counts(rna_file, "number_rna")
ribo <- read_selected_counts(ribo_file, "number_ribo")

# A transcript is considered covered only when its RNA count table contains a
# positive-count codon at the candidate codon start.  This avoids replacing a
# genuine zero-editing observation while rescuing transcript-ID mismatches.
covered_codons <- unique(rna[is.finite(number_rna) & number_rna > 0 & CDSpos_codonstart >
    0, .(FBtr, CDSpos_codonstart)])
setkey(covered_codons, FBtr, CDSpos_codonstart)
has_rna_coverage <- function(tx, pos) {
    if (is.na(tx) || tx == "" || is.na(pos) || pos < 1L)
        return(FALSE)
    codon_start <- pos - ((pos - 1L)%%3L)
    nrow(covered_codons[.(tx, codon_start), nomatch = 0L]) > 0L
}
edit[, new_tx_covered := mapply(has_rna_coverage, new_transcript, new_cds_pos)]

new_map <- edit[, .(raw_FBt_cdspos = new_handle, target_FBt_cdspos = new_handle, target_handle,
    count_match_route = "new_transcript")]
source_map <- edit[!new_tx_covered & !is.na(source_handle), .(raw_FBt_cdspos = source_handle,
    target_FBt_cdspos = new_handle, target_handle, count_match_route = "source_transcript_fallback")]
site_map <- rbindlist(list(new_map, source_map), use.names = TRUE)
site_map <- site_map[!duplicated(raw_FBt_cdspos)]

if (nrow(site_map) == 0L) stop("No valid CDS candidate sites after coordinate preparation")
candidate_tx <- unique(sub("_[0-9]+$", "", site_map$raw_FBt_cdspos))

dt <- merge(rna, ribo[, .(handle2, number_ribo)], by = "handle2", all = TRUE)
dt[is.na(dt)] <- 0
dt <- dt[CDSpos_codonstart > 0L & number_rna >= 1 & FBtr %in% candidate_tx]

compare_codons_group <- function(group) {
    if (nrow(group) < 2L)
        return(NULL)
    ans <- list()
    combs <- combn(nrow(group), 2L)
    for (i in seq_len(ncol(combs))) {
        a <- group[combs[1L, i]]
        b <- group[combs[2L, i]]
        aa <- strsplit(a$codon, "")[[1L]]
        bb <- strsplit(b$codon, "")[[1L]]
        if (length(aa) != length(bb))
            next
        dif <- which(aa != bb)
        if (!length(dif) || !all(aa[dif] %in% c("A", "G")) || !all(bb[dif] %in% c("A", "G")))
            next
        if ("A" %in% aa[dif] && "G" %in% bb[dif]) {
            ref <- a
            alt <- b
        } else if ("G" %in% aa[dif] && "A" %in% bb[dif]) {
            ref <- b
            alt <- a
        } else next
        ans[[length(ans) + 1L]] <- data.table(ref_codon = ref$codon, alt_codon = alt$codon,
            ref_number_rna = ref$number_rna, alt_number_rna = alt$number_rna, ref_number_ribo = ref$number_ribo,
            alt_number_ribo = alt$number_ribo)
    }
    if (!length(ans))
        return(NULL)
    rbindlist(ans)
}

setDT(dt)
results <- dt[, compare_codons_group(.SD), by = handle]

empty_result <- function() {
    data.table(handle = character(), ref_codon = character(), alt_codon = character(), ref_number_rna = numeric(),
        alt_number_rna = numeric(), ref_number_ribo = numeric(), alt_number_ribo = numeric(),
        edit_pos = character(), FBt_cdspos = character(), raw_FBt_cdspos = character(), count_match_route = character(),
        alt_ratio = numeric(), ref_ratio = numeric(), gmRC = numeric(), diff = numeric())
}

if (is.null(results) || !nrow(results) || !all(c("ref_codon", "alt_codon") %in% names(results))) {
    result <- empty_result()
} else {
    edit_pos <- function(ref_codon, alt_codon) {
        rr <- strsplit(ref_codon, "")[[1L]]
        aa <- strsplit(alt_codon, "")[[1L]]
        d <- which(rr != aa & rr == "A" & aa == "G")
        if (!length(d))
            "0" else paste(d, collapse = ", ")
    }
    results[, edit_pos := mapply(edit_pos, ref_codon, alt_codon)]
    results[, codon_start := as.integer(sub(".*_", "", handle))]
    results[, FBtr_id := sub("_[0-9]+$", "", handle)]
    results[, raw_FBt_cdspos := fifelse(edit_pos == "1", paste0(FBtr_id, "_", codon_start),
        fifelse(edit_pos == "2", paste0(FBtr_id, "_", codon_start + 1L), fifelse(edit_pos ==
            "3", paste0(FBtr_id, "_", codon_start + 2L), NA_character_)))]
    results[, c("codon_start", "FBtr_id") := NULL]
    idx <- match(results$raw_FBt_cdspos, site_map$raw_FBt_cdspos)
    keep <- !is.na(idx)
    if (!any(keep)) {
        result <- empty_result()
    } else {
        result <- copy(results[keep])
        idx <- idx[keep]
        result[, `:=`(FBt_cdspos = site_map$target_FBt_cdspos[idx], count_match_route = site_map$count_match_route[idx],
            raw_handle = handle, handle = site_map$target_handle[idx])]
        result[, alt_ratio := (alt_number_ribo + 1)/(alt_number_rna + 1)]
        result[, ref_ratio := (ref_number_ribo + 1)/(ref_number_rna + 1)]
        result[, gmRC := alt_ratio/ref_ratio]
        result[, diff := alt_ratio - ref_ratio]
    }
}

result1 <- result[(ref_number_rna + alt_number_rna) >= 1 | (ref_number_ribo + alt_number_ribo) >=
    1]
output_file <- paste0(output_prefix, "_df_alledit")
fwrite(result1, file = output_file, sep = "\t", quote = FALSE)
message("Matched editing codon pairs: ", nrow(result1))
message("Output written: ", output_file)
if (nrow(result1)) print(result1[, .N, by = count_match_route])
