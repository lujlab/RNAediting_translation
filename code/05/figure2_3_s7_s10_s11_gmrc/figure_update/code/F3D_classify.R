#!/usr/bin/env Rscript

# Classify the species-sex mean-gmRC codons into four AA groups.
#
# This script consumes the tables made by
# 00_build_species_sex_mean_gmRC_tables.R. It does not calculate a second
# sample-level mean. Every retained row remains one FBtr_codonpos and therefore
# becomes one point in the detailed violin plots.

suppressPackageStartupMessages({
    library(data.table)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
    getwd()
}
args <- commandArgs(trailingOnly = TRUE)

mean_table_dir <- if (length(args) >= 1L && nzchar(args[[1L]])) {
    args[[1L]]
} else {
    file.path(script_dir, "species_sex_mean_gmRC_tables")
}
output_dir <- if (length(args) >= 2L && nzchar(args[[2L]])) {
    args[[2L]]
} else {
    file.path(script_dir, "AA_gmRC_species_sex_syn_nonsyn_results")
}
edited_reference_dir <- if (length(args) >= 3L && nzchar(args[[3L]])) {
    args[[3L]]
} else {
    file.path(script_dir, "rebuilt_edited_references_3R", "split")
}
codon_table_file <- if (length(args) >= 4L && nzchar(args[[4L]])) {
    args[[4L]]
} else {
    "/path/to/RNAi/240719drosophila_check/Codon_table.new.txt"
}
rscu_file <- if (length(args) >= 5L && nzchar(args[[5L]])) {
    args[[5L]]
} else {
    "/path/to/RNAi/240719drosophila_check/RSCU_flybase.txt"
}
requested_units <- if (length(args) >= 6L && nzchar(args[[6L]])) {
    unique(trimws(strsplit(args[[6L]], ",", fixed = TRUE)[[1L]]))
} else {
    character()
}

unedit_A_only <- tolower(Sys.getenv("AA_UNEDIT_A_ONLY", unset = "false")) %in% c("1", "true",
    "yes", "y")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

m6a_reference_config <- Sys.getenv("M6A_REFERENCE_CONFIG", "")
if (!nzchar(m6a_reference_config) || !file.exists(m6a_reference_config)) {
    stop("Missing M6A_REFERENCE_CONFIG")
}
m6a_reference <- fread(m6a_reference_config)
stopifnot(all(c("species", "m6a_file") %in% names(m6a_reference)), !anyDuplicated(m6a_reference$species))
m6a_by_species <- setNames(lapply(m6a_reference$m6a_file, function(path) {
    if (!file.exists(path))
        stop("Missing m6A file: ", path)
    unique(scan(path, what = character(), quiet = TRUE))
}), m6a_reference$species)

min_edited_codons <- 5L
flank_codons <- 5L
group_levels <- c("Unedited_preferred", "Unedited_unpreferred", "Edited")
amino_acid_order <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

reference_files <- c(Dmel = "dmel-all-CDS-r6.12.fasta.edited.tab.split", Dsim = "dsim-all-CDS-r2.02.fasta.edited.tab.split",
    Dsec = "dsec-all-CDS-r1.3.fasta.edited.tab.split", Dere = "dere-all-CDS-r1.05.fasta.edited.tab.split")

manifest_file <- file.path(mean_table_dir, "mean_gmRC_table_manifest.tsv")
if (file.exists(manifest_file)) {
    manifest <- fread(manifest_file)
} else if (length(requested_units) > 0L) {
    # Permit a partial preview (for example, Dmel female/male while the full
    # nine-unit merge is still running) without fabricating a full manifest.
    preview_files <- file.path(mean_table_dir, paste0(requested_units, ".mean_gmRC.tsv.gz"))
    if (any(!file.exists(preview_files))) {
        stop("Missing requested preview mean-gmRC table(s): ", paste(preview_files[!file.exists(preview_files)],
            collapse = ", "))
    }
    preview_species <- sub("__.*$", "", requested_units)
    preview_sex <- sub("^.*__", "", requested_units)
    manifest <- data.table(species = preview_species, cohort = fifelse(preview_species ==
        "Dmel", "Dmel_WBg_ISO1", preview_species), sex = preview_sex, unit_id = requested_units,
        analysis_role = "main_8", n_input_replicates = NA_integer_, input_samples = NA_character_,
        mean_gmRC_file = normalizePath(preview_files, mustWork = TRUE))
} else {
    stop("Missing mean-gmRC manifest: ", manifest_file, "\nRun 00_build_species_sex_mean_gmRC_tables.R first.")
}
if (!file.exists(codon_table_file)) stop("Missing codon table: ", codon_table_file)
if (!file.exists(rscu_file)) stop("Missing RSCU table: ", rscu_file)

required_manifest_columns <- c("species", "cohort", "sex", "unit_id", "analysis_role", "n_input_replicates",
    "input_samples", "mean_gmRC_file")
missing_manifest_columns <- setdiff(required_manifest_columns, names(manifest))
if (length(missing_manifest_columns) > 0L) {
    stop("Missing manifest columns: ", paste(missing_manifest_columns, collapse = ", "))
}
if (length(requested_units) > 0L) {
    unknown_units <- setdiff(requested_units, manifest$unit_id)
    if (length(unknown_units) > 0L) {
        stop("Requested units are absent from manifest: ", paste(unknown_units, collapse = ", "))
    }
    manifest <- manifest[unit_id %chin% requested_units]
}
expected_main_units <- if (length(requested_units) > 0L) {
    uniqueN(manifest[analysis_role == "main_8", unit_id])
} else {
    8L
}
if (manifest[analysis_role == "main_8", uniqueN(unit_id)] != expected_main_units || expected_main_units ==
    0L) {
    stop("The selected manifest contains no main species-sex units")
}
missing_mean_tables <- manifest[!file.exists(mean_gmRC_file), mean_gmRC_file]
if (length(missing_mean_tables) > 0L) {
    stop("Missing mean-gmRC tables:\n", paste(missing_mean_tables, collapse = "\n"))
}

aa_table <- fread(codon_table_file, header = TRUE)
if (!all(c("Codon", "No") %in% names(aa_table))) {
    stop("Codon table must contain Codon and No columns")
}
aa_map <- unique(aa_table[, .(codon_key = toupper(trimws(as.character(Codon))), AA_translated = as.character(No))])

rscu <- fread(rscu_file, header = TRUE)
if (!all(c("codon", "category") %in% names(rscu))) {
    stop("RSCU table must contain codon and category columns")
}
preference_map <- unique(rscu[, .(codon_key = toupper(trimws(as.character(codon))), preference_category = tolower(trimws(as.character(category))))])

is_valid_a_to_g_change <- function(before, after) {
    before <- toupper(as.character(before))
    after <- toupper(as.character(after))
    valid_triplet <- !is.na(before) & !is.na(after) & nchar(before) == 3L & nchar(after) ==
        3L
    before_bases <- lapply(1:3, function(i) substr(before, i, i))
    after_bases <- lapply(1:3, function(i) substr(after, i, i))
    allowed <- Reduce(`&`, lapply(1:3, function(i) {
        before_bases[[i]] == after_bases[[i]] | (before_bases[[i]] == "A" & after_bases[[i]] ==
            "G")
    }))
    valid_triplet & before != after & allowed
}

edited_reference_cache <- new.env(parent = emptyenv())

read_edited_reference <- function(species_name) {
    if (exists(species_name, envir = edited_reference_cache, inherits = FALSE)) {
        return(get(species_name, envir = edited_reference_cache, inherits = FALSE))
    }
    filename <- reference_files[[species_name]]
    if (is.null(filename))
        stop("No edited reference configured for ", species_name)
    path <- file.path(edited_reference_dir, filename)
    if (!file.exists(path))
        stop("Missing edited reference: ", path)

    x <- fread(path, header = FALSE, select = 1:2, col.names = c("FBtr_codonpos", "codon_edited"))
    x[, `:=`(FBtr_codonpos = as.character(FBtr_codonpos), codon_edited = toupper(trimws(as.character(codon_edited))))]
    if (x[, anyDuplicated(FBtr_codonpos)] > 0L) {
        stop("Edited reference has duplicated FBtr_codonpos values: ", path)
    }
    assign(species_name, x, envir = edited_reference_cache)
    x[]
}

annotated_rows <- vector("list", nrow(manifest))
unit_qc_rows <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    info <- manifest[i]
    species_name <- info$species[[1L]]
    message("Classifying ", info$unit_id[[1L]])

    x <- fread(info$mean_gmRC_file[[1L]])
    # Input adapter: eligibility is already established by the new gmRC rules.
    # These metadata aliases do not introduce older RNA/RPF thresholds; the explicit species m6A filter is applied below.
    x[, `:=`(CDS_length = cds_coding_codons + 1L, codonpos_from_end = cds_coding_codons -
        codonpos + 1L, handle = paste(chr, pos, sep = "_"), m6A = NA_character_, m6A_mod = NA_character_,
        RNA_count = mRNA_sum, RPF_count = Asite_sum, RNA_count_threshold = NA_real_, RPF_count_threshold = NA_real_,
        valid_codon = is.finite(gmRC), n_replicates_total = n_libraries_expected, n_replicates_with_gmRC = n_gmRC_reps_finite,
        replicates_with_gmRC = sample_ids_finite, source_samples = info$input_samples[[1L]])]

    n_input_codons_before_m6a <- nrow(x)
    is_m6a <- !is.na(x$handle) & x$handle %chin% m6a_by_species[[species_name]]
    n_m6a_codons_removed <- sum(is_m6a)
    x <- x[!is_m6a]

    required_input_columns <- c("FBtr", "FBtr_codonpos", "codon", "AA", "gmRC", "chr", "pos",
        "codonpos", "CDS_length", "codonpos_from_end", "handle", "m6A", "m6A_mod", "RNA_count",
        "RPF_count", "RNA_count_threshold", "RPF_count_threshold", "valid_codon", "n_replicates_total",
        "n_replicates_with_gmRC", "replicates_with_gmRC", "source_samples")
    missing_input_columns <- setdiff(required_input_columns, names(x))
    if (length(missing_input_columns) > 0L) {
        stop("Missing columns in ", info$mean_gmRC_file[[1L]], ": ", paste(missing_input_columns,
            collapse = ", "))
    }
    if (any(x$valid_codon != TRUE, na.rm = TRUE) || anyNA(x$valid_codon)) {
        stop("The mean-gmRC table contains invalid codons in ", info$mean_gmRC_file[[1L]],
            "; rerun 00_build_species_sex_mean_gmRC_tables.R so filtering occurs ", "after replicate-count aggregation.")
    }

    setnames(x, c("codon", "AA"), c("codon_original", "AA_original_input"))
    x[, `:=`(FBtr = as.character(FBtr), FBtr_codonpos = as.character(FBtr_codonpos), codon_original = toupper(trimws(as.character(codon_original))),
        gmRC = suppressWarnings(as.numeric(gmRC)))]

    edited_reference <- read_edited_reference(species_name)
    x <- merge(x, edited_reference, by = "FBtr_codonpos", all.x = TRUE, sort = FALSE)
    x <- merge(x, aa_map[, .(codon_original = codon_key, AA_original = AA_translated)], by = "codon_original",
        all.x = TRUE, sort = FALSE)
    x <- merge(x, aa_map[, .(codon_edited = codon_key, AA_edited = AA_translated)], by = "codon_edited",
        all.x = TRUE, sort = FALSE)
    x <- merge(x, preference_map[, .(codon_edited = codon_key, preference_category)], by = "codon_edited",
        all.x = TRUE, sort = FALSE)

    x[, `:=`(species = info$species[[1L]], cohort = info$cohort[[1L]], sex = info$sex[[1L]],
        unit_id = info$unit_id[[1L]], analysis_role = info$analysis_role[[1L]], has_reference_change = !is.na(codon_original) &
            !is.na(codon_edited) & codon_original != codon_edited, valid_a_to_g_change = is_valid_a_to_g_change(codon_original,
            codon_edited))]

    # The old per-sample long tables had already removed unedited codons within
    # +/-5 codons of every edit. Because this version starts from the merged raw
    # codonTE tables, rebuild that exclusion explicitly so the biological group
    # definition is unchanged.
    canonical_flank_file <- file.path(edited_reference_dir, paste0(species_name, ".edited_codon_pm5_keys.tsv"))
    if (!file.exists(canonical_flank_file))
        stop("Missing canonical editing neighborhoods: ", canonical_flank_file)
    edit_neighborhood <- unique(fread(canonical_flank_file, select = "FBtr_codonpos")$FBtr_codonpos)
    x[, `:=`(is_edit_site = valid_a_to_g_change, is_edit_neighborhood = FBtr_codonpos %chin%
        edit_neighborhood)]

    # AA_original is used only to decide synonymous versus nonsynonymous.
    # The plotted/faceted amino acid and the preference class are both assigned
    # from the edited-reference codon, never from the pre-editing AA column.
    x[, group := fcase(!has_reference_change & !is_edit_neighborhood & preference_category ==
        "optimal", "Unedited_preferred", !has_reference_change & !is_edit_neighborhood &
        preference_category == "nonoptimal", "Unedited_unpreferred", has_reference_change &
        valid_a_to_g_change, "Edited", default = NA_character_)]
    x[, AA := AA_edited]

    unit_qc_rows[[i]] <- x[, .(species = first(species), cohort = first(cohort), sex = first(sex),
        unit_id = first(unit_id), analysis_role = first(analysis_role), n_input_codons_before_m6a = n_input_codons_before_m6a,
        n_m6a_codons_removed = n_m6a_codons_removed, n_input_codons = .N, n_finite_nonnegative_gmRC = sum(is.finite(gmRC) &
            gmRC >= 0), n_edited_reference_matched = sum(!is.na(codon_edited)), n_valid_A_to_G_changed_codons = sum(valid_a_to_g_change),
        n_invalid_changed_codons = sum(has_reference_change & !valid_a_to_g_change), n_edit_neighborhood_codons = sum(is_edit_neighborhood),
        n_AA_original_input_disagrees = sum(!is.na(AA_original_input) & !is.na(AA_original) &
            AA_original_input != AA_original), n_four_group_codons_before_AA_filter = sum(!is.na(group)))]

    annotated_rows[[i]] <- x[!is.na(group) & AA %chin% amino_acid_order & is.finite(gmRC) &
        gmRC >= 0]
}

annotated <- rbindlist(annotated_rows, use.names = TRUE, fill = TRUE)
if (nrow(annotated) == 0L) stop("No codons remained after annotation")

if (unedit_A_only) {
    # Keep every edited codon unchanged. For the two unedited groups, retain
    # only reference codons containing A; for an unedited row codon_original is
    # identical to codon_edited.
    annotated <- annotated[group == "Edited" | (group %chin% c("Unedited_preferred", "Unedited_unpreferred") &
        !is.na(codon_original) & grepl("A", codon_original, fixed = TRUE))]
    if (nrow(annotated) == 0L) {
        stop("No codons remained after the unedited-A-codon filter")
    }
}

# Preserve the established AA filter, now applied independently within every
# species-sex unit. Both unedited preference groups must exist and at least
# five edited FBtr_codonpos values must exist in total. One edited subtype is
# allowed to be absent.
group_counts <- annotated[, .(n_codons = uniqueN(FBtr_codonpos)), by = .(species, cohort,
    sex, unit_id, analysis_role, AA, group)]
group_counts_wide <- dcast(group_counts, species + cohort + sex + unit_id + analysis_role +
    AA ~ group, value.var = "n_codons", fill = 0L)
for (group_name in group_levels) {
    if (!group_name %in% names(group_counts_wide)) {
        group_counts_wide[, (group_name) := 0L]
    }
}
group_counts_wide[, n_edited_total := Edited]
group_counts_wide[, keep_AA := Unedited_preferred > 0L & Unedited_unpreferred > 0L & n_edited_total >=
    min_edited_codons]

keep_aa <- group_counts_wide[keep_AA == TRUE, .(unit_id, AA)]
codon_long <- merge(annotated, keep_aa, by = c("unit_id", "AA"), all = FALSE)
if (nrow(codon_long) == 0L) {
    stop("No unit-AA combination passed the AA filter")
}
codon_long[, `:=`(AA = factor(AA, levels = amino_acid_order), group = factor(group, levels = group_levels))]

unit_order <- c("Dmel__female", "Dmel__male", "Dsim__female", "Dsim__male", "Dsec__female",
    "Dsec__male", "Dere__female", "Dere__male", "Dmel_Adar5G1__male")
codon_long[, unit_order_value := match(unit_id, unit_order)]
setorder(codon_long, unit_order_value, AA, group, FBtr, FBtr_codonpos)
codon_long[, unit_order_value := NULL]

unit_medians <- codon_long[, .(median_gmRC = median(gmRC), n_codons = uniqueN(FBtr_codonpos),
    min_gmRC = min(gmRC), q1_gmRC = unname(quantile(gmRC, 0.25)), q3_gmRC = unname(quantile(gmRC,
        0.75)), max_gmRC = max(gmRC)), by = .(species, cohort, sex, unit_id, analysis_role,
    AA, group)]
unit_medians[, unit_order_value := match(unit_id, unit_order)]
setorder(unit_medians, unit_order_value, AA, group)
unit_medians[, unit_order_value := NULL]

unit_qc <- rbindlist(unit_qc_rows, use.names = TRUE)
unit_qc <- group_counts_wide[, .(n_AA_tested = .N, n_AA_retained = sum(keep_AA)), by = unit_id][unit_qc,
    on = "unit_id"]

codon_output <- file.path(output_dir, "AA_four_group_species_sex_codon_long.tsv.gz")
median_output <- file.path(output_dir, "AA_four_group_species_sex_unit_medians.tsv")

fwrite(codon_long, codon_output, sep = "\t", na = "NA", compress = "gzip", compressLevel = 1L)
fwrite(unit_medians, median_output, sep = "\t", na = "NA")
fwrite(group_counts_wide, file.path(output_dir, "AA_four_group_species_sex_filter_QC.tsv"),
    sep = "\t", na = "NA")
fwrite(unit_qc, file.path(output_dir, "AA_four_group_species_sex_unit_QC.tsv"), sep = "\t",
    na = "NA")
fwrite(data.table(parameter = c("unedit_background", "m6A_filter"), value = c(if (unedit_A_only) "only_reference_codons_containing_A" else "all_reference_codons",
    "exclude_exact_species_genomic_m6A_handles_before_classification")), file.path(output_dir,
    "analysis_parameters.tsv"), sep = "\t", na = "NA")

message("Done: ", codon_output)
message("Done: ", median_output)
