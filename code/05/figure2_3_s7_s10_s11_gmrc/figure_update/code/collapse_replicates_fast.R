# Equivalent implementation of section 11's two expensive aggregations.
# The input data.table is never modified. Group order is first appearance,
# as in raw_all[, ..., by = .(group_id, FBtr_codonpos)].
collapse_replicates_fast <- function(raw_all, expected_design) {
    stopifnot(data.table::is.data.table(raw_all), nrow(raw_all) > 0L)
    keys <- c("group_id", "FBtr_codonpos")
    key_dt <- unique(raw_all[, ..keys])
    gid <- key_dt[raw_all, on = keys, which = TRUE]
    ngroup <- nrow(key_dt)
    first_row <- match(seq_len(ngroup), gid)

    # uniqueN counts NA as a distinct value, including disagreement between
    # NA and a nonmissing annotation. Whole-column unique preserves that rule.
    annotation_map <- c(n_FBtr = "FBtr", n_codonpos = "codonpos", n_codon = "codon", n_AA = "AA",
        n_codon_edited = "codon_edited", n_codonlength = "codonlength", n_gene_id = "gene_id",
        n_structure = "secondary_structure")
    static_qc <- data.table::copy(key_dt)
    for (out_name in names(annotation_map)) {
        values <- raw_all[[annotation_map[[out_name]]]]
        if (out_name == "n_structure")
            values <- as.character(values)
        unique_values <- unique(data.table::data.table(.gid = gid, value = values))
        data.table::set(static_qc, j = out_name, value = tabulate(unique_values$.gid, nbins = ngroup))
    }
    static_bad <- static_qc[n_FBtr > 1L | n_codonpos > 1L | n_codon > 1L | n_AA > 1L | n_codon_edited >
        1L | n_codonlength > 1L | n_gene_id > 1L | n_structure > 1L]
    if (nrow(static_bad)) {
        return(list(static_qc = static_qc, static_bad = static_bad, group_mean = NULL))
    }
    rm(unique_values, values)

    finite_numeric <- function(x) {
        x <- as.numeric(x)
        x[!is.finite(x)] <- NA_real_
        x
    }
    work <- data.table::data.table(.gid = gid, mRNA = finite_numeric(raw_all$mRNA), asite_positive = is.finite(raw_all$Asite_median) &
        raw_all$Asite_median > 0, gmRC = finite_numeric(raw_all$gmRC), gmRC_finite = is.finite(as.numeric(raw_all$gmRC)),
        logTE = finite_numeric(raw_all$logTE), logTE_finite = is.finite(as.numeric(raw_all$logTE)),
        phyloP = finite_numeric(raw_all$codon_phyloP_mean), phyloP_finite = is.finite(as.numeric(raw_all$codon_phyloP_mean)),
        phyloP_covered = finite_numeric(raw_all$codon_phyloP_covered))
    # These simple direct aggregates are optimized by data.table's GForce.
    # No user-defined R callback executes once for every individual codon.
    summaries <- suppressWarnings(work[, .(mRNA_sum = sum(mRNA, na.rm = TRUE), Asite_positive_reps = sum(asite_positive),
        gmRC = mean(gmRC, na.rm = TRUE), n_gmRC_reps_finite = sum(gmRC_finite), logTE = mean(logTE,
            na.rm = TRUE), n_logTE_reps_finite = sum(logTE_finite), codon_phyloP_mean = mean(phyloP,
            na.rm = TRUE), .phyloP_n = sum(phyloP_finite), codon_phyloP_covered = min(phyloP_covered,
            na.rm = TRUE), gmRC_rep_min = min(gmRC, na.rm = TRUE), gmRC_rep_max = max(gmRC,
            na.rm = TRUE)), by = .gid])
    # GForce mean can overflow at extreme finite magnitudes (e.g. 1e308),
    # whereas base::mean handles those values. Preserve original semantics
    # with a base-R fallback only for affected, nonempty groups.
    mean_columns <- c(gmRC = "gmRC", logTE = "logTE", codon_phyloP_mean = "phyloP")
    n_columns <- c(gmRC = "n_gmRC_reps_finite", logTE = "n_logTE_reps_finite", codon_phyloP_mean = ".phyloP_n")
    for (column in names(mean_columns)) {
        bad <- which(!is.finite(summaries[[column]]) & summaries[[n_columns[[column]]]] >
            0L)
        if (length(bad)) {
            source_column <- mean_columns[[column]]
            affected <- summaries$.gid[bad]
            replacement <- work[.gid %in% affected, .(value = base::mean(as.numeric(get(source_column)),
                na.rm = TRUE)), by = .gid]
            data.table::set(summaries, i = bad, j = column, value = replacement$value[match(affected,
                replacement$.gid)])
        }
    }
    summaries[, .phyloP_n := NULL]
    # safe_mean/min/max return NA_real_ for empty finite subsets. In contrast,
    # GForce's corresponding empty values may be NaN or +/-Inf.
    for (column in c("gmRC", "logTE", "codon_phyloP_mean", "codon_phyloP_covered", "gmRC_rep_min",
        "gmRC_rep_max")) {
        v <- summaries[[column]]
        v[!is.finite(v)] <- NA_real_
        data.table::set(summaries, j = column, value = v)
    }
    rm(work)

    # The three biological designs have 2 or 3 replicate IDs. Collapse their
    # presence into a numeric bit set and format each observed combination
    # once. NA counts in uniqueN(replicate), but sort() omits it in the string.
    rep_values <- sort(unique(raw_all$replicate), na.last = TRUE)
    rep_pairs <- unique(data.table::data.table(.gid = gid, replicate = raw_all$replicate))
    if (length(rep_values) <= 52L) {
        rep_pairs[, .rep_bit := 2^(match(replicate, rep_values) - 1L)]
        rep_summary <- rep_pairs[, .(n_reps_present = .N, .mask = sum(.rep_bit)), by = .gid]
        masks <- unique(rep_summary$.mask)
        labels <- vapply(masks, function(mask) {
            included <- (floor(mask/2^(seq_along(rep_values) - 1L))%%2) == 1
            paste(sort(rep_values[included]), collapse = ",")
        }, character(1))
        rep_summary[, replicate_ids_present := labels[match(.mask, masks)]]
    } else {
        # General fallback is deliberately retained for unexpected designs with
        # more than 52 different replicate IDs, beyond exact double bit sets.
        rep_summary <- rep_pairs[, .(n_reps_present = .N, replicate_ids_present = paste(sort(replicate),
            collapse = ",")), by = .gid]
    }
    order_summary <- match(seq_len(ngroup), summaries$.gid)
    order_rep <- match(seq_len(ngroup), rep_summary$.gid)
    group_mean <- data.table::copy(key_dt)
    group_mean[, expected_n_reps := expected_design$expected_n_reps[match(group_id, expected_design$group_id)]]
    group_mean[, `:=`(n_reps_present = rep_summary$n_reps_present[order_rep], replicate_ids_present = rep_summary$replicate_ids_present[order_rep])]
    static_fields <- c("FBtr", "codonpos", "codon", "codon_edited", "AA", "codonlength",
        "gene_id")
    for (column in static_fields) {
        data.table::set(group_mean, j = column, value = raw_all[[column]][first_row])
    }
    group_mean[, secondary_structure_chr := as.character(raw_all$secondary_structure[first_row])]
    for (column in setdiff(names(summaries), ".gid")) {
        data.table::set(group_mean, j = column, value = summaries[[column]][order_summary])
    }
    group_mean[, `:=`(mRNA_mean = mRNA_sum/expected_n_reps, Asite_any_positive = Asite_positive_reps >
        0L)]
    data.table::setcolorder(group_mean, c(keys, "expected_n_reps", "n_reps_present", "replicate_ids_present",
        static_fields, "secondary_structure_chr", "mRNA_sum", "mRNA_mean", "Asite_any_positive",
        "Asite_positive_reps", "gmRC", "n_gmRC_reps_finite", "logTE", "n_logTE_reps_finite",
        "codon_phyloP_mean", "codon_phyloP_covered", "gmRC_rep_min", "gmRC_rep_max"))
    list(static_qc = static_qc, static_bad = static_bad, group_mean = group_mean)
}
