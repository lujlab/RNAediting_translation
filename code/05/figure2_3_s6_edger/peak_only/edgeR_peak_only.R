#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Missing R package 'data.table'.", call. = FALSE)
    }
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop(paste0("Missing Bioconductor package 'edgeR'. Install it once with:\n", "  Rscript -e 'if (!requireNamespace(\"BiocManager\", quietly=TRUE)) ",
            "install.packages(\"BiocManager\"); BiocManager::install(\"edgeR\")'"), call. = FALSE)
    }
})

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
edgeR_script_path <- if (length(script_arg)) {
    normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
    NA_character_
}
edgeR_script_md5 <- if (!is.na(edgeR_script_path)) {
    unname(tools::md5sum(edgeR_script_path)[[1L]])
} else {
    NA_character_
}
analysis_schema_version <- "edgeR_peak_only_CDSge150_terminal15_5_m6a_20261002"

parameter_value <- function(x) {
    if (is.null(x) || !length(x) || is.na(x[[1L]]))
        return("<NA>")
    if (is.numeric(x)) {
        return(format(x[[1L]], digits = 17L, scientific = FALSE, trim = TRUE))
    }
    as.character(x[[1L]])
}

make_parameter_table <- function(values) {
    answer <- data.table(parameter = names(values), value = vapply(values, parameter_value,
        character(1L)))
    setorder(answer, parameter)
    answer
}

parameter_differences <- function(previous, current) {
    if (!identical(names(previous), c("parameter", "value"))) {
        stop("Invalid parameter file columns; expected: parameter, value", call. = FALSE)
    }
    if (anyDuplicated(previous$parameter)) {
        stop("Duplicate keys in existing parameter file", call. = FALSE)
    }
    old_values <- setNames(previous$value, previous$parameter)
    new_values <- setNames(current$value, current$parameter)
    keys <- union(names(old_values), names(new_values))
    keys[vapply(keys, function(key) {
        old <- old_values[[key]]
        new <- new_values[[key]]
        is.null(old) || is.null(new) || !identical(old, new)
    }, logical(1L))]
}

usage <- function() {
    cat(paste0("Usage: edgeR_peak_only.R [options]\n\n", "Options:\n", "  --master-manifest FILE   Master table manifest [stalling_master_tables/master_table_manifest.tsv]\n",
        "  --output-dir DIR         Output directory [edgeR_stalling_results]\n", "  --rna-scale-manifest FILE  TSV with sample_id and rna_scale columns [required]\n",
        "  --focal-reference-manifest FILE  Species CDS lengths and mapped m6A references [required]\n",
        "  --groups CSV             Optional comma-separated group IDs\n", "  --flank-sizes CSV        Flank sizes to analyze [20,30]\n",
        "  --min-replicates INT     Minimum biological replicates [2]\n", "  --exclude-head-replicates  Exclude manifest head replicates [included by default]\n",
        "  --min-flank-n-per-rep INT    Minimum flank_n exposure per replicate [10]\n", "  --min-focal-rna-per-rep INT  Minimum focal RNA count per replicate [5]\n",
        "  --min-flank-rna-per-rep INT  Minimum flank RNA count per replicate [100]\n", "  --min-flank-n-sum INT    Optional fixed total flank_n override [disabled]\n",
        "  --min-focal-rna-sum INT  Optional fixed total focal RNA override [disabled]\n",
        "  --min-flank-rna-sum INT  Optional fixed total flank RNA override [disabled]\n",
        "  --min-focal-rpf-sum INT  Minimum focal RPF summed over replicates [0]\n", "  --min-transcript-median-rna FLOAT  Transcript focal-RNA median threshold [0; applied in all replicates]\n",
        "  --min-transcript-median-rpf FLOAT  Transcript focal-RPF median threshold [0; applied in all replicates]\n",
        "  --transcript-median-rule TEXT  all or any_both; any_both requires one replicate with both medians above threshold [any_both]\n",
        "  --normalization-method TEXT  TMM within RNA and RPF separately, or none [TMM]\n",
        "  --fdr FLOAT              One-sided BH FDR cutoff [0.05]\n", "  --min-stalling-log2fc FLOAT  Minimum positive stalling log2 fold change [0]\n",
        "  --sample-fraction FLOAT  Random fraction of eligible sites to fit [1]\n", "  --sample-seed INT        Reproducible sampling seed [20260805]\n",
        "  --run-label TEXT         Optional filename label, e.g. pilot10pct\n", "  --summary-file FILE      Per-run summary TSV path\n",
        "  --threads INT            data.table threads [4]\n", "  --overwrite              Replace existing results\n",
        "  --help                   Print this message\n"))
}

parse_args <- function(args) {
    out <- list(master_manifest = "stalling_master_tables/master_table_manifest.tsv", output_dir = "edgeR_stalling_results",
        rna_scale_manifest = NULL, focal_reference_manifest = NULL, groups = NULL, flank_sizes = "20,30",
        min_replicates = 2L, exclude_head_replicates = FALSE, min_flank_n_per_rep = 10L,
        min_focal_rna_per_rep = 5L, min_flank_rna_per_rep = 100L, min_flank_n_sum = NA_integer_,
        min_focal_rna_sum = NA_integer_, min_flank_rna_sum = NA_integer_, min_focal_rpf_sum = 0L,
        min_transcript_median_rna = 0, min_transcript_median_rpf = 0, transcript_median_rule = "any_both",
        normalization_method = "TMM", fdr = 0.05, min_stalling_log2fc = 0, sample_fraction = 1,
        sample_seed = 20260805L, run_label = NULL, summary_file = NULL, threads = 4L, overwrite = FALSE)
    value_options <- c("--master-manifest", "--output-dir", "--rna-scale-manifest", "--focal-reference-manifest",
        "--groups", "--flank-sizes", "--min-replicates", "--min-flank-n-per-rep", "--min-focal-rna-per-rep",
        "--min-flank-rna-per-rep", "--min-flank-n-sum", "--min-focal-rna-sum", "--min-flank-rna-sum",
        "--min-focal-rpf-sum", "--min-transcript-median-rna", "--min-transcript-median-rpf",
        "--transcript-median-rule", "--normalization-method", "--fdr", "--min-stalling-log2fc",
        "--threads", "--sample-fraction", "--sample-seed", "--run-label", "--summary-file")
    i <- 1L
    while (i <= length(args)) {
        arg <- args[[i]]
        if (arg == "--help") {
            usage()
            quit(save = "no", status = 0L)
        } else if (arg == "--exclude-head-replicates") {
            out$exclude_head_replicates <- TRUE
            i <- i + 1L
            next
        } else if (arg == "--overwrite") {
            out$overwrite <- TRUE
            i <- i + 1L
            next
        } else if (arg %in% value_options) {
            if (i == length(args))
                stop("Missing value after ", arg, call. = FALSE)
            value <- args[[i + 1L]]
            key <- gsub("-", "_", substring(arg, 3L), fixed = TRUE)
            out[[key]] <- value
            i <- i + 2L
            next
        } else {
            stop("Unknown option: ", arg, call. = FALSE)
        }
    }
    if (!is.null(out$groups))
        out$groups <- strsplit(out$groups, ",", fixed = TRUE)[[1L]]
    out$flank_sizes <- as.integer(strsplit(as.character(out$flank_sizes), ",", fixed = TRUE)[[1L]])
    integer_keys <- c("min_replicates", "min_flank_n_per_rep", "min_focal_rna_per_rep", "min_flank_rna_per_rep",
        "min_focal_rpf_sum", "sample_seed", "threads")
    optional_integer_keys <- c("min_flank_n_sum", "min_focal_rna_sum", "min_flank_rna_sum")
    for (key in integer_keys) out[[key]] <- as.integer(out[[key]])
    for (key in optional_integer_keys) out[[key]] <- as.integer(out[[key]])
    out$min_transcript_median_rna <- as.numeric(out$min_transcript_median_rna)
    out$min_transcript_median_rpf <- as.numeric(out$min_transcript_median_rpf)
    if (!out$transcript_median_rule %in% c("all", "any_both")) {
        stop("--transcript-median-rule must be 'all' or 'any_both'", call. = FALSE)
    }
    out$normalization_method <- toupper(out$normalization_method)
    if (!out$normalization_method %in% c("TMM", "NONE")) {
        stop("--normalization-method must be TMM or none", call. = FALSE)
    }
    out$fdr <- as.numeric(out$fdr)
    out$min_stalling_log2fc <- as.numeric(out$min_stalling_log2fc)
    out$sample_fraction <- as.numeric(out$sample_fraction)
    if (anyNA(out$flank_sizes) || any(out$flank_sizes < 1L))
        stop("Invalid --flank-sizes", call. = FALSE)
    if (any(vapply(out[integer_keys], function(x) is.na(x) || x < 0L, logical(1L)))) {
        stop("Integer thresholds must be non-negative", call. = FALSE)
    }
    if (any(vapply(out[optional_integer_keys], function(x) !is.na(x) && x < 0L, logical(1L)))) {
        stop("Optional fixed-sum thresholds must be non-negative", call. = FALSE)
    }
    if (is.na(out$fdr) || out$fdr <= 0 || out$fdr >= 1)
        stop("--fdr must be between 0 and 1", call. = FALSE)
    if (is.na(out$min_stalling_log2fc) || out$min_stalling_log2fc < 0) {
        stop("--min-stalling-log2fc must be non-negative", call. = FALSE)
    }
    if (!is.na(out$min_transcript_median_rna) && out$min_transcript_median_rna < 0) {
        stop("--min-transcript-median-rna must be non-negative", call. = FALSE)
    }
    if (!is.na(out$min_transcript_median_rpf) && out$min_transcript_median_rpf < 0) {
        stop("--min-transcript-median-rpf must be non-negative", call. = FALSE)
    }
    if (is.na(out$sample_fraction) || out$sample_fraction <= 0 || out$sample_fraction > 1) {
        stop("--sample-fraction must be in (0, 1]", call. = FALSE)
    }
    if (!is.null(out$run_label) && !grepl("^[A-Za-z0-9_.-]+$", out$run_label)) {
        stop("--run-label may contain only letters, digits, dot, underscore, and hyphen",
            call. = FALSE)
    }
    if (is.null(out$rna_scale_manifest) || !nzchar(out$rna_scale_manifest)) {
        stop("--rna-scale-manifest is required", call. = FALSE)
    }
    if (is.null(out$focal_reference_manifest))
        stop("--focal-reference-manifest is required", call. = FALSE)
    if (!is.null(out$summary_file) && !nzchar(out$summary_file))
        out$summary_file <- NULL
    out
}

load_rna_scale_manifest <- function(path) {
    if (!file.exists(path))
        stop("Missing RNA scale manifest: ", path, call. = FALSE)
    scale_table <- fread(path, colClasses = c(sample_id = "character"), showProgress = FALSE)
    required <- c("sample_id", "rna_scale")
    missing <- setdiff(required, names(scale_table))
    if (length(missing)) {
        stop("RNA scale manifest is missing columns: ", paste(missing, collapse = ", "),
            call. = FALSE)
    }
    if (anyNA(scale_table$sample_id) || any(!nzchar(trimws(scale_table$sample_id)))) {
        stop("RNA scale manifest contains an empty sample_id", call. = FALSE)
    }
    scale_table[, sample_id := trimws(sample_id)]
    if (anyDuplicated(scale_table$sample_id)) {
        stop("RNA scale manifest contains duplicate sample_id values", call. = FALSE)
    }
    scale_table[, rna_scale := as.numeric(rna_scale)]
    if (any(!is.finite(scale_table$rna_scale)) || any(scale_table$rna_scale <= 0)) {
        stop("All rna_scale values must be finite and positive", call. = FALSE)
    }
    scale_table
}

replicate_labels <- function(column_names) {
    hits <- grep("^rep[1-9][0-9]*__rna_count$", column_names, value = TRUE)
    reps <- sub("__rna_count$", "", hits)
    reps[order(as.integer(sub("^rep", "", reps)))]
}

manifest_labels <- function(value) {
    if (!length(value) || is.na(value[[1L]]) || !nzchar(trimws(value[[1L]]))) {
        return(character())
    }
    trimws(strsplit(as.character(value[[1L]]), ",", fixed = TRUE)[[1L]])
}

sum_columns <- function(dt, columns) {
    answer <- numeric(nrow(dt))
    for (column in columns) answer <- answer + as.numeric(dt[[column]])
    answer
}

estimate_assay_normalization <- function(dt, reps, assay, library_values, method) {
    count_suffix <- if (assay == "RNA")
        "rna_count" else "rpf_count"
    library_suffix <- if (assay == "RNA")
        "library_RNA_sum" else "library_RPF_sum"
    count_columns <- paste0(reps, "__", count_suffix)
    raw_library_sizes <- vapply(reps, function(rep) library_values[[paste0(rep, "__", library_suffix)]],
        numeric(1L))
    normalization_factors <- rep(1, length(reps))
    if (method == "TMM" && length(reps) >= 2L) {
        assay_counts <- vapply(count_columns, function(column) as.numeric(dt[[column]]),
            numeric(nrow(dt)))
        if (!is.matrix(assay_counts))
            assay_counts <- matrix(assay_counts, nrow = nrow(dt))
        normalization_dge <- edgeR::DGEList(counts = assay_counts, lib.size = raw_library_sizes)
        normalization_dge <- edgeR::calcNormFactors(normalization_dge, method = "TMM")
        normalization_factors <- normalization_dge$samples$norm.factors
        rm(assay_counts, normalization_dge)
        invisible(gc())
    }
    data.table(replicate = reps, assay = assay, raw_library_size = raw_library_sizes, normalization_method = method,
        normalization_factor = normalization_factors, effective_library_size = raw_library_sizes *
            normalization_factors)
}

transcript_median_filter <- function(dt, reps, min_rna, min_rpf, rule) {
    # Counts are represented as zero for no coverage. The rule is either strict
    # across all replicates, or requires one replicate with both RNA and RPF
    # medians above their thresholds.
    median_columns <- c(paste0(reps, "__rna_count"), paste0(reps, "__rpf_count"))
    tx_medians <- dt[, lapply(.SD, function(x) median(as.numeric(x), na.rm = TRUE)), by = transcript,
        .SDcols = median_columns]
    rna_pass <- matrix(TRUE, nrow = nrow(tx_medians), ncol = length(reps))
    rpf_pass <- matrix(TRUE, nrow = nrow(tx_medians), ncol = length(reps))
    if (!is.na(min_rna)) {
        for (j in seq_along(reps)) {
            column <- paste0(reps[[j]], "__rna_count")
            rna_pass[, j] <- is.finite(tx_medians[[column]]) & tx_medians[[column]] > min_rna
        }
    }
    if (!is.na(min_rpf)) {
        for (j in seq_along(reps)) {
            column <- paste0(reps[[j]], "__rpf_count")
            rpf_pass[, j] <- is.finite(tx_medians[[column]]) & tx_medians[[column]] > min_rpf
        }
    }
    tx_keep <- if (rule == "all") {
        rowSums(rna_pass) == ncol(rna_pass) & rowSums(rpf_pass) == ncol(rpf_pass)
    } else {
        rowSums(rna_pass & rpf_pass) >= 1L
    }
    tx_index <- match(as.character(dt$transcript), as.character(tx_medians$transcript))
    site_keep <- tx_keep[tx_index]
    site_keep[is.na(site_keep)] <- FALSE
    list(site_keep = site_keep, n_transcripts = nrow(tx_medians), n_transcripts_pass = sum(tx_keep),
        n_sites_pass = sum(site_keep))
}

analyze_master <- function(info, options) {
    group_id <- info$group_id[[1L]]
    flank_size <- as.integer(info$flank_size[[1L]])
    manifest_reps <- manifest_labels(info$replicates)
    manifest_sample_ids <- manifest_labels(info$sample_ids)
    if (length(manifest_sample_ids) != length(manifest_reps)) {
        stop(group_id, ": manifest replicate/sample count mismatch", call. = FALSE)
    }
    sample_by_rep <- setNames(manifest_sample_ids, manifest_reps)
    head_reps <- if ("head_replicates" %in% names(info)) {
        manifest_labels(info$head_replicates)
    } else {
        character()
    }
    analysis_reps <- if (options$exclude_head_replicates) {
        setdiff(manifest_reps, head_reps)
    } else {
        manifest_reps
    }
    analysis_sample_ids <- unname(sample_by_rep[analysis_reps])
    rna_scale_by_rep <- options$rna_scale_lookup[analysis_sample_ids]
    if (anyNA(rna_scale_by_rep)) {
        stop(group_id, ": samples missing from RNA scale manifest: ", paste(analysis_sample_ids[is.na(rna_scale_by_rep)],
            collapse = ", "), call. = FALSE)
    }
    names(rna_scale_by_rep) <- analysis_reps
    n_replicates_expected <- length(analysis_reps)
    threshold_for_group <- function(fixed_sum, per_rep) {
        if (!is.na(fixed_sum))
            fixed_sum else per_rep * n_replicates_expected
    }
    min_flank_n_sum_required <- threshold_for_group(options$min_flank_n_sum, options$min_flank_n_per_rep)
    min_focal_rna_sum_required <- threshold_for_group(options$min_focal_rna_sum, options$min_focal_rna_per_rep)
    min_flank_rna_sum_required <- threshold_for_group(options$min_flank_rna_sum, options$min_flank_rna_per_rep)
    input_path <- info$output_path[[1L]]
    output_stem <- sprintf("%s.flank%d", group_id, flank_size)
    if (!is.null(options$run_label))
        output_stem <- paste(output_stem, options$run_label, sep = ".")
    output_prefix <- file.path(options$output_dir, output_stem)
    all_output <- paste0(output_prefix, ".edgeR.tsv.gz")
    significant_output <- paste0(output_prefix, ".significant.tsv.gz")
    parameter_output <- paste0(output_prefix, ".analysis_parameters.tsv")
    normalization_output <- paste0(output_prefix, ".normalization_factors.tsv")
    model_diagnostic_output <- paste0(output_prefix, ".model_diagnostics.tsv")
    design_output <- paste0(output_prefix, ".design.tsv")
    observation_output <- paste0(output_prefix, ".observation_info.tsv")

    if (!file.exists(input_path))
        stop("Missing master table: ", input_path, call. = FALSE)
    input_info <- file.info(input_path)
    manifest_value <- function(column) {
        if (column %in% names(info))
            info[[column]][[1L]] else NA_character_
    }
    reference_manifest <- fread(options$focal_reference_manifest)
    reference <- reference_manifest[species == info$species[[1L]]]
    if (nrow(reference) != 1L)
        stop("Expected one focal reference for species ", info$species[[1L]])
    current_parameters <- make_parameter_table(list(focal_reference_manifest = normalizePath(options$focal_reference_manifest),
        cds_lengths_file = normalizePath(reference$cds_lengths), cds_lengths_md5 = unname(tools::md5sum(reference$cds_lengths)[[1L]]),
        m6a_map_file = normalizePath(reference$m6a_map), m6a_map_md5 = unname(tools::md5sum(reference$m6a_map)[[1L]]),
        min_cds_codons_excluding_stop = 150L, focal_position_rule = "16 <= codon_index_1based <= codonlength - 5",
        focal_m6a_excluded = TRUE, focal_edited_retained = TRUE, flank_mask_rule = "first 15, last 5 non-stop codons, m6A and known editing codons excluded",
        transcript_RNA_fraction_filter = "none", median_and_TMM_universe = "complete unscaled master rows, RNA scaling applied before TMM and medians",
        analysis_schema_version = analysis_schema_version, script_path = edgeR_script_path,
        script_md5 = edgeR_script_md5, R_version = as.character(getRversion()), edgeR_version = as.character(packageVersion("edgeR")),
        data_table_version = as.character(packageVersion("data.table")), group_id = group_id,
        flank_size = flank_size, n_replicates = n_replicates_expected, manifest_replicate_labels = paste(manifest_reps,
            collapse = ","), head_replicate_labels = paste(head_reps, collapse = ","), exclude_head_replicates = options$exclude_head_replicates,
        analyzed_replicate_labels = paste(analysis_reps, collapse = ","), sample_ids = manifest_value("sample_ids"),
        analyzed_sample_ids = paste(analysis_sample_ids, collapse = ","), rna_scale_manifest = normalizePath(options$rna_scale_manifest,
            mustWork = TRUE), rna_scale_manifest_md5 = options$rna_scale_manifest_md5, rna_scale_by_replicate = paste(paste0(analysis_reps,
            "=", format(rna_scale_by_rep, trim = TRUE)), collapse = ","), rna_scaled_fields = "rna_count,flank_rna_count,library_RNA_sum",
        rpf_scaled = FALSE, input_path = normalizePath(input_path, mustWork = TRUE), input_size_bytes = input_info$size[[1L]],
        input_mtime_utc = format(input_info$mtime[[1L]], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        input_rows_manifest = manifest_value("rows"), min_replicates = options$min_replicates,
        min_focal_rna_per_rep = options$min_focal_rna_per_rep, min_focal_rna_sum_override = options$min_focal_rna_sum,
        min_focal_rna_sum_required = min_focal_rna_sum_required, min_flank_rna_per_rep = options$min_flank_rna_per_rep,
        min_flank_rna_sum_override = options$min_flank_rna_sum, min_flank_rna_sum_required = min_flank_rna_sum_required,
        min_flank_n_per_rep = options$min_flank_n_per_rep, min_flank_n_sum_override = options$min_flank_n_sum,
        min_flank_n_sum_required = min_flank_n_sum_required, min_focal_rpf_sum = options$min_focal_rpf_sum,
        min_transcript_median_rna = options$min_transcript_median_rna, min_transcript_median_rpf = options$min_transcript_median_rpf,
        transcript_median_rule = options$transcript_median_rule, normalization_method = options$normalization_method,
        fdr = options$fdr, p_value_alternative_stalling = "greater", p_adjust_method = "BH primary and BY sensitivity within group and flank size",
        min_stalling_log2fc = options$min_stalling_log2fc, sample_fraction = options$sample_fraction,
        sample_seed = options$sample_seed, run_label = options$run_label, model_formula = "~ library + region + ribo_focal",
        tested_coefficient = "ribo_focal", normalization = paste0(options$normalization_method,
            " within RNA/RPF separately; explicit effective-library-size offset; flank adds log(flank_n)"),
        abundance_covariate = "AveLogCPM computed from the complete explicit offset", ql_fit = "edgeR glmQLFit robust=TRUE"))

    if (file.exists(all_output) && !options$overwrite) {
        expected_outputs <- c(all_output, significant_output, parameter_output, normalization_output,
            model_diagnostic_output, design_output, observation_output)
        missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
        if (length(missing_outputs)) {
            stop("Existing edgeR result is incomplete; missing companion output(s): ", paste(missing_outputs,
                collapse = ", "), ". Use --overwrite or a new output directory.", call. = FALSE)
        }
        if (!file.exists(parameter_output)) {
            stop("Existing edgeR result has no parameter record: ", parameter_output, ". Refusing to reuse an unverifiable result. Use --overwrite or a new output directory.",
                call. = FALSE)
        }
        previous_parameters <- fread(parameter_output, colClasses = "character", na.strings = NULL,
            showProgress = FALSE)
        differing_parameters <- parameter_differences(previous_parameters, current_parameters)
        if (length(differing_parameters)) {
            old_values <- setNames(previous_parameters$value, previous_parameters$parameter)
            new_values <- setNames(current_parameters$value, current_parameters$parameter)
            shown <- head(differing_parameters, 12L)
            difference_text <- vapply(shown, function(key) {
                old <- old_values[[key]]
                new <- new_values[[key]]
                if (is.null(old))
                  old <- "<missing>"
                if (is.null(new))
                  new <- "<missing>"
                sprintf("%s: existing=%s; current=%s", key, old, new)
            }, character(1L))
            stop("Existing edgeR result was generated with different input/code/parameters. ",
                "Refusing to reuse it. Differences: ", paste(difference_text, collapse = " | "),
                if (length(differing_parameters) > length(shown))
                  " | ..." else "", ". Use --overwrite or a new output directory.", call. = FALSE)
        }
        message("[skip existing] ", all_output)
        existing_columns <- "stalling"
        existing <- fread(all_output, select = existing_columns, integer64 = "double", showProgress = FALSE)
        existing_diagnostics <- fread(model_diagnostic_output, showProgress = FALSE)
        existing_summary <- data.table(group_id = group_id, flank_size = flank_size, status = "existing",
            n_replicates = existing_diagnostics$n_replicates[[1L]], n_input = existing_diagnostics$n_input[[1L]],
            n_eligible = existing_diagnostics$n_eligible[[1L]], n_tested = nrow(existing),
            n_significant = sum(existing$stalling), min_focal_rna_per_rep = options$min_focal_rna_per_rep,
            min_focal_rna_sum_required = min_focal_rna_sum_required, min_flank_rna_per_rep = options$min_flank_rna_per_rep,
            min_flank_rna_sum_required = min_flank_rna_sum_required, min_flank_n_per_rep = options$min_flank_n_per_rep,
            min_flank_n_sum_required = min_flank_n_sum_required, transcript_median_rule = options$transcript_median_rule,
            output_path = normalizePath(all_output), parameter_path = normalizePath(parameter_output),
            normalization_path = normalizePath(normalization_output), model_diagnostic_path = normalizePath(model_diagnostic_output),
            design_path = normalizePath(design_output), observation_info_path = normalizePath(observation_output))
        return(existing_summary)
    }

    header <- names(fread(input_path, nrows = 0L, integer64 = "double", showProgress = FALSE))
    header_reps <- replicate_labels(header)
    if (!setequal(header_reps, manifest_reps)) {
        stop(group_id, ": manifest/header replicate labels mismatch", call. = FALSE)
    }
    reps <- header_reps[header_reps %in% analysis_reps]
    if (!identical(reps, analysis_reps)) {
        stop(group_id, ": analyzed replicate order/labels mismatch", call. = FALSE)
    }
    if (length(reps) < options$min_replicates) {
        message("[skip] ", group_id, " flank", flank_size, ": only ", length(reps), " replicate(s)")
        return(data.table(group_id = group_id, flank_size = flank_size, status = "too_few_replicates",
            n_input = as.integer(info$rows[[1L]]), n_eligible = 0L, n_tested = 0L, n_significant = 0L,
            output_path = NA_character_))
    }

    input_meta_columns <- c("chrom", "genomic_codon_first_pos_1based", "transcript", "cds_block_index",
        "strand", "cds_codon_first_pos_1based", "codon_index_1based")
    meta_columns <- c(input_meta_columns, "site_id")
    count_suffixes <- c("rna_count", "rpf_count", "flank_rna_count", "flank_rpf_count", "flank_n")
    # site_id is reconstructed only after filtering to avoid loading millions of
    # unique strings from a multi-GB input table.
    data_columns <- c(input_meta_columns, unlist(lapply(reps, function(rep) paste0(rep, "__",
        count_suffixes))))
    library_columns <- unlist(lapply(reps, function(rep) {
        paste0(rep, "__", c("library_RNA_sum", "library_RPF_sum"))
    }))
    required <- c(data_columns, library_columns)
    missing <- setdiff(required, header)
    if (length(missing))
        stop(group_id, ": missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

    # Library sums can exceed 2^31. Reading them as character first avoids an
    # implicit integer64 representation (and a dependency on bit64) before the
    # explicit, exact conversion to IEEE-754 numeric values.
    library_row <- fread(input_path, select = library_columns, nrows = 1L, colClasses = "character",
        showProgress = FALSE)
    library_values <- vapply(library_columns, function(column) as.numeric(library_row[[column]][[1L]]),
        numeric(1L))
    if (any(!is.finite(library_values)) || any(library_values <= 0)) {
        stop(group_id, ": library sizes must be positive", call. = FALSE)
    }
    original_library_values <- library_values
    for (rep in reps) {
        rna_library_column <- paste0(rep, "__library_RNA_sum")
        library_values[[rna_library_column]] <- library_values[[rna_library_column]] * rna_scale_by_rep[[rep]]
    }

    message("[read] ", group_id, " flank", flank_size, " (", length(reps), " reps)")
    dt <- fread(input_path, select = data_columns, integer64 = "double", showProgress = TRUE)
    n_input <- nrow(dt)
    lengths <- fread(reference$cds_lengths)
    if (anyDuplicated(lengths$transcript))
        stop("Duplicate transcript length reference")
    dt[, codonlength := lengths$codonlength[match(transcript, lengths$transcript)]]
    if (anyNA(dt$codonlength))
        stop("Master transcript missing from CDS length reference")
    tx_check <- dt[, .(n = .N, max_pos = max(codon_index_1based), L = codonlength[1L]), by = transcript]
    if (any(tx_check$n != tx_check$L | tx_check$max_pos != tx_check$L))
        stop("Master does not contain the complete reference CDS codon universe")
    if (any(dt$cds_codon_first_pos_1based != 3 * dt$codon_index_1based - 2))
        stop("Inconsistent codon coordinates")
    m6a <- fread(reference$m6a_map)
    m6a <- unique(m6a[mapping_status == "mapped", .(transcript = selected_transcript, cds_codon_first_pos_1based = as.integer(codon_first_cds_pos_1based))])
    m6a_rows <- dt[m6a, on = .(transcript, cds_codon_first_pos_1based), which = TRUE, nomatch = NULL]
    focal_eligible <- dt$codonlength >= 150L & dt$codon_index_1based > 15L & dt$codon_index_1based <=
        dt$codonlength - 5L
    focal_eligible[m6a_rows] <- FALSE
    n_focal_structural <- sum(focal_eligible)
    message("[focal CDS/position/m6A] ", group_id, ": ", n_focal_structural, " / ", n_input)
    meta_columns <- c(meta_columns, "codonlength")
    rm(lengths, tx_check, m6a, m6a_rows)

    # The master tables contain focal single-codon RNA counts and flank sums of
    # the same per-codon RNA counts. Scaling both fields before filtering and
    # fitting is algebraically identical to scaling every input codon first and
    # then rebuilding the flank sums. RPF counts and flank exposures are intact.
    for (rep in reps) {
        scale <- rna_scale_by_rep[[rep]]
        for (suffix in c("rna_count", "flank_rna_count")) {
            column <- paste0(rep, "__", suffix)
            set(dt, j = column, value = as.numeric(dt[[column]]) * scale)
        }
    }
    message("[RNA scaling] ", group_id, ": ", paste(paste0(reps, "=", format(rna_scale_by_rep[reps],
        trim = TRUE)), collapse = ", "))

    message("[normalization] ", group_id, " flank", flank_size, ": ", options$normalization_method)
    normalization_table <- rbindlist(list(estimate_assay_normalization(dt, reps, "RNA", library_values,
        options$normalization_method), estimate_assay_normalization(dt, reps, "Ribo", library_values,
        options$normalization_method)))
    normalization_table[, sample_id := unname(sample_by_rep[replicate])]
    normalization_table[, rna_scale := fifelse(assay == "RNA", as.numeric(rna_scale_by_rep[as.character(replicate)]),
        1)]
    normalization_table[, original_library_size := vapply(seq_len(.N), function(i) {
        suffix <- if (assay[[i]] == "RNA")
            "library_RNA_sum" else "library_RPF_sum"
        original_library_values[[paste0(replicate[[i]], "__", suffix)]]
    }, numeric(1L))]
    setcolorder(normalization_table, c("replicate", "sample_id", "assay", "rna_scale", "original_library_size",
        "raw_library_size", "normalization_method", "normalization_factor", "effective_library_size"))
    normalization_table[, library_key := paste(replicate, assay, sep = "__")]
    normalization_factor_lookup <- setNames(normalization_table$normalization_factor, normalization_table$library_key)

    focal_rpf_columns <- paste0(reps, "__rpf_count")
    flank_rpf_columns <- paste0(reps, "__flank_rpf_count")
    focal_rna_columns <- paste0(reps, "__rna_count")
    flank_rna_columns <- paste0(reps, "__flank_rna_count")
    flank_n_columns <- paste0(reps, "__flank_n")

    focal_rpf_sum <- sum_columns(dt, focal_rpf_columns)
    flank_rpf_sum <- sum_columns(dt, flank_rpf_columns)
    focal_rna_sum <- sum_columns(dt, focal_rna_columns)
    flank_rna_sum <- sum_columns(dt, flank_rna_columns)
    min_observed_flank_n <- rep.int(Inf, n_input)
    flank_n_sum <- numeric(n_input)
    for (column in flank_n_columns) {
        min_observed_flank_n <- pmin(min_observed_flank_n, as.numeric(dt[[column]]))
        flank_n_sum <- flank_n_sum + as.numeric(dt[[column]])
    }
    message("[coverage thresholds] ", group_id, " flank", flank_size, " (", n_replicates_expected,
        " reps): focal_RNA_sum >= ", min_focal_rna_sum_required, "; flank_RNA_sum >= ", min_flank_rna_sum_required,
        "; flank_n_sum >= ", min_flank_n_sum_required)
    # RNA coverage is an outcome-independent eligibility filter. By default,
    # each total threshold is its per-replicate standard multiplied by the
    # number of biological replicates. A --min-*-sum option explicitly
    # overrides the corresponding scaled threshold. Every individual flank_n
    # must also be positive because it enters a logarithmic offset.
    keep <- (focal_eligible & min_observed_flank_n > 0 & flank_n_sum >= min_flank_n_sum_required &
        focal_rna_sum >= min_focal_rna_sum_required & flank_rna_sum >= min_flank_rna_sum_required &
        focal_rpf_sum >= options$min_focal_rpf_sum)
    if (!is.na(options$min_transcript_median_rna) || !is.na(options$min_transcript_median_rpf)) {
        tx_filter <- transcript_median_filter(dt, reps, options$min_transcript_median_rna,
            options$min_transcript_median_rpf, options$transcript_median_rule)
        keep <- keep & tx_filter$site_keep
        message("[transcript medians] ", group_id, " flank", flank_size, ": ", format(tx_filter$n_transcripts_pass,
            big.mark = ","), " / ", format(tx_filter$n_transcripts, big.mark = ","), " transcripts; ",
            format(tx_filter$n_sites_pass, big.mark = ","), " sites pass")
        rm(tx_filter)
    }
    if (!any(keep))
        stop(group_id, ": no sites passed the count filters", call. = FALSE)

    dt <- dt[keep]
    focal_rpf_sum <- focal_rpf_sum[keep]
    flank_rpf_sum <- flank_rpf_sum[keep]
    focal_rna_sum <- focal_rna_sum[keep]
    flank_rna_sum <- flank_rna_sum[keep]
    min_observed_flank_n <- min_observed_flank_n[keep]
    flank_n_sum <- flank_n_sum[keep]
    rm(keep, focal_eligible)
    n_eligible <- nrow(dt)

    if (options$sample_fraction < 1) {
        sample_n <- max(1L, floor(n_eligible * options$sample_fraction))
        set.seed(options$sample_seed)
        selected <- sort(sample.int(n_eligible, size = sample_n, replace = FALSE))
        dt <- dt[selected]
        focal_rpf_sum <- focal_rpf_sum[selected]
        flank_rpf_sum <- flank_rpf_sum[selected]
        focal_rna_sum <- focal_rna_sum[selected]
        flank_rna_sum <- flank_rna_sum[selected]
        min_observed_flank_n <- min_observed_flank_n[selected]
        flank_n_sum <- flank_n_sum[selected]
        rm(selected)
        message("[sample] ", group_id, " flank", flank_size, ": ", format(nrow(dt), big.mark = ","),
            " / ", format(n_eligible, big.mark = ","), " eligible sites (seed ", options$sample_seed,
            ")")
    }
    dt[, site_id := paste0(transcript, "_", codon_index_1based)]
    invisible(gc())

    observation_info <- rbindlist(lapply(reps, function(rep) {
        data.table(replicate = rep, sample_id = sample_by_rep[[rep]], assay = c("RNA", "RNA",
            "Ribo", "Ribo"), region = c("Flank", "Focal", "Flank", "Focal"), count_column = paste0(rep,
            "__", c("flank_rna_count", "rna_count", "flank_rpf_count", "rpf_count")), flank_n_column = c(paste0(rep,
            "__flank_n"), NA_character_, paste0(rep, "__flank_n"), NA_character_), library_size = c(library_values[[paste0(rep,
            "__library_RNA_sum")]], library_values[[paste0(rep, "__library_RNA_sum")]], library_values[[paste0(rep,
            "__library_RPF_sum")]], library_values[[paste0(rep, "__library_RPF_sum")]]),
            rna_scale = c(rna_scale_by_rep[[rep]], rna_scale_by_rep[[rep]], 1, 1))
    }))
    observation_info[, replicate := factor(replicate, levels = reps)]
    observation_info[, assay := factor(assay, levels = c("RNA", "Ribo"))]
    observation_info[, region := factor(region, levels = c("Flank", "Focal"))]
    observation_info[, library_key := paste(replicate, assay, sep = "__")]
    observation_info[, normalization_factor := normalization_factor_lookup[library_key]]
    if (any(!is.finite(observation_info$normalization_factor))) {
        stop("Missing or invalid normalization factor", call. = FALSE)
    }
    observation_info[, effective_library_size := library_size * normalization_factor]
    observation_info[, library := factor(library_key, levels = unique(library_key))]
    observation_info[, ribo_focal := as.integer(assay == "Ribo" & region == "Focal")]
    design <- model.matrix(~library + region + ribo_focal, data = observation_info)
    coefficient <- "ribo_focal"
    if (!coefficient %in% colnames(design))
        stop("Expected coefficient absent from design", call. = FALSE)
    if (qr(design)$rank != ncol(design))
        stop("Design matrix is not full rank", call. = FALSE)

    counts <- vapply(observation_info$count_column, function(column) as.numeric(dt[[column]]),
        numeric(nrow(dt)))
    if (!is.matrix(counts))
        counts <- matrix(counts, nrow = nrow(dt))
    colnames(counts) <- paste(observation_info$replicate, observation_info$assay, observation_info$region,
        sep = "_")

    offsets <- matrix(rep(log(observation_info$effective_library_size), each = nrow(dt)),
        nrow = nrow(dt), ncol = nrow(observation_info))
    flank_observations <- which(observation_info$region == "Flank")
    for (column_index in flank_observations) {
        offsets[, column_index] <- offsets[, column_index] + log(as.numeric(dt[[observation_info$flank_n_column[[column_index]]]]))
    }
    metadata <- dt[, ..meta_columns]
    rm(dt)
    invisible(gc())

    message("[fit] ", group_id, " flank", flank_size, ": ", format(nrow(counts), big.mark = ","),
        " sites")
    dge <- edgeR::DGEList(counts = counts, lib.size = observation_info$effective_library_size)
    dge$samples$norm.factors <- 1
    dge$offset <- offsets
    # glmQLFit.DGEList otherwise computes AveLogCPM only from DGEList library
    # sizes. Compute it explicitly from the same full offset used by the model,
    # including both TMM factors and each focal site's actual flank_n exposure.
    dge$AveLogCPM <- edgeR::aveLogCPM(dge$counts, offset = offsets)
    # Current edgeR QL pipeline estimates the overall NB dispersion directly in
    # glmQLFit when dispersion is NULL. Avoiding a preceding estimateDisp pass is
    # important for multi-million-row codon tables and does not change the model
    # or coefficient being tested.
    fit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, coef = coefficient)
    stats <- as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table)
    setnames(stats, c("PValue", "FDR"), c("PValue_two_sided", "FDR_two_sided"))
    stats[, PValue_greater := fifelse(logFC > 0, PValue_two_sided/2, fifelse(logFC < 0,
        1 - PValue_two_sided/2, 0.5))]
    stats[, `:=`(FDR_greater = p.adjust(PValue_greater, method = "BH"), FDR_greater_BY = p.adjust(PValue_greater,
        method = "BY"))]

    result <- cbind(metadata, data.table(focal_RNA_sum = focal_rna_sum, flank_RNA_sum = flank_rna_sum,
        focal_RPF_sum = focal_rpf_sum, flank_RPF_sum = flank_rpf_sum, flank_n_sum = flank_n_sum,
        min_flank_n = min_observed_flank_n), stats)
    # The scientific hypotheses are directional and pre-specified. For a
    # one-coefficient QL F test, combine the two-sided P value with the fitted
    # coefficient sign to obtain one-sided greater P values, then perform
    # BH correction across the complete tested codon universe for this group.
    result[, stalling := logFC > options$min_stalling_log2fc & FDR_greater <= options$fdr]
    result[, stalling_BY := logFC > options$min_stalling_log2fc & FDR_greater_BY <= options$fdr]
    dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(result, all_output, sep = "\t", quote = FALSE, compress = "gzip")
    fwrite(result[stalling == TRUE], significant_output, sep = "\t", quote = FALSE, compress = "gzip")
    normalization_table[, library_key := NULL]
    fwrite(normalization_table, normalization_output, sep = "\t", quote = FALSE)
    observation_export <- copy(observation_info)
    for (column in names(observation_export)) {
        if (is.factor(observation_export[[column]])) {
            set(observation_export, j = column, value = as.character(observation_export[[column]]))
        }
    }
    fwrite(observation_export, observation_output, sep = "\t", quote = FALSE)
    design_export <- as.data.table(design)
    design_export[, observation := colnames(counts)]
    setcolorder(design_export, c("observation", setdiff(names(design_export), "observation")))
    fwrite(design_export, design_output, sep = "\t", quote = FALSE)
    model_diagnostics <- data.table(group_id = group_id, flank_size = flank_size, n_replicates = n_replicates_expected,
        design_observations = nrow(design), design_columns = ncol(design), design_rank = qr(design)$rank,
        residual_df_design = nrow(design) - qr(design)$rank, n_input = n_input, n_eligible = n_eligible,
        n_after_focal_structure = n_focal_structural, n_tested = nrow(result), nb_dispersion_median = median(as.numeric(fit$dispersion),
            na.rm = TRUE), ql_prior_df_median = median(as.numeric(fit$df.prior), na.rm = TRUE),
        ql_posterior_variance_median = median(as.numeric(fit$s2.post), na.rm = TRUE), p_greater_le_0_01_fraction = mean(result$PValue_greater <=
            0.01), p_greater_le_0_05_fraction = mean(result$PValue_greater <= 0.05), bh_stalling_fraction = mean(result$stalling),
        by_stalling_fraction = mean(result$stalling_BY), median_logFC = median(result$logFC),
        focal_RPF_zero_fraction = mean(result$focal_RPF_sum == 0), min_logCPM = min(result$logCPM),
        median_logCPM = median(result$logCPM), max_logCPM = max(result$logCPM))
    fwrite(model_diagnostics, model_diagnostic_output, sep = "\t", quote = FALSE)
    fwrite(current_parameters, parameter_output, sep = "\t", quote = FALSE)
    n_significant <- sum(result$stalling)
    message("[done] ", group_id, " flank", flank_size, ": ", n_significant, " peak sites")
    summary_row <- data.table(group_id = group_id, flank_size = flank_size, status = "completed",
        n_replicates = n_replicates_expected, n_input = n_input, n_eligible = n_eligible,
        n_tested = nrow(result), n_significant = n_significant, min_focal_rna_per_rep = options$min_focal_rna_per_rep,
        min_focal_rna_sum_required = min_focal_rna_sum_required, min_flank_rna_per_rep = options$min_flank_rna_per_rep,
        min_flank_rna_sum_required = min_flank_rna_sum_required, min_flank_n_per_rep = options$min_flank_n_per_rep,
        min_flank_n_sum_required = min_flank_n_sum_required, transcript_median_rule = options$transcript_median_rule,
        output_path = normalizePath(all_output), parameter_path = normalizePath(parameter_output),
        normalization_path = normalizePath(normalization_output), model_diagnostic_path = normalizePath(model_diagnostic_output),
        design_path = normalizePath(design_output), observation_info_path = normalizePath(observation_output))
    summary_row
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
rna_scale_table <- load_rna_scale_manifest(options$rna_scale_manifest)
options$rna_scale_lookup <- setNames(rna_scale_table$rna_scale, rna_scale_table$sample_id)
options$rna_scale_manifest_md5 <- unname(tools::md5sum(options$rna_scale_manifest)[[1L]])
setDTthreads(options$threads)
dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
if (is.null(options$summary_file)) options$summary_file <- file.path(options$output_dir,
    "edgeR_run_summary.tsv")
dir.create(dirname(options$summary_file), recursive = TRUE, showWarnings = FALSE)
manifest <- fread(options$master_manifest, integer64 = "double")
required_manifest <- c("group_id", "flank_size", "n_replicates", "rows", "output_path")
if (length(setdiff(required_manifest, names(manifest)))) {
    stop("Invalid master manifest; missing: ", paste(setdiff(required_manifest, names(manifest)),
        collapse = ", "), call. = FALSE)
}
manifest <- manifest[flank_size %in% options$flank_sizes]
if (!is.null(options$groups)) manifest <- manifest[group_id %in% options$groups]
if (!nrow(manifest)) stop("No master tables selected", call. = FALSE)
setorder(manifest, group_id, flank_size)

summaries <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
    summaries[[i]] <- tryCatch(analyze_master(manifest[i], options), error = function(error) {
        message("[ERROR] ", manifest$group_id[[i]], " flank", manifest$flank_size[[i]], ": ",
            conditionMessage(error))
        data.table(group_id = manifest$group_id[[i]], flank_size = manifest$flank_size[[i]],
            status = paste0("error: ", conditionMessage(error)), n_input = NA_integer_, n_eligible = NA_integer_,
            n_tested = NA_integer_, n_significant = NA_integer_, output_path = NA_character_)
    })
    fwrite(rbindlist(summaries[seq_len(i)], fill = TRUE), options$summary_file, sep = "\t")
}
if (any(startsWith(vapply(summaries, function(x) x$status[[1L]], character(1L)), "error:"))) {
    quit(save = "no", status = 1L)
}
