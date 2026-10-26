#!/usr/bin/env Rscript

# Summarize edited/unedited edgeR v2 peak calls separately
# inside and outside CDS hairpin intervals. A codon is a hairpin site when its
# CDS interval
# [cds_codon_first_pos_1based, cds_codon_first_pos_1based + 2] overlaps a
# same-transcript hairpin interval. Overlapping hairpins are merged first.
#
# Within each hairpin stratum, rates use edgeR-tested sites in that stratum as
# denominator. The unedited background is restricted to codons containing at
# least one A, while mapped edited codons are retained unchanged. Sites filtered
# before edgeR are not included because they cannot have a peak call.

suppressPackageStartupMessages(library(data.table))
fread <- function(file, ...) {
    if (is.character(file) && length(file) == 1L && grepl("[.]gz$", file)) {
        data.table::fread(cmd = paste("gzip -dc --", shQuote(file)), ...)
    } else data.table::fread(file, ...)
}

args <- commandArgs(trailingOnly = TRUE)
cli_value <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx)) return(default)
    if (idx >= length(args)) stop("Missing value after ", flag, call. = FALSE)
    args[[idx + 1L]]
}
positive_int <- function(value, label) {
    answer <- suppressWarnings(as.integer(value))
    if (is.na(answer) || answer < 1L) {
        stop(label, " must be a positive integer", call. = FALSE)
    }
    answer
}

result_dir <- cli_value(
    "--result-dir", "edgeR_stalling_rigorous_v2/results_txmedian_any_both"
)
manifest_path <- cli_value(
    "--manifest", "stalling_master_tables/master_table_manifest.tsv"
)
annotation_dir <- cli_value("--annotation-dir", ".")
editing_dir <- cli_value(
    "--editing-dir",
    "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed"
)
hairpin_dir <- cli_value(
    "--hairpin-dir", "/data/dongxn/RNAediting/codex/02_stalling/3flank/hairpin"
)

codon_dir <- cli_value(
    "--codon-dir", "/path/to/RNA_editing_DATAprocessing/ref_genome"
)
output_dir <- cli_value(
    "--output-dir", file.path(result_dir, "edited_unedited_hairpin_summary")
)
run_label <- cli_value("--run-label", "rigorous_v2_txmedian_any_both")
jobs <- positive_int(cli_value("--jobs", Sys.getenv("SUMMARY_JOBS", "4")), "--jobs")
threads <- positive_int(
    cli_value("--threads", Sys.getenv("SUMMARY_THREADS", "2")), "--threads"
)
groups_arg <- cli_value("--groups", NULL)

for (path in c(result_dir, annotation_dir, editing_dir, hairpin_dir)) {
    if (!dir.exists(path)) stop("Missing directory: ", path)
}
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setDTthreads(threads)

species_specs <- data.table(
    species = c("dmel", "dere", "dsec", "dsim"),
    annotation_file = c(
        "dmelr6.12.gtf.cds.bed.ordered.1codon.longest.tsv",
        "derer1.05.gtf.cds.bed.ordered.1codon.longest.tsv",
        "dsecr1.3.gtf.cds.bed.ordered.1codon.longest.tsv",
        "dsimr2.02.gtf.cds.bed.ordered.1codon.longest.tsv"
    ),
    editing_file = c(
        "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
        "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
        "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
        "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
    ),
    hairpin_file = c(
        "Dmel.fa.T25Z3.2.merged.tsv",
        "Dere.fa.T25Z3.2.merged.tsv",
        "Dsec.fa.T25Z3.2.merged.tsv",
        "Dsim.fa.T25Z3.2.merged.tsv"
    ),
    codon_file = c(
        "dmel-all-CDS-r6.12.fasta.edited.tab.split",
        "dere-all-CDS-r1.05.fasta.edited.tab.split",
        "dsec-all-CDS-r1.3.fasta.edited.tab.split",
        "dsim-all-CDS-r2.02.fasta.edited.tab.split"
    )
)
species_specs[, annotation_path := file.path(annotation_dir, annotation_file)]
species_specs[, editing_path := file.path(editing_dir, editing_file)]
species_specs[, hairpin_path := file.path(hairpin_dir, hairpin_file)]
species_specs[, codon_path := file.path(codon_dir, codon_file)]

annotation_cols <- c(
    "chrom", "genomic_codon_first_pos_1based", "transcript",
    "cds_block_index", "strand", "cds_codon_first_pos_1based",
    "codon_index_1based"
)

load_codon_map <- function(spec) {
    x <- fread(
        spec$codon_path,
        header = FALSE,
        col.names = c("site", "codon"),
        showProgress = FALSE
    )
    x[, transcript := sub("_[0-9]+$", "", site)]
    x[, codon_index_1based := as.integer(sub(".*_", "", site))]
    unique(x[, .(
        transcript,
        codon_index_1based,
        codon = toupper(trimws(codon))
    )])
}

load_editing_map <- function(spec) {
    if (!file.exists(spec$annotation_path)) {
        stop("Missing annotation: ", spec$annotation_path)
    }
    if (!file.exists(spec$editing_path)) {
        stop("Missing editing file: ", spec$editing_path)
    }
    annotation <- fread(
        spec$annotation_path,
        header = FALSE, integer64 = "double",
        showProgress = FALSE
    )
    if (ncol(annotation) < length(annotation_cols)) {
        stop("Unexpected annotation columns in ", spec$annotation_path)
    }
    setnames(annotation, names(annotation)[seq_along(annotation_cols)], annotation_cols)
    annotation <- unique(annotation[, .(
        transcript = as.character(transcript),
        cds_codon_first_pos_1based = as.numeric(cds_codon_first_pos_1based),
        codon_index_1based = as.integer(codon_index_1based)
    )])
    annotation <- annotation[
        !is.na(transcript) & transcript != "" &
            is.finite(cds_codon_first_pos_1based) & !is.na(codon_index_1based)
    ]
    editing <- fread(
        spec$editing_path,
        select = c("new_transcript", "new_cds_pos_1based", "target_region"),
        integer64 = "double", showProgress = FALSE
    )
    input_rows <- nrow(editing)
    editing <- editing[
        target_region == "CDS" &
            !is.na(new_transcript) & new_transcript != "" &
            is.finite(as.numeric(new_cds_pos_1based))
    ]
    editing[, new_transcript := as.character(new_transcript)]
    editing[, new_cds_pos_1based := as.numeric(new_cds_pos_1based)]
    editing <- editing[new_cds_pos_1based >= 1]
    editing[, cds_codon_first_pos_1based :=
        1 + floor((new_cds_pos_1based - 1) / 3) * 3]
    editing <- unique(editing[, .(
        transcript = new_transcript, cds_codon_first_pos_1based
    )])
    mapped <- annotation[
        editing,
        on = .(transcript, cds_codon_first_pos_1based), nomatch = 0L
    ]
    mapped[, site_id := paste0(transcript, "_", codon_index_1based)]
    list(
        site_ids = unique(mapped$site_id),
        diagnostic = data.table(
            species = spec$species,
            editing_input_rows = input_rows,
            editing_cds_codon_rows = nrow(editing),
            editing_codon_rows_mapped = nrow(mapped),
            editing_site_ids_mapped = uniqueN(mapped$site_id),
            annotation_rows = nrow(annotation)
        )
    )
}

merge_hairpin_intervals <- function(spec) {
    if (!file.exists(spec$hairpin_path)) {
        stop("Missing hairpin file: ", spec$hairpin_path)
    }
    h <- fread(spec$hairpin_path, integer64 = "double", showProgress = FALSE)
    required <- c("transcript", "start", "end")
    missing <- setdiff(required, names(h))
    if (length(missing)) {
        stop(
            "Missing hairpin columns in ", spec$hairpin_path, ": ",
            paste(missing, collapse = ", ")
        )
    }
    h <- h[, .(
        transcript = as.character(transcript),
        start = as.numeric(start), end = as.numeric(end)
    )]
    h <- h[
        !is.na(transcript) & transcript != "" &
            is.finite(start) & is.finite(end) & start >= 1 & end >= start
    ]
    setorder(h, transcript, start, end)
    n_raw <- nrow(h)
    n_tx <- uniqueN(h$transcript)
    h[, running_end := cummax(end), by = transcript]
    h[, new_interval := is.na(shift(running_end)) |
        start > shift(running_end) + 1, by = transcript]
    h[, interval_group := cumsum(new_interval), by = transcript]
    h <- h[, .(
        start = min(start), end = max(end)
    ), by = .(transcript, interval_group)]
    h[, hairpin_interval_id := .I]
    list(
        intervals = h[, .(transcript, start, end, hairpin_interval_id)],
        diagnostic = data.table(
            species = spec$species, hairpin_file = spec$hairpin_path,
            hairpin_intervals_raw = n_raw,
            hairpin_intervals_merged = nrow(h),
            hairpin_transcripts = n_tx
        )
    )
}

message("Loading codon, editing and hairpin maps")

codon_maps <- setNames(
    lapply(seq_len(nrow(species_specs)), function(i) {
        load_codon_map(species_specs[i])
    }),
    species_specs$species
)
editing_maps <- setNames(
    lapply(seq_len(nrow(species_specs)), function(i) {
        message("  editing: ", species_specs$species[[i]])
        load_editing_map(species_specs[i])
    }),
    species_specs$species
)
hairpin_maps <- setNames(
    lapply(seq_len(nrow(species_specs)), function(i) {
        message("  hairpin: ", species_specs$species[[i]])
        merge_hairpin_intervals(species_specs[i])
    }),
    species_specs$species
)

manifest <- fread(manifest_path, integer64 = "double", showProgress = FALSE)
required_manifest <- c("group_id", "species", "genotype", "sex", "flank_size")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest)) {
    stop("Manifest missing columns: ", paste(missing_manifest, collapse = ", "))
}
manifest30 <- unique(manifest[
    flank_size == 30L,
    .(
        group_id = as.character(group_id), species = as.character(species),
        genotype = as.character(genotype), sex = as.character(sex),
        flank_size = as.integer(flank_size)
    )
])
parameter_pattern <- paste0(
    "\\.flank30\\.", run_label, "\\.analysis_parameters\\.tsv$"
)
parameter_files <- list.files(
    result_dir,
    pattern = parameter_pattern, full.names = TRUE
)
if (!length(parameter_files)) {
    stop("No parameter files found in ", result_dir)
}
group_from_parameter <- function(path) {
    stem <- sub("\\.analysis_parameters\\.tsv$", "", basename(path))
    sub(paste0("\\.flank30\\.", run_label, "$"), "", stem)
}
tasks <- data.table(
    parameter_path = parameter_files,
    group_id = vapply(parameter_files, group_from_parameter, character(1L))
)
tasks <- merge(tasks, manifest30, by = "group_id", all.x = TRUE, sort = FALSE)
if (anyNA(tasks$species)) {
    stop(
        "Could not map groups: ",
        paste(tasks[is.na(species), group_id], collapse = ", ")
    )
}
if (!is.null(groups_arg)) {
    requested <- trimws(strsplit(groups_arg, ",", fixed = TRUE)[[1L]])
    tasks <- tasks[group_id %chin% requested]
}
tasks[, edgeR_path := file.path(
    result_dir, paste0(group_id, ".flank30.", run_label, ".edgeR.tsv.gz")
)]
tasks[, complete_bundle := file.exists(edgeR_path)]
if (any(!tasks$complete_bundle)) {
    message(
        "Skipping incomplete bundles: ",
        paste(tasks[complete_bundle == FALSE, group_id], collapse = ", ")
    )
}
tasks <- tasks[complete_bundle == TRUE]
if (!nrow(tasks)) stop("No complete edgeR bundles remain")
setorder(tasks, species, genotype, sex, group_id)

parse_flag <- function(x) {
    if (is.logical(x)) return(x)
    toupper(as.character(x)) %chin% c("TRUE", "T", "1")
}

resolve_codon_column <- function(path) {
    header_names <- names(fread(
        path,
        nrows = 0L, integer64 = "double", showProgress = FALSE
    ))
    preferred <- c(
        "codon", "reference_codon", "ref_codon", "focal_codon",
        "codon_before", "codon_seq", "focal_codon_seq"
    )
    exact_hits <- preferred[preferred %chin% header_names]
    if (length(exact_hits) > 0L) return(exact_hits[[1L]])

    normalize_name <- function(x) {
        tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
    }
    header_keys <- normalize_name(header_names)
    preferred_keys <- normalize_name(preferred)
    normalized_hits <- match(preferred_keys, header_keys, nomatch = 0L)
    normalized_hits <- normalized_hits[normalized_hits > 0L]
    if (length(normalized_hits) > 0L) {
        return(header_names[[normalized_hits[[1L]]]])
    }

    stop(
        "Could not identify a codon-sequence column in ", path,
        ". Expected one of: ", paste(preferred, collapse = ", "),
        ". Available columns: ", paste(header_names, collapse = ", "),
        call. = FALSE
    )
}

safe_fraction <- function(numerator, denominator) {
    if (denominator > 0) numerator / denominator else NA_real_
}
fisher_greater <- function(
  edited_inside, edited_outside, unedited_inside, unedited_outside
) {
    if (edited_inside + edited_outside == 0L ||
        unedited_inside + unedited_outside == 0L) {
        return(c(odds_ratio = NA_real_, p_value = NA_real_))
    }
    tab <- matrix(
        c(edited_inside, edited_outside, unedited_inside, unedited_outside),
        nrow = 2L, byrow = TRUE
    )
    test <- stats::fisher.test(tab, alternative = "greater")
    c(
        odds_ratio = unname(as.numeric(test$estimate)),
        p_value = unname(as.numeric(test$p.value))
    )
}

find_hairpin_sites <- function(sites, hairpins) {
    if (!nrow(sites) || !nrow(hairpins)) return(character())
    x <- sites[, .(
        site_id = as.character(site_id),
        transcript = as.character(transcript),
        start = as.numeric(cds_codon_first_pos_1based),
        end = as.numeric(cds_codon_first_pos_1based) + 2
    )]
    x <- x[
        !is.na(site_id) & !is.na(transcript) &
            is.finite(start) & is.finite(end)
    ]
    x <- x[transcript %chin% hairpins$transcript]
    if (!nrow(x)) return(character())
    y <- copy(hairpins)
    setkey(y, transcript, start, end)
    overlaps <- foverlaps(
        x, y,
        by.x = c("transcript", "start", "end"),
        by.y = c("transcript", "start", "end"),
        type = "any", nomatch = 0L
    )
    site_column <- if ("i.site_id" %in% names(overlaps)) {
        "i.site_id"
    } else if ("site_id" %in% names(overlaps)) {
        "site_id"
    } else {
        stop("Could not identify site_id after hairpin overlap")
    }
    unique(as.character(overlaps[[site_column]]))
}

process_one <- function(i) {
    row <- tasks[i]
    species <- as.character(row$species)
    message("[start] ", row$group_id, " pid=", Sys.getpid())

    e <- fread(
        row$edgeR_path,
        select = c(
            "site_id",
            "transcript",
            "codon_index_1based",
            "cds_codon_first_pos_1based",
            "stalling"
        ),
        integer64 = "double",
        showProgress = FALSE
    )

    e <- merge(
        e,
        codon_maps[[species]],
        by = c("transcript", "codon_index_1based"),
        all.x = TRUE
    )

    if (anyNA(e$codon)) {
        stop(
            "Missing codon annotation in ",
            row$edgeR_path,
            ": ",
            sum(is.na(e$codon))
        )
    }

    e[, `:=`(
        site_id = as.character(site_id),
        transcript = as.character(transcript),
        cds_codon_first_pos_1based = as.numeric(cds_codon_first_pos_1based),
        codon = toupper(trimws(as.character(codon))),
        stalling = parse_flag(stalling)
    )]
    if (anyNA(e$site_id) || anyNA(e$stalling)) {
        stop("NA site_id or call flag in ", row$edgeR_path)
    }

    n_tested_before_A_filter <- nrow(e)
    hairpin_site_ids <- find_hairpin_sites(
        e, hairpin_maps[[species]]$intervals
    )
    e[, `:=`(
        edited = site_id %chin% editing_maps[[species]]$site_ids,
        hairpin = site_id %chin% hairpin_site_ids,
        codon_has_A = !is.na(codon) & grepl("A", codon, fixed = TRUE)
    )]

    n_unedited_nonA_removed <- sum(!e$edited & !e$codon_has_A)
    n_edited_without_A_retained <- sum(e$edited & !e$codon_has_A)

    # Only the unedited background is narrowed. All mapped edited codons are
    # retained; unedited codons must contain at least one A.
    e <- e[edited | codon_has_A]

    hairpin_in <- e[hairpin == TRUE]
    hairpin_out <- e[hairpin == FALSE]
    interval_qc <- hairpin_maps[[species]]$diagnostic

    summarize_region <- function(region_data, hairpin_status) {
        event_rows <- lapply("stalling", function(event) {
            inside <- region_data$stalling
            edited_total <- sum(region_data$edited)
            unedited_total <- sum(!region_data$edited)
            edited_inside <- sum(region_data$edited & inside)
            unedited_inside <- sum(!region_data$edited & inside)
            edited_outside <- edited_total - edited_inside
            unedited_outside <- unedited_total - unedited_inside
            fisher <- fisher_greater(
                edited_inside, edited_outside, unedited_inside, unedited_outside
            )
            data.table(
                species = species, library = as.character(row$group_id),
                genotype = as.character(row$genotype), sex = as.character(row$sex),
                flank_size = as.integer(row$flank_size),
                hairpin_status = hairpin_status,
                event = event,
                n_tested_all = nrow(e),
                n_region_tested = nrow(region_data),
                region_fraction_of_tested = safe_fraction(nrow(region_data), nrow(e)),
                edited_total_region_tested = edited_total,
                edited_inside_event = edited_inside,
                edited_outside_event = edited_outside,
                edited_event_rate = safe_fraction(edited_inside, edited_total),
                unedited_total_region_tested = unedited_total,
                unedited_inside_event = unedited_inside,
                unedited_outside_event = unedited_outside,
                unedited_event_rate = safe_fraction(unedited_inside, unedited_total),
                event_total_region = edited_inside + unedited_inside,
                event_edited_fraction = safe_fraction(
                    edited_inside, edited_inside + unedited_inside
                ),
                event_unedited_fraction = safe_fraction(
                    unedited_inside, edited_inside + unedited_inside
                ),
                event_rate_all_region_tested = safe_fraction(
                    edited_inside + unedited_inside, nrow(region_data)
                ),
                event_rate_all_tested = safe_fraction(
                    edited_inside + unedited_inside, nrow(e)
                ),
                event_odds_ratio_edited_vs_unedited = fisher[["odds_ratio"]],
                event_fisher_p_greater = fisher[["p_value"]],
                hairpin_intervals_raw = interval_qc$hairpin_intervals_raw,
                hairpin_intervals_merged = interval_qc$hairpin_intervals_merged,
                hairpin_transcripts = interval_qc$hairpin_transcripts,
                mapping_site_ids_available = length(editing_maps[[species]]$site_ids)
            )
        })
        rbindlist(event_rows, fill = TRUE)
    }

    counts_in <- summarize_region(hairpin_in, "in")
    counts_out <- summarize_region(hairpin_out, "out")

    diagnostic <- cbind(
        data.table(
            species = species, library = as.character(row$group_id),
            genotype = as.character(row$genotype), sex = as.character(row$sex),
            flank_size = as.integer(row$flank_size),
            n_tested_before_A_filter = n_tested_before_A_filter,
            n_unedited_nonA_removed = n_unedited_nonA_removed,
            n_edited_without_A_retained = n_edited_without_A_retained,
            codon_column_used = "fasta.edited.tab.split",
            n_tested_all = nrow(e),
            n_hairpin_in_tested = nrow(hairpin_in),
            n_hairpin_out_tested = nrow(hairpin_out),
            edited_total_hairpin_in_tested = sum(hairpin_in$edited),
            unedited_total_hairpin_in_tested = sum(!hairpin_in$edited),
            edited_total_hairpin_out_tested = sum(hairpin_out$edited),
            unedited_total_hairpin_out_tested = sum(!hairpin_out$edited),
            stalling_hairpin_in_total = sum(hairpin_in$stalling),
            # Original hairpin-only diagnostic names remain aliases of hairpin in.
            n_hairpin_tested = nrow(hairpin_in),
            edited_total_hairpin_tested = sum(hairpin_in$edited),
            unedited_total_hairpin_tested = sum(!hairpin_in$edited),
            stalling_hairpin_total = sum(hairpin_in$stalling),
            stalling_hairpin_out_total = sum(hairpin_out$stalling),
            hairpin_site_ids_with_any_overlap = length(hairpin_site_ids),
            edgeR_path = as.character(row$edgeR_path),
            parameter_path = as.character(row$parameter_path)
        ),
        interval_qc
    )
    message(
        "[done] ", row$group_id,
        " eligible tested=", format(nrow(e), big.mark = ","),
        "; hairpin in=", format(nrow(hairpin_in), big.mark = ","),
        "; hairpin out=", format(nrow(hairpin_out), big.mark = ","),
        "; removed unedited non-A=", format(n_unedited_nonA_removed, big.mark = ",")
    )
    list(
        counts_in = counts_in,
        counts_out = counts_out,
        diagnostic = diagnostic
    )
}

worker <- function(i) {
    tryCatch(
        process_one(i),
        error = function(error) list(
            error = conditionMessage(error),
            group_id = as.character(tasks[i, group_id])
        )
    )
}

message(
    "Complete samples=", nrow(tasks), ", jobs=", jobs,
    ", data.table threads/worker=", threads
)
jobs <- min(jobs, nrow(tasks))
results <- if (jobs > 1L) {
    parallel::mclapply(
        seq_len(nrow(tasks)), worker,
        mc.preschedule = FALSE, mc.cores = jobs
    )
} else {
    lapply(seq_len(nrow(tasks)), worker)
}
errors <- results[
    vapply(results, function(x) !is.null(x$error), logical(1L))
]
if (length(errors)) {
    stop(
        paste(
            vapply(
                errors, function(x) paste0(x$group_id, ": ", x$error),
                character(1L)
            ),
            collapse = "\n"
        ),
        call. = FALSE
    )
}

counts_in <- rbindlist(
    lapply(results, function(x) x$counts_in),
    fill = TRUE
)
counts_out <- rbindlist(
    lapply(results, function(x) x$counts_out),
    fill = TRUE
)
diagnostics <- rbindlist(
    lapply(results, function(x) x$diagnostic),
    fill = TRUE
)
setorder(counts_in, species, genotype, sex, library, event)
setorder(counts_out, species, genotype, sex, library, event)
setorder(diagnostics, species, genotype, sex, library)

make_pooled <- function(counts_data, group_cols, label) {
    by_cols <- c(group_cols, "hairpin_status", "event")
    pooled <- counts_data[, .(
        libraries = uniqueN(library),
        n_tested_all = sum(n_tested_all),
        n_region_tested = sum(n_region_tested),
        edited_total_region_tested = sum(edited_total_region_tested),
        edited_inside_event = sum(edited_inside_event),
        unedited_total_region_tested = sum(unedited_total_region_tested),
        unedited_inside_event = sum(unedited_inside_event),
        hairpin_intervals_raw = max(hairpin_intervals_raw),
        hairpin_intervals_merged = max(hairpin_intervals_merged),
        hairpin_transcripts = max(hairpin_transcripts),
        mapping_site_ids_available = max(mapping_site_ids_available)
    ), by = by_cols]
    pooled[, pooling := label]
    pooled[, c("edited_outside_event", "unedited_outside_event") := list(
        edited_total_region_tested - edited_inside_event,
        unedited_total_region_tested - unedited_inside_event
    )]
    pooled[, c(
        "region_fraction_of_tested", "edited_event_rate",
        "unedited_event_rate", "event_total_region"
    ) := list(
        n_region_tested / n_tested_all,
        edited_inside_event / edited_total_region_tested,
        unedited_inside_event / unedited_total_region_tested,
        edited_inside_event + unedited_inside_event
    )]
    pooled[, c(
        "event_edited_fraction", "event_unedited_fraction",
        "event_rate_all_region_tested", "event_rate_all_tested"
    ) := list(
        edited_inside_event / event_total_region,
        unedited_inside_event / event_total_region,
        event_total_region / n_region_tested,
        event_total_region / n_tested_all
    )]
    fisher_values <- lapply(seq_len(nrow(pooled)), function(i) {
        fisher_greater(
            pooled$edited_inside_event[[i]],
            pooled$edited_outside_event[[i]],
            pooled$unedited_inside_event[[i]],
            pooled$unedited_outside_event[[i]]
        )
    })
    pooled[, event_odds_ratio_edited_vs_unedited := vapply(
        fisher_values, function(x) x[["odds_ratio"]], numeric(1L)
    )]
    pooled[, event_fisher_p_greater := vapply(
        fisher_values, function(x) x[["p_value"]], numeric(1L)
    )]
    pooled
}

make_wide <- function(counts_data) {
    id_cols <- c(
        "species", "library", "genotype", "sex", "flank_size",
        "hairpin_status"
    )
    wide_metrics <- setdiff(names(counts_data), c(id_cols, "event"))
    wide_stalling <- counts_data[event == "stalling",
        c(id_cols, wide_metrics),
        with = FALSE
    ]
    setnames(wide_stalling, wide_metrics, paste0(wide_metrics, "_stalling"))
    answer <- wide_stalling
    setorder(answer, species, genotype, sex, library)
    answer
}

write_region_bundle <- function(counts_data, region_name) {
    pooled_species <- make_pooled(counts_data, "species", "species")
    pooled_species_sex <- make_pooled(
        counts_data, c("species", "sex"), "species_sex"
    )
    pooled_all <- make_pooled(counts_data, character(), "all_libraries")
    counts_wide <- make_wide(counts_data)

    prefix <- paste0(
        "edited_unedited_hairpin_", region_name, "_peaks"
    )
    fwrite(counts_data, file.path(
        output_dir, paste0(prefix, "_by_library.long.tsv")
    ), sep = "\t")
    fwrite(counts_wide, file.path(
        output_dir, paste0(prefix, "_by_library.tsv")
    ), sep = "\t")
    fwrite(pooled_species, file.path(
        output_dir, paste0(prefix, "_by_species.tsv")
    ), sep = "\t")
    fwrite(pooled_species_sex, file.path(
        output_dir, paste0(prefix, "_by_species_sex.tsv")
    ), sep = "\t")
    fwrite(pooled_all, file.path(
        output_dir, paste0(prefix, "_all_libraries.tsv")
    ), sep = "\t")

    list(
        counts_wide = counts_wide,
        pooled_species = pooled_species,
        pooled_species_sex = pooled_species_sex,
        pooled_all = pooled_all
    )
}

in_outputs <- write_region_bundle(counts_in, "in")
out_outputs <- write_region_bundle(counts_out, "out")

# Preserve the original hairpin-only filenames as aliases of the explicit
# hairpin-in result so existing downstream commands continue to work.
legacy_hairpin_names <- function(dt) {
    answer <- copy(dt)
    rename_map <- c(
        n_region_tested = "n_hairpin_tested",
        region_fraction_of_tested = "hairpin_fraction_of_tested",
        edited_total_region_tested = "edited_total_hairpin_tested",
        unedited_total_region_tested = "unedited_total_hairpin_tested",
        event_total_region = "event_total_hairpin",
        event_rate_all_region_tested = "event_rate_all_hairpin_tested"
    )
    for (old_name in names(rename_map)) {
        hit_names <- names(answer)[
            names(answer) == old_name |
                startsWith(names(answer), paste0(old_name, "_"))
        ]
        if (length(hit_names) > 0L) {
            replacement_names <- sub(
                paste0("^", old_name),
                unname(rename_map[[old_name]]),
                hit_names
            )
            setnames(answer, hit_names, replacement_names)
        }
    }
    if ("hairpin_status" %in% names(answer)) {
        answer[, hairpin_status := NULL]
    }
    answer
}

fwrite(legacy_hairpin_names(counts_in), file.path(
    output_dir, "edited_unedited_hairpin_peaks_by_library.long.tsv"
), sep = "\t")
fwrite(legacy_hairpin_names(in_outputs$counts_wide), file.path(
    output_dir, "edited_unedited_hairpin_peaks_by_library.tsv"
), sep = "\t")
fwrite(legacy_hairpin_names(in_outputs$pooled_species), file.path(
    output_dir, "edited_unedited_hairpin_peaks_by_species.tsv"
), sep = "\t")
fwrite(legacy_hairpin_names(in_outputs$pooled_species_sex), file.path(
    output_dir, "edited_unedited_hairpin_peaks_by_species_sex.tsv"
), sep = "\t")
fwrite(legacy_hairpin_names(in_outputs$pooled_all), file.path(
    output_dir, "edited_unedited_hairpin_peaks_all_libraries.tsv"
), sep = "\t")

fwrite(diagnostics, file.path(
    output_dir, "edited_unedited_hairpin_in_out_peaks_diagnostics.tsv"
), sep = "\t")
fwrite(diagnostics, file.path(
    output_dir, "edited_unedited_hairpin_peaks_diagnostics.tsv"
), sep = "\t")
fwrite(rbindlist(lapply(editing_maps, function(x) x$diagnostic), fill = TRUE),
    file.path(output_dir, "editing_map_diagnostics.tsv"),
    sep = "\t"
)
fwrite(rbindlist(lapply(hairpin_maps, function(x) x$diagnostic), fill = TRUE),
    file.path(output_dir, "hairpin_interval_diagnostics.tsv"),
    sep = "\t"
)

metadata <- data.table(
    parameter = c(
        "result_dir", "run_label", "n_samples", "jobs", "threads_per_worker",
        "manifest", "annotation_dir", "editing_dir", "hairpin_dir",
        "hairpin_coordinate_system", "hairpin_overlap_rule",
        "unedited_background", "edited_background", "output_strata",
        "event_denominator", "timestamp"
    ),
    value = c(
        result_dir, run_label, nrow(tasks), jobs, threads, manifest_path,
        annotation_dir, editing_dir, hairpin_dir, "CDS nucleotide coordinates",
        "codon CDS interval [codon_start, codon_start+2] overlaps hairpin interval",
        "edgeR-tested unedited codons containing at least one A",
        "all mapped edited codons retained",
        "hairpin in and hairpin out are summarized separately",
        "edgeR v2 tested sites within each hairpin stratum",
        as.character(Sys.time())
    )
)
fwrite(metadata, file.path(output_dir, "run_metadata.tsv"), sep = "\t")

invisible(NULL)

message(
    "Wrote hairpin-in and hairpin-out summary tables to ",
    normalizePath(output_dir)
)
print(rbindlist(list(
    in_outputs$pooled_all,
    out_outputs$pooled_all
), fill = TRUE))
