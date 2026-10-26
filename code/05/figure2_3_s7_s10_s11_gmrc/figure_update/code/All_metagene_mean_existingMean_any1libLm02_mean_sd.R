#!/usr/bin/env Rscript

# Metagene profiles from codon-level replicate means.
# Editing centers require RNA Lm >= 0.2 in at least one species-matched library.
# Profiles summarize sites after trimming the two lowest and two highest values.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
    sub("^--file=", "", script_arg[[1L]])
} else {
    file.path(getwd(), "Dmel_metagene_mean_existingMean_any1libLm02.R")
}
script_dir <- normalizePath(dirname(script_file), mustWork = TRUE)

resolve_output_path <- function(path) {
    if (grepl("^/", path)) path else file.path(script_dir, path)
}
outdir <- resolve_output_path(
    Sys.getenv("METAGENE_OUTDIR", unset = "metagene_gmRC")
)
mean_dir <- Sys.getenv("METAGENE_MEAN_DIR", file.path(script_dir, "species_sex_mean_gmRC_tables"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

mean_table_dir <- function() {
    if (!dir.exists(mean_dir)) {
        stop("Missing mean-gmRC directory: ", mean_dir)
    }
    mean_dir
}

alledit_dir <- file.path(
    script_dir, "..", "04_alledit", "result",
    "alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed_newdata",
    "results"
)
alledit_dir <- normalizePath(alledit_dir, mustWork = TRUE)

cdspos_dir <- file.path(
    script_dir, "..", "03_cdspos", "result",
    "cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed"
)
cdspos_dir <- normalizePath(Sys.getenv("METAGENE_CDSPOS_DIR", cdspos_dir), mustWork = TRUE)

trans_files <- c(
    Dmel = "/path/to/Proj2_allele_specific_translation/02_codon/dmel-r6.12_CDS_position",
    Dsim = "/path/to/Proj1_RNAediting/03_APEsite/dsim_r2.02_CDS_position",
    Dsec = "/path/to/Proj1_RNAediting/03_APEsite/dsec_r1.3_CDS_position",
    Dere = "/path/to/Proj1_RNAediting/03_APEsite/dere_r1.05_CDS_position"
)

m6a_files <- c(
    Dmel = "/path/to/Proj1_editing/codex/02_stalling/3flank/peak_find/m6A_list.handle",
    Dsim = "/path/to/Proj1_editing/codex/02_stalling/3flank/peak_find/m6A_list.droSim1.handle",
    Dsec = "/path/to/Proj1_editing/codex/02_stalling/3flank/peak_find/m6A_list.droSec1ToDm6.handle",
    Dere = "/path/to/Proj1_editing/codex/02_stalling/3flank/peak_find/m6A_list.droEre2ToDm6.handle"
)

reference_config <- Sys.getenv("METAGENE_REFERENCE_CONFIG", "")
if (nzchar(reference_config)) {
    references <- fread(reference_config)
    stopifnot(
        all(c("species", "trans_file", "m6a_file") %in% names(references)),
        !anyDuplicated(references$species),
        setequal(references$species, names(trans_files))
    )
    trans_files <- setNames(references$trans_file, references$species)
    m6a_files <- setNames(references$m6a_file, references$species)
}
if (any(!file.exists(m6a_files))) {
    stop("Missing m6A reference file(s):\n", paste(m6a_files[!file.exists(m6a_files)], collapse = "\n"))
}
m6a_by_species <- setNames(lapply(m6a_files, function(path) {
    unique(scan(path, what = character(), quiet = TRUE))
}), names(m6a_files))

cdspos_files <- c(
    Dmel = "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    Dsim = "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    Dsec = "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    Dere = "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
)

species_samples <- list(
    Dmel = c(
        "Dmel_ISO1_head_female_ribo", "Dmel_ISO1_head_male_ribo",
        "WBg_F1", "WBg_F2", "WBg_F3", "WBg_M1", "WBg_M2", "WBg_M3"
    ),
    Dsim = c(
        "Dsim_F1", "Dsim_F2", "Dsim_F3", "Dsim_female_head_ribo",
        "Dsim_M1", "Dsim_M2", "Dsim_M3", "Dsim_male_head_ribo"
    ),
    Dsec = c(
        "Dsec_F1", "Dsec_F2", "Dsec_F3",
        "Dsec_M1", "Dsec_M2", "Dsec_M3", "Dsec_male_head_ribo"
    ),
    Dere = c(
        "Dere_F1", "Dere_F2", "Dere_F3",
        "Dere_M1", "Dere_M2", "Dere_M3", "Dere_male_head_ribo"
    )
)


unit_info <- data.table(
    unit_id = c(
        "Dmel_ISO1__female", "Dmel_ISO1__male",
        "Dmel_WBg__female", "Dmel_WBg__male",
        "Dmel_ISO1_WBg__female", "Dmel_ISO1_WBg__male",
        "Dsim__female", "Dsim__male",
        "Dsec__female", "Dsec__male",
        "Dere__female", "Dere__male",
        "Dmel_Adar5G1__male"
    ),
    species = c(
        rep("Dmel", 6L), rep("Dsim", 2L), rep("Dsec", 2L),
        rep("Dere", 2L), "Dmel"
    ),
    cohort = c(
        "Dmel_ISO1", "Dmel_ISO1", "Dmel_WBg", "Dmel_WBg",
        "Dmel_ISO1_WBg", "Dmel_ISO1_WBg",
        "Dsim", "Dsim", "Dsec", "Dsec", "Dere", "Dere",
        "Dmel_Adar5G1"
    ),
    sex = c(
        "female", "male", "female", "male",
        "female", "male",
        "female", "male", "female", "male", "female", "male", "male"
    ),
    label = c(
        "Dmel ISO1, female", "Dmel ISO1, male",
        "Dmel WBg, female", "Dmel WBg, male",
        "Dmel ISO1 + WBg, female", "Dmel ISO1 + WBg, male",
        "Dsim, female", "Dsim, male",
        "Dsec, female", "Dsec, male",
        "Dere, female", "Dere, male",
        "Dmel Adar5G1, male"
    ),
    edit_species = c(
        rep("Dmel", 6L), "Dsim", "Dsim", "Dsec", "Dsec",
        "Dere", "Dere", "Dmel"
    )
)

# Additional analysis: merge Dmel ISO1 and WBg by retaining only the combined
# Dmel female/male units. Other species and the Adar5G1 male control are unchanged.
unit_info <- unit_info[unit_id %chin% c(
    "Dmel_ISO1_WBg__female", "Dmel_ISO1_WBg__male",
    "Dsim__female", "Dsim__male", "Dsec__female", "Dsec__male",
    "Dere__female", "Dere__male", "Dmel_Adar5G1__male"
)]
stopifnot(nrow(unit_info) == 9L)

# Use only the mean tables already produced upstream.  The combined Dmel
# ISO1+WBg units are the main Dmel female/male mean tables; every other unit
# has a mean-table filename matching its unit_id.
unit_info[, mean_file := file.path(
    mean_table_dir(),
    fcase(
        unit_id == "Dmel_ISO1_WBg__female", "Dmel__female.mean_gmRC.tsv.gz",
        unit_id == "Dmel_ISO1_WBg__male",   "Dmel__male.mean_gmRC.tsv.gz",
        default = paste0(unit_id, ".mean_gmRC.tsv.gz")
    )
)]
if (any(!file.exists(unit_info$mean_file))) {
    stop(
        "Missing precomputed mean-gmRC table(s):\n",
        paste(unit_info[!file.exists(mean_file), mean_file], collapse = "\n")
    )
}
desired_order <- c(
    "Dmel_Adar5G1__male",
    "Dmel_ISO1__female", "Dmel_ISO1__male",
    "Dmel_WBg__female", "Dmel_WBg__male",
    "Dmel_ISO1_WBg__female", "Dmel_ISO1_WBg__male",
    "Dsim__female", "Dsim__male", "Dsec__female", "Dsec__male",
    "Dere__female", "Dere__male"
)
unit_info[, unit_order := match(unit_id, desired_order)]
setorder(unit_info, unit_order)
unit_info[, unit_order := NULL]
fwrite(unit_info, file.path(outdir, "metagene_units.tsv"), sep = "\t")

read_transcript_lengths <- function(species) {
    file <- trans_files[[species]]
    if (!file.exists(file)) stop("Missing transcript-length file: ", file)
    x <- fread(
        file,
        header = FALSE,
        col.names = c("FBtr", "start", "end", "cds")
    )
    x[, `:=`(
        FBtr = as.character(FBtr),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end))
    )]
    x <- x[!is.na(FBtr) & FBtr != "" & is.finite(start) & is.finite(end)]
    x[, codonlength := (end + 1 - start) / 3]
    x[is.finite(codonlength), .(codonlength = max(codonlength)), by = FBtr]
}

read_candidate_sites <- function(species) {
    file <- file.path(cdspos_dir, cdspos_files[[species]])
    if (!file.exists(file)) stop("Missing cdspos file: ", file)
    x <- fread(file, select = c("new_transcript", "new_cds_pos_1based"))
    x[, `:=`(
        FBtr = as.character(new_transcript),
        cds_pos = suppressWarnings(as.integer(new_cds_pos_1based))
    )]
    x <- unique(x[!is.na(FBtr) & FBtr != "" & is.finite(cds_pos) & cds_pos >= 1L])
    x[, `:=`(
        site_key = paste(FBtr, cds_pos, sep = "_"),
        codon_pos = (cds_pos - 1L) %/% 3L + 1L,
        FBtr_codonpos = paste(
            FBtr, (cds_pos - 1L) %/% 3L + 1L,
            sep = "_"
        )
    )]
    x[]
}

read_lm_site_set <- function(sample_name, candidate_sites) {
    file <- file.path(alledit_dir, paste0(sample_name, ".Asite_df_alledit"))
    if (!file.exists(file)) stop("Missing A-site file: ", file)
    x <- fread(
        file,
        select = c("FBt_cdspos", "ref_number_rna", "alt_number_rna"),
        na.strings = c("", "NA", "na", "NaN", "nan")
    )
    x[, `:=`(
        FBt_cdspos = as.character(FBt_cdspos),
        ref_number_rna = suppressWarnings(as.numeric(ref_number_rna)),
        alt_number_rna = suppressWarnings(as.numeric(alt_number_rna))
    )]
    x <- x[
        FBt_cdspos %chin% candidate_sites &
            is.finite(ref_number_rna) & is.finite(alt_number_rna) &
            (ref_number_rna + alt_number_rna) > 0
    ]
    unique(x[alt_number_rna / (ref_number_rna + alt_number_rna) >= 0.2, FBt_cdspos])
}

candidate_by_species <- setNames(
    lapply(names(cdspos_files), read_candidate_sites),
    names(cdspos_files)
)
known_edit_by_species <- lapply(candidate_by_species, function(x) {
    unique(x$FBtr_codonpos)
})

any_edit_by_species <- lapply(names(candidate_by_species), function(species) {
    candidate <- candidate_by_species[[species]]
    candidate_sites <- unique(candidate$site_key)
    sample_sets <- lapply(
        species_samples[[species]],
        read_lm_site_set,
        candidate_sites = candidate_sites
    )
    any_sites <- unique(unlist(sample_sets, use.names = FALSE))
    any_codons <- unique(
        candidate[site_key %chin% any_sites, FBtr_codonpos]
    )
    message(
        species, ": ", length(candidate_sites), " candidate sites; ",
        length(any_sites), " any-library Lm>=0.2 sites; ",
        length(any_codons), " any-library edit codon centers"
    )
    any_codons
})
names(any_edit_by_species) <- names(candidate_by_species)

make_flank_set <- function(codon_keys, radius = 5L) {
    out <- unique(unlist(lapply(codon_keys, function(key) {
        pieces <- strsplit(key, "_", fixed = TRUE)[[1L]]
        if (length(pieces) != 2L) return(character())
        positions <- as.integer(pieces[[2L]]) + seq(-radius, radius)
        paste(pieces[[1L]], positions, sep = "_")
    })))
    out[!is.na(out) & out != ""]
}

read_mean_table <- function(file) {
    x <- fread(
        file,
        select = c("FBtr_codonpos", "FBtr", "codon", "AA", "gmRC", "chr", "pos"),
        na.strings = c("", "NA", "na", "NaN", "nan")
    )
    x[, `:=`(
        FBtr_codonpos = as.character(FBtr_codonpos),
        FBtr = as.character(FBtr),
        codon = toupper(as.character(codon)),
        AA = as.character(AA),
        gmRC = suppressWarnings(as.numeric(gmRC)),
        chr = as.character(chr),
        pos = suppressWarnings(as.numeric(pos))
    )]
    x <- x[!is.na(FBtr_codonpos) & FBtr_codonpos != ""]
    x[, codon_pos := suppressWarnings(as.integer(sub(".*_", "", FBtr_codonpos)))]
    x[, FBtr := fifelse(
        !is.na(FBtr) & FBtr != "",
        FBtr,
        sub("_.*", "", FBtr_codonpos)
    )]
    x[, handle := fifelse(
        !is.na(chr) & chr != "" & is.finite(pos),
        paste(chr, format(pos, scientific = FALSE, trim = TRUE), sep = "_"),
        NA_character_
    )]
    x <- unique(x, by = "FBtr_codonpos")
    x[]
}

trimmed_values <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) <= 4L) return(numeric())
    x <- sort(x)
    x[3L:(length(x) - 2L)]
}

summarize_vector <- function(x) {
    finite_x <- as.numeric(x)
    finite_x <- finite_x[is.finite(finite_x)]
    trimmed_x <- trimmed_values(finite_x)
    c(
        mean = if (length(trimmed_x)) mean(trimmed_x) else NA_real_,
        sd = if (length(trimmed_x) > 1L) sd(trimmed_x) else NA_real_,
        n_before_trim = length(finite_x),
        n_after_trim = length(trimmed_x)
    )
}

summarize_matrix <- function(x) {
    if (nrow(x) == 0L) {
        return(list(
            mean = rep(NA_real_, ncol(x)),
            sd = rep(NA_real_, ncol(x)),
            n_before_trim = rep(0L, ncol(x)),
            n_after_trim = rep(0L, ncol(x))
        ))
    }
    summary_mat <- vapply(
        seq_len(ncol(x)), function(j) summarize_vector(x[, j]), numeric(4L)
    )
    list(
        mean = as.numeric(summary_mat["mean", ]),
        sd = as.numeric(summary_mat["sd", ]),
        n_before_trim = as.integer(summary_mat["n_before_trim", ]),
        n_after_trim = as.integer(summary_mat["n_after_trim", ])
    )
}

process_unit <- function(info) {
    species <- info$species[[1L]]
    edit_species <- info$edit_species[[1L]]
    mean_dt <- read_mean_table(info$mean_file[[1L]])
    m6a_handles <- m6a_by_species[[species]]
    mean_dt[, is_m6a := !is.na(handle) & handle %chin% m6a_handles]
    n_mean_codons_before_m6a <- nrow(mean_dt)
    m6a_codon_keys <- unique(mean_dt[is_m6a == TRUE, FBtr_codonpos])
    n_m6a_codons_removed <- length(m6a_codon_keys)
    mean_dt <- mean_dt[is_m6a == FALSE]
    mean_dt[, is_m6a := NULL]
    length_dt <- read_transcript_lengths(species)
    mean_dt <- merge(mean_dt, length_dt, by = "FBtr", all = FALSE, sort = FALSE)
    mean_dt <- mean_dt[
        is.finite(gmRC) & is.finite(codonlength) &
            codonlength > 150 & codon_pos > 15 &
            (codonlength - codon_pos) > 5
    ]
    if (nrow(mean_dt) == 0L) {
        stop("No finite, edge-filtered mean-gmRC rows for ", info$unit_id[[1L]])
    }

    map_names <- mean_dt$FBtr_codonpos
    map_values <- mean_dt$gmRC
    keep_map <- !duplicated(map_names)
    map_names <- map_names[keep_map]
    map_values <- map_values[keep_map]

    known_edit <- known_edit_by_species[[edit_species]]
    known_edit_flank5 <- make_flank_set(known_edit, radius = 5L)
    edit_dt <- data.table(
        FBtr_codonpos = unique(any_edit_by_species[[edit_species]])
    )
    edit_dt[, `:=`(
        FBtr = sub("_.*", "", FBtr_codonpos),
        codon_pos = suppressWarnings(as.integer(sub(".*_", "", FBtr_codonpos)))
    )]
    n_edit_centers_before_m6a <- nrow(edit_dt)
    edit_dt <- edit_dt[!(FBtr_codonpos %chin% m6a_codon_keys)]
    n_m6a_edit_centers_removed <- n_edit_centers_before_m6a - nrow(edit_dt)
    edit_dt <- merge(
        edit_dt,
        unique(length_dt[, .(FBtr, codonlength)]),
        by = "FBtr", all.x = TRUE, sort = FALSE
    )
    edit_dt <- edit_dt[
        is.finite(codon_pos) & is.finite(codonlength) &
            codonlength > 150 & codon_pos > 15 &
            (codonlength - codon_pos) > 5
    ]
    if (nrow(edit_dt) == 0L) {
        stop("No eligible any-library edit centers for ", info$unit_id[[1L]])
    }

    offsets <- -100:100
    n_offsets <- length(offsets)
    process_edited_row <- function(row) {
        positions <- as.integer(row$codon_pos[[1L]]) + offsets
        identifiers <- paste0(row$FBtr[[1L]], "_", positions)
        hit_map <- match(identifiers, map_names)
        result <- rep(NA_real_, n_offsets)
        valid <- !is.na(hit_map) & !(identifiers %chin% known_edit)
        result[valid] <- map_values[hit_map[valid]]
        center <- match(
            paste0(row$FBtr[[1L]], "_", row$codon_pos[[1L]]), map_names
        )
        if (!is.na(center)) result[101L] <- map_values[[center]]
        result
    }
    edited_matrix <- matrix(
        unlist(lapply(seq_len(nrow(edit_dt)), function(i) {
            process_edited_row(edit_dt[i])
        })),
        nrow = nrow(edit_dt), byrow = TRUE, ncol = n_offsets
    )

    # Restrict only the control center (codon=0) to unedited reference codons
    # containing A. Once a center is selected, its surrounding profile is not
    # filtered by codon identity; neighboring codons are taken from the full
    # mean-gmRC table, subject to the existing known-edit +/-5 exclusion.
    control_center_dt <- mean_dt[
        is.finite(gmRC) &
            !is.na(codon) &
            !(FBtr_codonpos %chin% known_edit_flank5),
        .(control_positions = list(codon_pos)),
        by = FBtr
    ]
    if (control_center_A_only) {
        control_center_dt <- mean_dt[
            is.finite(gmRC) &
                !is.na(codon) &
                grepl("A", codon, fixed = TRUE) &
                !(FBtr_codonpos %chin% known_edit_flank5),
            .(control_positions = list(codon_pos)),
            by = FBtr
        ]
    }
    control_center_positions <- setNames(
        control_center_dt$control_positions,
        control_center_dt$FBtr
    )

    process_control_row <- function(row) {
        eligible_positions <- control_center_positions[[row$FBtr[[1L]]]]
        if (is.null(eligible_positions) || length(eligible_positions) == 0L) {
            return(matrix(NA_real_, nrow = 20L, ncol = n_offsets))
        }
        random_positions <- sample(eligible_positions, size = 20L, replace = TRUE)
        position_matrix <- outer(random_positions, offsets, "+")
        identifiers <- paste0(row$FBtr[[1L]], "_", as.vector(position_matrix))
        hit_map <- match(identifiers, map_names)
        result <- rep(NA_real_, length(identifiers))
        valid <- !is.na(hit_map) & !(identifiers %chin% known_edit_flank5)
        result[valid] <- map_values[hit_map[valid]]
        # outer() is flattened column-wise: all sampled controls for one
        # relative position are contiguous. Keep that order when reshaping.
        matrix(result, nrow = 20L, ncol = n_offsets)
    }
    control_matrix <- do.call(rbind, lapply(seq_len(nrow(edit_dt)), function(i) {
        process_control_row(edit_dt[i])
    }))

    edited_summary <- summarize_matrix(edited_matrix)
    control_summary <- summarize_matrix(control_matrix)
    result <- data.table(
        unit_id = info$unit_id[[1L]],
        species = species,
        cohort = info$cohort[[1L]],
        sex = info$sex[[1L]],
        label = info$label[[1L]],
        Position = offsets,
        edited = edited_summary$mean,
        control = control_summary$mean,
        edited_sd = edited_summary$sd,
        control_sd = control_summary$sd,
        edited_n_before_trim = edited_summary$n_before_trim,
        control_n_before_trim = control_summary$n_before_trim,
        edited_n_after_trim = edited_summary$n_after_trim,
        control_n_after_trim = control_summary$n_after_trim,
        n_edit_centers = nrow(edit_dt),
        n_mean_codons_before_m6a = n_mean_codons_before_m6a,
        n_m6a_codons_removed = n_m6a_codons_removed,
        n_m6a_edit_centers_removed = n_m6a_edit_centers_removed
    )
    rm(mean_dt, edited_matrix, control_matrix)
    invisible(gc())
    result
}

set.seed(18907675L)
existing_meta_file <- file.path(outdir, "metagene_gmRC_summary.tsv")
control_center_A_only <- tolower(
    Sys.getenv("METAGENE_CONTROL_CENTER_A_ONLY", unset = "false")
) %in% c("1", "true", "yes", "y")
control_filter_version <- paste(
    if (control_center_A_only) {
        "control_center=reference_codon_contains_A"
    } else {
        "control_center=all_reference_codons"
    },
    "flanks=unrestricted_except_known_edit_flank5",
    "m6A=excluded_from_centers_and_all_profile_positions",
    sep = ";"
)
control_filter_file <- file.path(outdir, "metagene_control_center_filter.txt")
cache_is_current <- file.exists(control_filter_file) &&
    identical(
        readLines(control_filter_file, warn = FALSE),
        control_filter_version
    )
existing_meta <- if (file.exists(existing_meta_file)) {
    fread(existing_meta_file, na.strings = c("", "NA", "NaN"))
} else {
    NULL
}
meta_results <- lapply(seq_len(nrow(unit_info)), function(i) {
    current_unit_id <- unit_info$unit_id[[i]]
    cached <- if (cache_is_current && !is.null(existing_meta)) {
        existing_meta[unit_id == current_unit_id]
    } else {
        NULL
    }
    cache_has_mean_sd <- !is.null(cached) && all(c(
        "edited_sd", "control_sd",
        "edited_n_after_trim", "control_n_after_trim"
    ) %in% names(cached))
    if (cache_has_mean_sd && nrow(cached) == length(-100:100)) {
        message(
            "Reusing completed metagene unit ", i, "/", nrow(unit_info), ": ",
            current_unit_id
        )
        cached
    } else {
        message(
            "Processing metagene unit ", i, "/", nrow(unit_info), ": ",
            current_unit_id
        )
        process_unit(unit_info[i])
    }
})
meta_dt <- rbindlist(meta_results, use.names = TRUE, fill = TRUE)
meta_dt[, unit_order := match(unit_id, unit_info$unit_id)]
setorder(meta_dt, unit_order, Position)
meta_dt[, unit_order := NULL]

fwrite(
    meta_dt,
    file.path(outdir, "metagene_gmRC_summary.tsv"),
    sep = "\t", na = "NA"
)
fwrite(
    unique(meta_dt[, .(
        unit_id, species, cohort, sex,
        n_mean_codons_before_m6a, n_m6a_codons_removed,
        n_m6a_edit_centers_removed, n_edit_centers
    )]),
    file.path(outdir, "metagene_m6A_filter_QC.tsv"),
    sep = "\t", na = "NA"
)
fwrite(
    meta_dt[, .(
        unit_id, species, cohort, sex, label, Position,
        edited_n_before_trim, control_n_before_trim,
        edited_n_after_trim, control_n_after_trim, n_edit_centers
    )],
    file.path(outdir, "metagene_gmRC_site_count_before_trim.tsv"),
    sep = "\t", na = "NA"
)
saveRDS(
    split(meta_dt, meta_dt$unit_id),
    file.path(outdir, "metagene_gmRC_meta.rds"),
    version = 2
)

plot_data <- melt(
    meta_dt,
    id.vars = c(
        "unit_id", "species", "cohort", "sex", "label", "Position",
        "edited_sd", "control_sd",
        "edited_n_before_trim", "control_n_before_trim",
        "edited_n_after_trim", "control_n_after_trim", "n_edit_centers",
        "n_mean_codons_before_m6a", "n_m6a_codons_removed",
        "n_m6a_edit_centers_removed"
    ),
    measure.vars = c("edited", "control"),
    variable.name = "profile",
    value.name = "gmRC"
)
plot_data[, profile := factor(
    profile,
    levels = c("edited", "control"),
    labels = c("Edited sites", "Control sites")
)]

format_count <- function(x) {
    if (length(x) == 0L || !is.finite(x)) return("NA")
    formatC(as.numeric(x), format = "d", big.mark = ",")
}

format_gmrc <- function(x) {
    if (length(x) == 0L || !is.finite(x)) return("NA")
    formatC(as.numeric(x), format = "f", digits = 3L)
}

make_unit_title <- function(info) {
    species_name <- as.character(info$species[[1L]])
    species_plotmath <- if (species_name %in% c("Dsim", "Dsec", "Dere")) {
        sub("^D", "", species_name)
    } else {
        species_name
    }
    cohort_name <- as.character(info$cohort[[1L]])
    cohort_label <- if (cohort_name == "Dmel_ISO1") {
        "ISO1"
    } else if (cohort_name == "Dmel_WBg") {
        "WBg"
    } else if (cohort_name == "Dmel_ISO1_WBg") {
        "ISO1 + WBg"
    } else if (cohort_name == "Dmel_Adar5G1") {
        "Adar5G1"
    } else {
        ""
    }
    sex_label <- if (as.character(info$sex[[1L]]) == "female") "F" else "M"
    title_tail <- paste(c(cohort_label, sex_label)[nzchar(c(cohort_label, sex_label))],
        collapse = " "
    )
    parse(text = paste0(
        "italic(", species_plotmath, ")~'", title_tail, "'"
    ))[[1L]]
}

plot_objects <- lapply(seq_len(nrow(unit_info)), function(i) {
    info <- unit_info[i]
    current_unit_id <- info$unit_id[[1L]]
    d <- plot_data[unit_id == current_unit_id]
    center <- meta_dt[unit_id == current_unit_id & Position == 0][1L]
    center_label <- paste0(
        "Edited: ", format_gmrc(center$edited), " \u00b1 ",
        format_gmrc(center$edited_sd),
        " (n=", format_count(center$edited_n_before_trim), ")\n",
        "bg: ", format_gmrc(center$control), " \u00b1 ",
        format_gmrc(center$control_sd)
    )
    ggplot(d, aes(x = Position, y = gmRC, color = profile)) +
        geom_vline(
            xintercept = 0, linetype = 2, linewidth = 0.45,
            color = "grey45"
        ) +
        geom_line(linewidth = 0.45, na.rm = TRUE) +
        scale_color_manual(values = c(
            "Edited sites" = "#C0392B", "Control sites" = "#2C7FB8"
        ), drop = FALSE) +
        scale_x_continuous(
            breaks = c(-100, -50, 0, 50, 100),
            labels = c("-100", "-50", "0", "50", "100"),
            expand = expansion(mult = c(0.01, 0.01))
        ) +
        scale_y_continuous(
            breaks = seq(0.2, 2.2, by = 0.4),
            expand = expansion(mult = c(0, 0))
        ) +
        annotate(
            "text",
            x = 0, y = 2.5, label = center_label,
            size = 2.64, lineheight = 0.9, hjust = 0.5, vjust = 0.5
        ) +
        coord_cartesian(ylim = c(0.2, 2.2), clip = "off") +
        labs(
            title = make_unit_title(info),
            x = "Codon position relative to edit",
            y = "gmRC",
            color = NULL
        ) +
        theme_classic(base_size = 7.5) +
        theme(
            plot.title = element_text(
                size = 8, face = "bold", hjust = 0.5,
                margin = margin(t = 0, r = 0, b = 20, l = 0)
            ),
            axis.title = element_text(size = 7.5),
            axis.text = element_text(size = 7.5, color = "black"),
            axis.text.x = element_text(size = 7.5),
            axis.ticks = element_line(linewidth = 0.3, color = "black"),
            axis.line = element_line(linewidth = 0.35, color = "black"),
            legend.position = "bottom",
            legend.text = element_text(size = 7.5),
            legend.key.width = grid::unit(0.65, "cm"),
            legend.key.height = grid::unit(0.3, "cm"),
            plot.margin = margin(6, 3, 4, 3)
        )
})
names(plot_objects) <- unit_info$unit_id

center_stats <- meta_dt[Position == 0, .(
    unit_id, species, cohort, sex, label, Position,
    edited_mean = edited,
    edited_sd,
    edited_n_sites = edited_n_before_trim,
    edited_n_after_trim,
    bg_mean = control,
    bg_sd = control_sd,
    bg_n_observations = control_n_before_trim,
    bg_n_after_trim = control_n_after_trim,
    n_edit_centers
)]
fwrite(
    center_stats,
    file.path(outdir, "metagene_gmRC_center_mean_sd.tsv"),
    sep = "\t", na = "NA"
)

metagene_pdf <- file.path(outdir, "metagene_gmRC_9_units.pdf")
pdf(
    metagene_pdf,
    width = 8.27, height = 11.69, onefile = TRUE, useDingbats = FALSE
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(
    y = grid::unit(0.5, "npc"),
    height = grid::unit(0.94, "npc"),
    layout = grid::grid.layout(
        6L, 3L,
        heights = grid::unit(c(0.18, 0.92, 0.92, 0.92, 0.92, 0.92), "null")
    ),
    gp = grid::gpar(fontfamily = "sans")
))
for (i in seq_along(plot_objects)) {
    row <- (i - 1L) %/% 3L + 2L
    col <- (i - 1L) %% 3L + 1L
    print(
        plot_objects[[i]],
        vp = grid::viewport(layout.pos.row = row, layout.pos.col = col)
    )
}
grid::popViewport()
dev.off()
file.copy(
    metagene_pdf,
    file.path(outdir, "metagene_gmRC_9_units_alias.pdf"),
    overwrite = TRUE
)
writeLines(control_filter_version, control_filter_file)

message("Done. Metagene outputs written to: ", normalizePath(outdir))
