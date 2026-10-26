#!/usr/bin/env Rscript

# ============================================================
# WT-only group-mean gmRC mixed-model analysis
# WBg F/M: 3 biological replicates per sex
#
# Key differences from the original per-replicate script:
#   1) Replicates are combined BEFORE coverage filtering/model fitting.
#   2) For each codon:
#        mRNA_sum >= 10 * n_replicates
#        at least one replicate has Asite_median > 0
#        gmRC = arithmetic mean across replicates
#        logTE = arithmetic mean across replicate-specific logTE values
#   3) WBg editing level = arithmetic mean of replicate-level editing
#      levels within the same sex.
#
# All other codon filters and Full/Matched model logic follow the supplied
# merged-TSV sites_final script.
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(lme4)
    library(lmerTest)
})

# ============================================================
# 1. Paths and parameters
# ============================================================

phyloP_dir <- paste0("/path/to/RNA_editing_DATAprocessing/", "05.psite/psite_2times/reproduce_old_flow/peak/",
    "codon_AA/phyloP_annotation")

output_dir <- Sys.getenv("GMRC_OUTPUT_DIR", unset = paste0("/path/to/RNA_editing_DATAprocessing_new/10.stalling/04_GLM/",
    "gmRC_groupMean_WBg_WT_only_m6A_excluded"))

trans_file <- paste0("/path/to/Proj2_allele_specific_translation/", "02_codon/dmel-r6.12_CDS_position")

codon_edit_file <- paste0("/path/to/RNA_editing_DATAprocessing/", "ref_genome/dmel-all-CDS-r6.12.fasta.edited.tab.split")

tr_gn_file <- paste0("/path/to/RNA_editing_DATAprocessing/", "05.psite/dmel_r6.12_txinfo.biotype.tsv")

logTE_file_WBg <- paste0("/path/to/RNA_editing_DATAprocessing/", "10.stalling/04_GLM/log2TE_matrix.tsv")

rscu_file <- paste0("/path/to/RNAi/", "240719drosophila_check/RSCU_flybase.txt")

structure_bed_file <- file.path("/path/to/RNA_editing_DATAprocessing_new/", "10.stalling/04_GLM/RNALfold_result",
    "Dmel.fa.T25Z3.2.merged.tsv")

edit_site_dir <- paste0("/path/to/Proj1_editing/code/", "01_identify_sites/sites_final")

# Only WBG editing is required in this version.
wbg_edit_file <- file.path(edit_site_dir, "WBG.candidate.AtoG.annotations.long.tsv.new")

# Exact genomic m6A handles. F3E and FS11 use D. melanogaster only.
m6a_reference_config <- Sys.getenv("M6A_REFERENCE_CONFIG", unset = "")
m6a_file <- Sys.getenv("GLM_M6A_FILE", unset = "")

# Input paths and configuration
path_env <- c(phyloP_dir = "GLM_PHYLOP_DIR", trans_file = "GLM_TRANS_FILE", codon_edit_file = "GLM_CODON_FILE",
    tr_gn_file = "GLM_TX_GENE_FILE", logTE_file_WBg = "GLM_TE_FILE", rscu_file = "GLM_RSCU_FILE",
    structure_bed_file = "GLM_STRUCTURE_FILE", wbg_edit_file = "GLM_EDIT_WBG")
for (path_name in names(path_env)) {
    value <- Sys.getenv(path_env[[path_name]], unset = "")
    if (nzchar(value))
        assign(path_name, value)
}
# END packaged path overrides

RUN_MATCHED <- TRUE
SAVE_MODELS <- FALSE
MAKE_PLOTS <- TRUE

MATCH_RATIO <- 20L
BASE_SEED <- 20260816L

# Group-level coverage criteria:
#   mRNA_sum >= MIN_MRNA_PER_REP * number_of_replicates
MIN_MRNA_PER_REP <- 10

# Expected biological-replicate design.
expected_design <- data.table(group_id = c("WBg_F_mean", "WBg_M_mean"), genotype = c("WBg",
    "WBg"), sex = c("F", "M"), expected_n_reps = c(3L, 3L), editing_source = c("WBg_F_mean_editing",
    "WBg_M_mean_editing"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "per_group"), recursive = TRUE, showWarnings = FALSE)

if (!nzchar(m6a_file)) {
    if (!nzchar(m6a_reference_config) || !file.exists(m6a_reference_config)) {
        stop("Missing M6A_REFERENCE_CONFIG or GLM_M6A_FILE")
    }
    m6a_reference <- fread(m6a_reference_config)
    if (!all(c("species", "m6a_file") %in% names(m6a_reference))) {
        stop("m6A reference config must contain species and m6a_file columns")
    }
    m6a_file <- m6a_reference[species == "Dmel", m6a_file][1L]
}
if (is.na(m6a_file) || !nzchar(m6a_file) || !file.exists(m6a_file)) {
    stop("Missing Dmel m6A handle file: ", m6a_file)
}
m6a_handles <- unique(scan(m6a_file, what = character(), quiet = TRUE))
if (!length(m6a_handles)) stop("Dmel m6A handle file is empty: ", m6a_file)

# ============================================================
# 2. Generic helpers
# ============================================================

normalize_sample_name <- function(x) {
    toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
}

resolve_unique_name <- function(aliases, candidates) {
    candidates <- unique(as.character(candidates))
    candidates <- candidates[!is.na(candidates) & candidates != ""]

    aliases <- unique(as.character(aliases))
    aliases <- aliases[!is.na(aliases) & aliases != ""]

    # Exact match first.
    exact_hits <- candidates[candidates %chin% aliases]
    exact_hits <- unique(exact_hits)

    if (length(exact_hits) == 1L) {
        return(exact_hits)
    }

    # Normalized match next.
    alias_keys <- unique(normalize_sample_name(aliases))
    candidate_keys <- normalize_sample_name(candidates)
    normalized_hits <- unique(candidates[candidate_keys %in% alias_keys])

    if (length(normalized_hits) == 1L) {
        return(normalized_hits)
    }

    return(NA_character_)
}

safe_mean <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok))
        return(NA_real_)
    mean(x[ok])
}

safe_sum <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok))
        return(0)
    sum(x[ok])
}

safe_min <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok))
        return(NA_real_)
    min(x[ok])
}

safe_max <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok))
        return(NA_real_)
    max(x[ok])
}

safe_zscore <- function(x) {
    x <- as.numeric(x)
    x_mean <- mean(x, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)

    if (!is.finite(x_sd) || x_sd == 0) {
        return(rep(NA_real_, length(x)))
    }

    (x - x_mean)/x_sd
}

parse_target_sample <- function(sample_name) {
    # WBg_F1 / WBg_M1 etc.
    if (grepl("^WBg_[FM][0-9]+$", sample_name)) {
        suffix <- sub("^WBg_", "", sample_name)
        sex <- sub("^([FM]).*$", "\\1", suffix)
        replicate <- suppressWarnings(as.integer(sub("^[FM]", "", suffix)))

        return(data.table(sample = sample_name, genotype = "WBg", sex = sex, replicate = replicate,
            group_id = paste0("WBg_", sex, "_mean")))
    }

    return(NULL)
}

get_logTE_aliases <- function(sample_name, genotype, sex, replicate) {
    unique(c(sample_name, paste0("WBG_", sex, replicate), paste0("WBg_", sex, replicate)))
}

resolve_logTE_column <- function(sample_name, genotype, sex, replicate, candidates) {
    aliases <- get_logTE_aliases(sample_name, genotype, sex, replicate)
    resolve_unique_name(aliases, candidates)
}

resolve_wbg_edit_lib <- function(sex, replicate, candidates) {
    aliases <- c(paste0("WBG_", sex, replicate), paste0("WBg_", sex, replicate), paste0(sex,
        replicate))

    hit <- resolve_unique_name(aliases, candidates)
    if (!is.na(hit))
        return(hit)

    candidate_keys <- normalize_sample_name(candidates)
    suffix_key <- normalize_sample_name(paste0(sex, replicate))
    suffix_hits <- unique(candidates[endsWith(candidate_keys, suffix_key)])

    if (length(suffix_hits) == 1L)
        return(suffix_hits)

    NA_character_
}

# ============================================================
# 3. RNA secondary-structure helpers
# ============================================================

read_secondary_structure_tsv <- function(path) {
    if (!file.exists(path)) {
        stop("Missing RNALfold merged TSV: ", path)
    }

    raw <- fread(path, header = TRUE)

    required <- c("transcript", "start", "end")
    missing <- setdiff(required, names(raw))

    if (length(missing) > 0L) {
        stop("RNALfold merged TSV is missing columns: ", paste(missing, collapse = ", "))
    }

    bed <- raw[, .(FBtr = trimws(as.character(transcript)), structure_start_1based = suppressWarnings(as.integer(start)),
        structure_end_1based = suppressWarnings(as.integer(end)))]

    invalid <- bed[is.na(FBtr) | FBtr == "" | is.na(structure_start_1based) | is.na(structure_end_1based) |
        structure_start_1based < 1L | structure_end_1based < structure_start_1based]

    if (nrow(invalid) > 0L) {
        print(head(invalid, 20L))
        stop("RNALfold merged TSV contains ", nrow(invalid), " invalid rows.")
    }

    n_length_mismatch <- 0L

    if ("length" %in% names(raw)) {
        reported_length <- suppressWarnings(as.integer(raw$length))
        calculated_length <- suppressWarnings(as.integer(raw$end) - as.integer(raw$start) +
            1L)

        bad <- is.na(reported_length) | is.na(calculated_length) | reported_length != calculated_length

        n_length_mismatch <- sum(bad)

        if (n_length_mismatch > 0L) {
            print(head(raw[bad, .(transcript, start, end, length)], 20L))
            stop("RNALfold merged TSV has ", n_length_mismatch, " rows where length != end-start+1.")
        }
    }

    bed <- unique(bed)
    setkey(bed, FBtr, structure_start_1based, structure_end_1based)

    attr(bed, "n_source_rows") <- nrow(raw)
    attr(bed, "n_length_mismatch") <- n_length_mismatch

    bed
}

annotate_secondary_structure <- function(data, structure_bed) {
    dt <- as.data.table(copy(data))

    dt[, codonpos := suppressWarnings(as.integer(codonpos))]

    if (any(is.na(dt$codonpos) | dt$codonpos < 1L)) {
        stop("Invalid codonpos before secondary-structure annotation.")
    }

    dt[, structure_query_id := .I]
    dt[, `:=`(codon_cds_start_1based = (codonpos - 1L) * 3L + 1L, codon_cds_end_1based = codonpos *
        3L)]

    q <- dt[, .(structure_query_id, FBtr = as.character(FBtr), codon_cds_start_1based, codon_cds_end_1based)]

    setkey(q, FBtr, codon_cds_start_1based, codon_cds_end_1based)

    hits <- foverlaps(x = q, y = structure_bed, by.x = c("FBtr", "codon_cds_start_1based",
        "codon_cds_end_1based"), by.y = c("FBtr", "structure_start_1based", "structure_end_1based"),
        type = "any", nomatch = 0L)

    inside_ids <- unique(hits$structure_query_id)

    dt[, secondary_structure := factor(fifelse(structure_query_id %in% inside_ids, "inside",
        "outside"), levels = c("outside", "inside"))]

    dt[, structure_query_id := NULL]
    dt
}

# ============================================================
# 4. Model-output helpers
# ============================================================

get_max_gradient <- function(model) {
    if (is.null(model@optinfo$derivs) || is.null(model@optinfo$derivs$gradient)) {
        return(NA_real_)
    }
    max(abs(model@optinfo$derivs$gradient))
}

get_convergence_message <- function(model) {
    x <- model@optinfo$conv$lme4$messages
    if (is.null(x))
        return("OK")
    paste(x, collapse = "; ")
}

get_optimizer_code <- function(model) {
    x <- model@optinfo$conv$opt
    if (is.null(x))
        return(NA_character_)
    paste(x, collapse = "; ")
}

extract_fixed_effects <- function(model, group_id, model_type, genotype, sex, n_replicates,
    editing_definition) {
    coef_matrix <- coef(summary(model))

    result <- data.table(group_id = group_id, genotype = genotype, sex = sex, n_replicates = n_replicates,
        editing_definition = editing_definition, model_type = model_type, term = rownames(coef_matrix),
        estimate = coef_matrix[, "Estimate"], std_error = coef_matrix[, "Std. Error"])

    result[, df := if ("df" %in% colnames(coef_matrix)) {
        coef_matrix[, "df"]
    } else {
        NA_real_
    }]

    result[, statistic := if ("t value" %in% colnames(coef_matrix)) {
        coef_matrix[, "t value"]
    } else {
        estimate/std_error
    }]

    result[, p_value := if ("Pr(>|t|)" %in% colnames(coef_matrix)) {
        coef_matrix[, "Pr(>|t|)"]
    } else {
        NA_real_
    }]

    result[, critical_value := ifelse(is.finite(df), qt(0.975, df = df), qnorm(0.975))]

    result[, `:=`(conf_low = estimate - critical_value * std_error, conf_high = estimate +
        critical_value * std_error)]

    result[, `:=`(effect_percent = 100 * (exp(estimate) - 1), lower_percent = 100 * (exp(conf_low) -
        1), upper_percent = 100 * (exp(conf_high) - 1))]

    result[, critical_value := NULL]
    result
}

extract_type2_anova <- function(model, group_id, model_type, genotype, sex, n_replicates,
    editing_definition) {
    x <- as.data.frame(lmerTest::anova(model, type = 2, ddf = "Satterthwaite"))

    x$term <- rownames(x)
    rownames(x) <- NULL

    x <- as.data.table(x)
    x[, `:=`(group_id = group_id, genotype = genotype, sex = sex, n_replicates = n_replicates,
        editing_definition = editing_definition, model_type = model_type)]

    setcolorder(x, c("group_id", "genotype", "sex", "n_replicates", "editing_definition",
        "model_type", "term", setdiff(names(x), c("group_id", "genotype", "sex", "n_replicates",
            "editing_definition", "model_type", "term"))))

    x
}

# ============================================================
# 5. Original reproducible 1:20 matching logic
# ============================================================

match_edited_controls <- function(data, ratio = 20L, seed = 20260816L, near_pool_factor = 5L,
    max_codon_distance = Inf, without_replacement = TRUE, verbose = FALSE) {
    dt <- as.data.table(copy(data))

    required_cols <- c("FBtr", "codon", "AA", "codonpos", "upstream_category", "downstream_category",
        "edit_status")

    missing_cols <- setdiff(required_cols, names(dt))

    if (length(missing_cols) > 0L) {
        stop("Missing matching variables: ", paste(missing_cols, collapse = ", "))
    }

    dt[, row_id := .I]

    dt[, `:=`(FBtr_match = as.character(FBtr), codon_match = toupper(as.character(codon)),
        AA_match = as.character(AA), upstream_match = fifelse(is.na(upstream_category), "__MISSING__",
            as.character(upstream_category)), downstream_match = fifelse(is.na(downstream_category),
            "__MISSING__", as.character(downstream_category)), edit_status_match = as.character(edit_status))]

    edited <- dt[edit_status_match == "edited"]
    controls <- dt[edit_status_match == "unedited"]

    if (nrow(edited) == 0L)
        stop("No edited/pseudo-edited sites were found.")
    if (nrow(controls) == 0L)
        stop("No unedited controls were found.")

    setkey(controls, FBtr_match)
    used_control <- rep(FALSE, nrow(dt))

    set.seed(seed)
    edited_order <- sample(seq_len(nrow(edited)))

    pair_list <- vector("list", nrow(edited))
    pair_counter <- 0L

    for (ii in seq_along(edited_order)) {
        case <- edited[edited_order[ii]]

        case_row_id <- case$row_id[1]
        case_fbtr <- case$FBtr_match[1]
        case_codon <- case$codon_match[1]
        case_AA <- case$AA_match[1]
        case_pos <- case$codonpos[1]
        case_up <- case$upstream_match[1]
        case_down <- case$downstream_match[1]

        cand <- controls[.(case_fbtr), nomatch = 0L]

        if (nrow(cand) == 0L)
            next

        if (without_replacement) {
            cand <- cand[!used_control[row_id]]
        }

        if (nrow(cand) == 0L)
            next

        cand[, `:=`(exact_codon = codon_match == case_codon, same_AA = AA_match == case_AA,
            same_upstream = upstream_match == case_up, same_downstream = downstream_match ==
                case_down, codon_distance = abs(codonpos - case_pos))]

        cand <- cand[(exact_codon | same_AA) & codon_distance <= max_codon_distance]

        if (nrow(cand) == 0L)
            next

        cand[, match_tier := fcase(exact_codon & same_upstream & same_downstream, 1L,
            exact_codon & (same_upstream | same_downstream), 2L, exact_codon, 3L, same_AA &
                same_upstream & same_downstream, 4L, same_AA & (same_upstream | same_downstream),
            5L, same_AA, 6L, default = 99L)]

        cand <- cand[match_tier < 99L]
        if (nrow(cand) == 0L)
            next

        selected_list <- list()
        selected_counter <- 0L
        n_needed <- ratio

        for (tier_i in 1:6) {
            if (n_needed <= 0L)
                break

            pool <- cand[match_tier == tier_i][order(codon_distance)]

            if (nrow(pool) == 0L)
                next

            pool_size <- min(nrow(pool), max(n_needed, n_needed * near_pool_factor))

            near_pool <- pool[seq_len(pool_size)]
            n_take <- min(n_needed, nrow(near_pool))

            sampled_index <- sample.int(n = nrow(near_pool), size = n_take, replace = FALSE)

            selected_counter <- selected_counter + 1L
            selected_list[[selected_counter]] <- near_pool[sampled_index]
            n_needed <- n_needed - n_take
        }

        if (length(selected_list) == 0L)
            next

        selected <- rbindlist(selected_list)

        if (without_replacement) {
            used_control[selected$row_id] <- TRUE
        }

        pair_counter <- pair_counter + 1L

        selected[, `:=`(match_id = sprintf("match_%05d", pair_counter), case_row_id = case_row_id,
            control_row_id = row_id)]

        pair_list[[pair_counter]] <- selected[, .(match_id, case_row_id, control_row_id,
            match_tier, codon_distance, exact_codon, same_AA, same_upstream, same_downstream)]

        if (verbose && ii%%50L == 0L) {
            message("Processed ", ii, "/", nrow(edited), " edited sites")
        }
    }

    if (pair_counter == 0L) {
        stop("No matched controls were found.")
    }

    pairs <- rbindlist(pair_list[seq_len(pair_counter)])

    case_map <- unique(pairs[, .(match_id, row_id = case_row_id, matched_role = "edited")])

    control_map <- pairs[, .(match_id, row_id = control_row_id, matched_role = "unedited",
        match_tier, codon_distance, exact_codon, same_AA, same_upstream, same_downstream)]

    matched_map <- rbindlist(list(case_map, control_map), fill = TRUE)

    matched_data <- merge(matched_map, dt, by = "row_id", all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)

    matched_data[, matched_role := factor(matched_role, levels = c("unedited", "edited"))]

    setorder(matched_data, match_id, matched_role)

    list(matched_data = matched_data, pairs = pairs, unmatched_case_row_ids = setdiff(edited$row_id,
        unique(pairs$case_row_id)), seed = seed, ratio = ratio)
}

# ============================================================
# 6. Locate target phyloP samples
# ============================================================

phyloP_files <- list.files(phyloP_dir, pattern = "\\.codonTE_AEP\\.AA\\.phyloP\\.txt$", full.names = TRUE,
    recursive = TRUE)

if (length(phyloP_files) == 0L) {
    stop("No phyloP files were found in: ", phyloP_dir)
}

sample_names <- sub("\\.codonTE_AEP\\.AA\\.phyloP\\.txt$", "", basename(phyloP_files))

parsed_list <- lapply(sample_names, parse_target_sample)
keep <- !vapply(parsed_list, is.null, logical(1L))

if (!any(keep)) {
    stop("No target WBg_F/M phyloP files were found.")
}

sample_map <- rbindlist(parsed_list[keep])
sample_map[, phyloP_file := phyloP_files[keep]]

# Select the configured replicate numbers.
sample_map <- sample_map[(group_id == "WBg_F_mean" & replicate %in% 1:3) | (group_id == "WBg_M_mean" &
    replicate %in% 1:3)]

setorder(sample_map, group_id, replicate)

design_check <- merge(expected_design, sample_map[, .(observed_n_reps = uniqueN(replicate),
    replicate_ids = paste(sort(unique(replicate)), collapse = ","), samples = paste(sample[order(replicate)],
        collapse = ",")), by = group_id], by = "group_id", all.x = TRUE)

design_check[is.na(observed_n_reps), observed_n_reps := 0L]

if (any(design_check$observed_n_reps != design_check$expected_n_reps)) {
    print(design_check)
    stop("Target replicate design is incomplete. Expected WBg F=3 and WBg M=3.")
}

fwrite(design_check, file.path(output_dir, "group_replicate_design.tsv"), sep = "\t")

cat("[INFO] Target biological-replicate design:\n")
print(design_check)

# ============================================================
# 7. Read common files
# ============================================================

if (!file.exists(wbg_edit_file)) {
    stop("Missing WBG editing file: ", wbg_edit_file)
}

trans_length <- fread(trans_file, header = FALSE)

if (ncol(trans_length) < 4L) {
    stop("trans_file must contain at least four columns.")
}

trans_length <- trans_length[, 1:4]
setnames(trans_length, c("FBtr", "start", "end", "cds"))
trans_length[, codonlength := (end + 1 - start)/3]
trans_length <- unique(trans_length[, .(FBtr, codonlength)], by = "FBtr")

codon_edit <- fread(codon_edit_file, header = FALSE)

if (ncol(codon_edit) < 2L) {
    stop("codon_edit_file must contain at least two columns.")
}

codon_edit <- codon_edit[, 1:2]
setnames(codon_edit, c("FBtr_codonpos", "codon_edited"))
codon_edit <- unique(codon_edit, by = "FBtr_codonpos")

tr_gn <- fread(tr_gn_file)

if (!all(c("tx_name", "gene_id") %in% names(tr_gn))) {
    stop("tr_gn_file must contain tx_name and gene_id.")
}

required_tx_columns <- c("tx_name", "tx_len", "cds_len", "utr5_len", "utr3_len")
if (!all(required_tx_columns %in% names(tr_gn))) stop("Transcript annotation lacks full-length/UTR fields")
if (anyDuplicated(tr_gn$tx_name)) stop("Duplicate transcript IDs in transcript annotation")
for (column in required_tx_columns[-1L]) set(tr_gn, j = column, value = as.numeric(tr_gn[[column]]))
coding_tx <- tr_gn[cds_len > 0]
if (any(!is.finite(coding_tx$tx_len) | coding_tx$tx_len != coding_tx$cds_len + coding_tx$utr5_len +
    coding_tx$utr3_len)) {
    stop("Full transcript length does not equal CDS + 5-prime UTR + 3-prime UTR")
}

logTE_WBg <- fread(logTE_file_WBg)

if (!"gene_id" %in% names(logTE_WBg)) {
    stop("WBg logTE matrix does not contain gene_id.")
}

logTE_candidates_WBg <- setdiff(names(logTE_WBg), "gene_id")

rscu_dmel_df <- fread(rscu_file)

if (!all(c("codon", "category") %in% names(rscu_dmel_df))) {
    stop("RSCU file must contain codon and category.")
}

category_map <- setNames(as.character(rscu_dmel_df$category), toupper(as.character(rscu_dmel_df$codon)))

structure_bed <- read_secondary_structure_tsv(structure_bed_file)

fwrite(data.table(structure_bed_file = structure_bed_file, input_format = "RNALfold_merged_tsv_with_header",
    coordinate_system = "1_based_closed_inclusive", codon_rule = "inside_if_any_codon_nucleotide_overlaps_any_interval",
    n_source_rows = as.integer(attr(structure_bed, "n_source_rows")), n_structure_intervals = nrow(structure_bed),
    n_structure_transcripts = uniqueN(structure_bed$FBtr), n_length_mismatch = as.integer(attr(structure_bed,
        "n_length_mismatch"))), file.path(output_dir, "secondary_structure_bed_summary.tsv"),
    sep = "\t")

# ============================================================
# 8. Read WBG editing table and construct sex-specific mean editing
# ============================================================

edit_all <- fread(wbg_edit_file)

required_edit <- c("handle", "lib", "refcount", "altcount", "new_transcript", "new_cds_pos_1based")

missing_edit <- setdiff(required_edit, names(edit_all))

if (length(missing_edit) > 0L) {
    stop("WBG editing file is missing columns: ", paste(missing_edit, collapse = ", "))
}

if ("target_region" %in% names(edit_all)) {
    # Preserve the original sites_final behavior:
    # explicit non-CDS is removed; NA target_region is retained.
    edit_all <- edit_all[is.na(target_region) | target_region == "CDS"]
}

edit_all[, `:=`(site_key = trimws(as.character(handle)), lib = trimws(as.character(lib)),
    FBtr = trimws(as.character(new_transcript)), cds_edit_pos = suppressWarnings(as.integer(new_cds_pos_1based)),
    refcount = suppressWarnings(as.numeric(refcount)), altcount = suppressWarnings(as.numeric(altcount)))]

edit_all <- edit_all[!is.na(site_key) & site_key != "" & !is.na(lib) & lib != "" & !is.na(FBtr) &
    FBtr != "" & !is.na(cds_edit_pos) & cds_edit_pos >= 1L]

edit_all <- unique(edit_all, by = c("lib", "site_key", "FBtr", "cds_edit_pos"))

edit_all[, codonpos := (cds_edit_pos - 1L)%/%3L + 1L]

available_wbg_libs <- sort(unique(edit_all$lib))

wbg_edit_lib_map <- rbindlist(lapply(c("F", "M"), function(sex_i) {
    rbindlist(lapply(1:3, function(rep_i) {
        data.table(sex = sex_i, replicate = rep_i, edit_lib = resolve_wbg_edit_lib(sex_i,
            rep_i, available_wbg_libs))
    }))
}))

if (anyNA(wbg_edit_lib_map$edit_lib)) {
    print(wbg_edit_lib_map)
    cat("[INFO] Available WBG lib names:\n")
    print(available_wbg_libs)
    stop("Could not map all WBG F1-F3/M1-M3 editing libraries.")
}

fwrite(wbg_edit_lib_map, file.path(output_dir, "wbg_edit_library_mapping.tsv"), sep = "\t")

build_wbg_mean_editing <- function(sex_i) {
    libs <- wbg_edit_lib_map[sex == sex_i, .(replicate, edit_lib)]

    rep_tables <- lapply(seq_len(nrow(libs)), function(i) {
        lib_i <- libs$edit_lib[i]
        rep_i <- libs$replicate[i]

        x <- edit_all[lib == lib_i, .(edit_refcount = sum(refcount, na.rm = TRUE), edit_altcount = sum(altcount,
            na.rm = TRUE), n_edit_sites = .N), by = .(FBtr, codonpos)]

        x[, editing_level_rep := fifelse(edit_refcount + edit_altcount > 0, edit_altcount/(edit_refcount +
            edit_altcount), NA_real_)]

        x[, `:=`(replicate = rep_i, edit_lib = lib_i)]

        x
    })

    long <- rbindlist(rep_tables, use.names = TRUE, fill = TRUE)

    # Create a complete codon x replicate grid for the union of WBG edit
    # codons. Missing replicate records are assigned editing_level = 0.
    # This makes 'mean editing level' a mean across all three biological
    # replicates rather than a mean only across replicates where a row exists.
    site_union <- unique(long[, .(FBtr, codonpos)])

    grid <- CJ(site_row = seq_len(nrow(site_union)), replicate = sort(unique(libs$replicate)))

    grid <- merge(grid, site_union[, site_row := .I], by = "site_row", all.x = TRUE)

    grid <- merge(grid, long[, .(FBtr, codonpos, replicate, editing_level_rep, edit_refcount,
        edit_altcount, n_edit_sites)], by = c("FBtr", "codonpos", "replicate"), all.x = TRUE)

    grid[, editing_level_for_mean := fcoalesce(editing_level_rep, 0)]

    mean_map <- grid[, .(editing_level_mean = mean(editing_level_for_mean), n_editing_reps_with_record = sum(!is.na(editing_level_rep)),
        n_editing_reps_positive = sum(!is.na(editing_level_rep) & editing_level_rep > 0),
        edit_refcount_sum = sum(edit_refcount, na.rm = TRUE), edit_altcount_sum = sum(edit_altcount,
            na.rm = TRUE), n_edit_sites_rep_sum = sum(n_edit_sites, na.rm = TRUE)), by = .(FBtr,
        codonpos)]

    mean_map[, editing_reference_sex := sex_i]
    mean_map[, edit_status_reference := "edited"]

    list(long = long, mean = mean_map)
}

wbg_edit_F <- build_wbg_mean_editing("F")
wbg_edit_M <- build_wbg_mean_editing("M")

fwrite(wbg_edit_F$mean, file.path(output_dir, "WBg_F_mean_editing_by_codon.tsv.gz"), sep = "\t")

fwrite(wbg_edit_M$mean, file.path(output_dir, "WBg_M_mean_editing_by_codon.tsv.gz"), sep = "\t")

# ============================================================
# 9. Resolve replicate-specific logTE columns
# ============================================================

sample_map[, logTE_source := "WBg"]

sample_map[, logTE_column := vapply(seq_len(.N), function(i) resolve_logTE_column(sample_name = sample[i],
    genotype = genotype[i], sex = sex[i], replicate = replicate[i], candidates = logTE_candidates_WBg),
    character(1))]

if (anyNA(sample_map$logTE_column)) {
    print(sample_map[is.na(logTE_column)])
    cat("[INFO] Available WBg logTE columns:\n")
    print(logTE_candidates_WBg)
    stop("Could not map one or more target WBg samples to logTE columns.")
}

fwrite(sample_map, file.path(output_dir, "group_sample_mapping.tsv"), sep = "\t")

# ============================================================
# 10. Read one raw replicate WITHOUT the old replicate-level
#     mRNA>=10 / Asite_median>0 filtering.
# ============================================================

m6a_filter_qc_list <- list()

read_raw_replicate <- function(phyloP_file, sample_name, group_id, genotype, sex, replicate,
    logTE_source, logTE_column) {
    cat("[READ] ", sample_name, " group=", group_id, " logTE=", logTE_column, "\n", sep = "")

    tmp <- fread(phyloP_file)

    if (!"FBtr_codonpos" %in% names(tmp) && all(c("FBtr", "codonpos") %in% names(tmp))) {
        tmp[, FBtr_codonpos := paste(FBtr, codonpos, sep = "_")]
    }

    required <- c("FBtr", "FBtr_codonpos", "codonpos", "codon", "AA", "gmRC", "mRNA", "Asite_median",
        "codon_phyloP_mean", "codon_phyloP_covered", "chr", "pos")

    missing <- setdiff(required, names(tmp))

    if (length(missing) > 0L) {
        stop(sample_name, " phyloP file is missing columns: ", paste(missing, collapse = ", "))
    }

    if ("codonlength" %in% names(tmp)) {
        tmp[, codonlength := NULL]
    }

    tmp[, `:=`(FBtr = as.character(FBtr), FBtr_codonpos = as.character(FBtr_codonpos), codonpos = suppressWarnings(as.integer(codonpos)),
        codon = toupper(as.character(codon)), AA = as.character(AA), gmRC = suppressWarnings(as.numeric(gmRC)),
        mRNA = suppressWarnings(as.numeric(mRNA)), Asite_median = suppressWarnings(as.numeric(Asite_median)),
        codon_phyloP_mean = suppressWarnings(as.numeric(codon_phyloP_mean)), codon_phyloP_covered = suppressWarnings(as.integer(codon_phyloP_covered)),
        chr = trimws(as.character(chr)), pos = suppressWarnings(as.integer(pos)))]

    tmp <- merge(tmp, trans_length, by = "FBtr", all = FALSE)

    # Apply the transcript-length and codon-position criteria from the
    # script. Only the mRNA and Asite filters are moved to the group level.
    tmp <- tmp[codonlength > 150 & codonpos > 15 & (codonlength - codonpos) > 5]

    # Preserve the original codon universe restriction.
    df <- merge(codon_edit, tmp, by = "FBtr_codonpos", all = FALSE)

    setDT(df)

    df <- annotate_secondary_structure(data = df, structure_bed = structure_bed)

    df[, gene_id := tr_gn$gene_id[match(FBtr, tr_gn$tx_name)]]

    logTE_table <- logTE_WBg

    df[, logTE := suppressWarnings(as.numeric(logTE_table[[logTE_column]][match(gene_id,
        logTE_table$gene_id)]))]

    df[, codon_edited := toupper(as.character(codon_edited))]

    df[, `:=`(sample = sample_name, group_id = group_id, genotype = genotype, sex = sex,
        replicate = as.integer(replicate))]

    # Use the same exact genomic-handle exclusion as the other updated gmRC
    # panels. Filtering here precedes replicate collapse, coverage filtering,
    # covariate construction, model fitting, and matched-control sampling.
    df[, genomic_handle := fifelse(!is.na(chr) & chr != "" & !is.na(pos), paste(chr, format(pos,
        scientific = FALSE, trim = TRUE), sep = "_"), NA_character_)]
    is_m6a <- !is.na(df$genomic_handle) & df$genomic_handle %chin% m6a_handles
    m6a_filter_qc_list[[sample_name]] <<- data.table(sample = sample_name, group_id = group_id,
        genotype = genotype, sex = sex, replicate = as.integer(replicate), n_before_m6a = nrow(df),
        n_m6a_removed = sum(is_m6a), n_after_m6a = sum(!is_m6a))
    df <- df[!is_m6a]

    # There should be one row per sample/codon after codon_edit was reduced
    # to one row per FBtr_codonpos.
    dup <- df[, .N, by = FBtr_codonpos][N > 1L]
    if (nrow(dup) > 0L) {
        print(head(dup, 20L))
        stop("Duplicate FBtr_codonpos rows remain in sample ", sample_name)
    }

    df[, .(sample, group_id, genotype, sex, replicate, FBtr, FBtr_codonpos, codonpos, codon,
        codon_edited, AA, gmRC, mRNA, Asite_median, codon_phyloP_mean, codon_phyloP_covered,
        codonlength, gene_id, logTE, secondary_structure)]
}

raw_rep_list <- lapply(seq_len(nrow(sample_map)), function(i) {
    read_raw_replicate(phyloP_file = sample_map$phyloP_file[i], sample_name = sample_map$sample[i],
        group_id = sample_map$group_id[i], genotype = sample_map$genotype[i], sex = sample_map$sex[i],
        replicate = sample_map$replicate[i], logTE_source = sample_map$logTE_source[i], logTE_column = sample_map$logTE_column[i])
})

fwrite(rbindlist(m6a_filter_qc_list, use.names = TRUE, fill = TRUE), file.path(output_dir,
    "m6A_filter_QC.tsv"), sep = "\t")

raw_all <- rbindlist(raw_rep_list, use.names = TRUE, fill = TRUE)

rm(raw_rep_list)
gc()

# ============================================================
# 11. Collapse biological replicates to one codon row per group
# ============================================================

# Check static annotations and collapse replicate observations using equivalent
# vectorized aggregation.
source(Sys.getenv("GLM_FAST_AGGREGATION_HELPER", unset = "/path/to/RNA_editing_DATAprocessing_new/code_upload/new_mRC_gmRC/figure_update/code/collapse_replicates_fast.R"))
message("Starting vectorized static-annotation QC and replicate collapse: ", nrow(raw_all),
    " rows")
collapse_result <- collapse_replicates_fast(raw_all, expected_design)
static_qc <- collapse_result$static_qc
static_bad <- collapse_result$static_bad
if (nrow(static_bad) > 0L) {
    fwrite(static_bad, file.path(output_dir, "ERROR_static_codon_annotation_conflicts.tsv"),
        sep = "\t")
    stop("Static codon annotations disagree across replicates. See ERROR_static_codon_annotation_conflicts.tsv.")
}
group_mean <- collapse_result$group_mean
rm(collapse_result)
message("Replicate collapse complete: ", nrow(group_mean), " group/codon observations")

group_mean <- merge(group_mean, expected_design, by = c("group_id", "expected_n_reps"), all.x = TRUE,
    sort = FALSE)

# QC before coverage filtering.
group_filter_qc <- group_mean[, .(n_codons_before_group_coverage = .N, n_codons_all_expected_rep_rows = sum(n_reps_present ==
    expected_n_reps), n_pass_mRNA_sum = sum(mRNA_sum >= MIN_MRNA_PER_REP * expected_n_reps),
    n_pass_Asite_any = sum(Asite_any_positive), n_pass_both_group_coverage = sum(mRNA_sum >=
        MIN_MRNA_PER_REP * expected_n_reps & Asite_any_positive)), by = .(group_id, genotype,
    sex, expected_n_reps, editing_source)]

# Apply group-level coverage filters:
#   old: every replicate independently mRNA>=10 and Asite_median>0
#   new: sum(mRNA)>=10*n and ANY replicate Asite_median>0
group_mean <- group_mean[mRNA_sum >= MIN_MRNA_PER_REP * expected_n_reps & Asite_any_positive]

# ============================================================
# 12. Construct codon covariates after group filtering
# ============================================================

group_mean[, category_before_chr := unname(category_map[toupper(codon)])]

setorder(group_mean, group_id, FBtr, codonpos)

group_mean[, `:=`(upstream_pos = shift(codonpos), downstream_pos = shift(codonpos, type = "lead"),
    upstream_category_chr = shift(category_before_chr), downstream_category_chr = shift(category_before_chr,
        type = "lead")), by = .(group_id, FBtr)]

group_mean[is.na(upstream_pos) | codonpos - upstream_pos != 1, upstream_category_chr :=
    NA_character_]

group_mean[is.na(downstream_pos) | downstream_pos - codonpos != 1, downstream_category_chr :=
    NA_character_]

group_mean[, c("upstream_pos", "downstream_pos") := NULL]

# ============================================================
# 13. Attach WBG sex-specific mean editing levels
# ============================================================

prepare_edit_map_for_join <- function(x) {
    copy(x)[, .(FBtr, codonpos, editing_level_mean, n_editing_reps_with_record, n_editing_reps_positive,
        edit_refcount_sum, edit_altcount_sum, n_edit_sites_rep_sum)]
}

edit_F_map <- prepare_edit_map_for_join(wbg_edit_F$mean)
edit_M_map <- prepare_edit_map_for_join(wbg_edit_M$mean)

group_list <- lapply(unique(group_mean$group_id), function(group_i) {
    x <- copy(group_mean[group_id == group_i])

    if (group_i == "WBg_F_mean") {
        edit_map <- edit_F_map
        definition <- "real_WBg_F_mean_editing"
    } else if (group_i == "WBg_M_mean") {
        edit_map <- edit_M_map
        definition <- "real_WBg_M_mean_editing"
    } else {
        stop("Unexpected group: ", group_i)
    }

    x[edit_map, on = .(FBtr, codonpos), `:=`(editing_level_mean = i.editing_level_mean, n_editing_reps_with_record = i.n_editing_reps_with_record,
        n_editing_reps_positive = i.n_editing_reps_positive, edit_refcount_sum = i.edit_refcount_sum,
        edit_altcount_sum = i.edit_altcount_sum, n_edit_sites_rep_sum = i.n_edit_sites_rep_sum)]

    # Preserve the original sites_final logic conceptually:
    # presence in the corresponding sex-specific WBg map = edited.
    x[, edit_status := factor(fifelse(is.na(editing_level_mean), "unedited", "edited"),
        levels = c("unedited", "edited"))]

    x[, editing_level_mean := fcoalesce(editing_level_mean, 0)]
    x[, editing_level_10pct := editing_level_mean/0.1]
    x[, editing_definition := definition]

    x
})

group_mean <- rbindlist(group_list, use.names = TRUE, fill = TRUE)

rm(group_list)
gc()

# ============================================================
# 14. Apply final codon and model-variable filters
# ============================================================

prepare_group_model_data <- function(x, group_id_i) {
    dt <- copy(x[group_id == group_id_i])

    if (nrow(dt) == 0L) {
        stop("No group-level codons available for ", group_id_i)
    }

    # Original script retained only reference codons containing A.
    dt <- dt[!is.na(codon) & grepl("A", toupper(codon), fixed = TRUE)]

    n_after_A <- nrow(dt)

    # Original phyloP completeness condition.
    dt <- dt[codon_phyloP_covered == 3]
    n_after_phyloP <- nrow(dt)

    dt[, secondary_structure := factor(secondary_structure_chr, levels = c("outside",
        "inside"))]

    dt[, category_before := factor(category_before_chr, levels = c("nonoptimal", "optimal"))]

    dt[, upstream_category := factor(upstream_category_chr, levels = c("nonoptimal", "optimal"))]

    dt[, downstream_category := factor(downstream_category_chr, levels = c("nonoptimal",
        "optimal"))]

    dt[, AA := factor(as.character(AA))]

    dt[, log_gmRC := log1p(gmRC)]
    # Full mature transcript length includes both UTRs; codonlength remains
    # exclusively the CDS-based positional/eligibility variable.
    dt[, transcript_length_nt := tr_gn$tx_len[match(FBtr, tr_gn$tx_name)]]
    if (any(!is.finite(dt$transcript_length_nt) | dt$transcript_length_nt <= 0)) {
        stop("Missing or invalid full transcript length in model group: ", group_id_i)
    }
    dt[, log_transcript_length := log(transcript_length_nt)]

    dt[, phyloP_z := safe_zscore(codon_phyloP_mean)]
    dt[, logTE_z := safe_zscore(logTE)]
    dt[, transcript_length_z := safe_zscore(log_transcript_length)]

    joint_vars <- c("log_gmRC", "edit_status", "editing_level_10pct", "category_before",
        "AA", "phyloP_z", "logTE_z", "transcript_length_z", "upstream_category", "downstream_category",
        "secondary_structure", "FBtr", "codon", "codonpos")

    dt <- dt[complete.cases(dt[, ..joint_vars])]

    n_after_complete <- nrow(dt)

    dt <- dt[is.finite(log_gmRC) & is.finite(phyloP_z) & is.finite(logTE_z) & is.finite(transcript_length_z)]

    n_after_finite <- nrow(dt)

    dt <- as.data.table(droplevels(as.data.frame(dt)))

    dt[, edit_status := factor(edit_status, levels = c("unedited", "edited"))]

    dt[, category_before := factor(category_before, levels = c("nonoptimal", "optimal"))]

    dt[, upstream_category := factor(upstream_category, levels = c("nonoptimal", "optimal"))]

    dt[, downstream_category := factor(downstream_category, levels = c("nonoptimal", "optimal"))]

    dt[, secondary_structure := factor(secondary_structure, levels = c("outside", "inside"))]

    dt[, AA := droplevels(factor(AA))]

    if (!"Ala" %in% levels(dt$AA)) {
        stop(group_id_i, ": Ala is absent after filtering, so a common AA reference ", "cannot be used.")
    }

    dt[, AA := relevel(AA, ref = "Ala")]

    n_edited <- sum(dt$edit_status == "edited")
    n_unedited <- sum(dt$edit_status == "unedited")

    if (n_edited < 2L) {
        stop(group_id_i, ": fewer than two edited/pseudo-edited codons.")
    }

    if (n_unedited < 2L) {
        stop(group_id_i, ": fewer than two unedited codons.")
    }

    mean_edit <- dt[edit_status == "edited", mean(editing_level_10pct, na.rm = TRUE)]

    sd_edit <- dt[edit_status == "edited", sd(editing_level_10pct, na.rm = TRUE)]

    if (!is.finite(sd_edit) || sd_edit == 0) {
        stop(group_id_i, ": mean editing level has no variation among edited/pseudo-edited ",
            "codons.")
    }

    dt[, editing_level_within_editSD := fifelse(edit_status == "edited", (editing_level_10pct -
        mean_edit)/sd_edit, 0)]

    factor_variables <- c("edit_status", "category_before", "AA", "upstream_category", "downstream_category",
        "secondary_structure")

    factor_level_counts <- vapply(dt[, ..factor_variables], function(z) length(unique(z)),
        integer(1))

    invalid_factors <- names(factor_level_counts[factor_level_counts < 2L])

    if (length(invalid_factors) > 0L) {
        stop(group_id_i, ": these factors have fewer than two levels: ", paste(invalid_factors,
            collapse = ", "))
    }

    qc <- data.table(group_id = group_id_i, n_after_group_coverage = nrow(x[group_id == group_id_i]),
        n_after_A_codon = n_after_A, n_after_phyloP_covered3 = n_after_phyloP, n_after_complete_cases = n_after_complete,
        n_after_finite = n_after_finite, n_final = nrow(dt), n_FBtr = uniqueN(dt$FBtr), n_edited = n_edited,
        n_unedited = n_unedited, n_structure_inside = sum(dt$secondary_structure == "inside"),
        n_structure_outside = sum(dt$secondary_structure == "outside"), mean_editing_level_10pct = mean_edit,
        sd_editing_level_10pct = sd_edit)

    list(data = dt, mean_edit = mean_edit, sd_edit = sd_edit, qc = qc)
}

prepared_groups <- lapply(expected_design$group_id, function(g) prepare_group_model_data(group_mean,
    g))
names(prepared_groups) <- expected_design$group_id

# Save the actual model-input codon tables for direct inspection.
for (g in names(prepared_groups)) {
    gdir <- file.path(output_dir, "per_group", g)
    dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

    fwrite(prepared_groups[[g]]$data, file.path(gdir, "group_mean_model_input.tsv.gz"), sep = "\t")
}

final_filter_qc <- rbindlist(lapply(prepared_groups, function(z) z$qc), use.names = TRUE,
    fill = TRUE)

group_filter_qc <- merge(group_filter_qc, final_filter_qc, by = "group_id", all = TRUE)

fwrite(group_filter_qc, file.path(output_dir, "group_mean_codon_filter_qc.tsv"), sep = "\t")

# ============================================================
# 15. Fit Full and Matched models
# ============================================================

full_formula <- log_gmRC ~ edit_status + editing_level_within_editSD + category_before +
    AA + phyloP_z + logTE_z + transcript_length_z + upstream_category + downstream_category +
    secondary_structure + (1 | FBtr)

matched_formula <- log_gmRC ~ edit_status_matched + category_before + AA + phyloP_z + logTE_z +
    transcript_length_z + upstream_category + downstream_category + secondary_structure +
    (1 | FBtr) + (1 | match_id)

all_fixed_list <- list()
all_anova_list <- list()
all_model_qc_list <- list()
failed_list <- list()

for (group_i in names(prepared_groups)) {
    meta <- expected_design[group_id == group_i]
    dt <- prepared_groups[[group_i]]$data

    gdir <- file.path(output_dir, "per_group", group_i)

    cat("\n==================================================\n", "[MODEL] ", group_i, " n=",
        nrow(dt), " edited/pseudo-edited=", sum(dt$edit_status == "edited"), " unedited=",
        sum(dt$edit_status == "unedited"), "\n", sep = "")

    editing_definition_i <- unique(dt$editing_definition)
    if (length(editing_definition_i) != 1L) {
        stop("Unexpected multiple editing definitions in ", group_i)
    }

    full_model <- tryCatch(lmerTest::lmer(full_formula, data = dt, REML = FALSE, na.action = na.omit,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))),
        error = function(e) {
            failed_list[[length(failed_list) + 1L]] <<- data.table(group_id = group_i, stage = "full_model",
                error = conditionMessage(e))
            NULL
        })

    if (is.null(full_model)) {
        next
    }

    full_fixed <- extract_fixed_effects(full_model, group_id = group_i, model_type = "GroupMean_Full",
        genotype = meta$genotype, sex = meta$sex, n_replicates = meta$expected_n_reps, editing_definition = editing_definition_i)

    full_anova <- tryCatch(extract_type2_anova(full_model, group_id = group_i, model_type = "GroupMean_Full",
        genotype = meta$genotype, sex = meta$sex, n_replicates = meta$expected_n_reps, editing_definition = editing_definition_i),
        error = function(e) NULL)

    all_fixed_list[[length(all_fixed_list) + 1L]] <- full_fixed

    if (!is.null(full_anova)) {
        all_anova_list[[length(all_anova_list) + 1L]] <- full_anova
    }

    fwrite(full_fixed, file.path(gdir, "fixed_effects_full.tsv"), sep = "\t")

    if (!is.null(full_anova)) {
        fwrite(full_anova, file.path(gdir, "type2_anova_full.tsv"), sep = "\t")
    }

    capture.output(summary(full_model), file = file.path(gdir, "summary_full.txt"))

    if (SAVE_MODELS) {
        saveRDS(full_model, file.path(gdir, "model_full.rds"))
    }

    fwrite(dt[, .N, by = .(edit_status, secondary_structure)], file.path(gdir, "secondary_structure_counts_full.tsv"),
        sep = "\t")

    model_qc <- data.table(group_id = group_i, genotype = meta$genotype, sex = meta$sex,
        n_replicates = meta$expected_n_reps, editing_definition = editing_definition_i, n_full = nrow(dt),
        n_FBtr_full = uniqueN(dt$FBtr), n_edited_full = sum(dt$edit_status == "edited"),
        n_unedited_full = sum(dt$edit_status == "unedited"), full_singular = lme4::isSingular(full_model,
            tol = 1e-05), full_max_gradient = get_max_gradient(full_model), full_optimizer_code = get_optimizer_code(full_model),
        full_convergence = get_convergence_message(full_model), matched_status = ifelse(RUN_MATCHED,
            "not_started", "not_requested"), n_matched = NA_integer_, n_match_groups = NA_integer_,
        n_FBtr_matched = NA_integer_, matched_singular = NA, matched_max_gradient = NA_real_,
        matched_optimizer_code = NA_character_, matched_convergence = NA_character_)

    if (RUN_MATCHED) {
        seed_i <- BASE_SEED + sum(utf8ToInt(group_i) * seq_along(utf8ToInt(group_i)))%%100000L

        matched_result <- tryCatch(match_edited_controls(data = dt, ratio = MATCH_RATIO,
            seed = seed_i, near_pool_factor = 5L, max_codon_distance = Inf, without_replacement = TRUE,
            verbose = FALSE), error = function(e) {
            failed_list[[length(failed_list) + 1L]] <<- data.table(group_id = group_i, stage = "matching",
                error = conditionMessage(e))
            NULL
        })

        if (!is.null(matched_result)) {
            matched_dat <- matched_result$matched_data
            matched_pairs <- matched_result$pairs

            matched_dat[, edit_status_matched := factor(matched_role, levels = c("unedited",
                "edited"))]

            matched_model <- tryCatch(lmerTest::lmer(matched_formula, data = matched_dat,
                REML = FALSE, na.action = na.omit, control = lme4::lmerControl(optimizer = "bobyqa",
                  optCtrl = list(maxfun = 2e+05))), error = function(e) {
                failed_list[[length(failed_list) + 1L]] <<- data.table(group_id = group_i,
                  stage = "matched_model", error = conditionMessage(e))
                NULL
            })

            if (!is.null(matched_model)) {
                matched_fixed <- extract_fixed_effects(matched_model, group_id = group_i,
                  model_type = "GroupMean_Matched", genotype = meta$genotype, sex = meta$sex,
                  n_replicates = meta$expected_n_reps, editing_definition = editing_definition_i)

                matched_anova <- tryCatch(extract_type2_anova(matched_model, group_id = group_i,
                  model_type = "GroupMean_Matched", genotype = meta$genotype, sex = meta$sex,
                  n_replicates = meta$expected_n_reps, editing_definition = editing_definition_i),
                  error = function(e) NULL)

                all_fixed_list[[length(all_fixed_list) + 1L]] <- matched_fixed

                if (!is.null(matched_anova)) {
                  all_anova_list[[length(all_anova_list) + 1L]] <- matched_anova
                }

                fwrite(matched_fixed, file.path(gdir, "fixed_effects_matched.tsv"), sep = "\t")

                if (!is.null(matched_anova)) {
                  fwrite(matched_anova, file.path(gdir, "type2_anova_matched.tsv"), sep = "\t")
                }

                fwrite(matched_pairs, file.path(gdir, "matched_pairs.tsv.gz"), sep = "\t")

                fwrite(data.table(unmatched_case_row_id = matched_result$unmatched_case_row_ids),
                  file.path(gdir, "unmatched_edited_case_row_ids.tsv"), sep = "\t")

                capture.output(summary(matched_model), file = file.path(gdir, "summary_matched.txt"))

                if (SAVE_MODELS) {
                  saveRDS(matched_model, file.path(gdir, "model_matched.rds"))
                }

                model_qc[, `:=`(matched_status = "success", n_matched = nrow(matched_dat),
                  n_match_groups = uniqueN(matched_dat$match_id), n_FBtr_matched = uniqueN(matched_dat$FBtr),
                  matched_singular = lme4::isSingular(matched_model, tol = 1e-05), matched_max_gradient = get_max_gradient(matched_model),
                  matched_optimizer_code = get_optimizer_code(matched_model), matched_convergence = get_convergence_message(matched_model))]
            } else {
                model_qc[, matched_status := "model_failed"]
            }
        } else {
            model_qc[, matched_status := "matching_failed"]
        }
    }

    all_model_qc_list[[length(all_model_qc_list) + 1L]] <- model_qc
}

if (length(all_fixed_list) == 0L) {
    stop("No group-mean model was successfully fitted.")
}

all_fixed <- rbindlist(all_fixed_list, use.names = TRUE, fill = TRUE)

all_anova <- if (length(all_anova_list) > 0L) {
    rbindlist(all_anova_list, use.names = TRUE, fill = TRUE)
} else {
    data.table()
}

all_model_qc <- if (length(all_model_qc_list) > 0L) {
    rbindlist(all_model_qc_list, use.names = TRUE, fill = TRUE)
} else {
    data.table()
}

fwrite(all_fixed, file.path(output_dir, "all_group_mean_fixed_effects.tsv"), sep = "\t")

if (nrow(all_anova) > 0L) {
    fwrite(all_anova, file.path(output_dir, "all_group_mean_type2_anova.tsv"), sep = "\t")
}

if (nrow(all_model_qc) > 0L) {
    fwrite(all_model_qc, file.path(output_dir, "all_group_mean_model_qc.tsv"), sep = "\t")
}

if (length(failed_list) > 0L) {
    fwrite(rbindlist(failed_list, use.names = TRUE, fill = TRUE), file.path(output_dir, "failed_groups.tsv"),
        sep = "\t")
}

# ============================================================
# 16. Export common and AA effects together
# ============================================================

term_label_map <- c(edit_statusedited = "Edited vs unedited", edit_status_matchededited = "Edited vs unedited",
    editing_level_within_editSD = "Editing level (+1 SD; edited)", category_beforeoptimal = "Current codon optimal",
    phyloP_z = "Codon phyloP (+1 SD)", logTE_z = "Translation efficiency (+1 SD)", transcript_length_z = "Transcript length (+1 SD)",
    upstream_categoryoptimal = "Upstream codon optimal", downstream_categoryoptimal = "Downstream codon optimal",
    secondary_structureinside = "Inside RNA secondary structure")

common_effects <- copy(all_fixed[term %in% names(term_label_map)])

common_effects[, variable := unname(term_label_map[term])]

fwrite(common_effects, file.path(output_dir, "all_group_mean_common_effects.tsv"), sep = "\t")

aa_effects <- copy(all_fixed[grepl("^AA", term)])
aa_effects[, amino_acid := sub("^AA", "", term)]

fwrite(aa_effects, file.path(output_dir, "all_group_mean_amino_acid_effects.tsv"), sep = "\t")

# ============================================================
# 17. Optional plots
# ============================================================

if (MAKE_PLOTS && requireNamespace("ggplot2", quietly = TRUE)) {
    plot_dat <- copy(common_effects)

    plot_dat[, group_id := factor(group_id, levels = c("WBg_F_mean", "WBg_M_mean"))]

    plot_dat[, model_type := factor(model_type, levels = c("GroupMean_Full", "GroupMean_Matched"))]

    p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = variable, y = effect_percent, ymin = lower_percent,
        ymax = upper_percent, shape = model_type)) + ggplot2::geom_hline(yintercept = 0,
        linetype = "dashed", linewidth = 0.45) + ggplot2::geom_errorbar(width = 0.12, position = ggplot2::position_dodge(width = 0.5)) +
        ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.5)) +
        ggplot2::facet_wrap(~group_id, scales = "free_y") + ggplot2::labs(x = NULL, y = "Estimated change in 1 + gmRC (%)",
        shape = "Model", title = "Group-mean gmRC models", subtitle = "WBg F/M: three biological replicates per sex") +
        ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "top",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), strip.text = ggplot2::element_text(face = "bold"))

    ggplot2::ggsave(file.path(output_dir, "all_group_mean_common_effects.pdf"), p, width = 14,
        height = 8.5, units = "in")
}

# ============================================================
# 18. Metadata and completion
# ============================================================

metadata <- data.table(parameter = c("analysis", "mRNA_filter", "Asite_filter", "m6A_filter",
    "m6A_reference", "gmRC_aggregation", "editing_level_aggregation", "WBg_edit_status",
    "length_covariate", "length_annotation", "CDS_length_filter", "five_prime_filter", "three_prime_filter",
    "codon_universe", "A_codon_filter", "phyloP_filter", "matched_ratio", "WBg_logTE_file",
    "structure_file", "timestamp"), value = c("genotype_sex_group_mean_before_model", "sum(mRNA across group replicates) >= 10 * n_replicates",
    "at least one replicate has Asite_median > 0", "exclude exact Dmel genomic m6A handles before replicate collapse and modeling",
    m6a_file, "arithmetic mean of replicate gmRC values", paste0("replicate-level alt/(alt+ref), then arithmetic mean across all ",
        "3 WBg sex-specific replicates; missing editing records counted as 0"), "sex-specific WBg editing map",
    "z-score of ln(tx_len in nt), including 5-prime and 3-prime UTRs", tr_gn_file, "codonlength > 150",
    "codonpos > 15", "(codonlength - codonpos) > 5", basename(codon_edit_file), "reference codon contains A",
    "codon_phyloP_covered == 3", as.character(MATCH_RATIO), logTE_file_WBg, structure_bed_file,
    as.character(Sys.time())))

fwrite(metadata, file.path(output_dir, "run_metadata.tsv"), sep = "\t")

invisible(NULL)

cat("\nWT-only group-mean analysis completed.\n", "Groups: WBg_F_mean, WBg_M_mean\n", "Coverage: mRNA_sum >= 10*n; any replicate Asite_median > 0\n",
    "m6A: exact Dmel genomic handles excluded before replicate collapse\n", "Output directory:\n",
    normalizePath(output_dir), "\n", sep = "")
