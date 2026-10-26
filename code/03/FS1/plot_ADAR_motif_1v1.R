#!/usr/bin/env Rscript

# Figure S1: ADAR sequence motifs at one-to-one editing sites.
#
# This is an adaptation of plot_ADAR_motif_two_versions.R.  It uses the
# cleaned project bundle rather than the former sites_candidate directory:
#
#   ISO1/WBG : 01_rebuild/.../candidate.final and WBG annotations
#   sim/sec/ere: 02_species_1v1/.../*.candidate.1v1.kept.final.tsv
#   quantitative table: 05_plot/data/candidate.editing_level.observed.long.tsv
#
# The primary output, current_1v1_kept, uses all retained one-to-one sites.
# current_1v1_min2 is written as a sensitivity/comparability analysis using
# the current *.1v1.min2.handles.txt lists.  Both use the species' own genome
# coordinates and species-specific best A>G annotations.
#
# Motif definition is unchanged from the old script:
#   - pool alt counts and depth by site within a group (or sex);
#   - high editing: pooled depth >= MIN_DEPTH and editing level >= MIN_EDITING;
#   - extract a three-base genomic motif by default (FLANK=1);
#   - reverse-complement sites whose genomic center is T so the edited center
#     is represented as A;
#   - draw all, coding, and noncoding high-editing motifs for ALL/F/M.
#
# Sequence logos are drawn with ggseqlogo and arranged with gridExtra.
#
# Usage:
#   Rscript plot_ADAR_motif_1v1.R [OUTPUT_DIR] [MIN_DEPTH] [MIN_EDITING]
#       [PLOT_MODE] [FLANK]
#
# Defaults are the project output directory, 50, 0.5, both, and 1.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggseqlogo)
    library(grid)
    library(gridExtra)
})

font_family <- "sans"


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
    getwd()
}

required_env <- function(name) {
    value <- Sys.getenv(name)
    if (!nzchar(value))
        stop("Set input directory: ", name)
    normalizePath(value, mustWork = TRUE)
}
project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = TRUE)
candidate_root <- required_env("FS1_CANDIDATE_DIR")
species_root <- required_env("FS1_SPECIES_DIR")
plot_data_dir <- required_env("FS1_PLOT_DATA_DIR")
reference_dir <- required_env("REF_DIR")

output_dir <- if (length(args) >= 1L && nzchar(args[[1L]])) {
    args[[1L]]
} else {
    file.path(project_root, "results/03/FS1")
}
min_depth <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 50
min_editing <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.5
plot_mode <- if (length(args) >= 4L) tolower(args[[4L]]) else "both"
flank <- if (length(args) >= 5L) as.integer(args[[5L]]) else 1L

if (!is.finite(min_depth) || min_depth < 1) {
    stop("MIN_DEPTH must be a positive number.")
}
if (!is.finite(min_editing) || min_editing < 0 || min_editing > 1) {
    stop("MIN_EDITING must be between 0 and 1.")
}
if (!plot_mode %chin% c("bits", "prob", "both")) {
    stop("PLOT_MODE must be bits, prob, or both.")
}
if (is.na(flank) || flank < 1L) {
    stop("FLANK must be a positive integer.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
reference_dir <- normalizePath(reference_dir, mustWork = TRUE)

group_order <- c("ISO1", "WBG", "sim", "sec", "ere")
sex_order <- c("ALL", "F", "M")
species_order <- c("sim", "sec", "ere")

group_display <- c(ISO1 = "iso-1", WBG = "WBG", sim = "D. simulans", sec = "D. sechellia",
    ere = "D. erecta")

group_genome <- c(ISO1 = "dmel", WBG = "dmel", sim = "dsim", sec = "dsec", ere = "dere")

genome_fasta <- c(dmel = file.path(reference_dir, "dmel-all-chromosome-r6.12.fasta"), dsim = file.path(reference_dir,
    "dsim-all-chromosome-r2.02.fasta"), dsec = file.path(reference_dir, "dsec-all-chromosome-r1.3.fasta"),
    dere = file.path(reference_dir, "dere-all-chromosome-r1.05.fasta"))

annotation_files <- c(ISO1 = file.path(candidate_root, "candidate.master.ISO1_annotation.tsv"),
    WBG = file.path(candidate_root, "candidate_lists_by_group/WBG.candidate.detected.annotated.tsv"),
    sim = file.path(species_root, "sim.candidate.1v1.kept.final.tsv"), sec = file.path(species_root,
        "sec.candidate.1v1.kept.final.tsv"), ere = file.path(species_root, "ere.candidate.1v1.kept.final.tsv"))

quant_file <- file.path(plot_data_dir, "candidate.editing_level.observed.long.tsv")

version_config <- list(current_1v1_kept = list(label = "Current 1v1 kept", species_source = "02_species_1v1/*.candidate.1v1.kept.final.tsv",
    species_files = c(sim = file.path(species_root, "sim.candidate.1v1.kept.final.tsv"),
        sec = file.path(species_root, "sec.candidate.1v1.kept.final.tsv"), ere = file.path(species_root,
            "ere.candidate.1v1.kept.final.tsv")), source_type = "table"), current_1v1_min2 = list(label = "Current 1v1 min2 sensitivity",
    species_source = "01_rebuild/site_lists/*.1v1.min2.handles.txt", species_files = c(sim = file.path(candidate_root,
        "site_lists/sim_overlap_candidate.1v1.min2.handles.txt"), sec = file.path(candidate_root,
        "site_lists/sec_overlap_candidate.1v1.min2.handles.txt"), ere = file.path(candidate_root,
        "site_lists/ere_overlap_candidate.1v1.min2.handles.txt")), source_type = "list"))

coding_pattern <- paste(c("missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
    "start_lost", "stop_retained_variant", "coding_sequence_variant", "protein_altering_variant",
    "inframe_insertion", "inframe_deletion", "frameshift_variant"), collapse = "|")


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

check_file <- function(path, label = path) {
    if (!file.exists(path) || is.na(file.info(path)$size)) {
        stop("Missing ", label, ": ", path)
    }
    normalizePath(path, mustWork = TRUE)
}

require_columns <- function(x, columns, label) {
    missing <- setdiff(columns, names(x))
    if (length(missing) > 0L) {
        stop(label, " is missing required column(s): ", paste(missing, collapse = ", "))
    }
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

first_nonempty <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x) & toupper(x) != "NA"]
    if (length(x) == 0L)
        NA_character_ else x[[1L]]
}

first_finite <- function(x) {
    x <- safe_numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L)
        NA_real_ else x[[1L]]
}

paste_unique <- function(x, collapse = ";") {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x) & toupper(x) != "NA"]
    if (length(x) == 0L)
        return(NA_character_)
    paste(sort(unique(x)), collapse = collapse)
}

safe_max <- function(x) {
    x <- safe_numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L)
        NA_real_ else max(x)
}

read_handles <- function(path) {
    path <- check_file(path, "site list")
    x <- fread(path, header = FALSE, showProgress = FALSE)
    if (nrow(x) == 0L || ncol(x) == 0L)
        return(character())
    h <- unique(trimws(as.character(x[[1L]])))
    h[!is.na(h) & nzchar(h)]
}

read_table_handles <- function(path) {
    x <- fread(check_file(path), showProgress = FALSE)
    require_columns(x, "handle", basename(path))
    h <- unique(trimws(as.character(x$handle)))
    h[!is.na(h) & nzchar(h)]
}

standardize_sex <- function(x, library) {
    x <- as.character(x)
    library <- as.character(library)
    answer <- rep(NA_character_, length(x))
    answer[tolower(x) %chin% c("f", "female")] <- "F"
    answer[tolower(x) %chin% c("m", "male")] <- "M"

    missing <- is.na(answer)
    answer[missing & grepl("(^|[_-])F([0-9]+|$)|female", library, ignore.case = TRUE)] <- "F"
    answer[missing & grepl("(^|[_-])M([0-9]+|$)|male", library, ignore.case = TRUE)] <- "M"
    answer
}

classify_region <- function(type) {
    type <- as.character(type)
    ifelse(!is.na(type) & grepl(coding_pattern, type, ignore.case = TRUE), "coding", "noncoding")
}


# -----------------------------------------------------------------------------
# Annotation and membership
# -----------------------------------------------------------------------------

read_group_annotation <- function(group_name) {
    path <- check_file(annotation_files[[group_name]], paste0(group_name, " annotation"))
    x <- fread(path, showProgress = FALSE)
    require_columns(x, c("handle", "type"), basename(path))

    if ("change2" %in% names(x)) {
        x <- x[trimws(as.character(change2)) == "A>G"]
    }
    if (nrow(x) == 0L) {
        stop("No A>G annotation remains in ", path)
    }

    x[, handle := as.character(handle)]
    x[, type := as.character(type)]
    x[, region_row := classify_region(type)]

    if (group_name %chin% species_order) {
        require_columns(x, c("chr", "pos"), basename(path))
        x[, motif_chr := as.character(chr)]
        x[, motif_pos := safe_numeric(pos)]
        if ("species_coord" %in% names(x)) {
            x[, motif_handle := as.character(species_coord)]
        } else {
            x[, motif_handle := paste(motif_chr, motif_pos, sep = "_")]
        }
        x[is.na(motif_handle) | !nzchar(motif_handle) | toupper(motif_handle) == "NA", motif_handle :=
            paste(motif_chr, motif_pos, sep = "_")]
        basis <- "species_own_best_AtoG_annotation"
    } else {
        chr_column <- if ("dmel_chr" %in% names(x))
            "dmel_chr" else "chr"
        pos_column <- if ("dmel_pos" %in% names(x))
            "dmel_pos" else "pos"
        require_columns(x, c(chr_column, pos_column), basename(path))
        x[, motif_chr := as.character(get(chr_column))]
        x[, motif_pos := safe_numeric(get(pos_column))]
        x[, motif_handle := handle]
        basis <- if (group_name == "ISO1") {
            "ISO1_authoritative_annotation"
        } else {
            "WBG_Dmel_best_AtoG_annotation"
        }
    }

    answer <- x[!is.na(handle) & nzchar(handle), .(motif_chr = first_nonempty(motif_chr),
        motif_pos = first_finite(motif_pos), motif_handle = first_nonempty(motif_handle),
        type = paste_unique(type), type_all = if ("type_all" %in% names(x)) {
            paste_unique(type_all)
        } else {
            paste_unique(type)
        }, region = if (any(region_row == "coding", na.rm = TRUE)) {
            "coding"
        } else {
            "noncoding"
        }, n_annotation_rows = .N), by = handle]
    answer[, group := group_name]
    answer[, genome_key := unname(group_genome[group_name])]
    answer[, annotation_basis := basis]
    answer[, annotation_file := path]

    setcolorder(answer, c("group", "genome_key", "handle", "motif_handle", "motif_chr", "motif_pos",
        "region", "type", "type_all", "n_annotation_rows", "annotation_basis", "annotation_file"))
    answer[]
}

annotation_map <- rbindlist(lapply(group_order, read_group_annotation), use.names = TRUE,
    fill = TRUE)

if (annotation_map[is.na(motif_chr) | !nzchar(motif_chr) | !is.finite(motif_pos) | motif_pos <
    1, .N] > 0L) {
    fwrite(annotation_map[is.na(motif_chr) | !nzchar(motif_chr) | !is.finite(motif_pos) |
        motif_pos < 1], file.path(output_dir, "annotation_rows_with_invalid_coordinates.tsv"),
        sep = "\t", quote = FALSE)
    warning("Some annotation rows have invalid motif coordinates.")
}

fwrite(annotation_map, file.path(output_dir, "motif_coordinate_annotation_map.tsv"), sep = "\t",
    quote = FALSE)


# -----------------------------------------------------------------------------
# Quantitative table
# -----------------------------------------------------------------------------

quant <- fread(check_file(quant_file, "editing-level long table"), showProgress = FALSE)
require_columns(quant, c("group", "library", "sex", "handle", "altcount"), basename(quant_file))
quant[, group := as.character(group)]
quant[, library := as.character(library)]
quant[, handle := as.character(handle)]
quant[, sex := standardize_sex(sex, library)]
quant[, altcount := safe_numeric(altcount)]
if (!"depth" %in% names(quant)) {
    require_columns(quant, "refcount", basename(quant_file))
    quant[, depth := safe_numeric(refcount) + altcount]
} else {
    quant[, depth := safe_numeric(depth)]
}

quant <- quant[group %chin% group_order & !is.na(library) & nzchar(library) & !is.na(handle) &
    nzchar(handle) & !is.na(sex) & sex %chin% c("F", "M") & is.finite(altcount) & is.finite(depth) &
    depth > 0]

# Match the current Figure 1 workflow: ISO1 is represented by the two head
# libraries, one female and one male.  AdarS/AdarG/Adar5G1 are not part of the
# Dmel motif groups.
quant <- quant[group != "ISO1" | library %chin% c("ISO1_head_F", "ISO1_head_M")]

# Defensive collapse only; the source table is expected to have one row per
# group/library/handle.  It prevents duplicate annotation rows from doubling
# a pooled count if an upstream table is regenerated.
quant <- quant[, .(altcount = safe_max(altcount), depth = safe_max(depth), n_quant_rows = .N),
    by = .(group, library, sex, handle)]

fwrite(quant, file.path(output_dir, "quantitative_motif_table.after_normalization.tsv"),
    sep = "\t", quote = FALSE)


# -----------------------------------------------------------------------------
# Current memberships
# -----------------------------------------------------------------------------

dmel_handles <- read_handles(file.path(candidate_root, "site_lists/candidate.final.handles.txt"))
wbg_handles <- read_table_handles(annotation_files[["WBG"]])

build_membership <- function(version_name, config) {
    fixed_sets <- list(ISO1 = intersect(dmel_handles, annotation_map[group == "ISO1", handle]),
        WBG = intersect(dmel_handles, wbg_handles))

    species_sets <- lapply(species_order, function(species) {
        if (config$source_type == "table") {
            handles <- read_table_handles(config$species_files[[species]])
        } else {
            handles <- read_handles(config$species_files[[species]])
        }
        handles
    })
    names(species_sets) <- species_order

    selected <- c(fixed_sets, species_sets)
    raw_counts <- data.table(group = names(selected), n_source_sites = vapply(selected, function(x) uniqueN(x),
        integer(1)))

    membership <- rbindlist(lapply(names(selected), function(group_name) {
        data.table(version = version_name, group = group_name, handle = sort(unique(selected[[group_name]])))
    }), use.names = TRUE)

    membership <- merge(membership, annotation_map[, .(group, handle)], by = c("group", "handle"),
        all = FALSE, sort = FALSE)
    selected_counts <- membership[, .(n_selected_sites = uniqueN(handle)), by = .(version,
        group)]
    selected_counts <- merge(data.table(version = version_name, group = group_order), selected_counts,
        by = c("version", "group"), all.x = TRUE, sort = FALSE)
    selected_counts <- merge(selected_counts, raw_counts, by = "group", all.x = TRUE, sort = FALSE)
    selected_counts[is.na(n_selected_sites), n_selected_sites := 0L]
    selected_counts[is.na(n_source_sites), n_source_sites := 0L]
    selected_counts[, n_without_annotation := n_source_sites - n_selected_sites]
    selected_counts[, group_order_index := match(group, group_order)]
    setorder(selected_counts, group_order_index)
    selected_counts[, group_order_index := NULL]
    list(membership = membership, counts = selected_counts)
}


# -----------------------------------------------------------------------------
# FASTA random-access extraction from .fai
# -----------------------------------------------------------------------------

revcomp_one <- function(sequence) {
    if (is.na(sequence))
        return(NA_character_)
    sequence <- toupper(sequence)
    paste0(rev(strsplit(chartr("ACGTN", "TGCAN", sequence), "", fixed = TRUE)[[1L]]), collapse = "")
}

read_fai <- function(fasta) {
    fai <- check_file(paste0(fasta, ".fai"), paste0(fasta, ".fai"))
    x <- fread(fai, header = FALSE, showProgress = FALSE)
    if (ncol(x) < 5L)
        stop("Invalid FASTA index: ", fai)
    x <- x[, 1:5]
    setnames(x, c("motif_chr", "chromosome_length", "byte_offset", "line_bases", "line_width"))
    x[, motif_chr := as.character(motif_chr)]
    for (j in c("chromosome_length", "byte_offset", "line_bases", "line_width")) {
        x[, (j) := safe_numeric(get(j))]
    }
    x
}

fai_tables <- lapply(genome_fasta, function(fasta) {
    fasta <- check_file(fasta, "reference FASTA")
    read_fai(fasta)
})

fetch_one_from_fai <- function(con, idx, chromosome, start_pos, width) {
    row <- idx[motif_chr == chromosome]
    if (nrow(row) != 1L || !is.finite(row$chromosome_length)) {
        return(NA_character_)
    }
    if (start_pos < 1 || start_pos + width - 1L > row$chromosome_length) {
        return(NA_character_)
    }

    first_zero <- start_pos - 1
    byte_start <- row$byte_offset + floor(first_zero/row$line_bases) * row$line_width + (first_zero%%row$line_bases)
    nbytes <- width + (ceiling(width/row$line_bases) + 2) * row$line_width + 4

    seek(con, where = as.double(byte_start), origin = "start", rw = "read")
    raw <- readBin(con, what = "raw", n = as.integer(nbytes))
    if (length(raw) == 0L)
        return(NA_character_)
    # Remove FASTA line breaks at the raw-byte level.  This avoids preserving
    # a newline when a requested three-base window crosses a FASTA line wrap.
    raw <- raw[!(raw %in% as.raw(c(10L, 13L)))]
    text <- rawToChar(raw)
    text <- toupper(substr(text, 1L, width))
    if (nchar(text) != width)
        NA_character_ else text
}

extract_one_genome <- function(coordinates, genome_key) {
    if (nrow(coordinates) == 0L)
        return(data.table())
    fasta <- genome_fasta[[genome_key]]
    idx <- fai_tables[[genome_key]]
    con <- file(fasta, open = "rb")
    on.exit(close(con), add = TRUE)

    x <- copy(coordinates)
    x[, motif_chr := as.character(motif_chr)]
    x[, motif_pos := safe_numeric(motif_pos)]
    x[, chromosome_length := idx$chromosome_length[match(motif_chr, idx$motif_chr)]]
    x[, faidx_region := ifelse(is.finite(motif_pos), paste0(motif_chr, ":", as.integer(motif_pos -
        flank), "-", as.integer(motif_pos + flank)), NA_character_)]
    x[, raw_motif := vapply(seq_len(.N), function(i) {
        if (!is.finite(motif_pos[[i]]))
            return(NA_character_)
        fetch_one_from_fai(con, idx, motif_chr[[i]], as.integer(motif_pos[[i]] - flank),
            2L * flank + 1L)
    }, character(1))]

    center_index <- flank + 1L
    motif_width <- 2L * flank + 1L
    x[, raw_motif := toupper(raw_motif)]
    x[, raw_center := substr(raw_motif, center_index, center_index)]
    x[, orientation := fifelse(raw_center == "A", "+", fifelse(raw_center == "T", "-",
        "invalid"))]
    x[, motif := raw_motif]
    x[orientation == "-", motif := vapply(raw_motif, revcomp_one, character(1))]
    x[, center := substr(motif, center_index, center_index)]
    x[, motif_valid := !is.na(motif) & nchar(motif) == motif_width & grepl("^[ACGT]+$",
        motif) & center == "A"]
    x[, invalid_reason := fcase(is.na(chromosome_length), "chromosome_not_in_fai", !is.finite(motif_pos),
        "invalid_position", motif_pos - flank < 1L | motif_pos + flank > chromosome_length,
        "out_of_bounds", is.na(raw_motif) | !nzchar(raw_motif), "sequence_not_extracted",
        nchar(raw_motif) != motif_width, "wrong_sequence_length", !grepl("^[ACGT]+$", raw_motif),
        "non_ACGT_sequence", !raw_center %chin% c("A", "T"), "center_not_A_or_T", center !=
            "A", "orientation_failed", default = "")]
    x[, .(genome_key, motif_chr, motif_pos, motif_handle, chromosome_length, faidx_region,
        raw_motif, raw_center, orientation, motif, center, motif_valid, invalid_reason)]
}

extract_motifs <- function(high_sites) {
    if (nrow(high_sites) == 0L)
        return(data.table())
    coordinates <- unique(high_sites[, .(genome_key, motif_chr, motif_pos, motif_handle)])
    rbindlist(lapply(names(genome_fasta), function(gk) {
        extract_one_genome(coordinates[coordinates$genome_key == gk], gk)
    }), use.names = TRUE, fill = TRUE)
}


# -----------------------------------------------------------------------------
# Sequence logos with ggseqlogo
# -----------------------------------------------------------------------------

dna_colors <- c(A = "#D62728", C = "#1F5CB8", G = "#F2A900", T = "#17823B")

dna_col_scheme <- make_col_scheme(chars = names(dna_colors), groups = names(dna_colors),
    cols = unname(dna_colors), name = "ADAR DNA")

make_empty_panel <- function(title_text, mode, show_y) {
    ymax <- if (mode == "bits")
        2 else 1

    p <- ggplot() + annotate("text", x = flank + 1L, y = ymax/2, label = "No sites", size = 3) +
        scale_x_continuous(limits = c(0.5, 2L * flank + 1.5), breaks = seq_len(2L * flank +
            1L), labels = seq_len(2L * flank + 1L) - (flank + 1L), expand = expansion(mult = c(0.02,
            0.02))) + scale_y_continuous(limits = c(0, ymax), breaks = if (mode == "bits") {
        c(0, 1, 2)
    } else {
        seq(0, 1, by = 0.25)
    }, expand = expansion(mult = c(0, 0))) + labs(title = paste0(title_text, "\n(n = 0)"),
        x = "Position relative to edited A", y = if (mode == "bits")
            "bits" else "probability") + theme_classic(base_size = 8, base_family = font_family) + theme(plot.title = element_text(hjust = 0.5,
        size = 7), axis.text = element_text(color = "black", size = 6), axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.35), plot.margin = margin(3, 3, 3, 3))

    if (!show_y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    p
}

make_logo_panel <- function(sequences, title_text, mode, show_y) {
    sequences <- toupper(as.character(sequences))
    sequences <- sequences[!is.na(sequences) & grepl("^[ACGT]+$", sequences)]

    if (length(sequences) == 0L) {
        return(make_empty_panel(title_text, mode, show_y))
    }

    widths <- unique(nchar(sequences))
    if (length(widths) != 1L) {
        stop("Logo sequences have different widths.")
    }

    width <- widths[[1L]]
    ymax <- if (mode == "bits")
        2 else 1
    logo_method <- if (mode == "bits")
        "bits" else "prob"

    p <- ggseqlogo(sequences, method = logo_method, col_scheme = dna_col_scheme) + scale_x_continuous(breaks = seq_len(width),
        labels = seq_len(width) - (flank + 1L), expand = expansion(mult = c(0.02, 0.02))) +
        scale_y_continuous(limits = c(0, ymax), breaks = if (mode == "bits") {
            c(0, 1, 2)
        } else {
            seq(0, 1, by = 0.25)
        }, expand = expansion(mult = c(0, 0))) + labs(title = paste0(title_text, "\n(n = ",
        length(sequences), ")"), x = "Position relative to edited A", y = if (mode == "bits")
        "bits" else "probability") + coord_cartesian(xlim = c(0.5, width + 0.5), ylim = c(0, ymax),
        clip = "off") + theme_classic(base_size = 8, base_family = font_family) + theme(plot.title = element_text(hjust = 0.5,
        size = 7), axis.text = element_text(color = "black", size = 6), axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.35), legend.position = "none", plot.margin = margin(3,
            3, 3, 3))

    if (!show_y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    p
}

make_logo_grid <- function(plot_sites, sex_scope, mode, version_label) {
    region_order <- c("All high-edited sites", "Coding high-edited sites", "Noncoding high-edited sites")

    panels <- vector("list", length(region_order) * length(group_order))

    panel_index <- 1L

    for (region_name in region_order) {
        for (group_name in group_order) {
            subset <- plot_sites[sex == sex_scope & region_class == region_name & group ==
                group_name]

            title <- paste(unname(group_display[group_name]), region_name, sep = "\n")

            panels[[panel_index]] <- make_logo_panel(subset$motif, title, mode, show_y = group_name ==
                group_order[[1L]])

            panel_index <- panel_index + 1L
        }
    }

    arrangeGrob(grobs = panels, ncol = length(group_order), top = textGrob(paste0("ADAR motif, ",
        version_label, ", ", sex_scope, " (", ifelse(mode == "bits", "bits", "probability"),
        ")"), gp = gpar(fontface = "bold", fontsize = 13, fontfamily = font_family)))
}

save_logo <- function(grob_object, pdf_path, png_path, width = 17, height = 9) {
    if (capabilities("cairo")) {
        cairo_pdf(pdf_path, width = width, height = height, family = font_family)
    } else {
        pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
    }

    grid.newpage()
    grid.draw(grob_object)
    invisible(dev.off())

    if (capabilities("cairo")) {
        png(png_path, width = width, height = height, units = "in", res = 300, type = "cairo")
    } else {
        png(png_path, width = width, height = height, units = "in", res = 300)
    }

    grid.newpage()
    grid.draw(grob_object)
    invisible(dev.off())
}


# -----------------------------------------------------------------------------
# Run one membership version
# -----------------------------------------------------------------------------

run_version <- function(version_name, config) {
    message("\n[version] ", version_name)
    version_dir <- file.path(output_dir, version_name)
    dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)

    membership_result <- build_membership(version_name, config)
    membership <- membership_result$membership
    selected_counts <- membership_result$counts
    fwrite(membership, file.path(version_dir, "selected_candidate_sites.tsv"), sep = "\t",
        quote = FALSE)
    fwrite(selected_counts, file.path(version_dir, "selected_candidate_site_counts.tsv"),
        sep = "\t", quote = FALSE)
    message("  selected sites: ", paste(selected_counts[, paste0(group, "=", n_selected_sites)],
        collapse = ", "))

    replicate_rows <- merge(quant, membership[, .(group, handle)], by = c("group", "handle"),
        all = FALSE, sort = FALSE)
    replicate_rows <- merge(replicate_rows, annotation_map, by = c("group", "handle"), all.x = TRUE,
        sort = FALSE)
    replicate_rows[, editing_level_library := altcount/depth]
    fwrite(replicate_rows, file.path(version_dir, "replicate_level_motif_input.tsv"), sep = "\t",
        quote = FALSE)

    group_pool <- replicate_rows[, .(altcount = sum(altcount, na.rm = TRUE), depth = sum(depth,
        na.rm = TRUE), n_replicates = uniqueN(library), motif_handle = first_nonempty(motif_handle),
        motif_chr = first_nonempty(motif_chr), motif_pos = first_finite(motif_pos), region = first_nonempty(region),
        type = first_nonempty(type), type_all = first_nonempty(type_all), annotation_basis = first_nonempty(annotation_basis)),
        by = .(group, genome_key, handle)]
    group_pool[, sex := "ALL"]
    group_pool[, editing_level := altcount/depth]

    sex_pool <- replicate_rows[, .(altcount = sum(altcount, na.rm = TRUE), depth = sum(depth,
        na.rm = TRUE), n_replicates = uniqueN(library), motif_handle = first_nonempty(motif_handle),
        motif_chr = first_nonempty(motif_chr), motif_pos = first_finite(motif_pos), region = first_nonempty(region),
        type = first_nonempty(type), type_all = first_nonempty(type_all), annotation_basis = first_nonempty(annotation_basis)),
        by = .(group, genome_key, sex, handle)]
    sex_pool[, editing_level := altcount/depth]

    high_sites <- rbindlist(list(group_pool[depth >= min_depth & editing_level >= min_editing],
        sex_pool[sex %chin% c("F", "M") & depth >= min_depth & editing_level >= min_editing]),
        use.names = TRUE, fill = TRUE)
    high_sites[, version := version_name]
    fwrite(high_sites, file.path(version_dir, "high_editing_sites.before_motif.tsv"), sep = "\t",
        quote = FALSE)

    motif_info <- extract_motifs(high_sites)
    high_motif <- merge(high_sites, motif_info, by = c("genome_key", "motif_chr", "motif_pos",
        "motif_handle"), all.x = TRUE, sort = FALSE)
    fwrite(high_motif, file.path(version_dir, "high_editing_sites.with_motif.tsv"), sep = "\t",
        quote = FALSE)
    fwrite(high_motif[is.na(motif_valid) | motif_valid == FALSE], file.path(version_dir,
        "invalid_or_unextracted_sites.tsv"), sep = "\t", quote = FALSE)

    valid_sites <- high_motif[motif_valid == TRUE]
    plot_sites <- rbindlist(list(copy(valid_sites)[, region_class := "All high-edited sites"],
        copy(valid_sites[region == "coding"])[, region_class := "Coding high-edited sites"],
        copy(valid_sites[region == "noncoding"])[, region_class := "Noncoding high-edited sites"]),
        use.names = TRUE, fill = TRUE)
    fwrite(plot_sites, file.path(version_dir, "motif_plot_sites.tsv"), sep = "\t", quote = FALSE)

    logo_counts <- plot_sites[, .(n_sites = uniqueN(handle), n_motif_rows = .N), by = .(sex,
        region_class, group)]
    logo_counts[, version := version_name]
    fwrite(logo_counts, file.path(version_dir, "motif_logo_site_counts.tsv"), sep = "\t",
        quote = FALSE)

    modes <- if (plot_mode == "both")
        c("bits", "prob") else plot_mode
    for (sex_scope in sex_order) {
        for (mode in modes) {
            logo <- make_logo_grid(plot_sites, sex_scope, mode, config$label)
            mode_tag <- if (mode == "bits")
                "bits" else "probability"
            prefix <- paste0("ADAR_motif.", version_name, ".", sex_scope, ".", mode_tag)
            save_logo(logo, file.path(version_dir, paste0(prefix, ".pdf")), file.path(version_dir,
                paste0(prefix, ".png")))
        }
    }

    high_counts <- data.table(group = group_order)
    high_counts <- merge(high_counts, high_sites[, .(n_high_sites_ALL = uniqueN(handle[sex ==
        "ALL"]), n_high_sites_F = uniqueN(handle[sex == "F"]), n_high_sites_M = uniqueN(handle[sex ==
        "M"])), by = group], by = "group", all.x = TRUE, sort = FALSE)
    valid_counts <- data.table(group = group_order)
    valid_counts <- merge(valid_counts, valid_sites[, .(n_valid_motifs_ALL = uniqueN(handle[sex ==
        "ALL"]), n_valid_motifs_F = uniqueN(handle[sex == "F"]), n_valid_motifs_M = uniqueN(handle[sex ==
        "M"])), by = group], by = "group", all.x = TRUE, sort = FALSE)
    for (x in list(high_counts, valid_counts)) {
        invisible(x)
    }
    for (column in c("n_high_sites_ALL", "n_high_sites_F", "n_high_sites_M")) {
        if (!column %in% names(high_counts))
            high_counts[, (column) := 0L]
        set(high_counts, which(is.na(high_counts[[column]])), column, 0L)
    }
    for (column in c("n_valid_motifs_ALL", "n_valid_motifs_F", "n_valid_motifs_M")) {
        if (!column %in% names(valid_counts))
            valid_counts[, (column) := 0L]
        set(valid_counts, which(is.na(valid_counts[[column]])), column, 0L)
    }

    result <- merge(selected_counts, high_counts, by = "group", all.x = TRUE, sort = FALSE)
    result <- merge(result, valid_counts, by = "group", all.x = TRUE, sort = FALSE)
    result[, version_label := config$label]
    result[, version := version_name]
    result[]
}


# -----------------------------------------------------------------------------
# Execute and write audit configuration
# -----------------------------------------------------------------------------

summaries <- rbindlist(lapply(names(version_config), function(version_name) {
    run_version(version_name, version_config[[version_name]])
}), use.names = TRUE, fill = TRUE)
summaries[, group := factor(group, levels = group_order)]
setorder(summaries, version, group)
summaries[, group := as.character(group)]
fwrite(summaries, file.path(output_dir, "motif_1v1_summary.tsv"), sep = "\t", quote = FALSE)

run_config <- data.table(item = c("project_root", "candidate_root", "species_root", "quantitative_table",
    "reference_dir", "output_dir", "minimum_pooled_depth", "minimum_pooled_editing_level",
    "plot_mode", "flank", "motif_width", "ISO1_libraries", "Dmel_mutant_groups_excluded",
    "groups", "primary_version", "sequence_extraction"), value = c(project_root, candidate_root,
    species_root, quant_file, reference_dir, output_dir, min_depth, min_editing, plot_mode,
    flank, 2L * flank + 1L, "ISO1_head_F;ISO1_head_M", "Adar5G1;AdarS;AdarG", paste(group_order,
        collapse = ";"), "current_1v1_kept", "R .fai random access; no samtools/ggseqlogo dependency"))
fwrite(run_config, file.path(output_dir, "motif_run_configuration.tsv"), sep = "\t", quote = FALSE)

message("\nDone.")
message("Summary: ", file.path(output_dir, "motif_1v1_summary.tsv"))
