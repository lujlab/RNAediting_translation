#!/usr/bin/env Rscript

# ============================================================
# Table 1 allele-specific RNA/Ribo analysis
#
# Analysis groups
# ----------------
# Dmel = ISO1 + WBg only.
# Four species are analyzed:
#   Dmel, Dsim, Dsec, Dere
#
# Each species is analyzed in:
#   (1) females only
#   (2) males only
#   (3) females + males combined
#
# Biological replicates are NEVER pooled at the read-count level.
# Independent "head" samples are retained as biological replicates.
#
# Statistical analyses
# --------------------
# Primary:
#   edgeR negative-binomial robust quasi-likelihood model
# Beta-binomial sensitivity analysis is not run in this two-method entry.
# Legacy/reference:
#   one-sided Fisher in both directions, all eligible replicates; joint directional BH
#
# Biological hypothesis
# ---------------------
# Does the edited/unedited allele ratio change from RNA-seq to Ribo-seq?
#
# effect_log2OR =
#   log2[(Ribo_I / Ribo_A) / (RNA_I / RNA_A)]
#
# effect > 0  : I>A, edited allele enriched in Ribo-seq
# effect < 0  : I<A, edited allele depleted in Ribo-seq
#
# A-site and P-site are tested independently.
# AP is constructed only after testing:
#   AP significant = A significant OR P significant
# A/P read counts are NEVER summed.
# ============================================================

suppressPackageStartupMessages({
    for (pkg in c("data.table", "ggplot2")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(
                "Missing R package '", pkg, "'. Install with: install.packages('", pkg, "')",
                call. = FALSE
            )
        }
    }

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop(
            paste0(
                "Missing Bioconductor package 'edgeR'. Install once with:\n",
                "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
                "  BiocManager::install('edgeR')"
            ),
            call. = FALSE
        )
    }
})

library(data.table)
library(ggplot2)

# ============================================================
# 0. Arguments
# ============================================================

parse_args <- function(args) {
    out <- list(
        input_dir = "/path/to/num/Table1_alledit",
        output_dir = NULL,
        sample_meta = NULL,
        threads = 8L,
        fdr = 0.05,

        # Common site-level eligibility.  A site must meet these coverage
        # thresholds in at least min_replicates independent replicates.
        min_rna_per_rep_equiv = 5,
        min_ribo_per_rep_equiv = 5,
        min_lm = 0.05,

        # Analyze groups/sites supported by at least 2 biological replicates.
        min_replicates = 2L,

        # Pseudocount is ONLY used to calculate descriptive Fisher log2OR.
        fisher_pc = 0.5
    )

    value_options <- c(
        "--input-dir", "--output-dir", "--sample-meta",
        "--threads", "--fdr",
        "--min-rna-per-rep-equiv", "--min-ribo-per-rep-equiv",
        "--min-lm", "--min-replicates", "--fisher-pc"
    )

    i <- 1L
    while (i <= length(args)) {
        a <- args[[i]]

        if (a == "--help") {
            cat(
                "Usage:\n",
                "  Rscript Table1_DmelISO1WBg_4species_edgeR_BetaBin_Fisher.R \\\n",
                "    --input-dir DIR \\\n",
                "    [--output-dir DIR] \\\n",
                "    [--sample-meta FILE] \\\n",
                "    [--threads 8]\n\n",
                "Optional sample_meta.tsv columns:\n",
                "  sample_id, species, sex, background, display_name\n",
                "species may be mel/sim/sec/ere or Dmel/Dsim/Dsec/Dere.\n",
                "background is only used to restrict Dmel to ISO1 + WBg.\n",
                sep = ""
            )
            quit(save = "no", status = 0L)
        }

        if (!a %in% value_options) {
            stop("Unknown option: ", a, call. = FALSE)
        }
        if (i == length(args)) stop("Missing value after ", a, call. = FALSE)

        key <- gsub("-", "_", sub("^--", "", a))
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
    }

    out$threads <- as.integer(out$threads)
    out$fdr <- as.numeric(out$fdr)
    out$min_rna_per_rep_equiv <- as.numeric(out$min_rna_per_rep_equiv)
    out$min_ribo_per_rep_equiv <- as.numeric(out$min_ribo_per_rep_equiv)
    out$min_lm <- as.numeric(out$min_lm)
    out$min_replicates <- as.integer(out$min_replicates)
    out$fisher_pc <- as.numeric(out$fisher_pc)

    if (is.null(out$output_dir)) {
        out$output_dir <- file.path(
            out$input_dir,
            "Table1_DmelISO1WBg_4species_edgeR_BetaBin_Fisher"
        )
    }

    if (is.null(out$sample_meta)) {
        out$sample_meta <- file.path(out$input_dir, "sample_info.tsv")
    }

    out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
setDTthreads(opt$threads)

message("[INFO] input_dir  = ", opt$input_dir)
message("[INFO] output_dir = ", opt$output_dir)
message("[INFO] Dmel = ISO1 + WBg only.")
message("[INFO] Head samples are retained as independent biological replicates.")
message("[INFO] Minimum biological replicates per analysis group = ", opt$min_replicates)

# ============================================================
# 1. Display order
# ============================================================

analysis_group_order <- c(
    "Dmel_F", "Dmel_M",
    "Dsim_F", "Dsim_M",
    "Dsec_F", "Dsec_M",
    "Dere_F", "Dere_M",
    "Dmel", "Dsim", "Dsec", "Dere"
)

analysis_group_labels <- c(
    Dmel_F = "Dmel F",
    Dmel_M = "Dmel M",
    Dsim_F = "Dsim F",
    Dsim_M = "Dsim M",
    Dsec_F = "Dsec F",
    Dsec_M = "Dsec M",
    Dere_F = "Dere F",
    Dere_M = "Dere M",
    Dmel = "Dmel",
    Dsim = "Dsim",
    Dsec = "Dsec",
    Dere = "Dere"
)

# ============================================================
# 2. Helpers: metadata and files
# ============================================================

first_non_na <- function(x) {
    if (is.character(x) || is.factor(x)) {
        x <- as.character(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (!length(x)) return(NA_character_)
        return(x[[1L]])
    }

    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    x[[1L]]
}

normalize_species <- function(x) {
    s <- tolower(as.character(x))
    out <- rep(NA_character_, length(s))

    out[s %in% c("mel", "dmel", "d. melanogaster", "dmelanogaster")] <- "mel"
    out[s %in% c("sim", "dsim", "d. simulans", "dsimulans")] <- "sim"
    out[s %in% c("sec", "dsec", "d. sechellia", "dsechellia")] <- "sec"
    out[s %in% c("ere", "dere", "d. erecta", "derecta")] <- "ere"

    out
}

species_display <- function(x) {
    map <- c(mel = "Dmel", sim = "Dsim", sec = "Dsec", ere = "Dere")
    unname(map[as.character(x)])
}

infer_position_mode <- function(file) {
    b <- basename(file)

    if (grepl("Psite", b, ignore.case = TRUE)) return("P")
    if (grepl("Asite", b, ignore.case = TRUE)) return("A")

    # Old files without explicit Asite in the name were A-site files.
    if (grepl("ribo\\.codon_df_alledit", b, ignore.case = TRUE)) return("A")

    NA_character_
}

infer_sample_id <- function(file) {
    b <- basename(file)

    b <- sub("\\.gz$", "", b, ignore.case = TRUE)
    b <- sub("\\.(txt|tsv|csv)$", "", b, ignore.case = TRUE)

    suffix_patterns <- c(
        "\\.ribo\\.codon_Psite_df_alledit$",
        "\\.ribo\\.codon_Asite_df_alledit$",
        "\\.ribo\\.codon_df_alledit$",
        "\\.codon_Psite_df_alledit$",
        "\\.codon_Asite_df_alledit$",
        "\\.Psite_df_alledit$",
        "\\.Asite_df_alledit$",
        "_Psite_df_alledit$",
        "_Asite_df_alledit$"
    )

    for (pat in suffix_patterns) {
        b <- sub(pat, "", b, ignore.case = TRUE)
    }

    b
}

infer_sex <- function(sample_id) {
    s <- tolower(sample_id)

    if (grepl("female", s)) return("F")
    if (grepl("male", s)) return("M")

    if (grepl("(^|[_.-])f[0-9]*($|[_.-])", s)) return("F")
    if (grepl("(^|[_.-])m[0-9]*($|[_.-])", s)) return("M")

    NA_character_
}

infer_species <- function(sample_id) {
    s <- tolower(sample_id)

    if (grepl("iso1|iso-1|wbg|w1118|adarg|adars|adar5g1|dmel|mel", s)) {
        return("mel")
    }
    if (grepl("dsim|(^|[_.-])sim([_.-]|$)", s)) return("sim")
    if (grepl("dsec|(^|[_.-])sec([_.-]|$)", s)) return("sec")
    if (grepl("dere|(^|[_.-])ere([_.-]|$)", s)) return("ere")

    NA_character_
}

infer_background <- function(sample_id) {
    s <- tolower(sample_id)

    if (grepl("iso1|iso-1", s)) return("ISO1")
    if (grepl("wbg", s)) return("WBg")
    if (grepl("w1118", s)) return("w1118")
    if (grepl("adarg", s)) return("AdarG")
    if (grepl("adars", s)) return("AdarS")
    if (grepl("adar5g1", s)) return("Adar5G1")

    sp <- infer_species(sample_id)
    if (!is.na(sp) && sp %in% c("sim", "sec", "ere")) return(species_display(sp))

    NA_character_
}

auto_metadata <- function(sample_ids) {
    rbindlist(lapply(sample_ids, function(sid) {
        data.table(
            sample_id = sid,
            species = infer_species(sid),
            sex = infer_sex(sid),
            background = infer_background(sid),
            display_name = sid
        )
    }))
}

# ============================================================
# 3. Helpers: codon annotation
# ============================================================

standardize_codon <- function(x) {
    x <- toupper(as.character(x))
    x <- gsub("U", "T", x)
    x <- gsub("[^ACGT]", "N", x)
    x[nchar(x) != 3] <- NA_character_
    x
}

genetic_code <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

translate_codon <- function(x) {
    codon <- standardize_codon(x)
    aa <- unname(genetic_code[codon])

    bad <- is.na(codon) | grepl("N", codon)
    aa[bad] <- NA_character_

    aa
}

infer_NS <- function(ref_codon, alt_codon, annotation_type = NA_character_) {
    ref_aa <- translate_codon(ref_codon)
    alt_aa <- translate_codon(alt_codon)

    out <- rep(NA_character_, length(ref_aa))

    ok <- !is.na(ref_aa) & !is.na(alt_aa)
    out[ok & ref_aa == alt_aa] <- "S"
    out[ok & ref_aa != alt_aa] <- "N"

    # Annotation fallback.
    miss <- is.na(out)
    if (any(miss)) {
        t <- tolower(as.character(annotation_type))

        out[
            miss &
                grepl("nonsyn|non-syn|missense|recoding", t)
        ] <- "N"

        out[
            miss &
                grepl("synonymous|\\bsyn\\b", t) &
                !grepl("nonsyn|non-syn|missense|recoding", t)
        ] <- "S"
    }

    out
}

format_NS <- function(total, nN, nS) {
    paste0(total, " (", nN, "/", nS, ")")
}

# ============================================================
# 4. Statistical helpers
# ============================================================

safe_fisher <- function(rna_ref, rna_alt, ribo_ref, ribo_alt) {
    x <- c(rna_ref, rna_alt, ribo_ref, ribo_alt)

    if (any(!is.finite(x)) || any(x < 0)) return(NA_real_)

    if ((rna_ref + rna_alt) <= 0 || (ribo_ref + ribo_alt) <= 0) {
        return(NA_real_)
    }

    tab <- matrix(
        c(
            ribo_alt, ribo_ref,
            rna_alt, rna_ref
        ),
        nrow = 2,
        byrow = TRUE
    )

    # Fisher requires integer counts; make R's original ties-to-even
    # rounding explicit for half-integer RNA counts. edgeR retains fractions.
    tab <- round(tab)
    tryCatch(
        fisher.test(tab)$p.value,
        error = function(e) NA_real_
    )
}

fisher_combine_p <- function(p) {
    p <- p[is.finite(p) & !is.na(p)]

    if (!length(p)) return(NA_real_)

    p <- pmin(pmax(p, 1e-300), 1)

    stat <- -2 * sum(log(p))

    pchisq(
        stat,
        df = 2 * length(p),
        lower.tail = FALSE
    )
}

# ============================================================
# 5. Detect A/P input files
# ============================================================

paths <- list.files(
    opt$input_dir,
    full.names = TRUE,
    recursive = FALSE
)

path_info <- file.info(paths)

files <- paths[
    !is.na(path_info$isdir) &
        path_info$isdir == FALSE &
        !is.na(path_info$size) &
        path_info$size > 0
]

file_info <- data.table(
    file = files,
    basename = basename(files)
)

file_info[, position := vapply(
    file,
    infer_position_mode,
    character(1L)
)]

file_info <- file_info[!is.na(position)]

file_info[, sample_id := vapply(
    file,
    infer_sample_id,
    character(1L)
)]

if (!nrow(file_info)) {
    stop("No A-site/P-site allele-count files detected in: ", opt$input_dir)
}

# When two filenames represent the same sample/position, keep only one.
# Head samples themselves are NOT removed.
file_info[, is_codon_file := grepl("codon_", basename, ignore.case = TRUE)]
file_info[, file_rank := fifelse(is_codon_file, 2L, 1L)]

setorder(
    file_info,
    sample_id,
    position,
    file_rank,
    basename
)

file_info <- file_info[
    ,
    .SD[1L],
    by = .(sample_id, position)
]

file_info[, c("is_codon_file", "file_rank") := NULL]

# ============================================================
# 6. Sample metadata
# ============================================================

if (file.exists(opt$sample_meta)) {
    sample_meta <- fread(opt$sample_meta)

    if (!"sample_id" %in% names(sample_meta)) {
        stop("sample_info.tsv must contain sample_id.")
    }

    if (!"species" %in% names(sample_meta)) {
        sample_meta[, species := vapply(
            sample_id,
            infer_species,
            character(1L)
        )]
    } else {
        sample_meta[, species := normalize_species(species)]
    }

    if (!"sex" %in% names(sample_meta)) {
        sample_meta[, sex := vapply(
            sample_id,
            infer_sex,
            character(1L)
        )]
    }

    if (!"background" %in% names(sample_meta)) {
        sample_meta[, background := vapply(
            sample_id,
            infer_background,
            character(1L)
        )]
    }

    if (!"display_name" %in% names(sample_meta)) {
        sample_meta[, display_name := sample_id]
    }
} else {
    sample_meta <- auto_metadata(
        unique(file_info$sample_id)
    )
}

file_info <- merge(
    file_info,
    sample_meta,
    by = "sample_id",
    all.x = TRUE
)

# Select the four study species.
file_info <- file_info[
    species %in% c("mel", "sim", "sec", "ere")
]

# Dmel = ISO1 + WBg only.
file_info <- file_info[
    species != "mel" |
        background %in% c("ISO1", "WBg")
]

# F/M must be known for both sex-specific and combined analyses.
file_info <- file_info[
    sex %in% c("F", "M")
]

if (!nrow(file_info)) {
    stop("No requested ISO1/WBg/Dsim/Dsec/Dere samples remained after filtering.")
}

fwrite(
    unique(
        file_info[
            ,
            .(
                sample_id,
                species,
                sex,
                background,
                display_name
            )
        ]
    ),
    file.path(opt$output_dir, "included_independent_samples.tsv"),
    sep = "\t"
)

message("[INFO] Included samples:")
print(
    unique(file_info[, .(sample_id, species, sex, background)])
)

# ============================================================
# 7. Read allele counts
# ============================================================

read_one_count_file <- function(
  file,
  sample_id,
  species,
  sex,
  background,
  position
) {
    x <- fread(file)
    setDT(x)

    count_cols <- c(
        "ref_number_rna",
        "alt_number_rna",
        "ref_number_ribo",
        "alt_number_ribo"
    )

    missing_count_cols <- setdiff(
        count_cols,
        names(x)
    )

    if (length(missing_count_cols)) {
        stop(
            basename(file),
            " missing count columns: ",
            paste(missing_count_cols, collapse = ", ")
        )
    }

    x[, (count_cols) := lapply(
        .SD,
        as.numeric
    ), .SDcols = count_cols]

    for (cc in count_cols) {
        x[is.na(get(cc)), (cc) := 0]
    }

    if (!"FBt_cdspos" %in% names(x)) x[, FBt_cdspos := NA_character_]
    if (!"handle" %in% names(x)) x[, handle := NA_character_]
    if (!"ref_codon" %in% names(x)) x[, ref_codon := NA_character_]
    if (!"alt_codon" %in% names(x)) x[, alt_codon := NA_character_]
    if (!"edit_pos" %in% names(x)) x[, edit_pos := NA_real_]
    if (!"type" %in% names(x)) x[, type := NA_character_]
    if (!"gene" %in% names(x)) x[, gene := NA_character_]

    x[, site_id := as.character(FBt_cdspos)]

    x[
        is.na(site_id) | !nzchar(site_id),
        site_id := as.character(handle)
    ]

    # Last-resort unique identifier.
    x[
        is.na(site_id) | !nzchar(site_id),
        site_id := paste0(sample_id, "__", position, "__row", .I)
    ]

    x[, NS_class := infer_NS(
        ref_codon,
        alt_codon,
        type
    )]

    out <- x[
        ,
        .(
            sample_id = sample_id,
            species = species,
            sex = sex,
            background = background,
            position = position,
            site_id = site_id,
            handle = as.character(handle),
            FBt_cdspos = as.character(FBt_cdspos),
            gene = as.character(gene),
            ref_codon = as.character(ref_codon),
            alt_codon = as.character(alt_codon),
            edit_pos = as.numeric(edit_pos),
            NS_class = NS_class,
            RNA_ref = ref_number_rna,
            RNA_alt = alt_number_rna,
            Ribo_ref = ref_number_ribo,
            Ribo_alt = alt_number_ribo
        )
    ]

    # Sum duplicate site rows within one independent sample.
    out <- out[
        ,
        .(
            species = first_non_na(species),
            sex = first_non_na(sex),
            background = first_non_na(background),
            handle = first_non_na(handle),
            FBt_cdspos = first_non_na(FBt_cdspos),
            gene = first_non_na(gene),
            ref_codon = first_non_na(ref_codon),
            alt_codon = first_non_na(alt_codon),
            edit_pos = first_non_na(edit_pos),
            NS_class = first_non_na(NS_class),
            RNA_ref = sum(RNA_ref, na.rm = TRUE),
            RNA_alt = sum(RNA_alt, na.rm = TRUE),
            Ribo_ref = sum(Ribo_ref, na.rm = TRUE),
            Ribo_alt = sum(Ribo_alt, na.rm = TRUE)
        ),
        by = .(
            sample_id,
            position,
            site_id
        )
    ]

    out[]
}

message("[INFO] Reading allele-count files ...")

count_raw <- rbindlist(
    lapply(seq_len(nrow(file_info)), function(i) {
        read_one_count_file(
            file = file_info$file[[i]],
            sample_id = file_info$sample_id[[i]],
            species = file_info$species[[i]],
            sex = file_info$sex[[i]],
            background = file_info$background[[i]],
            position = file_info$position[[i]]
        )
    }),
    use.names = TRUE,
    fill = TRUE
)

# ============================================================
# 8. Build species and sex analysis groups
# ============================================================

# Sex-specific copy.
count_sex <- copy(count_raw)

count_sex[, analysis_level := "sex_specific"]

count_sex[, analysis_group := paste0(
    species_display(species),
    "_",
    sex
)]

# Female + male combined copy.
# Raw counts are still kept separate for every biological replicate.
count_all <- copy(count_raw)

count_all[, analysis_level := "sex_combined"]

count_all[, analysis_group := species_display(species)]

count_long <- rbindlist(
    list(
        count_sex,
        count_all
    ),
    use.names = TRUE,
    fill = TRUE
)

count_long[, analysis_group := factor(
    analysis_group,
    levels = analysis_group_order
)]

count_long <- count_long[
    !is.na(analysis_group)
]

count_long[, analysis_group := as.character(analysis_group)]

# Sample composition for verification.
composition <- unique(
    count_long[
        ,
        .(
            analysis_group,
            analysis_level,
            sample_id,
            species,
            sex,
            background
        )
    ]
)

composition_summary <- composition[
    ,
    .(
        n_replicates = uniqueN(sample_id),
        samples = paste(
            sort(unique(sample_id)),
            collapse = ";"
        )
    ),
    by = .(
        analysis_group,
        analysis_level
    )
]

composition_summary[, analysis_group := factor(
    analysis_group,
    levels = analysis_group_order
)]

setorder(
    composition_summary,
    analysis_group
)

composition_summary[, analysis_group := as.character(analysis_group)]

fwrite(
    composition,
    file.path(opt$output_dir, "analysis_sample_composition.tsv"),
    sep = "\t"
)

fwrite(
    composition_summary,
    file.path(opt$output_dir, "analysis_group_replicate_counts.tsv"),
    sep = "\t"
)

message("[INFO] Analysis group replicate counts:")
print(composition_summary)

# ============================================================
# 9. Build matrices for one group + one ribosomal position
# ============================================================

build_group_matrices <- function(x) {
    reps <- sort(unique(x$sample_id))
    sites <- sort(unique(x$site_id))

    make_matrix <- function(value_col) {
        w <- dcast(
            x,
            site_id ~ sample_id,
            value.var = value_col,
            fun.aggregate = sum,
            fill = 0
        )

        w <- w[
            match(sites, site_id)
        ]

        missing_reps <- setdiff(
            reps,
            names(w)
        )

        for (r in missing_reps) {
            w[, (r) := 0]
        }

        mat <- as.matrix(
            w[, ..reps]
        )

        storage.mode(mat) <- "numeric"
        rownames(mat) <- sites

        mat
    }

    meta <- x[
        ,
        .(
            handle = first_non_na(handle),
            FBt_cdspos = first_non_na(FBt_cdspos),
            gene = first_non_na(gene),
            ref_codon = first_non_na(ref_codon),
            alt_codon = first_non_na(alt_codon),
            edit_pos = first_non_na(edit_pos),
            NS_class = first_non_na(NS_class)
        ),
        by = site_id
    ]

    meta <- meta[
        match(sites, site_id)
    ]

    list(
        reps = reps,
        sites = sites,
        meta = meta,
        RNA_ref = make_matrix("RNA_ref"),
        RNA_alt = make_matrix("RNA_alt"),
        Ribo_ref = make_matrix("Ribo_ref"),
        Ribo_alt = make_matrix("Ribo_alt")
    )
}

# ============================================================
# 10. Common test universe
# ============================================================

build_common_filter <- function(m) {
    nrep <- length(m$reps)

    RNA_total <- m$RNA_ref + m$RNA_alt
    Ribo_total <- m$Ribo_ref + m$Ribo_alt

    RNA_sum <- rowSums(RNA_total)
    Ribo_sum <- rowSums(Ribo_total)

    pooled_Lm <- fifelse(
        RNA_sum > 0,
        rowSums(m$RNA_alt) / RNA_sum,
        NA_real_
    )

    # Number of replicates with at least non-zero RNA and Ribo information.
    paired_nonzero_reps <- rowSums(
        RNA_total > 0 &
            Ribo_total > 0
    )

    # Number meeting the old coverage threshold in the same replicate.
    paired_coverage_reps <- rowSums(
        RNA_total >= opt$min_rna_per_rep_equiv &
            Ribo_total >= opt$min_ribo_per_rep_equiv
    )

    # Replicate support, rather than average coverage across every replicate,
    # defines the common universe.  The previous RNA_sum/Ribo_sum >=
    # threshold * nrep rule was unnecessarily stringent for large pooled
    # groups (for example, Dmel has eight replicates and therefore required
    # 40 RNA and 40 Ribo reads even when only two replicates were required).
    keep <- (
        nrep >= opt$min_replicates &
            paired_coverage_reps >= opt$min_replicates &
            is.finite(pooled_Lm) &
            pooled_Lm >= opt$min_lm
    )

    keep[is.na(keep)] <- FALSE

    data.table(
        site_id = m$sites,
        test_universe = keep,
        n_group_replicates = nrep,
        RNA_sum = RNA_sum,
        Ribo_sum = Ribo_sum,
        pooled_Lm = pooled_Lm,
        paired_nonzero_reps = paired_nonzero_reps,
        paired_coverage_reps = paired_coverage_reps
    )
}

subset_matrices <- function(m, keep) {
    idx <- which(keep)

    list(
        reps = m$reps,
        sites = m$sites[idx],
        meta = m$meta[idx],
        RNA_ref = m$RNA_ref[idx, , drop = FALSE],
        RNA_alt = m$RNA_alt[idx, , drop = FALSE],
        Ribo_ref = m$Ribo_ref[idx, , drop = FALSE],
        Ribo_alt = m$Ribo_alt[idx, , drop = FALSE]
    )
}

# ============================================================
# 11. Primary analysis: edgeR robust QL
# ============================================================

run_edgeR_QL <- function(
  m,
  analysis_group,
  analysis_level,
  position
) {
    nrep <- length(m$reps)

    base <- cbind(
        m$meta,
        data.table(
            analysis_group = analysis_group,
            analysis_level = analysis_level,
            position = position,
            method = "edgeR_QL",
            n_replicates = nrep,
            testable = FALSE,
            model_converged = FALSE,
            effect_log2OR = NA_real_,
            raw_p = NA_real_,
            q_value = NA_real_,
            call = "not_testable"
        )
    )

    if (
        nrep < opt$min_replicates ||
            !length(m$sites)
    ) {
        return(base)
    }

    # ----------------------------------------------------------
    # For each biological replicate there are four observations:
    # RNA_A, RNA_I, Ribo_A, Ribo_I.
    #
    # The model includes:
    # 1. one fixed effect for every replicate x assay library
    #    -> controls library depth / assay baseline
    # 2. one replicate-specific edited-allele baseline
    #    -> controls replicate-to-replicate RNA editing variation
    # 3. a common Ribo-specific edited-allele coefficient
    #    -> target effect:
    #       log[(Ribo_I/Ribo_A)/(RNA_I/RNA_A)]
    # ----------------------------------------------------------

    obs_list <- vector(
        "list",
        nrep
    )

    count_list <- vector(
        "list",
        nrep
    )

    for (j in seq_along(m$reps)) {
        rep_i <- m$reps[[j]]

        obs_list[[j]] <- data.table(
            replicate = rep_i,
            assay = c(
                "RNA",
                "RNA",
                "Ribo",
                "Ribo"
            ),
            allele = c(
                "A",
                "I",
                "A",
                "I"
            ),
            library = c(
                paste0(rep_i, "__RNA"),
                paste0(rep_i, "__RNA"),
                paste0(rep_i, "__Ribo"),
                paste0(rep_i, "__Ribo")
            )
        )

        count_list[[j]] <- cbind(
            m$RNA_ref[, j],
            m$RNA_alt[, j],
            m$Ribo_ref[, j],
            m$Ribo_alt[, j]
        )
    }

    obs <- rbindlist(
        obs_list
    )

    counts <- do.call(
        cbind,
        count_list
    )

    rownames(counts) <- m$sites

    colnames(counts) <- paste(
        obs$replicate,
        obs$assay,
        obs$allele,
        sep = "__"
    )

    obs[, replicate := factor(
        replicate,
        levels = m$reps
    )]

    obs[, assay := factor(
        assay,
        levels = c("RNA", "Ribo")
    )]

    obs[, allele := factor(
        allele,
        levels = c("A", "I")
    )]

    obs[, library := factor(
        library,
        levels = unique(library)
    )]

    # One intercept for each independent RNA/Ribo library.
    library_matrix <- model.matrix(
        ~ 0 + library,
        data = obs
    )

    # Replicate-specific baseline I/A ratio.
    rep_matrix <- model.matrix(
        ~ 0 + replicate,
        data = obs
    )

    edited_indicator <- as.numeric(
        obs$allele == "I"
    )

    replicate_edit_matrix <- (
        rep_matrix *
            edited_indicator
    )

    colnames(replicate_edit_matrix) <- paste0(
        "baseline_edit_",
        levels(obs$replicate)
    )

    # Target term.
    ribo_edited <- as.numeric(
        obs$assay == "Ribo" &
            obs$allele == "I"
    )

    design <- cbind(
        library_matrix,
        replicate_edit_matrix,
        ribo_edited = ribo_edited
    )

    design_rank <- qr(design)$rank
    residual_df <- nrow(design) - design_rank

    if (
        design_rank != ncol(design) ||
            residual_df < 1L
    ) {
        warning(
            analysis_group,
            " ",
            position,
            ": edgeR design has insufficient residual df."
        )
        return(base)
    }

    # No external library-size normalization is necessary for the
    # allele-ratio coefficient because each RNA/Ribo library has its
    # own fixed effect and A/I share the same library exposure.
    dge <- edgeR::DGEList(
        counts = counts,
        lib.size = rep(
            1,
            ncol(counts)
        )
    )

    dge$samples$norm.factors <- 1

    dge$offset <- matrix(
        0,
        nrow = nrow(counts),
        ncol = ncol(counts)
    )

    ans <- tryCatch(
        {
            dge <- edgeR::estimateDisp(
                dge,
                design = design,
                robust = TRUE
            )

            fit <- edgeR::glmQLFit(
                dge,
                design = design,
                robust = TRUE
            )

            test <- edgeR::glmQLFTest(
                fit,
                coef = "ribo_edited"
            )

            stats <- as.data.table(
                edgeR::topTags(
                    test,
                    n = Inf,
                    sort.by = "none"
                )$table
            )

            data.table(
                effect_log2OR = stats$logFC,
                raw_p = stats$PValue,
                q_value = p.adjust(
                    stats$PValue,
                    method = "BH"
                ),
                model_converged = TRUE,
                edgeR_common_dispersion = dge$common.dispersion
            )
        },
        error = function(e) {
            warning(
                analysis_group,
                " ",
                position,
                ": edgeR failed: ",
                conditionMessage(e)
            )
            NULL
        }
    )

    if (is.null(ans)) {
        return(base)
    }

    base[
        ,
        `:=`(
            testable = TRUE,
            model_converged = ans$model_converged,
            effect_log2OR = ans$effect_log2OR,
            raw_p = ans$raw_p,
            q_value = ans$q_value,
            edgeR_common_dispersion = ans$edgeR_common_dispersion
        )
    ]

    base[, call := fcase(
        q_value <= opt$fdr &
            effect_log2OR > 0,
        "I>A",
        q_value <= opt$fdr &
            effect_log2OR < 0,
        "I<A",
        default = "not_significant"
    )]

    base[]
}

# ============================================================
# 12. Sensitivity analysis: beta-binomial GLM
# ============================================================

fit_beta_binomial_site <- function(
  i,
  m
) {
    RNA_total <- (
        m$RNA_ref[i, ] +
            m$RNA_alt[i, ]
    )

    Ribo_total <- (
        m$Ribo_ref[i, ] +
            m$Ribo_alt[i, ]
    )

    # Beta-binomial requires defined denominators.
    usable <- (
        RNA_total > 0 &
            Ribo_total > 0
    )

    if (
        sum(usable) <
            opt$min_replicates
    ) {
        return(
            list(
                n_reps_used = sum(usable),
                effect_log2OR = NA_real_,
                raw_p = NA_real_,
                converged = FALSE
            )
        )
    }

    idx <- which(usable)

    d <- rbind(
        data.table(
            replicate = m$reps[idx],
            assay = "RNA",
            edited = as.numeric(
                m$RNA_alt[i, idx]
            ),
            unedited = as.numeric(
                m$RNA_ref[i, idx]
            )
        ),
        data.table(
            replicate = m$reps[idx],
            assay = "Ribo",
            edited = as.numeric(
                m$Ribo_alt[i, idx]
            ),
            unedited = as.numeric(
                m$Ribo_ref[i, idx]
            )
        )
    )

    d[, replicate := factor(
        replicate,
        levels = m$reps[idx]
    )]

    d[, assay := factor(
        assay,
        levels = c("RNA", "Ribo")
    )]

    ans <- tryCatch(
        {
            # Null: each replicate can have its own baseline editing fraction.
            fit0 <- glmmTMB::glmmTMB(
                cbind(
                    edited,
                    unedited
                ) ~ replicate,
                family = glmmTMB::betabinomial(
                    link = "logit"
                ),
                data = d
            )

            # Full: common RNA -> Ribo effect added.
            fit1 <- glmmTMB::glmmTMB(
                cbind(
                    edited,
                    unedited
                ) ~ replicate + assay,
                family = glmmTMB::betabinomial(
                    link = "logit"
                ),
                data = d
            )

            lrt <- suppressMessages(
                anova(
                    fit0,
                    fit1
                )
            )

            p_value <- as.numeric(
                lrt$`Pr(>Chisq)`[2L]
            )

            coefs <- summary(
                fit1
            )$coefficients$cond

            estimate <- if (
                "assayRibo" %in%
                    rownames(coefs)
            ) {
                coefs[
                    "assayRibo",
                    "Estimate"
                ]
            } else {
                NA_real_
            }

            converged <- (
                identical(
                    fit0$fit$convergence,
                    0L
                ) &&
                    identical(
                        fit1$fit$convergence,
                        0L
                    ) &&
                    isTRUE(
                        fit0$sdr$pdHess
                    ) &&
                    isTRUE(
                        fit1$sdr$pdHess
                    )
            )

            list(
                n_reps_used = sum(usable),

                # logit coefficient is natural-log OR.
                # Convert to log2 OR.
                effect_log2OR = estimate / log(2),
                raw_p = p_value,
                converged = converged
            )
        },
        error = function(e) {
            list(
                n_reps_used = sum(usable),
                effect_log2OR = NA_real_,
                raw_p = NA_real_,
                converged = FALSE
            )
        }
    )

    ans
}

run_beta_binomial <- function(
  m,
  analysis_group,
  analysis_level,
  position
) {
    nrep <- length(
        m$reps
    )

    base <- cbind(
        m$meta,
        data.table(
            analysis_group = analysis_group,
            analysis_level = analysis_level,
            position = position,
            method = "BetaBin_GLM",
            n_replicates = nrep,
            n_reps_used = NA_integer_,
            testable = FALSE,
            model_converged = FALSE,
            effect_log2OR = NA_real_,
            raw_p = NA_real_,
            q_value = NA_real_,
            call = "not_testable"
        )
    )

    if (
        nrep < opt$min_replicates ||
            !length(m$sites)
    ) {
        return(base)
    }

    indexes <- seq_along(
        m$sites
    )

    if (
        .Platform$OS.type == "unix" &&
            opt$threads > 1L
    ) {
        fits <- parallel::mclapply(
            indexes,
            fit_beta_binomial_site,
            m = m,
            mc.cores = opt$threads,
            mc.preschedule = TRUE
        )
    } else {
        fits <- lapply(
            indexes,
            fit_beta_binomial_site,
            m = m
        )
    }

    base[, n_reps_used := vapply(
        fits,
        function(z) z$n_reps_used,
        integer(1L)
    )]

    base[, effect_log2OR := vapply(
        fits,
        function(z) z$effect_log2OR,
        numeric(1L)
    )]

    base[, raw_p := vapply(
        fits,
        function(z) z$raw_p,
        numeric(1L)
    )]

    base[, model_converged := vapply(
        fits,
        function(z) z$converged,
        logical(1L)
    )]

    base[, testable := (
        model_converged &
            is.finite(raw_p) &
            is.finite(effect_log2OR)
    )]

    base[, q_value := NA_real_]

    base[
        testable == TRUE,
        q_value := p.adjust(
            raw_p,
            method = "BH"
        )
    ]

    base[, call := fcase(
        testable &
            q_value <= opt$fdr &
            effect_log2OR > 0,
        "I>A",
        testable &
            q_value <= opt$fdr &
            effect_log2OR < 0,
        "I<A",
        testable,
        "not_significant",
        default = "not_testable"
    )]

    base[]
}

# ============================================================
# 13. Original/reference analysis: Fisher
# ============================================================

# Directional Fisher reference: all eligible replicates contribute to BOTH tails.
# BH family = every testable site x 2 directions, within analysis_group x A/P.
# Counts are rounded only for Fisher, exactly as in the previous implementation.
fisher_replicate_audit <- list()
run_fisher_analysis <- function(m, analysis_group, analysis_level, position) {
    nrep <- length(m$reps)
    stats_list <- vector('list', length(m$sites))
    audit_list <- vector('list', length(m$sites))
    for (i in seq_along(m$sites)) {
        rna_ref <- m$RNA_ref[i, ]; rna_alt <- m$RNA_alt[i, ]
        ribo_ref <- m$Ribo_ref[i, ]; ribo_alt <- m$Ribo_alt[i, ]
        rna_total <- rna_ref + rna_alt
        Lm <- ifelse(rna_total > 0, rna_alt / rna_total, NA_real_)
        # Site eligibility was already decided by edgeR's common filter.
        # Consider EVERY group replicate; do not reapply per-replicate depth or Lm cuts.
        idx <- seq_along(m$reps)
        p_gt <- p_lt <- effects <- numeric()
        directions <- test_directions <- numeric()
        rep_rows <- list()
        for (j in idx) {
            # OR = (Ribo_I/Ribo_A)/(RNA_I/RNA_A); rows=assay, columns=I,A.
            tab <- round(rbind(c(ribo_alt[j], ribo_ref[j]), c(rna_alt[j], rna_ref[j])))
            reason <- if (any(rowSums(tab) == 0)) 'zero_assay_total_after_rounding' else if (any(colSums(tab) == 0)) 'zero_allele_margin_after_rounding' else 'informative'
            informative <- reason == 'informative'
            # Degenerate conditional tables contain no OR information. Keep them in
            # the audit, but not in the Fisher chi-square degrees of freedom.
            pv <- if (informative) c(
                fisher.test(tab, alternative = 'greater')$p.value,
                fisher.test(tab, alternative = 'less')$p.value
            ) else c(NA_real_, NA_real_)
            if (informative && !all(is.finite(pv))) stop('Non-finite one-sided Fisher P')
            # Direction counts describe unrounded allele ratios, without pseudocount bias.
            dr <- sign(ribo_alt[j] * rna_ref[j] - ribo_ref[j] * rna_alt[j])
            dt <- sign(tab[1, 1] * tab[2, 2] - tab[1, 2] * tab[2, 1])
            eff <- log2(((ribo_alt[j] + opt$fisher_pc) / (ribo_ref[j] + opt$fisher_pc)) /
                ((rna_alt[j] + opt$fisher_pc) / (rna_ref[j] + opt$fisher_pc)))
            used <- informative
            if (used) {
                p_gt <- c(p_gt, pv[1]); p_lt <- c(p_lt, pv[2])
                effects <- c(effects, eff); directions <- c(directions, dr)
                test_directions <- c(test_directions, dt)
            }
            rep_rows[[length(rep_rows) + 1L]] <- data.table(
                analysis_group = analysis_group, position = position, site_id = m$sites[i],
                sample_id = m$reps[j], inclusion_reason = reason, RNA_ref = rna_ref[j], RNA_alt = rna_alt[j],
                Ribo_ref = ribo_ref[j], Ribo_alt = ribo_alt[j],
                RNA_ref_Fisher = tab[2, 2], RNA_alt_Fisher = tab[2, 1],
                Ribo_ref_Fisher = tab[1, 2], Ribo_alt_Fisher = tab[1, 1],
                p_I_gt_A = pv[1], p_I_lt_A = pv[2], used_in_both_directions = used,
                direction_unrounded = dr, direction_Fisher_counts = dt,
                descriptive_log2OR_pc = eff
            )
        }
        audit_list[[i]] <- rbindlist(rep_rows)
        nused <- length(p_gt); ok <- TRUE
        # The >=2 repeat requirement is SITE coverage support, identical to edgeR.
        # With no informative integer table, retain this site in BH with P=1.
        npos <- sum(directions > 0); nneg <- sum(directions < 0); nzero <- sum(directions == 0)
        stats_list[[i]] <- data.table(
            n_reps_eligible = length(idx), n_reps_used = nused,
            n_reps_uninformative = length(idx) - nused, fisher_has_information = nused > 0,
            n_positive_reps = npos, n_negative_reps = nneg, n_neutral_reps = nzero,
            n_positive_reps_Fisher_counts = sum(test_directions > 0),
            n_negative_reps_Fisher_counts = sum(test_directions < 0),
            n_neutral_reps_Fisher_counts = sum(test_directions == 0),
            majority_direction = if (npos > nneg) 'I>A' else if (nneg > npos) 'I<A' else 'tie',
            majority_fraction = if (nused > 0) max(npos, nneg) / nused else NA_real_,
            effect_log2OR = if (nused > 0) median(effects) else NA_real_,
            raw_p_I_gt_A = if (nused > 0) fisher_combine_p(p_gt) else 1,
            raw_p_I_lt_A = if (nused > 0) fisher_combine_p(p_lt) else 1,
            testable = ok, model_converged = nused > 0
        )
    }
    fisher_replicate_audit[[paste(analysis_group, position, sep = '__')]] <<- rbindlist(audit_list, fill = TRUE)
    out <- cbind(
        m$meta, data.table(
            analysis_group = analysis_group,
            analysis_level = analysis_level, position = position, method = 'Fisher', n_replicates = nrep
        ),
        rbindlist(stats_list)
    )
    out[, `:=`(q_I_gt_A = NA_real_, q_I_lt_A = NA_real_)]
    eligible_sites <- which(out$testable)
    if (length(eligible_sites)) {
        n <- length(eligible_sites)
        q <- p.adjust(c(out$raw_p_I_gt_A[eligible_sites], out$raw_p_I_lt_A[eligible_sites]), method = 'BH')
        out[eligible_sites, `:=`(q_I_gt_A = q[seq_len(n)], q_I_lt_A = q[n + seq_len(n)])]
    }
    out[, `:=`(
        significant_I_gt_A = testable & !is.na(q_I_gt_A) & q_I_gt_A <= opt$fdr,
        significant_I_lt_A = testable & !is.na(q_I_lt_A) & q_I_lt_A <= opt$fdr,
        # Compatibility summaries only; decisions use the TWO directional q columns.
        raw_p = pmin(1, 2 * pmin(raw_p_I_gt_A, raw_p_I_lt_A)),
        q_value = pmin(q_I_gt_A, q_I_lt_A),
        direction_counts_position = position,
        fisher_BH_family = 'all_testable_sites_x_2_directions_within_group_position',
        fisher_q_semantics = 'minimum_directional_BH_q; see q_I_gt_A and q_I_lt_A'
    )]
    out[, call := fcase(
        significant_I_gt_A & significant_I_lt_A, 'direction_conflict',
        significant_I_gt_A, 'I>A', significant_I_lt_A, 'I<A',
        testable, 'not_significant',
        default = 'not_testable'
    )]
    out[, direction_consistency := fifelse(
        call == 'I>A', n_positive_reps / pmax(n_reps_used, 1),
        fifelse(call == 'I<A', n_negative_reps / pmax(n_reps_used, 1), NA_real_)
    )]
    out[]
}

# AP remains an A/P significance union, NOT a separately FDR-controlled test.
# Fisher conflicts (within a position or between A/P) get one exclusive row.
make_fisher_AP_union <- function(x) {
    x[,
        {
            z <- copy(.SD)
            gt <- any(z$significant_I_gt_A, na.rm = TRUE)
            lt <- any(z$significant_I_lt_A, na.rm = TRUE)
            z[, rank_q := fifelse(is.finite(q_value), q_value, Inf)]
            setorder(z, rank_q, position)
            rep <- copy(z[1L]); rep[, rank_q := NULL]
            minimum <- function(v) if (any(is.finite(v))) min(v[is.finite(v)]) else NA_real_
            rep[, `:=`(
                position = 'AP', testable = any(z$testable),
                model_converged = any(z$model_converged),
                call = if (gt && lt) 'direction_conflict' else if (gt) 'I>A' else if (lt) 'I<A' else if (any(z$testable)) 'not_significant' else 'not_testable',
                significant_I_gt_A = gt, significant_I_lt_A = lt,
                raw_p_I_gt_A = minimum(z$raw_p_I_gt_A), raw_p_I_lt_A = minimum(z$raw_p_I_lt_A),
                q_I_gt_A = minimum(z$q_I_gt_A), q_I_lt_A = minimum(z$q_I_lt_A),
                raw_p = NA_real_, q_value = minimum(z$q_value),
                AP_definition = 'A_significant_OR_P_significant; conflicts exclusive',
                AP_significant_modes = paste(z$position[z$significant_I_gt_A | z$significant_I_lt_A], collapse = ';'),
                AP_I_gt_A_modes = paste(z$position[z$significant_I_gt_A], collapse = ';'),
                AP_I_lt_A_modes = paste(z$position[z$significant_I_lt_A], collapse = ';'),
                AP_conflict = gt && lt,
                fisher_q_semantics = 'AP min of A/P directional q values; NOT AP-adjusted FDR'
            )]
            rep[, direction_consistency := fifelse(
                call == 'I>A', n_positive_reps / pmax(n_reps_used, 1),
                fifelse(call == 'I<A', n_negative_reps / pmax(n_reps_used, 1), NA_real_)
            )]
            # n_reps and majority fields refer to direction_counts_position, the representative A or P.
            rep[]
        },
        by = .(method, analysis_group, analysis_level, site_id)
    ]
}

# ============================================================
# 14. Run the configured groups
# ============================================================

group_info <- unique(
    count_long[
        ,
        .(
            analysis_group,
            analysis_level,
            position
        )
    ]
)

group_info[, analysis_group_order_num := match(
    analysis_group,
    analysis_group_order
)]

setorder(
    group_info,
    analysis_group_order_num,
    position
)

group_info[, analysis_group_order_num := NULL]

all_results_list <- list()
universe_list <- list()

out_index <- 0L

for (
    i in seq_len(
        nrow(group_info)
    )
) {
    group_i <- group_info$analysis_group[[i]]
    level_i <- group_info$analysis_level[[i]]
    pos_i <- group_info$position[[i]]

    message(
        "[ANALYSIS] ",
        group_i,
        " | ",
        pos_i
    )

    x <- count_long[
        analysis_group == group_i &
            position == pos_i
    ]

    m0 <- build_group_matrices(
        x
    )

    filter_dt <- build_common_filter(
        m0
    )

    universe_list[[
        length(universe_list) + 1L
    ]] <- cbind(
        data.table(
            analysis_group = group_i,
            analysis_level = level_i,
            position = pos_i
        ),
        filter_dt
    )

    if (
        length(m0$reps) <
            opt$min_replicates
    ) {
        warning(
            group_i,
            " has only ",
            length(m0$reps),
            " replicates; skipped."
        )
        next
    }

    m <- subset_matrices(
        m0,
        filter_dt$test_universe
    )

    if (!length(m$sites)) {
        warning(
            group_i,
            " ",
            pos_i,
            ": no sites passed common filtering."
        )
        next
    }

    edgeR_result <- run_edgeR_QL(
        m,
        group_i,
        level_i,
        pos_i
    )

    # Requested workflow: edgeR QL primary and Fisher reference only.
    beta_result <- NULL

    fisher_result <- run_fisher_analysis(
        m,
        group_i,
        level_i,
        pos_i
    )

    out_index <- out_index + 1L

    all_results_list[[out_index]] <- rbindlist(
        list(
            edgeR_result,
            beta_result,
            fisher_result
        ),
        use.names = TRUE,
        fill = TRUE
    )
}

if (!length(all_results_list)) {
    stop("No analysis group produced results.")
}

individual_AP_results <- rbindlist(
    all_results_list,
    use.names = TRUE,
    fill = TRUE
)

test_universe <- rbindlist(
    universe_list,
    use.names = TRUE,
    fill = TRUE
)

fwrite(
    test_universe,
    file.path(
        opt$output_dir,
        "common_test_universe.tsv.gz"
    ),
    sep = "\t"
)

# ============================================================
# 15. AP = A significant OR P significant
# ============================================================

make_AP_union <- function(x) {
    source_rows <- x[
        position %in% c("A", "P")
    ]

    if (!nrow(source_rows)) {
        return(data.table())
    }

    group_cols <- c(
        "method",
        "analysis_group",
        "analysis_level",
        "site_id"
    )

    out <- source_rows[
        ,
        {
            z <- copy(.SD)

            significant_directions <- unique(
                z[
                    call %in% c(
                        "I>A",
                        "I<A"
                    ),
                    call
                ]
            )

            significant_directions <- c(
                intersect(
                    "I>A",
                    significant_directions
                ),
                intersect(
                    "I<A",
                    significant_directions
                )
            )

            if (
                !length(
                    significant_directions
                )
            ) {
                rank_p <- fifelse(
                    is.finite(z$q_value),
                    z$q_value,
                    fifelse(
                        is.finite(z$raw_p),
                        z$raw_p,
                        Inf
                    )
                )

                z[, temp_rank := rank_p]

                setorder(
                    z,
                    temp_rank,
                    position
                )

                representative <- copy(
                    z[1L]
                )

                representative[
                    ,
                    `:=`(
                        position = "AP",
                        call = ifelse(
                            any(z$testable),
                            "not_significant",
                            "not_testable"
                        ),
                        AP_definition = "A_significant_OR_P_significant",
                        AP_significant_modes = NA_character_,
                        AP_conflict = FALSE
                    )
                ]

                representative[
                    ,
                    temp_rank := NULL
                ]

                representative[]
            } else {
                pieces <- lapply(
                    significant_directions,
                    function(direction_i) {
                        z2 <- copy(
                            z[
                                call == direction_i
                            ]
                        )

                        z2[, temp_rank := fifelse(
                            is.finite(q_value),
                            q_value,
                            fifelse(
                                is.finite(raw_p),
                                raw_p,
                                Inf
                            )
                        )]

                        setorder(
                            z2,
                            temp_rank,
                            position
                        )

                        representative <- copy(
                            z2[1L]
                        )

                        representative[
                            ,
                            `:=`(
                                position = "AP",
                                call = direction_i,
                                AP_definition = "A_significant_OR_P_significant",
                                AP_significant_modes = paste(
                                    sort(
                                        unique(
                                            z2$position
                                        )
                                    ),
                                    collapse = ";"
                                ),
                                AP_conflict = all(
                                    c(
                                        "I>A",
                                        "I<A"
                                    ) %in%
                                        significant_directions
                                )
                            )
                        ]

                        representative[
                            ,
                            temp_rank := NULL
                        ]

                        representative[]
                    }
                )

                rbindlist(
                    pieces,
                    use.names = TRUE,
                    fill = TRUE
                )
            }
        },
        by = group_cols
    ]

    out[]
}

AP_results <- rbindlist(list(
    make_AP_union(individual_AP_results[method != "Fisher"]),
    make_fisher_AP_union(individual_AP_results[method == "Fisher"])
), use.names = TRUE, fill = TRUE)

all_results <- rbindlist(
    list(
        individual_AP_results,
        AP_results
    ),
    use.names = TRUE,
    fill = TRUE
)

# ============================================================
# 16. Save site-level results
# ============================================================

fwrite(
    all_results[
        method == "edgeR_QL"
    ],
    file.path(
        opt$output_dir,
        "site_results.edgeR_QL.tsv.gz"
    ),
    sep = "\t"
)


fwrite(
    all_results[
        method == "Fisher"
    ],
    file.path(
        opt$output_dir,
        "site_results.Fisher.tsv.gz"
    ),
    sep = "\t"
)

fwrite(
    all_results,
    file.path(
        opt$output_dir,
        "site_results.ALL_methods.tsv.gz"
    ),
    sep = "\t"
)

# ============================================================
# 17. Table 1-style summary
# ============================================================

summary_long <- all_results[
    testable == TRUE,
    {
        total_ids <- unique(
            site_id
        )

        Igt_ids <- unique(
            site_id[
                call == "I>A"
            ]
        )

        Ilt_ids <- unique(
            site_id[
                call == "I<A"
            ]
        )

        .(
            Total_sites = length(
                total_ids
            ),
            Total_N = uniqueN(
                site_id[
                    NS_class == "N"
                ]
            ),
            Total_S = uniqueN(
                site_id[
                    NS_class == "S"
                ]
            ),
            I_gt_A = length(
                Igt_ids
            ),
            I_gt_A_N = uniqueN(
                site_id[
                    call == "I>A" &
                        NS_class == "N"
                ]
            ),
            I_gt_A_S = uniqueN(
                site_id[
                    call == "I>A" &
                        NS_class == "S"
                ]
            ),
            I_lt_A = length(
                Ilt_ids
            ),
            I_lt_A_N = uniqueN(
                site_id[
                    call == "I<A" &
                        NS_class == "N"
                ]
            ),
            Direction_conflict = uniqueN(site_id[call == "direction_conflict"]),
            I_lt_A_S = uniqueN(
                site_id[
                    call == "I<A" &
                        NS_class == "S"
                ]
            )
        )
    },
    by = .(
        method,
        analysis_group,
        analysis_level,
        position
    )
]

summary_long[, `:=`(
    Total_sites_N_S = format_NS(
        Total_sites,
        Total_N,
        Total_S
    ),
    I_gt_A_N_S = format_NS(
        I_gt_A,
        I_gt_A_N,
        I_gt_A_S
    ),
    I_lt_A_N_S = format_NS(
        I_lt_A,
        I_lt_A_N,
        I_lt_A_S
    )
)]

summary_long[, order_num := match(
    analysis_group,
    analysis_group_order
)]

setorder(
    summary_long,
    method,
    position,
    order_num
)

summary_long[, order_num := NULL]

fwrite(
    summary_long,
    file.path(
        opt$output_dir,
        "Table1_summary_long.tsv"
    ),
    sep = "\t"
)

summary_wide <- dcast(
    summary_long,
    method + analysis_group + analysis_level ~ position,
    value.var = c(
        "Total_sites_N_S",
        "I_gt_A_N_S",
        "I_lt_A_N_S"
    )
)

fwrite(
    summary_wide,
    file.path(
        opt$output_dir,
        "Table1_summary_wide.tsv"
    ),
    sep = "\t"
)

# Final manuscript-style Table 1 using the primary edgeR QL analysis.
final_sample_order <- data.table(
    analysis_group = c(
        "Dmel_F", "Dmel_M",
        "Dsim_F", "Dsim_M",
        "Dsec_F", "Dsec_M",
        "Dere_F", "Dere_M",
        "Dmel", "Dsim", "Dsec", "Dere"
    ),
    sample = c(
        "mel_F", "mel_M",
        "sim_F", "sim_M",
        "sec_F", "sec_M",
        "ere_F", "ere_M",
        "mel", "sim", "sec", "ere"
    ),
    final_order = seq_len(12L)
)

final_table1 <- dcast(
    summary_long[
        method == "edgeR_QL" &
            position %in% c("A", "AP")
    ],
    analysis_group ~ position,
    value.var = c(
        "Total_sites_N_S",
        "I_gt_A_N_S",
        "I_lt_A"
    )
)

final_table1 <- final_sample_order[
    final_table1,
    on = "analysis_group"
]

setorder(final_table1, final_order)

final_table1 <- final_table1[, .(
    sample,
    A_Total_sites_N_S = Total_sites_N_S_A,
    A_I_gt_A_N_S = I_gt_A_N_S_A,
    A_I_lt_A = I_lt_A_A,
    AP_Total_sites_N_S = Total_sites_N_S_AP,
    AP_I_gt_A_N_S = I_gt_A_N_S_AP,
    AP_I_lt_A = I_lt_A_AP
)]

fwrite(
    final_table1,
    file.path(
        opt$output_dir,
        "Table1_edgeR_final.tsv"
    ),
    sep = "\t"
)

# A display-ready version with a two-row column header.
final_table1_display_file <- file.path(
    opt$output_dir,
    "Table1_edgeR_final_two_row_header.tsv"
)

writeLines(
    c(
        paste(
            c(
                "sample",
                "A-site vs. mRNA", "", "",
                "A/P-site vs. mRNA", "", ""
            ),
            collapse = "\t"
        ),
        paste(
            c(
                "",
                "Total sites (N/S)", "I>A (N/S)", "I<A",
                "Total sites (N/S)", "I>A (N/S)", "I<A"
            ),
            collapse = "\t"
        )
    ),
    final_table1_display_file
)

fwrite(
    final_table1,
    final_table1_display_file,
    sep = "\t",
    col.names = FALSE,
    append = TRUE
)

# ============================================================
# 18. Two-method overlap for final figures: edgeR QL versus Fisher
# ============================================================

two_method_overlap_one <- function(
  x,
  group_i,
  position_i,
  direction_i
) {
    d <- x[
        analysis_group == group_i &
            position == position_i &
            call == direction_i
    ]

    edgeR_set <- unique(
        d[
            method == "edgeR_QL",
            site_id
        ]
    )

    fisher_set <- unique(
        d[
            method == "Fisher",
            site_id
        ]
    )

    edgeR_fisher <- intersect(
        edgeR_set,
        fisher_set
    )

    edgeR_only <- setdiff(
        edgeR_set,
        fisher_set
    )

    fisher_only <- setdiff(
        fisher_set,
        edgeR_set
    )

    data.table(
        analysis_group = group_i,
        position = position_i,
        direction = direction_i,
        Fisher_only = length(
            fisher_only
        ),
        edgeR_only = length(
            edgeR_only
        ),
        Fisher_edgeR = length(
            edgeR_fisher
        ),
        Fisher_total = length(
            fisher_set
        ),
        edgeR_total = length(
            edgeR_set
        ),
        Union_total = length(
            union(
                edgeR_set,
                fisher_set
            )
        )
    )
}

make_overlap_table <- function(
  position_i
) {
    groups_present <- analysis_group_order[
        analysis_group_order %in%
            unique(
                all_results$analysis_group
            )
    ]

    rbindlist(
        lapply(
            groups_present,
            function(group_i) {
                rbindlist(
                    list(
                        two_method_overlap_one(
                            all_results,
                            group_i,
                            position_i,
                            "I>A"
                        ),
                        two_method_overlap_one(
                            all_results,
                            group_i,
                            position_i,
                            "I<A"
                        )
                    ),
                    use.names = TRUE,
                    fill = TRUE
                )
            }
        ),
        use.names = TRUE,
        fill = TRUE
    )
}

# ============================================================
# 19. Detailed 3-class overlap plot
# ============================================================

plot_overlap_detailed <- function(
  overlap_dt,
  title_text,
  output_stem
) {
    class_cols <- c(
        "Fisher_only",
        "edgeR_only",
        "Fisher_edgeR"
    )

    pdat <- melt(
        overlap_dt,
        id.vars = c(
            "analysis_group",
            "position",
            "direction"
        ),
        measure.vars = class_cols,
        variable.name = "overlap_class",
        value.name = "N"
    )

    pdat[, analysis_group := factor(
        analysis_group,
        levels = analysis_group_order
    )]

    pdat[, direction := factor(
        direction,
        levels = c(
            "I>A",
            "I<A"
        )
    )]

    pdat[, overlap_class := factor(
        overlap_class,
        levels = class_cols
    )]

    pdat[, block_label := fifelse(
        N > 0,
        as.character(N),
        ""
    )]

    colors <- c(
        Fisher_only = "#E777A7",
        edgeR_only = "#6AACDD",
        Fisher_edgeR = "#D9D9D9"
    )

    labels <- c(
        Fisher_only = "Fisher only",
        edgeR_only = "edgeR only",
        Fisher_edgeR = "Shared by Fisher and edgeR"
    )

    p <- ggplot(
        pdat,
        aes(
            x = analysis_group,
            y = N,
            fill = overlap_class
        )
    ) +
        geom_col(
            width = 0.72
        ) +
        geom_text(
            aes(
                label = block_label
            ),
            position = position_stack(
                vjust = 0.5
            ),
            color = "black",
            size = 2.5
        ) +
        facet_grid(
            direction ~ .,
            scales = "free_y",
            space = "free_y"
        ) +
        scale_fill_manual(
            values = colors,
            labels = labels,
            drop = FALSE
        ) +
        scale_x_discrete(
            labels = analysis_group_labels
        ) +
        labs(
            x = NULL,
            y = "Sites number",
            fill = NULL,
            title = title_text
        ) +
        theme_classic(
            base_size = 10
        ) +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                face = "bold",
                size = 11
            ),
            strip.background = element_blank(),
            strip.text.y = element_text(
                face = "bold",
                size = 10
            ),
            axis.text.x = element_text(
                angle = 35,
                hjust = 1,
                size = 8
            ),
            axis.text.y = element_text(
                size = 8
            ),
            axis.title.y = element_text(
                size = 9
            ),
            legend.position = "right",
            legend.text = element_text(
                size = 8
            ),
            panel.spacing.y = grid::unit(
                0.7,
                "lines"
            )
        )

    # Visual separators:
    # Dmel / Dsim / Dsec / Dere sex-specific pairs,
    # then sex-combined species results.
    for (
        xint in c(
            2.5,
            4.5,
            6.5,
            8.5
        )
    ) {
        p <- p +
            geom_vline(
                xintercept = xint,
                linetype = "dashed",
                linewidth = 0.25,
                color = "grey65"
            )
    }

    ggsave(
        file.path(
            opt$output_dir,
            paste0(
                output_stem,
                ".pdf"
            )
        ),
        p,
        width = 12,
        height = 6.5
    )

    ggsave(
        file.path(
            opt$output_dir,
            paste0(
                output_stem,
                ".png"
            )
        ),
        p,
        width = 12,
        height = 6.5,
        dpi = 300
    )

    p
}

# ============================================================
# 20. Simplified overlap plot, closest to the supplied style
# ============================================================

plot_overlap_simple <- function(
  overlap_dt,
  title_text,
  output_stem
) {
    pdat <- copy(
        overlap_dt
    )

    pdat[, `:=`(
        One_method_only = (
            Fisher_only +
                edgeR_only
        ),
        Both_methods = Fisher_edgeR
    )]

    pdat <- melt(
        pdat,
        id.vars = c(
            "analysis_group",
            "position",
            "direction"
        ),
        measure.vars = c(
            "One_method_only",
            "Both_methods"
        ),
        variable.name = "support_class",
        value.name = "N"
    )

    pdat[, analysis_group := factor(
        analysis_group,
        levels = analysis_group_order
    )]

    pdat[, direction := factor(
        direction,
        levels = c(
            "I>A",
            "I<A"
        )
    )]

    pdat[, support_class := factor(
        support_class,
        levels = c(
            "One_method_only",
            "Both_methods"
        )
    )]

    colors <- c(
        One_method_only = "#E777A7",
        Both_methods = "#6AACDD"
    )

    labels <- c(
        One_method_only = "One method only",
        Both_methods = "Shared by Fisher and edgeR"
    )

    p <- ggplot(
        pdat,
        aes(
            x = analysis_group,
            y = N,
            fill = support_class
        )
    ) +
        geom_col(
            width = 0.72
        ) +
        facet_grid(
            direction ~ .,
            scales = "free_y",
            space = "free_y"
        ) +
        scale_fill_manual(
            values = colors,
            labels = labels,
            drop = FALSE
        ) +
        scale_x_discrete(
            labels = analysis_group_labels
        ) +
        labs(
            x = NULL,
            y = "Sites number",
            fill = NULL,
            title = title_text
        ) +
        theme_classic(
            base_size = 10
        ) +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                face = "bold",
                size = 11
            ),
            strip.background = element_blank(),
            strip.text.y = element_text(
                face = "bold",
                size = 10
            ),
            axis.text.x = element_text(
                angle = 35,
                hjust = 1,
                size = 8
            ),
            axis.text.y = element_text(
                size = 8
            ),
            axis.title.y = element_text(
                size = 9
            ),
            legend.position = "right",
            legend.text = element_text(
                size = 8
            )
        )

    for (
        xint in c(
            2.5,
            4.5,
            6.5,
            8.5
        )
    ) {
        p <- p +
            geom_vline(
                xintercept = xint,
                linetype = "dashed",
                linewidth = 0.25,
                color = "grey65"
            )
    }

    ggsave(
        file.path(
            opt$output_dir,
            paste0(
                output_stem,
                ".pdf"
            )
        ),
        p,
        width = 12,
        height = 6.5
    )

    ggsave(
        file.path(
            opt$output_dir,
            paste0(
                output_stem,
                ".png"
            )
        ),
        p,
        width = 12,
        height = 6.5,
        dpi = 300
    )

    p
}

# ============================================================
# 21. Generate A-site and AP overlap plots
# ============================================================

overlap_A <- make_overlap_table(
    "A"
)

fwrite(
    overlap_A,
    file.path(
        opt$output_dir,
        "overlap_2methods.Asite.tsv"
    ),
    sep = "\t"
)

plot_overlap_detailed(
    overlap_A,
    "A-site overlap between edgeR QL and Fisher",
    "overlap_2methods.Asite.detailed"
)

plot_overlap_simple(
    overlap_A,
    "A-site concordance between edgeR QL and Fisher",
    "overlap_2methods.Asite.simple"
)

overlap_AP <- make_overlap_table(
    "AP"
)

fwrite(
    overlap_AP,
    file.path(
        opt$output_dir,
        "overlap_2methods.AP_AorP.tsv"
    ),
    sep = "\t"
)

plot_overlap_detailed(
    overlap_AP,
    "A/P-site overlap between edgeR QL and Fisher",
    "overlap_2methods.AP_AorP.detailed"
)

plot_overlap_simple(
    overlap_AP,
    "A/P-site concordance between edgeR QL and Fisher",
    "overlap_2methods.AP_AorP.simple"
)

# ============================================================
# 22. Method diagnostics
# ============================================================

method_diagnostics <- all_results[
    position %in% c(
        "A",
        "P"
    ),
    .(
        n_sites = uniqueN(
            site_id
        ),
        n_testable = uniqueN(
            site_id[
                testable == TRUE
            ]
        ),
        n_converged = uniqueN(
            site_id[
                model_converged == TRUE
            ]
        ),
        n_I_gt_A = uniqueN(
            site_id[
                call == "I>A"
            ]
        ),
        n_I_lt_A = uniqueN(
            site_id[
                call == "I<A"
            ]
        ),
        median_effect_log2OR = median(
            effect_log2OR[
                testable == TRUE
            ],
            na.rm = TRUE
        )
    ),
    by = .(
        method,
        analysis_group,
        analysis_level,
        position
    )
]

method_diagnostics[, order_num := match(
    analysis_group,
    analysis_group_order
)]

setorder(
    method_diagnostics,
    method,
    position,
    order_num
)

method_diagnostics[, order_num := NULL]

fwrite(
    method_diagnostics,
    file.path(
        opt$output_dir,
        "method_diagnostics.tsv"
    ),
    sep = "\t"
)

# ============================================================
# 23. Save parameters
# ============================================================

parameter_table <- data.table(
    parameter = c(
        "Dmel_definition",
        "species",
        "analysis_groups",
        "head_samples",
        "min_replicates",
        "min_rna_per_rep_equiv",
        "min_ribo_per_rep_equiv",
        "min_lm",
        "common_universe_rule",
        "fdr",
        "primary_model",
        "sensitivity_model",
        "reference_model",
        "edgeR_design",
        "beta_binomial_design",
        "AP_definition"
    ),
    value = c(
        "ISO1 + WBg",
        "Dmel,Dsim,Dsec,Dere",
        "female,male,female+male",
        "included as independent biological replicates",
        opt$min_replicates,
        opt$min_rna_per_rep_equiv,
        opt$min_ribo_per_rep_equiv,
        opt$min_lm,
        paste0(
            "at least ", opt$min_replicates,
            " replicates each with RNA and Ribo coverage >= thresholds; ",
            "pooled Lm >= min_lm"
        ),
        opt$fdr,
        "edgeR negative-binomial robust quasi-likelihood",
        "not run (edgeR QL + Fisher workflow)",
        "replicate-wise one-sided Fisher; both directions; joint BH over sites x 2 directions",
        paste0(
            "replicate-assay library fixed effects + ",
            "replicate-specific baseline edited-allele effects + ",
            "common Ribo-specific edited-allele coefficient"
        ),
        "cbind(edited,unedited) ~ replicate + assay; beta-binomial LRT",
        "A significant OR P significant; A/P counts are not summed"
    )
)

fwrite(
    parameter_table,
    file.path(
        opt$output_dir,
        "analysis_parameters.tsv"
    ),
    sep = "\t"
)

# ============================================================
# 24. Completion
# ============================================================

message("")
message("[DONE] Results written to: ", opt$output_dir)

message("")
message("Main outputs:")
message("  analysis_group_replicate_counts.tsv")
message("  common_test_universe.tsv.gz")
message("  site_results.edgeR_QL.tsv.gz")

message("  site_results.Fisher.tsv.gz")
message("  site_results.ALL_methods.tsv.gz")
message("  Table1_summary_long.tsv")
message("  Table1_summary_wide.tsv")
message("  Table1_edgeR_final.tsv")
message("  Table1_edgeR_final_two_row_header.tsv")
message("  overlap_2methods.Asite.detailed.pdf")
message("  overlap_2methods.Asite.simple.pdf")
message("  overlap_2methods.AP_AorP.detailed.pdf")
message("  overlap_2methods.AP_AorP.simple.pdf")
message("  method_diagnostics.tsv")

invisible(NULL)

# Directional Fisher tables keep conflicts explicit and mutually exclusive.
fisher_summary <- all_results[method == 'Fisher' & testable == TRUE, .(
    Total_sites = uniqueN(site_id), I_gt_A = uniqueN(site_id[call == 'I>A']),
    I_lt_A = uniqueN(site_id[call == 'I<A']), Direction_conflict = uniqueN(site_id[call == 'direction_conflict']),
    Not_significant = uniqueN(site_id[call == 'not_significant'])
), by = .(analysis_group, position)]
fisher_summary[, order_num := match(analysis_group, analysis_group_order)]
setorder(fisher_summary, order_num, position); fisher_summary[, order_num := NULL]
fwrite(fisher_summary, file.path(opt$output_dir, 'Table1_Fisher_directional_summary.tsv'), sep = '\t')
fwrite(all_results[method == 'Fisher' & call == 'direction_conflict'],
    file.path(opt$output_dir, 'Fisher_direction_conflicts.tsv'),
    sep = '\t'
)
fwrite(rbindlist(fisher_replicate_audit, fill = TRUE),
    file.path(opt$output_dir, 'Fisher_replicate_directional_tests.tsv.gz'),
    sep = '\t'
)
fwrite(
    data.table(
        parameter = c('Fisher_directions', 'Fisher_BH_scope', 'Fisher_replicate_selection', 'Fisher_conflict_policy', 'Fisher_direction_counts', 'Fisher_AP_statistics'),
        value = c(
            'I>A and I<A, direct one-sided exact tests',
            '2 x number of testable sites, separately in each analysis group x A/P',
            'edgeR common site filter; all group replicates considered; only degenerate rounded tables omitted from combination',
            'both directions significant: exclusive direction_conflict, excluded from I>A and I<A counts',
            'raw unrounded cross-product sign; also report sign after Fisher rounding; no pseudocount',
            'union of A/P directional calls, no new AP hypothesis test; counts refer to representative position'
        )
    ),
    file.path(opt$output_dir, 'Fisher_directional_parameters.tsv'),
    sep = '\t'
)
