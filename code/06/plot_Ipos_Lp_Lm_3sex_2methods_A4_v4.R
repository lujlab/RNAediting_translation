#!/usr/bin/env Rscript

# =============================================================================
# I-position analysis of protein editing vs mRNA editing:
#
#   log2((Lp + eps) / (Lm + eps))
#
# Three sex strata:
#   1) Female
#   2) Male
#   3) Female + Male combined
#
# Each stratum contains TWO methods:
#   A) All sites pool
#      One point per species/site after within-species cohort collapsing.
#
#   B) Conserved sites mean
#      Only conserved site_key units present in >=2 species;
#      Lp and Lm are averaged across species BEFORE calculating the ratio.
#
# Combined-sex definition:
#   Female and Male are first averaged within the SAME species/site:
#       Lp_combined = mean(non-NA Lp_F, Lp_M)
#       Lm_combined = mean(non-NA Lm_F, Lm_M)
#   Then the all-sites / conserved-sites procedures are applied.
#
# Dmel WT:
#   only iso + WBg are included.
#   For a given Dmel site/sex, available iso/WBg Lp and Lm values are averaged
#   before the species/site unit is constructed.
#
# I-position:
#   single nucleotide edit:
#       ((new_cds_pos_1based - 1) %% 3) + 1 -> pos1/pos2/pos3
#   >1 distinct nucleotide positions represented in an AA-level site -> multi
#
# Conserved-site mapping:
#   new_handle -> species-specific CDS projection table ->
#   site_key, following the logic of the supplied reference Ipos script.
#
# Usage:
#   Rscript plot_Ipos_Lp_Lm_3sex_2methods_A4_v3.R \
#     <output_dir> \
#     <cdspos_dir> \
#     <dmel.min3_LmLr.tsv> \
#     <dsim.min3_LmLr.tsv> \
#     <dsec.min3_LmLr.tsv> \
#     <dere.min3_LmLr.tsv>
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
})

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 0. Arguments
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4L) {
    stop(
        paste0(
            "\nUsage:\n",
            "Rscript plot_Ipos_Lp_Lm_3sex_2methods_A4_v3.R \\\n",
            "  <output_dir> \\\n",
            "  <cdspos_dir> \\\n",
            "  <input1.min3_LmLr.tsv> [input2 ...]\n"
        )
    )
}

outdir <- args[[1L]]
cdspos_dir <- args[[2L]]
input_files <- args[-c(1L, 2L)]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(cdspos_dir)) {
    stop("Missing cdspos directory: ", cdspos_dir)
}

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste0("  ", missing_files, collapse = "\n")
    )
}

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------

species_order <- c("Dmel", "Dsim", "Dsec", "Dere")

group_levels <- c(
    "pos1",
    "pos2",
    "pos3",
    "multi"
)

group_labels <- c(
    pos1 = "Pos1",
    pos2 = "Pos2",
    pos3 = "Pos3",
    multi = "Multi-editing"
)

stratum_levels <- c(
    "Female",
    "Male",
    "Female + Male"
)

method_levels <- c(
    "All sites pool",
    "Conserved sites mean"
)

# SAME pseudovalue as the supplied Lr/Lm reference analysis.
eps <- 1e-6

# Reference-like palette, slightly refined for print.
group_cols <- c(
    pos1  = "#E13A9D",
    pos2  = "#D99BC7",
    pos3  = "#DCCDD9",
    multi = "#38B5CB"
)

background_col <- "#F2F2F2"
letter_col <- "#0B1DB3"

BASE_FAMILY <- "Arial"

# -----------------------------------------------------------------------------
# 2. Generic helpers
# -----------------------------------------------------------------------------

as_numeric_na <- function(x) {
    z <- trimws(as.character(x))

    z[
        tolower(z) %in%
            c(
                "",
                "na",
                "nan",
                "null",
                "."
            )
    ] <- NA_character_

    suppressWarnings(
        as.numeric(z)
    )
}

mean_finite <- function(x) {
    x <- as.numeric(x)
    x <- x[
        is.finite(x)
    ]

    if (length(x) == 0L) {
        NA_real_
    } else {
        mean(x)
    }
}

sd_finite <- function(x) {
    x <- as.numeric(x)
    x <- x[
        is.finite(x)
    ]

    if (length(x) < 2L) {
        NA_real_
    } else {
        stats::sd(x)
    }
}

first_nonempty <- function(x) {
    z <- as.character(x)

    z <- z[
        !is.na(z) &
            nzchar(z)
    ]

    if (length(z) == 0L) {
        NA_character_
    } else {
        z[[1L]]
    }
}

format_value <- function(
  x,
  digits = 2L
) {
    ifelse(
        is.finite(x),
        formatC(
            x,
            format = "f",
            digits = digits
        ),
        "NA"
    )
}

canonical <- function(x) {
    tolower(
        gsub(
            "[^A-Za-z0-9]",
            "",
            trimws(
                as.character(x)
            )
        )
    )
}

normalize_genotype <- function(x) {
    z <- canonical(x)

    fifelse(
        z %chin% c(
            "iso",
            "iso1",
            "dmeliso1"
        ),
        "iso",
        fifelse(
            z %chin% c(
                "wbg",
                "dmelwbg"
            ),
            "WBg",
            fifelse(
                z %chin% c(
                    "ag",
                    "adarg",
                    "dmelag"
                ),
                "AG",
                fifelse(
                    z %chin% c(
                        "as",
                        "adars",
                        "dmelas"
                    ),
                    "AS",
                    fifelse(
                        z %chin% c(
                            "sim",
                            "dsim",
                            "dsimulans"
                        ),
                        "sim",
                        fifelse(
                            z %chin% c(
                                "sec",
                                "dsec",
                                "dsechellia"
                            ),
                            "sec",
                            fifelse(
                                z %chin% c(
                                    "ere",
                                    "dere",
                                    "derecta"
                                ),
                                "ere",
                                as.character(x)
                            )
                        )
                    )
                )
            )
        )
    )
}

normalize_species <- function(
  x,
  genotype
) {
    xs <- canonical(x)
    g <- normalize_genotype(genotype)

    out <- rep(
        NA_character_,
        length(g)
    )

    out[
        g %chin%
            c(
                "iso",
                "WBg",
                "AG",
                "AS"
            )
    ] <- "Dmel"

    out[
        g == "sim"
    ] <- "Dsim"

    out[
        g == "sec"
    ] <- "Dsec"

    out[
        g == "ere"
    ] <- "Dere"

    out[
        is.na(out) &
            xs == "dmel"
    ] <- "Dmel"

    out[
        is.na(out) &
            xs == "dsim"
    ] <- "Dsim"

    out[
        is.na(out) &
            xs == "dsec"
    ] <- "Dsec"

    out[
        is.na(out) &
            xs == "dere"
    ] <- "Dere"

    out
}

split_tokens <- function(x) {
    z <- trimws(
        as.character(x)
    )

    if (
        is.na(z) ||
            !nzchar(z)
    ) {
        return(
            character()
        )
    }

    ans <- trimws(
        unlist(
            strsplit(
                z,
                "[;,|]+"
            )
        )
    )

    unique(
        ans[
            !is.na(ans) &
                nzchar(ans)
        ]
    )
}

split_positions <- function(x) {
    toks <- split_tokens(x)

    if (length(toks) == 0L) {
        return(
            integer()
        )
    }

    vals <- suppressWarnings(
        as.integer(toks)
    )

    vals <- vals[
        !is.na(vals)
    ]

    unique(vals)
}

classify_ipos <- function(x) {
    pos <- split_positions(x)

    if (length(pos) == 0L) {
        return(
            NA_character_
        )
    }

    if (length(pos) > 1L) {
        return(
            "multi"
        )
    }

    codon_pos <- (
        (pos[[1L]] - 1L) %%
            3L
    ) + 1L

    paste0(
        "pos",
        codon_pos
    )
}

# -----------------------------------------------------------------------------
# 3. Discover/load CDS projection tables for conserved site_key mapping
# -----------------------------------------------------------------------------

discover_one <- function(
  pattern,
  species_i
) {
    hits <- list.files(
        cdspos_dir,
        pattern = pattern,
        full.names = TRUE
    )

    if (length(hits) == 0L) {
        stop(
            "No CDS projection file found for ",
            species_i,
            " with pattern: ",
            pattern
        )
    }

    if (length(hits) > 1L) {
        stop(
            "Multiple CDS projection files found for ",
            species_i,
            ":\n",
            paste0(
                "  ",
                hits,
                collapse = "\n"
            )
        )
    }

    hits[[1L]]
}

cdspos_files <- c(
    Dmel = discover_one(
        "^ISO1.*cdspos_for_alledit\\.tsv$",
        "Dmel"
    ),
    Dsim = discover_one(
        "^sim.*cdspos_for_alledit\\.tsv$",
        "Dsim"
    ),
    Dsec = discover_one(
        "^sec.*cdspos_for_alledit\\.tsv$",
        "Dsec"
    ),
    Dere = discover_one(
        "^ere.*cdspos_for_alledit\\.tsv$",
        "Dere"
    )
)

read_site_key_map <- function(
  species_i,
  file_i
) {
    x <- fread(
        file_i,
        sep = "\t",
        na.strings = c(
            "",
            "NA",
            "NaN",
            "na"
        )
    )

    needed <- c(
        "new_transcript",
        "new_cds_pos_1based",
        "site_key"
    )

    miss <- setdiff(
        needed,
        names(x)
    )

    if (length(miss) > 0L) {
        stop(
            "Missing columns in ",
            file_i,
            ": ",
            paste(
                miss,
                collapse = ", "
            )
        )
    }

    cds_pos <- as_numeric_na(
        x$new_cds_pos_1based
    )

    out <- data.table(
        species = species_i,
        new_handle = paste(
            as.character(
                x$new_transcript
            ),
            as.integer(
                cds_pos
            ),
            sep = "_"
        ),
        site_key = as.character(
            x$site_key
        )
    )

    unique(
        out[
            !is.na(new_handle) &
                new_handle != "NA_NA" &
                !is.na(site_key) &
                site_key != ""
        ]
    )
}

site_key_map <- rbindlist(
    lapply(
        names(cdspos_files),
        function(sp) {
            read_site_key_map(
                sp,
                cdspos_files[[sp]]
            )
        }
    ),
    use.names = TRUE,
    fill = TRUE
)

setkey(
    site_key_map,
    species,
    new_handle
)

map_new_handle_to_site_key <- function(
  handles_i,
  species_i
) {
    ids <- split_tokens(
        handles_i
    )

    if (length(ids) == 0L) {
        return(
            NA_character_
        )
    }

    z <- site_key_map[
        .(
            species_i,
            ids
        ),
        on = .(
            species,
            new_handle
        ),
        nomatch = 0L,
        site_key
    ]

    z <- sort(
        unique(
            as.character(
                z[
                    !is.na(z) &
                        nzchar(z)
                ]
            )
        )
    )

    if (length(z) == 0L) {
        NA_character_
    } else {
        paste(
            z,
            collapse = ","
        )
    }
}

# -----------------------------------------------------------------------------
# 4. Read new min3 Lp/Lm files
# -----------------------------------------------------------------------------

dt_list <- lapply(
    input_files,
    function(f) {
        x <- fread(
            f,
            sep = "\t",
            header = TRUE,
            quote = "\"",
            na.strings = c(
                "",
                "NA",
                "NaN",
                "nan",
                "null"
            )
        )

        x[
            ,
            source_file := f
        ]

        x
    }
)

dt <- rbindlist(
    dt_list,
    fill = TRUE,
    use.names = TRUE
)

if (!"genotype" %in% names(dt)) {
    if (
        "group" %in%
            names(dt)
    ) {
        dt[
            ,
            genotype :=
                as.character(group)
        ]
    } else {
        stop(
            "Input tables contain neither `genotype` nor `group`."
        )
    }
}

required_cols <- c(
    "site",
    "genotype",
    "coordinate_species",
    "new_cds_pos_1based",
    "new_handle",
    "Lp_mean_F",
    "Lp_mean_M",
    "mean_Lm_F",
    "mean_Lm_M"
)

missing_cols <- setdiff(
    required_cols,
    names(dt)
)

if (
    length(missing_cols) > 0L
) {
    stop(
        "Missing required columns:\n  ",
        paste(
            missing_cols,
            collapse = "\n  "
        )
    )
}

dt[
    ,
    genotype :=
        normalize_genotype(
            genotype
        )
]

dt[
    ,
    species :=
        normalize_species(
            coordinate_species,
            genotype
        )
]

for (
    cc in c(
        "Lp_mean_F",
        "Lp_mean_M",
        "mean_Lm_F",
        "mean_Lm_M"
    )
) {
    set(
        dt,
        j = cc,
        value = as_numeric_na(
            dt[[cc]]
        )
    )
}

dt[
    ,
    ipos :=
        vapply(
            new_cds_pos_1based,
            classify_ipos,
            character(1L)
        )
]

dt[
    ,
    conserved_key :=
        vapply(
            seq_len(.N),
            function(i) {
                map_new_handle_to_site_key(
                    new_handle[[i]],
                    species[[i]]
                )
            },
            character(1L)
        )
]

# -----------------------------------------------------------------------------
# 5. WT filtering
#
# Dmel: iso + WBg only.
# Other species: their WT cohort.
# -----------------------------------------------------------------------------

dt_use <- dt[
    (
        species == "Dmel" &
            genotype %chin%
                c(
                    "iso",
                    "WBg"
                )
    ) |
        (
            species == "Dsim" &
                genotype == "sim"
        ) |
        (
            species == "Dsec" &
                genotype == "sec"
        ) |
        (
            species == "Dere" &
                genotype == "ere"
        )
]

if (
    nrow(dt_use) == 0L
) {
    stop(
        "No rows remain after WT species/genotype filtering."
    )
}

# -----------------------------------------------------------------------------
# 6. Wide F/M -> long
# -----------------------------------------------------------------------------

female_long <- dt_use[
    ,
    .(
        species,
        site,
        genotype,
        ipos,
        conserved_key,
        new_handle,
        Lp = Lp_mean_F,
        Lm = mean_Lm_F
    )
]

female_long[
    ,
    sex := "Female"
]

male_long <- dt_use[
    ,
    .(
        species,
        site,
        genotype,
        ipos,
        conserved_key,
        new_handle,
        Lp = Lp_mean_M,
        Lm = mean_Lm_M
    )
]

male_long[
    ,
    sex := "Male"
]

sex_long <- rbindlist(
    list(
        female_long,
        male_long
    ),
    use.names = TRUE,
    fill = TRUE
)

# -----------------------------------------------------------------------------
# 7. One species/site/sex unit
#
# For Dmel this collapses iso + WBg first.
# For other species there is normally one genotype row.
#
# IMPORTANT:
# Lp and Lm are averaged independently.
# -----------------------------------------------------------------------------

species_site_sex <- sex_long[
    ,
    .(
        Lp = mean_finite(Lp),
        Lm = mean_finite(Lm),
        n_genotype_rows = .N,
        genotypes_used = paste(
            sort(
                unique(
                    genotype[
                        !is.na(genotype)
                    ]
                )
            ),
            collapse = ","
        ),
        conserved_key =
            first_nonempty(
                conserved_key
            ),
        new_handle =
            first_nonempty(
                new_handle
            )
    ),
    by = .(
        species,
        sex,
        site,
        ipos
    )
]

# -----------------------------------------------------------------------------
# 8. Combined sex:
#
# Female + Male are averaged within SAME species/site first.
# This avoids counting the same biological editing site twice.
# -----------------------------------------------------------------------------

combined_site <- species_site_sex[
    ,
    .(
        Lp = mean_finite(Lp),
        Lm = mean_finite(Lm),
        n_sexes_available =
            sum(
                is.finite(Lp) |
                    is.finite(Lm)
            ),
        sexes_used = paste(
            sort(
                unique(
                    sex[
                        is.finite(Lp) |
                            is.finite(Lm)
                    ]
                )
            ),
            collapse = ","
        ),
        conserved_key =
            first_nonempty(
                conserved_key
            ),
        new_handle =
            first_nonempty(
                new_handle
            )
    ),
    by = .(
        species,
        site,
        ipos
    )
]

combined_site[
    ,
    sex := "Female + Male"
]

setcolorder(
    combined_site,
    c(
        "species",
        "sex",
        "site",
        "ipos",
        "conserved_key",
        "new_handle",
        "Lp",
        "Lm",
        "n_sexes_available",
        "sexes_used"
    )
)

# Harmonize columns.
species_site_sex[
    ,
    `:=`(
        n_sexes_available = NA_integer_,
        sexes_used = sex
    )
]

all_strata_species_site <- rbindlist(
    list(
        species_site_sex,
        combined_site
    ),
    fill = TRUE,
    use.names = TRUE
)

all_strata_species_site[
    ,
    sex :=
        factor(
            sex,
            levels = stratum_levels
        )
]

# -----------------------------------------------------------------------------
# 9. Method A: All sites pool
#
# One point per species/site for each sex stratum.
# -----------------------------------------------------------------------------

pooled_sites <- copy(
    all_strata_species_site[
        is.finite(Lp) &
            is.finite(Lm) &
            ipos %chin%
                group_levels
    ]
)

pooled_sites[
    ,
    lp_lm :=
        log2(
            (Lp + eps) /
                (Lm + eps)
        )
]

pooled_sites[
    ,
    method :=
        "All sites pool"
]

pooled_sites[
    ,
    unit_id :=
        paste(
            as.character(sex),
            species,
            site,
            sep = "__"
        )
]

pooled_sites <- pooled_sites[
    is.finite(lp_lm)
]

if (
    nrow(pooled_sites) == 0L
) {
    stop(
        "No All-sites-pool values remain."
    )
}

# -----------------------------------------------------------------------------
# 10. Method B: Conserved sites mean
#
# Same logic as the reference:
#   conserved_key must occur in >=2 species;
#   I-position group remains in aggregation key;
#   average Lp and Lm ACROSS SPECIES;
#   then calculate log2(Lp/Lm).
# -----------------------------------------------------------------------------

key_group_qc <- all_strata_species_site[
    !is.na(conserved_key) &
        conserved_key != "" &
        ipos %chin%
            group_levels,
    .(
        n_species =
            uniqueN(species),
        groups = paste(
            sort(
                unique(ipos)
            ),
            collapse = ","
        ),
        n_groups =
            uniqueN(ipos)
    ),
    by = .(
        sex,
        conserved_key
    )
]

group_conflicts <- key_group_qc[
    n_groups > 1L
]

conserved_sites <- all_strata_species_site[
    !is.na(conserved_key) &
        conserved_key != "" &
        ipos %chin%
            group_levels,
    .(
        Lp =
            mean_finite(Lp),
        Lm =
            mean_finite(Lm),
        n_species =
            uniqueN(species),
        species_present = paste(
            sort(
                unique(species)
            ),
            collapse = ","
        ),
        n_species_contributions =
            .N,
        source_species_sites = paste(
            sort(
                unique(
                    paste(
                        species,
                        site,
                        sep = "__"
                    )
                )
            ),
            collapse = ","
        )
    ),
    by = .(
        sex,
        conserved_key,
        ipos
    )
]

conserved_sites <- conserved_sites[
    n_species >= 2L &
        is.finite(Lp) &
        is.finite(Lm)
]

conserved_sites[
    ,
    lp_lm :=
        log2(
            (Lp + eps) /
                (Lm + eps)
        )
]

conserved_sites[
    ,
    method :=
        "Conserved sites mean"
]

conserved_sites[
    ,
    unit_id :=
        paste(
            as.character(sex),
            conserved_key,
            ipos,
            sep = "__"
        )
]

conserved_sites[
    ,
    conserved_group_conflict :=
        paste(
            as.character(sex),
            conserved_key,
            sep = "__"
        ) %chin%
        paste(
            as.character(
                group_conflicts$sex
            ),
            group_conflicts$conserved_key,
            sep = "__"
        )
]

conserved_sites <- conserved_sites[
    is.finite(lp_lm)
]

if (
    nrow(conserved_sites) == 0L
) {
    stop(
        "No conserved sites shared by >=2 species remain."
    )
}

# -----------------------------------------------------------------------------
# 11. Save method-level data
# -----------------------------------------------------------------------------

fwrite(
    site_key_map,
    file.path(
        outdir,
        "conserved_site_key_map.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    all_strata_species_site,
    file.path(
        outdir,
        "species_site_values_three_sex_strata.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    pooled_sites,
    file.path(
        outdir,
        "Lp_Lm_All_sites_pool_three_sex_strata.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    conserved_sites,
    file.path(
        outdir,
        "Lp_Lm_Conserved_sites_mean_three_sex_strata.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    key_group_qc,
    file.path(
        outdir,
        "conserved_site_group_QC.tsv"
    ),
    sep = "\t",
    na = "NA"
)

# -----------------------------------------------------------------------------
# 12. Pairwise Wilcoxon + BH
# -----------------------------------------------------------------------------

unpaired_wilcoxon_bh <- function(
  d,
  stratum_name,
  method_name
) {
    d <- d[
        is.finite(lp_lm) &
            ipos %chin%
                group_levels
    ]

    groups <- group_levels[
        group_levels %chin%
            unique(d$ipos)
    ]

    if (
        length(groups) < 2L
    ) {
        return(
            data.table(
                sex_stratum =
                    stratum_name,
                method =
                    method_name,
                group1 =
                    character(),
                group2 =
                    character(),
                statistic =
                    numeric(),
                p_value =
                    numeric(),
                p_adj =
                    numeric(),
                n_group1 =
                    integer(),
                n_group2 =
                    integer(),
                test =
                    character()
            )
        )
    }

    pairs <- combn(
        groups,
        2L,
        simplify = FALSE
    )

    ans <- rbindlist(
        lapply(
            pairs,
            function(pair) {
                x <- d[
                    ipos == pair[[1L]],
                    lp_lm
                ]

                y <- d[
                    ipos == pair[[2L]],
                    lp_lm
                ]

                test_i <- if (
                    length(x) >= 2L &&
                        length(y) >= 2L
                ) {
                    tryCatch(
                        wilcox.test(
                            x,
                            y,
                            paired = FALSE,
                            exact = FALSE
                        ),
                        error = function(e) NULL
                    )
                } else {
                    NULL
                }

                data.table(
                    sex_stratum =
                        stratum_name,
                    method =
                        method_name,
                    group1 =
                        pair[[1L]],
                    group2 =
                        pair[[2L]],
                    statistic =
                        if (
                            is.null(test_i)
                        ) {
                            NA_real_
                        } else {
                            unname(
                                as.numeric(
                                    test_i$statistic
                                )
                            )
                        },
                    p_value =
                        if (
                            is.null(test_i)
                        ) {
                            NA_real_
                        } else {
                            test_i$p.value
                        },
                    n_group1 =
                        length(x),
                    n_group2 =
                        length(y),
                    test =
                        "unpaired Wilcoxon rank-sum"
                )
            }
        ),
        fill = TRUE
    )

    ans[
        ,
        p_adj :=
            p.adjust(
                p_value,
                method = "BH"
            )
    ]

    ans[]
}

# -----------------------------------------------------------------------------
# 13. Compact-letter display
# -----------------------------------------------------------------------------

compact_letters <- function(
  groups,
  pairwise,
  alpha = 0.05
) {
    groups <- as.character(groups)

    n_groups <- length(groups)

    if (
        n_groups == 0L
    ) {
        return(
            data.table(
                ipos =
                    character(),
                cld =
                    character()
            )
        )
    }

    significant <- matrix(
        FALSE,
        nrow = n_groups,
        ncol = n_groups,
        dimnames = list(
            groups,
            groups
        )
    )

    if (
        nrow(pairwise) > 0L
    ) {
        for (
            i in seq_len(
                nrow(pairwise)
            )
        ) {
            g1 <- as.character(
                pairwise$group1[[i]]
            )

            g2 <- as.character(
                pairwise$group2[[i]]
            )

            if (
                g1 %in% groups &&
                    g2 %in% groups &&
                    is.finite(
                        pairwise$p_adj[[i]]
                    ) &&
                    pairwise$p_adj[[i]] <
                        alpha
            ) {
                significant[
                    g1,
                    g2
                ] <- TRUE

                significant[
                    g2,
                    g1
                ] <- TRUE
            }
        }
    }

    candidate_cliques <- list()

    if (
        n_groups >= 2L
    ) {
        for (
            size in 2:n_groups
        ) {
            candidate_cliques <- c(
                candidate_cliques,
                combn(
                    groups,
                    size,
                    simplify = FALSE
                )
            )
        }
    }

    candidate_cliques <- Filter(
        function(clique) {
            idx <- match(
                clique,
                groups
            )

            all(
                !significant[
                    idx,
                    idx,
                    drop = FALSE
                ]
            )
        },
        candidate_cliques
    )

    if (
        length(candidate_cliques) >
            0L
    ) {
        is_subsumed <- vapply(
            seq_along(
                candidate_cliques
            ),
            function(i) {
                any(
                    vapply(
                        seq_along(
                            candidate_cliques
                        ),
                        function(j) {
                            j != i &&
                                length(
                                    candidate_cliques[[j]]
                                ) >
                                    length(
                                        candidate_cliques[[i]]
                                    ) &&
                                all(
                                    candidate_cliques[[i]] %in%
                                        candidate_cliques[[j]]
                                )
                        },
                        logical(1)
                    )
                )
            },
            logical(1)
        )

        candidate_cliques <-
            candidate_cliques[
                !is_subsumed
            ]
    }

    covered <- unique(
        unlist(
            candidate_cliques,
            use.names = FALSE
        )
    )

    uncovered <- setdiff(
        groups,
        covered
    )

    if (
        length(uncovered) >
            0L
    ) {
        candidate_cliques <- c(
            candidate_cliques,
            as.list(
                uncovered
            )
        )
    }

    if (
        length(candidate_cliques) ==
            0L
    ) {
        candidate_cliques <-
            as.list(
                groups
            )
    }

    label_letters <- letters[
        seq_along(
            candidate_cliques
        )
    ]

    data.table(
        ipos =
            groups,
        cld =
            vapply(
                groups,
                function(g) {
                    paste(
                        label_letters[
                            vapply(
                                candidate_cliques,
                                function(clique) {
                                    g %in% clique
                                },
                                logical(1)
                            )
                        ],
                        collapse = ""
                    )
                },
                character(1)
            )
    )
}

# -----------------------------------------------------------------------------
# 14. Panel summary
# -----------------------------------------------------------------------------

make_summary <- function(
  d,
  stratum_name,
  method_name,
  pairwise
) {
    s <- d[
        is.finite(lp_lm),
        .(
            n_points =
                .N,
            mean_lp_lm =
                mean(lp_lm),
            sd_lp_lm =
                sd_finite(
                    lp_lm
                )
        ),
        by = ipos
    ]

    s <- merge(
        data.table(
            ipos =
                group_levels
        ),
        s,
        by = "ipos",
        all.x = TRUE,
        sort = FALSE
    )

    cld <- compact_letters(
        group_levels,
        pairwise
    )

    s <- merge(
        s,
        cld,
        by = "ipos",
        all.x = TRUE,
        sort = FALSE
    )

    s[
        ,
        `:=`(
            sex_stratum =
                stratum_name,
            method =
                method_name
        )
    ]

    s[
        ,
        ipos :=
            factor(
                ipos,
                levels =
                    group_levels
            )
    ]

    setorder(
        s,
        ipos
    )

    s[]
}

# -----------------------------------------------------------------------------
# 15. Build six statistical panels
# -----------------------------------------------------------------------------

panel_data <- list()
panel_pw <- list()
panel_summary <- list()

counter <- 0L

for (
    stratum_i in
    stratum_levels
) {
    d_pool <- pooled_sites[
        as.character(sex) ==
            stratum_i
    ]

    d_cons <- conserved_sites[
        as.character(sex) ==
            stratum_i
    ]

    for (
        method_i in
        method_levels
    ) {
        counter <- counter + 1L

        d_i <- if (
            method_i ==
                "All sites pool"
        ) {
            d_pool
        } else {
            d_cons
        }

        pw_i <- unpaired_wilcoxon_bh(
            d_i,
            stratum_i,
            method_i
        )

        sum_i <- make_summary(
            d_i,
            stratum_i,
            method_i,
            pw_i
        )

        key_i <- paste(
            stratum_i,
            method_i,
            sep = "__"
        )

        panel_data[[key_i]] <- d_i
        panel_pw[[key_i]] <- pw_i
        panel_summary[[key_i]] <- sum_i
    }
}

pairwise_all <- rbindlist(
    panel_pw,
    fill = TRUE
)

summary_all <- rbindlist(
    panel_summary,
    fill = TRUE
)

fwrite(
    pairwise_all,
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_Wilcoxon_BH.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    summary_all,
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_group_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

# -----------------------------------------------------------------------------
# 16. Reference-style plot function
# -----------------------------------------------------------------------------

make_panel <- function(
  d,
  summary_dt,
  title_text
) {
    d <- copy(
        d[
            is.finite(lp_lm) &
                ipos %chin%
                    group_levels
        ]
    )

    d[
        ,
        ipos :=
            factor(
                ipos,
                levels =
                    group_levels
            )
    ]

    if (
        nrow(d) == 0L
    ) {
        stop(
            "No finite values to plot for: ",
            title_text
        )
    }

    # Fixed y-axis range.
    y_min <- -20
    y_max <- 15
    y_range <- y_max - y_min

    # Match visual spacing of reference figure:
    # CLD at top, then mean<U+00B1>SD, then italic n.
    ann <- copy(
        summary_dt
    )

    ann[
        ,
        ipos :=
            factor(
                ipos,
                levels =
                    group_levels
            )
    ]

    ann[
        ,
        `:=`(
            y_sig = 14.0,
            y_mean_sd = 12.0,
            y_n = 10.0,
            n_label =
                paste0(
                    "italic(n)==",
                    n_points
                ),
            mean_sd_label =
                paste0(
                    format_value(
                        mean_lp_lm,
                        2L
                    ),
                    " <U+00B1> ",
                    format_value(
                        sd_lp_lm,
                        2L
                    )
                )
        )
    ]

    bracket_y <- -25.5

    bracket_text_y <- -28.0

    ggplot(
        d,
        aes(
            x = ipos,
            y = lp_lm,
            fill = ipos
        )
    ) +
        geom_hline(
            yintercept = 0,
            linetype = 2,
            linewidth = 0.45,
            color = "#BBBBBB"
        ) +
        geom_violin(
            width = 0.88,
            trim = TRUE,
            scale = "width",
            alpha = 0.98,
            color = NA
        ) +
        geom_boxplot(
            width = 0.15,
            outlier.shape = NA,
            fill = "white",
            color = "black",
            linewidth = 0.60
        ) +
        geom_text(
            data = ann,
            aes(
                x = ipos,
                y = y_sig,
                label = cld
            ),
            inherit.aes = FALSE,
            color = letter_col,
            size = 8.0 /
                ggplot2::.pt,
            family = BASE_FAMILY,
            fontface = "bold",
            vjust = 0
        ) +
        geom_text(
            data = ann,
            aes(
                x = ipos,
                y = y_mean_sd,
                label = mean_sd_label
            ),
            inherit.aes = FALSE,
            size = 7.5 /
                ggplot2::.pt,
            family = BASE_FAMILY,
            fontface = "plain",
            vjust = 0
        ) +
        geom_text(
            data = ann,
            aes(
                x = ipos,
                y = y_n,
                label = n_label
            ),
            inherit.aes = FALSE,
            parse = TRUE,
            size = 7.5 /
                ggplot2::.pt,
            family = BASE_FAMILY,
            fontface = "italic",
            vjust = 0
        ) +
        annotate(
            "segment",
            x = 0.82,
            xend = 3.18,
            y = bracket_y,
            yend = bracket_y,
            linewidth = 0.55,
            color = "black"
        ) +
        annotate(
            "text",
            x = 2,
            y = bracket_text_y,
            label = "Single-editing",
            size = 9.0 /
                ggplot2::.pt,
            family = BASE_FAMILY,
            fontface = "plain",
            vjust = 0
        ) +
        scale_fill_manual(
            values =
                group_cols,
            drop = FALSE
        ) +
        scale_x_discrete(
            drop = FALSE,
            labels =
                group_labels
        ) +
        scale_y_continuous(
            expand =
                expansion(
                    mult =
                        c(
                            0,
                            0
                        )
                )
        ) +
        coord_cartesian(
            ylim = c(-20, 15),
            clip = "off"
        ) +
        labs(
            title =
                title_text,
            x =
                NULL,
            y =
                expression(
                    log[2](
                        italic(L)[p] /
                            italic(L)[m]
                    )
                )
        ) +
        theme_classic(
            base_size = 10.5,
            base_family =
                BASE_FAMILY
        ) +
        theme(
            plot.background =
                element_rect(
                    fill =
                        background_col,
                    color = NA
                ),
            panel.background =
                element_rect(
                    fill =
                        background_col,
                    color = NA
                ),
            text =
                element_text(
                    family =
                        BASE_FAMILY,
                    color = "black"
                ),
            axis.title.y =
                element_text(
                    size = 9.5,
                    face = "italic",
                    margin =
                        margin(
                            r = 4
                        )
                ),
            axis.text.x =
                element_text(
                    size = 9.0,
                    face = "plain"
                ),
            axis.text.y =
                element_text(
                    size = 9.0,
                    face = "plain"
                ),
            axis.line =
                element_line(
                    linewidth = 0.60,
                    color = "black"
                ),
            axis.ticks =
                element_line(
                    linewidth = 0.60,
                    color = "black"
                ),
            axis.ticks.length =
                unit(
                    0.10,
                    "cm"
                ),
            plot.title =
                element_text(
                    size = 11.5,
                    face = "bold",
                    hjust = 0.5,
                    margin =
                        margin(
                            b = 5
                        )
                ),
            legend.position =
                "none",
            plot.margin =
                margin(
                    12,
                    10,
                    26,
                    12
                )
        )
}

# -----------------------------------------------------------------------------
# 17. Build six plots
# -----------------------------------------------------------------------------

plots <- list()

for (
    stratum_i in
    stratum_levels
) {
    for (
        method_i in
        method_levels
    ) {
        key_i <- paste(
            stratum_i,
            method_i,
            sep = "__"
        )

        plots[[key_i]] <-
            make_panel(
                panel_data[[key_i]],
                panel_summary[[key_i]],
                method_i
            )
    }
}

# -----------------------------------------------------------------------------
# 18. Main A4: 3 rows x 2 columns
#
# Row 1: Female
# Row 2: Male
# Row 3: Female + Male
#
# Col 1: All sites pool
# Col 2: Conserved sites mean
# -----------------------------------------------------------------------------

main_pdf <- file.path(
    outdir,
    "Ipos_Lp_over_Lm_3sex_x_2methods_A4_portrait.pdf"
)

grDevices::cairo_pdf(
    main_pdf,
    width = 8.27,
    height = 11.69,
    onefile = TRUE,
    family =
        BASE_FAMILY
)

grid::grid.newpage()

# Main body.
grid::pushViewport(
    grid::viewport(
        x = 0.52,
        y = 0.50,
        width =
            unit(
                0.94,
                "npc"
            ),
        height =
            unit(
                0.91,
                "npc"
            ),
        layout =
            grid::grid.layout(
                nrow = 3L,
                ncol = 2L,
                heights =
                    unit(
                        c(
                            1,
                            1,
                            1
                        ),
                        "null"
                    ),
                widths =
                    unit(
                        c(
                            1,
                            1
                        ),
                        "null"
                    )
            ),
        clip = "off"
    )
)

for (
    r in
    seq_along(
        stratum_levels
    )
) {
    stratum_i <-
        stratum_levels[[r]]

    # Left.
    key_left <- paste(
        stratum_i,
        "All sites pool",
        sep = "__"
    )

    grid::pushViewport(
        grid::viewport(
            layout.pos.row =
                r,
            layout.pos.col =
                1L,
            clip = "off"
        )
    )

    grid::grid.draw(
        ggplot2::ggplotGrob(
            plots[[key_left]]
        )
    )

    grid::popViewport()

    # Right.
    key_right <- paste(
        stratum_i,
        "Conserved sites mean",
        sep = "__"
    )

    grid::pushViewport(
        grid::viewport(
            layout.pos.row =
                r,
            layout.pos.col =
                2L,
            clip = "off"
        )
    )

    grid::grid.draw(
        ggplot2::ggplotGrob(
            plots[[key_right]]
        )
    )

    grid::popViewport()
}

grid::popViewport()

# Row labels placed on left edge.
row_y <- c(
    0.845,
    0.525,
    0.205
)

for (
    i in seq_along(
        stratum_levels
    )
) {
    grid::grid.text(
        stratum_levels[[i]],
        x =
            unit(
                0.012,
                "npc"
            ),
        y =
            unit(
                row_y[[i]],
                "npc"
            ),
        rot = 90,
        just = "center",
        gp =
            grid::gpar(
                fontfamily =
                    BASE_FAMILY,
                fontsize = 11,
                fontface = "bold"
            )
    )
}

grDevices::dev.off()

# -----------------------------------------------------------------------------
# 19. Additional analysis: exclude final sites with Lp == 0 or Lp == 1.
#
# The filter is applied to the final plotted site unit:
#   All sites pool      -> after species/site/sex collapsing
#   Conserved sites mean -> after cross-species Lp/Lm averaging
# This retains the strict interior 0 < Lp < 1 without changing the original
# analysis above.
# -----------------------------------------------------------------------------

keep_interior_lp <- function(d) {
    d[
        is.finite(Lp) &
            Lp > 0 &
            Lp < 1
    ]
}

pooled_sites_lp_interior <- keep_interior_lp(pooled_sites)
conserved_sites_lp_interior <- keep_interior_lp(conserved_sites)

lp_filter_summary <- rbindlist(
    lapply(
        stratum_levels,
        function(stratum_i) {
            rbindlist(
                lapply(
                    method_levels,
                    function(method_i) {
                        d_all <- if (method_i == "All sites pool") {
                            pooled_sites[as.character(sex) == stratum_i]
                        } else {
                            conserved_sites[as.character(sex) == stratum_i]
                        }
                        d_keep <- if (method_i == "All sites pool") {
                            pooled_sites_lp_interior[as.character(sex) == stratum_i]
                        } else {
                            conserved_sites_lp_interior[as.character(sex) == stratum_i]
                        }
                        data.table(
                            sex_stratum = stratum_i,
                            method = method_i,
                            n_original = nrow(d_all),
                            n_retained = nrow(d_keep),
                            n_excluded = nrow(d_all) - nrow(d_keep),
                            filter = "0 < final Lp < 1"
                        )
                    }
                ),
                fill = TRUE
            )
        }
    ),
    fill = TRUE
)

fwrite(
    pooled_sites_lp_interior,
    file.path(
        outdir,
        "Lp_Lm_All_sites_pool_three_sex_strata_Lp_between_0_and_1.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    conserved_sites_lp_interior,
    file.path(
        outdir,
        "Lp_Lm_Conserved_sites_mean_three_sex_strata_Lp_between_0_and_1.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    lp_filter_summary,
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_Lp_filter_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

panel_data_lp_interior <- list()
panel_pw_lp_interior <- list()
panel_summary_lp_interior <- list()

for (stratum_i in stratum_levels) {
    d_pool <- pooled_sites_lp_interior[
        as.character(sex) == stratum_i
    ]
    d_cons <- conserved_sites_lp_interior[
        as.character(sex) == stratum_i
    ]

    for (method_i in method_levels) {
        d_i <- if (method_i == "All sites pool") d_pool else d_cons
        pw_i <- unpaired_wilcoxon_bh(d_i, stratum_i, method_i)
        sum_i <- make_summary(d_i, stratum_i, method_i, pw_i)
        key_i <- paste(stratum_i, method_i, sep = "__")
        panel_data_lp_interior[[key_i]] <- d_i
        panel_pw_lp_interior[[key_i]] <- pw_i
        panel_summary_lp_interior[[key_i]] <- sum_i
    }
}

pairwise_lp_interior <- rbindlist(
    panel_pw_lp_interior,
    fill = TRUE
)

summary_lp_interior <- rbindlist(
    panel_summary_lp_interior,
    fill = TRUE
)

fwrite(
    pairwise_lp_interior,
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_Wilcoxon_BH_Lp_between_0_and_1.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    summary_lp_interior,
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_group_summary_Lp_between_0_and_1.tsv"
    ),
    sep = "\t",
    na = "NA"
)

plots_lp_interior <- list()

for (stratum_i in stratum_levels) {
    for (method_i in method_levels) {
        key_i <- paste(stratum_i, method_i, sep = "__")
        plots_lp_interior[[key_i]] <- make_panel(
            panel_data_lp_interior[[key_i]],
            panel_summary_lp_interior[[key_i]],
            paste0(method_i, " (0 < Lp < 1)")
        )
    }
}

write_a4_pdf <- function(plot_list, pdf_file) {
    grDevices::cairo_pdf(
        pdf_file,
        width = 8.27,
        height = 11.69,
        onefile = TRUE,
        family = BASE_FAMILY
    )

    grid::grid.newpage()
    grid::pushViewport(
        grid::viewport(
            x = 0.52,
            y = 0.50,
            width = unit(0.94, "npc"),
            height = unit(0.91, "npc"),
            layout = grid::grid.layout(
                nrow = 3L,
                ncol = 2L,
                heights = unit(c(1, 1, 1), "null"),
                widths = unit(c(1, 1), "null")
            ),
            clip = "off"
        )
    )

    for (r in seq_along(stratum_levels)) {
        stratum_i <- stratum_levels[[r]]
        for (c in seq_along(method_levels)) {
            method_i <- method_levels[[c]]
            key_i <- paste(stratum_i, method_i, sep = "__")
            grid::pushViewport(
                grid::viewport(
                    layout.pos.row = r,
                    layout.pos.col = c,
                    clip = "off"
                )
            )
            grid::grid.draw(ggplot2::ggplotGrob(plot_list[[key_i]]))
            grid::popViewport()
        }
    }

    grid::popViewport()

    row_y <- c(0.845, 0.525, 0.205)
    for (i in seq_along(stratum_levels)) {
        grid::grid.text(
            stratum_levels[[i]],
            x = unit(0.012, "npc"),
            y = unit(row_y[[i]], "npc"),
            rot = 90,
            just = "center",
            gp = grid::gpar(
                fontfamily = BASE_FAMILY,
                fontsize = 11,
                fontface = "bold"
            )
        )
    }

    grDevices::dev.off()
}

lp_interior_pdf <- file.path(
    outdir,
    "Ipos_Lp_over_Lm_3sex_x_2methods_A4_portrait_Lp_between_0_and_1.pdf"
)

write_a4_pdf(
    plots_lp_interior,
    lp_interior_pdf
)

# -----------------------------------------------------------------------------
# 20. Run summary
# -----------------------------------------------------------------------------

method_summary <- rbindlist(
    lapply(
        stratum_levels,
        function(stratum_i) {
            rbindlist(
                lapply(
                    method_levels,
                    function(method_i) {
                        key_i <- paste(
                            stratum_i,
                            method_i,
                            sep = "__"
                        )

                        d_i <- panel_data[[key_i]]

                        data.table(
                            sex_stratum =
                                stratum_i,
                            method =
                                method_i,
                            n_total =
                                nrow(d_i),
                            n_pos1 =
                                sum(
                                    d_i$ipos ==
                                        "pos1",
                                    na.rm = TRUE
                                ),
                            n_pos2 =
                                sum(
                                    d_i$ipos ==
                                        "pos2",
                                    na.rm = TRUE
                                ),
                            n_pos3 =
                                sum(
                                    d_i$ipos ==
                                        "pos3",
                                    na.rm = TRUE
                                ),
                            n_multi =
                                sum(
                                    d_i$ipos ==
                                        "multi",
                                    na.rm = TRUE
                                )
                        )
                    }
                ),
                fill = TRUE
            )
        }
    ),
    fill = TRUE
)

fwrite(
    method_summary,
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_run_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

# -----------------------------------------------------------------------------
# 21. Manifest
# -----------------------------------------------------------------------------

fwrite(
    data.table(
        item = c(
            "metric",
            "pseudovalue",
            "female_definition",
            "male_definition",
            "combined_sex_definition",
            "Dmel_definition",
            "all_sites_pool",
            "conserved_sites_mean",
            "conserved_filter",
            "statistics",
            "I_position"
        ),
        value = c(
            "log2((Lp + 1e-6)/(Lm + 1e-6))",
            "1e-6; same as supplied Lr/Lm reference code",
            "female Lp and Lm only",
            "male Lp and Lm only",
            "mean of available Female/Male Lp and Lm within same species/site before ratio",
            "iso + WBg averaged within species/site/sex; AG/AS excluded",
            "one point per species/site within stratum",
            "Lp and Lm averaged across species for conserved site_key + I-position, then ratio",
            "site_key present in >=2 species",
            "unpaired Wilcoxon rank-sum across I-position groups; BH within each panel",
            "pos1/pos2/pos3 from new CDS coordinate; >1 nucleotide position = multi"
        )
    ),
    file.path(
        outdir,
        "Lp_Lm_Ipos_three_sex_two_methods_manifest.tsv"
    ),
    sep = "\t",
    na = "NA"
)

message(
    "Done.\n",
    "Main 3x2 PDF: ",
    normalizePath(
        main_pdf,
        mustWork = FALSE
    ),
    "\nLp 0/1-excluded PDF: ",
    normalizePath(
        lp_interior_pdf,
        mustWork = FALSE
    ),
    "\nOutput directory: ",
    normalizePath(
        outdir,
        mustWork = FALSE
    )
)
