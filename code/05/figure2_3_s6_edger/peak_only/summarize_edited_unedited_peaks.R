#!/usr/bin/env Rscript

# summarize_edgeR_v2_edited_unedited_Acodons_codonmap_v2.R
# Add codon sequence back from *.fasta.edited.tab.split
# because edgeR output does not contain codon column.

suppressPackageStartupMessages(library(data.table))
fread <- function(file, ...) {
    if (is.character(file) && length(file) == 1L && grepl("[.]gz$", file)) {
        data.table::fread(cmd = paste("gzip -dc --", shQuote(file)), ...)
    } else data.table::fread(file, ...)
}

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(x, default = NULL) {
    i <- match(x, args)
    if (is.na(i))
        return(default)
    args[i + 1]
}

result_dir <- get_arg("--result-dir", "edgeR_stalling_rigorous_v2/results_txmedian_any_both_final")

manifest_path <- get_arg("--manifest", "stalling_master_tables/master_table_manifest.tsv")

codon_dir <- get_arg("--codon-dir", "/path/to/RNA_editing_DATAprocessing_new/ref_genome")

editing_dir <- get_arg("--editing-dir", "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed")

run_label <- get_arg("--run-label", "rigorous_v2_txmedian_any_both")

outdir <- get_arg("--output-dir", file.path(result_dir, "edited_unedited_summary_codonA"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species_codons <- list(dmel = "dmel-all-CDS-r6.12.fasta.edited.tab.split", dere = "dere-all-CDS-r1.05.fasta.edited.tab.split",
    dsec = "dsec-all-CDS-r1.3.fasta.edited.tab.split", dsim = "dsim-all-CDS-r2.02.fasta.edited.tab.split")

load_codon <- function(spec) {
    x <- fread(file.path(codon_dir, species_codons[[spec]]), header = FALSE, col.names = c("site",
        "codon"))
    x[, transcript := sub("_[0-9]+$", "", site)]
    x[, codon_index_1based := as.integer(sub(".*_", "", site))]
    x[, .(transcript, codon_index_1based, codon = toupper(codon))]
}

codon_maps <- lapply(names(species_codons), load_codon)
names(codon_maps) <- names(species_codons)


species_edit <- list(dmel = "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    dere = "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    dsec = "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv",
    dsim = "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv")

load_edit_sites <- function(spec) {

    f <- file.path(editing_dir, species_edit[[spec]])

    x <- fread(f, select = c("new_transcript", "new_cds_pos_1based", "target_region"))

    x <- x[target_region == "CDS"]

    x[, codon_index_1based := 1L + floor((new_cds_pos_1based - 1)/3)]

    x[, site_id := paste0(new_transcript, "_", codon_index_1based)]

    unique(x$site_id)
}

editing_maps <- lapply(names(species_edit), load_edit_sites)
names(editing_maps) <- names(species_edit)


manifest <- fread(manifest_path)

files <- list.files(result_dir, pattern = paste0("\\.flank30\\.", run_label, "\\.edgeR\\.tsv\\.gz$"),
    full.names = TRUE)


run_one <- function(f) {

    group <- sub("\\.flank30\\..*\\.edgeR\\.tsv\\.gz$", "", basename(f))

    m <- manifest[group_id == group][1]

    species <- as.character(m$species)

    e <- fread(f, select = c("site_id", "transcript", "codon_index_1based", "stalling"))

    e <- merge(e, codon_maps[[species]], by = c("transcript", "codon_index_1based"), all.x = TRUE)

    if (anyNA(e$codon))
        stop(group, ": missing codon=", sum(is.na(e$codon)))

    e[, edited := site_id %in% editing_maps[[species]]]
    e[, codon_has_A := grepl("A", codon, fixed = TRUE)]

    e <- e[edited | codon_has_A]

    ans <- rbindlist(lapply("stalling", function(ev) {

        call <- e[[ev]]

        ei <- sum(e$edited & call)
        eo <- sum(e$edited) - ei
        ui <- sum(!e$edited & call)
        uo <- sum(!e$edited) - ui

        ft <- fisher.test(matrix(c(ei, eo, ui, uo), nrow = 2), alternative = "greater")

        data.table(group_id = group, species = species, genotype = m$genotype, sex = m$sex,
            event = ev, n_tested = nrow(e), edited_total = sum(e$edited), edited_event = ei,
            edited_rate = ei/sum(e$edited), unedited_total = sum(!e$edited), unedited_event = ui,
            unedited_rate = ui/sum(!e$edited), odds_ratio = unname(ft$estimate), fisher_p = ft$p.value)
    }))

    ans
}


result <- rbindlist(lapply(files, run_one), fill = TRUE)

fwrite(result, file.path(outdir, "edited_unedited_peaks_by_library.tsv"), sep = "\t")

message("Finished: ", outdir)
