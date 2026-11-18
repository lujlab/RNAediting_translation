# Central configuration for the Drosophila RNA-editing / mRNA half-life analysis.
# HALF_LIFE_* variables configure paths. Figure S18 requires a minimum of one CDS editing site.
set.seed(20260827)

project_dir <- normalizePath(getwd(), mustWork = TRUE)
editing_file <- Sys.getenv("HALF_LIFE_EDITING_FILE", unset = "inputs/published/Duan2017_Graveley2011_embryo_CDS_editing_sites_long.tsv")
# Optional editing-experiment files define background sets only.
# Regression expression covariates come from the half-life study.
editing_expression_file <- NULL  # B2: gene plus normalized_count/TPM/FPKM
editing_expression_threshold <- 1
callable_gene_file <- NULL  # B3: gene plus optional callable flag
gtf_file <- Sys.getenv("HALF_LIFE_GTF", unset = "/path/to/RNA_editing_DATAprocessing_new/ref_genome/dmel-all-r6.12.gtf")
transcript_fasta <- Sys.getenv("HALF_LIFE_TRANSCRIPT_FASTA", unset = "/path/to/RNA_editing_DATAprocessing_new/ref_genome/dmel-all-transcript-r6.12.fasta")

# Thomsen et al. (2010) normalized probe-level expression time series and decay classes.
half_life_2010_file <- "inputs/published/Thomsenetal.normalised.classification_orig.txt"
# PLOS Biology 2023 S1/S2 Tables are bundled for offline preparation.
expression_2023_file <- "inputs/published/2023_S1_Table.xlsx"
half_life_2023_file <- "inputs/published/2023_S2_Table.xlsx"

out_dir <- Sys.getenv("HALF_LIFE_OUT_DIR", unset = "output/duan2017_embryo_cds_min1")
editing_min_depth <- 0L
editing_min_ratio <- 0
# Editing-site threshold for exploratory edited and unedited group comparisons.
# Genes with detected sites below the threshold are excluded.
min_editing_sites_per_gene <- as.integer(Sys.getenv("HALF_LIFE_MIN_EDITING_SITES", unset = "1"))
if (is.na(min_editing_sites_per_gene) || min_editing_sites_per_gene < 1L) {
    stop("HALF_LIFE_MIN_EDITING_SITES must be a positive integer")
}
max_matching_controls <- 5L
min_match_controls <- 1L
top_outlier_fraction <- 0.01
write_png <- FALSE

required_columns_2010 <- c("gene", "decay_class", "half_life_min")
required_columns_2023 <- c("transcript_id", "half_life_min")
