# Shared paths and deterministic analysis settings.
required_env <- function(name) {
    value <- Sys.getenv(name)
    if (!nzchar(value))
        stop("Required environment variable is unset: ", name)
    value
}
analysis_seed <- as.integer(Sys.getenv("ANALYSIS_SEED", "20260819"))
set.seed(analysis_seed)
base_dir <- required_env("FEATURECOUNTS_DIR")
count_out_dir <- required_env("RNA_OUT_DIR")
input_rds <- file.path(required_env("OUT_DIR"), "inputs", "rna_inputs.rds")
analysis_rds <- file.path(count_out_dir, "rna_analysis.rds")
dir.create(count_out_dir, recursive = TRUE, showWarnings = FALSE)
