# Requires data.table, jsonlite and PyMuPDF.
s18_existing_reference <- function(filename, configured, variable, base_dir) {
    explicit <- Sys.getenv(variable, unset = "")
    if (nzchar(explicit)) {
        if (!file.exists(explicit)) stop(variable, " does not exist: ", explicit)
        return(normalizePath(explicit))
    }
    if (length(configured) && file.exists(configured)) return(normalizePath(configured))
    current <- normalizePath(base_dir)
    for (i in seq_len(10L)) {
        candidate <- file.path(current, "ref_genome", filename)
        if (file.exists(candidate)) return(normalizePath(candidate))
        parent <- dirname(current)
        if (identical(parent, current)) break
        current <- parent
    }
    stop("Cannot locate ", filename, "; set ", variable, ".")
}

s18_command <- function(command, args) {
    status <- system2(command, vapply(args, shQuote, character(1)))
    if (!identical(as.integer(status), 0L)) stop("Command failed: ", command)
    invisible(status)
}

s18_render <- function(base_dir, data_dir, output_pdf) {
    python <- Sys.getenv("PDF_PYTHON", unset = Sys.getenv("PYTHON", unset = "python3"))
    s18_command(python, c(file.path(base_dir, "scripts", "figure_s18", "draw_editable_figure.py"),
        "--data-dir", normalizePath(data_dir), "--output-pdf", output_pdf))
}

s18_run <- function(base_dir) {
    base_dir <- normalizePath(base_dir)
    config <- new.env(parent = globalenv())
    previous <- setwd(base_dir)
    on.exit(setwd(previous), add = TRUE)
    sys.source(file.path(base_dir, "config.R"), envir = config)
    if (config$min_editing_sites_per_gene != 1L)
        stop("This Figure S18 uses at least one CDS editing site per measurement unit.")
    for (package in c("data.table", "jsonlite")) {
        if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
    }
    gtf <- s18_existing_reference("dmel-all-r6.12.gtf", config$gtf_file,
        "HALF_LIFE_GTF", base_dir)
    tx <- s18_existing_reference("dmel-all-transcript-r6.12.fasta", config$transcript_fasta,
        "HALF_LIFE_TRANSCRIPT_FASTA", base_dir)
    genome <- s18_existing_reference("dmel-all-chromosome-r6.12.fasta", "",
        "HALF_LIFE_GENOME_FASTA", base_dir)
    data_dir <- config$out_dir
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    data_dir <- normalizePath(data_dir)
    default_pdf <- file.path(base_dir, "..", "..", "..", "output",
        "Supplementary_Figure_18_r6_12_transcripts_editable.pdf")
    output_pdf <- Sys.getenv("HALF_LIFE_PDF", unset = default_pdf)
    output_pdf <- file.path(normalizePath(dirname(output_pdf), mustWork = FALSE), basename(output_pdf))
    helper <- file.path(base_dir, "scripts", "figure_s18")
    python <- Sys.getenv("PYTHON", unset = "python3")
    s18_command(python, c(file.path(helper, "build_data.py"),
        "--output-dir", data_dir, "--gtf", gtf, "--genome-fasta", genome,
        "--transcript-fasta", tx, "--editing-file", normalizePath(config$editing_file),
        "--half-life-file", normalizePath(config$half_life_2010_file)))
    s18_command(file.path(R.home("bin"), "Rscript"), c("--vanilla",
        file.path(helper, "calculate_statistics.R"), data_dir))
    s18_render(base_dir, data_dir, output_pdf)
    cat("[DONE]", output_pdf, "\n")
    invisible(list(data_dir = data_dir, output_pdf = output_pdf))
}
