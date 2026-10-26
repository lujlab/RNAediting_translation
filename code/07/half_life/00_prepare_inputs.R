source("config.R")
source("R/utils.R")
make_output_dirs(out_dir)
dir_create("inputs/published")
# Public input tables are bundled; never download or silently substitute them.
urls <- c(`2010_ArrayExpress_processed` = "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/580/E-MEXP-2580/Files/Thomsenetal.normalised.classification_orig.txt",
    `2023_PLOS_S1` = "https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3001956.s002",
    `2023_PLOS_S2` = "https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3001956.s003")
dests <- c(half_life_2010_file, expression_2023_file, half_life_2023_file)
missing_tables <- dests[!file.exists(dests)]
if (length(missing_tables)) stop("Missing bundled public tables: ", paste(missing_tables,
    collapse = ", "), ". Restore inputs/published from the release package.")
if (any(file.info(dests)$size <= 0)) stop("A bundled public input table is empty.")
download_status <- rep("bundled_local_file", length(dests))
input_path <- c(editing_file, or_default(editing_expression_file, "NOT_PROVIDED"), or_default(callable_gene_file,
    "NOT_PROVIDED"), gtf_file, transcript_fasta, half_life_2010_file, expression_2023_file,
    half_life_2023_file)
write_tsv(data.table(input = c("editing_file", "editing_expression_file_B2", "callable_gene_file_B3",
    "gtf_file", "transcript_fasta", "half_life_2010_file", "expression_2023_file", "half_life_2023_file"),
    path = input_path, exists = file.exists(input_path)), file.path(out_dir, "00_input_QC",
    "input_manifest.tsv"))
write_tsv(data.table(dataset = names(urls), url = unname(urls), destination = dests, status = download_status),
    file.path(out_dir, "00_input_QC", "download_status.tsv"))
message("Input preparation complete. 2010 table present: ", file.exists(half_life_2010_file),
    "; 2023 table present: ", file.exists(half_life_2023_file))
