#!/usr/bin/env Rscript
# Publish the existing edgeR table as the manuscript Table 1.
# Fisher statistics stay in the comparison outputs and overlap_2methods plots.
a <- commandArgs(TRUE)
if (length(a) != 1L) stop("Usage: export_edgeR_manuscript_Table1.R ANALYSIS_OUTPUT_DIR")
out <- a[[1]]
source <- file.path(out, "Table1_edgeR_final.tsv")
header <- file.path(out, "Table1_edgeR_final_two_row_header.tsv")
if (!file.exists(source) || !file.exists(header)) stop("Missing edgeR final tables: ", out)
x <- read.delim(source, check.names = FALSE)
expected <- c("sample", "A_Total_sites_N_S", "A_I_gt_A_N_S", "A_I_lt_A", "AP_Total_sites_N_S",
    "AP_I_gt_A_N_S", "AP_I_lt_A")
if (!identical(names(x), expected) || nrow(x) != 12L || anyDuplicated(x$sample)) stop("Unexpected edgeR Table 1 schema or groups")
if (!file.copy(source, file.path(out, "Table1_manuscript.tsv"), overwrite = TRUE)) stop("Cannot export Table 1")
if (!file.copy(header, file.path(out, "Table1_manuscript_two_row_header.tsv"), overwrite = TRUE)) stop("Cannot export display table")
write.table(data.frame(role = c("manuscript_Table1", "comparison_only"), method = c("edgeR_QL",
    "Fisher_directional"), source = c("Table1_edgeR_final.tsv", "overlap_2methods.*")), file.path(out,
    "Table1_method_roles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Manuscript Table 1 exported from edgeR; Fisher retained for comparison only.")
