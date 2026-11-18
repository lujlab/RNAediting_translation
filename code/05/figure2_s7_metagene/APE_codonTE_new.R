#!/usr/bin/env Rscript
# Historical single-library gmRC calculation used upstream of the selected PDF.
# Input: headerless 15-column whole-codon counts, in the order specified below.
# Output: tab-separated table on stdout (before codon/AA annotation).
# Retain transcripts with >75 codons having mRNA >=5; calculate gmRC only at
# those codons. The normalization denominator is arithmetic mean(mRC).
# Asite_median/Esite_median/Psite_median are separate metadata, NOT denominators.
args <- commandArgs(trailingOnly = TRUE)
library(data.table)

tmp <- fread(file = args[1], header = FALSE, col.names = c("chr", "pos", "FBtr", "no", "strand",
    "cdsposstart", "codonpos", "mRNA", "mRNA_median", "Asite", "Asite_median", "Esite", "Esite_median",
    "Psite", "Psite_median"))

result <- tmp[tmp[, .I[sum(mRNA >= 5) > 75], by = FBtr]$V1]

result[, `:=`(mRNA_median = median(mRNA[mRNA >= 5]), Asite_median = median(Asite[mRNA >=
    5]), Esite_median = median(Esite[mRNA >= 5]), Psite_median = median(Psite[mRNA >= 5])),
    by = FBtr]

hq <- result[mRNA >= 5]

hq[, mRC := Asite/mRNA]
hq[, gmRC := mRC/mean(mRC, na.rm = TRUE), by = FBtr]

result <- merge(result, hq[, .(FBtr, pos, gmRC)], by = c("FBtr", "pos"), all.x = TRUE)

fwrite(result, file = "", quote = FALSE, sep = "\t", col.names = TRUE)
