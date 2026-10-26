library(data.table)
library(dplyr)

# specify path
input_file <- commandArgs(trailingOnly = TRUE)[1]

# load data
data <- fread(input_file, sep = "\t")

# rename colnames
setnames(data, c("chr", "pos", "Ref", "Alt", "refcount", "altcount", "MQ", "type", "gene",
    "change", "ann"))

# calculate depth
data[, dp := refcount + altcount]

# filter
data <- data[nchar(Ref) == 1 & nchar(Alt) == 1]
data <- data[dp > 1]
data[, handle := paste(chr, pos, sep = "_")]

# define order for sites annotation
data[, class2 := 10]
data[grep("missense_variant", type), class2 := 1]
data[grep("synonymous_variant", type), class2 := 2]
data[grep("prime_UTR_variant", type), class2 := 3]
data[grep("non_coding_transcript_exon_variant", type), class2 := 4]
data[grep("intragenic_variant", type), class2 := 5]

# define strand
data[, change2 := gsub(".*[0-9]", "", change)]
data[, change := paste(Ref, Alt, sep = ">")]

data[, strand := "+"]
data[change != change2, strand := "-"]

data[, class := 1]
data[grep("A>G", change2), class := 0]

# sort sites
data <- data[order(handle, class, class2)]

# add 'ctl'
data[, ctl := "NA"]
data$ctl[2:nrow(data)] <- as.character(data$handle[1:(nrow(data) - 1)])

# filter
data <- data[ctl != handle]

# define P-value(probability that an observed editing site is derived from sequencing error)
eps0 <- 0.01

perr0 <- function(x) {
    p <- pbinom(x[2] - 1, x[1], prob = eps0, lower.tail = F)
    return(p)
}
data$pvalue <- apply(data[, c("dp", "altcount")], 1, function(x) {
    x = c(x[1], x[2])
    perr0(x)
})
data$pvalue <- as.numeric(data$pvalue)
data$fdr <- p.adjust(p = data$pvalue, method = "fdr")

# filter sites
filtered_data <- data[data$dp >= 5 & data$fdr < 0.05 & data$altcount/(data$refcount + data$altcount) >
    0.01 & data$MQ >= 10 & data$altcount >= 2]

# statistical results
table1 <- table(unique(filtered_data[, c("handle", "type", "change2")])[, c(2, 3)])
table2 <- table(unique(filtered_data[, c("handle", "change2")])[, c(2)])
#result <- AG[, .(min = min(ratio), max = max(ratio), mean = mean(ratio), median = median(ratio)), by = class2]

# output
AG <- filtered_data[change2 == "A>G"]
AG[, ratio := altcount/dp]
output_file_AG <- paste0(sub("\\.[^.]+$", "", input_file), "_AG_dp5mq10fdr0.05ratio0.01at2")
fwrite(AG, file = output_file_AG, sep = "\t", quote = FALSE)


AG2 <- data[dp >= 5 & MQ >= 10 & change2 == "A>G"]
output_file_AG2 <- paste0(sub("\\.[^.]+$", "", input_file), "_AG_dp5mq10")
fwrite(AG2, file = output_file_AG2, sep = "\t", quote = FALSE)
