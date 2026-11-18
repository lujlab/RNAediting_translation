#!/usr/bin/env Rscript
# A-site I>A vs remaining testable sites, using the supplied method-specific Table1.
# The plot compares N/S composition using a separate two-sided Fisher test.
a <- commandArgs(TRUE)
stopifnot(length(a) == 2L)
x <- read.delim(a[1], check.names = FALSE)
order <- paste(rep(c("mel", "sim", "sec", "ere"), each = 2), rep(c("F", "M"), 4), sep = "_")
x <- x[match(order, x$sample), ]
parse_ns <- function(v) {
    z <- regmatches(v, regexec("^([0-9]+) \\(([0-9]+)/([0-9]+)\\)$", v))[[1]]
    if (length(z) != 4)
        stop("Invalid N/S cell: ", v)
    n <- as.integer(z[2:4])
    if (n[1] != sum(n[2:3]))
        stop("Unknown N/S classes: ", v)
    n[2:3]
}
z <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    total <- parse_ns(x$A_Total_sites_N_S[i])
    ia <- parse_ns(x$A_I_gt_A_N_S[i])
    rem <- total - ia
    stopifnot(all(rem >= 0), ia[2] > 0, rem[2] > 0)
    p <- fisher.test(rbind(ia, rem), alternative = "two.sided")$p.value
    data.frame(sample = x$sample[i], IgtA_N = ia[1], IgtA_S = ia[2], Remaining_N = rem[1],
        Remaining_S = rem[2], IgtA_NS = ia[1]/ia[2], Remaining_NS = rem[1]/rem[2], Fisher_P = p,
        P_label = sprintf(" = %.3f", p))
}))
write.table(z, a[2], sep = "\t", row.names = FALSE, quote = FALSE)
