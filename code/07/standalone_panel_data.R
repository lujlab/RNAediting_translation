# Calculate F4D and S15D-G from experimental measurements.
out <- Sys.getenv("REPORTER_WORK_DIR", "results/.reporter_work")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
raw <- read.delim("data/07/dual_plot.txt", check.names = FALSE)
dose <- as.numeric(sub(" mM", "", raw$Inosine_concentration, fixed = TRUE))
groups <- list(1:3, 4:6, 7:9)
reps <- do.call(rbind, lapply(seq_along(groups), function(i) {
    data.frame(Panel = "F4D", Concentration = dose, Replicate = paste0("rep", i), Value = rowMeans(raw[paste0("well_",
        groups[[i]])]), Source = paste(paste0("well_", groups[[i]]), collapse = ";"))
}))
wb <- readRDS("data/07/Westernblot_and_luciferase.rds")
configs <- list(FS15D = c(4, 3), FS15E = c(2, 2), FS15F = c(3, 3), FS15G = c(1, 2))
order <- c("c_0mM", "c_0.1mM", "c_1mM", "c_2.5mM", "c_5mM", "c_10mM")
for (panel in names(configs)) {
    config <- configs[[panel]]
    d <- data.frame(wb[config[1]])[order, , drop = FALSE]
    stopifnot(ncol(d) == config[2] + 2, !anyNA(d))
    for (i in seq_len(config[2])) reps <- rbind(reps, data.frame(Panel = panel, Concentration = c(0,
        0.1, 1, 2.5, 5, 10), Replicate = paste0("rep", i), Value = d[[i]], Source = names(d)[i]))
}
summaries <- do.call(rbind, lapply(split(reps, reps$Panel), function(d) {
    do.call(rbind, lapply(sort(unique(d$Concentration)), function(x) {
        v <- d$Value[d$Concentration == x]
        data.frame(Panel = d$Panel[1], Concentration = x, N = length(v), Mean = mean(v),
            SD = sd(v))
    }))
}))
tests <- do.call(rbind, lapply(split(summaries, summaries$Panel), function(d) {
    r <- cor.test(d$Concentration, d$Mean, method = "spearman", alternative = "less", exact = TRUE)
    # With distinct, strictly decreasing means, exact lower-tail P = 1 / n!.
    stopifnot(unname(r$estimate) == -1, all(diff(d$Mean) < 0))
    data.frame(Panel = d$Panel[1], rho = unname(r$estimate), P = 1/factorial(nrow(d)), N_doses = nrow(d),
        Alternative = "less", Test = "Exact Spearman", R_P = r$p.value)
}))
stopifnot(all(summaries$N[summaries$Panel == "F4D"] == 3), abs(tests$P[tests$Panel == "F4D"] -
    2.50521083854417e-08) < 1e-20)
for (item in c("reps", "summaries", "tests")) write.table(get(item), file.path(out, paste0(item,
    ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
invisible(NULL)
