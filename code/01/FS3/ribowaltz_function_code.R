##-------Read length distribution  
rlength_distr <- function(data, sample, multisamples = "average", plot_style = "split", scale_factors = "auto",
    transcripts = NULL, length_range = NULL, cl = 100, colour = NULL) {

    # required packages
    if (!requireNamespace("data.table", quietly = TRUE))
        stop("package 'data.table' is required")
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("package 'ggplot2' is required")
    if (!requireNamespace("scales", quietly = TRUE))
        stop("package 'scales' is required")

    library(data.table)
    library(ggplot2)
    library(scales)

    if (class(data[[1]])[1] == "GRanges") {
        data_tmp <- list()
        for (i in unlist(sample)) {
            data_tmp[[i]] <- as.data.table(data[[i]])[, c("width", "strand") := NULL][,
                seqnames := as.character(seqnames)]
            setnames(data_tmp[[i]], c("seqnames", "start", "end"), c("transcript", "end5",
                "end3"))
        }
        data <- data_tmp
    }

    check_sample <- setdiff(unlist(sample), names(data))
    if (length(check_sample) != 0) {
        cat("\n")
        stop(sprintf("incorrect sample name(s): \"%s\" not found\n\n", paste(check_sample,
            collapse = ", ")))
    }

    if (length(sample) == 0) {
        cat("\n")
        stop("at least one sample name must be spcified\n\n")
    }

    if (is.numeric(scale_factors)) {
        if (!all(unlist(sample) %in% names(scale_factors))) {
            cat("\n")
            stop("scale factor for one or more sample is missing\n\n")
        }
    }

    if (!(multisamples %in% c("average", "independent"))) {
        cat("\n")
        warning("parameter multisamples must be either \"average\" or \"independent\".
            Set to default \"average\"\n",
            call. = FALSE)
        multisamples <- "average"
    }

    if (multisamples == "independent" & is.list(sample)) {
        cat("\n")
        warning("parameter multisamples is set to \"independent\" but parameter sample is a list:
            parameter multisamples will be coerced to default \"average\"\n",
            call. = FALSE)
        multisamples <- "average"
    }

    if (length(sample) != 2 & plot_style == "mirror") {
        cat("\n")
        warning("parameter plot_style is set to \"mirror\" but parameter sample is a list of dimension > 2:
            parameter plot_style will be coerced to default \"split\"\n",
            call. = FALSE)
        plot_style <- "split"
    }

    if (is.character(sample) & length(sample) == 1) {
        multisamples <- "independent"
        plot_style <- "split"
    }

    if (is.character(sample) & length(sample) > 1 & multisamples == "average") {
        sample <- list(Average = sample)
        plot_style <- "split"
        cat("\n")
        warning("Default name of averaged samples is \"Average\":
            consider to use a named list of one element to provide a meaningful plot title\n",
            call. = FALSE)
    }

    if (is.list(sample) & length(sample) == 1) {
        plot_style <- "split"
    }

    if (!(plot_style %in% c("mirror", "dodge", "split", "facet", "heatmap"))) {
        cat("\n")
        warning("parameter plot_style must be either \"split\", \"facet\", \"dodge\", \"mirror\" or \"heatmap\".
            Set to default \"split\"\n",
            call. = FALSE)
        plot_style <- "split"
    }

    #define color vector
    if ((length(colour) < length(sample)) & ((plot_style %in% c("dodge", "mirror")) | (plot_style %in%
        c("split", "facet") & length(colour) != 1))) {

        if (length(colour) != 0) {
            warning("Not enough colors specified:\n
            default ggplot color palette will be used\n",
                call. = FALSE)
        }

        default_gg_col <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]
        }
        colour <- default_gg_col(length(sample))

    } else {
        if (plot_style %in% c("split", "facet") & length(colour) == 1) {
            colour <- rep(colour, length(sample))
        }
    }

    #define length range taking into account all (unlisted) samples
    if (length(length_range) == 0) {
        for (samp in as.character(unlist(sample))) {
            if (length(transcripts) == 0) {
                dt <- data[[samp]]
            } else {
                dt <- data[[samp]][transcript %in% transcripts]
            }

            if (length(length_range) == 0) {
                length_range <- seq(quantile(dt$length, (1 - cl/100)/2), quantile(dt$length,
                  1 - (1 - cl/100)/2))
            } else {
                xmin <- min(min(length_range), quantile(dt$length, (1 - cl/100)/2))
                xmax <- max(max(length_range), quantile(dt$length, 1 - (1 - cl/100)/2))
                length_range <- seq(xmin, xmax)
            }
        }
    }

    xmin = min(length_range)
    xmax = max(length_range)

    # check if all samples have reads of the specified lengths
    # especially required if only one read length is passed
    if (length(length_range) != 0) {
        if (is.list(sample)) {
            samp_dt <- data.table(stack(sample))
            setnames(samp_dt, c("sample", "sample_l"))
        } else {
            samp_dt <- data.table(sample = sample, sample_l = sample)
        }

        for (samp in samp_dt$sample) {

            dt <- data[[samp]][cds_start != 0 & cds_stop != 0]

            if (length(transcripts) != 0) {
                dt <- dt[transcript %in% transcripts]
            }

            len_check <- unique(dt$length)
            if (sum(length_range %in% len_check) == 0) {
                cat("\n")
                warning(sprintf("\"%s\" doesn't contain any reads of the selected lengths: sample removed\n",
                  samp), call. = FALSE)
                sel_l_samp <- samp_dt[sample == samp, sample_l]
                if (is.list(sample)) {
                  sample[[sel_l_samp]] <- sample[[sel_l_samp]][sample[[sel_l_samp]] != samp]
                } else {
                  sample <- sample[sample != samp]
                }
            }
        }
    }

    if (is.null(unlist(sample))) {
        cat("\n")
        stop("none of the data tables listed in sample contains any reads of the specified lengths\n\n")
    }

    # compute count of reads of defined lengths and scale them 
    final_dt <- data.table()
    for (samp in as.character(unlist(sample))) {
        if (length(transcripts) == 0) {
            dt <- data[[samp]]
        } else {
            dt <- data[[samp]][transcript %in% transcripts]
        }

        setkey(dt, length)
        dist_dt <- dt[CJ(length_range), .(count = .N), by = .EACHI]

        #scaling/normalization
        if (is.character(scale_factors) & scale_factors[1] == "auto") {
            dist_dt[, scaled_count := (count/sum(count)) * 100]
            y_title <- "% reads"
        } else {
            y_title <- "# reads"
            if (is.numeric(scale_factors)) {
                dist_dt[, scaled_count := count * scale_factors[samp]]
            } else {
                dist_dt[, scaled_count := count]
            }
        }

        dist_dt[, tmp_samp := samp]

        final_dt <- rbind(final_dt, dist_dt)
    }

    output <- list()
    output[["count_dt"]] <- copy(final_dt[, c("tmp_samp", "length", "count", "scaled_count")])
    if (is.character(scale_factors) & scale_factors[1] == "auto") {
        output[["count_dt"]][, scaled_count := NULL]
    }
    setnames(output[["count_dt"]], "tmp_samp", "sample")

    # compute mean and se of samples if a list is provided
    if (is.list(sample)) {

        samp_dt <- data.table(stack(sample))
        setnames(samp_dt, c("tmp_samp", "sample"))

        # set names of samples as specified in parameter sample  
        final_dt <- merge.data.table(final_dt, samp_dt, sort = FALSE)[, tmp_samp := NULL]

        # compute mean and se
        plot_dt <- final_dt[, .(mean_scaled_count = mean(scaled_count), se_scaled_count = sd(scaled_count/sqrt(.N))),
            by = .(length, sample)]

        if (any(lengths(sample) != 1)) {
            output[["plot_dt"]] <- copy(plot_dt[, c("sample", "length", "mean_scaled_count",
                "se_scaled_count")])
            setnames(output[["plot_dt"]], c("length", "mean_scaled_count", "se_scaled_count"),
                c("x", "y", "y_se"))
        } else {
            output[["plot_dt"]] <- copy(final_dt[, c("sample", "length", "scaled_count")])
            setnames(output[["plot_dt"]], c("length", "scaled_count"), c("x", "y"))
        }

    } else {
        plot_dt <- final_dt[, sample := tmp_samp][, se_scaled_count := NA]
        setnames(plot_dt, "scaled_count", "mean_scaled_count")

        output[["plot_dt"]] <- copy(plot_dt[, c("sample", "length", "mean_scaled_count")])
        setnames(output[["plot_dt"]], c("length", "mean_scaled_count"), c("x", "y"))
    }

    plot_dt[, sample := factor(sample, levels = unique(sample))]

    oldw <- getOption("warn")
    options(warn = -1)

    # ============================
    # Plotting
    # ============================
    if (plot_style == "heatmap") {

        hm_dt <- copy(plot_dt)
        hm_dt[, length := as.numeric(length)]
        hm_dt[, sample := as.character(sample)]
        hm_dt[, value := mean_scaled_count]
        hm_dt[, sample := factor(sample, levels = unique(sample))]
        hm_dt[, length_f := factor(length, levels = sort(unique(length)))]

        output[["plot"]] <- ggplot(hm_dt, aes(x = length_f, y = sample, fill = value)) +
            geom_tile(color = "white", linewidth = 0.2) + labs(x = "Read length", y = NULL,
            fill = y_title) + theme_bw(base_size = 18) + theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text.y = element_text(size = 11)) +
            scale_fill_gradientn(colours = c("#2166ac", "#f7f7f7", "#b2182b"), values = rescale(c(min(hm_dt$value,
                na.rm = TRUE), median(hm_dt$value, na.rm = TRUE), max(hm_dt$value, na.rm = TRUE))))

        output[["heatmap_dt"]] <- hm_dt[, .(sample, length, value)]

    } else if (plot_style == "split") {

        i <- 0
        for (samp in unique(plot_dt$sample)) {
            # generate a plot for each sample and store it
            i <- i + 1
            sel_col = colour[i]

            sub_plot_dt <- plot_dt[sample == samp]
            plot <- ggplot(sub_plot_dt, aes(as.numeric(length), mean_scaled_count)) + geom_bar(stat = "identity",
                fill = sel_col) + geom_errorbar(aes(ymin = mean_scaled_count - se_scaled_count,
                ymax = mean_scaled_count + se_scaled_count), width = 0.35, linewidth = 1.1,
                na.rm = TRUE, color = sel_col, show.legend = FALSE) + labs(title = samp,
                x = "Read length", y = y_title) + theme_bw(base_size = 23) + scale_x_continuous(limits = c(xmin -
                0.5, xmax + 0.5), breaks = seq(xmin + ((xmin)%%2), xmax, by = max(c(1, floor((xmax -
                xmin)/7))))) + theme(panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5))

            output[[paste0("plot_", samp)]] <- plot
        }

    } else {
        if (plot_style == "mirror") {
            plot_dt[sample == unique(plot_dt$sample)[2], mean_scaled_count := -mean_scaled_count]
        }

        plot <- ggplot(plot_dt, aes(as.numeric(length), mean_scaled_count, fill = sample))

        if (plot_style %in% c("facet", "mirror")) {
            plot <- plot + geom_bar(stat = "identity") + geom_errorbar(aes(ymin = mean_scaled_count -
                se_scaled_count, ymax = mean_scaled_count + se_scaled_count, color = sample),
                width = 0.35, linewidth = 1.1, na.rm = TRUE)
            if (identical(plot_style, "mirror")) {
                plot <- plot + geom_hline(yintercept = 0, linetype = 2, color = "gray20")
            }
        } else {
            plot <- plot + geom_bar(stat = "identity", position = position_dodge(0.9)) +
                geom_errorbar(aes(ymin = mean_scaled_count - se_scaled_count, ymax = mean_scaled_count +
                  se_scaled_count, color = sample), width = 0.35, linewidth = 1.1, na.rm = TRUE,
                  position = position_dodge(0.9))
        }

        plot <- plot + labs(x = "Read length", y = y_title) + theme_bw(base_size = 23) +
            scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5), breaks = seq(xmin + ((xmin)%%2),
                xmax, by = max(c(1, floor((xmax - xmin)/7))))) + theme(panel.grid.minor.x = element_blank()) +
            scale_fill_manual(values = colour) + scale_color_manual(values = colour) + scale_y_continuous(labels = abs)

        if (data.table::uniqueN(colour) > 1 & plot_style != "facet") {
            plot <- plot + theme(legend.position = c(0.98, 1), legend.justification = c(1,
                1), legend.title = element_blank(), legend.background = element_blank())
        } else {
            plot <- plot + theme(legend.position = "none")
        }

        if (plot_style == "facet") {
            plot <- plot + facet_wrap(sample ~ ., ncol = ceiling(sqrt(length(sample)))) +
                theme(strip.background = element_blank())
        }

        output[["plot"]] <- plot
    }

    options(warn = oldw)
    return(output)
}

BiocManager::install("GenomicFeatures")
