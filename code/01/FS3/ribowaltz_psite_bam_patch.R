# ribowaltz_psite_bam_patch.R
# Convert psite-tagged transcriptome BAM files to riboWaltz-like data.table objects
# and plot P-site metaprofile / frame percentage with the same count logic as riboWaltz.
#
# Required R packages: data.table, ggplot2
# Required command-line program: samtools in PATH

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

.default_gg_col <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

.prepare_annotation_for_psite_bam <- function(annotation) {
    ann <- as.data.table(annotation)
    ann <- copy(ann)

    req <- c("transcript", "l_utr5", "l_cds", "l_utr3")
    miss <- setdiff(req, names(ann))
    if (length(miss) > 0) {
        stop(sprintf("annotation is missing column(s): %s. Use riboWaltz::create_annotation() output or a table with transcript/l_utr5/l_cds/l_utr3.",
            paste(miss, collapse = ", ")))
    }

    ann[, transcript := as.character(transcript)]
    ann[, `:=`(l_utr5 = as.integer(l_utr5), l_cds = as.integer(l_cds), l_utr3 = as.integer(l_utr3))]

    # riboWaltz-style coordinates on transcriptome:
    # start codon first nucleotide has psite_from_start == 0;
    # stop codon first nucleotide has psite_from_stop == 1.
    # If cds_start/cds_stop already exist, keep them; otherwise infer from lengths.
    if (!all(c("cds_start", "cds_stop") %in% names(ann))) {
        ann[, `:=`(cds_start = l_utr5 + 1L, cds_stop = l_utr5 + l_cds)]
    } else {
        ann[, `:=`(cds_start = as.integer(cds_start), cds_stop = as.integer(cds_stop))]
    }

    unique(ann[, .(transcript, l_utr5, l_cds, l_utr3, cds_start, cds_stop)])
}

# Read one psite-tagged BAM and return a named list containing one riboWaltz-like data.table.
# BAM optional tags expected:
#   PS:i:<offset>  P-site offset from SAM POS; for your example psite = POS + PS.
#   XF:i:<frame>   frame tag, 0/1/2.
#   NH:i:<n>       optional, used only when unique_only = TRUE.
read_psite_bam_as_ribowaltz <- function(bam, annotation, sample_name = NULL, samtools = "samtools",
    exclude_flag = 4L, unique_only = FALSE, psite_shift = 0L, check_frame = TRUE, showProgress = TRUE) {
    if (is.null(sample_name)) {
        sample_name <- tools::file_path_sans_ext(basename(bam))
    }

    ann <- .prepare_annotation_for_psite_bam(annotation)

    awk <- paste0("awk 'BEGIN{OFS=\"\\t\"} ", "{ps=\"NA\"; xf=\"NA\"; nh=\"NA\"; ", "for(i=12;i<=NF;i++){",
        "if($i ~ /^PS:i:/){split($i,a,\":\"); ps=a[3]} ", "else if($i ~ /^XF:i:/){split($i,b,\":\"); xf=b[3]} ",
        "else if($i ~ /^NH:i:/){split($i,c,\":\"); nh=c[3]}}; ", "print $3,$4,length($10),ps,xf,nh,$2}'")
    cmd <- paste(samtools, "view -F", as.integer(exclude_flag), shQuote(bam), "|", awk)

    dt <- fread(cmd = cmd, col.names = c("transcript", "end5", "length", "psite_offset",
        "frame_bam", "NH", "flag"), na.strings = c("NA", "*"), showProgress = showProgress)

    if (nrow(dt) == 0L) {
        stop("No aligned records were read from BAM after samtools filtering.")
    }

    dt[, transcript := as.character(transcript)]
    int_cols <- c("end5", "length", "psite_offset", "frame_bam", "NH", "flag")
    dt[, (int_cols) := lapply(.SD, as.integer), .SDcols = int_cols]

    n_no_ps <- dt[is.na(psite_offset), .N]
    if (n_no_ps > 0L) {
        warning(sprintf("Dropping %d BAM records without PS:i tag.", n_no_ps), call. = FALSE)
    }
    dt <- dt[!is.na(psite_offset)]

    if (unique_only) {
        n_before <- nrow(dt)
        dt <- dt[!is.na(NH) & NH == 1L]
        message(sprintf("unique_only=TRUE: kept %d/%d records with NH:i:1.", nrow(dt), n_before))
    }

    # In the SAM shown in the question, POS=13 and PS:i:13 gives P-site at transcript coordinate 26.
    # If your frame check shows a systematic one-nt mismatch, try psite_shift = -1L.
    dt[, `:=`(end3 = end5 + length - 1L, psite = end5 + psite_offset + as.integer(psite_shift))]

    n_before_join <- nrow(dt)
    dt <- merge(dt, ann, by = "transcript", all.x = FALSE, sort = FALSE)
    n_dropped <- n_before_join - nrow(dt)
    if (n_dropped > 0L) {
        warning(sprintf("Dropped %d records because their transcript IDs are absent from annotation.",
            n_dropped), call. = FALSE)
    }

    dt[, `:=`(psite_from_start = psite - cds_start, psite_from_stop = psite - cds_stop, psite_region = fifelse(psite <
        cds_start, "5utr", fifelse(psite <= cds_stop, "cds", "3utr")), frame_calc = (psite -
        cds_start)%%3L)]

    if (check_frame && any(!is.na(dt$frame_bam))) {
        n_xf <- dt[!is.na(frame_bam) & frame_bam %in% 0:2, .N]
        n_bad <- dt[!is.na(frame_bam) & frame_bam %in% 0:2 & frame_bam != frame_calc, .N]
        if (n_xf > 0L && n_bad > 0L) {
            warning(sprintf("XF frame differs from computed frame for %d/%d tagged records (%.2f%%). Check psite_shift or annotation convention.",
                n_bad, n_xf, 100 * n_bad/n_xf), call. = FALSE)
        }
    }

    setcolorder(dt, c("transcript", "end5", "end3", "length", "psite", "psite_offset", "frame_bam",
        "frame_calc", "psite_region", "psite_from_start", "psite_from_stop", "cds_start",
        "cds_stop", "l_utr5", "l_cds", "l_utr3", "NH", "flag"))

    out <- list(dt)
    names(out) <- sample_name
    out
}

# Read multiple BAMs and return one riboWaltz-like reads_list.
read_psite_bams_as_ribowaltz <- function(bams, annotation, sample_names = NULL, ...) {
    if (is.null(sample_names)) {
        sample_names <- names(bams)
        if (is.null(sample_names) || any(sample_names == "")) {
            sample_names <- tools::file_path_sans_ext(basename(bams))
        }
    }
    if (length(bams) != length(sample_names)) {
        stop("bams and sample_names must have the same length.")
    }

    out <- vector("list", length(bams))
    names(out) <- sample_names
    for (i in seq_along(bams)) {
        tmp <- read_psite_bam_as_ribowaltz(bam = bams[[i]], annotation = annotation, sample_name = sample_names[[i]],
            ...)
        out[[sample_names[[i]]]] <- tmp[[sample_names[[i]]]]
    }
    out
}

.make_sample_map <- function(data, sample, multisamples = "average") {
    if (length(sample) == 0L)
        stop("at least one sample name must be specified")

    if (is.character(sample)) {
        if (length(sample) == 1L || identical(multisamples, "independent")) {
            sample_map <- data.table(tmp_samp = as.character(sample), sample = as.character(sample))
        } else {
            sample_map <- data.table(tmp_samp = as.character(sample), sample = "Average")
        }
    } else if (is.list(sample)) {
        if (is.null(names(sample)) || any(names(sample) == "")) {
            names(sample) <- paste0("Sample", seq_along(sample))
        }
        sample_map <- rbindlist(lapply(names(sample), function(nm) {
            data.table(tmp_samp = as.character(sample[[nm]]), sample = nm)
        }))
    } else {
        stop("sample must be a character vector or a named list")
    }

    missing <- setdiff(sample_map$tmp_samp, names(data))
    if (length(missing) > 0L) {
        stop(sprintf("incorrect sample name(s): %s not found", paste(missing, collapse = ", ")))
    }
    sample_map
}

.compute_length_range <- function(data, actual_samples, transcripts, length_range = NULL,
    cl = 100) {
    if (length(length_range) != 0L)
        return(as.integer(length_range))
    if (!is.numeric(cl) || length(cl) != 1L || cl <= 0 || cl > 100)
        stop("cl must be in (0, 100]")

    xmin <- Inf
    xmax <- -Inf
    for (samp in actual_samples) {
        dt <- data[[samp]][transcript %in% transcripts]
        if (nrow(dt) == 0L)
            next
        lo <- floor(as.numeric(quantile(dt$length, (1 - cl/100)/2, na.rm = TRUE)))
        hi <- ceiling(as.numeric(quantile(dt$length, 1 - (1 - cl/100)/2, na.rm = TRUE)))
        xmin <- min(xmin, lo)
        xmax <- max(xmax, hi)
    }
    if (!is.finite(xmin) || !is.finite(xmax))
        stop("No reads left for selected transcripts; cannot define length_range.")
    seq(as.integer(xmin), as.integer(xmax))
}

# P-site metaprofile around start and stop codons.
# This mirrors riboWaltz metaprofile_psite counting after the BAM has been converted above.
metaprofile_psite_bam <- function(data, annotation, sample, multisamples = "average", plot_style = "split",
    scale_factors = "auto", transcripts = NULL, length_range = NULL, cl = 100, utr5l = 25,
    cdsl = 50, utr3l = 25, colour = NULL) {
    ann <- .prepare_annotation_for_psite_bam(annotation)
    sample_map <- .make_sample_map(data, sample, multisamples)
    actual_samples <- unique(sample_map$tmp_samp)
    plot_samples <- unique(sample_map$sample)

    if (!plot_style %in% c("split", "facet", "overlap", "mirror"))
        plot_style <- "split"
    if (length(plot_samples) == 1L)
        plot_style <- "split"
    if (identical(plot_style, "mirror") && length(plot_samples) != 2L)
        plot_style <- "split"

    if (length(colour) == 0L)
        colour <- .default_gg_col(length(plot_samples))
    if (length(colour) < length(plot_samples)) {
        colour <- rep(colour, length.out = length(plot_samples))
    } else {
        colour <- colour[seq_along(plot_samples)]
    }
    names(colour) <- plot_samples

    l_transcripts <- as.character(ann[l_utr5 >= utr5l & l_cds > 2 * (cdsl + 1) & l_utr3 >=
        utr3l, transcript])
    c_transcripts <- if (length(transcripts) == 0L)
        l_transcripts else intersect(l_transcripts, transcripts)
    if (length(c_transcripts) == 0L)
        stop("No transcripts passed the metaprofile annotation filters.")

    length_range <- .compute_length_range(data, actual_samples, c_transcripts, length_range,
        cl)

    if (is.numeric(scale_factors) && !all(actual_samples %in% names(scale_factors))) {
        stop("scale factor for one or more samples is missing")
    }

    y_title <- if (is.character(scale_factors) && scale_factors[1] == "auto")
        "P-site frequency" else "# P-sites"

    final_dt <- rbindlist(lapply(actual_samples, function(samp) {
        dt <- data[[samp]][transcript %in% c_transcripts & length %in% length_range]

        start_grid <- data.table(distance = seq(-utr5l, cdsl))
        start_count <- dt[psite_from_start >= -utr5l & psite_from_start <= cdsl, .(count = .N),
            by = .(distance = psite_from_start)]
        start_tab <- merge(start_grid, start_count, by = "distance", all.x = TRUE, sort = FALSE)
        start_tab[is.na(count), count := 0L]
        start_tab[, region := "start"]

        stop_grid <- data.table(distance = seq(-cdsl, utr3l))
        stop_count <- dt[psite_from_stop >= -cdsl & psite_from_stop <= utr3l, .(count = .N),
            by = .(distance = psite_from_stop)]
        stop_tab <- merge(stop_grid, stop_count, by = "distance", all.x = TRUE, sort = FALSE)
        stop_tab[is.na(count), count := 0L]
        stop_tab[, region := "stop"]

        samp_tab <- rbind(start_tab, stop_tab, use.names = TRUE)
        total <- sum(samp_tab$count)
        if (is.character(scale_factors) && scale_factors[1] == "auto") {
            samp_tab[, scaled_count := if (total > 0)
                count/total else 0]
        } else if (is.numeric(scale_factors)) {
            samp_tab[, scaled_count := count * scale_factors[samp]]
        } else {
            samp_tab[, scaled_count := count]
        }
        samp_tab[, tmp_samp := samp]
        samp_tab
    }), use.names = TRUE)

    final_dt[, region := factor(region, levels = c("start", "stop"), labels = c("Distance from start (nt)",
        "Distance from stop (nt)"))]

    output <- list()
    output[["count_dt"]] <- copy(final_dt[, .(sample = tmp_samp, region, distance, count,
        scaled_count)])
    if (is.character(scale_factors) && scale_factors[1] == "none")
        output[["count_dt"]][, scaled_count := NULL]

    plot_dt <- merge(final_dt, sample_map, by = "tmp_samp", sort = FALSE)
    plot_dt <- plot_dt[, .(mean_scaled_count = mean(scaled_count), se_scaled_count = if (.N >
        1L)
        sd(scaled_count)/sqrt(.N) else NA_real_), by = .(sample, region, distance)]

    output[["plot_dt"]] <- copy(plot_dt[, .(sample, region, x = distance, y = mean_scaled_count,
        y_se = se_scaled_count)])
    plot_dt[, sample := factor(sample, levels = plot_samples)]

    lines3nt <- data.table(region = rep(c("Distance from start (nt)", "Distance from stop (nt)"),
        times = c(length(seq(3, cdsl, 3)), length(seq(-2, -cdsl, -3)))), line = c(seq(3,
        cdsl, 3), rev(seq(-2, -cdsl, -3))))
    linered <- data.table(region = c("Distance from start (nt)", "Distance from stop (nt)"),
        line = c(0, 1))

    oldw <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = oldw), add = TRUE)

    if (identical(plot_style, "split")) {
        for (samp in plot_samples) {
            sub_plot_dt <- plot_dt[as.character(sample) == samp]
            sel_col <- colour[[samp]]
            p <- ggplot(sub_plot_dt, aes(distance, mean_scaled_count)) + geom_line(linewidth = 1.5,
                color = sel_col)
            if (!all(is.na(sub_plot_dt$se_scaled_count))) {
                p <- p + geom_ribbon(aes(ymin = mean_scaled_count - se_scaled_count, ymax = mean_scaled_count +
                  se_scaled_count), fill = sel_col, color = NA, alpha = 0.2, show.legend = FALSE)
            }
            p <- p + geom_vline(data = linered, aes(xintercept = line), inherit.aes = FALSE,
                linetype = 1, color = "red") + geom_vline(data = lines3nt, aes(xintercept = line),
                inherit.aes = FALSE, linetype = 3, color = "gray60") + theme_bw(base_size = 30) +
                facet_grid(. ~ region, switch = "x", scales = "free_x") + theme(panel.grid.major.x = element_blank(),
                panel.grid.minor = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
                strip.background = element_blank(), strip.placement = "outside") + labs(title = samp,
                y = y_title)
            output[[paste0("plot_", samp)]] <- p
        }
    } else {
        if (identical(plot_style, "mirror")) {
            mirror_samp <- plot_samples[[2]]
            plot_dt[as.character(sample) == mirror_samp, mean_scaled_count := -mean_scaled_count]
        }

        p <- ggplot(plot_dt, aes(distance, mean_scaled_count, color = sample, fill = sample)) +
            geom_line(linewidth = 1.5) + geom_ribbon(aes(ymin = mean_scaled_count - se_scaled_count,
            ymax = mean_scaled_count + se_scaled_count), color = NA, alpha = 0.2, show.legend = FALSE) +
            labs(y = y_title) + geom_vline(data = linered, aes(xintercept = line), inherit.aes = FALSE,
            linetype = 1, color = "red") + geom_vline(data = lines3nt, aes(xintercept = line),
            inherit.aes = FALSE, linetype = 3, color = "gray60") + theme_bw(base_size = 30) +
            theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                strip.placement = "outside") + scale_fill_manual(values = colour) + scale_color_manual(values = colour) +
            scale_y_continuous(labels = abs)

        if (identical(plot_style, "mirror"))
            p <- p + geom_hline(yintercept = 0, linetype = 2, color = "gray20")
        if (identical(plot_style, "facet")) {
            p <- p + facet_grid(sample ~ region, switch = "x", scales = "free_x")
        } else {
            p <- p + facet_grid(. ~ region, switch = "x", scales = "free_x")
        }
        if (length(unique(colour)) > 1L && !identical(plot_style, "facet")) {
            p <- p + theme(legend.position = c(0.98, 1), legend.justification = c(1, 1),
                legend.title = element_blank(), legend.background = element_blank())
        } else {
            p <- p + theme(legend.position = "none")
        }
        output[["plot"]] <- p
    }

    output
}

# Frame percentage plot.
# frame_source = 'XF' uses the BAM XF:i tag; frame_source = 'computed' mimics riboWaltz exactly.
frame_psite_bam <- function(data, annotation, sample, multisamples = "average", plot_style = "split",
    transcripts = NULL, region = "all", length_range = NULL, cl = 100, colour = NULL, frame_source = c("XF",
        "computed")) {
    frame_source <- match.arg(frame_source)
    ann <- .prepare_annotation_for_psite_bam(annotation)
    sample_map <- .make_sample_map(data, sample, multisamples)
    actual_samples <- unique(sample_map$tmp_samp)
    plot_samples <- unique(sample_map$sample)

    if (!plot_style %in% c("mirror", "dodge", "split", "facet"))
        plot_style <- "split"
    if (length(plot_samples) == 1L)
        plot_style <- "split"
    if (identical(plot_style, "mirror") && length(plot_samples) != 2L)
        plot_style <- "split"
    if (!region %in% c("all", "cds", "5utr", "3utr"))
        region <- "all"

    if (length(colour) == 0L)
        colour <- .default_gg_col(length(plot_samples))
    if (length(colour) < length(plot_samples)) {
        colour <- rep(colour, length.out = length(plot_samples))
    } else {
        colour <- colour[seq_along(plot_samples)]
    }
    names(colour) <- plot_samples

    l_transcripts <- as.character(ann[l_utr5 > 0 & l_cds > 0 & l_cds%%3 == 0 & l_utr3 > 0,
        transcript])
    c_transcripts <- if (length(transcripts) == 0L)
        l_transcripts else intersect(l_transcripts, transcripts)
    if (length(c_transcripts) == 0L)
        stop("No transcripts passed the frame annotation filters.")

    length_range <- .compute_length_range(data, actual_samples, c_transcripts, length_range,
        cl)
    target_regions <- if (identical(region, "all"))
        c("5utr", "cds", "3utr") else region

    final_dt <- rbindlist(lapply(actual_samples, function(samp) {
        frame_dt <- data[[samp]][length %in% length_range & transcript %in% c_transcripts]
        if (!identical(region, "all"))
            frame_dt <- frame_dt[psite_region == region]

        if (identical(frame_source, "XF")) {
            if (!"frame_bam" %in% names(frame_dt))
                stop("frame_source='XF' but frame_bam column is absent.")
            frame_dt[, frame := as.integer(frame_bam)]
        } else {
            frame_dt[, frame := psite_from_start%%3L]
        }
        frame_dt <- frame_dt[psite_region %in% target_regions & frame %in% 0:2]

        grid <- CJ(psite_region = target_regions, frame = 0:2, unique = TRUE)
        counts <- frame_dt[, .(count = .N), by = .(psite_region, frame)]
        tab <- merge(grid, counts, by = c("psite_region", "frame"), all.x = TRUE, sort = FALSE)
        tab[is.na(count), count := 0L]
        tab[, scaled_count := if (sum(count) > 0)
            (count/sum(count)) * 100 else 0, by = psite_region]
        tab[, tmp_samp := samp]
        tab
    }), use.names = TRUE)

    final_dt[, psite_region := factor(psite_region, levels = c("5utr", "cds", "3utr"),
        labels = c("5' UTR", "CDS", "3' UTR"))]
    final_dt <- final_dt[order(match(tmp_samp, actual_samples), psite_region, frame)]
    setnames(final_dt, "psite_region", "region")

    output <- list()
    output[["count_dt"]] <- copy(final_dt[, .(sample = tmp_samp, region, frame, count, scaled_count)])

    plot_dt <- merge(final_dt, sample_map, by = "tmp_samp", sort = FALSE)
    plot_dt <- plot_dt[, .(mean_scaled_count = mean(scaled_count), se_scaled_count = if (.N >
        1L)
        sd(scaled_count)/sqrt(.N) else NA_real_), by = .(sample, region, frame)]

    output[["plot_dt"]] <- copy(plot_dt[, .(sample, region, x = frame, y = mean_scaled_count,
        y_se = se_scaled_count)])
    plot_dt[, sample := factor(sample, levels = plot_samples)]

    oldw <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = oldw), add = TRUE)

    if (identical(plot_style, "split")) {
        for (samp in plot_samples) {
            sub_plot_dt <- plot_dt[as.character(sample) == samp]
            sel_col <- colour[[samp]]
            p <- ggplot(sub_plot_dt, aes(x = as.factor(frame), y = mean_scaled_count)) +
                geom_bar(stat = "identity", fill = sel_col, width = 0.8) + geom_errorbar(aes(ymin = mean_scaled_count -
                se_scaled_count, ymax = mean_scaled_count + se_scaled_count), width = 0.35,
                linewidth = 1.1, na.rm = TRUE, color = sel_col, show.legend = FALSE) + labs(title = samp,
                x = "Frame", y = "% P-sites") + theme_bw(base_size = 27.5) + theme(panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5))
            if (identical(region, "all")) {
                p <- p + facet_grid(. ~ region) + theme(strip.background = element_blank())
            }
            output[[paste0("plot_", samp)]] <- p
        }
    } else {
        if (identical(plot_style, "mirror")) {
            mirror_samp <- plot_samples[[2]]
            plot_dt[as.character(sample) == mirror_samp, mean_scaled_count := -mean_scaled_count]
        }

        p <- ggplot(plot_dt, aes(x = as.factor(frame), y = mean_scaled_count, fill = sample))

        if (plot_style %in% c("facet", "mirror")) {
            p <- p + geom_bar(stat = "identity", width = 0.8) + geom_errorbar(aes(ymin = mean_scaled_count -
                se_scaled_count, ymax = mean_scaled_count + se_scaled_count, color = sample),
                width = 0.35, linewidth = 1.1, na.rm = TRUE)
            if (identical(plot_style, "mirror"))
                p <- p + geom_hline(yintercept = 0, linetype = 2, color = "gray20")
        } else {
            p <- p + geom_bar(stat = "identity", position = position_dodge(0.9)) + geom_errorbar(aes(ymin = mean_scaled_count -
                se_scaled_count, ymax = mean_scaled_count + se_scaled_count, color = sample),
                width = 0.35, linewidth = 1.1, na.rm = TRUE, position = position_dodge(0.9))
        }

        p <- p + labs(x = "Frame", y = "% P-sites") + theme_bw(base_size = 27.5) + theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) + scale_fill_manual(values = colour) +
            scale_color_manual(values = colour) + scale_y_continuous(labels = abs)

        if (length(unique(colour)) > 1L && !identical(plot_style, "facet")) {
            p <- p + theme(legend.position = "right", legend.title = element_blank(), legend.background = element_blank())
        } else {
            p <- p + theme(legend.position = "none")
        }

        if (identical(region, "all")) {
            if (identical(plot_style, "facet")) {
                p <- p + facet_grid(sample ~ region)
            } else {
                p <- p + facet_grid(. ~ region)
            }
            p <- p + theme(strip.background = element_blank())
        } else {
            if (identical(plot_style, "facet")) {
                p <- p + facet_wrap(sample ~ ., ncol = ceiling(sqrt(length(plot_samples))))
            }
            p <- p + theme(strip.background = element_blank())
        }

        output[["plot"]] <- p
    }

    output
}
