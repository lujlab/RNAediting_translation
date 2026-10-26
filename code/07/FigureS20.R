library(data.table)
library(ggplot2)
library(ggpubr)
library(psych)
#-------------m6A
m6A_glori <- fread("/path/to/m6A/2023NB_GLORI_HEK293_m6A.txt")
if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    BiocManager::install("ChIPseeker")
}
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # hg38
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Convert sites to GRanges
m6A_gr <- makeGRangesFromDataFrame(
    m6A_glori,
    keep.extra.columns = TRUE,
    seqnames.field = "Chr",
    start.field = "Sites",
    end.field = "Sites",
    strand.field = "Strand"
)


print(m6A_gr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


anno <- annotatePeak(
    m6A_gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db",
    addFlankGeneInfo = TRUE,
    flankDistance = 5000,
    verbose = FALSE
)

print(anno)


anno_df <- as.data.frame(anno)
print(table(anno_df$annotation))

classify_region <- function(annotation) {
    if (grepl("Promoter", annotation, ignore.case = TRUE)) {
        return("Promoter")
    } else if (grepl("5' UTR|five prime UTR|UTR5", annotation, ignore.case = TRUE)) {
        return("UTR5")
    } else if (grepl("3' UTR|three prime UTR|UTR3", annotation, ignore.case = TRUE)) {
        return("UTR3")
    } else if (grepl("Exon|CDS|Coding", annotation, ignore.case = TRUE)) {
        return("CDS")
    } else if (grepl("Intron", annotation, ignore.case = TRUE)) {
        return("Intron")
    } else if (grepl("Intergenic", annotation, ignore.case = TRUE)) {
        return("Intergenic")
    } else {
        return("Other")
    }
}

anno_df <- anno_df %>%
    mutate(
        region_simple = sapply(annotation, classify_region)
    )

print(table(anno_df$region_simple))

m6A_glori$row_id <- 1:nrow(m6A_glori)
anno_df$row_id <- 1:nrow(anno_df)


combined_data <- merge(
    m6A_glori,
    anno_df[, c(
        "row_id", "annotation", "region_simple", "SYMBOL", "geneId",
        "distanceToTSS", "transcriptId"
    )],
    by = "row_id",
    all.x = TRUE
)

combined_data$row_id <- NULL

combined_data <- combined_data %>%
    mutate(
        gene_for_count = ifelse(is.na(SYMBOL) | SYMBOL == "", Gene, SYMBOL)
    )


gene_region_counts <- combined_data %>%
    group_by(gene_for_count, region_simple) %>%
    summarise(
        count = n(),
        .groups = "drop"
    ) %>%
    pivot_wider(
        id_cols = gene_for_count,
        names_from = region_simple,
        values_from = count,
        values_fill = 0
    )


all_regions <- c("CDS", "UTR3", "UTR5", "Promoter", "Intron", "Intergenic", "Other")
for (reg in all_regions) {
    if (!reg %in% names(gene_region_counts)) {
        gene_region_counts[[reg]] <- 0
    }
}

names(gene_region_counts)[names(gene_region_counts) == "NA"] <- "Unknown"
colnames(gene_region_counts) <- c("gene", "CDS_m6A", "UTR3_m6A", "UTR5_m6A", "Promoter_m6A", "Intron_m6A", "Intergenic_m6A", "Other_m6A", "Unknown_m6A")
print(head(gene_region_counts, 10))


# =======


combined_data <- fread("/path/to/RNAi/2026RNAediting/data/editing_halflife_combined", quote = FALSE)
HEK_hls <- combined_data[, c(2, 8)]

m6A_number_glori_data <- merge(gene_region_counts[, 1:4], HEK_hls, by = "gene")
head(m6A_number_glori_data)


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


m6A_long <- m6A_number_glori_data %>%
    pivot_longer(
        cols = c(CDS_m6A, UTR3_m6A, UTR5_m6A),
        names_to = "location",
        values_to = "m6A_count"
    ) %>%
    mutate(
        location = gsub("_m6A$", "", location),
        location = factor(location, levels = c("CDS", "UTR5", "UTR3"))
    )


m6A_long <- m6A_long %>%
    mutate(
        bin_m6A = case_when(
            m6A_count == 0 ~ "0",
            m6A_count >= 1 & m6A_count <= 3 ~ "1-3",
            m6A_count >= 4 & m6A_count <= 6 ~ "4-6",
            m6A_count >= 7 & m6A_count <= 9 ~ "7-9",
            m6A_count >= 10 & m6A_count <= 12 ~ "10-12",
            m6A_count > 12 ~ ">12"
        ),
        bin_m6A = factor(bin_m6A,
            levels = c("0", "1-3", "4-6", "7-9", "10-12", ">12"),
            ordered = TRUE
        )
    )


create_simple_stats_violin_m6A <- function(data) {
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(RColorBrewer)

    long_data <- data %>%
        pivot_longer(
            cols = c(CDS_m6A, UTR5_m6A, UTR3_m6A),
            names_to = "location",
            values_to = "m6A_count"
        ) %>%
        mutate(
            location = gsub("_m6A$", "", location),
            location = factor(location, levels = c("CDS", "UTR5", "UTR3")),
            bin_m6A = case_when(
                m6A_count == 0 ~ "0",
                m6A_count <= 3 ~ "1-3",
                m6A_count <= 6 ~ "4-6",
                m6A_count <= 9 ~ "7-9",
                m6A_count <= 12 ~ "10-12",
                TRUE ~ ">12"
            ),
            bin_m6A = factor(bin_m6A,
                levels = c("0", "1-3", "4-6", "7-9", "10-12", ">12")
            ),
            color_group = paste(location, bin_m6A, sep = "_")
        )

    median_hl <- median(data$half_life, na.rm = TRUE)

    get_gradient <- function(palette, n = 6) {
        colors <- RColorBrewer::brewer.pal(9, palette)
        return(colors[4:9])
    }

    cds_colors <- get_gradient("Oranges")
    utr5_colors <- get_gradient("Blues")
    utr3_colors <- get_gradient("Greens")

    all_colors <- c(
        setNames(cds_colors, paste("CDS", levels(long_data$bin_m6A), sep = "_")),
        setNames(utr5_colors, paste("UTR5", levels(long_data$bin_m6A), sep = "_")),
        setNames(utr3_colors, paste("UTR3", levels(long_data$bin_m6A), sep = "_"))
    )

    n_values <- long_data %>%
        group_by(location, bin_m6A) %>%
        summarise(n = n(), .groups = "drop")

    cor_results <- long_data %>%
        group_by(location) %>%
        summarise(
            spearman = cor.test(m6A_count, half_life, method = "spearman", exact = FALSE)$estimate,
            p_value = cor.test(m6A_count, half_life, method = "spearman", exact = FALSE)$p.value,
            .groups = "drop"
        ) %>%
        mutate(
            p_text = sprintf("p = %.2e", p_value),
            label = sprintf("<U+03C1> = %.3f\n%s", spearman, p_text)
        )

    p <- ggplot(long_data, aes(x = bin_m6A, y = half_life, fill = color_group)) +
        geom_violin(
            alpha = 0.8,
            position = position_dodge(0.9),
            width = 0.8
        ) +
        geom_boxplot(
            width = 0.15,
            alpha = 0.9,
            outlier.shape = NA,
            position = position_dodge(0.9)
        ) +
        facet_wrap(~location, nrow = 1) +
        coord_cartesian(ylim = c(0, 8)) +
        labs(
            x = "Number of m6A sites",
            y = "mRNA half-life (hours)",
            title = "mRNA Half-life by m6A Sites in Genomic Regions",
            subtitle = "Color gradients: CDS (orange), UTR5 (blue), UTR3 (green)"
        ) +
        geom_hline(
            yintercept = median_hl,
            linetype = "dashed",
            color = "red",
            size = 0.8
        ) +
        scale_fill_manual(values = all_colors) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "gray90", color = NA),
            legend.position = "none",
            panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1.2, "lines")
        )

    p <- p +
        geom_text(
            data = n_values,
            aes(x = bin_m6A, y = 7.8, label = paste0("n=", n)),
            size = 3.2,
            color = "gray30",
            position = position_dodge(width = 0.9),
            inherit.aes = FALSE
        )

    p <- p +
        geom_text(
            data = cor_results,
            aes(x = 3.5, y = 7.0, label = label),
            size = 3.5,
            fontface = "bold",
            hjust = 0.5,
            vjust = 0,
            inherit.aes = FALSE
        )
    p <- p +
        ggprism::theme_prism(
            base_fontface = "plain",
            base_size = 18,
            base_line_size = 0.8,
            axis_text_angle = 0
        )
    return(p)
}

simple_stats_plot_fixed <- create_simple_stats_violin_m6A(m6A_number_glori_data)
print(simple_stats_plot_fixed)


#--------------------------------------A-TO-I-------------------------------------
combined_data <- fread("/path/to/RNAi/2026RNAediting/data/editing_halflife_combined", quote = FALSE)

AI_number_data <- combined_data[, c(2, 4, 5, 6, 8)]
colnames(AI_number_data) <- c("gene", "CDS_edits", "UTR3_edits", "UTR5_edits", "half_life")
head(AI_number_data)

create_simple_stats_violin_editing <- function(data) {
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(RColorBrewer)

    long_data <- data %>%
        pivot_longer(
            cols = c(CDS_edits, UTR5_edits, UTR3_edits),
            names_to = "location",
            values_to = "edits_count"
        ) %>%
        mutate(
            location = gsub("_edits$", "", location),
            location = factor(location, levels = c("CDS", "UTR5", "UTR3")),
            bin_edits = case_when(
                edits_count == 0 ~ "0",
                edits_count <= 3 ~ "1-3",
                edits_count <= 6 ~ "4-6",
                edits_count <= 9 ~ "7-9",
                edits_count <= 12 ~ "10-12",
                TRUE ~ ">12"
            ),
            bin_edits = factor(bin_edits,
                levels = c("0", "1-3", "4-6", "7-9", "10-12", ">12")
            ),
            color_group = paste(location, bin_edits, sep = "_")
        )

    median_hl <- median(data$half_life, na.rm = TRUE)

    get_gradient <- function(palette, n = 6) {
        colors <- RColorBrewer::brewer.pal(9, palette)
        return(colors[4:9])
    }

    cds_colors <- get_gradient("Oranges")
    utr5_colors <- get_gradient("Blues")
    utr3_colors <- get_gradient("Greens")

    all_colors <- c(
        setNames(cds_colors, paste("CDS", levels(long_data$bin_edits), sep = "_")),
        setNames(utr5_colors, paste("UTR5", levels(long_data$bin_edits), sep = "_")),
        setNames(utr3_colors, paste("UTR3", levels(long_data$bin_edits), sep = "_"))
    )

    n_values <- long_data %>%
        group_by(location, bin_edits) %>%
        summarise(n = n(), .groups = "drop")

    cor_results <- long_data %>%
        group_by(location) %>%
        summarise(
            spearman = cor.test(edits_count, half_life, method = "spearman", exact = FALSE)$estimate,
            p_value = cor.test(edits_count, half_life, method = "spearman", exact = FALSE)$p.value,
            .groups = "drop"
        ) %>%
        mutate(
            p_text = sprintf("p = %.2e", p_value),
            label = sprintf("<U+03C1> = %.3f\n%s", spearman, p_text)
        )

    p <- ggplot(long_data, aes(x = bin_edits, y = half_life, fill = color_group)) +
        geom_violin(
            alpha = 0.8,
            position = position_dodge(0.9),
            width = 0.8
        ) +
        geom_boxplot(
            width = 0.15,
            alpha = 0.9,
            outlier.shape = NA,
            position = position_dodge(0.9)
        ) +
        facet_wrap(~location, nrow = 1) +
        coord_cartesian(ylim = c(0, 8)) +
        labs(
            x = "Number of editing sites",
            y = "mRNA half-life (hours)",
            title = "mRNA Half-life by editing Sites in Genomic Regions",
            subtitle = "Color gradients: CDS (orange), UTR5 (blue), UTR3 (green)"
        ) +
        geom_hline(
            yintercept = median_hl,
            linetype = "dashed",
            color = "red",
            size = 0.8
        ) +
        scale_fill_manual(values = all_colors) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "gray90", color = NA),
            legend.position = "none",
            panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1.2, "lines")
        )

    p <- p +
        geom_text(
            data = n_values,
            aes(x = bin_edits, y = 7.8, label = paste0("n=", n)),
            size = 3.2,
            color = "gray30",
            position = position_dodge(width = 0.9),
            inherit.aes = FALSE
        )

    p <- p +
        geom_text(
            data = cor_results,
            aes(x = 3.5, y = 7.0, label = label),
            size = 3.5,
            fontface = "bold",
            hjust = 0.5,
            vjust = 0,
            inherit.aes = FALSE
        )
    p <- p +
        ggprism::theme_prism(
            base_fontface = "plain",
            base_size = 18,
            base_line_size = 0.8,
            axis_text_angle = 0
        )
    return(p)
}

simple_stats_plot_editing <- create_simple_stats_violin_editing(AI_number_data)
print(simple_stats_plot_editing)


library(patchwork)
combined <- simple_stats_plot_fixed + simple_stats_plot_editing +
    plot_layout(
        guides = "collect",
        nrow = 2, # Two rows
        ncol = 1,
    )

combined
# 12 11
