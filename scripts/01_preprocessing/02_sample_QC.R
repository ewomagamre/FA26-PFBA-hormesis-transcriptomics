#!/usr/bin/env Rscript

###############################################################
# FA26 RNA-seq Gene Atlas Pipeline
# Script: 02_sample_QC.R
#
# Purpose:
# Sample-level QC and global transcriptomic structure for the
# 24-sample primary FA26 RNA-seq dataset.
#
# Primary samples:
#   Control, 2 h
#   Control, 4 h
#   PFBA 0.01 ug/g diet, 2 h
#   PFBA 0.01 ug/g diet, 4 h
#   PFBA 1 ug/g diet, 2 h
#   PFBA 1 ug/g diet, 4 h
#
# Control 0 h samples are excluded from the primary analysis.
#
# PCA:
#   500 most variable VST-transformed genes
#
# Heatmaps:
#   Euclidean distance on full VST matrix
#   Pearson correlation on full VST matrix
#
# Figure output:
#   PDF
#   TIFF, 600 dpi, LZW compression
###############################################################


#==============================================================
# 1. Load packages
#==============================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(grid)
})


#==============================================================
# 2. Directories
#==============================================================

object_dir      <- "objects"
main_figure_dir <- "figures/main"
supp_figure_dir <- "figures/supplementary"
qc_result_dir   <- "results/qc"

dir.create(
    main_figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    supp_figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    qc_result_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


#==============================================================
# 3. Remove previous QC figure versions
#==============================================================

old_figure_files <- c(
    file.path(
        main_figure_dir,
        "Figure1A_PCA.pdf"
    ),
    file.path(
        main_figure_dir,
        "Figure1A_PCA.tiff"
    ),
    file.path(
        main_figure_dir,
        "Figure1B_SampleDistanceHeatmap.pdf"
    ),
    file.path(
        main_figure_dir,
        "Figure1B_SampleDistanceHeatmap.tiff"
    ),
    file.path(
        supp_figure_dir,
        "FigureS1_SampleCorrelationHeatmap.pdf"
    ),
    file.path(
        supp_figure_dir,
        "FigureS1_SampleCorrelationHeatmap.tiff"
    ),
    file.path(
        supp_figure_dir,
        "FigureS2_LibrarySizeDistribution.pdf"
    ),
    file.path(
        supp_figure_dir,
        "FigureS2_LibrarySizeDistribution.tiff"
    )
)

unlink(old_figure_files)


#==============================================================
# 4. Input files
#==============================================================

dds_file <- file.path(
    object_dir,
    "FA26_dds_primary24.rds"
)

vsd_file <- file.path(
    object_dir,
    "FA26_vsd_primary24.rds"
)

metadata_file <- file.path(
    object_dir,
    "FA26_metadata_primary24.rds"
)

required_files <- c(
    dds_file,
    vsd_file,
    metadata_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        paste0(
            "Missing required files:\n",
            paste(
                missing_files,
                collapse = "\n"
            )
        )
    )
}


#==============================================================
# 5. Load objects
#==============================================================

dds      <- readRDS(dds_file)
vsd      <- readRDS(vsd_file)
metadata <- readRDS(metadata_file)

metadata <- as.data.frame(metadata)


#==============================================================
# 6. Validate sample structure
#==============================================================

if (is.null(rownames(metadata))) {
    stop(
        "Metadata must contain sample IDs as row names."
    )
}

if (!identical(
    colnames(dds),
    colnames(vsd)
)) {
    stop(
        "DDS and VST sample orders do not match."
    )
}

if (!all(
    colnames(dds) %in% rownames(metadata)
)) {
    stop(
        "Some DDS samples are absent from metadata."
    )
}

metadata <- metadata[
    colnames(dds),
    ,
    drop = FALSE
]

if (!identical(
    colnames(dds),
    rownames(metadata)
)) {
    stop(
        "Metadata could not be aligned to the expression objects."
    )
}

if (ncol(dds) != 24) {

    warning(
        paste(
            "Expected 24 primary samples but found",
            ncol(dds)
        )
    )
}


#==============================================================
# 7. Prepare metadata
#==============================================================

required_metadata_columns <- c(
    "SampleID",
    "Treatment",
    "Replicate",
    "Time(hrs)",
    "Condition"
)

missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {

    stop(
        paste0(
            "Metadata missing columns: ",
            paste(
                missing_metadata_columns,
                collapse = ", "
            )
        )
    )
}


# Raw treatment factor used for analysis and ggplot figures

metadata$Treatment <- factor(
    metadata$Treatment,
    levels = c(
        "Control",
        "PFBA_0.01",
        "PFBA_1"
    )
)


# Time factor

metadata$Time <- factor(
    as.character(
        metadata[["Time(hrs)"]]
    ),
    levels = c(
        "2",
        "4"
    )
)


# Condition order

metadata$Condition <- factor(
    metadata$Condition,
    levels = c(
        "Control_2h",
        "Control_4h",
        "PFBA_0.01_2h",
        "PFBA_0.01_4h",
        "PFBA_1_2h",
        "PFBA_1_4h"
    )
)


#==============================================================
# 8. Compact sample names for heatmaps
#
# C = Control
# L = Low PFBA, 0.01 ug/g diet
# H = High PFBA, 1 ug/g diet
#
# Examples:
#   C2_R1
#   L4_R3
#   H2_R4
#==============================================================

metadata$DisplayName <- ifelse(
    metadata$Treatment == "Control",
    paste0(
        "C",
        metadata$Time,
        "_R",
        metadata$Replicate
    ),
    ifelse(
        metadata$Treatment == "PFBA_0.01",
        paste0(
            "L",
            metadata$Time,
            "_R",
            metadata$Replicate
        ),
        paste0(
            "H",
            metadata$Time,
            "_R",
            metadata$Replicate
        )
    )
)

names(
    metadata$DisplayName
) <- rownames(metadata)


#==============================================================
# 9. Publication colors
#==============================================================

treatment_colors <- c(
    "Control"   = "#C73E3A",
    "PFBA_0.01" = "#2E8B57",
    "PFBA_1"    = "#3478B8"
)

time_colors <- c(
    "2" = "#BDBDBD",
    "4" = "#252525"
)


#==============================================================
# 10. Treatment annotation labels for heatmaps
#==============================================================

metadata$TreatmentPlot <- factor(
    metadata$Treatment,
    levels = c(
        "Control",
        "PFBA_0.01",
        "PFBA_1"
    ),
    labels = c(
        "Control",
        "Low PFBA (0.01 ug/g diet)",
        "High PFBA (1 ug/g diet)"
    )
)

annotation_data <- data.frame(
    Time = metadata$Time,
    Treatment = metadata$TreatmentPlot
)

rownames(annotation_data) <- rownames(metadata)

annotation_colors <- list(

    Time = c(
        "2" = "#BDBDBD",
        "4" = "#252525"
    ),

    Treatment = c(
        "Control" = "#C73E3A",
        "Low PFBA (0.01 ug/g diet)" = "#2E8B57",
        "High PFBA (1 ug/g diet)" = "#3478B8"
    )
)

sample_labels <- metadata$DisplayName

names(
    sample_labels
) <- rownames(metadata)


#==============================================================
# 11. Figure-saving functions
#==============================================================

save_ggplot <- function(
    plot_object,
    pdf_file,
    tiff_file,
    width,
    height
) {

    ggsave(
        filename = pdf_file,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        filename = tiff_file,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = 600,
        compression = "lzw"
    )
}


save_pheatmap <- function(
    heatmap_object,
    pdf_file,
    tiff_file,
    width,
    height
) {

    pdf(
        file = pdf_file,
        width = width,
        height = height,
        useDingbats = FALSE
    )

    grid.newpage()
    grid.draw(
        heatmap_object$gtable
    )

    dev.off()


    tiff(
        filename = tiff_file,
        width = width,
        height = height,
        units = "in",
        res = 600,
        compression = "lzw"
    )

    grid.newpage()
    grid.draw(
        heatmap_object$gtable
    )

    dev.off()
}


#==============================================================
# 12. VST matrix
#==============================================================

vst_matrix <- assay(vsd)


#==============================================================
# 13. PCA
#
# Match DESeq2 plotPCA behavior:
# use the 500 genes with highest variance.
#==============================================================

gene_variance <- apply(
    vst_matrix,
    1,
    var
)

ntop <- min(
    500,
    length(gene_variance)
)

top_gene_index <- order(
    gene_variance,
    decreasing = TRUE
)[
    seq_len(ntop)
]

pca_matrix <- vst_matrix[
    top_gene_index,
    ,
    drop = FALSE
]

pca_result <- prcomp(
    t(pca_matrix),
    center = TRUE,
    scale. = FALSE
)

percent_var <- round(
    100 *
        (
            pca_result$sdev^2 /
            sum(
                pca_result$sdev^2
            )
        ),
    digits = 1
)


pca_data <- data.frame(
    SampleID = rownames(
        pca_result$x
    ),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    stringsAsFactors = FALSE
)

pca_data$Treatment <- metadata[
    pca_data$SampleID,
    "Treatment"
]

pca_data$Time <- metadata[
    pca_data$SampleID,
    "Time"
]

pca_data$Condition <- metadata[
    pca_data$SampleID,
    "Condition"
]

pca_data$Replicate <- metadata[
    pca_data$SampleID,
    "Replicate"
]

pca_data$DisplayName <- metadata[
    pca_data$SampleID,
    "DisplayName"
]


write.csv(
    pca_data,
    file.path(
        qc_result_dir,
        "FA26_PCA_coordinates.csv"
    ),
    row.names = FALSE
)


#==============================================================
# 14. Publication PCA
#==============================================================

pca_plot <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        color = Treatment,
        shape = Time
    )
) +

    geom_hline(
        yintercept = 0,
        color = "grey82",
        linewidth = 0.45
    ) +

    geom_vline(
        xintercept = 0,
        color = "grey82",
        linewidth = 0.45
    ) +

    geom_point(
        size = 4.8,
        alpha = 0.95
    ) +

    scale_color_manual(
        name = "Treatment",
        values = treatment_colors,
        labels = c(
            "Control",
            "PFBA 0.01 ug/g diet",
            "PFBA 1 ug/g diet"
        ),
        drop = FALSE
    ) +

    scale_shape_manual(
        name = "Sampling time",
        values = c(
            "2" = 17,
            "4" = 15
        ),
        labels = c(
            "2 h",
            "4 h"
        ),
        drop = FALSE
    ) +

    labs(
        title = "Principal Component Analysis",
        x = paste0(
            "Principal component 1 (",
            percent_var[1],
            "%)"
        ),
        y = paste0(
            "Principal component 2 (",
            percent_var[2],
            "%)"
        )
    ) +

    guides(

        color = guide_legend(
            order = 1,
            override.aes = list(
                shape = 16,
                size = 4.5,
                alpha = 1
            )
        ),

        shape = guide_legend(
            order = 2,
            override.aes = list(
                color = "black",
                size = 4.5,
                alpha = 1
            )
        )
    ) +

    theme_classic(
        base_size = 13
    ) +

    theme(

        plot.title = element_text(
            size = 19,
            face = "bold",
            hjust = 0.5,
            margin = margin(
                b = 16
            )
        ),

        axis.title.x = element_text(
            size = 15,
            face = "bold",
            margin = margin(
                t = 10
            )
        ),

        axis.title.y = element_text(
            size = 15,
            face = "bold",
            margin = margin(
                r = 10
            )
        ),

        axis.text = element_text(
            size = 12,
            color = "black"
        ),

        axis.line = element_line(
            color = "black",
            linewidth = 0.8
        ),

        axis.ticks = element_line(
            color = "black",
            linewidth = 0.6
        ),

        legend.position = "right",

        legend.title = element_text(
            size = 13,
            face = "bold"
        ),

        legend.text = element_text(
            size = 11.5
        ),

        legend.key.height = unit(
            0.60,
            "cm"
        ),

        legend.spacing.y = unit(
            0.25,
            "cm"
        ),

        plot.margin = margin(
            t = 14,
            r = 18,
            b = 14,
            l = 14
        )
    )


save_ggplot(
    plot_object = pca_plot,
    pdf_file = file.path(
        main_figure_dir,
        "Figure1A_PCA.pdf"
    ),
    tiff_file = file.path(
        main_figure_dir,
        "Figure1A_PCA.tiff"
    ),
    width = 8.5,
    height = 6.5
)


#==============================================================
# 15. Sample-to-sample Euclidean distance
#==============================================================

sample_distance <- dist(
    t(vst_matrix),
    method = "euclidean"
)

sample_distance_matrix <- as.matrix(
    sample_distance
)


write.csv(
    sample_distance_matrix,
    file.path(
        qc_result_dir,
        "FA26_sample_distance_matrix.csv"
    ),
    row.names = TRUE
)


#==============================================================
# 16. Distance heatmap color palette
#
# Designed to reproduce the warm visual appearance of the
# previous FA26 heatmap:
#
# white/cream -> pale peach -> tan -> orange-brown -> dark brown
#==============================================================

distance_palette <- colorRampPalette(
    c(
        "#FFFDF9",
        "#FBE9D6",
        "#F4C99E",
        "#E9A567",
        "#D7823F",
        "#B75A27",
        "#873516",
        "#57200D"
    )
)(120)


#==============================================================
# 17. Sample-distance heatmap
#==============================================================

distance_heatmap <- pheatmap(

    sample_distance_matrix,

    color = distance_palette,

    clustering_distance_rows = sample_distance,
    clustering_distance_cols = sample_distance,

    clustering_method = "complete",

    annotation_row = annotation_data,
    annotation_col = annotation_data,

    annotation_colors = annotation_colors,

    labels_row = sample_labels[
        rownames(
            sample_distance_matrix
        )
    ],

    labels_col = sample_labels[
        colnames(
            sample_distance_matrix
        )
    ],

    border_color = NA,

    fontsize = 13,

    fontsize_row = 10.5,

    fontsize_col = 10.5,

    angle_col = 90,

    main = "Sample-to-Sample Distance",

    treeheight_row = 60,

    treeheight_col = 60,

    annotation_names_row = TRUE,

    annotation_names_col = TRUE,

    silent = TRUE
)


save_pheatmap(
    heatmap_object = distance_heatmap,
    pdf_file = file.path(
        main_figure_dir,
        "Figure1B_SampleDistanceHeatmap.pdf"
    ),
    tiff_file = file.path(
        main_figure_dir,
        "Figure1B_SampleDistanceHeatmap.tiff"
    ),
    width = 11.5,
    height = 9.5
)


#==============================================================
# 18. Sample Pearson correlation
#==============================================================

sample_correlation_matrix <- cor(
    vst_matrix,
    method = "pearson"
)


write.csv(
    sample_correlation_matrix,
    file.path(
        qc_result_dir,
        "FA26_sample_correlation_matrix.csv"
    ),
    row.names = TRUE
)


#==============================================================
# 19. Correlation heatmap
#==============================================================

correlation_palette <- colorRampPalette(
    c(
        "#313695",
        "#74ADD1",
        "#FFFFBF",
        "#F46D43",
        "#A50026"
    )
)(100)

correlation_breaks <- seq(
    min(
        sample_correlation_matrix
    ),
    1,
    length.out = 101
)


correlation_heatmap <- pheatmap(

    sample_correlation_matrix,

    color = correlation_palette,

    breaks = correlation_breaks,

    clustering_distance_rows = "correlation",

    clustering_distance_cols = "correlation",

    clustering_method = "complete",

    annotation_row = annotation_data,

    annotation_col = annotation_data,

    annotation_colors = annotation_colors,

    labels_row = sample_labels[
        rownames(
            sample_correlation_matrix
        )
    ],

    labels_col = sample_labels[
        colnames(
            sample_correlation_matrix
        )
    ],

    border_color = NA,

    fontsize = 12,

    fontsize_row = 9.5,

    fontsize_col = 9.5,

    angle_col = 90,

    main = "Sample-to-Sample Pearson Correlation",

    silent = TRUE
)


save_pheatmap(
    heatmap_object = correlation_heatmap,
    pdf_file = file.path(
        supp_figure_dir,
        "FigureS1_SampleCorrelationHeatmap.pdf"
    ),
    tiff_file = file.path(
        supp_figure_dir,
        "FigureS1_SampleCorrelationHeatmap.tiff"
    ),
    width = 11.5,
    height = 9.5
)


#==============================================================
# 20. Library size table
#==============================================================

library_sizes <- colSums(
    counts(dds)
)


library_data <- data.frame(

    SampleID = names(
        library_sizes
    ),

    LibrarySize = as.numeric(
        library_sizes
    ),

    stringsAsFactors = FALSE
)


library_data$LibrarySizeMillions <- (
    library_data$LibrarySize /
    1000000
)


library_data$Treatment <- metadata[
    library_data$SampleID,
    "Treatment"
]

library_data$Time <- metadata[
    library_data$SampleID,
    "Time"
]

library_data$Condition <- metadata[
    library_data$SampleID,
    "Condition"
]

library_data$Replicate <- metadata[
    library_data$SampleID,
    "Replicate"
]

library_data$DisplayName <- metadata[
    library_data$SampleID,
    "DisplayName"
]


library_data$DisplayName <- factor(
    library_data$DisplayName,
    levels = library_data$DisplayName[
        order(
            library_data$LibrarySizeMillions
        )
    ]
)


write.csv(
    library_data,
    file.path(
        qc_result_dir,
        "FA26_library_sizes.csv"
    ),
    row.names = FALSE
)


#==============================================================
# 21. Library-size figure
#==============================================================

median_library_size <- median(
    library_data$LibrarySizeMillions
)


library_plot <- ggplot(
    library_data,
    aes(
        x = DisplayName,
        y = LibrarySizeMillions,
        fill = Treatment
    )
) +

    geom_col(
        width = 0.78,
        color = "black",
        linewidth = 0.25
    ) +

    geom_hline(
        yintercept = median_library_size,
        linetype = "dashed",
        linewidth = 0.65,
        color = "grey35"
    ) +

    scale_fill_manual(
        name = "Treatment",
        values = treatment_colors,
        labels = c(
            "Control",
            "PFBA 0.01 ug/g diet",
            "PFBA 1 ug/g diet"
        ),
        drop = FALSE
    ) +

    labs(
        title = "RNA-seq Library Size Distribution",
        x = "Sample",
        y = "Library size (millions of reads)"
    ) +

    theme_classic(
        base_size = 12
    ) +

    theme(

        plot.title = element_text(
            size = 17,
            face = "bold",
            hjust = 0.5,
            margin = margin(
                b = 14
            )
        ),

        axis.title = element_text(
            size = 13.5,
            face = "bold"
        ),

        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 9,
            color = "black"
        ),

        axis.text.y = element_text(
            size = 10.5,
            color = "black"
        ),

        axis.line = element_line(
            linewidth = 0.7,
            color = "black"
        ),

        legend.position = "right",

        legend.title = element_text(
            face = "bold"
        )
    )


save_ggplot(
    plot_object = library_plot,
    pdf_file = file.path(
        supp_figure_dir,
        "FigureS2_LibrarySizeDistribution.pdf"
    ),
    tiff_file = file.path(
        supp_figure_dir,
        "FigureS2_LibrarySizeDistribution.tiff"
    ),
    width = 10,
    height = 6.5
)


#==============================================================
# 22. Save session information
#==============================================================

session_file <- file.path(
    qc_result_dir,
    "FA26_QC_sessionInfo.txt"
)

sink(
    session_file
)

print(
    sessionInfo()
)

sink()


#==============================================================
# 23. Completion report
#==============================================================

cat("\n")
cat("=========================================\n")
cat("Pipeline Step 02 completed successfully\n")
cat("=========================================\n")

cat(
    "Samples analyzed:",
    ncol(dds),
    "\n"
)

cat(
    "Genes analyzed  :",
    nrow(dds),
    "\n"
)

cat(
    "Genes used for PCA:",
    ntop,
    "\n"
)

cat("\n")

cat("PCA variance explained:\n")

cat(
    "  PC1:",
    percent_var[1],
    "%\n"
)

cat(
    "  PC2:",
    percent_var[2],
    "%\n"
)

cat("\n")

cat("Main figures:\n")
cat("  - Figure1A_PCA.pdf\n")
cat("  - Figure1A_PCA.tiff\n")
cat("  - Figure1B_SampleDistanceHeatmap.pdf\n")
cat("  - Figure1B_SampleDistanceHeatmap.tiff\n")

cat("\n")

cat("Supplementary figures:\n")
cat("  - FigureS1_SampleCorrelationHeatmap.pdf\n")
cat("  - FigureS1_SampleCorrelationHeatmap.tiff\n")
cat("  - FigureS2_LibrarySizeDistribution.pdf\n")
cat("  - FigureS2_LibrarySizeDistribution.tiff\n")

cat("\n")
cat("QC tables saved in results/qc/\n")
cat("=========================================\n")
