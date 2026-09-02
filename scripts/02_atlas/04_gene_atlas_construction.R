# ============================================================
# FA26 24-sample primary analysis
# Step 04 — Gene-level transcriptomic atlas
#
# Definitive feature set:
#   3,000 genes selected in Step 03
#
# PURPOSE:
#   Build the unsupervised gene-expression atlas and descriptive
#   highest/lowest mean-expression overlays.
#
# IMPORTANT:
#   No k-means clustering is performed in this step.
#   No differential-expression inference is performed.
# ============================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(uwot)
    library(patchwork)
})

cat("============================================================\n")
cat("FA26 Step 04: Gene-level transcriptomic atlas\n")
cat("============================================================\n\n")


# ============================================================
# 1. Parameters
# ============================================================

checkpoint_file <- "objects/FA26_step03_atlas_optimization.rds"

umap_neighbors <- 30
umap_min_dist <- 0.30
umap_metric <- "euclidean"
umap_seed <- 260725

result_dir <- "results/gene_atlas"
figure_main_dir <- file.path("figures", "main")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_main_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Load Step 03 checkpoint
# ============================================================

step03 <- readRDS(checkpoint_file)

vst_mat <- step03$selected_vst_matrix
metadata <- step03$metadata
selected_ids <- step03$selected_gene_ids

stopifnot(
    is.matrix(vst_mat),
    nrow(vst_mat) == 3000,
    ncol(vst_mat) == 24,
    length(selected_ids) == 3000,
    identical(rownames(vst_mat), selected_ids),
    identical(colnames(vst_mat), rownames(metadata))
)

cat("Genes:", nrow(vst_mat), "\n")
cat("Samples:", ncol(vst_mat), "\n\n")


# ============================================================
# 3. Experimental groups
# ============================================================

condition_levels <- c(
    "Control_2h",
    "Control_4h",
    "PFBA_0.01_2h",
    "PFBA_0.01_4h",
    "PFBA_1_2h",
    "PFBA_1_4h"
)

stopifnot(
    setequal(
        unique(as.character(metadata$Condition)),
        condition_levels
    )
)

metadata$Condition <- factor(
    metadata$Condition,
    levels = condition_levels
)


# ============================================================
# 4. Sharp manuscript palette
# ============================================================
#
# Retains visual identity of the earlier FA atlas.
# Dark green has been removed.
# ============================================================

condition_colors <- c(
    "Control_2h"   = "#E6A000",
    "Control_4h"   = "#D73027",
    "PFBA_0.01_2h" = "#008C95",
    "PFBA_0.01_4h" = "#0066CC",
    "PFBA_1_2h"    = "#7A00CC",
    "PFBA_1_4h"    = "#E6007E"
)

condition_labels <- c(
    "Control_2h"   = "Control, 2 h",
    "Control_4h"   = "Control, 4 h",
    "PFBA_0.01_2h" = "PFBA 0.01, 2 h",
    "PFBA_0.01_4h" = "PFBA 0.01, 4 h",
    "PFBA_1_2h"    = "PFBA 1, 2 h",
    "PFBA_1_4h"    = "PFBA 1, 4 h"
)


# ============================================================
# 5. Row-standardize gene expression
# ============================================================

row_zscore <- function(mat) {

    mu <- rowMeans(mat)
    s <- apply(mat, 1, stats::sd)

    s[!is.finite(s) | s == 0] <- 1

    out <- sweep(mat, 1, mu, "-")
    out <- sweep(out, 1, s, "/")

    out
}

z_mat <- row_zscore(vst_mat)

stopifnot(all(is.finite(z_mat)))

saveRDS(
    z_mat,
    "objects/FA26_atlas_3000_row_zscore_matrix.rds"
)


# ============================================================
# 6. Definitive UMAP
# ============================================================

cat("Constructing definitive 3,000-gene UMAP...\n")

set.seed(umap_seed)

embedding <- uwot::umap(
    z_mat,
    n_neighbors = umap_neighbors,
    min_dist = umap_min_dist,
    metric = umap_metric,
    n_components = 2,
    init = "spectral",
    n_threads = 1,
    verbose = FALSE
)

rownames(embedding) <- rownames(z_mat)
colnames(embedding) <- c("UMAP1", "UMAP2")

atlas_df <- data.frame(
    gene_id = rownames(embedding),
    UMAP1 = embedding[, 1],
    UMAP2 = embedding[, 2],
    stringsAsFactors = FALSE
)


# ============================================================
# 7. Group mean expression
# ============================================================

condition_means <- sapply(
    condition_levels,
    function(cond) {

        ids <- rownames(metadata)[
            metadata$Condition == cond
        ]

        rowMeans(
            vst_mat[, ids, drop = FALSE]
        )
    }
)

rownames(condition_means) <- rownames(vst_mat)


# ============================================================
# 8. Highest / lowest mean group
# ============================================================

high_idx <- max.col(
    condition_means,
    ties.method = "first"
)

low_idx <- max.col(
    -condition_means,
    ties.method = "first"
)

atlas_df$highest_mean_group <- factor(
    colnames(condition_means)[high_idx],
    levels = condition_levels
)

atlas_df$lowest_mean_group <- factor(
    colnames(condition_means)[low_idx],
    levels = condition_levels
)


# ============================================================
# 9. Add descriptive margins
# ============================================================

sorted_high <- t(
    apply(condition_means, 1, sort, decreasing = TRUE)
)

sorted_low <- t(
    apply(condition_means, 1, sort, decreasing = FALSE)
)

atlas_df$highest_mean_margin <-
    sorted_high[, 1] - sorted_high[, 2]

atlas_df$lowest_mean_margin <-
    sorted_low[, 2] - sorted_low[, 1]

atlas_df$condition_mean_range <- apply(
    condition_means,
    1,
    function(x) max(x) - min(x)
)


# ============================================================
# 10. Add six condition means
# ============================================================

mean_df <- data.frame(
    gene_id = rownames(condition_means),
    condition_means,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

atlas_df <- atlas_df %>%
    left_join(mean_df, by = "gene_id")

write_csv(
    atlas_df,
    file.path(
        result_dir,
        "FA26_gene_atlas_master_unannotated.csv"
    )
)

write_csv(
    atlas_df %>% select(gene_id, UMAP1, UMAP2),
    file.path(
        result_dir,
        "FA26_gene_atlas_coordinates.csv"
    )
)


# ============================================================
# 11. Manuscript theme
# ============================================================

atlas_theme <- theme_classic(base_size = 14) +
    theme(
        axis.title = element_text(
            size = 14,
            face = "plain"
        ),
        axis.text = element_text(
            size = 10.5,
            color = "black"
        ),
        axis.line = element_line(
            linewidth = 0.7,
            color = "black"
        ),
        axis.ticks = element_line(
            linewidth = 0.55,
            color = "black"
        ),
        plot.title = element_text(
            face = "bold",
            size = 16,
            hjust = 0.5,
            margin = margin(b = 8)
        ),
        legend.title = element_text(
            face = "plain",
            size = 11.5
        ),
        legend.text = element_text(
            face = "plain",
            size = 10.5
        ),
        legend.key.height = unit(0.50, "cm"),
        legend.spacing.y = unit(0.08, "cm"),
        plot.margin = margin(10, 10, 10, 10)
    )


# ============================================================
# 12. Figure 2A — definitive atlas structure
# ============================================================
#
# No artificial cluster colors.
# Each point is one gene.
# ============================================================

p_atlas <- ggplot(
    atlas_df,
    aes(x = UMAP1, y = UMAP2)
) +
    geom_point(
        size = 0.78,
        alpha = 0.68,
        color = "#303030"
    ) +
    labs(
        title = "Gene-level transcriptomic atlas",
        x = "UMAP 1",
        y = "UMAP 2"
    ) +
    atlas_theme

ggsave(
    file.path(
        figure_main_dir,
        "Figure2A_GeneAtlas.pdf"
    ),
    p_atlas,
    width = 7.5,
    height = 6.5
)

ggsave(
    file.path(
        figure_main_dir,
        "Figure2A_GeneAtlas.tiff"
    ),
    p_atlas,
    width = 7.5,
    height = 6.5,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# ============================================================
# ============================================================
# ============================================================
# 13. Figure 2B — Nature-style high / low expression atlas
# ============================================================

# One definitive UMAP geometry is shown in both panels.
#
# Panel A:
#   gene colored by experimental group with highest mean
#   VST expression.
#
# Panel B:
#   gene colored by experimental group with lowest mean
#   VST expression.
#
# Color is descriptive and does NOT represent differential
# expression significance.
# ============================================================


# ------------------------------------------------------------
# Sharp six-condition palette
# ------------------------------------------------------------

condition_colors <- c(
    "Control_2h"   = "#E69F00",
    "Control_4h"   = "#D73027",
    "PFBA_0.01_2h" = "#008B8B",
    "PFBA_0.01_4h" = "#0066CC",
    "PFBA_1_2h"    = "#7400C8",
    "PFBA_1_4h"    = "#EC008C"
)

condition_labels <- c(
    "Control_2h"   = "Control, 2 h",
    "Control_4h"   = "Control, 4 h",
    "PFBA_0.01_2h" = "PFBA 0.01, 2 h",
    "PFBA_0.01_4h" = "PFBA 0.01, 4 h",
    "PFBA_1_2h"    = "PFBA 1, 2 h",
    "PFBA_1_4h"    = "PFBA 1, 4 h"
)


# ------------------------------------------------------------
# Use identical axis limits in both panels
# ------------------------------------------------------------
#
# This ensures the geometry is directly comparable between
# highest- and lowest-expression maps.
# ------------------------------------------------------------

x_range <- range(atlas_df$UMAP1)
y_range <- range(atlas_df$UMAP2)

x_pad <- diff(x_range) * 0.04
y_pad <- diff(y_range) * 0.04

x_limits <- c(
    x_range[1] - x_pad,
    x_range[2] + x_pad
)

y_limits <- c(
    y_range[1] - y_pad,
    y_range[2] + y_pad
)


# ------------------------------------------------------------
# Nature-style atlas theme
# ------------------------------------------------------------

nature_atlas_theme <- theme_bw(base_size = 16) +
    theme(

        # Panel title
        plot.title = element_text(
            face = "bold",
            size = 20,
            hjust = 0.5,
            margin = margin(b = 7)
        ),

        # Panel tag
        plot.tag = element_text(
            face = "bold",
            size = 18,
            hjust = 0
        ),

        plot.tag.position = c(
            0.015,
            0.985
        ),

        # Axis titles
        axis.title.x = element_text(
            face = "bold",
            size = 18,
            margin = margin(t = 8)
        ),

        axis.title.y = element_text(
            face = "bold",
            size = 18,
            margin = margin(r = 8)
        ),

        # Numeric axis labels
        axis.text.x = element_text(
            size = 13.5,
            color = "black"
        ),

        axis.text.y = element_text(
            size = 13.5,
            color = "black"
        ),

        # Tick marks
        axis.ticks = element_line(
            linewidth = 0.8,
            color = "black"
        ),

        axis.ticks.length = unit(
            0.17,
            "cm"
        ),

        # Full strong frame
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 1.05
        ),

        # Very light major guides only
        panel.grid.major = element_line(
            color = "#E9E9E9",
            linewidth = 0.38
        ),

        panel.grid.minor = element_blank(),

        # Clean background
        panel.background = element_rect(
            fill = "white",
            color = NA
        ),

        plot.background = element_rect(
            fill = "white",
            color = NA
        ),

        # Legend hierarchy
        legend.title = element_text(
            face = "bold",
            size = 14.5,
            hjust = 0,
            margin = margin(b = 6)
        ),

        legend.text = element_text(
            face = "plain",
            size = 13.2,
            color = "black"
        ),

        legend.key = element_blank(),

        legend.key.height = unit(
            0.55,
            "cm"
        ),

        legend.key.width = unit(
            0.55,
            "cm"
        ),

        legend.spacing.y = unit(
            0.02,
            "cm"
        ),

        legend.margin = margin(
            0, 0, 0, 0
        ),

        # Keep plot area dominant
        plot.margin = margin(
            t = 5,
            r = 6,
            b = 5,
            l = 6
        )
    )


# ------------------------------------------------------------
# Panel A — Highest mean expression
# ------------------------------------------------------------

p_high <- ggplot(
    atlas_df,
    aes(
        x = UMAP1,
        y = UMAP2,
        color = highest_mean_group
    )
) +

    geom_point(
        size = 1.05,
        alpha = 0.98
    ) +

    scale_color_manual(
        values = condition_colors,
        labels = condition_labels,
        breaks = condition_levels,
        drop = FALSE,
        name = "Experimental group"
    ) +

    scale_x_continuous(
        limits = x_limits,
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = 0)
    ) +

    scale_y_continuous(
        limits = y_limits,
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = 0)
    ) +

    guides(
        color = guide_legend(
            override.aes = list(
                size = 4.6,
                alpha = 1
            ),
            title.position = "top",
            title.hjust = 0
        )
    ) +

    labs(
        tag = "A",
        title = "Highest mean expression",
        x = "UMAP1",
        y = "UMAP2"
    ) +

    nature_atlas_theme


# ------------------------------------------------------------
# Panel B — Lowest mean expression
# ------------------------------------------------------------

p_low <- ggplot(
    atlas_df,
    aes(
        x = UMAP1,
        y = UMAP2,
        color = lowest_mean_group
    )
) +

    geom_point(
        size = 1.05,
        alpha = 0.98
    ) +

    scale_color_manual(
        values = condition_colors,
        labels = condition_labels,
        breaks = condition_levels,
        drop = FALSE,
        name = "Experimental group"
    ) +

    scale_x_continuous(
        limits = x_limits,
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = 0)
    ) +

    scale_y_continuous(
        limits = y_limits,
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = 0)
    ) +

    guides(
        color = guide_legend(
            override.aes = list(
                size = 4.6,
                alpha = 1
            ),
            title.position = "top",
            title.hjust = 0
        )
    ) +

    labs(
        tag = "B",
        title = "Lowest mean expression",
        x = "UMAP1",
        y = "UMAP2"
    ) +

    nature_atlas_theme


# ------------------------------------------------------------
# Extract one shared legend
# ------------------------------------------------------------

legend_source <- p_high +
    theme(
        legend.position = "right"
    )

shared_legend <- cowplot::get_legend(
    legend_source
)


# Data panels do not draw legends

p_high_clean <- p_high +
    theme(
        legend.position = "none"
    )

p_low_clean <- p_low +
    theme(
        legend.position = "none"
    )


# ------------------------------------------------------------
# Final composition
# ------------------------------------------------------------
#
# The entire figure is intentionally wide.
# Both UMAPs receive large equal plotting regions.
# The center legend is compact and does not become a panel.
# ------------------------------------------------------------

legend_panel <- cowplot::ggdraw() +
    cowplot::draw_grob(
        shared_legend,
        x = 0,
        y = 0.22,
        width = 1,
        height = 0.56
    )

high_low <- cowplot::plot_grid(
    p_high_clean,
    legend_panel,
    p_low_clean,
    nrow = 1,
    rel_widths = c(
        1,
        0.27,
        1
    ),
    align = "h",
    axis = "tb"
)


# ------------------------------------------------------------
# Replace Figure 2B outputs
# ------------------------------------------------------------

ggsave(
    file.path(
        figure_main_dir,
        "Figure2B_GeneAtlas_HighLowExpression.pdf"
    ),
    high_low,
    width = 17.5,
    height = 7.6,
    units = "in",
    bg = "white"
)

ggsave(
    file.path(
        figure_main_dir,
        "Figure2B_GeneAtlas_HighLowExpression.tiff"
    ),
    high_low,
    width = 17.5,
    height = 7.6,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

cat("Final Nature-style Figure 2B generated.\n")

# 16. Step 04 checkpoint
# ============================================================

step04 <- list(
    step = "FA26 Step 04 - Gene-level transcriptomic atlas",
    selected_n = 3000L,
    gene_ids = selected_ids,
    vst_matrix = vst_mat,
    row_zscore_matrix = z_mat,
    metadata = metadata,
    umap_coordinates = atlas_df[, c("gene_id", "UMAP1", "UMAP2")],
    atlas_master = atlas_df,
    condition_mean_expression = condition_means,
    condition_colors = condition_colors,
    condition_labels = condition_labels,
    created = Sys.time()
)

saveRDS(
    step04,
    "objects/FA26_step04_gene_atlas.rds"
)


# ============================================================
# 17. Session info
# ============================================================

capture.output(
    sessionInfo(),
    file = file.path(
        result_dir,
        "FA26_step04_sessionInfo.txt"
    )
)


# ============================================================
# 18. Remove obsolete k-means outputs
# ============================================================

obsolete <- c(
    "figures/main/Figure2A_GeneAtlas_Clusters.pdf",
    "figures/main/Figure2A_GeneAtlas_Clusters.tiff",
    "figures/main/Figure2C_ClusterExpressionHeatmap.pdf",
    "figures/main/Figure2C_ClusterExpressionHeatmap.tiff",
    "figures/main/Figure2D_ClusterTemporalProfiles.pdf",
    "figures/main/Figure2D_ClusterTemporalProfiles.tiff",
    "figures/supplementary/FigureS5_GeneClusterOptimization.pdf",
    "figures/supplementary/FigureS5_GeneClusterOptimization.tiff"
)

for (f in obsolete) {
    if (file.exists(f)) file.remove(f)
}


# ============================================================
# 19. Final summary
# ============================================================

cat("\n============================================================\n")
cat("FA26 Step 04 complete\n")
cat("============================================================\n\n")

cat("Genes in atlas:", nrow(atlas_df), "\n")
cat("No k-means clustering performed.\n\n")

cat("Highest mean-expression counts:\n")
print(table(atlas_df$highest_mean_group))

cat("\nLowest mean-expression counts:\n")
print(table(atlas_df$lowest_mean_group))

cat("\nMain figures replaced:\n")
cat("  figures/main/Figure2A_GeneAtlas.pdf\n")
cat("  figures/main/Figure2A_GeneAtlas.tiff\n")
cat("  figures/main/Figure2B_GeneAtlas_HighLowExpression.pdf\n")
cat("  figures/main/Figure2B_GeneAtlas_HighLowExpression.tiff\n")

cat("\nCheckpoint:\n")
cat("  objects/FA26_step04_gene_atlas.rds\n")

cat("\nReady for eggNOG gene annotation.\n")
