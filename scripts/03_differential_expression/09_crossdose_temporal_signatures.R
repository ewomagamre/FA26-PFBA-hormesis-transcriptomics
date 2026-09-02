
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
})

cat("\n============================================================\n")
cat("FA26 STEP 09 — CROSS-DOSE TEMPORAL SIGNATURES\n")
cat("============================================================\n\n")

high_file <-
    "results/deseq2/signatures/FA26_PFBA1_2h_4h_all_genes_signature_master.csv"

low_file <-
    "results/deseq2/signatures_lowdose/FA26_PFBA001_all_genes_signature_master.csv"

outdir <-
    "results/deseq2/crossdose"

figdir <-
    "figures/main"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figdir,
    recursive = TRUE,
    showWarnings = FALSE
)

high <- read_csv(
    high_file,
    show_col_types = FALSE
)

low <- read_csv(
    low_file,
    show_col_types = FALSE
)

# ============================================================
# 1. Build common gene-level comparison table
# ============================================================

cross <- high %>%
    dplyr::select(
        gene_id,

        high_sig_2h = sig_2h,
        high_sig_4h = sig_4h,

        high_direction_2h = direction_2h,
        high_direction_4h = direction_4h,

        high_LFC_2h = shrunk_LFC_2h,
        high_LFC_4h = shrunk_LFC_4h,

        high_padj_2h = padj_2h,
        high_padj_4h = padj_4h,

        high_signature = signature_class,

        preferred_name,
        NCBI_product,
        eggNOG_description,
        zone
    ) %>%

    full_join(
        low %>%
            dplyr::select(
                gene_id,

                low_sig_2h = sig_2h,
                low_sig_4h = sig_4h,

                low_direction_2h =
                    direction_2h,

                low_direction_4h =
                    direction_4h,

                low_LFC_2h =
                    shrunk_LFC_2h,

                low_LFC_4h =
                    shrunk_LFC_4h,

                low_padj_2h =
                    padj_2h,

                low_padj_4h =
                    padj_4h,

                low_signature =
                    signature_class
            ),
        by = "gene_id"
    )

# ============================================================
# 2. Cross-dose categories
# ============================================================

cross <- cross %>%
    mutate(

        shared_2h =
            high_sig_2h &
            low_sig_2h,

        shared_4h =
            high_sig_4h &
            low_sig_4h,

        high2_low4 =
            high_sig_2h &
            low_sig_4h,

        high2_low4_same_direction =
            high2_low4 &
            high_direction_2h ==
                low_direction_4h,

        high2_low4_opposite_direction =
            high2_low4 &
            high_direction_2h !=
                low_direction_4h,

        common_persistent =
            high_sig_2h &
            high_sig_4h &
            low_sig_2h &
            low_sig_4h
    )

write_csv(
    cross,
    file.path(
        outdir,
        "FA26_crossdose_gene_master.csv"
    )
)

# ============================================================
# 3. Key overlap summary
# ============================================================

overlap_summary <- tibble(

    comparison = c(
        "Both doses significant at 2 h",
        "Both doses significant at 4 h",
        "High dose 2 h AND low dose 4 h",
        "High dose 2 h AND low dose 4 h, same direction",
        "High dose 2 h AND low dose 4 h, opposite direction",
        "Significant at both times in both doses"
    ),

    n_genes = c(

        sum(
            cross$shared_2h,
            na.rm = TRUE
        ),

        sum(
            cross$shared_4h,
            na.rm = TRUE
        ),

        sum(
            cross$high2_low4,
            na.rm = TRUE
        ),

        sum(
            cross$high2_low4_same_direction,
            na.rm = TRUE
        ),

        sum(
            cross$high2_low4_opposite_direction,
            na.rm = TRUE
        ),

        sum(
            cross$common_persistent,
            na.rm = TRUE
        )
    )
)

write_csv(
    overlap_summary,
    file.path(
        outdir,
        "FA26_crossdose_overlap_summary.csv"
    )
)

# ============================================================
# 4. DEG architecture table for Figure 5A
#
# Use clean condition-model results already reflected in
# Step 07/08 master tables.
# ============================================================

deg_counts <- bind_rows(

    tibble(
        dose = "PFBA 0.01 µg/g",
        time = "2 h",
        Up =
            sum(
                low$sig_2h &
                low$direction_2h == "Up",
                na.rm = TRUE
            ),
        Down =
            sum(
                low$sig_2h &
                low$direction_2h == "Down",
                na.rm = TRUE
            )
    ),

    tibble(
        dose = "PFBA 0.01 µg/g",
        time = "4 h",
        Up =
            sum(
                low$sig_4h &
                low$direction_4h == "Up",
                na.rm = TRUE
            ),
        Down =
            sum(
                low$sig_4h &
                low$direction_4h == "Down",
                na.rm = TRUE
            )
    ),

    tibble(
        dose = "PFBA 1 µg/g",
        time = "2 h",
        Up =
            sum(
                high$sig_2h &
                high$direction_2h == "Up",
                na.rm = TRUE
            ),
        Down =
            sum(
                high$sig_2h &
                high$direction_2h == "Down",
                na.rm = TRUE
            )
    ),

    tibble(
        dose = "PFBA 1 µg/g",
        time = "4 h",
        Up =
            sum(
                high$sig_4h &
                high$direction_4h == "Up",
                na.rm = TRUE
            ),
        Down =
            sum(
                high$sig_4h &
                high$direction_4h == "Down",
                na.rm = TRUE
            )
    )
)

deg_long <- deg_counts %>%

    pivot_longer(
        cols = c(
            Up,
            Down
        ),
        names_to = "direction",
        values_to = "n_genes"
    ) %>%

    mutate(
        signed_count =
            ifelse(
                direction == "Down",
                -n_genes,
                n_genes
            ),

        dose =
            factor(
                dose,
                levels = c(
                    "PFBA 0.01 µg/g",
                    "PFBA 1 µg/g"
                )
            ),

        time =
            factor(
                time,
                levels = c(
                    "2 h",
                    "4 h"
                )
            )
    )

# ============================================================
# 5. Figure 5A — DEG architecture
# ============================================================

pA <- ggplot(
    deg_long,
    aes(
        x = time,
        y = signed_count,
        fill = direction
    )
) +

    geom_hline(
        yintercept = 0,
        linewidth = 0.65
    ) +

    geom_col(
        width = 0.68
    ) +

    facet_wrap(
        ~dose,
        nrow = 1
    ) +

    scale_fill_manual(
        values = c(
            "Up" = "#B2182B",
            "Down" = "#2166AC"
        ),
        name = NULL
    ) +

    scale_y_continuous(
        labels = abs
    ) +

    labs(
        tag = "A",
        x = "Time after exposure",
        y = "Number of DEGs"
    ) +

    theme_classic(
        base_size = 16
    ) +

    theme(
        strip.background =
            element_blank(),

        strip.text =
            element_text(
                face = "bold",
                size = 17
            ),

        axis.title =
            element_text(
                face = "bold",
                size = 17
            ),

        axis.text =
            element_text(
                size = 13,
                color = "black"
            ),

        legend.position =
            "bottom",

        legend.text =
            element_text(
                size = 13
            ),

        plot.tag =
            element_text(
                face = "bold",
                size = 20
            )
    )

# ============================================================
# 6. Figure 5B — signature architecture by dose
# ============================================================

sig_counts <- bind_rows(

    high %>%
        filter(
            signature_class !=
                "Not significant"
        ) %>%
        count(
            signature_class,
            name = "n_genes"
        ) %>%
        mutate(
            dose =
                "PFBA 1 µg/g"
        ),

    low %>%
        filter(
            signature_class !=
                "Not significant"
        ) %>%
        count(
            signature_class,
            name = "n_genes"
        ) %>%
        mutate(
            dose =
                "PFBA 0.01 µg/g"
        )
)

sig_counts$signature_class <- factor(
    sig_counts$signature_class,
    levels = c(
        "Early-only",
        "Persistent",
        "Later-emerging"
    )
)

pB <- ggplot(
    sig_counts,
    aes(
        x = signature_class,
        y = n_genes,
        group = dose,
        color = dose
    )
) +

    geom_line(
        linewidth = 1.5
    ) +

    geom_point(
        size = 4.5
    ) +

    scale_color_manual(
        values = c(
            "PFBA 0.01 µg/g" = "#0072B2",
            "PFBA 1 µg/g" = "#7A00CC"
        ),
        name = NULL
    ) +

    labs(
        tag = "B",
        x = NULL,
        y = "Number of genes"
    ) +

    theme_classic(
        base_size = 16
    ) +

    theme(
        axis.title.y =
            element_text(
                face = "bold",
                size = 17
            ),

        axis.text =
            element_text(
                size = 13,
                color = "black"
            ),

        axis.text.x =
            element_text(
                angle = 20,
                hjust = 1
            ),

        legend.position =
            "bottom",

        legend.text =
            element_text(
                size = 13
            ),

        plot.tag =
            element_text(
                face = "bold",
                size = 20
            )
    )

# ============================================================
# 7. Assemble Figure 5
# ============================================================

fig5 <- pA / pB +

    plot_layout(
        heights = c(
            1,
            0.9
        )
    )

ggsave(
    file.path(
        figdir,
        "Figure5_DESeq2_TemporalSignatureArchitecture.pdf"
    ),
    fig5,
    width = 13.5,
    height = 9.5,
    units = "in",
    bg = "white"
)

ggsave(
    file.path(
        figdir,
        "Figure5_DESeq2_TemporalSignatureArchitecture.tiff"
    ),
    fig5,
    width = 13.5,
    height = 9.5,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

cat("\n============================================================\n")
cat("STEP 09 COMPLETE\n")
cat("============================================================\n\n")

cat("CROSS-DOSE OVERLAPS\n")
print(overlap_summary)

cat("\nDEG COUNTS\n")
print(deg_counts)

cat("\nSIGNATURE COUNTS\n")
print(sig_counts)

cat("\nFigure:\n")
cat("  figures/main/Figure5_DESeq2_TemporalSignatureArchitecture.pdf\n")
cat("  figures/main/Figure5_DESeq2_TemporalSignatureArchitecture.tiff\n")

