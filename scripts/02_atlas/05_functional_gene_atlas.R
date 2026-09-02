# ============================================================
# FA26 24-sample primary analysis
# Step 05A — Functional gene atlas
#
# Uses the fixed 3,000-gene UMAP from Step 04.
# No UMAP recomputation is performed.
#
# Current annotation source:
#   existing FA24 eggNOG annotations merged onto FA26 atlas
#
# When the missing-164 eggNOG run finishes, this script can be
# rerun after refreshing the master annotation table.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(patchwork)
})

cat("============================================================\n")
cat("FA26 Step 05A: Functional gene atlas\n")
cat("============================================================\n\n")


# ============================================================
# 1. Inputs / outputs
# ============================================================

input_rds <- "objects/FA26_gene_atlas_master_annotated.rds"

result_dir <- "results/functional_atlas"
figure_main_dir <- file.path("figures", "main")
figure_supp_dir <- file.path("figures", "supplementary")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_supp_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Load atlas
# ============================================================

atlas <- readRDS(input_rds)

stopifnot(
    nrow(atlas) == 3000,
    all(c("gene_id","UMAP1","UMAP2") %in% colnames(atlas))
)

cat("Atlas genes:", nrow(atlas), "\n")


# ============================================================
# 3. Annotation presence flags
# ============================================================

valid_annotation <- function(x) {
    !is.na(x) & x != "" & x != "-"
}

atlas <- atlas %>%
    mutate(
        has_description = valid_annotation(eggNOG_description),
        has_name = valid_annotation(preferred_name),
        has_GO = valid_annotation(GOs),
        has_KEGG_KO = valid_annotation(KEGG_ko),
        has_KEGG_pathway = valid_annotation(KEGG_Pathway),
        has_COG = valid_annotation(COG_category)
    )

annotation_summary <- tibble(
    annotation = c(
        "eggNOG description",
        "Preferred gene name",
        "GO",
        "KEGG KO",
        "KEGG pathway",
        "COG"
    ),
    n_genes = c(
        sum(atlas$has_description),
        sum(atlas$has_name),
        sum(atlas$has_GO),
        sum(atlas$has_KEGG_KO),
        sum(atlas$has_KEGG_pathway),
        sum(atlas$has_COG)
    )
) %>%
    mutate(
        percent = 100 * n_genes / nrow(atlas)
    )

write_csv(
    annotation_summary,
    file.path(result_dir, "FA26_annotation_coverage_summary.csv")
)


# ============================================================
# 4. Gene identity table
# ============================================================

gene_identity <- atlas %>%
    transmute(
        gene_id,
        transcript_id,
        protein_id,
        NCBI_product,
        preferred_name,
        eggNOG_description,
        COG_category,
        GOs,
        KEGG_ko,
        KEGG_Pathway,
        KEGG_Module,
        PFAMs
    )

write_csv(
    gene_identity,
    file.path(result_dir, "FA26_gene_functional_annotations.csv")
)


# ============================================================
# 5. COG parsing
# ============================================================
#
# A gene may contain multiple COG category letters.
# We retain one row per gene-category membership.
# ============================================================

cog_long <- atlas %>%
    filter(has_COG) %>%
    select(gene_id, UMAP1, UMAP2, COG_category) %>%
    mutate(
        COG_category = gsub("[^A-Z]", "", COG_category)
    ) %>%
    separate_rows(COG_category, sep = "") %>%
    filter(COG_category != "")

cog_key <- tribble(
    ~COG_category, ~COG_description,
    "J", "Translation, ribosomal structure and biogenesis",
    "A", "RNA processing and modification",
    "K", "Transcription",
    "L", "Replication, recombination and repair",
    "B", "Chromatin structure and dynamics",
    "D", "Cell cycle control and chromosome partitioning",
    "Y", "Nuclear structure",
    "V", "Defense mechanisms",
    "T", "Signal transduction mechanisms",
    "M", "Cell wall/membrane/envelope biogenesis",
    "N", "Cell motility",
    "Z", "Cytoskeleton",
    "W", "Extracellular structures",
    "U", "Intracellular trafficking, secretion and vesicular transport",
    "O", "Posttranslational modification, protein turnover and chaperones",
    "X", "Mobilome",
    "C", "Energy production and conversion",
    "G", "Carbohydrate transport and metabolism",
    "E", "Amino acid transport and metabolism",
    "F", "Nucleotide transport and metabolism",
    "H", "Coenzyme transport and metabolism",
    "I", "Lipid transport and metabolism",
    "P", "Inorganic ion transport and metabolism",
    "Q", "Secondary metabolites biosynthesis, transport and catabolism",
    "R", "General function prediction only",
    "S", "Function unknown"
)

cog_long <- cog_long %>%
    left_join(cog_key, by = "COG_category")

cog_summary <- cog_long %>%
    count(COG_category, COG_description, sort = TRUE, name = "n_genes")

write_csv(
    cog_summary,
    file.path(result_dir, "FA26_COG_summary.csv")
)

write_csv(
    cog_long,
    file.path(result_dir, "FA26_COG_gene_membership.csv")
)


# ============================================================
# 6. GO membership table
# ============================================================

go_long <- atlas %>%
    filter(has_GO) %>%
    select(gene_id, UMAP1, UMAP2, GOs) %>%
    separate_rows(GOs, sep = ",") %>%
    mutate(GO = trimws(GOs)) %>%
    select(-GOs) %>%
    filter(GO != "")

write_csv(
    go_long,
    file.path(result_dir, "FA26_GO_membership.csv")
)


# ============================================================
# 7. KEGG KO membership
# ============================================================

kegg_ko_long <- atlas %>%
    filter(has_KEGG_KO) %>%
    select(gene_id, UMAP1, UMAP2, KEGG_ko) %>%
    separate_rows(KEGG_ko, sep = ",") %>%
    mutate(
        KEGG_ko = trimws(KEGG_ko),
        KEGG_ko = sub("^ko:", "", KEGG_ko)
    ) %>%
    filter(KEGG_ko != "")

write_csv(
    kegg_ko_long,
    file.path(result_dir, "FA26_KEGG_KO_membership.csv")
)


# ============================================================
# 8. KEGG pathway membership
# ============================================================

kegg_path_long <- atlas %>%
    filter(has_KEGG_pathway) %>%
    select(gene_id, UMAP1, UMAP2, KEGG_Pathway) %>%
    separate_rows(KEGG_Pathway, sep = ",") %>%
    mutate(
        KEGG_Pathway = trimws(KEGG_Pathway)
    ) %>%
    filter(KEGG_Pathway != "")

kegg_path_summary <- kegg_path_long %>%
    count(KEGG_Pathway, sort = TRUE, name = "n_genes")

write_csv(
    kegg_path_long,
    file.path(result_dir, "FA26_KEGG_pathway_membership.csv")
)

write_csv(
    kegg_path_summary,
    file.path(result_dir, "FA26_KEGG_pathway_summary.csv")
)


# ============================================================
# 9. Nature-style plotting theme
# ============================================================

nature_theme <- theme_bw(base_size = 16) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 19,
            hjust = 0.5,
            margin = margin(b = 7)
        ),
        plot.subtitle = element_text(
            size = 13,
            hjust = 0.5,
            margin = margin(b = 8)
        ),
        axis.title = element_text(
            face = "bold",
            size = 17
        ),
        axis.text = element_text(
            size = 12.5,
            color = "black"
        ),
        axis.ticks = element_line(
            linewidth = 0.75,
            color = "black"
        ),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 1.05
        ),
        panel.grid.major = element_line(
            color = "#ECECEC",
            linewidth = 0.35
        ),
        panel.grid.minor = element_blank(),
        legend.title = element_text(
            face = "bold",
            size = 13.5
        ),
        legend.text = element_text(
            size = 12
        ),
        plot.margin = margin(6, 8, 6, 8)
    )


# ============================================================
# 10. Supplementary annotation coverage figure
# ============================================================

p_coverage <- ggplot(
    annotation_summary,
    aes(
        x = reorder(annotation, percent),
        y = percent
    )
) +
    geom_col(width = 0.72) +
    coord_flip() +
    geom_text(
        aes(
            label = paste0(
                n_genes,
                " (",
                sprintf("%.1f", percent),
                "%)"
            )
        ),
        hjust = -0.08,
        size = 4.6
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
        title = "Functional annotation coverage of the FA26 gene atlas",
        x = NULL,
        y = "Atlas genes annotated (%)"
    ) +
    nature_theme

ggsave(
    file.path(
        figure_supp_dir,
        "FigureS5_FunctionalAnnotationCoverage.pdf"
    ),
    p_coverage,
    width = 9,
    height = 5.8
)

ggsave(
    file.path(
        figure_supp_dir,
        "FigureS5_FunctionalAnnotationCoverage.tiff"
    ),
    p_coverage,
    width = 9,
    height = 5.8,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 11. COG-category atlas
# ============================================================
#
# For clarity, only genes with exactly ONE COG category are
# used in this categorical UMAP. Multi-category genes remain
# fully retained in the membership tables.
# ============================================================

cog_single <- atlas %>%
    filter(has_COG) %>%
    mutate(
        COG_clean = gsub("[^A-Z]", "", COG_category),
        n_cog_letters = nchar(COG_clean)
    ) %>%
    filter(n_cog_letters == 1) %>%
    transmute(
        gene_id,
        UMAP1,
        UMAP2,
        COG_category = COG_clean
    ) %>%
    left_join(cog_key, by = "COG_category")

top_cog_levels <- cog_single %>%
    count(COG_category, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(COG_category)

cog_plot_df <- cog_single %>%
    mutate(
        COG_plot = ifelse(
            COG_category %in% top_cog_levels,
            COG_category,
            "Other"
        )
    )

p_cog <- ggplot() +
    geom_point(
        data = atlas,
        aes(x = UMAP1, y = UMAP2),
        color = "#D9D9D9",
        size = 0.65,
        alpha = 0.65
    ) +
    geom_point(
        data = cog_plot_df,
        aes(
            x = UMAP1,
            y = UMAP2,
            color = COG_plot
        ),
        size = 0.95,
        alpha = 0.92
    ) +
    labs(
        title = "Functional organization of the FA26 transcriptomic atlas",
        subtitle = "Major single-category COG assignments mapped onto the fixed 3,000-gene UMAP",
        x = "UMAP1",
        y = "UMAP2",
        color = "COG category"
    ) +
    nature_theme

ggsave(
    file.path(
        figure_main_dir,
        "Figure3A_FunctionalAtlas_COG.pdf"
    ),
    p_cog,
    width = 9.5,
    height = 7.5
)

ggsave(
    file.path(
        figure_main_dir,
        "Figure3A_FunctionalAtlas_COG.tiff"
    ),
    p_cog,
    width = 9.5,
    height = 7.5,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 12. Top KEGG pathway atlas panels
# ============================================================
#
# Display top pathways by number of atlas genes.
# These are descriptive membership maps, not enrichment yet.
# ============================================================

top_pathways <- kegg_path_summary %>%
    slice_head(n = 6) %>%
    pull(KEGG_Pathway)

pathway_plot_df <- kegg_path_long %>%
    filter(KEGG_Pathway %in% top_pathways)

pathway_panels <- lapply(
    top_pathways,
    function(pathway_id) {

        members <- pathway_plot_df %>%
            filter(KEGG_Pathway == pathway_id)

        ggplot() +
            geom_point(
                data = atlas,
                aes(x = UMAP1, y = UMAP2),
                color = "#D9D9D9",
                size = 0.55,
                alpha = 0.55
            ) +
            geom_point(
                data = members,
                aes(x = UMAP1, y = UMAP2),
                size = 1.05,
                alpha = 0.95
            ) +
            labs(
                title = pathway_id,
                x = "UMAP1",
                y = "UMAP2"
            ) +
            nature_theme +
            theme(
                plot.title = element_text(
                    face = "bold",
                    size = 15,
                    hjust = 0.5
                )
            )
    }
)

pathway_figure <- wrap_plots(
    pathway_panels,
    ncol = 3
) +
    plot_annotation(
        title = "Spatial distribution of major KEGG pathway memberships across the FA26 gene atlas",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 19,
                hjust = 0.5
            )
        )
    )

ggsave(
    file.path(
        figure_main_dir,
        "Figure3B_KEGG_PathwayAtlas.pdf"
    ),
    pathway_figure,
    width = 15,
    height = 9.5
)

ggsave(
    file.path(
        figure_main_dir,
        "Figure3B_KEGG_PathwayAtlas.tiff"
    ),
    pathway_figure,
    width = 15,
    height = 9.5,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 13. Step 05A checkpoint
# ============================================================

step05 <- list(
    step = "FA26 Step 05A - Functional gene atlas",
    atlas = atlas,
    annotation_summary = annotation_summary,
    COG_membership = cog_long,
    GO_membership = go_long,
    KEGG_KO_membership = kegg_ko_long,
    KEGG_pathway_membership = kegg_path_long,
    KEGG_pathway_summary = kegg_path_summary,
    created = Sys.time()
)

saveRDS(
    step05,
    "objects/FA26_step05_functional_atlas.rds"
)

capture.output(
    sessionInfo(),
    file = file.path(
        result_dir,
        "FA26_step05_sessionInfo.txt"
    )
)


# ============================================================
# 14. Final report
# ============================================================

cat("\n============================================================\n")
cat("FA26 Step 05A complete\n")
cat("============================================================\n\n")

print(annotation_summary)

cat("\nTop 10 COG categories:\n")
print(head(cog_summary, 10))

cat("\nTop 10 KEGG pathway IDs:\n")
print(head(kegg_path_summary, 10))

cat("\nMain figures:\n")
cat("  figures/main/Figure3A_FunctionalAtlas_COG.pdf\n")
cat("  figures/main/Figure3A_FunctionalAtlas_COG.tiff\n")
cat("  figures/main/Figure3B_KEGG_PathwayAtlas.pdf\n")
cat("  figures/main/Figure3B_KEGG_PathwayAtlas.tiff\n")

cat("\nSupplementary:\n")
cat("  figures/supplementary/FigureS5_FunctionalAnnotationCoverage.pdf\n")
cat("  figures/supplementary/FigureS5_FunctionalAnnotationCoverage.tiff\n")

cat("\nNext: define biologically defensible modules and perform GO/KEGG enrichment.\n")
