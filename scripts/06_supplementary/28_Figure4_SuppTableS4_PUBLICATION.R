#!/usr/bin/env Rscript

# ============================================================
# FA26 — PUBLICATION SUPPLEMENTARY TABLE S4
#
# Figure 4 source data
#
# S4A — GO-BP enrichment
# S4B — Persistent high-dose genes
# S4C — Concordant shared 4-h genes
# S4D — Gene-level functional modules underlying Figure 4D
#
# IMPORTANT:
# No DESeq2 model, GO enrichment, or biological analysis is
# rerun here.
#
# Existing locked Figure 4 source tables are reformatted and
# Figure 4D gene membership is reconstructed using the exact
# module definitions used in the final Figure 4 analysis.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(openxlsx)
})

options(stringsAsFactors = FALSE)

read_csv <- function(file, ...) {
    read.csv(
        file,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        ...
    )
}

cat("\n============================================================\n")
cat("FIGURE 4 — PUBLICATION SUPPLEMENTARY TABLE S4\n")
cat("============================================================\n\n")

# ============================================================
# 1. INPUTS
# ============================================================

fig4dir <- "results/figure4"
finaldir <- file.path(fig4dir, "final")

fileA <- file.path(
    finaldir,
    "FA26_Figure4_GO_BP_temporal_ORA_significant_v2.csv"
)

fileB <- file.path(
    finaldir,
    "FA26_Figure4B_genes_v2.csv"
)

fileC <- file.path(
    finaldir,
    "FA26_Figure4C_concordant_shared_genes.csv"
)

fileD_summary <- file.path(
    finaldir,
    "FA26_Figure4D_module_counts_v2.csv"
)

stress_file <- file.path(
    fig4dir,
    "stress_2h",
    "FA26_stress_fullcontrast_high_vs_low_2h.csv"
)

adaptive_file <- file.path(
    fig4dir,
    "FA26_adaptive_metabolic_candidate_genes.csv"
)

annotation_file <- file.path(
    fig4dir,
    "FA26_Figure4_standardized_annotation.csv"
)

high2_file <- "results/deseq2/PFBA_1_vs_Control_2h.csv"
low2_file  <- "results/deseq2/PFBA_0.01_vs_Control_2h.csv"

required <- c(
    fileA,
    fileB,
    fileC,
    fileD_summary,
    stress_file,
    adaptive_file,
    annotation_file,
    high2_file,
    low2_file
)

missing <- required[!file.exists(required)]

if (length(missing) > 0) {
    cat("Missing source files:\n")
    print(missing)
    stop("Required Figure 4 source file(s) missing.")
}

A <- read_csv(fileA)
B <- read_csv(fileB)
C <- read_csv(fileC)

D_locked <- read_csv(
    fileD_summary
)

stress <- read_csv(
    stress_file
)

adaptive <- read_csv(
    adaptive_file
)

ann <- read_csv(
    annotation_file
)

high2 <- read_csv(
    high2_file
)

low2 <- read_csv(
    low2_file
)

cat("Source files loaded successfully.\n")

# ============================================================
# 2. S4A — GO-BP ENRICHMENT
#
# Keep all 270 significant enrichment results.
#
# Add exact Figure 4A display status:
# score = -log10(padj) * gene_ratio
# top 5 per displayed temporal program after removing
# duplicate descriptions.
# ============================================================

stopifnot(
    nrow(A) == 270
)

A$Figure4A_score <-
    -log10(
        pmax(
            A$padj,
            1e-300
        )
    ) *
    A$gene_ratio

A$Displayed_in_Figure4A <- FALSE

set_order <- c(
    "1 ug/g Early",
    "1 ug/g Persistent",
    "1 ug/g Emerging",
    "0.01 ug/g Emerging"
)

for (s in set_order) {

    idx <- which(
        A$gene_set == s &
        !is.na(A$padj) &
        A$padj < 0.05 &
        !is.na(A$Description) &
        A$Description != ""
    )

    if (length(idx) == 0)
        next

    z <- A[idx, , drop = FALSE]

    ord <- order(
        -z$Figure4A_score
    )

    z <- z[
        ord,
        ,
        drop = FALSE
    ]

    keep_unique <-
        !duplicated(
            z$Description
        )

    z <- z[
        keep_unique,
        ,
        drop = FALSE
    ]

    z <- head(
        z,
        5
    )

    keys <- paste(
        z$gene_set,
        z$GO_ID,
        z$Description,
        sep = "|||"
    )

    all_keys <- paste(
        A$gene_set,
        A$GO_ID,
        A$Description,
        sep = "|||"
    )

    A$Displayed_in_Figure4A[
        all_keys %in% keys
    ] <- TRUE
}

S4A <- A %>%
    select(
        gene_set,
        GO_ID,
        Description,
        overlap,
        set_annotated,
        GO_background,
        gene_ratio,
        pvalue,
        padj,
        Figure4A_score,
        Displayed_in_Figure4A,
        genes
    ) %>%
    arrange(
        gene_set,
        padj,
        desc(gene_ratio)
    )

cat(
    "✓ S4A:",
    nrow(S4A),
    "significant GO-BP results\n"
)

cat(
    "  Terms displayed in Figure 4A:",
    sum(S4A$Displayed_in_Figure4A),
    "\n"
)

# ============================================================
# 3. S4B — COMPLETE PERSISTENT 1 ug/g PROGRAM
#
# Biological definition:
#   significant at 1 ug/g vs time-matched Control at BOTH
#   2 h and 4 h, using condition-model DESeq2 Wald padj < 0.05.
#
# Effect-size stabilization:
#   apeglm-shrunken log2FC.
#
# The main heatmap displays 18 selected Up -> Up genes.
# S4B retains ALL 138 persistent genes so reviewers can inspect
# the complete persistent response.
# ============================================================

high_master_file <- paste0(
    "results/deseq2/signatures/",
    "FA26_PFBA1_2h_4h_all_genes_signature_master.csv"
)

B_display_file <- paste0(
    "results/figure4/final/",
    "FA26_Figure4B_persistent_genes_displayed.csv"
)

stopifnot(
    file.exists(high_master_file),
    file.exists(B_display_file)
)

high_master <- read.csv(
    high_master_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

B_display <- read.csv(
    B_display_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

stopifnot(
    nrow(high_master) == 11661,
    length(unique(high_master$gene_id)) == 11661,
    length(unique(B_display$gene_id)) == 18
)

# ------------------------------------------------------------
# Complete significant 2 h ∩ 4 h population
# ------------------------------------------------------------

S4B <- high_master %>%

    filter(
        sig_2h %in% TRUE,
        sig_4h %in% TRUE
    ) %>%

    mutate(

        Persistent_direction =
            case_when(

                direction_2h == "Up" &
                direction_4h == "Up" ~
                    "Up -> Up",

                direction_2h == "Down" &
                direction_4h == "Down" ~
                    "Down -> Down",

                direction_2h == "Up" &
                direction_4h == "Down" ~
                    "Up -> Down",

                direction_2h == "Down" &
                direction_4h == "Up" ~
                    "Down -> Up",

                TRUE ~
                    "Other"
            ),

        Concordant_persistent =
            direction_2h == direction_4h,

        Displayed_in_Figure4B =
            gene_id %in%
                B_display$gene_id
    )

# ------------------------------------------------------------
# Recover exact Figure 4B selection metrics
# driver_count and mean_effect
# ------------------------------------------------------------

B_metrics <- B_display %>%
    dplyr::select(
        gene_id,
        any_of(
            c(
                "driver_count",
                "mean_effect"
            )
        )
    ) %>%
    distinct(
        gene_id,
        .keep_all = TRUE
    )

S4B <- S4B %>%
    left_join(
        B_metrics,
        by = "gene_id"
    )

# For genes not selected for the figure, mean apeglm effect can
# still be calculated across the two time points.
S4B$Mean_apeglm_effect_2h_4h <-
    rowMeans(
        cbind(
            as.numeric(
                S4B$shrunk_LFC_2h
            ),
            as.numeric(
                S4B$shrunk_LFC_4h
            )
        ),
        na.rm = TRUE
    )

# ------------------------------------------------------------
# Reviewer-facing column order
# ------------------------------------------------------------

S4B <- S4B %>%
    transmute(

        gene_id,

        gene_name =
            preferred_name,

        NCBI_product,

        eggNOG_description,

        Persistent_direction,

        Concordant_persistent,

        Displayed_in_Figure4B,

        GO_driver_count =
            if (
                "driver_count" %in%
                    names(.)
            ) driver_count
            else NA_real_,

        Figure_selection_mean_effect =
            if (
                "mean_effect" %in%
                    names(.)
            ) mean_effect
            else NA_real_,

        Mean_apeglm_effect_2h_4h,

        Significant_2h =
            sig_2h,

        Direction_2h =
            direction_2h,

        DESeq2_log2FC_2h =
            LFC_2h,

        DESeq2_SE_2h =
            SE_2h,

        Wald_pvalue_2h =
            pvalue_2h,

        Wald_padj_2h =
            padj_2h,

        apeglm_log2FC_2h =
            shrunk_LFC_2h,

        apeglm_SE_2h =
            shrunk_SE_2h,

        Significant_4h =
            sig_4h,

        Direction_4h =
            direction_4h,

        DESeq2_log2FC_4h =
            LFC_4h,

        DESeq2_SE_4h =
            SE_4h,

        Wald_pvalue_4h =
            pvalue_4h,

        Wald_padj_4h =
            padj_4h,

        apeglm_log2FC_4h =
            shrunk_LFC_4h,

        apeglm_SE_4h =
            shrunk_SE_4h,

        Replicate_support_2h =
            support_2h,

        Replicate_support_4h =
            support_4h,

        Overall_support =
            overall_support,

        Atlas_zone =
            zone,

        UMAP1,
        UMAP2,

        COG_category,
        GO_terms = GOs,
        KEGG_KO = KEGG_ko,
        KEGG_pathway = KEGG_Pathway,
        PFAMs
    )

# ------------------------------------------------------------
# S4B validation
# ------------------------------------------------------------

stopifnot(
    nrow(S4B) == 138,

    sum(
        S4B$Persistent_direction ==
            "Up -> Up"
    ) == 123,

    sum(
        S4B$Persistent_direction ==
            "Down -> Down"
    ) == 12,

    sum(
        !S4B$Concordant_persistent
    ) == 3,

    sum(
        S4B$Displayed_in_Figure4B
    ) == 18,

    all(
        S4B$Persistent_direction[
            S4B$Displayed_in_Figure4B
        ] == "Up -> Up"
    )
)

cat(
    "✓ S4B:",
    nrow(S4B),
    "persistent 1 ug/g genes\n"
)

cat(
    "    Up -> Up:",
    sum(
        S4B$Persistent_direction ==
            "Up -> Up"
    ),
    "\n"
)

cat(
    "    Down -> Down:",
    sum(
        S4B$Persistent_direction ==
            "Down -> Down"
    ),
    "\n"
)

cat(
    "    Direction-switching:",
    sum(
        !S4B$Concordant_persistent
    ),
    "\n"
)

cat(
    "    Displayed in Figure 4B:",
    sum(
        S4B$Displayed_in_Figure4B
    ),
    "\n"
)

# ============================================================
# 4. S4C — COMPLETE CROSS-DOSE SHARED 4-h RESPONSE
#
# Biological definition:
#   significant relative to time-matched Control at 4 h for
#   BOTH 0.01 ug/g and 1 ug/g PFBA.
#
# The full intersection contains exactly 22 genes:
#   14 concordantly upregulated
#    8 concordantly downregulated
#    0 discordant
#
# All 22 are displayed in the heatmap.
# ============================================================

low_master_file <- paste0(
    "results/deseq2/signatures_lowdose/",
    "FA26_PFBA001_all_genes_signature_master.csv"
)

C_display_file <- paste0(
    "results/figure4/final/",
    "FA26_Figure4C_shared22_displayed.csv"
)

stopifnot(
    file.exists(low_master_file),
    file.exists(C_display_file)
)

low_master <- read.csv(
    low_master_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

C_display <- read.csv(
    C_display_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

stopifnot(
    nrow(low_master) == 11661,
    setequal(
        high_master$gene_id,
        low_master$gene_id
    ),
    length(unique(C_display$gene_id)) == 22
)

# High-dose 4-h significant genes
high4_shared <- high_master %>%
    filter(
        sig_4h %in% TRUE
    ) %>%
    transmute(

        gene_id,

        Preferred_name_high =
            preferred_name,

        NCBI_product_high =
            NCBI_product,

        eggNOG_description_high =
            eggNOG_description,

        High_direction_4h =
            direction_4h,

        High_DESeq2_log2FC_4h =
            LFC_4h,

        High_DESeq2_SE_4h =
            SE_4h,

        High_Wald_pvalue_4h =
            pvalue_4h,

        High_Wald_padj_4h =
            padj_4h,

        High_apeglm_log2FC_4h =
            shrunk_LFC_4h,

        High_apeglm_SE_4h =
            shrunk_SE_4h,

        Atlas_zone_high =
            zone
    )

# Low-dose 4-h significant genes
low4_shared <- low_master %>%
    filter(
        sig_4h %in% TRUE
    ) %>%
    transmute(

        gene_id,

        Preferred_name_low =
            preferred_name,

        NCBI_product_low =
            NCBI_product,

        eggNOG_description_low =
            eggNOG_description,

        Low_direction_4h =
            direction_4h,

        Low_DESeq2_log2FC_4h =
            LFC_4h,

        Low_DESeq2_SE_4h =
            SE_4h,

        Low_Wald_pvalue_4h =
            pvalue_4h,

        Low_Wald_padj_4h =
            padj_4h,

        Low_apeglm_log2FC_4h =
            shrunk_LFC_4h,

        Low_apeglm_SE_4h =
            shrunk_SE_4h,

        Atlas_zone_low =
            zone
    )

S4C <- inner_join(
    high4_shared,
    low4_shared,
    by = "gene_id"
) %>%

    mutate(

        Shared_direction =
            case_when(

                High_direction_4h == "Up" &
                Low_direction_4h == "Up" ~
                    "Shared Up",

                High_direction_4h == "Down" &
                Low_direction_4h == "Down" ~
                    "Shared Down",

                TRUE ~
                    "Discordant"
            ),

        Concordant_effect =
            High_direction_4h ==
            Low_direction_4h,

        Displayed_in_Figure4C =
            gene_id %in%
                C_display$gene_id,

        gene_name =
            ifelse(
                !is.na(
                    Preferred_name_high
                ) &
                Preferred_name_high != "" &
                Preferred_name_high != "-",
                Preferred_name_high,
                Preferred_name_low
            ),

        product =
            ifelse(
                !is.na(
                    NCBI_product_high
                ) &
                NCBI_product_high != "",
                NCBI_product_high,
                NCBI_product_low
            ),

        description =
            ifelse(
                !is.na(
                    eggNOG_description_high
                ) &
                eggNOG_description_high != "",
                eggNOG_description_high,
                eggNOG_description_low
            ),

        Atlas_zone =
            ifelse(
                !is.na(
                    Atlas_zone_high
                ) &
                Atlas_zone_high != "",
                Atlas_zone_high,
                Atlas_zone_low
            )
    ) %>%

    dplyr::select(
        gene_id,
        gene_name,
        product,
        description,
        Atlas_zone,
        Shared_direction,
        Concordant_effect,
        Displayed_in_Figure4C,
        High_direction_4h,
        High_DESeq2_log2FC_4h,
        High_DESeq2_SE_4h,
        High_Wald_pvalue_4h,
        High_Wald_padj_4h,
        High_apeglm_log2FC_4h,
        High_apeglm_SE_4h,
        Low_direction_4h,
        Low_DESeq2_log2FC_4h,
        Low_DESeq2_SE_4h,
        Low_Wald_pvalue_4h,
        Low_Wald_padj_4h,
        Low_apeglm_log2FC_4h,
        Low_apeglm_SE_4h
    )

# ------------------------------------------------------------
# S4C validation
# ------------------------------------------------------------

stopifnot(
    nrow(S4C) == 22,

    sum(
        S4C$Shared_direction ==
            "Shared Up"
    ) == 14,

    sum(
        S4C$Shared_direction ==
            "Shared Down"
    ) == 8,

    sum(
        !S4C$Concordant_effect
    ) == 0,

    sum(
        S4C$Displayed_in_Figure4C
    ) == 22,

    setequal(
        S4C$gene_id,
        C_display$gene_id
    )
)

cat(
    "✓ S4C:",
    nrow(S4C),
    "cross-dose shared 4-h genes\n"
)

cat(
    "    Shared Up:",
    sum(
        S4C$Shared_direction ==
            "Shared Up"
    ),
    "\n"
)

cat(
    "    Shared Down:",
    sum(
        S4C$Shared_direction ==
            "Shared Down"
    ),
    "\n"
)

cat(
    "    Discordant:",
    sum(
        !S4C$Concordant_effect
    ),
    "\n"
)

cat(
    "    Displayed in Figure 4C:",
    sum(
        S4C$Displayed_in_Figure4C
    ),
    "\n"
)

# ============================================================
# 5. S4D — EXACT FIGURE 4D MODULE DEFINITIONS
#
# From scripts/20_Figure4_refine_GO_labels_and_shared_genes.R:
#
# Canonical stress/redox:
#     all unique stress2 gene IDs
#
# Mitochondrial respiration:
#     adaptive screen_category == "Mitochondrial_Respiration"
#
# ATP/bioenergetics:
#     adaptive screen_category == "ATP_Bioenergetics"
#
# Figure 4D significance:
#     condition-model DESeq2 padj < 0.05 at 2 h
# ============================================================

stress_ids <- unique(
    stress$gene_id
)

mito_ids <- unique(
    adaptive$gene_id[
        adaptive$screen_category ==
            "Mitochondrial_Respiration"
    ]
)

atp_ids <- unique(
    adaptive$gene_id[
        adaptive$screen_category ==
            "ATP_Bioenergetics"
    ]
)

modules <- list(

    `Canonical stress/redox` =
        stress_ids,

    `Mitochondrial respiration` =
        mito_ids,

    `ATP/bioenergetics` =
        atp_ids
)

cat("\nFigure 4D module candidate sizes:\n")

for (m in names(modules)) {

    cat(
        "  ",
        m,
        ": ",
        length(modules[[m]]),
        " genes\n",
        sep = ""
    )
}

# ============================================================
# 6. GENE-LEVEL MODULE MEMBERSHIP
#
# A gene belonging to multiple biological modules legitimately
# appears once for each module because Figure 4D counts each
# module independently.
# ============================================================

module_membership <- do.call(
    rbind,
    lapply(
        names(modules),
        function(m) {

            ids <- modules[[m]]

            data.frame(
                Module = rep(
                    m,
                    length(ids)
                ),
                gene_id = ids,
                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(module_membership) <- NULL

# ============================================================
# 7. DESEQ2 2-h STATISTICS
#
# These are the exact DESeq2 tables whose padj values were used
# to determine Figure 4D significance.
# ============================================================

high_stats <- high2 %>%
    transmute(
        gene_id,

        High_baseMean_2h =
            baseMean,

        High_DESeq2_log2FC_2h =
            log2FoldChange,

        High_DESeq2_SE_2h =
            lfcSE,

        High_Wald_pvalue_2h =
            pvalue,

        High_Wald_padj_2h =
            padj,

        High_sig_2h =
            !is.na(padj) &
            padj < 0.05
    )

low_stats <- low2 %>%
    transmute(
        gene_id,

        Low_baseMean_2h =
            baseMean,

        Low_DESeq2_log2FC_2h =
            log2FoldChange,

        Low_DESeq2_SE_2h =
            lfcSE,

        Low_Wald_pvalue_2h =
            pvalue,

        Low_Wald_padj_2h =
            padj,

        Low_sig_2h =
            !is.na(padj) &
            padj < 0.05
    )

# ============================================================
# 8. STANDARDIZED GENE ANNOTATION
# ============================================================

ann <- ann[
    !duplicated(
        ann$gene_id
    ),
    ,
    drop = FALSE
]

# Helper: safely retrieve a possible annotation column
pick_column <- function(
    dat,
    candidates
) {

    hit <- candidates[
        candidates %in%
        names(dat)
    ]

    if (length(hit) == 0) {

        return(
            rep(
                NA_character_,
                nrow(dat)
            )
        )
    }

    as.character(
        dat[[hit[1]]]
    )
}

annotation <- data.frame(
    gene_id = ann$gene_id,

    gene_name =
        pick_column(
            ann,
            c(
                "gene_name",
                "preferred_name",
                "Preferred_name",
                "display_name"
            )
        ),

    product =
        pick_column(
            ann,
            c(
                "product",
                "NCBI_product",
                "Product"
            )
        ),

    description =
        pick_column(
            ann,
            c(
                "description",
                "eggNOG_description",
                "Description"
            )
        ),

    GO =
        pick_column(
            ann,
            c(
                "GO",
                "GOs",
                "GO_terms"
            )
        ),

    KEGG_KO =
        pick_column(
            ann,
            c(
                "KEGG_KO",
                "KEGG_ko"
            )
        ),

    KEGG_Pathway =
        pick_column(
            ann,
            c(
                "KEGG_Pathway",
                "KEGG_pathway"
            )
        ),

    PFAM =
        pick_column(
            ann,
            c(
                "PFAM",
                "PFAMs"
            )
        ),

    COG =
        pick_column(
            ann,
            c(
                "COG",
                "COG_category"
            )
        ),

    stringsAsFactors = FALSE
)

# ============================================================
# 9. MODULE-SPECIFIC INFORMATION
#
# Stress source has stress category + product.
# Adaptive source has adaptive screen category and shrunken
# effect sizes.
# ============================================================

stress_meta <- stress %>%
    transmute(
        gene_id,

        Stress_display_name =
            display_name,

        Stress_product =
            product,

        Stress_category =
            stress_category,

        Stress_Low_LFC_2h =
            Low_LFC_2h,

        Stress_Low_padj_2h =
            Low_padj_2h,

        Stress_High_LFC_2h =
            High_LFC_2h,

        Stress_High_padj_2h =
            High_padj_2h
    ) %>%
    distinct(
        gene_id,
        .keep_all = TRUE
    )

adaptive_meta <- adaptive %>%
    transmute(
        gene_id,

        Adaptive_gene_name =
            gene_name,

        Adaptive_display_name =
            display_name,

        Adaptive_product =
            product,

        Adaptive_description =
            description,

        Adaptive_screen =
            screen,

        Adaptive_screen_category =
            screen_category,

        High_apeglm_log2FC_2h =
            PFBA1_shrunk_LFC_2h,

        High_apeglm_padj_2h =
            PFBA1_padj_2h,

        Low_apeglm_log2FC_2h =
            PFBA001_shrunk_LFC_2h,

        Low_apeglm_padj_2h =
            PFBA001_padj_2h
    ) %>%
    distinct(
        gene_id,
        .keep_all = TRUE
    )

# ============================================================
# 10. BUILD REVIEWER-FACING S4D
# ============================================================

S4D <- module_membership %>%

    left_join(
        annotation,
        by = "gene_id"
    ) %>%

    left_join(
        stress_meta,
        by = "gene_id"
    ) %>%

    left_join(
        adaptive_meta,
        by = "gene_id"
    ) %>%

    left_join(
        low_stats,
        by = "gene_id"
    ) %>%

    left_join(
        high_stats,
        by = "gene_id"
    )

# ------------------------------------------------------------
# Improve human-readable gene naming.
#
# Prefer standardized annotation, then adaptive name/display,
# then stress display name, finally gene_id.
# ------------------------------------------------------------

usable_text <- function(x) {

    !is.na(x) &
    x != "" &
    x != "-"
}

S4D$Display_name <- S4D$gene_name

idx <- !usable_text(S4D$Display_name) &
    usable_text(S4D$Adaptive_gene_name)

S4D$Display_name[idx] <-
    S4D$Adaptive_gene_name[idx]

idx <- !usable_text(S4D$Display_name) &
    usable_text(S4D$Adaptive_display_name)

S4D$Display_name[idx] <-
    S4D$Adaptive_display_name[idx]

idx <- !usable_text(S4D$Display_name) &
    usable_text(S4D$Stress_display_name)

S4D$Display_name[idx] <-
    S4D$Stress_display_name[idx]

idx <- !usable_text(S4D$Display_name)

S4D$Display_name[idx] <-
    S4D$gene_id[idx]

# Prefer best available product
S4D$Final_product <- S4D$product

idx <- !usable_text(S4D$Final_product) &
    usable_text(S4D$Adaptive_product)

S4D$Final_product[idx] <-
    S4D$Adaptive_product[idx]

idx <- !usable_text(S4D$Final_product) &
    usable_text(S4D$Stress_product)

S4D$Final_product[idx] <-
    S4D$Stress_product[idx]

# Prefer best available description
S4D$Final_description <- S4D$description

idx <- !usable_text(S4D$Final_description) &
    usable_text(S4D$Adaptive_description)

S4D$Final_description[idx] <-
    S4D$Adaptive_description[idx]

# ------------------------------------------------------------
# Explicit Figure 4D response classification
# ------------------------------------------------------------

S4D$Figure4D_response <- case_when(

    S4D$Low_sig_2h &
    S4D$High_sig_2h ~
        "Both doses",

    S4D$Low_sig_2h &
    !S4D$High_sig_2h ~
        "0.01 µg/g only",

    !S4D$Low_sig_2h &
    S4D$High_sig_2h ~
        "1 µg/g only",

    TRUE ~
        "Neither significant"
)

# Reviewer-facing column order
S4D <- S4D %>%
    select(
        Module,

        gene_id,
        Display_name,
        Final_product,
        Final_description,

        Stress_category,
        Adaptive_screen_category,

        Low_sig_2h,
        Low_DESeq2_log2FC_2h,
        Low_DESeq2_SE_2h,
        Low_Wald_pvalue_2h,
        Low_Wald_padj_2h,

        High_sig_2h,
        High_DESeq2_log2FC_2h,
        High_DESeq2_SE_2h,
        High_Wald_pvalue_2h,
        High_Wald_padj_2h,

        Low_apeglm_log2FC_2h,
        High_apeglm_log2FC_2h,

        Figure4D_response,

        GO,
        KEGG_KO,
        KEGG_Pathway,
        PFAM,
        COG
    ) %>%

    arrange(
        factor(
            Module,
            levels = c(
                "Canonical stress/redox",
                "Mitochondrial respiration",
                "ATP/bioenergetics"
            )
        ),
        desc(High_sig_2h),
        High_Wald_padj_2h,
        gene_id
    )

cat(
    "✓ S4D gene-level module table:",
    nrow(S4D),
    "gene × module records\n"
)

# ============================================================
# 11. RECONSTRUCT FIGURE 4D FROM S4D
# ============================================================

D_check <- bind_rows(
    S4D %>%
        group_by(Module) %>%
        summarise(
            Dose = "0.01",
            Significant =
                sum(
                    Low_sig_2h,
                    na.rm = TRUE
                ),
            .groups = "drop"
        ),

    S4D %>%
        group_by(Module) %>%
        summarise(
            Dose = "1",
            Significant =
                sum(
                    High_sig_2h,
                    na.rm = TRUE
                ),
            .groups = "drop"
        )
) %>%

    arrange(
        Module,
        Dose
    )

D_locked2 <- D_locked %>%
    mutate(
        Module = as.character(Module),
        Dose = as.character(Dose),
        Significant = as.integer(Significant)
    ) %>%
    arrange(
        Module,
        Dose
    )

D_check <- D_check %>%
    mutate(
        Module = as.character(Module),
        Dose = as.character(Dose),
        Significant = as.integer(Significant)
    ) %>%
    arrange(
        Module,
        Dose
    )

cat("\n============================================================\n")
cat("FIGURE 4D RECONSTRUCTION\n")
cat("============================================================\n\n")

print(
    D_check,
    row.names = FALSE
)

stopifnot(
    nrow(D_check) ==
        nrow(D_locked2)
)

comparison <- merge(
    D_locked2,
    D_check,
    by = c(
        "Module",
        "Dose"
    ),
    suffixes = c(
        "_Locked",
        "_GeneLevel"
    ),
    all = TRUE
)

comparison$Status <- ifelse(
    comparison$Significant_Locked ==
        comparison$Significant_GeneLevel,
    "PASS",
    "FAIL"
)

cat("\nLocked Figure 4D versus S4D reconstruction:\n\n")

print(
    comparison,
    row.names = FALSE
)

stopifnot(
    all(
        comparison$Status ==
            "PASS"
    )
)

cat(
    "\n✓ Every Figure 4D bar is exactly reproduced by S4D gene-level records\n"
)

# ============================================================
# 12. S4A INTERNAL VALIDATION
# ============================================================

stopifnot(
    all(
        S4A$padj < 0.05
    )
)

stopifnot(
    all(
        S4A$overlap >= 2
    )
)

cat(
    "✓ All S4A entries satisfy locked significant GO-BP criteria\n"
)

# ============================================================
# 13. OUTPUT DIRECTORY
# ============================================================

outdir <- file.path(
    fig4dir,
    "supplementary_S4_PUBLICATION"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 14. WORKBOOK FORMATTING
# ============================================================

header_style <- createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE,
    border = "Bottom"
)

wrap_style <- createStyle(
    valign = "top",
    wrapText = TRUE
)

boolean_style <- createStyle(
    halign = "center"
)

add_publication_sheet <- function(
    wb,
    sheet,
    dat
) {

    addWorksheet(
        wb,
        sheet
    )

    writeData(
        wb,
        sheet,
        dat,
        withFilter = TRUE
    )

    addStyle(
        wb,
        sheet,
        header_style,
        rows = 1,
        cols = seq_len(ncol(dat)),
        gridExpand = TRUE
    )

    freezePane(
        wb,
        sheet,
        firstActiveRow = 2,
        firstActiveCol = 2
    )

    setRowHeights(
        wb,
        sheet,
        rows = 1,
        heights = 36
    )

    setColWidths(
        wb,
        sheet,
        cols = seq_len(ncol(dat)),
        widths = "auto"
    )

    long_cols <- grep(
        paste0(
            "Description|description|product|",
            "^genes$|^GO$|KEGG_Pathway|PFAM"
        ),
        names(dat),
        ignore.case = TRUE
    )

    if (length(long_cols) > 0) {

        for (j in long_cols) {

            setColWidths(
                wb,
                sheet,
                cols = j,
                widths = 42
            )

            if (nrow(dat) > 0) {

                addStyle(
                    wb,
                    sheet,
                    wrap_style,
                    rows = 2:(nrow(dat) + 1),
                    cols = j,
                    gridExpand = TRUE
                )
            }
        }
    }

    logical_cols <- which(
        vapply(
            dat,
            is.logical,
            logical(1)
        )
    )

    if (
        length(logical_cols) > 0 &&
        nrow(dat) > 0
    ) {

        addStyle(
            wb,
            sheet,
            boolean_style,
            rows = 2:(nrow(dat) + 1),
            cols = logical_cols,
            gridExpand = TRUE
        )
    }
}

# ============================================================
# 15. CREATE FINAL FOUR-SHEET WORKBOOK
# ============================================================

wb <- createWorkbook()

add_publication_sheet(
    wb,
    "S4A_GO_Enrichment",
    S4A
)

add_publication_sheet(
    wb,
    "S4B_Persistent_Genes",
    S4B
)

add_publication_sheet(
    wb,
    "S4C_Shared_4h_Genes",
    S4C
)

add_publication_sheet(
    wb,
    "S4D_Functional_Modules",
    S4D
)

outfile <- file.path(
    outdir,
    "Supplementary_Table_S4_Figure4_FINAL.xlsx"
)

saveWorkbook(
    wb,
    outfile,
    overwrite = TRUE
)

# ============================================================
# 16. FINAL AUDIT
# ============================================================

cat("\n============================================================\n")
cat("FINAL S4 PUBLICATION WORKBOOK\n")
cat("============================================================\n\n")

cat("Sheets: 4\n")
cat(
    "  S4A_GO_Enrichment       :",
    nrow(S4A),
    "rows\n"
)

cat(
    "  S4B_Persistent_Genes     :",
    nrow(S4B),
    "rows\n"
)

cat(
    "  S4C_Shared_4h_Genes      :",
    nrow(S4C),
    "rows\n"
)

cat(
    "  S4D_Functional_Modules   :",
    nrow(S4D),
    "gene × module rows\n"
)

cat("\nFigure 4D locked counts:\n")

for (i in seq_len(nrow(D_locked2))) {

    cat(
        sprintf(
            "  %-28s | %-4s µg/g | %d genes\n",
            D_locked2$Module[i],
            D_locked2$Dose[i],
            D_locked2$Significant[i]
        )
    )
}

cat("\nWorkbook:\n")
cat(outfile, "\n")

cat("\n============================================================\n")
cat("SUPPLEMENTARY TABLE S4 — FINAL PUBLICATION VERSION LOCKED\n")
cat("============================================================\n")
