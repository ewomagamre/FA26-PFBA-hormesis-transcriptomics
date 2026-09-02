#!/usr/bin/env Rscript

# ============================================================
# FA26 — PUBLICATION SUPPLEMENTARY TABLE S3
#
# FINAL STRUCTURE:
#   S3A_DEGs
#   S3B_Temporal_Signatures
#   S3C_Intersections
#   S3D_Definitions
#
# No DESeq2 analysis is rerun.
# All data derive from the previously validated signature
# master tables and Figure 3 intersection gene sets.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(openxlsx)
})

read_csv <- function(file, ...) {
    read.csv(
        file,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        ...
    )
}

cat("\n============================================================\n")
cat("FIGURE 3 — FINAL PUBLICATION SUPPLEMENTARY TABLE S3\n")
cat("============================================================\n\n")

# ============================================================
# 1. INPUTS
# ============================================================

low_file <- paste0(
    "results/deseq2/signatures_lowdose/",
    "FA26_PFBA001_all_genes_signature_master.csv"
)

high_file <- paste0(
    "results/deseq2/signatures/",
    "FA26_PFBA1_2h_4h_all_genes_signature_master.csv"
)

cross_file <- paste0(
    "results/deseq2/crossdose/",
    "FA26_crossdose_gene_master.csv"
)

figdir <- "results/figure3"

outdir <- file.path(
    figdir,
    "supplementary_S3_PUBLICATION"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

required <- c(
    low_file,
    high_file,
    cross_file,
    file.path(figdir, "Figure3C_crossdose_4h_up.csv"),
    file.path(figdir, "Figure3C_crossdose_4h_down.csv"),
    file.path(figdir, "Figure3C_high2h_low4h_up.csv"),
    file.path(figdir, "Figure3C_high2h_low4h_down.csv"),
    file.path(figdir, "Figure3C_sustained_1ug_up.csv"),
    file.path(figdir, "Figure3C_sustained_1ug_down.csv")
)

missing <- required[!file.exists(required)]

if (length(missing) > 0) {
    cat("Missing files:\n")
    print(missing)
    stop("Required source file(s) missing.")
}

low   <- read_csv(low_file)
high  <- read_csv(high_file)
cross <- read_csv(cross_file)

# ============================================================
# 2. SOURCE VALIDATION
# ============================================================

stopifnot(
    nrow(low) == 11661,
    nrow(high) == 11661,
    nrow(cross) == 11661
)

stopifnot(
    length(unique(low$gene_id)) == 11661,
    length(unique(high$gene_id)) == 11661,
    length(unique(cross$gene_id)) == 11661
)

stopifnot(
    setequal(low$gene_id, high$gene_id),
    setequal(low$gene_id, cross$gene_id)
)

cat("✓ Gene universe validated: 11,661 genes\n")

# ============================================================
# 3. S3A — GENE-LEVEL DEGs
# ============================================================

make_deg_table <- function(
    x,
    dose_label,
    time_label
) {

    if (time_label == "2 h") {

        x %>%
            filter(sig_2h) %>%
            transmute(
                Dose = dose_label,
                Time = time_label,

                gene_id,
                Preferred_name = preferred_name,
                NCBI_product,
                eggNOG_description,

                Direction = direction_2h,

                DESeq2_log2FC = LFC_2h,
                DESeq2_SE = SE_2h,
                Wald_pvalue = pvalue_2h,
                Wald_padj = padj_2h,

                apeglm_log2FC = shrunk_LFC_2h,
                apeglm_SE = shrunk_SE_2h,

                Replicate_support = support_2h,
                Overall_support = overall_support,

                Temporal_class = signature_class,
                Persistence_direction = persistence_direction,

                Atlas_zone = zone,
                UMAP1,
                UMAP2,

                COG_category,
                GO_terms = GOs,
                KEGG_KO = KEGG_ko,
                KEGG_pathway = KEGG_Pathway,
                PFAMs
            )

    } else {

        x %>%
            filter(sig_4h) %>%
            transmute(
                Dose = dose_label,
                Time = time_label,

                gene_id,
                Preferred_name = preferred_name,
                NCBI_product,
                eggNOG_description,

                Direction = direction_4h,

                DESeq2_log2FC = LFC_4h,
                DESeq2_SE = SE_4h,
                Wald_pvalue = pvalue_4h,
                Wald_padj = padj_4h,

                apeglm_log2FC = shrunk_LFC_4h,
                apeglm_SE = shrunk_SE_4h,

                Replicate_support = support_4h,
                Overall_support = overall_support,

                Temporal_class = signature_class,
                Persistence_direction = persistence_direction,

                Atlas_zone = zone,
                UMAP1,
                UMAP2,

                COG_category,
                GO_terms = GOs,
                KEGG_KO = KEGG_ko,
                KEGG_pathway = KEGG_Pathway,
                PFAMs
            )
    }
}

low2 <- make_deg_table(
    low,
    "0.01 µg/g",
    "2 h"
)

low4 <- make_deg_table(
    low,
    "0.01 µg/g",
    "4 h"
)

high2 <- make_deg_table(
    high,
    "1 µg/g",
    "2 h"
)

high4 <- make_deg_table(
    high,
    "1 µg/g",
    "4 h"
)

S3A <- bind_rows(
    low2,
    low4,
    high2,
    high4
)

stopifnot(
    nrow(low2) == 10,
    nrow(low4) == 587,
    nrow(high2) == 522,
    nrow(high4) == 1565,
    nrow(S3A) == 2684
)

stopifnot(
    sum(low2$Direction == "Up") == 10,
    sum(low2$Direction == "Down") == 0,

    sum(low4$Direction == "Up") == 50,
    sum(low4$Direction == "Down") == 537,

    sum(high2$Direction == "Up") == 431,
    sum(high2$Direction == "Down") == 91,

    sum(high4$Direction == "Up") == 437,
    sum(high4$Direction == "Down") == 1128
)

cat("✓ S3A validated: 2,684 gene × comparison records\n")

# ============================================================
# 4. S3B — TEMPORAL SIGNATURES
# ============================================================

make_temporal <- function(
    x,
    dose_label
) {

    x %>%
        filter(
            signature_class != "Not significant"
        ) %>%
        transmute(
            Dose = dose_label,

            gene_id,
            Preferred_name = preferred_name,
            NCBI_product,
            eggNOG_description,

            Temporal_class = signature_class,
            Persistence_direction = persistence_direction,

            Significant_2h = sig_2h,
            Direction_2h = direction_2h,
            DESeq2_log2FC_2h = LFC_2h,
            DESeq2_SE_2h = SE_2h,
            Wald_pvalue_2h = pvalue_2h,
            Wald_padj_2h = padj_2h,
            apeglm_log2FC_2h = shrunk_LFC_2h,
            apeglm_SE_2h = shrunk_SE_2h,

            Significant_4h = sig_4h,
            Direction_4h = direction_4h,
            DESeq2_log2FC_4h = LFC_4h,
            DESeq2_SE_4h = SE_4h,
            Wald_pvalue_4h = pvalue_4h,
            Wald_padj_4h = padj_4h,
            apeglm_log2FC_4h = shrunk_LFC_4h,
            apeglm_SE_4h = shrunk_SE_4h,

            Replicate_support_2h = support_2h,
            Replicate_support_4h = support_4h,
            Overall_support = overall_support,

            Atlas_zone = zone,
            UMAP1,
            UMAP2,

            COG_category,
            GO_terms = GOs,
            KEGG_KO = KEGG_ko,
            KEGG_pathway = KEGG_Pathway,
            PFAMs
        )
}

S3B <- bind_rows(
    make_temporal(
        low,
        "0.01 µg/g"
    ),
    make_temporal(
        high,
        "1 µg/g"
    )
)

stopifnot(
    sum(low$signature_class == "Early-only") == 9,
    sum(low$signature_class == "Persistent") == 1,
    sum(low$signature_class == "Later-emerging") == 586,

    sum(high$signature_class == "Early-only") == 384,
    sum(high$signature_class == "Persistent") == 138,
    sum(high$signature_class == "Later-emerging") == 1427,

    nrow(S3B) == 2545
)

cat("✓ S3B validated: 2,545 temporal-signature records\n")

# ============================================================
# 5. FIGURE 3C GENE SETS
# ============================================================

read_gene_set <- function(filename) {

    read_csv(
        file.path(figdir, filename)
    )$gene_id
}

sustained_up <- read_gene_set(
    "Figure3C_sustained_1ug_up.csv"
)

sustained_down <- read_gene_set(
    "Figure3C_sustained_1ug_down.csv"
)

cross4_up <- read_gene_set(
    "Figure3C_crossdose_4h_up.csv"
)

cross4_down <- read_gene_set(
    "Figure3C_crossdose_4h_down.csv"
)

high2low4_up <- read_gene_set(
    "Figure3C_high2h_low4h_up.csv"
)

high2low4_down <- read_gene_set(
    "Figure3C_high2h_low4h_down.csv"
)

stopifnot(
    length(sustained_up) == 123,
    length(sustained_down) == 12,

    length(cross4_up) == 14,
    length(cross4_down) == 8,

    length(high2low4_up) == 2,
    length(high2low4_down) == 0
)

# ============================================================
# 6. S3C — ONE ROW PER UNIQUE INTERSECTION GENE
# ============================================================

intersection_genes <- unique(
    c(
        sustained_up,
        sustained_down,
        cross4_up,
        cross4_down,
        high2low4_up,
        high2low4_down
    )
)

S3C <- data.frame(
    gene_id = intersection_genes,
    stringsAsFactors = FALSE
) %>%

    mutate(

        # ----------------------------------------------------
        # EXACT FIGURE 3C RELATIONSHIPS
        # ----------------------------------------------------

        Sustained_1ug_2to4h =
            gene_id %in%
            c(
                sustained_up,
                sustained_down
            ),

        Sustained_direction =
            case_when(
                gene_id %in% sustained_up ~
                    "Upregulated",

                gene_id %in% sustained_down ~
                    "Downregulated",

                TRUE ~
                    NA_character_
            ),

        CrossDose_Shared_4h =
            gene_id %in%
            c(
                cross4_up,
                cross4_down
            ),

        CrossDose_4h_direction =
            case_when(
                gene_id %in% cross4_up ~
                    "Upregulated",

                gene_id %in% cross4_down ~
                    "Downregulated",

                TRUE ~
                    NA_character_
            ),

        High1ug_2h_Low001ug_4h_Shared =
            gene_id %in%
            c(
                high2low4_up,
                high2low4_down
            ),

        High2h_Low4h_direction =
            case_when(
                gene_id %in% high2low4_up ~
                    "Upregulated",

                gene_id %in% high2low4_down ~
                    "Downregulated",

                TRUE ~
                    NA_character_
            )
    )

# ============================================================
# 7. ADD GENE NAMES / ANNOTATION
# ============================================================

annotation <- cross %>%
    transmute(
        gene_id,

        Preferred_name = preferred_name,
        NCBI_product,
        eggNOG_description,

        Atlas_zone = zone
    ) %>%
    distinct(
        gene_id,
        .keep_all = TRUE
    )

S3C <- S3C %>%
    left_join(
        annotation,
        by = "gene_id"
    )

# ============================================================
# 8. ADD HIGH-DOSE STATISTICS
# ============================================================

high_stats <- high %>%
    transmute(
        gene_id,

        High_sig_2h = sig_2h,
        High_direction_2h = direction_2h,
        High_DESeq2_log2FC_2h = LFC_2h,
        High_Wald_padj_2h = padj_2h,
        High_apeglm_log2FC_2h = shrunk_LFC_2h,

        High_sig_4h = sig_4h,
        High_direction_4h = direction_4h,
        High_DESeq2_log2FC_4h = LFC_4h,
        High_Wald_padj_4h = padj_4h,
        High_apeglm_log2FC_4h = shrunk_LFC_4h
    )

S3C <- S3C %>%
    left_join(
        high_stats,
        by = "gene_id"
    )

# ============================================================
# 9. ADD LOW-DOSE STATISTICS
# ============================================================

low_stats <- low %>%
    transmute(
        gene_id,

        Low_sig_2h = sig_2h,
        Low_direction_2h = direction_2h,
        Low_DESeq2_log2FC_2h = LFC_2h,
        Low_Wald_padj_2h = padj_2h,
        Low_apeglm_log2FC_2h = shrunk_LFC_2h,

        Low_sig_4h = sig_4h,
        Low_direction_4h = direction_4h,
        Low_DESeq2_log2FC_4h = LFC_4h,
        Low_Wald_padj_4h = padj_4h,
        Low_apeglm_log2FC_4h = shrunk_LFC_4h
    )

S3C <- S3C %>%
    left_join(
        low_stats,
        by = "gene_id"
    )

# Put the biologically useful columns first
S3C <- S3C %>%
    select(
        gene_id,
        Preferred_name,
        NCBI_product,
        eggNOG_description,
        Atlas_zone,

        Sustained_1ug_2to4h,
        Sustained_direction,

        CrossDose_Shared_4h,
        CrossDose_4h_direction,

        High1ug_2h_Low001ug_4h_Shared,
        High2h_Low4h_direction,

        everything()
    )

# ============================================================
# 10. VALIDATE THAT S3C RECONSTRUCTS FIGURE 3C
# ============================================================

stopifnot(

    sum(
        S3C$Sustained_1ug_2to4h &
        S3C$Sustained_direction ==
            "Upregulated",
        na.rm = TRUE
    ) == 123,

    sum(
        S3C$Sustained_1ug_2to4h &
        S3C$Sustained_direction ==
            "Downregulated",
        na.rm = TRUE
    ) == 12,

    sum(
        S3C$CrossDose_Shared_4h &
        S3C$CrossDose_4h_direction ==
            "Upregulated",
        na.rm = TRUE
    ) == 14,

    sum(
        S3C$CrossDose_Shared_4h &
        S3C$CrossDose_4h_direction ==
            "Downregulated",
        na.rm = TRUE
    ) == 8,

    sum(
        S3C$High1ug_2h_Low001ug_4h_Shared &
        S3C$High2h_Low4h_direction ==
            "Upregulated",
        na.rm = TRUE
    ) == 2,

    sum(
        S3C$High1ug_2h_Low001ug_4h_Shared &
        S3C$High2h_Low4h_direction ==
            "Downregulated",
        na.rm = TRUE
    ) == 0
)

cat("✓ S3C reproduces Figure 3C exactly\n")

cat("\nFigure 3C:\n")
cat("  Sustained 1 µg/g (2 → 4 h): 123 Up / 12 Down\n")
cat("  Cross-dose shared (4 h):     14 Up / 8 Down\n")
cat("  High 2 h ↔ Low 4 h:           2 Up / 0 Down\n")
cat("  Unique genes represented:    ", nrow(S3C), "\n")

# ============================================================
# 11. S3D — DEFINITIONS
# ============================================================

S3D <- data.frame(

    Parameter = c(
        "Gene universe",
        "DEG significance",
        "Direction",
        "DESeq2 log2FC",
        "apeglm log2FC",
        "Role of shrinkage",

        "Early-only",
        "Persistent",
        "Later-emerging",

        "Sustained 1 µg/g (2 → 4 h)",
        "Cross-dose shared (4 h)",
        "1 µg/g 2 h ↔ 0.01 µg/g 4 h",

        "Replicate support",
        "Atlas zone"
    ),

    Definition = c(
        "11,661 genes retained in the DESeq2 analysis.",

        paste(
            "Adjusted Wald P value (padj) < 0.05 from",
            "the unshrunken condition-model DESeq2 result."
        ),

        paste(
            "Up or Down assigned from the sign of the",
            "unshrunken DESeq2 log2 fold change among",
            "significant genes."
        ),

        paste(
            "Unshrunken condition-model DESeq2 log2 fold",
            "change corresponding to the Wald test."
        ),

        paste(
            "apeglm-shrunken log2 fold-change estimate",
            "obtained from a relevelled condition-model",
            "coefficient."
        ),

        paste(
            "Shrinkage was used to stabilize effect-size",
            "estimation, not to determine DEG significance."
        ),

        "Significant at 2 h but not at 4 h.",

        "Significant at both 2 h and 4 h.",

        "Not significant at 2 h but significant at 4 h.",

        paste(
            "Genes significant at both 2 h and 4 h under",
            "the 1 µg/g PFBA treatment, separated by",
            "concordant direction."
        ),

        paste(
            "Genes significant at 4 h for both PFBA doses,",
            "separated by concordant response direction."
        ),

        paste(
            "Genes significant at 1 µg/g PFBA at 2 h and",
            "0.01 µg/g PFBA at 4 h, separated by",
            "concordant response direction."
        ),

        paste(
            "Normalized-count replicate support;",
            "not an additional significance criterion."
        ),

        paste(
            "Atlas-zone assignment is available only for",
            "genes represented in the selected 3,000-gene",
            "transcriptomic atlas."
        )
    ),

    stringsAsFactors = FALSE
)

# ============================================================
# 12. WORKBOOK FORMATTING
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

add_sheet <- function(
    wb,
    sheet_name,
    dat
) {

    addWorksheet(
        wb,
        sheet_name
    )

    writeData(
        wb,
        sheet_name,
        dat,
        withFilter = TRUE
    )

    addStyle(
        wb,
        sheet_name,
        header_style,
        rows = 1,
        cols = seq_len(ncol(dat)),
        gridExpand = TRUE
    )

    freezePane(
        wb,
        sheet_name,
        firstActiveRow = 2,
        firstActiveCol = 2
    )

    setRowHeights(
        wb,
        sheet_name,
        rows = 1,
        heights = 35
    )

    setColWidths(
        wb,
        sheet_name,
        cols = seq_len(ncol(dat)),
        widths = "auto"
    )

    long_cols <- grep(
        paste(
            "product|description|GO_terms|",
            "KEGG_pathway|PFAM|Definition",
            sep = ""
        ),
        names(dat),
        ignore.case = TRUE
    )

    if (length(long_cols) > 0) {

        for (j in long_cols) {

            setColWidths(
                wb,
                sheet_name,
                cols = j,
                widths = 38
            )

            if (nrow(dat) > 0) {

                addStyle(
                    wb,
                    sheet_name,
                    wrap_style,
                    rows = 2:(nrow(dat) + 1),
                    cols = j,
                    gridExpand = TRUE
                )
            }
        }
    }
}

# ============================================================
# 13. WRITE FINAL FOUR-SHEET WORKBOOK
# ============================================================

wb <- createWorkbook()

add_sheet(
    wb,
    "S3A_DEGs",
    S3A
)

add_sheet(
    wb,
    "S3B_Temporal_Signatures",
    S3B
)

add_sheet(
    wb,
    "S3C_Intersections",
    S3C
)

add_sheet(
    wb,
    "S3D_Definitions",
    S3D
)

outfile <- file.path(
    outdir,
    "Supplementary_Table_S3_Figure3_FINAL.xlsx"
)

saveWorkbook(
    wb,
    outfile,
    overwrite = TRUE
)

# ============================================================
# 14. FINAL CONSOLE AUDIT
# ============================================================

cat("\n============================================================\n")
cat("FINAL S3 PUBLICATION WORKBOOK\n")
cat("============================================================\n\n")

cat("Sheets: 4\n")
cat("  S3A_DEGs                  :", nrow(S3A), "rows\n")
cat("  S3B_Temporal_Signatures   :", nrow(S3B), "rows\n")
cat("  S3C_Intersections         :", nrow(S3C), "unique genes\n")
cat("  S3D_Definitions           :", nrow(S3D), "rows\n\n")

cat("FIGURE 3A VERIFIED:\n")
cat("  Low 2 h  = 10 Up / 0 Down\n")
cat("  Low 4 h  = 50 Up / 537 Down\n")
cat("  High 2 h = 431 Up / 91 Down\n")
cat("  High 4 h = 437 Up / 1128 Down\n\n")

cat("FIGURE 3B VERIFIED:\n")
cat("  Low  = 9 Early / 1 Persistent / 586 Later\n")
cat("  High = 384 Early / 138 Persistent / 1427 Later\n\n")

cat("FIGURE 3C VERIFIED:\n")
cat("  Sustained high = 123 Up / 12 Down\n")
cat("  Cross-dose 4 h = 14 Up / 8 Down\n")
cat("  High2h ↔ Low4h = 2 Up / 0 Down\n\n")

cat("Workbook:\n")
cat(outfile, "\n")

cat("\n============================================================\n")
cat("SUPPLEMENTARY TABLE S3 — FINAL PUBLICATION VERSION LOCKED\n")
cat("============================================================\n")
