# ============================================================
# FA26 — Figure 2 Source Data
# Supplementary Table S2
#
# Outputs:
#   S2A.xlsx — atlas genes and expression data
#   S2B.xlsx — gene-level temporal dynamics
#   S2C.xlsx — global temporal movement summary
#   S2D.xlsx — atlas optimization
#   Combined S2 workbook with all four sheets
# ============================================================

suppressPackageStartupMessages({
    library(openxlsx)
})

cat("\n============================================================\n")
cat("FIGURE 2 — SUPPLEMENTARY TABLE S2\n")
cat("============================================================\n\n")

outdir <- "results/figure2/supplementary_S2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. INPUT FILES
# ============================================================

master_file <-
    "results/gene_atlas/FA26_gene_atlas_master_unannotated.csv"

variance_file <-
    "results/gene_variance/FA26_top3000_variable_genes.csv"

optimization_file <-
    "results/gene_variance/FA26_feature_selection_optimization_summary.csv"

stored_summary_file <-
    "results/figure2/Figure2_global_temporal_movement.csv"

stopifnot(
    file.exists(master_file),
    file.exists(variance_file),
    file.exists(optimization_file),
    file.exists(stored_summary_file)
)

master <- read.csv(master_file, check.names = FALSE)
variance <- read.csv(variance_file, check.names = FALSE)
optimization <- read.csv(optimization_file, check.names = FALSE)
stored_summary <- read.csv(stored_summary_file, check.names = FALSE)

stopifnot(nrow(master) == 3000)
stopifnot(nrow(variance) == 3000)

# ============================================================
# 2. S2A — ATLAS GENES
# ============================================================

S2A <- merge(
    variance[, c(
        "gene_id",
        "variance",
        "rank",
        "cumulative_variance",
        "fraction_of_genes"
    )],
    master,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

S2A <- S2A[order(S2A$rank), ]

stopifnot(nrow(S2A) == 3000)
stopifnot(!any(is.na(S2A$UMAP1)))
stopifnot(!any(is.na(S2A$UMAP2)))

# ============================================================
# 3. S2B — GENE-LEVEL TEMPORAL DYNAMICS
# ============================================================

S2B <- master[, c(
    "gene_id",
    "UMAP1",
    "UMAP2",
    "Control_2h",
    "Control_4h",
    "PFBA_0.01_2h",
    "PFBA_0.01_4h",
    "PFBA_1_2h",
    "PFBA_1_4h"
)]

S2B$DeltaVST_Control <-
    S2B$Control_4h - S2B$Control_2h

S2B$Abs_DeltaVST_Control <-
    abs(S2B$DeltaVST_Control)

S2B$DeltaVST_PFBA_0.01 <-
    S2B$PFBA_0.01_4h - S2B$PFBA_0.01_2h

S2B$Abs_DeltaVST_PFBA_0.01 <-
    abs(S2B$DeltaVST_PFBA_0.01)

S2B$DeltaVST_PFBA_1 <-
    S2B$PFBA_1_4h - S2B$PFBA_1_2h

S2B$Abs_DeltaVST_PFBA_1 <-
    abs(S2B$DeltaVST_PFBA_1)

# ============================================================
# 4. S2C — GLOBAL TEMPORAL MOVEMENT SUMMARY
# ============================================================

summarize_delta <- function(v, treatment) {

    data.frame(
        Treatment = treatment,
        Mean_DeltaVST = mean(v),
        Mean_Abs_DeltaVST = mean(abs(v)),
        Median_Abs_DeltaVST = median(abs(v)),
        IQR_Abs_DeltaVST = IQR(abs(v)),
        Percent_Abs_DeltaVST_gt_0.5 =
            100 * mean(abs(v) > 0.5),
        Percent_Abs_DeltaVST_gt_1 =
            100 * mean(abs(v) > 1),
        Percent_Abs_DeltaVST_gt_2 =
            100 * mean(abs(v) > 2)
    )
}

S2C <- rbind(
    summarize_delta(
        S2B$DeltaVST_Control,
        "Control"
    ),
    summarize_delta(
        S2B$DeltaVST_PFBA_0.01,
        "PFBA 0.01 ug/g"
    ),
    summarize_delta(
        S2B$DeltaVST_PFBA_1,
        "PFBA 1 ug/g"
    )
)

# ============================================================
# 5. VERIFY FIGURE 2 TEMPORAL SUMMARY
# ============================================================

check_cols <- c(
    "mean_delta",
    "mean_abs_delta",
    "median_abs_delta",
    "IQR_abs_delta",
    "pct_abs_gt_0.5",
    "pct_abs_gt_1",
    "pct_abs_gt_2"
)

calc_check <- data.frame(
    mean_delta = S2C$Mean_DeltaVST,
    mean_abs_delta = S2C$Mean_Abs_DeltaVST,
    median_abs_delta = S2C$Median_Abs_DeltaVST,
    IQR_abs_delta = S2C$IQR_Abs_DeltaVST,
    pct_abs_gt_0.5 = S2C$Percent_Abs_DeltaVST_gt_0.5,
    pct_abs_gt_1 = S2C$Percent_Abs_DeltaVST_gt_1,
    pct_abs_gt_2 = S2C$Percent_Abs_DeltaVST_gt_2
)

max_diff <- max(
    abs(
        as.matrix(calc_check[, check_cols]) -
        as.matrix(stored_summary[, check_cols])
    )
)

cat(
    "Maximum difference from stored Figure 2 summary:",
    format(max_diff, scientific = TRUE),
    "\n"
)

stopifnot(max_diff < 1e-10)

# ============================================================
# 6. S2D — FEATURE-SELECTION OPTIMIZATION
# ============================================================

S2D <- optimization

S2D$UMAP_n_neighbors <- 30
S2D$UMAP_min_dist <- 0.30
S2D$UMAP_metric <- "euclidean"
S2D$UMAP_seed <- 260725

# ============================================================
# 7. EXCEL FORMATTING FUNCTION
# ============================================================

headerStyle <- createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE,
    border = "Bottom"
)

save_single_table <- function(data, filename, sheetname) {

    wb <- createWorkbook()

    addWorksheet(wb, sheetname)

    writeData(
        wb,
        sheetname,
        data,
        headerStyle = headerStyle,
        withFilter = TRUE
    )

    freezePane(
        wb,
        sheetname,
        firstRow = TRUE
    )

    setColWidths(
        wb,
        sheetname,
        cols = 1:ncol(data),
        widths = "auto"
    )

    saveWorkbook(
        wb,
        filename,
        overwrite = TRUE
    )
}

# ============================================================
# 8. WRITE INDIVIDUAL XLSX FILES
# ============================================================

save_single_table(
    S2A,
    file.path(
        outdir,
        "Supplementary_Table_S2A_Atlas_Genes.xlsx"
    ),
    "S2A_Atlas_Genes"
)

save_single_table(
    S2B,
    file.path(
        outdir,
        "Supplementary_Table_S2B_Temporal_Dynamics.xlsx"
    ),
    "S2B_Temporal_Dynamics"
)

save_single_table(
    S2C,
    file.path(
        outdir,
        "Supplementary_Table_S2C_Global_Temporal_Summary.xlsx"
    ),
    "S2C_Global_Summary"
)

save_single_table(
    S2D,
    file.path(
        outdir,
        "Supplementary_Table_S2D_Atlas_Optimization.xlsx"
    ),
    "S2D_Atlas_Optimization"
)

# ============================================================
# 9. COMBINED SUPPLEMENTARY TABLE S2 WORKBOOK
# ============================================================

wb <- createWorkbook()

add_combined_sheet <- function(wb, sheet_name, dat) {

    addWorksheet(wb, sheet_name)

    writeData(
        wb,
        sheet_name,
        dat,
        headerStyle = headerStyle,
        withFilter = TRUE
    )

    freezePane(
        wb,
        sheet_name,
        firstRow = TRUE
    )

    setColWidths(
        wb,
        sheet_name,
        cols = 1:ncol(dat),
        widths = "auto"
    )
}

add_combined_sheet(
    wb,
    "S2A_Atlas_Genes",
    S2A
)

add_combined_sheet(
    wb,
    "S2B_Temporal_Dynamics",
    S2B
)

add_combined_sheet(
    wb,
    "S2C_Global_Summary",
    S2C
)

add_combined_sheet(
    wb,
    "S2D_Atlas_Optimization",
    S2D
)

combined_file <- file.path(
    outdir,
    "Supplementary_Table_S2_Figure2_Transcriptomic_Atlas.xlsx"
)

saveWorkbook(
    wb,
    combined_file,
    overwrite = TRUE
)

# ============================================================
# 10. FINAL AUDIT
# ============================================================

cat("\n============================================================\n")
cat("S2 TABLE DIMENSIONS\n")
cat("============================================================\n")

cat("S2A — Atlas genes:        ", nrow(S2A), "rows\n")
cat("S2B — Temporal dynamics:  ", nrow(S2B), "rows\n")
cat("S2C — Global summary:     ", nrow(S2C), "rows\n")
cat("S2D — Optimization:       ", nrow(S2D), "rows\n")

cat("\nSelected atlas solution:\n")

print(
    S2D[
        S2D$selected,
        c(
            "candidate_size",
            "fraction_retained_genes",
            "cumulative_variance_pct",
            "procrustes_r",
            "permutation_p",
            "selected"
        )
    ],
    row.names = FALSE
)

cat("\nTemporal movement summary:\n")
print(S2C, row.names = FALSE)

cat("\nFiles written to:\n")
cat(outdir, "\n\n")

print(
    list.files(
        outdir,
        pattern = "\\.xlsx$",
        full.names = TRUE
    )
)

cat("\n============================================================\n")
cat("SUPPLEMENTARY TABLE S2 COMPLETE\n")
cat("============================================================\n")
