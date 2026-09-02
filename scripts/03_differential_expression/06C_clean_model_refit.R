
suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(readr)
    library(tibble)
    library(tidyr)
})

cat("\n============================================================\n")
cat("FA26 STEP 06C — CLEAN EQUIVALENT MODEL REFIT\n")
cat("============================================================\n\n")

dir.create(
    "results/deseq2/clean_refit",
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 1. Load original filtered object only as source of
#    raw counts, metadata, and established size factors
# ============================================================

source_dds <- readRDS(
    "objects/FA26_dds_primary24.rds"
)

count_mat <- counts(
    source_dds,
    normalized = FALSE
)

meta <- as.data.frame(
    colData(source_dds)
)

sf <- sizeFactors(
    source_dds
)

cat("Counts:", dim(count_mat), "\n")
cat("Samples:", nrow(meta), "\n")


# ============================================================
# 2. Clean metadata
# ============================================================

meta$Treatment <- factor(
    as.character(meta$Treatment),
    levels = c(
        "Control",
        "PFBA_0.01",
        "PFBA_1"
    )
)

meta$Time <- factor(
    as.character(meta[["Time.hrs."]]),
    levels = c(
        "2",
        "4"
    )
)

meta$Condition <- factor(
    as.character(meta$Condition),
    levels = c(
        "Control_2h",
        "Control_4h",
        "PFBA_0.01_2h",
        "PFBA_0.01_4h",
        "PFBA_1_2h",
        "PFBA_1_4h"
    )
)

stopifnot(
    identical(
        colnames(count_mat),
        rownames(meta)
    )
)


# ============================================================
# 3. Fresh factorial model
# ============================================================

dds_fac <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = meta,
    design = ~ Treatment * Time
)

sizeFactors(dds_fac) <- sf

dds_fac <- DESeq(
    dds_fac,
    test = "Wald",
    fitType = "local",
    betaPrior = FALSE,
    quiet = FALSE
)

cat("\nFactorial coefficients:\n")
print(resultsNames(dds_fac))


# ============================================================
# 4. Fresh six-condition model
# ============================================================

dds_con <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = meta,
    design = ~ Condition
)

sizeFactors(dds_con) <- sf

dds_con <- DESeq(
    dds_con,
    test = "Wald",
    fitType = "local",
    betaPrior = FALSE,
    quiet = FALSE
)


# ============================================================
# 5. Factorial contrasts
# ============================================================

fac_low2 <- results(
    dds_fac,
    name = "Treatment_PFBA_0.01_vs_Control",
    alpha = 0.05
)

fac_high2 <- results(
    dds_fac,
    name = "Treatment_PFBA_1_vs_Control",
    alpha = 0.05
)

fac_low4 <- results(
    dds_fac,
    contrast = list(
        c(
            "Treatment_PFBA_0.01_vs_Control",
            "TreatmentPFBA_0.01.Time4"
        )
    ),
    alpha = 0.05
)

fac_high4 <- results(
    dds_fac,
    contrast = list(
        c(
            "Treatment_PFBA_1_vs_Control",
            "TreatmentPFBA_1.Time4"
        )
    ),
    alpha = 0.05
)


# ============================================================
# 6. Condition-model contrasts
# ============================================================

con_low2 <- results(
    dds_con,
    contrast = c(
        "Condition",
        "PFBA_0.01_2h",
        "Control_2h"
    ),
    alpha = 0.05
)

con_high2 <- results(
    dds_con,
    contrast = c(
        "Condition",
        "PFBA_1_2h",
        "Control_2h"
    ),
    alpha = 0.05
)

con_low4 <- results(
    dds_con,
    contrast = c(
        "Condition",
        "PFBA_0.01_4h",
        "Control_4h"
    ),
    alpha = 0.05
)

con_high4 <- results(
    dds_con,
    contrast = c(
        "Condition",
        "PFBA_1_4h",
        "Control_4h"
    ),
    alpha = 0.05
)


# ============================================================
# 7. Exact model-equivalence audit
# ============================================================

pairs <- list(

    PFBA_0.01_vs_Control_2h =
        list(fac_low2, con_low2),

    PFBA_1_vs_Control_2h =
        list(fac_high2, con_high2),

    PFBA_0.01_vs_Control_4h =
        list(fac_low4, con_low4),

    PFBA_1_vs_Control_4h =
        list(fac_high4, con_high4)
)


audit <- bind_rows(

    lapply(
        names(pairs),
        function(nm) {

            a <- as.data.frame(
                pairs[[nm]][[1]]
            ) %>%
                rownames_to_column(
                    "gene_id"
                )

            b <- as.data.frame(
                pairs[[nm]][[2]]
            ) %>%
                rownames_to_column(
                    "gene_id"
                )

            m <- a %>%
                dplyr::select(
                    gene_id,
                    fac_LFC = log2FoldChange,
                    fac_p = pvalue,
                    fac_padj = padj
                ) %>%
                inner_join(
                    b %>%
                        dplyr::select(
                            gene_id,
                            con_LFC = log2FoldChange,
                            con_p = pvalue,
                            con_padj = padj
                        ),
                    by = "gene_id"
                )

            write_csv(
                m %>%
                    mutate(
                        abs_LFC_diff =
                            abs(
                                fac_LFC -
                                con_LFC
                            )
                    ) %>%
                    arrange(
                        desc(
                            abs_LFC_diff
                        )
                    ),
                file.path(
                    "results/deseq2/clean_refit",
                    paste0(
                        nm,
                        "_model_equivalence.csv"
                    )
                )
            )

            tibble(
                comparison = nm,

                LFC_correlation =
                    cor(
                        m$fac_LFC,
                        m$con_LFC,
                        use = "complete.obs"
                    ),

                max_abs_LFC_difference =
                    max(
                        abs(
                            m$fac_LFC -
                            m$con_LFC
                        ),
                        na.rm = TRUE
                    ),

                factorial_sig =
                    sum(
                        !is.na(m$fac_padj) &
                        m$fac_padj < 0.05
                    ),

                condition_sig =
                    sum(
                        !is.na(m$con_padj) &
                        m$con_padj < 0.05
                    )
            )
        }
    )
)

write_csv(
    audit,
    "results/deseq2/clean_refit/FA26_clean_model_equivalence_summary.csv"
)


# ============================================================
# 8. Save the clean models
# ============================================================

saveRDS(
    dds_fac,
    "objects/FA26_dds_factorial_clean.rds"
)

saveRDS(
    dds_con,
    "objects/FA26_dds_condition_clean.rds"
)


# ============================================================
# 9. Examine group counts for extreme-discrepancy genes
# ============================================================

problem_genes <- c(
    "LOC118274061",
    "LOC118275628",
    "LOC118279874",
    "LOC118264537",
    "LOC118277132",
    "LOC118273857",
    "LOC118267647",
    "LOC118265473",
    "LOC118277187"
)

problem_genes <- intersect(
    problem_genes,
    rownames(count_mat)
)

problem_counts <- as.data.frame(
    count_mat[
        problem_genes,
        ,
        drop = FALSE
    ]
) %>%
    rownames_to_column(
        "gene_id"
    ) %>%
    pivot_longer(
        -gene_id,
        names_to = "SampleID",
        values_to = "count"
    ) %>%
    left_join(
        meta %>%
            dplyr::select(
                SampleID,
                Treatment,
                Time,
                Condition
            ),
        by = "SampleID"
    ) %>%
    group_by(
        gene_id,
        Condition
    ) %>%
    summarise(
        mean_count = mean(count),
        median_count = median(count),
        min_count = min(count),
        max_count = max(count),
        .groups = "drop"
    )

write_csv(
    problem_counts,
    "results/deseq2/clean_refit/FA26_extreme_LFC_gene_counts.csv"
)


cat("\n============================================================\n")
cat("CLEAN MODEL EQUIVALENCE\n")
cat("============================================================\n\n")

print(
    audit,
    n = nrow(audit),
    width = Inf
)

cat("\nEXTREME-LFC GENE GROUP COUNTS\n\n")

print(
    problem_counts,
    n = nrow(problem_counts),
    width = Inf
)

