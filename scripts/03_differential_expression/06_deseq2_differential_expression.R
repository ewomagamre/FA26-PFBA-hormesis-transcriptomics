
# ============================================================
# FA26 primary 24-sample analysis
# Step 06 — Differential expression with DESeq2
#
# Dataset:
#   11,661 filtered genes
#   24 libraries
#   Control, PFBA 0.01 ug/g, PFBA 1 ug/g
#   2 h and 4 h
#   n = 4 per treatment-time condition
#
# Model:
#   ~ Treatment * Time
#
# References:
#   Treatment = Control
#   Time      = 2 h
#
# This analysis is independent of the 3,000-gene UMAP atlas.
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(readr)
    library(tibble)
})

cat("\n============================================================\n")
cat("FA26 STEP 06 — DESeq2 DIFFERENTIAL EXPRESSION\n")
cat("============================================================\n\n")


# ============================================================
# 1. Paths
# ============================================================

input_file <- "objects/FA26_dds_primary24.rds"

result_dir <- "results/deseq2"
object_dir <- "objects"

dir.create(
    result_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    object_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 2. Load filtered DESeqDataSet
# ============================================================

dds <- readRDS(input_file)

cat("Input dimensions:\n")
print(dim(dds))

stopifnot(
    nrow(dds) == 11661,
    ncol(dds) == 24
)


# ============================================================
# 3. Construct clean factorial metadata
# ============================================================

# Explicit treatment factor
dds$Treatment <- factor(
    as.character(dds$Treatment),
    levels = c(
        "Control",
        "PFBA_0.01",
        "PFBA_1"
    )
)

# Explicit time factor from numeric Time(hrs)
dds$Time <- factor(
    as.character(
        colData(dds)[["Time(hrs)"]]
    ),
    levels = c(
        "2",
        "4"
    )
)

stopifnot(
    !any(is.na(dds$Treatment)),
    !any(is.na(dds$Time))
)

cat("\nTreatment levels:\n")
print(levels(dds$Treatment))

cat("\nTime levels:\n")
print(levels(dds$Time))

cat("\nSample counts by Treatment x Time:\n")
print(
    table(
        dds$Treatment,
        dds$Time
    )
)


# ============================================================
# 4. Set factorial design
# ============================================================

design(dds) <- ~ Treatment * Time

cat("\nFactorial design:\n")
print(design(dds))


# ============================================================
# 5. Fit Wald model
# ============================================================

set.seed(42)

dds_wald <- DESeq(
    dds,
    test = "Wald",
    quiet = FALSE
)

cat("\nDESeq2 coefficient names:\n")

rn <- resultsNames(dds_wald)
print(rn)


# Save fitted model immediately
saveRDS(
    dds_wald,
    file.path(
        object_dir,
        "FA26_dds_factorial_wald.rds"
    )
)


# ============================================================
# 6. Identify coefficient names robustly
# ============================================================

find_coef <- function(pattern) {

    hit <- grep(
        pattern,
        rn,
        value = TRUE
    )

    if (length(hit) != 1) {

        stop(
            paste0(
                "Expected one coefficient matching: ",
                pattern,
                "\nFound: ",
                paste(hit, collapse = ", ")
            )
        )
    }

    hit
}


coef_low_2h <- find_coef(
    "^Treatment_PFBA_0\\.01_vs_Control$"
)

coef_high_2h <- find_coef(
    "^Treatment_PFBA_1_vs_Control$"
)

coef_time_control <- find_coef(
    "^Time_4_vs_2$"
)


# Interaction coefficient names vary slightly between DESeq2
# versions, so identify them by both treatment and Time4.

interaction_names <- rn[
    grepl(
        "Treatment.*Time|Time.*Treatment",
        rn
    )
]

cat("\nInteraction coefficients:\n")
print(interaction_names)


coef_low_interaction <- interaction_names[
    grepl(
        "PFBA_0\\.01",
        interaction_names
    )
]

coef_high_interaction <- interaction_names[
    grepl(
        "PFBA_1",
        interaction_names
    )
]

if (length(coef_low_interaction) != 1) {
    stop("Could not uniquely identify low-dose interaction coefficient.")
}

if (length(coef_high_interaction) != 1) {
    stop("Could not uniquely identify high-dose interaction coefficient.")
}


cat("\nResolved coefficients:\n")

cat(
    "Low dose at 2 h:      ",
    coef_low_2h,
    "\n"
)

cat(
    "High dose at 2 h:     ",
    coef_high_2h,
    "\n"
)

cat(
    "Control 4 h vs 2 h:   ",
    coef_time_control,
    "\n"
)

cat(
    "Low-dose interaction: ",
    coef_low_interaction,
    "\n"
)

cat(
    "High-dose interaction:",
    coef_high_interaction,
    "\n"
)


# ============================================================
# 7. Helper to extract and save results
# ============================================================

extract_result <- function(
    dds_object,
    contrast_spec,
    comparison_name,
    alpha = 0.05
) {

    if (
        is.character(contrast_spec) &&
        length(contrast_spec) == 1
    ) {

        res <- results(
            dds_object,
            name = contrast_spec,
            alpha = alpha
        )

    } else {

        res <- results(
            dds_object,
            contrast = contrast_spec,
            alpha = alpha
        )
    }


    df <- as.data.frame(res) %>%
        rownames_to_column(
            "gene_id"
        ) %>%
        arrange(
            padj
        )

    df <- df %>%
        mutate(
            comparison = comparison_name,
            significant =
                !is.na(padj) &
                padj < alpha,
            direction =
                case_when(
                    significant &
                    log2FoldChange > 0 ~ "Up",
                    significant &
                    log2FoldChange < 0 ~ "Down",
                    TRUE ~ "Not significant"
                )
        )

    write_csv(
        df,
        file.path(
            result_dir,
            paste0(
                comparison_name,
                ".csv"
            )
        )
    )

    df
}


# ============================================================
# 8. Primary treatment-vs-control contrasts
# ============================================================

# ------------------------------------------------------------
# 8A. PFBA 0.01 vs Control at 2 h
# Main treatment coefficient because 2 h is reference time.
# ------------------------------------------------------------

res_low_2h <- extract_result(
    dds_wald,
    coef_low_2h,
    "PFBA_0.01_vs_Control_2h"
)


# ------------------------------------------------------------
# 8B. PFBA 1 vs Control at 2 h
# ------------------------------------------------------------

res_high_2h <- extract_result(
    dds_wald,
    coef_high_2h,
    "PFBA_1_vs_Control_2h"
)


# ------------------------------------------------------------
# 8C. PFBA 0.01 vs Control at 4 h
#
# treatment main effect + treatment:time interaction
# ------------------------------------------------------------

res_low_4h <- extract_result(
    dds_wald,
    list(
        c(
            coef_low_2h,
            coef_low_interaction
        )
    ),
    "PFBA_0.01_vs_Control_4h"
)


# ------------------------------------------------------------
# 8D. PFBA 1 vs Control at 4 h
# ------------------------------------------------------------

res_high_4h <- extract_result(
    dds_wald,
    list(
        c(
            coef_high_2h,
            coef_high_interaction
        )
    ),
    "PFBA_1_vs_Control_4h"
)


# ============================================================
# 9. Temporal contrasts within each treatment
# ============================================================

# ------------------------------------------------------------
# 9A. Control: 4 h vs 2 h
# ------------------------------------------------------------

res_control_time <- extract_result(
    dds_wald,
    coef_time_control,
    "Control_4h_vs_2h"
)


# ------------------------------------------------------------
# 9B. PFBA 0.01: 4 h vs 2 h
#
# time main effect + low-dose interaction
# ------------------------------------------------------------

res_low_time <- extract_result(
    dds_wald,
    list(
        c(
            coef_time_control,
            coef_low_interaction
        )
    ),
    "PFBA_0.01_4h_vs_2h"
)


# ------------------------------------------------------------
# 9C. PFBA 1: 4 h vs 2 h
# ------------------------------------------------------------

res_high_time <- extract_result(
    dds_wald,
    list(
        c(
            coef_time_control,
            coef_high_interaction
        )
    ),
    "PFBA_1_4h_vs_2h"
)


# ============================================================
# 10. Interaction contrasts
#
# Difference-in-differences:
#
# [PFBA 4h - PFBA 2h] -
# [Control 4h - Control 2h]
#
# These directly test whether temporal behavior differs
# between PFBA and Control.
# ============================================================

res_low_interaction <- extract_result(
    dds_wald,
    coef_low_interaction,
    "Interaction_PFBA_0.01_vs_Control"
)

res_high_interaction <- extract_result(
    dds_wald,
    coef_high_interaction,
    "Interaction_PFBA_1_vs_Control"
)


# ============================================================
# 11. Omnibus interaction LRT
#
# Full:
#   ~ Treatment * Time
#
# Reduced:
#   ~ Treatment + Time
#
# Tests whether ANY treatment-specific temporal response
# improves model fit for each gene.
# ============================================================

dds_lrt <- dds

design(dds_lrt) <- ~ Treatment * Time

dds_lrt <- DESeq(
    dds_lrt,
    test = "LRT",
    reduced = ~ Treatment + Time,
    quiet = FALSE
)

res_lrt <- results(
    dds_lrt,
    alpha = 0.05
)

res_lrt_df <- as.data.frame(
    res_lrt
) %>%
    rownames_to_column(
        "gene_id"
    ) %>%
    arrange(
        padj
    ) %>%
    mutate(
        comparison =
            "Omnibus_Treatment_Time_Interaction",
        significant =
            !is.na(padj) &
            padj < 0.05
    )

write_csv(
    res_lrt_df,
    file.path(
        result_dir,
        "Omnibus_Treatment_Time_Interaction.csv"
    )
)

saveRDS(
    dds_lrt,
    file.path(
        object_dir,
        "FA26_dds_factorial_interaction_LRT.rds"
    )
)


# ============================================================
# 12. Summarize DEG counts
#
# Significance:
#   BH-adjusted P < 0.05
#
# No LFC threshold imposed at this stage.
# ============================================================

contrast_results <- list(

    "PFBA_0.01_vs_Control_2h" =
        res_low_2h,

    "PFBA_1_vs_Control_2h" =
        res_high_2h,

    "PFBA_0.01_vs_Control_4h" =
        res_low_4h,

    "PFBA_1_vs_Control_4h" =
        res_high_4h,

    "Control_4h_vs_2h" =
        res_control_time,

    "PFBA_0.01_4h_vs_2h" =
        res_low_time,

    "PFBA_1_4h_vs_2h" =
        res_high_time,

    "Interaction_PFBA_0.01_vs_Control" =
        res_low_interaction,

    "Interaction_PFBA_1_vs_Control" =
        res_high_interaction
)


summary_table <- bind_rows(

    lapply(
        names(contrast_results),
        function(nm) {

            df <-
                contrast_results[[nm]]

            tibble(
                comparison = nm,

                tested_genes =
                    sum(
                        !is.na(
                            df$pvalue
                        )
                    ),

                padj_available =
                    sum(
                        !is.na(
                            df$padj
                        )
                    ),

                significant =
                    sum(
                        df$significant,
                        na.rm = TRUE
                    ),

                up =
                    sum(
                        df$direction == "Up",
                        na.rm = TRUE
                    ),

                down =
                    sum(
                        df$direction == "Down",
                        na.rm = TRUE
                    )
            )
        }
    )
)


summary_table <- bind_rows(

    summary_table,

    tibble(
        comparison =
            "Omnibus_Treatment_Time_Interaction",

        tested_genes =
            sum(
                !is.na(
                    res_lrt_df$pvalue
                )
            ),

        padj_available =
            sum(
                !is.na(
                    res_lrt_df$padj
                )
            ),

        significant =
            sum(
                res_lrt_df$significant,
                na.rm = TRUE
            ),

        up =
            NA_integer_,

        down =
            NA_integer_
    )
)


write_csv(
    summary_table,
    file.path(
        result_dir,
        "FA26_DESeq2_contrast_summary.csv"
    )
)


# ============================================================
# 13. Save model information
# ============================================================

model_info <- list(

    input_object =
        input_file,

    n_genes =
        nrow(dds_wald),

    n_samples =
        ncol(dds_wald),

    design =
        as.character(
            design(dds_wald)
        ),

    Treatment_levels =
        levels(
            dds_wald$Treatment
        ),

    Time_levels =
        levels(
            dds_wald$Time
        ),

    resultsNames =
        resultsNames(
            dds_wald
        ),

    alpha =
        0.05,

    created =
        Sys.time()
)


saveRDS(
    model_info,
    file.path(
        object_dir,
        "FA26_step06_DESeq2_model_info.rds"
    )
)


capture.output(
    sessionInfo(),
    file = file.path(
        result_dir,
        "FA26_step06_sessionInfo.txt"
    )
)


# ============================================================
# 14. Final report
# ============================================================

cat("\n============================================================\n")
cat("FA26 STEP 06 COMPLETE\n")
cat("============================================================\n\n")

cat("Design:\n")
print(design(dds_wald))

cat("\nReference treatment:\n")
cat(levels(dds_wald$Treatment)[1], "\n")

cat("\nReference time:\n")
cat(levels(dds_wald$Time)[1], "h\n")

cat("\nDEG summary — FDR < 0.05:\n\n")

print(
    summary_table,
    n = nrow(summary_table)
)

cat("\nOutputs:\n")
cat("  objects/FA26_dds_factorial_wald.rds\n")
cat("  objects/FA26_dds_factorial_interaction_LRT.rds\n")
cat("  objects/FA26_step06_DESeq2_model_info.rds\n")
cat("  results/deseq2/*.csv\n\n")

cat(
    "No log2FC threshold has been imposed.\n",
    "This is the initial statistical discovery layer.\n",
    sep = ""
)

