
# ============================================================
# FA26 Step 06B
# Validation of DESeq2 contrasts, interaction tests,
# effect sizes, and LFC shrinkage
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(readr)
    library(tibble)
})

dir.create("results/deseq2/validation",
           recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("FA26 STEP 06B — DESeq2 VALIDATION\n")
cat("============================================================\n")


# ============================================================
# 1. Load factorial Wald model
# ============================================================

dds <- readRDS("objects/FA26_dds_factorial_wald.rds")

cat("\nFactorial design:\n")
print(design(dds))

cat("\nCoefficient names:\n")
print(resultsNames(dds))


# ============================================================
# 2. Construct independent ~Condition model
#
# This provides a direct check of the four treatment-vs-control
# contrasts without relying on factorial coefficient arithmetic.
# ============================================================

dds_condition <- dds

dds_condition$Condition <- factor(
    as.character(dds_condition$Condition),
    levels = c(
        "Control_2h",
        "Control_4h",
        "PFBA_0.01_2h",
        "PFBA_0.01_4h",
        "PFBA_1_2h",
        "PFBA_1_4h"
    )
)

design(dds_condition) <- ~ Condition

dds_condition <- DESeq(
    dds_condition,
    test = "Wald",
    fitType = "local",
    quiet = FALSE
)

saveRDS(
    dds_condition,
    "objects/FA26_dds_condition_validation.rds"
)


# ============================================================
# 3. Extract direct Condition contrasts
# ============================================================

get_condition_result <- function(numerator, denominator) {

    as.data.frame(
        results(
            dds_condition,
            contrast = c(
                "Condition",
                numerator,
                denominator
            ),
            alpha = 0.05
        )
    ) %>%
        rownames_to_column("gene_id")
}


condition_results <- list(

    PFBA_0.01_vs_Control_2h =
        get_condition_result(
            "PFBA_0.01_2h",
            "Control_2h"
        ),

    PFBA_1_vs_Control_2h =
        get_condition_result(
            "PFBA_1_2h",
            "Control_2h"
        ),

    PFBA_0.01_vs_Control_4h =
        get_condition_result(
            "PFBA_0.01_4h",
            "Control_4h"
        ),

    PFBA_1_vs_Control_4h =
        get_condition_result(
            "PFBA_1_4h",
            "Control_4h"
        )
)


# ============================================================
# 4. Compare factorial vs Condition-model results
# ============================================================

validation_summary <- list()

for (nm in names(condition_results)) {

    factorial_file <- file.path(
        "results/deseq2",
        paste0(nm, ".csv")
    )

    fac <- read_csv(
        factorial_file,
        show_col_types = FALSE
    ) %>%
        select(
            gene_id,
            factorial_LFC = log2FoldChange,
            factorial_padj = padj
        )

    con <- condition_results[[nm]] %>%
        select(
            gene_id,
            condition_LFC = log2FoldChange,
            condition_padj = padj
        )

    z <- inner_join(
        fac,
        con,
        by = "gene_id"
    )

    lfc_cor <- cor(
        z$factorial_LFC,
        z$condition_LFC,
        use = "complete.obs",
        method = "pearson"
    )

    max_lfc_difference <- max(
        abs(
            z$factorial_LFC -
            z$condition_LFC
        ),
        na.rm = TRUE
    )

    fac_sig <-
        !is.na(z$factorial_padj) &
        z$factorial_padj < 0.05

    con_sig <-
        !is.na(z$condition_padj) &
        z$condition_padj < 0.05

    validation_summary[[nm]] <- tibble(

        comparison = nm,

        LFC_correlation = lfc_cor,

        max_abs_LFC_difference =
            max_lfc_difference,

        factorial_significant =
            sum(fac_sig),

        condition_significant =
            sum(con_sig),

        significant_in_both =
            sum(fac_sig & con_sig),

        factorial_only =
            sum(fac_sig & !con_sig),

        condition_only =
            sum(!fac_sig & con_sig)
    )

    write_csv(
        z,
        file.path(
            "results/deseq2/validation",
            paste0(
                nm,
                "_factorial_vs_condition.csv"
            )
        )
    )
}

validation_summary <- bind_rows(
    validation_summary
)

write_csv(
    validation_summary,
    "results/deseq2/validation/FA26_factorial_condition_validation_summary.csv"
)


# ============================================================
# 5. Effect-size summary for all Step 06 contrasts
#
# We are NOT imposing an LFC cutoff.
# This asks how many FDR-significant genes also exceed
# descriptive effect-size thresholds.
# ============================================================

files <- list.files(
    "results/deseq2",
    pattern = "\\.csv$",
    full.names = TRUE
)

files <- files[
    !grepl(
        "summary|Omnibus",
        basename(files)
    )
]

effect_summary <- list()

for (f in files) {

    z <- read_csv(
        f,
        show_col_types = FALSE
    )

    if (!all(
        c(
            "log2FoldChange",
            "padj"
        ) %in% colnames(z)
    )) next

    sig <-
        !is.na(z$padj) &
        z$padj < 0.05

    effect_summary[[basename(f)]] <- tibble(

        comparison =
            sub("\\.csv$", "", basename(f)),

        FDR_005 =
            sum(sig),

        FDR_005_absLFC_0.5 =
            sum(
                sig &
                abs(z$log2FoldChange) >= 0.5,
                na.rm = TRUE
            ),

        FDR_005_absLFC_1 =
            sum(
                sig &
                abs(z$log2FoldChange) >= 1,
                na.rm = TRUE
            ),

        median_absLFC_significant =
            if (sum(sig) > 0)
                median(
                    abs(z$log2FoldChange[sig]),
                    na.rm = TRUE
                )
            else NA_real_,

        max_absLFC_significant =
            if (sum(sig) > 0)
                max(
                    abs(z$log2FoldChange[sig]),
                    na.rm = TRUE
                )
            else NA_real_
    )
}

effect_summary <- bind_rows(
    effect_summary
)

write_csv(
    effect_summary,
    "results/deseq2/validation/FA26_effect_size_summary.csv"
)


# ============================================================
# 6. Interaction audit
# ============================================================

low_int <- as.data.frame(
    results(
        dds,
        name = "TreatmentPFBA_0.01.Time4",
        alpha = 0.05
    )
) %>%
    rownames_to_column("gene_id")

high_int <- as.data.frame(
    results(
        dds,
        name = "TreatmentPFBA_1.Time4",
        alpha = 0.05
    )
) %>%
    rownames_to_column("gene_id")


# ============================================================
# 7. Refit omnibus interaction LRT explicitly with fitType=local
#
# Full:    Treatment + Time + Treatment:Time
# Reduced: Treatment + Time
#
# This is a 2-df joint test of BOTH interaction coefficients.
# ============================================================

dds_lrt <- dds

design(dds_lrt) <- ~ Treatment * Time

dds_lrt <- DESeq(
    dds_lrt,
    test = "LRT",
    reduced = ~ Treatment + Time,
    fitType = "local",
    quiet = FALSE
)

lrt <- as.data.frame(
    results(
        dds_lrt,
        alpha = 0.05
    )
) %>%
    rownames_to_column("gene_id")


interaction_audit <- low_int %>%

    select(
        gene_id,
        low_LFC = log2FoldChange,
        low_p = pvalue,
        low_padj = padj
    ) %>%

    full_join(
        high_int %>%
            select(
                gene_id,
                high_LFC = log2FoldChange,
                high_p = pvalue,
                high_padj = padj
            ),
        by = "gene_id"
    ) %>%

    full_join(
        lrt %>%
            select(
                gene_id,
                LRT_stat = stat,
                LRT_p = pvalue,
                LRT_padj = padj
            ),
        by = "gene_id"
    ) %>%

    mutate(

        low_sig =
            !is.na(low_padj) &
            low_padj < 0.05,

        high_sig =
            !is.na(high_padj) &
            high_padj < 0.05,

        either_Wald_sig =
            low_sig | high_sig,

        LRT_sig =
            !is.na(LRT_padj) &
            LRT_padj < 0.05
    )


write_csv(
    interaction_audit,
    "results/deseq2/validation/FA26_interaction_Wald_LRT_audit.csv"
)


interaction_summary <- tibble(

    low_Wald_FDR05 =
        sum(
            interaction_audit$low_sig,
            na.rm = TRUE
        ),

    high_Wald_FDR05 =
        sum(
            interaction_audit$high_sig,
            na.rm = TRUE
        ),

    either_Wald_FDR05 =
        sum(
            interaction_audit$either_Wald_sig,
            na.rm = TRUE
        ),

    omnibus_LRT_FDR05 =
        sum(
            interaction_audit$LRT_sig,
            na.rm = TRUE
        ),

    Wald_and_LRT =
        sum(
            interaction_audit$either_Wald_sig &
            interaction_audit$LRT_sig,
            na.rm = TRUE
        ),

    Wald_only =
        sum(
            interaction_audit$either_Wald_sig &
            !interaction_audit$LRT_sig,
            na.rm = TRUE
        ),

    LRT_only =
        sum(
            !interaction_audit$either_Wald_sig &
            interaction_audit$LRT_sig,
            na.rm = TRUE
        )
)


write_csv(
    interaction_summary,
    "results/deseq2/validation/FA26_interaction_summary.csv"
)


# ============================================================
# 8. LFC shrinkage availability
#
# apeglm is preferred for individual coefficients.
# ashr is an alternative if installed.
#
# We only generate shrunk estimates for simple named
# coefficients here. Composite 4 h contrasts will be handled
# separately once the validation is complete.
# ============================================================

cat("\nLFC shrinkage packages:\n")

has_apeglm <- requireNamespace(
    "apeglm",
    quietly = TRUE
)

has_ashr <- requireNamespace(
    "ashr",
    quietly = TRUE
)

cat("apeglm:", has_apeglm, "\n")
cat("ashr:  ", has_ashr, "\n")


if (has_apeglm) {

    shrink_names <- c(
        "Treatment_PFBA_0.01_vs_Control",
        "Treatment_PFBA_1_vs_Control",
        "Time_4_vs_2",
        "TreatmentPFBA_0.01.Time4",
        "TreatmentPFBA_1.Time4"
    )

    for (coef_name in shrink_names) {

        shr <- lfcShrink(
            dds,
            coef = coef_name,
            type = "apeglm"
        )

        shr_df <- as.data.frame(
            shr
        ) %>%
            rownames_to_column(
                "gene_id"
            )

        write_csv(
            shr_df,
            file.path(
                "results/deseq2/validation",
                paste0(
                    "shrunk_",
                    coef_name,
                    ".csv"
                )
            )
        )
    }
}


# ============================================================
# 9. Final report
# ============================================================

cat("\n============================================================\n")
cat("STEP 06B VALIDATION COMPLETE\n")
cat("============================================================\n")

cat("\nFACTORIAL vs CONDITION MODEL:\n\n")
print(
    validation_summary,
    n = nrow(validation_summary)
)

cat("\nEFFECT-SIZE SUMMARY:\n\n")
print(
    effect_summary,
    n = nrow(effect_summary)
)

cat("\nINTERACTION AUDIT:\n\n")
print(
    interaction_summary
)

cat("\nInterpretation rules:\n")
cat("  FDR threshold = 0.05\n")
cat("  |LFC| 0.5 and 1.0 are descriptive checks only\n")
cat("  No effect-size threshold has yet been imposed\n")
cat("  LRT is a joint 2-df interaction test\n")
cat("  Wald tests evaluate individual interaction coefficients\n\n")

