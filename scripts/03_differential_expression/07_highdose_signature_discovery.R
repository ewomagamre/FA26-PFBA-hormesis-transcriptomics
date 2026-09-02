
# ============================================================
# FA26 STEP 07
# Early 1 ug/g PFBA transcriptional signature discovery
#
# Objectives
#   1. Define robust 1 ug/g vs Control DEGs at 2 h and 4 h
#   2. Classify:
#        Early-only
#        Persistent
#        Later-emerging
#   3. Preserve Up/Down direction
#   4. Calculate replicate-support metrics
#   5. Add shrunken LFC estimates
#   6. Integrate eggNOG/NCBI annotation
#   7. Integrate atlas-zone membership
#   8. Prepare signature sets for GO enrichment
#
# Statistical significance:
#   BH FDR < 0.05
#
# No gene is excluded simply because a condition contains zero
# counts. Sparse-expression metrics are used as interpretation
# flags rather than as a second significance test.
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(apeglm)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
})

cat("\n============================================================\n")
cat("FA26 STEP 07 — HIGH-DOSE SIGNATURE DISCOVERY\n")
cat("============================================================\n\n")

# ============================================================
# 1. Paths
# ============================================================

dds_file <-
    "objects/FA26_dds_condition_clean.rds"

annotation_file <-
    "objects/FA26_gene_atlas_master_annotated_zones_clockwise.rds"

outdir <-
    "results/deseq2/signatures"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 2. Load definitive clean condition model
# ============================================================

dds <- readRDS(dds_file)

cat("Dimensions:", dim(dds), "\n")
cat("Design:", deparse(design(dds)), "\n\n")

cat("Condition levels:\n")
print(levels(dds$Condition))

cat("\nCoefficient names:\n")
print(resultsNames(dds))

# ============================================================
# 3. Direct unshrunken DESeq2 results
#
# These are used for significance.
# ============================================================

res_2h <- results(
    dds,
    contrast = c(
        "Condition",
        "PFBA_1_2h",
        "Control_2h"
    ),
    alpha = 0.05
)

res_4h <- results(
    dds,
    contrast = c(
        "Condition",
        "PFBA_1_4h",
        "Control_4h"
    ),
    alpha = 0.05
)

res_2h_df <- as.data.frame(res_2h) %>%
    rownames_to_column("gene_id")

res_4h_df <- as.data.frame(res_4h) %>%
    rownames_to_column("gene_id")

# ============================================================
# 4. Shrunken LFC
#
# apeglm requires a coefficient rather than a general contrast.
#
# Because the condition model contains six levels, use
# relevelled copies so each requested comparison becomes a
# simple named coefficient.
#
# Significance remains based on the unshrunken Wald result.
# ============================================================

make_shrunken_comparison <- function(
    dds_input,
    numerator,
    denominator
) {

    d <- dds_input

    d$Condition <- relevel(
        factor(
            as.character(d$Condition)
        ),
        ref = denominator
    )

    design(d) <- ~ Condition

    # Refit coefficients using the already normalized counts.
    d <- nbinomWaldTest(
        d,
        betaPrior = FALSE,
        quiet = TRUE
    )

    rn <- resultsNames(d)

    expected <- paste0(
        "Condition_",
        make.names(numerator),
        "_vs_",
        make.names(denominator)
    )

    # DESeq2 may preserve punctuation differently, so find by
    # numerator and denominator if exact string differs.
    hit <- rn[
        grepl(
            gsub("\\.", "\\\\.", numerator),
            rn
        ) &
        grepl(
            gsub("\\.", "\\\\.", denominator),
            rn
        )
    ]

    if (length(hit) != 1) {

        # fallback: coefficient should be numerator vs reference
        hit <- rn[
            grepl(
                gsub("\\.", "\\\\.", numerator),
                rn
            )
        ]
    }

    if (length(hit) != 1) {
        cat("\nAvailable coefficient names:\n")
        print(rn)
        stop(
            paste(
                "Could not uniquely resolve coefficient for",
                numerator,
                "vs",
                denominator
            )
        )
    }

    cat(
        "Shrinking:",
        numerator,
        "vs",
        denominator,
        "using",
        hit,
        "\n"
    )

    shr <- lfcShrink(
        d,
        coef = hit,
        type = "apeglm"
    )

    as.data.frame(shr) %>%
        rownames_to_column("gene_id") %>%
        dplyr::select(
            gene_id,
            shrunk_LFC = log2FoldChange,
            shrunk_SE = lfcSE
        )
}

shr_2h <- make_shrunken_comparison(
    dds,
    "PFBA_1_2h",
    "Control_2h"
)

shr_4h <- make_shrunken_comparison(
    dds,
    "PFBA_1_4h",
    "Control_4h"
)

# ============================================================
# 5. Normalized-count replicate support
# ============================================================

norm_counts <- counts(
    dds,
    normalized = TRUE
)

meta <- as.data.frame(
    colData(dds)
)

stopifnot(
    identical(
        colnames(norm_counts),
        rownames(meta)
    )
)

group_support <- lapply(
    levels(meta$Condition),
    function(grp) {

        samples <- rownames(meta)[
            meta$Condition == grp
        ]

        mat <- norm_counts[
            ,
            samples,
            drop = FALSE
        ]

        tibble(
            gene_id = rownames(mat),

            Condition = grp,

            mean_norm_count =
                rowMeans(mat),

            median_norm_count =
                apply(
                    mat,
                    1,
                    median
                ),

            min_norm_count =
                apply(
                    mat,
                    1,
                    min
                ),

            max_norm_count =
                apply(
                    mat,
                    1,
                    max
                ),

            n_nonzero =
                rowSums(
                    mat > 0
                ),

            n_ge10 =
                rowSums(
                    mat >= 10
                )
        )
    }
) %>%
    bind_rows()

# ============================================================
# 6. Build support metrics for high-dose comparisons
# ============================================================

support_wide <- group_support %>%

    filter(
        Condition %in% c(
            "Control_2h",
            "PFBA_1_2h",
            "Control_4h",
            "PFBA_1_4h"
        )
    ) %>%

    pivot_wider(
        names_from = Condition,
        values_from = c(
            mean_norm_count,
            median_norm_count,
            min_norm_count,
            max_norm_count,
            n_nonzero,
            n_ge10
        ),
        names_sep = "__"
    )

# ============================================================
# 7. Merge 2 h + 4 h statistics
# ============================================================

master <- res_2h_df %>%

    transmute(
        gene_id,

        baseMean_2h =
            baseMean,

        LFC_2h =
            log2FoldChange,

        SE_2h =
            lfcSE,

        pvalue_2h =
            pvalue,

        padj_2h =
            padj,

        sig_2h =
            !is.na(padj) &
            padj < 0.05,

        direction_2h =
            case_when(
                sig_2h &
                LFC_2h > 0 ~ "Up",
                sig_2h &
                LFC_2h < 0 ~ "Down",
                TRUE ~ "NS"
            )
    ) %>%

    full_join(
        res_4h_df %>%

            transmute(
                gene_id,

                baseMean_4h =
                    baseMean,

                LFC_4h =
                    log2FoldChange,

                SE_4h =
                    lfcSE,

                pvalue_4h =
                    pvalue,

                padj_4h =
                    padj,

                sig_4h =
                    !is.na(padj) &
                    padj < 0.05,

                direction_4h =
                    case_when(
                        sig_4h &
                        LFC_4h > 0 ~ "Up",
                        sig_4h &
                        LFC_4h < 0 ~ "Down",
                        TRUE ~ "NS"
                    )
            ),

        by = "gene_id"
    ) %>%

    left_join(
        shr_2h %>%
            rename(
                shrunk_LFC_2h =
                    shrunk_LFC,
                shrunk_SE_2h =
                    shrunk_SE
            ),
        by = "gene_id"
    ) %>%

    left_join(
        shr_4h %>%
            rename(
                shrunk_LFC_4h =
                    shrunk_LFC,
                shrunk_SE_4h =
                    shrunk_SE
            ),
        by = "gene_id"
    ) %>%

    left_join(
        support_wide,
        by = "gene_id"
    )

# ============================================================
# 8. Signature classification
# ============================================================

master <- master %>%

    mutate(

        signature_class =
            case_when(

                sig_2h &
                sig_4h ~
                    "Persistent",

                sig_2h &
                !sig_4h ~
                    "Early-only",

                !sig_2h &
                sig_4h ~
                    "Later-emerging",

                TRUE ~
                    "Not significant"
            ),

        persistence_direction =
            case_when(

                signature_class ==
                    "Persistent" &
                direction_2h ==
                    direction_4h ~
                    paste0(
                        "Persistent ",
                        direction_2h
                    ),

                signature_class ==
                    "Persistent" &
                direction_2h !=
                    direction_4h ~
                    "Direction-switching",

                TRUE ~
                    NA_character_
            )
    )

# ============================================================
# 9. Replicate-support classification
#
# This is NOT a significance filter.
#
# High confidence:
#   all 4 replicates have nonzero normalized counts in either
#   comparison group at the relevant significant time.
#
# Moderate:
#   >=3/4 in both groups
#
# Sparse:
#   either group has <=2/4 nonzero replicates.
#
# We retain every DESeq2-significant gene regardless.
# ============================================================

support_time <- function(
    sig,
    control_n,
    treated_n
) {

    case_when(
        !sig ~ NA_character_,

        control_n == 4 &
        treated_n == 4 ~
            "High",

        control_n >= 3 &
        treated_n >= 3 ~
            "Moderate",

        TRUE ~
            "Sparse"
    )
}

master$support_2h <- support_time(
    master$sig_2h,
    master$n_nonzero__Control_2h,
    master$n_nonzero__PFBA_1_2h
)

master$support_4h <- support_time(
    master$sig_4h,
    master$n_nonzero__Control_4h,
    master$n_nonzero__PFBA_1_4h
)

master <- master %>%
    mutate(

        overall_support =
            case_when(

                signature_class ==
                    "Persistent" &
                support_2h == "High" &
                support_4h == "High" ~
                    "High",

                signature_class ==
                    "Persistent" &
                support_2h != "Sparse" &
                support_4h != "Sparse" ~
                    "Moderate",

                signature_class ==
                    "Early-only" ~
                    support_2h,

                signature_class ==
                    "Later-emerging" ~
                    support_4h,

                signature_class ==
                    "Persistent" ~
                    "Sparse",

                TRUE ~
                    NA_character_
            )
    )

# ============================================================
# 10. Add atlas zone + functional annotation where available
#
# Note:
# only the 3,000 atlas genes have atlas-zone assignment.
# Non-atlas DEGs remain fully retained.
# ============================================================

if (file.exists(annotation_file)) {

    atlas <- readRDS(
        annotation_file
    )

    annotation_cols <- intersect(
        c(
            "gene_id",
            "zone",
            "UMAP1",
            "UMAP2",
            "transcript_id",
            "protein_id",
            "NCBI_product",
            "eggNOG_description",
            "preferred_name",
            "COG_category",
            "GOs",
            "KEGG_ko",
            "KEGG_Pathway",
            "PFAMs"
        ),
        colnames(atlas)
    )

    annotation <- atlas %>%
        dplyr::select(
            all_of(
                annotation_cols
            )
        ) %>%
        distinct(
            gene_id,
            .keep_all = TRUE
        )

    master <- master %>%
        left_join(
            annotation,
            by = "gene_id"
        )
}

# ============================================================
# 11. Rank candidate genes
#
# Statistical status comes first.
# Then replicate support.
# Then shrunken effect size.
#
# Raw extreme MLE LFC is NOT used as the primary ranking
# variable.
# ============================================================

support_rank <- c(
    "High" = 1,
    "Moderate" = 2,
    "Sparse" = 3
)

master <- master %>%

    mutate(

        support_rank =
            unname(
                support_rank[
                    overall_support
                ]
            ),

        candidate_effect =
            case_when(

                signature_class ==
                    "Early-only" ~
                    abs(
                        shrunk_LFC_2h
                    ),

                signature_class ==
                    "Later-emerging" ~
                    abs(
                        shrunk_LFC_4h
                    ),

                signature_class ==
                    "Persistent" ~
                    pmin(
                        abs(
                            shrunk_LFC_2h
                        ),
                        abs(
                            shrunk_LFC_4h
                        )
                    ),

                TRUE ~
                    NA_real_
            ),

        candidate_padj =
            case_when(

                signature_class ==
                    "Early-only" ~
                    padj_2h,

                signature_class ==
                    "Later-emerging" ~
                    padj_4h,

                signature_class ==
                    "Persistent" ~
                    pmax(
                        padj_2h,
                        padj_4h,
                        na.rm = TRUE
                    ),

                TRUE ~
                    NA_real_
            )
    )

candidate_table <- master %>%

    filter(
        signature_class !=
            "Not significant"
    ) %>%

    arrange(
        signature_class,
        support_rank,
        candidate_padj,
        desc(
            candidate_effect
        )
    )

# ============================================================
# 12. Save complete outputs
# ============================================================

write_csv(
    master,
    file.path(
        outdir,
        "FA26_PFBA1_2h_4h_all_genes_signature_master.csv"
    )
)

write_csv(
    candidate_table,
    file.path(
        outdir,
        "FA26_PFBA1_candidate_signatures.csv"
    )
)

for (
    cls in
    c(
        "Early-only",
        "Persistent",
        "Later-emerging"
    )
) {

    z <- candidate_table %>%
        filter(
            signature_class ==
                cls
        )

    write_csv(
        z,
        file.path(
            outdir,
            paste0(
                "FA26_PFBA1_",
                gsub(
                    "-",
                    "_",
                    tolower(cls)
                ),
                "_genes.csv"
            )
        )
    )
}

# ============================================================
# 13. Signature summary
# ============================================================

signature_summary <- candidate_table %>%

    count(
        signature_class,
        direction_2h,
        direction_4h,
        name = "n_genes"
    ) %>%

    arrange(
        signature_class,
        desc(n_genes)
    )

support_summary <- candidate_table %>%

    count(
        signature_class,
        overall_support,
        name = "n_genes"
    )

zone_summary <- candidate_table %>%

    mutate(
        atlas_zone =
            ifelse(
                is.na(zone),
                "Outside top-3000 atlas",
                as.character(zone)
            )
    ) %>%

    count(
        signature_class,
        atlas_zone,
        name = "n_genes"
    )

write_csv(
    signature_summary,
    file.path(
        outdir,
        "FA26_PFBA1_signature_summary.csv"
    )
)

write_csv(
    support_summary,
    file.path(
        outdir,
        "FA26_PFBA1_signature_support_summary.csv"
    )
)

write_csv(
    zone_summary,
    file.path(
        outdir,
        "FA26_PFBA1_signature_zone_summary.csv"
    )
)

# ============================================================
# 14. Top robust annotated candidates
#
# Reporting table only.
# Does not alter the underlying signature sets.
# ============================================================

top_candidates <- candidate_table %>%

    filter(
        overall_support %in%
            c(
                "High",
                "Moderate"
            )
    ) %>%

    group_by(
        signature_class
    ) %>%

    arrange(
        candidate_padj,
        desc(
            candidate_effect
        ),
        .by_group = TRUE
    ) %>%

    slice_head(
        n = 25
    ) %>%

    ungroup()

write_csv(
    top_candidates,
    file.path(
        outdir,
        "FA26_PFBA1_top25_candidates_per_signature.csv"
    )
)

# ============================================================
# 15. Print terminal report
# ============================================================

cat("\n============================================================\n")
cat("STEP 07 COMPLETE\n")
cat("============================================================\n\n")

cat("SIGNATURE CLASSES\n\n")
print(signature_summary)

cat("\nREPLICATE SUPPORT\n\n")
print(support_summary)

cat("\nATLAS-ZONE DISTRIBUTION\n\n")
print(zone_summary)

cat("\nTOP ROBUST CANDIDATES\n\n")

display_cols <- intersect(
    c(
        "signature_class",
        "gene_id",
        "preferred_name",
        "NCBI_product",
        "eggNOG_description",
        "zone",
        "direction_2h",
        "shrunk_LFC_2h",
        "padj_2h",
        "direction_4h",
        "shrunk_LFC_4h",
        "padj_4h",
        "overall_support"
    ),
    colnames(
        top_candidates
    )
)

print(
    top_candidates %>%
        dplyr::select(
            all_of(display_cols)
        )
)

cat("\nOutputs:\n")
cat("  results/deseq2/signatures/\n\n")

cat(
    "Important: replicate support is an interpretation/ranking metric,\n",
    "not an additional differential-expression significance test.\n",
    sep = ""
)

