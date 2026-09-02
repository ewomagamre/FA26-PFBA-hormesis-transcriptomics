
suppressPackageStartupMessages({
    library(DESeq2)
    library(apeglm)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
})

cat("\n============================================================\n")
cat("FA26 STEP 08 — LOW-DOSE SIGNATURE DISCOVERY\n")
cat("============================================================\n\n")

dds_file <- "objects/FA26_dds_condition_clean.rds"
annotation_file <- "objects/FA26_gene_atlas_master_annotated_zones_clockwise.rds"
outdir <- "results/deseq2/signatures_lowdose"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(dds_file)

# ------------------------------------------------------------
# Direct DE results
# ------------------------------------------------------------

res_2h <- results(
    dds,
    contrast = c(
        "Condition",
        "PFBA_0.01_2h",
        "Control_2h"
    ),
    alpha = 0.05
)

res_4h <- results(
    dds,
    contrast = c(
        "Condition",
        "PFBA_0.01_4h",
        "Control_4h"
    ),
    alpha = 0.05
)

res_2h_df <- as.data.frame(res_2h) %>%
    rownames_to_column("gene_id")

res_4h_df <- as.data.frame(res_4h) %>%
    rownames_to_column("gene_id")

# ------------------------------------------------------------
# apeglm shrinkage using relevelled copies
# ------------------------------------------------------------

make_shrunken_comparison <- function(
    dds_input,
    numerator,
    denominator
) {

    d <- dds_input

    d$Condition <- relevel(
        factor(as.character(d$Condition)),
        ref = denominator
    )

    design(d) <- ~ Condition

    d <- nbinomWaldTest(
        d,
        betaPrior = FALSE,
        quiet = TRUE
    )

    rn <- resultsNames(d)

    hit <- rn[
        grepl(
            gsub("\\.", "\\\\.", numerator),
            rn
        )
    ]

    if (length(hit) != 1) {
        cat("\nAvailable coefficients:\n")
        print(rn)
        stop(
            paste(
                "Could not resolve coefficient for",
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
    "PFBA_0.01_2h",
    "Control_2h"
)

shr_4h <- make_shrunken_comparison(
    dds,
    "PFBA_0.01_4h",
    "Control_4h"
)

# ------------------------------------------------------------
# Replicate support from normalized counts
# ------------------------------------------------------------

norm_counts <- counts(
    dds,
    normalized = TRUE
)

meta <- as.data.frame(colData(dds))

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
                apply(mat, 1, median),

            n_nonzero =
                rowSums(mat > 0),

            n_ge10 =
                rowSums(mat >= 10)
        )
    }
) %>%
    bind_rows()

support_wide <- group_support %>%
    filter(
        Condition %in% c(
            "Control_2h",
            "PFBA_0.01_2h",
            "Control_4h",
            "PFBA_0.01_4h"
        )
    ) %>%
    pivot_wider(
        names_from = Condition,
        values_from = c(
            mean_norm_count,
            median_norm_count,
            n_nonzero,
            n_ge10
        ),
        names_sep = "__"
    )

# ------------------------------------------------------------
# Master table
# ------------------------------------------------------------

master <- res_2h_df %>%
    transmute(
        gene_id,

        baseMean_2h = baseMean,
        LFC_2h = log2FoldChange,
        SE_2h = lfcSE,
        pvalue_2h = pvalue,
        padj_2h = padj,

        sig_2h =
            !is.na(padj) &
            padj < 0.05,

        direction_2h =
            case_when(
                sig_2h & LFC_2h > 0 ~ "Up",
                sig_2h & LFC_2h < 0 ~ "Down",
                TRUE ~ "NS"
            )
    ) %>%

    full_join(
        res_4h_df %>%
            transmute(
                gene_id,

                baseMean_4h = baseMean,
                LFC_4h = log2FoldChange,
                SE_4h = lfcSE,
                pvalue_4h = pvalue,
                padj_4h = padj,

                sig_4h =
                    !is.na(padj) &
                    padj < 0.05,

                direction_4h =
                    case_when(
                        sig_4h & LFC_4h > 0 ~ "Up",
                        sig_4h & LFC_4h < 0 ~ "Down",
                        TRUE ~ "NS"
                    )
            ),
        by = "gene_id"
    ) %>%

    left_join(
        shr_2h %>%
            rename(
                shrunk_LFC_2h = shrunk_LFC,
                shrunk_SE_2h = shrunk_SE
            ),
        by = "gene_id"
    ) %>%

    left_join(
        shr_4h %>%
            rename(
                shrunk_LFC_4h = shrunk_LFC,
                shrunk_SE_4h = shrunk_SE
            ),
        by = "gene_id"
    ) %>%

    left_join(
        support_wide,
        by = "gene_id"
    )

# ------------------------------------------------------------
# Signature classes
# ------------------------------------------------------------

master <- master %>%
    mutate(
        signature_class =
            case_when(
                sig_2h & sig_4h ~ "Persistent",
                sig_2h & !sig_4h ~ "Early-only",
                !sig_2h & sig_4h ~ "Later-emerging",
                TRUE ~ "Not significant"
            ),

        persistence_direction =
            case_when(
                signature_class == "Persistent" &
                    direction_2h == direction_4h ~
                    paste0(
                        "Persistent ",
                        direction_2h
                    ),

                signature_class == "Persistent" &
                    direction_2h != direction_4h ~
                    "Direction-switching",

                TRUE ~ NA_character_
            )
    )

# ------------------------------------------------------------
# Replicate support
# ------------------------------------------------------------

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
    master$n_nonzero__PFBA_0.01_2h
)

master$support_4h <- support_time(
    master$sig_4h,
    master$n_nonzero__Control_4h,
    master$n_nonzero__PFBA_0.01_4h
)

master <- master %>%
    mutate(
        overall_support =
            case_when(

                signature_class == "Persistent" &
                    support_2h == "High" &
                    support_4h == "High" ~
                    "High",

                signature_class == "Persistent" &
                    support_2h != "Sparse" &
                    support_4h != "Sparse" ~
                    "Moderate",

                signature_class == "Early-only" ~
                    support_2h,

                signature_class == "Later-emerging" ~
                    support_4h,

                signature_class == "Persistent" ~
                    "Sparse",

                TRUE ~
                    NA_character_
            )
    )

# ------------------------------------------------------------
# Atlas annotation where available
# ------------------------------------------------------------

if (file.exists(annotation_file)) {

    atlas <- readRDS(annotation_file)

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
            all_of(annotation_cols)
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

# ------------------------------------------------------------
# Candidate ranking
# ------------------------------------------------------------

support_order <- c(
    "High" = 1,
    "Moderate" = 2,
    "Sparse" = 3
)

master <- master %>%
    mutate(
        support_rank =
            unname(
                support_order[
                    overall_support
                ]
            ),

        candidate_effect =
            case_when(

                signature_class == "Early-only" ~
                    abs(shrunk_LFC_2h),

                signature_class == "Later-emerging" ~
                    abs(shrunk_LFC_4h),

                signature_class == "Persistent" ~
                    pmin(
                        abs(shrunk_LFC_2h),
                        abs(shrunk_LFC_4h)
                    ),

                TRUE ~ NA_real_
            ),

        candidate_padj =
            case_when(

                signature_class == "Early-only" ~
                    padj_2h,

                signature_class == "Later-emerging" ~
                    padj_4h,

                signature_class == "Persistent" ~
                    pmax(
                        padj_2h,
                        padj_4h,
                        na.rm = TRUE
                    ),

                TRUE ~ NA_real_
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
        desc(candidate_effect)
    )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write_csv(
    master,
    file.path(
        outdir,
        "FA26_PFBA001_all_genes_signature_master.csv"
    )
)

write_csv(
    candidate_table,
    file.path(
        outdir,
        "FA26_PFBA001_candidate_signatures.csv"
    )
)

signature_summary <- candidate_table %>%
    count(
        signature_class,
        direction_2h,
        direction_4h,
        name = "n_genes"
    )

support_summary <- candidate_table %>%
    count(
        signature_class,
        overall_support,
        name = "n_genes"
    )

zone_summary <- candidate_table %>%
    mutate(
        atlas_status =
            ifelse(
                is.na(zone),
                "Not included in 3,000-gene atlas",
                as.character(zone)
            )
    ) %>%
    count(
        signature_class,
        atlas_status,
        name = "n_genes"
    )

write_csv(
    signature_summary,
    file.path(
        outdir,
        "FA26_PFBA001_signature_summary.csv"
    )
)

write_csv(
    support_summary,
    file.path(
        outdir,
        "FA26_PFBA001_support_summary.csv"
    )
)

write_csv(
    zone_summary,
    file.path(
        outdir,
        "FA26_PFBA001_atlas_distribution.csv"
    )
)

cat("\n============================================================\n")
cat("STEP 08 COMPLETE\n")
cat("============================================================\n\n")

cat("SIGNATURE CLASSES\n")
print(signature_summary)

cat("\nREPLICATE SUPPORT\n")
print(support_summary)

cat("\nATLAS DISTRIBUTION\n")
print(zone_summary)

