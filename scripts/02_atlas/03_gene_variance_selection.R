# ============================================================
# FA26 24-sample primary analysis
# Step 03 — Gene variance selection and atlas optimization
#
# Goal:
#   Determine the appropriate number of variable genes for the
#   gene-level transcriptomic atlas without pre-selecting 3,000.
#
# Candidate feature sets:
#   2,000 | 3,000 | 4,000 | 6,000 genes
#
# Evidence used:
#   1. Ranked gene variance
#   2. Cumulative variance captured
#   3. Data-derived cumulative-variance elbow
#   4. Gene-level UMAP topology
#   5. Procrustes topology agreement
#
# Selection rule:
#   Choose the smallest tested feature set at or beyond the
#   variance elbow that also shows strong topology agreement
#   with the broad 6,000-gene atlas.
#
# 3,000 genes is NOT assumed beforehand.
# ============================================================


# ------------------------------------------------------------
# 0. Required packages
# ------------------------------------------------------------

required_packages <- c(
    "DESeq2",
    "ggplot2",
    "dplyr",
    "readr",
    "uwot",
    "vegan",
    "patchwork"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0) {

    stop(
        paste0(
            "\nMissing required packages:\n",
            paste(missing_packages, collapse = ", "),
            "\n"
        )
    )
}

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(uwot)
    library(vegan)
    library(patchwork)
})


cat("============================================================\n")
cat("FA26 Step 03: Gene variance selection and atlas optimization\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Analysis parameters
# ------------------------------------------------------------

candidate_sizes <- c(
    2000,
    3000,
    4000,
    6000
)

reference_size <- max(candidate_sizes)

# Keep identical UMAP parameters for all candidate sets
umap_neighbors <- 30
umap_min_dist <- 0.30
umap_metric <- "euclidean"
umap_seed <- 260725

# Topology validation
n_permutations <- 999

# Criterion for strong topology agreement
procrustes_stability_threshold <- 0.90


# ------------------------------------------------------------
# 2. File locations
# ------------------------------------------------------------

input_vsd <- "objects/FA26_vsd_primary24.rds"

variance_dir <- "results/gene_variance"
umap_dir <- file.path(
    variance_dir,
    "candidate_umaps"
)

figure_dir <- file.path(
    variance_dir,
    "figures"
)

dir.create(
    variance_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    umap_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "objects",
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Load 24-sample VST object
# ------------------------------------------------------------

if (!file.exists(input_vsd)) {

    stop(
        paste(
            "Input file not found:",
            input_vsd
        )
    )
}

vsd24 <- readRDS(
    input_vsd
)

vst_mat <- assay(
    vsd24
)

cat(
    "Genes in VST matrix:",
    nrow(vst_mat),
    "\n"
)

cat(
    "Samples:",
    ncol(vst_mat),
    "\n\n"
)

if (nrow(vst_mat) < reference_size) {

    stop(
        paste0(
            "Only ",
            nrow(vst_mat),
            " genes are available, fewer than the 6,000-gene candidate."
        )
    )
}


# ------------------------------------------------------------
# 4. Calculate gene-wise variance
# ------------------------------------------------------------

gene_variance <- apply(
    vst_mat,
    1,
    stats::var
)

variance_df <- data.frame(
    gene_id = rownames(vst_mat),
    variance = as.numeric(gene_variance),
    stringsAsFactors = FALSE
) %>%
    arrange(
        desc(variance)
    ) %>%
    mutate(
        rank = seq_len(n()),
        cumulative_variance =
            cumsum(variance) /
            sum(variance) *
            100,
        fraction_of_genes =
            rank /
            n() *
            100
    )

write_csv(
    variance_df,
    file.path(
        variance_dir,
        "FA26_gene_variance_all.csv"
    )
)

cat("Gene-wise variance calculated.\n")
cat(
    "Total ranked genes:",
    nrow(variance_df),
    "\n\n"
)


# ------------------------------------------------------------
# 5. Variance summary
# ------------------------------------------------------------

cat("Variance summary:\n")

print(
    summary(
        variance_df$variance
    )
)

cat("\nVariance quantiles:\n")

variance_quantiles <- quantile(
    variance_df$variance,
    probs = c(
        0.50,
        0.60,
        0.70,
        0.75,
        0.80,
        0.85,
        0.90,
        0.95,
        0.975,
        0.99
    )
)

print(
    variance_quantiles
)


# ------------------------------------------------------------
# 6. Determine cumulative-variance elbow
# ------------------------------------------------------------
#
# Maximum perpendicular distance from the straight line
# connecting the beginning and end of the normalized cumulative
# variance curve.
#
# This estimates the point of diminishing return in added genes.
# ------------------------------------------------------------

x <- variance_df$rank
y <- variance_df$cumulative_variance

x_norm <- (
    x - min(x)
) / (
    max(x) - min(x)
)

y_norm <- (
    y - min(y)
) / (
    max(y) - min(y)
)

x1 <- x_norm[1]
y1 <- y_norm[1]

x2 <- x_norm[length(x_norm)]
y2 <- y_norm[length(y_norm)]

distance_to_line <- abs(
    (y2 - y1) * x_norm -
    (x2 - x1) * y_norm +
    x2 * y1 -
    y2 * x1
) /
sqrt(
    (y2 - y1)^2 +
    (x2 - x1)^2
)

elbow_index <- which.max(
    distance_to_line
)

elbow_rank <- variance_df$rank[
    elbow_index
]

elbow_variance <- variance_df$variance[
    elbow_index
]

elbow_cumvar <- variance_df$cumulative_variance[
    elbow_index
]

cat("\n============================================================\n")
cat("Cumulative variance elbow\n")
cat("============================================================\n")

cat(
    "Elbow gene rank:",
    elbow_rank,
    "\n"
)

cat(
    "Variance at elbow:",
    round(
        elbow_variance,
        4
    ),
    "\n"
)

cat(
    "Cumulative variance at elbow:",
    round(
        elbow_cumvar,
        2
    ),
    "%\n"
)


# ------------------------------------------------------------
# 7. Create candidate feature sets
# ------------------------------------------------------------

candidate_ids <- list()
candidate_matrices <- list()

candidate_summary <- data.frame(
    candidate_size = integer(),
    fraction_retained_genes = numeric(),
    minimum_variance = numeric(),
    median_variance = numeric(),
    maximum_variance = numeric(),
    cumulative_variance_pct = numeric(),
    stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Creating candidate feature sets\n")
cat("============================================================\n\n")


for (n_genes in candidate_sizes) {

    selected_df <- variance_df %>%
        slice_head(
            n = n_genes
        )

    selected_ids <- selected_df$gene_id

    selected_matrix <- vst_mat[
        selected_ids,
        ,
        drop = FALSE
    ]

    candidate_ids[[as.character(n_genes)]] <- selected_ids

    candidate_matrices[[as.character(n_genes)]] <- selected_matrix

    cumulative_pct <- selected_df$cumulative_variance[
        nrow(selected_df)
    ]

    write_csv(
        selected_df,
        file.path(
            variance_dir,
            paste0(
                "FA26_top",
                n_genes,
                "_variable_genes.csv"
            )
        )
    )

    saveRDS(
        selected_ids,
        file.path(
            "objects",
            paste0(
                "FA26_top",
                n_genes,
                "_variable_gene_ids.rds"
            )
        )
    )

    saveRDS(
        selected_matrix,
        file.path(
            "objects",
            paste0(
                "FA26_vst_top",
                n_genes,
                "_variable_genes.rds"
            )
        )
    )

    candidate_summary <- bind_rows(
        candidate_summary,
        data.frame(
            candidate_size =
                n_genes,

            fraction_retained_genes =
                n_genes /
                nrow(vst_mat) *
                100,

            minimum_variance =
                min(
                    selected_df$variance
                ),

            median_variance =
                median(
                    selected_df$variance
                ),

            maximum_variance =
                max(
                    selected_df$variance
                ),

            cumulative_variance_pct =
                cumulative_pct,

            stringsAsFactors = FALSE
        )
    )

    cat(
        "Top",
        n_genes,
        "genes | min variance =",
        round(
            min(selected_df$variance),
            4
        ),
        "| cumulative variance =",
        round(
            cumulative_pct,
            2
        ),
        "%\n"
    )
}


# ------------------------------------------------------------
# 8. Marginal information gain
# ------------------------------------------------------------

candidate_summary <- candidate_summary %>%
    arrange(
        candidate_size
    ) %>%
    mutate(
        added_genes =
            candidate_size -
            lag(candidate_size),

        marginal_variance_gain =
            cumulative_variance_pct -
            lag(cumulative_variance_pct),

        variance_gain_per_1000_genes =
            marginal_variance_gain /
            added_genes *
            1000
    )

write_csv(
    candidate_summary,
    file.path(
        variance_dir,
        "FA26_candidate_gene_sets_initial_summary.csv"
    )
)


# ------------------------------------------------------------
# 9. Row-wise gene expression pattern scaling
# ------------------------------------------------------------
#
# Genes are observations.
# Samples are features.
#
# Row z-scoring emphasizes the pattern of each gene across the
# 24 samples rather than absolute expression magnitude.
# ------------------------------------------------------------

row_zscore <- function(mat) {

    gene_means <- rowMeans(
        mat
    )

    gene_sds <- apply(
        mat,
        1,
        stats::sd
    )

    bad_sd <- (
        !is.finite(gene_sds) |
        gene_sds == 0
    )

    gene_sds[
        bad_sd
    ] <- 1

    scaled_mat <- sweep(
        mat,
        1,
        gene_means,
        FUN = "-"
    )

    scaled_mat <- sweep(
        scaled_mat,
        1,
        gene_sds,
        FUN = "/"
    )

    return(
        scaled_mat
    )
}


# ------------------------------------------------------------
# 10. Construct candidate gene-level UMAPs
# ------------------------------------------------------------

candidate_umaps <- list()

cat("\n============================================================\n")
cat("Constructing candidate UMAPs\n")
cat("============================================================\n\n")


for (n_genes in candidate_sizes) {

    cat(
        "UMAP:",
        n_genes,
        "genes\n"
    )

    candidate_matrix <- candidate_matrices[[as.character(n_genes)]]

    candidate_scaled <- row_zscore(
        candidate_matrix
    )

    set.seed(
        umap_seed
    )

    embedding <- uwot::umap(
        candidate_scaled,
        n_neighbors = umap_neighbors,
        min_dist = umap_min_dist,
        metric = umap_metric,
        n_components = 2,
        init = "spectral",
        n_threads = 1,
        verbose = FALSE
    )

    rownames(
        embedding
    ) <- rownames(
        candidate_scaled
    )

    colnames(
        embedding
    ) <- c(
        "UMAP1",
        "UMAP2"
    )

    embedding_df <- data.frame(
        gene_id =
            rownames(embedding),

        UMAP1 =
            embedding[, 1],

        UMAP2 =
            embedding[, 2],

        stringsAsFactors = FALSE
    )

    candidate_umaps[[as.character(n_genes)]] <- embedding_df

    write_csv(
        embedding_df,
        file.path(
            umap_dir,
            paste0(
                "FA26_gene_UMAP_top",
                n_genes,
                ".csv"
            )
        )
    )

    saveRDS(
        embedding_df,
        file.path(
            "objects",
            paste0(
                "FA26_gene_UMAP_top",
                n_genes,
                ".rds"
            )
        )
    )
}


# ------------------------------------------------------------
# 11. Procrustes comparison function
# ------------------------------------------------------------

compare_umaps <- function(
    smaller_n,
    larger_n,
    umap_list,
    permutations = 999
) {

    small_df <- umap_list[[as.character(smaller_n)]]

    large_df <- umap_list[[as.character(larger_n)]]

    shared_ids <- intersect(
        small_df$gene_id,
        large_df$gene_id
    )

    small_shared <- small_df[
        match(
            shared_ids,
            small_df$gene_id
        ),
        c(
            "UMAP1",
            "UMAP2"
        )
    ]

    large_shared <- large_df[
        match(
            shared_ids,
            large_df$gene_id
        ),
        c(
            "UMAP1",
            "UMAP2"
        )
    ]

    small_shared <- as.matrix(
        small_shared
    )

    large_shared <- as.matrix(
        large_shared
    )

    proc_fit <- vegan::procrustes(
        small_shared,
        large_shared,
        symmetric = TRUE
    )

    set.seed(
        umap_seed
    )

    proc_test <- vegan::protest(
        small_shared,
        large_shared,
        permutations = permutations
    )

    return(
        data.frame(
            comparison =
                paste0(
                    smaller_n,
                    "_vs_",
                    larger_n
                ),

            smaller_set =
                smaller_n,

            larger_set =
                larger_n,

            shared_genes =
                length(shared_ids),

            procrustes_r =
                as.numeric(
                    proc_test$t0
                ),

            permutation_p =
                as.numeric(
                    proc_test$signif
                ),

            procrustes_ss =
                as.numeric(
                    proc_fit$ss
                ),

            stringsAsFactors = FALSE
        )
    )
}


# ------------------------------------------------------------
# 12. Compare each candidate with 6K benchmark
# ------------------------------------------------------------

reference_results <- data.frame(
    comparison = character(),
    smaller_set = integer(),
    larger_set = integer(),
    shared_genes = integer(),
    procrustes_r = numeric(),
    permutation_p = numeric(),
    procrustes_ss = numeric(),
    stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Procrustes comparisons against 6,000-gene benchmark\n")
cat("============================================================\n\n")


for (
    n_genes in candidate_sizes[
        candidate_sizes < reference_size
    ]
) {

    result <- compare_umaps(
        smaller_n = n_genes,
        larger_n = reference_size,
        umap_list = candidate_umaps,
        permutations = n_permutations
    )

    reference_results <- bind_rows(
        reference_results,
        result
    )

    cat(
        n_genes,
        "vs",
        reference_size,
        "| r =",
        round(
            result$procrustes_r,
            3
        ),
        "| P =",
        result$permutation_p,
        "\n"
    )
}


# Add benchmark itself
reference_results <- bind_rows(
    reference_results,
    data.frame(
        comparison =
            paste0(
                reference_size,
                "_vs_",
                reference_size
            ),

        smaller_set =
            reference_size,

        larger_set =
            reference_size,

        shared_genes =
            reference_size,

        procrustes_r =
            1,

        permutation_p =
            NA_real_,

        procrustes_ss =
            0,

        stringsAsFactors = FALSE
    )
)

write_csv(
    reference_results,
    file.path(
        variance_dir,
        "FA26_Procrustes_vs_6000.csv"
    )
)


# ------------------------------------------------------------
# 13. Adjacent candidate comparisons
# ------------------------------------------------------------

adjacent_pairs <- list(
    c(
        2000,
        3000
    ),
    c(
        3000,
        4000
    ),
    c(
        4000,
        6000
    )
)

adjacent_results <- data.frame(
    comparison = character(),
    smaller_set = integer(),
    larger_set = integer(),
    shared_genes = integer(),
    procrustes_r = numeric(),
    permutation_p = numeric(),
    procrustes_ss = numeric(),
    stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Adjacent candidate topology comparisons\n")
cat("============================================================\n\n")


for (pair in adjacent_pairs) {

    result <- compare_umaps(
        smaller_n = pair[1],
        larger_n = pair[2],
        umap_list = candidate_umaps,
        permutations = n_permutations
    )

    adjacent_results <- bind_rows(
        adjacent_results,
        result
    )

    cat(
        pair[1],
        "vs",
        pair[2],
        "| r =",
        round(
            result$procrustes_r,
            3
        ),
        "| P =",
        result$permutation_p,
        "\n"
    )
}

write_csv(
    adjacent_results,
    file.path(
        variance_dir,
        "FA26_Procrustes_adjacent_candidates.csv"
    )
)


# ------------------------------------------------------------
# 14. Add topology statistics to candidate table
# ------------------------------------------------------------

candidate_summary <- candidate_summary %>%
    left_join(
        reference_results %>%
            select(
                smaller_set,
                procrustes_r,
                permutation_p,
                procrustes_ss
            ),
        by = c(
            "candidate_size" =
                "smaller_set"
        )
    )


# ------------------------------------------------------------
# 15. Data-driven feature selection
# ------------------------------------------------------------
#
# Candidate must satisfy:
#
# 1. candidate size >= data-derived variance elbow
# 2. Procrustes r >= predefined topology stability threshold
#
# Choose smallest tested candidate satisfying both.
#
# If no candidate satisfies topology threshold, choose smallest
# tested candidate at or beyond the variance elbow and flag it
# for manual review rather than pretending the criterion passed.
# ------------------------------------------------------------

candidate_summary <- candidate_summary %>%
    mutate(
        beyond_variance_elbow =
            candidate_size >= elbow_rank,

        topology_stable =
            procrustes_r >=
            procrustes_stability_threshold,

        selection_eligible =
            beyond_variance_elbow &
            topology_stable
    )


eligible_candidates <- candidate_summary %>%
    filter(
        selection_eligible
    )


if (nrow(eligible_candidates) > 0) {

    selected_n <- min(
        eligible_candidates$candidate_size
    )

    selection_status <-
        "PASS"

    selection_reason <- paste0(
        "Selected the smallest tested candidate at or beyond ",
        "the cumulative-variance elbow with Procrustes r >= ",
        procrustes_stability_threshold,
        " relative to the 6,000-gene benchmark."
    )

} else {

    above_elbow <- candidate_summary %>%
        filter(
            beyond_variance_elbow
        )

    if (nrow(above_elbow) > 0) {

        selected_n <- min(
            above_elbow$candidate_size
        )

        selection_status <-
            "REVIEW"

        selection_reason <- paste0(
            "No tested candidate satisfied the predefined ",
            "Procrustes stability threshold. The smallest ",
            "candidate at or beyond the variance elbow is ",
            selected_n,
            " genes and requires manual topology review."
        )

    } else {

        selected_n <- reference_size

        selection_status <-
            "REVIEW"

        selection_reason <- paste0(
            "The cumulative-variance elbow lies beyond all ",
            "tested candidate sizes. The 6,000-gene set was ",
            "retained for manual review."
        )
    }
}


candidate_summary <- candidate_summary %>%
    mutate(
        selected =
            candidate_size == selected_n
    )

write_csv(
    candidate_summary,
    file.path(
        variance_dir,
        "FA26_feature_selection_optimization_summary.csv"
    )
)


# ------------------------------------------------------------
# 16. Save selected feature set
# ------------------------------------------------------------

selected_ids <- 

selected_matrix <- 

selected_gene_table <- variance_df %>%
    filter(
        gene_id %in% selected_ids
    ) %>%
    arrange(
        rank
    )

saveRDS(
    selected_ids,
    "objects/FA26_atlas_selected_gene_ids.rds"
)

saveRDS(
    selected_matrix,
    "objects/FA26_atlas_selected_vst_matrix.rds"
)

write_csv(
    selected_gene_table,
    file.path(
        variance_dir,
        "FA26_atlas_selected_variable_genes.csv"
    )
)


# ------------------------------------------------------------
# ============================================================
# 17. Manuscript-grade figure styling
# ============================================================

# Publication figures live in centralized figure directories.
# Statistical tables and objects remain under results/objects.

supp_figure_dir <- "figures/supplementary"

dir.create(
    supp_figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# Consistent manuscript palette
# Red = selected feature threshold
# Gray = statistical elbow/reference
# Other colors distinguish data curves without competing
# with the selected-threshold emphasis.

accent_red <- "#B22222"
elbow_gray <- "#777777"
variance_blue <- "#356A92"
log_green <- "#176B2C"
cumulative_purple <- "#55218A"
point_charcoal <- "#303030"

publication_theme <- theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 17,
            margin = margin(b = 8)
        ),
        plot.subtitle = element_text(
            size = 12,
            margin = margin(b = 12)
        ),
        axis.title = element_text(
            size = 14
        ),
        axis.text = element_text(
            size = 11,
            color = "black"
        ),
        axis.line = element_line(
            linewidth = 0.7,
            color = "black"
        ),
        axis.ticks = element_line(
            linewidth = 0.6,
            color = "black"
        ),
        plot.margin = margin(12, 14, 12, 12)
    )


# ============================================================
# 18. Values used in feature-selection figure
# ============================================================

selected_row <- variance_df %>%
    filter(rank == selected_n)

selected_variance <- selected_row$variance
selected_cumvar <- selected_row$cumulative_variance
selected_fraction <- selected_n / nrow(variance_df) * 100

elbow_label <- paste0(
    "Elbow = ",
    format(elbow_rank, big.mark = ","),
    " genes (",
    round(elbow_cumvar, 1),
    "% cumulative variance)"
)

selection_label <- paste0(
    "Selected = ",
    format(selected_n, big.mark = ","),
    " genes (",
    round(selected_cumvar, 1),
    "% cumulative variance)"
)


# ============================================================
# 19. Panel A — Ranked variance
# ============================================================

pA <- ggplot(
    variance_df,
    aes(x = rank, y = variance)
) +
    geom_line(
        linewidth = 0.75,
        color = variance_blue
    ) +
    geom_vline(
        xintercept = elbow_rank,
        linetype = "dotted",
        linewidth = 0.8,
        color = elbow_gray
    ) +
    geom_vline(
        xintercept = selected_n,
        linetype = "dashed",
        linewidth = 0.9,
        color = accent_red
    ) +
    geom_point(
        data = selected_row,
        size = 3.2,
        color = accent_red
    ) +
    annotate(
        "text",
        x = selected_n + 350,
        y = selected_variance + 1.6,
        label = paste0(
            "Selected threshold\n",
            "Rank = ",
            format(selected_n, big.mark = ","),
            "\nVariance = ",
            round(selected_variance, 3)
        ),
        hjust = 0,
        size = 4.2,
        color = accent_red
    ) +
    labs(
        title = "Ranked gene variance",
        subtitle = elbow_label,
        x = "Gene rank",
        y = "Variance"
    ) +
    publication_theme


# ============================================================
# 20. Panel B — Log-ranked variance
# ============================================================

pB <- ggplot(
    variance_df,
    aes(x = rank, y = log10(variance))
) +
    geom_line(
        linewidth = 0.75,
        color = log_green
    ) +
    geom_vline(
        xintercept = elbow_rank,
        linetype = "dotted",
        linewidth = 0.8,
        color = elbow_gray
    ) +
    geom_vline(
        xintercept = selected_n,
        linetype = "dashed",
        linewidth = 0.9,
        color = accent_red
    ) +
    geom_point(
        data = selected_row,
        size = 3.2,
        color = accent_red
    ) +
    annotate(
        "text",
        x = selected_n + 400,
        y = log10(selected_variance) + 0.17,
        label = "Selected feature threshold",
        hjust = 0,
        size = 4.2,
        color = accent_red
    ) +
    labs(
        title = "Ranked gene variance on log scale",
        subtitle = "Gray dotted line = variance elbow; red dashed line = selected threshold",
        x = "Gene rank",
        y = expression(log[10](Variance))
    ) +
    publication_theme


# ============================================================
# 21. Panel C — Cumulative variance
# ============================================================

pC <- ggplot(
    variance_df,
    aes(x = rank, y = cumulative_variance)
) +
    geom_line(
        linewidth = 0.9,
        color = cumulative_purple
    ) +
    geom_vline(
        xintercept = elbow_rank,
        linetype = "dotted",
        linewidth = 0.8,
        color = elbow_gray
    ) +
    geom_vline(
        xintercept = selected_n,
        linetype = "dashed",
        linewidth = 0.9,
        color = accent_red
    ) +
    geom_hline(
        yintercept = selected_cumvar,
        linetype = "dotted",
        linewidth = 0.75,
        color = accent_red
    ) +
    geom_point(
        data = selected_row,
        size = 3.4,
        color = accent_red
    ) +
    annotate(
        "text",
        x = elbow_rank - 140,
        y = elbow_cumvar - 10,
        label = paste0(
            "Variance elbow\n",
            format(elbow_rank, big.mark = ","),
            " genes\n",
            round(elbow_cumvar, 1),
            "% variance"
        ),
        hjust = 1,
        size = 4.0,
        color = elbow_gray
    ) +
    annotate(
        "text",
        x = selected_n + 430,
        y = selected_cumvar - 10,
        label = paste0(
            "Selected threshold\n",
            format(selected_n, big.mark = ","),
            " genes\n",
            round(selected_fraction, 1),
            "% of retained genes\n",
            round(selected_cumvar, 1),
            "% cumulative variance"
        ),
        hjust = 0,
        size = 4.1,
        color = accent_red
    ) +
    labs(
        title = "Cumulative variance captured",
        x = "Gene rank",
        y = "Cumulative variance (%)"
    ) +
    publication_theme


# ============================================================
# 22. Assemble Supplementary feature-selection figure
# ============================================================

optimization_figure <- (pA | pB) / pC +
    patchwork::plot_annotation(
        title = "Optimization of feature selection for construction of the FA26 gene-level transcriptomic atlas",
        subtitle = paste0(
            elbow_label,
            "; ",
            selection_label,
            "."
        ),
        tag_levels = "A",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 19,
                hjust = 0.5
            ),
            plot.subtitle = element_text(
                size = 12,
                hjust = 0.5,
                margin = margin(b = 10)
            )
        )
    )

ggsave(
    file.path(
        supp_figure_dir,
        "FigureS3_FeatureSelectionOptimization.pdf"
    ),
    optimization_figure,
    width = 14,
    height = 10,
    units = "in"
)

ggsave(
    file.path(
        supp_figure_dir,
        "FigureS3_FeatureSelectionOptimization.tiff"
    ),
    optimization_figure,
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 23. Candidate UMAP plotting function
# ============================================================

make_umap_plot <- function(
    n_genes,
    embedding_df,
    summary_df
) {

    stats_row <- summary_df %>%
        filter(candidate_size == n_genes)

    if (n_genes == reference_size) {

        topology_text <- "Maximal 6,000-gene benchmark"

    } else {

        topology_text <- paste0(
            "Procrustes r = ",
            sprintf("%.3f", stats_row$procrustes_r),
            "; permutation P = ",
            format.pval(
                stats_row$permutation_p,
                digits = 2,
                eps = 0.001
            )
        )
    }

    title_text <- paste0(
        format(n_genes, big.mark = ","),
        " genes"
    )

    border_color <- if (n_genes == selected_n) {
        accent_red
    } else {
        "black"
    }

    border_width <- if (n_genes == selected_n) {
        1.1
    } else {
        0.65
    }

    ggplot(
        embedding_df,
        aes(x = UMAP1, y = UMAP2)
    ) +
        geom_point(
            size = 0.62,
            alpha = 0.52,
            color = point_charcoal
        ) +
        labs(
            title = title_text,
            subtitle = topology_text,
            x = "UMAP 1",
            y = "UMAP 2"
        ) +
        theme_classic(base_size = 13) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 16,
                margin = margin(b = 3)
            ),
            plot.subtitle = element_text(
                size = 11,
                margin = margin(b = 8)
            ),
            axis.title = element_text(size = 13),
            axis.text = element_text(
                size = 10,
                color = "black"
            ),
            axis.line = element_line(
                linewidth = border_width,
                color = border_color
            ),
            axis.ticks = element_line(
                linewidth = 0.55,
                color = "black"
            ),
            plot.margin = margin(10, 12, 10, 10)
        )
}


# ============================================================
# 24. Four-panel topology robustness figure
# ============================================================

p2000 <- make_umap_plot(
    2000,
    candidate_umaps[["2000"]],
    candidate_summary
)

p3000 <- make_umap_plot(
    3000,
    candidate_umaps[["3000"]],
    candidate_summary
)

p4000 <- make_umap_plot(
    4000,
    candidate_umaps[["4000"]],
    candidate_summary
)

p6000 <- make_umap_plot(
    6000,
    candidate_umaps[["6000"]],
    candidate_summary
)

robustness_figure <- (
    p2000 | p3000
) / (
    p4000 | p6000
) +
    patchwork::plot_annotation(
        title = "Robustness of FA26 transcriptomic atlas topology across feature-selection thresholds",
        subtitle = paste0(
            "Identical UMAP parameters across candidate sets; ",
            "symmetric Procrustes analysis with ",
            n_permutations,
            " permutations. The selected threshold is indicated by the accent axis."
        ),
        tag_levels = "A",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 19,
                hjust = 0.5
            ),
            plot.subtitle = element_text(
                size = 12,
                hjust = 0.5,
                margin = margin(b = 10)
            )
        )
    )

ggsave(
    file.path(
        supp_figure_dir,
        "FigureS4_FeatureSelectionRobustness.pdf"
    ),
    robustness_figure,
    width = 14,
    height = 10,
    units = "in"
)

ggsave(
    file.path(
        supp_figure_dir,
        "FigureS4_FeatureSelectionRobustness.tiff"
    ),
    robustness_figure,
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 25. Write selection report
# ============================================================

report_file <- file.path(
    variance_dir,
    "FA26_feature_selection_decision.txt"
)

sink(report_file)

cat("FA26 gene-atlas feature-selection optimization\n")
cat("==============================================\n\n")

cat("Retained genes:", nrow(vst_mat), "\n")
cat("Samples:", ncol(vst_mat), "\n\n")

cat("Variance elbow rank:", elbow_rank, "\n")
cat("Variance at elbow:", round(elbow_variance, 4), "\n")
cat(
    "Cumulative variance at elbow:",
    round(elbow_cumvar, 2),
    "%\n\n"
)

cat("Candidate optimization table:\n\n")
print(candidate_summary)

cat("\nSelected candidate:", selected_n, "genes\n")
cat("Selection status:", selection_status, "\n")
cat(
    "Cumulative variance captured:",
    round(selected_cumvar, 2),
    "%\n"
)
cat(
    "Fraction of retained genes:",
    round(selected_fraction, 2),
    "%\n\n"
)
cat("Selection rationale:\n")
cat(selection_reason, "\n")

sink()


# ============================================================
# 26. Save session information
# ============================================================

capture.output(
    sessionInfo(),
    file = file.path(
        variance_dir,
        "FA26_step03_sessionInfo.txt"
    )
)


# ============================================================
# 27. Final terminal summary
# ============================================================

cat("\n============================================================\n")
cat("FA26 Step 03 optimization complete\n")
cat("============================================================\n\n")

cat("Variance elbow rank:", elbow_rank, "\n")
cat(
    "Variance elbow cumulative variance:",
    round(elbow_cumvar, 2),
    "%\n\n"
)

cat("Candidate optimization results:\n\n")

print(
    candidate_summary %>%
        select(
            candidate_size,
            cumulative_variance_pct,
            marginal_variance_gain,
            variance_gain_per_1000_genes,
            procrustes_r,
            permutation_p,
            beyond_variance_elbow,
            topology_stable,
            selected
        )
)

cat("\n------------------------------------------------------------\n")
cat("DATA-DRIVEN CANDIDATE:", selected_n, "genes\n")
cat("Selection status:", selection_status, "\n")
cat(
    "Cumulative variance captured:",
    round(selected_cumvar, 2),
    "%\n"
)

cat("\nRationale:\n")
cat(selection_reason, "\n")

cat("\nManuscript supplementary figures:\n")
cat("  figures/supplementary/FigureS3_FeatureSelectionOptimization.pdf\n")
cat("  figures/supplementary/FigureS3_FeatureSelectionOptimization.tiff\n")
cat("  figures/supplementary/FigureS4_FeatureSelectionRobustness.pdf\n")
cat("  figures/supplementary/FigureS4_FeatureSelectionRobustness.tiff\n")

cat("\nStep 03 complete pending final biological review of 3K versus 4K.\n")
