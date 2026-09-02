
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
    library(GO.db)
    library(AnnotationDbi)
    library(ggplot2)
    library(patchwork)
})

cat("\n============================================================\n")
cat("FA26 STEP 11 — DEG GO BIOLOGICAL PROCESS ENRICHMENT\n")
cat("============================================================\n\n")

annotation_file <- paste0(
    "results/annotation/full_DESeq2/",
    "FA26_DESeq2_11661_gene_annotation_master.csv"
)

go_file <- paste0(
    "results/annotation/full_DESeq2/",
    "FA26_DESeq2_GO_membership.csv"
)

high_file <- paste0(
    "results/deseq2/signatures/",
    "FA26_PFBA1_candidate_signatures.csv"
)

low_file <- paste0(
    "results/deseq2/signatures_lowdose/",
    "FA26_PFBA001_candidate_signatures.csv"
)

outdir <- "results/deseq2/GO_enrichment"
figdir <- "figures/main"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 1. Load full annotation universe
# ============================================================

annotation <- read_csv(
    annotation_file,
    show_col_types = FALSE
)

go_mem <- read_csv(
    go_file,
    show_col_types = FALSE
)

# Use genes with >=1 GO annotation as the GO-testable universe.
go_universe <- unique(
    go_mem$gene_id
)

cat(
    "DESeq2 genes:",
    nrow(annotation),
    "\n"
)

cat(
    "Genes with GO annotation:",
    length(go_universe),
    "\n"
)

# ============================================================
# 2. GO Biological Process ontology
# ============================================================

go_lookup <- AnnotationDbi::select(
    GO.db,
    keys = unique(
        go_mem$GO_ID
    ),
    keytype = "GOID",
    columns = c(
        "TERM",
        "ONTOLOGY"
    )
) %>%
    distinct(
        GOID,
        .keep_all = TRUE
    )

go_bp <- go_mem %>%
    left_join(
        go_lookup,
        by = c(
            "GO_ID" = "GOID"
        )
    ) %>%
    filter(
        ONTOLOGY == "BP"
    ) %>%
    dplyr::select(
        gene_id,
        GO_ID,
        TERM
    ) %>%
    distinct()

bp_universe <- unique(
    go_bp$gene_id
)

cat(
    "Genes with GO Biological Process annotation:",
    length(bp_universe),
    "\n"
)

# ============================================================
# 3. Load signatures
# ============================================================

high <- read_csv(
    high_file,
    show_col_types = FALSE
)

low <- read_csv(
    low_file,
    show_col_types = FALSE
)

# ============================================================
# 4. Define biologically interpretable gene sets
#
# Formal signature sets remain based on DESeq2 significance.
# Robust low-dose subset is a sensitivity analysis.
# ============================================================

gene_sets <- list(

    High_EarlyOnly_All =
        high %>%
        filter(
            signature_class == "Early-only"
        ) %>%
        pull(gene_id),

    High_EarlyOnly_Up =
        high %>%
        filter(
            signature_class == "Early-only",
            direction_2h == "Up"
        ) %>%
        pull(gene_id),

    High_EarlyOnly_Down =
        high %>%
        filter(
            signature_class == "Early-only",
            direction_2h == "Down"
        ) %>%
        pull(gene_id),

    High_Persistent_All =
        high %>%
        filter(
            signature_class == "Persistent"
        ) %>%
        pull(gene_id),

    High_Persistent_UpUp =
        high %>%
        filter(
            signature_class == "Persistent",
            direction_2h == "Up",
            direction_4h == "Up"
        ) %>%
        pull(gene_id),

    High_LaterEmerging_All =
        high %>%
        filter(
            signature_class == "Later-emerging"
        ) %>%
        pull(gene_id),

    High_LaterEmerging_Up =
        high %>%
        filter(
            signature_class == "Later-emerging",
            direction_4h == "Up"
        ) %>%
        pull(gene_id),

    High_LaterEmerging_Down =
        high %>%
        filter(
            signature_class == "Later-emerging",
            direction_4h == "Down"
        ) %>%
        pull(gene_id),

    Low_EarlyOnly_All =
        low %>%
        filter(
            signature_class == "Early-only"
        ) %>%
        pull(gene_id),

    Low_LaterEmerging_All =
        low %>%
        filter(
            signature_class == "Later-emerging"
        ) %>%
        pull(gene_id),

    Low_LaterEmerging_Up =
        low %>%
        filter(
            signature_class == "Later-emerging",
            direction_4h == "Up"
        ) %>%
        pull(gene_id),

    Low_LaterEmerging_Down =
        low %>%
        filter(
            signature_class == "Later-emerging",
            direction_4h == "Down"
        ) %>%
        pull(gene_id),

    Low_LaterEmerging_Robust =
        low %>%
        filter(
            signature_class == "Later-emerging",
            overall_support %in%
                c(
                    "High",
                    "Moderate"
                )
        ) %>%
        pull(gene_id),

    Low_LaterEmerging_RobustDown =
        low %>%
        filter(
            signature_class == "Later-emerging",
            direction_4h == "Down",
            overall_support %in%
                c(
                    "High",
                    "Moderate"
                )
        ) %>%
        pull(gene_id)
)

# ============================================================
# 5. ORA function
# ============================================================

run_ORA <- function(
    genes,
    set_name,
    min_term_genes = 10,
    min_set_hits = 4
) {

    genes <- unique(
        intersect(
            genes,
            bp_universe
        )
    )

    background <- bp_universe

    set_size <- length(
        genes
    )

    universe_size <- length(
        background
    )

    if (set_size < min_set_hits) {

        return(
            tibble()
        )
    }

    terms <- go_bp %>%
        group_by(
            GO_ID,
            TERM
        ) %>%
        summarise(
            term_genes =
                n_distinct(
                    gene_id
                ),
            .groups = "drop"
        ) %>%
        filter(
            term_genes >=
                min_term_genes
        )

    results_list <- lapply(
        seq_len(
            nrow(terms)
        ),
        function(i) {

            go_id <- terms$GO_ID[i]
            term_name <- terms$TERM[i]

            members <- go_bp %>%
                filter(
                    GO_ID == go_id
                ) %>%
                pull(
                    gene_id
                ) %>%
                unique()

            members <- intersect(
                members,
                background
            )

            a <- length(
                intersect(
                    genes,
                    members
                )
            )

            if (a < min_set_hits) {
                return(NULL)
            }

            b <- set_size - a
            c <- length(members) - a

            d <- universe_size -
                a -
                b -
                c

            mat <- matrix(
                c(
                    a,
                    b,
                    c,
                    d
                ),
                nrow = 2,
                byrow = TRUE
            )

            ft <- fisher.test(
                mat,
                alternative = "greater"
            )

            zone_fraction <-
                a /
                set_size

            background_fraction <-
                length(members) /
                universe_size

            tibble(
                gene_set =
                    set_name,

                GO_ID =
                    go_id,

                TERM =
                    term_name,

                set_hits =
                    a,

                set_size =
                    set_size,

                term_genes =
                    length(
                        members
                    ),

                universe_size =
                    universe_size,

                set_fraction =
                    zone_fraction,

                background_fraction =
                    background_fraction,

                fold_enrichment =
                    zone_fraction /
                    background_fraction,

                odds_ratio =
                    unname(
                        ft$estimate
                    ),

                p_value =
                    ft$p.value
            )
        }
    )

    out <- bind_rows(
        results_list
    )

    if (nrow(out) == 0) {
        return(out)
    }

    out %>%
        mutate(
            padj =
                p.adjust(
                    p_value,
                    method = "BH"
                )
        ) %>%
        arrange(
            padj,
            desc(
                fold_enrichment
            )
        )
}

# ============================================================
# 6. Run all sets
# ============================================================

all_enrichment <- bind_rows(

    lapply(
        names(gene_sets),
        function(nm) {

            cat(
                "Testing:",
                nm,
                "| input genes:",
                length(
                    unique(
                        gene_sets[[nm]]
                    )
                ),
                "\n"
            )

            run_ORA(
                gene_sets[[nm]],
                nm
            )
        }
    )
)

write_csv(
    all_enrichment,
    file.path(
        outdir,
        "FA26_DEG_GO_BP_enrichment_all.csv"
    )
)

significant <- all_enrichment %>%
    filter(
        padj < 0.05,
        fold_enrichment > 1
    )

write_csv(
    significant,
    file.path(
        outdir,
        "FA26_DEG_GO_BP_enrichment_significant.csv"
    )
)

# ============================================================
# 7. Sensitivity summary
# ============================================================

enrichment_summary <- significant %>%
    count(
        gene_set,
        name = "significant_BP_terms"
    )

set_summary <- tibble(
    gene_set =
        names(
            gene_sets
        ),

    total_signature_genes =
        sapply(
            gene_sets,
            function(g)
                length(
                    unique(g)
                )
        ),

    GO_BP_testable_genes =
        sapply(
            gene_sets,
            function(g)
                length(
                    intersect(
                        unique(g),
                        bp_universe
                    )
                )
        )
) %>%
    left_join(
        enrichment_summary,
        by = "gene_set"
    ) %>%
    mutate(
        significant_BP_terms =
            replace_na(
                significant_BP_terms,
                0L
            )
    )

write_csv(
    set_summary,
    file.path(
        outdir,
        "FA26_DEG_GO_BP_set_summary.csv"
    )
)

# ============================================================
# 8. Reduce redundancy for visualization
#
# Jaccard >= 0.70 -> retain higher-ranked term.
# ============================================================

reduce_terms <- function(
    df,
    max_terms = 8,
    cutoff = 0.70
) {

    if (nrow(df) == 0) {
        return(df)
    }

    df <- df %>%
        arrange(
            padj,
            desc(
                fold_enrichment
            )
        )

    retained <- character()

    for (
        go_i in df$GO_ID
    ) {

        genes_i <- go_bp %>%
            filter(
                GO_ID == go_i
            ) %>%
            pull(
                gene_id
            ) %>%
            unique()

        redundant <- FALSE

        if (length(retained) > 0) {

            for (
                go_j in retained
            ) {

                genes_j <- go_bp %>%
                    filter(
                        GO_ID == go_j
                    ) %>%
                    pull(
                        gene_id
                    ) %>%
                    unique()

                union_n <- length(
                    union(
                        genes_i,
                        genes_j
                    )
                )

                if (union_n == 0) {
                    next
                }

                jac <- length(
                    intersect(
                        genes_i,
                        genes_j
                    )
                ) /
                    union_n

                if (jac >= cutoff) {

                    redundant <- TRUE
                    break
                }
            }
        }

        if (!redundant) {

            retained <- c(
                retained,
                go_i
            )
        }

        if (
            length(retained) >=
                max_terms
        ) {
            break
        }
    }

    df %>%
        filter(
            GO_ID %in%
                retained
        )
}

# Main manuscript-level sets
main_sets <- c(
    "High_EarlyOnly_Up",
    "High_Persistent_UpUp",
    "High_LaterEmerging_Down",
    "Low_LaterEmerging_RobustDown"
)

plot_df <- bind_rows(

    lapply(
        main_sets,
        function(s) {

            reduce_terms(
                significant %>%
                    filter(
                        gene_set == s
                    ),
                max_terms = 7
            )
        }
    )
)

write_csv(
    plot_df,
    file.path(
        outdir,
        "FA26_Figure6_GO_terms_used.csv"
    )
)

# ============================================================
# 9. Figure 6
# ============================================================

plot_names <- c(

    "High_EarlyOnly_Up" =
        "1 µg/g: early activation",

    "High_Persistent_UpUp" =
        "1 µg/g: persistent activation",

    "High_LaterEmerging_Down" =
        "1 µg/g: later repression",

    "Low_LaterEmerging_RobustDown" =
        "0.01 µg/g: later repression\n(high/moderate support)"
)

plot_df <- plot_df %>%
    mutate(
        program =
            factor(
                gene_set,
                levels =
                    main_sets,
                labels =
                    unname(
                        plot_names[
                            main_sets
                        ]
                    )
            ),

        minus_log10_FDR =
            -log10(
                padj
            )
    )

if (nrow(plot_df) > 0) {

    fig6 <- ggplot(
        plot_df,
        aes(
            x =
                fold_enrichment,

            y =
                reorder(
                    TERM,
                    fold_enrichment
                ),

            size =
                set_hits,

            color =
                minus_log10_FDR
        )
    ) +

        geom_point(
            alpha = 0.95
        ) +

        facet_wrap(
            ~program,
            ncol = 2,
            scales = "free_y"
        ) +

        scale_size_continuous(
            name =
                "Genes",
            range =
                c(
                    4,
                    9
                )
        ) +

        scale_color_viridis_c(
            name =
                expression(
                    -log[10](
                        FDR
                    )
                ),
            option =
                "D",
            direction =
                1
        ) +

        labs(
            x =
                "Fold enrichment",
            y =
                NULL
        ) +

        theme_bw(
            base_size = 15
        ) +

        theme(
            strip.background =
                element_blank(),

            strip.text =
                element_text(
                    face = "bold",
                    size = 15
                ),

            axis.title.x =
                element_text(
                    face = "bold",
                    size = 17
                ),

            axis.text.x =
                element_text(
                    size = 12.5,
                    color = "black"
                ),

            axis.text.y =
                element_text(
                    size = 12,
                    color = "black"
                ),

            panel.grid.major.y =
                element_blank(),

            panel.grid.minor =
                element_blank(),

            legend.title =
                element_text(
                    face = "bold",
                    size = 12.5
                ),

            legend.text =
                element_text(
                    size = 11.5
                ),

            legend.position =
                "bottom"
        )

    ggsave(
        file.path(
            figdir,
            "Figure6_DEG_GO_BiologicalProcesses.pdf"
        ),
        fig6,
        width = 14,
        height = 10,
        units = "in",
        bg = "white"
    )

    ggsave(
        file.path(
            figdir,
            "Figure6_DEG_GO_BiologicalProcesses.tiff"
        ),
        fig6,
        width = 14,
        height = 10,
        units = "in",
        dpi = 600,
        compression = "lzw",
        bg = "white"
    )
}

cat("\n============================================================\n")
cat("STEP 11 COMPLETE\n")
cat("============================================================\n\n")

cat("GENE-SET SUMMARY\n\n")
print(set_summary)

cat("\nSIGNIFICANT GO BP TERMS BY SET\n\n")
print(enrichment_summary)

cat("\nFIGURE 6 TERMS\n\n")

print(
    plot_df %>%
        dplyr::select(
            program,
            GO_ID,
            TERM,
            set_hits,
            fold_enrichment,
            padj
        )
)

