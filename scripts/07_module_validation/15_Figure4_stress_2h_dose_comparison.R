
# ============================================================
# FA26 Figure 4
# Canonical stress response at 2 h:
# 0.01 ug/g PFBA versus 1 ug/g PFBA
# ============================================================

options(stringsAsFactors = FALSE)

library(ggplot2)

cat("\n============================================================\n")
cat("FIGURE 4 — 2 h CANONICAL STRESS COMPARISON\n")
cat("0.01 ug/g versus 1 ug/g PFBA\n")
cat("============================================================\n\n")

outdir <- "results/figure4/stress_2h"
figdir <- "figures/working"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. READ EXISTING STRESS SCREEN
# ------------------------------------------------------------

stress <- read.csv(
    "results/figure4/FA26_stress_candidate_genes.csv",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("Rows in stress table:", nrow(stress), "\n")
cat("Unique stress-associated genes:",
    length(unique(stress$gene_id)), "\n\n")

# ------------------------------------------------------------
# 2. COLLAPSE MULTIPLE STRESS-CATEGORY ASSIGNMENTS PER GENE
#
# A gene may appear more than once because it matches multiple
# stress modules. For the dose comparison we want ONE row/gene.
# ------------------------------------------------------------

collapse_text <- function(x) {
    x <- unique(x[!is.na(x) & x != "" & x != "-"])
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = "; ")
}

first_nonmissing <- function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0) return(NA)
    x2[1]
}

gene_ids <- unique(stress$gene_id)

gene <- do.call(
    rbind,
    lapply(
        gene_ids,
        function(g) {

            x <- stress[stress$gene_id == g, , drop = FALSE]

            data.frame(
                gene_id = g,

                display_name =
                    first_nonmissing(x$display_name),

                product =
                    first_nonmissing(x$product),

                stress_category =
                    collapse_text(x$screen_category),

                High_LFC_2h =
                    first_nonmissing(x$PFBA1_shrunk_LFC_2h),

                High_padj_2h =
                    first_nonmissing(x$PFBA1_padj_2h),

                High_sig_2h =
                    any(x$PFBA1_sig_2h %in% TRUE),

                High_direction_2h =
                    first_nonmissing(x$PFBA1_direction_2h),

                Low_LFC_2h =
                    first_nonmissing(x$PFBA001_shrunk_LFC_2h),

                Low_padj_2h =
                    first_nonmissing(x$PFBA001_padj_2h),

                Low_sig_2h =
                    any(x$PFBA001_sig_2h %in% TRUE),

                Low_direction_2h =
                    first_nonmissing(x$PFBA001_direction_2h),

                stringsAsFactors = FALSE
            )
        }
    )
)

# Ensure numeric
gene$High_LFC_2h  <- as.numeric(gene$High_LFC_2h)
gene$High_padj_2h <- as.numeric(gene$High_padj_2h)

gene$Low_LFC_2h   <- as.numeric(gene$Low_LFC_2h)
gene$Low_padj_2h  <- as.numeric(gene$Low_padj_2h)

# ------------------------------------------------------------
# 3. CLASSIFY 2-h DOSE RESPONSE
# ------------------------------------------------------------

gene$response_2h <- "Neither significant"

gene$response_2h[
    gene$High_sig_2h &
    !gene$Low_sig_2h
] <- "1 ug/g only"

gene$response_2h[
    !gene$High_sig_2h &
    gene$Low_sig_2h
] <- "0.01 ug/g only"

gene$response_2h[
    gene$High_sig_2h &
    gene$Low_sig_2h
] <- "Both doses"

# Same/opposite direction among shared genes
gene$shared_direction <- NA_character_

both <- gene$High_sig_2h & gene$Low_sig_2h

gene$shared_direction[both] <-
    ifelse(
        sign(gene$High_LFC_2h[both]) ==
        sign(gene$Low_LFC_2h[both]),
        "Same direction",
        "Opposite direction"
    )

# ------------------------------------------------------------
# 4. BASIC COUNTS
# ------------------------------------------------------------

n_high <- sum(gene$High_sig_2h)
n_low  <- sum(gene$Low_sig_2h)
n_both <- sum(gene$High_sig_2h & gene$Low_sig_2h)

n_high_only <- sum(
    gene$High_sig_2h &
    !gene$Low_sig_2h
)

n_low_only <- sum(
    !gene$High_sig_2h &
    gene$Low_sig_2h
)

cat("2-h CANONICAL STRESS RESPONSE\n")
cat("------------------------------------------------------------\n")
cat("Significant at 1 ug/g      :", n_high, "\n")
cat("Significant at 0.01 ug/g   :", n_low, "\n")
cat("Significant at both doses  :", n_both, "\n")
cat("1 ug/g only                :", n_high_only, "\n")
cat("0.01 ug/g only             :", n_low_only, "\n\n")

cat("RESPONSE CLASSIFICATION\n")
cat("------------------------------------------------------------\n")
print(table(gene$response_2h))

# ------------------------------------------------------------
# 5. LIST LOW-DOSE 2-h STRESS GENES
# ------------------------------------------------------------

low2 <- gene[
    gene$Low_sig_2h,
    ,
    drop = FALSE
]

low2 <- low2[
    order(low2$Low_padj_2h),
    ,
    drop = FALSE
]

cat("\n============================================================\n")
cat("STRESS-ASSOCIATED GENES SIGNIFICANT AT 0.01 ug/g — 2 h\n")
cat("============================================================\n\n")

if (nrow(low2) == 0) {

    cat("NONE\n")

} else {

    print(
        low2[
            ,
            c(
                "gene_id",
                "display_name",
                "stress_category",
                "Low_LFC_2h",
                "Low_padj_2h",
                "High_LFC_2h",
                "High_padj_2h",
                "High_sig_2h",
                "response_2h"
            )
        ],
        row.names = FALSE
    )
}

# ------------------------------------------------------------
# 6. LIST HIGH-DOSE 2-h STRESS GENES
# ------------------------------------------------------------

high2 <- gene[
    gene$High_sig_2h,
    ,
    drop = FALSE
]

high2 <- high2[
    order(high2$High_padj_2h),
    ,
    drop = FALSE
]

cat("\n============================================================\n")
cat("STRESS-ASSOCIATED GENES SIGNIFICANT AT 1 ug/g — 2 h\n")
cat("============================================================\n\n")

print(
    high2[
        ,
        c(
            "gene_id",
            "display_name",
            "stress_category",
            "High_LFC_2h",
            "High_padj_2h",
            "Low_LFC_2h",
            "Low_padj_2h",
            "Low_sig_2h",
            "response_2h"
        )
    ],
    row.names = FALSE
)

# ------------------------------------------------------------
# 7. SHARED 2-h STRESS GENES
# ------------------------------------------------------------

shared2 <- gene[
    gene$High_sig_2h &
    gene$Low_sig_2h,
    ,
    drop = FALSE
]

cat("\n============================================================\n")
cat("SHARED 2-h STRESS GENES\n")
cat("============================================================\n\n")

if (nrow(shared2) == 0) {

    cat("NONE\n")

} else {

    print(
        shared2[
            ,
            c(
                "gene_id",
                "display_name",
                "stress_category",
                "Low_LFC_2h",
                "Low_padj_2h",
                "High_LFC_2h",
                "High_padj_2h",
                "shared_direction"
            )
        ],
        row.names = FALSE
    )
}

# ------------------------------------------------------------
# 8. CATEGORY-LEVEL SUMMARY
#
# Use unique gene-category combinations so genes are not
# double-counted within the same category.
# ------------------------------------------------------------

category_pairs <- unique(
    stress[
        ,
        c(
            "gene_id",
            "screen_category"
        )
    ]
)

category_summary <- do.call(
    rbind,
    lapply(
        unique(category_pairs$screen_category),
        function(cat) {

            ids <- category_pairs$gene_id[
                category_pairs$screen_category == cat
            ]

            x <- gene[gene$gene_id %in% ids, ]

            data.frame(
                stress_category = cat,
                candidate_genes = nrow(x),

                High_sig_2h =
                    sum(x$High_sig_2h),

                Low_sig_2h =
                    sum(x$Low_sig_2h),

                Shared_sig_2h =
                    sum(
                        x$High_sig_2h &
                        x$Low_sig_2h
                    ),

                High_only_2h =
                    sum(
                        x$High_sig_2h &
                        !x$Low_sig_2h
                    ),

                Low_only_2h =
                    sum(
                        !x$High_sig_2h &
                        x$Low_sig_2h
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

category_summary <- category_summary[
    order(
        -category_summary$High_sig_2h,
        -category_summary$Low_sig_2h
    ),
]

cat("\n============================================================\n")
cat("STRESS MODULE SUMMARY — 2 h\n")
cat("============================================================\n\n")

print(category_summary, row.names = FALSE)

# ------------------------------------------------------------
# 9. SAVE TABLES
# ------------------------------------------------------------

write.csv(
    gene,
    file.path(
        outdir,
        "FA26_stress_genelevel_high_vs_low_2h.csv"
    ),
    row.names = FALSE
)

write.csv(
    low2,
    file.path(
        outdir,
        "FA26_lowPFBA_stress_genes_2h.csv"
    ),
    row.names = FALSE
)

write.csv(
    high2,
    file.path(
        outdir,
        "FA26_highPFBA_stress_genes_2h.csv"
    ),
    row.names = FALSE
)

write.csv(
    shared2,
    file.path(
        outdir,
        "FA26_shared_stress_genes_2h.csv"
    ),
    row.names = FALSE
)

write.csv(
    category_summary,
    file.path(
        outdir,
        "FA26_stress_module_summary_2h.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 10. EXPLORATORY SCATTER
#
# Same stress-associated genes:
# x = low-dose response
# y = stimulatory-dose response
# ------------------------------------------------------------

plotdat <- gene[
    !is.na(gene$Low_LFC_2h) &
    !is.na(gene$High_LFC_2h),
    ,
    drop = FALSE
]

# Label significant genes only.
plotdat$label <- ""

plotdat$label[
    plotdat$High_sig_2h |
    plotdat$Low_sig_2h
] <- plotdat$display_name[
    plotdat$High_sig_2h |
    plotdat$Low_sig_2h
]

p <- ggplot(
    plotdat,
    aes(
        x = Low_LFC_2h,
        y = High_LFC_2h
    )
) +
    geom_hline(
        yintercept = 0,
        linewidth = 0.35,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linewidth = 0.35,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linewidth = 0.45,
        linetype = "dotted"
    ) +
    geom_point(
        aes(shape = response_2h),
        size = 3.0,
        alpha = 0.85
    ) +
    labs(
        x = expression(
            "0.01 " * mu * "g g"^{-1} *
            " PFBA, 2 h (shrunk log"[2] * "FC)"
        ),
        y = expression(
            "1 " * mu * "g g"^{-1} *
            " PFBA, 2 h (shrunk log"[2] * "FC)"
        ),
        shape = "2-h significance",
        title = "Canonical stress-associated transcription at 2 h"
    ) +
    theme_classic(base_size = 13) +
    theme(
        plot.title = element_text(
            size = 14,
            face = "bold"
        ),
        axis.title = element_text(
            size = 13,
            face = "bold"
        ),
        axis.text = element_text(
            size = 11
        ),
        legend.title = element_text(
            size = 11,
            face = "bold"
        ),
        legend.text = element_text(
            size = 10
        )
    )

# ggrepel only if installed
if (requireNamespace("ggrepel", quietly = TRUE)) {

    p <- p +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 3.4,
            max.overlaps = 30,
            box.padding = 0.4,
            point.padding = 0.25,
            min.segment.length = 0
        )
}

ggsave(
    file.path(
        figdir,
        "Figure4_stress_2h_dose_comparison_exploratory.pdf"
    ),
    p,
    width = 7.2,
    height = 6.3,
    units = "in"
)

ggsave(
    file.path(
        figdir,
        "Figure4_stress_2h_dose_comparison_exploratory.tiff"
    ),
    p,
    width = 7.2,
    height = 6.3,
    units = "in",
    dpi = 600,
    compression = "lzw"
)

cat("\n============================================================\n")
cat("2-h STRESS COMPARISON COMPLETE\n")
cat("============================================================\n\n")

cat("Tables saved to:\n")
cat(outdir, "\n\n")

cat("Exploratory figure:\n")
cat(
    "figures/working/Figure4_stress_2h_dose_comparison_exploratory.pdf\n"
)
cat(
    "figures/working/Figure4_stress_2h_dose_comparison_exploratory.tiff\n\n"
)

