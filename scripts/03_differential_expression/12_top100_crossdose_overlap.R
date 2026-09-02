
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
})

cat("\n============================================================\n")
cat("FA26 TOP-100 CROSS-DOSE / TEMPORAL OVERLAP\n")
cat("============================================================\n\n")

high <- read_csv(
    "results/deseq2/signatures/FA26_PFBA1_candidate_signatures.csv",
    show_col_types = FALSE
)

low <- read_csv(
    "results/deseq2/signatures_lowdose/FA26_PFBA001_candidate_signatures.csv",
    show_col_types = FALSE
)

annotation <- read_csv(
    paste0(
        "results/annotation/full_DESeq2/",
        "FA26_DESeq2_11661_gene_annotation_master.csv"
    ),
    show_col_types = FALSE
)

outdir <- "results/deseq2/top_genes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Extract top N significant UP genes
# ============================================================

top_up <- function(df, time, n = 100) {

    if (time == "2h") {

        z <- df %>%
            filter(
                sig_2h,
                direction_2h == "Up",
                !is.na(shrunk_LFC_2h)
            ) %>%
            arrange(
                desc(shrunk_LFC_2h),
                padj_2h
            ) %>%
            slice_head(n = n) %>%
            transmute(
                gene_id,
                shrunk_LFC = shrunk_LFC_2h,
                raw_LFC = LFC_2h,
                padj = padj_2h,
                support = support_2h
            )

    } else {

        z <- df %>%
            filter(
                sig_4h,
                direction_4h == "Up",
                !is.na(shrunk_LFC_4h)
            ) %>%
            arrange(
                desc(shrunk_LFC_4h),
                padj_4h
            ) %>%
            slice_head(n = n) %>%
            transmute(
                gene_id,
                shrunk_LFC = shrunk_LFC_4h,
                raw_LFC = LFC_4h,
                padj = padj_4h,
                support = support_4h
            )
    }

    z %>%
        mutate(
            rank = row_number()
        ) %>%
        left_join(
            annotation %>%
                dplyr::select(
                    gene_id,
                    preferred_name,
                    NCBI_product,
                    eggNOG_description,
                    GOs
                ),
            by = "gene_id"
        )
}

H2 <- top_up(high, "2h", 100)
H4 <- top_up(high, "4h", 100)
L2 <- top_up(low,  "2h", 100)
L4 <- top_up(low,  "4h", 100)

cat("Top-set sizes:\n")
cat("1 ug/g 2 h:   ", nrow(H2), "\n")
cat("1 ug/g 4 h:   ", nrow(H4), "\n")
cat("0.01 ug/g 2 h:", nrow(L2), "\n")
cat("0.01 ug/g 4 h:", nrow(L4), "\n\n")

# ============================================================
# 2. Pairwise overlap function
# ============================================================

pair_overlap <- function(a, b, name_a, name_b) {

    ids <- intersect(
        a$gene_id,
        b$gene_id
    )

    tibble(
        comparison =
            paste0(
                name_a,
                " vs ",
                name_b
            ),

        n_A = nrow(a),
        n_B = nrow(b),

        shared =
            length(ids),

        pct_A_shared =
            round(
                100 *
                length(ids) /
                nrow(a),
                1
            ),

        pct_B_shared =
            round(
                100 *
                length(ids) /
                nrow(b),
                1
            )
    )
}

overlap_summary <- bind_rows(

    pair_overlap(
        H2, L2,
        "1 ug/g 2h",
        "0.01 ug/g 2h"
    ),

    pair_overlap(
        H4, L4,
        "1 ug/g 4h",
        "0.01 ug/g 4h"
    ),

    pair_overlap(
        H2, H4,
        "1 ug/g 2h",
        "1 ug/g 4h"
    ),

    pair_overlap(
        L2, L4,
        "0.01 ug/g 2h",
        "0.01 ug/g 4h"
    ),

    pair_overlap(
        H2, L4,
        "1 ug/g 2h",
        "0.01 ug/g 4h"
    ),

    pair_overlap(
        L2, H4,
        "0.01 ug/g 2h",
        "1 ug/g 4h"
    )
)

# ============================================================
# 3. Long membership table
# ============================================================

membership <- bind_rows(

    H2 %>%
        transmute(
            gene_id,
            group = "PFBA1_2h",
            rank,
            shrunk_LFC,
            padj
        ),

    H4 %>%
        transmute(
            gene_id,
            group = "PFBA1_4h",
            rank,
            shrunk_LFC,
            padj
        ),

    L2 %>%
        transmute(
            gene_id,
            group = "PFBA001_2h",
            rank,
            shrunk_LFC,
            padj
        ),

    L4 %>%
        transmute(
            gene_id,
            group = "PFBA001_4h",
            rank,
            shrunk_LFC,
            padj
        )
)

multi_membership <- membership %>%
    group_by(
        gene_id
    ) %>%
    summarise(
        n_top100_sets =
            n_distinct(group),

        groups =
            paste(
                sort(
                    unique(group)
                ),
                collapse = "; "
            ),

        best_rank =
            min(rank),

        .groups = "drop"
    ) %>%

    filter(
        n_top100_sets >= 2
    ) %>%

    left_join(
        annotation %>%
            dplyr::select(
                gene_id,
                preferred_name,
                NCBI_product,
                eggNOG_description,
                GOs
            ),
        by = "gene_id"
    ) %>%

    arrange(
        desc(n_top100_sets),
        best_rank
    )

# ============================================================
# 4. Detailed pairwise shared tables
# ============================================================

make_shared_table <- function(
    a,
    b,
    prefix_a,
    prefix_b
) {

    ids <- intersect(
        a$gene_id,
        b$gene_id
    )

    if (length(ids) == 0) {
        return(tibble())
    }

    a %>%
        filter(
            gene_id %in% ids
        ) %>%
        dplyr::select(
            gene_id,
            rank_A = rank,
            LFC_A = shrunk_LFC,
            padj_A = padj
        ) %>%

        inner_join(
            b %>%
                filter(
                    gene_id %in% ids
                ) %>%
                dplyr::select(
                    gene_id,
                    rank_B = rank,
                    LFC_B = shrunk_LFC,
                    padj_B = padj
                ),
            by = "gene_id"
        ) %>%

        mutate(
            comparison_A = prefix_a,
            comparison_B = prefix_b,
            same_direction =
                sign(LFC_A) ==
                sign(LFC_B)
        ) %>%

        left_join(
            annotation %>%
                dplyr::select(
                    gene_id,
                    preferred_name,
                    NCBI_product,
                    eggNOG_description
                ),
            by = "gene_id"
        ) %>%

        arrange(
            pmax(
                rank_A,
                rank_B
            )
        )
}

shared_H2_L2 <- make_shared_table(
    H2, L2,
    "PFBA1_2h",
    "PFBA001_2h"
)

shared_H4_L4 <- make_shared_table(
    H4, L4,
    "PFBA1_4h",
    "PFBA001_4h"
)

shared_H2_L4 <- make_shared_table(
    H2, L4,
    "PFBA1_2h",
    "PFBA001_4h"
)

shared_H2_H4 <- make_shared_table(
    H2, H4,
    "PFBA1_2h",
    "PFBA1_4h"
)

shared_L2_L4 <- make_shared_table(
    L2, L4,
    "PFBA001_2h",
    "PFBA001_4h"
)

# ============================================================
# 5. Save
# ============================================================

write_csv(
    overlap_summary,
    file.path(
        outdir,
        "FA26_top100_overlap_summary.csv"
    )
)

write_csv(
    multi_membership,
    file.path(
        outdir,
        "FA26_top100_genes_in_multiple_sets.csv"
    )
)

write_csv(
    shared_H2_L2,
    file.path(
        outdir,
        "Top100_shared_PFBA1_vs_PFBA001_2h.csv"
    )
)

write_csv(
    shared_H4_L4,
    file.path(
        outdir,
        "Top100_shared_PFBA1_vs_PFBA001_4h.csv"
    )
)

write_csv(
    shared_H2_L4,
    file.path(
        outdir,
        "Top100_shared_PFBA1_2h_vs_PFBA001_4h.csv"
    )
)

write_csv(
    shared_H2_H4,
    file.path(
        outdir,
        "Top100_shared_PFBA1_2h_vs_4h.csv"
    )
)

write_csv(
    shared_L2_L4,
    file.path(
        outdir,
        "Top100_shared_PFBA001_2h_vs_4h.csv"
    )
)

# Save top-100 lists too
write_csv(
    H2,
    file.path(
        outdir,
        "PFBA1_top100_up_2h.csv"
    )
)

write_csv(
    H4,
    file.path(
        outdir,
        "PFBA1_top100_up_4h.csv"
    )
)

write_csv(
    L2,
    file.path(
        outdir,
        "PFBA001_top100_up_2h.csv"
    )
)

write_csv(
    L4,
    file.path(
        outdir,
        "PFBA001_top100_up_4h.csv"
    )
)

# ============================================================
# 6. Terminal report
# ============================================================

cat("\n============================================================\n")
cat("PAIRWISE TOP-100 OVERLAPS\n")
cat("============================================================\n\n")

print(
    overlap_summary,
    width = Inf
)

cat("\n============================================================\n")
cat("GENES PRESENT IN >=2 TOP-100 SETS\n")
cat("============================================================\n\n")

cat(
    "Number of recurrent genes:",
    nrow(multi_membership),
    "\n\n"
)

print(
    multi_membership %>%
        dplyr::select(
            gene_id,
            preferred_name,
            NCBI_product,
            n_top100_sets,
            groups,
            best_rank
        ),
    n = Inf,
    width = Inf
)

cat("\n============================================================\n")
cat("HIGH DOSE 2h vs LOW DOSE 4h SHARED TOP-100 GENES\n")
cat("============================================================\n\n")

print(
    shared_H2_L4,
    n = Inf,
    width = Inf
)

