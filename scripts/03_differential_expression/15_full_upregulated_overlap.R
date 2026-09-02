
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
})

cat("\n============================================================\n")
cat("FA26 FULL SIGNIFICANT UPREGULATED OVERLAP ANALYSIS\n")
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
    "results/annotation/full_DESeq2/FA26_DESeq2_11661_gene_annotation_master.csv",
    show_col_types = FALSE
)

outdir <- "results/deseq2/upregulated_overlap"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Extract ALL significant upregulated genes
# ============================================================

get_up <- function(dat, time) {

    if (time == "2h") {

        dat %>%
            filter(
                sig_2h %in% TRUE,
                direction_2h == "Up",
                !is.na(shrunk_LFC_2h)
            ) %>%
            arrange(
                desc(shrunk_LFC_2h),
                padj_2h
            ) %>%
            mutate(
                rank = row_number(),
                shrunk_LFC = shrunk_LFC_2h,
                raw_LFC = LFC_2h,
                padj = padj_2h,
                support = support_2h
            )

    } else {

        dat %>%
            filter(
                sig_4h %in% TRUE,
                direction_4h == "Up",
                !is.na(shrunk_LFC_4h)
            ) %>%
            arrange(
                desc(shrunk_LFC_4h),
                padj_4h
            ) %>%
            mutate(
                rank = row_number(),
                shrunk_LFC = shrunk_LFC_4h,
                raw_LFC = LFC_4h,
                padj = padj_4h,
                support = support_4h
            )
    }
}

H2 <- get_up(high, "2h")
H4 <- get_up(high, "4h")
L2 <- get_up(low,  "2h")
L4 <- get_up(low,  "4h")

cat("SIGNIFICANT UPREGULATED SET SIZES\n")
cat("1 ug/g 2 h:    ", nrow(H2), "\n")
cat("1 ug/g 4 h:    ", nrow(H4), "\n")
cat("0.01 ug/g 2 h: ", nrow(L2), "\n")
cat("0.01 ug/g 4 h: ", nrow(L4), "\n\n")

# ============================================================
# 2. Add top-100 membership
# ============================================================

H2 <- H2 %>% mutate(top100 = rank <= 100)
H4 <- H4 %>% mutate(top100 = rank <= 100)
L2 <- L2 %>% mutate(top100 = rank <= 100)
L4 <- L4 %>% mutate(top100 = rank <= 100)

# ============================================================
# 3. Pairwise overlap summary
# ============================================================

overlap_summary <- function(A, B, label) {

    ids <- intersect(
        A$gene_id,
        B$gene_id
    )

    tibble(
        comparison = label,
        n_A = nrow(A),
        n_B = nrow(B),
        shared = length(ids),
        pct_A_shared =
            ifelse(
                nrow(A) > 0,
                round(
                    100 * length(ids) / nrow(A),
                    1
                ),
                NA_real_
            ),
        pct_B_shared =
            ifelse(
                nrow(B) > 0,
                round(
                    100 * length(ids) / nrow(B),
                    1
                ),
                NA_real_
            )
    )
}

summary <- bind_rows(

    overlap_summary(
        H2, L2,
        "1 ug/g 2h vs 0.01 ug/g 2h"
    ),

    overlap_summary(
        H4, L4,
        "1 ug/g 4h vs 0.01 ug/g 4h"
    ),

    overlap_summary(
        H2, H4,
        "1 ug/g 2h vs 1 ug/g 4h"
    ),

    overlap_summary(
        L2, L4,
        "0.01 ug/g 2h vs 0.01 ug/g 4h"
    ),

    overlap_summary(
        H2, L4,
        "1 ug/g 2h vs 0.01 ug/g 4h"
    ),

    overlap_summary(
        L2, H4,
        "0.01 ug/g 2h vs 1 ug/g 4h"
    )
)

write_csv(
    summary,
    file.path(
        outdir,
        "FA26_all_significant_up_overlap_summary.csv"
    )
)

# ============================================================
# 4. Annotated shared-gene table builder
# ============================================================

make_shared <- function(A, B, suffixA, suffixB) {

    ids <- intersect(
        A$gene_id,
        B$gene_id
    )

    if (length(ids) == 0) {
        return(tibble())
    }

    ann_cols <- c(
        "gene_id",
        "preferred_name",
        "NCBI_product",
        "eggNOG_description",
        "GOs",
        "zone"
    )

    ann_cols <- ann_cols[
        ann_cols %in% names(A)
    ]

    Aout <- A %>%
        filter(
            gene_id %in% ids
        ) %>%
        select(
            all_of(ann_cols),
            rank,
            shrunk_LFC,
            raw_LFC,
            padj,
            support,
            top100
        )

    names(Aout)[names(Aout) == "rank"] <-
        paste0("rank_", suffixA)

    names(Aout)[names(Aout) == "shrunk_LFC"] <-
        paste0("shrunk_LFC_", suffixA)

    names(Aout)[names(Aout) == "raw_LFC"] <-
        paste0("raw_LFC_", suffixA)

    names(Aout)[names(Aout) == "padj"] <-
        paste0("padj_", suffixA)

    names(Aout)[names(Aout) == "support"] <-
        paste0("support_", suffixA)

    names(Aout)[names(Aout) == "top100"] <-
        paste0("top100_", suffixA)

    Bout <- B %>%
        filter(
            gene_id %in% ids
        ) %>%
        select(
            gene_id,
            rank,
            shrunk_LFC,
            raw_LFC,
            padj,
            support,
            top100
        )

    names(Bout)[names(Bout) == "rank"] <-
        paste0("rank_", suffixB)

    names(Bout)[names(Bout) == "shrunk_LFC"] <-
        paste0("shrunk_LFC_", suffixB)

    names(Bout)[names(Bout) == "raw_LFC"] <-
        paste0("raw_LFC_", suffixB)

    names(Bout)[names(Bout) == "padj"] <-
        paste0("padj_", suffixB)

    names(Bout)[names(Bout) == "support"] <-
        paste0("support_", suffixB)

    names(Bout)[names(Bout) == "top100"] <-
        paste0("top100_", suffixB)

    out <- left_join(
        Aout,
        Bout,
        by = "gene_id"
    )

    rankA <- paste0(
        "rank_",
        suffixA
    )

    rankB <- paste0(
        "rank_",
        suffixB
    )

    out %>%
        mutate(
            best_rank =
                pmin(
                    .data[[rankA]],
                    .data[[rankB]],
                    na.rm = TRUE
                ),

            worst_rank =
                pmax(
                    .data[[rankA]],
                    .data[[rankB]],
                    na.rm = TRUE
                )
        ) %>%
        arrange(
            worst_rank,
            best_rank
        )
}

# ============================================================
# 5. Build all biologically useful intersections
# ============================================================

shared_2h_crossdose <- make_shared(
    H2, L2,
    "PFBA1_2h",
    "PFBA001_2h"
)

shared_4h_crossdose <- make_shared(
    H4, L4,
    "PFBA1_4h",
    "PFBA001_4h"
)

persistent_highdose <- make_shared(
    H2, H4,
    "PFBA1_2h",
    "PFBA1_4h"
)

persistent_lowdose <- make_shared(
    L2, L4,
    "PFBA001_2h",
    "PFBA001_4h"
)

high2_low4 <- make_shared(
    H2, L4,
    "PFBA1_2h",
    "PFBA001_4h"
)

low2_high4 <- make_shared(
    L2, H4,
    "PFBA001_2h",
    "PFBA1_4h"
)

# ============================================================
# 6. Save intersection tables
# ============================================================

write_csv(
    shared_2h_crossdose,
    file.path(
        outdir,
        "FA26_all_up_shared_crossdose_2h.csv"
    )
)

write_csv(
    shared_4h_crossdose,
    file.path(
        outdir,
        "FA26_all_up_shared_crossdose_4h.csv"
    )
)

write_csv(
    persistent_highdose,
    file.path(
        outdir,
        "FA26_all_up_shared_PFBA1_2h_4h.csv"
    )
)

write_csv(
    persistent_lowdose,
    file.path(
        outdir,
        "FA26_all_up_shared_PFBA001_2h_4h.csv"
    )
)

write_csv(
    high2_low4,
    file.path(
        outdir,
        "FA26_all_up_shared_PFBA1_2h_PFBA001_4h.csv"
    )
)

write_csv(
    low2_high4,
    file.path(
        outdir,
        "FA26_all_up_shared_PFBA001_2h_PFBA1_4h.csv"
    )
)

# ============================================================
# 7. Genes recurring in 2+ of all four upregulated sets
# ============================================================

membership <- bind_rows(

    H2 %>%
        transmute(
            gene_id,
            set = "PFBA1_2h",
            rank,
            shrunk_LFC,
            padj
        ),

    H4 %>%
        transmute(
            gene_id,
            set = "PFBA1_4h",
            rank,
            shrunk_LFC,
            padj
        ),

    L2 %>%
        transmute(
            gene_id,
            set = "PFBA001_2h",
            rank,
            shrunk_LFC,
            padj
        ),

    L4 %>%
        transmute(
            gene_id,
            set = "PFBA001_4h",
            rank,
            shrunk_LFC,
            padj
        )
)

recurrent <- membership %>%
    group_by(
        gene_id
    ) %>%
    summarise(
        n_sets =
            n_distinct(set),

        sets =
            paste(
                sort(
                    unique(set)
                ),
                collapse = "; "
            ),

        best_rank =
            min(rank),

        .groups = "drop"
    ) %>%

    filter(
        n_sets >= 2
    ) %>%

    left_join(
        annotation %>%
            select(
                gene_id,
                preferred_name,
                NCBI_product,
                eggNOG_description,
                GOs
            ),
        by = "gene_id"
    ) %>%

    arrange(
        desc(n_sets),
        best_rank
    )

write_csv(
    recurrent,
    file.path(
        outdir,
        "FA26_all_significant_up_recurrent_genes.csv"
    )
)

# ============================================================
# 8. Terminal report
# ============================================================

cat("\n============================================================\n")
cat("ALL SIGNIFICANT UPREGULATED OVERLAPS\n")
cat("============================================================\n\n")

print(
    summary,
    n = Inf,
    width = Inf
)

cat("\n============================================================\n")
cat("CROSS-DOSE SHARED UPREGULATED GENES — 4 h\n")
cat("============================================================\n\n")

cat(
    "Genes:",
    nrow(shared_4h_crossdose),
    "\n\n"
)

print(
    shared_4h_crossdose,
    n = Inf,
    width = Inf
)

cat("\n============================================================\n")
cat("PERSISTENT 1 ug/g UPREGULATED GENES — 2 h AND 4 h\n")
cat("============================================================\n\n")

cat(
    "Genes:",
    nrow(persistent_highdose),
    "\n\n"
)

print(
    persistent_highdose,
    n = Inf,
    width = Inf
)

cat("\n============================================================\n")
cat("1 ug/g 2 h vs 0.01 ug/g 4 h UPREGULATED INTERSECTION\n")
cat("============================================================\n\n")

cat(
    "Genes:",
    nrow(high2_low4),
    "\n\n"
)

print(
    high2_low4,
    n = Inf,
    width = Inf
)

cat("\n============================================================\n")
cat("GENES SIGNIFICANTLY UPREGULATED IN >=2 SETS\n")
cat("============================================================\n\n")

cat(
    "Recurrent genes:",
    nrow(recurrent),
    "\n\n"
)

print(
    recurrent %>%
        select(
            gene_id,
            preferred_name,
            NCBI_product,
            n_sets,
            sets,
            best_rank
        ),
    n = Inf,
    width = Inf
)

cat("\n============================================================\n")
cat("STEP 15 COMPLETE\n")
cat("============================================================\n")

