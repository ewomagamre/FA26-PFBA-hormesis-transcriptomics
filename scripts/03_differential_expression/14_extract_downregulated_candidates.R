
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

cat("\n============================================================\n")
cat("FA26 DOWNREGULATED CANDIDATE EXTRACTION\n")
cat("============================================================\n\n")

indir <- "results/deseq2/downregulated_overlap"
outdir <- "results/deseq2/candidate_genes"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 1. Load the two definitive downregulated intersection sets
# ============================================================

crossdose <- read_csv(
    file.path(
        indir,
        "FA26_all_down_shared_crossdose_4h.csv"
    ),
    show_col_types = FALSE
)

persistent <- read_csv(
    file.path(
        indir,
        "FA26_all_down_shared_PFBA1_2h_4h.csv"
    ),
    show_col_types = FALSE
)

cat("Cross-dose downregulated at 4 h:", nrow(crossdose), "\n")
cat("Persistent high-dose downregulated:", nrow(persistent), "\n")

# ============================================================
# 2. Check overlap between the 8 and 12 gene sets
# ============================================================

shared_between_classes <- intersect(
    crossdose$gene_id,
    persistent$gene_id
)

cat(
    "Genes present in BOTH downregulated candidate classes:",
    length(shared_between_classes),
    "\n\n"
)

# ============================================================
# 3. Add candidate class
# ============================================================

crossdose <- crossdose %>%
    mutate(
        candidate_class =
            "Cross-dose downregulated at 4 h"
    )

persistent <- persistent %>%
    mutate(
        candidate_class =
            "Persistent 1 ug/g downregulated"
    )

# ============================================================
# 4. Save individual candidate lists
# ============================================================

write_csv(
    crossdose,
    file.path(
        outdir,
        "FA26_crossdose_downregulated_4h_candidates.csv"
    )
)

write_csv(
    persistent,
    file.path(
        outdir,
        "FA26_PFBA1_persistent_downregulated_candidates.csv"
    )
)

# ============================================================
# 5. Build combined nonredundant candidate table
# ============================================================

combined <- bind_rows(
    crossdose %>%
        select(
            gene_id,
            preferred_name,
            NCBI_product,
            eggNOG_description,
            GOs,
            zone,
            candidate_class
        ),

    persistent %>%
        select(
            gene_id,
            preferred_name,
            NCBI_product,
            eggNOG_description,
            GOs,
            zone,
            candidate_class
        )
) %>%
    group_by(
        gene_id,
        preferred_name,
        NCBI_product,
        eggNOG_description,
        GOs,
        zone
    ) %>%
    summarise(
        candidate_class =
            paste(
                unique(candidate_class),
                collapse = "; "
            ),
        .groups = "drop"
    )

write_csv(
    combined,
    file.path(
        outdir,
        "FA26_downregulated_candidate_master.csv"
    )
)

# ============================================================
# 6. Save genes belonging to BOTH classes
# ============================================================

if (length(shared_between_classes) > 0) {

    dual_class <- combined %>%
        filter(
            gene_id %in%
                shared_between_classes
        )

} else {

    dual_class <- tibble(
        gene_id = character()
    )
}

write_csv(
    dual_class,
    file.path(
        outdir,
        "FA26_downregulated_candidates_in_both_classes.csv"
    )
)

# ============================================================
# 7. Print cross-dose candidates
# ============================================================

cat("\n============================================================\n")
cat("8 CROSS-DOSE DOWNREGULATED GENES — 4 h\n")
cat("============================================================\n\n")

cross_cols <- intersect(
    c(
        "gene_id",
        "preferred_name",
        "NCBI_product",
        "eggNOG_description",
        "zone",
        "shrunk_LFC_PFBA1_4h",
        "padj_PFBA1_4h",
        "support_PFBA1_4h",
        "shrunk_LFC_PFBA001_4h",
        "padj_PFBA001_4h",
        "support_PFBA001_4h"
    ),
    names(crossdose)
)

print(
    crossdose %>%
        select(
            all_of(cross_cols)
        ),
    n = Inf,
    width = Inf
)

# ============================================================
# 8. Print persistent high-dose candidates
# ============================================================

cat("\n============================================================\n")
cat("12 PERSISTENT 1 ug/g DOWNREGULATED GENES — 2 h AND 4 h\n")
cat("============================================================\n\n")

persistent_cols <- intersect(
    c(
        "gene_id",
        "preferred_name",
        "NCBI_product",
        "eggNOG_description",
        "zone",
        "shrunk_LFC_PFBA1_2h",
        "padj_PFBA1_2h",
        "support_PFBA1_2h",
        "shrunk_LFC_PFBA1_4h",
        "padj_PFBA1_4h",
        "support_PFBA1_4h"
    ),
    names(persistent)
)

print(
    persistent %>%
        select(
            all_of(persistent_cols)
        ),
    n = Inf,
    width = Inf
)

# ============================================================
# 9. Print overlap between candidate classes
# ============================================================

cat("\n============================================================\n")
cat("GENES BELONGING TO BOTH DOWNREGULATED CLASSES\n")
cat("============================================================\n\n")

if (length(shared_between_classes) == 0) {

    cat("No genes occur in both candidate classes.\n")

} else {

    print(
        combined %>%
            filter(
                gene_id %in%
                    shared_between_classes
            ),
        n = Inf,
        width = Inf
    )
}

cat("\n============================================================\n")
cat("CANDIDATE SUMMARY\n")
cat("============================================================\n")

cat(
    "Cross-dose 4 h downregulated:",
    nrow(crossdose),
    "\n"
)

cat(
    "Persistent high-dose downregulated:",
    nrow(persistent),
    "\n"
)

cat(
    "Overlap between classes:",
    length(shared_between_classes),
    "\n"
)

cat(
    "Unique downregulated candidate genes:",
    nrow(combined),
    "\n"
)

cat("\nOutputs:\n")
cat("  results/deseq2/candidate_genes/FA26_crossdose_downregulated_4h_candidates.csv\n")
cat("  results/deseq2/candidate_genes/FA26_PFBA1_persistent_downregulated_candidates.csv\n")
cat("  results/deseq2/candidate_genes/FA26_downregulated_candidate_master.csv\n")
cat("  results/deseq2/candidate_genes/FA26_downregulated_candidates_in_both_classes.csv\n")


