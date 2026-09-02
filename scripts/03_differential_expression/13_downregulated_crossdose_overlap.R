
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

cat("\n============================================================\n")
cat("FA26 DOWNREGULATED CROSS-DOSE / TEMPORAL OVERLAP\n")
cat("============================================================\n\n")

high <- read_csv(
  "results/deseq2/signatures/FA26_PFBA1_candidate_signatures.csv",
  show_col_types = FALSE
)

low <- read_csv(
  "results/deseq2/signatures_lowdose/FA26_PFBA001_candidate_signatures.csv",
  show_col_types = FALSE
)

outdir <- "results/deseq2/downregulated_overlap"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

get_down <- function(dat, time) {

  if (time == "2h") {

    dat %>%
      filter(
        sig_2h %in% TRUE,
        direction_2h == "Down",
        !is.na(shrunk_LFC_2h)
      ) %>%
      arrange(shrunk_LFC_2h) %>%
      mutate(
        rank = row_number(),
        shrunk_LFC = shrunk_LFC_2h,
        raw_LFC = LFC_2h,
        padj = padj_2h
      )

  } else {

    dat %>%
      filter(
        sig_4h %in% TRUE,
        direction_4h == "Down",
        !is.na(shrunk_LFC_4h)
      ) %>%
      arrange(shrunk_LFC_4h) %>%
      mutate(
        rank = row_number(),
        shrunk_LFC = shrunk_LFC_4h,
        raw_LFC = LFC_4h,
        padj = padj_4h
      )
  }
}

high2_all <- get_down(high, "2h")
high4_all <- get_down(high, "4h")
low2_all  <- get_down(low,  "2h")
low4_all  <- get_down(low,  "4h")

cat("SIGNIFICANT DOWNREGULATED SET SIZES\n")
cat("1 ug/g 2 h:    ", nrow(high2_all), "\n")
cat("1 ug/g 4 h:    ", nrow(high4_all), "\n")
cat("0.01 ug/g 2 h: ", nrow(low2_all),  "\n")
cat("0.01 ug/g 4 h: ", nrow(low4_all),  "\n\n")

high2 <- high2_all %>% slice_head(n = 100)
high4 <- high4_all %>% slice_head(n = 100)
low2  <- low2_all  %>% slice_head(n = 100)
low4  <- low4_all  %>% slice_head(n = 100)

overlap_summary <- function(A, B, label) {

  shared <- intersect(A$gene_id, B$gene_id)

  tibble(
    comparison = label,
    n_A = nrow(A),
    n_B = nrow(B),
    shared = length(shared),
    pct_A_shared = ifelse(
      nrow(A) > 0,
      round(100 * length(shared) / nrow(A), 1),
      NA_real_
    ),
    pct_B_shared = ifelse(
      nrow(B) > 0,
      round(100 * length(shared) / nrow(B), 1),
      NA_real_
    )
  )
}

top100_summary <- bind_rows(

  overlap_summary(
    high2, high4,
    "1 ug/g 2h vs 1 ug/g 4h"
  ),

  overlap_summary(
    high4, low4,
    "1 ug/g 4h vs 0.01 ug/g 4h"
  ),

  overlap_summary(
    high2, low4,
    "1 ug/g 2h vs 0.01 ug/g 4h"
  )
)

full_summary <- bind_rows(

  overlap_summary(
    high2_all, high4_all,
    "1 ug/g 2h vs 1 ug/g 4h"
  ),

  overlap_summary(
    high4_all, low4_all,
    "1 ug/g 4h vs 0.01 ug/g 4h"
  ),

  overlap_summary(
    high2_all, low4_all,
    "1 ug/g 2h vs 0.01 ug/g 4h"
  )
)

make_shared <- function(A, B, suffixA, suffixB) {

  ids <- intersect(A$gene_id, B$gene_id)

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
    filter(gene_id %in% ids) %>%
    select(
      all_of(ann_cols),
      rank,
      shrunk_LFC,
      raw_LFC,
      padj,
      overall_support
    )

  names(Aout)[names(Aout) == "rank"] <-
    paste0("rank_", suffixA)

  names(Aout)[names(Aout) == "shrunk_LFC"] <-
    paste0("shrunk_LFC_", suffixA)

  names(Aout)[names(Aout) == "raw_LFC"] <-
    paste0("raw_LFC_", suffixA)

  names(Aout)[names(Aout) == "padj"] <-
    paste0("padj_", suffixA)

  names(Aout)[names(Aout) == "overall_support"] <-
    paste0("support_", suffixA)

  Bout <- B %>%
    filter(gene_id %in% ids) %>%
    select(
      gene_id,
      rank,
      shrunk_LFC,
      raw_LFC,
      padj,
      overall_support
    )

  names(Bout)[names(Bout) == "rank"] <-
    paste0("rank_", suffixB)

  names(Bout)[names(Bout) == "shrunk_LFC"] <-
    paste0("shrunk_LFC_", suffixB)

  names(Bout)[names(Bout) == "raw_LFC"] <-
    paste0("raw_LFC_", suffixB)

  names(Bout)[names(Bout) == "padj"] <-
    paste0("padj_", suffixB)

  names(Bout)[names(Bout) == "overall_support"] <-
    paste0("support_", suffixB)

  left_join(
    Aout,
    Bout,
    by = "gene_id"
  )
}

top_high_temporal <- make_shared(
  high2, high4,
  "PFBA1_2h",
  "PFBA1_4h"
)

top_crossdose_4h <- make_shared(
  high4, low4,
  "PFBA1_4h",
  "PFBA001_4h"
)

top_high2_low4 <- make_shared(
  high2, low4,
  "PFBA1_2h",
  "PFBA001_4h"
)

full_high_temporal <- make_shared(
  high2_all, high4_all,
  "PFBA1_2h",
  "PFBA1_4h"
)

full_crossdose_4h <- make_shared(
  high4_all, low4_all,
  "PFBA1_4h",
  "PFBA001_4h"
)

full_high2_low4 <- make_shared(
  high2_all, low4_all,
  "PFBA1_2h",
  "PFBA001_4h"
)

write_csv(
  top100_summary,
  file.path(outdir, "FA26_top100_down_overlap_summary.csv")
)

write_csv(
  full_summary,
  file.path(outdir, "FA26_all_significant_down_overlap_summary.csv")
)

write_csv(
  top_high_temporal,
  file.path(outdir, "FA26_top100_down_shared_PFBA1_2h_4h.csv")
)

write_csv(
  top_crossdose_4h,
  file.path(outdir, "FA26_top100_down_shared_crossdose_4h.csv")
)

write_csv(
  top_high2_low4,
  file.path(outdir, "FA26_top100_down_shared_PFBA1_2h_PFBA001_4h.csv")
)

write_csv(
  full_high_temporal,
  file.path(outdir, "FA26_all_down_shared_PFBA1_2h_4h.csv")
)

write_csv(
  full_crossdose_4h,
  file.path(outdir, "FA26_all_down_shared_crossdose_4h.csv")
)

write_csv(
  full_high2_low4,
  file.path(outdir, "FA26_all_down_shared_PFBA1_2h_PFBA001_4h.csv")
)

cat("\n============================================================\n")
cat("TOP-100 DOWNREGULATED OVERLAPS\n")
cat("============================================================\n\n")

print(top100_summary)

cat("\n============================================================\n")
cat("ALL SIGNIFICANT DOWNREGULATED GENE OVERLAPS\n")
cat("============================================================\n\n")

print(full_summary)

cat("\n============================================================\n")
cat("TOP-100 DOWN — HIGH DOSE PERSISTENT (2 h AND 4 h)\n")
cat("============================================================\n\n")

if (nrow(top_high_temporal) > 0) {
  print(top_high_temporal)
} else {
  cat("No shared genes.\n")
}

cat("\n============================================================\n")
cat("TOP-100 DOWN — CROSS-DOSE SHARED AT 4 h\n")
cat("============================================================\n\n")

if (nrow(top_crossdose_4h) > 0) {
  print(top_crossdose_4h)
} else {
  cat("No shared genes.\n")
}

cat("\n============================================================\n")
cat("TOP-100 DOWN — 1 ug/g 2 h vs 0.01 ug/g 4 h\n")
cat("============================================================\n\n")

if (nrow(top_high2_low4) > 0) {
  print(top_high2_low4)
} else {
  cat("No shared genes.\n")
}

cat("\n============================================================\n")
cat("STEP 13 COMPLETE\n")
cat("============================================================\n")

