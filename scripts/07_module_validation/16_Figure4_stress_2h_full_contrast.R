
options(stringsAsFactors = FALSE)

library(ggplot2)

outdir <- "results/figure4/stress_2h"
figdir <- "figures/working"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

stress <- read.csv(
  "results/figure4/FA26_stress_candidate_genes.csv",
  check.names = FALSE
)

low <- read.csv(
  "results/deseq2/PFBA_0.01_vs_Control_2h.csv",
  check.names = FALSE
)

high <- read.csv(
  "results/deseq2/PFBA_1_vs_Control_2h.csv",
  check.names = FALSE
)

cat("LOW columns:\n")
print(names(low))

cat("\nHIGH columns:\n")
print(names(high))

# normalize gene ID column
find_id <- function(x) {
  hit <- names(x)[tolower(names(x)) %in% c("gene_id","geneid","gene","id")]
  if (length(hit) == 0) stop("Could not identify gene ID column.")
  hit[1]
}

names(low)[names(low) == find_id(low)] <- "gene_id"
names(high)[names(high) == find_id(high)] <- "gene_id"

# identify DESeq2 columns
find_col <- function(x, candidates) {
  hit <- intersect(candidates, names(x))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

lfc_low  <- find_col(low,  c("log2FoldChange","LFC","shrunk_LFC"))
padj_low <- find_col(low,  c("padj","adj.P.Val","FDR"))

lfc_high  <- find_col(high, c("log2FoldChange","LFC","shrunk_LFC"))
padj_high <- find_col(high, c("padj","adj.P.Val","FDR"))

if (is.null(lfc_low) || is.null(padj_low) ||
    is.null(lfc_high) || is.null(padj_high)) {
  stop("Could not identify LFC/padj columns.")
}

low2 <- low[, c("gene_id", lfc_low, padj_low)]
high2 <- high[, c("gene_id", lfc_high, padj_high)]

names(low2) <- c("gene_id","Low_LFC_2h","Low_padj_2h")
names(high2) <- c("gene_id","High_LFC_2h","High_padj_2h")

# one row per stress gene
collapse_text <- function(x) {
  x <- unique(x[!is.na(x) & x != "" & x != "-"])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

meta <- do.call(
  rbind,
  lapply(unique(stress$gene_id), function(g) {
    x <- stress[stress$gene_id == g, , drop = FALSE]

    data.frame(
      gene_id = g,
      display_name = {
        y <- x$display_name[!is.na(x$display_name) & x$display_name != ""]
        if (length(y)) y[1] else g
      },
      product = {
        y <- x$product[!is.na(x$product) & x$product != ""]
        if (length(y)) y[1] else NA_character_
      },
      stress_category = collapse_text(x$screen_category),
      stringsAsFactors = FALSE
    )
  })
)

dat <- merge(meta, low2, by = "gene_id", all.x = TRUE)
dat <- merge(dat, high2, by = "gene_id", all.x = TRUE)

dat$Low_sig_2h  <- !is.na(dat$Low_padj_2h)  & dat$Low_padj_2h  < 0.05
dat$High_sig_2h <- !is.na(dat$High_padj_2h) & dat$High_padj_2h < 0.05

dat$response_2h <- "Neither significant"
dat$response_2h[ dat$High_sig_2h & !dat$Low_sig_2h] <- "1 ug/g only"
dat$response_2h[!dat$High_sig_2h &  dat$Low_sig_2h] <- "0.01 ug/g only"
dat$response_2h[ dat$High_sig_2h &  dat$Low_sig_2h] <- "Both doses"

cat("\n============================================================\n")
cat("FULL-CONTRAST 2 h STRESS RESPONSE\n")
cat("============================================================\n\n")

cat("Stress genes tested:", nrow(dat), "\n")
cat("Significant at 1 ug/g    :", sum(dat$High_sig_2h), "\n")
cat("Significant at 0.01 ug/g :", sum(dat$Low_sig_2h), "\n")
cat("Shared significant      :", sum(dat$High_sig_2h & dat$Low_sig_2h), "\n\n")

cat("Response classes:\n")
print(table(dat$response_2h))

# high-significant genes with actual low-dose estimates
highsig <- dat[dat$High_sig_2h, ]
highsig <- highsig[order(highsig$High_padj_2h), ]

cat("\n============================================================\n")
cat("1 ug/g STRESS GENES WITH 0.01 ug/g EFFECT ESTIMATES\n")
cat("============================================================\n\n")

print(
  highsig[, c(
    "gene_id","display_name","stress_category",
    "High_LFC_2h","High_padj_2h",
    "Low_LFC_2h","Low_padj_2h",
    "Low_sig_2h","response_2h"
  )],
  row.names = FALSE
)

# low-significant genes, if any
lowsig <- dat[dat$Low_sig_2h, ]

cat("\n============================================================\n")
cat("0.01 ug/g STRESS GENES AT 2 h\n")
cat("============================================================\n\n")

if (nrow(lowsig) == 0) {
  cat("NONE\n")
} else {
  print(
    lowsig[, c(
      "gene_id","display_name","stress_category",
      "Low_LFC_2h","Low_padj_2h",
      "High_LFC_2h","High_padj_2h"
    )],
    row.names = FALSE
  )
}

write.csv(
  dat,
  file.path(outdir, "FA26_stress_fullcontrast_high_vs_low_2h.csv"),
  row.names = FALSE
)

write.csv(
  highsig,
  file.path(outdir, "FA26_high_stress_genes_with_low_effects_2h.csv"),
  row.names = FALSE
)

# scatter
plotdat <- dat[
  !is.na(dat$Low_LFC_2h) &
  !is.na(dat$High_LFC_2h),
]

plotdat$label <- ifelse(
  plotdat$High_sig_2h | plotdat$Low_sig_2h,
  plotdat$display_name,
  ""
)

p <- ggplot(plotdat, aes(Low_LFC_2h, High_LFC_2h)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.45) +
  geom_point(aes(shape = response_2h), size = 3, alpha = 0.85) +
  labs(
    x = expression("0.01 " * mu * "g g"^{-1} * " PFBA, 2 h log"[2] * "FC"),
    y = expression("1 " * mu * "g g"^{-1} * " PFBA, 2 h log"[2] * "FC"),
    shape = "2 h significance",
    title = "Canonical stress-associated transcription at 2 h"
  ) +
  theme_classic(base_size = 13)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 3.2,
      max.overlaps = 30,
      box.padding = 0.4,
      point.padding = 0.2
    )
}

ggsave(
  file.path(figdir, "Figure4_stress_2h_fullcontrast_exploratory.pdf"),
  p,
  width = 7.2,
  height = 6.3
)

ggsave(
  file.path(figdir, "Figure4_stress_2h_fullcontrast_exploratory.tiff"),
  p,
  width = 7.2,
  height = 6.3,
  dpi = 600,
  compression = "lzw"
)

cat("\nSaved full comparison and exploratory plot.\n")

