###############################################################
# FA26 RNA-seq Gene Atlas Pipeline
# Script: 01_preprocessing_filtering.R
#
# Purpose:
# Prepare the primary 24-sample dataset by excluding the
# pre-exposure benchmark (G0) samples.
#
# Author: Wayne Omagamre
# Project: FA26 RNA-seq
###############################################################

#-----------------------------#
# Load libraries
#-----------------------------#

library(DESeq2)

#-----------------------------#
# Read input objects
#-----------------------------#

dds27      <- readRDS("../dds_filtered.rds")
vsd27      <- readRDS("../vsd.rds")
metadata27 <- readRDS("../metadata_clean.rds")

#-----------------------------#
# Initial summary
#-----------------------------#

cat("\n=========================================\n")
cat("Original filtered dataset\n")
cat("=========================================\n")

cat("Genes   :", nrow(dds27), "\n")
cat("Samples :", ncol(dds27), "\n\n")

print(table(metadata27$Condition))

#-----------------------------#
# Remove pre-exposure benchmark (G0)
#-----------------------------#

keep_samples <- metadata27$Condition != "Control_0h"

metadata24 <- metadata27[keep_samples, ]

dds24 <- dds27[, rownames(metadata24)]

vsd24 <- vsd27[, rownames(metadata24)]

#-----------------------------#
# Verify
#-----------------------------#

stopifnot(
    all(colnames(dds24) == rownames(metadata24)),
    all(colnames(vsd24) == rownames(metadata24))
)

cat("\n=========================================\n")
cat("Primary analysis dataset\n")
cat("=========================================\n")

cat("Genes   :", nrow(dds24), "\n")
cat("Samples :", ncol(dds24), "\n\n")

print(table(metadata24$Condition))

#-----------------------------#
# Save objects
#-----------------------------#

dir.create("objects", showWarnings = FALSE)

saveRDS(
    dds24,
    "objects/FA26_dds_primary24.rds"
)

saveRDS(
    vsd24,
    "objects/FA26_vsd_primary24.rds"
)

saveRDS(
    metadata24,
    "objects/FA26_metadata_primary24.rds"
)

cat("\nObjects successfully saved.\n")
cat("Pipeline Step 01 complete.\n")
