#!/usr/bin/env Rscript

# ============================================================
# FA26 PFBA Hormesis Transcriptomics
# Step 04B: Topology-defined atlas regions
#
# Purpose:
#   Partition the fixed 3,000-gene UMAP atlas into three
#   topology-defined regions using hierarchical clustering
#   of the two-dimensional UMAP coordinates.
#
# Method:
#   Distance: Euclidean
#   Linkage:  single
#   Number of regions: 3
#
# This script reproduces the final atlas-region assignments
# used in the manuscript and downstream analyses.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
})

input_file <-
    "results/gene_atlas/FA26_gene_atlas_coordinates.csv"

output_file <-
    "results/gene_atlas/FA26_atlas3000_topology_zones_reproduced.csv"

reference_file <-
    "results/gene_atlas/FA26_atlas3000_topology_zones.csv"

# ------------------------------------------------------------
# 1. Load fixed 3,000-gene atlas coordinates
# ------------------------------------------------------------

atlas <- read.csv(
    input_file,
    stringsAsFactors = FALSE
)

required_cols <- c("gene_id", "UMAP1", "UMAP2")

if (!all(required_cols %in% names(atlas))) {
    stop(
        "Input atlas must contain: ",
        paste(required_cols, collapse = ", ")
    )
}

if (nrow(atlas) != 3000) {
    stop(
        "Expected 3,000 atlas genes; found ",
        nrow(atlas)
    )
}

# ------------------------------------------------------------
# 2. Hierarchical clustering of UMAP coordinates
# ------------------------------------------------------------

xy <- as.matrix(
    atlas[, c("UMAP1", "UMAP2")]
)

d <- dist(
    xy,
    method = "euclidean"
)

hc <- hclust(
    d,
    method = "single"
)

raw_cluster <- cutree(
    hc,
    k = 3
)

# ------------------------------------------------------------
# 3. Deterministic region labeling
#
# Label regions by spatial position so labels do not depend
# on arbitrary cluster numbering:
#
# Zone 1 = region with largest UMAP1 centroid
# Zone 2 = remaining region with largest UMAP2 centroid
# Zone 3 = remaining region
# ------------------------------------------------------------

tmp <- atlas %>%
    mutate(raw_cluster = raw_cluster) %>%
    group_by(raw_cluster) %>%
    summarise(
        centroid_UMAP1 = mean(UMAP1),
        centroid_UMAP2 = mean(UMAP2),
        .groups = "drop"
    )

zone1_cluster <-
    tmp$raw_cluster[which.max(tmp$centroid_UMAP1)]

remaining <-
    tmp %>%
    filter(raw_cluster != zone1_cluster)

zone2_cluster <-
    remaining$raw_cluster[
        which.max(remaining$centroid_UMAP2)
    ]

zone3_cluster <-
    remaining$raw_cluster[
        remaining$raw_cluster != zone2_cluster
    ]

zone_map <- c(
    setNames("Zone 1", zone1_cluster),
    setNames("Zone 2", zone2_cluster),
    setNames("Zone 3", zone3_cluster)
)

atlas$zone <-
    unname(zone_map[as.character(raw_cluster)])

# ------------------------------------------------------------
# 4. Save reproducible assignments
# ------------------------------------------------------------

out <- atlas %>%
    select(
        gene_id,
        UMAP1,
        UMAP2,
        zone
    )

write.csv(
    out,
    output_file,
    row.names = FALSE
)

# ------------------------------------------------------------
# 5. Report region summary
# ------------------------------------------------------------

cat("\nTopology-defined atlas regions:\n")
print(table(out$zone))

centroids <- out %>%
    group_by(zone) %>%
    summarise(
        n_genes = n(),
        centroid_UMAP1 = mean(UMAP1),
        centroid_UMAP2 = mean(UMAP2),
        .groups = "drop"
    )

cat("\nRegion centroids:\n")
print(centroids)

# ------------------------------------------------------------
# 6. Validate against archived final assignments
# ------------------------------------------------------------

if (file.exists(reference_file)) {

    reference <- read.csv(
        reference_file,
        stringsAsFactors = FALSE
    )

    validation <- out %>%
        select(gene_id, reproduced_zone = zone) %>%
        inner_join(
            reference %>%
                select(gene_id, archived_zone = zone),
            by = "gene_id"
        )

    agreement <-
        mean(
            validation$reproduced_zone ==
            validation$archived_zone
        )

    cat(
        "\nAgreement with archived final assignments:",
        sprintf("%.2f%%", agreement * 100),
        "\n"
    )

    if (
        nrow(validation) != 3000 ||
        agreement != 1
    ) {
        stop(
            "Reproduced topology zones do not exactly match ",
            "the archived final assignments."
        )
    }

    cat(
        "VALIDATION PASSED: all 3,000 gene assignments ",
        "were reproduced exactly.\n"
    )
}

cat("\nStep 04B complete.\n")
