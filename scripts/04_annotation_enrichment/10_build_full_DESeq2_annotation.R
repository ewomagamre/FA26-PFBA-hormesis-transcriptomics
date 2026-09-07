
suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
})

cat("\n============================================================\n")
cat("FA26 STEP 10 — FULL DESeq2 FUNCTIONAL ANNOTATION\n")
cat("============================================================\n\n")

# ============================================================
# 1. Paths
# ============================================================

dds_file <-
    "objects/FA26_dds_condition_clean.rds"

gff_file <-
    "../fa_genome/FA_annotation.gff"

protein_fasta <-
    "../fa_genome/FA_proteins_clean.fa"

eggnog_file <-
    "../fa_genome/FA_annotation.emapper.annotations"

outdir <-
    "results/annotation/full_DESeq2"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# 2. DESeq2 universe
# ============================================================

dds <- readRDS(dds_file)

universe <- tibble(
    gene_id = rownames(dds)
)

cat("DESeq2 universe:", nrow(universe), "genes\n")

stopifnot(
    nrow(universe) == 11661
)

# ============================================================
# 3. Parse protein FASTA
#
# Current FA26 cleaned FASTA headers correspond to transcript IDs.
# ============================================================

read_fasta_lengths <- function(f) {

    lines <- readLines(f)

    headers <- grep(
        "^>",
        lines
    )

    ends <- c(
        headers[-1] - 1,
        length(lines)
    )

    records <- lapply(
        seq_along(headers),
        function(i) {

            header <- sub(
                "^>",
                "",
                lines[
                    headers[i]
                ]
            )

            transcript_id <- strsplit(
                header,
                "\\s+"
            )[[1]][1]

            seq_lines <- lines[
                (headers[i] + 1):
                ends[i]
            ]

            seq <- paste0(
                seq_lines,
                collapse = ""
            )

            # Count amino-acid characters while ignoring formatting dots.
            aa_length <- nchar(
                gsub(
                    "[^A-Za-z*]",
                    "",
                    seq
                )
            )

            tibble(
                transcript_id =
                    transcript_id,
                aa_length =
                    aa_length
            )
        }
    )

    bind_rows(records)
}

protein_lengths <- read_fasta_lengths(
    protein_fasta
)

cat(
    "Protein FASTA records:",
    nrow(protein_lengths),
    "\n"
)

# ============================================================
# 4. Parse GFF CDS records
#
# CDS lines link:
#   gene
#   transcript_id
#   protein_id
#   product
# ============================================================

gff <- read.delim(
    gff_file,
    comment.char = "#",
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
)

colnames(gff) <- c(
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
)

extract_attr <- function(x, key) {

    pattern <- paste0(
        "(?:^|;)",
        key,
        "=([^;]+)"
    )

    out <- sub(
        paste0(
            ".*",
            pattern,
            ".*"
        ),
        "\\1",
        x
    )

    missing <- !grepl(
        pattern,
        x
    )

    out[missing] <- NA_character_

    out
}

cds <- gff %>%
    filter(
        type == "CDS"
    ) %>%
    transmute(
        gene_id =
            extract_attr(
                attributes,
                "gene"
            ),

        transcript_id =
            extract_attr(
                attributes,
                "Parent"
            ),

        protein_id =
            extract_attr(
                attributes,
                "protein_id"
            ),

        NCBI_product =
            extract_attr(
                attributes,
                "product"
            )
    ) %>%
    mutate(
        transcript_id =
            sub(
                "^rna-",
                "",
                transcript_id
            )
    ) %>%
    filter(
        !is.na(gene_id),
        !is.na(transcript_id)
    ) %>%
    distinct()

cat(
    "Unique CDS gene-transcript mappings:",
    nrow(cds),
    "\n"
)

# ============================================================
# 5. Representative transcript per gene
#
# Choose longest available protein transcript.
# ============================================================

mapping <- cds %>%
    left_join(
        protein_lengths,
        by = "transcript_id"
    ) %>%
    filter(
        gene_id %in%
            universe$gene_id
    )

representative <- mapping %>%
    group_by(
        gene_id
    ) %>%
    arrange(
        desc(
            aa_length
        ),
        transcript_id,
        .by_group = TRUE
    ) %>%
    slice_head(
        n = 1
    ) %>%
    ungroup()

cat(
    "DESeq2 genes with representative CDS/protein:",
    nrow(representative),
    "\n"
)

# ============================================================
# 6. Parse full previous eggNOG output
# ============================================================

egg_lines <- readLines(
    eggnog_file
)

header_index <- grep(
    "^#query\t",
    egg_lines
)

if (length(header_index) != 1) {
    stop(
        "Could not identify eggNOG #query header."
    )
}

header <- strsplit(
    sub(
        "^#",
        "",
        egg_lines[
            header_index
        ]
    ),
    "\t"
)[[1]]

data_lines <- egg_lines[
    !grepl(
        "^#",
        egg_lines
    ) &
    nzchar(
        egg_lines
    )
]

eggnog <- read.delim(
    text = paste(
        data_lines,
        collapse = "\n"
    ),
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (ncol(eggnog) < length(header)) {

    stop(
        paste(
            "eggNOG field count mismatch:",
            ncol(eggnog),
            "columns vs",
            length(header),
            "header fields"
        )
    )
}

colnames(eggnog)[
    seq_along(header)
] <- header

eggnog <- eggnog %>%
    rename(
        transcript_id =
            query,

        eggNOG_description =
            Description,

        preferred_name =
            Preferred_name
    )

cat(
    "Full eggNOG transcript annotations:",
    nrow(eggnog),
    "\n"
)

# ============================================================
# 7. Merge representative genes with eggNOG
# ============================================================

annotation <- universe %>%

    left_join(
        representative %>%
            dplyr::select(
                gene_id,
                transcript_id,
                protein_id,
                NCBI_product,
                aa_length
            ),
        by = "gene_id"
    ) %>%

    left_join(
        eggnog,
        by = "transcript_id"
    )

# ============================================================
# 8. Standardize missing annotation
# ============================================================

annotation <- annotation %>%
    mutate(
        has_representative_protein =
            !is.na(
                transcript_id
            ),

        has_eggnog =
            !is.na(
                seed_ortholog
            ) &
            seed_ortholog != "-",

        has_description =
            !is.na(
                eggNOG_description
            ) &
            eggNOG_description != "-" &
            eggNOG_description != "",

        has_GO =
            !is.na(
                GOs
            ) &
            GOs != "-" &
            GOs != "",

        has_KEGG =
            !is.na(
                KEGG_ko
            ) &
            KEGG_ko != "-" &
            KEGG_ko != ""
    )

# ============================================================
# 9. Coverage summary
# ============================================================

coverage <- tibble(

    metric = c(
        "DESeq2 genes",
        "Representative protein mapping",
        "Existing eggNOG match",
        "eggNOG description",
        "GO annotation",
        "KEGG KO annotation"
    ),

    n_genes = c(
        nrow(annotation),

        sum(
            annotation$has_representative_protein
        ),

        sum(
            annotation$has_eggnog
        ),

        sum(
            annotation$has_description
        ),

        sum(
            annotation$has_GO
        ),

        sum(
            annotation$has_KEGG
        )
    )
) %>%

    mutate(
        percent =
            round(
                100 *
                n_genes /
                nrow(annotation),
                1
            )
    )

# ============================================================
# 10. Missing eggNOG representative transcripts
# ============================================================

missing_eggnog <- annotation %>%
    filter(
        has_representative_protein,
        !has_eggnog
    ) %>%
    dplyr::select(
        gene_id,
        transcript_id,
        protein_id,
        NCBI_product,
        aa_length
    )

# ============================================================
# 11. Save
# ============================================================

write_csv(
    annotation,
    file.path(
        outdir,
        "FA26_DESeq2_11661_gene_annotation_master.csv"
    )
)

saveRDS(
    annotation,
    file.path(
        outdir,
        "FA26_DESeq2_11661_gene_annotation_master.rds"
    )
)

write_csv(
    coverage,
    file.path(
        outdir,
        "FA26_DESeq2_annotation_coverage.csv"
    )
)

write_csv(
    missing_eggnog,
    file.path(
        outdir,
        "FA26_DESeq2_missing_existing_eggnog.csv"
    )
)

# ============================================================
# 12. GO membership table for full DESeq2 universe
# ============================================================

go_membership <- annotation %>%
    dplyr::select(
        gene_id,
        GOs
    ) %>%
    filter(
        !is.na(GOs),
        GOs != "",
        GOs != "-"
    ) %>%
    separate_rows(
        GOs,
        sep = ","
    ) %>%
    mutate(
        GO_ID =
            trimws(
                GOs
            )
    ) %>%
    filter(
        grepl(
            "^GO:[0-9]+$",
            GO_ID
        )
    ) %>%
    dplyr::select(
        gene_id,
        GO_ID
    ) %>%
    distinct()

write_csv(
    go_membership,
    file.path(
        outdir,
        "FA26_DESeq2_GO_membership.csv"
    )
)

cat("\n============================================================\n")
cat("STEP 10 COMPLETE\n")
cat("============================================================\n\n")

cat("ANNOTATION COVERAGE\n\n")
print(coverage)

cat(
    "\nRepresentative transcripts missing from existing eggNOG:",
    nrow(missing_eggnog),
    "\n"
)

cat(
    "Unique gene-GO memberships:",
    nrow(go_membership),
    "\n"
)

