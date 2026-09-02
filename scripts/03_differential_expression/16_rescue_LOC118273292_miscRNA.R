
suppressPackageStartupMessages({
    library(Biostrings)
    library(dplyr)
    library(readr)
})

cat("\n============================================================\n")
cat("FA26 RESCUE — LOC118273292 / XR_007705569.1\n")
cat("============================================================\n\n")

genome_file <- "../fa_genome/FA_GCF_genomic.fna"
gff_file <- "../fa_genome/FA_annotation.gff"

outdir <- "results/annotation/candidate_rescue/LOC118273292"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

target_gene <- "LOC118273292"
target_tx <- "XR_007705569.1"

# ------------------------------------------------------------
# 1. Load genome
# ------------------------------------------------------------

genome <- readDNAStringSet(genome_file)

# Standardize FASTA headers to accession only
names(genome) <- sub(" .*", "", names(genome))

cat("Genome sequences:", length(genome), "\n")

# ------------------------------------------------------------
# 2. Read GFF as plain table
# ------------------------------------------------------------

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

target <- gff %>%
    filter(
        grepl(target_gene, attributes) |
        grepl(target_tx, attributes)
    )

exons <- target %>%
    filter(
        type == "exon",
        grepl(target_tx, attributes)
    )

if (nrow(exons) == 0) {
    stop("No exons found for target transcript.")
}

chr <- unique(exons$seqid)
strand_val <- unique(exons$strand)

if (length(chr) != 1) {
    stop("Unexpected multiple chromosomes for target transcript.")
}

if (length(strand_val) != 1) {
    stop("Unexpected multiple strand values for target transcript.")
}

cat("Chromosome:", chr, "\n")
cat("Strand:", strand_val, "\n")
cat("Exons:", nrow(exons), "\n\n")

# ------------------------------------------------------------
# 3. Order exons in transcript direction
# ------------------------------------------------------------

if (strand_val == "+") {

    exons <- exons %>%
        arrange(start)

} else {

    exons <- exons %>%
        arrange(desc(start))
}

# ------------------------------------------------------------
# 4. Match chromosome name
# ------------------------------------------------------------

if (!(chr %in% names(genome))) {

    cat("\nGenome FASTA headers:\n")
    print(head(names(genome), 20))

    stop(
        paste(
            "Chromosome",
            chr,
            "not found in genome FASTA."
        )
    )
}

chr_seq <- genome[[chr]]

# ------------------------------------------------------------
# 5. Extract exon sequences
# ------------------------------------------------------------

exon_seqs <- DNAStringSet(
    lapply(
        seq_len(nrow(exons)),
        function(i) {

            s <- subseq(
                chr_seq,
                start = exons$start[i],
                end = exons$end[i]
            )

            if (strand_val == "-") {
                s <- reverseComplement(s)
            }

            s
        }
    )
)

names(exon_seqs) <- paste0(
    "exon_",
    seq_along(exon_seqs),
    "|",
    exons$start,
    "-",
    exons$end
)

transcript <- DNAString(
    paste0(
        as.character(exon_seqs),
        collapse = ""
    )
)

# ------------------------------------------------------------
# 6. Basic transcript properties
# ------------------------------------------------------------

tx_length <- length(transcript)

base_freq <- letterFrequency(
    transcript,
    letters = c("A", "C", "G", "T"),
    as.prob = TRUE
)

gc_pct <- 100 * (
    base_freq["G"] +
    base_freq["C"]
)

cat("Spliced transcript length:", tx_length, "nt\n")
cat("GC content:", round(gc_pct, 2), "%\n\n")

# ------------------------------------------------------------
# 7. Save RNA-equivalent DNA transcript FASTA
# ------------------------------------------------------------

tx_fasta <- DNAStringSet(transcript)

names(tx_fasta) <- paste0(
    target_gene,
    "|",
    target_tx,
    "|misc_RNA"
)

writeXStringSet(
    tx_fasta,
    filepath = file.path(
        outdir,
        "LOC118273292_XR_007705569.1_spliced.fa"
    )
)

writeXStringSet(
    exon_seqs,
    filepath = file.path(
        outdir,
        "LOC118273292_exons.fa"
    )
)

# ------------------------------------------------------------
# 8. ORF finder
#    Finds ATG -> first downstream in-frame stop
# ------------------------------------------------------------

find_orfs_one_frame <- function(dna, frame, strand_label) {

    seq_char <- as.character(dna)

    starts_pos <- seq(
        from = frame,
        to = nchar(seq_char) - 2,
        by = 3
    )

    codons <- substring(
        seq_char,
        starts_pos,
        starts_pos + 2
    )

    atg_idx <- which(codons == "ATG")
    stop_idx <- which(
        codons %in% c(
            "TAA",
            "TAG",
            "TGA"
        )
    )

    if (
        length(atg_idx) == 0 ||
        length(stop_idx) == 0
    ) {
        return(tibble())
    }

    out <- list()
    k <- 1

    for (s in atg_idx) {

        downstream_stops <- stop_idx[
            stop_idx > s
        ]

        if (length(downstream_stops) == 0) {
            next
        }

        e <- downstream_stops[1]

        nt_start <- starts_pos[s]
        nt_end <- starts_pos[e] + 2

        orf_seq <- substring(
            seq_char,
            nt_start,
            nt_end
        )

        aa <- translate(
            DNAString(orf_seq),
            if.fuzzy.codon = "X"
        )

        out[[k]] <- tibble(
            strand = strand_label,
            frame = frame,
            nt_start = nt_start,
            nt_end = nt_end,
            nt_length = nchar(orf_seq),
            aa_length = length(aa),
            aa_sequence = as.character(aa)
        )

        k <- k + 1
    }

    bind_rows(out)
}

orf_fwd <- bind_rows(
    lapply(
        1:3,
        function(f) {
            find_orfs_one_frame(
                transcript,
                f,
                "+"
            )
        }
    )
)

rc <- reverseComplement(transcript)

orf_rev <- bind_rows(
    lapply(
        1:3,
        function(f) {
            find_orfs_one_frame(
                rc,
                f,
                "-"
            )
        }
    )
)

orfs <- bind_rows(
    orf_fwd,
    orf_rev
) %>%
    arrange(
        desc(aa_length)
    )

write_csv(
    orfs,
    file.path(
        outdir,
        "LOC118273292_ORF_scan.csv"
    )
)

# ------------------------------------------------------------
# 9. Save top ORFs
# ------------------------------------------------------------

if (nrow(orfs) > 0) {

    top_orfs <- orfs %>%
        slice_head(
            n = min(
                10,
                nrow(orfs)
            )
        )

    aa_set <- AAStringSet(
        top_orfs$aa_sequence
    )

    names(aa_set) <- paste0(
        target_gene,
        "_ORF",
        seq_len(length(aa_set)),
        "|strand=",
        top_orfs$strand,
        "|frame=",
        top_orfs$frame,
        "|aa=",
        top_orfs$aa_length
    )

    writeXStringSet(
        aa_set,
        filepath = file.path(
            outdir,
            "LOC118273292_top_ORFs.faa"
        )
    )
}

# ------------------------------------------------------------
# 10. Save exon coordinates
# ------------------------------------------------------------

write_csv(
    exons %>%
        select(
            seqid,
            type,
            start,
            end,
            strand,
            attributes
        ),
    file.path(
        outdir,
        "LOC118273292_exon_coordinates.csv"
    )
)

# ------------------------------------------------------------
# 11. Report
# ------------------------------------------------------------

cat("============================================================\n")
cat("TOP ORFs\n")
cat("============================================================\n\n")

if (nrow(orfs) == 0) {

    cat("No canonical ATG-to-stop ORFs found.\n")

} else {

    print(
        orfs %>%
            select(
                strand,
                frame,
                nt_start,
                nt_end,
                nt_length,
                aa_length
            ) %>%
            slice_head(
                n = 15
            ),
        n = Inf,
        width = Inf
    )
}

cat("\nOutputs:\n")
cat("  ", outdir, "/LOC118273292_XR_007705569.1_spliced.fa\n", sep="")
cat("  ", outdir, "/LOC118273292_exons.fa\n", sep="")
cat("  ", outdir, "/LOC118273292_ORF_scan.csv\n", sep="")
cat("  ", outdir, "/LOC118273292_top_ORFs.faa\n", sep="")
cat("  ", outdir, "/LOC118273292_exon_coordinates.csv\n", sep="")

cat("\n============================================================\n")
cat("STEP 16 COMPLETE\n")
cat("============================================================\n")

