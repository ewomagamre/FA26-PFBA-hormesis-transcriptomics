
# ============================================================
# FA26 RNA-seq
# Figure 4 — Functional architecture and gene-level discovery
#
# Goals:
#   1. Reconstruct temporal DEG classes
#   2. Reconstruct shared/cross-dose gene sets
#   3. Save complete membership table
#   4. Locate the most complete annotation resource available
#   5. Build gene-level evidence tables
#   6. Characterize the 22 shared 4-h genes
#   7. Screen early/persistent programs for:
#        a. canonical stress-response features
#        b. metabolic/bioenergetic/adaptive features
#
# NOTE:
# This does NOT replace unbiased GO/KEGG enrichment.
# Targeted stress/adaptation screening is an additional analysis.
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("FIGURE 4 — FUNCTIONAL ARCHITECTURE\n")
cat("============================================================\n\n")

dir.create("results/figure4", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. INPUT FILES
# ============================================================

high_file <- "results/deseq2/signatures/FA26_PFBA1_candidate_signatures.csv"
low_file  <- "results/deseq2/signatures_lowdose/FA26_PFBA001_candidate_signatures.csv"

if (!file.exists(high_file)) stop("Missing: ", high_file)
if (!file.exists(low_file))  stop("Missing: ", low_file)

high <- read.csv(high_file, check.names = FALSE)
low  <- read.csv(low_file,  check.names = FALSE)

cat("High-dose signature table:", nrow(high), "genes\n")
cat("Low-dose signature table :", nrow(low),  "genes\n\n")

# ============================================================
# 2. HELPERS
# ============================================================

sig <- function(x) !is.na(x) & x

clean_ids <- function(x) {
    unique(x[!is.na(x) & x != ""])
}

collapse_unique <- function(x) {
    x <- unique(x[!is.na(x) & x != "" & x != "-"])
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = "; ")
}

# ============================================================
# 3. BUILD TEMPORAL GENE SETS
# ============================================================

# --------------------------
# HIGH: 1 ug/g PFBA
# --------------------------

H2 <- clean_ids(high$gene_id[sig(high$sig_2h)])
H4 <- clean_ids(high$gene_id[sig(high$sig_4h)])

H_early      <- setdiff(H2, H4)
H_persistent <- intersect(H2, H4)
H_emerging   <- setdiff(H4, H2)

H_upup <- clean_ids(
    high$gene_id[
        sig(high$sig_2h) &
        sig(high$sig_4h) &
        high$direction_2h == "Up" &
        high$direction_4h == "Up"
    ]
)

H_downdown <- clean_ids(
    high$gene_id[
        sig(high$sig_2h) &
        sig(high$sig_4h) &
        high$direction_2h == "Down" &
        high$direction_4h == "Down"
    ]
)

H_updown <- clean_ids(
    high$gene_id[
        sig(high$sig_2h) &
        sig(high$sig_4h) &
        high$direction_2h == "Up" &
        high$direction_4h == "Down"
    ]
)

H_downup <- clean_ids(
    high$gene_id[
        sig(high$sig_2h) &
        sig(high$sig_4h) &
        high$direction_2h == "Down" &
        high$direction_4h == "Up"
    ]
)

H_switch <- union(H_updown, H_downup)

# --------------------------
# LOW: 0.01 ug/g PFBA
# --------------------------

L2 <- clean_ids(low$gene_id[sig(low$sig_2h)])
L4 <- clean_ids(low$gene_id[sig(low$sig_4h)])

L_early      <- setdiff(L2, L4)
L_persistent <- intersect(L2, L4)
L_emerging   <- setdiff(L4, L2)

# ============================================================
# 4. DIRECTIONAL CROSS-DOSE SETS
# ============================================================

get_dir <- function(dat, time, direction) {

    sig_col <- paste0("sig_", time)
    dir_col <- paste0("direction_", time)

    clean_ids(
        dat$gene_id[
            sig(dat[[sig_col]]) &
            dat[[dir_col]] == direction
        ]
    )
}

H2_up   <- get_dir(high, "2h", "Up")
H2_down <- get_dir(high, "2h", "Down")
H4_up   <- get_dir(high, "4h", "Up")
H4_down <- get_dir(high, "4h", "Down")

L2_up   <- get_dir(low, "2h", "Up")
L2_down <- get_dir(low, "2h", "Down")
L4_up   <- get_dir(low, "4h", "Up")
L4_down <- get_dir(low, "4h", "Down")

shared4_up   <- intersect(H4_up,   L4_up)
shared4_down <- intersect(H4_down, L4_down)

lead_up   <- intersect(H2_up,   L4_up)
lead_down <- intersect(H2_down, L4_down)

# ============================================================
# 5. MASTER SET LIST
# ============================================================

sets <- list(

    High_EarlyOnly             = H_early,

    High_Persistent_All        = H_persistent,
    High_Persistent_UpUp       = H_upup,
    High_Persistent_DownDown   = H_downdown,
    High_Persistent_UpDown     = H_updown,
    High_Persistent_DownUp     = H_downup,
    High_Persistent_Switch     = H_switch,

    High_LaterEmerging         = H_emerging,

    Low_EarlyOnly              = L_early,
    Low_Persistent             = L_persistent,
    Low_LaterEmerging          = L_emerging,

    Shared4h_Up                = shared4_up,
    Shared4h_Down              = shared4_down,

    High2h_Low4h_Up            = lead_up,
    High2h_Low4h_Down          = lead_down
)

cat("GENE SET SIZES\n")
cat("------------------------------------------------------------\n")

for (nm in names(sets)) {
    cat(sprintf("%-30s %5d\n", nm, length(sets[[nm]])))
}

# ============================================================
# 6. SAVE MEMBERSHIP TABLE
#    Correctly handles zero-length sets
# ============================================================

membership_list <- lapply(
    names(sets),
    function(nm) {

        ids <- sets[[nm]]

        if (length(ids) == 0)
            return(NULL)

        data.frame(
            gene_id  = ids,
            gene_set = rep(nm, length(ids)),
            stringsAsFactors = FALSE
        )
    }
)

membership <- do.call(rbind, membership_list)

write.csv(
    membership,
    "results/figure4/FA26_Figure4_gene_set_membership.csv",
    row.names = FALSE
)

saveRDS(
    sets,
    "results/figure4/FA26_Figure4_gene_sets.rds"
)

# ============================================================
# 7. SEARCH FOR FULL ANNOTATION RESOURCES
#
# We deliberately do NOT assume the 3,000-gene atlas is the
# correct universe for enrichment.
# ============================================================

cat("\n============================================================\n")
cat("SEARCHING FOR FULL ANNOTATION RESOURCE\n")
cat("============================================================\n\n")

search_dirs <- c(
    "results",
    "objects",
    "annotation",
    "annotations"
)

search_dirs <- search_dirs[dir.exists(search_dirs)]

candidate_files <- unique(unlist(
    lapply(
        search_dirs,
        function(d) {
            list.files(
                d,
                recursive = TRUE,
                full.names = TRUE,
                pattern = "\\.(csv|tsv|txt|rds)$",
                ignore.case = TRUE
            )
        }
    )
))

# prioritize annotation-like filenames
candidate_files <- candidate_files[
    grepl(
        "annot|eggnog|gene.*master|master.*gene|functional|gene.*info",
        basename(candidate_files),
        ignore.case = TRUE
    )
]

cat("Candidate annotation files found:", length(candidate_files), "\n\n")

# ============================================================
# 8. READER FOR CANDIDATE FILES
# ============================================================

read_candidate <- function(f) {

    tryCatch({

        if (grepl("\\.rds$", f, ignore.case = TRUE)) {

            x <- readRDS(f)

            if (!is.data.frame(x))
                x <- tryCatch(as.data.frame(x), error = function(e) NULL)

        } else if (grepl("\\.tsv$|\\.txt$", f, ignore.case = TRUE)) {

            x <- read.delim(
                f,
                check.names = FALSE,
                stringsAsFactors = FALSE,
                comment.char = "",
                quote = ""
            )

        } else {

            x <- read.csv(
                f,
                check.names = FALSE,
                stringsAsFactors = FALSE
            )
        }

        if (is.null(x) || nrow(x) == 0)
            return(NULL)

        # Normalize possible gene-ID column names
        id_candidates <- names(x)[
            tolower(names(x)) %in%
                c(
                    "gene_id",
                    "geneid",
                    "gene",
                    "query",
                    "query_name",
                    "#query"
                )
        ]

        if (length(id_candidates) == 0)
            return(NULL)

        if (!"gene_id" %in% names(x)) {
            names(x)[names(x) == id_candidates[1]] <- "gene_id"
        }

        x$gene_id <- as.character(x$gene_id)

        x

    }, error = function(e) NULL)
}

# ============================================================
# 9. SCORE ANNOTATION CANDIDATES
# ============================================================

candidate_summary <- data.frame()

annotation_objects <- list()

for (f in candidate_files) {

    x <- read_candidate(f)

    if (is.null(x))
        next

    n_ids <- length(unique(x$gene_id[
        !is.na(x$gene_id) & x$gene_id != ""
    ]))

    useful_cols <- sum(
        grepl(
            "product|preferred|description|GO|KEGG|PFAM|COG",
            names(x),
            ignore.case = TRUE
        )
    )

    candidate_summary <- rbind(
        candidate_summary,
        data.frame(
            file = f,
            rows = nrow(x),
            unique_gene_ids = n_ids,
            annotation_columns = useful_cols
        )
    )

    annotation_objects[[f]] <- x
}

if (nrow(candidate_summary) > 0) {

    candidate_summary <- candidate_summary[
        order(
            -candidate_summary$unique_gene_ids,
            -candidate_summary$annotation_columns
        ),
    ]

    write.csv(
        candidate_summary,
        "results/figure4/annotation_candidate_summary.csv",
        row.names = FALSE
    )

    cat("TOP ANNOTATION CANDIDATES\n")
    cat("------------------------------------------------------------\n")
    print(head(candidate_summary, 15), row.names = FALSE)

} else {

    cat("No annotation-like candidate containing a recognizable gene ID column was found.\n")
}

# ============================================================
# 10. SELECT BEST ANNOTATION TABLE
# ============================================================

ann <- NULL
ann_source <- NA_character_

if (nrow(candidate_summary) > 0) {

    ann_source <- candidate_summary$file[1]
    ann <- annotation_objects[[ann_source]]

    cat("\nSelected annotation resource:\n")
    cat(ann_source, "\n")
    cat("Rows:", nrow(ann), "\n")
    cat("Unique gene IDs:",
        length(unique(ann$gene_id)),
        "\n")
}

# Fall back to 3k atlas ONLY for gene-name display,
# never as the enrichment universe.

atlas_file <- "objects/FA26_gene_atlas_master_annotated_zones_clockwise.rds"

if (is.null(ann) && file.exists(atlas_file)) {

    cat("\nWARNING:\n")
    cat("No full annotation resource identified.\n")
    cat("Using 3,000-gene atlas ONLY for exploratory gene labeling.\n")
    cat("Do NOT use this as enrichment universe.\n\n")

    ann <- as.data.frame(readRDS(atlas_file))
    ann_source <- atlas_file
}

# ============================================================
# 11. BUILD STANDARDIZED DISPLAY ANNOTATION
# ============================================================

if (!is.null(ann)) {

    ann <- ann[!duplicated(ann$gene_id), ]

    # Find likely columns
    find_col <- function(patterns) {

        hit <- names(ann)[
            Reduce(
                `|`,
                lapply(
                    patterns,
                    function(p)
                        grepl(p, names(ann), ignore.case = TRUE)
                )
            )
        ]

        if (length(hit) == 0)
            return(NULL)

        hit[1]
    }

    product_col <- find_col(c("^NCBI_product$", "product"))
    name_col    <- find_col(c("^preferred_name$", "gene_name", "symbol"))
    desc_col    <- find_col(c("eggNOG_description", "description"))
    go_col      <- find_col(c("^GOs$", "^GO$", "go_terms"))
    ko_col      <- find_col(c("^KEGG_ko$", "kegg.*ko"))
    path_col    <- find_col(c("^KEGG_Pathway$", "kegg.*path"))
    pfam_col    <- find_col(c("PFAM"))
    cog_col     <- find_col(c("COG"))

    get_field <- function(col) {
        if (is.null(col))
            rep(NA_character_, nrow(ann))
        else
            as.character(ann[[col]])
    }

    ann_std <- data.frame(
        gene_id = ann$gene_id,
        gene_name = get_field(name_col),
        product = get_field(product_col),
        description = get_field(desc_col),
        GO = get_field(go_col),
        KEGG_KO = get_field(ko_col),
        KEGG_Pathway = get_field(path_col),
        PFAM = get_field(pfam_col),
        COG = get_field(cog_col),
        stringsAsFactors = FALSE
    )

    # Convert "-" to missing
    for (j in seq_along(ann_std)) {
        if (is.character(ann_std[[j]])) {
            ann_std[[j]][ann_std[[j]] == "-"] <- NA
        }
    }

    # Main-figure display label:
    # reliable gene name first; LOC ID when unnamed
    ann_std$display_name <- ann_std$gene_name

    bad_name <- is.na(ann_std$display_name) |
                ann_std$display_name == "" |
                ann_std$display_name == "-"

    ann_std$display_name[bad_name] <-
        ann_std$gene_id[bad_name]

    write.csv(
        ann_std,
        "results/figure4/FA26_Figure4_standardized_annotation.csv",
        row.names = FALSE
    )

} else {

    ann_std <- NULL
}

# ============================================================
# 12. BUILD COMBINED DE EVIDENCE TABLE
# ============================================================

high_cols <- intersect(
    c(
        "gene_id",
        "shrunk_LFC_2h",
        "padj_2h",
        "direction_2h",
        "sig_2h",
        "shrunk_LFC_4h",
        "padj_4h",
        "direction_4h",
        "sig_4h"
    ),
    names(high)
)

low_cols <- intersect(
    c(
        "gene_id",
        "shrunk_LFC_2h",
        "padj_2h",
        "direction_2h",
        "sig_2h",
        "shrunk_LFC_4h",
        "padj_4h",
        "direction_4h",
        "sig_4h"
    ),
    names(low)
)

high_de <- high[, high_cols, drop = FALSE]
low_de  <- low[,  low_cols,  drop = FALSE]

names(high_de)[names(high_de) != "gene_id"] <-
    paste0("PFBA1_", names(high_de)[names(high_de) != "gene_id"])

names(low_de)[names(low_de) != "gene_id"] <-
    paste0("PFBA001_", names(low_de)[names(low_de) != "gene_id"])

all_de <- merge(
    high_de,
    low_de,
    by = "gene_id",
    all = TRUE
)

if (!is.null(ann_std)) {
    evidence <- merge(
        ann_std,
        all_de,
        by = "gene_id",
        all.y = TRUE
    )
} else {
    evidence <- all_de
}

# ============================================================
# 13. ASSIGN TEMPORAL / SHARED CLASS LABELS
# ============================================================

evidence$High_temporal_class <- NA_character_

evidence$High_temporal_class[
    evidence$gene_id %in% H_early
] <- "Early-only"

evidence$High_temporal_class[
    evidence$gene_id %in% H_persistent
] <- "Persistent"

evidence$High_temporal_class[
    evidence$gene_id %in% H_emerging
] <- "Later-emerging"


evidence$High_persistent_direction <- NA_character_

evidence$High_persistent_direction[
    evidence$gene_id %in% H_upup
] <- "Up -> Up"

evidence$High_persistent_direction[
    evidence$gene_id %in% H_downdown
] <- "Down -> Down"

evidence$High_persistent_direction[
    evidence$gene_id %in% H_updown
] <- "Up -> Down"

evidence$High_persistent_direction[
    evidence$gene_id %in% H_downup
] <- "Down -> Up"


evidence$Low_temporal_class <- NA_character_

evidence$Low_temporal_class[
    evidence$gene_id %in% L_early
] <- "Early-only"

evidence$Low_temporal_class[
    evidence$gene_id %in% L_persistent
] <- "Persistent"

evidence$Low_temporal_class[
    evidence$gene_id %in% L_emerging
] <- "Later-emerging"


evidence$shared_4h <- NA_character_

evidence$shared_4h[
    evidence$gene_id %in% shared4_up
] <- "Shared up"

evidence$shared_4h[
    evidence$gene_id %in% shared4_down
] <- "Shared down"


evidence$high2h_low4h <- NA_character_

evidence$high2h_low4h[
    evidence$gene_id %in% lead_up
] <- "Shared up"

evidence$high2h_low4h[
    evidence$gene_id %in% lead_down
] <- "Shared down"

write.csv(
    evidence,
    "results/figure4/FA26_Figure4_master_gene_evidence.csv",
    row.names = FALSE
)

# ============================================================
# 14. SPECIFIC GENE TABLES FOR FIGURE 4
# ============================================================

write_subset <- function(ids, filename) {

    x <- evidence[evidence$gene_id %in% ids, , drop = FALSE]

    write.csv(
        x,
        file.path("results/figure4", filename),
        row.names = FALSE
    )

    x
}

tab_shared_up <- write_subset(
    shared4_up,
    "FA26_shared4h_up_14genes.csv"
)

tab_shared_down <- write_subset(
    shared4_down,
    "FA26_shared4h_down_8genes.csv"
)

tab_lead <- write_subset(
    union(lead_up, lead_down),
    "FA26_high2h_low4h_shared_genes.csv"
)

tab_early <- write_subset(
    H_early,
    "FA26_PFBA1_earlyOnly_384genes.csv"
)

tab_persist_up <- write_subset(
    H_upup,
    "FA26_PFBA1_persistent_UpUp_123genes.csv"
)

tab_persist_down <- write_subset(
    H_downdown,
    "FA26_PFBA1_persistent_DownDown_12genes.csv"
)

tab_emerging <- write_subset(
    H_emerging,
    "FA26_PFBA1_laterEmerging_1427genes.csv"
)

tab_low_emerging <- write_subset(
    L_emerging,
    "FA26_PFBA001_laterEmerging_586genes.csv"
)

# ============================================================
# 15. TARGETED BIOLOGICAL SCREEN
#
# IMPORTANT:
# This is NOT enrichment analysis.
# It asks whether genes annotated with canonical stress or
# adaptive/metabolic terminology are represented in each set.
# ============================================================

cat("\n============================================================\n")
cat("TARGETED STRESS / ADAPTIVE PROGRAM SCREEN\n")
cat("============================================================\n")

if (!is.null(ann_std)) {

    search_text <- apply(
        ann_std[
            ,
            c(
                "gene_name",
                "product",
                "description",
                "GO",
                "KEGG_Pathway",
                "PFAM",
                "COG"
            ),
            drop = FALSE
        ],
        1,
        function(z) {
            paste(
                z[!is.na(z)],
                collapse = " | "
            )
        }
    )

    search_text <- tolower(search_text)

    # --------------------------------------------------------
    # Canonical stress-response concepts
    # --------------------------------------------------------

    stress_patterns <- list(

        HeatShock_Chaperone =
            "heat shock|hsp[0-9]|chaperon|dna[jk]|groel|protein folding",

        OxidativeStress_ROS =
            "oxidative stress|reactive oxygen|superoxide|peroxidase|catalase|thioredoxin|peroxiredoxin",

        Glutathione_Redox =
            "glutathione|glutaredoxin|redox",

        Xenobiotic_Detoxification =
            "cytochrome p450|p450|glutathione s-transferase|gst|xenobiotic|detox",

        DNA_Damage_Repair =
            "dna damage|dna repair|mismatch repair|nucleotide excision|base excision|double.strand break",

        ER_UPR_Proteostasis =
            "unfolded protein|endoplasmic reticulum stress|er stress|proteostasis",

        Autophagy_Apoptosis =
            "autophagy|apoptosis|programmed cell death|caspase"
    )

    # --------------------------------------------------------
    # Metabolic/adaptive concepts
    # --------------------------------------------------------

    adaptive_patterns <- list(

        Mitochondrial_Respiration =
            "mitochond|oxidative phosphorylation|respiratory chain|cytochrome c oxidase|complex i|complex iv",

        ATP_Bioenergetics =
            "atp synthase|atp synthesis|proton gradient|energy metabolism",

        Lipid_Metabolism =
            "fatty acid|lipid metabolism|fatty acid binding|beta oxidation|acyl",

        Carbohydrate_Metabolism =
            "glycolysis|glucose|carbohydrate|tca cycle|citric acid|pyruvate",

        AminoAcid_ProteinTurnover =
            "amino acid|peptidase|protease|protein turnover",

        Transcription_Regulation =
            "transcription factor|transcription activator|transcriptional regulation",

        Signaling =
            "signal transduction|signaling|kinase|phosphatase|second messenger",

        Transport_Homeostasis =
            "transporter|transport|ion homeostasis|solute carrier",

        Growth_Development =
            "growth|development|morphogenesis|ecdys|eclosion|molting|cell proliferation"
    )

    screen_pattern <- function(pattern_list) {

        out <- list()

        for (nm in names(pattern_list)) {

            hit <- grepl(
                pattern_list[[nm]],
                search_text,
                ignore.case = TRUE
            )

            out[[nm]] <- ann_std$gene_id[hit]
        }

        out
    }

    stress_hits   <- screen_pattern(stress_patterns)
    adaptive_hits <- screen_pattern(adaptive_patterns)

    # gene sets central to current Figure 4 plan
    focal_sets <- list(
        High_EarlyOnly = H_early,
        High_Persistent_UpUp = H_upup,
        High_LaterEmerging = H_emerging,
        Low_LaterEmerging = L_emerging,
        Shared4h_Up = shared4_up,
        Shared4h_Down = shared4_down
    )

    summarize_screen <- function(hit_sets, screen_name) {

        rows <- list()

        for (gene_set_name in names(focal_sets)) {

            g <- focal_sets[[gene_set_name]]

            for (category in names(hit_sets)) {

                overlap <- intersect(
                    g,
                    hit_sets[[category]]
                )

                rows[[length(rows) + 1]] <- data.frame(
                    screen = screen_name,
                    gene_set = gene_set_name,
                    category = category,
                    set_size = length(g),
                    matching_genes = length(overlap),
                    percent_of_set =
                        ifelse(
                            length(g) == 0,
                            NA,
                            100 * length(overlap) / length(g)
                        ),
                    gene_ids = paste(overlap, collapse = ";"),
                    stringsAsFactors = FALSE
                )
            }
        }

        do.call(rbind, rows)
    }

    stress_summary <- summarize_screen(
        stress_hits,
        "Canonical stress"
    )

    adaptive_summary <- summarize_screen(
        adaptive_hits,
        "Metabolic/adaptive"
    )

    targeted_summary <- rbind(
        stress_summary,
        adaptive_summary
    )

    write.csv(
        targeted_summary,
        "results/figure4/FA26_targeted_stress_adaptive_screen.csv",
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Gene-level matched records
    # --------------------------------------------------------

    make_hit_table <- function(hit_sets, screen_name) {

        out <- list()

        for (category in names(hit_sets)) {

            ids <- hit_sets[[category]]

            if (length(ids) == 0)
                next

            x <- evidence[
                evidence$gene_id %in% ids,
                ,
                drop = FALSE
            ]

            if (nrow(x) == 0)
                next

            x$screen <- screen_name
            x$screen_category <- category

            out[[length(out) + 1]] <- x
        }

        if (length(out) == 0)
            return(NULL)

        do.call(rbind, out)
    }

    stress_gene_table <- make_hit_table(
        stress_hits,
        "Canonical stress"
    )

    adaptive_gene_table <- make_hit_table(
        adaptive_hits,
        "Metabolic/adaptive"
    )

    if (!is.null(stress_gene_table)) {

        write.csv(
            stress_gene_table,
            "results/figure4/FA26_stress_candidate_genes.csv",
            row.names = FALSE
        )
    }

    if (!is.null(adaptive_gene_table)) {

        write.csv(
            adaptive_gene_table,
            "results/figure4/FA26_adaptive_metabolic_candidate_genes.csv",
            row.names = FALSE
        )
    }

    # --------------------------------------------------------
    # Print focal summary
    # --------------------------------------------------------

    cat("\n1 ug/g EARLY-ONLY — CANONICAL STRESS SCREEN\n")
    cat("------------------------------------------------------------\n")

    print(
        stress_summary[
            stress_summary$gene_set == "High_EarlyOnly",
            c("category", "matching_genes", "percent_of_set")
        ],
        row.names = FALSE
    )

    cat("\n1 ug/g PERSISTENT UP->UP — CANONICAL STRESS SCREEN\n")
    cat("------------------------------------------------------------\n")

    print(
        stress_summary[
            stress_summary$gene_set == "High_Persistent_UpUp",
            c("category", "matching_genes", "percent_of_set")
        ],
        row.names = FALSE
    )

    cat("\n1 ug/g EARLY-ONLY — METABOLIC / ADAPTIVE SCREEN\n")
    cat("------------------------------------------------------------\n")

    print(
        adaptive_summary[
            adaptive_summary$gene_set == "High_EarlyOnly",
            c("category", "matching_genes", "percent_of_set")
        ],
        row.names = FALSE
    )

    cat("\n1 ug/g PERSISTENT UP->UP — METABOLIC / ADAPTIVE SCREEN\n")
    cat("------------------------------------------------------------\n")

    print(
        adaptive_summary[
            adaptive_summary$gene_set == "High_Persistent_UpUp",
            c("category", "matching_genes", "percent_of_set")
        ],
        row.names = FALSE
    )
}

# ============================================================
# 16. PRINT ALL SHARED 4-H GENES
# ============================================================

cat("\n============================================================\n")
cat("SHARED 4-h GENES — MAIN FIGURE CANDIDATES\n")
cat("============================================================\n\n")

shared_all <- evidence[
    evidence$gene_id %in% union(shared4_up, shared4_down),
    ,
    drop = FALSE
]

if ("display_name" %in% names(shared_all)) {

    show_cols <- intersect(
        c(
            "gene_id",
            "display_name",
            "product",
            "description",
            "shared_4h",
            "PFBA001_shrunk_LFC_4h",
            "PFBA001_padj_4h",
            "PFBA1_shrunk_LFC_4h",
            "PFBA1_padj_4h"
        ),
        names(shared_all)
    )

} else {

    show_cols <- intersect(
        c(
            "gene_id",
            "shared_4h",
            "PFBA001_shrunk_LFC_4h",
            "PFBA001_padj_4h",
            "PFBA1_shrunk_LFC_4h",
            "PFBA1_padj_4h"
        ),
        names(shared_all)
    )
}

print(
    shared_all[, show_cols, drop = FALSE],
    row.names = FALSE
)

write.csv(
    shared_all,
    "results/figure4/FA26_shared4h_all_22genes.csv",
    row.names = FALSE
)

# ============================================================
# 17. PRINT LEAD/LAG GENES
# ============================================================

cat("\n============================================================\n")
cat("HIGH 2 h -> LOW 4 h SHARED GENES\n")
cat("============================================================\n\n")

lead_all <- evidence[
    evidence$gene_id %in% union(lead_up, lead_down),
    ,
    drop = FALSE
]

print(
    lead_all[
        ,
        intersect(
            c(
                "gene_id",
                "display_name",
                "product",
                "description",
                "PFBA1_shrunk_LFC_2h",
                "PFBA1_padj_2h",
                "PFBA001_shrunk_LFC_4h",
                "PFBA001_padj_4h"
            ),
            names(lead_all)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)

# ============================================================
# 18. SAVE ANALYSIS OBJECT
# ============================================================

saveRDS(
    list(
        sets = sets,
        evidence = evidence,
        annotation_source = ann_source
    ),
    "results/figure4/FA26_Figure4_analysis_foundation.rds"
)

cat("\n============================================================\n")
cat("FIGURE 4 FOUNDATION COMPLETE\n")
cat("============================================================\n\n")

cat("Key outputs:\n")
cat("  results/figure4/FA26_Figure4_gene_sets.rds\n")
cat("  results/figure4/FA26_Figure4_gene_set_membership.csv\n")
cat("  results/figure4/FA26_Figure4_master_gene_evidence.csv\n")
cat("  results/figure4/FA26_shared4h_all_22genes.csv\n")
cat("  results/figure4/FA26_PFBA1_earlyOnly_384genes.csv\n")
cat("  results/figure4/FA26_PFBA1_persistent_UpUp_123genes.csv\n")
cat("  results/figure4/FA26_PFBA1_laterEmerging_1427genes.csv\n")
cat("  results/figure4/FA26_PFBA001_laterEmerging_586genes.csv\n")
cat("  results/figure4/FA26_targeted_stress_adaptive_screen.csv\n")
cat("  results/figure4/annotation_candidate_summary.csv\n")
cat("  results/figure4/FA26_Figure4_analysis_foundation.rds\n\n")

cat("Next: unbiased GO/KEGG enrichment + Figure 4 panel construction.\n\n")

