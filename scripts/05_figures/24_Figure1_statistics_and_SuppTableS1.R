options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("FIGURE 1 — PHENOTYPE STATISTICS + SUPPLEMENTARY TABLE S1\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

infile <- "results/phenotype/FA26_weight_change.csv"

outdir <- "results/phenotype/supplementary_S1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

xlsx_file <- file.path(
    outdir,
    "Supplementary_Table_S1_Figure1_Phenotype.xlsx"
)

if (!file.exists(infile)) {
    stop("Missing input file: ", infile)
}

# ============================================================
# 2. LOAD VERIFIED RAW DATA
# ============================================================

dat <- read.csv(
    infile,
    check.names = FALSE
)

required_cols <- c(
    "Treatment",
    "Time_hr",
    "Replicate",
    "Weight_change"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
    stop(
        "Missing required column(s): ",
        paste(missing_cols, collapse = ", ")
    )
}

# Expected treatment order
treatment_levels <- c(
    "Control",
    "PFBA 0.01",
    "PFBA 1"
)

dat$Treatment <- factor(
    dat$Treatment,
    levels = treatment_levels
)

# ============================================================
# 3. TIME-MATCHED CONTROL NORMALIZATION
# ============================================================

control_means <- aggregate(
    Weight_change ~ Time_hr,
    data = dat[dat$Treatment == "Control", ],
    FUN = mean
)

names(control_means)[2] <- "Time_matched_control_mean"

dat <- merge(
    dat,
    control_means,
    by = "Time_hr",
    all.x = TRUE,
    sort = FALSE
)

dat$Percent_of_time_matched_control <-
    100 *
    dat$Weight_change /
    dat$Time_matched_control_mean

dat$Treatment <- factor(
    dat$Treatment,
    levels = treatment_levels
)

dat <- dat[
    order(
        dat$Time_hr,
        dat$Treatment,
        dat$Replicate
    ),
]

# Clean treatment labels for manuscript tables
pretty_treatment <- function(x) {
    x <- as.character(x)

    x[x == "PFBA 0.01"] <- "PFBA 0.01 ug/g"
    x[x == "PFBA 1"]    <- "PFBA 1 ug/g"

    x
}

# ============================================================
# 4. S1A — INDIVIDUAL PHENOTYPE DATA
# ============================================================

S1A <- data.frame(
    Treatment =
        pretty_treatment(dat$Treatment),

    Time_h =
        dat$Time_hr,

    Replicate =
        dat$Replicate,

    Weight_change_raw =
        dat$Weight_change,

    Time_matched_control_mean =
        dat$Time_matched_control_mean,

    Percent_of_time_matched_control =
        dat$Percent_of_time_matched_control
)

# ============================================================
# 5. S1B — GROUP SUMMARY + FIGURE 1 DATA
# ============================================================

groups <- split(
    dat,
    interaction(
        dat$Time_hr,
        dat$Treatment,
        drop = TRUE
    )
)

summary_list <- lapply(
    groups,
    function(z) {

        n <- nrow(z)

        raw_mean <- mean(z$Weight_change)
        raw_sd   <- sd(z$Weight_change)
        raw_sem  <- raw_sd / sqrt(n)

        pct_mean <-
            mean(
                z$Percent_of_time_matched_control
            )

        pct_sd <-
            sd(
                z$Percent_of_time_matched_control
            )

        pct_sem <-
            pct_sd / sqrt(n)

        data.frame(
            Time_h = z$Time_hr[1],
            Treatment =
                pretty_treatment(
                    z$Treatment[1]
                ),
            n = n,

            Raw_mean = raw_mean,
            Raw_SD = raw_sd,
            Raw_SEM = raw_sem,

            Time_matched_control_mean =
                z$Time_matched_control_mean[1],

            Percent_of_time_matched_control =
                pct_mean,

            Percent_SD =
                pct_sd,

            Percent_SEM =
                pct_sem
        )
    }
)

S1B <- do.call(
    rbind,
    summary_list
)

rownames(S1B) <- NULL

S1B$Treatment_order <- match(
    S1B$Treatment,
    c(
        "Control",
        "PFBA 0.01 ug/g",
        "PFBA 1 ug/g"
    )
)

S1B <- S1B[
    order(
        S1B$Time_h,
        S1B$Treatment_order
    ),
]

S1B$Treatment_order <- NULL

# Mark how values appear in Figure 1
S1B$Figure_1_display <- ifelse(
    S1B$Treatment == "Control",
    "100% reference line",
    "Plotted bar"
)

# ============================================================
# 6. S1C — TWO-WAY ANOVA ON RAW DATA
# ============================================================

anova_dat <- dat

anova_dat$Treatment <- factor(
    anova_dat$Treatment,
    levels = treatment_levels
)

anova_dat$Time_factor <- factor(
    anova_dat$Time_hr
)

fit <- aov(
    Weight_change ~
        Treatment *
        Time_factor,
    data = anova_dat
)

aov_tab <- summary(fit)[[1]]

S1C <- data.frame(
    Term = c(
        "Treatment",
        "Time",
        "Treatment x Time",
        "Residuals"
    ),

    Df =
        aov_tab[, "Df"],

    Sum_Sq =
        aov_tab[, "Sum Sq"],

    Mean_Sq =
        aov_tab[, "Mean Sq"],

    F_value =
        c(
            aov_tab[1:3, "F value"],
            NA
        ),

    P_value =
        c(
            aov_tab[1:3, "Pr(>F)"],
            NA
        )
)

S1C$Significance <- ifelse(
    is.na(S1C$P_value),
    "",
    ifelse(
        S1C$P_value < 0.001,
        "***",
        ifelse(
            S1C$P_value < 0.01,
            "**",
            ifelse(
                S1C$P_value < 0.05,
                "*",
                "ns"
            )
        )
    )
)

# ============================================================
# 7. S1D — PLANNED PAIRWISE WELCH TESTS
# ============================================================

times <- sort(
    unique(dat$Time_hr)
)

doses <- c(
    "PFBA 0.01",
    "PFBA 1"
)

pairwise_list <- list()
k <- 1

for (tt in times) {

    ctrl <-
        dat$Weight_change[
            dat$Time_hr == tt &
            dat$Treatment == "Control"
        ]

    tmp <- list()

    for (dose in doses) {

        trt <-
            dat$Weight_change[
                dat$Time_hr == tt &
                dat$Treatment == dose
            ]

        test <-
            t.test(
                trt,
                ctrl,
                alternative = "two.sided",
                var.equal = FALSE
            )

        tmp[[dose]] <- data.frame(
            Time_h = tt,

            Comparison =
                paste0(
                    pretty_treatment(dose),
                    " vs Control"
                ),

            Test =
                "Two-sided Welch two-sample t-test",

            Raw_P =
                test$p.value
        )
    }

    tmp_df <- do.call(
        rbind,
        tmp
    )

    # Definitive manuscript procedure:
    # Holm correction across the two planned
    # PFBA-vs-Control comparisons within each time point.
    tmp_df$Holm_adjusted_P <-
        p.adjust(
            tmp_df$Raw_P,
            method = "holm"
        )

    tmp_df$Significance <-
        ifelse(
            tmp_df$Holm_adjusted_P < 0.001,
            "***",
            ifelse(
                tmp_df$Holm_adjusted_P < 0.01,
                "**",
                ifelse(
                    tmp_df$Holm_adjusted_P < 0.05,
                    "*",
                    "ns"
                )
            )
        )

    pairwise_list[[k]] <- tmp_df
    k <- k + 1
}

S1D <- do.call(
    rbind,
    pairwise_list
)

rownames(S1D) <- NULL

# ============================================================
# 8. WRITE FOUR AUDIT CSVs
# ============================================================

write.csv(
    S1A,
    file.path(
        outdir,
        "Supplementary_Table_S1A_individual_phenotype_data.csv"
    ),
    row.names = FALSE
)

write.csv(
    S1B,
    file.path(
        outdir,
        "Supplementary_Table_S1B_Figure1_summary_data.csv"
    ),
    row.names = FALSE
)

write.csv(
    S1C,
    file.path(
        outdir,
        "Supplementary_Table_S1C_two_way_ANOVA.csv"
    ),
    row.names = FALSE
)

write.csv(
    S1D,
    file.path(
        outdir,
        "Supplementary_Table_S1D_pairwise_Welch_Holm.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 9. CREATE EXCEL WORKBOOK
# ============================================================

if (!requireNamespace(
        "openxlsx",
        quietly = TRUE
    )) {

    cat(
        "\nWARNING:\n",
        "Package 'openxlsx' is not installed.\n",
        "The four CSV files were created successfully,\n",
        "but the combined XLSX workbook could not be generated.\n"
    )

} else {

    library(openxlsx)

    wb <- createWorkbook()

    # --------------------------------------------------------
    # Styles
    # --------------------------------------------------------

    title_style <- createStyle(
        fontSize = 14,
        fontColour = "#FFFFFF",
        fgFill = "#1F4E78",
        textDecoration = "bold",
        halign = "left",
        valign = "center"
    )

    header_style <- createStyle(
        fontSize = 11,
        textDecoration = "bold",
        fgFill = "#D9EAF7",
        border = "Bottom",
        halign = "center",
        valign = "center",
        wrapText = TRUE
    )

    note_style <- createStyle(
        fontSize = 10,
        textDecoration = "italic",
        fgFill = "#FFF2CC",
        wrapText = TRUE,
        valign = "top"
    )

    # --------------------------------------------------------
    # README
    # --------------------------------------------------------

    addWorksheet(
        wb,
        "README"
    )

    writeData(
        wb,
        "README",
        "Supplementary Table S1 — Early weight-change phenotype and statistical analyses",
        startRow = 1,
        startCol = 1
    )

    mergeCells(
        wb,
        "README",
        cols = 1:2,
        rows = 1
    )

    addStyle(
        wb,
        "README",
        title_style,
        rows = 1,
        cols = 1:2,
        gridExpand = TRUE
    )

    readme <- data.frame(
        Item = c(
            "Figure supported",
            "Experimental design",
            "Primary outcome",
            "Figure normalization",
            "Error bars",
            "Inferential statistics",
            "Overall model",
            "Planned comparisons",
            "Multiplicity correction",
            "Interpretive note",
            "S1A",
            "S1B",
            "S1C",
            "S1D"
        ),

        Description = c(
            "Main Figure 1.",

            paste0(
                "Control, 0.01 ug/g PFBA, and 1 ug/g PFBA; ",
                "2, 4, and 6 h; n = 4 per treatment x time group."
            ),

            paste0(
                "Individual raw weight-change measurements. ",
                "Original measurement unit should be inserted ",
                "once confirmed from the phenotype record."
            ),

            paste0(
                "For Figure 1, PFBA group means were expressed ",
                "relative to the mean of the time-matched control: ",
                "100 x treatment mean / time-matched control mean. ",
                "The dashed reference line represents 100%."
            ),

            "SEM.",

            paste0(
                "Inferential tests were conducted on the ",
                "untransformed individual weight-change measurements, ",
                "not on the percent-control visualization."
            ),

            paste0(
                "Two-way ANOVA with Treatment, Time, and ",
                "Treatment x Time."
            ),

            paste0(
                "At each time point, each PFBA dose was compared ",
                "with control using a two-sided Welch two-sample t-test."
            ),

            paste0(
                "Holm adjustment was applied within each time point ",
                "across the two planned PFBA-vs-Control comparisons."
            ),

            paste0(
                "Percent-of-control values become large at later ",
                "time points because the corresponding raw control ",
                "means are small. Raw observations are therefore ",
                "provided in S1A for full transparency."
            ),

            "Individual phenotype observations and Figure 1 normalization.",

            "Group means, SD, SEM, and percent-of-time-matched-control values.",

            "Two-way ANOVA performed on raw weight-change observations.",

            "Planned Welch comparisons with raw and Holm-adjusted P values."
        )
    )

    writeData(
        wb,
        "README",
        readme,
        startRow = 3,
        startCol = 1,
        headerStyle = header_style
    )

    setColWidths(
        wb,
        "README",
        cols = 1,
        widths = 24
    )

    setColWidths(
        wb,
        "README",
        cols = 2,
        widths = 85
    )

    setRowHeights(
        wb,
        "README",
        rows = 3:(nrow(readme) + 3),
        heights = "auto"
    )

    freezePane(
        wb,
        "README",
        firstActiveRow = 4
    )

    # --------------------------------------------------------
    # Helper to add data sheet
    # --------------------------------------------------------

    add_table_sheet <- function(
        wb,
        sheet,
        title,
        data
    ) {

        addWorksheet(
            wb,
            sheet
        )

        writeData(
            wb,
            sheet,
            title,
            startRow = 1,
            startCol = 1
        )

        mergeCells(
            wb,
            sheet,
            cols = 1:ncol(data),
            rows = 1
        )

        addStyle(
            wb,
            sheet,
            title_style,
            rows = 1,
            cols = 1:ncol(data),
            gridExpand = TRUE
        )

        writeData(
            wb,
            sheet,
            data,
            startRow = 3,
            startCol = 1,
            headerStyle = header_style
        )

        freezePane(
            wb,
            sheet,
            firstActiveRow = 4
        )

        setColWidths(
            wb,
            sheet,
            cols = 1:ncol(data),
            widths = "auto"
        )
    }

    add_table_sheet(
        wb,
        "S1A_Raw_Data",
        "S1A — Individual phenotype observations",
        S1A
    )

    add_table_sheet(
        wb,
        "S1B_Summary",
        "S1B — Figure 1 group summary and normalization",
        S1B
    )

    add_table_sheet(
        wb,
        "S1C_ANOVA",
        "S1C — Two-way ANOVA of raw weight-change measurements",
        S1C
    )

    add_table_sheet(
        wb,
        "S1D_Pairwise",
        "S1D — Planned PFBA-vs-Control Welch tests with Holm adjustment",
        S1D
    )

    saveWorkbook(
        wb,
        xlsx_file,
        overwrite = TRUE
    )
}

# ============================================================
# 10. AUDIT OUTPUT TO TERMINAL
# ============================================================

cat("\n============================================================\n")
cat("S1A — INDIVIDUAL DATA\n")
cat("============================================================\n")
print(S1A, row.names = FALSE)

cat("\n============================================================\n")
cat("S1B — SUMMARY / FIGURE 1 DATA\n")
cat("============================================================\n")
print(S1B, row.names = FALSE)

cat("\n============================================================\n")
cat("S1C — TWO-WAY ANOVA\n")
cat("============================================================\n")
print(S1C, row.names = FALSE)

cat("\n============================================================\n")
cat("S1D — PLANNED WELCH + HOLM COMPARISONS\n")
cat("============================================================\n")
print(S1D, digits = 10, row.names = FALSE)

cat("\n============================================================\n")
cat("KEY STATISTICAL RESULTS\n")
cat("============================================================\n")

cat(
    "Treatment effect: F(",
    S1C$Df[S1C$Term == "Treatment"],
    ",",
    S1C$Df[S1C$Term == "Residuals"],
    ") = ",
    round(
        S1C$F_value[
            S1C$Term == "Treatment"
        ],
        3
    ),
    ", P = ",
    signif(
        S1C$P_value[
            S1C$Term == "Treatment"
        ],
        4
    ),
    "\n",
    sep = ""
)

cat(
    "Time effect: P = ",
    signif(
        S1C$P_value[
            S1C$Term == "Time"
        ],
        4
    ),
    "\n",
    sep = ""
)

cat(
    "Treatment x Time: P = ",
    signif(
        S1C$P_value[
            S1C$Term == "Treatment x Time"
        ],
        4
    ),
    "\n",
    sep = ""
)

cat("\nFiles created in:\n")
cat(outdir, "\n\n")

print(
    list.files(
        outdir,
        full.names = TRUE
    )
)

cat("\nDONE.\n")
