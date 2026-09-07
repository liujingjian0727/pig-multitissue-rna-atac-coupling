#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

# ============================================================
# STEP 08B1c
# FINAL AME SUMMARY
#
# MEME Suite 5.5.8
# JASPAR2024 CORE vertebrates non-redundant
#
# Primary known-motif significance:
#     AME E-value < 0.05
#
# adj_p-value is retained but is NOT treated as the
# database-wide primary significance criterion.
#
# Cerebellum:
# MA1539.1 / NR2F6 reproducibly produces no AME Fisher/max
# result and is recorded as unavailable rather than p=1.
# ============================================================


TISSUES <- c(
    "Muscle",
    "Spleen",
    "Liver",
    "Cerebellum"
)


ROOT <- "08B1b_AME_all_motifs"

OUTDIR <- "08B1c_AME_final_summary"

MOTIF_DB <- paste0(
    "/datadisk2/liujingjian_data/ATAC-seq/Sus/",
    "footprints/footprinting/motif_databases/JASPAR/",
    "JASPAR2024_CORE_vertebrates_non-redundant_v2.meme"
)


dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


cat(
    "============================================================\n"
)

cat(
    "STEP 08B1c: FINAL AME SUMMARY\n"
)

cat(
    "============================================================\n\n"
)


# ============================================================
# READ JASPAR MOTIF IDS
# ============================================================

db_lines <- readLines(
    MOTIF_DB
)


motif_lines <- grep(
    "^MOTIF ",
    db_lines,
    value = TRUE
)


parse_motif <- function(x) {

    z <- strsplit(
        x,
        "[[:space:]]+"
    )[[1]]


    data.frame(
        motif_ID = ifelse(
            length(z) >= 2,
            z[2],
            NA
        ),

        motif_alt_ID_DB = ifelse(
            length(z) >= 3,
            paste(
                z[3:length(z)],
                collapse = " "
            ),
            NA
        ),

        stringsAsFactors = FALSE
    )
}


db <- do.call(
    rbind,
    lapply(
        motif_lines,
        parse_motif
    )
)


if (nrow(db) != 879) {

    stop(
        paste(
            "Expected 879 JASPAR motifs, observed",
            nrow(db)
        )
    )
}


if (anyDuplicated(
    db$motif_ID
)) {

    stop(
        "Duplicated motif_ID found in JASPAR database."
    )
}


cat(
    "JASPAR motifs:",
    nrow(db),
    "\n\n"
)


# ============================================================
# READ AME TABLE
# ============================================================

read_ame <- function(
    tissue
) {

    file <- file.path(
        ROOT,
        tissue,
        "ame.tsv"
    )


    if (!file.exists(file)) {

        stop(
            paste(
                "Missing AME file:",
                file
            )
        )
    }


    x <- read.delim(
        file,
        comment.char = "#",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )


    required <- c(
        "rank",
        "motif_ID",
        "motif_alt_ID",
        "p-value",
        "adj_p-value",
        "E-value",
        "pos",
        "neg",
        "TP",
        "%TP",
        "FP",
        "%FP"
    )


    missing_columns <- setdiff(
        required,
        colnames(x)
    )


    if (length(
        missing_columns
    ) > 0) {

        stop(
            paste(
                tissue,
                "missing columns:",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )
    }


    numeric_columns <- c(
        "rank",
        "p-value",
        "adj_p-value",
        "E-value",
        "tests",
        "FASTA_max",
        "pos",
        "neg",
        "PWM_min",
        "TP",
        "%TP",
        "FP",
        "%FP"
    )


    for (cc in intersect(
        numeric_columns,
        colnames(x)
    )) {

        x[[cc]] <- suppressWarnings(
            as.numeric(
                x[[cc]]
            )
        )
    }


    x$Tissue <- tissue


    x
}


# ============================================================
# PROCESS EACH TISSUE
# ============================================================

all_results <- list()

missing_results <- list()

qc_results <- list()


for (tissue in TISSUES) {

    cat(
        "Processing:",
        tissue,
        "\n"
    )


    x <- read_ame(
        tissue
    )


    # --------------------------------------------------------
    # Basic motif-ID QC
    # --------------------------------------------------------

    unknown_ids <- setdiff(
        x$motif_ID,
        db$motif_ID
    )


    if (length(
        unknown_ids
    ) > 0) {

        stop(
            paste(
                tissue,
                "contains motif IDs absent from JASPAR database:",
                paste(
                    unknown_ids,
                    collapse = ","
                )
            )
        )
    }


    missing_ids <- setdiff(
        db$motif_ID,
        x$motif_ID
    )


    if (length(
        missing_ids
    ) > 0) {

        missing_table <- db[
            db$motif_ID %in% missing_ids,
            ,
            drop = FALSE
        ]


        missing_table$Tissue <- tissue

        missing_table$Reason <- (
            "No_AME_Fisher_max_result_reported"
        )


        missing_results[[tissue]] <- missing_table

    } else {

        missing_results[[tissue]] <- NULL
    }


    # --------------------------------------------------------
    # Effect sizes
    # --------------------------------------------------------

    x$Foreground_fraction <- (
        x$TP /
        x$pos
    )


    x$Background_fraction <- (
        x$FP /
        x$neg
    )


    x$Foreground_percent_calc <- (
        100 *
        x$Foreground_fraction
    )


    x$Background_percent_calc <- (
        100 *
        x$Background_fraction
    )


    x$Absolute_percent_difference <- (
        x$Foreground_percent_calc -
        x$Background_percent_calc
    )


    # Haldane-Anscombe corrected prevalence ratio
    x$Prevalence_ratio_HA <- (
        (x$TP + 0.5) /
        (x$pos + 1)
    ) /
    (
        (x$FP + 0.5) /
        (x$neg + 1)
    )


    # Haldane-Anscombe corrected odds ratio
    a <- x$TP
    b <- x$pos - x$TP
    c <- x$FP
    d <- x$neg - x$FP


    x$Odds_ratio_HA <- (
        (a + 0.5) *
        (d + 0.5)
    ) /
    (
        (b + 0.5) *
        (c + 0.5)
    )


    x$log2_Odds_ratio_HA <- log2(
        x$Odds_ratio_HA
    )


    # --------------------------------------------------------
    # AME significance categories
    # --------------------------------------------------------

    x$AME_E_lt_0.05 <- (
        x$`E-value` < 0.05
    )


    x$AME_E_lt_1 <- (
        x$`E-value` < 1
    )


    x$AME_E_lt_10 <- (
        x$`E-value` < 10
    )


    x$AME_adjP_lt_0.05 <- (
        x$`adj_p-value` < 0.05
    )


    x$AME_nominalP_lt_0.05 <- (
        x$`p-value` < 0.05
    )


    x$Primary_status <- ifelse(
        x$AME_E_lt_0.05,
        "Database_wide_significant",
        "Not_database_wide_significant"
    )


    # --------------------------------------------------------
    # Sort using native AME E-value first
    # --------------------------------------------------------

    x <- x[
        order(
            x$`E-value`,
            x$`adj_p-value`,
            x$`p-value`,
            -x$Odds_ratio_HA
        ),
        ,
        drop = FALSE
    ]


    x$Summary_rank <- seq_len(
        nrow(x)
    )


    # --------------------------------------------------------
    # Reorder important columns first
    # --------------------------------------------------------

    first_columns <- c(
        "Tissue",
        "Summary_rank",
        "rank",
        "motif_ID",
        "motif_alt_ID",
        "consensus",

        "p-value",
        "adj_p-value",
        "E-value",

        "pos",
        "neg",
        "TP",
        "%TP",
        "FP",
        "%FP",

        "Foreground_percent_calc",
        "Background_percent_calc",
        "Absolute_percent_difference",

        "Prevalence_ratio_HA",
        "Odds_ratio_HA",
        "log2_Odds_ratio_HA",

        "AME_E_lt_0.05",
        "AME_E_lt_1",
        "AME_E_lt_10",
        "AME_adjP_lt_0.05",
        "AME_nominalP_lt_0.05",
        "Primary_status"
    )


    remaining_columns <- setdiff(
        colnames(x),
        first_columns
    )


    x <- x[
        ,
        c(
            first_columns[
                first_columns %in% colnames(x)
            ],
            remaining_columns
        ),
        drop = FALSE
    ]


    # --------------------------------------------------------
    # Full result
    # --------------------------------------------------------

    write.table(
        x,
        file = file.path(
            OUTDIR,
            paste0(
                tissue,
                ".AME_all_results.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )


    # --------------------------------------------------------
    # Top20
    # --------------------------------------------------------

    write.table(
        head(
            x,
            20
        ),
        file = file.path(
            OUTDIR,
            paste0(
                tissue,
                ".AME_top20.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )


    # --------------------------------------------------------
    # Threshold tables
    # --------------------------------------------------------

    write.table(
        x[
            x$AME_E_lt_0.05,
            ,
            drop = FALSE
        ],
        file = file.path(
            OUTDIR,
            paste0(
                tissue,
                ".AME_E_lt_0.05.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )


    write.table(
        x[
            x$AME_E_lt_1,
            ,
            drop = FALSE
        ],
        file = file.path(
            OUTDIR,
            paste0(
                tissue,
                ".AME_E_lt_1.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )


    write.table(
        x[
            x$AME_E_lt_10,
            ,
            drop = FALSE
        ],
        file = file.path(
            OUTDIR,
            paste0(
                tissue,
                ".AME_E_lt_10.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )


    # --------------------------------------------------------
    # QC summary
    # --------------------------------------------------------

    best <- x[1, ]


    qc <- data.frame(

        Tissue = tissue,

        JASPAR_input_motifs =
            nrow(db),

        AME_reported_motifs =
            nrow(x),

        Missing_AME_results =
            length(missing_ids),

        Missing_motif_ID =
            ifelse(
                length(missing_ids) == 0,
                "None",
                paste(
                    missing_ids,
                    collapse = ";"
                )
            ),

        E_lt_0.05 =
            sum(
                x$AME_E_lt_0.05,
                na.rm = TRUE
            ),

        E_lt_1 =
            sum(
                x$AME_E_lt_1,
                na.rm = TRUE
            ),

        E_lt_10 =
            sum(
                x$AME_E_lt_10,
                na.rm = TRUE
            ),

        AdjP_lt_0.05 =
            sum(
                x$AME_adjP_lt_0.05,
                na.rm = TRUE
            ),

        NominalP_lt_0.05 =
            sum(
                x$AME_nominalP_lt_0.05,
                na.rm = TRUE
            ),

        Best_motif_ID =
            best$motif_ID,

        Best_motif_alt_ID =
            best$motif_alt_ID,

        Best_p_value =
            best$`p-value`,

        Best_adj_p_value =
            best$`adj_p-value`,

        Best_E_value =
            best$`E-value`,

        Best_FG_percent =
            best$Foreground_percent_calc,

        Best_BG_percent =
            best$Background_percent_calc,

        Best_absolute_percent_difference =
            best$Absolute_percent_difference,

        Best_OR_HA =
            best$Odds_ratio_HA,

        stringsAsFactors = FALSE
    )


    qc_results[[tissue]] <- qc


    all_results[[tissue]] <- x
}


# ============================================================
# COMBINE QC
# ============================================================

qc_table <- do.call(
    rbind,
    qc_results
)


write.table(
    qc_table,
    file.path(
        OUTDIR,
        "08B1c_AME_QC_and_significance_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# COMBINE ALL RESULTS
# ============================================================

combined <- do.call(
    rbind,
    all_results
)


write.table(
    combined,
    file.path(
        OUTDIR,
        "08B1c_all_tissues_AME_results.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# MISSING MOTIF RESULTS
# ============================================================

valid_missing <- Filter(
    Negate(is.null),
    missing_results
)


if (length(
    valid_missing
) > 0) {

    missing_table <- do.call(
        rbind,
        valid_missing
    )


    missing_table <- missing_table[
        ,
        c(
            "Tissue",
            "motif_ID",
            "motif_alt_ID_DB",
            "Reason"
        )
    ]


} else {

    missing_table <- data.frame(
        Tissue = character(),
        motif_ID = character(),
        motif_alt_ID_DB = character(),
        Reason = character()
    )
}


write.table(
    missing_table,
    file.path(
        OUTDIR,
        "08B1c_AME_missing_motif_results.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# FINAL LOG
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "STEP 08B1c COMPLETED\n"
)

cat(
    "============================================================\n\n"
)


print(
    qc_table,
    row.names = FALSE
)


cat(
    "\nMissing AME motif results:\n"
)


if (nrow(
    missing_table
) == 0) {

    cat(
        "None\n"
    )

} else {

    print(
        missing_table,
        row.names = FALSE
    )
}


cat(
    "\nPrimary known-motif significance criterion:\n"
)

cat(
    "AME E-value < 0.05\n"
)


cat(
    "\nIMPORTANT:\n"
)

cat(
    "Effect sizes use the motif-specific optimized AME threshold.\n"
)

cat(
    "They are descriptive complements to the native AME statistic.\n"
)

cat(
    "Missing AME results are NOT imputed as p=1.\n"
)

cat(
    "============================================================\n"
)

