#!/usr/bin/env Rscript

###############################################################################
# STEP12B
# Mammalian Genome manuscript numeric consistency audit
#
# PURPOSE
#   1. Recompute key frozen values from Step03-Step09 source tables.
#   2. Search frozen source files for promoter/permutation/motif statistics.
#   3. Check that manuscript numbers occur in the correct biological context.
#   4. Report PASS / WARNING / MISSING without modifying any biological result.
#
# USAGE
#
#   Rscript 12B_Mammalian_Genome_numeric_consistency_audit.R \
#       Mammalian_Genome_submission_manuscript.txt
#
# or, if pandoc is installed:
#
#   Rscript 12B_Mammalian_Genome_numeric_consistency_audit.R \
#       Mammalian_Genome_submission_manuscript.docx
#
###############################################################################

options(
    stringsAsFactors = FALSE,
    scipen = 999,
    width = 200
)

ARGS <- commandArgs(trailingOnly = TRUE)

if (length(ARGS) < 1) {
    stop(
        "\nUsage:\n",
        "Rscript 12B_Mammalian_Genome_numeric_consistency_audit.R ",
        "Mammalian_Genome_submission_manuscript.txt\n",
        call. = FALSE
    )
}

MANUSCRIPT_FILE <- ARGS[1]

if (!file.exists(MANUSCRIPT_FILE)) {
    stop(
        "Manuscript file does not exist: ",
        MANUSCRIPT_FILE,
        call. = FALSE
    )
}

ROOT <- normalizePath(".", mustWork = TRUE)

OUTDIR <- file.path(
    ROOT,
    "12B_Mammalian_Genome_numeric_consistency_audit"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

LOGFILE <- file.path(
    OUTDIR,
    "Step12B_numeric_consistency_audit.log"
)

zz <- file(LOGFILE, open = "wt")

sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")


###############################################################################
# GENERAL HELPERS
###############################################################################

cat(
    "====================================================================\n",
    "STEP12B FINAL v4 — MAMMALIAN GENOME NUMERIC CONSISTENCY AUDIT\n",
    "====================================================================\n",
    "Root       : ", ROOT, "\n",
    "Manuscript : ", normalizePath(MANUSCRIPT_FILE), "\n",
    "Start time : ", format(Sys.time()), "\n",
    "====================================================================\n\n",
    sep = ""
)


die <- function(...) {
    stop(paste0(...), call. = FALSE)
}


normalize_text <- function(x) {

    x <- enc2utf8(x)

    # Normalize Unicode minus/dashes
    x <- gsub("\u2212", "-", x, fixed = TRUE)
    x <- gsub("\u2013", "-", x, fixed = TRUE)
    x <- gsub("\u2014", "-", x, fixed = TRUE)

    # Normalize multiplication symbol
    x <- gsub("\u00d7", "x", x, fixed = TRUE)

    # Normalize non-breaking space
    x <- gsub("\u00a0", " ", x, fixed = TRUE)

    # Normalize curly apostrophes
    x <- gsub("\u2018", "'", x, fixed = TRUE)
    x <- gsub("\u2019", "'", x, fixed = TRUE)

    x
}


read_manuscript <- function(file) {

    ext <- tolower(
        tools::file_ext(file)
    )

    if (ext %in% c("txt", "md", "markdown")) {

        x <- readLines(
            file,
            warn = FALSE,
            encoding = "UTF-8"
        )

    } else if (ext == "docx") {

        pandoc <- Sys.which("pandoc")

        if (pandoc == "") {
            die(
                "\nDOCX manuscript supplied but pandoc is not installed.\n",
                "Either install pandoc or save the manuscript as UTF-8 TXT.\n"
            )
        }

        x <- system2(
            pandoc,
            args = c(
                file,
                "-t",
                "plain"
            ),
            stdout = TRUE,
            stderr = FALSE
        )

    } else {

        die(
            "Unsupported manuscript format: .",
            ext,
            "\nUse TXT, MD, or DOCX."
        )
    }

    normalize_text(x)
}


make_paragraphs <- function(lines) {

    blank <- trimws(lines) == ""

    grp <- cumsum(
        c(
            TRUE,
            diff(blank) != 0
        )
    )

    result <- list()

    counter <- 0

    for (g in unique(grp)) {

        idx <- which(grp == g)

        if (all(blank[idx])) {
            next
        }

        counter <- counter + 1

        result[[counter]] <- data.frame(
            Paragraph = counter,
            Start_line = min(idx),
            End_line = max(idx),
            Text = paste(
                trimws(lines[idx]),
                collapse = " "
            ),
            stringsAsFactors = FALSE
        )
    }

    do.call(
        rbind,
        result
    )
}


extract_numbers <- function(x) {

    pattern <- paste0(
        "(?<![A-Za-z])",
        "[-+]?",
        "[0-9][0-9,]*",
        "(?:\\.[0-9]+)?",
        "(?:[eE][-+]?[0-9]+)?",
        "%?"
    )

    z <- regmatches(
        x,
        gregexpr(
            pattern,
            x,
            perl = TRUE
        )
    )[[1]]

    if (
        length(z) == 1 &&
        identical(z, character(0))
    ) {
        return("")
    }

    paste(
        unique(z),
        collapse = "; "
    )
}


safe_read_table <- function(file) {

    if (!file.exists(file)) {
        return(NULL)
    }

    tryCatch(
        read.delim(
            file,
            header = TRUE,
            check.names = FALSE,
            stringsAsFactors = FALSE
        ),
        error = function(e) {
            cat(
                "WARNING: failed reading ",
                file,
                "\n",
                sep = ""
            )

            NULL
        }
    )
}


bool_value <- function(x) {

    if (is.logical(x)) {
        return(x)
    }

    toupper(as.character(x)) %in%
        c(
            "TRUE",
            "T",
            "1",
            "YES",
            "Y"
        )
}


first_existing <- function(x) {

    x <- x[file.exists(x)]

    if (length(x) == 0) {
        return(NA_character_)
    }

    normalizePath(
        x[1],
        mustWork = TRUE
    )
}


detect_column <- function(
    df,
    patterns,
    numeric_only = FALSE
) {

    nm <- colnames(df)

    hits <- unique(
        unlist(
            lapply(
                patterns,
                function(p) {
                    grep(
                        p,
                        nm,
                        ignore.case = TRUE,
                        value = TRUE,
                        perl = TRUE
                    )
                }
            )
        )
    )

    if (numeric_only) {

        hits <- hits[
            sapply(
                df[hits],
                function(x) {
                    is.numeric(x) ||
                        all(
                            suppressWarnings(
                                !is.na(
                                    as.numeric(
                                        as.character(x)
                                    )
                                )
                            ) |
                            is.na(x)
                        )
                }
            )
        ]
    }

    if (length(hits) == 0) {
        return(NA_character_)
    }

    hits[1]
}


fmt_num <- function(x, digits = 6) {

    if (is.na(x)) {
        return("NA")
    }

    format(
        x,
        digits = digits,
        scientific = FALSE,
        trim = TRUE
    )
}


###############################################################################
# READ MANUSCRIPT
###############################################################################

cat(
    "====================================================================\n",
    "1. READING MANUSCRIPT\n",
    "====================================================================\n",
    sep = ""
)

MANUSCRIPT_LINES <- read_manuscript(
    MANUSCRIPT_FILE
)

PARAGRAPHS <- make_paragraphs(
    MANUSCRIPT_LINES
)

cat(
    "Manuscript lines      : ",
    length(MANUSCRIPT_LINES),
    "\n",
    "Manuscript paragraphs : ",
    nrow(PARAGRAPHS),
    "\n\n",
    sep = ""
)

write.table(
    data.frame(
        Line = seq_along(MANUSCRIPT_LINES),
        Text = MANUSCRIPT_LINES
    ),
    file.path(
        OUTDIR,
        "00_manuscript_numbered_lines.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


###############################################################################
# CANONICAL FROZEN RESULTS
###############################################################################

cat(
    "====================================================================\n",
    "2. BUILDING CANONICAL FROZEN METRIC TABLE\n",
    "====================================================================\n",
    sep = ""
)


METRICS <- data.frame(
    ID = character(),
    Section = character(),
    Metric = character(),
    Expected_display = character(),
    Expected_regex = character(),
    Context_regex = character(),
    Importance = character(),
    stringsAsFactors = FALSE
)


add_metric <- function(
    id,
    section,
    metric,
    expected,
    regex,
    context,
    importance = "PRIMARY"
) {

    METRICS <<- rbind(
        METRICS,
        data.frame(
            ID = id,
            Section = section,
            Metric = metric,
            Expected_display = expected,
            Expected_regex = regex,
            Context_regex = context,
            Importance = importance,
            stringsAsFactors = FALSE
        )
    )
}


###############################################################################
# CORE RESULTS
###############################################################################

add_metric(
    "RNA_LRT",
    "Core",
    "RNA genes with tissue effect",
    "18,870",
    "\\b18,?870\\b",
    "(RNA|gene).*(18,?870|tissue|FDR)|18,?870.*(RNA|gene)"
)

add_metric(
    "ATAC_LRT",
    "Core",
    "ATAC peaks with tissue effect",
    "176,917",
    "\\b176,?917\\b",
    "(ATAC|peak|accessible).*(176,?917|tissue|FDR)|176,?917.*(ATAC|peak)"
)

add_metric(
    "RNA_HC",
    "Core",
    "High-confidence tissue-specific RNA genes",
    "5,036",
    "\\b5,?036\\b",
    "(5,?036.*(gene|RNA|high-confidence))|((gene|RNA|high-confidence).*5,?036)"
)

add_metric(
    "ATAC_HC",
    "Core",
    "High-confidence tissue-specific ATAC peaks",
    "40,545",
    "\\b40,?545\\b",
    "(40,?545.*(ATAC|peak|accessible|high-confidence))|((ATAC|peak|accessible|high-confidence).*40,?545)"
)

add_metric(
    "CONCORDANT_PAIRS",
    "Core",
    "Concordant HC peak-gene associations",
    "4,674",
    "\\b4,?674\\b",
    "(4,?674.*(concordant|association|peak-gene))|((concordant|association|peak-gene).*4,?674)"
)

add_metric(
    "CONCORDANT_GENES",
    "Core",
    "Unique genes in 4,674 concordant associations",
    "1,417",
    "\\b1,?417\\b",
    "(1,?417.*gene)|(gene.*1,?417)"
)

add_metric(
    "POSITIVE_ALL_THREE",
    "Core",
    "Positive in P348/P350/mean",
    "3,916",
    "\\b3,?916\\b",
    "(3,?916.*(positive|correlation|coupl))|((positive|correlation|coupl).*3,?916)"
)

add_metric(
    "STRONG_PAIRS",
    "Core",
    "Strongly coupled pairs",
    "702",
    "\\b702\\b",
    "(702.*(strong|coupled|association|candidate))|((strong|coupled|association|candidate).*702)"
)

add_metric(
    "STRONG_GENES",
    "Core",
    "Unique genes in strong set",
    "450",
    "\\b450\\b",
    "(450.*gene)|(gene.*450)"
)

add_metric(
    "VERY_STRONG_PAIRS",
    "Core",
    "Very strongly coupled pairs",
    "226",
    "\\b226\\b",
    "(226.*very strong)|(very strong.*226)"
)

add_metric(
    "VERY_STRONG_GENES",
    "Core",
    "Unique genes in very strong set",
    "176",
    "\\b176\\b",
    "(176.*gene)|(gene.*176)"
)


###############################################################################
# PROMOTER / ARCHITECTURE
###############################################################################

add_metric(
    "ALL_ASSIGNED",
    "Promoter",
    "All UROPA-assigned pairs",
    "168,999",
    "\\b168,?999\\b",
    "(168,?999.*(UROPA|assigned|association|peak-gene))|((UROPA|assigned|association|peak-gene).*168,?999)"
)

add_metric(
    "ALL_PROMOTER",
    "Promoter",
    "Promoter pairs in all assigned set",
    "16,677",
    "\\b16,?677\\b",
    "(16,?677.*promoter)|(promoter.*16,?677)"
)

add_metric(
    "STRONG_PROMOTER",
    "Promoter",
    "Promoter pairs in 702 strong set",
    "121",
    "\\b121\\b",
    "(121.*promoter)|(promoter.*121)"
)

add_metric(
    "STRONG_PROMOTER_PCT",
    "Promoter",
    "Promoter percentage among strong pairs",
    "17.24%",
    "17\\.2(?:3[0-9]*|4)%?",
    "(17\\.2(?:3[0-9]*|4).*promoter)|(promoter.*17\\.2(?:3[0-9]*|4))"
)

add_metric(
    "PROMOTER_OR_ALL",
    "Promoter",
    "Strong vs all-non-strong promoter OR",
    "1.909",
    "\\b1\\.90(?:8[0-9]*|9)\\b",
    "(1\\.90(?:8[0-9]*|9).*(OR|odds|promoter))|((OR|odds|promoter).*1\\.90(?:8[0-9]*|9))"
)

add_metric(
    "PROMOTER_FDR_ALL",
    "Promoter",
    "Strong vs all-non-strong promoter FDR",
    "6.04e-9",
    "(6\\.04[eE]-?0?9)|(6\\.04\\s*x\\s*10\\^?-?9)",
    "(6\\.04.*(FDR|promoter))|((FDR|promoter).*6\\.04)"
)

add_metric(
    "PROMOTER_OR_HC",
    "Promoter",
    "Strong vs HC-non-strong promoter OR",
    "1.752",
    "\\b1\\.75(?:1[0-9]*|2)\\b",
    "(1\\.75(?:1[0-9]*|2).*(OR|odds|promoter))|((OR|odds|promoter).*1\\.75(?:1[0-9]*|2))"
)

add_metric(
    "PROMOTER_FDR_HC",
    "Promoter",
    "Strong vs HC-non-strong promoter FDR",
    "1.03e-5",
    "(1\\.03[eE]-?0?5)|(1\\.03\\s*x\\s*10\\^?-?5)",
    "(1\\.03.*(FDR|promoter))|((FDR|promoter).*1\\.03)"
)

add_metric(
    "TSS_ALL_MEDIAN",
    "Promoter",
    "All-pair median absolute TSS distance",
    "33,378 bp",
    "\\b33,?378\\b",
    "(33,?378.*TSS)|(TSS.*33,?378)"
)

add_metric(
    "TSS_HC_MEDIAN",
    "Promoter",
    "HC median absolute TSS distance",
    "17,248.5 bp",
    "17,?248\\.5",
    "(17,?248\\.5.*TSS)|(TSS.*17,?248\\.5)"
)

add_metric(
    "TSS_STRONG_MEDIAN",
    "Promoter",
    "Strong median absolute TSS distance",
    "11,549 bp",
    "\\b11,?549\\b",
    "(11,?549.*TSS)|(TSS.*11,?549)"
)


###############################################################################
# PERMUTATION
###############################################################################

add_metric(
    "PERM_STRONG_NULL",
    "Permutation",
    "Strong null mean",
    "588.39",
    "\\b588\\.39\\b",
    "(588\\.39.*(null|permutation|expect))|((null|permutation|expect).*588\\.39)"
)

add_metric(
    "PERM_STRONG_FOLD",
    "Permutation",
    "Strong fold enrichment",
    "1.193",
    "\\b1\\.193\\b",
    "(1\\.193.*(fold|enrichment))|((fold|enrichment).*1\\.193)"
)

add_metric(
    "PERM_STRONG_Z",
    "Permutation",
    "Strong permutation Z",
    "5.69",
    "\\b5\\.69\\b",
    "(5\\.69.*(Z|permutation))|((Z|permutation).*5\\.69)"
)

add_metric(
    "PERM_VERY_NULL",
    "Permutation",
    "Very-strong null mean",
    "167.93",
    "\\b167\\.93\\b",
    "(167\\.93.*(null|permutation|expect))|((null|permutation|expect).*167\\.93)"
)

add_metric(
    "PERM_VERY_FOLD",
    "Permutation",
    "Very-strong fold enrichment",
    "1.346",
    "\\b1\\.346\\b",
    "(1\\.346.*(fold|enrichment))|((fold|enrichment).*1\\.346)"
)

add_metric(
    "PERM_VERY_Z",
    "Permutation",
    "Very-strong permutation Z",
    "5.03",
    "\\b5\\.03\\b",
    "(5\\.03.*(Z|permutation))|((Z|permutation).*5\\.03)"
)

add_metric(
    "PERM_P",
    "Permutation",
    "Permutation empirical P",
    "Pperm < 1e-4",
    "(P[_ ]?perm[^0-9]*<[^0-9]*(1\\s*[eE]-?0?4|10\\^?-?4))|(<\\s*10\\^?-?4)",
    "(P[_ ]?perm|permutation)"
)


###############################################################################
# MOTIF — GLOBAL
###############################################################################

add_metric(
    "FIMO_LINKS",
    "Motif",
    "FIMO q<=0.05 motif-peak links",
    "612",
    "\\b612\\b",
    "(612.*(FIMO|motif|link))|((FIMO|motif|link).*612)"
)

add_metric(
    "FIMO_MUSCLE_PEAKS",
    "Motif",
    "Muscle FIMO q<=0.05 unique peaks",
    "191",
    "\\b191\\b",
    "(191.*(Muscle|muscle).*(peak|FIMO))|((Muscle|muscle).*(peak|FIMO).*191)"
)

add_metric(
    "FIMO_SPLEEN_PEAKS",
    "Motif",
    "Spleen FIMO q<=0.05 unique peaks",
    "163",
    "\\b163\\b",
    "(163.*(Spleen|spleen).*(peak|FIMO))|((Spleen|spleen).*(peak|FIMO).*163)"
)

add_metric(
    "FIMO_LIVER_PEAKS",
    "Motif",
    "Liver FIMO q<=0.05 unique peaks",
    "26",
    "\\b26\\b",
    "(26.*(Liver|liver).*(peak|FIMO))|((Liver|liver).*(peak|FIMO).*26)"
)

add_metric(
    "STREME_EVALUABLE",
    "Motif",
    "STREME motifs with evaluable holdout E-values",
    "9",
    "\\b9\\b",
    "(9.*(evaluable|holdout|STREME))|((evaluable|holdout|STREME).*9)"
)

add_metric(
    "CEREBELLUM_TRAINING_ONLY",
    "Motif",
    "Cerebellum STREME training-only motifs",
    "5",
    "\\b5\\b",
    "(5.*(cerebell|training|STREME))|((cerebell|training|STREME).*5)"
)


###############################################################################
# MUSCLE KLF9
###############################################################################

add_metric(
    "KLF9_QPEAKS",
    "Motif",
    "Muscle KLF9 q-positive peaks",
    "121",
    "\\b121\\b",
    "(121.*KLF9)|(KLF9.*121)"
)

add_metric(
    "KLF9_PROM_POS",
    "Motif",
    "KLF9-positive promoter peaks",
    "30/121",
    "30\\s*(?:/|of)\\s*121",
    "(KLF9.*promoter)|(promoter.*KLF9)"
)

add_metric(
    "KLF9_PROM_POS_PCT",
    "Motif",
    "KLF9-positive promoter percentage",
    "24.79%",
    "24\\.79%?",
    "(KLF9.*24\\.79)|(24\\.79.*KLF9)"
)

add_metric(
    "KLF9_PROM_NEG",
    "Motif",
    "KLF9-negative promoter peaks",
    "26/223",
    "26\\s*(?:/|of)\\s*223",
    "(KLF9.*promoter)|(promoter.*KLF9)"
)

add_metric(
    "KLF9_PROM_NEG_PCT",
    "Motif",
    "KLF9-negative promoter percentage",
    "11.66%",
    "11\\.66%?",
    "(KLF9.*11\\.66)|(11\\.66.*KLF9)"
)

add_metric(
    "KLF9_OR",
    "Motif",
    "KLF9 promoter odds ratio",
    "2.491",
    "\\b2\\.491\\b",
    "(KLF9.*2\\.491)|(2\\.491.*KLF9)"
)

add_metric(
    "KLF9_PROM_FDR",
    "Motif",
    "KLF9 promoter FDR",
    "0.00651",
    "\\b0\\.00651\\b",
    "(KLF9.*0\\.00651)|(0\\.00651.*KLF9)"
)

add_metric(
    "KLF9_TSS_POS",
    "Motif",
    "KLF9-positive median TSS distance",
    "4,199 bp",
    "\\b4,?199\\b",
    "(KLF9.*4,?199)|(4,?199.*KLF9)"
)

add_metric(
    "KLF9_TSS_NEG",
    "Motif",
    "KLF9-negative median TSS distance",
    "15,777 bp",
    "\\b15,?777\\b",
    "(KLF9.*15,?777)|(15,?777.*KLF9)"
)

add_metric(
    "KLF9_TSS_FDR",
    "Motif",
    "KLF9 TSS-distance FDR",
    "3.19e-4",
    "(3\\.19[eE]-?0?4)|(3\\.19\\s*x\\s*10\\^?-?4)",
    "(KLF9.*3\\.19)|(3\\.19.*KLF9)"
)


###############################################################################
# MUSCLE SP-LIKE STREME3
###############################################################################

add_metric(
    "SP_QPEAKS",
    "Motif",
    "SP-like q-positive peaks",
    "138",
    "\\b138\\b",
    "(138.*(SP-like|STREME))|((SP-like|STREME).*138)"
)

add_metric(
    "SP_PROM_POS",
    "Motif",
    "SP-like positive promoter peaks",
    "40/138",
    "40\\s*(?:/|of)\\s*138",
    "(SP-like|STREME).*promoter|promoter.*(SP-like|STREME)"
)

add_metric(
    "SP_PROM_POS_PCT",
    "Motif",
    "SP-like positive promoter percentage",
    "28.99%",
    "28\\.99%?",
    "(SP-like|STREME).*28\\.99|28\\.99.*(SP-like|STREME)"
)

add_metric(
    "SP_PROM_NEG",
    "Motif",
    "SP-like negative promoter peaks",
    "16/206",
    "16\\s*(?:/|of)\\s*206",
    "(SP-like|STREME).*promoter|promoter.*(SP-like|STREME)"
)

add_metric(
    "SP_PROM_NEG_PCT",
    "Motif",
    "SP-like negative promoter percentage",
    "7.77%",
    "7\\.77%?",
    "(SP-like|STREME).*7\\.77|7\\.77.*(SP-like|STREME)"
)

add_metric(
    "SP_OR",
    "Motif",
    "SP-like promoter odds ratio",
    "4.823",
    "\\b4\\.823\\b",
    "(SP-like|STREME).*4\\.823|4\\.823.*(SP-like|STREME)"
)

add_metric(
    "SP_PROM_FDR",
    "Motif",
    "SP-like promoter FDR",
    "2.42e-6",
    "(2\\.42[eE]-?0?6)|(2\\.42\\s*x\\s*10\\^?-?6)",
    "(SP-like|STREME).*2\\.42|2\\.42.*(SP-like|STREME)"
)

add_metric(
    "SP_TSS_POS",
    "Motif",
    "SP-like positive median TSS distance",
    "4,246.5 bp",
    "4,?246\\.5",
    "(SP-like|STREME).*4,?246\\.5|4,?246\\.5.*(SP-like|STREME)"
)

add_metric(
    "SP_TSS_NEG",
    "Motif",
    "SP-like negative median TSS distance",
    "16,942 bp",
    "\\b16,?942\\b",
    "(SP-like|STREME).*16,?942|16,?942.*(SP-like|STREME)"
)

add_metric(
    "SP_TSS_FDR",
    "Motif",
    "SP-like TSS-distance FDR",
    "6.84e-7",
    "(6\\.84[eE]-?0?7)|(6\\.84\\s*x\\s*10\\^?-?7)",
    "(SP-like|STREME).*6\\.84|6\\.84.*(SP-like|STREME)"
)


###############################################################################
# SPLEEN OVERLAP
###############################################################################

add_metric(
    "SPLEEN_SHARED_PEAKS",
    "Motif",
    "Shared spleen STREME1/STREME2 peaks",
    "104",
    "\\b104\\b",
    "(104.*(spleen|STREME|shared))|((spleen|STREME|shared).*104)"
)

add_metric(
    "SPLEEN_PEAK_JACCARD",
    "Motif",
    "Spleen peak-level Jaccard",
    "0.65",
    "\\b0\\.65(?:0+)?\\b",
    "(0\\.65.*(Jaccard|spleen))|((Jaccard|spleen).*0\\.65)"
)

add_metric(
    "SPLEEN_SHARED_GENES",
    "Motif",
    "Shared spleen genes",
    "74",
    "\\b74\\b",
    "(74.*(gene|spleen|shared))|((gene|spleen|shared).*74)"
)

add_metric(
    "SPLEEN_GENE_JACCARD",
    "Motif",
    "Spleen gene-level Jaccard",
    "0.72549",
    "0\\.72549[0-9]*",
    "(0\\.72549.*(Jaccard|gene|spleen))|((Jaccard|gene|spleen).*0\\.72549)"
)


write.table(
    METRICS,
    file.path(
        OUTDIR,
        "01_canonical_frozen_metrics.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "Canonical metrics registered: ",
    nrow(METRICS),
    "\n\n",
    sep = ""
)


###############################################################################
# SOURCE TABLE RECOMPUTATION — CORE RESULTS
###############################################################################

cat(
    "====================================================================\n",
    "3. RECOMPUTING CORE VALUES FROM FROZEN SOURCE TABLES\n",
    "====================================================================\n",
    sep = ""
)


source_checks <- data.frame(
    ID = character(),
    Source_file = character(),
    Recomputed_value = character(),
    Expected_value = character(),
    Status = character(),
    Note = character(),
    stringsAsFactors = FALSE
)


add_source_check <- function(
    id,
    file,
    actual,
    expected,
    status,
    note = ""
) {

    source_checks <<- rbind(
        source_checks,
        data.frame(
            ID = id,
            Source_file = file,
            Recomputed_value = as.character(actual),
            Expected_value = as.character(expected),
            Status = status,
            Note = note,
            stringsAsFactors = FALSE
        )
    )
}


###############################################################################
# Step03 RNA / ATAC
###############################################################################

RNA_STEP03 <- file.path(
    ROOT,
    "03_tissue_specificity",
    "RNA",
    "RNA_tissue_specificity_all.tsv"
)

ATAC_STEP03 <- file.path(
    ROOT,
    "03_tissue_specificity",
    "ATAC",
    "ATAC_tissue_specificity_all.tsv"
)


recompute_step03 <- function(
    file,
    assay,
    expected_lrt,
    expected_hc
) {

    x <- safe_read_table(file)

    if (is.null(x)) {

        add_source_check(
            paste0(assay, "_LRT"),
            file,
            NA,
            expected_lrt,
            "MISSING_SOURCE"
        )

        add_source_check(
            paste0(assay, "_HC"),
            file,
            NA,
            expected_hc,
            "MISSING_SOURCE"
        )

        return(NULL)
    }


    required <- c(
        "LRT_padj",
        "High_confidence_tissue_specific"
    )

    if (!all(required %in% colnames(x))) {

        add_source_check(
            paste0(assay, "_LRT"),
            file,
            NA,
            expected_lrt,
            "COLUMN_MISSING",
            paste(
                setdiff(required, colnames(x)),
                collapse = ", "
            )
        )

        return(NULL)
    }


    lrt_n <- sum(
        !is.na(x$LRT_padj) &
        x$LRT_padj < 0.01
    )

    hc <- bool_value(
        x$High_confidence_tissue_specific
    )

    hc_n <- sum(
        hc,
        na.rm = TRUE
    )


    add_source_check(
        paste0(assay, "_LRT"),
        file,
        lrt_n,
        expected_lrt,
        ifelse(
            lrt_n == expected_lrt,
            "PASS",
            "FAIL"
        )
    )

    add_source_check(
        paste0(assay, "_HC"),
        file,
        hc_n,
        expected_hc,
        ifelse(
            hc_n == expected_hc,
            "PASS",
            "FAIL"
        )
    )


    if (
        "Max_tissue" %in% colnames(x)
    ) {

        tt <- table(
            x$Max_tissue[hc]
        )

        out <- data.frame(
            Assay = assay,
            Tissue = names(tt),
            HC_N = as.integer(tt),
            stringsAsFactors = FALSE
        )

        write.table(
            out,
            file.path(
                OUTDIR,
                paste0(
                    "02_",
                    assay,
                    "_HC_tissue_counts_recomputed.tsv"
                )
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }


    x
}


rna03 <- recompute_step03(
    RNA_STEP03,
    "RNA",
    18870,
    5036
)

atac03 <- recompute_step03(
    ATAC_STEP03,
    "ATAC",
    176917,
    40545
)


###############################################################################
# Concordant 4,674 table
###############################################################################

CONCORDANT_FILE <- file.path(
    ROOT,
    "04_RNA_ATAC_peak_gene_master",
    "04_concordant_HC_RNA_ATAC_peak_gene_pairs.tsv"
)

con <- safe_read_table(
    CONCORDANT_FILE
)

if (!is.null(con)) {

    n_con <- nrow(con)

    gene_col <- detect_column(
        con,
        c(
            "^gene_id$",
            "gene.*id",
            "RNA.*gene"
        )
    )


    add_source_check(
        "CONCORDANT_PAIRS",
        CONCORDANT_FILE,
        n_con,
        4674,
        ifelse(
            n_con == 4674,
            "PASS",
            "FAIL"
        )
    )


    if (!is.na(gene_col)) {

        n_gene <- length(
            unique(
                con[[gene_col]][
                    !is.na(con[[gene_col]]) &
                    con[[gene_col]] != ""
                ]
            )
        )

        add_source_check(
            "CONCORDANT_GENES",
            CONCORDANT_FILE,
            n_gene,
            1417,
            ifelse(
                n_gene == 1417,
                "PASS",
                "FAIL"
            ),
            paste0(
                "gene column=",
                gene_col
            )
        )
    }


    tissue_col <- detect_column(
        con,
        c(
            "^Tissue$",
            "tissue"
        )
    )

    if (!is.na(tissue_col)) {

        ct <- as.data.frame(
            table(
                con[[tissue_col]]
            )
        )

        colnames(ct) <- c(
            "Tissue",
            "Concordant_pair_N"
        )

        write.table(
            ct,
            file.path(
                OUTDIR,
                "03_concordant_pair_tissue_counts_recomputed.tsv"
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }

} else {

    add_source_check(
        "CONCORDANT_PAIRS",
        CONCORDANT_FILE,
        NA,
        4674,
        "MISSING_SOURCE"
    )
}


###############################################################################
# Coupling 3,916 / 702 / 226 / gene counts
###############################################################################

COUPLING_FILE <- file.path(
    ROOT,
    "05_RNA_ATAC_quantitative_coupling",
    "05_concordant_HC_peak_gene_correlations.tsv"
)

coupling <- safe_read_table(
    COUPLING_FILE
)


if (!is.null(coupling)) {

    c348 <- detect_column(
        coupling,
        c(
            "rho.*348",
            "348.*rho",
            "cor.*348",
            "348.*cor"
        ),
        numeric_only = TRUE
    )

    c350 <- detect_column(
        coupling,
        c(
            "rho.*350",
            "350.*rho",
            "cor.*350",
            "350.*cor"
        ),
        numeric_only = TRUE
    )

    cmean <- detect_column(
        coupling,
        c(
            "rho.*mean",
            "mean.*rho",
            "cor.*mean",
            "mean.*cor"
        ),
        numeric_only = TRUE
    )


    cat(
        "Detected coupling columns:\n",
        "P348 = ", c348, "\n",
        "P350 = ", c350, "\n",
        "Mean = ", cmean, "\n\n",
        sep = ""
    )


    if (
        !any(
            is.na(
                c(
                    c348,
                    c350,
                    cmean
                )
            )
        )
    ) {

        r348 <- as.numeric(
            coupling[[c348]]
        )

        r350 <- as.numeric(
            coupling[[c350]]
        )

        rmean <- as.numeric(
            coupling[[cmean]]
        )

        valid <- is.finite(r348) &
            is.finite(r350) &
            is.finite(rmean)


        positive <- valid &
            r348 > 0 &
            r350 > 0 &
            rmean > 0


        strong <- valid &
            r348 >= 0.50 &
            r350 >= 0.50 &
            rmean >= 0.70


        very <- valid &
            r348 >= 0.70 &
            r350 >= 0.70 &
            rmean >= 0.80


        add_source_check(
            "POSITIVE_ALL_THREE",
            COUPLING_FILE,
            sum(positive),
            3916,
            ifelse(
                sum(positive) == 3916,
                "PASS",
                "FAIL"
            )
        )


        add_source_check(
            "STRONG_PAIRS",
            COUPLING_FILE,
            sum(strong),
            702,
            ifelse(
                sum(strong) == 702,
                "PASS",
                "FAIL"
            )
        )


        add_source_check(
            "VERY_STRONG_PAIRS",
            COUPLING_FILE,
            sum(very),
            226,
            ifelse(
                sum(very) == 226,
                "PASS",
                "FAIL"
            )
        )


        gene_col <- detect_column(
            coupling,
            c(
                "^gene_id$",
                "gene.*id"
            )
        )


        if (!is.na(gene_col)) {

            strong_genes <- length(
                unique(
                    coupling[[gene_col]][strong]
                )
            )

            very_genes <- length(
                unique(
                    coupling[[gene_col]][very]
                )
            )


            add_source_check(
                "STRONG_GENES",
                COUPLING_FILE,
                strong_genes,
                450,
                ifelse(
                    strong_genes == 450,
                    "PASS",
                    "FAIL"
                )
            )


            add_source_check(
                "VERY_STRONG_GENES",
                COUPLING_FILE,
                very_genes,
                176,
                ifelse(
                    very_genes == 176,
                    "PASS",
                    "FAIL"
                )
            )
        }


        tissue_col <- detect_column(
            coupling,
            c(
                "^Tissue$",
                "tissue"
            )
        )


        if (!is.na(tissue_col)) {

            strong_tissue <- as.data.frame(
                table(
                    coupling[[tissue_col]][strong]
                )
            )

            colnames(strong_tissue) <- c(
                "Tissue",
                "Strong_pair_N"
            )

            write.table(
                strong_tissue,
                file.path(
                    OUTDIR,
                    "04_strong_pair_tissue_counts_recomputed.tsv"
                ),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE
            )
        }

    } else {

        add_source_check(
            "STRONG_PAIRS",
            COUPLING_FILE,
            NA,
            702,
            "COLUMN_MISSING",
            "Could not detect P348/P350/mean rho columns"
        )
    }

} else {

    add_source_check(
        "STRONG_PAIRS",
        COUPLING_FILE,
        NA,
        702,
        "MISSING_SOURCE"
    )
}


###############################################################################
# PROMOTER FIRST TABLES
###############################################################################

ALL_PROM_FILE <- file.path(
    ROOT,
    "07_promoter_first_classification",
    "07A_all_UROPA_assigned_promoter_first.tsv"
)

HC_PROM_FILE <- file.path(
    ROOT,
    "07_promoter_first_enrichment",
    "07B_4674_HC_promoter_first.tsv"
)

STRONG_PROM_FILE <- file.path(
    ROOT,
    "07_promoter_first_enrichment",
    "07B_702_strong_promoter_first.tsv"
)



standard_class <- function(x) {

    raw <- tolower(
        trimws(
            as.character(x)
        )
    )

    z <- gsub(
        "[^a-z0-9]+",
        "",
        raw
    )

    out <- rep(
        NA_character_,
        length(z)
    )

    out[
        grepl("^promoter", z)
    ] <- "Promoter"

    out[
        grepl(
            "^(5utr|fiveprimeutr|fiveutr)",
            z
        )
    ] <- "5UTR"

    out[
        grepl(
            "^(3utr|threeprimeutr|threeutr)",
            z
        )
    ] <- "3UTR"

    out[
        grepl("^exon", z)
    ] <- "Exon"

    out[
        grepl("^intron", z)
    ] <- "Intron"

    out[
        grepl("^distal", z)
    ] <- "Distal"

    out
}


detect_promoter_class <- function(df) {

    nm <- colnames(df)

    if ("Pair_specific_class" %in% nm) {

        cat(
            "Official Step07 classification column: ",
            "Pair_specific_class\n",
            sep = ""
        )

        return(
            "Pair_specific_class"
        )
    }

    stop(
        paste0(
            "Pair_specific_class is missing. ",
            "Step07 promoter statistics cannot be audited safely."
        ),
        call. = FALSE
    )
}

detect_tss_col <- function(df) {

    cand <- colnames(df)[
        grepl(
            "TSS.*distance|distance.*TSS",
            colnames(df),
            ignore.case = TRUE
        )
    ]

    if (length(cand) == 0) {
        return(NA_character_)
    }

    abs_cand <- cand[
        grepl(
            "abs|absolute",
            cand,
            ignore.case = TRUE
        )
    ]

    if (length(abs_cand) > 0) {
        return(abs_cand[1])
    }

    cand[1]
}


read_promoter_table <- function(
    file,
    label
) {

    x <- safe_read_table(file)

    if (is.null(x)) {
        return(NULL)
    }

    class_col <- detect_promoter_class(
        x
    )

    tss_col <- detect_tss_col(
        x
    )

    cat(
        label,
        " promoter class column = ",
        class_col,
        "\n",
        label,
        " TSS distance column    = ",
        tss_col,
        "\n",
        sep = ""
    )


    classes <- if (!is.na(class_col)) {
        standard_class(
            x[[class_col]]
        )
    } else {
        rep(
            NA_character_,
            nrow(x)
        )
    }


    tss <- if (!is.na(tss_col)) {

        abs(
            suppressWarnings(
                as.numeric(
                    x[[tss_col]]
                )
            )
        )

    } else {

        rep(
            NA_real_,
            nrow(x)
        )
    }


    list(
        data = x,
        classes = classes,
        tss = tss,
        class_col = class_col,
        tss_col = tss_col
    )
}


all_prom <- read_promoter_table(
    ALL_PROM_FILE,
    "ALL"
)

hc_prom <- read_promoter_table(
    HC_PROM_FILE,
    "HC"
)

strong_prom <- read_promoter_table(
    STRONG_PROM_FILE,
    "STRONG"
)


if (!is.null(all_prom)) {

    add_source_check(
        "ALL_ASSIGNED",
        ALL_PROM_FILE,
        nrow(all_prom$data),
        168999,
        ifelse(
            nrow(all_prom$data) == 168999,
            "PASS",
            "FAIL"
        )
    )

    np <- sum(
        all_prom$classes == "Promoter",
        na.rm = TRUE
    )

    add_source_check(
        "ALL_PROMOTER",
        ALL_PROM_FILE,
        np,
        16677,
        ifelse(
            np == 16677,
            "PASS",
            "FAIL"
        )
    )

    med <- median(
        all_prom$tss[
            is.finite(all_prom$tss)
        ],
        na.rm = TRUE
    )

    add_source_check(
        "TSS_ALL_MEDIAN",
        ALL_PROM_FILE,
        med,
        33378,
        ifelse(
            isTRUE(
                all.equal(
                    med,
                    33378,
                    tolerance = 1e-8
                )
            ),
            "PASS",
            "FAIL"
        )
    )
}


if (!is.null(hc_prom)) {

    med <- median(
        hc_prom$tss[
            is.finite(hc_prom$tss)
        ],
        na.rm = TRUE
    )

    add_source_check(
        "TSS_HC_MEDIAN",
        HC_PROM_FILE,
        med,
        17248.5,
        ifelse(
            isTRUE(
                all.equal(
                    med,
                    17248.5,
                    tolerance = 1e-8
                )
            ),
            "PASS",
            "FAIL"
        )
    )
}


if (!is.null(strong_prom)) {

    np <- sum(
        strong_prom$classes == "Promoter",
        na.rm = TRUE
    )

    add_source_check(
        "STRONG_PROMOTER",
        STRONG_PROM_FILE,
        np,
        121,
        ifelse(
            np == 121,
            "PASS",
            "FAIL"
        )
    )

    pct <- 100 * np /
        nrow(
            strong_prom$data
        )

    add_source_check(
        "STRONG_PROMOTER_PCT",
        STRONG_PROM_FILE,
        sprintf("%.4f%%", pct),
        "17.2365%",
        ifelse(
            abs(
                pct -
                17.2364672364672
            ) < 0.001,
            "PASS",
            "FAIL"
        )
    )

    med <- median(
        strong_prom$tss[
            is.finite(strong_prom$tss)
        ],
        na.rm = TRUE
    )

    add_source_check(
        "TSS_STRONG_MEDIAN",
        STRONG_PROM_FILE,
        med,
        11549,
        ifelse(
            isTRUE(
                all.equal(
                    med,
                    11549,
                    tolerance = 1e-8
                )
            ),
            "PASS",
            "FAIL"
        )
    )
}


###############################################################################
# RECOMPUTE PROMOTER ODDS RATIOS IF TABLES ARE AVAILABLE
###############################################################################

if (
    !is.null(all_prom) &&
    !is.null(hc_prom) &&
    !is.null(strong_prom)
) {

    strong_total <- nrow(
        strong_prom$data
    )

    strong_p <- sum(
        strong_prom$classes == "Promoter",
        na.rm = TRUE
    )

    strong_nonp <- strong_total -
        strong_p


    all_total <- nrow(
        all_prom$data
    )

    all_p <- sum(
        all_prom$classes == "Promoter",
        na.rm = TRUE
    )

    all_nonstrong_p <- all_p -
        strong_p

    all_nonstrong_total <-
        all_total -
        strong_total

    all_nonstrong_nonp <-
        all_nonstrong_total -
        all_nonstrong_p


    ft_all <- fisher.test(
        matrix(
            c(
                strong_p,
                strong_nonp,
                all_nonstrong_p,
                all_nonstrong_nonp
            ),
            nrow = 2,
            byrow = TRUE
        )
    )


    add_source_check(
        "PROMOTER_OR_ALL",
        paste(
            STRONG_PROM_FILE,
            ALL_PROM_FILE,
            sep = " + "
        ),
        unname(ft_all$estimate),
        1.9089,
        ifelse(
            abs(
                unname(ft_all$estimate) -
                1.9089
            ) < 0.01,
            "PASS",
            "FAIL"
        ),
        "Recomputed Fisher OR using disjoint all-non-strong reference"
    )


    hc_total <- nrow(
        hc_prom$data
    )

    hc_p <- sum(
        hc_prom$classes == "Promoter",
        na.rm = TRUE
    )

    hc_nonstrong_p <-
        hc_p -
        strong_p

    hc_nonstrong_total <-
        hc_total -
        strong_total

    hc_nonstrong_nonp <-
        hc_nonstrong_total -
        hc_nonstrong_p


    ft_hc <- fisher.test(
        matrix(
            c(
                strong_p,
                strong_nonp,
                hc_nonstrong_p,
                hc_nonstrong_nonp
            ),
            nrow = 2,
            byrow = TRUE
        )
    )


    add_source_check(
        "PROMOTER_OR_HC",
        paste(
            STRONG_PROM_FILE,
            HC_PROM_FILE,
            sep = " + "
        ),
        unname(ft_hc$estimate),
        1.7517,
        ifelse(
            abs(
                unname(ft_hc$estimate) -
                1.7517
            ) < 0.01,
            "PASS",
            "FAIL"
        ),
        "Recomputed Fisher OR using disjoint HC-non-strong reference"
    )
}


###############################################################################
# WRITE SOURCE RECOMPUTATION
###############################################################################

write.table(
    source_checks,
    file.path(
        OUTDIR,
        "05_source_recomputed_core_metrics.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\nSource recomputation summary:\n")

print(
    source_checks,
    row.names = FALSE
)

cat("\n")


###############################################################################
# SOURCE FILE NUMERIC SEARCH
#
# Used especially for:
#   FDR values
#   permutation summaries
#   motif architecture statistics
#
###############################################################################

cat(
    "====================================================================\n",
    "4. SEARCHING FROZEN SOURCE FILES FOR FORMAL NUMERIC STATISTICS\n",
    "====================================================================\n",
    sep = ""
)


top_dirs <- list.dirs(
    ROOT,
    full.names = TRUE,
    recursive = FALSE
)

SOURCE_DIRS <- top_dirs[
    grepl(
        paste0(
            "(/|^)",
            "(",
            "05B_|",
            "06_|",
            "07_|",
            "08",
            "|09B_|",
            "09C_",
            ")"
        ),
        top_dirs
    )
]


source_files <- unique(
    unlist(
        lapply(
            SOURCE_DIRS,
            function(d) {
                list.files(
                    d,
                    recursive = TRUE,
                    full.names = TRUE,
                    pattern = "\\.(tsv|txt|csv|log)$",
                    ignore.case = TRUE
                )
            }
        )
    )
)


source_files <- source_files[
    file.exists(source_files)
]


if (length(source_files) > 0) {

    sizes <- file.info(
        source_files
    )$size /
        1024^2

    # Large pair-level files are already handled via explicit recomputation.
    # Numeric token scanning focuses on summary/statistics files.
    source_files_scan <- source_files[
        is.finite(sizes) &
        sizes <= 25
    ]

} else {

    source_files_scan <- character(0)
}


write.table(
    data.frame(
        Source_file = source_files_scan,
        Size_MB = if (
            length(source_files_scan) > 0
        ) {
            round(
                file.info(
                    source_files_scan
                )$size /
                    1024^2,
                3
            )
        } else {
            numeric(0)
        }
    ),
    file.path(
        OUTDIR,
        "06_source_files_scanned.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


source_hits <- data.frame(
    ID = character(),
    Expected = character(),
    Source_file = character(),
    Line_number = integer(),
    Source_line = character(),
    stringsAsFactors = FALSE
)


scan_metric_source <- function(metric_row) {

    if (length(source_files_scan) == 0) {
        return(NULL)
    }

    regex <- metric_row$Expected_regex

    for (f in source_files_scan) {

        x <- tryCatch(
            readLines(
                f,
                warn = FALSE
            ),
            error = function(e) {
                character(0)
            }
        )

        if (length(x) == 0) {
            next
        }

        x <- normalize_text(x)

        hit <- grep(
            regex,
            x,
            perl = TRUE,
            ignore.case = TRUE
        )

        if (length(hit) > 0) {

            # Limit to first 20 lines per metric/file.
            hit <- head(
                hit,
                20
            )

            for (h in hit) {

                source_hits <<- rbind(
                    source_hits,
                    data.frame(
                        ID = metric_row$ID,
                        Expected =
                            metric_row$Expected_display,
                        Source_file =
                            sub(
                                paste0(
                                    "^",
                                    ROOT,
                                    "/?"
                                ),
                                "",
                                f
                            ),
                        Line_number = h,
                        Source_line = x[h],
                        stringsAsFactors = FALSE
                    )
                )
            }
        }
    }
}


# Source scan is particularly useful for promoter/permutation/motif values.
scan_ids <- METRICS$ID[
    METRICS$Section %in%
        c(
            "Promoter",
            "Permutation",
            "Motif"
        )
]


for (id in scan_ids) {

    m <- METRICS[
        METRICS$ID == id,
        ,
        drop = FALSE
    ]

    scan_metric_source(
        m[1, ]
    )
}


write.table(
    source_hits,
    file.path(
        OUTDIR,
        "07_frozen_source_numeric_hits.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


###############################################################################
# MANUSCRIPT CONTEXT AUDIT
###############################################################################

cat(
    "====================================================================\n",
    "5. AUDITING MANUSCRIPT NUMBERS AND CONTEXT\n",
    "====================================================================\n",
    sep = ""
)


manuscript_audit <- data.frame(
    ID = character(),
    Section = character(),
    Metric = character(),
    Expected = character(),
    Global_number_hits = integer(),
    Context_paragraph_hits = integer(),
    Correct_context_hits = integer(),
    Status = character(),
    Paragraphs = character(),
    Line_ranges = character(),
    Context_numbers = character(),
    Context_preview = character(),
    stringsAsFactors = FALSE
)


for (i in seq_len(nrow(METRICS))) {

    m <- METRICS[i, ]

    global_hits <- grep(
        m$Expected_regex,
        MANUSCRIPT_LINES,
        perl = TRUE,
        ignore.case = TRUE
    )


    context_idx <- grep(
        m$Context_regex,
        PARAGRAPHS$Text,
        perl = TRUE,
        ignore.case = TRUE
    )


    correct_idx <- context_idx[
        grepl(
            m$Expected_regex,
            PARAGRAPHS$Text[
                context_idx
            ],
            perl = TRUE,
            ignore.case = TRUE
        )
    ]


    if (length(correct_idx) >= 1) {

        status <- "PASS"

    } else if (
        length(context_idx) >= 1 &&
        length(global_hits) >= 1
    ) {

        status <- "WARNING_NUMBER_OUTSIDE_CONTEXT"

    } else if (
        length(context_idx) >= 1 &&
        length(global_hits) == 0
    ) {

        status <- "WARNING_CONTEXT_VALUE_MISSING"

    } else if (
        length(context_idx) == 0 &&
        length(global_hits) >= 1
    ) {

        status <- "WARNING_NUMBER_WITHOUT_EXPECTED_CONTEXT"

    } else {

        status <- "MISSING"
    }


    inspect_idx <- if (
        length(correct_idx) > 0
    ) {
        correct_idx
    } else {
        context_idx
    }


    paragraphs_str <- if (
        length(inspect_idx) > 0
    ) {
        paste(
            PARAGRAPHS$Paragraph[
                inspect_idx
            ],
            collapse = ";"
        )
    } else {
        ""
    }


    line_ranges <- if (
        length(inspect_idx) > 0
    ) {
        paste(
            paste0(
                PARAGRAPHS$Start_line[
                    inspect_idx
                ],
                "-",
                PARAGRAPHS$End_line[
                    inspect_idx
                ]
            ),
            collapse = ";"
        )
    } else {
        ""
    }


    nums <- if (
        length(inspect_idx) > 0
    ) {

        paste(
            unique(
                unlist(
                    lapply(
                        PARAGRAPHS$Text[
                            inspect_idx
                        ],
                        extract_numbers
                    )
                )
            ),
            collapse = " | "
        )

    } else {
        ""
    }


    previews <- if (
        length(inspect_idx) > 0
    ) {

        z <- PARAGRAPHS$Text[
            head(
                inspect_idx,
                3
            )
        ]

        z <- substr(
            z,
            1,
            500
        )

        paste(
            z,
            collapse = " || "
        )

    } else {
        ""
    }


    manuscript_audit <- rbind(
        manuscript_audit,
        data.frame(
            ID = m$ID,
            Section = m$Section,
            Metric = m$Metric,
            Expected =
                m$Expected_display,
            Global_number_hits =
                length(global_hits),
            Context_paragraph_hits =
                length(context_idx),
            Correct_context_hits =
                length(correct_idx),
            Status = status,
            Paragraphs =
                paragraphs_str,
            Line_ranges =
                line_ranges,
            Context_numbers =
                nums,
            Context_preview =
                previews,
            stringsAsFactors = FALSE
        )
    )
}


write.table(
    manuscript_audit,
    file.path(
        OUTDIR,
        "08_manuscript_context_numeric_audit.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


###############################################################################
# CORE TISSUE COUNT CHECKS IN MANUSCRIPT
###############################################################################

cat(
    "====================================================================\n",
    "6. AUDITING TISSUE-SPECIFIC COUNT VECTORS\n",
    "====================================================================\n",
    sep = ""
)


TISSUES <- c(
    "Adipose",
    "Cerebellum",
    "Cortex",
    "Hypothalamus",
    "Liver",
    "Lung",
    "Muscle",
    "Spleen"
)


expected_tissue <- rbind(
    data.frame(
        Set = "RNA_HC",
        Tissue = TISSUES,
        Expected = c(
            408,
            573,
            397,
            244,
            1182,
            610,
            733,
            889
        )
    ),

    data.frame(
        Set = "ATAC_HC",
        Tissue = TISSUES,
        Expected = c(
            47,
            4643,
            2689,
            785,
            5071,
            725,
            20853,
            5732
        )
    ),

    data.frame(
        Set = "Concordant",
        Tissue = TISSUES,
        Expected = c(
            4,
            197,
            66,
            16,
            481,
            38,
            3065,
            807
        )
    ),

    data.frame(
        Set = "Strong",
        Tissue = TISSUES,
        Expected = c(
            0,
            47,
            15,
            10,
            54,
            6,
            344,
            226
        )
    )
)


whole_text <- paste(
    MANUSCRIPT_LINES,
    collapse = " "
)


tissue_audit <- expected_tissue


tissue_audit$Number_present <- mapply(
    function(tissue, value) {

        number_regex <- if (value >= 1000) {

            z <- format(
                value,
                big.mark = ",",
                scientific = FALSE,
                trim = TRUE
            )

            paste0(
                "(",
                value,
                "|",
                gsub(
                    ",",
                    ",?",
                    z,
                    fixed = TRUE
                ),
                ")"
            )

        } else {

            paste0(
                "\\b",
                value,
                "\\b"
            )
        }


        pattern <- paste0(
            "(",
            tissue,
            ".{0,80}",
            number_regex,
            ")|(",
            number_regex,
            ".{0,80}",
            tissue,
            ")"
        )

        grepl(
            pattern,
            whole_text,
            ignore.case = TRUE,
            perl = TRUE
        )

    },
    tissue_audit$Tissue,
    tissue_audit$Expected
)


tissue_audit$Status <- ifelse(
    tissue_audit$Number_present,
    "PASS",
    "NOT_FOUND_NEAR_TISSUE"
)


write.table(
    tissue_audit,
    file.path(
        OUTDIR,
        "09_tissue_count_manuscript_audit.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


###############################################################################
# NUMERIC INVENTORY
###############################################################################

cat(
    "====================================================================\n",
    "7. BUILDING MANUSCRIPT NUMERIC INVENTORY\n",
    "====================================================================\n",
    sep = ""
)


number_pattern <- paste0(
    "(?<![A-Za-z])",
    "[-+]?",
    "[0-9][0-9,]*",
    "(?:\\.[0-9]+)?",
    "(?:[eE][-+]?[0-9]+)?",
    "%?"
)


inventory <- data.frame(
    Number = character(),
    Line = integer(),
    Context = character(),
    stringsAsFactors = FALSE
)


for (
    i in seq_along(
        MANUSCRIPT_LINES
    )
) {

    z <- regmatches(
        MANUSCRIPT_LINES[i],
        gregexpr(
            number_pattern,
            MANUSCRIPT_LINES[i],
            perl = TRUE
        )
    )[[1]]


    if (
        length(z) == 1 &&
        identical(
            z,
            character(0)
        )
    ) {
        next
    }


    for (n in z) {

        inventory <- rbind(
            inventory,
            data.frame(
                Number = n,
                Line = i,
                Context = MANUSCRIPT_LINES[i],
                stringsAsFactors = FALSE
            )
        )
    }
}


write.table(
    inventory,
    file.path(
        OUTDIR,
        "10_manuscript_numeric_inventory.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


if (nrow(inventory) > 0) {

    freq <- sort(
        table(
            inventory$Number
        ),
        decreasing = TRUE
    )

    freq_df <- data.frame(
        Number = names(freq),
        Frequency = as.integer(freq),
        stringsAsFactors = FALSE
    )

    write.table(
        freq_df,
        file.path(
            OUTDIR,
            "11_manuscript_numeric_frequency.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


###############################################################################
# CROSS-AUDIT FINAL SUMMARY
###############################################################################

cat(
    "====================================================================\n",
    "8. FINAL CROSS-AUDIT SUMMARY\n",
    "====================================================================\n",
    sep = ""
)


final <- merge(
    METRICS[
        ,
        c(
            "ID",
            "Section",
            "Metric",
            "Expected_display",
            "Importance"
        )
    ],
    manuscript_audit[
        ,
        c(
            "ID",
            "Status",
            "Line_ranges",
            "Context_numbers"
        )
    ],
    by = "ID",
    all.x = TRUE
)


colnames(final)[
    colnames(final) == "Status"
] <- "Manuscript_status"


source_summary <- aggregate(
    Status ~ ID,
    data = source_checks,
    FUN = function(x) {

        if (any(x == "FAIL")) {
            "FAIL"
        } else if (
            any(
                x %in%
                c(
                    "MISSING_SOURCE",
                    "COLUMN_MISSING"
                )
            )
        ) {
            "WARNING"
        } else if (
            all(x == "PASS")
        ) {
            "PASS"
        } else {
            paste(
                unique(x),
                collapse = ";"
            )
        }
    }
)


final <- merge(
    final,
    source_summary,
    by = "ID",
    all.x = TRUE
)


colnames(final)[
    colnames(final) == "Status"
] <- "Source_recompute_status"


source_hit_counts <- if (
    nrow(source_hits) > 0
) {

    aggregate(
        Source_file ~ ID,
        data = source_hits,
        FUN = function(x) {
            length(
                unique(x)
            )
        }
    )

} else {

    data.frame(
        ID = character(),
        Source_file = integer()
    )
}


colnames(source_hit_counts) <- c(
    "ID",
    "Source_files_with_numeric_hit"
)


final <- merge(
    final,
    source_hit_counts,
    by = "ID",
    all.x = TRUE
)


final$Source_files_with_numeric_hit[
    is.na(
        final$Source_files_with_numeric_hit
    )
] <- 0



# ==================================================================
# FINAL STATUS LOGIC — NA SAFE
# ==================================================================

final$Source_recompute_status[
    is.na(final$Source_recompute_status)
] <- "NOT_RECOMPUTED"

final$Manuscript_status[
    is.na(final$Manuscript_status)
] <- "NOT_AUDITED"

final$Source_files_with_numeric_hit[
    is.na(final$Source_files_with_numeric_hit)
] <- 0


final$Final_status <- rep(
    "REVIEW",
    nrow(final)
)


# ----------------------------------------------------------
# Genuine FAIL:
# source recomputation contradicts frozen value.
# ----------------------------------------------------------

source_fail <- (
    final$Source_recompute_status == "FAIL"
)


# ----------------------------------------------------------
# Missing manuscript metric.
#
# It is kept as REVIEW rather than automatic FAIL because
# some registered metrics may appropriately appear only in
# supplementary material or may not need to be stated in
# the main manuscript.
# ----------------------------------------------------------

manuscript_missing <- (
    final$Manuscript_status == "MISSING"
)


# ----------------------------------------------------------
# PASS:
# correct manuscript context AND either:
#   a) directly recomputed from frozen source, or
#   b) found in a frozen result file.
# ----------------------------------------------------------

pass_idx <- (
    final$Manuscript_status == "PASS"
) &
(
    final$Source_recompute_status == "PASS" |
    final$Source_files_with_numeric_hit > 0
)


final$Final_status[
    pass_idx
] <- "PASS"


# ----------------------------------------------------------
# Direct source contradiction always overrides PASS.
# ----------------------------------------------------------

final$Final_status[
    source_fail
] <- "FAIL"


# ----------------------------------------------------------
# Manuscript missing but no source contradiction = REVIEW.
# ----------------------------------------------------------

final$Final_status[
    manuscript_missing &
    !source_fail
] <- "REVIEW"

write.table(
    final,
    file.path(
        OUTDIR,
        "12_FINAL_cross_audit_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


###############################################################################
# HUMAN REVIEW FILE
###############################################################################

review <- final[
    final$Final_status != "PASS",
    ,
    drop = FALSE
]


write.table(
    review,
    file.path(
        OUTDIR,
        "13_ITEMS_REQUIRING_MANUAL_REVIEW.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


###############################################################################
# FINAL COUNTS
###############################################################################


n_pass <- sum(
    final$Final_status == "PASS",
    na.rm = TRUE
)

n_review <- sum(
    final$Final_status == "REVIEW",
    na.rm = TRUE
)

n_fail <- sum(
    final$Final_status == "FAIL",
    na.rm = TRUE
)

cat("\n")
cat(
    "====================================================================\n",
    "FINAL AUDIT RESULT\n",
    "====================================================================\n",
    "Canonical metrics audited : ",
    nrow(final),
    "\n",
    "PASS                      : ",
    n_pass,
    "\n",
    "REVIEW                    : ",
    n_review,
    "\n",
    "FAIL                      : ",
    n_fail,
    "\n",
    "====================================================================\n",
    sep = ""
)


if (
    n_fail == 0 &&
    n_review == 0
) {

    cat(
        "OVERALL STATUS: PASS\n",
        "All registered numeric claims are internally consistent.\n"
    )

} else if (
    n_fail == 0
) {

    cat(
        "OVERALL STATUS: PASS WITH MANUAL REVIEW\n",
        "No direct contradiction detected, but some items require review.\n"
    )

} else {

    cat(
        "OVERALL STATUS: FAIL — CORRECTION REQUIRED\n",
        "At least one manuscript/source inconsistency was detected.\n"
    )
}


cat(
    "\nPrimary manual-review file:\n",
    file.path(
        OUTDIR,
        "13_ITEMS_REQUIRING_MANUAL_REVIEW.tsv"
    ),
    "\n\nFinal cross-audit table:\n",
    file.path(
        OUTDIR,
        "12_FINAL_cross_audit_summary.tsv"
    ),
    "\n",
    sep = ""
)


###############################################################################
# README
###############################################################################

README <- c(
    "STEP12B — Mammalian Genome manuscript numeric consistency audit",
    "",
    "This audit is read-only.",
    "No biological result file is modified.",
    "",
    "Primary outputs:",
    "",
    "01_canonical_frozen_metrics.tsv",
    "  Canonical frozen values used for manuscript checking.",
    "",
    "05_source_recomputed_core_metrics.tsv",
    "  Metrics recomputed directly from frozen Step03-Step07 source tables.",
    "",
    "07_frozen_source_numeric_hits.tsv",
    "  Promoter/permutation/motif statistics located in frozen source files.",
    "",
    "08_manuscript_context_numeric_audit.tsv",
    "  Checks whether each number occurs in the correct manuscript context.",
    "",
    "09_tissue_count_manuscript_audit.tsv",
    "  Tissue-by-tissue count consistency check.",
    "",
    "10_manuscript_numeric_inventory.tsv",
    "  Every numeric token detected in the manuscript.",
    "",
    "12_FINAL_cross_audit_summary.tsv",
    "  Main combined audit result.",
    "",
    "13_ITEMS_REQUIRING_MANUAL_REVIEW.tsv",
    "  The only file that normally requires manual inspection.",
    "",
    "Interpretation:",
    "PASS   = expected metric supported and correctly represented.",
    "REVIEW = number/context/source requires manual confirmation.",
    "FAIL   = direct inconsistency or missing critical result.",
    "",
    "Important:",
    "A REVIEW result does not automatically mean the manuscript is wrong.",
    "Some values may appear only in figures, legends, or supplementary files."
)


writeLines(
    README,
    file.path(
        OUTDIR,
        "Step12B_README.txt"
    )
)


cat(
    "\nEnd time: ",
    format(Sys.time()),
    "\n",
    sep = ""
)


sink(type = "message")
sink(type = "output")
close(zz)

