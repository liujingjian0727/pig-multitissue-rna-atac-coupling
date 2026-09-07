#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# STEP 07B
# PROMOTER-FIRST GENOMIC CLASS ENRICHMENT
#
# PRIMARY annotation:
#   Pair_specific_class
#
# This describes each ATAC peak relative specifically to its
# UROPA-associated gene.
#
# Integration levels:
#   All UROPA-assigned associations : 168999
#   Concordant HC associations      :   4674
#   Strong associations             :    702
#
# Promoter definition inherited from Step07A:
#   TSS -2000 bp to +500 bp, strand-aware
#
# Priority:
#   Promoter_TSSproximal
#   > 5UTR
#   > 3UTR
#   > Exon
#   > Intron
#   > Distal_to_associated_gene
# ============================================================


# ============================================================
# INPUT
# ============================================================

all_file <- paste0(
    "07_promoter_first_classification/",
    "07A_all_UROPA_assigned_promoter_first.tsv"
)

hc_file <- paste0(
    "05_RNA_ATAC_quantitative_coupling/",
    "05_concordant_HC_peak_gene_correlations.tsv"
)

strong_file <- paste0(
    "05B_same_tissue_permutation_v2/",
    "05B_observed_strong_05_05_07_pairs.tsv"
)

outdir <- "07_promoter_first_enrichment"

dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
)


# ============================================================
# FILE CHECK
# ============================================================

required_files <- c(
    all_file,
    hc_file,
    strong_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing input file(s):\n",
        paste(
            missing_files,
            collapse = "\n"
        )
    )
}


cat("============================================================\n")
cat("STEP 07B: PROMOTER-FIRST GENOMIC CLASS ENRICHMENT\n")
cat("============================================================\n\n")


# ============================================================
# READ
# ============================================================

all <- read.delim(
    all_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

hc <- read.delim(
    hc_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

strong <- read.delim(
    strong_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


cat("Input:\n")
cat("All UROPA assigned :", nrow(all), "\n")
cat("Concordant HC      :", nrow(hc), "\n")
cat("Strong             :", nrow(strong), "\n\n")


# ============================================================
# EXPECTED NUMBER QC
# ============================================================

if (nrow(all) != 168999) {
    warning(
        "Expected 168999 all UROPA-assigned associations."
    )
}

if (nrow(hc) != 4674) {
    warning(
        "Expected 4674 concordant HC associations."
    )
}

if (nrow(strong) != 702) {
    warning(
        "Expected 702 strong associations."
    )
}


# ============================================================
# REQUIRED COLUMN QC
# ============================================================

required_all_cols <- c(
    "SAF_peak_id",
    "gene_id",
    "Pair_specific_class",
    "Genomic_class",
    "Peak_to_UROPA_TSS_min_abs_distance_bp",
    "Genomic_class_gene_matches_UROPA_gene"
)

missing_all_cols <- setdiff(
    required_all_cols,
    colnames(all)
)

if (length(missing_all_cols) > 0) {

    stop(
        "Step07A table missing columns: ",
        paste(
            missing_all_cols,
            collapse = ", "
        )
    )
}


required_subset_cols <- c(
    "SAF_peak_id",
    "gene_id"
)

for (nm in c("HC", "Strong")) {

    z <- if (nm == "HC") hc else strong

    missing_cols <- setdiff(
        required_subset_cols,
        colnames(z)
    )

    if (length(missing_cols) > 0) {

        stop(
            nm,
            " table missing columns: ",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }
}


# ============================================================
# UNIQUE PEAK QC
# ============================================================

if (anyDuplicated(all$SAF_peak_id)) {
    stop(
        "Duplicated SAF_peak_id in Step07A all table."
    )
}

if (anyDuplicated(hc$SAF_peak_id)) {
    stop(
        "Duplicated SAF_peak_id in HC table."
    )
}

if (anyDuplicated(strong$SAF_peak_id)) {
    stop(
        "Duplicated SAF_peak_id in strong table."
    )
}


# ============================================================
# MAP HC / STRONG TO THE ANNOTATED ALL TABLE
# ============================================================

hc_idx <- match(
    hc$SAF_peak_id,
    all$SAF_peak_id
)

strong_idx <- match(
    strong$SAF_peak_id,
    all$SAF_peak_id
)


if (anyNA(hc_idx)) {

    stop(
        sum(is.na(hc_idx)),
        " HC peaks were not found in Step07A."
    )
}

if (anyNA(strong_idx)) {

    stop(
        sum(is.na(strong_idx)),
        " strong peaks were not found in Step07A."
    )
}


ann_hc <- all[
    hc_idx,
    ,
    drop = FALSE
]

ann_strong <- all[
    strong_idx,
    ,
    drop = FALSE
]


# ============================================================
# PEAK-GENE PAIR QC
# ============================================================

hc_gene_match <- (
    ann_hc$gene_id ==
    hc$gene_id
)

strong_gene_match <- (
    ann_strong$gene_id ==
    strong$gene_id
)


if (!all(hc_gene_match)) {

    stop(
        "HC peak IDs matched, but associated gene IDs differed."
    )
}

if (!all(strong_gene_match)) {

    stop(
        "Strong peak IDs matched, but associated gene IDs differed."
    )
}


cat("Subset mapping QC:\n")
cat("HC peak-gene pairs matched     :", sum(hc_gene_match), "/", nrow(hc), "\n")
cat("Strong peak-gene pairs matched :", sum(strong_gene_match), "/", nrow(strong), "\n\n")


# ============================================================
# WRITE FINAL ANNOTATED SUBSETS
# ============================================================

write.table(
    ann_hc,
    file.path(
        outdir,
        "07B_4674_HC_promoter_first.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


write.table(
    ann_strong,
    file.path(
        outdir,
        "07B_702_strong_promoter_first.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# CLASS LEVELS
# ============================================================

pair_levels <- c(
    "Promoter_TSSproximal",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Distal_to_associated_gene",
    "Gene_not_in_GTF"
)


global_levels <- c(
    "Promoter_TSSproximal",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Distal_intergenic"
)


# ============================================================
# CHECK FOR UNEXPECTED CLASSES
# ============================================================

unexpected_pair <- setdiff(
    unique(all$Pair_specific_class),
    pair_levels
)

if (length(unexpected_pair) > 0) {

    stop(
        "Unexpected Pair_specific_class value(s): ",
        paste(
            unexpected_pair,
            collapse = ", "
        )
    )
}


unexpected_global <- setdiff(
    unique(all$Genomic_class),
    global_levels
)

if (length(unexpected_global) > 0) {

    stop(
        "Unexpected Genomic_class value(s): ",
        paste(
            unexpected_global,
            collapse = ", "
        )
    )
}


# ============================================================
# FACTORIZE
# ============================================================

all$Pair_specific_class <- factor(
    all$Pair_specific_class,
    levels = pair_levels
)

ann_hc$Pair_specific_class <- factor(
    ann_hc$Pair_specific_class,
    levels = pair_levels
)

ann_strong$Pair_specific_class <- factor(
    ann_strong$Pair_specific_class,
    levels = pair_levels
)


all$Genomic_class <- factor(
    all$Genomic_class,
    levels = global_levels
)

ann_hc$Genomic_class <- factor(
    ann_hc$Genomic_class,
    levels = global_levels
)

ann_strong$Genomic_class <- factor(
    ann_strong$Genomic_class,
    levels = global_levels
)


# ============================================================
# BASIC COMPOSITION FUNCTION
# ============================================================

make_summary <- function(
    x,
    class_column,
    level_name,
    class_levels
) {

    tab <- table(
        factor(
            x[[class_column]],
            levels = class_levels
        )
    )

    result <- data.frame(
        Genomic_class =
            names(tab),

        Number =
            as.integer(tab),

        Level =
            level_name,

        Total =
            nrow(x),

        stringsAsFactors = FALSE
    )


    result$Percent <- (
        100 *
        result$Number /
        result$Total
    )


    result
}


# ============================================================
# PAIR-SPECIFIC COMPOSITION
# ============================================================

pair_summary <- rbind(

    make_summary(
        all,
        "Pair_specific_class",
        "All_UROPA_assigned",
        pair_levels
    ),

    make_summary(
        ann_hc,
        "Pair_specific_class",
        "Concordant_HC",
        pair_levels
    ),

    make_summary(
        ann_strong,
        "Pair_specific_class",
        "Strong",
        pair_levels
    )

)


write.table(
    pair_summary,
    file.path(
        outdir,
        "07B_pair_specific_class_composition.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# GLOBAL GENOMIC COMPOSITION
# ============================================================

global_summary <- rbind(

    make_summary(
        all,
        "Genomic_class",
        "All_UROPA_assigned",
        global_levels
    ),

    make_summary(
        ann_hc,
        "Genomic_class",
        "Concordant_HC",
        global_levels
    ),

    make_summary(
        ann_strong,
        "Genomic_class",
        "Strong",
        global_levels
    )

)


write.table(
    global_summary,
    file.path(
        outdir,
        "07B_global_genomic_class_composition.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PAIR-SPECIFIC WIDE COUNT TABLE
# ============================================================

pair_count_wide <- xtabs(
    Number ~ Genomic_class + Level,
    data = pair_summary
)


write.table(
    cbind(
        Genomic_class =
            rownames(pair_count_wide),

        as.data.frame.matrix(
            pair_count_wide
        )
    ),
    file.path(
        outdir,
        "07B_pair_specific_counts_wide.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PAIR-SPECIFIC WIDE PERCENT TABLE
# ============================================================

pair_percent_wide <- xtabs(
    Percent ~ Genomic_class + Level,
    data = pair_summary
)


write.table(
    cbind(
        Genomic_class =
            rownames(pair_percent_wide),

        as.data.frame.matrix(
            pair_percent_wide
        )
    ),
    file.path(
        outdir,
        "07B_pair_specific_percent_wide.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# FISHER ENRICHMENT
#
# IMPORTANT:
#
# selected is nested inside background.
#
# Therefore Fisher table is:
#
#                        target class   other classes
# selected
# background-selected
#
# This avoids treating overlapping datasets as independent.
# ============================================================

count_class <- function(
    x,
    class_column,
    target_class
) {

    sum(
        x[[class_column]] ==
        target_class,
        na.rm = TRUE
    )
}


run_enrichment <- function(
    selected,
    background,
    class_column,
    classes,
    comparison_name
) {

    results <- lapply(
        classes,
        function(target_class) {

            selected_class <- count_class(
                selected,
                class_column,
                target_class
            )


            selected_other <- (
                nrow(selected) -
                selected_class
            )


            background_class_total <- count_class(
                background,
                class_column,
                target_class
            )


            complement_n <- (
                nrow(background) -
                nrow(selected)
            )


            complement_class <- (
                background_class_total -
                selected_class
            )


            complement_other <- (
                complement_n -
                complement_class
            )


            if (
                complement_class < 0 ||
                complement_other < 0
            ) {

                stop(
                    "Invalid nested comparison: ",
                    comparison_name,
                    " / ",
                    target_class
                )
            }


            selected_fraction <- (
                selected_class /
                nrow(selected)
            )


            background_fraction <- (
                background_class_total /
                nrow(background)
            )


            fold_enrichment <- if (
                background_fraction > 0
            ) {

                selected_fraction /
                    background_fraction

            } else {

                NA_real_
            }


            # If this class does not exist in the background,
            # there is no meaningful enrichment test.

            if (
                background_class_total == 0
            ) {

                odds_ratio <- NA_real_
                p_value <- 1

            } else {

                fisher_matrix <- matrix(
                    c(
                        selected_class,
                        selected_other,
                        complement_class,
                        complement_other
                    ),
                    nrow = 2,
                    byrow = TRUE
                )


                ft <- fisher.test(
                    fisher_matrix,
                    alternative = "two.sided"
                )


                odds_ratio <- if (
                    length(ft$estimate) > 0
                ) {

                    unname(
                        ft$estimate
                    )

                } else {

                    NA_real_
                }


                p_value <- (
                    ft$p.value
                )
            }


            data.frame(
                Comparison =
                    comparison_name,

                Genomic_class =
                    target_class,

                Selected_N =
                    nrow(selected),

                Selected_class_N =
                    selected_class,

                Selected_percent =
                    100 *
                    selected_fraction,

                Background_N =
                    nrow(background),

                Background_class_N =
                    background_class_total,

                Background_percent =
                    100 *
                    background_fraction,

                Fold_enrichment =
                    fold_enrichment,

                Odds_ratio =
                    odds_ratio,

                P_value =
                    p_value,

                stringsAsFactors = FALSE
            )
        }
    )


    results <- do.call(
        rbind,
        results
    )


    results$FDR <- p.adjust(
        results$P_value,
        method = "BH"
    )


    results
}


# ============================================================
# PRIMARY: PAIR-SPECIFIC ENRICHMENT
# ============================================================

pair_enrichment_hc_all <- run_enrichment(
    ann_hc,
    all,
    "Pair_specific_class",
    pair_levels,
    "Concordant_HC_vs_All"
)


pair_enrichment_strong_all <- run_enrichment(
    ann_strong,
    all,
    "Pair_specific_class",
    pair_levels,
    "Strong_vs_All"
)


pair_enrichment_strong_hc <- run_enrichment(
    ann_strong,
    ann_hc,
    "Pair_specific_class",
    pair_levels,
    "Strong_vs_Concordant_HC"
)


pair_enrichment <- rbind(
    pair_enrichment_hc_all,
    pair_enrichment_strong_all,
    pair_enrichment_strong_hc
)


write.table(
    pair_enrichment,
    file.path(
        outdir,
        "07B_pair_specific_enrichment_tests.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# SECONDARY: GLOBAL GENOMIC ENRICHMENT
# ============================================================

global_enrichment <- rbind(

    run_enrichment(
        ann_hc,
        all,
        "Genomic_class",
        global_levels,
        "Concordant_HC_vs_All"
    ),

    run_enrichment(
        ann_strong,
        all,
        "Genomic_class",
        global_levels,
        "Strong_vs_All"
    ),

    run_enrichment(
        ann_strong,
        ann_hc,
        "Genomic_class",
        global_levels,
        "Strong_vs_Concordant_HC"
    )

)


write.table(
    global_enrichment,
    file.path(
        outdir,
        "07B_global_genomic_enrichment_tests.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# PROMOTER vs NON-PROMOTER SUMMARY
# ============================================================

make_promoter_summary <- function(
    x,
    level_name
) {

    promoter <- (
        x$Pair_specific_class ==
        "Promoter_TSSproximal"
    )


    promoter_n <- sum(
        promoter,
        na.rm = TRUE
    )


    data.frame(
        Level =
            level_name,

        N =
            nrow(x),

        Promoter_N =
            promoter_n,

        Promoter_percent =
            100 *
            promoter_n /
            nrow(x),

        Nonpromoter_N =
            nrow(x) -
            promoter_n,

        Nonpromoter_percent =
            100 *
            (
                nrow(x) -
                promoter_n
            ) /
            nrow(x),

        stringsAsFactors =
            FALSE
    )
}


promoter_summary <- rbind(

    make_promoter_summary(
        all,
        "All_UROPA_assigned"
    ),

    make_promoter_summary(
        ann_hc,
        "Concordant_HC"
    ),

    make_promoter_summary(
        ann_strong,
        "Strong"
    )

)


write.table(
    promoter_summary,
    file.path(
        outdir,
        "07B_promoter_vs_nonpromoter_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# TRUE REPRESENTATIVE-TSS DISTANCE SUMMARY
# ============================================================

summarize_tss_distance <- function(
    x,
    level_name
) {

    d <- suppressWarnings(
        as.numeric(
            x$Peak_to_UROPA_TSS_min_abs_distance_bp
        )
    )


    d <- d[
        is.finite(d)
    ]


    if (length(d) == 0) {

        return(
            data.frame(
                Level =
                    level_name,

                N_with_TSS_distance =
                    0,

                Median_TSS_distance_bp =
                    NA,

                Q25_TSS_distance_bp =
                    NA,

                Q75_TSS_distance_bp =
                    NA,

                Maximum_TSS_distance_bp =
                    NA
            )
        )
    }


    data.frame(
        Level =
            level_name,

        N_with_TSS_distance =
            length(d),

        Median_TSS_distance_bp =
            median(d),

        Q25_TSS_distance_bp =
            as.numeric(
                quantile(
                    d,
                    0.25,
                    names = FALSE
                )
            ),

        Q75_TSS_distance_bp =
            as.numeric(
                quantile(
                    d,
                    0.75,
                    names = FALSE
                )
            ),

        Maximum_TSS_distance_bp =
            max(d)
    )
}


tss_summary <- rbind(

    summarize_tss_distance(
        all,
        "All_UROPA_assigned"
    ),

    summarize_tss_distance(
        ann_hc,
        "Concordant_HC"
    ),

    summarize_tss_distance(
        ann_strong,
        "Strong"
    )

)


write.table(
    tss_summary,
    file.path(
        outdir,
        "07B_TSS_distance_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 702 STRONG:
# TISSUE x PAIR-SPECIFIC CLASS
# ============================================================

tissues <- c(
    "Adipose",
    "Cerebellum",
    "Cortex",
    "Hypothalamus",
    "Liver",
    "Lung",
    "Muscle",
    "Spleen"
)


if (
    !"ATAC_Max_tissue" %in%
    colnames(ann_strong)
) {

    stop(
        "ATAC_Max_tissue column is missing from Step07A table."
    )
}


tissue_class <- as.data.frame(
    table(
        factor(
            ann_strong$ATAC_Max_tissue,
            levels = tissues
        ),

        factor(
            ann_strong$Pair_specific_class,
            levels = pair_levels
        )
    )
)


colnames(
    tissue_class
) <- c(
    "Tissue",
    "Genomic_class",
    "Number"
)


tissue_totals <- tapply(
    tissue_class$Number,
    tissue_class$Tissue,
    sum
)


tissue_class$Percent_within_tissue <- NA_real_


for (
    i in seq_len(
        nrow(tissue_class)
    )
) {

    tissue_name <- as.character(
        tissue_class$Tissue[i]
    )


    denominator <- tissue_totals[
        tissue_name
    ]


    if (
        !is.na(denominator) &&
        denominator > 0
    ) {

        tissue_class$Percent_within_tissue[i] <-
            100 *
            tissue_class$Number[i] /
            denominator
    }
}


write.table(
    tissue_class,
    file.path(
        outdir,
        "07B_702_tissue_by_pair_specific_class.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 702 STRONG:
# CLASS-SPECIFIC COUNTS + UNIQUE GENES
# ============================================================

strong_class_gene_summary <- do.call(
    rbind,
    lapply(
        pair_levels,
        function(cc) {

            z <- ann_strong[
                ann_strong$Pair_specific_class ==
                    cc,
                ,
                drop = FALSE
            ]


            data.frame(
                Genomic_class =
                    cc,

                N_peak_gene_pairs =
                    nrow(z),

                Unique_peaks =
                    length(
                        unique(
                            z$SAF_peak_id
                        )
                    ),

                Unique_genes =
                    length(
                        unique(
                            z$gene_id
                        )
                    ),

                stringsAsFactors =
                    FALSE
            )
        }
    )
)


write.table(
    strong_class_gene_summary,
    file.path(
        outdir,
        "07B_702_class_peak_gene_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# GLOBAL CLASS vs UROPA-GENE MATCH QC
# ============================================================

parse_boolean <- function(x) {

    as.character(x) %in%
        c(
            "TRUE",
            "True",
            "true",
            "1"
        )
}


make_gene_match_qc <- function(
    x,
    level_name
) {

    raw <- as.character(
        x$Genomic_class_gene_matches_UROPA_gene
    )


    valid <- (
        !is.na(raw) &
        raw != "" &
        raw != "NA"
    )


    flag <- parse_boolean(
        raw
    )


    n_valid <- sum(
        valid
    )


    n_match <- sum(
        flag &
        valid
    )


    data.frame(
        Level =
            level_name,

        N =
            nrow(x),

        Global_nonintergenic =
            sum(
                x$Genomic_class !=
                    "Distal_intergenic",
                na.rm = TRUE
            ),

        Valid_gene_match_flag =
            n_valid,

        Class_gene_matches_UROPA_gene =
            n_match,

        Match_percent =
            ifelse(
                n_valid > 0,
                100 *
                    n_match /
                    n_valid,
                NA_real_
            ),

        stringsAsFactors =
            FALSE
    )
}


gene_match_qc <- rbind(

    make_gene_match_qc(
        all,
        "All_UROPA_assigned"
    ),

    make_gene_match_qc(
        ann_hc,
        "Concordant_HC"
    ),

    make_gene_match_qc(
        ann_strong,
        "Strong"
    )

)


write.table(
    gene_match_qc,
    file.path(
        outdir,
        "07B_global_class_gene_match_QC.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# PLOT 1
# PROMOTER-FIRST COMPOSITION
# ============================================================

pair_summary$Level <- factor(
    pair_summary$Level,
    levels = c(
        "All_UROPA_assigned",
        "Concordant_HC",
        "Strong"
    )
)


pair_summary$Genomic_class <- factor(
    pair_summary$Genomic_class,
    levels = pair_levels
)


p1 <- ggplot(
    pair_summary,
    aes(
        x = Level,
        y = Percent,
        fill = Genomic_class
    )
) +
    geom_col() +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title =
            "Promoter-first genomic classification",

        subtitle =
            paste0(
                "Pair-specific genomic context relative ",
                "to the UROPA-associated gene"
            ),

        x = NULL,

        y = "Percentage",

        fill = "Genomic class"
    )


ggsave(
    file.path(
        outdir,
        "07B_pair_specific_class_composition.pdf"
    ),
    p1,
    width = 9,
    height = 6
)


# ============================================================
# PLOT 2
# FOLD ENRICHMENT
# ============================================================

plot_enrichment <- pair_enrichment[
    pair_enrichment$Comparison %in%
        c(
            "Concordant_HC_vs_All",
            "Strong_vs_All"
        ) &
        is.finite(
            pair_enrichment$Fold_enrichment
        ),
    ,
    drop = FALSE
]


plot_enrichment$Genomic_class <- factor(
    plot_enrichment$Genomic_class,
    levels = pair_levels
)


p2 <- ggplot(
    plot_enrichment,
    aes(
        x = Genomic_class,
        y = Fold_enrichment,
        shape = Comparison,
        group = Comparison
    )
) +
    geom_hline(
        yintercept = 1,
        linetype = 2
    ) +
    geom_point(
        size = 3
    ) +
    theme_bw(
        base_size = 11
    ) +
    theme(
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    ) +
    labs(
        title =
            "Promoter-first genomic class enrichment",

        subtitle =
            "Pair-specific classification",

        x = NULL,

        y = "Fold enrichment"
    )


ggsave(
    file.path(
        outdir,
        "07B_pair_specific_fold_enrichment.pdf"
    ),
    p2,
    width = 9,
    height = 5.5
)


# ============================================================
# PLOT 3
# PROMOTER PROPORTION
# ============================================================

promoter_summary$Level <- factor(
    promoter_summary$Level,
    levels = c(
        "All_UROPA_assigned",
        "Concordant_HC",
        "Strong"
    )
)


p3 <- ggplot(
    promoter_summary,
    aes(
        x = Level,
        y = Promoter_percent
    )
) +
    geom_col() +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title =
            "Promoter-associated peak-gene associations",

        subtitle =
            "Promoter = representative TSS -2 kb to +500 bp",

        x = NULL,

        y = "Promoter-associated pairs (%)"
    )


ggsave(
    file.path(
        outdir,
        "07B_promoter_percentage.pdf"
    ),
    p3,
    width = 6.5,
    height = 5
)


# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("STEP 07B COMPLETED\n")
cat("============================================================\n\n")


cat("PAIR-SPECIFIC PROMOTER-FIRST COMPOSITION:\n\n")

print(
    pair_summary,
    row.names = FALSE
)


cat("\nPROMOTER VS NON-PROMOTER:\n\n")

print(
    promoter_summary,
    row.names = FALSE
)


cat("\nPAIR-SPECIFIC ENRICHMENT:\n\n")

print(
    pair_enrichment[
        ,
        c(
            "Comparison",
            "Genomic_class",
            "Selected_percent",
            "Background_percent",
            "Fold_enrichment",
            "Odds_ratio",
            "P_value",
            "FDR"
        )
    ],
    row.names = FALSE
)


cat("\nTRUE TSS DISTANCE SUMMARY:\n\n")

print(
    tss_summary,
    row.names = FALSE
)


cat("\n702 STRONG PAIRS BY GENOMIC CLASS:\n\n")

print(
    strong_class_gene_summary,
    row.names = FALSE
)


cat("\nGLOBAL CLASS / UROPA-GENE MATCH QC:\n\n")

print(
    gene_match_qc,
    row.names = FALSE
)


cat("\nIMPORTANT:\n")

cat(
    "Primary annotation = Pair_specific_class.\n"
)

cat(
    "Promoter_TSSproximal = representative TSS -2000 bp to +500 bp.\n"
)

cat(
    "Distal_to_associated_gene means the peak does not overlap the\n",
    "promoter/UTR/exon/intron of its associated UROPA gene.\n",
    sep = ""
)

cat(
    "It does NOT necessarily mean genome-wide intergenic.\n"
)

cat(
    "Fisher enrichment tests use selected vs background-minus-selected.\n"
)


cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("============================================================\n")


sessionInfo()

