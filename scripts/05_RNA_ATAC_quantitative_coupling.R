#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# ============================================================
# STEP 05
# RNA-ATAC quantitative coupling across 8 tissues
#
# Main statistics:
#   rho_P348  : 8 tissues in P348
#   rho_P350  : 8 tissues in P350
#   rho_mean  : 8 tissue means across P348/P350
#   rho_all16 : all 16 matched samples (DESCRIPTIVE ONLY)
#
# Correlation:
#   Spearman
#
# Input values:
#   LINEAR TMM.TPM
#
# IMPORTANT:
#   rho_all16 is NOT treated as 16 independent biological
#   replicates.
# ============================================================


# ============================================================
# INPUT
# ============================================================

master_file <- paste0(
    "04_RNA_ATAC_peak_gene_master/",
    "04_peak_gene_RNA_ATAC_master.tsv"
)

atac_tpm_file <- "pig.TMM.TPM.matrix"
rna_tpm_file  <- "sus.TMM.TPM.matrix"

metadata_file <- "sample_pair_metadata.tsv"

outdir <- "05_RNA_ATAC_quantitative_coupling"

dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

read_expression_matrix <- function(file) {

    x <- read.delim(
        file,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
    )

    x <- as.matrix(x)

    storage.mode(x) <- "numeric"

    return(x)
}


as_flag <- function(x) {

    x %in% c(
        TRUE,
        "TRUE",
        "True",
        "true",
        1,
        "1"
    )
}


# ------------------------------------------------------------
# Fast row ranking
# ------------------------------------------------------------

rank_rows <- function(x) {

    if (
        requireNamespace(
            "matrixStats",
            quietly = TRUE
        )
    ) {

        return(
            matrixStats::rowRanks(
                x,
                ties.method = "average",
                useNames = FALSE
            )
        )

    } else {

        message(
            "Package matrixStats is not available. ",
            "Using base R apply(); this will be slower."
        )

        return(
            t(
                apply(
                    x,
                    1,
                    rank,
                    ties.method = "average"
                )
            )
        )
    }
}


# ------------------------------------------------------------
# Vectorized row-wise Spearman correlation
# ------------------------------------------------------------

row_spearman <- function(x, y, label = "") {

    if (!all(dim(x) == dim(y))) {
        stop(
            "Dimension mismatch in row_spearman: ",
            label
        )
    }

    cat(
        "  Calculating ",
        label,
        " (",
        nrow(x),
        " pairs x ",
        ncol(x),
        " observations)...\n",
        sep = ""
    )

    rx <- rank_rows(x)
    ry <- rank_rows(y)

    mx <- rowMeans(rx)
    my <- rowMeans(ry)

    xc <- rx - mx
    yc <- ry - my

    numerator <- rowSums(
        xc * yc
    )

    denominator <- sqrt(
        rowSums(xc^2) *
        rowSums(yc^2)
    )

    rho <- numerator / denominator

    rho[
        denominator == 0 |
        !is.finite(denominator)
    ] <- NA_real_

    rm(
        rx,
        ry,
        xc,
        yc,
        numerator,
        denominator
    )

    invisible(gc())

    return(rho)
}


safe_median <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    median(x)
}


safe_quantile <- function(x, p) {

    x <- x[is.finite(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    as.numeric(
        quantile(
            x,
            probs = p,
            names = FALSE
        )
    )
}


# ============================================================
# READ INPUT
# ============================================================

cat(
    "============================================================\n"
)
cat(
    "STEP 05: RNA-ATAC QUANTITATIVE COUPLING\n"
)
cat(
    "============================================================\n\n"
)

cat("Reading master table...\n")

master <- read.delim(
    master_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("Reading metadata...\n")

meta <- read.delim(
    metadata_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("Reading ATAC TMM.TPM...\n")

atac <- read_expression_matrix(
    atac_tpm_file
)

cat("Reading RNA TMM.TPM...\n")

rna <- read_expression_matrix(
    rna_tpm_file
)


cat("\nInput dimensions:\n")

cat(
    "Master:",
    nrow(master),
    "rows\n"
)

cat(
    "ATAC:",
    nrow(atac),
    "features x",
    ncol(atac),
    "samples\n"
)

cat(
    "RNA :",
    nrow(rna),
    "genes x",
    ncol(rna),
    "samples\n"
)

cat(
    "Metadata:",
    nrow(meta),
    "paired sample rows\n\n"
)


# ============================================================
# BASIC QC
# ============================================================

required_meta <- c(
    "Tissue",
    "Animal",
    "RNA_SRR",
    "ATAC_SRR"
)

if (!all(required_meta %in% colnames(meta))) {

    stop(
        "Required metadata columns are missing."
    )
}


required_master <- c(
    "SAF_peak_id",
    "BED_peak_id",
    "gene_id",
    "Has_gene_annotation",
    "Gene_present_in_RNA",
    "ATAC_Max_tissue",
    "RNA_Max_tissue",
    "ATAC_High_confidence_tissue_specific",
    "RNA_High_confidence_tissue_specific",
    "Concordant_HC_tissue_specific_pair"
)

if (!all(required_master %in% colnames(master))) {

    missing_cols <- setdiff(
        required_master,
        colnames(master)
    )

    stop(
        "Required master columns missing: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


if (anyDuplicated(rownames(atac))) {
    stop("Duplicated ATAC feature IDs.")
}

if (anyDuplicated(rownames(rna))) {
    stop("Duplicated RNA gene IDs.")
}


if (
    !all(
        meta$ATAC_SRR %in% colnames(atac)
    )
) {

    stop(
        "Some ATAC metadata samples are absent from ATAC matrix."
    )
}


if (
    !all(
        meta$RNA_SRR %in% colnames(rna)
    )
) {

    stop(
        "Some RNA metadata samples are absent from RNA matrix."
    )
}


if (anyNA(atac)) {
    stop("ATAC TMM.TPM matrix contains NA.")
}

if (anyNA(rna)) {
    stop("RNA TMM.TPM matrix contains NA.")
}


# ============================================================
# FIX TISSUE / ANIMAL ORDER
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

animals <- c(
    "P348",
    "P350"
)


meta_key <- paste(
    meta$Tissue,
    meta$Animal,
    sep = "__"
)


expected_key <- as.vector(
    outer(
        tissues,
        animals,
        paste,
        sep = "__"
    )
)


if (
    !all(
        expected_key %in% meta_key
    )
) {

    stop(
        "Not all 8 tissue x 2 animal combinations ",
        "are present in metadata."
    )
}


# ------------------------------------------------------------
# Tissue order for P348
# ------------------------------------------------------------

idx_P348 <- match(
    paste(
        tissues,
        "P348",
        sep = "__"
    ),
    meta_key
)


# ------------------------------------------------------------
# Tissue order for P350
# ------------------------------------------------------------

idx_P350 <- match(
    paste(
        tissues,
        "P350",
        sep = "__"
    ),
    meta_key
)


if (
    anyNA(idx_P348) ||
    anyNA(idx_P350)
) {

    stop(
        "Failed to construct animal-specific tissue order."
    )
}


cat(
    "PASS: metadata contains complete ",
    "2 animals x 8 tissues paired design.\n"
)


# ============================================================
# SELECT ALL UROPA-ASSIGNED ASSOCIATIONS
# ============================================================

assigned <- (
    as_flag(
        master$Has_gene_annotation
    ) &
    as_flag(
        master$Gene_present_in_RNA
    ) &
    !is.na(master$gene_id) &
    master$gene_id != ""
)


assoc <- master[
    assigned,
    ,
    drop = FALSE
]


cat(
    "\nUROPA-assigned peak-gene associations:",
    nrow(assoc),
    "\n"
)


# Should equal current Step04 result
if (nrow(assoc) != 168999) {

    warning(
        "Expected 168999 associations from Step04, but found ",
        nrow(assoc),
        ". Analysis will continue using the actual count."
    )
}


# ============================================================
# MAP PEAKS AND GENES TO TMM.TPM MATRICES
# ============================================================

atac_idx <- match(
    assoc$SAF_peak_id,
    rownames(atac)
)

rna_idx <- match(
    assoc$gene_id,
    rownames(rna)
)


if (anyNA(atac_idx)) {

    stop(
        sum(is.na(atac_idx)),
        " master peaks are missing from ATAC TMM.TPM matrix."
    )
}


if (anyNA(rna_idx)) {

    stop(
        sum(is.na(rna_idx)),
        " master genes are missing from RNA TMM.TPM matrix."
    )
}


cat(
    "PASS: all assigned peaks and genes matched ",
    "their TMM.TPM matrices.\n"
)


# ============================================================
# BUILD MATCHED EXPRESSION / ACCESSIBILITY MATRICES
#
# IMPORTANT:
# RNA and ATAC rows here correspond exactly to the same
# peak-gene association.
# Duplicate gene IDs are deliberately retained.
# ============================================================

cat(
    "\nConstructing 168999 matched peak-gene profiles...\n"
)


atac_assoc <- atac[
    atac_idx,
    meta$ATAC_SRR,
    drop = FALSE
]


rna_assoc <- rna[
    rna_idx,
    meta$RNA_SRR,
    drop = FALSE
]


# ============================================================
# P348: 8 tissues
# ============================================================

ATAC_P348 <- atac_assoc[
    ,
    idx_P348,
    drop = FALSE
]

RNA_P348 <- rna_assoc[
    ,
    idx_P348,
    drop = FALSE
]


colnames(ATAC_P348) <- tissues
colnames(RNA_P348)  <- tissues


# ============================================================
# P350: 8 tissues
# ============================================================

ATAC_P350 <- atac_assoc[
    ,
    idx_P350,
    drop = FALSE
]

RNA_P350 <- rna_assoc[
    ,
    idx_P350,
    drop = FALSE
]


colnames(ATAC_P350) <- tissues
colnames(RNA_P350)  <- tissues


# ============================================================
# TISSUE MEAN: mean of P348 and P350 in LINEAR TMM.TPM
# ============================================================

ATAC_mean <- (
    ATAC_P348 +
    ATAC_P350
) / 2


RNA_mean <- (
    RNA_P348 +
    RNA_P350
) / 2


colnames(ATAC_mean) <- tissues
colnames(RNA_mean)  <- tissues


# ============================================================
# ALL 16 MATCHED SAMPLES
#
# This is only descriptive.
# It is NOT interpreted as n=16 independent animals.
# ============================================================

ATAC_all16 <- cbind(
    ATAC_P348,
    ATAC_P350
)

RNA_all16 <- cbind(
    RNA_P348,
    RNA_P350
)


# ============================================================
# CALCULATE SPEARMAN CORRELATIONS
# ============================================================

cat(
    "\n============================================================\n"
)
cat(
    "Calculating row-wise Spearman correlations\n"
)
cat(
    "============================================================\n"
)


rho_P348 <- row_spearman(
    ATAC_P348,
    RNA_P348,
    "rho_P348"
)


rho_P350 <- row_spearman(
    ATAC_P350,
    RNA_P350,
    "rho_P350"
)


rho_mean <- row_spearman(
    ATAC_mean,
    RNA_mean,
    "rho_mean"
)


rho_all16 <- row_spearman(
    ATAC_all16,
    RNA_all16,
    "rho_all16"
)


# ============================================================
# BUILD OUTPUT TABLE
# ============================================================

key_columns <- c(
    "SAF_peak_id",
    "BED_peak_id",
    "peak_chr",
    "BED_start",
    "BED_end",

    "gene_id",

    "UROPA_feature",
    "UROPA_relative_location",
    "UROPA_distance",

    "ATAC_Tau",
    "ATAC_Max_SPM",
    "ATAC_Max_tissue",
    "ATAC_Max_mean_TMM_TPM",
    "ATAC_Dominance_ratio",
    "ATAC_High_confidence_tissue_specific",

    "RNA_Tau",
    "RNA_Max_SPM",
    "RNA_Max_tissue",
    "RNA_Max_mean_TMM_TPM",
    "RNA_Dominance_ratio",
    "RNA_High_confidence_tissue_specific",

    "ATAC_RNA_same_max_tissue",
    "Concordant_HC_tissue_specific_pair"
)


key_columns <- intersect(
    key_columns,
    colnames(assoc)
)


result <- assoc[
    ,
    key_columns,
    drop = FALSE
]


result$rho_P348 <- rho_P348
result$rho_P350 <- rho_P350
result$rho_mean <- rho_mean
result$rho_all16 <- rho_all16


# ============================================================
# DERIVED CORRELATION CONSISTENCY METRICS
# ============================================================

both_animal_valid <- (
    is.finite(result$rho_P348) &
    is.finite(result$rho_P350)
)


all_three_valid <- (
    both_animal_valid &
    is.finite(result$rho_mean)
)


result$Both_animals_valid <- both_animal_valid

result$All_three_primary_rho_valid <-
    all_three_valid


result$Same_sign_between_animals <- (
    both_animal_valid &
    (
        result$rho_P348 *
        result$rho_P350
    ) > 0
)


result$Positive_in_both_animals <- (
    both_animal_valid &
    result$rho_P348 > 0 &
    result$rho_P350 > 0
)


result$Positive_P348_P350_and_mean <- (
    all_three_valid &
    result$rho_P348 > 0 &
    result$rho_P350 > 0 &
    result$rho_mean > 0
)


result$Mean_animal_rho <- NA_real_

result$Mean_animal_rho[
    both_animal_valid
] <- (
    result$rho_P348[
        both_animal_valid
    ] +
    result$rho_P350[
        both_animal_valid
    ]
) / 2


result$Min_animal_rho <- NA_real_

result$Min_animal_rho[
    both_animal_valid
] <- pmin(
    result$rho_P348[
        both_animal_valid
    ],
    result$rho_P350[
        both_animal_valid
    ]
)


result$Animal_rho_absolute_difference <-
    NA_real_

result$Animal_rho_absolute_difference[
    both_animal_valid
] <- abs(
    result$rho_P348[
        both_animal_valid
    ] -
    result$rho_P350[
        both_animal_valid
    ]
)


# ============================================================
# WRITE ALL 168999 ASSOCIATIONS
# ============================================================

all_file <- file.path(
    outdir,
    "05_all_UROPA_assigned_peak_gene_correlations.tsv"
)


write.table(
    result,
    all_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# EXTRACT THE PRE-DEFINED 4674 CONCORDANT HC ASSOCIATIONS
#
# IMPORTANT:
# Do not recalculate Tau/SPM thresholds here.
# Membership comes directly from Step04.
# ============================================================

hc_flag <- as_flag(
    result$Concordant_HC_tissue_specific_pair
)


hc <- result[
    hc_flag,
    ,
    drop = FALSE
]


cat(
    "\nConcordant HC peak-gene associations:",
    nrow(hc),
    "\n"
)


if (nrow(hc) != 4674) {

    warning(
        "Expected 4674 concordant HC associations ",
        "from Step04, but found ",
        nrow(hc),
        "."
    )
}


hc_file <- file.path(
    outdir,
    "05_concordant_HC_peak_gene_correlations.tsv"
)


write.table(
    hc,
    hc_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# INVALID / CONSTANT PROFILE QC
# ============================================================

invalid <- result[
    !result$All_three_primary_rho_valid,
    ,
    drop = FALSE
]


write.table(
    invalid,
    file.path(
        outdir,
        "05_associations_with_invalid_primary_rho.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# SUMMARY FUNCTION
# ============================================================

make_summary <- function(x, dataset_name) {

    valid348 <- is.finite(
        x$rho_P348
    )

    valid350 <- is.finite(
        x$rho_P350
    )

    validmean <- is.finite(
        x$rho_mean
    )

    valid16 <- is.finite(
        x$rho_all16
    )

    both <- (
        valid348 &
        valid350
    )

    all3 <- (
        both &
        validmean
    )

    agreement <- suppressWarnings(
        cor(
            x$rho_P348,
            x$rho_P350,
            method = "spearman",
            use = "complete.obs"
        )
    )

    data.frame(
        Dataset = dataset_name,

        N_associations =
            nrow(x),

        Valid_rho_P348 =
            sum(valid348),

        Valid_rho_P350 =
            sum(valid350),

        Valid_rho_mean =
            sum(validmean),

        Valid_rho_all16 =
            sum(valid16),

        Valid_all_three =
            sum(all3),

        Median_rho_P348 =
            safe_median(
                x$rho_P348
            ),

        Median_rho_P350 =
            safe_median(
                x$rho_P350
            ),

        Median_rho_mean =
            safe_median(
                x$rho_mean
            ),

        Median_rho_all16 =
            safe_median(
                x$rho_all16
            ),

        Q25_rho_mean =
            safe_quantile(
                x$rho_mean,
                0.25
            ),

        Q75_rho_mean =
            safe_quantile(
                x$rho_mean,
                0.75
            ),

        Positive_both_animals =
            sum(
                both &
                x$rho_P348 > 0 &
                x$rho_P350 > 0
            ),

        Positive_P348_P350_mean =
            sum(
                all3 &
                x$rho_P348 > 0 &
                x$rho_P350 > 0 &
                x$rho_mean > 0
            ),

        Same_sign_animals =
            sum(
                both &
                (
                    x$rho_P348 *
                    x$rho_P350
                ) > 0
            ),

        Spearman_agreement_rho_P348_vs_P350 =
            agreement,

        stringsAsFactors = FALSE
    )
}


summary_all <- make_summary(
    result,
    "All_UROPA_assigned"
)


summary_hc <- make_summary(
    hc,
    "Concordant_HC"
)


summary_table <- rbind(
    summary_all,
    summary_hc
)


write.table(
    summary_table,
    file.path(
        outdir,
        "05_correlation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# EXPLORATORY EFFECT-SIZE SENSITIVITY
#
# These are NOT final thresholds.
# They are only provided to inspect the distribution before
# choosing a final high-confidence coupling criterion.
# ============================================================

make_sensitivity <- function(x, dataset_name) {

    configs <- data.frame(
        Criterion = c(
            "Positive_all_three",
            "Animal_ge_0.3_mean_ge_0.5",
            "Animal_ge_0.5_mean_ge_0.7",
            "Animal_ge_0.7_mean_ge_0.8"
        ),

        Animal_threshold = c(
            0,
            0.3,
            0.5,
            0.7
        ),

        Mean_threshold = c(
            0,
            0.5,
            0.7,
            0.8
        ),

        stringsAsFactors = FALSE
    )


    out <- lapply(
        seq_len(
            nrow(configs)
        ),
        function(i) {

            a <- configs$Animal_threshold[i]
            m <- configs$Mean_threshold[i]


            valid <- (
                is.finite(
                    x$rho_P348
                ) &
                is.finite(
                    x$rho_P350
                ) &
                is.finite(
                    x$rho_mean
                )
            )


            if (a == 0 && m == 0) {

                pass <- (
                    valid &
                    x$rho_P348 > 0 &
                    x$rho_P350 > 0 &
                    x$rho_mean > 0
                )

            } else {

                pass <- (
                    valid &
                    x$rho_P348 >= a &
                    x$rho_P350 >= a &
                    x$rho_mean >= m
                )
            }


            data.frame(
                Dataset =
                    dataset_name,

                Criterion =
                    configs$Criterion[i],

                Animal_rho_threshold =
                    a,

                Tissue_mean_rho_threshold =
                    m,

                Number =
                    sum(pass),

                Fraction_percent =
                    round(
                        100 *
                        sum(pass) /
                        nrow(x),
                        3
                    ),

                Unique_genes =
                    length(
                        unique(
                            x$gene_id[
                                pass &
                                !is.na(
                                    x$gene_id
                                )
                            ]
                        )
                    ),

                stringsAsFactors = FALSE
            )
        }
    )


    do.call(
        rbind,
        out
    )
}


sensitivity <- rbind(
    make_sensitivity(
        result,
        "All_UROPA_assigned"
    ),

    make_sensitivity(
        hc,
        "Concordant_HC"
    )
)


write.table(
    sensitivity,
    file.path(
        outdir,
        "05_correlation_threshold_sensitivity.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# CONCORDANT HC: TISSUE-LEVEL SUMMARY
# ============================================================

hc_tissue_summary <- do.call(
    rbind,
    lapply(
        tissues,
        function(tissue) {

            z <- hc[
                hc$ATAC_Max_tissue ==
                    tissue,
                ,
                drop = FALSE
            ]


            valid <- (
                is.finite(
                    z$rho_P348
                ) &
                is.finite(
                    z$rho_P350
                ) &
                is.finite(
                    z$rho_mean
                )
            )


            positive <- (
                valid &
                z$rho_P348 > 0 &
                z$rho_P350 > 0 &
                z$rho_mean > 0
            )


            moderate <- (
                valid &
                z$rho_P348 >= 0.5 &
                z$rho_P350 >= 0.5 &
                z$rho_mean >= 0.7
            )


            data.frame(
                Tissue =
                    tissue,

                HC_peak_gene_pairs =
                    nrow(z),

                Unique_genes =
                    length(
                        unique(
                            z$gene_id
                        )
                    ),

                Valid_all_three_rho =
                    sum(valid),

                Median_rho_P348 =
                    safe_median(
                        z$rho_P348
                    ),

                Median_rho_P350 =
                    safe_median(
                        z$rho_P350
                    ),

                Median_rho_mean =
                    safe_median(
                        z$rho_mean
                    ),

                Median_rho_all16 =
                    safe_median(
                        z$rho_all16
                    ),

                Positive_all_three =
                    sum(positive),

                Positive_all_three_percent =
                    ifelse(
                        nrow(z) > 0,
                        round(
                            100 *
                            sum(positive) /
                            nrow(z),
                            3
                        ),
                        NA
                    ),

                Exploratory_rho05_05_07 =
                    sum(moderate),

                Exploratory_rho05_05_07_unique_genes =
                    length(
                        unique(
                            z$gene_id[
                                moderate
                            ]
                        )
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)


write.table(
    hc_tissue_summary,
    file.path(
        outdir,
        "05_concordant_HC_tissue_coupling_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# OPTIONAL CORRELATION DISTRIBUTION PLOTS
# ============================================================

if (
    requireNamespace(
        "ggplot2",
        quietly = TRUE
    )
) {

    suppressPackageStartupMessages(
        library(ggplot2)
    )


    make_long <- function(x, dataset) {

        rbind(
            data.frame(
                Dataset = dataset,
                Metric = "rho_P348",
                rho = x$rho_P348
            ),

            data.frame(
                Dataset = dataset,
                Metric = "rho_P350",
                rho = x$rho_P350
            ),

            data.frame(
                Dataset = dataset,
                Metric = "rho_mean",
                rho = x$rho_mean
            ),

            data.frame(
                Dataset = dataset,
                Metric = "rho_all16",
                rho = x$rho_all16
            )
        )
    }


    plot_all <- make_long(
        result,
        "All assigned"
    )

    p1 <- ggplot(
        plot_all[
            is.finite(
                plot_all$rho
            ),
        ],
        aes(
            x = rho
        )
    ) +
        geom_histogram(
            bins = 50
        ) +
        facet_wrap(
            ~ Metric,
            ncol = 2
        ) +
        theme_bw(
            base_size = 11
        ) +
        labs(
            title = paste0(
                "RNA-ATAC coupling: ",
                "all UROPA-assigned peak-gene associations"
            ),
            x = "Spearman rho",
            y = "Number of associations"
        )


    ggsave(
        file.path(
            outdir,
            "05_all_assigned_rho_distributions.pdf"
        ),
        p1,
        width = 9,
        height = 7
    )


    plot_hc <- make_long(
        hc,
        "Concordant HC"
    )

    p2 <- ggplot(
        plot_hc[
            is.finite(
                plot_hc$rho
            ),
        ],
        aes(
            x = rho
        )
    ) +
        geom_histogram(
            bins = 40
        ) +
        facet_wrap(
            ~ Metric,
            ncol = 2
        ) +
        theme_bw(
            base_size = 11
        ) +
        labs(
            title = paste0(
                "RNA-ATAC coupling: ",
                "concordant tissue-specific peak-gene associations"
            ),
            x = "Spearman rho",
            y = "Number of associations"
        )


    ggsave(
        file.path(
            outdir,
            "05_concordant_HC_rho_distributions.pdf"
        ),
        p2,
        width = 9,
        height = 7
    )


    scatter_data <- hc[
        is.finite(
            hc$rho_P348
        ) &
        is.finite(
            hc$rho_P350
        ),
        ,
        drop = FALSE
    ]


    p3 <- ggplot(
        scatter_data,
        aes(
            x = rho_P348,
            y = rho_P350
        )
    ) +
        geom_point(
            alpha = 0.35,
            size = 1
        ) +
        geom_hline(
            yintercept = 0,
            linetype = 2
        ) +
        geom_vline(
            xintercept = 0,
            linetype = 2
        ) +
        coord_equal(
            xlim = c(-1, 1),
            ylim = c(-1, 1)
        ) +
        theme_bw(
            base_size = 11
        ) +
        labs(
            title = paste0(
                "Cross-animal consistency of ",
                "RNA-ATAC coupling"
            ),
            x = expression(
                rho[P348]
            ),
            y = expression(
                rho[P350]
            )
        )


    ggsave(
        file.path(
            outdir,
            "05_concordant_HC_rho_P348_vs_P350.pdf"
        ),
        p3,
        width = 6.5,
        height = 6
    )
}


# ============================================================
# FINAL REPORT
# ============================================================

cat(
    "\n============================================================\n"
)
cat(
    "STEP 05 COMPLETED\n"
)
cat(
    "============================================================\n\n"
)


cat(
    "All assigned peak-gene associations:",
    nrow(result),
    "\n"
)

cat(
    "Concordant HC associations:",
    nrow(hc),
    "\n\n"
)


cat(
    "Correlation summary:\n"
)

print(
    summary_table,
    row.names = FALSE
)


cat(
    "\nExploratory threshold sensitivity:\n"
)

print(
    sensitivity,
    row.names = FALSE
)


cat(
    "\nConcordant HC tissue summary:\n"
)

print(
    hc_tissue_summary,
    row.names = FALSE
)


cat(
    "\nIMPORTANT INTERPRETATION:\n"
)

cat(
    "1. rho_P348 and rho_P350 are the two independent-animal\n",
    "   cross-tissue coupling estimates.\n",
    sep = ""
)

cat(
    "2. rho_mean uses 8 tissue means and is the principal\n",
    "   cross-tissue summary coefficient.\n",
    sep = ""
)

cat(
    "3. rho_all16 is descriptive only and must NOT be interpreted\n",
    "   as evidence from 16 independent biological replicates.\n",
    sep = ""
)

cat(
    "4. No correlation threshold is finalized in this step.\n",
    "   Threshold tables are sensitivity analyses only.\n",
    sep = ""
)

cat(
    "5. Spearman correlations were calculated from LINEAR\n",
    "   TMM.TPM values. No log transformation was applied.\n",
    sep = ""
)

cat(
    "\nOutput directory:\n",
    outdir,
    "\n"
)

cat(
    "============================================================\n"
)

sessionInfo()

