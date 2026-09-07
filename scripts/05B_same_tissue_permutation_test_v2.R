#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(matrixStats)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# STEP 05B
# SAME-TISSUE PERMUTATION TEST
#
# Null hypothesis:
# RNA-ATAC coupling of observed UROPA-linked peak-gene pairs
# is no stronger than random peak-gene pairing within the same
# tissue-specific class.
#
# Permutation preserves:
#   1. tissue identity
#   2. number of peaks per tissue
#   3. multiset / frequency of RNA genes per tissue
#   4. ATAC and RNA tissue-specific profiles
#
# It destroys only:
#   observed peak -> UROPA gene pairing
# ============================================================


# ============================================================
# SETTINGS
# ============================================================

n_perm <- as.integer(
    Sys.getenv(
        "NPERM",
        "1000"
    )
)

seed <- 20260830

set.seed(seed)

# Numerical tolerance for threshold comparisons.
# Prevents theoretically identical Spearman coefficients
# such as 0.5 or 0.7 from being classified differently
# because of ~1e-16 floating-point error.
eps <- 1e-12

cat("============================================================\n")
cat("STEP 05B: SAME-TISSUE PERMUTATION TEST\n")
cat("============================================================\n")
cat("Number of permutations :", n_perm, "\n")
cat("Random seed            :", seed, "\n\n")


# ============================================================
# INPUT
# ============================================================

hc_file <- paste0(
    "05_RNA_ATAC_quantitative_coupling/",
    "05_concordant_HC_peak_gene_correlations.tsv"
)

atac_file <- "pig.TMM.TPM.matrix"
rna_file  <- "sus.TMM.TPM.matrix"

meta_file <- "sample_pair_metadata.tsv"

outdir <- "05B_same_tissue_permutation_v2"

dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
)


# ============================================================
# FUNCTIONS
# ============================================================

read_matrix <- function(file) {

    x <- read.delim(
        file,
        header = TRUE,
        row.names = 1,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    x <- as.matrix(x)
    storage.mode(x) <- "numeric"

    x
}


# ------------------------------------------------------------
# Row ranks -> centered and unit-length vectors
#
# Pearson correlation between these standardized rank vectors
# equals Spearman correlation.
# ------------------------------------------------------------

rank_standardize_rows <- function(x) {

    r <- matrixStats::rowRanks(
        x,
        ties.method = "average",
        useNames = FALSE
    )

    r <- sweep(
        r,
        1,
        rowMeans(r),
        "-"
    )

    norm <- sqrt(
        rowSums(r^2)
    )

    valid <- (
        is.finite(norm) &
        norm > 0
    )

    r[valid, ] <- sweep(
        r[valid, , drop = FALSE],
        1,
        norm[valid],
        "/"
    )

    r[!valid, ] <- NA_real_

    r
}


row_cor_standardized <- function(x, y) {

    rho <- rowSums(
        x * y
    )

    rho[
        !is.finite(rho)
    ] <- NA_real_

    rho
}


safe_median <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    median(x)
}


empirical_p_upper <- function(
    observed,
    null
) {

    null <- null[
        is.finite(null)
    ]

    if (
        !is.finite(observed) ||
        length(null) == 0
    ) {
        return(NA_real_)
    }

    (
        1 +
        sum(
            null >= observed
        )
    ) /
    (
        length(null) + 1
    )
}


summarize_stat <- function(
    statistic,
    observed,
    null
) {

    null <- null[
        is.finite(null)
    ]

    null_mean <- mean(null)
    null_sd   <- sd(null)

    z <- if (
        is.finite(null_sd) &&
        null_sd > 0
    ) {

        (
            observed -
            null_mean
        ) /
        null_sd

    } else {

        NA_real_
    }

    fold <- if (
        is.finite(null_mean) &&
        null_mean != 0
    ) {

        observed /
        null_mean

    } else {

        NA_real_
    }

    data.frame(
        Statistic = statistic,
        Observed = observed,
        Null_mean = null_mean,
        Null_median = median(null),
        Null_SD = null_sd,
        Null_Q025 = as.numeric(
            quantile(
                null,
                0.025,
                names = FALSE
            )
        ),
        Null_Q975 = as.numeric(
            quantile(
                null,
                0.975,
                names = FALSE
            )
        ),
        Observed_over_null_mean = fold,
        Z_score = z,
        Empirical_P_upper =
            empirical_p_upper(
                observed,
                null
            ),
        stringsAsFactors = FALSE
    )
}


# ============================================================
# READ DATA
# ============================================================

cat("Reading concordant HC associations...\n")

hc <- read.delim(
    hc_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("Reading ATAC matrix...\n")
atac <- read_matrix(
    atac_file
)

cat("Reading RNA matrix...\n")
rna <- read_matrix(
    rna_file
)

cat("Reading metadata...\n")

meta <- read.delim(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


cat("\nInput:\n")
cat("HC peak-gene pairs :", nrow(hc), "\n")
cat("Unique genes       :", length(unique(hc$gene_id)), "\n")
cat("ATAC matrix        :", nrow(atac), "x", ncol(atac), "\n")
cat("RNA matrix         :", nrow(rna), "x", ncol(rna), "\n\n")


if (nrow(hc) != 4674) {

    warning(
        "Expected 4674 concordant HC pairs, found ",
        nrow(hc)
    )
}


# ============================================================
# FIX TISSUE ORDER
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

meta_key <- paste(
    meta$Tissue,
    meta$Animal,
    sep = "__"
)

idx348 <- match(
    paste(
        tissues,
        "P348",
        sep = "__"
    ),
    meta_key
)

idx350 <- match(
    paste(
        tissues,
        "P350",
        sep = "__"
    ),
    meta_key
)

if (
    anyNA(idx348) ||
    anyNA(idx350)
) {

    stop(
        "Incomplete tissue x animal metadata."
    )
}


# ============================================================
# MAP ASSOCIATIONS TO MATRICES
# ============================================================

atac_idx <- match(
    hc$SAF_peak_id,
    rownames(atac)
)

rna_idx <- match(
    hc$gene_id,
    rownames(rna)
)

if (anyNA(atac_idx)) {
    stop("HC peak missing from ATAC matrix.")
}

if (anyNA(rna_idx)) {
    stop("HC gene missing from RNA matrix.")
}


ATAC348 <- atac[
    atac_idx,
    meta$ATAC_SRR[idx348],
    drop = FALSE
]

ATAC350 <- atac[
    atac_idx,
    meta$ATAC_SRR[idx350],
    drop = FALSE
]

RNA348 <- rna[
    rna_idx,
    meta$RNA_SRR[idx348],
    drop = FALSE
]

RNA350 <- rna[
    rna_idx,
    meta$RNA_SRR[idx350],
    drop = FALSE
]

ATACmean <- (
    ATAC348 +
    ATAC350
) / 2

RNAmean <- (
    RNA348 +
    RNA350
) / 2


# ============================================================
# PRECOMPUTE STANDARDIZED RANK PROFILES
# ============================================================

cat(
    "Precomputing rank-standardized profiles...\n"
)

A348 <- rank_standardize_rows(
    ATAC348
)

A350 <- rank_standardize_rows(
    ATAC350
)

Amean <- rank_standardize_rows(
    ATACmean
)

R348 <- rank_standardize_rows(
    RNA348
)

R350 <- rank_standardize_rows(
    RNA350
)

Rmean <- rank_standardize_rows(
    RNAmean
)


# ============================================================
# OBSERVED CORRELATIONS
# ============================================================

obs348_recomputed <- row_cor_standardized(
    A348,
    R348
)

obs350_recomputed <- row_cor_standardized(
    A350,
    R350
)

obsmean_recomputed <- row_cor_standardized(
    Amean,
    Rmean
)

# Canonical observed correlations are inherited from Step05.
# Recomputed values are used only for QC.
obs348 <- hc$rho_P348
obs350 <- hc$rho_P350
obsmean <- hc$rho_mean


# Verify against Step05
# ============================================================

check348 <- max(
    abs(
        obs348_recomputed -
        hc$rho_P348
    ),
    na.rm = TRUE
)

check350 <- max(
    abs(
        obs350_recomputed -
        hc$rho_P350
    ),
    na.rm = TRUE
)

checkmean <- max(
    abs(
        obsmean_recomputed -
        hc$rho_mean
    ),
    na.rm = TRUE
)

cat("\nCorrelation reconstruction QC:\n")
cat("Max difference rho_P348 :", check348, "\n")
cat("Max difference rho_P350 :", check350, "\n")
cat("Max difference rho_mean :", checkmean, "\n")

if (
    check348 > 1e-10 ||
    check350 > 1e-10 ||
    checkmean > 1e-10
) {

    warning(
        "Recomputed correlations differ slightly ",
        "from Step05."
    )

} else {

    cat(
        "PASS: Step05 correlations reproduced.\n"
    )
}


# ============================================================
# OBSERVED STATISTICS
# ============================================================

valid_obs <- (
    is.finite(obs348) &
    is.finite(obs350) &
    is.finite(obsmean)
)

positive_obs <- (
    valid_obs &
    obs348 > eps &
    obs350 > eps &
    obsmean > eps
)

moderate_obs <- (
    valid_obs &
    obs348 >= 0.3 - eps &
    obs350 >= 0.3 - eps &
    obsmean >= 0.5 - eps
)

strong_obs <- (
    valid_obs &
    obs348 >= 0.5 - eps &
    obs350 >= 0.5 - eps &
    obsmean >= 0.7 - eps
)

very_strong_obs <- (
    valid_obs &
    obs348 >= 0.7 - eps &
    obs350 >= 0.7 - eps &
    obsmean >= 0.8 - eps
)


observed_statistics <- c(
    Median_rho_P348 =
        safe_median(obs348),

    Median_rho_P350 =
        safe_median(obs350),

    Median_rho_mean =
        safe_median(obsmean),

    Positive_all_three_count =
        sum(positive_obs),

    Moderate_03_03_05_count =
        sum(moderate_obs),

    Strong_05_05_07_count =
        sum(strong_obs),

    Very_strong_07_07_08_count =
        sum(very_strong_obs)
)


cat("\nObserved statistics:\n")
print(observed_statistics)


# ============================================================
# TISSUE GROUPS
# ============================================================

pair_tissue <- hc$ATAC_Max_tissue

if (
    any(
        pair_tissue !=
        hc$RNA_Max_tissue
    )
) {

    stop(
        "Concordant HC table contains mismatched tissues."
    )
}


group_indices <- lapply(
    tissues,
    function(t) {
        which(
            pair_tissue == t
        )
    }
)

names(group_indices) <- tissues


cat("\nPair number per tissue:\n")

for (t in tissues) {

    cat(
        sprintf(
            "%-15s %d\n",
            t,
            length(
                group_indices[[t]]
            )
        )
    )
}


# ============================================================
# NULL STORAGE
# ============================================================

null_global <- data.frame(
    Permutation = seq_len(n_perm),

    Median_rho_P348 =
        NA_real_,

    Median_rho_P350 =
        NA_real_,

    Median_rho_mean =
        NA_real_,

    Positive_all_three_count =
        NA_integer_,

    Moderate_03_03_05_count =
        NA_integer_,

    Strong_05_05_07_count =
        NA_integer_,

    Very_strong_07_07_08_count =
        NA_integer_
)


null_tissue_list <- vector(
    "list",
    n_perm * length(tissues)
)

list_counter <- 0L


# ============================================================
# OPTIONAL PER-PAIR EMPIRICAL NULL FOR rho_mean
#
# For every observed peak, count how often a random same-tissue
# gene produces rho_mean >= observed rho_mean.
# ============================================================

pair_exceed_mean <- integer(
    nrow(hc)
)


# ============================================================
# PERMUTATION LOOP
# ============================================================

cat("\n============================================================\n")
cat("Running permutations...\n")
cat("============================================================\n")

n_pairs <- nrow(hc)


for (b in seq_len(n_perm)) {

    # --------------------------------------------------------
    # Randomly permute RNA association rows WITHIN each tissue
    # --------------------------------------------------------

    perm_idx <- seq_len(
        n_pairs
    )

    for (t in tissues) {

        idx <- group_indices[[t]]

        if (length(idx) > 1) {

            perm_idx[idx] <- sample(
                idx,
                length(idx),
                replace = FALSE
            )
        }
    }


    # --------------------------------------------------------
    # Permuted correlations
    # --------------------------------------------------------

    r348 <- row_cor_standardized(
        A348,
        R348[
            perm_idx,
            ,
            drop = FALSE
        ]
    )

    r350 <- row_cor_standardized(
        A350,
        R350[
            perm_idx,
            ,
            drop = FALSE
        ]
    )

    rmean <- row_cor_standardized(
        Amean,
        Rmean[
            perm_idx,
            ,
            drop = FALSE
        ]
    )


    valid <- (
        is.finite(r348) &
        is.finite(r350) &
        is.finite(rmean)
    )


    positive <- (
        valid &
        r348 > eps &
        r350 > eps &
        rmean > eps
    )


    moderate <- (
        valid &
        r348 >= 0.3 - eps &
        r350 >= 0.3 - eps &
        rmean >= 0.5 - eps
    )


    strong <- (
        valid &
        r348 >= 0.5 - eps &
        r350 >= 0.5 - eps &
        rmean >= 0.7 - eps
    )


    very_strong <- (
        valid &
        r348 >= 0.7 - eps &
        r350 >= 0.7 - eps &
        rmean >= 0.8 - eps
    )


    # --------------------------------------------------------
    # Global null
    # --------------------------------------------------------

    null_global$Median_rho_P348[b] <-
        safe_median(r348)

    null_global$Median_rho_P350[b] <-
        safe_median(r350)

    null_global$Median_rho_mean[b] <-
        safe_median(rmean)

    null_global$Positive_all_three_count[b] <-
        sum(positive)

    null_global$Moderate_03_03_05_count[b] <-
        sum(moderate)

    null_global$Strong_05_05_07_count[b] <-
        sum(strong)

    null_global$Very_strong_07_07_08_count[b] <-
        sum(very_strong)


    # --------------------------------------------------------
    # Per-pair empirical rho_mean null
    # --------------------------------------------------------

    pair_exceed_mean <- (
        pair_exceed_mean +
        as.integer(
            is.finite(rmean) &
            rmean >= obsmean
        )
    )


    # --------------------------------------------------------
    # Tissue-specific null
    # --------------------------------------------------------

    for (t in tissues) {

        idx <- group_indices[[t]]

        list_counter <- list_counter + 1L

        if (length(idx) == 0) {
            next
        }

        null_tissue_list[[list_counter]] <-
            data.frame(
                Permutation = b,
                Tissue = t,
                N_pairs = length(idx),

                Median_rho_mean =
                    safe_median(
                        rmean[idx]
                    ),

                Positive_all_three_count =
                    sum(
                        positive[idx]
                    ),

                Strong_05_05_07_count =
                    sum(
                        strong[idx]
                    ),

                stringsAsFactors = FALSE
            )
    }


    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------

    if (
        b %% 100 == 0 ||
        b == n_perm
    ) {

        cat(
            "Completed:",
            b,
            "/",
            n_perm,
            "\n"
        )
    }
}


# ============================================================
# COMBINE TISSUE NULL
# ============================================================

null_tissue <- do.call(
    rbind,
    null_tissue_list[
        !vapply(
            null_tissue_list,
            is.null,
            logical(1)
        )
    ]
)


# ============================================================
# WRITE NULL DISTRIBUTIONS
# ============================================================

write.table(
    null_global,
    file.path(
        outdir,
        "05B_global_permutation_null.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


write.table(
    null_tissue,
    file.path(
        outdir,
        "05B_tissue_permutation_null.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# GLOBAL PERMUTATION SUMMARY
# ============================================================

global_summary <- rbind(

    summarize_stat(
        "Median_rho_P348",
        observed_statistics[
            "Median_rho_P348"
        ],
        null_global$Median_rho_P348
    ),

    summarize_stat(
        "Median_rho_P350",
        observed_statistics[
            "Median_rho_P350"
        ],
        null_global$Median_rho_P350
    ),

    summarize_stat(
        "Median_rho_mean",
        observed_statistics[
            "Median_rho_mean"
        ],
        null_global$Median_rho_mean
    ),

    summarize_stat(
        "Positive_all_three_count",
        observed_statistics[
            "Positive_all_three_count"
        ],
        null_global$Positive_all_three_count
    ),

    summarize_stat(
        "Moderate_03_03_05_count",
        observed_statistics[
            "Moderate_03_03_05_count"
        ],
        null_global$Moderate_03_03_05_count
    ),

    summarize_stat(
        "Strong_05_05_07_count",
        observed_statistics[
            "Strong_05_05_07_count"
        ],
        null_global$Strong_05_05_07_count
    ),

    summarize_stat(
        "Very_strong_07_07_08_count",
        observed_statistics[
            "Very_strong_07_07_08_count"
        ],
        null_global$Very_strong_07_07_08_count
    )
)


write.table(
    global_summary,
    file.path(
        outdir,
        "05B_global_permutation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PER-TISSUE OBSERVED + NULL SUMMARY
# ============================================================

tissue_summary_list <- lapply(
    tissues,
    function(t) {

        idx <- group_indices[[t]]

        null_t <- null_tissue[
            null_tissue$Tissue == t,
            ,
            drop = FALSE
        ]


        obs_med <- safe_median(
            obsmean[idx]
        )


        obs_positive <- sum(
            positive_obs[idx]
        )


        obs_strong <- sum(
            strong_obs[idx]
        )


        data.frame(
            Tissue = t,

            N_pairs =
                length(idx),

            Observed_median_rho_mean =
                obs_med,

            Null_mean_median_rho_mean =
                mean(
                    null_t$Median_rho_mean
                ),

            Null_Q025_median_rho_mean =
                as.numeric(
                    quantile(
                        null_t$Median_rho_mean,
                        0.025,
                        names = FALSE
                    )
                ),

            Null_Q975_median_rho_mean =
                as.numeric(
                    quantile(
                        null_t$Median_rho_mean,
                        0.975,
                        names = FALSE
                    )
                ),

            P_median_rho_mean =
                empirical_p_upper(
                    obs_med,
                    null_t$Median_rho_mean
                ),

            Observed_positive_all_three =
                obs_positive,

            Null_mean_positive_all_three =
                mean(
                    null_t$Positive_all_three_count
                ),

            P_positive_all_three =
                empirical_p_upper(
                    obs_positive,
                    null_t$Positive_all_three_count
                ),

            Observed_strong_05_05_07 =
                obs_strong,

            Null_mean_strong_05_05_07 =
                mean(
                    null_t$Strong_05_05_07_count
                ),

            Strong_fold_enrichment =
                ifelse(
                    mean(
                        null_t$Strong_05_05_07_count
                    ) > 0,

                    obs_strong /
                    mean(
                        null_t$Strong_05_05_07_count
                    ),

                    NA
                ),

            P_strong_05_05_07 =
                empirical_p_upper(
                    obs_strong,
                    null_t$Strong_05_05_07_count
                ),

            stringsAsFactors = FALSE
        )
    }
)


tissue_summary <- do.call(
    rbind,
    tissue_summary_list
)


# BH adjustment for the 8 tissue-specific tests
tissue_summary$FDR_median_rho_mean <-
    p.adjust(
        tissue_summary$P_median_rho_mean,
        method = "BH"
    )

tissue_summary$FDR_positive_all_three <-
    p.adjust(
        tissue_summary$P_positive_all_three,
        method = "BH"
    )

tissue_summary$FDR_strong_05_05_07 <-
    p.adjust(
        tissue_summary$P_strong_05_05_07,
        method = "BH"
    )


write.table(
    tissue_summary,
    file.path(
        outdir,
        "05B_tissue_permutation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PER-PAIR EMPIRICAL P FOR rho_mean
# ============================================================

hc$Permutation_p_rho_mean <- (
    pair_exceed_mean + 1
) / (
    n_perm + 1
)

hc$Permutation_FDR_rho_mean <-
    p.adjust(
        hc$Permutation_p_rho_mean,
        method = "BH"
    )


hc$Observed_strong_05_05_07 <-
    strong_obs


write.table(
    hc,
    file.path(
        outdir,
        "05B_concordant_HC_pairs_with_permutation_p.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# OUTPUT THE 702-STYLE STRONG SET
# ============================================================

strong_pairs <- hc[
    strong_obs,
    ,
    drop = FALSE
]


write.table(
    strong_pairs,
    file.path(
        outdir,
        "05B_observed_strong_05_05_07_pairs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# PLOTS
# ============================================================

plot_long <- rbind(

    data.frame(
        Statistic =
            "Median rho(mean)",
        Null_value =
            null_global$Median_rho_mean,
        Observed =
            observed_statistics[
                "Median_rho_mean"
            ]
    ),

    data.frame(
        Statistic =
            "Positive all three",
        Null_value =
            null_global$Positive_all_three_count,
        Observed =
            observed_statistics[
                "Positive_all_three_count"
            ]
    ),

    data.frame(
        Statistic =
            "Strong 0.5/0.5/0.7",
        Null_value =
            null_global$Strong_05_05_07_count,
        Observed =
            observed_statistics[
                "Strong_05_05_07_count"
            ]
    )
)


p1 <- ggplot(
    plot_long,
    aes(
        x = Null_value
    )
) +
    geom_histogram(
        bins = 40
    ) +
    geom_vline(
        aes(
            xintercept = Observed
        ),
        linetype = 2,
        linewidth = 0.8
    ) +
    facet_wrap(
        ~ Statistic,
        scales = "free",
        ncol = 1
    ) +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title = paste0(
            "Same-tissue permutation null (",
            n_perm,
            " permutations)"
        ),
        x = "Null statistic",
        y = "Number of permutations"
    )


ggsave(
    file.path(
        outdir,
        "05B_global_permutation_null.pdf"
    ),
    p1,
    width = 7,
    height = 8
)


# ------------------------------------------------------------
# Tissue strong-pair observed vs null expectation
# ------------------------------------------------------------

plot_tissue <- tissue_summary

plot_tissue$Tissue <- factor(
    plot_tissue$Tissue,
    levels = tissues
)


p2 <- ggplot(
    plot_tissue,
    aes(
        x = Tissue
    )
) +
    geom_col(
        aes(
            y = Observed_strong_05_05_07
        ),
        alpha = 0.5
    ) +
    geom_point(
        aes(
            y = Null_mean_strong_05_05_07
        ),
        size = 2.5
    ) +
    theme_bw(
        base_size = 11
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    ) +
    labs(
        title =
            "Strong RNA-ATAC coupling: observed vs same-tissue null",
        x = NULL,
        y = "Number of strong peak-gene pairs"
    )


ggsave(
    file.path(
        outdir,
        "05B_tissue_strong_pair_enrichment.pdf"
    ),
    p2,
    width = 8,
    height = 5.5
)


# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("STEP 05B COMPLETED\n")
cat("============================================================\n\n")

cat("GLOBAL PERMUTATION RESULTS:\n\n")

print(
    global_summary,
    row.names = FALSE
)


cat("\nTISSUE-SPECIFIC RESULTS:\n\n")

print(
    tissue_summary,
    row.names = FALSE
)


cat("\nObserved strong 0.5/0.5/0.7 pairs:",
    nrow(strong_pairs),
    "\n"
)


cat("\nInterpretation:\n")
cat(
    "Empirical P = (1 + number of null statistics >= observed) / ",
    "(Nperm + 1)\n",
    sep = ""
)

cat(
    "The null preserves tissue identity and gene multiplicity, ",
    "but breaks the observed UROPA peak-gene association.\n"
)

cat(
    "Per-pair permutation P/FDR for rho_mean is exploratory; ",
    "the principal inference is the aggregate enrichment test.\n"
)

cat(
    "For very small tissue groups (especially Adipose), ",
    "tissue-specific empirical P-values should be interpreted ",
    "cautiously.\n"
)

cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("============================================================\n")

sessionInfo()

