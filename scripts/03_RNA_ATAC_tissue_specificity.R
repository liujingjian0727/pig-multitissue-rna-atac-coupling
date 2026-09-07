#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# INPUT
# ============================================================

rna_tpm_file  <- "sus.TMM.TPM.matrix"
atac_tpm_file <- "pig.TMM.TPM.matrix"
meta_file     <- "sample_pair_metadata.tsv"

rna_lrt_file <- paste0(
    "02_multi_tissue_effect/RNA/",
    "RNA_multi_tissue_LRT_all.tsv"
)

atac_lrt_file <- paste0(
    "02_multi_tissue_effect/ATAC/",
    "ATAC_multi_tissue_LRT_all.tsv"
)

outdir <- "03_tissue_specificity"

dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "RNA"), showWarnings = FALSE)
dir.create(file.path(outdir, "ATAC"), showWarnings = FALSE)

tissue_levels <- c(
    "Adipose",
    "Cerebellum",
    "Cortex",
    "Hypothalamus",
    "Liver",
    "Lung",
    "Muscle",
    "Spleen"
)

PRIMARY_FDR <- 0.01
PRIMARY_TAU <- 0.80

RNA_MIN_MAX  <- 1
ATAC_MIN_MAX <- 1

cat("============================================================\n")
cat("STEP 03: RNA / ATAC TISSUE SPECIFICITY\n")
cat("============================================================\n")
cat("Primary definition:\n")
cat("LRT FDR <", PRIMARY_FDR, "\n")
cat("Tau >=", PRIMARY_TAU, "\n")
cat("RNA maximum mean TMM.TPM >=", RNA_MIN_MAX, "\n")
cat("ATAC maximum mean TMM.TPM >=", ATAC_MIN_MAX, "\n\n")

# ============================================================
# READ MATRIX
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

    return(x)
}

meta <- read.delim(
    meta_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

meta$Tissue <- factor(
    meta$Tissue,
    levels = tissue_levels
)

rna_tpm  <- read_matrix(rna_tpm_file)
atac_tpm <- read_matrix(atac_tpm_file)

rna_lrt <- read.delim(
    rna_lrt_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

atac_lrt <- read.delim(
    atac_lrt_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

# ============================================================
# TAU
# ============================================================

calculate_tau <- function(x) {

    max_x <- apply(x, 1, max)

    tau <- rep(NA_real_, nrow(x))

    valid <- max_x > 0 & is.finite(max_x)

    scaled <- sweep(
        x[valid, , drop = FALSE],
        1,
        max_x[valid],
        "/"
    )

    tau[valid] <- rowSums(
        1 - scaled
    ) / (ncol(x) - 1)

    tau
}

# ============================================================
# MAIN FUNCTION
# ============================================================

run_specificity <- function(
    matrix,
    metadata,
    sample_column,
    lrt,
    assay,
    min_max,
    assay_dir
) {

    cat("\n============================================================\n")
    cat("Processing:", assay, "\n")
    cat("============================================================\n")

    expected_samples <- metadata[[sample_column]]

    if (!all(expected_samples %in% colnames(matrix))) {
        stop(assay, ": matrix samples do not match metadata.")
    }

    matrix <- matrix[
        ,
        expected_samples,
        drop = FALSE
    ]

    # --------------------------------------------------------
    # Tissue mean profile
    # --------------------------------------------------------

    tissue_mean <- sapply(
        tissue_levels,
        function(tissue) {

            samples <- metadata[
                metadata$Tissue == tissue,
                sample_column
            ]

            rowMeans(
                matrix[, samples, drop = FALSE]
            )
        }
    )

    colnames(tissue_mean) <- tissue_levels

    write.table(
        tissue_mean,
        file.path(
            assay_dir,
            paste0(assay, "_8_tissue_mean_TMM_TPM.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    # --------------------------------------------------------
    # Individual tissue profiles
    # --------------------------------------------------------

    get_animal_profile <- function(animal) {

        profile <- sapply(
            tissue_levels,
            function(tissue) {

                sample <- metadata[
                    metadata$Animal == animal &
                    metadata$Tissue == tissue,
                    sample_column
                ]

                if (length(sample) != 1) {
                    stop(
                        "Expected exactly one sample for ",
                        animal, " ", tissue
                    )
                }

                matrix[, sample]
            }
        )

        colnames(profile) <- tissue_levels

        profile
    }

    profile_P348 <- get_animal_profile("P348")
    profile_P350 <- get_animal_profile("P350")

    # --------------------------------------------------------
    # Tau
    # --------------------------------------------------------

    tau_mean <- calculate_tau(tissue_mean)
    tau_P348 <- calculate_tau(profile_P348)
    tau_P350 <- calculate_tau(profile_P350)

    # --------------------------------------------------------
    # Maximum tissue
    # --------------------------------------------------------

    max_index <- max.col(
        tissue_mean,
        ties.method = "first"
    )

    max_tissue <- tissue_levels[max_index]

    max_value <- tissue_mean[
        cbind(
            seq_len(nrow(tissue_mean)),
            max_index
        )
    ]

    second_value <- apply(
        tissue_mean,
        1,
        function(x) {
            sx <- sort(
                x,
                decreasing = TRUE
            )
            sx[2]
        }
    )

    sum_value <- rowSums(tissue_mean)

    max_fraction <- ifelse(
        sum_value > 0,
        max_value / sum_value,
        NA
    )

    dominance_ratio <- (
        max_value + 1e-6
    ) / (
        second_value + 1e-6
    )

    # --------------------------------------------------------
    # P348 / P350 maximum tissue
    # --------------------------------------------------------

    max_P348 <- tissue_levels[
        max.col(
            profile_P348,
            ties.method = "first"
        )
    ]

    max_P350 <- tissue_levels[
        max.col(
            profile_P350,
            ties.method = "first"
        )
    ]

    replicate_same_max <- (
        max_P348 == max_P350 &
        max_P348 == max_tissue
    )

    # --------------------------------------------------------
    # Add LRT
    # --------------------------------------------------------

    lrt_index <- match(
        rownames(tissue_mean),
        lrt$Feature_ID
    )

    lrt_pvalue <- lrt$pvalue[lrt_index]
    lrt_padj   <- lrt$padj[lrt_index]

    result <- data.frame(
        Feature_ID = rownames(tissue_mean),

        Tau = tau_mean,
        Tau_P348 = tau_P348,
        Tau_P350 = tau_P350,

        Max_tissue = max_tissue,
        Max_mean_TMM_TPM = max_value,
        Second_mean_TMM_TPM = second_value,
        Dominance_ratio = dominance_ratio,
        Max_fraction = max_fraction,

        P348_max_tissue = max_P348,
        P350_max_tissue = max_P350,
        Replicate_same_max = replicate_same_max,

        LRT_pvalue = lrt_pvalue,
        LRT_padj = lrt_padj,

        stringsAsFactors = FALSE
    )

    # --------------------------------------------------------
    # Classification
    # --------------------------------------------------------

    result$LRT_FDR001 <- (
        !is.na(result$LRT_padj) &
        result$LRT_padj < PRIMARY_FDR
    )

    result$Primary_tissue_specific <- (
        result$LRT_FDR001 &
        !is.na(result$Tau) &
        result$Tau >= PRIMARY_TAU &
        result$Max_mean_TMM_TPM >= min_max
    )

    result$High_confidence_tissue_specific <- (
        result$Primary_tissue_specific &
        result$Replicate_same_max
    )

    result <- result[
        order(
            -result$Tau,
            -result$Max_mean_TMM_TPM
        ),
    ]

    # --------------------------------------------------------
    # All metrics
    # --------------------------------------------------------

    write.table(
        result,
        file.path(
            assay_dir,
            paste0(
                assay,
                "_tissue_specificity_all.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Primary set
    # --------------------------------------------------------

    primary <- result[
        result$Primary_tissue_specific,
    ]

    write.table(
        primary,
        file.path(
            assay_dir,
            paste0(
                assay,
                "_primary_tissue_specific.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # High-confidence set
    # --------------------------------------------------------

    high_conf <- result[
        result$High_confidence_tissue_specific,
    ]

    write.table(
        high_conf,
        file.path(
            assay_dir,
            paste0(
                assay,
                "_high_confidence_tissue_specific.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Counts by tissue
    # --------------------------------------------------------

    primary_counts <- as.data.frame(
        table(
            factor(
                primary$Max_tissue,
                levels = tissue_levels
            )
        )
    )

    colnames(primary_counts) <- c(
        "Tissue",
        "Primary_specific_features"
    )

    high_counts <- as.data.frame(
        table(
            factor(
                high_conf$Max_tissue,
                levels = tissue_levels
            )
        )
    )

    colnames(high_counts) <- c(
        "Tissue",
        "High_confidence_features"
    )

    tissue_summary <- merge(
        primary_counts,
        high_counts,
        by = "Tissue",
        all = TRUE
    )

    tissue_summary$Tissue <- factor(
        tissue_summary$Tissue,
        levels = tissue_levels
    )

    tissue_summary <- tissue_summary[
        order(tissue_summary$Tissue),
    ]

    write.table(
        tissue_summary,
        file.path(
            assay_dir,
            paste0(
                assay,
                "_tissue_specific_counts.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Sensitivity analysis
    # --------------------------------------------------------

    sensitivity <- list()
    k <- 1

    for (tau_cut in c(0.70, 0.80, 0.90)) {

        for (signal_cut in c(0.5, 1, 2)) {

            x <- (
                result$LRT_FDR001 &
                !is.na(result$Tau) &
                result$Tau >= tau_cut &
                result$Max_mean_TMM_TPM >= signal_cut
            )

            sensitivity[[k]] <- data.frame(
                Tau_threshold = tau_cut,
                Max_signal_threshold = signal_cut,
                Number = sum(x),
                Number_replicate_consistent =
                    sum(x & result$Replicate_same_max)
            )

            k <- k + 1
        }
    }

    sensitivity <- do.call(
        rbind,
        sensitivity
    )

    write.table(
        sensitivity,
        file.path(
            assay_dir,
            paste0(
                assay,
                "_threshold_sensitivity.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Tau distribution
    # --------------------------------------------------------

    plot_data <- result[
        result$LRT_FDR001 &
        result$Max_mean_TMM_TPM >= min_max &
        !is.na(result$Tau),
    ]

    p_tau <- ggplot(
        plot_data,
        aes(x = Tau)
    ) +
        geom_histogram(
            bins = 50
        ) +
        geom_vline(
            xintercept = PRIMARY_TAU,
            linetype = 2
        ) +
        theme_bw(base_size = 12) +
        labs(
            title = paste0(
                assay,
                " tissue-specificity (Tau)"
            ),
            x = "Tau",
            y = "Number of features"
        )

    ggsave(
        file.path(
            assay_dir,
            paste0(
                assay,
                "_Tau_distribution.pdf"
            )
        ),
        p_tau,
        width = 7,
        height = 5
    )

    # --------------------------------------------------------
    # Tissue-specific count barplot
    # --------------------------------------------------------

    plot_counts <- tissue_summary

    p_count <- ggplot(
        plot_counts,
        aes(
            x = Tissue,
            y = High_confidence_features
        )
    ) +
        geom_col() +
        coord_flip() +
        theme_bw(base_size = 12) +
        labs(
            title = paste0(
                assay,
                " high-confidence tissue-specific features"
            ),
            x = NULL,
            y = "Number of features"
        )

    ggsave(
        file.path(
            assay_dir,
            paste0(
                assay,
                "_high_confidence_counts.pdf"
            )
        ),
        p_count,
        width = 7,
        height = 5
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary_table <- data.frame(
        Metric = c(
            "Total_features",
            "Features_in_LRT",
            "LRT_FDR_lt_0.01",
            "Tau_ge_0.80",
            "Primary_tissue_specific",
            "High_confidence_tissue_specific",
            "High_confidence_fraction_of_primary"
        ),
        Value = c(
            nrow(result),
            sum(!is.na(result$LRT_padj)),
            sum(result$LRT_FDR001),
            sum(
                !is.na(result$Tau) &
                result$Tau >= 0.80
            ),
            nrow(primary),
            nrow(high_conf),
            ifelse(
                nrow(primary) > 0,
                round(
                    100 *
                    nrow(high_conf) /
                    nrow(primary),
                    3
                ),
                NA
            )
        )
    )

    write.table(
        summary_table,
        file.path(
            assay_dir,
            paste0(
                assay,
                "_tissue_specificity_summary.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n", assay, " summary:\n", sep = "")
    print(summary_table)

    cat("\n", assay, " tissue counts:\n", sep = "")
    print(tissue_summary)

    return(
        list(
            result = result,
            primary = primary,
            high = high_conf,
            summary = summary_table,
            tissue_summary = tissue_summary
        )
    )
}

# ============================================================
# RNA
# ============================================================

rna_result <- run_specificity(
    matrix = rna_tpm,
    metadata = meta,
    sample_column = "RNA_SRR",
    lrt = rna_lrt,
    assay = "RNA",
    min_max = RNA_MIN_MAX,
    assay_dir = file.path(
        outdir,
        "RNA"
    )
)

# ============================================================
# ATAC
# ============================================================

atac_result <- run_specificity(
    matrix = atac_tpm,
    metadata = meta,
    sample_column = "ATAC_SRR",
    lrt = atac_lrt,
    assay = "ATAC",
    min_max = ATAC_MIN_MAX,
    assay_dir = file.path(
        outdir,
        "ATAC"
    )
)

# ============================================================
# COMBINED SUMMARY
# ============================================================

combined <- data.frame(
    Assay = c("RNA", "ATAC"),
    Primary_specific = c(
        nrow(rna_result$primary),
        nrow(atac_result$primary)
    ),
    High_confidence_specific = c(
        nrow(rna_result$high),
        nrow(atac_result$high)
    )
)

write.table(
    combined,
    file.path(
        outdir,
        "RNA_ATAC_tissue_specificity_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n============================================================\n")
cat("STEP 03 COMPLETED\n")
cat("============================================================\n")

print(combined)

cat("\nPrimary:\n")
cat(
    "LRT FDR < 0.01 + Tau >= 0.80 + max mean TMM.TPM >= 1\n"
)

cat("\nHigh confidence:\n")
cat(
    "Primary + same maximum tissue in P348 and P350\n"
)

cat("============================================================\n")

sessionInfo()
