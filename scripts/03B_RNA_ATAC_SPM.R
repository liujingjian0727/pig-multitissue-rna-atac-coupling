#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# INPUT
# ============================================================

rna_file  <- "sus.TMM.TPM.matrix"
atac_file <- "pig.TMM.TPM.matrix"
meta_file <- "sample_pair_metadata.tsv"

rna_tau_file <- paste0(
    "03_tissue_specificity/RNA/",
    "RNA_tissue_specificity_all.tsv"
)

atac_tau_file <- paste0(
    "03_tissue_specificity/ATAC/",
    "ATAC_tissue_specificity_all.tsv"
)

outdir <- "03B_SPM_analysis"

dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "RNA"), showWarnings = FALSE)
dir.create(file.path(outdir, "ATAC"), showWarnings = FALSE)

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

cat("============================================================\n")
cat("STEP 03B: RNA / ATAC SPM ANALYSIS\n")
cat("SPM calculated from LINEAR TMM.TPM values\n")
cat("NO log transformation\n")
cat("============================================================\n\n")

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

calculate_spm <- function(x) {

    denom <- sqrt(rowSums(x^2))

    result <- matrix(
        NA_real_,
        nrow = nrow(x),
        ncol = ncol(x),
        dimnames = dimnames(x)
    )

    valid <- denom > 0 & is.finite(denom)

    result[valid, ] <- sweep(
        x[valid, , drop = FALSE],
        1,
        denom[valid],
        "/"
    )

    result
}

# ============================================================
# READ
# ============================================================

meta <- read.delim(
    meta_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

rna <- read_matrix(rna_file)
atac <- read_matrix(atac_file)

rna_tau <- read.delim(
    rna_tau_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

atac_tau <- read.delim(
    atac_tau_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

# ============================================================
# MAIN
# ============================================================

run_spm <- function(
    mat,
    metadata,
    sample_column,
    tau_result,
    assay,
    output_dir
) {

    cat("\n============================================================\n")
    cat("Processing:", assay, "\n")
    cat("============================================================\n")

    samples_expected <- metadata[[sample_column]]

    if (!all(samples_expected %in% colnames(mat))) {
        stop(assay, ": sample mismatch.")
    }

    mat <- mat[, samples_expected, drop = FALSE]

    # --------------------------------------------------------
    # 8-tissue mean profile
    # --------------------------------------------------------

    tissue_mean <- sapply(
        tissues,
        function(tissue) {

            ss <- metadata[
                metadata$Tissue == tissue,
                sample_column
            ]

            rowMeans(
                mat[, ss, drop = FALSE]
            )
        }
    )

    colnames(tissue_mean) <- tissues

    # --------------------------------------------------------
    # SPM
    # --------------------------------------------------------

    spm <- calculate_spm(tissue_mean)

    colnames(spm) <- paste0(
        "SPM_",
        tissues
    )

    # --------------------------------------------------------
    # maximum SPM
    # --------------------------------------------------------

    max_spm <- apply(
        spm,
        1,
        function(x) {
            if (all(is.na(x))) return(NA_real_)
            max(x, na.rm = TRUE)
        }
    )

    max_index <- apply(
        spm,
        1,
        function(x) {
            if (all(is.na(x))) return(NA_integer_)
            which.max(x)
        }
    )

    max_tissue <- rep(
        NA_character_,
        length(max_index)
    )

    valid <- !is.na(max_index)

    max_tissue[valid] <- tissues[
        max_index[valid]
    ]

    # --------------------------------------------------------
    # second-highest SPM
    # --------------------------------------------------------

    second_spm <- apply(
        spm,
        1,
        function(x) {

            x <- x[is.finite(x)]

            if (length(x) < 2) {
                return(NA_real_)
            }

            sort(
                x,
                decreasing = TRUE
            )[2]
        }
    )

    spm_margin <- max_spm - second_spm

    spm_ratio <- (
        max_spm + 1e-8
    ) / (
        second_spm + 1e-8
    )

    # --------------------------------------------------------
    # combine
    # --------------------------------------------------------

    result <- data.frame(
        Feature_ID = rownames(spm),
        spm,
        Max_SPM = max_spm,
        Max_SPM_tissue = max_tissue,
        Second_SPM = second_spm,
        SPM_margin = spm_margin,
        SPM_ratio = spm_ratio,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    idx <- match(
        result$Feature_ID,
        tau_result$Feature_ID
    )

    result$Tau <- tau_result$Tau[idx]

    result$Tau_max_tissue <-
        tau_result$Max_tissue[idx]

    result$Max_mean_TMM_TPM <-
        tau_result$Max_mean_TMM_TPM[idx]

    result$LRT_padj <-
        tau_result$LRT_padj[idx]

    result$Replicate_same_max <-
        tau_result$Replicate_same_max[idx]

    # mathematically this should be TRUE
    result$SPM_Tau_same_tissue <- (
        result$Max_SPM_tissue ==
        result$Tau_max_tissue
    )

    # --------------------------------------------------------
    # Combined candidate definitions
    # --------------------------------------------------------

    result$SPM_ge_050 <- (
        !is.na(result$Max_SPM) &
        result$Max_SPM >= 0.50
    )

    result$SPM_ge_060 <- (
        !is.na(result$Max_SPM) &
        result$Max_SPM >= 0.60
    )

    result$SPM_ge_070 <- (
        !is.na(result$Max_SPM) &
        result$Max_SPM >= 0.70
    )

    result$Tau080_SPM050 <- (
        !is.na(result$Tau) &
        result$Tau >= 0.80 &
        result$SPM_ge_050
    )

    result$High_confidence_TauSPM <- (
        !is.na(result$LRT_padj) &
        result$LRT_padj < 0.01 &
        !is.na(result$Tau) &
        result$Tau >= 0.80 &
        result$SPM_ge_050 &
        result$Max_mean_TMM_TPM >= 1 &
        result$Replicate_same_max
    )

    # --------------------------------------------------------
    # Output full SPM matrix
    # --------------------------------------------------------

    write.table(
        result,
        file.path(
            output_dir,
            paste0(
                assay,
                "_SPM_all.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # High-confidence Tau + SPM
    # --------------------------------------------------------

    hc <- result[
        result$High_confidence_TauSPM,
    ]

    write.table(
        hc,
        file.path(
            output_dir,
            paste0(
                assay,
                "_Tau080_SPM050_high_confidence.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # threshold sensitivity
    # --------------------------------------------------------

    threshold_result <- do.call(
        rbind,
        lapply(
            c(0.50, 0.60, 0.70, 0.80, 0.90),
            function(cutoff) {

                candidate <- (
                    !is.na(result$LRT_padj) &
                    result$LRT_padj < 0.01 &
                    !is.na(result$Tau) &
                    result$Tau >= 0.80 &
                    !is.na(result$Max_SPM) &
                    result$Max_SPM >= cutoff &
                    result$Max_mean_TMM_TPM >= 1
                )

                data.frame(
                    SPM_threshold = cutoff,
                    Number = sum(candidate),
                    Number_replicate_consistent =
                        sum(
                            candidate &
                            result$Replicate_same_max
                        )
                )
            }
        )
    )

    write.table(
        threshold_result,
        file.path(
            output_dir,
            paste0(
                assay,
                "_SPM_threshold_sensitivity.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # tissue counts
    # --------------------------------------------------------

    counts <- as.data.frame(
        table(
            factor(
                hc$Max_SPM_tissue,
                levels = tissues
            )
        )
    )

    colnames(counts) <- c(
        "Tissue",
        "High_confidence_Tau080_SPM050"
    )

    write.table(
        counts,
        file.path(
            output_dir,
            paste0(
                assay,
                "_Tau080_SPM050_tissue_counts.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Max SPM distribution
    # --------------------------------------------------------

    plot_data <- result[
        !is.na(result$Max_SPM),
    ]

    p <- ggplot(
        plot_data,
        aes(x = Max_SPM)
    ) +
        geom_histogram(
            bins = 60
        ) +
        geom_vline(
            xintercept = 0.50,
            linetype = 2
        ) +
        geom_vline(
            xintercept = 0.70,
            linetype = 3
        ) +
        theme_bw(base_size = 12) +
        labs(
            title = paste0(
                assay,
                " maximum SPM distribution"
            ),
            x = "Maximum SPM",
            y = "Number of features"
        )

    ggsave(
        file.path(
            output_dir,
            paste0(
                assay,
                "_Max_SPM_distribution.pdf"
            )
        ),
        p,
        width = 7,
        height = 5
    )

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    summary <- data.frame(
        Metric = c(
            "Total_features",
            "Valid_SPM",
            "Max_SPM_ge_0.50",
            "Max_SPM_ge_0.60",
            "Max_SPM_ge_0.70",
            "Tau_ge_0.80_and_SPM_ge_0.50",
            "High_confidence_Tau080_SPM050",
            "SPM_Tau_same_max_tissue_percent"
        ),
        Value = c(
            nrow(result),
            sum(!is.na(result$Max_SPM)),
            sum(result$SPM_ge_050, na.rm = TRUE),
            sum(result$SPM_ge_060, na.rm = TRUE),
            sum(result$SPM_ge_070, na.rm = TRUE),
            sum(result$Tau080_SPM050, na.rm = TRUE),
            nrow(hc),
            round(
                100 *
                mean(
                    result$SPM_Tau_same_tissue,
                    na.rm = TRUE
                ),
                3
            )
        )
    )

    write.table(
        summary,
        file.path(
            output_dir,
            paste0(
                assay,
                "_SPM_summary.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n", assay, " SPM summary:\n", sep="")
    print(summary)

    cat("\n", assay, " SPM sensitivity:\n", sep="")
    print(threshold_result)

    cat("\n", assay, " tissue counts:\n", sep="")
    print(counts)

    invisible(
        list(
            result = result,
            high_confidence = hc,
            summary = summary
        )
    )
}

# ============================================================
# RNA
# ============================================================

rna_res <- run_spm(
    mat = rna,
    metadata = meta,
    sample_column = "RNA_SRR",
    tau_result = rna_tau,
    assay = "RNA",
    output_dir = file.path(
        outdir,
        "RNA"
    )
)

# ============================================================
# ATAC
# ============================================================

atac_res <- run_spm(
    mat = atac,
    metadata = meta,
    sample_column = "ATAC_SRR",
    tau_result = atac_tau,
    assay = "ATAC",
    output_dir = file.path(
        outdir,
        "ATAC"
    )
)

cat("\n============================================================\n")
cat("STEP 03B COMPLETED\n")
cat("============================================================\n")

sessionInfo()
