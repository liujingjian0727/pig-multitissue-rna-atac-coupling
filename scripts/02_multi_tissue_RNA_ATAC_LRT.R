#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

options(stringsAsFactors = FALSE)

# ============================================================
# INPUT
# ============================================================

rna_file  <- "sus.counts.matrix"
atac_file <- "pig.counts.matrix"
meta_file <- "sample_pair_metadata.tsv"

outdir <- "02_multi_tissue_effect"
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "RNA"), showWarnings = FALSE)
dir.create(file.path(outdir, "ATAC"), showWarnings = FALSE)

cat("============================================================\n")
cat("STEP 02: MULTI-TISSUE RNA / ATAC GLOBAL EFFECT ANALYSIS\n")
cat("Model: ~ Animal + Tissue\n")
cat("LRT reduced model: ~ Animal\n")
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

    return(x)
}

# ============================================================
# READ METADATA
# ============================================================

meta <- read.delim(
    meta_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

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

meta$Tissue <- factor(meta$Tissue, levels = tissue_levels)
meta$Animal <- factor(meta$Animal)

# RNA metadata
rna_meta <- data.frame(
    Animal = meta$Animal,
    Tissue = meta$Tissue,
    row.names = meta$RNA_SRR
)

# ATAC metadata
atac_meta <- data.frame(
    Animal = meta$Animal,
    Tissue = meta$Tissue,
    row.names = meta$ATAC_SRR
)

# ============================================================
# READ MATRICES
# ============================================================

cat("Reading RNA counts...\n")
rna <- read_matrix(rna_file)

cat("Reading ATAC counts...\n")
atac <- read_matrix(atac_file)

cat("\nRNA :", nrow(rna), "features x", ncol(rna), "samples\n")
cat("ATAC:", nrow(atac), "features x", ncol(atac), "samples\n")

# ============================================================
# CHECK SAMPLE ORDER
# ============================================================

if (!identical(colnames(rna), rownames(rna_meta))) {
    stop("RNA sample order does not match metadata.")
}

if (!identical(colnames(atac), rownames(atac_meta))) {
    stop("ATAC sample order does not match metadata.")
}

cat("\nPASS: sample matching confirmed.\n")

# ============================================================
# MAIN ANALYSIS FUNCTION
# ============================================================

run_lrt <- function(
    count_matrix,
    metadata,
    assay_name,
    output_dir,
    heatmap_n
) {

    cat("\n============================================================\n")
    cat("Processing:", assay_name, "\n")
    cat("============================================================\n")

    # --------------------------------------------------------
    # PREFILTER
    # Keep features with count >=10 in >=2 samples
    # --------------------------------------------------------

    keep <- rowSums(count_matrix >= 10) >= 2

    counts_filtered <- count_matrix[keep, , drop = FALSE]

    cat("Features before filtering:", nrow(count_matrix), "\n")
    cat("Features after filtering :", nrow(counts_filtered), "\n")

    # --------------------------------------------------------
    # DESeq2 dataset
    # --------------------------------------------------------

    dds <- DESeqDataSetFromMatrix(
        countData = round(counts_filtered),
        colData = metadata,
        design = ~ Animal + Tissue
    )

    # --------------------------------------------------------
    # LRT
    # Full:    Animal + Tissue
    # Reduced: Animal
    # --------------------------------------------------------

    cat("Running DESeq2 LRT...\n")

    dds <- DESeq(
        dds,
        test = "LRT",
        reduced = ~ Animal,
        quiet = TRUE
    )

    res <- results(
        dds,
        alpha = 0.05,
        independentFiltering = TRUE
    )

    res_df <- as.data.frame(res)

    res_df$Feature_ID <- rownames(res_df)

    # LRT log2FC column is not a global tissue effect size.
    # Therefore keep it in raw result but do not use it for
    # biological interpretation of the omnibus tissue test.

    res_df <- res_df[
        ,
        c(
            "Feature_ID",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        )
    ]

    res_df <- res_df[
        order(
            is.na(res_df$padj),
            res_df$padj,
            res_df$pvalue
        ),
    ]

    write.table(
        res_df,
        file.path(
            output_dir,
            paste0(assay_name, "_multi_tissue_LRT_all.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Significant tissue-dependent features
    # --------------------------------------------------------

    sig <- res_df[
        !is.na(res_df$padj) &
        res_df$padj < 0.05,
    ]

    write.table(
        sig,
        file.path(
            output_dir,
            paste0(assay_name, "_multi_tissue_LRT_FDR005.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # More stringent set
    sig001 <- res_df[
        !is.na(res_df$padj) &
        res_df$padj < 0.01,
    ]

    write.table(
        sig001,
        file.path(
            output_dir,
            paste0(assay_name, "_multi_tissue_LRT_FDR001.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # --------------------------------------------------------
    # Normalized counts
    # --------------------------------------------------------

    norm_counts <- counts(dds, normalized = TRUE)

    write.table(
        norm_counts,
        file.path(
            output_dir,
            paste0(assay_name, "_DESeq2_normalized_counts.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    # --------------------------------------------------------
    # VST for visualization
    # --------------------------------------------------------

    cat("Running VST for tissue-pattern visualization...\n")

    vsd <- vst(dds, blind = FALSE)

    vst_mat <- assay(vsd)

    # --------------------------------------------------------
    # Average P348/P350 within tissue
    # --------------------------------------------------------

    tissue_mean <- sapply(
        tissue_levels,
        function(tissue) {

            samples <- rownames(metadata)[
                metadata$Tissue == tissue
            ]

            rowMeans(
                vst_mat[, samples, drop = FALSE]
            )
        }
    )

    colnames(tissue_mean) <- tissue_levels

    write.table(
        tissue_mean,
        file.path(
            output_dir,
            paste0(assay_name, "_tissue_mean_VST.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    # --------------------------------------------------------
    # Top significant features heatmap
    # --------------------------------------------------------

    sig_ids <- sig$Feature_ID

    sig_ids <- sig_ids[
        sig_ids %in% rownames(tissue_mean)
    ]

    n_use <- min(
        heatmap_n,
        length(sig_ids)
    )

    if (n_use >= 2) {

        top_ids <- sig_ids[seq_len(n_use)]

        hm <- tissue_mean[
            top_ids,
            ,
            drop = FALSE
        ]

        # Row z-score
        hm_z <- t(scale(t(hm)))

        # Remove zero-variance / NA rows
        hm_z <- hm_z[
            complete.cases(hm_z),
            ,
            drop = FALSE
        ]

        pdf(
            file.path(
                output_dir,
                paste0(
                    assay_name,
                    "_top",
                    nrow(hm_z),
                    "_tissue_effect_heatmap.pdf"
                )
            ),
            width = 8,
            height = 10
        )

        pheatmap(
            hm_z,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            border_color = NA,
            scale = "none",
            main = paste0(
                assay_name,
                ": top tissue-dependent features"
            )
        )

        dev.off()
    }

    # --------------------------------------------------------
    # P-value / FDR summary
    # --------------------------------------------------------

    summary_table <- data.frame(
        Metric = c(
            "Input_features",
            "Filtered_features",
            "Tested_features",
            "FDR_lt_0.05",
            "FDR_lt_0.01",
            "FDR_lt_0.001",
            "Percent_FDR_lt_0.05"
        ),
        Value = c(
            nrow(count_matrix),
            nrow(counts_filtered),
            nrow(res_df),
            sum(!is.na(res_df$padj) & res_df$padj < 0.05),
            sum(!is.na(res_df$padj) & res_df$padj < 0.01),
            sum(!is.na(res_df$padj) & res_df$padj < 0.001),
            round(
                100 *
                sum(!is.na(res_df$padj) & res_df$padj < 0.05) /
                nrow(res_df),
                3
            )
        )
    )

    write.table(
        summary_table,
        file.path(
            output_dir,
            paste0(assay_name, "_LRT_summary.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n", assay_name, " LRT summary:\n", sep = "")
    print(summary_table)

    return(
        list(
            dds = dds,
            result = res_df,
            significant = sig,
            tissue_mean = tissue_mean,
            summary = summary_table
        )
    )
}

# ============================================================
# RNA
# ============================================================

rna_result <- run_lrt(
    count_matrix = rna,
    metadata = rna_meta,
    assay_name = "RNA",
    output_dir = file.path(outdir, "RNA"),
    heatmap_n = 200
)

# ============================================================
# ATAC
# ============================================================

atac_result <- run_lrt(
    count_matrix = atac,
    metadata = atac_meta,
    assay_name = "ATAC",
    output_dir = file.path(outdir, "ATAC"),
    heatmap_n = 500
)

# ============================================================
# Combined summary
# ============================================================

combined <- data.frame(
    Assay = c("RNA", "ATAC"),
    Input_features = c(
        nrow(rna),
        nrow(atac)
    ),
    Filtered_features = c(
        nrow(rna_result$result),
        nrow(atac_result$result)
    ),
    Tissue_dependent_FDR005 = c(
        nrow(rna_result$significant),
        nrow(atac_result$significant)
    ),
    Tissue_dependent_percent = c(
        round(
            100 *
            nrow(rna_result$significant) /
            nrow(rna_result$result),
            3
        ),
        round(
            100 *
            nrow(atac_result$significant) /
            nrow(atac_result$result),
            3
        )
    )
)

write.table(
    combined,
    file.path(outdir, "RNA_ATAC_multi_tissue_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n============================================================\n")
cat("STEP 02 COMPLETED\n")
cat("============================================================\n")

print(combined)

cat("\nImportant interpretation:\n")
cat(
    "The LRT tests whether Tissue contributes to expression/accessibility\n",
    "after controlling for Animal (P348/P350).\n",
    "The LRT log2FoldChange column is NOT an omnibus tissue effect size\n",
    "and should not be used to identify the tissue of maximum activity.\n",
    "Tissue specificity will be determined separately in Step 03.\n",
    sep = ""
)

cat("\n============================================================\n")

sessionInfo()
