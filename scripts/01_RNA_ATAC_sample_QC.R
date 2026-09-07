#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 0. Input
# ============================================================

rna_count_file <- "sus.counts.matrix"
atac_count_file <- "pig.counts.matrix"
rna_tpm_file <- "sus.TMM.TPM.matrix"
atac_tpm_file <- "pig.TMM.TPM.matrix"
meta_file <- "sample_pair_metadata.tsv"

outdir <- "01_sample_QC"

dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "RNA"), showWarnings = FALSE)
dir.create(file.path(outdir, "ATAC"), showWarnings = FALSE)

cat("============================================================\n")
cat("Pig paired RNA-seq / ATAC-seq sample QC\n")
cat("============================================================\n\n")

# ============================================================
# 1. Read data
# ============================================================

read_matrix <- function(file) {

    x <- read.delim(
        file,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        sep = "\t",
        quote = "",
        comment.char = ""
    )

    x <- as.matrix(x)

    storage.mode(x) <- "numeric"

    return(x)
}

cat("Reading metadata...\n")

meta <- read.delim(
    meta_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

cat("Reading RNA counts...\n")
rna_counts <- read_matrix(rna_count_file)

cat("Reading ATAC counts...\n")
atac_counts <- read_matrix(atac_count_file)

cat("Reading RNA TMM.TPM...\n")
rna_tpm <- read_matrix(rna_tpm_file)

cat("Reading ATAC TMM.TPM...\n")
atac_tpm <- read_matrix(atac_tpm_file)

cat("\nInput dimensions:\n")
cat("RNA counts :", nrow(rna_counts), "features x",
    ncol(rna_counts), "samples\n")
cat("RNA TPM    :", nrow(rna_tpm), "features x",
    ncol(rna_tpm), "samples\n")
cat("ATAC counts:", nrow(atac_counts), "features x",
    ncol(atac_counts), "samples\n")
cat("ATAC TPM   :", nrow(atac_tpm), "features x",
    ncol(atac_tpm), "samples\n")

# ============================================================
# 2. Basic sanity checks
# ============================================================

stopifnot(nrow(meta) == 16)

if (!identical(colnames(rna_counts), meta$RNA_SRR)) {
    stop("RNA count sample order does not match metadata.")
}

if (!identical(colnames(atac_counts), meta$ATAC_SRR)) {
    stop("ATAC count sample order does not match metadata.")
}

if (!identical(colnames(rna_tpm), meta$RNA_SRR)) {
    stop("RNA TPM sample order does not match metadata.")
}

if (!identical(colnames(atac_tpm), meta$ATAC_SRR)) {
    stop("ATAC TPM sample order does not match metadata.")
}

if (!identical(rownames(rna_counts), rownames(rna_tpm))) {
    stop("RNA count and TPM feature IDs/order are different.")
}

if (!identical(rownames(atac_counts), rownames(atac_tpm))) {
    stop("ATAC count and TPM feature IDs/order are different.")
}

if (anyNA(rna_counts) || anyNA(atac_counts)) {
    stop("NA detected in count matrices.")
}

if (any(rna_counts < 0) || any(atac_counts < 0)) {
    stop("Negative counts detected.")
}

if (any(abs(rna_counts - round(rna_counts)) > 1e-8)) {
    stop("RNA count matrix contains non-integer values.")
}

if (any(abs(atac_counts - round(atac_counts)) > 1e-8)) {
    stop("ATAC count matrix contains non-integer values.")
}

cat("\nPASS: matrix structure and metadata matching are correct.\n")

# ============================================================
# 3. Create metadata for each assay
# ============================================================

rna_meta <- data.frame(
    SRR = meta$RNA_SRR,
    Tissue = factor(
        meta$Tissue,
        levels = c(
            "Adipose",
            "Cerebellum",
            "Cortex",
            "Hypothalamus",
            "Liver",
            "Lung",
            "Muscle",
            "Spleen"
        )
    ),
    Animal = factor(meta$Animal),
    Label = paste(meta$Tissue, meta$Animal, sep = "_"),
    row.names = meta$RNA_SRR,
    check.names = FALSE
)

atac_meta <- data.frame(
    SRR = meta$ATAC_SRR,
    Tissue = factor(
        meta$Tissue,
        levels = c(
            "Adipose",
            "Cerebellum",
            "Cortex",
            "Hypothalamus",
            "Liver",
            "Lung",
            "Muscle",
            "Spleen"
        )
    ),
    Animal = factor(meta$Animal),
    Label = paste(meta$Tissue, meta$Animal, sep = "_"),
    row.names = meta$ATAC_SRR,
    check.names = FALSE
)

write.table(
    rna_meta,
    file.path(outdir, "RNA", "RNA_metadata_used.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

write.table(
    atac_meta,
    file.path(outdir, "ATAC", "ATAC_metadata_used.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

# ============================================================
# 4. General QC function
# ============================================================

run_qc <- function(counts, metadata, assay_name, assay_dir) {

    cat("\n============================================================\n")
    cat("Processing:", assay_name, "\n")
    cat("============================================================\n")

    # --------------------------------------------------------
    # Library size
    # --------------------------------------------------------

    lib_size <- colSums(counts)

    lib_table <- data.frame(
        Sample = colnames(counts),
        Tissue = metadata[colnames(counts), "Tissue"],
        Animal = metadata[colnames(counts), "Animal"],
        Library_size = lib_size,
        Million_counts = lib_size / 1e6
    )

    write.table(
        lib_table,
        file.path(assay_dir, paste0(assay_name, "_library_size.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    p_lib <- ggplot(
        lib_table,
        aes(
            x = reorder(paste(Tissue, Animal, sep = "_"), Million_counts),
            y = Million_counts,
            fill = Tissue
        )
    ) +
        geom_col() +
        coord_flip() +
        labs(
            x = NULL,
            y = "Library size (million counts)",
            title = paste0(assay_name, " library size")
        ) +
        theme_bw(base_size = 12) +
        theme(
            legend.position = "none",
            panel.grid.major.y = element_blank()
        )

    ggsave(
        file.path(assay_dir, paste0(assay_name, "_library_size.pdf")),
        p_lib,
        width = 7,
        height = 6
    )

    # --------------------------------------------------------
    # Filtering
    # count >=10 in >=2 samples
    # --------------------------------------------------------

    keep <- rowSums(counts >= 10) >= 2

    counts_filtered <- counts[keep, , drop = FALSE]

    filtering <- data.frame(
        Metric = c(
            "Features_before_filtering",
            "Features_after_filtering",
            "Features_removed",
            "Retention_percent"
        ),
        Value = c(
            nrow(counts),
            nrow(counts_filtered),
            nrow(counts) - nrow(counts_filtered),
            round(100 * nrow(counts_filtered) / nrow(counts), 3)
        )
    )

    write.table(
        filtering,
        file.path(assay_dir, paste0(assay_name, "_filtering_summary.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("Features before filtering:", nrow(counts), "\n")
    cat("Features after filtering :", nrow(counts_filtered), "\n")

    # --------------------------------------------------------
    # DESeq2 VST
    # --------------------------------------------------------

    count_integer <- round(counts_filtered)

    dds <- DESeqDataSetFromMatrix(
        countData = count_integer,
        colData = metadata,
        design = ~ Animal + Tissue
    )

    cat("Running variance stabilizing transformation...\n")

    vsd <- vst(dds, blind = TRUE)

    vst_mat <- assay(vsd)

    write.table(
        vst_mat,
        file.path(assay_dir, paste0(assay_name, "_VST_matrix.tsv")),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    # --------------------------------------------------------
    # PCA
    # --------------------------------------------------------

    pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

    variance <- (pca$sdev^2) / sum(pca$sdev^2)

    pca_table <- data.frame(
        Sample = rownames(pca$x),
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        Tissue = metadata[rownames(pca$x), "Tissue"],
        Animal = metadata[rownames(pca$x), "Animal"],
        Label = metadata[rownames(pca$x), "Label"]
    )

    write.table(
        pca_table,
        file.path(assay_dir, paste0(assay_name, "_PCA_coordinates.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    variance_table <- data.frame(
        PC = paste0("PC", seq_along(variance)),
        Variance_percent = 100 * variance
    )

    write.table(
        variance_table,
        file.path(assay_dir, paste0(assay_name, "_PCA_variance.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    p_pca <- ggplot(
        pca_table,
        aes(
            x = PC1,
            y = PC2,
            color = Tissue,
            shape = Animal
        )
    ) +
        geom_point(size = 4) +
        geom_text(
            aes(label = Label),
            size = 3,
            vjust = -0.8,
            show.legend = FALSE
        ) +
        labs(
            x = sprintf("PC1 (%.2f%%)", 100 * variance[1]),
            y = sprintf("PC2 (%.2f%%)", 100 * variance[2]),
            title = paste0(assay_name, " PCA")
        ) +
        theme_bw(base_size = 12)

    ggsave(
        file.path(assay_dir, paste0(assay_name, "_PCA.pdf")),
        p_pca,
        width = 8,
        height = 7
    )

    # --------------------------------------------------------
    # Correlation matrix
    # --------------------------------------------------------

    sample_cor <- cor(
        vst_mat,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    write.table(
        sample_cor,
        file.path(assay_dir, paste0(assay_name, "_sample_correlation.tsv")),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    annotation_col <- data.frame(
        Tissue = metadata[colnames(sample_cor), "Tissue"],
        Animal = metadata[colnames(sample_cor), "Animal"]
    )

    rownames(annotation_col) <- metadata[colnames(sample_cor), "Label"]

    cor_plot <- sample_cor

    new_labels <- metadata[colnames(sample_cor), "Label"]

    colnames(cor_plot) <- new_labels
    rownames(cor_plot) <- new_labels

    pdf(
        file.path(
            assay_dir,
            paste0(assay_name, "_sample_correlation_heatmap.pdf")
        ),
        width = 10,
        height = 9
    )

    pheatmap(
        cor_plot,
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "correlation",
        annotation_col = annotation_col,
        border_color = NA,
        fontsize = 8,
        main = paste0(assay_name, " sample correlation")
    )

    dev.off()

    # --------------------------------------------------------
    # Hierarchical clustering
    # --------------------------------------------------------

    sample_dist <- dist(t(vst_mat), method = "euclidean")

    hc <- hclust(sample_dist, method = "complete")

    hc$labels <- metadata[hc$labels, "Label"]

    pdf(
        file.path(
            assay_dir,
            paste0(assay_name, "_hierarchical_clustering.pdf")
        ),
        width = 10,
        height = 6
    )

    plot(
        hc,
        main = paste0(assay_name, " hierarchical clustering"),
        xlab = "",
        sub = "",
        cex = 0.8
    )

    dev.off()

    # --------------------------------------------------------
    # P348 vs P350 replicate correlation within each tissue
    # --------------------------------------------------------

    tissues <- levels(metadata$Tissue)

    replicate_results <- list()

    for (tissue in tissues) {

        samples <- rownames(metadata)[metadata$Tissue == tissue]

        if (length(samples) != 2) next

        x <- vst_mat[, samples[1]]
        y <- vst_mat[, samples[2]]

        replicate_results[[tissue]] <- data.frame(
            Tissue = tissue,
            Sample1 = samples[1],
            Sample2 = samples[2],
            Pearson = cor(x, y, method = "pearson"),
            Spearman = cor(x, y, method = "spearman"),
            N_features = length(x)
        )
    }

    replicate_table <- do.call(rbind, replicate_results)

    write.table(
        replicate_table,
        file.path(
            assay_dir,
            paste0(assay_name, "_P348_P350_replicate_correlations.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n", assay_name, " PCA variance:\n", sep = "")
    cat(
        sprintf(
            "PC1 = %.2f%%, PC2 = %.2f%%\n",
            100 * variance[1],
            100 * variance[2]
        )
    )

    cat("\nP348/P350 replicate correlations:\n")
    print(replicate_table)

    return(
        list(
            filtered_counts = counts_filtered,
            vst = vst_mat,
            pca = pca_table,
            correlations = replicate_table
        )
    )
}

# ============================================================
# 5. RNA QC
# ============================================================

rna_result <- run_qc(
    counts = rna_counts,
    metadata = rna_meta,
    assay_name = "RNA",
    assay_dir = file.path(outdir, "RNA")
)

# ============================================================
# 6. ATAC QC
# ============================================================

atac_result <- run_qc(
    counts = atac_counts,
    metadata = atac_meta,
    assay_name = "ATAC",
    assay_dir = file.path(outdir, "ATAC")
)

# ============================================================
# 7. Final summary
# ============================================================

cat("\n============================================================\n")
cat("STEP 01 SAMPLE QC COMPLETED\n")
cat("============================================================\n")

cat("\nRNA:\n")
cat("Original genes :", nrow(rna_counts), "\n")
cat("Filtered genes :", nrow(rna_result$filtered_counts), "\n")

cat("\nATAC:\n")
cat("Original peaks :", nrow(atac_counts), "\n")
cat("Filtered peaks :", nrow(atac_result$filtered_counts), "\n")

cat("\nOutput directory:", outdir, "\n")
cat("============================================================\n")

sessionInfo()
