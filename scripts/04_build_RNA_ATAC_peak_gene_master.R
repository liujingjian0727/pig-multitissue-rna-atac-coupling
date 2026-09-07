#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# ============================================================
# INPUT
# ============================================================

atac_spec_file <- paste0(
    "03_tissue_specificity/ATAC/",
    "ATAC_tissue_specificity_all.tsv"
)

rna_spec_file <- paste0(
    "03_tissue_specificity/RNA/",
    "RNA_tissue_specificity_all.tsv"
)

atac_spm_file <- paste0(
    "03B_SPM_analysis/ATAC/",
    "ATAC_SPM_all.tsv"
)

rna_spm_file <- paste0(
    "03B_SPM_analysis/RNA/",
    "RNA_SPM_all.tsv"
)

bridge_file <- "04_SAF_to_BED_peak_coordinate_bridge.tsv"

uropa_file <- "pig_8tissues_master_peaks_finalhits.txt"

outdir <- "04_RNA_ATAC_peak_gene_master"
dir.create(outdir, showWarnings = FALSE)

cat("============================================================\n")
cat("STEP 04: BUILD RNA-ATAC PEAK-GENE MASTER TABLE\n")
cat("============================================================\n\n")

# ============================================================
# READ
# ============================================================

atac <- read.delim(
    atac_spec_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

rna <- read.delim(
    rna_spec_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

atac_spm <- read.delim(
    atac_spm_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

rna_spm <- read.delim(
    rna_spm_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

bridge <- read.delim(
    bridge_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

uropa <- read.delim(
    uropa_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("ATAC specificity :", nrow(atac), "\n")
cat("RNA specificity  :", nrow(rna), "\n")
cat("ATAC SPM         :", nrow(atac_spm), "\n")
cat("RNA SPM          :", nrow(rna_spm), "\n")
cat("Bridge           :", nrow(bridge), "\n")
cat("UROPA            :", nrow(uropa), "\n\n")

# ============================================================
# UNIQUENESS QC
# ============================================================

stopifnot(!anyDuplicated(atac$Feature_ID))
stopifnot(!anyDuplicated(rna$Feature_ID))
stopifnot(!anyDuplicated(atac_spm$Feature_ID))
stopifnot(!anyDuplicated(rna_spm$Feature_ID))
stopifnot(!anyDuplicated(bridge$SAF_peak_id))
stopifnot(!anyDuplicated(bridge$BED_peak_id))
stopifnot(!anyDuplicated(uropa$peak_id))

cat("PASS: all primary identifiers are unique.\n")

# ============================================================
# 1. ATAC -> SAF/BED bridge
# ============================================================

bridge_idx <- match(
    atac$Feature_ID,
    bridge$SAF_peak_id
)

if (anyNA(bridge_idx)) {
    stop("ATAC features missing from SAF/BED bridge.")
}

master <- data.frame(
    SAF_peak_id = atac$Feature_ID,
    BED_peak_id = bridge$BED_peak_id[bridge_idx],

    peak_chr = bridge$Chr[bridge_idx],
    SAF_start = bridge$SAF_start[bridge_idx],
    SAF_end = bridge$SAF_end[bridge_idx],
    BED_start = bridge$BED_start[bridge_idx],
    BED_end = bridge$BED_end[bridge_idx],

    stringsAsFactors = FALSE
)

# ============================================================
# 2. Add ATAC tissue-specificity
# IMPORTANT:
# use stored Step03 classification directly.
# Do NOT re-apply Tau >= 0.8 after reading TSV.
# ============================================================

master$ATAC_Tau <- atac$Tau
master$ATAC_Tau_P348 <- atac$Tau_P348
master$ATAC_Tau_P350 <- atac$Tau_P350

master$ATAC_Max_tissue <- atac$Max_tissue
master$ATAC_Max_mean_TMM_TPM <- atac$Max_mean_TMM_TPM
master$ATAC_Second_mean_TMM_TPM <- atac$Second_mean_TMM_TPM
master$ATAC_Dominance_ratio <- atac$Dominance_ratio
master$ATAC_Max_fraction <- atac$Max_fraction

master$ATAC_P348_max_tissue <- atac$P348_max_tissue
master$ATAC_P350_max_tissue <- atac$P350_max_tissue
master$ATAC_Replicate_same_max <- atac$Replicate_same_max

master$ATAC_LRT_padj <- atac$LRT_padj

master$ATAC_Primary_tissue_specific <-
    atac$Primary_tissue_specific

master$ATAC_High_confidence_tissue_specific <-
    atac$High_confidence_tissue_specific

# ============================================================
# 3. Add ATAC SPM
# ============================================================

aidx <- match(
    master$SAF_peak_id,
    atac_spm$Feature_ID
)

if (anyNA(aidx)) {
    stop("ATAC peaks missing from ATAC SPM table.")
}

spm_cols <- paste0(
    "SPM_",
    c(
        "Adipose",
        "Cerebellum",
        "Cortex",
        "Hypothalamus",
        "Liver",
        "Lung",
        "Muscle",
        "Spleen"
    )
)

for (cc in spm_cols) {

    master[[paste0("ATAC_", cc)]] <-
        atac_spm[[cc]][aidx]
}

master$ATAC_Max_SPM <-
    atac_spm$Max_SPM[aidx]

master$ATAC_Second_SPM <-
    atac_spm$Second_SPM[aidx]

master$ATAC_SPM_margin <-
    atac_spm$SPM_margin[aidx]

master$ATAC_SPM_ratio <-
    atac_spm$SPM_ratio[aidx]

# ============================================================
# 4. Add UROPA annotation using BED coordinates
# ============================================================

uidx <- match(
    master$BED_peak_id,
    uropa$peak_id
)

if (anyNA(uidx)) {

    cat(
        "Missing UROPA annotations:",
        sum(is.na(uidx)),
        "\n"
    )

    stop("BED peak IDs are missing from UROPA.")
}

cat(
    "PASS: all",
    nrow(master),
    "BED peaks matched UROPA exactly.\n"
)

# Preserve available UROPA fields
uropa_fields <- c(
    "feature",
    "feat_start",
    "feat_end",
    "feat_strand",
    "feat_anchor",
    "distance",
    "relative_location",
    "feat_ovl_peak",
    "peak_ovl_feat",
    "gene_id",
    "name"
)

for (cc in uropa_fields) {

    if (cc %in% colnames(uropa)) {

        master[[paste0("UROPA_", cc)]] <-
            uropa[[cc]][uidx]
    }
}

# ============================================================
# 5. Standardize gene ID
# ============================================================

master$gene_id <- master$UROPA_gene_id

master$gene_id[
    is.na(master$gene_id) |
    master$gene_id == "" |
    master$gene_id == "NA" |
    master$gene_id == "."
] <- NA

master$Has_gene_annotation <- !is.na(master$gene_id)

# ============================================================
# 6. Link to RNA table
# ============================================================

ridx <- match(
    master$gene_id,
    rna$Feature_ID
)

master$Gene_present_in_RNA <- !is.na(ridx)

master$RNA_Tau <- rna$Tau[ridx]
master$RNA_Tau_P348 <- rna$Tau_P348[ridx]
master$RNA_Tau_P350 <- rna$Tau_P350[ridx]

master$RNA_Max_tissue <- rna$Max_tissue[ridx]

master$RNA_Max_mean_TMM_TPM <-
    rna$Max_mean_TMM_TPM[ridx]

master$RNA_Dominance_ratio <-
    rna$Dominance_ratio[ridx]

master$RNA_Max_fraction <-
    rna$Max_fraction[ridx]

master$RNA_Replicate_same_max <-
    rna$Replicate_same_max[ridx]

master$RNA_LRT_padj <-
    rna$LRT_padj[ridx]

master$RNA_Primary_tissue_specific <-
    rna$Primary_tissue_specific[ridx]

master$RNA_High_confidence_tissue_specific <-
    rna$High_confidence_tissue_specific[ridx]

# ============================================================
# 7. Add RNA SPM
# ============================================================

rsidx <- match(
    master$gene_id,
    rna_spm$Feature_ID
)

for (cc in spm_cols) {

    master[[paste0("RNA_", cc)]] <-
        rna_spm[[cc]][rsidx]
}

master$RNA_Max_SPM <-
    rna_spm$Max_SPM[rsidx]

master$RNA_Second_SPM <-
    rna_spm$Second_SPM[rsidx]

master$RNA_SPM_margin <-
    rna_spm$SPM_margin[rsidx]

master$RNA_SPM_ratio <-
    rna_spm$SPM_ratio[rsidx]

# ============================================================
# 8. RNA-ATAC concordance
# ============================================================

master$ATAC_RNA_same_max_tissue <- (
    !is.na(master$RNA_Max_tissue) &
    master$ATAC_Max_tissue ==
    master$RNA_Max_tissue
)

master$ATAC_RNA_same_max_tissue[
    is.na(master$ATAC_RNA_same_max_tissue)
] <- FALSE

master$Concordant_HC_tissue_specific_pair <- (
    master$ATAC_High_confidence_tissue_specific == TRUE &
    master$RNA_High_confidence_tissue_specific == TRUE &
    master$ATAC_RNA_same_max_tissue == TRUE
)

master$Concordant_HC_tissue_specific_pair[
    is.na(master$Concordant_HC_tissue_specific_pair)
] <- FALSE

# ============================================================
# 9. WRITE FULL MASTER TABLE
# ============================================================

write.table(
    master,
    file.path(
        outdir,
        "04_peak_gene_RNA_ATAC_master.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# 10. High-confidence ATAC subset
# ============================================================

hc_atac <- master[
    master$ATAC_High_confidence_tissue_specific == TRUE,
]

write.table(
    hc_atac,
    file.path(
        outdir,
        "04_HC_ATAC_peak_gene_master.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# 11. Concordant HC RNA-ATAC pairs
# ============================================================

concordant <- master[
    master$Concordant_HC_tissue_specific_pair == TRUE,
]

write.table(
    concordant,
    file.path(
        outdir,
        "04_concordant_HC_RNA_ATAC_peak_gene_pairs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# 12. Overall QC summary
# ============================================================

summary <- data.frame(
    Metric = c(
        "Total_ATAC_peaks",
        "UROPA_gene_assigned",
        "UROPA_no_gene",
        "Gene_present_in_RNA",
        "HC_ATAC_peaks",
        "HC_ATAC_with_gene",
        "HC_ATAC_gene_present_in_RNA",
        "HC_RNA_genes",
        "Concordant_HC_peak_gene_pairs",
        "Concordant_HC_unique_genes"
    ),
    Value = c(
        nrow(master),

        sum(
            master$Has_gene_annotation,
            na.rm = TRUE
        ),

        sum(
            !master$Has_gene_annotation,
            na.rm = TRUE
        ),

        sum(
            master$Gene_present_in_RNA,
            na.rm = TRUE
        ),

        nrow(hc_atac),

        sum(
            hc_atac$Has_gene_annotation,
            na.rm = TRUE
        ),

        sum(
            hc_atac$Gene_present_in_RNA,
            na.rm = TRUE
        ),

        sum(
            rna$High_confidence_tissue_specific == TRUE,
            na.rm = TRUE
        ),

        nrow(concordant),

        length(
            unique(
                concordant$gene_id[
                    !is.na(concordant$gene_id)
                ]
            )
        )
    )
)

write.table(
    summary,
    file.path(
        outdir,
        "04_master_QC_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 13. Tissue-level summary
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

tissue_summary <- do.call(
    rbind,
    lapply(
        tissues,
        function(tissue) {

            z <- (
                master$ATAC_High_confidence_tissue_specific == TRUE &
                master$ATAC_Max_tissue == tissue
            )

            z[is.na(z)] <- FALSE

            cpair <- (
                z &
                master$Concordant_HC_tissue_specific_pair == TRUE
            )

            cpair[is.na(cpair)] <- FALSE

            data.frame(
                Tissue = tissue,

                HC_ATAC_peaks =
                    sum(z),

                HC_ATAC_with_gene =
                    sum(
                        z &
                        master$Has_gene_annotation
                    ),

                HC_ATAC_gene_in_RNA =
                    sum(
                        z &
                        master$Gene_present_in_RNA
                    ),

                Concordant_HC_peak_gene_pairs =
                    sum(cpair),

                Concordant_HC_unique_genes =
                    length(
                        unique(
                            master$gene_id[
                                cpair &
                                !is.na(master$gene_id)
                            ]
                        )
                    )
            )
        }
    )
)

write.table(
    tissue_summary,
    file.path(
        outdir,
        "04_tissue_level_integration_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("STEP 04 MASTER TABLE COMPLETED\n")
cat("============================================================\n\n")

print(summary)

cat("\nTissue-level summary:\n")
print(tissue_summary)

cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("\nIMPORTANT:\n")
cat(
    "ATAC high-confidence membership is inherited directly\n",
    "from Step03 and was NOT recalculated from Tau after TSV\n",
    "serialization.\n",
    sep = ""
)

cat("============================================================\n")
