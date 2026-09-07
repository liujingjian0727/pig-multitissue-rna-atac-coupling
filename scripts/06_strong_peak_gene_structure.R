#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# STEP 06
# STRUCTURAL CHARACTERIZATION OF 702 STRONGLY COUPLED
# TISSUE-SPECIFIC RNA-ATAC PEAK-GENE ASSOCIATIONS
#
# Strong criterion inherited from Step05/05B:
# rho_P348 >= 0.5
# rho_P350 >= 0.5
# rho_mean >= 0.7
#
# IMPORTANT:
# Membership is NOT recalculated here.
# ============================================================


# ============================================================
# INPUT
# ============================================================

input_file <- paste0(
    "05B_same_tissue_permutation_v2/",
    "05B_observed_strong_05_05_07_pairs.tsv"
)

outdir <- "06_strong_peak_gene_structure"

dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
)


# ============================================================
# READ
# ============================================================

x <- read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("============================================================\n")
cat("STEP 06: STRONG PEAK-GENE STRUCTURAL ANALYSIS\n")
cat("============================================================\n\n")

cat("Input strong associations :", nrow(x), "\n")

if (nrow(x) != 702) {
    warning(
        "Expected 702 strong associations, but found ",
        nrow(x)
    )
}


# ============================================================
# REQUIRED COLUMNS
# ============================================================

required <- c(
    "SAF_peak_id",
    "BED_peak_id",
    "gene_id",
    "ATAC_Max_tissue",
    "RNA_Max_tissue",
    "UROPA_relative_location",
    "UROPA_distance",
    "rho_P348",
    "rho_P350",
    "rho_mean"
)

missing_cols <- setdiff(
    required,
    colnames(x)
)

if (length(missing_cols) > 0) {
    stop(
        "Missing required columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


# ============================================================
# BASIC QC
# ============================================================

cat("\nBasic QC:\n")

cat(
    "Unique SAF peaks :",
    length(unique(x$SAF_peak_id)),
    "\n"
)

cat(
    "Unique BED peaks :",
    length(unique(x$BED_peak_id)),
    "\n"
)

cat(
    "Unique genes     :",
    length(unique(x$gene_id)),
    "\n"
)

cat(
    "Duplicated SAF peak IDs:",
    sum(duplicated(x$SAF_peak_id)),
    "\n"
)

cat(
    "Duplicated BED peak IDs:",
    sum(duplicated(x$BED_peak_id)),
    "\n"
)

same_tissue <- (
    x$ATAC_Max_tissue ==
    x$RNA_Max_tissue
)

cat(
    "ATAC/RNA same tissue:",
    sum(same_tissue, na.rm = TRUE),
    "/",
    nrow(x),
    "\n"
)

if (!all(same_tissue)) {
    stop(
        "Strong input unexpectedly contains ",
        "ATAC/RNA tissue mismatches."
    )
}


# ============================================================
# STANDARDIZE UROPA FIELDS
# ============================================================

x$UROPA_relative_location[
    is.na(x$UROPA_relative_location) |
    x$UROPA_relative_location == ""
] <- "NA"

x$UROPA_distance_numeric <- suppressWarnings(
    as.numeric(
        x$UROPA_distance
    )
)

x$UROPA_abs_distance <- abs(
    x$UROPA_distance_numeric
)


# ============================================================
# STRUCTURAL CLASS
#
# Current UROPA annotation is gene-level.
#
# Gene-overlapping:
#   PeakInsideFeature
#   FeatureInsidePeak
#   OverlapStart
#   OverlapEnd
#
# Flanking:
#   Upstream / Downstream
#
# Distances below refer to UROPA feature-anchor distance,
# NOT necessarily distance to TSS.
# ============================================================

overlap_classes <- c(
    "PeakInsideFeature",
    "FeatureInsidePeak",
    "OverlapStart",
    "OverlapEnd"
)

flanking_classes <- c(
    "Upstream",
    "Downstream"
)

x$Structural_class <- "Other_or_NA"

x$Structural_class[
    x$UROPA_relative_location %in%
        overlap_classes
] <- "Gene_overlapping"


is_flanking <- (
    x$UROPA_relative_location %in%
        flanking_classes &
    is.finite(x$UROPA_abs_distance)
)


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance <= 2000
] <- "Flanking_0_2kb"


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance > 2000 &
    x$UROPA_abs_distance <= 20000
] <- "Flanking_2_20kb"


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance > 20000 &
    x$UROPA_abs_distance <= 100000
] <- "Flanking_20_100kb"


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance > 100000
] <- "Flanking_gt100kb"


class_levels <- c(
    "Gene_overlapping",
    "Flanking_0_2kb",
    "Flanking_2_20kb",
    "Flanking_20_100kb",
    "Flanking_gt100kb",
    "Other_or_NA"
)

x$Structural_class <- factor(
    x$Structural_class,
    levels = class_levels
)


# ============================================================
# WRITE ANNOTATED 702-PAIR MASTER
# ============================================================

write.table(
    x,
    file.path(
        outdir,
        "06_702_strong_pairs_with_structure.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 1. OVERALL SUMMARY
# ============================================================

overall_summary <- data.frame(
    Metric = c(
        "Strong_peak_gene_pairs",
        "Unique_SAF_peaks",
        "Unique_BED_peaks",
        "Unique_genes",
        "Mean_peaks_per_gene",
        "Median_peaks_per_gene",
        "Maximum_peaks_per_gene"
    ),

    Value = c(
        nrow(x),
        length(unique(x$SAF_peak_id)),
        length(unique(x$BED_peak_id)),
        length(unique(x$gene_id)),
        NA,
        NA,
        NA
    )
)


gene_peak_counts <- as.data.frame(
    table(
        x$gene_id
    ),
    stringsAsFactors = FALSE
)

colnames(gene_peak_counts) <- c(
    "gene_id",
    "N_strong_peaks"
)

gene_peak_counts$N_strong_peaks <-
    as.integer(
        gene_peak_counts$N_strong_peaks
    )


overall_summary$Value[
    overall_summary$Metric ==
        "Mean_peaks_per_gene"
] <- mean(
    gene_peak_counts$N_strong_peaks
)


overall_summary$Value[
    overall_summary$Metric ==
        "Median_peaks_per_gene"
] <- median(
    gene_peak_counts$N_strong_peaks
)


overall_summary$Value[
    overall_summary$Metric ==
        "Maximum_peaks_per_gene"
] <- max(
    gene_peak_counts$N_strong_peaks
)


write.table(
    overall_summary,
    file.path(
        outdir,
        "06_overall_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 2. GENE-LEVEL PEAK MULTIPLICITY
# ============================================================

# Add preferred tissue for each gene.
# Strong associations are already same-tissue,
# but theoretically the same gene should normally belong to one
# tissue-specific category.

gene_tissue <- tapply(
    x$ATAC_Max_tissue,
    x$gene_id,
    function(z) {
        paste(
            sort(unique(z)),
            collapse = ";"
        )
    }
)

gene_peak_counts$Tissue <-
    gene_tissue[
        gene_peak_counts$gene_id
    ]


# Mean / maximum coupling per gene

gene_median_rho <- tapply(
    x$rho_mean,
    x$gene_id,
    median,
    na.rm = TRUE
)

gene_max_rho <- tapply(
    x$rho_mean,
    x$gene_id,
    max,
    na.rm = TRUE
)

gene_peak_counts$Median_rho_mean <-
    gene_median_rho[
        gene_peak_counts$gene_id
    ]

gene_peak_counts$Max_rho_mean <-
    gene_max_rho[
        gene_peak_counts$gene_id
    ]


# Count structural classes per gene

for (cc in class_levels) {

    temp <- tapply(
        x$Structural_class == cc,
        x$gene_id,
        sum,
        na.rm = TRUE
    )

    gene_peak_counts[[paste0("N_", cc)]] <- temp[
        gene_peak_counts$gene_id
    ]
}


gene_peak_counts <- gene_peak_counts[
    order(
        -gene_peak_counts$N_strong_peaks,
        -gene_peak_counts$Max_rho_mean,
        gene_peak_counts$gene_id
    ),
]


write.table(
    gene_peak_counts,
    file.path(
        outdir,
        "06_gene_strong_peak_multiplicity.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# Peak-number bins per gene

gene_peak_counts$Peak_number_class <- cut(
    gene_peak_counts$N_strong_peaks,
    breaks = c(
        0,
        1,
        2,
        3,
        5,
        10,
        Inf
    ),
    labels = c(
        "1",
        "2",
        "3",
        "4-5",
        "6-10",
        ">10"
    ),
    right = TRUE
)


multiplicity_summary <- as.data.frame(
    table(
        gene_peak_counts$Peak_number_class
    )
)

colnames(multiplicity_summary) <- c(
    "Strong_peaks_per_gene",
    "Number_of_genes"
)

multiplicity_summary$Percent <- round(
    100 *
    multiplicity_summary$Number_of_genes /
    sum(
        multiplicity_summary$Number_of_genes
    ),
    3
)


write.table(
    multiplicity_summary,
    file.path(
        outdir,
        "06_gene_peak_multiplicity_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 3. UROPA RELATIVE LOCATION
# ============================================================

relative_summary <- as.data.frame(
    table(
        x$UROPA_relative_location
    )
)

colnames(relative_summary) <- c(
    "UROPA_relative_location",
    "Number_of_pairs"
)

relative_summary$Percent <- round(
    100 *
    relative_summary$Number_of_pairs /
    nrow(x),
    3
)

relative_summary <- relative_summary[
    order(
        -relative_summary$Number_of_pairs
    ),
]


write.table(
    relative_summary,
    file.path(
        outdir,
        "06_UROPA_relative_location_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 4. STRUCTURAL CLASS SUMMARY
# ============================================================

structure_summary <- as.data.frame(
    table(
        x$Structural_class
    )
)

colnames(structure_summary) <- c(
    "Structural_class",
    "Number_of_pairs"
)

structure_summary$Percent <- round(
    100 *
    structure_summary$Number_of_pairs /
    nrow(x),
    3
)


write.table(
    structure_summary,
    file.path(
        outdir,
        "06_structural_class_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 5. UPSTREAM / DOWNSTREAM DISTANCE SUMMARY
# ============================================================

flanking <- x[
    x$UROPA_relative_location %in%
        flanking_classes &
    is.finite(
        x$UROPA_abs_distance
    ),
    ,
    drop = FALSE
]


distance_summary <- do.call(
    rbind,
    lapply(
        c(
            "All_flanking",
            "Upstream",
            "Downstream"
        ),
        function(group) {

            if (group == "All_flanking") {

                z <- flanking$UROPA_abs_distance

            } else {

                z <- flanking$UROPA_abs_distance[
                    flanking$UROPA_relative_location ==
                        group
                ]
            }

            if (length(z) == 0) {

                return(
                    data.frame(
                        Group = group,
                        N = 0,
                        Median_distance = NA,
                        Q25_distance = NA,
                        Q75_distance = NA,
                        Min_distance = NA,
                        Max_distance = NA
                    )
                )
            }

            data.frame(
                Group = group,
                N = length(z),
                Median_distance =
                    median(z),

                Q25_distance =
                    as.numeric(
                        quantile(
                            z,
                            0.25,
                            names = FALSE
                        )
                    ),

                Q75_distance =
                    as.numeric(
                        quantile(
                            z,
                            0.75,
                            names = FALSE
                        )
                    ),

                Min_distance =
                    min(z),

                Max_distance =
                    max(z)
            )
        }
    )
)


write.table(
    distance_summary,
    file.path(
        outdir,
        "06_flanking_distance_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 6. TISSUE x STRUCTURE
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


tissue_structure <- as.data.frame(
    table(
        factor(
            x$ATAC_Max_tissue,
            levels = tissues
        ),
        x$Structural_class
    )
)

colnames(tissue_structure) <- c(
    "Tissue",
    "Structural_class",
    "Number_of_pairs"
)


tissue_total <- tapply(
    tissue_structure$Number_of_pairs,
    tissue_structure$Tissue,
    sum
)


tissue_structure$Percent_within_tissue <-
    round(
        100 *
        tissue_structure$Number_of_pairs /
        tissue_total[
            as.character(
                tissue_structure$Tissue
            )
        ],
        3
    )


write.table(
    tissue_structure,
    file.path(
        outdir,
        "06_tissue_by_structural_class.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 7. COUPLING BY STRUCTURAL CLASS
# ============================================================

safe_summary <- function(z) {

    z <- z[
        is.finite(z)
    ]

    if (length(z) == 0) {

        return(
            c(
                N = 0,
                Median = NA,
                Q25 = NA,
                Q75 = NA,
                Mean = NA
            )
        )
    }

    c(
        N = length(z),

        Median =
            median(z),

        Q25 =
            as.numeric(
                quantile(
                    z,
                    0.25,
                    names = FALSE
                )
            ),

        Q75 =
            as.numeric(
                quantile(
                    z,
                    0.75,
                    names = FALSE
                )
            ),

        Mean =
            mean(z)
    )
}


coupling_structure <- do.call(
    rbind,
    lapply(
        class_levels,
        function(cc) {

            z <- x[
                x$Structural_class == cc,
                ,
                drop = FALSE
            ]

            s348 <- safe_summary(
                z$rho_P348
            )

            s350 <- safe_summary(
                z$rho_P350
            )

            sm <- safe_summary(
                z$rho_mean
            )

            data.frame(
                Structural_class = cc,
                N_pairs = nrow(z),

                Median_rho_P348 =
                    s348["Median"],

                Median_rho_P350 =
                    s350["Median"],

                Median_rho_mean =
                    sm["Median"],

                Q25_rho_mean =
                    sm["Q25"],

                Q75_rho_mean =
                    sm["Q75"],

                Mean_rho_mean =
                    sm["Mean"]
            )
        }
    )
)


write.table(
    coupling_structure,
    file.path(
        outdir,
        "06_coupling_by_structural_class.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 8. COUPLING BY TISSUE
# ============================================================

coupling_tissue <- do.call(
    rbind,
    lapply(
        tissues,
        function(tissue) {

            z <- x[
                x$ATAC_Max_tissue ==
                    tissue,
                ,
                drop = FALSE
            ]

            sm <- safe_summary(
                z$rho_mean
            )

            data.frame(
                Tissue = tissue,
                N_pairs = nrow(z),
                Unique_genes =
                    length(
                        unique(
                            z$gene_id
                        )
                    ),

                Median_rho_mean =
                    sm["Median"],

                Q25_rho_mean =
                    sm["Q25"],

                Q75_rho_mean =
                    sm["Q75"]
            )
        }
    )
)


write.table(
    coupling_tissue,
    file.path(
        outdir,
        "06_coupling_by_tissue.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 9. TOP GENES BY NUMBER OF STRONG PEAKS
# ============================================================

top_genes <- head(
    gene_peak_counts,
    50
)


write.table(
    top_genes,
    file.path(
        outdir,
        "06_top50_genes_by_strong_peak_number.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 10. TOP PAIRS BY rho_mean
# ============================================================

top_pairs <- x[
    order(
        -x$rho_mean,
        -x$rho_P348,
        -x$rho_P350
    ),
]


write.table(
    head(
        top_pairs,
        100
    ),
    file.path(
        outdir,
        "06_top100_strong_pairs_by_rho_mean.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 11. PLOTS
# ============================================================

# ------------------------------------------------------------
# Relative-location bar plot
# ------------------------------------------------------------

p1 <- ggplot(
    relative_summary,
    aes(
        x = reorder(
            UROPA_relative_location,
            Number_of_pairs
        ),
        y = Number_of_pairs
    )
) +
    geom_col() +
    coord_flip() +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title =
            "UROPA relative locations of 702 strong peak-gene associations",
        x = NULL,
        y = "Number of peak-gene associations"
    )


ggsave(
    file.path(
        outdir,
        "06_UROPA_relative_location.pdf"
    ),
    p1,
    width = 7,
    height = 5
)


# ------------------------------------------------------------
# Structural-class bar plot
# ------------------------------------------------------------

p2 <- ggplot(
    structure_summary,
    aes(
        x = Structural_class,
        y = Number_of_pairs
    )
) +
    geom_col() +
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
            "Structural distribution of strongly coupled peak-gene associations",
        x = NULL,
        y = "Number of peak-gene associations"
    )


ggsave(
    file.path(
        outdir,
        "06_structural_class_counts.pdf"
    ),
    p2,
    width = 8,
    height = 5.5
)


# ------------------------------------------------------------
# Gene peak multiplicity
# ------------------------------------------------------------

p3 <- ggplot(
    multiplicity_summary,
    aes(
        x = Strong_peaks_per_gene,
        y = Number_of_genes
    )
) +
    geom_col() +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title =
            "Number of strongly coupled ATAC peaks per gene",
        x = "Strong ATAC peaks per gene",
        y = "Number of genes"
    )


ggsave(
    file.path(
        outdir,
        "06_gene_peak_multiplicity.pdf"
    ),
    p3,
    width = 6.5,
    height = 5
)


# ------------------------------------------------------------
# Tissue x structural class
# ------------------------------------------------------------

p4 <- ggplot(
    tissue_structure,
    aes(
        x = Tissue,
        y = Percent_within_tissue,
        fill = Structural_class
    )
) +
    geom_col() +
    theme_bw(
        base_size = 11
    ) +
    theme(
        axis.text.x =
            element_text(
                angle = 45,
                hjust = 1
            )
    ) +
    labs(
        title =
            "Structural composition of strong peak-gene associations by tissue",
        x = NULL,
        y = "Percentage within tissue",
        fill = "Structural class"
    )


ggsave(
    file.path(
        outdir,
        "06_tissue_structural_composition.pdf"
    ),
    p4,
    width = 9,
    height = 6
)


# ------------------------------------------------------------
# rho_mean by structural class
# ------------------------------------------------------------

plot_data <- x[
    is.finite(
        x$rho_mean
    ),
    ,
    drop = FALSE
]


p5 <- ggplot(
    plot_data,
    aes(
        x = Structural_class,
        y = rho_mean
    )
) +
    geom_boxplot(
        outlier.shape = NA
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
            "RNA-ATAC coupling by genomic structural class",
        x = NULL,
        y = expression(
            rho[mean]
        )
    )


ggsave(
    file.path(
        outdir,
        "06_rho_mean_by_structural_class.pdf"
    ),
    p5,
    width = 8,
    height = 5.5
)


# ------------------------------------------------------------
# Flanking distance distribution
# log10 only for visualization
# ------------------------------------------------------------

if (nrow(flanking) > 0) {

    flanking$log10_distance_plus1 <-
        log10(
            flanking$UROPA_abs_distance + 1
        )

    p6 <- ggplot(
        flanking,
        aes(
            x = log10_distance_plus1
        )
    ) +
        geom_histogram(
            bins = 50
        ) +
        facet_wrap(
            ~ UROPA_relative_location,
            ncol = 1
        ) +
        theme_bw(
            base_size = 11
        ) +
        labs(
            title =
                "Flanking distance distribution of strong peak-gene associations",
            subtitle =
                "Distance refers to UROPA feature-anchor distance, not necessarily TSS distance",
            x = "log10(|UROPA distance| + 1)",
            y = "Number of peak-gene associations"
        )

    ggsave(
        file.path(
            outdir,
            "06_flanking_distance_distribution.pdf"
        ),
        p6,
        width = 7,
        height = 7
    )
}


# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("STEP 06 COMPLETED\n")
cat("============================================================\n\n")

cat("OVERALL SUMMARY:\n")
print(
    overall_summary,
    row.names = FALSE
)

cat("\nGENE PEAK MULTIPLICITY:\n")
print(
    multiplicity_summary,
    row.names = FALSE
)

cat("\nUROPA RELATIVE LOCATION:\n")
print(
    relative_summary,
    row.names = FALSE
)

cat("\nSTRUCTURAL CLASS:\n")
print(
    structure_summary,
    row.names = FALSE
)

cat("\nFLANKING DISTANCE SUMMARY:\n")
print(
    distance_summary,
    row.names = FALSE
)

cat("\nCOUPLING BY STRUCTURAL CLASS:\n")
print(
    coupling_structure,
    row.names = FALSE
)

cat("\nCOUPLING BY TISSUE:\n")
print(
    coupling_tissue,
    row.names = FALSE
)

cat("\nTOP 20 GENES BY NUMBER OF STRONG PEAKS:\n")
print(
    head(
        gene_peak_counts,
        20
    ),
    row.names = FALSE
)

cat("\nIMPORTANT:\n")
cat(
    "UROPA distance in the current gene-level annotation is\n",
    "not assumed to be distance-to-TSS. Therefore structural\n",
    "classes use the neutral terms gene-overlapping and flanking.\n",
    sep = ""
)

cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("============================================================\n")

sessionInfo()

EOFcat > 06_strong_peak_gene_structure.R <<'EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# STEP 06
# STRUCTURAL CHARACTERIZATION OF 702 STRONGLY COUPLED
# TISSUE-SPECIFIC RNA-ATAC PEAK-GENE ASSOCIATIONS
#
# Strong criterion inherited from Step05/05B:
# rho_P348 >= 0.5
# rho_P350 >= 0.5
# rho_mean >= 0.7
#
# IMPORTANT:
# Membership is NOT recalculated here.
# ============================================================


# ============================================================
# INPUT
# ============================================================

input_file <- paste0(
    "05B_same_tissue_permutation_v2/",
    "05B_observed_strong_05_05_07_pairs.tsv"
)

outdir <- "06_strong_peak_gene_structure"

dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
)


# ============================================================
# READ
# ============================================================

x <- read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("============================================================\n")
cat("STEP 06: STRONG PEAK-GENE STRUCTURAL ANALYSIS\n")
cat("============================================================\n\n")

cat("Input strong associations :", nrow(x), "\n")

if (nrow(x) != 702) {
    warning(
        "Expected 702 strong associations, but found ",
        nrow(x)
    )
}


# ============================================================
# REQUIRED COLUMNS
# ============================================================

required <- c(
    "SAF_peak_id",
    "BED_peak_id",
    "gene_id",
    "ATAC_Max_tissue",
    "RNA_Max_tissue",
    "UROPA_relative_location",
    "UROPA_distance",
    "rho_P348",
    "rho_P350",
    "rho_mean"
)

missing_cols <- setdiff(
    required,
    colnames(x)
)

if (length(missing_cols) > 0) {
    stop(
        "Missing required columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


# ============================================================
# BASIC QC
# ============================================================

cat("\nBasic QC:\n")

cat(
    "Unique SAF peaks :",
    length(unique(x$SAF_peak_id)),
    "\n"
)

cat(
    "Unique BED peaks :",
    length(unique(x$BED_peak_id)),
    "\n"
)

cat(
    "Unique genes     :",
    length(unique(x$gene_id)),
    "\n"
)

cat(
    "Duplicated SAF peak IDs:",
    sum(duplicated(x$SAF_peak_id)),
    "\n"
)

cat(
    "Duplicated BED peak IDs:",
    sum(duplicated(x$BED_peak_id)),
    "\n"
)

same_tissue <- (
    x$ATAC_Max_tissue ==
    x$RNA_Max_tissue
)

cat(
    "ATAC/RNA same tissue:",
    sum(same_tissue, na.rm = TRUE),
    "/",
    nrow(x),
    "\n"
)

if (!all(same_tissue)) {
    stop(
        "Strong input unexpectedly contains ",
        "ATAC/RNA tissue mismatches."
    )
}


# ============================================================
# STANDARDIZE UROPA FIELDS
# ============================================================

x$UROPA_relative_location[
    is.na(x$UROPA_relative_location) |
    x$UROPA_relative_location == ""
] <- "NA"

x$UROPA_distance_numeric <- suppressWarnings(
    as.numeric(
        x$UROPA_distance
    )
)

x$UROPA_abs_distance <- abs(
    x$UROPA_distance_numeric
)


# ============================================================
# STRUCTURAL CLASS
#
# Current UROPA annotation is gene-level.
#
# Gene-overlapping:
#   PeakInsideFeature
#   FeatureInsidePeak
#   OverlapStart
#   OverlapEnd
#
# Flanking:
#   Upstream / Downstream
#
# Distances below refer to UROPA feature-anchor distance,
# NOT necessarily distance to TSS.
# ============================================================

overlap_classes <- c(
    "PeakInsideFeature",
    "FeatureInsidePeak",
    "OverlapStart",
    "OverlapEnd"
)

flanking_classes <- c(
    "Upstream",
    "Downstream"
)

x$Structural_class <- "Other_or_NA"

x$Structural_class[
    x$UROPA_relative_location %in%
        overlap_classes
] <- "Gene_overlapping"


is_flanking <- (
    x$UROPA_relative_location %in%
        flanking_classes &
    is.finite(x$UROPA_abs_distance)
)


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance <= 2000
] <- "Flanking_0_2kb"


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance > 2000 &
    x$UROPA_abs_distance <= 20000
] <- "Flanking_2_20kb"


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance > 20000 &
    x$UROPA_abs_distance <= 100000
] <- "Flanking_20_100kb"


x$Structural_class[
    is_flanking &
    x$UROPA_abs_distance > 100000
] <- "Flanking_gt100kb"


class_levels <- c(
    "Gene_overlapping",
    "Flanking_0_2kb",
    "Flanking_2_20kb",
    "Flanking_20_100kb",
    "Flanking_gt100kb",
    "Other_or_NA"
)

x$Structural_class <- factor(
    x$Structural_class,
    levels = class_levels
)


# ============================================================
# WRITE ANNOTATED 702-PAIR MASTER
# ============================================================

write.table(
    x,
    file.path(
        outdir,
        "06_702_strong_pairs_with_structure.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 1. OVERALL SUMMARY
# ============================================================

overall_summary <- data.frame(
    Metric = c(
        "Strong_peak_gene_pairs",
        "Unique_SAF_peaks",
        "Unique_BED_peaks",
        "Unique_genes",
        "Mean_peaks_per_gene",
        "Median_peaks_per_gene",
        "Maximum_peaks_per_gene"
    ),

    Value = c(
        nrow(x),
        length(unique(x$SAF_peak_id)),
        length(unique(x$BED_peak_id)),
        length(unique(x$gene_id)),
        NA,
        NA,
        NA
    )
)


gene_peak_counts <- as.data.frame(
    table(
        x$gene_id
    ),
    stringsAsFactors = FALSE
)

colnames(gene_peak_counts) <- c(
    "gene_id",
    "N_strong_peaks"
)

gene_peak_counts$N_strong_peaks <-
    as.integer(
        gene_peak_counts$N_strong_peaks
    )


overall_summary$Value[
    overall_summary$Metric ==
        "Mean_peaks_per_gene"
] <- mean(
    gene_peak_counts$N_strong_peaks
)


overall_summary$Value[
    overall_summary$Metric ==
        "Median_peaks_per_gene"
] <- median(
    gene_peak_counts$N_strong_peaks
)


overall_summary$Value[
    overall_summary$Metric ==
        "Maximum_peaks_per_gene"
] <- max(
    gene_peak_counts$N_strong_peaks
)


write.table(
    overall_summary,
    file.path(
        outdir,
        "06_overall_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 2. GENE-LEVEL PEAK MULTIPLICITY
# ============================================================

# Add preferred tissue for each gene.
# Strong associations are already same-tissue,
# but theoretically the same gene should normally belong to one
# tissue-specific category.

gene_tissue <- tapply(
    x$ATAC_Max_tissue,
    x$gene_id,
    function(z) {
        paste(
            sort(unique(z)),
            collapse = ";"
        )
    }
)

gene_peak_counts$Tissue <-
    gene_tissue[
        gene_peak_counts$gene_id
    ]


# Mean / maximum coupling per gene

gene_median_rho <- tapply(
    x$rho_mean,
    x$gene_id,
    median,
    na.rm = TRUE
)

gene_max_rho <- tapply(
    x$rho_mean,
    x$gene_id,
    max,
    na.rm = TRUE
)

gene_peak_counts$Median_rho_mean <-
    gene_median_rho[
        gene_peak_counts$gene_id
    ]

gene_peak_counts$Max_rho_mean <-
    gene_max_rho[
        gene_peak_counts$gene_id
    ]


# Count structural classes per gene

for (cc in class_levels) {

    temp <- tapply(
        x$Structural_class == cc,
        x$gene_id,
        sum,
        na.rm = TRUE
    )

    gene_peak_counts[[paste0("N_", cc)]] <- temp[
        gene_peak_counts$gene_id
    ]
}


gene_peak_counts <- gene_peak_counts[
    order(
        -gene_peak_counts$N_strong_peaks,
        -gene_peak_counts$Max_rho_mean,
        gene_peak_counts$gene_id
    ),
]


write.table(
    gene_peak_counts,
    file.path(
        outdir,
        "06_gene_strong_peak_multiplicity.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# Peak-number bins per gene

gene_peak_counts$Peak_number_class <- cut(
    gene_peak_counts$N_strong_peaks,
    breaks = c(
        0,
        1,
        2,
        3,
        5,
        10,
        Inf
    ),
    labels = c(
        "1",
        "2",
        "3",
        "4-5",
        "6-10",
        ">10"
    ),
    right = TRUE
)


multiplicity_summary <- as.data.frame(
    table(
        gene_peak_counts$Peak_number_class
    )
)

colnames(multiplicity_summary) <- c(
    "Strong_peaks_per_gene",
    "Number_of_genes"
)

multiplicity_summary$Percent <- round(
    100 *
    multiplicity_summary$Number_of_genes /
    sum(
        multiplicity_summary$Number_of_genes
    ),
    3
)


write.table(
    multiplicity_summary,
    file.path(
        outdir,
        "06_gene_peak_multiplicity_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 3. UROPA RELATIVE LOCATION
# ============================================================

relative_summary <- as.data.frame(
    table(
        x$UROPA_relative_location
    )
)

colnames(relative_summary) <- c(
    "UROPA_relative_location",
    "Number_of_pairs"
)

relative_summary$Percent <- round(
    100 *
    relative_summary$Number_of_pairs /
    nrow(x),
    3
)

relative_summary <- relative_summary[
    order(
        -relative_summary$Number_of_pairs
    ),
]


write.table(
    relative_summary,
    file.path(
        outdir,
        "06_UROPA_relative_location_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 4. STRUCTURAL CLASS SUMMARY
# ============================================================

structure_summary <- as.data.frame(
    table(
        x$Structural_class
    )
)

colnames(structure_summary) <- c(
    "Structural_class",
    "Number_of_pairs"
)

structure_summary$Percent <- round(
    100 *
    structure_summary$Number_of_pairs /
    nrow(x),
    3
)


write.table(
    structure_summary,
    file.path(
        outdir,
        "06_structural_class_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 5. UPSTREAM / DOWNSTREAM DISTANCE SUMMARY
# ============================================================

flanking <- x[
    x$UROPA_relative_location %in%
        flanking_classes &
    is.finite(
        x$UROPA_abs_distance
    ),
    ,
    drop = FALSE
]


distance_summary <- do.call(
    rbind,
    lapply(
        c(
            "All_flanking",
            "Upstream",
            "Downstream"
        ),
        function(group) {

            if (group == "All_flanking") {

                z <- flanking$UROPA_abs_distance

            } else {

                z <- flanking$UROPA_abs_distance[
                    flanking$UROPA_relative_location ==
                        group
                ]
            }

            if (length(z) == 0) {

                return(
                    data.frame(
                        Group = group,
                        N = 0,
                        Median_distance = NA,
                        Q25_distance = NA,
                        Q75_distance = NA,
                        Min_distance = NA,
                        Max_distance = NA
                    )
                )
            }

            data.frame(
                Group = group,
                N = length(z),
                Median_distance =
                    median(z),

                Q25_distance =
                    as.numeric(
                        quantile(
                            z,
                            0.25,
                            names = FALSE
                        )
                    ),

                Q75_distance =
                    as.numeric(
                        quantile(
                            z,
                            0.75,
                            names = FALSE
                        )
                    ),

                Min_distance =
                    min(z),

                Max_distance =
                    max(z)
            )
        }
    )
)


write.table(
    distance_summary,
    file.path(
        outdir,
        "06_flanking_distance_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 6. TISSUE x STRUCTURE
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


tissue_structure <- as.data.frame(
    table(
        factor(
            x$ATAC_Max_tissue,
            levels = tissues
        ),
        x$Structural_class
    )
)

colnames(tissue_structure) <- c(
    "Tissue",
    "Structural_class",
    "Number_of_pairs"
)


tissue_total <- tapply(
    tissue_structure$Number_of_pairs,
    tissue_structure$Tissue,
    sum
)


tissue_structure$Percent_within_tissue <-
    round(
        100 *
        tissue_structure$Number_of_pairs /
        tissue_total[
            as.character(
                tissue_structure$Tissue
            )
        ],
        3
    )


write.table(
    tissue_structure,
    file.path(
        outdir,
        "06_tissue_by_structural_class.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 7. COUPLING BY STRUCTURAL CLASS
# ============================================================

safe_summary <- function(z) {

    z <- z[
        is.finite(z)
    ]

    if (length(z) == 0) {

        return(
            c(
                N = 0,
                Median = NA,
                Q25 = NA,
                Q75 = NA,
                Mean = NA
            )
        )
    }

    c(
        N = length(z),

        Median =
            median(z),

        Q25 =
            as.numeric(
                quantile(
                    z,
                    0.25,
                    names = FALSE
                )
            ),

        Q75 =
            as.numeric(
                quantile(
                    z,
                    0.75,
                    names = FALSE
                )
            ),

        Mean =
            mean(z)
    )
}


coupling_structure <- do.call(
    rbind,
    lapply(
        class_levels,
        function(cc) {

            z <- x[
                x$Structural_class == cc,
                ,
                drop = FALSE
            ]

            s348 <- safe_summary(
                z$rho_P348
            )

            s350 <- safe_summary(
                z$rho_P350
            )

            sm <- safe_summary(
                z$rho_mean
            )

            data.frame(
                Structural_class = cc,
                N_pairs = nrow(z),

                Median_rho_P348 =
                    s348["Median"],

                Median_rho_P350 =
                    s350["Median"],

                Median_rho_mean =
                    sm["Median"],

                Q25_rho_mean =
                    sm["Q25"],

                Q75_rho_mean =
                    sm["Q75"],

                Mean_rho_mean =
                    sm["Mean"]
            )
        }
    )
)


write.table(
    coupling_structure,
    file.path(
        outdir,
        "06_coupling_by_structural_class.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 8. COUPLING BY TISSUE
# ============================================================

coupling_tissue <- do.call(
    rbind,
    lapply(
        tissues,
        function(tissue) {

            z <- x[
                x$ATAC_Max_tissue ==
                    tissue,
                ,
                drop = FALSE
            ]

            sm <- safe_summary(
                z$rho_mean
            )

            data.frame(
                Tissue = tissue,
                N_pairs = nrow(z),
                Unique_genes =
                    length(
                        unique(
                            z$gene_id
                        )
                    ),

                Median_rho_mean =
                    sm["Median"],

                Q25_rho_mean =
                    sm["Q25"],

                Q75_rho_mean =
                    sm["Q75"]
            )
        }
    )
)


write.table(
    coupling_tissue,
    file.path(
        outdir,
        "06_coupling_by_tissue.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 9. TOP GENES BY NUMBER OF STRONG PEAKS
# ============================================================

top_genes <- head(
    gene_peak_counts,
    50
)


write.table(
    top_genes,
    file.path(
        outdir,
        "06_top50_genes_by_strong_peak_number.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 10. TOP PAIRS BY rho_mean
# ============================================================

top_pairs <- x[
    order(
        -x$rho_mean,
        -x$rho_P348,
        -x$rho_P350
    ),
]


write.table(
    head(
        top_pairs,
        100
    ),
    file.path(
        outdir,
        "06_top100_strong_pairs_by_rho_mean.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)


# ============================================================
# 11. PLOTS
# ============================================================

# ------------------------------------------------------------
# Relative-location bar plot
# ------------------------------------------------------------

p1 <- ggplot(
    relative_summary,
    aes(
        x = reorder(
            UROPA_relative_location,
            Number_of_pairs
        ),
        y = Number_of_pairs
    )
) +
    geom_col() +
    coord_flip() +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title =
            "UROPA relative locations of 702 strong peak-gene associations",
        x = NULL,
        y = "Number of peak-gene associations"
    )


ggsave(
    file.path(
        outdir,
        "06_UROPA_relative_location.pdf"
    ),
    p1,
    width = 7,
    height = 5
)


# ------------------------------------------------------------
# Structural-class bar plot
# ------------------------------------------------------------

p2 <- ggplot(
    structure_summary,
    aes(
        x = Structural_class,
        y = Number_of_pairs
    )
) +
    geom_col() +
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
            "Structural distribution of strongly coupled peak-gene associations",
        x = NULL,
        y = "Number of peak-gene associations"
    )


ggsave(
    file.path(
        outdir,
        "06_structural_class_counts.pdf"
    ),
    p2,
    width = 8,
    height = 5.5
)


# ------------------------------------------------------------
# Gene peak multiplicity
# ------------------------------------------------------------

p3 <- ggplot(
    multiplicity_summary,
    aes(
        x = Strong_peaks_per_gene,
        y = Number_of_genes
    )
) +
    geom_col() +
    theme_bw(
        base_size = 11
    ) +
    labs(
        title =
            "Number of strongly coupled ATAC peaks per gene",
        x = "Strong ATAC peaks per gene",
        y = "Number of genes"
    )


ggsave(
    file.path(
        outdir,
        "06_gene_peak_multiplicity.pdf"
    ),
    p3,
    width = 6.5,
    height = 5
)


# ------------------------------------------------------------
# Tissue x structural class
# ------------------------------------------------------------

p4 <- ggplot(
    tissue_structure,
    aes(
        x = Tissue,
        y = Percent_within_tissue,
        fill = Structural_class
    )
) +
    geom_col() +
    theme_bw(
        base_size = 11
    ) +
    theme(
        axis.text.x =
            element_text(
                angle = 45,
                hjust = 1
            )
    ) +
    labs(
        title =
            "Structural composition of strong peak-gene associations by tissue",
        x = NULL,
        y = "Percentage within tissue",
        fill = "Structural class"
    )


ggsave(
    file.path(
        outdir,
        "06_tissue_structural_composition.pdf"
    ),
    p4,
    width = 9,
    height = 6
)


# ------------------------------------------------------------
# rho_mean by structural class
# ------------------------------------------------------------

plot_data <- x[
    is.finite(
        x$rho_mean
    ),
    ,
    drop = FALSE
]


p5 <- ggplot(
    plot_data,
    aes(
        x = Structural_class,
        y = rho_mean
    )
) +
    geom_boxplot(
        outlier.shape = NA
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
            "RNA-ATAC coupling by genomic structural class",
        x = NULL,
        y = expression(
            rho[mean]
        )
    )


ggsave(
    file.path(
        outdir,
        "06_rho_mean_by_structural_class.pdf"
    ),
    p5,
    width = 8,
    height = 5.5
)


# ------------------------------------------------------------
# Flanking distance distribution
# log10 only for visualization
# ------------------------------------------------------------

if (nrow(flanking) > 0) {

    flanking$log10_distance_plus1 <-
        log10(
            flanking$UROPA_abs_distance + 1
        )

    p6 <- ggplot(
        flanking,
        aes(
            x = log10_distance_plus1
        )
    ) +
        geom_histogram(
            bins = 50
        ) +
        facet_wrap(
            ~ UROPA_relative_location,
            ncol = 1
        ) +
        theme_bw(
            base_size = 11
        ) +
        labs(
            title =
                "Flanking distance distribution of strong peak-gene associations",
            subtitle =
                "Distance refers to UROPA feature-anchor distance, not necessarily TSS distance",
            x = "log10(|UROPA distance| + 1)",
            y = "Number of peak-gene associations"
        )

    ggsave(
        file.path(
            outdir,
            "06_flanking_distance_distribution.pdf"
        ),
        p6,
        width = 7,
        height = 7
    )
}


# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("STEP 06 COMPLETED\n")
cat("============================================================\n\n")

cat("OVERALL SUMMARY:\n")
print(
    overall_summary,
    row.names = FALSE
)

cat("\nGENE PEAK MULTIPLICITY:\n")
print(
    multiplicity_summary,
    row.names = FALSE
)

cat("\nUROPA RELATIVE LOCATION:\n")
print(
    relative_summary,
    row.names = FALSE
)

cat("\nSTRUCTURAL CLASS:\n")
print(
    structure_summary,
    row.names = FALSE
)

cat("\nFLANKING DISTANCE SUMMARY:\n")
print(
    distance_summary,
    row.names = FALSE
)

cat("\nCOUPLING BY STRUCTURAL CLASS:\n")
print(
    coupling_structure,
    row.names = FALSE
)

cat("\nCOUPLING BY TISSUE:\n")
print(
    coupling_tissue,
    row.names = FALSE
)

cat("\nTOP 20 GENES BY NUMBER OF STRONG PEAKS:\n")
print(
    head(
        gene_peak_counts,
        20
    ),
    row.names = FALSE
)

cat("\nIMPORTANT:\n")
cat(
    "UROPA distance in the current gene-level annotation is\n",
    "not assumed to be distance-to-TSS. Therefore structural\n",
    "classes use the neutral terms gene-overlapping and flanking.\n",
    sep = ""
)

cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("============================================================\n")

sessionInfo()

