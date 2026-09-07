#!/usr/bin/env Rscript

# ============================================================
# STEP 10B v2
# Publication-ready Figure 1 + Figure 2
#
# This script ONLY redraws finalized Step01-Step03 results.
# No PCA, LRT, Tau, SPM, or HC definition is recalculated.
#
# Major v2 changes:
#   1. Figure1 PCA panels share legends.
#   2. Figure1 replicate correlations use dots only.
#   3. Figure1 study-design panel is more compact.
#   4. Figure2 B/C display tissues top-to-bottom:
#      Adipose -> Cerebellum -> Cortex -> Hypothalamus
#      -> Liver -> Lung -> Muscle -> Spleen.
#   5. Figure2 top-row widths improved.
#   6. Tau=0.80 labels made horizontal.
#   7. Tau-SPM merge/non-finite QC explicitly reports
#      why ATAC has 203225 plotted features.
#
# Output:
#   PDF vector
#   TIFF 600 dpi LZW
#
# ============================================================


# ============================================================
# 0. Package setup
# ============================================================

options(timeout = 600)

pkgs <- c(
    "ggplot2",
    "dplyr",
    "readr",
    "tidyr",
    "scales",
    "patchwork"
)

for (pkg in pkgs) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

        message("Installing missing package: ", pkg)

        install.packages(
            pkg,
            repos = "https://cloud.r-project.org"
        )
    }
}


suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(tidyr)
    library(scales)
    library(patchwork)
})


# ============================================================
# 1. Input files
# ============================================================

META_FILE <- "sample_pair_metadata.tsv"


RNA_PCA_FILE <-
    "01_sample_QC/RNA/RNA_PCA_coordinates.tsv"

RNA_PCA_VAR_FILE <-
    "01_sample_QC/RNA/RNA_PCA_variance.tsv"

RNA_COR_FILE <-
    "01_sample_QC/RNA/RNA_P348_P350_replicate_correlations.tsv"


ATAC_PCA_FILE <-
    "01_sample_QC/ATAC/ATAC_PCA_coordinates.tsv"

ATAC_PCA_VAR_FILE <-
    "01_sample_QC/ATAC/ATAC_PCA_variance.tsv"

ATAC_COR_FILE <-
    "01_sample_QC/ATAC/ATAC_P348_P350_replicate_correlations.tsv"


RNA_LRT_ALL <-
    "02_multi_tissue_effect/RNA/RNA_multi_tissue_LRT_all.tsv"

RNA_LRT_FDR001 <-
    "02_multi_tissue_effect/RNA/RNA_multi_tissue_LRT_FDR001.tsv"

ATAC_LRT_ALL <-
    "02_multi_tissue_effect/ATAC/ATAC_multi_tissue_LRT_all.tsv"

ATAC_LRT_FDR001 <-
    "02_multi_tissue_effect/ATAC/ATAC_multi_tissue_LRT_FDR001.tsv"


RNA_HC_FILE <-
    "03_tissue_specificity/RNA/RNA_high_confidence_tissue_specific.tsv"

ATAC_HC_FILE <-
    "03_tissue_specificity/ATAC/ATAC_high_confidence_tissue_specific.tsv"


RNA_TAU_ALL <-
    "03_tissue_specificity/RNA/RNA_tissue_specificity_all.tsv"

ATAC_TAU_ALL <-
    "03_tissue_specificity/ATAC/ATAC_tissue_specificity_all.tsv"


RNA_SPM_ALL <-
    "03B_SPM_analysis/RNA/RNA_SPM_all.tsv"

ATAC_SPM_ALL <-
    "03B_SPM_analysis/ATAC/ATAC_SPM_all.tsv"


# ============================================================
# 2. Output directory
# ============================================================

OUTDIR <-
    "10_publication_figures/Step10B_Figure1_2_v2"

FIG1DIR <- file.path(
    OUTDIR,
    "Figure1"
)

FIG2DIR <- file.path(
    OUTDIR,
    "Figure2"
)

dir.create(
    FIG1DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    FIG2DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. Fixed publication settings
# ============================================================

TISSUE_LEVELS <- c(
    "Adipose",
    "Cerebellum",
    "Cortex",
    "Hypothalamus",
    "Liver",
    "Lung",
    "Muscle",
    "Spleen"
)


TISSUE_COLORS <- c(
    "Adipose"      = "#D89000",
    "Cerebellum"   = "#7CAE00",
    "Cortex"       = "#00BA38",
    "Hypothalamus" = "#00C08B",
    "Liver"        = "#00BFC4",
    "Lung"         = "#619CFF",
    "Muscle"       = "#C77CFF",
    "Spleen"       = "#F564E3"
)


OMIC_COLORS <- c(
    "RNA"  = "#377EB8",
    "ATAC" = "#E41A1C"
)


ANIMAL_SHAPES <- c(
    "P348" = 16,
    "P350" = 17
)


theme_publication <- function(
    base_size = 8.5
) {

    theme_classic(
        base_size = base_size,
        base_family = "sans"
    ) +

        theme(

            plot.title =
                element_text(
                    size = base_size + 0.5,
                    face = "bold",
                    hjust = 0
                ),

            plot.subtitle =
                element_text(
                    size = base_size - 0.5
                ),

            axis.title =
                element_text(
                    size = base_size
                ),

            axis.text =
                element_text(
                    size = base_size - 1
                ),

            legend.title =
                element_text(
                    size = base_size - 0.5,
                    face = "bold"
                ),

            legend.text =
                element_text(
                    size = base_size - 1
                ),

            strip.text =
                element_text(
                    size = base_size,
                    face = "bold"
                ),

            plot.margin =
                margin(
                    5,
                    5,
                    5,
                    5
                )
        )
}


# ============================================================
# 4. Generic helpers
# ============================================================

required_files <- c(
    META_FILE,

    RNA_PCA_FILE,
    RNA_PCA_VAR_FILE,
    RNA_COR_FILE,

    ATAC_PCA_FILE,
    ATAC_PCA_VAR_FILE,
    ATAC_COR_FILE,

    RNA_LRT_ALL,
    RNA_LRT_FDR001,

    ATAC_LRT_ALL,
    ATAC_LRT_FDR001,

    RNA_HC_FILE,
    ATAC_HC_FILE,

    RNA_TAU_ALL,
    ATAC_TAU_ALL,

    RNA_SPM_ALL,
    ATAC_SPM_ALL
)


for (path in required_files) {

    if (!file.exists(path)) {

        stop(
            "Required file not found: ",
            path
        )
    }
}


find_col <- function(
    df,
    exact = character(),
    regex = character(),
    required = TRUE,
    description = ""
) {

    nms <- names(df)

    for (x in exact) {

        hit <- which(
            tolower(nms) ==
                tolower(x)
        )

        if (length(hit) > 0) {

            return(
                nms[hit[1]]
            )
        }
    }


    for (pattern in regex) {

        hit <- grep(
            pattern,
            nms,
            ignore.case = TRUE
        )

        if (length(hit) > 0) {

            return(
                nms[hit[1]]
            )
        }
    }


    if (required) {

        stop(
            "Cannot detect column: ",
            description,
            "\nAvailable columns:\n",
            paste(
                nms,
                collapse = ", "
            )
        )
    }


    NULL
}


save_panel <- function(
    plot,
    prefix,
    width,
    height
) {

    ggsave(
        paste0(
            prefix,
            ".pdf"
        ),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = cairo_pdf
    )


    ggsave(
        paste0(
            prefix,
            "_600dpi.tiff"
        ),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = 600,
        compression = "lzw"
    )
}


# ============================================================
# 5. Metadata
# ============================================================

meta <- read_tsv(
    META_FILE,
    show_col_types = FALSE
)


required_meta <- c(
    "Tissue",
    "Animal",
    "RNA_SRR",
    "ATAC_SRR"
)


missing_meta <- setdiff(
    required_meta,
    names(meta)
)


if (length(missing_meta) > 0) {

    stop(
        "Metadata missing columns: ",
        paste(
            missing_meta,
            collapse = ", "
        )
    )
}


meta <- meta %>%

    mutate(

        Tissue =
            factor(
                Tissue,
                levels = TISSUE_LEVELS
            ),

        Animal =
            factor(
                Animal,
                levels = c(
                    "P348",
                    "P350"
                )
            )
    )


if (nrow(meta) != 16) {

    stop(
        "Expected 16 paired metadata rows; observed ",
        nrow(meta)
    )
}


# ============================================================
# 6. FIGURE 1A — study design
# ============================================================

design_long <- meta %>%

    select(
        Tissue,
        Animal,
        RNA_SRR,
        ATAC_SRR
    ) %>%

    pivot_longer(

        cols = c(
            RNA_SRR,
            ATAC_SRR
        ),

        names_to = "Assay",
        values_to = "SRR"
    ) %>%

    mutate(

        Assay =
            recode(
                Assay,
                RNA_SRR = "RNA",
                ATAC_SRR = "ATAC"
            ),

        Assay =
            factor(
                Assay,
                levels = c(
                    "RNA",
                    "ATAC"
                )
            ),

        Animal_y =
            as.numeric(
                Animal
            ),

        Assay_offset =
            ifelse(
                Assay == "RNA",
                -0.11,
                0.11
            ),

        y =
            Animal_y +
            Assay_offset
    )


p1A <- ggplot(
    design_long,
    aes(
        x = Tissue,
        y = y
    )
) +

    geom_segment(

        data = meta,

        aes(
            x = Tissue,
            xend = Tissue,
            y = as.numeric(Animal) - 0.12,
            yend = as.numeric(Animal) + 0.12
        ),

        inherit.aes = FALSE,

        linewidth = 0.35
    ) +

    geom_point(

        aes(
            shape = Assay,
            fill = Assay
        ),

        size = 2.9,
        stroke = 0.4
    ) +

    scale_shape_manual(
        values = c(
            "RNA" = 21,
            "ATAC" = 24
        )
    ) +

    scale_fill_manual(
        values = OMIC_COLORS
    ) +

    scale_y_continuous(

        breaks = c(
            1,
            2
        ),

        labels = c(
            "P348",
            "P350"
        ),

        limits = c(
            0.68,
            2.32
        )
    ) +

    scale_x_discrete(
        drop = FALSE
    ) +

    labs(
        x = NULL,
        y = "Animal",
        shape = "Assay",
        fill = "Assay"
    ) +

    theme_publication() +

    theme(

        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            ),

        legend.position =
            "top"
    )


save_panel(
    p1A,
    file.path(
        FIG1DIR,
        "Figure1A_study_design_v2"
    ),
    width = 7.0,
    height = 2.0
)


# ============================================================
# 7. PCA helpers
# ============================================================

extract_pca_variance <- function(path) {

    x <- read_tsv(
        path,
        show_col_types = FALSE
    )


    candidate <- find_col(

        x,

        exact = c(
            "variance_percent",
            "variance_explained",
            "percent_variance",
            "Variance"
        ),

        regex = c(
            "variance",
            "percent",
            "explained"
        ),

        required = FALSE
    )


    if (!is.null(candidate)) {

        values <- suppressWarnings(
            as.numeric(
                x[[candidate]]
            )
        )

    } else {

        numeric_cols <- names(x)[
            vapply(
                x,
                is.numeric,
                logical(1)
            )
        ]


        if (length(numeric_cols) == 0) {

            return(
                c(
                    NA_real_,
                    NA_real_
                )
            )
        }


        values <- as.numeric(
            x[[numeric_cols[1]]]
        )
    }


    values <- values[
        is.finite(
            values
        )
    ]


    if (length(values) < 2) {

        return(
            c(
                NA_real_,
                NA_real_
            )
        )
    }


    if (
        max(
            values,
            na.rm = TRUE
        ) <= 1
    ) {

        values <- values * 100
    }


    values[1:2]
}


prepare_pca <- function(
    pca_file,
    variance_file,
    modality
) {

    pca <- read_tsv(
        pca_file,
        show_col_types = FALSE
    )


    sample_col <- find_col(

        pca,

        exact = c(
            "sampleID",
            "Sample",
            "sample",
            "SampleID"
        ),

        regex = c(
            "^sample",
            "sample.*id"
        ),

        description =
            paste0(
                modality,
                " PCA sample ID"
            )
    )


    pc1_col <- find_col(
        pca,
        exact = "PC1",
        regex = "^PC1$",
        description =
            paste0(
                modality,
                " PC1"
            )
    )


    pc2_col <- find_col(
        pca,
        exact = "PC2",
        regex = "^PC2$",
        description =
            paste0(
                modality,
                " PC2"
            )
    )


    pca2 <- pca %>%

        transmute(

            Sample =
                as.character(
                    .data[[sample_col]]
                ),

            PC1 =
                as.numeric(
                    .data[[pc1_col]]
                ),

            PC2 =
                as.numeric(
                    .data[[pc2_col]]
                )
        )


    modality_meta <- if (
        modality == "RNA"
    ) {

        meta %>%

            transmute(
                Sample =
                    as.character(
                        RNA_SRR
                    ),
                Tissue,
                Animal
            )

    } else {

        meta %>%

            transmute(
                Sample =
                    as.character(
                        ATAC_SRR
                    ),
                Tissue,
                Animal
            )
    }


    pca2 <- pca2 %>%

        left_join(
            modality_meta,
            by = "Sample"
        )


    if (
        any(
            is.na(
                pca2$Tissue
            )
        )
    ) {

        stop(
            modality,
            " PCA samples cannot all be mapped to metadata."
        )
    }


    variance <- extract_pca_variance(
        variance_file
    )


    list(
        data = pca2,
        PC1 = variance[1],
        PC2 = variance[2]
    )
}


rna_pca <- prepare_pca(
    RNA_PCA_FILE,
    RNA_PCA_VAR_FILE,
    "RNA"
)


atac_pca <- prepare_pca(
    ATAC_PCA_FILE,
    ATAC_PCA_VAR_FILE,
    "ATAC"
)


# ============================================================
# 8. FIGURE 1B/C — PCA
# ============================================================

make_pca_plot <- function(obj) {

    xlab <- if (
        is.finite(
            obj$PC1
        )
    ) {

        sprintf(
            "PC1 (%.2f%%)",
            obj$PC1
        )

    } else {

        "PC1"
    }


    ylab <- if (
        is.finite(
            obj$PC2
        )
    ) {

        sprintf(
            "PC2 (%.2f%%)",
            obj$PC2
        )

    } else {

        "PC2"
    }


    ggplot(
        obj$data,
        aes(
            x = PC1,
            y = PC2,
            color = Tissue,
            shape = Animal
        )
    ) +

        geom_point(
            size = 2.8,
            alpha = 0.95
        ) +

        scale_color_manual(
            values = TISSUE_COLORS,
            drop = FALSE
        ) +

        scale_shape_manual(
            values = ANIMAL_SHAPES
        ) +

        labs(
            x = xlab,
            y = ylab,
            color = "Tissue",
            shape = "Animal"
        ) +

        theme_publication()
}


p1B <- make_pca_plot(
    rna_pca
)


p1C <- make_pca_plot(
    atac_pca
)


save_panel(
    p1B,
    file.path(
        FIG1DIR,
        "Figure1B_RNA_PCA_v2"
    ),
    width = 3.5,
    height = 3.1
)


save_panel(
    p1C,
    file.path(
        FIG1DIR,
        "Figure1C_ATAC_PCA_v2"
    ),
    width = 3.5,
    height = 3.1
)


# ============================================================
# 9. Replicate correlation
# ============================================================

prepare_replicate_cor <- function(
    path,
    modality
) {

    x <- read_tsv(
        path,
        show_col_types = FALSE
    )


    tissue_col <- find_col(

        x,

        exact = c(
            "Tissue",
            "tissue"
        ),

        regex = "tissue",

        description =
            paste0(
                modality,
                " tissue"
            )
    )


    cor_col <- find_col(

        x,

        exact = c(
            "Pearson",
            "Pearson_r",
            "cor",
            "Correlation",
            "correlation"
        ),

        regex = c(
            "pearson",
            "^cor$",
            "correlation"
        ),

        description =
            paste0(
                modality,
                " replicate correlation"
            )
    )


    x %>%

        transmute(

            Tissue =
                factor(
                    .data[[tissue_col]],
                    levels = TISSUE_LEVELS
                ),

            Correlation =
                as.numeric(
                    .data[[cor_col]]
                )
        )
}


rna_cor <- prepare_replicate_cor(
    RNA_COR_FILE,
    "RNA"
)


atac_cor <- prepare_replicate_cor(
    ATAC_COR_FILE,
    "ATAC"
)


if (
    nrow(rna_cor) != 8 ||
    nrow(atac_cor) != 8
) {

    stop(
        "Expected 8 tissue-level replicate correlations."
    )
}


make_cor_plot <- function(dat) {

    mean_r <- mean(
        dat$Correlation,
        na.rm = TRUE
    )


    ggplot(
        dat,
        aes(
            x = Tissue,
            y = Correlation,
            color = Tissue
        )
    ) +

        geom_point(
            size = 2.9
        ) +

        geom_hline(
            yintercept = mean_r,
            linetype = "dashed",
            linewidth = 0.45
        ) +

        scale_color_manual(
            values = TISSUE_COLORS,
            guide = "none"
        ) +

        scale_y_continuous(

            limits = c(
                0.80,
                1.005
            ),

            breaks = seq(
                0.80,
                1.00,
                0.05
            )
        ) +

        labs(
            x = NULL,
            y = "P348–P350 Pearson correlation",
            subtitle =
                sprintf(
                    "Mean Pearson r = %.3f",
                    mean_r
                )
        ) +

        theme_publication() +

        theme(

            axis.text.x =
                element_text(
                    angle = 35,
                    hjust = 1
                ),

            plot.subtitle =
                element_text(
                    hjust = 0.5,
                    size = 7.5
                )
        )
}


p1D <- make_cor_plot(
    rna_cor
)


p1E <- make_cor_plot(
    atac_cor
)


save_panel(
    p1D,
    file.path(
        FIG1DIR,
        "Figure1D_RNA_replicate_correlation_v2"
    ),
    width = 3.5,
    height = 2.7
)


save_panel(
    p1E,
    file.path(
        FIG1DIR,
        "Figure1E_ATAC_replicate_correlation_v2"
    ),
    width = 3.5,
    height = 2.7
)


# ============================================================
# 10. Assemble FIGURE 1
# ============================================================

pca_row <- (

    p1B +
    p1C +

    plot_layout(
        guides = "collect"
    )

) &

    theme(
        legend.position = "right"
    )


cor_row <- (
    p1D +
    p1E
)


Figure1 <- (

    p1A /

    pca_row /

    cor_row

) +

    plot_annotation(

        tag_levels = "A",

        theme = theme(

            plot.tag =
                element_text(
                    size = 11,
                    face = "bold",
                    family = "sans"
                )
        )
    ) +

    plot_layout(
        heights = c(
            0.58,
            1,
            0.82
        )
    )


ggsave(
    file.path(
        FIG1DIR,
        "Figure1_FINAL_v2.pdf"
    ),
    Figure1,
    width = 7.1,
    height = 7.0,
    device = cairo_pdf
)


ggsave(
    file.path(
        FIG1DIR,
        "Figure1_FINAL_v2_600dpi.tiff"
    ),
    Figure1,
    width = 7.1,
    height = 7.0,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 11. FIGURE 2A — LRT
# ============================================================

rna_lrt_all <- read_tsv(
    RNA_LRT_ALL,
    show_col_types = FALSE
)

rna_lrt_001 <- read_tsv(
    RNA_LRT_FDR001,
    show_col_types = FALSE
)

atac_lrt_all <- read_tsv(
    ATAC_LRT_ALL,
    show_col_types = FALSE
)

atac_lrt_001 <- read_tsv(
    ATAC_LRT_FDR001,
    show_col_types = FALSE
)


lrt_summary <- tibble(

    Omic =
        factor(
            c(
                "RNA",
                "ATAC"
            ),
            levels = c(
                "RNA",
                "ATAC"
            )
        ),

    Tested =
        c(
            nrow(
                rna_lrt_all
            ),
            nrow(
                atac_lrt_all
            )
        ),

    FDR001 =
        c(
            nrow(
                rna_lrt_001
            ),
            nrow(
                atac_lrt_001
            )
        )
) %>%

    mutate(

        Percent =
            100 *
            FDR001 /
            Tested,

        Label =
            paste0(
                comma(
                    FDR001
                ),
                " / ",
                comma(
                    Tested
                ),
                "\n",
                sprintf(
                    "%.1f%%",
                    Percent
                )
            )
    )


expected_lrt <- tibble(

    Omic =
        factor(
            c(
                "RNA",
                "ATAC"
            ),
            levels = c(
                "RNA",
                "ATAC"
            )
        ),

    Tested_expected =
        c(
            22240L,
            203210L
        ),

    FDR001_expected =
        c(
            18870L,
            176917L
        )
)


lrt_qc <- lrt_summary %>%

    left_join(
        expected_lrt,
        by = "Omic"
    ) %>%

    mutate(

        Status =
            ifelse(
                Tested ==
                    Tested_expected &
                FDR001 ==
                    FDR001_expected,
                "PASS",
                "FAIL"
            )
    )


if (
    any(
        lrt_qc$Status != "PASS"
    )
) {

    stop(
        "LRT frozen-result QC failed."
    )
}


write_tsv(
    lrt_qc,
    file.path(
        FIG2DIR,
        "Figure2A_LRT_source_QC_v2.tsv"
    )
)


p2A <- ggplot(
    lrt_summary,
    aes(
        x = Omic,
        y = Percent,
        fill = Omic
    )
) +

    geom_col(
        width = 0.60
    ) +

    geom_text(

        aes(
            label = Label
        ),

        vjust = -0.35,
        size = 2.55,
        lineheight = 0.95
    ) +

    scale_fill_manual(
        values = OMIC_COLORS,
        guide = "none"
    ) +

    scale_y_continuous(

        limits = c(
            0,
            100
        ),

        expand = expansion(
            mult = c(
                0,
                0.09
            )
        )
    ) +

    labs(
        x = NULL,
        y = "Genes/peaks with tissue effect (%)",
        subtitle = "DESeq2 LRT, BH FDR < 0.01"
    ) +

    theme_publication() +

    theme(
        plot.subtitle =
            element_text(
                hjust = 0.5
            )
    )


save_panel(
    p2A,
    file.path(
        FIG2DIR,
        "Figure2A_global_tissue_LRT_v2"
    ),
    width = 2.8,
    height = 3.0
)


# ============================================================
# 12. FIGURE 2B/C — HC counts
# ============================================================

prepare_hc_counts <- function(
    path,
    modality
) {

    x <- read_tsv(
        path,
        show_col_types = FALSE
    )


    tissue_col <- find_col(

        x,

        exact = c(
            "Max_tissue",
            "max_tissue",
            "Tissue",
            "tissue"
        ),

        regex = c(
            "max.*tissue",
            "^tissue$"
        ),

        description =
            paste0(
                modality,
                " HC tissue"
            )
    )


    x %>%

        count(
            Tissue =
                as.character(
                    .data[[tissue_col]]
                ),
            name = "N"
        ) %>%

        complete(
            Tissue = TISSUE_LEVELS,
            fill = list(
                N = 0L
            )
        ) %>%

        mutate(

            Tissue =
                factor(
                    Tissue,
                    levels = TISSUE_LEVELS
                ),

            Tissue_plot =
                factor(
                    as.character(
                        Tissue
                    ),
                    levels =
                        rev(
                            TISSUE_LEVELS
                        )
                ),

            Omic =
                modality
        )
}


rna_hc_counts <- prepare_hc_counts(
    RNA_HC_FILE,
    "RNA"
)


atac_hc_counts <- prepare_hc_counts(
    ATAC_HC_FILE,
    "ATAC"
)


expected_rna_hc <- c(
    Adipose = 408,
    Cerebellum = 573,
    Cortex = 397,
    Hypothalamus = 244,
    Liver = 1182,
    Lung = 610,
    Muscle = 733,
    Spleen = 889
)


expected_atac_hc <- c(
    Adipose = 47,
    Cerebellum = 4643,
    Cortex = 2689,
    Hypothalamus = 785,
    Liver = 5071,
    Lung = 725,
    Muscle = 20853,
    Spleen = 5732
)


check_hc <- function(
    dat,
    expected,
    modality
) {

    observed <- setNames(
        dat$N,
        as.character(
            dat$Tissue
        )
    )


    pass <- all(
        observed[
            names(
                expected
            )
        ] ==
            expected
    )


    tibble(

        Omic =
            modality,

        Observed_total =
            sum(
                dat$N
            ),

        Expected_total =
            sum(
                expected
            ),

        Per_tissue_status =
            ifelse(
                pass,
                "PASS",
                "FAIL"
            )
    )
}


hc_qc <- bind_rows(

    check_hc(
        rna_hc_counts,
        expected_rna_hc,
        "RNA"
    ),

    check_hc(
        atac_hc_counts,
        expected_atac_hc,
        "ATAC"
    )
)


if (
    any(
        hc_qc$Per_tissue_status != "PASS"
    )
) {

    stop(
        "HC frozen-result QC failed."
    )
}


write_tsv(
    hc_qc,
    file.path(
        FIG2DIR,
        "Figure2BC_HC_source_QC_v2.tsv"
    )
)


make_hc_count_plot <- function(
    dat,
    ylab
) {

    ggplot(
        dat,
        aes(
            x = Tissue_plot,
            y = N,
            fill = Tissue
        )
    ) +

        geom_col(
            width = 0.72
        ) +

        geom_text(

            aes(
                label =
                    comma(
                        N
                    )
            ),

            hjust = -0.10,
            size = 2.45
        ) +

        coord_flip(
            clip = "off"
        ) +

        scale_fill_manual(
            values = TISSUE_COLORS,
            guide = "none",
            drop = FALSE
        ) +

        scale_y_continuous(

            labels = comma,

            expand =
                expansion(
                    mult = c(
                        0,
                        0.20
                    )
                )
        ) +

        labs(
            x = NULL,
            y = ylab
        ) +

        theme_publication()
}


p2B <- make_hc_count_plot(
    rna_hc_counts,
    "HC tissue-specific RNA genes"
)


p2C <- make_hc_count_plot(
    atac_hc_counts,
    "HC tissue-specific ATAC peaks"
)


save_panel(
    p2B,
    file.path(
        FIG2DIR,
        "Figure2B_RNA_HC_counts_v2"
    ),
    width = 3.5,
    height = 3.0
)


save_panel(
    p2C,
    file.path(
        FIG2DIR,
        "Figure2C_ATAC_HC_counts_v2"
    ),
    width = 3.8,
    height = 3.0
)


# ============================================================
# 13. Tau-SPM source QC and preparation
# ============================================================

prepare_tau_spm <- function(
    tau_file,
    spm_file,
    modality
) {

    tau <- read_tsv(
        tau_file,
        show_col_types = FALSE
    )


    spm <- read_tsv(
        spm_file,
        show_col_types = FALSE
    )


    required_tau <- c(
        "Feature_ID",
        "Tau"
    )


    required_spm <- c(
        "Feature_ID",
        "Max_SPM"
    )


    if (
        !all(
            required_tau %in%
                names(
                    tau
                )
        )
    ) {

        stop(
            modality,
            " Tau file missing Feature_ID/Tau."
        )
    }


    if (
        !all(
            required_spm %in%
                names(
                    spm
                )
        )
    ) {

        stop(
            modality,
            " SPM file missing Feature_ID/Max_SPM."
        )
    }


    tau2 <- tau %>%

        transmute(

            Feature_ID =
                as.character(
                    Feature_ID
                ),

            Tau =
                as.numeric(
                    Tau
                )
        )


    spm2 <- spm %>%

        transmute(

            Feature_ID =
                as.character(
                    Feature_ID
                ),

            Max_SPM =
                as.numeric(
                    Max_SPM
                )
        )


    tau_dup <- sum(
        duplicated(
            tau2$Feature_ID
        )
    )


    spm_dup <- sum(
        duplicated(
            spm2$Feature_ID
        )
    )


    if (
        tau_dup > 0 ||
        spm_dup > 0
    ) {

        stop(
            modality,
            " Tau/SPM Feature_ID contains duplicate IDs."
        )
    }


    tau_only <- setdiff(
        tau2$Feature_ID,
        spm2$Feature_ID
    )


    spm_only <- setdiff(
        spm2$Feature_ID,
        tau2$Feature_ID
    )


    joined <- inner_join(
        tau2,
        spm2,
        by = "Feature_ID"
    )


    nonfinite <- joined %>%

        filter(
            !is.finite(
                Tau
            ) |
            !is.finite(
                Max_SPM
            )
        ) %>%

        mutate(
            Issue =
                case_when(

                    !is.finite(Tau) &
                    !is.finite(Max_SPM) ~
                        "Nonfinite_Tau_and_MaxSPM",

                    !is.finite(Tau) ~
                        "Nonfinite_Tau",

                    !is.finite(Max_SPM) ~
                        "Nonfinite_MaxSPM",

                    TRUE ~
                        "Unknown"
                )
        )


    valid <- joined %>%

        filter(
            is.finite(
                Tau
            ),
            is.finite(
                Max_SPM
            )
        )


    id_issue_rows <- bind_rows(

        tibble(
            Feature_ID = tau_only,
            Tau = NA_real_,
            Max_SPM = NA_real_,
            Issue = "Present_in_Tau_only"
        ),

        tibble(
            Feature_ID = spm_only,
            Tau = NA_real_,
            Max_SPM = NA_real_,
            Issue = "Present_in_SPM_only"
        ),

        nonfinite %>%
            select(
                Feature_ID,
                Tau,
                Max_SPM,
                Issue
            )
    )


    qc <- tibble(

        Omic =
            modality,

        Tau_source_rows =
            nrow(
                tau2
            ),

        SPM_source_rows =
            nrow(
                spm2
            ),

        Tau_unique_IDs =
            n_distinct(
                tau2$Feature_ID
            ),

        SPM_unique_IDs =
            n_distinct(
                spm2$Feature_ID
            ),

        Tau_duplicate_IDs =
            tau_dup,

        SPM_duplicate_IDs =
            spm_dup,

        Tau_only_IDs =
            length(
                tau_only
            ),

        SPM_only_IDs =
            length(
                spm_only
            ),

        Inner_join_rows =
            nrow(
                joined
            ),

        Nonfinite_rows =
            nrow(
                nonfinite
            ),

        Final_plot_rows =
            nrow(
                valid
            )
    )


    list(
        data = valid,
        qc = qc,
        issues = id_issue_rows
    )
}


rna_tau_spm_obj <- prepare_tau_spm(
    RNA_TAU_ALL,
    RNA_SPM_ALL,
    "RNA"
)


atac_tau_spm_obj <- prepare_tau_spm(
    ATAC_TAU_ALL,
    ATAC_SPM_ALL,
    "ATAC"
)


rna_tau_spm <- rna_tau_spm_obj$data
atac_tau_spm <- atac_tau_spm_obj$data


tau_spm_merge_qc <- bind_rows(
    rna_tau_spm_obj$qc,
    atac_tau_spm_obj$qc
)


write_tsv(
    tau_spm_merge_qc,
    file.path(
        FIG2DIR,
        "Figure2DE_Tau_SPM_merge_QC_v2.tsv"
    )
)


write_tsv(
    rna_tau_spm_obj$issues,
    file.path(
        FIG2DIR,
        "Figure2D_RNA_Tau_SPM_excluded_features.tsv"
    )
)


write_tsv(
    atac_tau_spm_obj$issues,
    file.path(
        FIG2DIR,
        "Figure2E_ATAC_Tau_SPM_excluded_features.tsv"
    )
)


# ============================================================
# 14. FIGURE 2D/E — Tau-SPM landscapes
# ============================================================

make_tau_spm_plot <- function(dat) {

    ggplot(
        dat,
        aes(
            x = Tau,
            y = Max_SPM
        )
    ) +

        geom_bin_2d(
            bins = 48
        ) +

        geom_vline(

            xintercept = 0.80,

            linetype = "dashed",

            linewidth = 0.45
        ) +

        annotate(

            "text",

            x = 0.815,

            y = 0.97,

            label = "\u03C4 = 0.80",

            hjust = 0,

            vjust = 1,

            size = 2.5
        ) +

        scale_fill_viridis_c(

            option = "C",

            trans = "log10",

            name =
                "Features\nper bin"
        ) +

        scale_x_continuous(

            breaks = seq(
                0,
                1,
                0.2
            )
        ) +

        scale_y_continuous(

            breaks = seq(
                0,
                1,
                0.2
            )
        ) +

        coord_cartesian(

            xlim = c(
                0,
                1
            ),

            ylim = c(
                0,
                1
            ),

            expand = FALSE
        ) +

        labs(

            x =
                expression(
                    "Tissue specificity (" *
                    tau *
                    ")"
                ),

            y =
                "Maximum SPM"
        ) +

        theme_publication() +

        theme(
            legend.position =
                "right"
        )
}


p2D <- make_tau_spm_plot(
    rna_tau_spm
)


p2E <- make_tau_spm_plot(
    atac_tau_spm
)


save_panel(
    p2D,
    file.path(
        FIG2DIR,
        "Figure2D_RNA_Tau_SPM_landscape_v2"
    ),
    width = 3.5,
    height = 3.1
)


save_panel(
    p2E,
    file.path(
        FIG2DIR,
        "Figure2E_ATAC_Tau_SPM_landscape_v2"
    ),
    width = 3.5,
    height = 3.1
)


# ============================================================
# 15. Assemble FIGURE 2
# ============================================================

top_row <- (

    p2A +
    p2B +
    p2C +

    plot_layout(
        widths = c(
            1.10,
            1.35,
            1.55
        )
    )
)


bottom_row <- (
    p2D +
    p2E
)


Figure2 <- (

    top_row /
    bottom_row

) +

    plot_annotation(

        tag_levels = "A",

        theme = theme(

            plot.tag =
                element_text(
                    size = 11,
                    face = "bold",
                    family = "sans"
                )
        )
    ) +

    plot_layout(
        heights = c(
            1,
            1.08
        )
    )


ggsave(
    file.path(
        FIG2DIR,
        "Figure2_FINAL_v2.pdf"
    ),
    Figure2,
    width = 7.1,
    height = 6.3,
    device = cairo_pdf
)


ggsave(
    file.path(
        FIG2DIR,
        "Figure2_FINAL_v2_600dpi.tiff"
    ),
    Figure2,
    width = 7.1,
    height = 6.3,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 16. Figure 1 source QC
# ============================================================

figure1_qc <- tibble(

    Metric = c(

        "Metadata_rows",
        "Unique_tissues",
        "Unique_animals",

        "RNA_PCA_samples",
        "ATAC_PCA_samples",

        "RNA_replicate_correlation_tissues",
        "ATAC_replicate_correlation_tissues",

        "RNA_mean_P348_P350_correlation",
        "ATAC_mean_P348_P350_correlation",

        "RNA_PC1_variance_percent",
        "RNA_PC2_variance_percent",

        "ATAC_PC1_variance_percent",
        "ATAC_PC2_variance_percent"
    ),

    Value = c(

        nrow(meta),
        n_distinct(meta$Tissue),
        n_distinct(meta$Animal),

        nrow(rna_pca$data),
        nrow(atac_pca$data),

        nrow(rna_cor),
        nrow(atac_cor),

        mean(
            rna_cor$Correlation
        ),

        mean(
            atac_cor$Correlation
        ),

        rna_pca$PC1,
        rna_pca$PC2,

        atac_pca$PC1,
        atac_pca$PC2
    )
)


write_tsv(
    figure1_qc,
    file.path(
        FIG1DIR,
        "Figure1_source_data_summary_v2.tsv"
    )
)


# ============================================================
# 17. Figure 2 source summary
# ============================================================

figure2_summary <- bind_rows(

    lrt_summary %>%

        transmute(

            Section =
                "LRT",

            Omic =
                as.character(
                    Omic
                ),

            Metric =
                "FDR001",

            Value =
                FDR001,

            Denominator =
                Tested,

            Percent =
                Percent
        ),


    rna_hc_counts %>%

        transmute(

            Section =
                "HC_counts",

            Omic =
                "RNA",

            Metric =
                as.character(
                    Tissue
                ),

            Value =
                N,

            Denominator =
                NA_real_,

            Percent =
                NA_real_
        ),


    atac_hc_counts %>%

        transmute(

            Section =
                "HC_counts",

            Omic =
                "ATAC",

            Metric =
                as.character(
                    Tissue
                ),

            Value =
                N,

            Denominator =
                NA_real_,

            Percent =
                NA_real_
        )
)


write_tsv(
    figure2_summary,
    file.path(
        FIG2DIR,
        "Figure2_source_data_summary_v2.tsv"
    )
)


# ============================================================
# 18. Caption notes
# ============================================================

caption_lines <- c(

    "FIGURE 1",

    "A. Paired multi-tissue study design comprising P348 and P350, with matched RNA-seq and ATAC-seq data from eight tissues.",

    sprintf(
        "B. RNA PCA based on the finalized Step01 DESeq2 VST matrix. PC1 = %.2f%% and PC2 = %.2f%%.",
        rna_pca$PC1,
        rna_pca$PC2
    ),

    sprintf(
        "C. ATAC PCA based on the finalized Step01 DESeq2 VST matrix. PC1 = %.2f%% and PC2 = %.2f%%.",
        atac_pca$PC1,
        atac_pca$PC2
    ),

    sprintf(
        "D. Within-tissue P348-P350 RNA reproducibility. Mean Pearson r = %.3f.",
        mean(
            rna_cor$Correlation
        )
    ),

    sprintf(
        "E. Within-tissue P348-P350 ATAC reproducibility. Mean Pearson r = %.3f.",
        mean(
            atac_cor$Correlation
        )
    ),

    "",

    "FIGURE 2",

    sprintf(
        paste0(
            "A. Global tissue effects assessed using DESeq2 likelihood-ratio tests ",
            "controlling for animal identity. At BH FDR < 0.01, ",
            "%s/%s RNA genes (%.1f%%) and %s/%s ATAC peaks (%.1f%%) ",
            "showed significant tissue effects."
        ),
        comma(
            lrt_summary$FDR001[
                lrt_summary$Omic == "RNA"
            ]
        ),
        comma(
            lrt_summary$Tested[
                lrt_summary$Omic == "RNA"
            ]
        ),
        lrt_summary$Percent[
            lrt_summary$Omic == "RNA"
        ],
        comma(
            lrt_summary$FDR001[
                lrt_summary$Omic == "ATAC"
            ]
        ),
        comma(
            lrt_summary$Tested[
                lrt_summary$Omic == "ATAC"
            ]
        ),
        lrt_summary$Percent[
            lrt_summary$Omic == "ATAC"
        ]
    ),

    sprintf(
        paste0(
            "B-C. Finalized high-confidence tissue-specific sets comprised ",
            "%s RNA genes and %s ATAC peaks."
        ),
        comma(
            sum(
                rna_hc_counts$N
            )
        ),
        comma(
            sum(
                atac_hc_counts$N
            )
        )
    ),

    paste0(
        "D-E. Joint Tau-SPM landscapes. Tau describes across-tissue specificity, ",
        "whereas maximum SPM provides a continuous description of tissue preference. ",
        "The dashed line marks Tau = 0.80. SPM is not used here to redefine ",
        "high-confidence membership."
    )
)


writeLines(
    caption_lines,
    file.path(
        OUTDIR,
        "Figure1_2_caption_notes_v2.txt"
    )
)


# ============================================================
# 19. Final file QC
# ============================================================

output_files <- c(

    file.path(
        FIG1DIR,
        "Figure1_FINAL_v2.pdf"
    ),

    file.path(
        FIG1DIR,
        "Figure1_FINAL_v2_600dpi.tiff"
    ),

    file.path(
        FIG2DIR,
        "Figure2_FINAL_v2.pdf"
    ),

    file.path(
        FIG2DIR,
        "Figure2_FINAL_v2_600dpi.tiff"
    )
)


output_status <- tibble(

    File =
        output_files,

    Exists =
        file.exists(
            output_files
        ),

    Size_bytes =
        ifelse(
            file.exists(
                output_files
            ),
            file.info(
                output_files
            )$size,
            NA_real_
        )
)


write_tsv(
    output_status,
    file.path(
        OUTDIR,
        "10B_v2_output_file_QC.tsv"
    )
)


if (
    any(
        !output_status$Exists
    )
) {

    stop(
        "Some final Figure1/2 v2 files were not created."
    )
}


# ============================================================
# 20. Session info
# ============================================================

capture.output(
    sessionInfo(),
    file =
        file.path(
            OUTDIR,
            "10B_v2_sessionInfo.txt"
        )
)


# ============================================================
# 21. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 10B v2 COMPLETED\n")
cat("============================================================\n")


cat("\nFIGURE 1\n")

cat(
    "Metadata rows                    : ",
    nrow(meta),
    "\n",
    sep = ""
)

cat(
    "RNA PCA samples                  : ",
    nrow(rna_pca$data),
    "\n",
    sep = ""
)

cat(
    "ATAC PCA samples                 : ",
    nrow(atac_pca$data),
    "\n",
    sep = ""
)

cat(
    "RNA PC1 / PC2                    : ",
    sprintf(
        "%.2f%% / %.2f%%",
        rna_pca$PC1,
        rna_pca$PC2
    ),
    "\n",
    sep = ""
)

cat(
    "ATAC PC1 / PC2                   : ",
    sprintf(
        "%.2f%% / %.2f%%",
        atac_pca$PC1,
        atac_pca$PC2
    ),
    "\n",
    sep = ""
)

cat(
    "Mean RNA P348-P350 correlation   : ",
    sprintf(
        "%.4f",
        mean(
            rna_cor$Correlation
        )
    ),
    "\n",
    sep = ""
)

cat(
    "Mean ATAC P348-P350 correlation  : ",
    sprintf(
        "%.4f",
        mean(
            atac_cor$Correlation
        )
    ),
    "\n",
    sep = ""
)


cat("\nFIGURE 2\n")

print(
    lrt_qc,
    n = Inf,
    width = Inf
)


cat(
    "RNA HC features                  : ",
    sum(
        rna_hc_counts$N
    ),
    "\n",
    sep = ""
)

cat(
    "ATAC HC features                 : ",
    sum(
        atac_hc_counts$N
    ),
    "\n",
    sep = ""
)


cat("\nTau-SPM merge QC:\n")

print(
    tau_spm_merge_qc,
    n = Inf,
    width = Inf
)


cat("\nATAC excluded Tau-SPM features:\n")

if (
    nrow(
        atac_tau_spm_obj$issues
    ) == 0
) {

    cat("None\n")

} else {

    print(
        atac_tau_spm_obj$issues,
        n = Inf,
        width = Inf
    )
}


cat("\nFinal files:\n")

print(
    output_status,
    n = Inf,
    width = Inf
)


cat("\nIMPORTANT:\n")

cat(
    "No biological threshold or upstream result was redefined.\n"
)

cat(
    "HC membership comes directly from finalized Step03 files.\n"
)

cat(
    "Tau-SPM exclusions are explicitly documented.\n"
)

cat("============================================================\n")

