suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

cat("\n")
cat("============================================================\n")
cat("STEP10F3B3c — SUPPLEMENTARY FIGURE S3 FINAL v3\n")
cat("Publication-label/layout refinement only\n")
cat("Frozen S3-v2 plotting data only\n")
cat("No biological analysis rerun\n")
cat("============================================================\n\n")


# ============================================================
# 1. Authoritative frozen S3-v2 plotting sources
# ============================================================

V2DIR <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3B_Supplementary_Figures/",
    "FigureS3_genomic_architecture_v2"
)

F_A <- file.path(
    V2DIR,
    "FigureS3A_UROPA_relative_location.tsv"
)

F_B <- file.path(
    V2DIR,
    "FigureS3B_structural_classes.tsv"
)

F_C <- file.path(
    V2DIR,
    "FigureS3C_structural_enrichment.tsv"
)

F_D <- file.path(
    V2DIR,
    "FigureS3D_coupling_by_structure.tsv"
)

F_E <- file.path(
    V2DIR,
    "FigureS3E_flanking_distance.tsv"
)

F_D_ZERO <- file.path(
    V2DIR,
    "FigureS3D_zero_count_classes_omitted_from_plot.tsv"
)

F_D_QC <- file.path(
    V2DIR,
    "FigureS3D_zero_count_omission_QC.tsv"
)

required_files <- c(
    F_A,
    F_B,
    F_C,
    F_D,
    F_E,
    F_D_ZERO,
    F_D_QC
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing frozen S3-v2 source file(s):\n",
        paste(
            missing_files,
            collapse = "\n"
        )
    )
}


# ============================================================
# 2. Output
# ============================================================

OUTDIR <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3B_Supplementary_Figures/",
    "FigureS3_genomic_architecture_FINAL_v3"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. Read frozen v2 plotting tables
# ============================================================

rel <- read_tsv(
    F_A,
    show_col_types = FALSE
)

struct <- read_tsv(
    F_B,
    show_col_types = FALSE
)

enrich <- read_tsv(
    F_C,
    show_col_types = FALSE
)

coupling <- read_tsv(
    F_D,
    show_col_types = FALSE
)

flank <- read_tsv(
    F_E,
    show_col_types = FALSE
)

zero_classes <- read_tsv(
    F_D_ZERO,
    show_col_types = FALSE
)

zero_qc <- read_tsv(
    F_D_QC,
    show_col_types = FALSE
)


# ============================================================
# 4. Frozen source QC
# ============================================================

qc <- tibble(
    Metric = c(
        "S3A_pair_sum",
        "S3B_pair_sum",
        "S3C_test_rows",
        "S3D_plotted_pair_sum",
        "S3D_zero_count_classes",
        "S3D_removed_observations",
        "S3E_all_flanking_N"
    ),

    Value = c(
        sum(
            rel$Number_of_pairs,
            na.rm = TRUE
        ),

        sum(
            struct$Number_of_pairs,
            na.rm = TRUE
        ),

        nrow(enrich),

        sum(
            coupling$N_pairs,
            na.rm = TRUE
        ),

        nrow(zero_classes),

        sum(
            zero_classes$N_pairs,
            na.rm = TRUE
        ),

        flank %>%
            filter(
                Group == "All_flanking"
            ) %>%
            pull(N) %>%
            first()
    ),

    Expected = c(
        702,
        702,
        18,
        702,
        2,
        0,
        221
    )
) %>%

    mutate(
        Status = if_else(
            Value == Expected,
            "PASS",
            "FAIL"
        )
    )


write_tsv(
    qc,
    file.path(
        OUTDIR,
        "FigureS3_FINAL_v3_frozen_source_QC.tsv"
    )
)

if (any(qc$Status != "PASS")) {

    print(qc)

    stop(
        "Frozen S3-v2 QC FAILED. ",
        "FINAL v3 was not generated."
    )
}


# ============================================================
# 5. Theme
# ============================================================

theme_pub <- theme_classic(
    base_size = 10
) +

    theme(
        text = element_text(
            family = "sans",
            colour = "black"
        ),

        axis.text = element_text(
            size = 8.4,
            colour = "black"
        ),

        axis.title = element_text(
            size = 9.5
        ),

        plot.title = element_text(
            size = 10.5,
            face = "bold",
            hjust = 0
        ),

        plot.subtitle = element_text(
            size = 8.2,
            colour = "grey25",
            hjust = 0
        ),

        plot.tag = element_text(
            face = "bold",
            size = 15,
            colour = "black"
        ),

        plot.tag.position = c(
            0.005,
            0.995
        ),

        legend.title = element_text(
            size = 8.5
        ),

        legend.text = element_text(
            size = 8
        ),

        plot.margin = margin(
            8,
            9,
            8,
            10
        )
    )


# ============================================================
# 6. PANEL A
# Same data, only additional upper headroom
# ============================================================

relative_order <- c(
    "PeakInsideFeature",
    "Upstream",
    "OverlapStart",
    "Downstream",
    "OverlapEnd"
)

rel_plot <- rel %>%

    mutate(
        UROPA_relative_location =
            factor(
                UROPA_relative_location,
                levels =
                    relative_order
            )
    )


panelA <- ggplot(
    rel_plot,
    aes(
        x = UROPA_relative_location,
        y = Number_of_pairs
    )
) +

    geom_col(
        width = 0.72,
        fill = "#4C78A8"
    ) +

    geom_text(
        aes(
            label = paste0(
                Number_of_pairs,
                "\n(",
                sprintf(
                    "%.1f",
                    Percent
                ),
                "%)"
            )
        ),

        vjust = -0.25,
        size = 2.75
    ) +

    # FINAL v3:
    # additional upper headroom for 395 / 56.3% label
    scale_y_continuous(
        expand = expansion(
            mult = c(
                0,
                0.21
            )
        )
    ) +

    labs(
        tag = "A",

        title =
            "UROPA relative-location categories",

        subtitle =
            "702 strongly coupled tissue-specific peak-gene candidates",

        x = NULL,

        y =
            "Peak-gene pairs"
    ) +

    theme_pub +

    theme(
        axis.text.x =
            element_text(
                angle = 25,
                hjust = 1
            )
    )


# ============================================================
# 7. PANEL B
# Fixed biologically interpretable structural order
# ============================================================

structural_order_top_to_bottom <- c(
    "Gene-overlapping",
    "Flanking 0–2 kb",
    "Flanking 2–20 kb",
    "Flanking 20–100 kb",
    "Flanking >100 kb",
    "Other or NA"
)


struct_plot <- struct %>%

    mutate(
        Class_label = as.character(
            Class_label
        ),

        # factor levels reversed because first factor level
        # appears at the bottom of a ggplot discrete y-axis
        Class_label = factor(
            Class_label,
            levels = rev(
                structural_order_top_to_bottom
            )
        )
    )


panelB <- ggplot(
    struct_plot,
    aes(
        y = Class_label,
        x = Percent
    )
) +

    geom_col(
        width = 0.65,
        fill = "#72A06A"
    ) +

    geom_text(
        aes(
            label = sprintf(
                "%d (%.1f%%)",
                Number_of_pairs,
                Percent
            )
        ),

        hjust = -0.08,
        size = 2.8
    ) +

    scale_x_continuous(
        limits = c(
            0,
            max(
                struct_plot$Percent,
                na.rm = TRUE
            ) * 1.25
        ),

        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    labs(
        tag = "B",

        title =
            "Coarse peak-gene structural classes",

        subtitle =
            "Distance/overlap-based Step06 definition",

        x =
            "Percentage of strong pairs",

        y = NULL
    ) +

    theme_pub


# ============================================================
# 8. PANEL C
# Clean comparison labels + explicit undefined cells
# ============================================================

comparison_display <- c(

    "Concordant_HC_vs_All" =
        "Concordant HC vs All",

    "Strong_vs_All" =
        "Strong vs All",

    "Strong_vs_Concordant_HC" =
        "Strong vs Concordant HC"
)


comparison_order <- c(
    "Concordant HC vs All",
    "Strong vs All",
    "Strong vs Concordant HC"
)


enrich_plot <- enrich %>%

    mutate(
        Class_label =
            as.character(
                Class_label
            ),

        Comparison_raw =
            as.character(
                Comparison
            ),

        Comparison_label =
            recode(
                Comparison_raw,
                !!!comparison_display
            ),

        Comparison_label =
            factor(
                Comparison_label,
                levels =
                    comparison_order
            ),

        Class_label =
            factor(
                Class_label,
                levels = rev(
                    structural_order_top_to_bottom
                )
            ),

        # Use frozen values if available.
        # Undefined cells are kept explicitly rather than dropped.
        Undefined =
            !is.finite(log2_fold) |
            is.na(log2_fold),

        FDR_label =
            case_when(
                Undefined ~ "",

                FDR < 0.001 ~ "***",

                FDR < 0.01 ~ "**",

                FDR < 0.05 ~ "*",

                TRUE ~ ""
            ),

        Undefined_label =
            case_when(

                # If both compared counts are explicitly zero,
                # state N=0.
                "Selected_class_N" %in%
                    names(.) &
                "Background_class_N" %in%
                    names(.) &
                Undefined &
                Selected_class_N == 0 &
                Background_class_N == 0
                    ~ "N=0",

                Undefined ~ "NA",

                TRUE ~ ""
            )
    )


# If dplyr parsing of the conditional column references above
# encounters an older version issue, rebuild Undefined_label
# safely here.
if (!"Undefined_label" %in% colnames(enrich_plot)) {

    enrich_plot$Undefined_label <-
        ifelse(
            enrich_plot$Undefined,
            "NA",
            ""
        )
}


# ------------------------------------------------------------
# Split defined and undefined cells
# ------------------------------------------------------------

enrich_defined <- enrich_plot %>%
    filter(
        !Undefined
    )

enrich_undefined <- enrich_plot %>%
    filter(
        Undefined
    )


panelC <- ggplot() +

    # Defined enrichment cells
    geom_tile(
        data = enrich_defined,

        aes(
            x = Comparison_label,
            y = Class_label,
            fill = log2_fold
        ),

        colour = "white",
        linewidth = 0.55
    ) +

    # Explicit undefined/NA cells
    geom_tile(
        data = enrich_undefined,

        aes(
            x = Comparison_label,
            y = Class_label
        ),

        fill = "grey70",
        colour = "white",
        linewidth = 0.55
    ) +

    # Significance stars only for defined cells
    geom_text(
        data = enrich_defined,

        aes(
            x = Comparison_label,
            y = Class_label,
            label = FDR_label
        ),

        size = 3.4,
        fontface = "bold"
    ) +

    # Explicit NA/N=0 labels
    geom_text(
        data = enrich_undefined,

        aes(
            x = Comparison_label,
            y = Class_label,
            label = Undefined_label
        ),

        size = 2.65,
        colour = "grey20",
        fontface = "bold"
    ) +

    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        name = "log2 fold"
    ) +

    labs(
        tag = "C",

        title =
            "Structural enrichment sensitivity",

        subtitle = paste0(
            "* FDR<0.05; ** FDR<0.01; *** FDR<0.001; ",
            "grey cells indicate undefined enrichment"
        ),

        x = NULL,
        y = NULL
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        text = element_text(
            family = "sans",
            colour = "black"
        ),

        panel.grid =
            element_blank(),

        axis.text.x =
            element_text(
                angle = 25,
                hjust = 1,
                size = 8.0,
                colour = "black"
            ),

        axis.text.y =
            element_text(
                size = 8.2,
                colour = "black"
            ),

        plot.title =
            element_text(
                size = 10.5,
                face = "bold",
                hjust = 0
            ),

        plot.subtitle =
            element_text(
                size = 8.0,
                colour = "grey25",
                hjust = 0
            ),

        plot.tag =
            element_text(
                face = "bold",
                size = 15,
                colour = "black"
            ),

        plot.tag.position = c(
            0.005,
            0.995
        ),

        legend.title =
            element_text(
                size = 8.3
            ),

        legend.text =
            element_text(
                size = 8
            ),

        plot.margin = margin(
            8,
            9,
            8,
            10
        )
    )


# ============================================================
# 9. PANEL D
# Keep v2 scientific visualization unchanged
# ============================================================

coupling_plot <- coupling %>%

    mutate(
        Class_label =
            as.character(
                Class_label
            ),

        Class_label = factor(
            Class_label,
            levels =
                unique(
                    Class_label
                )
        )
    )


panelD <- ggplot(
    coupling_plot,
    aes(
        y = Class_label
    )
) +

    geom_segment(
        aes(
            x = Q25_rho_mean,
            xend = Q75_rho_mean,
            yend = Class_label
        ),

        linewidth = 0.8,
        colour = "grey45"
    ) +

    geom_point(
        aes(
            x = Median_rho_mean,
            size = N_pairs
        ),

        shape = 21,
        fill = "#E69F00",
        colour = "black",
        stroke = 0.35
    ) +

    scale_size_continuous(
        range = c(
            2.2,
            5.0
        )
    ) +

    labs(
        tag = "D",

        title =
            "RNA-ATAC coupling by structural class",

        subtitle =
            "Point: median rho; line: interquartile range; N=0 classes omitted",

        x =
            "Mean-profile Spearman correlation",

        y = NULL,

        size =
            "Pairs"
    ) +

    theme_pub


# ============================================================
# 10. PANEL E
# Keep v2 scientific visualization unchanged
# ============================================================

flank_plot <- flank %>%

    mutate(
        Group_label =
            as.character(
                Group_label
            ),

        Group_label = factor(
            Group_label,
            levels =
                unique(
                    Group_label
                )
        )
    )


panelE <- ggplot(
    flank_plot,
    aes(
        y = Group_label
    )
) +

    geom_segment(
        aes(
            x = Q25_kb,
            xend = Q75_kb,
            yend = Group_label
        ),

        linewidth = 0.85,
        colour = "grey45"
    ) +

    geom_point(
        aes(
            x = Median_kb
        ),

        shape = 21,
        size = 3.2,
        fill = "#CC6677",
        colour = "black",
        stroke = 0.35
    ) +

    geom_text(
        aes(
            x = Q75_kb,
            label = paste0(
                "n=",
                N
            )
        ),

        hjust = -0.25,
        size = 2.8
    ) +

    scale_x_continuous(
        expand = expansion(
            mult = c(
                0.04,
                0.18
            )
        )
    ) +

    labs(
        tag = "E",

        title =
            "Flanking-distance characteristics",

        subtitle =
            "Median and interquartile range",

        x =
            "Peak-gene distance (kb)",

        y = NULL
    ) +

    theme_pub


# ============================================================
# 11. Assemble — same overall layout
# ============================================================

top_row <-
    panelA |
    panelB

middle_row <-
    panelC |
    panelD

final_plot <-
    top_row /
    middle_row /
    panelE +

    plot_layout(
        heights = c(
            1.0,
            1.15,
            0.78
        )
    )


# ============================================================
# 12. Export FINAL v3
# ============================================================

PDF_OUT <- file.path(
    OUTDIR,
    "Supplementary_Figure_S3_FINAL_v3.pdf"
)

TIFF_OUT <- file.path(
    OUTDIR,
    "Supplementary_Figure_S3_FINAL_v3_600dpi.tiff"
)


ggsave(
    filename = PDF_OUT,
    plot = final_plot,
    device = cairo_pdf,
    width = 11.4,
    height = 11.0,
    units = "in"
)

ggsave(
    filename = TIFF_OUT,
    plot = final_plot,
    device = "tiff",
    width = 11.4,
    height = 11.0,
    units = "in",
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 13. Export exact v3 plotting tables
# ============================================================

write_tsv(
    rel_plot,
    file.path(
        OUTDIR,
        "FigureS3A_FINAL_v3.tsv"
    )
)

write_tsv(
    struct_plot,
    file.path(
        OUTDIR,
        "FigureS3B_FINAL_v3.tsv"
    )
)

write_tsv(
    enrich_plot,
    file.path(
        OUTDIR,
        "FigureS3C_FINAL_v3.tsv"
    )
)

write_tsv(
    coupling_plot,
    file.path(
        OUTDIR,
        "FigureS3D_FINAL_v3.tsv"
    )
)

write_tsv(
    flank_plot,
    file.path(
        OUTDIR,
        "FigureS3E_FINAL_v3.tsv"
    )
)


# ============================================================
# 14. Explicit v3 change log
# ============================================================

change_log <- tibble(

    Panel = c(
        "S3A",
        "S3B",
        "S3C",
        "S3D",
        "S3E"
    ),

    FINAL_v3_change = c(
        "Increased upper y-axis headroom only",
        "Reordered classes by biological distance",
        paste0(
            "Removed internal underscores from comparison labels; ",
            "undefined cells explicitly grey-labelled"
        ),
        "No scientific or layout change",
        "No scientific or layout change"
    ),

    Biological_values_changed = c(
        "NO",
        "NO",
        "NO",
        "NO",
        "NO"
    )
)


write_tsv(
    change_log,
    file.path(
        OUTDIR,
        "FigureS3_FINAL_v3_change_log.tsv"
    )
)


# ============================================================
# 15. Undefined-cell audit
# ============================================================

undefined_audit <- enrich_plot %>%

    filter(
        Undefined
    ) %>%

    select(
        Comparison_raw,
        Comparison_label,
        Structural_class,
        Class_label,
        everything()
    )


write_tsv(
    undefined_audit,
    file.path(
        OUTDIR,
        "FigureS3C_undefined_cell_audit.tsv"
    )
)


# ============================================================
# 16. Final file QC
# ============================================================

file_qc <- tibble(

    File = c(
        PDF_OUT,
        TIFF_OUT
    ),

    Exists = file.exists(
        c(
            PDF_OUT,
            TIFF_OUT
        )
    )
) %>%

    mutate(
        Size_bytes =
            file.info(File)$size,

        Status =
            if_else(
                Exists &
                Size_bytes > 0,
                "PASS",
                "FAIL"
            )
    )


write_tsv(
    file_qc,
    file.path(
        OUTDIR,
        "FigureS3_FINAL_v3_file_QC.tsv"
    )
)


# ============================================================
# 17. Caption — FINAL v3 wording
# ============================================================

caption <- paste0(
    "Supplementary Figure S3. Complementary genomic architecture ",
    "of strongly coupled tissue-specific peak-gene candidates. ",
    "(A) UROPA relative-location categories among the 702 strongly ",
    "coupled peak-gene candidates. ",
    "(B) Coarse overlap- and distance-based structural classes ",
    "defined in Step06 and displayed in increasing distance from ",
    "gene-overlapping to distal categories. ",
    "(C) Structural-class enrichment across the frozen selection ",
    "comparisons. Asterisks indicate Benjamini-Hochberg-adjusted ",
    "significance (* FDR<0.05, ** FDR<0.01, *** FDR<0.001); ",
    "grey cells indicate comparisons for which enrichment was ",
    "undefined in the frozen source results. ",
    "(D) Mean-profile RNA-ATAC coupling according to structural ",
    "class; points represent median Spearman correlations and ",
    "horizontal lines represent interquartile ranges. Structural ",
    "classes with zero strong peak-gene candidates were retained ",
    "in the frozen source tables but omitted from this quantitative ",
    "coupling panel because correlation statistics were undefined. ",
    "(E) Distance characteristics of flanking peak-gene classes. ",
    "These analyses provide a complementary sensitivity analysis ",
    "using the earlier UROPA/distance-based structural framework; ",
    "the primary promoter-first, strand-aware TSS classification ",
    "is presented in Figure 4."
)


writeLines(
    caption,
    file.path(
        OUTDIR,
        "Supplementary_Figure_S3_FINAL_v3_caption.txt"
    )
)


# ============================================================
# 18. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10F3B3c COMPLETED\n")
cat("============================================================\n\n")

cat("Frozen source QC:\n")
print(qc)

cat("\nFINAL v3 change log:\n")
print(
    change_log,
    width = Inf
)

cat("\nUndefined Panel-C cells:\n")
cat(
    nrow(
        undefined_audit
    ),
    "\n"
)

if (nrow(undefined_audit) > 0) {

    print(
        undefined_audit %>%
            select(
                Comparison_label,
                Class_label,
                Undefined_label
            ),
        width = Inf
    )
}

cat("\nOutput QC:\n")
print(file_qc)

cat("\nIMPORTANT:\n")
cat("- Input is frozen S3-v2 plotting data only.\n")
cat("- Step06 and Step06B were NOT rerun.\n")
cat("- Strong candidates remain 702.\n")
cat("- Peak-gene assignments were NOT changed.\n")
cat("- Structural definitions were NOT changed.\n")
cat("- Panel A changed spacing only.\n")
cat("- Panel B changed display order only.\n")
cat("- Panel C changed display labels/NA annotation only.\n")
cat("- Panels D and E retain the v2 scientific visualization.\n")
cat("- Biological results modified: 0.\n")

if (
    all(qc$Status == "PASS") &&
    all(file_qc$Status == "PASS")
) {

    cat("\nSTEP10F3B3c STATUS: PASS\n")
    cat(
        "Supplementary Figure S3 FINAL v3 ",
        "is ready for final visual lock.\n"
    )

} else {

    cat("\nSTEP10F3B3c STATUS: CHECK\n")
}

cat("============================================================\n")
