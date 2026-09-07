suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

cat("\n")
cat("============================================================\n")
cat("STEP10F3B2b — SUPPLEMENTARY FIGURE S2 FINAL v2\n")
cat("Layout refinement only\n")
cat("No permutation rerun; no biological result modification\n")
cat("============================================================\n\n")


# ============================================================
# 1. Authoritative frozen sources
# ============================================================

SRC <- "05B_same_tissue_permutation_v2"

global_null_file <-
    file.path(SRC, "05B_global_permutation_null.tsv")

global_summary_file <-
    file.path(SRC, "05B_global_permutation_summary.tsv")

tissue_null_file <-
    file.path(SRC, "05B_tissue_permutation_null.tsv")

tissue_summary_file <-
    file.path(SRC, "05B_tissue_permutation_summary.tsv")

required_files <- c(
    global_null_file,
    global_summary_file,
    tissue_null_file,
    tissue_summary_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop(
        "Missing frozen Step05B-v2 file(s):\n",
        paste(missing_files, collapse = "\n")
    )
}


# ============================================================
# 2. Output
# ============================================================

OUTDIR <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3B_Supplementary_Figures/",
    "FigureS2_same_tissue_permutation_v2"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. Read frozen data
# ============================================================

global_null <- read_tsv(
    global_null_file,
    show_col_types = FALSE
)

global_summary <- read_tsv(
    global_summary_file,
    show_col_types = FALSE
)

tissue_null <- read_tsv(
    tissue_null_file,
    show_col_types = FALSE
)

tissue_summary <- read_tsv(
    tissue_summary_file,
    show_col_types = FALSE
)


# ============================================================
# 4. Frozen-result QC
# ============================================================

get_observed <- function(stat_name) {

    z <- global_summary %>%
        filter(Statistic == stat_name)

    if (nrow(z) != 1) {
        stop(
            "Cannot uniquely identify statistic: ",
            stat_name
        )
    }

    z$Observed[[1]]
}


obs_median <- get_observed(
    "Median_rho_mean"
)

obs_positive <- get_observed(
    "Positive_all_three_count"
)

obs_moderate <- get_observed(
    "Moderate_03_03_05_count"
)

obs_strong <- get_observed(
    "Strong_05_05_07_count"
)

obs_verystrong <- get_observed(
    "Very_strong_07_07_08_count"
)


frozen_qc <- tibble(
    Metric = c(
        "Global_permutation_N",
        "Tissue_permutation_rows",
        "Tissue_summary_rows",
        "Observed_positive_count",
        "Observed_moderate_count",
        "Observed_strong_count",
        "Observed_very_strong_count"
    ),

    Value = c(
        nrow(global_null),
        nrow(tissue_null),
        nrow(tissue_summary),
        obs_positive,
        obs_moderate,
        obs_strong,
        obs_verystrong
    ),

    Expected = c(
        10000,
        80000,
        8,
        3916,
        1778,
        702,
        226
    )
) %>%
    mutate(
        Status = if_else(
            abs(Value - Expected) < 1e-10,
            "PASS",
            "FAIL"
        )
    )

write_tsv(
    frozen_qc,
    file.path(
        OUTDIR,
        "FigureS2_v2_frozen_source_QC.tsv"
    )
)

if (any(frozen_qc$Status != "PASS")) {
    print(frozen_qc)
    stop(
        "Frozen-source QC FAILED. ",
        "S2 v2 will not be generated."
    )
}


# ============================================================
# 5. Tissue order
# ============================================================

tissue_order <- c(
    "Adipose",
    "Cerebellum",
    "Cortex",
    "Hypothalamus",
    "Liver",
    "Lung",
    "Muscle",
    "Spleen"
)

tissue_summary <- tissue_summary %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )

tissue_null <- tissue_null %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


# ============================================================
# 6. Publication theme
# ============================================================

theme_pub <- theme_classic(base_size = 10) +
    theme(
        text = element_text(
            family = "sans",
            colour = "black"
        ),

        axis.title = element_text(
            size = 10
        ),

        axis.text = element_text(
            size = 8.5,
            colour = "black"
        ),

        plot.title = element_text(
            size = 10.5,
            face = "bold",
            hjust = 0
        ),

        plot.subtitle = element_text(
            size = 8.4,
            colour = "grey25",
            hjust = 0
        ),

        legend.title = element_text(
            size = 8.5
        ),

        legend.text = element_text(
            size = 8
        ),

        # Same tag geometry for all four panels
        plot.tag = element_text(
            face = "bold",
            size = 15,
            colour = "black"
        ),

        plot.tag.position = c(
            0.005,
            0.995
        ),

        plot.margin = margin(
            8, 9, 8, 10
        )
    )


# ============================================================
# 7. PANEL A
# Exact discrete-frequency distribution
# ============================================================

rho_summary <- global_summary %>%
    filter(
        Statistic == "Median_rho_mean"
    )

if (nrow(rho_summary) != 1) {
    stop(
        "Median_rho_mean global summary row not unique."
    )
}

# The statistic is discrete because correlation is based on
# eight tissues. Count actual frozen values instead of
# imposing arbitrary histogram bins.
panelA_freq <- global_null %>%
    transmute(
        Median_rho_exact = round(
            Median_rho_mean,
            digits = 12
        )
    ) %>%
    count(
        Median_rho_exact,
        name = "Permutation_count"
    ) %>%
    arrange(
        Median_rho_exact
    )


unique_x <- sort(
    unique(
        panelA_freq$Median_rho_exact
    )
)

dx <- diff(unique_x)
dx <- dx[dx > 0]

if (length(dx) > 0) {

    bar_width <- min(dx) * 0.72

} else {

    bar_width <- 0.001
}


panelA <- ggplot(
    panelA_freq,
    aes(
        x = Median_rho_exact,
        y = Permutation_count
    )
) +

    geom_col(
        width = bar_width,
        fill = "grey78",
        colour = "grey55",
        linewidth = 0.25
    ) +

    geom_vline(
        xintercept =
            rho_summary$Null_mean,
        linetype = "dashed",
        linewidth = 0.6,
        colour = "grey25"
    ) +

    geom_vline(
        xintercept =
            rho_summary$Observed,
        linewidth = 0.9,
        colour = "#B2182B"
    ) +

    annotate(
        "text",
        x = rho_summary$Observed,
        y = Inf,
        label = sprintf(
            "Observed = %.3f",
            rho_summary$Observed
        ),
        hjust = 1.05,
        vjust = 1.45,
        size = 3.0,
        colour = "#B2182B"
    ) +

    annotate(
        "text",
        x = rho_summary$Null_mean,
        y = Inf,
        label = sprintf(
            "Null mean = %.3f",
            rho_summary$Null_mean
        ),
        hjust = -0.05,
        vjust = 3.0,
        size = 2.8,
        colour = "grey25"
    ) +

    labs(
        tag = "A",
        title =
            "Global cross-tissue correlation",
        subtitle =
            "Exact same-tissue permutation frequencies (10,000 permutations)",
        x =
            "Median Spearman correlation",
        y =
            "Permutation count"
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0, 0.08)
        )
    ) +

    theme_pub


# ============================================================
# 8. PANEL B
# Threshold-dependent enrichment
# ============================================================

threshold_patterns <- c(
    Positive =
        "Positive_all_three_count",

    Moderate =
        "Moderate_03_03_05_count",

    Strong =
        "Strong_05_05_07_count",

    `Very strong` =
        "Very_strong_07_07_08_count"
)


threshold_data <- bind_rows(
    lapply(
        names(threshold_patterns),
        function(label) {

            statistic <-
                threshold_patterns[[label]]

            z <- global_summary %>%
                filter(
                    Statistic == statistic
                )

            if (nrow(z) != 1) {
                stop(
                    "Cannot identify global statistic: ",
                    statistic
                )
            }

            tibble(
                Threshold = label,
                Observed = z$Observed,
                Null_mean = z$Null_mean,
                Fold =
                    z$Observed /
                    z$Null_mean,
                Empirical_P =
                    z$Empirical_P_upper
            )
        }
    )
) %>%
    mutate(
        Threshold = factor(
            Threshold,
            levels = c(
                "Positive",
                "Moderate",
                "Strong",
                "Very strong"
            )
        )
    )


panelB <- ggplot(
    threshold_data,
    aes(
        x = Threshold,
        y = Fold
    )
) +

    geom_hline(
        yintercept = 1,
        linetype = "dashed",
        linewidth = 0.5,
        colour = "grey50"
    ) +

    geom_segment(
        aes(
            xend = Threshold,
            y = 1,
            yend = Fold
        ),
        linewidth = 0.7,
        colour = "grey45"
    ) +

    geom_point(
        size = 3.0,
        shape = 21,
        fill = "#2166AC",
        colour = "black",
        stroke = 0.35
    ) +

    geom_text(
        aes(
            label = sprintf(
                "%.2fx",
                Fold
            )
        ),
        vjust = -0.85,
        size = 3
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0.03, 0.16)
        )
    ) +

    labs(
        tag = "B",
        title =
            "Threshold-dependent enrichment",
        subtitle =
            "Observed relative to same-tissue permutation null",
        x = NULL,
        y =
            "Observed / null mean"
    ) +

    theme_pub +

    theme(
        axis.text.x = element_text(
            angle = 20,
            hjust = 1
        )
    )


# ============================================================
# 9. PANEL C
# Tissue-specific Strong counts
# ============================================================

tissue_null_strong <- tissue_null %>%
    group_by(Tissue) %>%
    summarise(
        Null_mean =
            mean(
                Strong_05_05_07_count,
                na.rm = TRUE
            ),

        Null_Q025 =
            quantile(
                Strong_05_05_07_count,
                0.025,
                na.rm = TRUE
            ),

        Null_Q975 =
            quantile(
                Strong_05_05_07_count,
                0.975,
                na.rm = TRUE
            ),

        .groups = "drop"
    )


panelC_data <- tissue_summary %>%
    select(
        Tissue,
        N_pairs,

        Observed =
            Observed_strong_05_05_07,

        Strong_fold_enrichment,

        FDR =
            FDR_strong_05_05_07
    ) %>%

    left_join(
        tissue_null_strong,
        by = "Tissue"
    )


panelC <- ggplot(
    panelC_data,
    aes(
        y = Tissue
    )
) +

    # Horizontal 95% null interval.
    # geom_segment avoids the deprecated geom_errorbarh().
    geom_segment(
        aes(
            x = Null_Q025,
            xend = Null_Q975,
            yend = Tissue
        ),
        linewidth = 0.75,
        colour = "grey45"
    ) +

    # Small caps at interval boundaries
    geom_point(
        aes(
            x = Null_Q025
        ),
        shape = 124,
        size = 3.4,
        colour = "grey45"
    ) +

    geom_point(
        aes(
            x = Null_Q975
        ),
        shape = 124,
        size = 3.4,
        colour = "grey45"
    ) +

    # Null mean
    geom_point(
        aes(
            x = Null_mean
        ),
        shape = 21,
        size = 2.7,
        fill = "white",
        colour = "grey20",
        stroke = 0.6
    ) +

    # Observed
    geom_point(
        aes(
            x = Observed
        ),
        shape = 21,
        size = 3.15,
        fill = "#B2182B",
        colour = "black",
        stroke = 0.35
    ) +

    labs(
        tag = "C",
        title =
            "Tissue-specific Strong-pair enrichment",
        subtitle =
            "Grey lines: null 95% interval; open circles: null mean; red: observed",
        x =
            "Strong peak-gene pair count",
        y = NULL
    ) +

    theme_pub +

    theme(
        panel.grid.major.x =
            element_line(
                colour = "grey92",
                linewidth = 0.3
            )
    )


# ============================================================
# 10. PANEL D
# Tissue-level permutation FDR
# ============================================================

panelD_data <- tissue_summary %>%
    select(
        Tissue,

        `Median rho` =
            FDR_median_rho_mean,

        `Positive in all three` =
            FDR_positive_all_three,

        `Strong coupling` =
            FDR_strong_05_05_07
    ) %>%

    pivot_longer(
        cols = -Tissue,
        names_to = "Statistic",
        values_to = "FDR"
    ) %>%

    mutate(
        Statistic = factor(
            Statistic,
            levels = c(
                "Median rho",
                "Positive in all three",
                "Strong coupling"
            )
        ),

        minus_log10_FDR =
            -log10(
                pmax(
                    FDR,
                    1e-4
                )
            ),

        label = case_when(
            is.na(FDR) ~ "NA",
            FDR < 0.001 ~ "<0.001",
            TRUE ~ sprintf(
                "%.3f",
                FDR
            )
        )
    )


panelD <- ggplot(
    panelD_data,
    aes(
        x = Statistic,
        y = Tissue,
        fill = minus_log10_FDR
    )
) +

    geom_tile(
        colour = "white",
        linewidth = 0.7
    ) +

    geom_text(
        aes(
            label = label
        ),
        size = 2.75
    ) +

    scale_fill_gradient(
        low = "white",
        high = "#2166AC",
        name = "-log10(FDR)"
    ) +

    labs(
        tag = "D",
        title =
            "Tissue-level permutation significance",
        subtitle =
            "Numbers inside cells are BH-adjusted permutation FDR",
        x = NULL,
        y = NULL
    ) +

    theme_minimal(base_size = 10) +

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
                size = 8.2,
                colour = "black"
            ),

        axis.text.y =
            element_text(
                size = 8.5,
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
                size = 8.4,
                colour = "grey25",
                hjust = 0
            ),

        legend.title =
            element_text(
                size = 8.5
            ),

        legend.text =
            element_text(
                size = 8
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

        plot.margin =
            margin(
                8, 9, 8, 10
            )
    )


# ============================================================
# 11. Assemble unchanged 2 x 2 layout
# ============================================================

final_plot <-
    (
        panelA |
        panelB
    ) /
    (
        panelC |
        panelD
    ) +

    plot_layout(
        heights = c(
            1.0,
            1.12
        ),

        widths = c(
            1.0,
            1.0
        )
    )


# ============================================================
# 12. Export FINAL v2
# ============================================================

pdf_file <- file.path(
    OUTDIR,
    "Supplementary_Figure_S2_FINAL_v2.pdf"
)

tiff_file <- file.path(
    OUTDIR,
    "Supplementary_Figure_S2_FINAL_v2_600dpi.tiff"
)


ggsave(
    pdf_file,
    final_plot,
    device = cairo_pdf,
    width = 11.2,
    height = 8.4,
    units = "in"
)

ggsave(
    tiff_file,
    final_plot,
    device = "tiff",
    width = 11.2,
    height = 8.4,
    units = "in",
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 13. Export exact plotting datasets
# ============================================================

write_tsv(
    panelA_freq,
    file.path(
        OUTDIR,
        "FigureS2A_exact_permutation_frequency.tsv"
    )
)

write_tsv(
    threshold_data,
    file.path(
        OUTDIR,
        "FigureS2B_threshold_enrichment.tsv"
    )
)

write_tsv(
    panelC_data,
    file.path(
        OUTDIR,
        "FigureS2C_tissue_strong_null_interval.tsv"
    )
)

write_tsv(
    panelD_data,
    file.path(
        OUTDIR,
        "FigureS2D_tissue_FDR.tsv"
    )
)


# ============================================================
# 14. Change log / biological lock
# ============================================================

change_log <- tibble(
    Item = c(
        "Panel_A_visualization",
        "Panel_letters",
        "Panel_C_subtitle",
        "Deprecated_geom_errorbarh",
        "Permutation_rerun",
        "Coupling_recalculated",
        "Thresholds_modified",
        "Biological_results_modified"
    ),

    Previous = c(
        "Histogram with arbitrary bins",
        "Mixed positioning",
        "Grey: permutation 95% interval; red: observed",
        "Used",
        "No",
        "No",
        "No",
        "No"
    ),

    FINAL_v2 = c(
        "Exact frozen-value frequency bars",
        "Unified A-D positioning",
        "Grey lines: null 95% interval; open circles: null mean; red: observed",
        "Removed",
        "No",
        "No",
        "No",
        "No"
    )
)

write_tsv(
    change_log,
    file.path(
        OUTDIR,
        "FigureS2_FINAL_v2_change_log.tsv"
    )
)


# ============================================================
# 15. File QC
# ============================================================

files <- c(
    pdf_file,
    tiff_file
)

file_qc <- tibble(
    File = files,
    Exists = file.exists(files),
    Size_bytes =
        file.info(files)$size,
    Status =
        if_else(
            file.exists(files) &
            file.info(files)$size > 0,
            "PASS",
            "FAIL"
        )
)

write_tsv(
    file_qc,
    file.path(
        OUTDIR,
        "FigureS2_FINAL_v2_file_QC.tsv"
    )
)


# ============================================================
# 16. FINAL report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10F3B2b COMPLETED\n")
cat("============================================================\n\n")

cat("Frozen source QC:\n")
print(frozen_qc)

cat("\nFINAL v2 modifications:\n")
print(change_log)

cat("\nOutput files:\n")
print(file_qc)

cat("\nPanel A exact null values:\n")
cat(
    "Unique median-rho values : ",
    nrow(panelA_freq),
    "\n",
    sep = ""
)

cat(
    "Total permutations       : ",
    sum(panelA_freq$Permutation_count),
    "\n",
    sep = ""
)

cat("\nIMPORTANT:\n")
cat("- Same authoritative Step05B-v2 source was used.\n")
cat("- No permutation was rerun.\n")
cat("- No coupling statistic was recalculated.\n")
cat("- No candidate threshold was changed.\n")
cat("- Panel A changed visualization only: histogram -> exact-value frequency bars.\n")
cat("- Panel letters A-D now use identical placement rules.\n")
cat("- Panel C now explicitly identifies the null mean.\n")
cat("- Biological results modified: 0.\n")

if (
    all(frozen_qc$Status == "PASS") &&
    all(file_qc$Status == "PASS")
) {

    cat("\nSTEP10F3B2b STATUS: PASS\n")
    cat("Supplementary Figure S2 FINAL v2 is ready for visual inspection.\n")

} else {

    cat("\nSTEP10F3B2b STATUS: CHECK\n")
}

cat("============================================================\n")
