suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

cat("\n")
cat("============================================================\n")
cat("STEP10F3B4D — SUPPLEMENTARY FIGURE S4 FINAL v3\n")
cat("Visualization refinement only\n")
cat("Frozen S4-FIX plotting tables only\n")
cat("No AME/STREME/Tomtom rerun\n")
cat("============================================================\n\n")


# ============================================================
# 1. Frozen S4-FIX plotting sources
# ============================================================

V1DIR <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3B_Supplementary_Figures/",
    "FigureS4_motif_analysis_QC_FIX"
)

F_A <- file.path(
    V1DIR,
    "FigureS4A_background_balance.tsv"
)

F_B <- file.path(
    V1DIR,
    "FigureS4B_genomic_class_matching.tsv"
)

F_C <- file.path(
    V1DIR,
    "FigureS4C_AME_supportive_ranking.tsv"
)

F_D <- file.path(
    V1DIR,
    "FigureS4D_STREME_formal_evidence.tsv"
)

F_E <- file.path(
    V1DIR,
    "FigureS4E_Tomtom_crossrun.tsv"
)

F_QC <- file.path(
    V1DIR,
    "FigureS4_frozen_source_QC.tsv"
)

required_files <- c(
    F_A,
    F_B,
    F_C,
    F_D,
    F_E,
    F_QC
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {
    stop(
        "Missing frozen S4 source file(s):\n",
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
    "FigureS4_motif_analysis_QC_FINAL_v3"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. Read frozen plotting tables
# ============================================================

bg_long <- read_tsv(
    F_A,
    show_col_types = FALSE
)

class_qc <- read_tsv(
    F_B,
    show_col_types = FALSE
)

ame_best <- read_tsv(
    F_C,
    show_col_types = FALSE
)

streme_source <- read_tsv(
    F_D,
    show_col_types = FALSE
)

tomtom <- read_tsv(
    F_E,
    show_col_types = FALSE
)

previous_qc <- read_tsv(
    F_QC,
    show_col_types = FALSE
)


# ============================================================
# 4. Tissue order
# ============================================================

tissue_order <- c(
    "Cerebellum",
    "Liver",
    "Muscle",
    "Spleen"
)

bg_long <- bg_long %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )

class_qc <- class_qc %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )

ame_best <- ame_best %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )

tomtom <- tomtom %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


# ============================================================
# 5. Frozen biological QC
# ============================================================

required_previous_pass <- all(
    previous_qc$Status == "PASS"
)

qc <- tibble(
    Metric = c(
        "Previous_S4_frozen_QC_all_PASS",
        "PanelA_rows",
        "PanelB_cells",
        "PanelC_AME_tissues",
        "PanelD_reported_STREME_motifs",
        "PanelE_Tomtom_rows",
        "AME_formal_E_lt_0.05_hits",
        "STREME_holdout_E_lt_0.05_hits",
        "Biological_results_modified"
    ),

    Value = c(
        as.integer(required_previous_pass),
        nrow(bg_long),
        nrow(class_qc),
        nrow(ame_best),
        nrow(streme_source),
        nrow(tomtom),
        sum(
            ame_best$E_value < 0.05,
            na.rm = TRUE
        ),
        sum(
            streme_source$Formal_E_lt_0_05,
            na.rm = TRUE
        ),
        0
    ),

    Expected = c(
        1,
        8,
        24,
        4,
        14,
        9,
        0,
        0,
        0
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
        "FigureS4_FINAL_v2_frozen_source_QC.tsv"
    )
)

if (any(qc$Status != "PASS")) {
    print(qc)
    stop(
        "Frozen S4 source QC failed. ",
        "FINAL v2 was not generated."
    )
}


# ============================================================
# 6. Theme
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
            size = 8.3,
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
            size = 8.0,
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
            size = 8.3
        ),

        legend.text = element_text(
            size = 8
        ),

        plot.margin = margin(
            8, 9, 8, 10
        )
    )


# ============================================================
# 7. PANEL A
# Unchanged
# ============================================================

panelA <- ggplot(
    bg_long,
    aes(
        x = SMD,
        y = Tissue,
        shape = Metric
    )
) +

    geom_rect(
        inherit.aes = FALSE,
        xmin = -0.10,
        xmax = 0.10,
        ymin = -Inf,
        ymax = Inf,
        fill = "grey95"
    ) +

    geom_vline(
        xintercept = c(
            -0.10,
            0.10
        ),
        linetype = "dashed",
        linewidth = 0.45,
        colour = "grey45"
    ) +

    geom_vline(
        xintercept = 0,
        linewidth = 0.35,
        colour = "grey65"
    ) +

    geom_point(
        size = 3.1,
        stroke = 0.75
    ) +

    scale_shape_manual(
        values = c(
            "GC content" = 16,
            "log2 peak length" = 17
        )
    ) +

    labs(
        tag = "A",

        title =
            "Matched-background covariate balance",

        subtitle =
            "Absolute standardized mean difference <0.10 in all tissues",

        x =
            "Standardized mean difference",

        y = NULL,

        shape = NULL
    ) +

    theme_pub


# ============================================================
# 8. PANEL B — FINAL v2
# Uniform cells; remove misleading colour gradient
# ============================================================

if (!all(class_qc$Difference == 0)) {
    stop(
        "Panel B contains non-zero FG-BG differences. ",
        "Uniform zero-difference visualization is not appropriate."
    )
}


panelB <- ggplot(
    class_qc,
    aes(
        x = Tissue,
        y = Genomic_class
    )
) +

    geom_tile(
        fill = "grey94",
        colour = "white",
        linewidth = 0.75
    ) +

    geom_text(
        aes(
            label = Difference
        ),
        size = 2.9,
        colour = "black"
    ) +

    labs(
        tag = "B",

        title =
            "Exact genomic-class matching",

        subtitle =
            "Exact matching; all foreground-background differences = 0",

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

        panel.grid = element_blank(),

        axis.text.x = element_text(
            angle = 20,
            hjust = 1,
            size = 8.3,
            colour = "black"
        ),

        axis.text.y = element_text(
            size = 8.0,
            colour = "black"
        ),

        plot.title = element_text(
            size = 10.5,
            face = "bold",
            hjust = 0
        ),

        plot.subtitle = element_text(
            size = 8.0,
            colour = "grey25",
            hjust = 0
        ),

        plot.tag = element_text(
            face = "bold",
            size = 15
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
# 9. PANEL C — FINAL v2
# Explicit descriptive threshold annotation
# ============================================================

ame_threshold <- -log10(
    0.05
)

panelC <- ggplot(
    ame_best,
    aes(
        x = Tissue,
        y = minus_log10_adjP
    )
) +

    geom_hline(
        yintercept = ame_threshold,
        linetype = "dashed",
        linewidth = 0.5,
        colour = "grey45"
    ) +

    geom_col(
        width = 0.62,
        fill = "#7A5195"
    ) +

    geom_text(
        aes(
            label = paste0(
                motif_alt_ID,
                "\n",
                motif_ID
            )
        ),
        vjust = -0.35,
        size = 2.65
    ) +

    geom_text(
        aes(
            label = paste0(
                "E=",
                format(
                    E_value,
                    digits = 2,
                    scientific = TRUE
                )
            )
        ),
        y = 0.04,
        vjust = 0,
        size = 2.3,
        colour = "grey25"
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(
                0,
                0.21
            )
        )
    ) +

    labs(
        tag = "C",

        title =
            "Supportive AME ranking",

        subtitle =
            "Dashed line: descriptive adj.P=0.05; formal AME criterion: E<0.05 (0 hits)",

        x = NULL,

        y =
            "-log10(AME adjusted P)"
    ) +

    theme_pub +

    theme(
        axis.text.x = element_text(
            angle = 18,
            hjust = 1
        )
    )


# ============================================================
# 10. PANEL D — FINAL v2
# Explicit NR cells
# ============================================================

# Keep the 14 reported motifs exactly as frozen.
streme_reported <- streme_source %>%
    mutate(
        Tissue = as.character(
            Tissue
        ),
        Motif_number = as.integer(
            Motif_number
        )
    )


# Build complete 4 x 5 display grid.
streme_grid <- expand_grid(
    Tissue = tissue_order,
    Motif_number = 1:5
) %>%

    left_join(
        streme_reported %>%
            select(
                Tissue,
                Motif_number,
                Statistic_type,
                E_value,
                S_score,
                Formal_E_lt_0_05
            ),

        by = c(
            "Tissue",
            "Motif_number"
        )
    ) %>%

    mutate(
        Report_status = case_when(

            is.na(
                Statistic_type
            ) ~
                "Not reported",

            Statistic_type ==
                "Holdout E-value" ~
                "Holdout E-value",

            Statistic_type ==
                "Training score only" ~
                "Training score only",

            TRUE ~
                "Other"
        ),

        Value_label = case_when(

            Report_status ==
                "Not reported" ~
                "NR",

            Report_status ==
                "Holdout E-value" ~
                paste0(
                    "E=",
                    format(
                        E_value,
                        digits = 2,
                        scientific = TRUE
                    )
                ),

            Report_status ==
                "Training score only" ~
                paste0(
                    "S=",
                    format(
                        S_score,
                        digits = 2,
                        scientific = TRUE
                    )
                ),

            TRUE ~
                ""
        ),

        Tissue = factor(
            Tissue,
            levels = tissue_order
        ),

        Motif_row = factor(
            paste0(
                "STREME-",
                Motif_number
            ),

            levels = rev(
                paste0(
                    "STREME-",
                    1:5
                )
            )
        )
    )


nr_count <- sum(
    streme_grid$Report_status ==
        "Not reported"
)

if (nr_count != 6) {
    stop(
        "Expected 6 NR cells in the 4x5 STREME display grid; observed ",
        nr_count
    )
}


panelD <- ggplot(
    streme_grid,
    aes(
        x = Tissue,
        y = Motif_row,
        fill = Report_status
    )
) +

    geom_tile(
        colour = "white",
        linewidth = 0.7,
        width = 0.92,
        height = 0.92
    ) +

    geom_text(
        aes(
            label = Value_label
        ),
        size = 2.35,
        colour = "black"
    ) +

    scale_fill_manual(
        values = c(
            "Holdout E-value" =
                "#E6AB02",

            "Training score only" =
                "grey72",

            "Not reported" =
                "grey95",

            "Other" =
                "white"
        ),

        breaks = c(
            "Holdout E-value",
            "Training score only",
            "Not reported"
        ),

        labels = c(
            "Holdout E-value",
            "Training score only",
            "Not reported (NR)"
        ),

        name =
            "Reported statistic"
    ) +

    labs(
        tag = "D",

        title =
            "Default STREME formal-evidence status",

        subtitle = paste0(
            "14 motifs reported; 9 had holdout E-values (0/9 E<0.05); ",
            "5 Cerebellum motifs had training scores only"
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

        panel.grid = element_blank(),

        axis.text.x = element_text(
            angle = 18,
            hjust = 1,
            size = 8.2,
            colour = "black"
        ),

        axis.text.y = element_text(
            size = 8.0,
            colour = "black"
        ),

        plot.title = element_text(
            size = 10.5,
            face = "bold",
            hjust = 0
        ),

        plot.subtitle = element_text(
            size = 7.7,
            colour = "grey25",
            hjust = 0
        ),

        plot.tag = element_text(
            face = "bold",
            size = 15
        ),

        plot.tag.position = c(
            0.005,
            0.995
        ),

        legend.title = element_text(
            size = 8
        ),

        legend.text = element_text(
            size = 7.8
        ),

        plot.margin = margin(
            8, 9, 8, 10
        )
    )


# ============================================================
# 11. PANEL E
# Unchanged
# ============================================================

tomtom_summary <- tomtom %>%
    group_by(
        Tissue
    ) %>%
    summarise(
        Tomtom_matches =
            n(),

        q_lt_0_05 =
            sum(
                q_value < 0.05,
                na.rm = TRUE
            ),

        Best_q =
            min(
                q_value,
                na.rm = TRUE
            ),

        .groups = "drop"
    ) %>%

    mutate(
        Best_minus_log10_q =
            -log10(
                pmax(
                    Best_q,
                    1e-300
                )
            )
    )


panelE <- ggplot(
    tomtom,
    aes(
        x = Tissue,
        y = minus_log10_q
    )
) +

    geom_hline(
        yintercept =
            -log10(
                0.05
            ),
        linetype = "dashed",
        linewidth = 0.5,
        colour = "grey45"
    ) +

    geom_jitter(
        width = 0.15,
        height = 0,
        size = 2.6,
        shape = 21,
        fill = "#56B4E9",
        colour = "black",
        stroke = 0.3
    ) +

    geom_text(
        data = tomtom_summary,

        aes(
            x = Tissue,
            y = Best_minus_log10_q,

            label = paste0(
                q_lt_0_05,
                "/",
                Tomtom_matches,
                " q<0.05"
            )
        ),

        inherit.aes = FALSE,
        vjust = -0.75,
        size = 2.65
    ) +

    scale_y_continuous(
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
            "Cross-run STREME motif reproducibility",

        subtitle = paste0(
            "Default versus no-holdout sensitivity runs; ",
            "supportive sequence similarity, not independent replication"
        ),

        x = NULL,

        y =
            "-log10(Tomtom q)"
    ) +

    theme_pub


# ============================================================
# 12. Assemble
# ============================================================

final_plot <-
    (
        panelA |
        panelB
    ) /

    (
        panelC |
        panelD
    ) /

    panelE +

    plot_layout(
        heights = c(
            1.0,
            1.12,
            0.82
        )
    )


# ============================================================
# 13. Export FINAL v2
# ============================================================

PDF_OUT <- file.path(
    OUTDIR,
    "Supplementary_Figure_S4_FINAL_v3.pdf"
)

TIFF_OUT <- file.path(
    OUTDIR,
    "Supplementary_Figure_S4_FINAL_v3_600dpi.tiff"
)


ggsave(
    filename = PDF_OUT,
    plot = final_plot,
    device = cairo_pdf,
    width = 11.4,
    height = 11.2,
    units = "in"
)

ggsave(
    filename = TIFF_OUT,
    plot = final_plot,
    device = "tiff",
    width = 11.4,
    height = 11.2,
    units = "in",
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 14. Export exact v2 plotting tables
# ============================================================

write_tsv(
    bg_long,
    file.path(
        OUTDIR,
        "FigureS4A_FINAL_v2.tsv"
    )
)

write_tsv(
    class_qc,
    file.path(
        OUTDIR,
        "FigureS4B_FINAL_v2.tsv"
    )
)

write_tsv(
    ame_best,
    file.path(
        OUTDIR,
        "FigureS4C_FINAL_v2.tsv"
    )
)

write_tsv(
    streme_grid,
    file.path(
        OUTDIR,
        "FigureS4D_FINAL_v2_with_NR.tsv"
    )
)

write_tsv(
    tomtom,
    file.path(
        OUTDIR,
        "FigureS4E_FINAL_v2.tsv"
    )
)


# ============================================================
# 15. v2 change log
# ============================================================

change_log <- tibble(
    Panel = c(
        "S4A",
        "S4B",
        "S4C",
        "S4D",
        "S4E"
    ),

    FINAL_v2_change = c(
        "No change",
        paste0(
            "Removed diverging colour scale; ",
            "all exact-match cells shown uniformly"
        ),
        paste0(
            "Added explicit annotation that adj.P=0.05 ",
            "is descriptive and not the formal AME criterion"
        ),
        paste0(
            "Added six explicit NR cells for motif ranks ",
            "not reported by the corresponding default STREME run"
        ),
        "No change"
    ),

    Biological_values_changed =
        "NO"
)


write_tsv(
    change_log,
    file.path(
        OUTDIR,
        "FigureS4_FINAL_v2_change_log.tsv"
    )
)


# ============================================================
# 16. Panel-D QC
# ============================================================

panelD_qc <- tibble(
    Metric = c(
        "Display_grid_cells",
        "Reported_motif_cells",
        "NR_cells",
        "Holdout_E_cells",
        "Training_score_cells",
        "Formal_E_lt_0.05_cells",
        "Reported_motifs_changed"
    ),

    Value = c(
        nrow(streme_grid),

        sum(
            streme_grid$Report_status !=
                "Not reported"
        ),

        sum(
            streme_grid$Report_status ==
                "Not reported"
        ),

        sum(
            streme_grid$Report_status ==
                "Holdout E-value"
        ),

        sum(
            streme_grid$Report_status ==
                "Training score only"
        ),

        sum(
            streme_grid$Formal_E_lt_0_05,
            na.rm = TRUE
        ),

        0
    ),

    Expected = c(
        20,
        14,
        6,
        9,
        5,
        0,
        0
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
    panelD_qc,
    file.path(
        OUTDIR,
        "FigureS4D_FINAL_v2_NR_QC.tsv"
    )
)

if (any(panelD_qc$Status != "PASS")) {
    print(panelD_qc)
    stop(
        "Panel-D FINAL-v2 NR QC failed."
    )
}


# ============================================================
# 17. Updated caption
# ============================================================

caption <- paste0(
    "Supplementary Figure S4. Quality control and evidence hierarchy ",
    "for motif analyses of strongly coupled tissue-specific peak-gene ",
    "candidates. ",
    "(A) Standardized mean differences in GC content and log2 peak ",
    "length between foreground peaks and their one-to-one same-tissue ",
    "matched high-confidence non-strong background peaks. Dashed lines ",
    "indicate an absolute standardized mean difference of 0.10. ",
    "(B) Exact pair-specific genomic-class matching between foreground ",
    "and matched-background sets. All foreground-background count ",
    "differences were zero. ",
    "(C) Best-ranked AME motif for each tissue from the complete JASPAR ",
    "scan. Bars show within-motif adjusted P values for descriptive ",
    "candidate ranking. The dashed adjusted-P threshold is descriptive ",
    "and was not the formal significance criterion; none of the motifs ",
    "met the predefined formal AME criterion of E<0.05. ",
    "(D) Statistical-evidence status of the 14 motifs reported by the ",
    "default STREME runs. Nine motifs from Liver, Muscle and Spleen had ",
    "holdout-based E values and none met E<0.05. The five Cerebellum ",
    "motifs were reported with training-set scores (S) rather than ",
    "holdout E values and therefore were not interpreted as formally ",
    "significant. NR indicates a motif rank that was not reported by ",
    "the corresponding default STREME run. ",
    "(E) Tomtom sequence-similarity comparisons between default STREME ",
    "motifs and motifs recovered in the no-holdout sensitivity runs. ",
    "This comparison represents supportive sequence-level reproducibility ",
    "rather than independent replication, transcription-factor occupancy, ",
    "or causal regulatory validation. Motifs carried forward to downstream ",
    "analyses were therefore treated as prioritized exploratory sequence ",
    "features."
)


writeLines(
    caption,
    file.path(
        OUTDIR,
        "Supplementary_Figure_S4_FINAL_v2_caption.txt"
    )
)


# ============================================================
# 18. Final file QC
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

        Status = if_else(
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
        "FigureS4_FINAL_v2_file_QC.tsv"
    )
)


# ============================================================
# 19. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10F3B4D COMPLETED\n")
cat("============================================================\n\n")

cat("Frozen biological QC:\n")
print(
    qc,
    width = Inf
)

cat("\nPanel-D NR QC:\n")
print(
    panelD_qc,
    width = Inf
)

cat("\nFINAL v2 change log:\n")
print(
    change_log,
    width = Inf
)

cat("\nOutput QC:\n")
print(
    file_qc
)

cat("\nIMPORTANT:\n")
cat("- Input consists only of frozen S4-FIX plotting tables.\n")
cat("- No AME analysis was rerun.\n")
cat("- No STREME analysis was rerun.\n")
cat("- No Tomtom analysis was rerun.\n")
cat("- Panel B changed visual encoding only.\n")
cat("- Panel C changed threshold labeling only.\n")
cat("- Panel D added NR placeholders only.\n")
cat("- 14 reported STREME motifs remain unchanged.\n")
cat("- 9 holdout-E motifs remain unchanged.\n")
cat("- 5 score-only motifs remain unchanged.\n")
cat("- Formal AME E<0.05 hits remain 0.\n")
cat("- Formal STREME holdout E<0.05 hits remain 0.\n")
cat("- Biological results modified: 0.\n")

if (
    all(qc$Status == "PASS") &&
    all(panelD_qc$Status == "PASS") &&
    all(file_qc$Status == "PASS")
) {

    cat("\nSTEP10F3B4D STATUS: PASS\n")
    cat(
        "Supplementary Figure S4 FINAL v3 ",
        "is ready for final visual lock.\n"
    )

} else {

    cat("\nSTEP10F3B4C STATUS: CHECK\n")
}

cat("============================================================\n")
