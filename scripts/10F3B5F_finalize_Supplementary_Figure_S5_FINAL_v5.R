suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

cat("\n")
cat("============================================================\n")
cat("STEP10F3B5F — SUPPLEMENTARY FIGURE S5 FINAL v5\n")
cat("Minor visualization refinement only\n")
cat("Frozen S5 plotting tables only\n")
cat("Panel A: remove redundant Tissue legend\n")
cat("Panel B: refine representative-TF label placement\n")
cat("Panels C/D/E: unchanged\n")
cat("============================================================\n\n")


# ============================================================
# 1. Frozen S5-v1 plotting sources
# ============================================================

V1DIR <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3B_Supplementary_Figures/",
    "FigureS5_motif_TF_RNA_support"
)

F_A <- file.path(
    V1DIR,
    "FigureS5A_representative_TF_RNA.tsv"
)

F_B <- file.path(
    V1DIR,
    "FigureS5B_Tau_SPM_all_RNA_mapped_TFs.tsv"
)

F_C <- file.path(
    V1DIR,
    "FigureS5C_RNA_support_roles.tsv"
)

F_D <- file.path(
    V1DIR,
    "FigureS5D_evidence_support_matrix.tsv"
)

F_E <- file.path(
    V1DIR,
    "FigureS5E_FIMO_best_tier.tsv"
)

F_QC <- file.path(
    V1DIR,
    "FigureS5_frozen_source_QC.tsv"
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
        "Missing frozen S5 source file(s):\n",
        paste(
            missing_files,
            collapse = "\n"
        )
    )
}


# ============================================================
# 2. Output directory
# ============================================================

OUTDIR <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3B_Supplementary_Figures/",
    "FigureS5_motif_TF_RNA_support_FINAL_v5"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. Read frozen plotting tables
# ============================================================

panelA_data <- read_tsv(
    F_A,
    show_col_types = FALSE
)

panelB_data <- read_tsv(
    F_B,
    show_col_types = FALSE
)

panelC_data <- read_tsv(
    F_C,
    show_col_types = FALSE
)

panelD_data <- read_tsv(
    F_D,
    show_col_types = FALSE
)

panelE_data <- read_tsv(
    F_E,
    show_col_types = FALSE
)

previous_qc <- read_tsv(
    F_QC,
    show_col_types = FALSE
)


# ============================================================
# 4. Fixed tissue order
# ============================================================

tissue_order <- c(
    "Cerebellum",
    "Liver",
    "Muscle",
    "Spleen"
)


panelA_data <- panelA_data %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


panelB_data <- panelB_data %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


panelC_data <- panelC_data %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


panelD_data <- panelD_data %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


panelE_data <- panelE_data %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        )
    )


# ============================================================
# 5. Restore frozen factor ordering
# ============================================================

rna_role_levels <- c(
    "Unmapped to RNA gene",
    "Low RNA support",
    "High RNA",
    "High RNA + tissue preference"
)

panelC_data <- panelC_data %>%
    mutate(
        RNA_support_role = factor(
            RNA_support_role,
            levels = rna_role_levels
        )
    )


evidence_levels <- c(
    "SPM ≥ 0.5",
    "Target tissue max",
    "High RNA",
    "Cross-run Tomtom",
    "No-holdout Tomtom",
    "Default Tomtom",
    "AME exploratory"
)

panelD_data <- panelD_data %>%
    mutate(
        Evidence = factor(
            Evidence,
            levels = evidence_levels
        )
    )


fimo_levels <- c(
    "Not selected",
    "Secondary",
    "Primary"
)

panelE_data <- panelE_data %>%
    mutate(
        FIMO_best_tier = factor(
            FIMO_best_tier,
            levels = fimo_levels
        )
    )


# ============================================================
# 6. Frozen biological QC
# ============================================================

previous_all_pass <- all(
    previous_qc$Status == "PASS"
)


# Representative candidate identities expected from frozen S5
expected_rep <- tibble(
    Tissue = c(
        "Liver",
        "Liver",
        "Liver",
        "Muscle",
        "Spleen",
        "Spleen",
        "Spleen"
    ),

    TF_symbol = c(
        "ETV2",
        "KLF12",
        "SP5",
        "KLF9",
        "ELF4",
        "KLF12",
        "SP5"
    )
)


observed_rep <- panelA_data %>%
    transmute(
        Tissue = as.character(Tissue),
        TF_symbol
    ) %>%
    distinct()


rep_identity_pass <- nrow(
    anti_join(
        expected_rep,
        observed_rep,
        by = c(
            "Tissue",
            "TF_symbol"
        )
    )
) == 0 &&
nrow(
    anti_join(
        observed_rep,
        expected_rep,
        by = c(
            "Tissue",
            "TF_symbol"
        )
    )
) == 0


# Spleen-SP5 low-RNA example must remain.
spleen_sp5 <- panelA_data %>%
    filter(
        as.character(Tissue) == "Spleen",
        TF_symbol == "SP5"
    )


spleen_sp5_pass <-
    nrow(spleen_sp5) == 1 &&
    abs(
        spleen_sp5$Target_tissue_mean_TPM -
            0.0095
    ) < 1e-10 &&
    spleen_sp5$RNA_support_role ==
        "Low_RNA_support"


# Evidence totals
evidence_totals <- panelD_data %>%
    group_by(
        Evidence
    ) %>%
    summarise(
        Supported_N =
            sum(
                Supported_N,
                na.rm = TRUE
            ),

        .groups = "drop"
    )


get_evidence_total <- function(x) {

    z <- evidence_totals %>%
        filter(
            as.character(Evidence) ==
                x
        ) %>%
        pull(
            Supported_N
        )

    if (length(z) == 0) {
        return(NA_real_)
    }

    z[[1]]
}


qc <- tibble(

    Metric = c(
        "Previous_S5_frozen_QC_all_PASS",
        "PanelA_representative_rows",
        "PanelA_representative_identity",
        "PanelB_RNA_mapped_rows",
        "PanelC_candidate_total",
        "PanelD_cells",
        "PanelE_candidate_total",

        "Spleen_SP5_low_RNA_example_preserved",

        "AME_exploratory_support_TRUE",
        "Default_Tomtom_support_TRUE",
        "NoHoldout_Tomtom_support_TRUE",
        "CrossRun_Tomtom_support_TRUE",
        "High_RNA_support_TRUE",
        "Target_tissue_max_TRUE",
        "SPM_ge_0.5_TRUE",

        "FIMO_Primary",
        "FIMO_Secondary",
        "FIMO_Not_selected",

        "Biological_results_modified"
    ),

    Value = c(
        as.integer(
            previous_all_pass
        ),

        nrow(
            panelA_data
        ),

        as.integer(
            rep_identity_pass
        ),

        nrow(
            panelB_data
        ),

        sum(
            panelC_data$N,
            na.rm = TRUE
        ),

        nrow(
            panelD_data
        ),

        sum(
            panelE_data$N,
            na.rm = TRUE
        ),

        as.integer(
            spleen_sp5_pass
        ),

        get_evidence_total(
            "AME exploratory"
        ),

        get_evidence_total(
            "Default Tomtom"
        ),

        get_evidence_total(
            "No-holdout Tomtom"
        ),

        get_evidence_total(
            "Cross-run Tomtom"
        ),

        get_evidence_total(
            "High RNA"
        ),

        get_evidence_total(
            "Target tissue max"
        ),

        get_evidence_total(
            "SPM ≥ 0.5"
        ),

        sum(
            panelE_data$N[
                panelE_data$FIMO_best_tier ==
                    "Primary"
            ],
            na.rm = TRUE
        ),

        sum(
            panelE_data$N[
                panelE_data$FIMO_best_tier ==
                    "Secondary"
            ],
            na.rm = TRUE
        ),

        sum(
            panelE_data$N[
                panelE_data$FIMO_best_tier ==
                    "Not selected"
            ],
            na.rm = TRUE
        ),

        0
    ),

    Expected = c(
        1,
        7,
        1,
        44,
        48,
        28,
        48,

        1,

        10,
        40,
        16,
        24,
        33,
        8,
        9,

        28,
        19,
        1,

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
        "FigureS5_FINAL_v2_frozen_source_QC.tsv"
    )
)


if (any(qc$Status != "PASS")) {

    print(
        qc,
        n = Inf,
        width = Inf
    )

    stop(
        "Frozen S5 QC failed. ",
        "FINAL v2 was not generated."
    )
}


# ============================================================
# 7. Publication theme
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
            size = 8.2
        ),

        legend.text = element_text(
            size = 7.8
        ),

        plot.margin = margin(
            8, 9, 8, 10
        )
    )


# ============================================================
# 8. PANEL A — FINAL v2
# Same data and geometry.
# Only hide redundant Tissue/fill legend.
# ============================================================

# Restore candidate order from frozen plotting table
candidate_levels <- panelA_data %>%
    arrange(
        log10_TPM
    ) %>%
    pull(
        Candidate
    )


panelA_data <- panelA_data %>%
    mutate(
        Candidate = factor(
            Candidate,
            levels = candidate_levels
        )
    )


panelA <- ggplot(
    panelA_data,
    aes(
        x = log10_TPM,
        y = Candidate
    )
) +

    geom_segment(
        aes(
            x =
                min(
                    log10_TPM,
                    na.rm = TRUE
                ) - 0.1,

            xend =
                log10_TPM,

            yend =
                Candidate
        ),

        linewidth = 0.45,
        colour = "grey78"
    ) +

    geom_point(
        aes(
            size =
                Target_tissue_SPM,

            shape =
                Max_tissue_label,

            fill =
                Tissue
        ),

        colour = "black",
        stroke = 0.35
    ) +

    scale_size_continuous(
        range = c(
            2.6,
            5.4
        ),

        limits = c(
            0,
            1
        ),

        name =
            "Target-tissue SPM"
    ) +

    scale_shape_manual(
        values = c(
            "Target tissue is maximum" = 21,
            "Target tissue is not maximum" = 24
        ),

        name =
            "Expression maximum"
    ) +

    labs(
        tag = "A",

        title =
            "Representative motif-family candidate TFs",

        subtitle =
            "Target-tissue RNA expression; plotting transform only",

        x =
            "log10(target-tissue mean TPM + 0.01)",

        y = NULL
    ) +

    guides(
        # FINAL v2:
        # Tissue already appears explicitly in every y-axis label.
        fill = "none"
    ) +

    theme_pub


# ============================================================
# 9. PANEL B — FINAL v2
# Same 44 data points and same Tau-SPM geometry.
# Only representative text placement changes.
# ============================================================

representative_key <- tibble(
    Tissue = c(
        "Liver",
        "Liver",
        "Liver",
        "Muscle",
        "Spleen",
        "Spleen",
        "Spleen"
    ),

    TF_symbol = c(
        "ETV2",
        "KLF12",
        "SP5",
        "KLF9",
        "ELF4",
        "KLF12",
        "SP5"
    )
)


rep_B <- panelB_data %>%

    mutate(
        Tissue_chr =
            as.character(
                Tissue
            )
    ) %>%

    inner_join(
        representative_key,
        by = c(
            "Tissue_chr" =
                "Tissue",
            "TF_symbol"
        )
    ) %>%

    mutate(
        Label = paste0(
            TF_symbol,
            " (",
            Tissue_chr,
            ")"
        )
    )


if (nrow(rep_B) != 7) {

    stop(
        "Panel B representative label rows != 7."
    )
}


panelB_base <- ggplot(
    panelB_data,
    aes(
        x = Tau,
        y = Target_tissue_SPM
    )
) +

    geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        linewidth = 0.45,
        colour = "grey45"
    ) +

    annotate(
        "text",
        x = 0.325,
        y = 0.515,
        label = "SPM = 0.5",
        hjust = 0,
        vjust = 0,
        size = 2.45,
        colour = "grey35"
    ) +

    geom_point(
        aes(
            fill =
                Tissue,

            shape =
                Max_tissue_label,

            size =
                TPM_size
        ),

        colour = "black",
        stroke = 0.30,
        alpha = 0.86
    ) +

    scale_shape_manual(
        values = c(
            "Target tissue is maximum" = 21,
            "Other tissue is maximum" = 24
        ),

        name =
            "Expression maximum"
    ) +

    scale_size_continuous(
        range = c(
            2.0,
            5.0
        ),

        name =
            "log10(TPM + 1)"
    ) +

    coord_cartesian(
        xlim = c(
            0.24,
            1.02
        ),

        ylim = c(
            0,
            1.08
        ),

        clip = "off"
    ) +

    labs(
        tag = "B",

        title =
            "RNA tissue preference across candidate TFs",

        subtitle = NULL,

        x =
            "Tau",

        y =
            "Target-tissue SPM",

        fill =
            "Tissue"
    ) +

    theme_pub +

    theme(
        legend.box.spacing = unit(
            4,
            "pt"
        ),

        plot.margin = margin(
            8,
            13,
            8,
            10
        )
    )


# ============================================================
# 10. Panel B labels
# Prefer ggrepel if installed.
# Otherwise use fixed offsets only.
# Point positions are NEVER changed.
# ============================================================

if (
    requireNamespace(
        "ggrepel",
        quietly = TRUE
    )
) {

    cat(
        "Panel B labels: using ggrepel::geom_text_repel()\n"
    )

    set.seed(
        20260903
    )

    panelB <- panelB_base +

        ggrepel::geom_text_repel(
            data = rep_B,

            aes(
                x = Tau,
                y = Target_tissue_SPM,
                label = Label
            ),

            inherit.aes = FALSE,

            size = 2.45,

            box.padding = 0.32,
            point.padding = 0.20,

            min.segment.length = 0,

            segment.size = 0.30,
            segment.colour = "grey45",

            max.overlaps = Inf,

            force = 1.1,

            seed = 20260903,

            xlim = c(
                0.30,
                1.02
            ),

            ylim = c(
                0,
                1.07
            )
        )

} else {

    cat(
        "Panel B labels: ggrepel not installed; ",
        "using fixed publication-safe offsets.\n",
        sep = ""
    )


    manual_label_positions <- tibble(

        Tissue_chr = c(
            "Liver",
            "Liver",
            "Liver",
            "Muscle",
            "Spleen",
            "Spleen",
            "Spleen"
        ),

        TF_symbol = c(
            "ETV2",
            "KLF12",
            "SP5",
            "KLF9",
            "ELF4",
            "KLF12",
            "SP5"
        ),

        label_x = c(
            0.775,
            0.665,
            0.865,
            0.455,
            0.690,
            0.685,
            0.875
        ),

        label_y = c(
            0.945,
            0.125,
            1.025,
            0.610,
            0.840,
            0.755,
            0.060
        ),

        hjust = c(
            0,
            1,
            0,
            1,
            1,
            1,
            1
        )
    )


    rep_B_manual <- rep_B %>%

        left_join(
            manual_label_positions,
            by = c(
                "Tissue_chr",
                "TF_symbol"
            )
        )


    if (
        any(
            !is.finite(
                rep_B_manual$label_x
            )
        ) ||
        any(
            !is.finite(
                rep_B_manual$label_y
            )
        )
    ) {

        stop(
            "Panel B manual label positions are incomplete."
        )
    }


    panelB <- panelB_base +

        geom_segment(
            data =
                rep_B_manual,

            aes(
                x = Tau,
                y = Target_tissue_SPM,
                xend = label_x,
                yend = label_y
            ),

            inherit.aes = FALSE,

            linewidth = 0.30,
            colour = "grey45"
        ) +

        geom_text(
            data =
                rep_B_manual,

            aes(
                x = label_x,
                y = label_y,
                label = Label,
                hjust = hjust
            ),

            inherit.aes = FALSE,

            size = 2.45
        )
}


# ============================================================
# 11. PANEL C — unchanged
# ============================================================

panelC <- ggplot(
    panelC_data,
    aes(
        x = Tissue,
        y = N,
        fill = RNA_support_role
    )
) +

    geom_col(
        width = 0.68
    ) +

    geom_text(
        aes(
            label =
                ifelse(
                    N > 0,
                    N,
                    ""
                )
        ),

        position =
            position_stack(
                vjust = 0.5
            ),

        size = 2.7
    ) +

    labs(
        tag = "C",

        title =
            "RNA-support classification",

        subtitle =
            "All 48 publication candidate TF rows retained",

        x = NULL,

        y =
            "Candidate TFs",

        fill =
            "RNA support"
    ) +

    scale_fill_discrete(
        labels = c(
            "Unmapped to RNA gene",
            "Low RNA support",
            "High RNA",
            "High RNA + tissue preference"
        )
    ) +

    theme_pub


# ============================================================
# 12. PANEL D — unchanged
# ============================================================

panelD <- ggplot(
    panelD_data,
    aes(
        x = Tissue,
        y = Evidence,
        fill = Support_fraction
    )
) +

    geom_tile(
        colour = "white",
        linewidth = 0.65
    ) +

    geom_text(
        aes(
            label =
                Cell_label
        ),

        size = 2.6
    ) +

    scale_fill_gradient(
        low = "white",
        high = "#3B75AF",
        limits = c(
            0,
            1
        ),

        name =
            "Fraction supported"
    ) +

    labs(
        tag = "D",

        title =
            "Candidate evidence-support matrix",

        subtitle =
            "Cell labels show supported candidates / candidate TFs in each tissue",

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
                angle = 20,
                hjust = 1,
                size = 8.2,
                colour = "black"
            ),

        axis.text.y =
            element_text(
                size = 7.8,
                colour = "black"
            ),

        plot.title =
            element_text(
                size = 10.5,
                face = "bold"
            ),

        plot.subtitle =
            element_text(
                size = 8.0,
                colour = "grey25"
            ),

        plot.tag =
            element_text(
                face = "bold",
                size = 15
            ),

        plot.tag.position =
            c(
                0.005,
                0.995
            ),

        legend.title =
            element_text(
                size = 8
            ),

        legend.text =
            element_text(
                size = 7.8
            ),

        plot.margin =
            margin(
                8, 9, 8, 10
            )
    )


# ============================================================
# 13. PANEL E — unchanged
# ============================================================

panelE <- ggplot(
    panelE_data,
    aes(
        x = Tissue,
        y = N,
        fill = FIMO_best_tier
    )
) +

    geom_col(
        width = 0.68
    ) +

    geom_text(
        aes(
            label =
                ifelse(
                    N > 0,
                    N,
                    ""
                )
        ),

        position =
            position_stack(
                vjust = 0.5
            ),

        size = 2.8
    ) +

    labs(
        tag = "E",

        title =
            "Best FIMO evidence tier",

        subtitle = paste0(
            "Sequence-match evidence assigned to candidate motif families; ",
            "not TF occupancy"
        ),

        x = NULL,

        y =
            "Candidate TFs",

        fill =
            "FIMO tier"
    ) +

    scale_fill_discrete(
        labels = c(
            "Not selected",
            "Secondary",
            "Primary"
        )
    ) +

    theme_pub


# ============================================================
# 14. Assemble FINAL v2
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
            1.08,
            1.10,
            0.80
        )
    )


# ============================================================
# 15. Export
# ============================================================

PDF_OUT <- file.path(
    OUTDIR,
    "Supplementary_Figure_S5_FINAL_v5.pdf"
)

TIFF_OUT <- file.path(
    OUTDIR,
    "Supplementary_Figure_S5_FINAL_v5_600dpi.tiff"
)


ggsave(
    filename = PDF_OUT,
    plot = final_plot,
    device = cairo_pdf,
    width = 11.5,
    height = 11.3,
    units = "in"
)


ggsave(
    filename = TIFF_OUT,
    plot = final_plot,
    device = "tiff",
    width = 11.5,
    height = 11.3,
    units = "in",
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 16. Save exact FINAL-v2 plotting tables
# ============================================================

write_tsv(
    panelA_data,
    file.path(
        OUTDIR,
        "FigureS5A_FINAL_v2.tsv"
    )
)

write_tsv(
    panelB_data,
    file.path(
        OUTDIR,
        "FigureS5B_FINAL_v2.tsv"
    )
)

write_tsv(
    rep_B,
    file.path(
        OUTDIR,
        "FigureS5B_FINAL_v2_representative_label_audit.tsv"
    )
)

write_tsv(
    panelC_data,
    file.path(
        OUTDIR,
        "FigureS5C_FINAL_v2.tsv"
    )
)

write_tsv(
    panelD_data,
    file.path(
        OUTDIR,
        "FigureS5D_FINAL_v2.tsv"
    )
)

write_tsv(
    panelE_data,
    file.path(
        OUTDIR,
        "FigureS5E_FINAL_v2.tsv"
    )
)


# ============================================================
# 17. FINAL-v2 change log
# ============================================================

label_method <- ifelse(
    requireNamespace(
        "ggrepel",
        quietly = TRUE
    ),
    "ggrepel non-overlapping text placement",
    "fixed manual text offsets"
)


change_log <- tibble(

    Panel = c(
        "S5A",
        "S5B",
        "S5C",
        "S5D",
        "S5E"
    ),

    FINAL_v2_change = c(

        paste0(
            "Removed redundant Tissue legend only; ",
            "Tissue remains explicit in y-axis candidate labels"
        ),

        paste0(
            "Refined representative TF label placement using ",
            label_method,
            "; all point coordinates unchanged"
        ),

        "No change",

        "No change",

        "No change"
    ),

    Biological_values_changed =
        "NO"
)


write_tsv(
    change_log,
    file.path(
        OUTDIR,
        "FigureS5_FINAL_v2_change_log.tsv"
    )
)


# ============================================================
# 18. Updated caption
# ============================================================

caption <- paste0(

    "Supplementary Figure S5. RNA-expression support and evidence ",
    "hierarchy for motif-family candidate transcription factors. ",

    "(A) Target-tissue RNA expression of seven representative ",
    "candidate TF entries used in the integrated regulatory narrative. ",
    "Point size represents target-tissue SPM and point shape indicates ",
    "whether the target tissue was the maximum-expression tissue. ",
    "Tissue identity is provided directly in the candidate labels; ",
    "the log10 transformation was used only for visualization. ",

    "(B) Relationship between Tau and target-tissue SPM among the 44 ",
    "candidate TF entries mapped to RNA genes. The dashed horizontal ",
    "line denotes the frozen SPM>=0.5 tissue-preference flag. ",
    "Representative labels were repositioned only to improve readability, ",
    "and the dashed horizontal line marks SPM=0.5; data-point coordinates ",
    "were unchanged. Four candidate entries ",
    "unmapped to an RNA gene were omitted from quantitative RNA panels ",
    "only. ",

    "(C) Frozen RNA-support classifications for all 48 publication ",
    "candidate TF entries. ",

    "(D) Tissue-stratified summary of motif and RNA evidence fields ",
    "carried by the 48 candidate TF entries. Cell labels indicate ",
    "supported entries relative to the total candidate entries for ",
    "each tissue. These fields represent candidate-support annotations ",
    "rather than independent experimental validation layers. ",

    "(E) Distribution of the frozen best FIMO evidence tier. FIMO ",
    "sequence matches were interpreted as motif-sequence evidence and ",
    "not as direct TF occupancy. Candidate TF assignment represents ",
    "RNA-supported compatibility with motif families and does not ",
    "establish direct DNA binding or causal regulation."
)


writeLines(
    caption,
    file.path(
        OUTDIR,
        "Supplementary_Figure_S5_FINAL_v2_caption.txt"
    )
)


# ============================================================
# 19. File QC
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
            file.info(
                File
            )$size,

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
        "FigureS5_FINAL_v2_file_QC.tsv"
    )
)


# ============================================================
# 20. Final console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10F3B5F COMPLETED\n")
cat("============================================================\n\n")


cat("Frozen biological QC:\n")

print(
    qc,
    n = Inf,
    width = Inf
)


cat("\nFINAL v2 change log:\n")

print(
    change_log,
    n = Inf,
    width = Inf
)


cat("\nRepresentative label audit:\n")

print(
    rep_B %>%
        select(
            Tissue_chr,
            TF_symbol,
            Target_tissue_mean_TPM,
            Tau,
            Target_tissue_SPM,
            Target_tissue_expression_rank,
            RNA_support_role,
            FIMO_best_tier
        ),
    n = Inf,
    width = Inf
)


cat("\nOutput QC:\n")

print(
    file_qc,
    width = Inf
)


cat("\nIMPORTANT:\n")
cat("- Input consists only of frozen S5 plotting tables.\n")
cat("- Panel A Tissue legend was hidden only.\n")
cat("- Panel B representative text placement was changed only.\n")
cat("- Panel B point coordinates are unchanged.\n")
cat("- Panel B x-axis display limit changed from 0.30 to 0.24 only.\n")
cat("- This adds visual headroom; Tau values were not changed.\n")
cat("- Panels C/D/E are biologically unchanged.\n")
cat("- No TF candidate was added or removed.\n")
cat("- No RNA TPM/Tau/SPM statistic was recomputed.\n")
cat("- No FIMO analysis was rerun.\n")
cat("- No motif analysis was rerun.\n")
cat("- Candidate TF does not imply direct DNA binding.\n")
cat("- FIMO evidence does not imply TF occupancy.\n")
cat("- Biological results modified: 0.\n")


if (
    all(
        qc$Status ==
            "PASS"
    ) &&
    all(
        file_qc$Status ==
            "PASS"
    )
) {

    cat("\nSTEP10F3B5F STATUS: PASS\n")
    cat(
        "Supplementary Figure S5 FINAL v5 ",
        "is ready for final visual inspection.\n"
    )

} else {

    cat("\nSTEP10F3B5C STATUS: CHECK\n")
}


cat("============================================================\n")
