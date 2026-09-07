#!/usr/bin/env Rscript

# ============================================================
# STEP10C-v3 FINAL PATCH
#
# Only polish:
#   Figure3A
#   Figure4B
#
# All other panels are inherited from Step10C-v2.
# No statistical analysis is changed.
# ============================================================


# ============================================================
# 0. Rebuild frozen v2 plotting objects
# ============================================================

V2_SCRIPT <- "10C_v2_publication_polish_Figure3_Figure4.R"

if (!file.exists(V2_SCRIPT)) {
    stop(
        "Cannot find v2 script: ",
        V2_SCRIPT
    )
}

cat("============================================================\n")
cat("STEP10C-v3 FINAL PATCH\n")
cat("============================================================\n")
cat("Loading frozen Step10C-v2 plotting objects...\n")

source(
    V2_SCRIPT,
    local = FALSE
)


# ============================================================
# 1. New output directory
# ============================================================

OUTDIR_V3 <-
    "10_publication_figures/Step10C_Figure3_4_v3"

FIG3DIR_V3 <- file.path(
    OUTDIR_V3,
    "Figure3"
)

FIG4DIR_V3 <- file.path(
    OUTDIR_V3,
    "Figure4"
)

dir.create(
    FIG3DIR_V3,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    FIG4DIR_V3,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 2. FIGURE 3A v3
#
# Force labels onto multiple lines to avoid crowding.
# ============================================================

flow_nodes_v3 <- tibble(

    x = c(
        1,
        2,
        3
    ),

    N = c(
        168999,
        4674,
        702
    ),

    Stage = c(
        "All\nassigned",
        "Concordant\nHC",
        "Strong"
    ),

    Stage_key = c(
        "All assigned",
        "Concordant HC",
        "Strong"
    )
)


flow_arrows_v3 <- tibble(

    x = c(
        1.28,
        2.28
    ),

    xend = c(
        1.72,
        2.72
    ),

    label = c(
        "2.77%",
        "15.02%"
    )
)


p3A_v3 <- ggplot() +

    # Large counts
    geom_text(

        data = flow_nodes_v3,

        aes(
            x = x,
            y = 1.10,
            label = comma(N),
            color = Stage_key
        ),

        size = 4.0,
        fontface = "bold",
        lineheight = 0.95
    ) +

    # Stage names forced onto 1–2 lines
    geom_text(

        data = flow_nodes_v3,

        aes(
            x = x,
            y = 0.87,
            label = Stage
        ),

        size = 2.30,
        lineheight = 0.90
    ) +

    # Arrows
    geom_segment(

        data = flow_arrows_v3,

        aes(
            x = x,
            xend = xend,
            y = 1.025,
            yend = 1.025
        ),

        linewidth = 0.68,

        arrow =
            grid::arrow(
                length =
                    grid::unit(
                        0.085,
                        "inches"
                    ),
                type = "closed"
            )
    ) +

    # Selection percentages
    geom_text(

        data = flow_arrows_v3,

        aes(
            x =
                (
                    x +
                    xend
                ) / 2,

            y = 0.73,

            label = label
        ),

        size = 2.40
    ) +

    scale_color_manual(

        values = c(
            "All assigned" = "#737373",
            "Concordant HC" = "#377EB8",
            "Strong" = "#E41A1C"
        ),

        guide = "none"
    ) +

    coord_cartesian(

        xlim = c(
            0.53,
            3.47
        ),

        ylim = c(
            0.63,
            1.22
        ),

        clip = "off"
    ) +

    theme_void(
        base_family = "sans"
    )


# ============================================================
# 3. FIGURE 4B v3
#
# FDR moved to left labels.
# Right side displays OR only.
# ============================================================

promoter_tests_v3 <- tibble(

    Comparison = factor(

        c(
            "All non-Strong\nFDR = 6.0×10⁻⁹",
            "HC non-Strong\nFDR = 1.0×10⁻⁵"
        ),

        levels = rev(
            c(
                "All non-Strong\nFDR = 6.0×10⁻⁹",
                "HC non-Strong\nFDR = 1.0×10⁻⁵"
            )
        )
    ),

    OR = c(
        1.9088698783783362,
        1.7516790765064936
    ),

    CI_low = c(
        1.55501317924468,
        1.3930242446456977
    ),

    CI_high = c(
        2.3271189731953967,
        2.1918987049824867
    )
)


p4B_v3 <- ggplot(
    promoter_tests_v3,
    aes(
        y = Comparison
    )
) +

    geom_vline(

        xintercept = 1,

        linetype = "dashed",

        linewidth = 0.42
    ) +

    # Horizontal 95% CI
    geom_errorbar(

        aes(
            xmin = CI_low,
            xmax = CI_high,
            x = OR
        ),

        orientation = "y",

        width = 0.14,

        linewidth = 0.62
    ) +

    geom_point(

        aes(
            x = OR
        ),

        size = 3.05,

        color = "#E64B35"
    ) +

    # OR only on right side
    geom_text(

        aes(
            x = CI_high + 0.08,
            label =
                sprintf(
                    "OR = %.2f",
                    OR
                )
        ),

        hjust = 0,

        size = 2.40
    ) +

    scale_x_continuous(

        limits = c(
            0.85,
            2.95
        ),

        breaks = c(
            1,
            1.5,
            2,
            2.5
        )
    ) +

    labs(
        x = "Promoter enrichment odds ratio",
        y = NULL
    ) +

    theme_pub() +

    theme(

        axis.text.y =
            element_text(
                size = 7.0,
                lineheight = 0.88
            ),

        plot.margin =
            margin(
                5,
                8,
                5,
                5
            )
    )


# ============================================================
# 4. Save revised individual panels
# ============================================================

save_panel(
    p3A_v3,
    file.path(
        FIG3DIR_V3,
        "Figure3A_flow_v3"
    ),
    width = 3.3,
    height = 2.2
)


save_panel(
    p4B_v3,
    file.path(
        FIG4DIR_V3,
        "Figure4B_promoter_forest_v3"
    ),
    width = 3.8,
    height = 2.6
)


# ============================================================
# 5. Reassemble Figure 3 v3
#
# B-F unchanged from v2.
# ============================================================

top3_v3 <- (

    p3A_v3 |
    p3B |
    p3C

) +

    plot_layout(
        widths = c(
            1.12,
            1.00,
            1.00
        )
    )


bottom3_v3 <- (

    p3D |
    p3E |
    p3F

) +

    plot_layout(
        widths = c(
            1.05,
            0.90,
            1.00
        )
    )


Figure3_v3 <- (

    top3_v3 /
    bottom3_v3

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
            1
        )
    )


# ============================================================
# 6. Reassemble Figure 4 v3
#
# A,C,D,E unchanged from v2.
# ============================================================

top4_v3 <- (

    p4A |
    p4B_v3

) +

    plot_layout(
        widths = c(
            1.00,
            1.22
        )
    )


bottom4_v3 <- (
    p4D |
    p4E
)


Figure4_v3 <- (

    top4_v3 /
    p4C /
    bottom4_v3

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
            1.00,
            0.70,
            1.00
        )
    )


# ============================================================
# 7. Save final Figure 3 v3
# ============================================================

ggsave(
    file.path(
        FIG3DIR_V3,
        "Figure3_FINAL_v3.pdf"
    ),
    Figure3_v3,
    width = 7.1,
    height = 6.0,
    device = cairo_pdf
)


ggsave(
    file.path(
        FIG3DIR_V3,
        "Figure3_FINAL_v3_600dpi.tiff"
    ),
    Figure3_v3,
    width = 7.1,
    height = 6.0,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 8. Save final Figure 4 v3
# ============================================================

ggsave(
    file.path(
        FIG4DIR_V3,
        "Figure4_FINAL_v3.pdf"
    ),
    Figure4_v3,
    width = 7.1,
    height = 7.0,
    device = cairo_pdf
)


ggsave(
    file.path(
        FIG4DIR_V3,
        "Figure4_FINAL_v3_600dpi.tiff"
    ),
    Figure4_v3,
    width = 7.1,
    height = 7.0,
    dpi = 600,
    compression = "lzw"
)


# ============================================================
# 9. Frozen-value QC
# ============================================================

qc <- tibble(

    Metric = c(
        "All_UROPA_pairs",
        "Concordant_HC_pairs",
        "Official_Strong_pairs",
        "HC_animal_agreement",
        "Promoter_OR_All_nonStrong",
        "Promoter_OR_HC_nonStrong",
        "All_TSS_median",
        "HC_TSS_median",
        "Strong_TSS_median",
        "Strong_unique_genes"
    ),

    Value = c(
        nrow(all_cor),
        nrow(hc_cor),
        nrow(strong),
        agreement_hc,
        1.9088698783783362,
        1.7516790765064936,
        tss_summary$Median_bp[
            tss_summary$Group ==
                "All assigned"
        ],
        tss_summary$Median_bp[
            tss_summary$Group ==
                "Concordant HC"
        ],
        tss_summary$Median_bp[
            tss_summary$Group ==
                "Strong"
        ],
        sum(
            mult_summary$Gene_N
        )
    ),

    Expected = c(
        168999,
        4674,
        702,
        0.3851,
        1.909,
        1.752,
        33378,
        17248.5,
        11549,
        450
    )
)


qc <- qc %>%

    mutate(

        Status =
            case_when(

                Metric ==
                    "HC_animal_agreement" &
                    abs(
                        Value -
                        Expected
                    ) <
                    0.001 ~
                    "PASS",

                grepl(
                    "Promoter_OR",
                    Metric
                ) &
                    abs(
                        Value -
                        Expected
                    ) <
                    0.01 ~
                    "PASS",

                !grepl(
                    "Promoter_OR",
                    Metric
                ) &
                    Metric !=
                    "HC_animal_agreement" &
                    abs(
                        Value -
                        Expected
                    ) <
                    0.01 ~
                    "PASS",

                TRUE ~
                    "FAIL"
            )
    )


write_tsv(
    qc,
    file.path(
        OUTDIR_V3,
        "10C_v3_QC_summary.tsv"
    )
)


if (
    any(
        qc$Status != "PASS"
    )
) {

    stop(
        "Step10C-v3 frozen-result QC failed."
    )
}


# ============================================================
# 10. Output file QC
# ============================================================

final_files <- c(

    file.path(
        FIG3DIR_V3,
        "Figure3_FINAL_v3.pdf"
    ),

    file.path(
        FIG3DIR_V3,
        "Figure3_FINAL_v3_600dpi.tiff"
    ),

    file.path(
        FIG4DIR_V3,
        "Figure4_FINAL_v3.pdf"
    ),

    file.path(
        FIG4DIR_V3,
        "Figure4_FINAL_v3_600dpi.tiff"
    )
)


file_qc <- tibble(

    File = final_files,

    Exists =
        file.exists(
            final_files
        ),

    Size_bytes =
        file.info(
            final_files
        )$size
)


write_tsv(
    file_qc,
    file.path(
        OUTDIR_V3,
        "10C_v3_output_QC.tsv"
    )
)


if (
    any(
        !file_qc$Exists
    )
) {

    stop(
        "Some Step10C-v3 final files were not created."
    )
}


# ============================================================
# 11. Notes
# ============================================================

notes <- c(

    "STEP10C-v3 FINAL PUBLICATION PATCH",

    "",

    "Figure3A:",
    "Counts retained: 168999 -> 4674 -> 702.",
    "Stage labels were forced onto separate lines.",
    "Selection percentages remain 2.77% and 15.02%.",

    "",

    "Figure4B:",
    "Strong is the case group.",
    "Reference groups are disjoint non-Strong sets.",
    "Strong vs All non-Strong: OR = 1.909, FDR = 6.04e-9.",
    "Strong vs HC non-Strong: OR = 1.752, FDR = 1.03e-5.",
    "FDR annotations were moved into the left comparison labels.",

    "",

    "No biological or statistical result was recalculated or redefined."
)


writeLines(
    notes,
    file.path(
        OUTDIR_V3,
        "10C_v3_patch_notes.txt"
    )
)


capture.output(
    sessionInfo(),
    file =
        file.path(
            OUTDIR_V3,
            "10C_v3_sessionInfo.txt"
        )
)


# ============================================================
# 12. Final console output
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10C-v3 COMPLETED\n")
cat("============================================================\n")

cat("\nPATCHED PANELS:\n")
cat("Figure3A : multiline stage labels\n")
cat("Figure4B : FDR moved to left labels; OR only on right\n")

cat("\nFROZEN RESULTS:\n")
cat("All assigned       : 168999\n")
cat("Concordant HC      : 4674\n")
cat("Strong             : 702\n")
cat("HC agreement       : 0.3851\n")
cat("Strong vs All OR   : 1.909\n")
cat("Strong vs HC OR    : 1.752\n")

cat("\nQC:\n")

print(
    qc,
    n = Inf,
    width = Inf
)

cat("\nFINAL FILES:\n")

print(
    file_qc,
    n = Inf,
    width = Inf
)

cat("\nNo statistical result was changed.\n")
cat("============================================================\n")

