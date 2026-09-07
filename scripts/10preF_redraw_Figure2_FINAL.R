#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(scales)
    library(patchwork)
})

cat("\n============================================================\n")
cat("STEP10-preF — FIGURE 2 FINAL REDRAW\n")
cat("Layout-only; frozen biological results unchanged\n")
cat("============================================================\n\n")

outdir <- "10_publication_figures/Step10preF_Figure2_FINAL"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. FROZEN RESULTS
# ============================================================

lrt <- tibble(
    Assay  = factor(c("RNA", "ATAC"), levels = c("RNA", "ATAC")),
    Sig    = c(18870, 176917),
    Tested = c(22240, 203210)
) %>%
    mutate(
        Percent = 100 * Sig / Tested,
        Label = paste0(
            comma(Sig),
            " / ",
            comma(Tested),
            "\n",
            sprintf("%.1f%%", Percent)
        )
    )

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

rna_hc <- tibble(
    Tissue = tissue_order,
    Count = c(
        408,
        573,
        397,
        244,
        1182,
        610,
        733,
        889
    )
)

atac_hc <- tibble(
    Tissue = tissue_order,
    Count = c(
        47,
        4643,
        2689,
        785,
        5071,
        725,
        20853,
        5732
    )
)

# ------------------------------------------------------------
# Fixed tissue palette
# ------------------------------------------------------------

tissue_colors <- c(
    Adipose      = "#D89000",
    Cerebellum   = "#7CAE00",
    Cortex       = "#00BA38",
    Hypothalamus = "#00C08B",
    Liver        = "#00BFC4",
    Lung         = "#619CFF",
    Muscle       = "#C77CFF",
    Spleen       = "#F564E3"
)

# ============================================================
# 2. ROBUSTLY FIND FROZEN STEP03 TAU/SPM TABLES
# ============================================================

canonical <- function(x) {
    tolower(gsub("[^a-zA-Z0-9]", "", x))
}

has_tau_spm <- function(f) {

    z <- tryCatch(
        suppressMessages(
            read_tsv(
                f,
                n_max = 0,
                show_col_types = FALSE,
                progress = FALSE
            )
        ),
        error = function(e) NULL
    )

    if (is.null(z)) {
        return(FALSE)
    }

    nm <- canonical(names(z))

    has_tau <- any(grepl("tau$", nm))
    has_spm <- any(grepl("maxspm$", nm))

    has_tau && has_spm
}

count_file_lines <- function(f) {

    z <- tryCatch(
        system(
            paste(
                "wc -l",
                shQuote(f)
            ),
            intern = TRUE
        ),
        error = function(e) NA_character_
    )

    if (
        length(z) == 0 ||
        is.na(z[1])
    ) {
        return(NA_real_)
    }

    suppressWarnings(
        as.numeric(
            strsplit(
                trimws(z[1]),
                "[[:space:]]+"
            )[[1]][1]
        )
    )
}

find_specificity_table <- function(
    expected_rows,
    assay
) {

    files <- list.files(
        ".",
        pattern = "\\.tsv$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Prefer Step03 / specificity files first.
    preferred <- files[
        grepl(
            "03|specific|tau|spm",
            files,
            ignore.case = TRUE
        )
    ]

    if (length(preferred) == 0) {
        preferred <- files
    }

    header_ok <- vapply(
        preferred,
        has_tau_spm,
        logical(1)
    )

    candidates <- preferred[header_ok]

    if (length(candidates) == 0) {
        stop(
            "No Tau/MaxSPM candidate TSV files were found."
        )
    }

    nlines <- vapply(
        candidates,
        count_file_lines,
        numeric(1)
    )

    exact <- candidates[
        !is.na(nlines) &
        nlines == expected_rows + 1
    ]

    if (length(exact) == 0) {

        cat("\nCandidate Tau/SPM tables:\n")

        print(
            tibble(
                file = candidates,
                lines = nlines
            )
        )

        stop(
            paste0(
                "Could not find ",
                assay,
                " specificity table with exactly ",
                expected_rows,
                " data rows."
            )
        )
    }

    # Prefer assay name if multiple exact tables exist.
    assay_match <- exact[
        grepl(
            assay,
            exact,
            ignore.case = TRUE
        )
    ]

    if (length(assay_match) > 0) {
        exact <- assay_match
    }

    # Prefer tissue-specificity/Tau wording.
    score <- integer(length(exact))

    score <- score +
        3L * grepl(
            "tissue.*specific|specificity",
            exact,
            ignore.case = TRUE
        )

    score <- score +
        2L * grepl(
            "tau",
            exact,
            ignore.case = TRUE
        )

    exact <- exact[
        order(
            score,
            decreasing = TRUE
        )
    ]

    chosen <- exact[1]

    dat <- suppressMessages(
        read_tsv(
            chosen,
            show_col_types = FALSE,
            progress = FALSE
        )
    )

    if (nrow(dat) != expected_rows) {
        stop(
            "Unexpected feature count after reading: ",
            chosen
        )
    }

    list(
        file = chosen,
        data = dat
    )
}

rna_obj <- find_specificity_table(
    expected_rows = 35682,
    assay = "RNA"
)

atac_obj <- find_specificity_table(
    expected_rows = 203226,
    assay = "ATAC"
)

cat("Frozen RNA specificity source : ",
    rna_obj$file, "\n", sep = "")

cat("Frozen ATAC specificity source: ",
    atac_obj$file, "\n\n", sep = "")

# ============================================================
# 3. EXTRACT TAU / MAX SPM
# ============================================================

get_tau_spm <- function(dat) {

    nm <- names(dat)
    cn <- canonical(nm)

    # Prefer exact Tau if available.
    tau_candidates <- which(
        cn == "tau"
    )

    if (length(tau_candidates) == 0) {
        tau_candidates <- which(
            grepl("tau$", cn)
        )
    }

    if (length(tau_candidates) == 0) {
        stop("Tau column not found.")
    }

    tau_candidates <- tau_candidates[
        order(
            nchar(cn[tau_candidates])
        )
    ]

    tau_col <- nm[
        tau_candidates[1]
    ]

    # Prefer exact MaxSPM.
    spm_candidates <- which(
        cn == "maxspm"
    )

    if (length(spm_candidates) == 0) {
        spm_candidates <- which(
            grepl("maxspm$", cn)
        )
    }

    if (length(spm_candidates) == 0) {
        stop("Max_SPM column not found.")
    }

    spm_candidates <- spm_candidates[
        order(
            nchar(cn[spm_candidates])
        )
    ]

    spm_col <- nm[
        spm_candidates[1]
    ]

    tibble(
        Tau = suppressWarnings(
            as.numeric(dat[[tau_col]])
        ),
        MaxSPM = suppressWarnings(
            as.numeric(dat[[spm_col]])
        )
    )
}

rna_spec <- get_tau_spm(
    rna_obj$data
)

atac_spec <- get_tau_spm(
    atac_obj$data
)

rna_plot <- rna_spec %>%
    filter(
        is.finite(Tau),
        is.finite(MaxSPM)
    )

atac_plot <- atac_spec %>%
    filter(
        is.finite(Tau),
        is.finite(MaxSPM)
    )

# ============================================================
# 4. COMMON THEME
# ============================================================

base_theme <- theme_classic(
    base_size = 10
) +
    theme(
        axis.title = element_text(
            size = 10
        ),
        axis.text = element_text(
            size = 8.5,
            color = "black"
        ),
        plot.title = element_text(
            size = 10,
            face = "bold",
            hjust = 0.5,
            margin = margin(
                b = 7
            )
        ),
        plot.tag = element_text(
            size = 14,
            face = "bold"
        ),
        plot.tag.position = c(
            0.01,
            0.99
        ),
        plot.margin = margin(
            t = 12,
            r = 12,
            b = 8,
            l = 8
        )
    )

# ============================================================
# 5. PANEL A
#
# KEY FIXES:
# - two-line title
# - wider x expansion
# - labels split over 2 lines
# - extra top space
# ============================================================

pA <- ggplot(
    lrt,
    aes(
        x = Assay,
        y = Percent,
        fill = Assay
    )
) +
    geom_col(
        width = 0.58
    ) +
    geom_text(
        aes(
            label = Label
        ),
        y = c(
            91.0,
            94.0
        ),
        size = 3.0,
        lineheight = 0.90,
        color = "black"
    ) +
    scale_fill_manual(
        values = c(
            RNA = "#3F83B9",
            ATAC = "#F2181B"
        ),
        guide = "none"
    ) +
    scale_y_continuous(
        breaks = c(
            0,
            25,
            50,
            75,
            100
        ),
        limits = c(
            0,
            112
        ),
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +
    scale_x_discrete(
        expand = expansion(
            mult = c(
                0.35,
                0.35
            )
        )
    ) +
    labs(
        title = "DESeq2 LRT\nBH FDR < 0.01",
        x = NULL,
        y = "Genes/peaks with tissue effect (%)",
        tag = "A"
    ) +
    base_theme +
    theme(
        plot.title = element_text(
            size = 9.8,
            face = "bold",
            hjust = 0.5,
            lineheight = 0.95,
            margin = margin(
                b = 8
            )
        ),
        axis.text.x = element_text(
            size = 9
        ),
        plot.margin = margin(
            14,
            16,
            8,
            12
        )
    )

# ============================================================
# 6. PANEL B — RNA HC COUNTS
# ============================================================

rna_hc_plot <- rna_hc %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = rev(tissue_order)
        )
    )

pB <- ggplot(
    rna_hc_plot,
    aes(
        x = Count,
        y = Tissue,
        fill = Tissue
    )
) +
    geom_col(
        width = 0.70
    ) +
    geom_text(
        aes(
            label = comma(Count)
        ),
        hjust = -0.10,
        size = 3.0
    ) +
    scale_fill_manual(
        values = tissue_colors,
        guide = "none"
    ) +
    scale_x_continuous(
        breaks = c(
            0,
            500,
            1000
        ),
        labels = comma,
        limits = c(
            0,
            1350
        ),
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +
    labs(
        x = "HC tissue-specific RNA genes",
        y = NULL,
        tag = "B"
    ) +
    base_theme +
    theme(
        axis.text.y = element_text(
            size = 8
        ),
        axis.title.x = element_text(
            margin = margin(
                t = 6
            )
        ),
        plot.margin = margin(
            14,
            12,
            8,
            10
        )
    )

# ============================================================
# 7. PANEL C — ATAC HC COUNTS
#
# MAIN FIX:
# use compact 0 / 5k / 10k / ... labels
# ============================================================

atac_hc_plot <- atac_hc %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = rev(tissue_order)
        )
    )

pC <- ggplot(
    atac_hc_plot,
    aes(
        x = Count,
        y = Tissue,
        fill = Tissue
    )
) +
    geom_col(
        width = 0.70
    ) +
    geom_text(
        aes(
            label = comma(Count)
        ),
        hjust = -0.10,
        size = 3.0
    ) +
    scale_fill_manual(
        values = tissue_colors,
        guide = "none"
    ) +
    scale_x_continuous(
        breaks = c(
            0,
            5000,
            10000,
            15000,
            20000,
            25000
        ),
        labels = c(
            "0",
            "5k",
            "10k",
            "15k",
            "20k",
            "25k"
        ),
        limits = c(
            0,
            25500
        ),
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +
    labs(
        x = "HC tissue-specific ATAC peaks",
        y = NULL,
        tag = "C"
    ) +
    base_theme +
    theme(
        axis.text.y = element_text(
            size = 8
        ),
        axis.text.x = element_text(
            size = 8
        ),
        axis.title.x = element_text(
            margin = margin(
                t = 6
            )
        ),
        plot.margin = margin(
            14,
            18,
            8,
            10
        )
    )

# ============================================================
# 8. PANELS D/E — FROZEN TAU-SPM DISTRIBUTIONS
# ============================================================

make_specificity_plot <- function(
    dat,
    tag
) {

    ggplot(
        dat,
        aes(
            x = Tau,
            y = MaxSPM
        )
    ) +
        geom_bin_2d(
            bins = 65
        ) +
        geom_vline(
            xintercept = 0.80,
            linetype = 2,
            linewidth = 0.45
        ) +
        annotate(
            "text",
            x = 0.805,
            y = 0.965,
            label = "tau = 0.80",
            hjust = 0,
            vjust = 1,
            size = 3.0
        ) +
        scale_fill_viridis_c(
            option = "plasma",
            trans = "log10",
            breaks = c(
                1,
                10,
                100,
                1000
            ),
            labels = c(
                "1",
                "10",
                "100",
                "1000"
            ),
            name = "Features\nper bin"
        ) +
        scale_x_continuous(
            limits = c(
                0,
                1
            ),
            breaks = seq(
                0,
                1,
                0.2
            ),
            expand = c(
                0,
                0
            )
        ) +
        scale_y_continuous(
            limits = c(
                0,
                1
            ),
            breaks = seq(
                0,
                1,
                0.2
            ),
            expand = c(
                0,
                0
            )
        ) +
        labs(
            x = "Tissue specificity (tau)",
            y = "Maximum SPM",
            tag = tag
        ) +
        base_theme +
        theme(
            legend.title = element_text(
                size = 8.5,
                face = "bold"
            ),
            legend.text = element_text(
                size = 7.5
            ),
            legend.key.height = grid::unit(
                0.30,
                "in"
            ),
            plot.margin = margin(
                12,
                14,
                8,
                10
            )
        )
}

pD <- make_specificity_plot(
    rna_plot,
    "D"
)

pE <- make_specificity_plot(
    atac_plot,
    "E"
)

# ============================================================
# 9. FINAL LAYOUT
#
# Give Panel C substantially more width.
# ============================================================

top_row <- pA + pB + pC +
    plot_layout(
        widths = c(
            1.10,
            1.25,
            1.80
        )
    )

bottom_row <- pD + pE +
    plot_layout(
        widths = c(
            1,
            1
        ),
        guides = "collect"
    ) &
    theme(
        legend.position = "right"
    )

final_plot <- top_row / bottom_row +
    plot_layout(
        heights = c(
            0.88,
            1.12
        )
    )

# ============================================================
# 10. QC
# ============================================================

qc <- tibble(
    Metric = c(
        "RNA_LRT_significant",
        "RNA_LRT_tested",
        "ATAC_LRT_significant",
        "ATAC_LRT_tested",
        "RNA_HC_total",
        "ATAC_HC_total",
        "RNA_specificity_source_rows",
        "ATAC_specificity_source_rows",
        "RNA_finite_Tau_SPM_rows",
        "ATAC_finite_Tau_SPM_rows"
    ),
    Value = c(
        lrt$Sig[1],
        lrt$Tested[1],
        lrt$Sig[2],
        lrt$Tested[2],
        sum(rna_hc$Count),
        sum(atac_hc$Count),
        nrow(rna_obj$data),
        nrow(atac_obj$data),
        nrow(rna_plot),
        nrow(atac_plot)
    ),
    Expected = c(
        18870,
        22240,
        176917,
        203210,
        5036,
        40545,
        35682,
        203226,
        31063,
        203225
    )
) %>%
    mutate(
        Difference = Value - Expected,
        Status = if_else(
            Difference == 0,
            "PASS",
            "FAIL"
        )
    )

write_tsv(
    qc,
    file.path(
        outdir,
        "Figure2_FINAL_QC.tsv"
    )
)

if (any(qc$Status != "PASS")) {

    print(qc)

    stop(
        "Figure2 frozen-result QC failed."
    )
}

# ============================================================
# 11. EXPORT
# ============================================================

pdf_file <- file.path(
    outdir,
    "Figure2_FINAL.pdf"
)

tiff_file <- file.path(
    outdir,
    "Figure2_FINAL_600dpi.tiff"
)

ggsave(
    filename = pdf_file,
    plot = final_plot,
    device = cairo_pdf,
    width = 12.4,
    height = 8.3,
    units = "in"
)

ggsave(
    filename = tiff_file,
    plot = final_plot,
    device = "tiff",
    width = 12.4,
    height = 8.3,
    units = "in",
    dpi = 600,
    compression = "lzw"
)

cat("\n============================================================\n")
cat("STEP10-preF FIGURE2 FINAL COMPLETED\n")
cat("============================================================\n\n")

cat("Frozen specificity sources:\n")
cat("RNA : ", rna_obj$file, "\n", sep = "")
cat("ATAC: ", atac_obj$file, "\n\n", sep = "")

cat("QC:\n")
print(qc)

cat("\nFinal files:\n")
cat(pdf_file, "\n")
cat(tiff_file, "\n")

cat("\nLAYOUT-ONLY FIXES:\n")
cat("- Panel A title split into two lines.\n")
cat("- Panel A count/percentage labels separated and centered.\n")
cat("- Panel A/B panel-letter collision removed.\n")
cat("- Panel C enlarged substantially.\n")
cat("- Panel C x-axis uses 0, 5k, 10k, 15k, 20k, 25k.\n")
cat("- D/E retain frozen Tau/SPM data and tau = 0.80 threshold.\n")
cat("- No LRT, Tau, SPM, HC membership or threshold rerun.\n")
cat("============================================================\n")

