#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

cat("============================================================\n")
cat("STEP10-preF — FIGURE5 FINAL MICRO-PATCH\n")
cat("Panel A + Panel C + ASCII minus only\n")
cat("Frozen biological results unchanged\n")
cat("============================================================\n\n")

# ============================================================
# 1. Load the already validated v3-FIX plotting objects
# ============================================================

src <- "10preF_redraw_Figure5_v3_FIX.R"

if (!file.exists(src)) {
    stop("Cannot find: ", src)
}

env <- new.env(parent = globalenv())

sys.source(
    src,
    envir = env
)

need <- c(
    "pA_fix",
    "pB_fix",
    "pC_fix",
    "pD_fix",
    "pE_fix"
)

miss <- need[
    !vapply(
        need,
        exists,
        logical(1),
        envir = env,
        inherits = FALSE
    )
]

if (length(miss) > 0) {
    stop(
        "Missing plotting objects: ",
        paste(miss, collapse = ", ")
    )
}

pA <- get("pA_fix", envir = env)
pB <- get("pB_fix", envir = env)
pC <- get("pC_fix", envir = env)
pD <- get("pD_fix", envir = env)
pE <- get("pE_fix", envir = env)


# ============================================================
# 2. PANEL A
# Slightly shrink workflow geom_text labels.
# Final layout below also gives Panel A more vertical height.
# ============================================================

for (i in seq_along(pA$layers)) {

    if (inherits(
        pA$layers[[i]]$geom,
        "GeomText"
    )) {

        old_size <- pA$layers[[i]]$aes_params$size

        if (!is.null(old_size) &&
            is.numeric(old_size)) {

            pA$layers[[i]]$aes_params$size <-
                old_size * 0.93
        }
    }
}

pA <- pA +
    theme(
        plot.margin = margin(
            t = 8,
            r = 8,
            b = 12,
            l = 8
        )
    )


# ============================================================
# 3. PANEL C
# Remove old x-scale and give right-side labels more space.
# ============================================================

if (length(pC$scales$scales) > 0) {

    keep_scale <- vapply(
        pC$scales$scales,
        function(s) {

            aes <- s$aesthetics

            if (is.null(aes)) {
                return(TRUE)
            }

            !any(aes %in% "x")
        },
        logical(1)
    )

    pC$scales$scales <-
        pC$scales$scales[keep_scale]
}

pC <- pC +
    scale_x_continuous(
        limits = c(
            0,
            0.88
        ),
        breaks = seq(
            0,
            0.8,
            by = 0.2
        ),
        expand = expansion(
            mult = c(
                0.015,
                0.04
            )
        )
    ) +
    coord_cartesian(
        clip = "off"
    ) +
    theme(
        plot.margin = margin(
            t = 12,
            r = 34,
            b = 8,
            l = 20
        )
    )


# ============================================================
# 4. Replace Unicode minus with ordinary ASCII "-"
# Works recursively for ggplot / patchwork objects.
# ============================================================

unicode_minus <- intToUtf8(8722)

ascii_minus <- function(x) {

    if (is.character(x)) {

        return(
            gsub(
                unicode_minus,
                "-",
                x,
                fixed = TRUE
            )
        )
    }

    x
}


sanitize_plot_text <- function(p) {

    # ggplot labels
    # IMPORTANT:
    # ggplot2 >= 4 keeps p$labels as a dedicated <ggplot2::labels>
    # object. Modify individual elements instead of replacing the
    # whole object with lapply(), which would convert it to a list.
    if (!is.null(p$labels)) {

        label_fields <- c(
            "title",
            "subtitle",
            "caption",
            "x",
            "y",
            "tag"
        )

        for (nm in label_fields) {

            val <- tryCatch(
                p$labels[[nm]],
                error = function(e) NULL
            )

            if (is.character(val)) {

                p$labels[[nm]] <-
                    ascii_minus(val)
            }
        }
    }

    # main data
    if (is.data.frame(p$data)) {

        for (nm in names(p$data)) {

            if (is.character(p$data[[nm]])) {

                p$data[[nm]] <-
                    ascii_minus(
                        p$data[[nm]]
                    )
            }
        }
    }

    # layer data / fixed labels
    if (!is.null(p$layers)) {

        for (i in seq_along(p$layers)) {

            ld <- p$layers[[i]]$data

            if (is.data.frame(ld)) {

                for (nm in names(ld)) {

                    if (is.character(ld[[nm]])) {

                        ld[[nm]] <-
                            ascii_minus(
                                ld[[nm]]
                            )
                    }
                }

                p$layers[[i]]$data <- ld
            }

            if (!is.null(
                p$layers[[i]]$aes_params$label
            )) {

                p$layers[[i]]$aes_params$label <-
                    ascii_minus(
                        p$layers[[i]]$aes_params$label
                    )
            }
        }
    }

    # patchwork children
    if (!is.null(p$patches$plots)) {

        p$patches$plots <- lapply(
            p$patches$plots,
            sanitize_plot_text
        )
    }

    # patchwork annotation
    if (!is.null(p$patches$annotation)) {

        for (nm in names(
            p$patches$annotation
        )) {

            p$patches$annotation[[nm]] <-
                ascii_minus(
                    p$patches$annotation[[nm]]
                )
        }
    }

    p
}

pA <- sanitize_plot_text(pA)
pB <- sanitize_plot_text(pB)
pC <- sanitize_plot_text(pC)
pD <- sanitize_plot_text(pD)
pE <- sanitize_plot_text(pE)


# ============================================================
# 5. FINAL LAYOUT
#
# Only meaningful layout change:
# Panel A gets a little more height.
# ============================================================

pCD <- pC + pD +
    plot_layout(
        nrow = 1,
        widths = c(
            0.95,
            1.20
        )
    )

Figure5_FINAL <- (
    pA /
    pB /
    pCD /
    pE
) +
    plot_layout(
        heights = c(
            1.18,
            1.55,
            1.55,
            1.90
        )
    ) &
    theme(
        plot.background = element_rect(
            fill = "white",
            colour = NA
        )
    )


# ============================================================
# 6. OUTPUT
# ============================================================

outdir <- paste0(
    "10_publication_figures/",
    "Step10preF_Figure5_FINAL"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

pdf_file <- file.path(
    outdir,
    "Figure5_FINAL.pdf"
)

tiff_file <- file.path(
    outdir,
    "Figure5_FINAL_600dpi.tiff"
)

ggsave(
    pdf_file,
    Figure5_FINAL,
    width = 12.0,
    height = 9.9,
    units = "in",
    device = cairo_pdf,
    bg = "white"
)

ggsave(
    tiff_file,
    Figure5_FINAL,
    width = 12.0,
    height = 9.9,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)


# ============================================================
# 7. Preserve authoritative frozen QC
# ============================================================

qc_source <- paste0(
    "10_publication_figures/",
    "Step10preF_Figure5_v3_FIX/",
    "10preF_v3_FIX_Figure5_QC.tsv"
)

qc_target <- file.path(
    outdir,
    "Figure5_FINAL_QC.tsv"
)

if (!file.exists(qc_source)) {
    stop(
        "Cannot find frozen QC: ",
        qc_source
    )
}

file.copy(
    qc_source,
    qc_target,
    overwrite = TRUE
)

qc <- read.delim(
    qc_target,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!all(qc$Status == "PASS")) {
    stop("Frozen QC contains non-PASS status.")
}


# ============================================================
# 8. Final checks
# ============================================================

if (!file.exists(pdf_file)) {
    stop("PDF not generated.")
}

if (!file.exists(tiff_file)) {
    stop("TIFF not generated.")
}

cat("\n============================================================\n")
cat("FIGURE5 FINAL MICRO-PATCH COMPLETED\n")
cat("============================================================\n")

cat("Frozen QC       : ",
    sum(qc$Status == "PASS"),
    "/",
    nrow(qc),
    " PASS\n",
    sep = ""
)

cat("\nChanges only:\n")
cat("- Panel A: + vertical room; geom_text size x0.93\n")
cat("- Panel C: x-axis extended to 0.88; clipping disabled\n")
cat("- Unicode minus converted to ASCII '-'\n")
cat("- No frozen biological result changed\n")
cat("- No motif/FIMO/overlap/enrichment analysis rerun\n")

cat("\nFinal files:\n")
cat(pdf_file, "\n")
cat(tiff_file, "\n")
cat(qc_target, "\n")

cat("============================================================\n")
