#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("STEP10F3A1 — BUILD SUPPLEMENTARY FIGURE S1\n")
cat("ATAC-seq quality control\n")
cat("Frozen QC figures only; publication assembly only\n")
cat("============================================================\n\n")

# ============================================================
# 0. Packages
# ============================================================

required_pkgs <- c("magick", "pdftools")

missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
    stop(
        "Missing R package(s): ",
        paste(missing_pkgs, collapse = ", "),
        "\nInstall them before running Step10F3A1."
    )
}

library(magick)
library(pdftools)

# ============================================================
# 1. Authoritative frozen sources
# ============================================================

src <- list(

    A = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/fragments_distribution/",
        "pig_ATAC_fragment_size.pdf"
    ),

    B = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/plotCorrelation/",
        "pig_sample_correlation.pdf"
    ),

    C = paste0(
        "/datadisk2/liujingjian_data/ATAC-seq/Sus/",
        "ATAC-qc/TSS_enrichment_plot.pdf"
    ),

    D1 = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR54_58_tss_heatmap.pdf"
    ),

    D2 = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR61_66_tss_heatmap.pdf"
    ),

    D3 = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR69_74_tss_heatmap.pdf"
    ),

    D4 = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR77_82_tss_heatmap.pdf"
    )
)

# Companion source files: provenance only, not re-analysed here
companion <- list(

    A_raw = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/fragments_distribution/",
        "pig_ATAC_fragment_size_raw.tsv"
    ),

    A_summary = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/fragments_distribution/",
        "pig_ATAC_fragment_size_summary.tsv"
    ),

    B_correlation_table = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/plotCorrelation/",
        "pig_sample_correlation.tab"
    ),

    B_summary = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/plotCorrelation/",
        "bw_summary.gz"
    ),

    D1_matrix = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR54_58_tss_matrix.gz"
    ),

    D2_matrix = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR61_66_tss_matrix.gz"
    ),

    D3_matrix = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR69_74_tss_matrix.gz"
    ),

    D4_matrix = paste0(
        "/datadisk2/liujingjian_data/cattle_and_sus/",
        "sus_atac_analysis/bam2bigwig/reference_point_plot/",
        "SRR77_82_tss_matrix.gz"
    )
)

# ============================================================
# 2. Output directory
# ============================================================

outdir <- paste0(
    "10_publication_figures/",
    "Step10F_final_consistency/",
    "Step10F3A_ATAC_QC"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pdf_out <- file.path(
    outdir,
    "Supplementary_Figure_S1.pdf"
)

tiff_out <- file.path(
    outdir,
    "Supplementary_Figure_S1_600dpi.tiff"
)

preview_out <- file.path(
    outdir,
    "Supplementary_Figure_S1_preview.png"
)

manifest_out <- file.path(
    outdir,
    "10F3A1_S1_source_manifest.tsv"
)

companion_out <- file.path(
    outdir,
    "10F3A1_S1_companion_source_manifest.tsv"
)

qc_out <- file.path(
    outdir,
    "10F3A1_S1_QC.tsv"
)

caption_out <- file.path(
    outdir,
    "Supplementary_Figure_S1_caption.md"
)

# ============================================================
# 3. Frozen source preflight
# ============================================================

source_paths <- unlist(src, use.names = TRUE)

missing_sources <- source_paths[!file.exists(source_paths)]

if (length(missing_sources) > 0) {

    cat("Missing frozen source files:\n")
    cat(paste0("  ", missing_sources, "\n"))

    stop("STEP10F3A1 stopped because frozen source PDF(s) are missing.")
}

zero_sources <- source_paths[file.info(source_paths)$size <= 0]

if (length(zero_sources) > 0) {
    stop(
        "Zero-byte frozen source file(s):\n",
        paste(zero_sources, collapse = "\n")
    )
}

pdf_meta <- lapply(source_paths, function(x) {

    inf <- pdftools::pdf_info(x)

    data.frame(
        Pages = inf$pages,
        Width_pt = inf$width,
        Height_pt = inf$height,
        stringsAsFactors = FALSE
    )
})

pdf_meta <- do.call(rbind, pdf_meta)
rownames(pdf_meta) <- NULL

if (any(pdf_meta$Pages != 1)) {
    stop(
        "Every frozen source PDF must contain exactly one page. ",
        "At least one source does not."
    )
}

source_md5_before <- unname(tools::md5sum(source_paths))

cat("Frozen source preflight: PASS\n")
cat("Frozen source PDFs       :", length(source_paths), "\n\n")

# ============================================================
# 4. Source manifest
# ============================================================

panel_names <- c(
    "S1A",
    "S1B",
    "S1C",
    "S1D1",
    "S1D2",
    "S1D3",
    "S1D4"
)

panel_contents <- c(
    "ATAC fragment-size distribution",
    "ATAC sample correlation heatmap",
    "TSS enrichment profile",
    "TSS-centered heatmap: Adipose + Cerebellum",
    "TSS-centered heatmap: Cortex + Hypothalamus",
    "TSS-centered heatmap: Liver + Lung",
    "TSS-centered heatmap: Muscle + Spleen"
)

manifest <- data.frame(
    Panel = panel_names,
    Content = panel_contents,
    Source = unname(source_paths),
    Exists = file.exists(source_paths),
    Pages = pdf_meta$Pages,
    Width_pt = pdf_meta$Width_pt,
    Height_pt = pdf_meta$Height_pt,
    Size_bytes = file.info(source_paths)$size,
    MD5_before = source_md5_before,
    Frozen_status = "FROZEN",
    Step10F3A1_action = "READ_AND_ASSEMBLE_ONLY",
    stringsAsFactors = FALSE
)

write.table(
    manifest,
    manifest_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Companion provenance manifest

companion_paths <- unlist(companion, use.names = TRUE)

companion_manifest <- data.frame(
    Source_ID = names(companion_paths),
    File = unname(companion_paths),
    Exists = file.exists(companion_paths),
    Size_bytes = ifelse(
        file.exists(companion_paths),
        file.info(companion_paths)$size,
        NA
    ),
    MD5 = NA_character_,
    Step10F3A1_action = "PROVENANCE_ONLY_NOT_ANALYSED",
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(companion_manifest))) {

    f <- companion_manifest$File[i]

    if (file.exists(f)) {
        companion_manifest$MD5[i] <- unname(tools::md5sum(f))
    }
}

write.table(
    companion_manifest,
    companion_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. Helper functions
# ============================================================

read_frozen_pdf <- function(path, density = 450) {

    cat("Rendering frozen PDF:\n  ", path, "\n", sep = "")

    img <- magick::image_read_pdf(
        path,
        density = density
    )

    if (length(img) != 1) {
        img <- img[1]
    }

    # Crop only external blank margin.
    # This does not alter any data layer.
    img <- magick::image_trim(
        img,
        fuzz = 4
    )

    img
}


fit_image_to_box <- function(
    img,
    width,
    height,
    top_band = 115,
    side_pad = 45,
    bottom_pad = 35
) {

    info <- magick::image_info(img)

    available_w <- width - 2 * side_pad
    available_h <- height - top_band - bottom_pad

    scale_factor <- min(
        available_w / info$width,
        available_h / info$height
    )

    new_w <- max(1, floor(info$width * scale_factor))
    new_h <- max(1, floor(info$height * scale_factor))

    img2 <- magick::image_resize(
        img,
        geometry = sprintf("%dx%d!", new_w, new_h)
    )

    xoff <- floor((width - new_w) / 2)

    yoff <- top_band +
        floor((available_h - new_h) / 2)

    canvas <- magick::image_blank(
        width = width,
        height = height,
        color = "white"
    )

    magick::image_composite(
        canvas,
        img2,
        offset = sprintf("+%d+%d", xoff, yoff)
    )
}


add_panel_letter <- function(
    img,
    letter,
    x = 35,
    y = 22,
    size = 76
) {

    magick::image_annotate(
        img,
        text = letter,
        gravity = "northwest",
        location = sprintf("+%d+%d", x, y),
        size = size,
        weight = 700,
        color = "black"
    )
}


add_subpanel_title <- function(
    img,
    title,
    y = 23,
    size = 52
) {

    magick::image_annotate(
        img,
        text = title,
        gravity = "north",
        location = sprintf("+0+%d", y),
        size = size,
        weight = 600,
        color = "black"
    )
}

# ============================================================
# 6. Render frozen source PDFs
# ============================================================

# Densities are rendering choices only.
# They do not change any source value or analysis result.

img_A <- read_frozen_pdf(
    src$A,
    density = 500
)

img_B <- read_frozen_pdf(
    src$B,
    density = 450
)

img_C <- read_frozen_pdf(
    src$C,
    density = 600
)

img_D1 <- read_frozen_pdf(
    src$D1,
    density = 450
)

img_D2 <- read_frozen_pdf(
    src$D2,
    density = 450
)

img_D3 <- read_frozen_pdf(
    src$D3,
    density = 450
)

img_D4 <- read_frozen_pdf(
    src$D4,
    density = 450
)

# ============================================================
# 7. Final publication geometry
# ============================================================

# Final output:
# width  = 7200 px = 12.0 in at 600 dpi
# height = 7560 px = 12.6 in at 600 dpi

FINAL_W <- 7200
FINAL_H <- 7560

GAP <- 80

TOP_H <- 2500
C_H <- 1500
D_H <- 3400

AB_W <- floor((FINAL_W - GAP) / 2)

# ------------------------------------------------------------
# S1A
# ------------------------------------------------------------

panel_A <- fit_image_to_box(
    img_A,
    width = AB_W,
    height = TOP_H,
    top_band = 120,
    side_pad = 45,
    bottom_pad = 35
)

panel_A <- add_panel_letter(
    panel_A,
    "A"
)

# ------------------------------------------------------------
# S1B
# ------------------------------------------------------------

panel_B <- fit_image_to_box(
    img_B,
    width = AB_W,
    height = TOP_H,
    top_band = 120,
    side_pad = 45,
    bottom_pad = 35
)

panel_B <- add_panel_letter(
    panel_B,
    "B"
)

# ------------------------------------------------------------
# S1C
# ------------------------------------------------------------

panel_C <- fit_image_to_box(
    img_C,
    width = FINAL_W,
    height = C_H,
    top_band = 115,
    side_pad = 55,
    bottom_pad = 25
)

panel_C <- add_panel_letter(
    panel_C,
    "C"
)

# ============================================================
# 8. S1D — four TSS-centered heatmap blocks
# ============================================================

D_HEADER <- 100
D_GAP <- 80

D_CELL_W <- floor((FINAL_W - D_GAP) / 2)
D_CELL_H <- 1600

make_D_cell <- function(img, title) {

    p <- fit_image_to_box(
        img,
        width = D_CELL_W,
        height = D_CELL_H,
        top_band = 105,
        side_pad = 35,
        bottom_pad = 25
    )

    add_subpanel_title(
        p,
        title = title,
        y = 20,
        size = 52
    )
}

cell_D1 <- make_D_cell(
    img_D1,
    "Adipose + Cerebellum"
)

cell_D2 <- make_D_cell(
    img_D2,
    "Cortex + Hypothalamus"
)

cell_D3 <- make_D_cell(
    img_D3,
    "Liver + Lung"
)

cell_D4 <- make_D_cell(
    img_D4,
    "Muscle + Spleen"
)

panel_D <- magick::image_blank(
    width = FINAL_W,
    height = D_H,
    color = "white"
)

panel_D <- add_panel_letter(
    panel_D,
    "D",
    x = 35,
    y = 15,
    size = 76
)

D_row1_y <- D_HEADER
D_row2_y <- D_HEADER + D_CELL_H + D_GAP

panel_D <- magick::image_composite(
    panel_D,
    cell_D1,
    offset = sprintf("+0+%d", D_row1_y)
)

panel_D <- magick::image_composite(
    panel_D,
    cell_D2,
    offset = sprintf("+%d+%d", D_CELL_W + D_GAP, D_row1_y)
)

panel_D <- magick::image_composite(
    panel_D,
    cell_D3,
    offset = sprintf("+0+%d", D_row2_y)
)

panel_D <- magick::image_composite(
    panel_D,
    cell_D4,
    offset = sprintf("+%d+%d", D_CELL_W + D_GAP, D_row2_y)
)

# ============================================================
# 9. Assemble Supplementary Figure S1
# ============================================================

final_canvas <- magick::image_blank(
    width = FINAL_W,
    height = FINAL_H,
    color = "white"
)

# A
final_canvas <- magick::image_composite(
    final_canvas,
    panel_A,
    offset = "+0+0"
)

# B
final_canvas <- magick::image_composite(
    final_canvas,
    panel_B,
    offset = sprintf("+%d+0", AB_W + GAP)
)

# C
C_Y <- TOP_H + GAP

final_canvas <- magick::image_composite(
    final_canvas,
    panel_C,
    offset = sprintf("+0+%d", C_Y)
)

# D
D_Y <- TOP_H + GAP + C_H + GAP

final_canvas <- magick::image_composite(
    final_canvas,
    panel_D,
    offset = sprintf("+0+%d", D_Y)
)

# ============================================================
# 10. Export publication files
# ============================================================

cat("\nWriting final PDF...\n")

magick::image_write(
    final_canvas,
    path = pdf_out,
    format = "pdf",
    density = "600x600"
)

cat("Writing final 600-dpi TIFF...\n")

# Try LZW compression first.
tryCatch(
    {
        magick::image_write(
            final_canvas,
            path = tiff_out,
            format = "tiff",
            density = "600x600",
            defines = c(
                "tiff:compression" = "lzw"
            )
        )
    },
    error = function(e) {

        warning(
            "LZW write option was not accepted by this magick build. ",
            "Writing TIFF without explicit compression."
        )

        magick::image_write(
            final_canvas,
            path = tiff_out,
            format = "tiff",
            density = "600x600"
        )
    }
)

cat("Writing preview PNG...\n")

preview <- magick::image_resize(
    final_canvas,
    geometry = "1800"
)

magick::image_write(
    preview,
    path = preview_out,
    format = "png"
)

# ============================================================
# 11. Verify frozen sources were unchanged
# ============================================================

source_md5_after <- unname(
    tools::md5sum(source_paths)
)

source_identity <- source_md5_before == source_md5_after

manifest$MD5_after <- source_md5_after
manifest$MD5_identity <- ifelse(
    source_identity,
    "PASS",
    "FAIL"
)

write.table(
    manifest,
    manifest_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 12. Output technical QC
# ============================================================

pdf_exists <- file.exists(pdf_out)
tiff_exists <- file.exists(tiff_out)
preview_exists <- file.exists(preview_out)

pdf_pages <- NA_integer_

if (pdf_exists) {
    pdf_pages <- pdftools::pdf_info(pdf_out)$pages
}

final_info <- magick::image_info(final_canvas)

# Try to verify TIFF dimensions and density using ImageMagick CLI.
tiff_width <- FINAL_W
tiff_height <- FINAL_H
tiff_xdpi <- 600
tiff_ydpi <- 600

identify_bin <- Sys.which("identify")

if (nzchar(identify_bin) && tiff_exists) {

    identify_cmd <- sprintf(
        "%s -format '%%w\\t%%h\\t%%x\\t%%y\\t%%U' %s",
        shQuote(identify_bin),
        shQuote(tiff_out)
    )

    identify_result <- tryCatch(
        system(identify_cmd, intern = TRUE),
        error = function(e) character()
    )

    if (length(identify_result) == 1) {

        fields <- strsplit(
            identify_result,
            "\t",
            fixed = TRUE
        )[[1]]

        if (length(fields) >= 5) {

            tiff_width <- suppressWarnings(
                as.numeric(fields[1])
            )

            tiff_height <- suppressWarnings(
                as.numeric(fields[2])
            )

            tiff_xdpi <- suppressWarnings(
                as.numeric(fields[3])
            )

            tiff_ydpi <- suppressWarnings(
                as.numeric(fields[4])
            )

            # If ImageMagick reports PixelsPerCentimeter,
            # convert to dpi.
            if (
                length(fields) >= 5 &&
                grepl(
                    "PixelsPerCentimeter",
                    fields[5],
                    fixed = TRUE
                )
            ) {

                tiff_xdpi <- tiff_xdpi * 2.54
                tiff_ydpi <- tiff_ydpi * 2.54
            }
        }
    }
}

qc <- data.frame(

    Metric = c(
        "Frozen_source_PDFs",
        "Frozen_source_PDFs_existing",
        "Frozen_source_PDFs_single_page",
        "Frozen_source_MD5_failures",
        "Panel_A_present",
        "Panel_B_present",
        "Panel_C_present",
        "Panel_D_heatmap_blocks",
        "Final_width_px",
        "Final_height_px",
        "Final_width_in_at_600dpi",
        "Final_height_in_at_600dpi",
        "Final_PDF_exists",
        "Final_PDF_pages",
        "Final_TIFF_exists",
        "Final_TIFF_Xdpi",
        "Final_TIFF_Ydpi",
        "Preview_PNG_exists",
        "ATAC_QC_recomputed",
        "Biological_results_modified"
    ),

    Value = c(
        length(source_paths),
        sum(file.exists(source_paths)),
        sum(pdf_meta$Pages == 1),
        sum(!source_identity),
        1,
        1,
        1,
        4,
        tiff_width,
        tiff_height,
        round(tiff_width / 600, 3),
        round(tiff_height / 600, 3),
        as.integer(pdf_exists),
        pdf_pages,
        as.integer(tiff_exists),
        round(tiff_xdpi, 2),
        round(tiff_ydpi, 2),
        as.integer(preview_exists),
        0,
        0
    ),

    Expected = c(
        7,
        7,
        7,
        0,
        1,
        1,
        1,
        4,
        7200,
        7560,
        12.000,
        12.600,
        1,
        1,
        1,
        600,
        600,
        1,
        0,
        0
    ),

    stringsAsFactors = FALSE
)

tol <- c(
    rep(0, 10),
    0.01,
    0.01,
    0,
    0,
    0,
    1,
    1,
    0,
    0,
    0
)

qc$Difference <- abs(
    suppressWarnings(
        as.numeric(qc$Value) -
        as.numeric(qc$Expected)
    )
)

qc$Status <- ifelse(
    qc$Difference <= tol,
    "PASS",
    "FAIL"
)

write.table(
    qc,
    qc_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 13. Supplementary Figure S1 caption
# ============================================================

caption <- c(
    "# Supplementary Figure S1",
    "",
    "## Supplementary Figure S1. ATAC-seq quality control across the 16 pig tissue samples.",
    "",
    "(A) ATAC-seq fragment-size distribution across samples.",
    "(B) Pairwise sample-correlation heatmap derived from genome-wide ATAC-seq signal profiles.",
    "(C) Aggregate transcription-start-site (TSS) enrichment profile.",
    "(D) TSS-centered ATAC-seq signal heatmaps for the eight tissues, displayed as four paired tissue blocks: Adipose and Cerebellum, Cortex and Hypothalamus, Liver and Lung, and Muscle and Spleen.",
    "These panels summarize previously generated ATAC-seq quality-control outputs. Step10F3A1 performed publication assembly only; fragment-size, correlation, TSS-enrichment, and TSS-centered signal analyses were not recomputed."
)

writeLines(
    caption,
    caption_out
)

# ============================================================
# 14. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10F3A1 COMPLETED\n")
cat("============================================================\n\n")

cat("Final panel structure:\n")
cat("S1A : ATAC fragment-size distribution\n")
cat("S1B : 16-sample ATAC correlation heatmap\n")
cat("S1C : aggregate TSS enrichment profile\n")
cat("S1D : four TSS-centered heatmap blocks\n\n")

cat("S1D blocks:\n")
cat("  D1 Adipose + Cerebellum\n")
cat("  D2 Cortex + Hypothalamus\n")
cat("  D3 Liver + Lung\n")
cat("  D4 Muscle + Spleen\n\n")

cat("Final files:\n")
cat(pdf_out, "\n")
cat(tiff_out, "\n")
cat(preview_out, "\n")
cat(manifest_out, "\n")
cat(companion_out, "\n")
cat(qc_out, "\n")
cat(caption_out, "\n\n")

cat("QC:\n")
print(qc, row.names = FALSE)

cat("\n")

if (all(qc$Status == "PASS")) {

    cat("STEP10F3A1 STATUS: PASS\n")
    cat("Supplementary Figure S1 is technically assembled.\n")

} else {

    cat("STEP10F3A1 STATUS: CHECK\n")
    cat("At least one technical QC item requires inspection.\n")
}

cat("\nIMPORTANT:\n")
cat("- No ATAC QC metric was recalculated.\n")
cat("- No BAM, bigWig, matrix, or QC table was modified.\n")
cat("- Source PDF MD5 values were checked before and after assembly.\n")
cat("- PDF rendering/cropping is publication layout only.\n")
cat("- Panel D contains four visualization blocks but remains one main panel.\n")
cat("- Biological results remain frozen.\n")
cat("============================================================\n")

