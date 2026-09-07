#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
    library(grid)
})

cat("============================================================\n")
cat("STEP10E-v3 FINAL — PUBLICATION POLISH FIGURE 6\n")
cat("Rebuild visualization objects from frozen Step09C Top5\n")
cat("============================================================\n\n")


# ============================================================
# 0. Input / output
# ============================================================

TOP5_FILE <- paste0(
    "09C_publication_candidate_modules/",
    "09C5_top5_primary_supporting_figure_candidates.tsv"
)

STEP10E_QC <- paste0(
    "10_publication_figures/Step10E_Figure6/",
    "10E7_Figure6_overall_QC_summary.tsv"
)

OUTDIR <- "10_publication_figures/Step10E_Figure6_FINAL"

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    file.path(OUTDIR, "Figure6A_Muscle"),
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    file.path(OUTDIR, "Figure6B_Spleen"),
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    file.path(OUTDIR, "Figure6C_Liver"),
    recursive = TRUE,
    showWarnings = FALSE
)


if (!file.exists(TOP5_FILE)) {
    stop("Missing authoritative Step09C Top5 file.")
}

if (!file.exists(STEP10E_QC)) {
    stop("Missing original Step10E QC file.")
}


# ============================================================
# 1. Read authoritative frozen data
# ============================================================

top5 <- read_tsv(
    TOP5_FILE,
    col_types = cols(.default = col_character()),
    na = c("", "NA"),
    show_col_types = FALSE,
    progress = FALSE
)

old_qc <- read_tsv(
    STEP10E_QC,
    col_types = cols(.default = col_character()),
    na = c("", "NA"),
    show_col_types = FALSE,
    progress = FALSE
)


if (nrow(top5) != 30) {
    stop(
        paste0(
            "Authoritative Step09C Top5 contains ",
            nrow(top5),
            " rows, expected 30."
        )
    )
}

if (!all(old_qc$Status == "PASS")) {
    print(old_qc)
    stop("Original Step10E QC is not completely PASS.")
}


# ============================================================
# 2. Required columns
# ============================================================

required_cols <- c(
    "SAF_peak_id",
    "gene_id",
    "gene_name",
    "Tissue",
    "Step07_Pair_specific_class",
    "Evidence_tier",
    "Module_ID",
    "Module_role",
    "Module_motif_ID",
    "Module_motif_name",
    "Module_motif_family_label",
    "Module_representative_candidate_TFs",
    "Publication_named_gene_rank_within_module"
)

missing_cols <- setdiff(
    required_cols,
    names(top5)
)

if (length(missing_cols) > 0) {

    cat("Missing columns:\n")
    cat(
        paste0(
            "  ",
            missing_cols,
            collapse = "\n"
        ),
        "\n"
    )

    stop("Step09C Top5 structure incompatible.")
}


# ============================================================
# 3. Frozen main-figure module definitions
# ============================================================

module_meta <- tribble(

    ~Panel,
    ~Tissue,
    ~Module_ID,
    ~Module_order,
    ~Expected_role,
    ~Motif_display,
    ~TF_display,
    ~Required_TFs,
    ~Module_color,
    ~Module_y,

    # --------------------------------------------------------
    # Muscle
    # --------------------------------------------------------

    "A",
    "Muscle",
    "Muscle__MA1107.3",
    1,
    "Primary_structural",
    "KLF9 motif\nMA1107.3",
    "KLF9\ncandidate TF",
    "KLF9",
    "#0072B2",
    3.55,

    "A",
    "Muscle",
    "Muscle__3-CCTCCGCCCCTG",
    2,
    "Supporting_q005_structural",
    "SP-like STREME3\nCCTCCGCCCCTG",
    "SP1 / SP4\nfamily candidates",
    "SP1;SP4",
    "#D55E00",
    5.55,


    # --------------------------------------------------------
    # Spleen
    # --------------------------------------------------------

    "B",
    "Spleen",
    "Spleen__1-MGGGGSAGGAGCMG",
    1,
    "Primary_integrated",
    "GC-rich STREME1\nMGGGGSAGGAGCMG",
    "ZNF148\nfamily candidate",
    "ZNF148",
    "#009E73",
    7.20,

    "B",
    "Spleen",
    "Spleen__2-GCCCAGCCCMGCCCM",
    2,
    "Primary_integrated",
    "GC-rich STREME2\nGCCCAGCCCMGCCCM",
    "KLF12 / SP1 / SP3\nSP4 / ZBED4\nfamily candidates",
    "KLF12;SP1;SP3;SP4;ZBED4",
    "#CC79A7",
    4.65,

    "B",
    "Spleen",
    "Spleen__MA0742.2",
    3,
    "Primary_integrated",
    "KLF12 motif\nMA0742.2",
    "KLF12\ncandidate TF",
    "KLF12",
    "#E69F00",
    2.10,


    # --------------------------------------------------------
    # Liver
    # --------------------------------------------------------

    "C",
    "Liver",
    "Liver__3-AGGGGGAGGGGRGGR",
    1,
    "Supporting_q005",
    "GC-rich STREME3\nAGGGGGAGGGGRGGR",
    "SP5\nfamily candidate",
    "SP5",
    "#7B3294",
    3.00
)


# ============================================================
# 4. Exact frozen representative gene sets
# ============================================================

expected_gene_sets <- list(

    "Muscle__MA1107.3" =
        c(
            "SLC2A4",
            "DOK5",
            "PHKG1",
            "KBTBD12",
            "TPM2"
        ),

    "Muscle__3-CCTCCGCCCCTG" =
        c(
            "SLC2A4",
            "SIX4",
            "FHL3",
            "DOK5",
            "PHKG1"
        ),

    "Spleen__1-MGGGGSAGGAGCMG" =
        c(
            "IKZF1",
            "MZB1",
            "GIMAP6",
            "CD22",
            "SASH3"
        ),

    "Spleen__2-GCCCAGCCCMGCCCM" =
        c(
            "IKZF1",
            "LYN",
            "MZB1",
            "GIMAP6",
            "CD22"
        ),

    "Spleen__MA0742.2" =
        c(
            "IKZF1",
            "MZB1",
            "CD22",
            "IL2RG",
            "CCR7"
        ),

    "Liver__3-AGGGGGAGGGGRGGR" =
        c(
            "CDO1",
            "CREB3L3",
            "HNF4A",
            "ARG1",
            "RIPK3"
        )
)


# ============================================================
# 5. Publication display gene order
#
# Only plotting order.
# No membership changes.
# ============================================================

panel_gene_order <- list(

    "A" = c(
        "SIX4",
        "FHL3",
        "KBTBD12",
        "TPM2",
        "SLC2A4",
        "DOK5",
        "PHKG1"
    ),

    "B" = c(
        "LYN",
        "IL2RG",
        "CCR7",
        "SASH3",
        "IKZF1",
        "MZB1",
        "CD22",
        "GIMAP6"
    ),

    "C" = c(
        "CDO1",
        "CREB3L3",
        "HNF4A",
        "ARG1",
        "RIPK3"
    )
)


# ============================================================
# 6. Select exactly six frozen modules
# ============================================================

selected <- top5 %>%

    filter(
        Module_ID %in%
            module_meta$Module_ID
    ) %>%

    inner_join(
        module_meta,
        by = c(
            "Module_ID",
            "Tissue"
        )
    ) %>%

    mutate(
        Publication_named_gene_rank_within_module =
            suppressWarnings(
                as.integer(
                    Publication_named_gene_rank_within_module
                )
            )
    ) %>%

    arrange(
        Panel,
        Module_order,
        Publication_named_gene_rank_within_module
    )


if (nrow(selected) != 30) {
    stop(
        paste0(
            "Selected Figure6 frozen rows = ",
            nrow(selected),
            ", expected 30."
        )
    )
}


# ============================================================
# 7. Frozen result QC
# ============================================================

# ------------------------------------------------------------
# 7A. Exactly 5 representatives/module
# ------------------------------------------------------------

module_n_qc <- selected %>%

    count(
        Module_ID,
        name = "Observed_N"
    ) %>%

    mutate(
        Expected_N = 5L,
        Status =
            if_else(
                Observed_N == Expected_N,
                "PASS",
                "FAIL"
            )
    )


if (any(module_n_qc$Status != "PASS")) {
    print(module_n_qc)
    stop("Top5 module-count QC failed.")
}


# ------------------------------------------------------------
# 7B. Exact frozen gene sets
# ------------------------------------------------------------

gene_qc_list <- list()

for (mid in names(expected_gene_sets)) {

    observed <- selected %>%
        filter(
            Module_ID == mid
        ) %>%
        pull(
            gene_name
        ) %>%
        unique() %>%
        sort()

    expected <- sort(
        expected_gene_sets[[mid]]
    )

    gene_qc_list[[mid]] <- tibble(
        Module_ID = mid,
        Observed_N = length(observed),
        Expected_N = length(expected),
        Observed_genes =
            paste(
                observed,
                collapse = ";"
            ),
        Expected_genes =
            paste(
                expected,
                collapse = ";"
            ),
        Status =
            ifelse(
                setequal(
                    observed,
                    expected
                ),
                "PASS",
                "FAIL"
            )
    )
}


gene_set_qc <- bind_rows(
    gene_qc_list
)


if (any(gene_set_qc$Status != "PASS")) {
    print(gene_set_qc)
    stop("Frozen representative gene-set QC failed.")
}


# ------------------------------------------------------------
# 7C. Frozen module roles
# ------------------------------------------------------------

role_qc <- selected %>%

    distinct(
        Module_ID,
        Module_role,
        Expected_role
    ) %>%

    mutate(
        Status =
            if_else(
                Module_role ==
                    Expected_role,
                "PASS",
                "FAIL"
            )
    )


if (any(role_qc$Status != "PASS")) {
    print(role_qc)
    stop("Frozen module-role QC failed.")
}


# ------------------------------------------------------------
# 7D. Candidate TF source support
# ------------------------------------------------------------

contains_all <- function(
    observed,
    required
) {

    if (
        is.na(observed) ||
        is.na(required)
    ) {
        return(FALSE)
    }

    observed_vec <-
        str_split(
            observed,
            fixed(";")
        )[[1]]

    required_vec <-
        str_split(
            required,
            fixed(";")
        )[[1]]

    all(
        required_vec %in%
            observed_vec
    )
}


tf_qc <- selected %>%

    group_by(
        Module_ID
    ) %>%

    summarise(
        Candidate_TFs =
            first(
                Module_representative_candidate_TFs
            ),
        Required_TFs =
            first(
                Required_TFs
            ),
        .groups = "drop"
    )


tf_qc$Status <- mapply(
    contains_all,
    tf_qc$Candidate_TFs,
    tf_qc$Required_TFs
)


tf_qc <- tf_qc %>%

    mutate(
        Status =
            if_else(
                Status,
                "PASS",
                "FAIL"
            )
    )


if (any(tf_qc$Status != "PASS")) {
    print(tf_qc)
    stop("Frozen candidate-TF source QC failed.")
}


# ============================================================
# 8. Network construction
#
# Visualization reconstruction only.
# ============================================================

build_panel <- function(
    panel_code
) {

    df <- selected %>%
        filter(
            Panel == panel_code
        )

    meta <- module_meta %>%
        filter(
            Panel == panel_code
        )

    tissue <- unique(
        df$Tissue
    )

    if (length(tissue) != 1) {
        stop("More than one tissue found in panel.")
    }


    # ========================================================
    # Gene nodes
    # ========================================================

    gene_order <-
        panel_gene_order[[panel_code]]


    gene_nodes <- df %>%

        distinct(
            Tissue,
            gene_id,
            gene_name
        ) %>%

        left_join(

            df %>%
                distinct(
                    Module_ID,
                    gene_id
                ) %>%
                count(
                    gene_id,
                    name = "Module_count"
                ),

            by = "gene_id"
        ) %>%

        mutate(
            Gene_order =
                match(
                    gene_name,
                    gene_order
                )
        )


    if (any(is.na(gene_nodes$Gene_order))) {

        print(gene_nodes)

        stop(
            paste0(
                "Gene display-order QC failed for panel ",
                panel_code
            )
        )
    }


    gene_nodes <- gene_nodes %>%

        arrange(
            Gene_order
        ) %>%

        mutate(
            x = 5.15,

            y =
                length(gene_order) -
                Gene_order +
                1,

            Shared_gene =
                Module_count > 1,

            node_id =
                paste0(
                    "Gene:",
                    Tissue,
                    ":",
                    gene_id
                ),

            node_type =
                "Associated_gene",

            label =
                gene_name,

            Module_ID =
                NA_character_,

            SAF_peak_id =
                NA_character_,

            Pair_specific_class =
                NA_character_,

            Evidence_tier =
                NA_character_
        )


    # ========================================================
    # Peak nodes
    # ========================================================

    peak_nodes <- df %>%

        distinct(
            Tissue,
            SAF_peak_id,
            gene_id,
            gene_name,
            Step07_Pair_specific_class,
            Evidence_tier
        ) %>%

        left_join(

            gene_nodes %>%
                select(
                    gene_id,
                    gene_name,
                    gene_y = y
                ),

            by = c(
                "gene_id",
                "gene_name"
            )
        ) %>%

        group_by(
            gene_id,
            gene_name
        ) %>%

        arrange(
            SAF_peak_id,
            .by_group = TRUE
        ) %>%

        mutate(
            Peak_N =
                n(),

            Peak_index =
                row_number(),

            Peak_offset =
                (
                    Peak_index -
                    (Peak_N + 1) / 2
                ) * 0.18,

            x = 3.65,

            y =
                gene_y +
                Peak_offset,

            node_id =
                paste0(
                    "Peak:",
                    Tissue,
                    ":",
                    SAF_peak_id
                ),

            node_type =
                "ATAC_peak",

            label =
                SAF_peak_id,

            Module_ID =
                NA_character_,

            Pair_specific_class =
                Step07_Pair_specific_class,

            Shared_gene =
                NA
        ) %>%

        ungroup()


    # ========================================================
    # Motif nodes
    # ========================================================

    motif_nodes <- meta %>%

        transmute(
            Panel,
            Tissue,
            Module_ID,
            Module_order,
            Module_color,

            x = 2.10,

            y = Module_y,

            node_id =
                paste0(
                    "Motif:",
                    Module_ID
                ),

            node_type =
                "Motif",

            label =
                Motif_display,

            gene_id =
                NA_character_,

            SAF_peak_id =
                NA_character_,

            Pair_specific_class =
                NA_character_,

            Evidence_tier =
                NA_character_,

            Shared_gene =
                NA
        )


    # ========================================================
    # Candidate TF annotation nodes
    # ========================================================

    tf_nodes <- meta %>%

        transmute(
            Panel,
            Tissue,
            Module_ID,
            Module_order,
            Module_color,

            x = 0.25,

            y = Module_y,

            node_id =
                paste0(
                    "TFAnnotation:",
                    Module_ID
                ),

            node_type =
                "Candidate_TF_annotation",

            label =
                TF_display,

            gene_id =
                NA_character_,

            SAF_peak_id =
                NA_character_,

            Pair_specific_class =
                NA_character_,

            Evidence_tier =
                NA_character_,

            Shared_gene =
                NA
        )


    # ========================================================
    # TF -> motif edges
    # ========================================================

    tf_motif_edges <- tf_nodes %>%

        select(
            Module_ID,
            source = node_id,
            x,
            y
        ) %>%

        left_join(

            motif_nodes %>%
                select(
                    Module_ID,
                    target = node_id,
                    xend = x,
                    yend = y
                ),

            by = "Module_ID"
        ) %>%

        mutate(
            edge_type =
                "candidate_TF_annotation_to_motif"
        )


    # ========================================================
    # Motif -> peak edges
    # ========================================================

    motif_peak_edges <- df %>%

        distinct(
            Module_ID,
            SAF_peak_id
        ) %>%

        left_join(

            motif_nodes %>%
                select(
                    Module_ID,
                    source = node_id,
                    x,
                    y
                ),

            by = "Module_ID"
        ) %>%

        left_join(

            peak_nodes %>%
                select(
                    SAF_peak_id,
                    target = node_id,
                    xend = x,
                    yend = y
                ),

            by = "SAF_peak_id"
        ) %>%

        mutate(
            edge_type =
                "motif_to_peak"
        )


    # ========================================================
    # Peak -> gene edges
    # ========================================================

    peak_gene_edges <- peak_nodes %>%

        left_join(

            gene_nodes %>%
                select(
                    gene_id,
                    gene_name,
                    target = node_id,
                    gene_x = x,
                    gene_y2 = y
                ),

            by = c(
                "gene_id",
                "gene_name"
            )
        ) %>%

        transmute(
            Module_ID =
                NA_character_,

            source =
                node_id,

            target,

            x,
            y,

            xend =
                gene_x,

            yend =
                gene_y2,

            edge_type =
                "peak_to_associated_gene"
        )


    # ========================================================
    # Combined nodes
    # ========================================================

    nodes <- bind_rows(

        tf_nodes %>%
            select(
                Panel,
                Tissue,
                node_id,
                node_type,
                label,
                Module_ID,
                x,
                y,
                gene_id,
                SAF_peak_id,
                Pair_specific_class,
                Evidence_tier,
                Shared_gene
            ),

        motif_nodes %>%
            select(
                Panel,
                Tissue,
                node_id,
                node_type,
                label,
                Module_ID,
                x,
                y,
                gene_id,
                SAF_peak_id,
                Pair_specific_class,
                Evidence_tier,
                Shared_gene
            ),

        peak_nodes %>%
            transmute(
                Panel =
                    panel_code,

                Tissue,

                node_id,
                node_type,
                label,
                Module_ID,

                x,
                y,

                gene_id,
                SAF_peak_id,
                Pair_specific_class,
                Evidence_tier,
                Shared_gene
            ),

        gene_nodes %>%
            transmute(
                Panel =
                    panel_code,

                Tissue,

                node_id,
                node_type,
                label,
                Module_ID,

                x,
                y,

                gene_id,
                SAF_peak_id,
                Pair_specific_class,
                Evidence_tier,
                Shared_gene
            )
    )


    # ========================================================
    # Combined edges
    # ========================================================

    edges <- bind_rows(
        tf_motif_edges,
        motif_peak_edges,
        peak_gene_edges
    ) %>%

        mutate(
            Panel =
                panel_code,

            Tissue =
                tissue,

            edge_id =
                sprintf(
                    "%s_E%03d",
                    panel_code,
                    row_number()
                )
        ) %>%

        select(
            Panel,
            Tissue,
            edge_id,
            source,
            target,
            edge_type,
            Module_ID,
            x,
            y,
            xend,
            yend
        )


    # ========================================================
    # Summary
    # ========================================================

    summary <- tibble(
        Panel =
            panel_code,

        Tissue =
            tissue,

        Module_N =
            n_distinct(
                df$Module_ID
            ),

        Selected_pair_rows =
            nrow(df),

        Unique_selected_peaks =
            n_distinct(
                df$SAF_peak_id
            ),

        Unique_selected_genes =
            n_distinct(
                df$gene_id
            ),

        Shared_gene_N =
            sum(
                gene_nodes$Shared_gene
            ),

        Promoter_peak_N =
            sum(
                peak_nodes$Pair_specific_class ==
                    "Promoter_TSSproximal",
                na.rm = TRUE
            ),

        Nonpromoter_peak_N =
            sum(
                peak_nodes$Pair_specific_class !=
                    "Promoter_TSSproximal" |
                    is.na(
                        peak_nodes$Pair_specific_class
                    )
            ),

        Plot_node_N =
            nrow(nodes),

        Plot_edge_N =
            nrow(edges)
    )


    list(
        df =
            df,

        nodes =
            nodes,

        edges =
            edges,

        summary =
            summary
    )
}


A <- build_panel("A")
B <- build_panel("B")
C <- build_panel("C")


# ============================================================
# 9. Frozen network count QC
# ============================================================

summary_all <- bind_rows(
    A$summary,
    B$summary,
    C$summary
)


expected_summary <- tribble(

    ~Panel,
    ~Expected_rows,
    ~Expected_peaks,
    ~Expected_genes,
    ~Expected_shared,
    ~Expected_nodes,
    ~Expected_edges,

    "A", 10L, 7L, 7L, 3L, 18L, 19L,
    "B", 15L, 8L, 8L, 4L, 22L, 26L,
    "C",  5L, 5L, 5L, 0L, 12L, 11L
)


summary_qc <- summary_all %>%

    left_join(
        expected_summary,
        by = "Panel"
    ) %>%

    mutate(
        Rows_status =
            if_else(
                Selected_pair_rows ==
                    Expected_rows,
                "PASS",
                "FAIL"
            ),

        Peaks_status =
            if_else(
                Unique_selected_peaks ==
                    Expected_peaks,
                "PASS",
                "FAIL"
            ),

        Genes_status =
            if_else(
                Unique_selected_genes ==
                    Expected_genes,
                "PASS",
                "FAIL"
            ),

        Shared_status =
            if_else(
                Shared_gene_N ==
                    Expected_shared,
                "PASS",
                "FAIL"
            ),

        Nodes_status =
            if_else(
                Plot_node_N ==
                    Expected_nodes,
                "PASS",
                "FAIL"
            ),

        Edges_status =
            if_else(
                Plot_edge_N ==
                    Expected_edges,
                "PASS",
                "FAIL"
            )
    )


status_cols <- grep(
    "_status$",
    names(summary_qc),
    value = TRUE
)


if (
    any(
        unlist(
            summary_qc[
                status_cols
            ]
        ) != "PASS"
    )
) {

    print(
        summary_qc,
        width = Inf
    )

    stop(
        "Frozen Figure6 reconstruction QC failed."
    )
}


nodes_all <- bind_rows(
    A$nodes,
    B$nodes,
    C$nodes
)

edges_all <- bind_rows(
    A$edges,
    B$edges,
    C$edges
)


if (nrow(nodes_all) != 52) {
    stop(
        paste0(
            "Rebuilt nodes = ",
            nrow(nodes_all),
            ", expected 52."
        )
    )
}

if (nrow(edges_all) != 56) {
    stop(
        paste0(
            "Rebuilt edges = ",
            nrow(edges_all),
            ", expected 56."
        )
    )
}


# ============================================================
# 10. Colors
# ============================================================

MODULE_COLORS <- setNames(
    module_meta$Module_color,
    module_meta$Module_ID
)


# ============================================================
# 11. Plotting function
# ============================================================

draw_panel <- function(
    object,
    title,
    subtitle
) {

    panel_code <-
        unique(
            object$nodes$Panel
        )

    nodes <-
        object$nodes

    edges <-
        object$edges


    tf_nodes <- nodes %>%
        filter(
            node_type ==
                "Candidate_TF_annotation"
        )


    motif_nodes <- nodes %>%
        filter(
            node_type ==
                "Motif"
        )


    peak_nodes <- nodes %>%
        filter(
            node_type ==
                "ATAC_peak"
        )


    gene_nodes <- nodes %>%
        filter(
            node_type ==
                "Associated_gene"
        )


    tf_edges <- edges %>%
        filter(
            edge_type ==
                "candidate_TF_annotation_to_motif"
        )


    motif_edges <- edges %>%
        filter(
            edge_type ==
                "motif_to_peak"
        )


    gene_edges <- edges %>%
        filter(
            edge_type ==
                "peak_to_associated_gene"
        )


    promoter_peaks <- peak_nodes %>%
        filter(
            Pair_specific_class ==
                "Promoter_TSSproximal"
        )


    other_peaks <- peak_nodes %>%
        filter(
            is.na(
                Pair_specific_class
            ) |
            Pair_specific_class !=
                "Promoter_TSSproximal"
        )


    shared_genes <- gene_nodes %>%
        filter(
            !is.na(
                Shared_gene
            ) &
            Shared_gene
        )


    nonshared_genes <- gene_nodes %>%
        filter(
            is.na(
                Shared_gene
            ) |
            !Shared_gene
        )


    max_y <-
        max(
            nodes$y,
            na.rm = TRUE
        ) + 0.90


    tf_size <-
        ifelse(
            panel_code == "B",
            1.82,
            1.98
        )


    motif_size <-
        ifelse(
            panel_code == "B",
            1.82,
            1.98
        )


    gene_size <-
        ifelse(
            panel_code == "B",
            2.35,
            2.48
        )


    ggplot() +

        # ----------------------------------------------------
        # Peak -> gene
        # ----------------------------------------------------

        geom_segment(
            data =
                gene_edges,

            aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend
            ),

            linewidth =
                0.48,

            colour =
                "#686868",

            alpha =
                0.82,

            lineend =
                "round"
        ) +


        # ----------------------------------------------------
        # Motif -> peak
        # ----------------------------------------------------

        geom_segment(
            data =
                motif_edges,

            aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend,
                colour = Module_ID
            ),

            linewidth =
                0.58,

            alpha =
                0.70,

            lineend =
                "round"
        ) +


        # ----------------------------------------------------
        # Candidate TF annotation -> motif
        # ----------------------------------------------------

        geom_segment(
            data =
                tf_edges,

            aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend,
                colour = Module_ID
            ),

            linewidth =
                0.62,

            linetype =
                "22",

            alpha =
                0.92,

            lineend =
                "round"
        ) +


        # ----------------------------------------------------
        # Non-promoter peaks
        # ----------------------------------------------------

        geom_point(
            data =
                other_peaks,

            aes(
                x = x,
                y = y
            ),

            shape =
                22,

            size =
                2.45,

            stroke =
                0.48,

            fill =
                "white",

            colour =
                "#404040"
        ) +


        # ----------------------------------------------------
        # Promoter/TSS-proximal peaks
        # ----------------------------------------------------

        geom_point(
            data =
                promoter_peaks,

            aes(
                x = x,
                y = y
            ),

            shape =
                22,

            size =
                2.45,

            stroke =
                0.48,

            fill =
                "#484848",

            colour =
                "#303030"
        ) +


        # ----------------------------------------------------
        # Candidate TF
        # ----------------------------------------------------

        geom_label(
            data =
                tf_nodes,

            aes(
                x = x,
                y = y,
                label = label,
                fill = Module_ID
            ),

            colour =
                "white",

            size =
                tf_size,

            lineheight =
                0.90,

            linewidth =
                0.28,

            label.padding =
                unit(
                    0.085,
                    "lines"
                )
        ) +


        # ----------------------------------------------------
        # Motifs
        # ----------------------------------------------------

        geom_label(
            data =
                motif_nodes,

            aes(
                x = x,
                y = y,
                label = label,
                fill = Module_ID
            ),

            colour =
                "white",

            size =
                motif_size,

            lineheight =
                0.90,

            linewidth =
                0.28,

            label.padding =
                unit(
                    0.085,
                    "lines"
                )
        ) +


        # ----------------------------------------------------
        # Non-shared genes
        # ----------------------------------------------------

        geom_label(
            data =
                nonshared_genes,

            aes(
                x = x,
                y = y,
                label = label
            ),

            fill =
                "white",

            colour =
                "#202020",

            size =
                gene_size,

            linewidth =
                0.34,

            label.padding =
                unit(
                    0.085,
                    "lines"
                )
        ) +


        # ----------------------------------------------------
        # Shared genes
        # ----------------------------------------------------

        geom_label(
            data =
                shared_genes,

            aes(
                x = x,
                y = y,
                label = label
            ),

            fill =
                "#FFF2CC",

            colour =
                "#202020",

            fontface =
                "bold",

            size =
                gene_size,

            linewidth =
                0.40,

            label.padding =
                unit(
                    0.085,
                    "lines"
                )
        ) +


        # ----------------------------------------------------
        # Column headings
        # ----------------------------------------------------

        annotate(
            "text",

            x = c(
                0.25,
                2.10,
                3.65,
                5.15
            ),

            y =
                max_y,

            label = c(
                "Candidate TF",
                "Motif",
                "ATAC peak",
                "Associated gene"
            ),

            fontface =
                "bold",

            colour =
                "#303030",

            size =
                2.50
        ) +


        scale_colour_manual(
            values =
                MODULE_COLORS,

            guide =
                "none"
        ) +


        scale_fill_manual(
            values =
                MODULE_COLORS,

            guide =
                "none"
        ) +


        coord_cartesian(
            xlim = c(
                -0.25,
                5.70
            ),

            ylim = c(
                0.30,
                max_y + 0.12
            ),

            clip =
                "off"
        ) +


        labs(
            title =
                paste0(
                    panel_code,
                    "  ",
                    title
                ),

            subtitle =
                subtitle
        ) +


        theme_void(
            base_family =
                "sans"
        ) +


        theme(
            plot.title =
                element_text(
                    size =
                        9.4,

                    face =
                        "bold",

                    hjust =
                        0,

                    margin =
                        margin(
                            b = 2
                        )
                ),

            plot.subtitle =
                element_text(
                    size =
                        7.2,

                    colour =
                        "#555555",

                    hjust =
                        0,

                    lineheight =
                        1.0,

                    margin =
                        margin(
                            b = 6
                        )
                ),

            plot.margin =
                margin(
                    8,
                    18,
                    8,
                    18
                )
        )
}


# ============================================================
# 12. Draw panels
# ============================================================

pA <- draw_panel(
    A,
    "Muscle: promoter-centered candidate network",
    "KLF9 primary structural + SP-like supporting structural module"
)


pB <- draw_panel(
    B,
    "Spleen: overlapping GC-rich candidate network",
    "STREME1, STREME2 and KLF12 primary integrated modules"
)


pC <- draw_panel(
    C,
    "Liver: supporting GC-rich candidate network",
    "SP5-like family annotation within the supporting q <= 0.05 module"
)


# ============================================================
# 13. Save individual panels
# ============================================================

ggsave(
    file.path(
        OUTDIR,
        "Figure6A_Muscle",
        "Figure6A_Muscle_candidate_network_v3.pdf"
    ),
    pA,
    width = 7.8,
    height = 5.0,
    device = cairo_pdf
)


ggsave(
    file.path(
        OUTDIR,
        "Figure6B_Spleen",
        "Figure6B_Spleen_candidate_network_v3.pdf"
    ),
    pB,
    width = 8.4,
    height = 5.4,
    device = cairo_pdf
)


ggsave(
    file.path(
        OUTDIR,
        "Figure6C_Liver",
        "Figure6C_Liver_candidate_network_v3.pdf"
    ),
    pC,
    width = 8.0,
    height = 3.35,
    device = cairo_pdf
)


# ============================================================
# 14. Publication layout
#
# Spleen gets more horizontal space.
# Liver is compressed vertically.
# ============================================================

top_row <- (
    pA |
    pB
) +
    plot_layout(
        widths = c(
            1.00,
            1.20
        )
    )


Figure6_v2 <- (
    top_row /
    pC
) +
    plot_layout(
        heights = c(
            1.00,
            0.50
        )
    ) +
    plot_annotation(

        title =
            "Integrated candidate TF-motif-peak-gene networks",

        subtitle =
            paste0(
                "Representative frozen Step09C modules; ",
                "links summarize multi-evidence candidate relationships"
            ),

        caption = paste0(
            "Dashed: candidate TF-motif annotation; colored: motif-FIMO-positive peak; ",
            "gray: Strong coupled peak-gene association.\n",
            "Filled squares: promoter/TSS-proximal peaks; ",
            "shaded genes: shared across selected modules."
        ),

        theme =
            theme(
                plot.title =
                    element_text(
                        family =
                            "sans",
                        size =
                            11.5,
                        face =
                            "bold",
                        hjust =
                            0
                    ),

                plot.subtitle =
                    element_text(
                        family =
                            "sans",
                        size =
                            7.7,
                        colour =
                            "#4D4D4D",
                        hjust =
                            0,
                        margin =
                            margin(
                                b = 5
                            )
                    ),

                plot.caption =
                    element_text(
                        family =
                            "sans",
                        size =
                            6.3,
                        colour =
                            "#555555",
                        hjust =
                            0,
                        lineheight =
                            1.0,
                        margin =
                            margin(
                                t = 6
                            )
                    ),

                plot.margin =
                    margin(
                        8,
                        10,
                        8,
                        10
                    )
            )
    )


# ============================================================
# 15. Save FINAL
# ============================================================

FINAL_PDF <- file.path(
    OUTDIR,
    "Figure6_FINAL.pdf"
)

FINAL_TIFF <- file.path(
    OUTDIR,
    "Figure6_FINAL_600dpi.tiff"
)


ggsave(
    FINAL_PDF,
    Figure6_v2,
    width = 8.5,
    height = 7.0,
    units = "in",
    device = cairo_pdf
)


ggsave(
    FINAL_TIFF,
    Figure6_v2,
    width = 8.5,
    height = 7.0,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw"
)


# ============================================================
# 16. Export visualization tables
#
# IMPORTANT:
# Replace literal newline characters in label fields with "\\n"
# before TSV writing so they can safely be read again later.
# ============================================================

nodes_export <- nodes_all %>%
    mutate(
        label =
            str_replace_all(
                label,
                "\n",
                "\\\\n"
            )
    )


write_tsv(
    nodes_export,
    file.path(
        OUTDIR,
        "10E_v3_Figure6_plot_nodes.tsv"
    )
)


write_tsv(
    edges_all,
    file.path(
        OUTDIR,
        "10E_v3_Figure6_plot_edges.tsv"
    )
)


write_tsv(
    summary_all,
    file.path(
        OUTDIR,
        "10E_v3_panel_summary.tsv"
    )
)


write_tsv(
    module_n_qc,
    file.path(
        OUTDIR,
        "10E_v3_module_count_QC.tsv"
    )
)


write_tsv(
    gene_set_qc,
    file.path(
        OUTDIR,
        "10E_v3_gene_set_QC.tsv"
    )
)


write_tsv(
    tf_qc,
    file.path(
        OUTDIR,
        "10E_v3_candidate_TF_QC.tsv"
    )
)


# ============================================================
# 17. Final QC summary
# ============================================================

qc <- tribble(

    ~Metric,
    ~Value,
    ~Expected,

    "Authoritative_Step09C_Top5_rows",
    nrow(top5),
    30,

    "Figure6_selected_rows",
    nrow(selected),
    30,

    "Figure6_selected_modules",
    n_distinct(selected$Module_ID),
    6,

    "Figure6_total_plot_nodes",
    nrow(nodes_all),
    52,

    "Figure6_total_plot_edges",
    nrow(edges_all),
    56,

    "Figure6A_nodes",
    nrow(A$nodes),
    18,

    "Figure6A_edges",
    nrow(A$edges),
    19,

    "Figure6B_nodes",
    nrow(B$nodes),
    22,

    "Figure6B_edges",
    nrow(B$edges),
    26,

    "Figure6C_nodes",
    nrow(C$nodes),
    12,

    "Figure6C_edges",
    nrow(C$edges),
    11,

    "Figure6A_unique_genes",
    A$summary$Unique_selected_genes,
    7,

    "Figure6A_shared_genes",
    A$summary$Shared_gene_N,
    3,

    "Figure6B_unique_genes",
    B$summary$Unique_selected_genes,
    8,

    "Figure6B_shared_genes",
    B$summary$Shared_gene_N,
    4,

    "Figure6C_unique_genes",
    C$summary$Unique_selected_genes,
    5,

    "Figure6C_shared_genes",
    C$summary$Shared_gene_N,
    0,

    "Frozen_gene_set_failures",
    sum(
        gene_set_qc$Status !=
            "PASS"
    ),
    0,

    "Frozen_module_role_failures",
    sum(
        role_qc$Status !=
            "PASS"
    ),
    0,

    "Frozen_candidate_TF_failures",
    sum(
        tf_qc$Status !=
            "PASS"
    ),
    0
) %>%

    mutate(
        Difference =
            as.numeric(Value) -
            as.numeric(Expected),

        Status =
            if_else(
                as.numeric(Value) ==
                    as.numeric(Expected),
                "PASS",
                "FAIL"
            )
    )


if (any(qc$Status != "PASS")) {

    print(
        qc,
        n = Inf,
        width = Inf
    )

    stop(
        "Step10E-v2 final QC failed."
    )
}


write_tsv(
    qc,
    file.path(
        OUTDIR,
        "10E_v3_QC_summary.tsv"
    )
)


# ============================================================
# 18. Output file QC
# ============================================================

output_files <- tibble(
    File = c(
        FINAL_PDF,
        FINAL_TIFF,

        file.path(
            OUTDIR,
            "Figure6A_Muscle",
            "Figure6A_Muscle_candidate_network_v3.pdf"
        ),

        file.path(
            OUTDIR,
            "Figure6B_Spleen",
            "Figure6B_Spleen_candidate_network_v3.pdf"
        ),

        file.path(
            OUTDIR,
            "Figure6C_Liver",
            "Figure6C_Liver_candidate_network_v3.pdf"
        )
    )
) %>%

    mutate(
        Exists =
            file.exists(
                File
            ),

        Size_bytes =
            file.info(
                File
            )$size
    )


if (
    any(
        !output_files$Exists
    ) ||
    any(
        is.na(
            output_files$Size_bytes
        )
    ) ||
    any(
        output_files$Size_bytes <= 0
    )
) {

    print(output_files)

    stop(
        "Output-file QC failed."
    )
}


write_tsv(
    output_files,
    file.path(
        OUTDIR,
        "10E_v3_output_file_QC.tsv"
    )
)


# ============================================================
# 19. Notes
# ============================================================

notes <- c(

    "STEP10E-v3 FINAL PUBLICATION POLISH",

    "",

    "Authoritative input:",
    "09C5_top5_primary_supporting_figure_candidates.tsv",

    "",

    "No Step08 analysis was rerun.",
    "No Step09 biological analysis was rerun.",
    "No motif significance was recalculated.",
    "No FIMO threshold was changed.",
    "No RNA-ATAC correlation was recalculated.",
    "No candidate peak, gene, motif, TF, module, node or edge was added.",

    "",

    "Visualization reconstruction:",
    "52 plot nodes and 56 plot edges were reconstructed directly from the same frozen 30 Step09C representative rows.",

    "",

    "Publication polish:",
    "Spleen modules were assigned separate vertical levels.",
    "Spleen receives more horizontal space.",
    "Liver receives less vertical space.",
    "The main caption was shortened.",

    "",

    "Interpretation:",
    "Candidate TF-motif links remain candidate annotations.",
    "Motif-peak links remain FIMO sequence matches.",
    "Peak-gene links remain Strong coupled UROPA-associated candidates.",
    "No edge represents validated TF occupancy or causal regulation."
)


writeLines(
    notes,
    file.path(
        OUTDIR,
        "10E_v3_publication_polish_notes.txt"
    )
)


# ============================================================
# 20. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10E-v3 FINAL COMPLETED\n")
cat("============================================================\n\n")


cat("Frozen source:\n")
cat("Step09C Top5 rows   : ", nrow(top5), "\n", sep = "")
cat("Selected modules    : ", n_distinct(selected$Module_ID), "\n", sep = "")
cat("Selected rows       : ", nrow(selected), "\n", sep = "")
cat("Rebuilt plot nodes  : ", nrow(nodes_all), "\n", sep = "")
cat("Rebuilt plot edges  : ", nrow(edges_all), "\n", sep = "")


cat("\nPanel summary:\n")

print(
    summary_all,
    n = Inf,
    width = Inf
)


cat("\nFrozen gene-set QC:\n")

print(
    gene_set_qc %>%
        select(
            Module_ID,
            Observed_N,
            Expected_N,
            Status
        ),
    n = Inf,
    width = Inf
)


cat("\nCandidate TF QC:\n")

print(
    tf_qc,
    n = Inf,
    width = Inf
)


cat("\nFinal QC:\n")

print(
    qc,
    n = Inf,
    width = Inf
)


cat("\nFinal files:\n")

print(
    output_files,
    n = Inf,
    width = Inf
)


cat("\nPUBLICATION POLISH ONLY:\n")
cat("No frozen biological result changed.\n")
cat("No analysis threshold changed.\n")
cat("Spleen module positions separated vertically.\n")
cat("Liver panel compressed.\n")
cat("Node labels exported with escaped newlines for safe TSV reuse.\n")

cat("============================================================\n")

