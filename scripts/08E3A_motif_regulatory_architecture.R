#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# STEP 08E3A
#
# Motif-positive vs same-tissue motif-negative strong peaks
#
# Analyses:
#   1. Promoter_TSSproximal enrichment
#   2. Full pair-specific genomic-class enrichment
#   3. True-TSS distance comparison
#   4. <=2 kb absolute TSS proximity enrichment
#
# Evidence levels:
#   A. Reported_p1e-4
#   B. High_confidence_q005
#
# IMPORTANT:
#   Background is motif-specific and tissue-matched.
#   TF-expanded Step08E2D/E is NOT used for inference.
# ============================================================

EPS <- 1e-12

coupling_file <- paste0(
    "05_RNA_ATAC_quantitative_coupling/",
    "05_concordant_HC_peak_gene_correlations.tsv"
)

step07_file <- paste0(
    "07_promoter_first_classification/",
    "07A_all_UROPA_assigned_promoter_first.tsv"
)

motif_file <- paste0(
    "08E2_motif_peak_gene_TF_network/",
    "08E2A_all_reported_motif_peak_gene_links.tsv"
)

outdir <- "08E3A_motif_regulatory_architecture"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# Helper functions
# ============================================================

bool_value <- function(x) {
    tolower(as.character(x)) %in% c(
        "true", "t", "1", "yes"
    )
}

safe_div <- function(x, y) {
    ifelse(
        y == 0,
        NA_real_,
        x / y
    )
}

haldane_or <- function(a, b, c, d) {

    ((a + 0.5) * (d + 0.5)) /
        ((b + 0.5) * (c + 0.5))
}

fisher_binary <- function(
    positive_status,
    negative_status
) {

    a <- sum(positive_status, na.rm = TRUE)
    b <- sum(!positive_status, na.rm = TRUE)

    c <- sum(negative_status, na.rm = TRUE)
    d <- sum(!negative_status, na.rm = TRUE)

    mat <- matrix(
        c(a, b, c, d),
        nrow = 2,
        byrow = TRUE
    )

    ft <- fisher.test(
        mat,
        alternative = "two.sided"
    )

    fisher_or <- unname(ft$estimate)

    if (length(fisher_or) == 0) {
        fisher_or <- NA_real_
    }

    tibble(
        Positive_yes = a,
        Positive_no = b,
        Negative_yes = c,
        Negative_no = d,

        Positive_fraction =
            safe_div(a, a + b),

        Negative_fraction =
            safe_div(c, c + d),

        Fold_enrichment =
            safe_div(
                safe_div(a, a + b),
                safe_div(c, c + d)
            ),

        Fisher_OR = fisher_or,

        Haldane_OR =
            haldane_or(a, b, c, d),

        P_value = ft$p.value
    )
}

significance_label <- function(fdr) {

    case_when(
        is.na(fdr) ~ "",
        fdr < 0.001 ~ "***",
        fdr < 0.01  ~ "**",
        fdr < 0.05  ~ "*",
        TRUE ~ ""
    )
}

# ============================================================
# 1. Read Step05 and reconstruct official 702 strong pairs
# ============================================================

cat("Reading Step05 coupling table...\n")

coupling <- read_tsv(
    coupling_file,
    show_col_types = FALSE
)

required_coupling <- c(
    "SAF_peak_id",
    "BED_peak_id",
    "gene_id",
    "ATAC_Max_tissue",
    "rho_P348",
    "rho_P350",
    "rho_mean"
)

missing <- setdiff(
    required_coupling,
    names(coupling)
)

if (length(missing) > 0) {
    stop(
        "Missing Step05 columns: ",
        paste(missing, collapse = ", ")
    )
}

strong <- coupling %>%
    mutate(
        rho_P348 = as.numeric(rho_P348),
        rho_P350 = as.numeric(rho_P350),
        rho_mean = as.numeric(rho_mean)
    ) %>%
    filter(
        rho_P348 >= 0.5 - EPS,
        rho_P350 >= 0.5 - EPS,
        rho_mean >= 0.7 - EPS
    )

cat(
    "Strong pairs reconstructed:",
    nrow(strong),
    "\n"
)

if (nrow(strong) != 702) {
    stop(
        "ERROR: expected 702 strong pairs, observed ",
        nrow(strong)
    )
}

if (n_distinct(strong$SAF_peak_id) != 702) {
    stop(
        "ERROR: strong SAF_peak_id is not unique."
    )
}

# ============================================================
# 2. Join Step07 pair-specific annotation
# ============================================================

cat("Reading Step07 promoter-first annotation...\n")

step07 <- read_tsv(
    step07_file,
    show_col_types = FALSE
)

required_step07 <- c(
    "SAF_peak_id",
    "gene_id",
    "Pair_specific_class",
    "Peak_to_UROPA_TSS_min_abs_distance_bp"
)

missing <- setdiff(
    required_step07,
    names(step07)
)

if (length(missing) > 0) {
    stop(
        "Missing Step07 columns: ",
        paste(missing, collapse = ", ")
    )
}

step07_keep <- step07 %>%
    select(
        SAF_peak_id,
        gene_id,
        Pair_specific_class,
        Peak_to_UROPA_TSS_min_abs_distance_bp,
        any_of(c(
            "Peak_center_to_UROPA_TSS_signed_bp",
            "UROPA_gene_TSS_1based",
            "Genomic_class"
        ))
    )

strong_annot <- strong %>%
    left_join(
        step07_keep,
        by = c(
            "SAF_peak_id",
            "gene_id"
        )
    ) %>%
    mutate(
        Peak_to_UROPA_TSS_min_abs_distance_bp =
            as.numeric(
                Peak_to_UROPA_TSS_min_abs_distance_bp
            )
    )

if (
    any(is.na(strong_annot$Pair_specific_class))
) {
    stop(
        "ERROR: some strong pairs have no ",
        "Pair_specific_class after Step07 join."
    )
}

if (nrow(strong_annot) != 702) {
    stop(
        "ERROR: Step07 join changed strong-pair row count."
    )
}

cat(
    "PASS: all 702 strong pairs have Step07 annotation.\n"
)

# ============================================================
# 3. Read Step08E2 motif-positive links
# ============================================================

cat("Reading Step08E2 motif-positive links...\n")

motif_links <- read_tsv(
    motif_file,
    show_col_types = FALSE
)

required_motif <- c(
    "Tissue",
    "FIMO_tier",
    "Motif_source",
    "Motif_family_ID",
    "Motif_family_label",
    "Motif_ID",
    "Motif_alt_ID",
    "SAF_peak_id",
    "FIMO_q_le_0.05"
)

missing <- setdiff(
    required_motif,
    names(motif_links)
)

if (length(missing) > 0) {
    stop(
        "Missing Step08E2 columns: ",
        paste(missing, collapse = ", ")
    )
}

motif_links <- motif_links %>%
    mutate(
        FIMO_q_le_0.05 =
            bool_value(FIMO_q_le_0.05)
    )

cat(
    "Reported motif-peak links:",
    nrow(motif_links),
    "\n"
)

if (nrow(motif_links) != 902) {
    stop(
        "ERROR: expected 902 motif-peak links."
    )
}

# ============================================================
# 4. Motif metadata
# ============================================================

motif_meta <- motif_links %>%
    distinct(
        Tissue,
        FIMO_tier,
        Motif_source,
        Motif_family_ID,
        Motif_family_label,
        Motif_ID,
        Motif_alt_ID
    ) %>%
    arrange(
        Tissue,
        FIMO_tier,
        Motif_source,
        Motif_ID
    )

cat(
    "Selected motifs:",
    nrow(motif_meta),
    "\n"
)

# ============================================================
# Fixed pair-specific genomic classes
# ============================================================

class_levels <- c(
    "Promoter_TSSproximal",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Distal_to_associated_gene",
    "Gene_not_in_GTF"
)

observed_classes <- unique(
    strong_annot$Pair_specific_class
)

unexpected_classes <- setdiff(
    observed_classes,
    class_levels
)

if (length(unexpected_classes) > 0) {

    warning(
        "Unexpected Pair_specific_class values: ",
        paste(
            unexpected_classes,
            collapse = ", "
        )
    )

    class_levels <- c(
        class_levels,
        unexpected_classes
    )
}

# ============================================================
# 5. Build motif-specific positive/negative sets
# ============================================================

evidence_levels <- c(
    "Reported_p1e-4",
    "High_confidence_q005"
)

set_rows <- list()
promoter_rows <- list()
class_rows <- list()
tss_rows <- list()
tss2kb_rows <- list()
composition_rows <- list()

idx_set <- 1
idx_prom <- 1
idx_class <- 1
idx_tss <- 1
idx_tss2 <- 1
idx_comp <- 1

for (i in seq_len(nrow(motif_meta))) {

    meta <- motif_meta[i, ]

    tissue <- meta$Tissue
    motif_id <- meta$Motif_ID

    tissue_strong <- strong_annot %>%
        filter(
            ATAC_Max_tissue == tissue
        )

    if (nrow(tissue_strong) == 0) {
        stop(
            "No strong background for tissue ",
            tissue
        )
    }

    motif_hits <- motif_links %>%
        filter(
            Tissue == tissue,
            Motif_ID == motif_id
        )

    for (evidence in evidence_levels) {

        if (evidence == "Reported_p1e-4") {

            pos_ids <- unique(
                motif_hits$SAF_peak_id
            )

        } else {

            pos_ids <- unique(
                motif_hits$SAF_peak_id[
                    motif_hits$FIMO_q_le_0.05
                ]
            )
        }

        positive <- tissue_strong %>%
            filter(
                SAF_peak_id %in% pos_ids
            )

        negative <- tissue_strong %>%
            filter(
                !SAF_peak_id %in% pos_ids
            )

        n_pos <- nrow(positive)
        n_neg <- nrow(negative)

        if (
            n_pos + n_neg !=
            nrow(tissue_strong)
        ) {
            stop(
                "Positive/negative partition failure: ",
                tissue, " ", motif_id
            )
        }

        low_power <- (
            n_pos < 10 ||
            n_neg < 10
        )

        set_rows[[idx_set]] <- tibble(
            Tissue = tissue,
            FIMO_tier = meta$FIMO_tier,
            Motif_source = meta$Motif_source,
            Motif_family_ID =
                meta$Motif_family_ID,
            Motif_family_label =
                meta$Motif_family_label,
            Motif_ID = motif_id,
            Motif_alt_ID =
                meta$Motif_alt_ID,
            Evidence_level = evidence,
            Same_tissue_strong_N =
                nrow(tissue_strong),
            Motif_positive_N = n_pos,
            Motif_negative_N = n_neg,
            Motif_positive_fraction =
                n_pos / nrow(tissue_strong),
            Low_power_Npos_or_Nneg_lt10 =
                low_power
        )

        idx_set <- idx_set + 1

        # ====================================================
        # 5A. Promoter enrichment
        # ====================================================

        if (
            n_pos > 0 &&
            n_neg > 0
        ) {

            prom_test <- fisher_binary(
                positive$Pair_specific_class ==
                    "Promoter_TSSproximal",

                negative$Pair_specific_class ==
                    "Promoter_TSSproximal"
            )

            promoter_rows[[idx_prom]] <- bind_cols(
                tibble(
                    Tissue = tissue,
                    FIMO_tier =
                        meta$FIMO_tier,
                    Motif_source =
                        meta$Motif_source,
                    Motif_family_ID =
                        meta$Motif_family_ID,
                    Motif_family_label =
                        meta$Motif_family_label,
                    Motif_ID = motif_id,
                    Motif_alt_ID =
                        meta$Motif_alt_ID,
                    Evidence_level =
                        evidence,
                    Feature =
                        "Promoter_TSSproximal",
                    Motif_positive_N =
                        n_pos,
                    Motif_negative_N =
                        n_neg,
                    Low_power_Npos_or_Nneg_lt10 =
                        low_power
                ),
                prom_test
            )

            idx_prom <- idx_prom + 1

            # ================================================
            # 5B. TSS <= 2 kb enrichment
            # ================================================

            pos_tss_valid <- !is.na(
                positive$
                    Peak_to_UROPA_TSS_min_abs_distance_bp
            )

            neg_tss_valid <- !is.na(
                negative$
                    Peak_to_UROPA_TSS_min_abs_distance_bp
            )

            tss2_test <- fisher_binary(
                positive$
                    Peak_to_UROPA_TSS_min_abs_distance_bp[
                        pos_tss_valid
                    ] <= 2000,

                negative$
                    Peak_to_UROPA_TSS_min_abs_distance_bp[
                        neg_tss_valid
                    ] <= 2000
            )

            tss2kb_rows[[idx_tss2]] <- bind_cols(
                tibble(
                    Tissue = tissue,
                    FIMO_tier =
                        meta$FIMO_tier,
                    Motif_source =
                        meta$Motif_source,
                    Motif_family_ID =
                        meta$Motif_family_ID,
                    Motif_family_label =
                        meta$Motif_family_label,
                    Motif_ID = motif_id,
                    Motif_alt_ID =
                        meta$Motif_alt_ID,
                    Evidence_level =
                        evidence,
                    Feature =
                        "Absolute_TSS_distance_le_2kb",
                    Motif_positive_N =
                        n_pos,
                    Motif_negative_N =
                        n_neg,
                    Positive_valid_TSS_N =
                        sum(pos_tss_valid),
                    Negative_valid_TSS_N =
                        sum(neg_tss_valid),
                    Low_power_Npos_or_Nneg_lt10 =
                        low_power
                ),
                tss2_test
            )

            idx_tss2 <- idx_tss2 + 1

            # ================================================
            # 5C. Continuous absolute TSS distance
            # ================================================

            pos_dist <- positive$
                Peak_to_UROPA_TSS_min_abs_distance_bp

            neg_dist <- negative$
                Peak_to_UROPA_TSS_min_abs_distance_bp

            pos_dist <- pos_dist[
                !is.na(pos_dist)
            ]

            neg_dist <- neg_dist[
                !is.na(neg_dist)
            ]

            wilcox_p <- NA_real_
            wilcox_W <- NA_real_

            if (
                length(pos_dist) > 0 &&
                length(neg_dist) > 0
            ) {

                wt <- suppressWarnings(
                    wilcox.test(
                        pos_dist,
                        neg_dist,
                        alternative = "two.sided",
                        exact = FALSE
                    )
                )

                wilcox_p <- wt$p.value
                wilcox_W <- unname(
                    wt$statistic
                )
            }

            tss_rows[[idx_tss]] <- tibble(
                Tissue = tissue,
                FIMO_tier =
                    meta$FIMO_tier,
                Motif_source =
                    meta$Motif_source,
                Motif_family_ID =
                    meta$Motif_family_ID,
                Motif_family_label =
                    meta$Motif_family_label,
                Motif_ID = motif_id,
                Motif_alt_ID =
                    meta$Motif_alt_ID,
                Evidence_level =
                    evidence,

                Positive_valid_TSS_N =
                    length(pos_dist),

                Negative_valid_TSS_N =
                    length(neg_dist),

                Positive_median_abs_TSS_bp =
                    ifelse(
                        length(pos_dist) > 0,
                        median(pos_dist),
                        NA_real_
                    ),

                Negative_median_abs_TSS_bp =
                    ifelse(
                        length(neg_dist) > 0,
                        median(neg_dist),
                        NA_real_
                    ),

                Positive_Q25_abs_TSS_bp =
                    ifelse(
                        length(pos_dist) > 0,
                        quantile(
                            pos_dist,
                            0.25,
                            names = FALSE
                        ),
                        NA_real_
                    ),

                Positive_Q75_abs_TSS_bp =
                    ifelse(
                        length(pos_dist) > 0,
                        quantile(
                            pos_dist,
                            0.75,
                            names = FALSE
                        ),
                        NA_real_
                    ),

                Negative_Q25_abs_TSS_bp =
                    ifelse(
                        length(neg_dist) > 0,
                        quantile(
                            neg_dist,
                            0.25,
                            names = FALSE
                        ),
                        NA_real_
                    ),

                Negative_Q75_abs_TSS_bp =
                    ifelse(
                        length(neg_dist) > 0,
                        quantile(
                            neg_dist,
                            0.75,
                            names = FALSE
                        ),
                        NA_real_
                    ),

                Median_distance_ratio_Pos_over_Neg =
                    ifelse(
                        length(pos_dist) > 0 &&
                        length(neg_dist) > 0 &&
                        median(neg_dist) > 0,

                        median(pos_dist) /
                            median(neg_dist),

                        NA_real_
                    ),

                Wilcoxon_W =
                    wilcox_W,

                P_value =
                    wilcox_p,

                Low_power_Npos_or_Nneg_lt10 =
                    low_power
            )

            idx_tss <- idx_tss + 1

            # ================================================
            # 5D. Full pair-specific class enrichment
            # ================================================

            for (class_name in class_levels) {

                class_test <- fisher_binary(
                    positive$Pair_specific_class ==
                        class_name,

                    negative$Pair_specific_class ==
                        class_name
                )

                class_rows[[idx_class]] <- bind_cols(
                    tibble(
                        Tissue = tissue,
                        FIMO_tier =
                            meta$FIMO_tier,
                        Motif_source =
                            meta$Motif_source,
                        Motif_family_ID =
                            meta$Motif_family_ID,
                        Motif_family_label =
                            meta$Motif_family_label,
                        Motif_ID =
                            motif_id,
                        Motif_alt_ID =
                            meta$Motif_alt_ID,
                        Evidence_level =
                            evidence,
                        Pair_specific_class =
                            class_name,
                        Motif_positive_N =
                            n_pos,
                        Motif_negative_N =
                            n_neg,
                        Low_power_Npos_or_Nneg_lt10 =
                            low_power
                    ),
                    class_test
                )

                idx_class <- idx_class + 1
            }
        }

        # ====================================================
        # 5E. Composition table
        # ====================================================

        for (
            set_name in c(
                "Motif_positive",
                "Motif_negative"
            )
        ) {

            current <- if (
                set_name == "Motif_positive"
            ) {
                positive
            } else {
                negative
            }

            denom <- nrow(current)

            for (class_name in class_levels) {

                nn <- sum(
                    current$Pair_specific_class ==
                        class_name,
                    na.rm = TRUE
                )

                composition_rows[[idx_comp]] <- tibble(
                    Tissue = tissue,
                    FIMO_tier =
                        meta$FIMO_tier,
                    Motif_source =
                        meta$Motif_source,
                    Motif_family_ID =
                        meta$Motif_family_ID,
                    Motif_family_label =
                        meta$Motif_family_label,
                    Motif_ID = motif_id,
                    Motif_alt_ID =
                        meta$Motif_alt_ID,
                    Evidence_level =
                        evidence,
                    Set = set_name,
                    Pair_specific_class =
                        class_name,
                    N = nn,
                    Total_N = denom,
                    Fraction =
                        ifelse(
                            denom > 0,
                            nn / denom,
                            NA_real_
                        )
                )

                idx_comp <- idx_comp + 1
            }
        }
    }
}

# ============================================================
# 6. Combine results
# ============================================================

set_sizes <- bind_rows(set_rows)

promoter_enrichment <-
    bind_rows(promoter_rows)

class_enrichment <-
    bind_rows(class_rows)

tss_distance <-
    bind_rows(tss_rows)

tss2kb_enrichment <-
    bind_rows(tss2kb_rows)

class_composition <-
    bind_rows(composition_rows)

# ============================================================
# 7. BH adjustment
#
# FDR is adjusted separately by evidence level and test family.
# ============================================================

if (nrow(promoter_enrichment) > 0) {

    promoter_enrichment <-
        promoter_enrichment %>%
        group_by(Evidence_level) %>%
        mutate(
            FDR =
                p.adjust(
                    P_value,
                    method = "BH"
                ),
            Significance =
                significance_label(FDR),
            log2_Haldane_OR =
                log2(Haldane_OR)
        ) %>%
        ungroup()
}

if (nrow(class_enrichment) > 0) {

    class_enrichment <-
        class_enrichment %>%
        group_by(Evidence_level) %>%
        mutate(
            FDR =
                p.adjust(
                    P_value,
                    method = "BH"
                ),
            Significance =
                significance_label(FDR),
            log2_Haldane_OR =
                log2(Haldane_OR)
        ) %>%
        ungroup()
}

if (nrow(tss_distance) > 0) {

    tss_distance <-
        tss_distance %>%
        group_by(Evidence_level) %>%
        mutate(
            FDR =
                p.adjust(
                    P_value,
                    method = "BH"
                ),
            Significance =
                significance_label(FDR)
        ) %>%
        ungroup()
}

if (nrow(tss2kb_enrichment) > 0) {

    tss2kb_enrichment <-
        tss2kb_enrichment %>%
        group_by(Evidence_level) %>%
        mutate(
            FDR =
                p.adjust(
                    P_value,
                    method = "BH"
                ),
            Significance =
                significance_label(FDR),
            log2_Haldane_OR =
                log2(Haldane_OR)
        ) %>%
        ungroup()
}

# ============================================================
# 8. Write tables
# ============================================================

write_tsv(
    set_sizes,
    file.path(
        outdir,
        "08E3A1_motif_positive_negative_set_sizes.tsv"
    )
)

write_tsv(
    promoter_enrichment,
    file.path(
        outdir,
        "08E3A2_promoter_enrichment.tsv"
    )
)

write_tsv(
    class_enrichment,
    file.path(
        outdir,
        "08E3A3_pair_specific_class_enrichment.tsv"
    )
)

write_tsv(
    tss_distance,
    file.path(
        outdir,
        "08E3A4_TSS_distance_comparison.tsv"
    )
)

write_tsv(
    tss2kb_enrichment,
    file.path(
        outdir,
        "08E3A5_TSS_2kb_enrichment.tsv"
    )
)

write_tsv(
    class_composition,
    file.path(
        outdir,
        "08E3A6_pair_specific_class_composition.tsv"
    )
)

# ============================================================
# 9. Build plotting data
# ============================================================

motif_labels <- motif_meta %>%
    transmute(
        Tissue,
        Motif_ID,
        Plot_label = paste0(
            Tissue,
            " | ",
            Motif_alt_ID,
            "\n",
            Motif_ID
        )
    )

# ============================================================
# 10. Promoter enrichment forest plot
# ============================================================

if (nrow(promoter_enrichment) > 0) {

    plot_prom <- promoter_enrichment %>%
        left_join(
            motif_labels,
            by = c(
                "Tissue",
                "Motif_ID"
            )
        ) %>%
        mutate(
            Plot_label =
                factor(
                    Plot_label,
                    levels = rev(
                        unique(
                            Plot_label
                        )
                    )
                )
        )

    p <- ggplot(
        plot_prom,
        aes(
            x = log2_Haldane_OR,
            y = Plot_label
        )
    ) +
        geom_vline(
            xintercept = 0,
            linetype = 2
        ) +
        geom_point(
            aes(
                shape =
                    FDR < 0.05
            ),
            size = 2.5
        ) +
        facet_wrap(
            ~ Evidence_level,
            scales = "free_y"
        ) +
        labs(
            x = "log2 odds ratio: promoter enrichment",
            y = NULL,
            shape = "FDR < 0.05"
        ) +
        theme_bw(base_size = 11)

    ggsave(
        file.path(
            outdir,
            "08E3A7_promoter_enrichment_forest.pdf"
        ),
        p,
        width = 10,
        height = 7
    )
}

# ============================================================
# 11. Genomic-class enrichment heatmap
# ============================================================

if (nrow(class_enrichment) > 0) {

    plot_class <- class_enrichment %>%
        left_join(
            motif_labels,
            by = c(
                "Tissue",
                "Motif_ID"
            )
        ) %>%
        mutate(
            Pair_specific_class =
                factor(
                    Pair_specific_class,
                    levels =
                        class_levels
                )
        )

    p <- ggplot(
        plot_class,
        aes(
            x = Pair_specific_class,
            y = Plot_label,
            fill = log2_Haldane_OR
        )
    ) +
        geom_tile() +
        geom_text(
            aes(
                label = Significance
            ),
            size = 3
        ) +
        facet_wrap(
            ~ Evidence_level,
            scales = "free_y"
        ) +
        labs(
            x = "Pair-specific genomic class",
            y = NULL,
            fill = "log2 OR"
        ) +
        theme_bw(base_size = 10) +
        theme(
            axis.text.x =
                element_text(
                    angle = 45,
                    hjust = 1
                )
        )

    ggsave(
        file.path(
            outdir,
            "08E3A8_genomic_class_enrichment_heatmap.pdf"
        ),
        p,
        width = 12,
        height = 8
    )
}

# ============================================================
# 12. Pair-specific class composition plot
# ============================================================

if (nrow(class_composition) > 0) {

    plot_comp <- class_composition %>%
        filter(
            Evidence_level ==
                "Reported_p1e-4"
        ) %>%
        left_join(
            motif_labels,
            by = c(
                "Tissue",
                "Motif_ID"
            )
        ) %>%
        mutate(
            Pair_specific_class =
                factor(
                    Pair_specific_class,
                    levels = class_levels
                )
        )

    p <- ggplot(
        plot_comp,
        aes(
            x = Set,
            y = Fraction,
            fill = Pair_specific_class
        )
    ) +
        geom_col(
            position = "fill"
        ) +
        facet_wrap(
            ~ Plot_label,
            scales = "free_x"
        ) +
        labs(
            x = NULL,
            y = "Proportion",
            fill = "Pair-specific class"
        ) +
        theme_bw(base_size = 9) +
        theme(
            axis.text.x =
                element_text(
                    angle = 45,
                    hjust = 1
                )
        )

    ggsave(
        file.path(
            outdir,
            "08E3A9_genomic_class_composition.pdf"
        ),
        p,
        width = 12,
        height = 10
    )
}

# ============================================================
# 13. TSS-distance visualization
#
# Reconstruct long plotting data.
# Only Reported_p1e-4 is shown to avoid duplicated facets.
# ============================================================

tss_plot_list <- list()
idx <- 1

for (i in seq_len(nrow(motif_meta))) {

    meta <- motif_meta[i, ]

    tissue <- meta$Tissue
    motif_id <- meta$Motif_ID

    tissue_strong <- strong_annot %>%
        filter(
            ATAC_Max_tissue ==
                tissue
        )

    pos_ids <- motif_links %>%
        filter(
            Tissue == tissue,
            Motif_ID == motif_id
        ) %>%
        pull(SAF_peak_id) %>%
        unique()

    tmp <- tissue_strong %>%
        transmute(
            Tissue = tissue,
            Motif_ID = motif_id,
            Motif_alt_ID =
                meta$Motif_alt_ID,
            Set = ifelse(
                SAF_peak_id %in% pos_ids,
                "Motif_positive",
                "Motif_negative"
            ),
            Abs_TSS_distance_bp =
                Peak_to_UROPA_TSS_min_abs_distance_bp
        ) %>%
        filter(
            !is.na(
                Abs_TSS_distance_bp
            )
        )

    tss_plot_list[[idx]] <- tmp
    idx <- idx + 1
}

tss_plot_data <- bind_rows(
    tss_plot_list
) %>%
    left_join(
        motif_labels,
        by = c(
            "Tissue",
            "Motif_ID"
        )
    ) %>%
    mutate(
        log10_TSS_distance_plus1 =
            log10(
                Abs_TSS_distance_bp + 1
            )
    )

write_tsv(
    tss_plot_data,
    file.path(
        outdir,
        "08E3A10_TSS_distance_plot_data.tsv"
    )
)

p <- ggplot(
    tss_plot_data,
    aes(
        x = Set,
        y = log10_TSS_distance_plus1
    )
) +
    geom_boxplot(
        outlier.shape = NA
    ) +
    facet_wrap(
        ~ Plot_label,
        scales = "free_y"
    ) +
    labs(
        x = NULL,
        y = "log10(|peak-to-TSS distance| + 1 bp)"
    ) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.x =
            element_text(
                angle = 45,
                hjust = 1
            )
    )

ggsave(
    file.path(
        outdir,
        "08E3A11_TSS_distance_boxplot.pdf"
    ),
    p,
    width = 12,
    height = 10
)

# ============================================================
# 14. QC and summary
# ============================================================

reported_sets <- set_sizes %>%
    filter(
        Evidence_level ==
            "Reported_p1e-4"
    )

q_sets <- set_sizes %>%
    filter(
        Evidence_level ==
            "High_confidence_q005"
    )

qc <- tibble(
    Metric = c(
        "Official_strong_pairs",
        "Strong_pairs_with_Step07_annotation",
        "Selected_motifs",
        "Reported_motif_peak_links_input",
        "Reported_set_definitions",
        "Q005_set_definitions",
        "Reported_positive_peak_sum",
        "Q005_positive_peak_sum",
        "Promoter_tests_reported",
        "Promoter_tests_q005",
        "Pair_class_tests_reported",
        "TSS_distance_tests_reported"
    ),

    Value = c(
        nrow(strong),
        nrow(strong_annot),
        nrow(motif_meta),
        nrow(motif_links),
        nrow(reported_sets),
        nrow(q_sets),
        sum(
            reported_sets$
                Motif_positive_N
        ),
        sum(
            q_sets$
                Motif_positive_N
        ),
        sum(
            promoter_enrichment$
                Evidence_level ==
                "Reported_p1e-4"
        ),
        sum(
            promoter_enrichment$
                Evidence_level ==
                "High_confidence_q005"
        ),
        sum(
            class_enrichment$
                Evidence_level ==
                "Reported_p1e-4"
        ),
        sum(
            tss_distance$
                Evidence_level ==
                "Reported_p1e-4"
        )
    ),

    Expected_or_note = c(
        "702",
        "702",
        "11",
        "902",
        "11",
        "11",
        "902",
        "612",
        "11",
        "motifs with q-positive peaks",
        "11 x genomic classes",
        "11"
    )
)

write_tsv(
    qc,
    file.path(
        outdir,
        "08E3A12_QC_summary.tsv"
    )
)

# ============================================================
# 15. Console result
# ============================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("STEP 08E3A COMPLETED\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\nInput/QC:\n")
print(qc, n = Inf)

cat("\nPromoter enrichment: Reported p<=1e-4\n")

print(
    promoter_enrichment %>%
        filter(
            Evidence_level ==
                "Reported_p1e-4"
        ) %>%
        select(
            Tissue,
            FIMO_tier,
            Motif_ID,
            Motif_alt_ID,
            Motif_positive_N,
            Positive_fraction,
            Negative_fraction,
            Fold_enrichment,
            Haldane_OR,
            P_value,
            FDR,
            Significance
        ) %>%
        arrange(FDR),
    n = Inf
)

cat("\nTSS distance comparison: Reported p<=1e-4\n")

print(
    tss_distance %>%
        filter(
            Evidence_level ==
                "Reported_p1e-4"
        ) %>%
        select(
            Tissue,
            Motif_ID,
            Motif_alt_ID,
            Positive_median_abs_TSS_bp,
            Negative_median_abs_TSS_bp,
            Median_distance_ratio_Pos_over_Neg,
            P_value,
            FDR,
            Significance
        ) %>%
        arrange(FDR),
    n = Inf
)

cat("\nSignificant genomic-class enrichments/depletions FDR<0.05:\n")

sig_classes <- class_enrichment %>%
    filter(
        FDR < 0.05
    ) %>%
    select(
        Tissue,
        Evidence_level,
        Motif_ID,
        Motif_alt_ID,
        Pair_specific_class,
        Positive_fraction,
        Negative_fraction,
        Fold_enrichment,
        Haldane_OR,
        P_value,
        FDR,
        Significance
    ) %>%
    arrange(
        Evidence_level,
        FDR
    )

if (nrow(sig_classes) == 0) {

    cat(
        "No pair-specific class test reached FDR < 0.05.\n"
    )

} else {

    print(
        sig_classes,
        n = Inf
    )
}

cat("\nMain outputs:\n")

cat(
    file.path(
        outdir,
        "08E3A1_motif_positive_negative_set_sizes.tsv"
    ),
    "\n"
)

cat(
    file.path(
        outdir,
        "08E3A2_promoter_enrichment.tsv"
    ),
    "\n"
)

cat(
    file.path(
        outdir,
        "08E3A3_pair_specific_class_enrichment.tsv"
    ),
    "\n"
)

cat(
    file.path(
        outdir,
        "08E3A4_TSS_distance_comparison.tsv"
    ),
    "\n"
)

cat(
    file.path(
        outdir,
        "08E3A5_TSS_2kb_enrichment.tsv"
    ),
    "\n"
)

cat(
    file.path(
        outdir,
        "08E3A12_QC_summary.tsv"
    ),
    "\n"
)

cat(paste(rep("=", 80), collapse = ""), "\n")

