suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
})

cat("\n")
cat("============================================================\n")
cat("STEP10F3B5A — MOTIF-TF RNA SUPPORT SOURCE LOCK\n")
cat("Frozen Step08C/08D/08E2/09C results only\n")
cat("Inspection and category audit only\n")
cat("No candidate TF reselection\n")
cat("============================================================\n\n")


# ============================================================
# 1. Authoritative frozen sources
# ============================================================

F_MASTER <- paste0(
    "08D_motif_family_TF_master/",
    "08D2_publication_TF_candidate_master.tsv"
)

F_PRIORITY <- paste0(
    "08C_TF_RNA_integration/",
    "08C5_priority_TF_candidates.tsv"
)

F_TPM <- paste0(
    "08C_TF_RNA_integration/",
    "08C6_candidate_TF_tissue_mean_TPM.tsv"
)

F_NETWORK <- paste0(
    "08E2_motif_peak_gene_TF_network/",
    "08E2E_high_RNA_TF_motif_peak_gene_master.tsv"
)

F_MODULE <- paste0(
    "09C_publication_candidate_modules/",
    "09C6_module_candidate_summary.tsv"
)


required_files <- c(
    F_MASTER,
    F_PRIORITY,
    F_TPM,
    F_NETWORK,
    F_MODULE
)


missing_files <- required_files[
    !file.exists(required_files)
]


if (length(missing_files) > 0) {

    stop(
        "Missing frozen source file(s):\n",
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
    "FigureS5_motif_TF_RNA_support_source_lock"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. Read frozen tables
# ============================================================

master <- read_tsv(
    F_MASTER,
    show_col_types = FALSE
)

priority <- read_tsv(
    F_PRIORITY,
    show_col_types = FALSE
)

tpm <- read_tsv(
    F_TPM,
    show_col_types = FALSE
)

network <- read_tsv(
    F_NETWORK,
    show_col_types = FALSE
)

module <- read_tsv(
    F_MODULE,
    show_col_types = FALSE
)


# ============================================================
# 4. Required columns in authoritative master
# ============================================================

required_master_cols <- c(
    "Tissue",
    "TF_symbol",
    "gene_id",
    "Evidence_sources",
    "AME_exploratory_support",
    "Default_Tomtom_support",
    "NoHoldout_Tomtom_support",
    "CrossRun_reproducible_Tomtom_support",
    "Best_AME_p",
    "Best_AME_adj_p",
    "Best_AME_E",
    "Best_Tomtom_q",
    "Target_tissue_mean_TPM",
    "Target_tissue_P348_TPM",
    "Target_tissue_P350_TPM",
    "High_RNA_expression_support",
    "Tau",
    "Target_tissue_SPM",
    "Target_tissue_expression_rank",
    "Max_expression_tissue",
    "Target_is_max_expression_tissue",
    "SPM_ge_0.5",
    "RNA_support_role",
    "FIMO_best_tier",
    "Interpretation"
)


missing_master_cols <- setdiff(
    required_master_cols,
    colnames(master)
)


if (length(missing_master_cols) > 0) {

    stop(
        "08D2 master missing required column(s): ",
        paste(
            missing_master_cols,
            collapse = ", "
        )
    )
}


# ============================================================
# 5. Frozen-table row-count QC
# ============================================================

row_qc <- tibble(
    Source = c(
        "08D2_publication_TF_candidate_master",
        "08C5_priority_TF_candidates",
        "08C6_candidate_TF_tissue_mean_TPM",
        "08E2E_high_RNA_TF_motif_peak_gene_master",
        "09C6_module_candidate_summary"
    ),

    Observed_rows = c(
        nrow(master),
        nrow(priority),
        nrow(tpm),
        nrow(network),
        nrow(module)
    ),

    Expected_rows = c(
        48,
        41,
        36,
        3178,
        11
    )
) %>%

    mutate(
        Status = if_else(
            Observed_rows ==
                Expected_rows,

            "PASS",
            "FAIL"
        )
    )


write_tsv(
    row_qc,
    file.path(
        OUTDIR,
        "10F3B5A_frozen_source_row_QC.tsv"
    )
)


if (any(row_qc$Status != "PASS")) {

    print(
        row_qc,
        width = Inf
    )

    stop(
        "Frozen source row-count QC failed. ",
        "Do not build Supplementary Figure S5 yet."
    )
}


# ============================================================
# 6. Tissue audit
# ============================================================

tissue_order <- c(
    "Cerebellum",
    "Liver",
    "Muscle",
    "Spleen"
)


master_tissue_counts <- master %>%

    count(
        Tissue,
        name = "Candidate_TF_N"
    ) %>%

    arrange(
        match(
            Tissue,
            tissue_order
        )
    )


write_tsv(
    master_tissue_counts,
    file.path(
        OUTDIR,
        "10F3B5A_candidate_TF_counts_by_tissue.tsv"
    )
)


unexpected_tissues <- setdiff(
    unique(
        master$Tissue
    ),
    tissue_order
)


tissue_qc <- tibble(
    Metric = c(
        "Authoritative_master_rows",
        "Distinct_target_tissues",
        "Unexpected_target_tissues"
    ),

    Value = c(
        nrow(master),
        n_distinct(
            master$Tissue
        ),
        length(
            unexpected_tissues
        )
    ),

    Expected = c(
        48,
        4,
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
    tissue_qc,
    file.path(
        OUTDIR,
        "10F3B5A_tissue_QC.tsv"
    )
)


# ============================================================
# 7. Key uniqueness audit
# ============================================================

duplicate_key <- master %>%

    count(
        Tissue,
        TF_symbol,
        name = "N"
    ) %>%

    filter(
        N > 1
    )


write_tsv(
    duplicate_key,
    file.path(
        OUTDIR,
        "10F3B5A_duplicate_Tissue_TF_keys.tsv"
    )
)


# Do not fail automatically:
# repeated rows could reflect source-defined evidence structure.
# We only audit them before plotting.


# ============================================================
# 8. Numeric completeness
# ============================================================

numeric_qc <- tibble(
    Metric = c(
        "Target_tissue_mean_TPM_finite",
        "Target_tissue_P348_TPM_finite",
        "Target_tissue_P350_TPM_finite",
        "Tau_finite",
        "Target_tissue_SPM_finite",
        "Expression_rank_finite"
    ),

    Finite_N = c(
        sum(
            is.finite(
                master$Target_tissue_mean_TPM
            )
        ),

        sum(
            is.finite(
                master$Target_tissue_P348_TPM
            )
        ),

        sum(
            is.finite(
                master$Target_tissue_P350_TPM
            )
        ),

        sum(
            is.finite(
                master$Tau
            )
        ),

        sum(
            is.finite(
                master$Target_tissue_SPM
            )
        ),

        sum(
            is.finite(
                master$Target_tissue_expression_rank
            )
        )
    ),

    Total_N = nrow(
        master
    )
) %>%

    mutate(
        Missing_N =
            Total_N -
            Finite_N,

        Status =
            if_else(
                Missing_N == 0,
                "PASS",
                "CHECK"
            )
    )


write_tsv(
    numeric_qc,
    file.path(
        OUTDIR,
        "10F3B5A_numeric_completeness_QC.tsv"
    )
)


# ============================================================
# 9. Robust binary parser for audit only
# ============================================================

as_flag <- function(x) {

    z <- tolower(
        trimws(
            as.character(
                x
            )
        )
    )

    case_when(
        z %in% c(
            "true",
            "t",
            "1",
            "yes",
            "y",
            "pass"
        ) ~ TRUE,

        z %in% c(
            "false",
            "f",
            "0",
            "no",
            "n",
            "fail"
        ) ~ FALSE,

        TRUE ~ NA
    )
}


flag_columns <- c(
    "AME_exploratory_support",
    "Default_Tomtom_support",
    "NoHoldout_Tomtom_support",
    "CrossRun_reproducible_Tomtom_support",
    "High_RNA_expression_support",
    "Target_is_max_expression_tissue",
    "SPM_ge_0.5"
)


flag_long <- master %>%

    select(
        Tissue,
        TF_symbol,
        all_of(
            flag_columns
        )
    ) %>%

    pivot_longer(
        cols = all_of(
            flag_columns
        ),

        names_to = "Evidence",
        values_to = "Raw_value"
    ) %>%

    mutate(
        Parsed_flag =
            as_flag(
                Raw_value
            )
    )


flag_value_audit <- flag_long %>%

    count(
        Evidence,
        Raw_value,
        Parsed_flag,
        name = "N"
    ) %>%

    arrange(
        Evidence,
        desc(N)
    )


write_tsv(
    flag_value_audit,
    file.path(
        OUTDIR,
        "10F3B5A_boolean_field_value_audit.tsv"
    )
)


flag_summary <- flag_long %>%

    group_by(
        Tissue,
        Evidence
    ) %>%

    summarise(
        Candidate_N =
            n(),

        TRUE_N =
            sum(
                Parsed_flag %in% TRUE,
                na.rm = TRUE
            ),

        FALSE_N =
            sum(
                Parsed_flag %in% FALSE,
                na.rm = TRUE
            ),

        Unparsed_N =
            sum(
                is.na(
                    Parsed_flag
                )
            ),

        .groups = "drop"
    )


write_tsv(
    flag_summary,
    file.path(
        OUTDIR,
        "10F3B5A_binary_evidence_counts_by_tissue.tsv"
    )
)


# ============================================================
# 10. Exact categorical values
# ============================================================

categorical_columns <- c(
    "RNA_support_role",
    "FIMO_best_tier",
    "Interpretation",
    "Evidence_sources",
    "Max_expression_tissue"
)


category_audit <- bind_rows(
    lapply(
        categorical_columns,
        function(col) {

            master %>%

                transmute(
                    Field = col,

                    Value =
                        as.character(
                            .data[[col]]
                        )
                ) %>%

                count(
                    Field,
                    Value,
                    name = "N"
                )
        }
    )
) %>%

    arrange(
        Field,
        desc(N),
        Value
    )


write_tsv(
    category_audit,
    file.path(
        OUTDIR,
        "10F3B5A_exact_category_value_audit.tsv"
    )
)


# ============================================================
# 11. RNA-support-role × tissue
# ============================================================

rna_role_by_tissue <- master %>%

    count(
        Tissue,
        RNA_support_role,
        name = "N"
    ) %>%

    arrange(
        match(
            Tissue,
            tissue_order
        ),
        RNA_support_role
    )


write_tsv(
    rna_role_by_tissue,
    file.path(
        OUTDIR,
        "10F3B5A_RNA_support_role_by_tissue.tsv"
    )
)


# ============================================================
# 12. FIMO tier × tissue
# ============================================================

fimo_by_tissue <- master %>%

    count(
        Tissue,
        FIMO_best_tier,
        name = "N"
    ) %>%

    arrange(
        match(
            Tissue,
            tissue_order
        ),
        FIMO_best_tier
    )


write_tsv(
    fimo_by_tissue,
    file.path(
        OUTDIR,
        "10F3B5A_FIMO_best_tier_by_tissue.tsv"
    )
)


# ============================================================
# 13. Representative frozen TF audit
# Not a new selection; only verify whether known manuscript
# candidates exist in the authoritative master.
# ============================================================

representative_symbols <- c(
    "KLF9",
    "KLF12",
    "ELF4",
    "ETV2",
    "SP5"
)


representative_audit <- master %>%

    filter(
        TF_symbol %in%
            representative_symbols
    ) %>%

    select(
        Tissue,
        TF_symbol,
        gene_id,
        Target_tissue_mean_TPM,
        Tau,
        Target_tissue_SPM,
        Target_tissue_expression_rank,
        Max_expression_tissue,
        Target_is_max_expression_tissue,
        High_RNA_expression_support,
        RNA_support_role,
        FIMO_best_tier,
        Interpretation
    ) %>%

    arrange(
        match(
            Tissue,
            tissue_order
        ),
        TF_symbol
    )


write_tsv(
    representative_audit,
    file.path(
        OUTDIR,
        "10F3B5A_representative_TF_audit.tsv"
    )
)


# ============================================================
# 14. Source manifest
# ============================================================

manifest <- tibble(
    Role = c(
        "Primary_S5_source",
        "Priority_TF_crosscheck",
        "Eight_tissue_TPM_crosscheck",
        "High_RNA_network_crosscheck",
        "Publication_module_crosscheck"
    ),

    Frozen_source = c(
        F_MASTER,
        F_PRIORITY,
        F_TPM,
        F_NETWORK,
        F_MODULE
    ),

    Rows = c(
        nrow(master),
        nrow(priority),
        nrow(tpm),
        nrow(network),
        nrow(module)
    ),

    Used_for_reselection = "NO"
)


write_tsv(
    manifest,
    file.path(
        OUTDIR,
        "10F3B5A_source_manifest.tsv"
    )
)


# ============================================================
# 15. Human-readable preview
# ============================================================

preview_file <- file.path(
    OUTDIR,
    "10F3B5A_source_category_preview.txt"
)


sink(
    preview_file
)

cat("============================================================\n")
cat("STEP10F3B5A — S5 SOURCE/CATEGORY PREVIEW\n")
cat("============================================================\n\n")


cat("---- FROZEN SOURCE ROW QC ----\n")
print(
    row_qc,
    n = Inf,
    width = Inf
)


cat("\n---- CANDIDATE TF COUNTS BY TISSUE ----\n")
print(
    master_tissue_counts,
    n = Inf,
    width = Inf
)


cat("\n---- NUMERIC COMPLETENESS ----\n")
print(
    numeric_qc,
    n = Inf,
    width = Inf
)


cat("\n---- EXACT CATEGORY VALUES ----\n")
print(
    category_audit,
    n = Inf,
    width = Inf
)


cat("\n---- BOOLEAN FIELD VALUES ----\n")
print(
    flag_value_audit,
    n = Inf,
    width = Inf
)


cat("\n---- RNA SUPPORT ROLE BY TISSUE ----\n")
print(
    rna_role_by_tissue,
    n = Inf,
    width = Inf
)


cat("\n---- FIMO BEST TIER BY TISSUE ----\n")
print(
    fimo_by_tissue,
    n = Inf,
    width = Inf
)


cat("\n---- REPRESENTATIVE TF AUDIT ----\n")
print(
    representative_audit,
    n = Inf,
    width = Inf
)


cat("\n---- DUPLICATE TISSUE+TF KEYS ----\n")

if (nrow(duplicate_key) == 0) {

    cat("NONE\n")

} else {

    print(
        duplicate_key,
        n = Inf,
        width = Inf
    )
}


cat("\n============================================================\n")
cat("No candidate TF was reselected.\n")
cat("No RNA expression statistic was recalculated.\n")
cat("No FIMO analysis was rerun.\n")
cat("No motif analysis was rerun.\n")
cat("Biological results modified: 0.\n")
cat("============================================================\n")


sink()


# ============================================================
# 16. Console summary
# ============================================================

cat("Frozen source row QC:\n")
print(
    row_qc,
    width = Inf
)


cat("\nCandidate TF counts by tissue:\n")
print(
    master_tissue_counts,
    width = Inf
)


cat("\nExact RNA_support_role values:\n")

print(
    category_audit %>%
        filter(
            Field ==
                "RNA_support_role"
        ),
    n = Inf,
    width = Inf
)


cat("\nExact FIMO_best_tier values:\n")

print(
    category_audit %>%
        filter(
            Field ==
                "FIMO_best_tier"
        ),
    n = Inf,
    width = Inf
)


cat("\nBoolean-field parsing audit:\n")
print(
    flag_value_audit,
    n = Inf,
    width = Inf
)


cat("\nRepresentative TF audit:\n")
print(
    representative_audit,
    n = Inf,
    width = Inf
)


cat("\nOutputs:\n")
cat(
    preview_file,
    "\n"
)


cat("\nIMPORTANT:\n")
cat("- 08D2 is the authoritative source for S5.\n")
cat("- Other tables are crosschecks only.\n")
cat("- Exact category values were audited before plotting.\n")
cat("- No TF candidate was added or removed.\n")
cat("- No RNA statistic was recomputed.\n")
cat("- No FIMO/motif analysis was rerun.\n")
cat("- Biological results modified: 0.\n")


if (
    all(row_qc$Status == "PASS") &&
    all(tissue_qc$Status == "PASS")
) {

    cat("\nSTEP10F3B5A STATUS: PASS\n")
    cat(
        "S5 frozen sources and categorical values are ready for figure construction.\n"
    )

} else {

    cat("\nSTEP10F3B5A STATUS: CHECK\n")
}

cat("============================================================\n")
