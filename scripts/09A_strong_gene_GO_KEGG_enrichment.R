#!/usr/bin/env Rscript

# ============================================================
# STEP 09A
#
# Tissue-specific functional enrichment of strongly coupled
# RNA-ATAC peak-gene associations.
#
# Foreground:
#   Unique genes from strong peak-gene pairs in each tissue.
#
# Background:
#   Unique genes from ALL concordant-HC peak-gene pairs
#   in the SAME tissue.
#
# Strong definition:
#   rho_P348 >= 0.5
#   rho_P350 >= 0.5
#   rho_mean >= 0.7
#
# ID mapping:
#   ENSSSCG
#     -> Sus_longest.gtf gene_name
#     -> org.Ss.eg.db SYMBOL
#     -> ENTREZID
#
# GO:
#   ENTREZID + org.Ss.eg.db
#
# KEGG:
#   ENTREZID + Sus scrofa organism = "ssc"
#
# Primary enrichment significance:
#   within-tissue BH adjusted P < 0.05
#
# Additional sensitivity:
#   global BH across tissue gene sets
#
# ============================================================
#
# SOFTWARE INSTALLATION
#
# The script automatically installs missing packages into:
#
#   ./R_lib_08E3B
#
# This reuses the project-local Bioconductor library already
# established during Step08E3B.
#
# Alternative conda installation if needed:
#
# conda install -y \
#   -c conda-forge \
#   -c bioconda \
#   r-base \
#   r-biocmanager \
#   r-dplyr \
#   r-readr \
#   r-ggplot2 \
#   bioconductor-clusterprofiler \
#   bioconductor-org.ss.eg.db \
#   bioconductor-annotationdbi \
#   bioconductor-enrichplot \
#   bioconductor-keggrest
#
# ============================================================


# ============================================================
# 0. Configuration
# ============================================================

options(timeout = 600)

EPS <- 1e-12

COUPLING_FILE <- paste0(
    "05_RNA_ATAC_quantitative_coupling/",
    "05_concordant_HC_peak_gene_correlations.tsv"
)

GTF_FILE <- "Sus_longest.gtf"

OUTDIR <- "09A_strong_gene_functional_enrichment"

LOCAL_LIB <- file.path(
    getwd(),
    "R_lib_08E3B"
)

dir.create(
    LOCAL_LIB,
    recursive = TRUE,
    showWarnings = FALSE
)

.libPaths(
    c(
        LOCAL_LIB,
        .libPaths()
    )
)


# ============================================================
# 1. Automatic package installation
# ============================================================

cran_packages <- c(
    "dplyr",
    "readr",
    "ggplot2"
)

for (pkg in cran_packages) {

    if (
        !requireNamespace(
            pkg,
            quietly = TRUE
        )
    ) {

        message(
            "Installing CRAN package: ",
            pkg
        )

        install.packages(
            pkg,
            repos = "https://cloud.r-project.org",
            lib = LOCAL_LIB
        )
    }
}


if (
    !requireNamespace(
        "BiocManager",
        quietly = TRUE
    )
) {

    message(
        "Installing BiocManager..."
    )

    install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org",
        lib = LOCAL_LIB
    )
}


bioc_packages <- c(
    "clusterProfiler",
    "org.Ss.eg.db",
    "AnnotationDbi",
    "enrichplot",
    "KEGGREST"
)

for (pkg in bioc_packages) {

    if (
        !requireNamespace(
            pkg,
            quietly = TRUE
        )
    ) {

        message(
            "Installing Bioconductor package: ",
            pkg
        )

        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE,
            lib = LOCAL_LIB
        )
    }
}


# ============================================================
# 2. Load packages
# ============================================================

suppressPackageStartupMessages({

    library(dplyr)
    library(readr)
    library(ggplot2)

    library(clusterProfiler)
    library(org.Ss.eg.db)
    library(AnnotationDbi)
    library(enrichplot)
    library(KEGGREST)
})


# ============================================================
# 3. Output directories
# ============================================================

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

GENESET_DIR <- file.path(
    OUTDIR,
    "09A_gene_sets"
)

TARGET_DIR <- file.path(
    GENESET_DIR,
    "strong_targets"
)

UNIVERSE_DIR <- file.path(
    GENESET_DIR,
    "same_tissue_concordant_HC_universes"
)

MAPPING_DIR <- file.path(
    OUTDIR,
    "09A_ID_mapping"
)

TISSUE_DIR <- file.path(
    OUTDIR,
    "09A_per_tissue_enrichment"
)

dir.create(
    TARGET_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    UNIVERSE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    MAPPING_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    TISSUE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 4. Analysis settings
# ============================================================

GO_ONTOLOGIES <- c(
    "BP",
    "MF",
    "CC"
)

PADJ_CUTOFF <- 0.05

MIN_GS_SIZE <- 5

MAX_GS_SIZE <- 500

KEGG_ORGANISM <- "ssc"


# ============================================================
# 5. Expected counts from completed upstream analysis
#
# Used only as QC; no biological result is generated from
# these hard-coded values.
# ============================================================

expected_counts <- tibble::tribble(

    ~Tissue,        ~Concordant_pair_N, ~Concordant_gene_N, ~Strong_pair_N, ~Strong_gene_N,

    "Adipose",      4L,                  4L,                  0L,             0L,
    "Cerebellum",   197L,                124L,                47L,            41L,
    "Cortex",       66L,                 51L,                 15L,            14L,
    "Hypothalamus", 16L,                 15L,                 10L,            10L,
    "Liver",        481L,                299L,                54L,            47L,
    "Lung",         38L,                 36L,                 6L,             6L,
    "Muscle",       3065L,               539L,                344L,           202L,
    "Spleen",       807L,                349L,                226L,           130L
)


# ============================================================
# 6. Helper functions
# ============================================================

safe_filename <- function(x) {

    gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        x
    )
}


write_status <- function(
    path,
    status_text
) {

    write_tsv(
        tibble(
            Status = status_text
        ),
        path
    )
}


safe_df <- function(x) {

    if (is.null(x)) {
        return(data.frame())
    }

    as.data.frame(x)
}


extract_attr <- function(
    x,
    field
) {

    pattern <- paste0(
        field,
        ' "([^"]+)"'
    )

    m <- regexec(
        pattern,
        x
    )

    z <- regmatches(
        x,
        m
    )

    vapply(
        z,
        function(y) {

            if (length(y) >= 2) {
                y[2]
            } else {
                NA_character_
            }

        },
        FUN.VALUE = character(1)
    )
}


size_flag <- function(n) {

    case_when(

        n >= 30 ~
            "Adequate_ge30",

        n >= 10 ~
            "Limited_10_29",

        n > 0 ~
            "Very_small_lt10",

        TRUE ~
            "No_strong_genes"
    )
}


mapping_flag <- function(x) {

    case_when(

        is.na(x) ~
            "NA",

        x >= 0.80 ~
            "PASS_high",

        x >= 0.60 ~
            "PASS_moderate",

        TRUE ~
            "LOW_mapping"
    )
}


# ============================================================
# 7. Software / OrgDb QC
# ============================================================

cat("============================================================\n")
cat("STEP 09A\n")
cat("Strong RNA-ATAC gene functional enrichment\n")
cat("============================================================\n")

cat(
    "R version              : ",
    R.version.string,
    "\n",
    sep = ""
)

cat(
    "Bioconductor version   : ",
    as.character(
        BiocManager::version()
    ),
    "\n",
    sep = ""
)

cat(
    "clusterProfiler        : ",
    as.character(
        packageVersion(
            "clusterProfiler"
        )
    ),
    "\n",
    sep = ""
)

cat(
    "org.Ss.eg.db           : ",
    as.character(
        packageVersion(
            "org.Ss.eg.db"
        )
    ),
    "\n",
    sep = ""
)


kt <- keytypes(
    org.Ss.eg.db
)

if (
    !"SYMBOL" %in% kt ||
    !"ENTREZID" %in% kt
) {

    stop(
        "ERROR: org.Ss.eg.db lacks SYMBOL or ENTREZID keytype."
    )
}

cat(
    "OrgDb SYMBOL           : TRUE\n"
)

cat(
    "OrgDb ENTREZID         : TRUE\n"
)

cat(
    "OrgDb ENSEMBL          : ",
    "ENSEMBL" %in% kt,
    "\n",
    sep = ""
)


# ============================================================
# 8. Read all concordant-HC peak-gene pairs
# ============================================================

cat("\nReading Step05 concordant-HC coupling table...\n")

coupling <- read_tsv(
    COUPLING_FILE,
    show_col_types = FALSE
)


required_cols <- c(
    "SAF_peak_id",
    "gene_id",
    "ATAC_Max_tissue",
    "RNA_Max_tissue",
    "rho_P348",
    "rho_P350",
    "rho_mean"
)


missing_cols <- setdiff(
    required_cols,
    names(coupling)
)


if (
    length(missing_cols) > 0
) {

    stop(
        "Missing Step05 columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


if (
    nrow(coupling) != 4674
) {

    stop(
        "ERROR: expected 4674 concordant-HC pairs, observed ",
        nrow(coupling)
    )
}


coupling <- coupling %>%

    mutate(

        Tissue =
            ATAC_Max_tissue,

        rho_P348 =
            as.numeric(
                rho_P348
            ),

        rho_P350 =
            as.numeric(
                rho_P350
            ),

        rho_mean =
            as.numeric(
                rho_mean
            )
    )


if (
    any(
        coupling$ATAC_Max_tissue !=
        coupling$RNA_Max_tissue
    )
) {

    stop(
        "ERROR: some concordant-HC rows have different ",
        "ATAC and RNA max tissues."
    )
}


cat(
    "Concordant-HC peak-gene pairs: ",
    nrow(coupling),
    "\n",
    sep = ""
)


# ============================================================
# 9. Reconstruct official strong 702 pairs
# ============================================================

strong <- coupling %>%

    filter(

        rho_P348 >=
            0.5 - EPS,

        rho_P350 >=
            0.5 - EPS,

        rho_mean >=
            0.7 - EPS
    )


cat(
    "Strong peak-gene pairs: ",
    nrow(strong),
    "\n",
    sep = ""
)


if (
    nrow(strong) != 702
) {

    stop(
        "ERROR: expected 702 strong pairs, observed ",
        nrow(strong)
    )
}


# ============================================================
# 10. Tissue pair/gene count QC
# ============================================================

actual_counts <- lapply(

    expected_counts$Tissue,

    function(tissue) {

        bg <- coupling %>%
            filter(
                Tissue == tissue
            )

        fg <- strong %>%
            filter(
                Tissue == tissue
            )


        tibble(

            Tissue =
                tissue,

            Concordant_pair_N =
                nrow(bg),

            Concordant_gene_N =
                n_distinct(
                    bg$gene_id
                ),

            Strong_pair_N =
                nrow(fg),

            Strong_gene_N =
                n_distinct(
                    fg$gene_id
                )
        )
    }

) %>%
    bind_rows()


count_qc <- expected_counts %>%

    rename_with(
        ~paste0(
            .x,
            "_expected"
        ),
        -Tissue
    ) %>%

    left_join(
        actual_counts,
        by = "Tissue"
    ) %>%

    mutate(

        Pair_gene_count_QC =
            ifelse(

                Concordant_pair_N_expected ==
                    Concordant_pair_N &

                Concordant_gene_N_expected ==
                    Concordant_gene_N &

                Strong_pair_N_expected ==
                    Strong_pair_N &

                Strong_gene_N_expected ==
                    Strong_gene_N,

                "PASS",
                "FAIL"
            ),

        Strong_gene_size_flag =
            size_flag(
                Strong_gene_N
            )
    )


write_tsv(
    count_qc,
    file.path(
        OUTDIR,
        "09A1_pair_gene_count_QC.tsv"
    )
)


cat("\nTissue pair/gene counts:\n")

print(
    count_qc,
    n = Inf,
    width = Inf
)


if (
    any(
        count_qc$Pair_gene_count_QC !=
            "PASS"
    )
) {

    stop(
        "ERROR: tissue-level upstream counts do not match ",
        "the completed Step05 results."
    )
}


# ============================================================
# 11. Generate tissue-specific target and universe gene lists
# ============================================================

cat("\nGenerating tissue gene sets...\n")


gene_set_manifest <- list()


for (
    i in seq_len(
        nrow(actual_counts)
    )
) {

    tissue <- actual_counts$Tissue[i]


    universe_genes <- coupling %>%

        filter(
            Tissue == tissue
        ) %>%

        pull(
            gene_id
        ) %>%

        unique() %>%

        sort()


    target_genes <- strong %>%

        filter(
            Tissue == tissue
        ) %>%

        pull(
            gene_id
        ) %>%

        unique() %>%

        sort()


    if (
        length(target_genes) > 0 &&
        !all(
            target_genes %in%
                universe_genes
        )
    ) {

        stop(
            tissue,
            ": strong genes are not a subset of ",
            "same-tissue concordant-HC universe."
        )
    }


    target_file <- file.path(
        TARGET_DIR,
        paste0(
            tissue,
            ".strong_genes.txt"
        )
    )


    universe_file <- file.path(
        UNIVERSE_DIR,
        paste0(
            tissue,
            ".concordant_HC_gene_universe.txt"
        )
    )


    writeLines(
        target_genes,
        target_file
    )


    writeLines(
        universe_genes,
        universe_file
    )


    gene_set_manifest[[i]] <- tibble(

        Tissue =
            tissue,

        Strong_pair_N =
            actual_counts$Strong_pair_N[i],

        Strong_gene_N =
            length(target_genes),

        Concordant_HC_pair_N =
            actual_counts$Concordant_pair_N[i],

        Universe_gene_N =
            length(universe_genes),

        Strong_gene_fraction_of_universe =
            ifelse(
                length(universe_genes) > 0,
                length(target_genes) /
                    length(universe_genes),
                NA_real_
            ),

        Input_size_flag =
            size_flag(
                length(target_genes)
            ),

        Analysis_status =
            ifelse(
                length(target_genes) == 0,
                "SKIP_no_strong_genes",
                "ANALYZE"
            ),

        Target_gene_file =
            target_file,

        Universe_gene_file =
            universe_file
    )
}


gene_set_manifest <- bind_rows(
    gene_set_manifest
)


write_tsv(
    gene_set_manifest,
    file.path(
        OUTDIR,
        "09A2_tissue_gene_set_manifest.tsv"
    )
)


cat("\nGene-set manifest:\n")

print(
    gene_set_manifest,
    n = Inf,
    width = Inf
)


# ============================================================
# 12. Read representative GTF:
#     ENSSSCG -> gene_name
# ============================================================

cat(
    "\nReading representative GTF for ",
    "ENSSSCG -> gene_name mapping...\n",
    sep = ""
)


gtf_lines <- readLines(
    GTF_FILE,
    warn = FALSE
)


gtf_lines <- gtf_lines[
    !grepl(
        "^#",
        gtf_lines
    )
]


gtf_parts <- strsplit(
    gtf_lines,
    "\t",
    fixed = TRUE
)


is_gene <- vapply(

    gtf_parts,

    function(x) {

        length(x) >= 9 &&
        x[3] == "gene"
    },

    FUN.VALUE =
        logical(1)
)


gene_parts <- gtf_parts[
    is_gene
]


gene_attrs <- vapply(

    gene_parts,

    function(x) {
        x[9]
    },

    FUN.VALUE =
        character(1)
)


gtf_gene_ids <- extract_attr(
    gene_attrs,
    "gene_id"
)


gtf_gene_names <- extract_attr(
    gene_attrs,
    "gene_name"
)


gtf_map <- tibble(

    ENSEMBL_gene_id =
        gtf_gene_ids,

    SYMBOL =
        gtf_gene_names
) %>%

    filter(
        !is.na(
            ENSEMBL_gene_id
        )
    ) %>%

    distinct(
        ENSEMBL_gene_id,
        .keep_all = TRUE
    )


cat(
    "Unique GTF genes       : ",
    nrow(gtf_map),
    "\n",
    sep = ""
)


cat(
    "Genes with gene_name   : ",
    sum(
        !is.na(
            gtf_map$SYMBOL
        ) &
        gtf_map$SYMBOL != ""
    ),
    "\n",
    sep = ""
)


cat(
    "GTF symbol coverage    : ",
    sprintf(
        "%.2f%%",
        100 *
        mean(
            !is.na(
                gtf_map$SYMBOL
            ) &
            gtf_map$SYMBOL != ""
        )
    ),
    "\n",
    sep = ""
)


write_tsv(
    gtf_map,
    file.path(
        MAPPING_DIR,
        "09A3A_GTF_ENSSSCG_to_SYMBOL.tsv"
    )
)


# ============================================================
# 13. Build SYMBOL -> ENTREZ map
# ============================================================

cat(
    "\nBuilding org.Ss.eg.db SYMBOL -> ENTREZID map...\n"
)


org_symbols <- keys(
    org.Ss.eg.db,
    keytype = "SYMBOL"
)


symbol_entrez <- suppressMessages(

    AnnotationDbi::select(

        org.Ss.eg.db,

        keys =
            org_symbols,

        keytype =
            "SYMBOL",

        columns =
            c(
                "SYMBOL",
                "ENTREZID"
            )
    )
) %>%

    filter(
        !is.na(
            SYMBOL
        ),
        !is.na(
            ENTREZID
        )
    ) %>%

    distinct(
        SYMBOL,
        ENTREZID
    )


write_tsv(
    symbol_entrez,
    file.path(
        MAPPING_DIR,
        "09A3B_orgSs_SYMBOL_to_ENTREZID.tsv"
    )
)


# ============================================================
# 14. Master mapping
#
# Many-to-many SYMBOL mappings are retained intentionally
# and quantified below.
# ============================================================

master_map <- suppressWarnings(

    left_join(
        gtf_map,
        symbol_entrez,
        by = "SYMBOL"
    )
)


write_tsv(
    master_map,
    file.path(
        MAPPING_DIR,
        "09A3C_master_ENSSSCG_SYMBOL_ENTREZ.tsv"
    )
)


cat(
    "Master mapping rows    : ",
    nrow(master_map),
    "\n",
    sep = ""
)


cat(
    "Mapped ENSSSCG genes   : ",
    n_distinct(
        master_map$ENSEMBL_gene_id[
            !is.na(
                master_map$ENTREZID
            )
        ]
    ),
    "\n",
    sep = ""
)


# ============================================================
# 15. Global mapping ambiguity QC
# ============================================================

ens_multi <- master_map %>%

    filter(
        !is.na(
            ENTREZID
        )
    ) %>%

    distinct(
        ENSEMBL_gene_id,
        ENTREZID
    ) %>%

    count(
        ENSEMBL_gene_id,
        name = "ENTREZ_N"
    ) %>%

    filter(
        ENTREZ_N > 1
    )


symbol_multi <- master_map %>%

    filter(
        !is.na(
            SYMBOL
        ),
        !is.na(
            ENTREZID
        )
    ) %>%

    distinct(
        SYMBOL,
        ENTREZID
    ) %>%

    count(
        SYMBOL,
        name = "ENTREZ_N"
    ) %>%

    filter(
        ENTREZ_N > 1
    )


entrez_multi <- master_map %>%

    filter(
        !is.na(
            ENTREZID
        )
    ) %>%

    distinct(
        ENSEMBL_gene_id,
        ENTREZID
    ) %>%

    count(
        ENTREZID,
        name = "ENSEMBL_N"
    ) %>%

    filter(
        ENSEMBL_N > 1
    )


global_mapping_ambiguity <- tibble(

    Metric = c(
        "Unique_GTF_ENSSSCG",
        "Mapped_ENSSSCG",
        "ENSSSCG_with_multiple_ENTREZ",
        "SYMBOL_with_multiple_ENTREZ",
        "ENTREZ_with_multiple_ENSSSCG"
    ),

    Value = c(
        nrow(gtf_map),

        n_distinct(
            master_map$ENSEMBL_gene_id[
                !is.na(
                    master_map$ENTREZID
                )
            ]
        ),

        nrow(
            ens_multi
        ),

        nrow(
            symbol_multi
        ),

        nrow(
            entrez_multi
        )
    )
)


write_tsv(
    global_mapping_ambiguity,
    file.path(
        MAPPING_DIR,
        "09A3D_global_mapping_ambiguity_QC.tsv"
    )
)


write_tsv(
    ens_multi,
    file.path(
        MAPPING_DIR,
        "09A3E_ENSSSCG_multiple_ENTREZ.tsv"
    )
)


# ============================================================
# 16. Analyze each tissue
# ============================================================

mapping_qc_list <- list()

enrichment_summary_list <- list()

go_all_list <- list()

kegg_all_list <- list()

go_index <- 1
kegg_index <- 1
summary_index <- 1
mapping_index <- 1


for (
    tissue in gene_set_manifest$Tissue
) {

    meta <- gene_set_manifest %>%
        filter(
            Tissue == tissue
        )


    cat("\n")
    cat("============================================================\n")
    cat(
        "Processing tissue: ",
        tissue,
        "\n",
        sep = ""
    )
    cat("============================================================\n")


    tissue_out <- file.path(
        TISSUE_DIR,
        tissue
    )


    dir.create(
        tissue_out,
        recursive = TRUE,
        showWarnings = FALSE
    )


    target_ens <- readLines(
        meta$Target_gene_file,
        warn = FALSE
    )


    universe_ens <- readLines(
        meta$Universe_gene_file,
        warn = FALSE
    )


    target_ens <- unique(
        target_ens[
            target_ens != ""
        ]
    )


    universe_ens <- unique(
        universe_ens[
            universe_ens != ""
        ]
    )


    # --------------------------------------------------------
    # Adipose has no strong genes
    # --------------------------------------------------------

    if (
        length(target_ens) == 0
    ) {

        cat(
            "No strong genes. Enrichment skipped.\n"
        )


        mapping_qc_list[[mapping_index]] <- tibble(

            Tissue =
                tissue,

            Input_target_ENSSSCG_N =
                0,

            Input_universe_ENSSSCG_N =
                length(
                    universe_ens
                ),

            Target_with_SYMBOL_N =
                0,

            Universe_with_SYMBOL_N =
                NA_integer_,

            Target_unique_ENTREZ_N =
                0,

            Universe_unique_ENTREZ_N =
                NA_integer_,

            Target_ENTREZ_fraction =
                NA_real_,

            Universe_ENTREZ_fraction =
                NA_real_,

            Target_ambiguous_ENSSSCG_N =
                0,

            Universe_ambiguous_ENSSSCG_N =
                NA_integer_,

            Target_mapping_QC =
                "SKIP",

            Universe_mapping_QC =
                "SKIP"
        )


        mapping_index <-
            mapping_index + 1


        enrichment_summary_list[[summary_index]] <- tibble(

            Tissue =
                tissue,

            Strong_gene_N =
                0,

            Universe_gene_N =
                length(
                    universe_ens
                ),

            Mapped_target_ENTREZ_N =
                0,

            Mapped_universe_ENTREZ_N =
                NA_integer_,

            Input_size_flag =
                "No_strong_genes",

            GO_BP_FDR005_N =
                0,

            GO_MF_FDR005_N =
                0,

            GO_CC_FDR005_N =
                0,

            KEGG_FDR005_N =
                0,

            Analysis_status =
                "SKIPPED_no_strong_genes"
        )


        summary_index <-
            summary_index + 1


        next
    }


    # --------------------------------------------------------
    # ENSSSCG -> SYMBOL -> ENTREZID mapping
    # --------------------------------------------------------

    target_mapping <- tibble(
        ENSEMBL_gene_id =
            target_ens
    ) %>%

        left_join(
            master_map,
            by = "ENSEMBL_gene_id"
        )


    universe_mapping <- tibble(
        ENSEMBL_gene_id =
            universe_ens
    ) %>%

        left_join(
            master_map,
            by = "ENSEMBL_gene_id"
        )


    write_tsv(
        target_mapping,
        file.path(
            tissue_out,
            paste0(
                tissue,
                ".strong_gene_ID_mapping.tsv"
            )
        )
    )


    write_tsv(
        universe_mapping,
        file.path(
            tissue_out,
            paste0(
                tissue,
                ".universe_gene_ID_mapping.tsv"
            )
        )
    )


    target_symbol_n <- n_distinct(
        target_mapping$ENSEMBL_gene_id[
            !is.na(
                target_mapping$SYMBOL
            ) &
            target_mapping$SYMBOL != ""
        ]
    )


    universe_symbol_n <- n_distinct(
        universe_mapping$ENSEMBL_gene_id[
            !is.na(
                universe_mapping$SYMBOL
            ) &
            universe_mapping$SYMBOL != ""
        ]
    )


    target_entrez <- unique(
        target_mapping$ENTREZID[
            !is.na(
                target_mapping$ENTREZID
            )
        ]
    )


    universe_entrez <- unique(
        universe_mapping$ENTREZID[
            !is.na(
                universe_mapping$ENTREZID
            )
        ]
    )


    # Target is forced into the mapped same-tissue universe.
    target_entrez <- intersect(
        target_entrez,
        universe_entrez
    )


    target_ambig_n <- target_mapping %>%

        filter(
            !is.na(
                ENTREZID
            )
        ) %>%

        distinct(
            ENSEMBL_gene_id,
            ENTREZID
        ) %>%

        count(
            ENSEMBL_gene_id
        ) %>%

        filter(
            n > 1
        ) %>%

        nrow()


    universe_ambig_n <- universe_mapping %>%

        filter(
            !is.na(
                ENTREZID
            )
        ) %>%

        distinct(
            ENSEMBL_gene_id,
            ENTREZID
        ) %>%

        count(
            ENSEMBL_gene_id
        ) %>%

        filter(
            n > 1
        ) %>%

        nrow()


    target_map_fraction <-
        length(target_entrez) /
        length(target_ens)


    universe_map_fraction <-
        length(universe_entrez) /
        length(universe_ens)


    mapping_qc_list[[mapping_index]] <- tibble(

        Tissue =
            tissue,

        Input_target_ENSSSCG_N =
            length(
                target_ens
            ),

        Input_universe_ENSSSCG_N =
            length(
                universe_ens
            ),

        Target_with_SYMBOL_N =
            target_symbol_n,

        Universe_with_SYMBOL_N =
            universe_symbol_n,

        Target_unique_ENTREZ_N =
            length(
                target_entrez
            ),

        Universe_unique_ENTREZ_N =
            length(
                universe_entrez
            ),

        Target_SYMBOL_fraction =
            target_symbol_n /
            length(
                target_ens
            ),

        Universe_SYMBOL_fraction =
            universe_symbol_n /
            length(
                universe_ens
            ),

        Target_ENTREZ_fraction =
            target_map_fraction,

        Universe_ENTREZ_fraction =
            universe_map_fraction,

        Target_ambiguous_ENSSSCG_N =
            target_ambig_n,

        Universe_ambiguous_ENSSSCG_N =
            universe_ambig_n,

        Target_ambiguous_fraction =
            target_ambig_n /
            length(
                target_ens
            ),

        Universe_ambiguous_fraction =
            universe_ambig_n /
            length(
                universe_ens
            ),

        Target_mapping_QC =
            mapping_flag(
                target_map_fraction
            ),

        Universe_mapping_QC =
            mapping_flag(
                universe_map_fraction
            )
    )


    mapping_index <-
        mapping_index + 1


    cat(
        "Target ENSSSCG          : ",
        length(target_ens),
        "\n",
        sep = ""
    )

    cat(
        "Target unique ENTREZ    : ",
        length(target_entrez),
        "\n",
        sep = ""
    )

    cat(
        "Universe ENSSSCG        : ",
        length(universe_ens),
        "\n",
        sep = ""
    )

    cat(
        "Universe unique ENTREZ  : ",
        length(universe_entrez),
        "\n",
        sep = ""
    )

    cat(
        "Target mapping          : ",
        sprintf(
            "%.2f%%",
            100 *
            target_map_fraction
        ),
        "\n",
        sep = ""
    )

    cat(
        "Universe mapping        : ",
        sprintf(
            "%.2f%%",
            100 *
            universe_map_fraction
        ),
        "\n",
        sep = ""
    )


    # ========================================================
    # 17. GO enrichment
    # ========================================================

    go_sig_counts <- c(
        BP = 0,
        MF = 0,
        CC = 0
    )


    for (
        ont in GO_ONTOLOGIES
    ) {

        cat(
            "Running GO ",
            ont,
            "...\n",
            sep = ""
        )


        if (
            length(
                target_entrez
            ) < 3 ||
            length(
                universe_entrez
            ) < 10
        ) {

            write_status(
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_all.tsv"
                    )
                ),
                "Insufficient_mapped_genes"
            )

            write_status(
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005.tsv"
                    )
                ),
                "Insufficient_mapped_genes"
            )

            next
        }


        ego <- tryCatch(

            enrichGO(

                gene =
                    target_entrez,

                universe =
                    universe_entrez,

                OrgDb =
                    org.Ss.eg.db,

                keyType =
                    "ENTREZID",

                ont =
                    ont,

                pAdjustMethod =
                    "BH",

                pvalueCutoff =
                    1,

                qvalueCutoff =
                    1,

                minGSSize =
                    MIN_GS_SIZE,

                maxGSSize =
                    MAX_GS_SIZE,

                readable =
                    TRUE
            ),

            error = function(e) {

                message(
                    tissue,
                    " GO ",
                    ont,
                    " ERROR: ",
                    conditionMessage(e)
                )

                NULL
            }
        )


        go_df <- safe_df(
            ego
        )


        if (
            nrow(
                go_df
            ) == 0
        ) {

            write_status(
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_all.tsv"
                    )
                ),
                "No_enrichment_results"
            )

            write_status(
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005.tsv"
                    )
                ),
                "No_FDR005_terms"
            )

            next
        }


        go_df <- go_df %>%

            mutate(

                Tissue =
                    tissue,

                Strong_gene_input_N =
                    length(
                        target_ens
                    ),

                Universe_gene_input_N =
                    length(
                        universe_ens
                    ),

                Mapped_target_ENTREZ_N =
                    length(
                        target_entrez
                    ),

                Mapped_universe_ENTREZ_N =
                    length(
                        universe_entrez
                    ),

                Input_size_flag =
                    size_flag(
                        length(
                            target_ens
                        )
                    ),

                Ontology =
                    ont
            )


        go_sig <- go_df %>%

            filter(
                !is.na(
                    p.adjust
                ),
                p.adjust <
                    PADJ_CUTOFF
            )


        go_sig_counts[ont] <-
            nrow(
                go_sig
            )


        write_tsv(
            go_df,
            file.path(
                tissue_out,
                paste0(
                    "GO_",
                    ont,
                    "_all.tsv"
                )
            )
        )


        if (
            nrow(
                go_sig
            ) > 0
        ) {

            write_tsv(
                go_sig,
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005.tsv"
                    )
                )
            )


            pdf(
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005_dotplot.pdf"
                    )
                ),
                width = 9,
                height = 6
            )


            print(
                enrichplot::dotplot(
                    ego,
                    showCategory =
                        min(
                            20,
                            nrow(
                                go_sig
                            )
                        )
                ) +
                    ggtitle(
                        paste0(
                            tissue,
                            " strong genes | GO ",
                            ont
                        )
                    )
            )


            dev.off()


        } else {

            write_status(
                file.path(
                    tissue_out,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005.tsv"
                    )
                ),
                "No_FDR005_terms"
            )
        }


        go_all_list[[go_index]] <-
            go_df

        go_index <-
            go_index + 1
    }


    # ========================================================
    # 18. KEGG enrichment
    # ========================================================

    cat(
        "Running KEGG...\n"
    )


    kegg_sig_n <- 0

    kegg_status <- "PASS"


    if (
        length(
            target_entrez
        ) < 3 ||
        length(
            universe_entrez
        ) < 10
    ) {

        kegg_status <-
            "Insufficient_mapped_genes"

        write_status(
            file.path(
                tissue_out,
                "KEGG_all.tsv"
            ),
            kegg_status
        )

        write_status(
            file.path(
                tissue_out,
                "KEGG_FDR005.tsv"
            ),
            kegg_status
        )


    } else {

        ekegg <- tryCatch(

            enrichKEGG(

                gene =
                    target_entrez,

                universe =
                    universe_entrez,

                organism =
                    KEGG_ORGANISM,

                keyType =
                    "ncbi-geneid",

                pAdjustMethod =
                    "BH",

                pvalueCutoff =
                    1,

                qvalueCutoff =
                    1,

                minGSSize =
                    MIN_GS_SIZE,

                maxGSSize =
                    MAX_GS_SIZE
            ),

            error = function(e) {

                message(
                    tissue,
                    " KEGG ERROR: ",
                    conditionMessage(e)
                )

                NULL
            }
        )


        kegg_df <- safe_df(
            ekegg
        )


        if (
            nrow(
                kegg_df
            ) == 0
        ) {

            kegg_status <-
                "No_results_or_KEGG_unavailable"

            write_status(
                file.path(
                    tissue_out,
                    "KEGG_all.tsv"
                ),
                kegg_status
            )

            write_status(
                file.path(
                    tissue_out,
                    "KEGG_FDR005.tsv"
                ),
                "No_FDR005_terms"
            )


        } else {

            kegg_df <- kegg_df %>%

                mutate(

                    Tissue =
                        tissue,

                    Strong_gene_input_N =
                        length(
                            target_ens
                        ),

                    Universe_gene_input_N =
                        length(
                            universe_ens
                        ),

                    Mapped_target_ENTREZ_N =
                        length(
                            target_entrez
                        ),

                    Mapped_universe_ENTREZ_N =
                        length(
                            universe_entrez
                        ),

                    Input_size_flag =
                        size_flag(
                            length(
                                target_ens
                            )
                        )
                )


            kegg_sig <- kegg_df %>%

                filter(
                    !is.na(
                        p.adjust
                    ),
                    p.adjust <
                        PADJ_CUTOFF
            )


            kegg_sig_n <-
                nrow(
                    kegg_sig
                )


            write_tsv(
                kegg_df,
                file.path(
                    tissue_out,
                    "KEGG_all.tsv"
                )
            )


            if (
                nrow(
                    kegg_sig
                ) > 0
            ) {

                write_tsv(
                    kegg_sig,
                    file.path(
                        tissue_out,
                        "KEGG_FDR005.tsv"
                    )
                )


                pdf(
                    file.path(
                        tissue_out,
                        "KEGG_FDR005_dotplot.pdf"
                    ),
                    width = 9,
                    height = 6
                )


                print(
                    enrichplot::dotplot(
                        ekegg,
                        showCategory =
                            min(
                                20,
                                nrow(
                                    kegg_sig
                                )
                            )
                    ) +
                        ggtitle(
                            paste0(
                                tissue,
                                " strong genes | KEGG"
                            )
                        )
                )


                dev.off()


            } else {

                write_status(
                    file.path(
                        tissue_out,
                        "KEGG_FDR005.tsv"
                    ),
                    "No_FDR005_terms"
                )
            }


            kegg_all_list[[kegg_index]] <-
                kegg_df

            kegg_index <-
                kegg_index + 1
        }
    }


    # ========================================================
    # 19. Tissue enrichment summary
    # ========================================================

    enrichment_summary_list[[summary_index]] <- tibble(

        Tissue =
            tissue,

        Strong_pair_N =
            meta$Strong_pair_N,

        Strong_gene_N =
            length(
                target_ens
            ),

        Concordant_HC_pair_N =
            meta$Concordant_HC_pair_N,

        Universe_gene_N =
            length(
                universe_ens
            ),

        Mapped_target_ENTREZ_N =
            length(
                target_entrez
            ),

        Mapped_universe_ENTREZ_N =
            length(
                universe_entrez
            ),

        Input_size_flag =
            size_flag(
                length(
                    target_ens
                )
            ),

        GO_BP_FDR005_N =
            go_sig_counts["BP"],

        GO_MF_FDR005_N =
            go_sig_counts["MF"],

        GO_CC_FDR005_N =
            go_sig_counts["CC"],

        KEGG_FDR005_N =
            kegg_sig_n,

        KEGG_status =
            kegg_status,

        Analysis_status =
            "COMPLETED"
    )


    summary_index <-
        summary_index + 1
}


# ============================================================
# 20. Mapping QC table
# ============================================================

mapping_qc <- bind_rows(
    mapping_qc_list
)


write_tsv(
    mapping_qc,
    file.path(
        OUTDIR,
        "09A4_ID_mapping_QC.tsv"
    )
)


# ============================================================
# 21. Combined GO results
# ============================================================

if (
    length(
        go_all_list
    ) > 0
) {

    go_combined <- bind_rows(
        go_all_list
    )


    # Additional sensitivity:
    # BH across all tissue gene sets within each ontology.
    #
    # Primary significance remains the original
    # clusterProfiler within-tissue p.adjust.
    go_combined <- go_combined %>%

        group_by(
            Ontology
        ) %>%

        mutate(

            Global_BH_across_tissues =
                p.adjust(
                    pvalue,
                    method = "BH"
                )
        ) %>%

        ungroup()


    write_tsv(
        go_combined,
        file.path(
            OUTDIR,
            "09A5_GO_all_tissues_combined.tsv"
        )
    )


    go_sig_combined <- go_combined %>%

        filter(
            !is.na(
                p.adjust
            ),
            p.adjust <
                PADJ_CUTOFF
        )


    write_tsv(
        go_sig_combined,
        file.path(
            OUTDIR,
            "09A6_GO_FDR005_combined.tsv"
        )
    )


    # Top nominal terms for review only.
    go_top <- go_combined %>%

        group_by(
            Tissue,
            Ontology
        ) %>%

        arrange(
            pvalue,
            .by_group = TRUE
        ) %>%

        slice_head(
            n = 10
        ) %>%

        ungroup()


    write_tsv(
        go_top,
        file.path(
            OUTDIR,
            "09A7_GO_top10_nominal_for_review.tsv"
        )
    )
}


# ============================================================
# 22. Combined KEGG results
# ============================================================

if (
    length(
        kegg_all_list
    ) > 0
) {

    kegg_combined <- bind_rows(
        kegg_all_list
    )


    kegg_combined <- kegg_combined %>%

        mutate(

            Global_BH_across_tissues =
                p.adjust(
                    pvalue,
                    method = "BH"
                )
        )


    write_tsv(
        kegg_combined,
        file.path(
            OUTDIR,
            "09A8_KEGG_all_tissues_combined.tsv"
        )
    )


    kegg_sig_combined <- kegg_combined %>%

        filter(
            !is.na(
                p.adjust
            ),
            p.adjust <
                PADJ_CUTOFF
        )


    write_tsv(
        kegg_sig_combined,
        file.path(
            OUTDIR,
            "09A9_KEGG_FDR005_combined.tsv"
        )
    )


    kegg_top <- kegg_combined %>%

        group_by(
            Tissue
        ) %>%

        arrange(
            pvalue,
            .by_group = TRUE
        ) %>%

        slice_head(
            n = 10
        ) %>%

        ungroup()


    write_tsv(
        kegg_top,
        file.path(
            OUTDIR,
            "09A10_KEGG_top10_nominal_for_review.tsv"
        )
    )
}


# ============================================================
# 23. Final enrichment summary
# ============================================================

enrichment_summary <- bind_rows(
    enrichment_summary_list
)


write_tsv(
    enrichment_summary,
    file.path(
        OUTDIR,
        "09A11_enrichment_summary.tsv"
    )
)


# ============================================================
# 24. Final integrated QC
# ============================================================

overall_qc <- tibble(

    Metric = c(

        "Concordant_HC_peak_gene_pairs",

        "Official_strong_peak_gene_pairs",

        "Tissues_total",

        "Tissues_with_strong_genes",

        "Adipose_strong_genes",

        "Pair_gene_count_reconstruction",

        "OrgDb_SYMBOL_available",

        "OrgDb_ENTREZID_available"
    ),

    Value = c(

        as.character(
            nrow(
                coupling
            )
        ),

        as.character(
            nrow(
                strong
            )
        ),

        as.character(
            nrow(
                expected_counts
            )
        ),

        as.character(
            sum(
                actual_counts$Strong_gene_N >
                    0
            )
        ),

        as.character(
            actual_counts$Strong_gene_N[
                actual_counts$Tissue ==
                    "Adipose"
            ]
        ),

        ifelse(
            all(
                count_qc$Pair_gene_count_QC ==
                    "PASS"
            ),
            "PASS",
            "FAIL"
        ),

        as.character(
            "SYMBOL" %in%
                kt
        ),

        as.character(
            "ENTREZID" %in%
                kt
        )
    ),

    Expected = c(

        "4674",

        "702",

        "8",

        "7",

        "0",

        "PASS",

        "TRUE",

        "TRUE"
    )
)


overall_qc <- overall_qc %>%

    mutate(

        Status =
            ifelse(
                Value ==
                    Expected,
                "PASS",
                "CHECK"
            )
    )


write_tsv(
    overall_qc,
    file.path(
        OUTDIR,
        "09A12_overall_QC_summary.tsv"
    )
)


# ============================================================
# 25. Session information
# ============================================================

capture.output(
    sessionInfo(),
    file =
        file.path(
            OUTDIR,
            "09A13_sessionInfo.txt"
        )
)


# ============================================================
# 26. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 09A COMPLETED\n")
cat("============================================================\n")


cat("\nOverall QC:\n")

print(
    overall_qc,
    n = Inf,
    width = Inf
)


cat("\nID mapping QC:\n")

print(
    mapping_qc,
    n = Inf,
    width = Inf
)


cat("\nFunctional enrichment summary:\n")

print(
    enrichment_summary,
    n = Inf,
    width = Inf
)


cat("\nInterpretation rules:\n")

cat(
    "1. Foreground = unique strong-coupled genes within tissue.\n"
)

cat(
    "2. Background = unique concordant-HC genes from SAME tissue.\n"
)

cat(
    "3. Primary significance = within-tissue BH p.adjust < 0.05.\n"
)

cat(
    "4. Global_BH_across_tissues is an additional sensitivity statistic.\n"
)

cat(
    "5. Very_small_lt10 tissues are exploratory because of low power.\n"
)

cat(
    "6. Strong peak-gene associations remain candidate genomic associations, not validated regulation.\n"
)


cat("\nMain outputs:\n")

cat(
    file.path(
        OUTDIR,
        "09A1_pair_gene_count_QC.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "09A2_tissue_gene_set_manifest.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "09A4_ID_mapping_QC.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "09A6_GO_FDR005_combined.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "09A9_KEGG_FDR005_combined.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "09A11_enrichment_summary.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "09A12_overall_QC_summary.tsv"
    ),
    "\n"
)

cat("============================================================\n")

