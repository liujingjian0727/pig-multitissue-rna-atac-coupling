#!/usr/bin/env Rscript

# ============================================================
# STEP 08E3B2-v2
#
# GO / KEGG enrichment for motif-associated gene sets
# using SAME-TISSUE STRONG-GENE UNIVERSES.
#
# ID strategy:
#
# ENSSSCG Ensembl gene ID
#        ↓
# Sus_longest.gtf gene_name
#        ↓
# org.Ss.eg.db SYMBOL
#        ↓
# ENTREZID
#
# GO:
#   ENTREZID + org.Ss.eg.db
#
# KEGG:
#   ENTREZID + organism = ssc
#
# ============================================================
#
# SOFTWARE INSTALLATION
#
# Missing R/Bioconductor packages are installed automatically.
#
# Alternative conda installation:
#
# conda install -y \
#   -c conda-forge \
#   -c bioconda \
#   r-base \
#   r-biocmanager \
#   bioconductor-clusterprofiler \
#   bioconductor-org.ss.eg.db \
#   bioconductor-annotationdbi \
#   bioconductor-enrichplot \
#   bioconductor-keggrest \
#   r-dplyr \
#   r-readr \
#   r-ggplot2
#
# ============================================================


# ============================================================
# 0. Local library and automatic installation
# ============================================================

local_lib <- file.path(
    getwd(),
    "R_lib_08E3B"
)

dir.create(
    local_lib,
    recursive = TRUE,
    showWarnings = FALSE
)

.libPaths(
    c(
        local_lib,
        .libPaths()
    )
)


cran_packages <- c(
    "dplyr",
    "readr",
    "ggplot2"
)


for (pkg in cran_packages) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

        message(
            "Installing CRAN package: ",
            pkg
        )

        install.packages(
            pkg,
            repos = "https://cloud.r-project.org",
            lib = local_lib
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
        lib = local_lib
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
            lib = local_lib
        )
    }
}


# ============================================================
# 1. Load packages
# ============================================================

suppressPackageStartupMessages({

    library(clusterProfiler)
    library(org.Ss.eg.db)
    library(AnnotationDbi)

    library(dplyr)
    library(readr)
    library(ggplot2)
})


# ============================================================
# 2. Configuration
# ============================================================

MANIFEST <- paste0(
    "08E3B_gene_sets/",
    "08E3B1_gene_set_manifest.tsv"
)

GTF_FILE <- "Sus_longest.gtf"

OUTDIR <- "08E3B_GO_KEGG_v2"

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


GO_ONTOLOGIES <- c(
    "BP",
    "MF",
    "CC"
)

KEGG_ORGANISM <- "ssc"

PADJ_CUTOFF <- 0.05

MIN_GS_SIZE <- 5

MAX_GS_SIZE <- 500


# ============================================================
# 3. Helper functions
# ============================================================

safe_read_gene_list <- function(path) {

    x <- readLines(
        path,
        warn = FALSE
    )

    x <- trimws(x)

    x <- x[
        x != ""
    ]

    unique(x)
}


safe_filename <- function(x) {

    gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        x
    )
}


write_status_file <- function(
    path,
    status
) {

    write_tsv(
        tibble(
            Status = status
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


# ============================================================
# 4. Package / OrgDb QC
# ============================================================

cat("============================================================\n")
cat("STEP 08E3B2-v2\n")
cat("============================================================\n")

cat(
    "R version:",
    R.version.string,
    "\n"
)

cat(
    "Bioconductor:",
    as.character(
        BiocManager::version()
    ),
    "\n"
)

cat(
    "org.Ss.eg.db:",
    as.character(
        packageVersion(
            "org.Ss.eg.db"
        )
    ),
    "\n"
)

cat(
    "clusterProfiler:",
    as.character(
        packageVersion(
            "clusterProfiler"
        )
    ),
    "\n"
)


db_keytypes <- keytypes(
    org.Ss.eg.db
)

cat("\norg.Ss.eg.db keytypes:\n")
cat(
    paste(
        db_keytypes,
        collapse = ", "
    ),
    "\n"
)


if (
    !"SYMBOL" %in% db_keytypes
) {

    stop(
        "ERROR: org.Ss.eg.db does not provide SYMBOL keytype."
    )
}


if (
    !"ENTREZID" %in% db_keytypes
) {

    stop(
        "ERROR: org.Ss.eg.db does not provide ENTREZID keytype."
    )
}


if (
    "ENSEMBL" %in% db_keytypes
) {

    cat(
        "NOTE: ENSEMBL is available, ",
        "but v2 intentionally uses GTF SYMBOL -> ENTREZID ",
        "for consistency.\n",
        sep = ""
    )

} else {

    cat(
        "PASS: ENSEMBL keytype is unavailable; ",
        "using GTF gene_name -> SYMBOL -> ENTREZID.\n",
        sep = ""
    )
}


# ============================================================
# 5. Read representative GTF:
#    Ensembl gene ID -> gene_name
# ============================================================

cat("\nReading GTF gene ID -> gene_name mapping...\n")


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
        FUN.VALUE =
            character(1)
    )
}


gtf_fields <- strsplit(
    gtf_lines,
    "\t",
    fixed = TRUE
)


gtf_attr <- vapply(
    gtf_fields,
    function(x) {

        if (length(x) >= 9) {
            x[9]
        } else {
            NA_character_
        }

    },
    FUN.VALUE =
        character(1)
)


gtf_gene_id <- extract_attr(
    gtf_attr,
    "gene_id"
)


gtf_gene_name <- extract_attr(
    gtf_attr,
    "gene_name"
)


gtf_map <- tibble(
    ENSEMBL_gene_id =
        gtf_gene_id,

    gene_name =
        gtf_gene_name
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
    "Unique GTF gene IDs:",
    nrow(gtf_map),
    "\n"
)


cat(
    "GTF genes with gene_name:",
    sum(
        !is.na(gtf_map$gene_name) &
        gtf_map$gene_name != ""
    ),
    "\n"
)


cat(
    "GTF gene_name coverage:",
    sprintf(
        "%.2f%%",
        100 *
        mean(
            !is.na(gtf_map$gene_name) &
            gtf_map$gene_name != ""
        )
    ),
    "\n"
)


write_tsv(
    gtf_map,
    file.path(
        OUTDIR,
        "08E3B2A_GTF_ENSEMBL_to_gene_name.tsv"
    )
)


# ============================================================
# 6. SYMBOL -> ENTREZID mapping from org.Ss.eg.db
# ============================================================

cat(
    "\nBuilding SYMBOL -> ENTREZID mapping from org.Ss.eg.db...\n"
)


all_symbols <- keys(
    org.Ss.eg.db,
    keytype = "SYMBOL"
)


symbol_entrez <- suppressMessages(

    AnnotationDbi::select(

        org.Ss.eg.db,

        keys =
            all_symbols,

        keytype =
            "SYMBOL",

        columns =
            c(
                "SYMBOL",
                "ENTREZID"
            )
    )
)


symbol_entrez <- symbol_entrez %>%

    filter(
        !is.na(SYMBOL),
        !is.na(ENTREZID)
    ) %>%

    distinct(
        SYMBOL,
        ENTREZID
    )


cat(
    "OrgDb SYMBOL->ENTREZ mappings:",
    nrow(symbol_entrez),
    "\n"
)


write_tsv(
    symbol_entrez,
    file.path(
        OUTDIR,
        "08E3B2B_orgSs_SYMBOL_to_ENTREZID.tsv"
    )
)


# ============================================================
# 7. Build master ENSSSCG -> SYMBOL -> ENTREZ map
# ============================================================

master_map <- gtf_map %>%

    rename(
        SYMBOL =
            gene_name
    ) %>%

    left_join(
        symbol_entrez,
        by = "SYMBOL"
    )


cat(
    "\nMaster GTF genes:",
    nrow(master_map),
    "\n"
)


cat(
    "Master genes with ENTREZID:",
    sum(
        !is.na(
            master_map$ENTREZID
        )
    ),
    "\n"
)


cat(
    "Overall ENSSSCG -> ENTREZ coverage:",
    sprintf(
        "%.2f%%",
        100 *
        mean(
            !is.na(
                master_map$ENTREZID
            )
        )
    ),
    "\n"
)


write_tsv(
    master_map,
    file.path(
        OUTDIR,
        "08E3B2C_master_ENSEMBL_SYMBOL_ENTREZ_mapping.tsv"
    )
)


# ============================================================
# 8. Read manifest
# ============================================================

manifest <- read_tsv(
    MANIFEST,
    show_col_types = FALSE
)


cat(
    "\nGene sets:",
    nrow(manifest),
    "\n"
)


if (
    nrow(manifest) != 8
) {

    stop(
        "ERROR: expected 8 frozen gene sets."
    )
}


# ============================================================
# 9. Containers
# ============================================================

mapping_qc_list <- list()

summary_list <- list()

go_all_list <- list()

kegg_all_list <- list()

go_index <- 1
kegg_index <- 1


# ============================================================
# 10. Analyze each frozen gene set
# ============================================================

for (
    i in seq_len(
        nrow(manifest)
    )
) {

    meta <- manifest[i, ]


    set_id <- meta$Set_ID

    tissue <- meta$Tissue

    motif_id <- meta$Motif_ID


    cat("\n")
    cat("============================================================\n")

    cat(
        "Processing:",
        set_id,
        "\n"
    )

    cat(
        "Tissue:",
        tissue,
        "\n"
    )

    cat(
        "Motif:",
        motif_id,
        "\n"
    )


    set_dir <- file.path(
        OUTDIR,
        safe_filename(
            set_id
        )
    )


    dir.create(
        set_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )


    # ========================================================
    # Read target and universe ENSSSCG IDs
    # ========================================================

    target_ens <- safe_read_gene_list(
        meta$Gene_list_txt
    )


    universe_ens <- safe_read_gene_list(
        meta$Universe_txt
    )


    if (
        !all(
            target_ens %in%
                universe_ens
        )
    ) {

        stop(
            set_id,
            ": target is not a subset ",
            "of same-tissue universe."
        )
    }


    # ========================================================
    # Map target
    # ========================================================

    target_mapping <- tibble(
        ENSEMBL_gene_id =
            target_ens
    ) %>%

        left_join(
            master_map,
            by =
                "ENSEMBL_gene_id"
        )


    universe_mapping <- tibble(
        ENSEMBL_gene_id =
            universe_ens
    ) %>%

        left_join(
            master_map,
            by =
                "ENSEMBL_gene_id"
        )


    write_tsv(
        target_mapping,
        file.path(
            set_dir,
            "target_ENSEMBL_SYMBOL_ENTREZ.tsv"
        )
    )


    write_tsv(
        universe_mapping,
        file.path(
            set_dir,
            "universe_ENSEMBL_SYMBOL_ENTREZ.tsv"
        )
    )


    # ========================================================
    # Mapping counts
    # ========================================================

    target_symbol_n <- sum(
        !is.na(
            target_mapping$SYMBOL
        ) &
        target_mapping$SYMBOL != ""
    )


    universe_symbol_n <- sum(
        !is.na(
            universe_mapping$SYMBOL
        ) &
        universe_mapping$SYMBOL != ""
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


    # Important:
    # target Entrez universe must be exactly matched to the
    # successfully mapped same-tissue universe.
    target_entrez <- intersect(
        target_entrez,
        universe_entrez
    )


    cat(
        "Target ENSEMBL:",
        length(target_ens),
        "\n"
    )

    cat(
        "Target with SYMBOL:",
        target_symbol_n,
        "\n"
    )

    cat(
        "Target unique ENTREZ:",
        length(target_entrez),
        "\n"
    )

    cat(
        "Universe ENSEMBL:",
        length(universe_ens),
        "\n"
    )

    cat(
        "Universe unique ENTREZ:",
        length(universe_entrez),
        "\n"
    )


    mapping_qc_list[[i]] <- tibble(

        Set_ID =
            set_id,

        Analysis_tier =
            meta$Analysis_tier,

        Evidence_level =
            meta$Evidence_level,

        Tissue =
            tissue,

        Motif_ID =
            motif_id,

        Input_target_ENSEMBL_N =
            length(target_ens),

        Input_universe_ENSEMBL_N =
            length(universe_ens),

        Target_with_GTF_symbol_N =
            target_symbol_n,

        Universe_with_GTF_symbol_N =
            universe_symbol_n,

        Target_unique_ENTREZ_N =
            length(target_entrez),

        Universe_unique_ENTREZ_N =
            length(universe_entrez),

        Target_symbol_fraction =
            target_symbol_n /
            length(target_ens),

        Universe_symbol_fraction =
            universe_symbol_n /
            length(universe_ens),

        Target_ENTREZ_fraction =
            length(target_entrez) /
            length(target_ens),

        Universe_ENTREZ_fraction =
            length(universe_entrez) /
            length(universe_ens)
    )


    # ========================================================
    # GO enrichment
    #
    # GO now uses ENTREZID, matching KEGG.
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
            length(target_entrez) < 3 ||
            length(universe_entrez) < 10
        ) {

            write_status_file(
                file.path(
                    set_dir,
                    paste0(
                        "GO_",
                        ont,
                        "_all.tsv"
                    )
                ),
                "Insufficient_mapped_ENTREZ_genes"
            )

            write_status_file(
                file.path(
                    set_dir,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005.tsv"
                    )
                ),
                "Insufficient_mapped_ENTREZ_genes"
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
                    "GO ",
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
            nrow(go_df) == 0
        ) {

            write_status_file(
                file.path(
                    set_dir,
                    paste0(
                        "GO_",
                        ont,
                        "_all.tsv"
                    )
                ),
                "No_enrichment_results"
            )

            write_status_file(
                file.path(
                    set_dir,
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

                Set_ID =
                    set_id,

                Analysis_tier =
                    meta$Analysis_tier,

                Evidence_level =
                    meta$Evidence_level,

                Tissue =
                    tissue,

                Motif_ID =
                    motif_id,

                Motif_name =
                    meta$Motif_name,

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
                set_dir,
                paste0(
                    "GO_",
                    ont,
                    "_all.tsv"
                )
            )
        )


        if (
            nrow(go_sig) > 0
        ) {

            write_tsv(
                go_sig,
                file.path(
                    set_dir,
                    paste0(
                        "GO_",
                        ont,
                        "_FDR005.tsv"
                    )
                )
            )


            pdf(
                file.path(
                    set_dir,
                    paste0(
                        "GO_",
                        ont,
                        "_dotplot.pdf"
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
                            nrow(go_sig)
                        )
                ) +
                    ggtitle(
                        paste0(
                            set_id,
                            " | GO ",
                            ont
                        )
                    )
            )


            dev.off()


        } else {

            write_status_file(
                file.path(
                    set_dir,
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
    # KEGG enrichment
    # ========================================================

    cat(
        "Running KEGG...\n"
    )


    kegg_sig_n <- 0

    kegg_status <- "PASS"


    if (
        length(target_entrez) < 3 ||
        length(universe_entrez) < 10
    ) {

        kegg_status <-
            "Insufficient_mapped_ENTREZ_genes"

        write_status_file(
            file.path(
                set_dir,
                "KEGG_all.tsv"
            ),
            kegg_status
        )

        write_status_file(
            file.path(
                set_dir,
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
                    "KEGG ERROR for ",
                    set_id,
                    ": ",
                    conditionMessage(e)
                )

                NULL
            }
        )


        kegg_df <- safe_df(
            ekegg
        )


        if (
            nrow(kegg_df) == 0
        ) {

            kegg_status <-
                "No_results_or_KEGG_unavailable"

            write_status_file(
                file.path(
                    set_dir,
                    "KEGG_all.tsv"
                ),
                kegg_status
            )

            write_status_file(
                file.path(
                    set_dir,
                    "KEGG_FDR005.tsv"
                ),
                "No_FDR005_terms"
            )


        } else {

            kegg_df <- kegg_df %>%

                mutate(

                    Set_ID =
                        set_id,

                    Analysis_tier =
                        meta$Analysis_tier,

                    Evidence_level =
                        meta$Evidence_level,

                    Tissue =
                        tissue,

                    Motif_ID =
                        motif_id,

                    Motif_name =
                        meta$Motif_name
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
                    set_dir,
                    "KEGG_all.tsv"
                )
            )


            if (
                nrow(kegg_sig) > 0
            ) {

                write_tsv(
                    kegg_sig,
                    file.path(
                        set_dir,
                        "KEGG_FDR005.tsv"
                    )
                )


                pdf(
                    file.path(
                        set_dir,
                        "KEGG_dotplot.pdf"
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
                                nrow(kegg_sig)
                            )
                    ) +
                        ggtitle(
                            paste0(
                                set_id,
                                " | KEGG"
                            )
                        )
                )


                dev.off()


            } else {

                write_status_file(
                    file.path(
                        set_dir,
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
    # Per-set summary
    # ========================================================

    summary_list[[i]] <- tibble(

        Set_ID =
            set_id,

        Analysis_tier =
            meta$Analysis_tier,

        Evidence_level =
            meta$Evidence_level,

        Tissue =
            tissue,

        Motif_ID =
            motif_id,

        Motif_name =
            meta$Motif_name,

        Input_gene_N =
            length(target_ens),

        Universe_gene_N =
            length(universe_ens),

        Mapped_target_ENTREZ_N =
            length(target_entrez),

        Mapped_universe_ENTREZ_N =
            length(universe_entrez),

        GO_BP_FDR005_N =
            go_sig_counts["BP"],

        GO_MF_FDR005_N =
            go_sig_counts["MF"],

        GO_CC_FDR005_N =
            go_sig_counts["CC"],

        KEGG_FDR005_N =
            kegg_sig_n,

        KEGG_status =
            kegg_status
    )
}


# ============================================================
# 11. Mapping QC
# ============================================================

mapping_qc <- bind_rows(
    mapping_qc_list
)


write_tsv(
    mapping_qc,
    file.path(
        OUTDIR,
        "08E3B2_ID_mapping_QC.tsv"
    )
)


# ============================================================
# 12. Combined GO
# ============================================================

if (
    length(
        go_all_list
    ) > 0
) {

    go_combined <- bind_rows(
        go_all_list
    )


    # Additional cross-gene-set sensitivity statistic.
    # Primary FDR remains clusterProfiler p.adjust per set.
    go_combined <- go_combined %>%

        group_by(
            Ontology
        ) %>%

        mutate(
            Global_BH_across_gene_sets =
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
            "08E3B3_GO_all_gene_sets_combined.tsv"
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
            "08E3B4_GO_FDR005_combined.tsv"
        )
    )
}


# ============================================================
# 13. Combined KEGG
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
            Global_BH_across_gene_sets =
                p.adjust(
                    pvalue,
                    method = "BH"
                )
        )


    write_tsv(
        kegg_combined,
        file.path(
            OUTDIR,
            "08E3B5_KEGG_all_gene_sets_combined.tsv"
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
            "08E3B6_KEGG_FDR005_combined.tsv"
        )
    )
}


# ============================================================
# 14. Final enrichment summary
# ============================================================

summary_df <- bind_rows(
    summary_list
)


write_tsv(
    summary_df,
    file.path(
        OUTDIR,
        "08E3B7_enrichment_summary.tsv"
    )
)


# ============================================================
# 15. Mapping warning classification
# ============================================================

mapping_qc <- mapping_qc %>%

    mutate(

        Target_mapping_QC =
            case_when(

                Target_ENTREZ_fraction >= 0.80 ~
                    "PASS_high",

                Target_ENTREZ_fraction >= 0.60 ~
                    "PASS_moderate",

                TRUE ~
                    "LOW_mapping"
            ),

        Universe_mapping_QC =
            case_when(

                Universe_ENTREZ_fraction >= 0.80 ~
                    "PASS_high",

                Universe_ENTREZ_fraction >= 0.60 ~
                    "PASS_moderate",

                TRUE ~
                    "LOW_mapping"
            )
    )


write_tsv(
    mapping_qc,
    file.path(
        OUTDIR,
        "08E3B2_ID_mapping_QC.tsv"
    )
)


# ============================================================
# 16. Session info
# ============================================================

capture.output(
    sessionInfo(),
    file =
        file.path(
            OUTDIR,
            "08E3B8_sessionInfo.txt"
        )
)


# ============================================================
# 17. Console output
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 08E3B2-v2 COMPLETED\n")
cat("============================================================\n")


cat("\nID mapping QC:\n")

print(
    mapping_qc,
    n = Inf,
    width = Inf
)


cat("\nEnrichment summary:\n")

print(
    summary_df,
    n = Inf,
    width = Inf
)


cat(
    "\nIMPORTANT:\n"
)

cat(
    "Background = mapped unique genes from the SAME-TISSUE ",
    "strong peak-gene universe.\n",
    sep = ""
)

cat(
    "GO and KEGG both use ENTREZID after ",
    "GTF gene_name -> org.Ss.eg.db SYMBOL mapping.\n",
    sep = ""
)

cat(
    "Primary significance = within-set BH p.adjust < 0.05.\n"
)

cat(
    "Global_BH_across_gene_sets is sensitivity only.\n"
)


cat("\nMain outputs:\n")

cat(
    file.path(
        OUTDIR,
        "08E3B2_ID_mapping_QC.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "08E3B4_GO_FDR005_combined.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "08E3B6_KEGG_FDR005_combined.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUTDIR,
        "08E3B7_enrichment_summary.tsv"
    ),
    "\n"
)

cat("============================================================\n")

