#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

map_file <- paste0(
    "08E3B_GO_KEGG_v2/",
    "08E3B2C_master_ENSEMBL_SYMBOL_ENTREZ_mapping.tsv"
)

outdir <- "08E3B_GO_KEGG_v2"

x <- read_tsv(
    map_file,
    show_col_types = FALSE
)

cat("============================================================\n")
cat("STEP 08E3B2C ID MAPPING AMBIGUITY QC\n")
cat("============================================================\n")

cat("Mapping rows             :", nrow(x), "\n")
cat("Unique ENSEMBL gene IDs  :", n_distinct(x$ENSEMBL_gene_id), "\n")
cat(
    "Mapped ENSEMBL gene IDs  :",
    n_distinct(
        x$ENSEMBL_gene_id[
            !is.na(x$ENTREZID)
        ]
    ),
    "\n"
)

# ENSEMBL -> multiple ENTREZ
ens_multi <- x %>%
    filter(
        !is.na(ENTREZID)
    ) %>%
    distinct(
        ENSEMBL_gene_id,
        SYMBOL,
        ENTREZID
    ) %>%
    count(
        ENSEMBL_gene_id,
        name = "ENTREZ_N"
    ) %>%
    filter(
        ENTREZ_N > 1
    ) %>%
    arrange(
        desc(ENTREZ_N)
    )

# SYMBOL -> multiple ENTREZ
symbol_multi <- x %>%
    filter(
        !is.na(SYMBOL),
        !is.na(ENTREZID)
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
    ) %>%
    arrange(
        desc(ENTREZ_N)
    )

# multiple ENSEMBL -> same ENTREZ
entrez_multi <- x %>%
    filter(
        !is.na(ENTREZID)
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
    ) %>%
    arrange(
        desc(ENSEMBL_N)
    )

cat(
    "ENSEMBL IDs -> >1 ENTREZ :",
    nrow(ens_multi),
    "\n"
)

cat(
    "SYMBOLs -> >1 ENTREZ     :",
    nrow(symbol_multi),
    "\n"
)

cat(
    "ENTREZ -> >1 ENSEMBL     :",
    nrow(entrez_multi),
    "\n"
)

write_tsv(
    ens_multi,
    file.path(
        outdir,
        "08E3B2C_ENSEMBL_multiple_ENTREZ.tsv"
    )
)

write_tsv(
    symbol_multi,
    file.path(
        outdir,
        "08E3B2C_SYMBOL_multiple_ENTREZ.tsv"
    )
)

write_tsv(
    entrez_multi,
    file.path(
        outdir,
        "08E3B2C_ENTREZ_multiple_ENSEMBL.tsv"
    )
)

# ------------------------------------------------------------
# Check actual 8 target/universe mapping files
# ------------------------------------------------------------

dirs <- list.dirs(
    outdir,
    recursive = FALSE,
    full.names = TRUE
)

qc <- list()
ii <- 1

for (d in dirs) {

    tf <- file.path(
        d,
        "target_ENSEMBL_SYMBOL_ENTREZ.tsv"
    )

    uf <- file.path(
        d,
        "universe_ENSEMBL_SYMBOL_ENTREZ.tsv"
    )

    if (
        !file.exists(tf) ||
        !file.exists(uf)
    ) {
        next
    }

    target <- read_tsv(
        tf,
        show_col_types = FALSE
    )

    universe <- read_tsv(
        uf,
        show_col_types = FALSE
    )

    target_ambig <- target %>%
        filter(
            !is.na(ENTREZID)
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
        )

    universe_ambig <- universe %>%
        filter(
            !is.na(ENTREZID)
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
        )

    qc[[ii]] <- tibble(
        Set_ID = basename(d),

        Target_input_ENSSSCG_N =
            n_distinct(
                target$ENSEMBL_gene_id
            ),

        Target_ambiguous_ENSSSCG_N =
            nrow(target_ambig),

        Target_ambiguous_fraction =
            nrow(target_ambig) /
            n_distinct(
                target$ENSEMBL_gene_id
            ),

        Universe_input_ENSSSCG_N =
            n_distinct(
                universe$ENSEMBL_gene_id
            ),

        Universe_ambiguous_ENSSSCG_N =
            nrow(universe_ambig),

        Universe_ambiguous_fraction =
            nrow(universe_ambig) /
            n_distinct(
                universe$ENSEMBL_gene_id
            )
    )

    ii <- ii + 1
}

qc <- bind_rows(qc)

write_tsv(
    qc,
    file.path(
        outdir,
        "08E3B2C_selected_sets_mapping_ambiguity_QC.tsv"
    )
)

cat("\nSelected gene-set QC:\n")
print(
    qc,
    n = Inf,
    width = Inf
)

cat("\nInterpretation guideline:\n")
cat("<1% ambiguous ENSSSCG: negligible\n")
cat("1-5%: minor, document\n")
cat(">5%: inspect before final enrichment interpretation\n")

cat("============================================================\n")
cat("STEP 08E3B2C COMPLETED\n")
cat("============================================================\n")
