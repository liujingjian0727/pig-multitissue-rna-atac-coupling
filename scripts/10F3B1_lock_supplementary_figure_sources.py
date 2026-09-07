#!/usr/bin/env python3

from pathlib import Path
import csv
import os
import re

OUTDIR = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F3B_supplementary_inventory/"
    "10F3B1_source_lock"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Supplementary figure source groups
# ============================================================

GROUPS = {

    "FigureS2_permutation": {
        "roots": [
            "05B_same_tissue_permutation_v2",
        ],
        "keywords": [
            "global",
            "tissue",
            "summary",
            "permutation",
            "threshold",
            "null",
            "strong",
            "observed",
        ],
    },

    "FigureS3_genomic_architecture": {
        "roots": [
            "06_strong_peak_gene_structure",
            "06B_structural_enrichment",
            "07_promoter_first_enrichment",
        ],
        "keywords": [
            "summary",
            "structure",
            "class",
            "enrichment",
            "relative",
            "multiplicity",
            "distance",
            "tss",
            "fisher",
        ],
    },

    "FigureS4_motif_QC": {
        "roots": [
            "08A_matched_motif_background_v2",
            "08B1_AME",
            "08B2_STREME",
        ],
        "keywords": [
            "summary",
            "qc",
            "match",
            "background",
            "ame",
            "streme",
            "tomtom",
            "motif",
            "significant",
            "result",
        ],
    },

    "FigureS5_motif_TF_support": {
        "roots": [
            "08C_TF_RNA_integration",
            "08D_motif_family_TF_master",
            "08E1c_FIMO",
            "08E2_motif_peak_gene_TF_network",
            "08E3A_motif_regulatory_architecture",
            "09C_publication_candidate_modules",
        ],
        "keywords": [
            "summary",
            "master",
            "candidate",
            "tf",
            "rna",
            "fimo",
            "motif",
            "module",
            "architecture",
            "network",
        ],
    },
}


# ============================================================
# Helpers
# ============================================================

allowed_suffixes = {
    ".tsv",
    ".csv",
    ".txt",
}

def row_col_count(path):

    suffix = path.suffix.lower()

    try:

        with open(
            path,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as f:

            first = f.readline().rstrip("\n\r")

        if not first:
            return 0, 0, ""

        if suffix == ".csv":
            delim = ","
        else:
            delim = "\t"

        header = first.split(delim)
        col_n = len(header)

        with open(
            path,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as f:

            n = sum(1 for _ in f)

        row_n = max(0, n - 1)

        return row_n, col_n, first

    except Exception:

        return -1, -1, ""


def keyword_score(filename, keywords):

    name = filename.lower()

    score = 0
    hits = []

    for k in keywords:

        if k.lower() in name:
            score += 1
            hits.append(k)

    return score, ",".join(hits)


# ============================================================
# Inventory
# ============================================================

records = []

for group, cfg in GROUPS.items():

    for root_string in cfg["roots"]:

        root = Path(root_string)

        if not root.exists():
            continue

        for path in root.rglob("*"):

            if not path.is_file():
                continue

            if path.suffix.lower() not in allowed_suffixes:
                continue

            score, hits = keyword_score(
                path.name,
                cfg["keywords"]
            )

            # Keep all tables, but score the likely source tables.
            rows, cols, header = row_col_count(path)

            records.append({
                "Supplementary_figure": group,
                "Root": str(root),
                "File": str(path),
                "Filename": path.name,
                "Keyword_score": score,
                "Keyword_hits": hits,
                "Rows_excluding_header": rows,
                "Columns": cols,
                "Size_bytes": path.stat().st_size,
                "Header": header,
            })


# ============================================================
# Sort
# ============================================================

records.sort(
    key=lambda x: (
        x["Supplementary_figure"],
        -x["Keyword_score"],
        x["Filename"],
    )
)


# ============================================================
# Full source-lock inventory
# ============================================================

full_file = OUTDIR / "10F3B1_all_candidate_source_tables.tsv"

fields = [
    "Supplementary_figure",
    "Root",
    "File",
    "Filename",
    "Keyword_score",
    "Keyword_hits",
    "Rows_excluding_header",
    "Columns",
    "Size_bytes",
    "Header",
]

with open(
    full_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    writer.writeheader()
    writer.writerows(records)


# ============================================================
# Top candidate tables
# ============================================================

top_records = []

for group in GROUPS:

    subset = [
        r for r in records
        if r["Supplementary_figure"] == group
    ]

    subset.sort(
        key=lambda x: (
            -x["Keyword_score"],
            x["Filename"]
        )
    )

    # Top 20 per supplementary figure.
    top_records.extend(
        subset[:20]
    )

top_file = OUTDIR / "10F3B1_top_candidate_source_tables.tsv"

with open(
    top_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    writer.writeheader()
    writer.writerows(top_records)


# ============================================================
# Group summary
# ============================================================

summary_file = OUTDIR / "10F3B1_source_group_summary.tsv"

with open(
    summary_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.writer(
        f,
        delimiter="\t"
    )

    writer.writerow([
        "Supplementary_figure",
        "Candidate_table_N",
        "Top_scored_file",
        "Top_score",
    ])

    for group in GROUPS:

        subset = [
            r for r in records
            if r["Supplementary_figure"] == group
        ]

        if subset:

            subset.sort(
                key=lambda x: (
                    -x["Keyword_score"],
                    x["Filename"]
                )
            )

            top = subset[0]

            writer.writerow([
                group,
                len(subset),
                top["File"],
                top["Keyword_score"],
            ])

        else:

            writer.writerow([
                group,
                0,
                "NA",
                "NA",
            ])


# ============================================================
# Human-readable preview
# ============================================================

preview_file = OUTDIR / "10F3B1_top_candidate_preview.txt"

with open(
    preview_file,
    "w",
    encoding="utf-8"
) as out:

    out.write(
        "============================================================\n"
    )
    out.write(
        "STEP10F3B1 — SUPPLEMENTARY FIGURE SOURCE LOCK PREVIEW\n"
    )
    out.write(
        "============================================================\n\n"
    )

    for group in GROUPS:

        out.write(
            f"==================== {group} ====================\n"
        )

        subset = [
            r for r in records
            if r["Supplementary_figure"] == group
        ]

        subset.sort(
            key=lambda x: (
                -x["Keyword_score"],
                x["Filename"]
            )
        )

        for r in subset[:12]:

            out.write(
                f"\nFILE : {r['File']}\n"
            )
            out.write(
                f"SCORE: {r['Keyword_score']}\n"
            )
            out.write(
                f"HITS : {r['Keyword_hits']}\n"
            )
            out.write(
                f"ROWS : {r['Rows_excluding_header']}\n"
            )
            out.write(
                f"COLS : {r['Columns']}\n"
            )
            out.write(
                f"HEADER:\n{r['Header']}\n"
            )

        out.write("\n")


# ============================================================
# Final console report
# ============================================================

print(
    "============================================================"
)
print(
    "STEP10F3B1 — SUPPLEMENTARY FIGURE SOURCE LOCK"
)
print(
    "Inspection only; no result modified"
)
print(
    "============================================================"
)

for group in GROUPS:

    subset = [
        r for r in records
        if r["Supplementary_figure"] == group
    ]

    print(
        f"{group:35s}: {len(subset):4d} candidate tables"
    )

print()
print(
    "Full inventory:"
)
print(full_file)

print()
print(
    "Top candidates:"
)
print(top_file)

print()
print(
    "Preview:"
)
print(preview_file)

print()
print(
    "No biological analysis was rerun."
)
print(
    "No frozen result was modified."
)
print(
    "No supplementary figure was drawn."
)
print(
    "============================================================"
)
