#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import re

print("=" * 60)
print("STEP10F3B4A — LOCK MOTIF-QC AUTHORITATIVE SOURCES")
print("Exact basename/header classification only")
print("No motif analysis rerun")
print("=" * 60)

OUTDIR = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F3B_Supplementary_Figures/"
    "FigureS4_motif_QC_source_lock"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ============================================================
# 1. Known authoritative roots
# ============================================================

roots = []

for p in [
    Path("08B1_AME"),
    Path("08B2_STREME"),
]:
    if p.exists():
        roots.append(p)

# fallback only for actual Step08B motif-analysis directories
for p in Path(".").glob("08B*AME*"):
    if p.is_dir() and p not in roots:
        roots.append(p)

for p in Path(".").glob("08B*STREME*"):
    if p.is_dir() and p not in roots:
        roots.append(p)

roots = sorted(
    set(roots),
    key=lambda x: str(x)
)

BACKGROUND_SUMMARY = Path(
    "08A_matched_motif_background_v2/"
    "08A_v2_matching_QC_summary.tsv"
)

TF_EVIDENCE_LONG = Path(
    "08C_TF_RNA_integration/"
    "08C2_candidate_motif_TF_evidence_long.tsv"
)

TF_MASTER = Path(
    "08D_motif_family_TF_master/"
    "08D2_publication_TF_candidate_master.tsv"
)


# ============================================================
# 2. Helpers
# ============================================================

def md5sum(path):

    h = hashlib.md5()

    with open(path, "rb") as f:

        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(chunk)

    return h.hexdigest()


def first_nonempty_line(path):

    try:

        with open(
            path,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as f:

            for line in f:

                line = line.rstrip("\r\n")

                if line.strip():
                    return line

    except Exception:
        pass

    return ""


def count_rows_cols(path):

    header = first_nonempty_line(path)

    if not header:
        return 0, 0, ""

    delimiter = "\t"

    if "\t" not in header and "," in header:
        delimiter = ","

    cols = len(
        header.split(delimiter)
    )

    n = 0

    try:

        with open(
            path,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as f:

            for line in f:

                if line.strip():
                    n += 1

    except Exception:

        return -1, cols, header

    return max(0, n - 1), cols, header


# ============================================================
# 3. Strict MEME-result classification
# ============================================================

def classify(path):

    name = path.name.lower()

    # --------------------------------------------------------
    # AME
    # --------------------------------------------------------

    if name in {
        "ame.tsv",
        "ame.txt",
        "ame_results.tsv",
        "ame_results.txt"
    }:
        return "AME"

    if re.fullmatch(
        r"ame([._-].*)?\.(tsv|txt)",
        name
    ):
        return "AME"

    # --------------------------------------------------------
    # STREME
    # --------------------------------------------------------

    if name in {
        "streme.tsv",
        "streme.txt",
        "streme_results.tsv",
        "streme_results.txt"
    }:
        return "STREME"

    if re.fullmatch(
        r"streme([._-].*)?\.(tsv|txt)",
        name
    ):
        return "STREME"

    # --------------------------------------------------------
    # TOMTOM
    # --------------------------------------------------------

    if name in {
        "tomtom.tsv",
        "tomtom.txt",
        "tomtom_results.tsv",
        "tomtom_results.txt"
    }:
        return "Tomtom"

    if re.fullmatch(
        r"tomtom([._-].*)?\.(tsv|txt)",
        name
    ):
        return "Tomtom"

    return None


# ============================================================
# 4. Discover exact MEME outputs
# ============================================================

records = []

for root in roots:

    for path in root.rglob("*"):

        if not path.is_file():
            continue

        if path.suffix.lower() not in {
            ".tsv",
            ".txt"
        }:
            continue

        category = classify(path)

        if category is None:
            continue

        rows, cols, header = count_rows_cols(
            path
        )

        records.append({
            "Category": category,
            "Root": str(root),
            "File": str(path),
            "Rows": rows,
            "Columns": cols,
            "Size_bytes": path.stat().st_size,
            "MD5": md5sum(path),
            "Header": header,
        })


records.sort(
    key=lambda x: (
        x["Category"],
        x["File"]
    )
)


# ============================================================
# 5. Source inventory
# ============================================================

inventory_file = OUTDIR / (
    "10F3B4A_exact_MEME_source_inventory.tsv"
)

fields = [
    "Category",
    "Root",
    "File",
    "Rows",
    "Columns",
    "Size_bytes",
    "MD5",
    "Header",
]

with open(
    inventory_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(records)


# ============================================================
# 6. Auxiliary authoritative sources
# ============================================================

auxiliary = [
    (
        "Matched_background_QC",
        BACKGROUND_SUMMARY
    ),
    (
        "Integrated_motif_TF_evidence",
        TF_EVIDENCE_LONG
    ),
    (
        "Publication_TF_master",
        TF_MASTER
    ),
]

aux_records = []

for category, path in auxiliary:

    if path.exists():

        rows, cols, header = count_rows_cols(
            path
        )

        aux_records.append({
            "Category": category,
            "File": str(path),
            "Exists": "TRUE",
            "Rows": rows,
            "Columns": cols,
            "Size_bytes": path.stat().st_size,
            "MD5": md5sum(path),
            "Header": header,
        })

    else:

        aux_records.append({
            "Category": category,
            "File": str(path),
            "Exists": "FALSE",
            "Rows": "NA",
            "Columns": "NA",
            "Size_bytes": "NA",
            "MD5": "NA",
            "Header": "NA",
        })


aux_file = OUTDIR / (
    "10F3B4A_auxiliary_source_inventory.tsv"
)

with open(
    aux_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Category",
            "File",
            "Exists",
            "Rows",
            "Columns",
            "Size_bytes",
            "MD5",
            "Header",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(aux_records)


# ============================================================
# 7. Category summary
# ============================================================

categories = [
    "AME",
    "STREME",
    "Tomtom"
]

summary_rows = []

for category in categories:

    subset = [
        x for x in records
        if x["Category"] == category
    ]

    summary_rows.append({
        "Category": category,
        "File_N": len(subset),
        "Total_result_rows": sum(
            max(0, int(x["Rows"]))
            for x in subset
        )
    })


summary_file = OUTDIR / (
    "10F3B4A_exact_MEME_category_summary.tsv"
)

with open(
    summary_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Category",
            "File_N",
            "Total_result_rows",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary_rows)


# ============================================================
# 8. Background-QC validation
# ============================================================

qc_rows = []


# Background summary
if BACKGROUND_SUMMARY.exists():

    bg_rows, bg_cols, bg_header = count_rows_cols(
        BACKGROUND_SUMMARY
    )

    qc_rows.append({
        "Metric":
            "Matched_background_summary_exists",
        "Value": 1,
        "Expected": 1,
        "Status": "PASS"
    })

    qc_rows.append({
        "Metric":
            "Matched_background_tissues",
        "Value": bg_rows,
        "Expected": 4,
        "Status":
            "PASS"
            if bg_rows == 4
            else "FAIL"
    })

else:

    qc_rows.append({
        "Metric":
            "Matched_background_summary_exists",
        "Value": 0,
        "Expected": 1,
        "Status": "FAIL"
    })


# Integrated candidate evidence
qc_rows.append({
    "Metric":
        "Integrated_motif_TF_evidence_exists",
    "Value":
        int(TF_EVIDENCE_LONG.exists()),
    "Expected": 1,
    "Status":
        "PASS"
        if TF_EVIDENCE_LONG.exists()
        else "FAIL"
})


# Publication TF master
qc_rows.append({
    "Metric":
        "Publication_TF_master_exists",
    "Value":
        int(TF_MASTER.exists()),
    "Expected": 1,
    "Status":
        "PASS"
        if TF_MASTER.exists()
        else "FAIL"
})


# Actual MEME outputs
for category in categories:

    n = sum(
        1 for x in records
        if x["Category"] == category
    )

    qc_rows.append({
        "Metric":
            f"{category}_exact_result_files",
        "Value": n,
        "Expected":
            ">0",
        "Status":
            "PASS"
            if n > 0
            else "CHECK"
    })


qc_file = OUTDIR / (
    "10F3B4A_source_lock_QC.tsv"
)

with open(
    qc_file,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Metric",
            "Value",
            "Expected",
            "Status",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(qc_rows)


# ============================================================
# 9. Human-readable preview
# ============================================================

preview = OUTDIR / (
    "10F3B4A_exact_MEME_source_preview.txt"
)

with open(
    preview,
    "w",
    encoding="utf-8"
) as out:

    out.write(
        "=" * 60 + "\n"
    )

    out.write(
        "STEP10F3B4A — EXACT MEME SOURCE PREVIEW\n"
    )

    out.write(
        "=" * 60 + "\n\n"
    )

    out.write(
        "AUTHORITATIVE ROOTS\n"
    )

    for root in roots:
        out.write(
            f"{root}\n"
        )

    out.write("\n")

    for category in categories:

        out.write(
            "=" * 20 +
            f" {category} " +
            "=" * 20 +
            "\n"
        )

        subset = [
            x for x in records
            if x["Category"] == category
        ]

        if not subset:

            out.write(
                "NO EXACT RESULT FILE FOUND\n\n"
            )

            continue

        for rec in subset:

            out.write(
                f"\nFILE   : {rec['File']}\n"
            )

            out.write(
                f"ROWS   : {rec['Rows']}\n"
            )

            out.write(
                f"COLUMNS: {rec['Columns']}\n"
            )

            out.write(
                f"HEADER : {rec['Header']}\n"
            )

            # First three data rows
            try:

                with open(
                    rec["File"],
                    "r",
                    encoding="utf-8",
                    errors="replace"
                ) as f:

                    lines = [
                        x.rstrip("\r\n")
                        for x in f
                        if x.strip()
                    ]

                out.write(
                    "FIRST DATA ROWS:\n"
                )

                for line in lines[1:4]:
                    out.write(
                        line + "\n"
                    )

            except Exception:
                pass

        out.write("\n")

    out.write(
        "=" * 20 +
        " AUXILIARY " +
        "=" * 20 +
        "\n\n"
    )

    for rec in aux_records:

        out.write(
            f"{rec['Category']}\n"
        )

        out.write(
            f"FILE  : {rec['File']}\n"
        )

        out.write(
            f"EXISTS: {rec['Exists']}\n"
        )

        out.write(
            f"ROWS  : {rec['Rows']}\n"
        )

        out.write(
            f"HEADER: {rec['Header']}\n\n"
        )


# ============================================================
# 10. Console
# ============================================================

print()
print("Authoritative roots:")

for r in roots:
    print(" ", r)

print()
print("Exact MEME result files:")

for row in summary_rows:

    print(
        f"{row['Category']:10s}"
        f" files={row['File_N']:3d}"
        f" result_rows={row['Total_result_rows']:5d}"
    )

print()
print("Source-lock QC:")

for x in qc_rows:

    print(
        f"{x['Metric']:40s}"
        f" {str(x['Value']):>8s}"
        f" {str(x['Status']):>6s}"
    )

print()
print("Inventory:")
print(inventory_file)

print()
print("Auxiliary sources:")
print(aux_file)

print()
print("Preview:")
print(preview)

print()
print("IMPORTANT:")
print("- Classification uses exact MEME-result filenames.")
print("- 'same' can no longer be misclassified as AME.")
print("- No AME/STREME/Tomtom analysis was rerun.")
print("- No candidate motif was reselected.")
print("- No biological result was modified.")
print("=" * 60)

