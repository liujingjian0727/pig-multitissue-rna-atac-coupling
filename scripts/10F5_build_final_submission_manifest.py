#!/usr/bin/env python3

from pathlib import Path
import hashlib
import csv
import re
import subprocess
import sys

try:
    from PIL import Image
except Exception as e:
    raise SystemExit(
        "ERROR: Pillow is required.\n"
        f"{e}"
    )

print("=" * 96)
print("STEP10F5 — FINAL MAIN + SUPPLEMENTARY SUBMISSION MANIFEST")
print("Figures 1–6 + Supplementary Figures S1–S5 + Supplementary Tables S1–S8")
print("No figure redraw")
print("No table reconstruction")
print("No biological analysis rerun")
print("=" * 96)

ROOT = Path(".")

OUTDIR = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F5_final_submission_manifest"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================================
# 1. FINAL FROZEN ASSETS
# ======================================================================

ASSETS = [

# ----------------------------------------------------------------------
# MAIN FIGURES
# ----------------------------------------------------------------------

{
    "ID": "Figure1",
    "Category": "Main_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10B_Figure1_2_v2/Figure1/"
        "Figure1_FINAL_v2.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10B_Figure1_2_v2/Figure1/"
        "Figure1_FINAL_v2_600dpi.tiff",
},

{
    "ID": "Figure2",
    "Category": "Main_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10preF_Figure2_FINAL/"
        "Figure2_FINAL.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10preF_Figure2_FINAL/"
        "Figure2_FINAL_600dpi.tiff",
},

{
    "ID": "Figure3",
    "Category": "Main_Figure",
    "Panel_count": 6,
    "Panel_labels": "A-F",
    "PDF":
        "10_publication_figures/"
        "Step10C_Figure3_4_v3/Figure3/"
        "Figure3_FINAL_v3.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10C_Figure3_4_v3/Figure3/"
        "Figure3_FINAL_v3_600dpi.tiff",
},

{
    "ID": "Figure4",
    "Category": "Main_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10C_Figure3_4_v3/Figure4/"
        "Figure4_FINAL_v3.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10C_Figure3_4_v3/Figure4/"
        "Figure4_FINAL_v3_600dpi.tiff",
},

{
    "ID": "Figure5",
    "Category": "Main_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10preF_Figure5_FINAL/"
        "Figure5_FINAL.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10preF_Figure5_FINAL/"
        "Figure5_FINAL_600dpi.tiff",
},

{
    "ID": "Figure6",
    "Category": "Main_Figure",
    "Panel_count": 3,
    "Panel_labels": "A-C",
    "PDF":
        "10_publication_figures/"
        "Step10E_Figure6_FINAL/"
        "Figure6_FINAL.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10E_Figure6_FINAL/"
        "Figure6_FINAL_600dpi.tiff",
},


# ----------------------------------------------------------------------
# SUPPLEMENTARY FIGURES
# ----------------------------------------------------------------------

{
    "ID": "FigureS1",
    "Category": "Supplementary_Figure",
    "Panel_count": 6,
    "Panel_labels": "A,B,C1-C4",
    "PDF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS1_ATAC_QC_FINAL_v2/"
        "Supplementary_Figure_S1_FINAL_v2.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS1_ATAC_QC_FINAL_v2/"
        "Supplementary_Figure_S1_FINAL_v2_600dpi.tiff",
},

{
    "ID": "FigureS2",
    "Category": "Supplementary_Figure",
    "Panel_count": 4,
    "Panel_labels": "A-D",
    "PDF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS2_same_tissue_permutation_v2/"
        "Supplementary_Figure_S2_FINAL_v2.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS2_same_tissue_permutation_v2/"
        "Supplementary_Figure_S2_FINAL_v2_600dpi.tiff",
},

{
    "ID": "FigureS3",
    "Category": "Supplementary_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS3_genomic_architecture_FINAL_v3/"
        "Supplementary_Figure_S3_FINAL_v3.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS3_genomic_architecture_FINAL_v3/"
        "Supplementary_Figure_S3_FINAL_v3_600dpi.tiff",
},

{
    "ID": "FigureS4",
    "Category": "Supplementary_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS4_motif_analysis_QC_FINAL_v3/"
        "Supplementary_Figure_S4_FINAL_v3.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS4_motif_analysis_QC_FINAL_v3/"
        "Supplementary_Figure_S4_FINAL_v3_600dpi.tiff",
},

{
    "ID": "FigureS5",
    "Category": "Supplementary_Figure",
    "Panel_count": 5,
    "Panel_labels": "A-E",
    "PDF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS5_motif_TF_RNA_support_FINAL_v5/"
        "Supplementary_Figure_S5_FINAL_v5.pdf",
    "TIFF":
        "10_publication_figures/"
        "Step10F_final_consistency/"
        "Step10F3B_Supplementary_Figures/"
        "FigureS5_motif_TF_RNA_support_FINAL_v5/"
        "Supplementary_Figure_S5_FINAL_v5_600dpi.tiff",
},
]


# ----------------------------------------------------------------------
# SUPPLEMENTARY TABLES
# ----------------------------------------------------------------------

TABLE_DIR = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F4_Supplementary_Tables/"
    "Step10F4B_final_workbooks"
)

TABLES = [
    {
        "ID": f"TableS{i}",
        "Category": "Supplementary_Table",
        "XLSX":
            str(
                TABLE_DIR /
                f"Supplementary_Table_S{i}.xlsx"
            )
    }
    for i in range(1, 9)
]


# ======================================================================
# 2. LOCK FILES
# ======================================================================

SUPP_FIG_LOCK = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F3C_S1_S5_consistency_lock/"
    "10F3C_S1_S5_overall_status.tsv"
)

SUPP_TABLE_LOCK = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F4_Supplementary_Tables/"
    "Step10F4C_consistency_lock/"
    "10F4C_overall_status.tsv"
)

SUPP_TABLE_LOCK_MARKER = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F4_Supplementary_Tables/"
    "Step10F4C_consistency_lock/"
    "SUPPLEMENTARY_TABLES_S1_S8_LOCKED.txt"
)


# ======================================================================
# 3. HELPERS
# ======================================================================

def md5(path):

    h = hashlib.md5()

    with open(path, "rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def run(cmd):

    return subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False
    )


def pdf_pages(path):

    x = run(
        [
            "pdfinfo",
            str(path)
        ]
    )

    if x.returncode != 0:
        return None

    m = re.search(
        r"^Pages:\s+(\d+)",
        x.stdout,
        re.MULTILINE
    )

    if not m:
        return None

    return int(
        m.group(1)
    )


def pdf_size_points(path):

    x = run(
        [
            "pdfinfo",
            str(path)
        ]
    )

    if x.returncode != 0:
        return None, None

    m = re.search(
        r"Page size:\s+"
        r"([0-9.]+)\s+x\s+([0-9.]+)\s+pts",
        x.stdout
    )

    if not m:
        return None, None

    return (
        float(m.group(1)),
        float(m.group(2))
    )


def tiff_info(path):

    with Image.open(path) as im:

        width, height = im.size

        dpi = im.info.get(
            "dpi",
            (None, None)
        )

        if isinstance(
            dpi,
            tuple
        ):
            dx, dy = dpi
        else:
            dx = dy = None

    return (
        width,
        height,
        dx,
        dy
    )


def dpi_pass(x):

    try:
        return (
            abs(
                float(x) - 600
            )
            <= 5
        )
    except Exception:
        return False


def read_status_tsv(path):

    if not path.exists():
        return {}

    d = {}

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as f:

        reader = csv.reader(
            f,
            delimiter="\t"
        )

        header = next(
            reader,
            None
        )

        for row in reader:

            if len(row) >= 2:
                d[
                    row[0]
                ] = row[1]

    return d


# ======================================================================
# 4. VERIFY LOCKS
# ======================================================================

fig_lock_data = read_status_tsv(
    SUPP_FIG_LOCK
)

table_lock_data = read_status_tsv(
    SUPP_TABLE_LOCK
)

fig_lock_pass = (
    fig_lock_data.get(
        "S1_S5_consistency_lock"
    ) == "PASS"
)

table_lock_pass = (
    table_lock_data.get(
        "Supplementary_Tables_S1_S8_consistency_lock"
    ) == "PASS"
)

table_marker_pass = (
    SUPP_TABLE_LOCK_MARKER.exists()
)


# ======================================================================
# 5. FIGURE QC
# ======================================================================

figure_rows = []

hard_fail = False

for x in ASSETS:

    pdf = ROOT / x["PDF"]
    tiff = ROOT / x["TIFF"]

    pdf_exists = pdf.exists()
    tiff_exists = tiff.exists()

    pages = (
        pdf_pages(pdf)
        if pdf_exists
        else None
    )

    pw, ph = (
        pdf_size_points(pdf)
        if pdf_exists
        else (
            None,
            None
        )
    )

    if tiff_exists:

        try:
            (
                tw,
                th,
                dx,
                dy
            ) = tiff_info(
                tiff
            )

        except Exception:

            tw = th = dx = dy = None

    else:

        tw = th = dx = dy = None


    pdf_status = (
        "PASS"
        if pdf_exists
        else "FAIL"
    )

    page_status = (
        "PASS"
        if pages == 1
        else "FAIL"
    )

    tiff_status = (
        "PASS"
        if tiff_exists
        else "FAIL"
    )

    dpi_status = (
        "PASS"
        if (
            dpi_pass(dx)
            and
            dpi_pass(dy)
        )
        else "FAIL"
    )


    for status in [
        pdf_status,
        page_status,
        tiff_status,
        dpi_status
    ]:

        if status != "PASS":
            hard_fail = True


    figure_rows.append({
        "ID":
            x["ID"],

        "Category":
            x["Category"],

        "Panel_count":
            x["Panel_count"],

        "Panel_labels":
            x["Panel_labels"],

        "PDF":
            str(pdf),

        "PDF_exists":
            int(pdf_exists),

        "PDF_pages":
            pages
            if pages is not None
            else "NA",

        "PDF_width_pt":
            pw
            if pw is not None
            else "NA",

        "PDF_height_pt":
            ph
            if ph is not None
            else "NA",

        "PDF_MD5":
            md5(pdf)
            if pdf_exists
            else "",

        "TIFF":
            str(tiff),

        "TIFF_exists":
            int(tiff_exists),

        "TIFF_width_px":
            tw
            if tw is not None
            else "NA",

        "TIFF_height_px":
            th
            if th is not None
            else "NA",

        "TIFF_dpi_x":
            dx
            if dx is not None
            else "NA",

        "TIFF_dpi_y":
            dy
            if dy is not None
            else "NA",

        "TIFF_MD5":
            md5(tiff)
            if tiff_exists
            else "",

        "PDF_status":
            pdf_status,

        "Single_page_status":
            page_status,

        "TIFF_status":
            tiff_status,

        "600dpi_status":
            dpi_status
    })


# ======================================================================
# 6. TABLE QC
# ======================================================================

table_rows = []

for x in TABLES:

    p = ROOT / x["XLSX"]

    exists = p.exists()

    status = (
        "PASS"
        if exists
        else "FAIL"
    )

    if status != "PASS":
        hard_fail = True

    table_rows.append({
        "ID":
            x["ID"],

        "Category":
            x["Category"],

        "Workbook":
            str(p),

        "Exists":
            int(exists),

        "Size_bytes":
            (
                p.stat().st_size
                if exists
                else "NA"
            ),

        "MD5":
            (
                md5(p)
                if exists
                else ""
            ),

        "Status":
            status
    })


# ======================================================================
# 7. MAIN FIGURE PANEL COUNT CHECK
# ======================================================================

main_panel_count = sum(
    x["Panel_count"]
    for x in ASSETS
    if x["Category"]
    == "Main_Figure"
)

main_panel_status = (
    "PASS"
    if main_panel_count == 29
    else "FAIL"
)

if main_panel_status != "PASS":
    hard_fail = True


# ======================================================================
# 8. LOCK CHECK
# ======================================================================

lock_checks = [

    {
        "Check":
            "Supplementary_Figures_S1_S5_lock",
        "Observed":
            "PASS"
            if fig_lock_pass
            else "FAIL",
        "Expected":
            "PASS",
        "Status":
            "PASS"
            if fig_lock_pass
            else "FAIL"
    },

    {
        "Check":
            "Supplementary_Tables_S1_S8_lock",
        "Observed":
            "PASS"
            if table_lock_pass
            else "FAIL",
        "Expected":
            "PASS",
        "Status":
            "PASS"
            if table_lock_pass
            else "FAIL"
    },

    {
        "Check":
            "Supplementary_Table_lock_marker",
        "Observed":
            int(
                table_marker_pass
            ),
        "Expected":
            1,
        "Status":
            "PASS"
            if table_marker_pass
            else "FAIL"
    },

    {
        "Check":
            "Main_Figure_panel_count",
        "Observed":
            main_panel_count,
        "Expected":
            29,
        "Status":
            main_panel_status
    }
]


for x in lock_checks:

    if x["Status"] != "PASS":
        hard_fail = True


# ======================================================================
# 9. WRITE MANIFESTS
# ======================================================================

figure_file = (
    OUTDIR /
    "10F5_final_figure_manifest.tsv"
)

with open(
    figure_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=
            list(
                figure_rows[0].keys()
            ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        figure_rows
    )


table_file = (
    OUTDIR /
    "10F5_final_table_manifest.tsv"
)

with open(
    table_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=
            list(
                table_rows[0].keys()
            ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        table_rows
    )


lock_file = (
    OUTDIR /
    "10F5_lock_status.tsv"
)

with open(
    lock_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Check",
            "Observed",
            "Expected",
            "Status"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        lock_checks
    )


# ======================================================================
# 10. UNIFIED SUBMISSION FILE MANIFEST
# ======================================================================

all_file = (
    OUTDIR /
    "10F5_UNIFIED_submission_file_manifest.tsv"
)

with open(
    all_file,
    "w",
    newline=""
) as f:

    fields = [
        "Asset_ID",
        "Category",
        "Format",
        "Path",
        "MD5",
        "Status"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()


    for x in figure_rows:

        writer.writerow({
            "Asset_ID":
                x["ID"],
            "Category":
                x["Category"],
            "Format":
                "PDF",
            "Path":
                x["PDF"],
            "MD5":
                x["PDF_MD5"],
            "Status":
                x["PDF_status"]
        })

        writer.writerow({
            "Asset_ID":
                x["ID"],
            "Category":
                x["Category"],
            "Format":
                "TIFF_600dpi",
            "Path":
                x["TIFF"],
            "MD5":
                x["TIFF_MD5"],
            "Status":
                (
                    "PASS"
                    if (
                        x["TIFF_status"]
                        == "PASS"
                        and
                        x["600dpi_status"]
                        == "PASS"
                    )
                    else "FAIL"
                )
        })


    for x in table_rows:

        writer.writerow({
            "Asset_ID":
                x["ID"],
            "Category":
                x["Category"],
            "Format":
                "XLSX",
            "Path":
                x["Workbook"],
            "MD5":
                x["MD5"],
            "Status":
                x["Status"]
        })


# ======================================================================
# 11. OVERALL STATUS
# ======================================================================

figure_pass = sum(
    (
        x["PDF_status"] == "PASS"
        and
        x["Single_page_status"] == "PASS"
        and
        x["TIFF_status"] == "PASS"
        and
        x["600dpi_status"] == "PASS"
    )
    for x in figure_rows
)

table_pass = sum(
    x["Status"] == "PASS"
    for x in table_rows
)

overall = (
    "FAIL"
    if hard_fail
    else "PASS"
)


status_file = (
    OUTDIR /
    "10F5_overall_status.tsv"
)

with open(
    status_file,
    "w"
) as f:

    f.write(
        "Metric\tValue\n"
    )

    f.write(
        "Main_figures_expected\t6\n"
    )

    f.write(
        "Supplementary_figures_expected\t5\n"
    )

    f.write(
        "Figure_pairs_verified\t"
        f"{figure_pass}\n"
    )

    f.write(
        "Figure_pairs_expected\t11\n"
    )

    f.write(
        "Supplementary_tables_verified\t"
        f"{table_pass}\n"
    )

    f.write(
        "Supplementary_tables_expected\t8\n"
    )

    f.write(
        "Main_figure_total_panels\t"
        f"{main_panel_count}\n"
    )

    f.write(
        "Supplementary_figures_lock\t"
        f"{'PASS' if fig_lock_pass else 'FAIL'}\n"
    )

    f.write(
        "Supplementary_tables_lock\t"
        f"{'PASS' if table_lock_pass else 'FAIL'}\n"
    )

    f.write(
        "Final_submission_manifest\t"
        f"{overall}\n"
    )

    f.write(
        "Biological_results_modified\t0\n"
    )


# ======================================================================
# 12. FINAL LOCK MARKER
# ======================================================================

marker = (
    OUTDIR /
    "FINAL_FIGURE_TABLE_SUBMISSION_ASSETS_LOCKED.txt"
)

if overall == "PASS":

    with open(
        marker,
        "w"
    ) as f:

        f.write(
            "FINAL FIGURE/TABLE SUBMISSION ASSETS ARE LOCKED.\n"
        )

        f.write(
            "Main Figures: 1-6\n"
        )

        f.write(
            "Supplementary Figures: S1-S5\n"
        )

        f.write(
            "Supplementary Tables: S1-S8\n"
        )

        f.write(
            "Main figure panels: 29\n"
        )

        f.write(
            "Figure PDF/TIFF pairs verified: 11/11\n"
        )

        f.write(
            "Supplementary Tables verified: 8/8\n"
        )

        f.write(
            "Biological results modified: 0\n"
        )

        f.write(
            "Step10F5 final submission manifest: PASS\n"
        )

else:

    if marker.exists():
        marker.unlink()


# ======================================================================
# 13. CONSOLE
# ======================================================================

print()
print("=" * 96)
print("FIGURE MANIFEST QC")
print("=" * 96)

for x in figure_rows:

    overall_figure = (
        "PASS"
        if (
            x["PDF_status"]
            == "PASS"
            and
            x["Single_page_status"]
            == "PASS"
            and
            x["TIFF_status"]
            == "PASS"
            and
            x["600dpi_status"]
            == "PASS"
        )
        else "FAIL"
    )

    print(
        f"{x['ID']:<9} "
        f"panels={str(x['Panel_count']):>2}  "
        f"PDF={x['PDF_status']}  "
        f"page={x['Single_page_status']}  "
        f"TIFF={x['TIFF_status']}  "
        f"600dpi={x['600dpi_status']}  "
        f"{overall_figure}"
    )


print()
print("=" * 96)
print("SUPPLEMENTARY TABLE MANIFEST QC")
print("=" * 96)

for x in table_rows:

    print(
        f"{x['ID']:<9} "
        f"{x['Status']}"
    )


print()
print("=" * 96)
print("FINAL LOCK CHECKS")
print("=" * 96)

for x in lock_checks:

    print(
        f"{x['Check']:<42} "
        f"{str(x['Observed']):>6} / "
        f"{str(x['Expected']):<6} "
        f"{x['Status']}"
    )


print()
print("=" * 96)
print(
    f"STEP10F5 STATUS: {overall}"
)

if overall == "PASS":

    print(
        "Main Figures 1–6, Supplementary Figures S1–S5, "
        "and Supplementary Tables S1–S8 are FINAL and LOCKED."
    )

    print(
        "Ready for manuscript-level Results/Methods/Discussion integration."
    )

else:

    print(
        "Submission assets are NOT fully locked."
    )

    print(
        "Resolve failed checks before manuscript finalization."
    )

print(
    "Biological results modified: 0"
)

print("=" * 96)

print("\nOutputs:")
print(figure_file)
print(table_file)
print(lock_file)
print(all_file)
print(status_file)

if overall == "PASS":
    print(marker)
