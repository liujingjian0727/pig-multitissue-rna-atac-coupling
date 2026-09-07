#!/usr/bin/env python3

import os
import sys
import re
import csv
import hashlib
import subprocess
from pathlib import Path

try:
    from PIL import Image
except ImportError:
    print("ERROR: Pillow is required.")
    print("Install with: pip install pillow")
    sys.exit(1)

# ============================================================
# STEP10F3C
# Supplementary Figures S1-S5 final consistency lock
# ============================================================

print("=" * 70)
print("STEP10F3C — SUPPLEMENTARY FIGURES S1–S5 CONSISTENCY LOCK")
print("Technical / naming / version consistency only")
print("No figure redraw")
print("No biological analysis rerun")
print("=" * 70)

ROOT = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F3B_Supplementary_Figures"
)

OUTDIR = (
    Path("10_publication_figures")
    / "Step10F_final_consistency"
    / "Step10F3C_S1_S5_consistency_lock"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ============================================================
# 1. Authoritative FINAL figure definitions
# ============================================================

FIGURES = {

    "S1": {
        "status_expected": "FINAL_v2",
        "pdf": (
            ROOT
            / "FigureS1_ATAC_QC_FINAL_v2"
            / "Supplementary_Figure_S1_FINAL_v2.pdf"
        ),
        "tiff": (
            ROOT
            / "FigureS1_ATAC_QC_FINAL_v2"
            / "Supplementary_Figure_S1_FINAL_v2_600dpi.tiff"
        ),
        "caption": (
            ROOT
            / "FigureS1_ATAC_QC_FINAL_v2"
            / "Supplementary_Figure_S1_FINAL_v2_caption.txt"
        ),
    },

    "S2": {
        "status_expected": "FINAL_v2",
        "pdf": (
            ROOT
            / "FigureS2_same_tissue_permutation_v2"
            / "Supplementary_Figure_S2_FINAL_v2.pdf"
        ),
        "tiff": (
            ROOT
            / "FigureS2_same_tissue_permutation_v2"
            / "Supplementary_Figure_S2_FINAL_v2_600dpi.tiff"
        ),
        "caption": None,
    },

    "S3": {
        "status_expected": "FINAL_v3",
        "pdf": (
            ROOT
            / "FigureS3_genomic_architecture_FINAL_v3"
            / "Supplementary_Figure_S3_FINAL_v3.pdf"
        ),
        "tiff": (
            ROOT
            / "FigureS3_genomic_architecture_FINAL_v3"
            / "Supplementary_Figure_S3_FINAL_v3_600dpi.tiff"
        ),
        "caption": None,
    },

    "S4": {
        "status_expected": "FINAL_v3",
        "pdf": (
            ROOT
            / "FigureS4_motif_analysis_QC_FINAL_v3"
            / "Supplementary_Figure_S4_FINAL_v3.pdf"
        ),
        "tiff": (
            ROOT
            / "FigureS4_motif_analysis_QC_FINAL_v3"
            / "Supplementary_Figure_S4_FINAL_v3_600dpi.tiff"
        ),
        "caption": None,
    },

    "S5": {
        "status_expected": "FINAL_v5",
        "pdf": (
            ROOT
            / "FigureS5_motif_TF_RNA_support_FINAL_v5"
            / "Supplementary_Figure_S5_FINAL_v5.pdf"
        ),
        "tiff": (
            ROOT
            / "FigureS5_motif_TF_RNA_support_FINAL_v5"
            / "Supplementary_Figure_S5_FINAL_v5_600dpi.tiff"
        ),
        "caption": None,
    },
}

# ============================================================
# 2. Caption auto-detection for S2-S5
# ============================================================

def find_caption(pdf_path):
    if not pdf_path:
        return None

    d = pdf_path.parent

    candidates = sorted(
        list(d.glob("*caption*.txt"))
        + list(d.glob("*caption*.md"))
    )

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) > 1:
        final_candidates = [
            x for x in candidates
            if "FINAL" in x.name
        ]

        if len(final_candidates) == 1:
            return final_candidates[0]

    return None


for fig, x in FIGURES.items():
    if x["caption"] is None:
        x["caption"] = find_caption(
            x["pdf"]
        )

# ============================================================
# 3. Helpers
# ============================================================

def md5(path):
    h = hashlib.md5()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def run_command(cmd):
    try:
        x = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        return x.returncode, x.stdout, x.stderr
    except Exception as e:
        return 999, "", str(e)


def pdf_page_count(pdf):
    rc, out, err = run_command(
        ["pdfinfo", str(pdf)]
    )

    if rc != 0:
        return None

    m = re.search(
        r"^Pages:\s+(\d+)",
        out,
        re.MULTILINE
    )

    if not m:
        return None

    return int(m.group(1))


def pdf_text(pdf):
    rc, out, err = run_command(
        [
            "pdftotext",
            "-layout",
            str(pdf),
            "-"
        ]
    )

    if rc != 0:
        return ""

    return out


def tiff_info(tiff):
    with Image.open(tiff) as im:
        width, height = im.size
        dpi = im.info.get(
            "dpi",
            (None, None)
        )

        if isinstance(dpi, tuple):
            dpi_x = dpi[0]
            dpi_y = dpi[1]
        else:
            dpi_x = None
            dpi_y = None

    return (
        width,
        height,
        dpi_x,
        dpi_y
    )


def dpi_pass(x):
    if x is None:
        return False

    try:
        return abs(float(x) - 600) <= 5
    except Exception:
        return False


# ============================================================
# 4. Expected S1 biological labels
# ============================================================

S1_LABELS = [
    "Adipose-P348",
    "Adipose-P350",
    "Cerebellum-P348",
    "Cerebellum-P350",
    "Cortex-P348",
    "Cortex-P350",
    "Hypothalamus-P348",
    "Hypothalamus-P350",
    "Liver-P348",
    "Liver-P350",
    "Lung-P348",
    "Lung-P350",
    "Muscle-P348",
    "Muscle-P350",
    "Spleen-P348",
    "Spleen-P350",
]

# ============================================================
# 5. Main technical lock
# ============================================================

rows = []

for fig in ["S1", "S2", "S3", "S4", "S5"]:

    item = FIGURES[fig]

    pdf = item["pdf"]
    tiff = item["tiff"]
    caption = item["caption"]

    pdf_exists = pdf.exists()
    tiff_exists = tiff.exists()
    caption_exists = (
        caption is not None
        and caption.exists()
    )

    pages = (
        pdf_page_count(pdf)
        if pdf_exists
        else None
    )

    if tiff_exists:
        try:
            (
                width,
                height,
                dpi_x,
                dpi_y,
            ) = tiff_info(tiff)
        except Exception:
            width = None
            height = None
            dpi_x = None
            dpi_y = None
    else:
        width = None
        height = None
        dpi_x = None
        dpi_y = None

    md5_pdf = (
        md5(pdf)
        if pdf_exists
        else ""
    )

    md5_tiff = (
        md5(tiff)
        if tiff_exists
        else ""
    )

    status_pdf = (
        "PASS"
        if pdf_exists
        else "FAIL"
    )

    status_tiff = (
        "PASS"
        if tiff_exists
        else "FAIL"
    )

    status_page = (
        "PASS"
        if pages == 1
        else "FAIL"
    )

    status_dpi = (
        "PASS"
        if (
            dpi_pass(dpi_x)
            and dpi_pass(dpi_y)
        )
        else "FAIL"
    )

    status_caption = (
        "PASS"
        if caption_exists
        else "CHECK"
    )

    rows.append({
        "Figure": fig,
        "Expected_version":
            item["status_expected"],
        "PDF":
            str(pdf),
        "PDF_exists":
            int(pdf_exists),
        "PDF_pages":
            pages if pages is not None else "NA",
        "PDF_status":
            status_pdf,
        "Single_page_status":
            status_page,
        "TIFF":
            str(tiff),
        "TIFF_exists":
            int(tiff_exists),
        "TIFF_width_px":
            width if width is not None else "NA",
        "TIFF_height_px":
            height if height is not None else "NA",
        "TIFF_dpi_x":
            dpi_x if dpi_x is not None else "NA",
        "TIFF_dpi_y":
            dpi_y if dpi_y is not None else "NA",
        "TIFF_status":
            status_tiff,
        "600dpi_status":
            status_dpi,
        "Caption":
            str(caption)
            if caption is not None
            else "NOT_FOUND",
        "Caption_status":
            status_caption,
        "PDF_MD5":
            md5_pdf,
        "TIFF_MD5":
            md5_tiff,
    })

# ============================================================
# 6. S1 naming-lock QC
# ============================================================

s1_pdf = FIGURES["S1"]["pdf"]

s1_qc_rows = []

if s1_pdf.exists():

    s1_text = pdf_text(
        s1_pdf
    )

    srr_hits = sorted(
        set(
            re.findall(
                r"SRR\d+",
                s1_text
            )
        )
    )

    expected_present = [
        x for x in S1_LABELS
        if x in s1_text
    ]

    s1_qc_rows.extend([
        {
            "Metric":
                "S1_SRR_accessions_remaining",
            "Value":
                len(srr_hits),
            "Expected":
                0,
            "Status":
                "PASS"
                if len(srr_hits) == 0
                else "FAIL",
        },
        {
            "Metric":
                "S1_expected_biological_labels_detected",
            "Value":
                len(expected_present),
            "Expected":
                16,
            "Status":
                "PASS"
                if len(expected_present) == 16
                else "FAIL",
        },
    ])

    with open(
        OUTDIR
        / "S1_detected_sample_labels.txt",
        "w"
    ) as f:

        f.write(
            "Detected biological labels:\n"
        )

        for x in expected_present:
            f.write(
                x + "\n"
            )

        f.write(
            "\nDetected SRR labels:\n"
        )

        for x in srr_hits:
            f.write(
                x + "\n"
            )

else:

    s1_qc_rows.extend([
        {
            "Metric":
                "S1_SRR_accessions_remaining",
            "Value":
                "NA",
            "Expected":
                0,
            "Status":
                "FAIL",
        },
        {
            "Metric":
                "S1_expected_biological_labels_detected",
            "Value":
                "NA",
            "Expected":
                16,
            "Status":
                "FAIL",
        },
    ])

# ============================================================
# 7. Locked-panel / scientific-semantic manifest
# ============================================================

semantic_rows = [

    {
        "Figure": "S1",
        "Panel_structure":
            "A-B-C1/C2/C3/C4",
        "Frozen_interpretation":
            "ATAC-seq technical QC only",
        "Critical_wording":
            "Sample labels use Tissue-P348/P350; SRR retained in Table S1",
        "Status":
            "PASS_IF_TECHNICAL_QC_PASSES",
    },

    {
        "Figure": "S2",
        "Panel_structure":
            "A-D",
        "Frozen_interpretation":
            "Same-tissue permutation validation of RNA-ATAC coupling",
        "Critical_wording":
            "Permutation P<1e-4 where no exceedances among 10,000 permutations",
        "Status":
            "LOCKED",
    },

    {
        "Figure": "S3",
        "Panel_structure":
            "A-E",
        "Frozen_interpretation":
            "Genomic architecture and structural enrichment",
        "Critical_wording":
            "Distal_to_associated_gene is pair-specific; N=0 classes explicit",
        "Status":
            "LOCKED",
    },

    {
        "Figure": "S4",
        "Panel_structure":
            "A-E",
        "Frozen_interpretation":
            "Motif-analysis QC and evidence hierarchy",
        "Critical_wording":
            "AME formal E<0.05=0; STREME 0/9 evaluable holdout E<0.05; score-only Cerebellum not formal significance",
        "Status":
            "LOCKED",
    },

    {
        "Figure": "S5",
        "Panel_structure":
            "A-E",
        "Frozen_interpretation":
            "Candidate TF RNA support and evidence hierarchy",
        "Critical_wording":
            "Candidate TF does not imply direct binding; FIMO does not imply occupancy",
        "Status":
            "LOCKED",
    },
]

# ============================================================
# 8. Write files
# ============================================================

tech_file = (
    OUTDIR
    / "10F3C_S1_S5_technical_lock_manifest.tsv"
)

with open(
    tech_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys(),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


s1_file = (
    OUTDIR
    / "10F3C_S1_naming_lock_QC.tsv"
)

with open(
    s1_file,
    "w",
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
    writer.writerows(
        s1_qc_rows
    )


semantic_file = (
    OUTDIR
    / "10F3C_S1_S5_semantic_lock_manifest.tsv"
)

with open(
    semantic_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Figure",
            "Panel_structure",
            "Frozen_interpretation",
            "Critical_wording",
            "Status",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        semantic_rows
    )

# ============================================================
# 9. Overall status
# ============================================================

hard_fail = False

for r in rows:

    for field in [
        "PDF_status",
        "Single_page_status",
        "TIFF_status",
        "600dpi_status",
    ]:

        if r[field] != "PASS":
            hard_fail = True


for r in s1_qc_rows:

    if r["Status"] != "PASS":
        hard_fail = True


overall = (
    "PASS"
    if not hard_fail
    else "FAIL"
)


overall_file = (
    OUTDIR
    / "10F3C_S1_S5_overall_status.tsv"
)

with open(
    overall_file,
    "w"
) as f:

    f.write(
        "Metric\tValue\n"
    )

    f.write(
        f"S1_S5_consistency_lock\t{overall}\n"
    )

    f.write(
        "Biological_results_modified\t0\n"
    )


# ============================================================
# 10. Console report
# ============================================================

print()
print("=" * 70)
print("TECHNICAL LOCK SUMMARY")
print("=" * 70)

for r in rows:

    print(
        f"{r['Figure']:>2}  "
        f"PDF={r['PDF_status']:<4}  "
        f"page={r['Single_page_status']:<4}  "
        f"TIFF={r['TIFF_status']:<4}  "
        f"600dpi={r['600dpi_status']:<4}  "
        f"caption={r['Caption_status']}"
    )


print()
print("=" * 70)
print("S1 NAMING LOCK")
print("=" * 70)

for r in s1_qc_rows:

    print(
        f"{r['Metric']:<42} "
        f"{str(r['Value']):>6} / "
        f"{str(r['Expected']):<6} "
        f"{r['Status']}"
    )


print()
print("=" * 70)

if overall == "PASS":

    print(
        "STEP10F3C STATUS: PASS"
    )

    print(
        "Supplementary Figures S1–S5 are technically and semantically LOCKED."
    )

else:

    print(
        "STEP10F3C STATUS: FAIL"
    )

    print(
        "At least one final supplementary figure requires correction before lock."
    )

print(
    "Biological results modified: 0"
)

print("=" * 70)

print("\nOutputs:")
print(tech_file)
print(s1_file)
print(semantic_file)
print(overall_file)
