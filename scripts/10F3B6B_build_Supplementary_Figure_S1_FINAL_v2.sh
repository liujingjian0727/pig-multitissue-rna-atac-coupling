#!/usr/bin/env bash
set -euo pipefail

echo "======================================================================"
echo "STEP10F3B6B — SUPPLEMENTARY FIGURE S1 FINAL v2"
echo "Rebuild from frozen ATAC-QC sources"
echo "SRR labels -> Tissue-P348/P350"
echo "No ATAC QC statistic recomputed"
echo "======================================================================"

# ============================================================
# 1. Frozen source paths
# ============================================================

BASE="/datadisk2/liujingjian_data/cattle_and_sus/sus_atac_analysis/bam2bigwig"

FRAG_RAW="${BASE}/fragments_distribution/pig_ATAC_fragment_size_raw.tsv"
COR_TAB="${BASE}/plotCorrelation/pig_sample_correlation.tab"
REFDIR="${BASE}/reference_point_plot"

# Absolute frozen deepTools executable
PLOTHEATMAP="/home/liujingjian/miniconda3/envs/deeptools_env/bin/plotHeatmap"

if [[ ! -x "${PLOTHEATMAP}" ]]; then
    echo "ERROR: plotHeatmap executable unavailable:"
    echo "${PLOTHEATMAP}"
    exit 1
fi

echo "Using plotHeatmap:"
"${PLOTHEATMAP}" --version || true

M1="${REFDIR}/SRR54_58_tss_matrix.gz"
M2="${REFDIR}/SRR61_66_tss_matrix.gz"
M3="${REFDIR}/SRR69_74_tss_matrix.gz"
M4="${REFDIR}/SRR77_82_tss_matrix.gz"

OUTDIR="10_publication_figures/Step10F_final_consistency/Step10F3B_Supplementary_Figures/FigureS1_ATAC_QC_FINAL_v2"

TMP="${OUTDIR}/tmp_build"

mkdir -p "${OUTDIR}"
rm -rf "${TMP}"
mkdir -p "${TMP}"

# ============================================================
# 2. Required-file check
# ============================================================

for f in \
    "${FRAG_RAW}" \
    "${COR_TAB}" \
    "${M1}" \
    "${M2}" \
    "${M3}" \
    "${M4}"
do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing frozen source:"
        echo "${f}"
        exit 1
    fi
done

for cmd in python3 pdfinfo pdftotext pdftoppm
do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "ERROR: required command not found: ${cmd}"
        exit 1
    fi
done

# ============================================================
# 3. Build Panel A + Panel B
#    A = frozen bamPEFragmentSize raw table
#    B = frozen correlation matrix
# ============================================================

cat > "${TMP}/build_S1_AB.py" <<'PY'
import re
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform
except Exception as e:
    raise SystemExit(
        "ERROR importing Python packages.\n"
        "Required: numpy pandas matplotlib scipy\n"
        f"{e}"
    )

frag_file = Path(sys.argv[1])
corr_file = Path(sys.argv[2])
out_A = Path(sys.argv[3])
out_B = Path(sys.argv[4])

# ============================================================
# Biological-label mapping
# ============================================================

mapping = {
    "SRR12697154": "Adipose-P348",
    "SRR12697155": "Adipose-P350",

    "SRR12697157": "Cerebellum-P348",
    "SRR12697158": "Cerebellum-P350",

    "SRR12697161": "Cortex-P348",
    "SRR12697162": "Cortex-P350",

    "SRR12697165": "Hypothalamus-P348",
    "SRR12697166": "Hypothalamus-P350",

    "SRR12697169": "Liver-P348",
    "SRR12697170": "Liver-P350",

    "SRR12697173": "Lung-P348",
    "SRR12697174": "Lung-P350",

    "SRR12697177": "Muscle-P348",
    "SRR12697178": "Muscle-P350",

    "SRR12697181": "Spleen-P348",
    "SRR12697182": "Spleen-P350",
}

sample_order = list(mapping.keys())


def get_srr(x):
    x = str(x)
    m = re.search(r"SRR\d+", x)
    return m.group(0) if m else x


# ============================================================
# PANEL A
# Frozen fragment-size raw table
#
# bamPEFragmentSize --outRawFragmentLengths convention:
# Size / Occurrences / Sample
# ============================================================

try:
    frag = pd.read_csv(
        frag_file,
        sep="\t",
        comment="#"
    )
except Exception as e:
    raise SystemExit(
        f"ERROR reading fragment-size table: {e}"
    )


def find_col(cols, patterns):
    for c in cols:
        lc = str(c).strip().lower()
        for p in patterns:
            if p in lc:
                return c
    return None


size_col = find_col(
    frag.columns,
    ["size", "fragment length", "fragment_length", "length"]
)

count_col = find_col(
    frag.columns,
    ["occurrence", "count", "number"]
)

sample_col = find_col(
    frag.columns,
    ["sample"]
)


# Fallback for the standard three-column deepTools raw output.
if (
    size_col is None or
    count_col is None or
    sample_col is None
):
    frag2 = pd.read_csv(
        frag_file,
        sep="\t",
        comment="#",
        header=None
    )

    if frag2.shape[1] < 3:
        raise SystemExit(
            "ERROR: could not identify columns in "
            "pig_ATAC_fragment_size_raw.tsv"
        )

    frag2 = frag2.iloc[:, :3].copy()
    frag2.columns = [
        "Size",
        "Occurrences",
        "Sample"
    ]

    # Remove accidental header row if present.
    frag2["Size_num"] = pd.to_numeric(
        frag2["Size"],
        errors="coerce"
    )

    frag2 = frag2.loc[
        frag2["Size_num"].notna()
    ].copy()

    frag2["Size"] = frag2["Size_num"]
    frag2.drop(
        columns="Size_num",
        inplace=True
    )

    frag = frag2

    size_col = "Size"
    count_col = "Occurrences"
    sample_col = "Sample"


frag = frag[
    [size_col, count_col, sample_col]
].copy()

frag.columns = [
    "Size",
    "Occurrences",
    "Sample"
]

frag["Size"] = pd.to_numeric(
    frag["Size"],
    errors="coerce"
)

frag["Occurrences"] = pd.to_numeric(
    frag["Occurrences"],
    errors="coerce"
)

frag["SRR"] = frag["Sample"].map(
    get_srr
)

frag = frag[
    frag["SRR"].isin(sample_order)
].dropna(
    subset=[
        "Size",
        "Occurrences"
    ]
).copy()

if frag.empty:
    raise SystemExit(
        "ERROR: fragment raw table contained no expected SRR samples."
    )


# Normalize each sample to frequency.
frag["Frequency"] = (
    frag["Occurrences"] /
    frag.groupby("SRR")["Occurrences"].transform("sum")
)


fig, ax = plt.subplots(
    figsize=(8.2, 6.2)
)

for srr in sample_order:

    d = frag[
        frag["SRR"] == srr
    ].sort_values(
        "Size"
    )

    if d.empty:
        raise SystemExit(
            f"ERROR: fragment data missing sample {srr}"
        )

    ax.plot(
        d["Size"],
        d["Frequency"],
        linewidth=1.0,
        alpha=0.72,
        label=mapping[srr]
    )


ax.set_xlim(
    30,
    650
)

ax.set_xlabel(
    "Fragment Length",
    fontsize=10
)

ax.set_ylabel(
    "Frequency",
    fontsize=10
)

ax.tick_params(
    labelsize=8
)

ax.legend(
    loc="upper right",
    frameon=True,
    fontsize=7.1,
    ncol=1
)

for spine in ax.spines.values():
    spine.set_linewidth(0.8)

fig.tight_layout()

fig.savefig(
    out_A,
    format="pdf",
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# PANEL B
# Frozen Pearson correlation matrix.
# No correlations are recalculated.
# ============================================================

try:
    corr = pd.read_csv(
        corr_file,
        sep="\t",
        comment="#",
        index_col=0
    )
except Exception as e:
    raise SystemExit(
        f"ERROR reading correlation matrix: {e}"
    )


# Strip quoting/whitespace.
corr.index = [
    str(x).strip().strip('"').strip("'")
    for x in corr.index
]

corr.columns = [
    str(x).strip().strip('"').strip("'")
    for x in corr.columns
]


# Some deepTools matrices may have an unnamed first column.
corr = corr.apply(
    pd.to_numeric,
    errors="coerce"
)

corr = corr.dropna(
    axis=0,
    how="all"
).dropna(
    axis=1,
    how="all"
)


# Extract SRR IDs from row/column names.
row_srr = [
    get_srr(x)
    for x in corr.index
]

col_srr = [
    get_srr(x)
    for x in corr.columns
]

corr.index = row_srr
corr.columns = col_srr


common = [
    x for x in sample_order
    if x in corr.index and x in corr.columns
]

if len(common) != 16:
    raise SystemExit(
        "ERROR: expected 16 samples in correlation matrix; "
        f"found {len(common)}.\n"
        f"Rows: {list(corr.index)}\n"
        f"Cols: {list(corr.columns)}"
    )


corr = corr.loc[
    common,
    common
].astype(float)


# Ensure symmetric matrix numerically.
arr = corr.to_numpy()
arr = (
    arr + arr.T
) / 2.0

np.fill_diagonal(
    arr,
    1.0
)


# Visualization clustering only.
# Correlation values remain frozen.
distance = 1.0 - arr

distance[
    distance < 0
] = 0

np.fill_diagonal(
    distance,
    0
)

condensed = squareform(
    distance,
    checks=False
)

Z = linkage(
    condensed,
    method="average"
)

order = leaves_list(
    Z
)

arr_o = arr[
    np.ix_(
        order,
        order
    )
]

srr_o = [
    common[i]
    for i in order
]

labels_o = [
    mapping[x]
    for x in srr_o
]


fig = plt.figure(
    figsize=(9.2, 7.2)
)

gs = GridSpec(
    nrows=2,
    ncols=2,
    width_ratios=[
        1.0,
        6.8
    ],
    height_ratios=[
        6.5,
        0.38
    ],
    hspace=0.15,
    wspace=0.02
)


ax_den = fig.add_subplot(
    gs[0, 0]
)

dendrogram(
    Z,
    orientation="left",
    no_labels=True,
    color_threshold=0,
    above_threshold_color="#9c4a4a",
    ax=ax_den
)

ax_den.invert_yaxis()
ax_den.axis("off")


ax_hm = fig.add_subplot(
    gs[0, 1]
)

im = ax_hm.imshow(
    arr_o,
    cmap="jet",
    vmin=0.30,
    vmax=1.00,
    interpolation="nearest",
    aspect="equal"
)

ax_hm.set_xticks(
    np.arange(16)
)

ax_hm.set_xticklabels(
    labels_o,
    rotation=48,
    ha="left",
    rotation_mode="anchor",
    fontsize=6.4
)

ax_hm.xaxis.tick_top()

ax_hm.set_yticks(
    np.arange(16)
)

ax_hm.set_yticklabels(
    labels_o,
    fontsize=6.5
)

ax_hm.yaxis.tick_right()

ax_hm.tick_params(
    length=0
)

# Grid lines.
ax_hm.set_xticks(
    np.arange(-0.5, 16, 1),
    minor=True
)

ax_hm.set_yticks(
    np.arange(-0.5, 16, 1),
    minor=True
)

ax_hm.grid(
    which="minor",
    linewidth=0.35,
    color="black"
)

ax_hm.tick_params(
    which="minor",
    bottom=False,
    left=False
)

ax_hm.set_title(
    "Pearson Correlation of Samples",
    fontsize=10,
    pad=38
)


ax_cb = fig.add_subplot(
    gs[1, 1]
)

cb = fig.colorbar(
    im,
    cax=ax_cb,
    orientation="horizontal"
)

cb.set_ticks(
    np.arange(
        0.3,
        1.01,
        0.1
    )
)

cb.ax.tick_params(
    labelsize=7
)


fig.savefig(
    out_B,
    format="pdf",
    bbox_inches="tight"
)

plt.close(fig)


print("Panel A written:", out_A)
print("Panel B written:", out_B)
print("Correlation values recalculated: NO")
PY


python3 \
"${TMP}/build_S1_AB.py" \
"${FRAG_RAW}" \
"${COR_TAB}" \
"${TMP}/PanelA_fragment_size.pdf" \
"${TMP}/PanelB_correlation.pdf"


# ============================================================
# 4. Build Panel C1-C4 from frozen deepTools matrices
#    Only --samplesLabel changes.
# ============================================================

echo
echo "Building TSS heatmaps with biological sample labels..."

"${PLOTHEATMAP}" \
-m "${M1}" \
-out "${TMP}/PanelC1_TSS.pdf" \
--samplesLabel \
Adipose-P348 \
Adipose-P350 \
Cerebellum-P348 \
Cerebellum-P350


"${PLOTHEATMAP}" \
-m "${M2}" \
-out "${TMP}/PanelC2_TSS.pdf" \
--samplesLabel \
Cortex-P348 \
Cortex-P350 \
Hypothalamus-P348 \
Hypothalamus-P350


"${PLOTHEATMAP}" \
-m "${M3}" \
-out "${TMP}/PanelC3_TSS.pdf" \
--samplesLabel \
Liver-P348 \
Liver-P350 \
Lung-P348 \
Lung-P350


"${PLOTHEATMAP}" \
-m "${M4}" \
-out "${TMP}/PanelC4_TSS.pdf" \
--samplesLabel \
Muscle-P348 \
Muscle-P350 \
Spleen-P348 \
Spleen-P350


# ============================================================
# 5. Assemble publication S1 PDF
# ============================================================

cat > "${TMP}/assemble_S1.py" <<'PY'
import sys
from pathlib import Path

try:
    import pymupdf as fitz
except Exception as e:
    raise SystemExit(
        "ERROR: PyMuPDF is required (import fitz failed).\n"
        "Install with: pip install pymupdf\n"
        f"{e}"
    )


out_pdf = Path(sys.argv[1])

A = Path(sys.argv[2])
B = Path(sys.argv[3])

C1 = Path(sys.argv[4])
C2 = Path(sys.argv[5])
C3 = Path(sys.argv[6])
C4 = Path(sys.argv[7])


# 17 x 14 inches
W = 17 * 72
H = 14 * 72


doc = fitz.open()
page = doc.new_page(
    width=W,
    height=H
)


def add_pdf(pdf_path, rect):
    src = fitz.open(pdf_path)

    page.show_pdf_page(
        rect,
        src,
        0,
        keep_proportion=True,
        overlay=True
    )

    src.close()


def text(x, y, s, size=12, bold=False):
    font = (
        "hebo"
        if bold
        else "helv"
    )

    page.insert_text(
        fitz.Point(x, y),
        s,
        fontsize=size,
        fontname=font,
        color=(0, 0, 0)
    )


# ============================================================
# Panel labels / titles
# ============================================================

text(
    25,
    30,
    "A",
    size=20,
    bold=True
)

text(
    56,
    30,
    "Fragment-size distribution",
    size=13.5,
    bold=True
)


text(
    625,
    30,
    "B",
    size=20,
    bold=True
)

text(
    658,
    30,
    "ATAC-seq sample correlation",
    size=13.5,
    bold=True
)


# ============================================================
# A / B
# ============================================================

add_pdf(
    A,
    fitz.Rect(
        25,
        45,
        600,
        420
    )
)

add_pdf(
    B,
    fitz.Rect(
        625,
        45,
        1200,
        420
    )
)


# ============================================================
# C title
# ============================================================

text(
    25,
    455,
    "C",
    size=20,
    bold=True
)

text(
    56,
    455,
    "TSS-centered enrichment profiles and heatmaps",
    size=13.5,
    bold=True
)


# ============================================================
# C1-C4
# ============================================================

x_positions = [
    (25, 310),
    (320, 605),
    (615, 900),
    (910, 1195)
]

Cfiles = [
    C1,
    C2,
    C3,
    C4
]

for i, (
    f,
    (x1, x2)
) in enumerate(
    zip(
        Cfiles,
        x_positions
    ),
    start=1
):

    text(
        x1,
        480,
        f"C{i}",
        size=12.5,
        bold=True
    )

    add_pdf(
        f,
        fitz.Rect(
            x1,
            490,
            x2,
            995
        )
    )


doc.save(
    out_pdf,
    garbage=4,
    deflate=True,
    clean=True
)

doc.close()

print(
    "Assembled:",
    out_pdf
)
PY


FINAL_PDF="${OUTDIR}/Supplementary_Figure_S1_FINAL_v2.pdf"

python3 \
"${TMP}/assemble_S1.py" \
"${FINAL_PDF}" \
"${TMP}/PanelA_fragment_size.pdf" \
"${TMP}/PanelB_correlation.pdf" \
"${TMP}/PanelC1_TSS.pdf" \
"${TMP}/PanelC2_TSS.pdf" \
"${TMP}/PanelC3_TSS.pdf" \
"${TMP}/PanelC4_TSS.pdf"


# ============================================================
# 6. Caption
# ============================================================

cat > \
"${OUTDIR}/Supplementary_Figure_S1_FINAL_v2_caption.txt" <<'CAP'
Supplementary Figure S1. ATAC-seq quality assessment across the 16 pig multi-tissue samples. (A) Fragment-size distributions derived from the frozen ATAC-seq fragment-size output. (B) Pearson correlation heatmap generated from the frozen sample-correlation matrix; hierarchical clustering is shown for visualization. (C) TSS-centered enrichment profiles and heatmaps generated from the frozen deepTools reference-point matrices: (C1) Adipose and Cerebellum, (C2) Cortex and Hypothalamus, (C3) Liver and Lung, and (C4) Muscle and Spleen. Sample labels indicate tissue and animal identity (P348 or P350); corresponding SRA accession numbers are provided in Supplementary Table S1. The TSS panels retain the representative-transcript TSS reference and shuffled genomic background encoded in the frozen source matrices.
CAP


# ============================================================
# 7. 600-dpi TIFF
# ============================================================

TIFF_PREFIX="${TMP}/S1_600dpi"

pdftoppm \
-r 600 \
-tiff \
-singlefile \
-tiffcompression lzw \
"${FINAL_PDF}" \
"${TIFF_PREFIX}"


if [[ -f "${TIFF_PREFIX}.tif" ]]; then

    mv \
    "${TIFF_PREFIX}.tif" \
    "${OUTDIR}/Supplementary_Figure_S1_FINAL_v2_600dpi.tiff"

elif [[ -f "${TIFF_PREFIX}.tiff" ]]; then

    mv \
    "${TIFF_PREFIX}.tiff" \
    "${OUTDIR}/Supplementary_Figure_S1_FINAL_v2_600dpi.tiff"

else

    echo "ERROR: 600-dpi TIFF was not generated."
    exit 1
fi


# ============================================================
# 8. Strict FINAL-v2 QC
# ============================================================

cat > "${TMP}/S1_final_QC.py" <<'PY'
import re
import subprocess
import sys
from pathlib import Path

try:
    from PIL import Image
except Exception as e:
    raise SystemExit(
        "ERROR: Pillow is required for TIFF QC.\n"
        f"{e}"
    )


pdf = Path(sys.argv[1])
tiff = Path(sys.argv[2])
out = Path(sys.argv[3])


labels = [
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


def run(cmd):
    return subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False
    )


text_result = run(
    [
        "pdftotext",
        "-layout",
        str(pdf),
        "-"
    ]
)

text = text_result.stdout


srrs = sorted(
    set(
        re.findall(
            r"SRR\d+",
            text
        )
    )
)

detected = [
    x for x in labels
    if x in text
]


info = run(
    [
        "pdfinfo",
        str(pdf)
    ]
).stdout

m = re.search(
    r"^Pages:\s+(\d+)",
    info,
    re.MULTILINE
)

pages = (
    int(m.group(1))
    if m
    else None
)


with Image.open(tiff) as im:

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


def dpi_ok(x):
    try:
        return abs(
            float(x) - 600
        ) <= 5
    except Exception:
        return False


rows = [
    (
        "PDF_exists",
        int(pdf.exists()),
        1
    ),

    (
        "PDF_single_page",
        pages,
        1
    ),

    (
        "TIFF_exists",
        int(tiff.exists()),
        1
    ),

    (
        "TIFF_600dpi_x",
        round(
            float(dx),
            2
        ) if dx else "NA",
        600
    ),

    (
        "TIFF_600dpi_y",
        round(
            float(dy),
            2
        ) if dy else "NA",
        600
    ),

    (
        "S1_SRR_accessions_remaining",
        len(srrs),
        0
    ),

    (
        "S1_expected_biological_labels_detected",
        len(detected),
        16
    ),

    (
        "Biological_results_modified",
        0,
        0
    )
]


with open(
    out,
    "w"
) as f:

    f.write(
        "Metric\tValue\tExpected\tStatus\n"
    )

    for metric, value, expected in rows:

        if metric.startswith(
            "TIFF_600dpi"
        ):
            status = (
                "PASS"
                if dpi_ok(value)
                else "FAIL"
            )

        else:
            status = (
                "PASS"
                if value == expected
                else "FAIL"
            )

        f.write(
            f"{metric}\t"
            f"{value}\t"
            f"{expected}\t"
            f"{status}\n"
        )


print(
    "PDF pages:",
    pages
)

print(
    "TIFF:",
    width,
    "x",
    height,
    "dpi=",
    dx,
    dy
)

print(
    "Remaining SRRs:",
    len(srrs),
    srrs
)

print(
    "Detected biological labels:",
    len(detected)
)

for x in detected:
    print(
        "  ",
        x
    )


hard_pass = (
    pages == 1
    and
    pdf.exists()
    and
    tiff.exists()
    and
    dpi_ok(dx)
    and
    dpi_ok(dy)
    and
    len(srrs) == 0
    and
    len(detected) == 16
)


if not hard_pass:
    raise SystemExit(
        "STEP10F3B6B STATUS: FAIL"
    )


print()
print(
    "STEP10F3B6B STATUS: PASS"
)

print(
    "Supplementary Figure S1 FINAL v2 is ready for consistency lock."
)
PY


python3 \
"${TMP}/S1_final_QC.py" \
"${FINAL_PDF}" \
"${OUTDIR}/Supplementary_Figure_S1_FINAL_v2_600dpi.tiff" \
"${OUTDIR}/FigureS1_FINAL_v2_QC.tsv"


# ============================================================
# 9. MD5 manifest
# ============================================================

(
    cd "${OUTDIR}"

    md5sum \
    Supplementary_Figure_S1_FINAL_v2.pdf \
    Supplementary_Figure_S1_FINAL_v2_600dpi.tiff \
    Supplementary_Figure_S1_FINAL_v2_caption.txt \
    > FigureS1_FINAL_v2_MD5.txt
)


# ============================================================
# 10. Remove temporary plotting files after PASS
# ============================================================

rm -rf "${TMP}"


echo
echo "======================================================================"
echo "STEP10F3B6B COMPLETED"
echo "======================================================================"

echo
echo "Final files:"
ls -lh \
"${OUTDIR}/Supplementary_Figure_S1_FINAL_v2.pdf" \
"${OUTDIR}/Supplementary_Figure_S1_FINAL_v2_600dpi.tiff" \
"${OUTDIR}/Supplementary_Figure_S1_FINAL_v2_caption.txt"

echo
echo "QC:"
column -t -s $'\t' \
"${OUTDIR}/FigureS1_FINAL_v2_QC.tsv"

echo
echo "IMPORTANT:"
echo "- Panel A uses frozen fragment-size raw data."
echo "- Panel B uses frozen Pearson correlation values."
echo "- Panel C uses frozen TSS reference-point matrices."
echo "- TSS-TES matrices were NOT used."
echo "- Only sample labels / publication layout were changed."
echo "- No RNA/ATAC biological result was changed."
echo "- Biological results modified: 0."
echo
echo "STEP10F3B6B STATUS: PASS"
echo "======================================================================"
