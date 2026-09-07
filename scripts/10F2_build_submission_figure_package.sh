#!/usr/bin/env bash
set -euo pipefail

ROOT="10_publication_figures/Step10F_final_consistency"
F1DIR="${ROOT}/Step10F1"
OUT="${ROOT}/Step10F2_submission_package"

PDFDIR="${OUT}/01_Main_Figures_PDF"
TIFFDIR="${OUT}/02_Main_Figures_TIFF_600dpi"
DOCDIR="${OUT}/03_Figure_Legends_and_QC"

mkdir -p "${PDFDIR}" "${TIFFDIR}" "${DOCDIR}"

echo "============================================================"
echo "STEP10F2 — BUILD SUBMISSION-READY FIGURE PACKAGE"
echo "Copy-only from locked Step10F1 figures"
echo "No redraw; no conversion; no biological analysis"
echo "============================================================"

# ------------------------------------------------------------
# Authoritative locked sources
# ------------------------------------------------------------

PDF1="10_publication_figures/Step10B_Figure1_2_v2/Figure1/Figure1_FINAL_v2.pdf"
TIFF1="10_publication_figures/Step10B_Figure1_2_v2/Figure1/Figure1_FINAL_v2_600dpi.tiff"

PDF2="10_publication_figures/Step10preF_Figure2_FINAL/Figure2_FINAL.pdf"
TIFF2="10_publication_figures/Step10preF_Figure2_FINAL/Figure2_FINAL_600dpi.tiff"

PDF3="10_publication_figures/Step10C_Figure3_4_v3/Figure3/Figure3_FINAL_v3.pdf"
TIFF3="10_publication_figures/Step10C_Figure3_4_v3/Figure3/Figure3_FINAL_v3_600dpi.tiff"

PDF4="10_publication_figures/Step10C_Figure3_4_v3/Figure4/Figure4_FINAL_v3.pdf"
TIFF4="10_publication_figures/Step10C_Figure3_4_v3/Figure4/Figure4_FINAL_v3_600dpi.tiff"

PDF5="10_publication_figures/Step10preF_Figure5_FINAL/Figure5_FINAL.pdf"
TIFF5="10_publication_figures/Step10preF_Figure5_FINAL/Figure5_FINAL_600dpi.tiff"

PDF6="10_publication_figures/Step10E_Figure6_FINAL/Figure6_FINAL.pdf"
TIFF6="10_publication_figures/Step10E_Figure6_FINAL/Figure6_FINAL_600dpi.tiff"

# ------------------------------------------------------------
# Safety check
# ------------------------------------------------------------

for f in \
    "$PDF1" "$TIFF1" \
    "$PDF2" "$TIFF2" \
    "$PDF3" "$TIFF3" \
    "$PDF4" "$TIFF4" \
    "$PDF5" "$TIFF5" \
    "$PDF6" "$TIFF6"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty locked source:"
        echo "$f"
        exit 1
    fi
done

# ------------------------------------------------------------
# Copy only
# ------------------------------------------------------------

cp -p "$PDF1"  "${PDFDIR}/Figure1.pdf"
cp -p "$PDF2"  "${PDFDIR}/Figure2.pdf"
cp -p "$PDF3"  "${PDFDIR}/Figure3.pdf"
cp -p "$PDF4"  "${PDFDIR}/Figure4.pdf"
cp -p "$PDF5"  "${PDFDIR}/Figure5.pdf"
cp -p "$PDF6"  "${PDFDIR}/Figure6.pdf"

cp -p "$TIFF1" "${TIFFDIR}/Figure1_600dpi.tiff"
cp -p "$TIFF2" "${TIFFDIR}/Figure2_600dpi.tiff"
cp -p "$TIFF3" "${TIFFDIR}/Figure3_600dpi.tiff"
cp -p "$TIFF4" "${TIFFDIR}/Figure4_600dpi.tiff"
cp -p "$TIFF5" "${TIFFDIR}/Figure5_600dpi.tiff"
cp -p "$TIFF6" "${TIFFDIR}/Figure6_600dpi.tiff"

# ------------------------------------------------------------
# Copy frozen captions / manifests / QC
# ------------------------------------------------------------

cp -p \
"${F1DIR}/10F1_main_figure_captions.md" \
"${DOCDIR}/Main_Figure_Legends.md"

cp -p \
"${F1DIR}/10F1_main_figure_lock_manifest.tsv" \
"${DOCDIR}/Main_Figure_Lock_Manifest.tsv"

cp -p \
"${F1DIR}/10F1_panel_semantic_manifest.tsv" \
"${DOCDIR}/Main_Figure_Panel_Manifest.tsv"

cp -p \
"${F1DIR}/10F1_terminology_lock.tsv" \
"${DOCDIR}/Main_Figure_Terminology_Lock.tsv"

cp -p \
"${F1DIR}/10F1_technical_lock_QC.tsv" \
"${DOCDIR}/Main_Figure_Technical_QC.tsv"

cp -p \
"${F1DIR}/10F1_overall_status.tsv" \
"${DOCDIR}/Step10F1_Overall_Status.tsv"

# ------------------------------------------------------------
# Source -> submission checksum identity
# ------------------------------------------------------------

CHECKSUM="${OUT}/10F2_source_vs_submission_MD5.tsv"

printf \
"Figure\tFormat\tSource_MD5\tSubmission_MD5\tIdentical\n" \
> "$CHECKSUM"

check_pair () {

    fig="$1"
    format="$2"
    source="$3"
    target="$4"

    src_md5=$(md5sum "$source" | awk '{print $1}')
    dst_md5=$(md5sum "$target" | awk '{print $1}')

    if [[ "$src_md5" == "$dst_md5" ]]; then
        identical="PASS"
    else
        identical="FAIL"
    fi

    printf \
    "%s\t%s\t%s\t%s\t%s\n" \
    "$fig" "$format" "$src_md5" "$dst_md5" "$identical" \
    >> "$CHECKSUM"
}

check_pair Figure1 PDF  "$PDF1"  "${PDFDIR}/Figure1.pdf"
check_pair Figure1 TIFF "$TIFF1" "${TIFFDIR}/Figure1_600dpi.tiff"

check_pair Figure2 PDF  "$PDF2"  "${PDFDIR}/Figure2.pdf"
check_pair Figure2 TIFF "$TIFF2" "${TIFFDIR}/Figure2_600dpi.tiff"

check_pair Figure3 PDF  "$PDF3"  "${PDFDIR}/Figure3.pdf"
check_pair Figure3 TIFF "$TIFF3" "${TIFFDIR}/Figure3_600dpi.tiff"

check_pair Figure4 PDF  "$PDF4"  "${PDFDIR}/Figure4.pdf"
check_pair Figure4 TIFF "$TIFF4" "${TIFFDIR}/Figure4_600dpi.tiff"

check_pair Figure5 PDF  "$PDF5"  "${PDFDIR}/Figure5.pdf"
check_pair Figure5 TIFF "$TIFF5" "${TIFFDIR}/Figure5_600dpi.tiff"

check_pair Figure6 PDF  "$PDF6"  "${PDFDIR}/Figure6.pdf"
check_pair Figure6 TIFF "$TIFF6" "${TIFFDIR}/Figure6_600dpi.tiff"

# ------------------------------------------------------------
# Final technical QC of copied package
# ------------------------------------------------------------

QC="${OUT}/10F2_submission_package_QC.tsv"

printf \
"Figure\tPDF_exists\tPDF_pages\tTIFF_exists\tWidth_px\tHeight_px\tX_dpi\tY_dpi\tFonts_embedded\tMD5_identity\tStatus\n" \
> "$QC"

for fig in 1 2 3 4 5 6
do

    pdf="${PDFDIR}/Figure${fig}.pdf"
    tif="${TIFFDIR}/Figure${fig}_600dpi.tiff"

    pdf_exists="FALSE"
    tif_exists="FALSE"
    pages="NA"
    width="NA"
    height="NA"
    xdpi="NA"
    ydpi="NA"
    fonts="FALSE"

    [[ -s "$pdf" ]] && pdf_exists="TRUE"
    [[ -s "$tif" ]] && tif_exists="TRUE"

    if [[ "$pdf_exists" == "TRUE" ]]; then

        pages=$(
            pdfinfo "$pdf" |
            awk '/^Pages:/ {print $2}'
        )

        if pdffonts "$pdf" |
           awk 'NR>2 && NF>0 {
               if ($(NF-4)!="yes") bad=1
           }
           END {exit bad}'
        then
            fonts="TRUE"
        fi
    fi

    if [[ "$tif_exists" == "TRUE" ]]; then

        read -r width height xdpi ydpi <<< "$(
            identify -format '%w %h %x %y' "$tif" |
            sed 's/ PixelsPerInch//g'
        )"
    fi

    md5status=$(
        awk -F '\t' \
        -v f="Figure${fig}" \
        '$1==f {
            if ($5!="PASS") bad=1
        }
        END {
            if (bad) print "FAIL";
            else print "PASS"
        }' \
        "$CHECKSUM"
    )

    status="PASS"

    [[ "$pdf_exists" == "TRUE" ]] || status="FAIL"
    [[ "$pages" == "1" ]] || status="FAIL"
    [[ "$tif_exists" == "TRUE" ]] || status="FAIL"
    [[ "$fonts" == "TRUE" ]] || status="FAIL"
    [[ "$md5status" == "PASS" ]] || status="FAIL"

    awk -v x="$xdpi" -v y="$ydpi" \
    'BEGIN {
        if (!(x>=599 && x<=601 && y>=599 && y<=601))
            exit 1
    }' || status="FAIL"

    printf \
    "Figure%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$fig" \
    "$pdf_exists" \
    "$pages" \
    "$tif_exists" \
    "$width" \
    "$height" \
    "$xdpi" \
    "$ydpi" \
    "$fonts" \
    "$md5status" \
    "$status" \
    >> "$QC"
done

# ------------------------------------------------------------
# File manifest
# ------------------------------------------------------------

MANIFEST="${OUT}/10F2_submission_file_manifest.tsv"

printf \
"Category\tFilename\tSize_bytes\tMD5\n" \
> "$MANIFEST"

for f in \
    "${PDFDIR}"/* \
    "${TIFFDIR}"/* \
    "${DOCDIR}"/*
do

    printf \
    "%s\t%s\t%s\t%s\n" \
    "$(basename "$(dirname "$f")")" \
    "$(basename "$f")" \
    "$(stat -c '%s' "$f")" \
    "$(md5sum "$f" | awk '{print $1}')" \
    >> "$MANIFEST"
done

# ------------------------------------------------------------
# README
# ------------------------------------------------------------

cat > "${OUT}/README_submission_figures.txt" <<'TXT'
Main Figure Submission Package
==============================

Contents

01_Main_Figures_PDF
    Figure1.pdf
    Figure2.pdf
    Figure3.pdf
    Figure4.pdf
    Figure5.pdf
    Figure6.pdf

02_Main_Figures_TIFF_600dpi
    Figure1_600dpi.tiff
    Figure2_600dpi.tiff
    Figure3_600dpi.tiff
    Figure4_600dpi.tiff
    Figure5_600dpi.tiff
    Figure6_600dpi.tiff

03_Figure_Legends_and_QC
    Main_Figure_Legends.md
    Main_Figure_Lock_Manifest.tsv
    Main_Figure_Panel_Manifest.tsv
    Main_Figure_Terminology_Lock.tsv
    Main_Figure_Technical_QC.tsv
    Step10F1_Overall_Status.tsv

All files were copied directly from the locked Step10F1
authoritative figures.

No figure was redrawn.
No image was rescaled.
No PDF was converted.
No TIFF was converted.
No biological analysis was rerun.
No analysis threshold was changed.

MD5 identity between locked source figures and submission
copies is recorded in:

10F2_source_vs_submission_MD5.tsv
TXT

# ------------------------------------------------------------
# Overall status
# ------------------------------------------------------------

FAIL_N=$(
    awk -F '\t' \
    'NR>1 && $11!="PASS" {n++}
     END {print n+0}' \
    "$QC"
)

MD5_FAIL_N=$(
    awk -F '\t' \
    'NR>1 && $5!="PASS" {n++}
     END {print n+0}' \
    "$CHECKSUM"
)

OVERALL="${OUT}/10F2_overall_status.tsv"

printf \
"Metric\tValue\tExpected\tStatus\n" \
> "$OVERALL"

printf \
"Main_PDF_files\t6\t6\tPASS\n" \
>> "$OVERALL"

printf \
"Main_TIFF_files\t6\t6\tPASS\n" \
>> "$OVERALL"

printf \
"Source_submission_MD5_failures\t%s\t0\t%s\n" \
"$MD5_FAIL_N" \
"$([[ "$MD5_FAIL_N" -eq 0 ]] && echo PASS || echo FAIL)" \
>> "$OVERALL"

printf \
"Technical_QC_failures\t%s\t0\t%s\n" \
"$FAIL_N" \
"$([[ "$FAIL_N" -eq 0 ]] && echo PASS || echo FAIL)" \
>> "$OVERALL"

printf \
"Figures_redrawn_in_Step10F2\t0\t0\tPASS\n" \
>> "$OVERALL"

printf \
"Biological_results_modified\t0\t0\tPASS\n" \
>> "$OVERALL"

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

echo
echo "============================================================"
echo "STEP10F2 COMPLETED"
echo "============================================================"

echo
echo "Submission package QC:"
column -t -s $'\t' "$QC"

echo
echo "Source -> submission MD5 identity:"
column -t -s $'\t' "$CHECKSUM"

echo
echo "Overall:"
column -t -s $'\t' "$OVERALL"

echo
echo "Submission directory:"
echo "$OUT"

echo
echo "Main PDFs:"
ls -lh "${PDFDIR}"

echo
echo "Main TIFFs:"
ls -lh "${TIFFDIR}"

echo
echo "IMPORTANT:"
echo "- Figure1-Figure6 remain LOCKED."
echo "- Files were copied only."
echo "- No redraw or conversion occurred."
echo "- MD5 identity verifies unchanged figure content."

if [[ "$FAIL_N" -ne 0 || "$MD5_FAIL_N" -ne 0 ]]; then
    echo
    echo "STEP10F2 STATUS: FAIL"
    exit 1
fi

echo
echo "STEP10F2 STATUS: PASS"
echo "Main figures are submission-package ready."
echo "============================================================"
