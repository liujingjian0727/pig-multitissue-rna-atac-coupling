#!/usr/bin/env bash
set -euo pipefail

OUTDIR="10_publication_figures/Step10F_final_consistency/Step10F3A_ATAC_QC"
mkdir -p "${OUTDIR}"

REPORT="${OUTDIR}/10F3A0_ATAC_QC_source_inventory.txt"
TABLE="${OUTDIR}/10F3A0_ATAC_QC_candidate_files.tsv"

echo "============================================================" | tee "${REPORT}"
echo "STEP10F3A0 — ATAC QC SOURCE INVENTORY" | tee -a "${REPORT}"
echo "Inspection only; no analysis and no redraw" | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"
echo | tee -a "${REPORT}"

# ------------------------------------------------------------
# Search roots
# ------------------------------------------------------------

SEARCH_ROOTS=(
    "."
    "/datadisk2/liujingjian_data/ATAC-seq/Sus"
)

printf "Category\tFile\tExtension\tSize_bytes\n" > "${TABLE}"

collect_files () {

    local category="$1"
    shift

    for root in "${SEARCH_ROOTS[@]}"
    do
        [[ -d "$root" ]] || continue

        find "$root" \
            -type f \
            \( "$@" \) \
            2>/dev/null |
        sort -u |
        while IFS= read -r f
        do
            [[ -f "$f" ]] || continue

            ext="${f##*.}"
            size=$(stat -c '%s' "$f" 2>/dev/null || echo "NA")

            printf "%s\t%s\t%s\t%s\n" \
                "$category" "$f" "$ext" "$size" \
                >> "${TABLE}"
        done
    done
}

# ------------------------------------------------------------
# A. Fragment-size / insert-size QC
# ------------------------------------------------------------

collect_files "Fragment_size" \
    -iname '*fragment*size*.pdf' -o \
    -iname '*fragment*size*.png' -o \
    -iname '*fragment*size*.tif' -o \
    -iname '*fragment*size*.tiff' -o \
    -iname '*insert*size*.pdf' -o \
    -iname '*insert*size*.png' -o \
    -iname '*insert*size*.tif' -o \
    -iname '*insert*size*.tiff'

# ------------------------------------------------------------
# B. ATAC sample-correlation heatmap
# ------------------------------------------------------------

collect_files "ATAC_correlation" \
    -iname '*ATAC*correlation*.pdf' -o \
    -iname '*ATAC*correlation*.png' -o \
    -iname '*ATAC*correlation*.tif' -o \
    -iname '*ATAC*correlation*.tiff' -o \
    -iname '*ATAC*corr*.pdf' -o \
    -iname '*ATAC*corr*.png' -o \
    -iname '*correlation*heatmap*.pdf' -o \
    -iname '*correlation*heatmap*.png' -o \
    -iname '*sample*correlation*.pdf' -o \
    -iname '*sample*correlation*.png'

# ------------------------------------------------------------
# C. TSS enrichment profile
# ------------------------------------------------------------

collect_files "TSS_profile" \
    -iname '*TSS*profile*.pdf' -o \
    -iname '*TSS*profile*.png' -o \
    -iname '*TSS*profile*.tif' -o \
    -iname '*TSS*profile*.tiff' -o \
    -iname '*TSS*enrichment*.pdf' -o \
    -iname '*TSS*enrichment*.png' -o \
    -iname '*TSS*enrichment*.tif' -o \
    -iname '*TSS*enrichment*.tiff' -o \
    -iname '*TSS*only*.pdf' -o \
    -iname '*TSS*only*.png'

# ------------------------------------------------------------
# D. TSS-centered heatmap / shuffled control
# ------------------------------------------------------------

collect_files "TSS_heatmap" \
    -iname '*TSS*heatmap*.pdf' -o \
    -iname '*TSS*heatmap*.png' -o \
    -iname '*TSS*heatmap*.tif' -o \
    -iname '*TSS*heatmap*.tiff' -o \
    -iname '*TSS*center*.pdf' -o \
    -iname '*TSS*center*.png' -o \
    -iname '*shuffle*.pdf' -o \
    -iname '*shuffle*.png' -o \
    -iname '*random*TSS*.pdf' -o \
    -iname '*random*TSS*.png'

# ------------------------------------------------------------
# Related QC tables that may be useful for caption/QC
# ------------------------------------------------------------

collect_files "ATAC_QC_table" \
    -iname '*fragment*summary*.tsv' -o \
    -iname '*fragment*summary*.txt' -o \
    -iname '*TSS*summary*.tsv' -o \
    -iname '*TSS*summary*.txt' -o \
    -iname '*FRiP*.tsv' -o \
    -iname '*FRiP*.txt' -o \
    -iname '*filter*metric*.tsv' -o \
    -iname '*filter*metric*.txt' -o \
    -iname '*ATAC*QC*.tsv' -o \
    -iname '*ATAC*QC*.txt'

# ------------------------------------------------------------
# Remove exact duplicate rows
# ------------------------------------------------------------

{
    head -n 1 "${TABLE}"
    tail -n +2 "${TABLE}" | sort -u
} > "${TABLE}.tmp"

mv "${TABLE}.tmp" "${TABLE}"

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

echo "Candidate counts by category:" | tee -a "${REPORT}"
awk -F'\t' '
NR > 1 {n[$1]++}
END {
    for (k in n)
        printf "%-20s %d\n", k, n[k]
}' "${TABLE}" | sort | tee -a "${REPORT}"

echo | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"
echo "CANDIDATE FILES" | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"

column -t -s $'\t' "${TABLE}" | tee -a "${REPORT}"

# ------------------------------------------------------------
# PDF metadata
# ------------------------------------------------------------

PDFINFO="${OUTDIR}/10F3A0_candidate_PDF_metadata.tsv"

printf "Category\tFile\tPages\tWidth_pt\tHeight_pt\n" > "${PDFINFO}"

tail -n +2 "${TABLE}" |
while IFS=$'\t' read -r category file ext size
do
    case "${ext,,}" in
        pdf)
            pages=$(pdfinfo "$file" 2>/dev/null |
                    awk -F':' '/^Pages/{gsub(/^[ \t]+/,"",$2); print $2}')

            dims=$(pdfinfo "$file" 2>/dev/null |
                   awk '
                   /^Page size:/ {
                       for(i=1;i<=NF;i++){
                           if($i=="x"){
                               print $(i-1),$(i+1)
                               exit
                           }
                       }
                   }')

            width=$(echo "$dims" | awk '{print $1}')
            height=$(echo "$dims" | awk '{print $2}')

            printf "%s\t%s\t%s\t%s\t%s\n" \
                "$category" "$file" \
                "${pages:-NA}" \
                "${width:-NA}" \
                "${height:-NA}" \
                >> "${PDFINFO}"
            ;;
    esac
done

# ------------------------------------------------------------
# Raster metadata
# ------------------------------------------------------------

RASTER="${OUTDIR}/10F3A0_candidate_raster_metadata.tsv"

printf "Category\tFile\tWidth_px\tHeight_px\tX_resolution\tY_resolution\n" \
    > "${RASTER}"

tail -n +2 "${TABLE}" |
while IFS=$'\t' read -r category file ext size
do
    case "${ext,,}" in
        png|tif|tiff|jpg|jpeg)

            values=$(identify -format '%w %h %x %y' "$file" \
                     2>/dev/null || true)

            if [[ -n "$values" ]]
            then
                printf "%s\t%s\t%s\n" \
                    "$category" "$file" "$values" \
                    >> "${RASTER}"
            fi
            ;;
    esac
done

echo | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"
echo "STEP10F3A0 COMPLETED" | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"
echo | tee -a "${REPORT}"

echo "Source inventory:" | tee -a "${REPORT}"
echo "${TABLE}" | tee -a "${REPORT}"

echo "PDF metadata:" | tee -a "${REPORT}"
echo "${PDFINFO}" | tee -a "${REPORT}"

echo "Raster metadata:" | tee -a "${REPORT}"
echo "${RASTER}" | tee -a "${REPORT}"

echo | tee -a "${REPORT}"
echo "No source figure was modified." | tee -a "${REPORT}"
echo "No ATAC QC analysis was rerun." | tee -a "${REPORT}"
echo "No biological result was modified." | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"
