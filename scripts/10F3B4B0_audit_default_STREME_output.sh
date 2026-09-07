#!/usr/bin/env bash
set -euo pipefail

OUTDIR="10_publication_figures/Step10F_final_consistency/Step10F3B_Supplementary_Figures/FigureS4_motif_analysis_QC"
mkdir -p "$OUTDIR"

OUT="${OUTDIR}/10F3B4B0_default_STREME_output_audit.txt"
COUNT="${OUTDIR}/10F3B4B0_default_STREME_motif_counts.tsv"

printf "Tissue\tMOTIF_lines\n" > "$COUNT"

{
echo "============================================================"
echo "STEP10F3B4B0 — AUDIT DEFAULT STREME OUTPUT"
echo "No analysis rerun"
echo "============================================================"
echo
} > "$OUT"


for tissue in Cerebellum Liver Muscle Spleen
do

    f="08B2_STREME_de_novo/${tissue}/streme.txt"

    echo "============================================================" >> "$OUT"
    echo "TISSUE: ${tissue}" >> "$OUT"
    echo "FILE  : ${f}" >> "$OUT"
    echo "============================================================" >> "$OUT"

    if [[ ! -f "$f" ]]
    then
        echo "FILE MISSING" >> "$OUT"
        printf "%s\tNA\n" "$tissue" >> "$COUNT"
        echo >> "$OUT"
        continue
    fi

    # Count genuine motif declarations
    n=$(
        grep -c '^MOTIF[[:space:]]' "$f" || true
    )

    printf "%s\t%s\n" \
        "$tissue" "$n" \
        >> "$COUNT"

    echo "MOTIF declaration count: $n" >> "$OUT"
    echo >> "$OUT"

    echo "---- MOTIF DECLARATIONS ----" >> "$OUT"

    grep -n \
        '^MOTIF[[:space:]]' \
        "$f" \
        >> "$OUT" || true

    echo >> "$OUT"

    echo "---- MOTIF + MATRIX / P / E CONTEXT ----" >> "$OUT"

    # Show lines that are likely to contain the per-motif
    # statistics without assuming a parser yet.
    grep -n -E \
        '^MOTIF[[:space:]]|letter-probability matrix|[[:space:]]P=|[[:space:]]E=|p-value|E-value|holdout|test|train' \
        "$f" \
        >> "$OUT" || true

    echo >> "$OUT"

    echo "---- FIRST 180 LINES ----" >> "$OUT"

    sed -n '1,180p' \
        "$f" \
        >> "$OUT"

    echo >> "$OUT"

done


{
echo "============================================================"
echo "MOTIF COUNT SUMMARY"
echo "============================================================"
column -t -s $'\t' "$COUNT"
echo
echo "TOTAL MOTIFS:"
awk -F'\t' '
NR > 1 && $2 ~ /^[0-9]+$/ {
    s += $2
}
END {
    print s+0
}' "$COUNT"
echo
echo "No STREME analysis was rerun."
echo "No result was modified."
echo "============================================================"
} >> "$OUT"


cat "$OUT"
