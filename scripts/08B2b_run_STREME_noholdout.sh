#!/usr/bin/env bash
set -euo pipefail

INPUT_ROOT="08A_matched_motif_background_v2"
OUTROOT="08B2b_STREME_noholdout"

TISSUES=(
    Muscle
    Spleen
    Liver
    Cerebellum
)

mkdir -p "$OUTROOT"

SUMMARY="${OUTROOT}/08B2b_STREME_run_summary.tsv"

printf \
"Tissue\tForeground_N\tBackground_N\tStatus\tMotif_N\n" \
> "$SUMMARY"


for tissue in "${TISSUES[@]}"
do

    echo
    echo "============================================================"
    echo "Tissue: $tissue"
    echo "============================================================"

    FG="${INPUT_ROOT}/${tissue}/${tissue}.strong.foreground.fa"
    BG="${INPUT_ROOT}/${tissue}/${tissue}.matched.background.fa"
    OUT="${OUTROOT}/${tissue}"

    NFG=$(grep -c '^>' "$FG")
    NBG=$(grep -c '^>' "$BG")

    echo "Foreground = $NFG"
    echo "Background = $NBG"

    rm -rf "$OUT"

    streme \
        --oc "$OUT" \
        --p "$FG" \
        --n "$BG" \
        --dna \
        --objfun de \
        --minw 6 \
        --maxw 20 \
        --hofract 0 \
        --evalue \
        --thresh 0.05

    MOTIF_N=$(
        grep -c '^MOTIF ' \
        "${OUT}/streme.txt" \
        || true
    )

    printf \
    "%s\t%s\t%s\tPASS\t%s\n" \
    "$tissue" \
    "$NFG" \
    "$NBG" \
    "$MOTIF_N" \
    >> "$SUMMARY"

done


echo
echo "============================================================"
echo "STEP 08B2b COMPLETED"
echo "============================================================"

column -t -s $'\t' "$SUMMARY"

