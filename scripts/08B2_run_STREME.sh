#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STEP 08B2
# STREME DE NOVO MOTIF DISCOVERY
#
# MEME Suite 5.5.8
#
# Primary:
#   strong tissue-specific RNA-ATAC coupled peaks
#
# Control:
#   same-tissue HC non-strong matched peaks
#
# Matching already controls:
#   tissue
#   genomic class
#   GC
#   peak length
#
# Primary de novo significance:
#   STREME E-value < 0.05
# ============================================================


INPUT_ROOT="08A_matched_motif_background_v2"
OUTROOT="08B2_STREME_de_novo"

TISSUES=(
    Muscle
    Spleen
    Liver
    Cerebellum
)


echo "============================================================"
echo "STEP 08B2: STREME DE NOVO MOTIF DISCOVERY"
echo "============================================================"
echo
echo "STREME version:"
streme --version
echo


mkdir -p "$OUTROOT"


SUMMARY="${OUTROOT}/08B2_STREME_run_summary.tsv"

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


    if [[ ! -s "$FG" ]]
    then
        echo "ERROR: foreground missing:"
        echo "$FG"
        exit 1
    fi


    if [[ ! -s "$BG" ]]
    then
        echo "ERROR: background missing:"
        echo "$BG"
        exit 1
    fi


    NFG=$(grep -c '^>' "$FG")
    NBG=$(grep -c '^>' "$BG")


    echo "Foreground N : $NFG"
    echo "Background N : $NBG"


    if [[ "$NFG" -ne "$NBG" ]]
    then
        echo "ERROR: foreground/background sequence number mismatch."
        exit 1
    fi


    rm -rf "$OUT"


    streme \
        --oc "$OUT" \
        --p "$FG" \
        --n "$BG" \
        --dna \
        --objfun de \
        --minw 6 \
        --maxw 20 \
        --evalue \
        --thresh 0.05


    if [[ ! -s "${OUT}/streme.txt" ]]
    then

        STATUS="FAIL"
        MOTIF_N=0

    else

        STATUS="PASS"

        MOTIF_N=$(
            grep -c '^MOTIF ' \
            "${OUT}/streme.txt" \
            || true
        )

    fi


    printf \
    "%s\t%s\t%s\t%s\t%s\n" \
    "$tissue" \
    "$NFG" \
    "$NBG" \
    "$STATUS" \
    "$MOTIF_N" \
    >> "$SUMMARY"


    echo
    echo "$tissue STREME status : $STATUS"
    echo "$tissue motif count   : $MOTIF_N"

done


echo
echo "============================================================"
echo "STEP 08B2 COMPLETED"
echo "============================================================"

column -t -s $'\t' "$SUMMARY"

