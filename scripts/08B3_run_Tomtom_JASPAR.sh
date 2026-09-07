#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STEP 08B3
# TOMTOM ANNOTATION OF STREME DE NOVO MOTIFS
#
# Query:
#   STREME motifs
#
# Target:
#   JASPAR2024 CORE Vertebrates non-redundant
#
# Primary match significance:
#   Tomtom q-value < 0.05
# ============================================================


MOTIF_DB="/datadisk2/liujingjian_data/ATAC-seq/Sus/footprints/footprinting/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_v2.meme"

INPUT_ROOT="08B2_STREME_de_novo"
OUTROOT="08B3_Tomtom_JASPAR2024"

TISSUES=(
    Muscle
    Spleen
    Liver
    Cerebellum
)


echo "============================================================"
echo "STEP 08B3: TOMTOM JASPAR2024 ANNOTATION"
echo "============================================================"
echo
echo "Tomtom version:"
tomtom --version
echo


if [[ ! -s "$MOTIF_DB" ]]
then
    echo "ERROR: JASPAR motif database missing."
    exit 1
fi


mkdir -p "$OUTROOT"


SUMMARY="${OUTROOT}/08B3_Tomtom_run_summary.tsv"

printf \
"Tissue\tQuery_motif_N\tStatus\tTomtom_match_rows\n" \
> "$SUMMARY"


for tissue in "${TISSUES[@]}"
do

    echo
    echo "============================================================"
    echo "Tissue: $tissue"
    echo "============================================================"


    QUERY="${INPUT_ROOT}/${tissue}/streme.txt"

    OUT="${OUTROOT}/${tissue}"


    if [[ ! -s "$QUERY" ]]
    then
        echo "WARNING: no STREME query file:"
        echo "$QUERY"

        printf \
        "%s\t0\tSKIPPED\t0\n" \
        "$tissue" \
        >> "$SUMMARY"

        continue
    fi


    QUERY_N=$(grep -c '^MOTIF ' "$QUERY" || true)


    echo "Query motifs: $QUERY_N"


    if [[ "$QUERY_N" -eq 0 ]]
    then

        printf \
        "%s\t0\tNO_STREME_MOTIFS\t0\n" \
        "$tissue" \
        >> "$SUMMARY"

        continue
    fi


    rm -rf "$OUT"


    tomtom \
        -oc "$OUT" \
        -dist pearson \
        -thresh 0.05 \
        "$QUERY" \
        "$MOTIF_DB"


    if [[ ! -s "${OUT}/tomtom.tsv" ]]
    then

        STATUS="FAIL"
        MATCH_N=0

    else

        STATUS="PASS"

        MATCH_N=$(
            awk -F $'\t' '
                !/^#/ &&
                $1 != "Query_ID" &&
                NF > 1 {
                    n++
                }
                END {
                    print n+0
                }
            ' \
            "${OUT}/tomtom.tsv"
        )

    fi


    printf \
    "%s\t%s\t%s\t%s\n" \
    "$tissue" \
    "$QUERY_N" \
    "$STATUS" \
    "$MATCH_N" \
    >> "$SUMMARY"


    echo "Tomtom status     : $STATUS"
    echo "Significant rows  : $MATCH_N"

done


echo
echo "============================================================"
echo "STEP 08B3 COMPLETED"
echo "============================================================"

column -t -s $'\t' "$SUMMARY"

