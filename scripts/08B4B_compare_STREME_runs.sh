#!/usr/bin/env bash
set -euo pipefail

DEFAULT_ROOT="08B2_STREME_de_novo"
NOHOLD_ROOT="08B2b_STREME_noholdout"
OUTROOT="08B4B_STREME_reproducibility"

mkdir -p "$OUTROOT"

SUMMARY="${OUTROOT}/08B4B_reproducibility_summary.tsv"

printf \
"Tissue\tDefault_motifs\tNoHoldout_motifs\tTomtom_q005_rows\tDefault_motifs_reproduced\n" \
> "$SUMMARY"


for tissue in Muscle Spleen Liver Cerebellum
do

    echo
    echo "============================================================"
    echo "Tissue: $tissue"
    echo "============================================================"

    QUERY="${DEFAULT_ROOT}/${tissue}/streme.txt"
    TARGET="${NOHOLD_ROOT}/${tissue}/streme.txt"
    OUT="${OUTROOT}/${tissue}"

    NQ=$(grep -c '^MOTIF ' "$QUERY" || true)
    NT=$(grep -c '^MOTIF ' "$TARGET" || true)

    rm -rf "$OUT"

    tomtom \
        -oc "$OUT" \
        -dist pearson \
        -thresh 0.05 \
        "$QUERY" \
        "$TARGET"

    NMATCH=$(
        awk -F $'\t' '
            !/^#/ &&
            $1 != "Query_ID" &&
            NF > 1 {
                n++
            }
            END {
                print n+0
            }
        ' "${OUT}/tomtom.tsv"
    )

    NREPRO=$(
        awk -F $'\t' '
            !/^#/ &&
            $1 != "Query_ID" &&
            NF > 1 {
                seen[$1]=1
            }
            END {
                n=0
                for (x in seen) n++
                print n
            }
        ' "${OUT}/tomtom.tsv"
    )

    printf \
    "%s\t%s\t%s\t%s\t%s\n" \
    "$tissue" \
    "$NQ" \
    "$NT" \
    "$NMATCH" \
    "$NREPRO" \
    >> "$SUMMARY"

done


echo
echo "============================================================"
echo "STEP 08B4B COMPLETED"
echo "============================================================"

column -t -s $'\t' "$SUMMARY"

