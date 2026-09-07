#!/usr/bin/env bash
set -euo pipefail

MOTIF_DB="/datadisk2/liujingjian_data/ATAC-seq/Sus/footprints/footprinting/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_v2.meme"

INPUT_ROOT="08B2b_STREME_noholdout"
OUTROOT="08B3b_Tomtom_noholdout_JASPAR2024"

mkdir -p "$OUTROOT"

printf \
"Tissue\tQuery_motifs\tTomtom_q005_rows\tAnnotated_query_motifs\n" \
> "${OUTROOT}/08B3b_summary.tsv"

for tissue in Muscle Spleen Liver Cerebellum
do

    echo
    echo "============================================================"
    echo "Tissue: $tissue"
    echo "============================================================"

    QUERY="${INPUT_ROOT}/${tissue}/streme.txt"
    OUT="${OUTROOT}/${tissue}"

    rm -rf "$OUT"

    tomtom \
        -oc "$OUT" \
        -dist pearson \
        -thresh 0.05 \
        "$QUERY" \
        "$MOTIF_DB"

    NQUERY=$(grep -c '^MOTIF ' "$QUERY" || true)

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

    NANNOT=$(
        awk -F $'\t' '
            !/^#/ &&
            $1 != "Query_ID" &&
            NF > 1 {
                seen[$1]=1
            }
            END {
                for (x in seen) n++
                print n+0
            }
        ' "${OUT}/tomtom.tsv"
    )

    printf \
    "%s\t%s\t%s\t%s\n" \
    "$tissue" \
    "$NQUERY" \
    "$NMATCH" \
    "$NANNOT" \
    >> "${OUTROOT}/08B3b_summary.tsv"

done

echo
echo "============================================================"
echo "STEP 08B3b COMPLETED"
echo "============================================================"

column -t -s $'\t' \
"${OUTROOT}/08B3b_summary.tsv"

