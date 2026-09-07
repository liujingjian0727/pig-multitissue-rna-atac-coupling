#!/usr/bin/env bash
set -euo pipefail

MOTIF_DB="/datadisk2/liujingjian_data/ATAC-seq/Sus/footprints/footprinting/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_v2.meme"

INPUT_ROOT="08A_matched_motif_background_v2"
OUTROOT="08B1b_AME_all_motifs"

mkdir -p "$OUTROOT"

TISSUES=(
    Muscle
    Spleen
    Liver
    Cerebellum
)

echo "============================================================"
echo "STEP 08B-1b: AME COMPLETE MOTIF OUTPUT"
echo "============================================================"
echo "AME version: $(ame --version)"
echo "Motifs: $(grep -c '^MOTIF ' "$MOTIF_DB")"
echo

for tissue in "${TISSUES[@]}"
do

    FG="${INPUT_ROOT}/${tissue}/${tissue}.strong.foreground.fa"
    BG="${INPUT_ROOT}/${tissue}/${tissue}.matched.background.fa"
    OUT="${OUTROOT}/${tissue}"

    NFG=$(grep -c '^>' "$FG")
    NBG=$(grep -c '^>' "$BG")

    echo "============================================================"
    echo "$tissue : ${NFG} vs ${NBG}"
    echo "============================================================"

    rm -rf "$OUT"

    ame \
        --oc "$OUT" \
        --control "$BG" \
        --method fisher \
        --scoring max \
        --evalue-report-threshold 1000 \
        "$FG" \
        "$MOTIF_DB"

    if [[ ! -s "${OUT}/ame.tsv" ]]
    then
        echo "ERROR: ${OUT}/ame.tsv missing"
        exit 1
    fi

done

echo
echo "============================================================"
echo "STEP 08B-1b COMPLETED"
echo "============================================================"
