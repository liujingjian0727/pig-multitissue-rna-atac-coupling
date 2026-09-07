#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STEP 08E1c
# Final FIMO scan using unique peak IDs
#
# Reported sites:
#   FIMO p <= 1e-4
#
# High-confidence sites:
#   FIMO q <= 0.05
# ============================================================

MOTIF_ROOT="08D_motif_family_TF_master/FIMO_motif_sets"
SEQ_ROOT="08E1b_unique_peak_FASTA"
OUT_ROOT="08E1c_FIMO_unique_peak_ids"

mkdir -p "$OUT_ROOT"

SUMMARY="${OUT_ROOT}/08E1c_FIMO_run_summary.tsv"
PER_MOTIF="${OUT_ROOT}/08E1c_FIMO_per_motif_summary.tsv"

printf \
"Tissue\tTier\tMotif_source\tMotif_N\tForeground_N\tSites_p_le_1e-4\tSites_q_le_0.05\tPeaks_p_le_1e-4\tPeaks_q_le_0.05\tPeak_fraction_p_le_1e-4\tPeak_fraction_q_le_0.05\tStatus\n" \
> "$SUMMARY"

printf \
"Tissue\tTier\tMotif_source\tMotif_ID\tMotif_alt_ID\tSites_p_le_1e-4\tSites_q_le_0.05\tPeaks_p_le_1e-4\tPeaks_q_le_0.05\tForeground_N\tPeak_fraction_p_le_1e-4\tPeak_fraction_q_le_0.05\n" \
> "$PER_MOTIF"

for tissue in Muscle Spleen Liver Cerebellum
do

    FASTA="${SEQ_ROOT}/${tissue}/${tissue}.strong.foreground.unique.fa"

    if [[ ! -s "$FASTA" ]]
    then
        echo "ERROR: missing FASTA: $FASTA"
        exit 1
    fi

    FG_N=$(grep -c '^>' "$FASTA")

    echo
    echo "============================================================"
    echo "Tissue: $tissue"
    echo "Foreground sequences: $FG_N"
    echo "============================================================"

    for tier in primary secondary
    do
        for source in known denovo
        do

            MOTIF_FILE="${MOTIF_ROOT}/${tissue}/${tissue}.${source}.${tier}.meme"

            [[ -s "$MOTIF_FILE" ]] || continue

            MOTIF_N=$(grep -c '^MOTIF ' "$MOTIF_FILE")

            OUT="${OUT_ROOT}/${tissue}/${tier}/${source}"

            # Critical fix:
            # FIMO requires the complete parent path to exist.
            mkdir -p "$OUT"

            echo
            echo "------------------------------------------------------------"
            echo "$tissue | $tier | $source | motifs=$MOTIF_N"
            echo "Motif file: $MOTIF_FILE"
            echo "Sequence file: $FASTA"
            echo "Output: $OUT"
            echo "------------------------------------------------------------"

            fimo \
                --oc "$OUT" \
                --thresh 1e-4 \
                "$MOTIF_FILE" \
                "$FASTA"

            FIMO="${OUT}/fimo.tsv"

            if [[ ! -s "$FIMO" ]]
            then
                printf \
                "%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\tNO_OUTPUT\n" \
                "$tissue" "$tier" "$source" "$MOTIF_N" "$FG_N" \
                >> "$SUMMARY"

                continue
            fi

            # ------------------------------------------------
            # Run-level statistics
            # ------------------------------------------------

            STATS=$(
                awk -F '\t' '
                BEGIN {
                    sites  = 0
                    qsites = 0
                    qcol   = 0
                    seqcol = 0
                }

                NR == 1 {
                    for (i=1; i<=NF; i++) {
                        if ($i=="q-value") qcol=i
                        if ($i=="sequence_name") seqcol=i
                    }
                    next
                }

                /^#/ { next }

                NF > 1 {
                    sites++

                    if (seqcol > 0 && $seqcol != "") {
                        seq_all[$seqcol] = 1
                    }

                    if (qcol > 0) {
                        qv=$qcol

                        if (qv != "" && qv != "NA" && (qv+0) <= 0.05) {
                            qsites++

                            if (seqcol > 0 && $seqcol != "") {
                                seq_q[$seqcol] = 1
                            }
                        }
                    }
                }

                END {
                    na=0
                    nq=0

                    for (x in seq_all) na++
                    for (x in seq_q) nq++

                    print sites "\t" qsites "\t" na "\t" nq
                }
                ' "$FIMO"
            )

            SITES=$(printf '%s\n' "$STATS" | cut -f1)
            QSITES=$(printf '%s\n' "$STATS" | cut -f2)
            PEAKS=$(printf '%s\n' "$STATS" | cut -f3)
            QPEAKS=$(printf '%s\n' "$STATS" | cut -f4)

            FRAC=$(
                awk -v x="$PEAKS" -v n="$FG_N" \
                'BEGIN {
                    if (n > 0) printf "%.6f", x/n
                    else print "NA"
                }'
            )

            QFRAC=$(
                awk -v x="$QPEAKS" -v n="$FG_N" \
                'BEGIN {
                    if (n > 0) printf "%.6f", x/n
                    else print "NA"
                }'
            )

            printf \
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n" \
            "$tissue" \
            "$tier" \
            "$source" \
            "$MOTIF_N" \
            "$FG_N" \
            "$SITES" \
            "$QSITES" \
            "$PEAKS" \
            "$QPEAKS" \
            "$FRAC" \
            "$QFRAC" \
            >> "$SUMMARY"

            # ------------------------------------------------
            # Per-motif statistics
            #
            # Avoid multi-line if() syntax because some awk
            # versions are sensitive to line breaks.
            # ------------------------------------------------

            awk -F '\t' \
                -v OFS='\t' \
                -v tissue="$tissue" \
                -v tier="$tier" \
                -v source="$source" \
                -v fg="$FG_N" '
            BEGIN {
                mid=0
                maid=0
                sid=0
                qid=0
            }

            NR == 1 {
                for (i=1; i<=NF; i++) {
                    if ($i=="motif_id") mid=i
                    if ($i=="motif_alt_id") maid=i
                    if ($i=="sequence_name") sid=i
                    if ($i=="q-value") qid=i
                }
                next
            }

            /^#/ { next }

            NF > 1 {
                m=$mid
                alt=$maid
                key=m SUBSEP alt

                motif[key]=m
                altname[key]=alt

                sites[key]++

                peakkey=key SUBSEP $sid
                peak_all[peakkey]=1

                qv=$qid

                if (qv != "" && qv != "NA" && (qv+0) <= 0.05) {
                    qsites[key]++
                    peak_q[peakkey]=1
                }
            }

            END {
                for (k in motif) {

                    np=0
                    nqp=0

                    for (x in peak_all) {
                        split(x,a,SUBSEP)
                        testkey=a[1] SUBSEP a[2]

                        if (testkey == k) {
                            np++
                        }
                    }

                    for (x in peak_q) {
                        split(x,a,SUBSEP)
                        testkey=a[1] SUBSEP a[2]

                        if (testkey == k) {
                            nqp++
                        }
                    }

                    if (fg > 0) {
                        frac=np/fg
                        qfrac=nqp/fg
                    } else {
                        frac=0
                        qfrac=0
                    }

                    print \
                        tissue,
                        tier,
                        source,
                        motif[k],
                        altname[k],
                        sites[k]+0,
                        qsites[k]+0,
                        np,
                        nqp,
                        fg,
                        sprintf("%.6f",frac),
                        sprintf("%.6f",qfrac)
                }
            }
            ' "$FIMO" \
            >> "$PER_MOTIF"

        done
    done
done

echo
echo "============================================================"
echo "STEP 08E1c COMPLETED"
echo "============================================================"

echo
echo "Run-level summary:"
column -t -s $'\t' "$SUMMARY"

echo
echo "Per-motif summary:"
column -t -s $'\t' "$PER_MOTIF"

echo
echo "Main outputs:"
echo "$SUMMARY"
echo "$PER_MOTIF"

echo "============================================================"
