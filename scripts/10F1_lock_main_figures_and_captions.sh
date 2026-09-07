#!/usr/bin/env bash
set -euo pipefail

OUTDIR="10_publication_figures/Step10F_final_consistency/Step10F1"
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "STEP10F1 — MAIN FIGURE VISUAL / SEMANTIC LOCK"
echo "Figure1–Figure6"
echo "No biological analysis and no figure redraw"
echo "============================================================"

# ============================================================
# 1. Authoritative final figure sources
# ============================================================

cat > "${OUTDIR}/10F1_main_figure_lock_manifest.tsv" <<'TSV'
FigurePanelsPDFTIFF600dpiLayout_statusBiological_statusStep10F1_status
Figure1A-E10_publication_figures/Step10B_Figure1_2_v2/Figure1/Figure1_FINAL_v2.pdf10_publication_figures/Step10B_Figure1_2_v2/Figure1/Figure1_FINAL_v2_600dpi.tiffLOCKEDFROZENPASS
Figure2A-E10_publication_figures/Step10preF_Figure2_FINAL/Figure2_FINAL.pdf10_publication_figures/Step10preF_Figure2_FINAL/Figure2_FINAL_600dpi.tiffLOCKEDFROZENPASS
Figure3A-F10_publication_figures/Step10C_Figure3_4_v3/Figure3/Figure3_FINAL_v3.pdf10_publication_figures/Step10C_Figure3_4_v3/Figure3/Figure3_FINAL_v3_600dpi.tiffLOCKEDFROZENPASS
Figure4A-E10_publication_figures/Step10C_Figure3_4_v3/Figure4/Figure4_FINAL_v3.pdf10_publication_figures/Step10C_Figure3_4_v3/Figure4/Figure4_FINAL_v3_600dpi.tiffLOCKEDFROZENPASS
Figure5A-E10_publication_figures/Step10preF_Figure5_FINAL/Figure5_FINAL.pdf10_publication_figures/Step10preF_Figure5_FINAL/Figure5_FINAL_600dpi.tiffLOCKEDFROZENPASS
Figure6A-C10_publication_figures/Step10E_Figure6_FINAL/Figure6_FINAL.pdf10_publication_figures/Step10E_Figure6_FINAL/Figure6_FINAL_600dpi.tiffLOCKEDFROZENPASS
TSV


# ============================================================
# 2. Panel semantic manifest
# ============================================================

cat > "${OUTDIR}/10F1_panel_semantic_manifest.tsv" <<'TSV'
FigurePanelPanel_content
Figure1APaired RNA-seq and ATAC-seq experimental design across eight tissues and two animals
Figure1BRNA-seq PCA
Figure1CATAC-seq PCA
Figure1DWithin-tissue RNA replicate correlation
Figure1EWithin-tissue ATAC replicate correlation
Figure2AGlobal tissue effects detected by DESeq2 LRT
Figure2BHigh-confidence tissue-specific RNA genes
Figure2CHigh-confidence tissue-specific ATAC peaks
Figure2DRNA tissue-specificity landscape using Tau and maximum SPM
Figure2EATAC tissue-specificity landscape using Tau and maximum SPM
Figure3APeak-gene candidate-selection flow
Figure3BDistribution of cross-tissue RNA-ATAC correlations
Figure3CAnimal-to-animal agreement of peak-gene correlations
Figure3DTissue distribution of strongly coupled peak-gene pairs
Figure3ESame-tissue permutation test
Figure3FCoupling-threshold sensitivity
Figure4APair-specific genomic-class composition
Figure4BPromoter enrichment odds ratios
Figure4CPeak-to-associated-gene TSS distance
Figure4DGenomic classes of strongly coupled peak-gene pairs
Figure4ENumber of strong ATAC peaks per associated gene
Figure5AMotif-analysis workflow using matched HC backgrounds
Figure5BMuscle promoter/TSS architecture of prioritized motif-positive peaks
Figure5COverlap between Muscle KLF9 and SP-like motif programs
Figure5DOverlap among Spleen GC-rich/KLF12 motif programs
Figure5EArchitecture-enrichment comparison across prioritized Muscle and Spleen motif programs
Figure6ARepresentative Muscle promoter-centered candidate network
Figure6BRepresentative Spleen overlapping GC-rich candidate network
Figure6CRepresentative Liver supporting GC-rich candidate network
TSV


# ============================================================
# 3. Final English captions
# ============================================================

cat > "${OUTDIR}/10F1_main_figure_captions.md" <<'CAPTIONS'
# Main figure legends

## Figure 1. Paired multi-tissue RNA-seq and ATAC-seq dataset and sample-level quality assessment.

(A) Experimental design comprising paired RNA-seq and ATAC-seq profiles from eight tissues in two pigs, P348 and P350.  
(B) Principal-component analysis of RNA-seq samples after count filtering and variance-stabilizing transformation.  
(C) Principal-component analysis of ATAC-seq samples using the corresponding filtered peak-count matrix.  
(D) Pearson correlations between P348 and P350 RNA-seq profiles within each tissue.  
(E) Pearson correlations between P348 and P350 ATAC-seq profiles within each tissue. Dashed horizontal lines indicate the mean within-tissue correlation across the eight tissues. Tissue colors are used consistently across panels.

## Figure 2. Global tissue effects and tissue specificity of RNA expression and chromatin accessibility.

(A) Numbers and percentages of RNA genes and ATAC peaks showing a significant global tissue effect in DESeq2 likelihood-ratio tests after controlling for animal identity (Benjamini-Hochberg FDR < 0.01).  
(B) Numbers of high-confidence tissue-specific RNA genes assigned to each tissue.  
(C) Numbers of high-confidence tissue-specific ATAC peaks assigned to each tissue. High-confidence tissue specificity required a significant global tissue effect, Tau >= 0.80, sufficient maximum mean signal, and concordant maximum-tissue assignment in P348 and P350.  
(D,E) Relationship between Tau tissue-specificity scores and maximum specificity-measure (SPM) values for RNA genes (D) and ATAC peaks (E). Only features with finite Tau and SPM values are displayed; all-zero features for which specificity statistics are mathematically undefined are omitted from these visualization panels only. The dashed line marks Tau = 0.80. Color intensity represents feature density.

## Figure 3. Quantitative coupling between chromatin accessibility and expression across tissues.

(A) Selection of UROPA-associated ATAC peak-gene pairs from all gene-assigned peaks to concordant high-confidence tissue-specific pairs and finally to strongly coupled candidate pairs.  
(B) Distributions of Spearman correlations between ATAC accessibility and expression across eight tissues for all assigned pairs, concordant high-confidence pairs, and strongly coupled pairs.  
(C) Agreement between P348- and P350-specific cross-tissue Spearman correlations among concordant high-confidence peak-gene pairs.  
(D) Tissue distribution of the 702 strongly coupled peak-gene candidate pairs.  
(E) Comparison of the observed number of strongly coupled pairs with the same-tissue permutation null distribution. Peak-to-gene assignments were permuted within tissue while preserving tissue composition and the RNA-profile pool.  
(F) Observed-to-permutation-null ratios across increasingly stringent coupling thresholds. Permutation significance is reported using empirical permutation probabilities. Correlation and permutation analyses describe coordinated accessibility-expression variation and do not establish causal regulation.

## Figure 4. Genomic architecture of coupled ATAC peak-gene associations.

(A) Pair-specific genomic-class composition of all UROPA-assigned peak-gene pairs, concordant high-confidence pairs, and strongly coupled pairs using a promoter-first annotation hierarchy. Promoters were defined relative to the representative-transcript TSS.  
(B) Promoter enrichment of strongly coupled pairs relative to disjoint non-Strong comparison sets from all assigned pairs and concordant high-confidence pairs. Points show odds ratios from Fisher's exact tests.  
(C) Absolute distance between ATAC peaks and the TSS of their associated genes across the three candidate sets.  
(D) Pair-specific genomic-class distribution of the 702 strongly coupled peak-gene pairs.  
(E) Number of strongly coupled ATAC peaks associated with each of the 450 unique genes represented among the 702 strong pairs. These genomic associations identify candidate regulatory architecture but do not by themselves establish direct regulatory relationships.

## Figure 5. Motif-associated regulatory architecture of strongly coupled peak-gene candidates.

(A) Motif-analysis framework. Strong tissue-specific peak-gene pairs were compared with one-to-one matched concordant high-confidence non-Strong backgrounds from the same tissue, matched by genomic class, GC content, and peak length. Known-motif prioritization, de novo motif discovery/annotation, and FIMO sequence localization were integrated with genomic architecture.  
(B) Promoter/TSS-proximal representation and median absolute TSS distance for motif-positive versus motif-negative strongly coupled Muscle peaks for the KLF9 motif and the prioritized SP-like STREME3 motif.  
(C) Peak- and gene-level overlap between the Muscle KLF9 and SP-like motif-associated candidate programs.  
(D) Peak- and gene-level overlap among the prioritized Spleen GC-rich 1, GC-rich 2, and KLF12-associated programs.  
(E) Promoter- and TSS-centering enrichment ratios for prioritized Muscle and Spleen motif-associated programs. Filled points denote FDR < 0.05 for the corresponding architecture comparison. FIMO q <= 0.05 indicates high-confidence sequence matches and should not be interpreted as TF-binding validation. Known- and de novo-motif analyses did not yield motifs meeting the primary AME E < 0.05 or independent STREME holdout E < 0.05 criteria; therefore, these panels represent prioritized sequence-architecture candidates rather than validated TF occupancy.

## Figure 6. Representative integrated candidate TF-motif-peak-gene networks.

(A) Representative Muscle promoter-centered candidate network integrating KLF9- and SP-like motif-associated modules.  
(B) Representative Spleen network illustrating the overlapping GC-rich and KLF12-associated candidate architecture.  
(C) Representative Liver supporting GC-rich candidate network. Candidate TF-to-motif edges represent motif-family annotation, motif-to-peak edges represent FIMO sequence matches, and peak-to-gene edges represent strongly coupled UROPA-associated peak-gene candidates. Filled peak symbols indicate promoter/TSS-proximal peaks, and highlighted associated genes are shared across selected modules where applicable. The networks summarize integrated candidate relationships and do not demonstrate TF occupancy or causal transcriptional regulation.
CAPTIONS


# ============================================================
# 4. Terminology lock
# ============================================================

cat > "${OUTDIR}/10F1_terminology_lock.tsv" <<'TSV'
ConceptLocked_termAvoid
702 pairsstrongly coupled peak-gene candidate pairsvalidated regulatory pairs
4674 pairsconcordant high-confidence tissue-specific peak-gene associationsregulatory interactions
Peak-gene relationUROPA-associated peak-gene candidatedirect target
FIMOFIMO sequence matchTF binding
TF-motif relationcandidate TF-to-motif annotationTF occupancy
Motif familySP-like or GC-rich KLF/SP/ZNF-like motif familyexact TF identity
Promoterpromoter/TSS-proximalpromoter-bound TF
TSSrepresentative-transcript TSStissue-specific dominant TSS
TauTau tissue-specificity scoreexpression significance
SPMspecificity measure (SPM)independent validation of Tau
Permutationsame-tissue permutation nullrandom unrestricted permutation
TSV


# ============================================================
# 5. Technical QC of authoritative files
# ============================================================

printf \
"Figure\tPDF_exists\tPDF_pages\tTIFF_exists\tTIFF_Xdpi\tTIFF_Ydpi\tFonts_embedded\tStatus\n" \
> "${OUTDIR}/10F1_technical_lock_QC.tsv"

tail -n +2 "${OUTDIR}/10F1_main_figure_lock_manifest.tsv" |
while IFS=$'\t' read -r fig panels pdf tif layout bio lock
do
    pdf_exists=FALSE
    tif_exists=FALSE
    pages=NA
    xdpi=NA
    ydpi=NA
    fonts=NA
    status=PASS

    if [[ -f "$pdf" ]]; then
        pdf_exists=TRUE
        pages=$(pdfinfo "$pdf" | awk '/^Pages:/ {print $2}')
    else
        status=FAIL
    fi

    if [[ -f "$tif" ]]; then
        tif_exists=TRUE
        read -r xdpi ydpi <<< "$(
            identify -format '%x %y' "$tif" |
            sed 's/ PixelsPerInch//g'
        )"
    else
        status=FAIL
    fi

    if [[ "$pdf_exists" == "TRUE" ]]; then
        if pdffonts "$pdf" |
           awk 'NR>2 && NF>0 {if($6!="yes") bad=1} END{exit bad ? 1 : 0}'
        then
            fonts=TRUE
        else
            fonts=FALSE
            status=FAIL
        fi
    fi

    if [[ "$pages" != "1" ]]; then
        status=FAIL
    fi

    if [[ "$xdpi" != "NA" ]]; then
        dpi_ok=$(awk -v x="$xdpi" -v y="$ydpi" \
            'BEGIN{
                if(x>=599 && x<=601 && y>=599 && y<=601) print "TRUE";
                else print "FALSE"
            }')
        if [[ "$dpi_ok" != "TRUE" ]]; then
            status=FAIL
        fi
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$fig" "$pdf_exists" "$pages" "$tif_exists" \
        "$xdpi" "$ydpi" "$fonts" "$status" \
        >> "${OUTDIR}/10F1_technical_lock_QC.tsv"
done


# ============================================================
# 6. Overall lock status
# ============================================================

fail_n=$(
    awk -F'\t' '
    NR>1 && $8!="PASS" {n++}
    END{print n+0}
    ' "${OUTDIR}/10F1_technical_lock_QC.tsv"
)

caption_n=$(
    grep -c '^## Figure [1-6]\.' \
    "${OUTDIR}/10F1_main_figure_captions.md"
)

panel_n=$(
    awk 'END{print NR-1}' \
    "${OUTDIR}/10F1_panel_semantic_manifest.tsv"
)

cat > "${OUTDIR}/10F1_overall_status.tsv" <<TSV
MetricValueExpectedStatus
Authoritative_main_figures66PASS
Final_caption_blocks${caption_n}6$([[ "$caption_n" -eq 6 ]] && echo PASS || echo FAIL)
Panel_manifest_rows${panel_n}29$([[ "$panel_n" -eq 29 ]] && echo PASS || echo FAIL)
Technical_QC_failures${fail_n}0$([[ "$fail_n" -eq 0 ]] && echo PASS || echo FAIL)
Biological_results_modified00PASS
Figures_redrawn_in_Step10F100PASS
TSV


echo
echo "============================================================"
echo "STEP10F1 COMPLETED"
echo "============================================================"

echo
echo "Authoritative figure lock:"
column -t -s $'\t' \
    "${OUTDIR}/10F1_main_figure_lock_manifest.tsv"

echo
echo "Technical QC:"
column -t -s $'\t' \
    "${OUTDIR}/10F1_technical_lock_QC.tsv"

echo
echo "Overall status:"
column -t -s $'\t' \
    "${OUTDIR}/10F1_overall_status.tsv"

echo
echo "Caption file:"
echo "${OUTDIR}/10F1_main_figure_captions.md"

echo
echo "IMPORTANT:"
echo "- Figure1-Figure6 are treated as visually locked."
echo "- Scientific Unicode symbols are allowed."
echo "- Different aspect ratios are retained intentionally."
echo "- No PDF or TIFF was modified."
echo "- No biological analysis was rerun."
echo "- No submission package was created."
echo "- Step10F2 must not start until Step10F1 QC is PASS."
echo "============================================================"

if [[ "$fail_n" -ne 0 || "$caption_n" -ne 6 || "$panel_n" -ne 29 ]]; then
    exit 1
fi
