#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import sys

print("=" * 84)
print("STEP10F4A2 — FINAL SUPPLEMENTARY TABLE COMPONENT SOURCE LOCK")
print("Exact frozen sources only")
print("No biological analysis rerun")
print("No statistic recomputed")
print("No Excel workbook generated")
print("=" * 84)

ROOT = Path(".")

OUTDIR = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F4_Supplementary_Tables/"
    "Step10F4A2_final_source_lock"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================================
# Helpers
# ======================================================================

def md5(path):

    h = hashlib.md5()

    with open(path, "rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def count_table_rows(path):

    if path.suffix.lower() not in {
        ".tsv",
        ".csv"
    }:
        return None

    delimiter = (
        "\t"
        if path.suffix.lower() == ".tsv"
        else ","
    )

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        lines = [
            x.rstrip("\n\r")
            for x in f
            if x.strip()
            and not x.startswith("#")
        ]

    if not lines:
        return 0

    return len(lines) - 1


def get_header(path):

    if path.suffix.lower() not in {
        ".tsv",
        ".csv"
    }:
        return ""

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for line in f:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            return line.rstrip(
                "\n\r"
            )

    return ""


# ======================================================================
# Frozen source specification
# ======================================================================

SOURCES = [

# ----------------------------------------------------------------------
# TABLE S1
# ----------------------------------------------------------------------

{
    "Table": "S1",
    "Sheet": "Sample_metadata",
    "Role": "PRIMARY",
    "Path": "sample_pair_metadata.tsv",
    "Expected_rows": 16,
    "Description":
        "Paired RNA-seq and ATAC-seq biological sample metadata."
},


# ----------------------------------------------------------------------
# TABLE S2
# QC + multi-tissue effect + tissue specificity + SPM
# ----------------------------------------------------------------------

{
    "Table": "S2",
    "Sheet": "RNA_filtering",
    "Role": "SUMMARY",
    "Path":
        "01_sample_QC/RNA/RNA_filtering_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "RNA-seq filtering summary."
},

{
    "Table": "S2",
    "Sheet": "ATAC_filtering",
    "Role": "SUMMARY",
    "Path":
        "01_sample_QC/ATAC/ATAC_filtering_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "ATAC-seq filtering summary."
},

{
    "Table": "S2",
    "Sheet": "Multi_tissue_LRT",
    "Role": "SUMMARY",
    "Path":
        "02_multi_tissue_effect/"
        "RNA_ATAC_multi_tissue_summary.tsv",
    "Expected_rows": 2,
    "Description":
        "RNA/ATAC multi-tissue DESeq2 LRT summary."
},

{
    "Table": "S2",
    "Sheet": "Tissue_specificity",
    "Role": "SUMMARY",
    "Path":
        "03_tissue_specificity/"
        "RNA_ATAC_tissue_specificity_summary.tsv",
    "Expected_rows": 2,
    "Description":
        "RNA/ATAC tissue-specificity summary."
},

{
    "Table": "S2",
    "Sheet": "RNA_tissue_counts",
    "Role": "SUMMARY",
    "Path":
        "03_tissue_specificity/RNA/"
        "RNA_tissue_specific_counts.tsv",
    "Expected_rows": 8,
    "Description":
        "RNA high-confidence tissue-specific gene counts."
},

{
    "Table": "S2",
    "Sheet": "ATAC_tissue_counts",
    "Role": "SUMMARY",
    "Path":
        "03_tissue_specificity/ATAC/"
        "ATAC_tissue_specific_counts.tsv",
    "Expected_rows": 8,
    "Description":
        "ATAC high-confidence tissue-specific peak counts."
},

{
    "Table": "S2",
    "Sheet": "RNA_SPM_summary",
    "Role": "SUPPORTING",
    "Path":
        "03B_SPM_analysis/RNA/"
        "RNA_SPM_summary.tsv",
    "Expected_rows": 8,
    "Description":
        "RNA SPM tissue-preference summary."
},

{
    "Table": "S2",
    "Sheet": "ATAC_SPM_summary",
    "Role": "SUPPORTING",
    "Path":
        "03B_SPM_analysis/ATAC/"
        "ATAC_SPM_summary.tsv",
    "Expected_rows": 8,
    "Description":
        "ATAC SPM tissue-preference summary."
},

{
    "Table": "S2",
    "Sheet": "RNA_threshold_sensitivity",
    "Role": "SUPPORTING",
    "Path":
        "03B_SPM_analysis/RNA/"
        "RNA_SPM_threshold_sensitivity.tsv",
    "Expected_rows": 5,
    "Description":
        "RNA tissue-specificity threshold sensitivity."
},

{
    "Table": "S2",
    "Sheet": "ATAC_threshold_sensitivity",
    "Role": "SUPPORTING",
    "Path":
        "03B_SPM_analysis/ATAC/"
        "ATAC_SPM_threshold_sensitivity.tsv",
    "Expected_rows": 5,
    "Description":
        "ATAC tissue-specificity threshold sensitivity."
},


# ----------------------------------------------------------------------
# TABLE S3
# ----------------------------------------------------------------------

{
    "Table": "S3",
    "Sheet": "Strong_peak_gene_pairs",
    "Role": "PRIMARY",
    "Path":
        "09B_702_pair_multi_evidence_master/"
        "09B1_702_strong_pair_multi_evidence_master.tsv",
    "Expected_rows": 702,
    "Description":
        "Frozen 702 strongly coupled tissue-specific peak-gene associations."
},


# ----------------------------------------------------------------------
# TABLE S4
# Structure + true promoter-first + TSS architecture
# ----------------------------------------------------------------------

{
    "Table": "S4",
    "Sheet": "Strong_structure",
    "Role": "PRIMARY",
    "Path":
        "06_strong_peak_gene_structure/"
        "06_702_strong_pairs_with_structure.tsv",
    "Expected_rows": 702,
    "Description":
        "Frozen structural annotation for the 702 strong peak-gene pairs."
},

{
    "Table": "S4",
    "Sheet": "UROPA_location_summary",
    "Role": "SUMMARY",
    "Path":
        "06_strong_peak_gene_structure/"
        "06_UROPA_relative_location_summary.tsv",
    "Expected_rows": 5,
    "Description":
        "UROPA relative-location composition."
},

{
    "Table": "S4",
    "Sheet": "Structural_class_summary",
    "Role": "SUMMARY",
    "Path":
        "06_strong_peak_gene_structure/"
        "06_structural_class_summary.tsv",
    "Expected_rows": 6,
    "Description":
        "Pair-specific structural-class composition."
},

{
    "Table": "S4",
    "Sheet": "Structural_enrichment",
    "Role": "SUMMARY",
    "Path":
        "06B_structural_enrichment/"
        "06B_structural_enrichment_tests.tsv",
    "Expected_rows": 18,
    "Description":
        "Frozen structural enrichment tests."
},

{
    "Table": "S4",
    "Sheet": "Promoter_first_702",
    "Role": "PRIMARY",
    "Path":
        "07_promoter_first_enrichment/"
        "07B_702_strong_promoter_first.tsv",
    "Expected_rows": 702,
    "Description":
        "True promoter-first and TSS annotation of the 702 strong pairs."
},

{
    "Table": "S4",
    "Sheet": "Promoter_class_summary",
    "Role": "SUMMARY",
    "Path":
        "07_promoter_first_enrichment/"
        "07B_702_class_peak_gene_summary.tsv",
    "Expected_rows": 7,
    "Description":
        "Promoter-first genomic-class composition of the 702 pairs."
},

{
    "Table": "S4",
    "Sheet": "TSS_distance_summary",
    "Role": "SUMMARY",
    "Path":
        "07_promoter_first_enrichment/"
        "07B_TSS_distance_summary.tsv",
    "Expected_rows": 3,
    "Description":
        "TSS-distance summary for All, HC and Strong cohorts."
},

{
    "Table": "S4",
    "Sheet": "Promoter_summary",
    "Role": "SUMMARY",
    "Path":
        "07_promoter_first_enrichment/"
        "07B_promoter_vs_nonpromoter_summary.tsv",
    "Expected_rows": 3,
    "Description":
        "Promoter versus non-promoter summary."
},

{
    "Table": "S4",
    "Sheet": "Global_enrichment",
    "Role": "SUMMARY",
    "Path":
        "07_promoter_first_enrichment/"
        "07B_global_genomic_enrichment_tests.tsv",
    "Expected_rows": 18,
    "Description":
        "Global promoter-first genomic enrichment tests."
},

{
    "Table": "S4",
    "Sheet": "Pair_specific_enrichment",
    "Role": "SUMMARY",
    "Path":
        "07_promoter_first_enrichment/"
        "07B_pair_specific_enrichment_tests.tsv",
    "Expected_rows": 21,
    "Description":
        "Pair-specific promoter-first enrichment tests."
},


# ----------------------------------------------------------------------
# TABLE S5
# Motif-analysis evidence and QC
# ----------------------------------------------------------------------

{
    "Table": "S5",
    "Sheet": "Matched_background_QC",
    "Role": "PRIMARY_QC",
    "Path":
        "08A_matched_motif_background_v2/"
        "08A_v2_matching_QC_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "GC/length/genomic-class matched-background QC."
},

{
    "Table": "S5",
    "Sheet": "AME_QC",
    "Role": "SUMMARY",
    "Path":
        "08B1c_AME_final_summary/"
        "08B1c_AME_QC_and_significance_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "Final AME QC and formal significance summary."
},

{
    "Table": "S5",
    "Sheet": "AME_all_results",
    "Role": "SUPPORTING",
    "Path":
        "08B1c_AME_final_summary/"
        "08B1c_all_tissues_AME_results.tsv",
    "Expected_rows": 3515,
    "Description":
        "Final cross-tissue AME motif result table."
},

{
    "Table": "S5",
    "Sheet": "STREME_default_summary",
    "Role": "SUMMARY",
    "Path":
        "08B2_STREME_de_novo/"
        "08B2_STREME_run_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "Default STREME run summary."
},

{
    "Table": "S5",
    "Sheet": "STREME_noholdout_summary",
    "Role": "SUPPORTING",
    "Path":
        "08B2b_STREME_noholdout/"
        "08B2b_STREME_run_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "No-holdout STREME sensitivity run summary."
},

{
    "Table": "S5",
    "Sheet": "Tomtom_repro_summary",
    "Role": "SUMMARY",
    "Path":
        "08B4B_STREME_reproducibility/"
        "08B4B_reproducibility_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "Default versus no-holdout STREME Tomtom reproducibility summary."
},


# Raw STREME text sources.
# Existence is locked, but row count is not used as a scientific criterion.

{
    "Table": "S5",
    "Sheet": "STREME_default_Cerebellum",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2_STREME_de_novo/Cerebellum/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen default STREME Cerebellum text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_default_Liver",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2_STREME_de_novo/Liver/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen default STREME Liver text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_default_Muscle",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2_STREME_de_novo/Muscle/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen default STREME Muscle text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_default_Spleen",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2_STREME_de_novo/Spleen/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen default STREME Spleen text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_noholdout_Cerebellum",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2b_STREME_noholdout/Cerebellum/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen no-holdout STREME Cerebellum text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_noholdout_Liver",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2b_STREME_noholdout/Liver/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen no-holdout STREME Liver text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_noholdout_Muscle",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2b_STREME_noholdout/Muscle/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen no-holdout STREME Muscle text result."
},

{
    "Table": "S5",
    "Sheet": "STREME_noholdout_Spleen",
    "Role": "RAW_TEXT_SOURCE",
    "Path":
        "08B2b_STREME_noholdout/Spleen/streme.txt",
    "Expected_rows": None,
    "Description":
        "Frozen no-holdout STREME Spleen text result."
},


# Raw Tomtom result sources: 5 + 1 + 1 + 2 = 9 true rows.

{
    "Table": "S5",
    "Sheet": "Tomtom_Cerebellum",
    "Role": "RAW_RESULT",
    "Path":
        "08B4B_STREME_reproducibility/"
        "Cerebellum/tomtom.tsv",
    "Expected_rows": 5,
    "Description":
        "Cerebellum cross-run Tomtom result."
},

{
    "Table": "S5",
    "Sheet": "Tomtom_Liver",
    "Role": "RAW_RESULT",
    "Path":
        "08B4B_STREME_reproducibility/"
        "Liver/tomtom.tsv",
    "Expected_rows": 1,
    "Description":
        "Liver cross-run Tomtom result."
},

{
    "Table": "S5",
    "Sheet": "Tomtom_Muscle",
    "Role": "RAW_RESULT",
    "Path":
        "08B4B_STREME_reproducibility/"
        "Muscle/tomtom.tsv",
    "Expected_rows": 1,
    "Description":
        "Muscle cross-run Tomtom result."
},

{
    "Table": "S5",
    "Sheet": "Tomtom_Spleen",
    "Role": "RAW_RESULT",
    "Path":
        "08B4B_STREME_reproducibility/"
        "Spleen/tomtom.tsv",
    "Expected_rows": 2,
    "Description":
        "Spleen cross-run Tomtom result."
},


# ----------------------------------------------------------------------
# TABLE S6
# FIMO sequence-match evidence
# ----------------------------------------------------------------------

{
    "Table": "S6",
    "Sheet": "FIMO_q005_links",
    "Role": "PRIMARY",
    "Path":
        "08E2_motif_peak_gene_TF_network/"
        "08E2B_high_confidence_q005_motif_peak_gene_links.tsv",
    "Expected_rows": 612,
    "Description":
        "High-confidence FIMO q<=0.05 motif-peak-gene sequence-match links."
},

{
    "Table": "S6",
    "Sheet": "Motif_target_summary",
    "Role": "SUMMARY",
    "Path":
        "08E2_motif_peak_gene_TF_network/"
        "08E2C_motif_level_target_summary.tsv",
    "Expected_rows": 11,
    "Description":
        "Motif-level target summary."
},

{
    "Table": "S6",
    "Sheet": "Tissue_network_summary",
    "Role": "SUMMARY",
    "Path":
        "08E2_motif_peak_gene_TF_network/"
        "08E2F_tissue_network_summary.tsv",
    "Expected_rows": 4,
    "Description":
        "Tissue-level motif network summary."
},

{
    "Table": "S6",
    "Sheet": "High_RNA_integrated",
    "Role": "SUPPORTING",
    "Path":
        "08E2_motif_peak_gene_TF_network/"
        "08E2E_high_RNA_TF_motif_peak_gene_master.tsv",
    "Expected_rows": 3178,
    "Description":
        "High-RNA TF motif-peak-gene integrated supporting master."
},

{
    "Table": "S6",
    "Sheet": "FIMO_QC",
    "Role": "QC",
    "Path":
        "08E2_motif_peak_gene_TF_network/"
        "08E2G_QC_summary.tsv",
    "Expected_rows": 9,
    "Description":
        "Frozen motif-network QC summary."
},


# ----------------------------------------------------------------------
# TABLE S7
# ----------------------------------------------------------------------

{
    "Table": "S7",
    "Sheet": "Candidate_TF_master",
    "Role": "PRIMARY",
    "Path":
        "08D_motif_family_TF_master/"
        "08D2_publication_TF_candidate_master.tsv",
    "Expected_rows": 48,
    "Description":
        "Frozen publication candidate TF RNA-support master."
},


# ----------------------------------------------------------------------
# TABLE S8
# Candidate modules
# ----------------------------------------------------------------------

{
    "Table": "S8",
    "Sheet": "Module_summary",
    "Role": "PRIMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C6_module_candidate_summary.tsv",
    "Expected_rows": 11,
    "Description":
        "Candidate regulatory module summary."
},

{
    "Table": "S8",
    "Sheet": "Module_definitions",
    "Role": "SUMMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C1_module_definition.tsv",
    "Expected_rows": 11,
    "Description":
        "Frozen module definitions."
},

{
    "Table": "S8",
    "Sheet": "Module_peak_gene_pairs",
    "Role": "PRIMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C2_all_module_peak_gene_candidates.tsv",
    "Expected_rows": 787,
    "Description":
        "All module-specific peak-gene candidate rows."
},

{
    "Table": "S8",
    "Sheet": "Gene_ranking",
    "Role": "PRIMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C3_gene_level_candidate_ranking.tsv",
    "Expected_rows": 596,
    "Description":
        "Module-specific gene-level candidate ranking."
},

{
    "Table": "S8",
    "Sheet": "Top10_representatives",
    "Role": "SUMMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C4_top10_representative_genes_per_module.tsv",
    "Expected_rows": 89,
    "Description":
        "Top representative genes per module."
},

{
    "Table": "S8",
    "Sheet": "Top5_primary_supporting",
    "Role": "SUMMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C5_top5_primary_supporting_figure_candidates.tsv",
    "Expected_rows": 30,
    "Description":
        "Top primary/supporting figure candidates."
},

{
    "Table": "S8",
    "Sheet": "Within_tissue_overlap",
    "Role": "SUMMARY",
    "Path":
        "09C_publication_candidate_modules/"
        "09C7_within_tissue_module_overlap.tsv",
    "Expected_rows": 11,
    "Description":
        "Within-tissue module overlap."
},

{
    "Table": "S8",
    "Sheet": "Representative_edges",
    "Role": "NETWORK",
    "Path":
        "09C_publication_candidate_modules/"
        "09C9_representative_figure_network_edges.tsv",
    "Expected_rows": 119,
    "Description":
        "Representative figure network edges."
},

{
    "Table": "S8",
    "Sheet": "Representative_nodes",
    "Role": "NETWORK",
    "Path":
        "09C_publication_candidate_modules/"
        "09C9_representative_figure_network_nodes.tsv",
    "Expected_rows": 73,
    "Description":
        "Representative figure network nodes."
},

]


# ======================================================================
# Validation
# ======================================================================

rows = []

hard_fail = False

for spec in SOURCES:

    path = ROOT / spec["Path"]

    exists = path.exists()

    observed_rows = (
        count_table_rows(path)
        if exists
        else None
    )

    expected = spec["Expected_rows"]

    if not exists:

        status = "FAIL"

    elif expected is None:

        status = "PASS"

    elif observed_rows == expected:

        status = "PASS"

    else:

        status = "FAIL"

    if status != "PASS":
        hard_fail = True

    rows.append({

        "Table":
            spec["Table"],

        "Sheet":
            spec["Sheet"],

        "Role":
            spec["Role"],

        "Path":
            spec["Path"],

        "Expected_rows":
            (
                expected
                if expected is not None
                else "EXISTENCE_ONLY"
            ),

        "Observed_rows":
            (
                observed_rows
                if observed_rows is not None
                else "NA"
            ),

        "Status":
            status,

        "MD5":
            (
                md5(path)
                if exists
                else ""
            ),

        "Header_preview":
            (
                get_header(path)[:600]
                if exists
                else ""
            ),

        "Description":
            spec["Description"]
    })


# ======================================================================
# Cross-source biological-count checks
# ======================================================================

checks = []


def add_check(metric, observed, expected):

    global hard_fail

    status = (
        "PASS"
        if observed == expected
        else "FAIL"
    )

    if status != "PASS":
        hard_fail = True

    checks.append({
        "Metric":
            metric,
        "Observed":
            observed,
        "Expected":
            expected,
        "Status":
            status
    })


# Core frozen counts.
add_check(
    "S1_sample_pairs",
    count_table_rows(
        ROOT / "sample_pair_metadata.tsv"
    ),
    16
)

add_check(
    "S3_strong_pairs",
    count_table_rows(
        ROOT /
        "09B_702_pair_multi_evidence_master/"
        "09B1_702_strong_pair_multi_evidence_master.tsv"
    ),
    702
)

add_check(
    "S4_promoter_first_strong_pairs",
    count_table_rows(
        ROOT /
        "07_promoter_first_enrichment/"
        "07B_702_strong_promoter_first.tsv"
    ),
    702
)

add_check(
    "S6_FIMO_q005_links",
    count_table_rows(
        ROOT /
        "08E2_motif_peak_gene_TF_network/"
        "08E2B_high_confidence_q005_motif_peak_gene_links.tsv"
    ),
    612
)

add_check(
    "S7_candidate_TFs",
    count_table_rows(
        ROOT /
        "08D_motif_family_TF_master/"
        "08D2_publication_TF_candidate_master.tsv"
    ),
    48
)

add_check(
    "S8_module_pairs",
    count_table_rows(
        ROOT /
        "09C_publication_candidate_modules/"
        "09C2_all_module_peak_gene_candidates.tsv"
    ),
    787
)

add_check(
    "S8_module_genes",
    count_table_rows(
        ROOT /
        "09C_publication_candidate_modules/"
        "09C3_gene_level_candidate_ranking.tsv"
    ),
    596
)

add_check(
    "S8_top10_representatives",
    count_table_rows(
        ROOT /
        "09C_publication_candidate_modules/"
        "09C4_top10_representative_genes_per_module.tsv"
    ),
    89
)

add_check(
    "S8_top5_primary_supporting",
    count_table_rows(
        ROOT /
        "09C_publication_candidate_modules/"
        "09C5_top5_primary_supporting_figure_candidates.tsv"
    ),
    30
)

add_check(
    "S8_representative_network_edges",
    count_table_rows(
        ROOT /
        "09C_publication_candidate_modules/"
        "09C9_representative_figure_network_edges.tsv"
    ),
    119
)

add_check(
    "S8_representative_network_nodes",
    count_table_rows(
        ROOT /
        "09C_publication_candidate_modules/"
        "09C9_representative_figure_network_nodes.tsv"
    ),
    73
)


# Tomtom true raw rows = 9
tomtom_total = sum(
    count_table_rows(
        ROOT / x
    )
    for x in [
        "08B4B_STREME_reproducibility/Cerebellum/tomtom.tsv",
        "08B4B_STREME_reproducibility/Liver/tomtom.tsv",
        "08B4B_STREME_reproducibility/Muscle/tomtom.tsv",
        "08B4B_STREME_reproducibility/Spleen/tomtom.tsv"
    ]
)

add_check(
    "S5_Tomtom_raw_result_rows",
    tomtom_total,
    9
)


# ======================================================================
# Write source manifest
# ======================================================================

manifest = (
    OUTDIR /
    "10F4A2_final_component_source_lock.tsv"
)

with open(
    manifest,
    "w",
    newline=""
) as f:

    fields = [
        "Table",
        "Sheet",
        "Role",
        "Path",
        "Expected_rows",
        "Observed_rows",
        "Status",
        "MD5",
        "Description",
        "Header_preview"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for row in rows:
        writer.writerow(row)


# ======================================================================
# Crosscheck file
# ======================================================================

check_file = (
    OUTDIR /
    "10F4A2_cross_source_QC.tsv"
)

with open(
    check_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Metric",
            "Observed",
            "Expected",
            "Status"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(checks)


# ======================================================================
# Workbook plan
# ======================================================================

plan_file = (
    OUTDIR /
    "10F4A2_final_workbook_plan.tsv"
)

with open(
    plan_file,
    "w"
) as f:

    f.write(
        "Table\tFinal_title\tPrimary_content\n"
    )

    f.write(
        "S1\t"
        "Paired multi-tissue sample metadata\t"
        "16 RNA-seq/ATAC-seq paired biological samples\n"
    )

    f.write(
        "S2\t"
        "RNA-seq and ATAC-seq QC, tissue effects and specificity\t"
        "QC + LRT + Tau/SPM summaries\n"
    )

    f.write(
        "S3\t"
        "Strongly coupled tissue-specific peak-gene associations\t"
        "702 frozen strong peak-gene pairs\n"
    )

    f.write(
        "S4\t"
        "Genomic architecture and promoter-first annotation\t"
        "Structural + true promoter/TSS architecture of 702 pairs\n"
    )

    f.write(
        "S5\t"
        "Motif analysis and matched-background evidence\t"
        "Matched background + AME + STREME + Tomtom\n"
    )

    f.write(
        "S6\t"
        "High-confidence motif-peak-gene sequence-match evidence\t"
        "Primary 612 FIMO q<=0.05 links + supporting TF integration\n"
    )

    f.write(
        "S7\t"
        "Candidate transcription-factor RNA support\t"
        "48 publication candidate TF records\n"
    )

    f.write(
        "S8\t"
        "Candidate regulatory modules and representative network support\t"
        "Module summary, 787 pairs, 596 genes and representative networks\n"
    )


# ======================================================================
# Overall
# ======================================================================

overall = (
    "PASS"
    if not hard_fail
    else "FAIL"
)

status_file = (
    OUTDIR /
    "10F4A2_overall_status.tsv"
)

with open(
    status_file,
    "w"
) as f:

    f.write(
        "Metric\tValue\n"
    )

    f.write(
        f"Final_component_source_lock\t{overall}\n"
    )

    f.write(
        f"Locked_source_files\t{len(rows)}\n"
    )

    f.write(
        "Biological_results_modified\t0\n"
    )


# ======================================================================
# Console
# ======================================================================

print()
print("=" * 84)
print("FINAL SOURCE LOCK")
print("=" * 84)

for table in [
    "S1",
    "S2",
    "S3",
    "S4",
    "S5",
    "S6",
    "S7",
    "S8"
]:

    subset = [
        x for x in rows
        if x["Table"] == table
    ]

    passed = sum(
        x["Status"] == "PASS"
        for x in subset
    )

    print(
        f"{table}: "
        f"{passed}/{len(subset)} sources PASS"
    )


print()
print("=" * 84)
print("CORE CROSS-SOURCE QC")
print("=" * 84)

for x in checks:

    print(
        f"{x['Metric']:<42} "
        f"{str(x['Observed']):>6} / "
        f"{str(x['Expected']):<6} "
        f"{x['Status']}"
    )


print()
print("=" * 84)

print(
    f"STEP10F4A2 STATUS: {overall}"
)

print(
    "Biological results modified: 0"
)

if overall == "PASS":

    print(
        "Supplementary Table S1-S8 source selection is FINAL and LOCKED."
    )

    print(
        "Ready for Step10F4B workbook construction."
    )

else:

    print(
        "Do NOT build Supplementary Tables until failed sources are resolved."
    )

print("=" * 84)

print("\nOutputs:")
print(manifest)
print(check_file)
print(plan_file)
print(status_file)
