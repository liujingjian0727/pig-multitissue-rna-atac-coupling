#!/usr/bin/env python3

# ============================================================
# STEP 09B
#
# Build the final 702-row strong RNA-ATAC peak-gene
# multi-evidence master table.
#
# One row = one official strong peak-associated gene pair.
#
# NO new statistical test is performed here.
# This step only integrates evidence already established in:
#
# Step05   RNA-ATAC quantitative coupling
# Step07   promoter-first / true-TSS annotation
# Step08D  motif family / TF candidate prioritization
# Step08E2 FIMO motif localization
# Step08E3A motif-family regulatory architecture
# Step09A  tissue-level functional enrichment context
#
# ------------------------------------------------------------
# Evidence tier:
#
# Tier0:
#   Strong coupling only
#
# Tier1:
#   Strong + reported FIMO motif (p <= 1e-4)
#
# Tier2:
#   Strong + FIMO q <= 0.05
#
# Tier3:
#   Tier2 + primary prioritized motif
#
# Tier4:
#   Tier3 + family-level structural support:
#     promoter enrichment FDR < 0.05
#       OR
#     significantly shorter TSS distance FDR < 0.05
#
# IMPORTANT:
#   Evidence tiers are NOT statistical significance tiers.
#
# IMPORTANT:
#   Step08 motif analysis was performed only in:
#     Muscle
#     Spleen
#     Liver
#     Cerebellum
#
#   Therefore Tier0 is subdivided into:
#     motif evaluated but no selected FIMO motif site
#     motif not evaluated in Step08
#
# ------------------------------------------------------------
# Peak-gene terminology:
#
# Peak-associated genes originate from the UROPA genomic
# association framework and are NOT experimentally validated
# regulatory targets.
#
# Candidate TF assignments, especially for de novo motifs,
# remain family-level candidates rather than demonstrated
# TF occupancy.
# ============================================================


import os
import re
import csv
import math
from collections import defaultdict, Counter


# ============================================================
# 0. Input files
# ============================================================

STEP05 = (
    "05_RNA_ATAC_quantitative_coupling/"
    "05_concordant_HC_peak_gene_correlations.tsv"
)

STEP07 = (
    "07_promoter_first_classification/"
    "07A_all_UROPA_assigned_promoter_first.tsv"
)

STEP08D1 = (
    "08D_motif_family_TF_master/"
    "08D1_publication_motif_family_master.tsv"
)

STEP08D2 = (
    "08D_motif_family_TF_master/"
    "08D2_publication_TF_candidate_master.tsv"
)

STEP08E2A = (
    "08E2_motif_peak_gene_TF_network/"
    "08E2A_all_reported_motif_peak_gene_links.tsv"
)

STEP08E2B = (
    "08E2_motif_peak_gene_TF_network/"
    "08E2B_high_confidence_q005_motif_peak_gene_links.tsv"
)

STEP08E2F = (
    "08E2_motif_peak_gene_TF_network/"
    "08E2F_tissue_network_summary.tsv"
)

STEP08E3A_PROMOTER = (
    "08E3A_motif_regulatory_architecture/"
    "08E3A2_promoter_enrichment.tsv"
)

STEP08E3A_CLASS = (
    "08E3A_motif_regulatory_architecture/"
    "08E3A3_pair_specific_class_enrichment.tsv"
)

STEP08E3A_TSS = (
    "08E3A_motif_regulatory_architecture/"
    "08E3A4_TSS_distance_comparison.tsv"
)

STEP09A = (
    "09A_strong_gene_functional_enrichment/"
    "09A11_enrichment_summary.tsv"
)

STEP09A_MANIFEST = (
    "09A_strong_gene_functional_enrichment/"
    "09A2_tissue_gene_set_manifest.tsv"
)

GTF = "Sus_longest.gtf"

OUTDIR = "09B_702_pair_multi_evidence_master"

os.makedirs(
    OUTDIR,
    exist_ok=True
)

EPS = 1e-12


# ============================================================
# 1. Helper functions
# ============================================================

def read_tsv(path):

    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Required file not found: {path}"
        )

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


def write_tsv(path, rows, fields):

    with open(
        path,
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:

            writer.writerow({
                field: row.get(field, "")
                for field in fields
            })


def as_float(x):

    if x is None:
        return None

    x = str(x).strip()

    if x == "":
        return None

    try:
        return float(x)

    except Exception:
        return None


def as_int(x):

    v = as_float(x)

    if v is None:
        return None

    return int(v)


def bool_value(x):

    if x is None:
        return False

    return str(x).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y"
    }


def safe_fraction(x):

    if x is None:
        return ""

    return f"{x:.8g}"


def split_multi(x):

    if x is None:
        return []

    x = str(x).strip()

    if x == "":
        return []

    out = []

    for item in re.split(
        r"[;,]",
        x
    ):

        item = item.strip()

        if item:
            out.append(item)

    return out


def unique_sorted(values):

    vals = set()

    for v in values:

        if v is None:
            continue

        v = str(v).strip()

        if v == "":
            continue

        vals.add(v)

    return sorted(vals)


def join_values(values):

    return ";".join(
        unique_sorted(values)
    )


def join_multi_fields(rows, field):

    vals = []

    for row in rows:

        vals.extend(
            split_multi(
                row.get(
                    field,
                    ""
                )
            )
        )

    return join_values(vals)


def first_nonempty(row, fields):

    for field in fields:

        if field in row:

            value = str(
                row.get(
                    field,
                    ""
                )
            ).strip()

            if value != "":
                return value

    return ""


def find_column(
    header,
    candidates,
    required=False
):

    for x in candidates:

        if x in header:
            return x

    if required:

        raise RuntimeError(
            "Could not find required column. "
            f"Tried: {candidates}"
        )

    return None


# ============================================================
# 2. Console header
# ============================================================

print("=" * 80)
print("STEP 09B")
print("702-pair RNA-ATAC multi-evidence master table")
print("=" * 80)

print(
    "No additional Python packages are required."
)


# ============================================================
# 3. Read Step05 and reconstruct official 702 strong pairs
# ============================================================

print("\nReading Step05 coupling table...")

step05 = read_tsv(
    STEP05
)

if len(step05) != 4674:

    raise RuntimeError(
        "Expected 4674 concordant-HC pairs in Step05, "
        f"observed {len(step05)}"
    )


required_step05 = [
    "SAF_peak_id",
    "BED_peak_id",
    "gene_id",
    "ATAC_Max_tissue",
    "RNA_Max_tissue",
    "rho_P348",
    "rho_P350",
    "rho_mean"
]


for col in required_step05:

    if col not in step05[0]:

        raise RuntimeError(
            f"Missing Step05 column: {col}"
        )


strong = []

for row in step05:

    r348 = as_float(
        row.get("rho_P348")
    )

    r350 = as_float(
        row.get("rho_P350")
    )

    rmean = as_float(
        row.get("rho_mean")
    )

    if None in (
        r348,
        r350,
        rmean
    ):
        continue

    if (
        r348 >= 0.5 - EPS
        and
        r350 >= 0.5 - EPS
        and
        rmean >= 0.7 - EPS
    ):

        row2 = dict(row)

        row2["Tissue"] = row[
            "ATAC_Max_tissue"
        ]

        strong.append(
            row2
        )


print(
    "Official strong pairs:",
    len(strong)
)


if len(strong) != 702:

    raise RuntimeError(
        "Expected exactly 702 strong pairs."
    )


strong_peak_ids = [
    x["SAF_peak_id"]
    for x in strong
]


if len(set(strong_peak_ids)) != 702:

    raise RuntimeError(
        "Strong SAF_peak_id is not unique."
    )


strong_keys = {
    (
        x["Tissue"],
        x["SAF_peak_id"]
    )
    for x in strong
}


# ============================================================
# 4. Strong peak multiplicity per associated gene
# ============================================================

gene_tissue_peak_count = Counter()
gene_all_peak_count = Counter()

for row in strong:

    gene_tissue_peak_count[
        (
            row["Tissue"],
            row["gene_id"]
        )
    ] += 1

    gene_all_peak_count[
        row["gene_id"]
    ] += 1


# ============================================================
# 5. Gene names from representative GTF
# ============================================================

print("\nReading representative GTF gene names...")

gene_name = {}

gene_id_re = re.compile(
    r'gene_id "([^"]+)"'
)

gene_name_re = re.compile(
    r'gene_name "([^"]+)"'
)


with open(
    GTF,
    "r",
    encoding="utf-8"
) as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip(
            "\n"
        ).split("\t")

        if len(fields) < 9:
            continue

        attrs = fields[8]

        gid_m = gene_id_re.search(
            attrs
        )

        if not gid_m:
            continue

        gid = gid_m.group(1)

        gname_m = gene_name_re.search(
            attrs
        )

        gname = (
            gname_m.group(1)
            if gname_m
            else ""
        )

        if gid not in gene_name:

            gene_name[gid] = gname


print(
    "GTF gene IDs:",
    len(gene_name)
)


# ============================================================
# 6. Step07 promoter-first / true-TSS annotation
# ============================================================

print("\nReading Step07 promoter-first annotation...")

step07 = read_tsv(
    STEP07
)


if (
    "SAF_peak_id" not in step07[0]
    or
    "gene_id" not in step07[0]
):

    raise RuntimeError(
        "Step07 requires SAF_peak_id and gene_id."
    )


step07_map = {}

for row in step07:

    key = (
        row["SAF_peak_id"],
        row["gene_id"]
    )

    if key not in step07_map:

        step07_map[
            key
        ] = row


step07_fields = [
    "Genomic_class",
    "Genomic_class_gene_id",
    "Genomic_class_transcript_id",
    "Genomic_class_overlap_bp",
    "Genomic_class_gene_matches_UROPA_gene",
    "Pair_specific_class",
    "Pair_class_transcript_id",
    "Pair_class_overlap_bp",
    "UROPA_gene_TSS_1based",
    "Peak_to_UROPA_TSS_min_abs_distance_bp",
    "Peak_center_to_UROPA_TSS_signed_bp"
]


step07_missing = 0

for row in strong:

    key = (
        row["SAF_peak_id"],
        row["gene_id"]
    )

    if key not in step07_map:

        step07_missing += 1


print(
    "Strong pairs without Step07 annotation:",
    step07_missing
)


if step07_missing != 0:

    raise RuntimeError(
        "Some strong pairs are missing Step07 annotation."
    )


# ============================================================
# 7. Read Step08D motif-family dictionary
# ============================================================

print("\nReading Step08D motif-family master...")

d1 = read_tsv(
    STEP08D1
)


required_d1 = [
    "Tissue",
    "Motif_family_ID",
    "Motif_family_label",
    "Representative_motif_ID",
    "Representative_motif_name",
    "Representative_source",
    "Consensus",
    "Candidate_TFs",
    "High_RNA_candidate_TFs",
    "Tissue_preferred_candidate_TFs",
    "FIMO_tier"
]


for col in required_d1:

    if col not in d1[0]:

        raise RuntimeError(
            f"Missing Step08D1 column: {col}"
        )


motif_dictionary = {}

for row in d1:

    key = (
        row["Tissue"],
        row["Representative_motif_ID"]
    )

    motif_dictionary[
        key
    ] = row


print(
    "Motif families in Step08D:",
    len(d1)
)


# ============================================================
# 8. Read Step08E2 motif-positive links
# ============================================================

print("\nReading Step08E2 motif localization...")

reported_links = read_tsv(
    STEP08E2A
)

q005_links = read_tsv(
    STEP08E2B
)


for data_name, data in [
    ("reported", reported_links),
    ("q005", q005_links)
]:

    for col in [
        "Tissue",
        "Motif_ID",
        "SAF_peak_id"
    ]:

        if col not in data[0]:

            raise RuntimeError(
                f"Step08E2 {data_name} links "
                f"missing column: {col}"
            )


print(
    "Reported motif-peak-gene links:",
    len(reported_links)
)

print(
    "q<=0.05 motif-peak-gene links:",
    len(q005_links)
)


if len(reported_links) != 902:

    raise RuntimeError(
        "Expected 902 reported motif-peak-gene links."
    )


if len(q005_links) != 612:

    raise RuntimeError(
        "Expected 612 q<=0.05 motif-peak-gene links."
    )


# ------------------------------------------------------------
# Validate all motif links map to official strong peaks.
# ------------------------------------------------------------

unmapped_reported = []

for row in reported_links:

    key = (
        row["Tissue"],
        row["SAF_peak_id"]
    )

    if key not in strong_keys:

        unmapped_reported.append(
            key
        )


unmapped_q005 = []

for row in q005_links:

    key = (
        row["Tissue"],
        row["SAF_peak_id"]
    )

    if key not in strong_keys:

        unmapped_q005.append(
            key
        )


print(
    "Reported motif links not mapping to strong pairs:",
    len(unmapped_reported)
)

print(
    "q005 motif links not mapping to strong pairs:",
    len(unmapped_q005)
)


if unmapped_reported or unmapped_q005:

    raise RuntimeError(
        "Step08E2 contains motif links outside "
        "the official strong-pair set."
    )


# ------------------------------------------------------------
# q005 must be a subset of reported links.
# ------------------------------------------------------------

reported_link_keys = {
    (
        row["Tissue"],
        row["SAF_peak_id"],
        row["Motif_ID"]
    )
    for row in reported_links
}


q005_link_keys = {
    (
        row["Tissue"],
        row["SAF_peak_id"],
        row["Motif_ID"]
    )
    for row in q005_links
}


if not q005_link_keys.issubset(
    reported_link_keys
):

    raise RuntimeError(
        "q005 motif links are not a subset "
        "of reported motif links."
    )


# ============================================================
# 9. Step08E2 tissue scope
# ============================================================

e2f = read_tsv(
    STEP08E2F
)

motif_scope_tissues = {
    row["Tissue"]
    for row in e2f
}


print(
    "Step08 motif-analysis tissues:",
    ", ".join(
        sorted(
            motif_scope_tissues
        )
    )
)


expected_reported_unique = sum(
    as_int(
        row.get(
            "Unique_FIMO_positive_peaks"
        )
    ) or 0
    for row in e2f
)


expected_q005_unique = sum(
    as_int(
        row.get(
            "Unique_Q005_peaks"
        )
    ) or 0
    for row in e2f
)


observed_reported_unique = len({
    (
        row["Tissue"],
        row["SAF_peak_id"]
    )
    for row in reported_links
})


observed_q005_unique = len({
    (
        row["Tissue"],
        row["SAF_peak_id"]
    )
    for row in q005_links
})


print(
    "Unique reported motif-positive strong peaks:",
    observed_reported_unique
)

print(
    "Unique q005 motif-positive strong peaks:",
    observed_q005_unique
)


if (
    expected_reported_unique
    and
    observed_reported_unique !=
        expected_reported_unique
):

    raise RuntimeError(
        "Reported unique-peak count disagrees "
        "with Step08E2F."
    )


if (
    expected_q005_unique
    and
    observed_q005_unique !=
        expected_q005_unique
):

    raise RuntimeError(
        "q005 unique-peak count disagrees "
        "with Step08E2F."
    )


# ============================================================
# 10. Index motif links per strong peak
# ============================================================

reported_by_peak = defaultdict(
    list
)

q005_by_peak = defaultdict(
    list
)


for row in reported_links:

    reported_by_peak[
        (
            row["Tissue"],
            row["SAF_peak_id"]
        )
    ].append(
        row
    )


for row in q005_links:

    q005_by_peak[
        (
            row["Tissue"],
            row["SAF_peak_id"]
        )
    ].append(
        row
    )


# ============================================================
# 11. Read Step08E3A architecture results
# ============================================================

print(
    "\nReading Step08E3A motif-family architecture evidence..."
)

promoter_rows = read_tsv(
    STEP08E3A_PROMOTER
)

class_rows = read_tsv(
    STEP08E3A_CLASS
)

tss_rows = read_tsv(
    STEP08E3A_TSS
)


promoter_map = {}
tss_map = {}
distal_map = {}


for row in promoter_rows:

    key = (
        row["Tissue"],
        row["Motif_ID"],
        row["Evidence_level"]
    )

    pfrac = as_float(
        row.get(
            "Positive_fraction"
        )
    )

    nfrac = as_float(
        row.get(
            "Negative_fraction"
        )
    )

    fdr = as_float(
        row.get(
            "FDR"
        )
    )

    promoter_map[
        key
    ] = {
        "FDR":
            fdr,

        "Positive_fraction":
            pfrac,

        "Negative_fraction":
            nfrac,

        "OR":
            as_float(
                row.get(
                    "Fisher_OR"
                )
            ),

        "enriched":
            (
                fdr is not None
                and
                fdr < 0.05
                and
                pfrac is not None
                and
                nfrac is not None
                and
                pfrac > nfrac
            )
    }


for row in tss_rows:

    key = (
        row["Tissue"],
        row["Motif_ID"],
        row["Evidence_level"]
    )

    pos_med = as_float(
        row.get(
            "Positive_median_abs_TSS_bp"
        )
    )

    neg_med = as_float(
        row.get(
            "Negative_median_abs_TSS_bp"
        )
    )

    fdr = as_float(
        row.get(
            "FDR"
        )
    )

    tss_map[
        key
    ] = {
        "FDR":
            fdr,

        "Positive_median":
            pos_med,

        "Negative_median":
            neg_med,

        "closer":
            (
                fdr is not None
                and
                fdr < 0.05
                and
                pos_med is not None
                and
                neg_med is not None
                and
                pos_med < neg_med
            )
    }


for row in class_rows:

    if (
        row.get(
            "Pair_specific_class"
        )
        !=
        "Distal_to_associated_gene"
    ):
        continue

    key = (
        row["Tissue"],
        row["Motif_ID"],
        row["Evidence_level"]
    )

    pfrac = as_float(
        row.get(
            "Positive_fraction"
        )
    )

    nfrac = as_float(
        row.get(
            "Negative_fraction"
        )
    )

    fdr = as_float(
        row.get(
            "FDR"
        )
    )

    distal_map[
        key
    ] = {
        "FDR":
            fdr,

        "Positive_fraction":
            pfrac,

        "Negative_fraction":
            nfrac,

        "depleted":
            (
                fdr is not None
                and
                fdr < 0.05
                and
                pfrac is not None
                and
                nfrac is not None
                and
                pfrac < nfrac
            )
    }


# ============================================================
# 12. Helper: motif metadata + structural evidence
# ============================================================

def get_motif_info(
    tissue,
    motif_id,
    evidence_level
):

    d = motif_dictionary.get(
        (
            tissue,
            motif_id
        ),
        {}
    )


    p = promoter_map.get(
        (
            tissue,
            motif_id,
            evidence_level
        ),
        {}
    )


    t = tss_map.get(
        (
            tissue,
            motif_id,
            evidence_level
        ),
        {}
    )


    dist = distal_map.get(
        (
            tissue,
            motif_id,
            evidence_level
        ),
        {}
    )


    promoter_support = bool(
        p.get(
            "enriched",
            False
        )
    )


    tss_support = bool(
        t.get(
            "closer",
            False
        )
    )


    distal_depletion = bool(
        dist.get(
            "depleted",
            False
        )
    )


    return {

        "Motif_ID":
            motif_id,

        "Motif_family_ID":
            d.get(
                "Motif_family_ID",
                ""
            ),

        "Motif_family_label":
            d.get(
                "Motif_family_label",
                ""
            ),

        "Motif_source":
            d.get(
                "Representative_source",
                ""
            ),

        "Motif_name":
            d.get(
                "Representative_motif_name",
                ""
            ),

        "Consensus":
            d.get(
                "Consensus",
                ""
            ),

        "FIMO_tier":
            d.get(
                "FIMO_tier",
                ""
            ),

        "Candidate_TFs":
            d.get(
                "Candidate_TFs",
                ""
            ),

        "High_RNA_candidate_TFs":
            d.get(
                "High_RNA_candidate_TFs",
                ""
            ),

        "Tissue_preferred_candidate_TFs":
            d.get(
                "Tissue_preferred_candidate_TFs",
                ""
            ),

        "AME_exploratory_support":
            d.get(
                "AME_exploratory_support",
                ""
            ),

        "AME_formal_E_lt_0.05":
            d.get(
                "AME_formal_E_lt_0.05",
                ""
            ),

        "Default_STREME_formal_E_lt_0.05":
            d.get(
                "Default_STREME_formal_E_lt_0.05",
                ""
            ),

        "CrossRun_reproduced":
            d.get(
                "CrossRun_reproduced",
                ""
            ),

        "Promoter_support":
            promoter_support,

        "Promoter_FDR":
            p.get(
                "FDR"
            ),

        "TSS_support":
            tss_support,

        "TSS_FDR":
            t.get(
                "FDR"
            ),

        "Distal_depletion":
            distal_depletion,

        "Distal_FDR":
            dist.get(
                "FDR"
            ),

        "Structural_support":
            (
                promoter_support
                or
                tss_support
            )
    }


# ============================================================
# 13. Validate motif dictionary coverage
# ============================================================

used_motifs = {
    (
        row["Tissue"],
        row["Motif_ID"]
    )
    for row in reported_links
}


missing_motif_dictionary = [
    x
    for x in sorted(
        used_motifs
    )
    if x not in motif_dictionary
]


print(
    "Selected FIMO motifs used:",
    len(
        used_motifs
    )
)

print(
    "Used motifs missing Step08D metadata:",
    len(
        missing_motif_dictionary
    )
)


if missing_motif_dictionary:

    raise RuntimeError(
        "Some Step08E2 motifs cannot be mapped "
        "to Step08D1."
    )


# ============================================================
# 14. Step09A tissue functional context
# ============================================================

print("\nReading Step09A tissue functional context...")

step09a = read_tsv(
    STEP09A
)

step09a_manifest = read_tsv(
    STEP09A_MANIFEST
)


step09a_map = {
    row["Tissue"]:
        row
    for row in step09a
}


step09a_manifest_map = {
    row["Tissue"]:
        row
    for row in step09a_manifest
}


# ============================================================
# 15. Build master rows
# ============================================================

print("\nBuilding 702-row multi-evidence master...")

master_rows = []


for base in strong:

    tissue = base[
        "Tissue"
    ]

    saf_peak = base[
        "SAF_peak_id"
    ]

    gid = base[
        "gene_id"
    ]

    key_peak = (
        tissue,
        saf_peak
    )


    # --------------------------------------------------------
    # Step07 pair annotation
    # --------------------------------------------------------

    s07 = step07_map[
        (
            saf_peak,
            gid
        )
    ]


    pair_class = s07.get(
        "Pair_specific_class",
        ""
    )


    abs_tss = as_float(
        s07.get(
            "Peak_to_UROPA_TSS_min_abs_distance_bp"
        )
    )


    is_promoter = (
        pair_class ==
        "Promoter_TSSproximal"
    )


    abs_tss_le_2kb = (
        abs_tss is not None
        and
        abs_tss <= 2000
    )


    # --------------------------------------------------------
    # Motif links
    # --------------------------------------------------------

    reported_peak_links = reported_by_peak.get(
        key_peak,
        []
    )


    q005_peak_links = q005_by_peak.get(
        key_peak,
        []
    )


    reported_motif_ids = unique_sorted(
        row["Motif_ID"]
        for row in reported_peak_links
    )


    q005_motif_ids = unique_sorted(
        row["Motif_ID"]
        for row in q005_peak_links
    )


    reported_info = [
        get_motif_info(
            tissue,
            motif,
            "Reported_p1e-4"
        )
        for motif in reported_motif_ids
    ]


    q005_info = [
        get_motif_info(
            tissue,
            motif,
            "High_confidence_q005"
        )
        for motif in q005_motif_ids
    ]


    primary_q005_info = [
        x
        for x in q005_info
        if str(
            x.get(
                "FIMO_tier",
                ""
            )
        ).lower()
        ==
        "primary"
    ]


    # --------------------------------------------------------
    # Structural motif support
    # --------------------------------------------------------

    reported_structural = [
        x
        for x in reported_info
        if x[
            "Structural_support"
        ]
    ]


    q005_structural = [
        x
        for x in q005_info
        if x[
            "Structural_support"
        ]
    ]


    primary_q005_structural = [
        x
        for x in primary_q005_info
        if x[
            "Structural_support"
        ]
    ]


    # --------------------------------------------------------
    # Evidence tier
    # --------------------------------------------------------

    if primary_q005_structural:

        evidence_tier = "Tier4"

        evidence_tier_numeric = 4

        evidence_tier_definition = (
            "Strong + primary q005 motif + "
            "family-level promoter/TSS structural support"
        )


    elif primary_q005_info:

        evidence_tier = "Tier3"

        evidence_tier_numeric = 3

        evidence_tier_definition = (
            "Strong + primary prioritized motif "
            "with FIMO q<=0.05"
        )


    elif q005_info:

        evidence_tier = "Tier2"

        evidence_tier_numeric = 2

        evidence_tier_definition = (
            "Strong + selected motif with "
            "FIMO q<=0.05"
        )


    elif reported_info:

        evidence_tier = "Tier1"

        evidence_tier_numeric = 1

        evidence_tier_definition = (
            "Strong + selected motif with "
            "FIMO p<=1e-4"
        )


    else:

        evidence_tier = "Tier0"

        evidence_tier_numeric = 0

        evidence_tier_definition = (
            "Strong coupling only"
        )


    # --------------------------------------------------------
    # Motif scope / Tier0 subtype
    # --------------------------------------------------------

    if tissue in motif_scope_tissues:

        motif_analysis_status = (
            "Evaluated_in_Step08"
        )

    else:

        motif_analysis_status = (
            "Not_evaluated_in_Step08"
        )


    if evidence_tier != "Tier0":

        tier0_subtype = ""

    elif tissue in motif_scope_tissues:

        tier0_subtype = (
            "Motif_evaluated_no_selected_reported_site"
        )

    else:

        tier0_subtype = (
            "Motif_not_evaluated_in_Step08"
        )


    # --------------------------------------------------------
    # Candidate TF aggregation
    # --------------------------------------------------------

    all_candidate_tfs = []

    all_high_rna_tfs = []

    all_preferred_tfs = []

    q005_candidate_tfs = []

    q005_high_rna_tfs = []

    q005_preferred_tfs = []

    primary_q005_candidate_tfs = []

    primary_q005_high_rna_tfs = []

    primary_q005_preferred_tfs = []


    for info in reported_info:

        all_candidate_tfs.extend(
            split_multi(
                info[
                    "Candidate_TFs"
                ]
            )
        )

        all_high_rna_tfs.extend(
            split_multi(
                info[
                    "High_RNA_candidate_TFs"
                ]
            )
        )

        all_preferred_tfs.extend(
            split_multi(
                info[
                    "Tissue_preferred_candidate_TFs"
                ]
            )
        )


    for info in q005_info:

        q005_candidate_tfs.extend(
            split_multi(
                info[
                    "Candidate_TFs"
                ]
            )
        )

        q005_high_rna_tfs.extend(
            split_multi(
                info[
                    "High_RNA_candidate_TFs"
                ]
            )
        )

        q005_preferred_tfs.extend(
            split_multi(
                info[
                    "Tissue_preferred_candidate_TFs"
                ]
            )
        )


    for info in primary_q005_info:

        primary_q005_candidate_tfs.extend(
            split_multi(
                info[
                    "Candidate_TFs"
                ]
            )
        )

        primary_q005_high_rna_tfs.extend(
            split_multi(
                info[
                    "High_RNA_candidate_TFs"
                ]
            )
        )

        primary_q005_preferred_tfs.extend(
            split_multi(
                info[
                    "Tissue_preferred_candidate_TFs"
                ]
            )
        )


    # --------------------------------------------------------
    # Step09A tissue functional context
    # --------------------------------------------------------

    fctx = step09a_map.get(
        tissue,
        {}
    )

    fmanifest = step09a_manifest_map.get(
        tissue,
        {}
    )


    # --------------------------------------------------------
    # Build output row
    # --------------------------------------------------------

    out = dict(
        base
    )


    out["gene_name"] = gene_name.get(
        gid,
        ""
    )


    out[
        "Strong_peak_N_for_gene_within_tissue"
    ] = gene_tissue_peak_count[
        (
            tissue,
            gid
        )
    ]


    out[
        "Strong_peak_N_for_gene_all_tissues"
    ] = gene_all_peak_count[
        gid
    ]


    # --------------------------------------------------------
    # Step07 prefixed columns
    # --------------------------------------------------------

    for field in step07_fields:

        out[
            "Step07_" + field
        ] = s07.get(
            field,
            ""
        )


    out[
        "Pair_is_Promoter_TSSproximal"
    ] = str(
        is_promoter
    )


    out[
        "Pair_abs_TSS_distance_le_2kb"
    ] = str(
        abs_tss_le_2kb
    )


    # --------------------------------------------------------
    # Motif-analysis scope
    # --------------------------------------------------------

    out[
        "Motif_analysis_status"
    ] = motif_analysis_status


    out[
        "Evidence_tier"
    ] = evidence_tier


    out[
        "Evidence_tier_numeric"
    ] = evidence_tier_numeric


    out[
        "Evidence_tier_definition"
    ] = evidence_tier_definition


    out[
        "Tier0_subtype"
    ] = tier0_subtype


    # --------------------------------------------------------
    # Reported motif evidence
    # --------------------------------------------------------

    out[
        "Reported_motif_N"
    ] = len(
        reported_info
    )


    out[
        "Reported_motif_IDs"
    ] = join_values(
        x[
            "Motif_ID"
        ]
        for x in reported_info
    )


    out[
        "Reported_motif_family_IDs"
    ] = join_values(
        x[
            "Motif_family_ID"
        ]
        for x in reported_info
    )


    out[
        "Reported_motif_family_labels"
    ] = join_values(
        x[
            "Motif_family_label"
        ]
        for x in reported_info
    )


    out[
        "Reported_motif_sources"
    ] = join_values(
        x[
            "Motif_source"
        ]
        for x in reported_info
    )


    out[
        "Reported_FIMO_tiers"
    ] = join_values(
        x[
            "FIMO_tier"
        ]
        for x in reported_info
    )


    # --------------------------------------------------------
    # q005 motif evidence
    # --------------------------------------------------------

    out[
        "Q005_motif_N"
    ] = len(
        q005_info
    )


    out[
        "Q005_motif_IDs"
    ] = join_values(
        x[
            "Motif_ID"
        ]
        for x in q005_info
    )


    out[
        "Q005_motif_family_IDs"
    ] = join_values(
        x[
            "Motif_family_ID"
        ]
        for x in q005_info
    )


    out[
        "Q005_motif_family_labels"
    ] = join_values(
        x[
            "Motif_family_label"
        ]
        for x in q005_info
    )


    out[
        "Q005_motif_sources"
    ] = join_values(
        x[
            "Motif_source"
        ]
        for x in q005_info
    )


    # --------------------------------------------------------
    # Primary q005 motif evidence
    # --------------------------------------------------------

    out[
        "Primary_q005_motif_N"
    ] = len(
        primary_q005_info
    )


    out[
        "Primary_q005_motif_IDs"
    ] = join_values(
        x[
            "Motif_ID"
        ]
        for x in primary_q005_info
    )


    out[
        "Primary_q005_motif_family_labels"
    ] = join_values(
        x[
            "Motif_family_label"
        ]
        for x in primary_q005_info
    )


    # --------------------------------------------------------
    # Candidate TF layers
    # --------------------------------------------------------

    out[
        "Reported_candidate_TFs"
    ] = join_values(
        all_candidate_tfs
    )


    out[
        "Reported_high_RNA_candidate_TFs"
    ] = join_values(
        all_high_rna_tfs
    )


    out[
        "Reported_tissue_preferred_candidate_TFs"
    ] = join_values(
        all_preferred_tfs
    )


    out[
        "Q005_candidate_TFs"
    ] = join_values(
        q005_candidate_tfs
    )


    out[
        "Q005_high_RNA_candidate_TFs"
    ] = join_values(
        q005_high_rna_tfs
    )


    out[
        "Q005_tissue_preferred_candidate_TFs"
    ] = join_values(
        q005_preferred_tfs
    )


    out[
        "Primary_q005_candidate_TFs"
    ] = join_values(
        primary_q005_candidate_tfs
    )


    out[
        "Primary_q005_high_RNA_candidate_TFs"
    ] = join_values(
        primary_q005_high_rna_tfs
    )


    out[
        "Primary_q005_tissue_preferred_candidate_TFs"
    ] = join_values(
        primary_q005_preferred_tfs
    )


    # --------------------------------------------------------
    # Family-level structural support
    # --------------------------------------------------------

    out[
        "Reported_structurally_supported_motif_IDs"
    ] = join_values(
        x[
            "Motif_ID"
        ]
        for x in reported_structural
    )


    out[
        "Q005_structurally_supported_motif_IDs"
    ] = join_values(
        x[
            "Motif_ID"
        ]
        for x in q005_structural
    )


    out[
        "Primary_q005_structurally_supported_motif_IDs"
    ] = join_values(
        x[
            "Motif_ID"
        ]
        for x in primary_q005_structural
    )


    out[
        "Any_reported_family_promoter_enrichment"
    ] = str(
        any(
            x[
                "Promoter_support"
            ]
            for x in reported_info
        )
    )


    out[
        "Any_q005_family_promoter_enrichment"
    ] = str(
        any(
            x[
                "Promoter_support"
            ]
            for x in q005_info
        )
    )


    out[
        "Any_reported_family_TSS_proximity"
    ] = str(
        any(
            x[
                "TSS_support"
            ]
            for x in reported_info
        )
    )


    out[
        "Any_q005_family_TSS_proximity"
    ] = str(
        any(
            x[
                "TSS_support"
            ]
            for x in q005_info
        )
    )


    out[
        "Any_reported_family_distal_depletion"
    ] = str(
        any(
            x[
                "Distal_depletion"
            ]
            for x in reported_info
        )
    )


    out[
        "Any_q005_family_distal_depletion"
    ] = str(
        any(
            x[
                "Distal_depletion"
            ]
            for x in q005_info
        )
    )


    out[
        "Primary_q005_has_structural_support"
    ] = str(
        len(
            primary_q005_structural
        ) > 0
    )


    # --------------------------------------------------------
    # Formal motif significance reminder
    # --------------------------------------------------------

    out[
        "Any_reported_AME_formal_E_lt_0.05"
    ] = str(
        any(
            bool_value(
                x[
                    "AME_formal_E_lt_0.05"
                ]
            )
            for x in reported_info
        )
    )


    out[
        "Any_reported_STREME_holdout_formal_E_lt_0.05"
    ] = str(
        any(
            bool_value(
                x[
                    "Default_STREME_formal_E_lt_0.05"
                ]
            )
            for x in reported_info
        )
    )


    # --------------------------------------------------------
    # Step09A tissue context
    # --------------------------------------------------------

    out[
        "Step09A_tissue_Strong_gene_N"
    ] = first_nonempty(
        fctx,
        [
            "Strong_gene_N"
        ]
    )


    out[
        "Step09A_tissue_Universe_gene_N"
    ] = first_nonempty(
        fctx,
        [
            "Universe_gene_N"
        ]
    )


    out[
        "Step09A_tissue_GO_BP_FDR005_N"
    ] = first_nonempty(
        fctx,
        [
            "GO_BP_FDR005_N"
        ]
    )


    out[
        "Step09A_tissue_GO_MF_FDR005_N"
    ] = first_nonempty(
        fctx,
        [
            "GO_MF_FDR005_N"
        ]
    )


    out[
        "Step09A_tissue_GO_CC_FDR005_N"
    ] = first_nonempty(
        fctx,
        [
            "GO_CC_FDR005_N"
        ]
    )


    out[
        "Step09A_tissue_KEGG_FDR005_N"
    ] = first_nonempty(
        fctx,
        [
            "KEGG_FDR005_N"
        ]
    )


    out[
        "Step09A_input_size_flag"
    ] = first_nonempty(
        fctx,
        [
            "Input_size_flag"
        ]
    )


    out[
        "Step09A_strong_gene_fraction_of_universe"
    ] = first_nonempty(
        fmanifest,
        [
            "Strong_gene_fraction_of_universe"
        ]
    )


    # --------------------------------------------------------
    # Candidate publication flags
    # --------------------------------------------------------

    out[
        "Publication_core_candidate"
    ] = str(
        evidence_tier_numeric >= 3
    )


    out[
        "Publication_structural_candidate"
    ] = str(
        evidence_tier_numeric >= 4
    )


    master_rows.append(
        out
    )


# ============================================================
# 16. Final row count / uniqueness QC
# ============================================================

if len(master_rows) != 702:

    raise RuntimeError(
        "Final master table does not contain 702 rows."
    )


if len({
    row["SAF_peak_id"]
    for row in master_rows
}) != 702:

    raise RuntimeError(
        "Final master SAF_peak_id is not unique."
    )


# ============================================================
# 17. Output column order
# ============================================================

step05_fields = list(
    step05[0].keys()
)


new_fields = [

    "Tissue",
    "gene_name",

    "Strong_peak_N_for_gene_within_tissue",
    "Strong_peak_N_for_gene_all_tissues",

    # Step07
    "Step07_Genomic_class",
    "Step07_Genomic_class_gene_id",
    "Step07_Genomic_class_transcript_id",
    "Step07_Genomic_class_overlap_bp",
    "Step07_Genomic_class_gene_matches_UROPA_gene",
    "Step07_Pair_specific_class",
    "Step07_Pair_class_transcript_id",
    "Step07_Pair_class_overlap_bp",
    "Step07_UROPA_gene_TSS_1based",
    "Step07_Peak_to_UROPA_TSS_min_abs_distance_bp",
    "Step07_Peak_center_to_UROPA_TSS_signed_bp",

    "Pair_is_Promoter_TSSproximal",
    "Pair_abs_TSS_distance_le_2kb",

    # Evidence tier
    "Motif_analysis_status",
    "Evidence_tier",
    "Evidence_tier_numeric",
    "Evidence_tier_definition",
    "Tier0_subtype",

    # Reported motif
    "Reported_motif_N",
    "Reported_motif_IDs",
    "Reported_motif_family_IDs",
    "Reported_motif_family_labels",
    "Reported_motif_sources",
    "Reported_FIMO_tiers",

    # q005 motif
    "Q005_motif_N",
    "Q005_motif_IDs",
    "Q005_motif_family_IDs",
    "Q005_motif_family_labels",
    "Q005_motif_sources",

    # Primary q005
    "Primary_q005_motif_N",
    "Primary_q005_motif_IDs",
    "Primary_q005_motif_family_labels",

    # Candidate TF
    "Reported_candidate_TFs",
    "Reported_high_RNA_candidate_TFs",
    "Reported_tissue_preferred_candidate_TFs",

    "Q005_candidate_TFs",
    "Q005_high_RNA_candidate_TFs",
    "Q005_tissue_preferred_candidate_TFs",

    "Primary_q005_candidate_TFs",
    "Primary_q005_high_RNA_candidate_TFs",
    "Primary_q005_tissue_preferred_candidate_TFs",

    # Structural motif evidence
    "Reported_structurally_supported_motif_IDs",
    "Q005_structurally_supported_motif_IDs",
    "Primary_q005_structurally_supported_motif_IDs",

    "Any_reported_family_promoter_enrichment",
    "Any_q005_family_promoter_enrichment",
    "Any_reported_family_TSS_proximity",
    "Any_q005_family_TSS_proximity",
    "Any_reported_family_distal_depletion",
    "Any_q005_family_distal_depletion",
    "Primary_q005_has_structural_support",

    # Formal motif statistics reminder
    "Any_reported_AME_formal_E_lt_0.05",
    "Any_reported_STREME_holdout_formal_E_lt_0.05",

    # Step09A tissue function context
    "Step09A_tissue_Strong_gene_N",
    "Step09A_tissue_Universe_gene_N",
    "Step09A_tissue_GO_BP_FDR005_N",
    "Step09A_tissue_GO_MF_FDR005_N",
    "Step09A_tissue_GO_CC_FDR005_N",
    "Step09A_tissue_KEGG_FDR005_N",
    "Step09A_input_size_flag",
    "Step09A_strong_gene_fraction_of_universe",

    # Publication use
    "Publication_core_candidate",
    "Publication_structural_candidate"
]


# Avoid duplicate field names.
final_fields = []

for field in (
    step05_fields +
    new_fields
):

    if field not in final_fields:

        final_fields.append(
            field
        )


master_file = os.path.join(
    OUTDIR,
    "09B1_702_strong_pair_multi_evidence_master.tsv"
)


write_tsv(
    master_file,
    master_rows,
    final_fields
)


# ============================================================
# 18. Evidence-tier summary
# ============================================================

tier_counter = Counter(
    row["Evidence_tier"]
    for row in master_rows
)


tier_order = [
    "Tier4",
    "Tier3",
    "Tier2",
    "Tier1",
    "Tier0"
]


tier_summary = []

for tier in tier_order:

    rows = [
        row
        for row in master_rows
        if row[
            "Evidence_tier"
        ] == tier
    ]

    tier_summary.append({

        "Evidence_tier":
            tier,

        "Pair_N":
            len(rows),

        "Percent_of_702":
            len(rows) /
            702 * 100,

        "Unique_gene_N":
            len({
                row["gene_id"]
                for row in rows
            }),

        "Definition":
            (
                rows[0][
                    "Evidence_tier_definition"
                ]
                if rows
                else ""
            )
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "09B2_evidence_tier_summary.tsv"
    ),
    tier_summary,
    [
        "Evidence_tier",
        "Pair_N",
        "Percent_of_702",
        "Unique_gene_N",
        "Definition"
    ]
)


# ============================================================
# 19. Tissue-level evidence summary
# ============================================================

tissues = sorted({
    row["Tissue"]
    for row in master_rows
})


tissue_summary = []


for tissue in tissues:

    rows = [
        row
        for row in master_rows
        if row[
            "Tissue"
        ] == tissue
    ]


    tier_counts = Counter(
        row["Evidence_tier"]
        for row in rows
    )


    reported_positive = sum(
        int(
            row[
                "Reported_motif_N"
            ]
        ) > 0
        for row in rows
    )


    q005_positive = sum(
        int(
            row[
                "Q005_motif_N"
            ]
        ) > 0
        for row in rows
    )


    primary_q005 = sum(
        int(
            row[
                "Primary_q005_motif_N"
            ]
        ) > 0
        for row in rows
    )


    promoter_pair_n = sum(
        row[
            "Pair_is_Promoter_TSSproximal"
        ] == "True"
        for row in rows
    )


    tissue_summary.append({

        "Tissue":
            tissue,

        "Strong_pair_N":
            len(rows),

        "Strong_unique_gene_N":
            len({
                row["gene_id"]
                for row in rows
            }),

        "Motif_analysis_status":
            rows[0][
                "Motif_analysis_status"
            ],

        "Reported_motif_positive_pair_N":
            reported_positive,

        "Q005_motif_positive_pair_N":
            q005_positive,

        "Primary_q005_pair_N":
            primary_q005,

        "Tier4_pair_N":
            tier_counts["Tier4"],

        "Tier3_pair_N":
            tier_counts["Tier3"],

        "Tier2_pair_N":
            tier_counts["Tier2"],

        "Tier1_pair_N":
            tier_counts["Tier1"],

        "Tier0_pair_N":
            tier_counts["Tier0"],

        "Pair_specific_promoter_N":
            promoter_pair_n,

        "Pair_specific_promoter_fraction":
            promoter_pair_n /
            len(rows),

        "Step09A_GO_BP_FDR005_N":
            rows[0].get(
                "Step09A_tissue_GO_BP_FDR005_N",
                ""
            ),

        "Step09A_KEGG_FDR005_N":
            rows[0].get(
                "Step09A_tissue_KEGG_FDR005_N",
                ""
            )
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "09B3_tissue_evidence_summary.tsv"
    ),
    tissue_summary,
    [
        "Tissue",
        "Strong_pair_N",
        "Strong_unique_gene_N",
        "Motif_analysis_status",
        "Reported_motif_positive_pair_N",
        "Q005_motif_positive_pair_N",
        "Primary_q005_pair_N",
        "Tier4_pair_N",
        "Tier3_pair_N",
        "Tier2_pair_N",
        "Tier1_pair_N",
        "Tier0_pair_N",
        "Pair_specific_promoter_N",
        "Pair_specific_promoter_fraction",
        "Step09A_GO_BP_FDR005_N",
        "Step09A_KEGG_FDR005_N"
    ]
)


# ============================================================
# 20. Used motif-family dictionary
# ============================================================

used_motif_rows = []


for tissue, motif_id in sorted(
    used_motifs
):

    d = motif_dictionary[
        (
            tissue,
            motif_id
        )
    ]


    r_arch = get_motif_info(
        tissue,
        motif_id,
        "Reported_p1e-4"
    )


    q_arch = get_motif_info(
        tissue,
        motif_id,
        "High_confidence_q005"
    )


    n_reported = sum(
        1
        for row in reported_links
        if (
            row["Tissue"] == tissue
            and
            row["Motif_ID"] == motif_id
        )
    )


    n_q005 = sum(
        1
        for row in q005_links
        if (
            row["Tissue"] == tissue
            and
            row["Motif_ID"] == motif_id
        )
    )


    out = dict(
        d
    )


    out[
        "Reported_motif_peak_link_N"
    ] = n_reported


    out[
        "Q005_motif_peak_link_N"
    ] = n_q005


    out[
        "Reported_promoter_enrichment"
    ] = str(
        r_arch[
            "Promoter_support"
        ]
    )


    out[
        "Reported_promoter_FDR"
    ] = (
        ""
        if r_arch[
            "Promoter_FDR"
        ] is None
        else
        r_arch[
            "Promoter_FDR"
        ]
    )


    out[
        "Reported_TSS_proximity"
    ] = str(
        r_arch[
            "TSS_support"
        ]
    )


    out[
        "Reported_TSS_FDR"
    ] = (
        ""
        if r_arch[
            "TSS_FDR"
        ] is None
        else
        r_arch[
            "TSS_FDR"
        ]
    )


    out[
        "Reported_distal_depletion"
    ] = str(
        r_arch[
            "Distal_depletion"
        ]
    )


    out[
        "Q005_promoter_enrichment"
    ] = str(
        q_arch[
            "Promoter_support"
        ]
    )


    out[
        "Q005_promoter_FDR"
    ] = (
        ""
        if q_arch[
            "Promoter_FDR"
        ] is None
        else
        q_arch[
            "Promoter_FDR"
        ]
    )


    out[
        "Q005_TSS_proximity"
    ] = str(
        q_arch[
            "TSS_support"
        ]
    )


    out[
        "Q005_TSS_FDR"
    ] = (
        ""
        if q_arch[
            "TSS_FDR"
        ] is None
        else
        q_arch[
            "TSS_FDR"
        ]
    )


    out[
        "Q005_distal_depletion"
    ] = str(
        q_arch[
            "Distal_depletion"
        ]
    )


    used_motif_rows.append(
        out
    )


motif_dictionary_fields = list(
    d1[0].keys()
) + [
    "Reported_motif_peak_link_N",
    "Q005_motif_peak_link_N",
    "Reported_promoter_enrichment",
    "Reported_promoter_FDR",
    "Reported_TSS_proximity",
    "Reported_TSS_FDR",
    "Reported_distal_depletion",
    "Q005_promoter_enrichment",
    "Q005_promoter_FDR",
    "Q005_TSS_proximity",
    "Q005_TSS_FDR",
    "Q005_distal_depletion"
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09B4_used_motif_family_dictionary.tsv"
    ),
    used_motif_rows,
    motif_dictionary_fields
)


# ============================================================
# 21. Tier3+ publication candidate table
# ============================================================

tier3plus = [
    row
    for row in master_rows
    if int(
        row[
            "Evidence_tier_numeric"
        ]
    ) >= 3
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09B5_Tier3_Tier4_publication_candidate_pairs.tsv"
    ),
    tier3plus,
    final_fields
)


tier4_only = [
    row
    for row in master_rows
    if row[
        "Evidence_tier"
    ] == "Tier4"
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09B6_Tier4_structurally_supported_candidate_pairs.tsv"
    ),
    tier4_only,
    final_fields
)


# ============================================================
# 22. QC summary
# ============================================================

expected_motif_scope_strong = sum(
    1
    for row in strong
    if row["Tissue"] in motif_scope_tissues
)


observed_motif_scope_strong = sum(
    1
    for row in master_rows
    if row[
        "Tissue"
    ] in motif_scope_tissues
)


qc_rows = [

    {
        "Metric":
            "Concordant_HC_input_pairs",
        "Value":
            len(step05),
        "Expected":
            4674,
        "Status":
            "PASS"
            if len(step05) == 4674
            else "FAIL"
    },

    {
        "Metric":
            "Official_strong_pairs",
        "Value":
            len(strong),
        "Expected":
            702,
        "Status":
            "PASS"
            if len(strong) == 702
            else "FAIL"
    },

    {
        "Metric":
            "Final_master_rows",
        "Value":
            len(master_rows),
        "Expected":
            702,
        "Status":
            "PASS"
            if len(master_rows) == 702
            else "FAIL"
    },

    {
        "Metric":
            "Final_unique_SAF_peak_ids",
        "Value":
            len({
                row["SAF_peak_id"]
                for row in master_rows
            }),
        "Expected":
            702,
        "Status":
            "PASS"
            if len({
                row["SAF_peak_id"]
                for row in master_rows
            }) == 702
            else "FAIL"
    },

    {
        "Metric":
            "Step07_missing_strong_pairs",
        "Value":
            step07_missing,
        "Expected":
            0,
        "Status":
            "PASS"
            if step07_missing == 0
            else "FAIL"
    },

    {
        "Metric":
            "Reported_motif_links",
        "Value":
            len(reported_links),
        "Expected":
            902,
        "Status":
            "PASS"
            if len(reported_links) == 902
            else "FAIL"
    },

    {
        "Metric":
            "Q005_motif_links",
        "Value":
            len(q005_links),
        "Expected":
            612,
        "Status":
            "PASS"
            if len(q005_links) == 612
            else "FAIL"
    },

    {
        "Metric":
            "Reported_unique_motif_positive_peaks",
        "Value":
            observed_reported_unique,
        "Expected":
            expected_reported_unique,
        "Status":
            "PASS"
            if (
                observed_reported_unique ==
                expected_reported_unique
            )
            else "FAIL"
    },

    {
        "Metric":
            "Q005_unique_motif_positive_peaks",
        "Value":
            observed_q005_unique,
        "Expected":
            expected_q005_unique,
        "Status":
            "PASS"
            if (
                observed_q005_unique ==
                expected_q005_unique
            )
            else "FAIL"
    },

    {
        "Metric":
            "Strong_pairs_in_motif_analysis_tissues",
        "Value":
            observed_motif_scope_strong,
        "Expected":
            expected_motif_scope_strong,
        "Status":
            "PASS"
            if (
                observed_motif_scope_strong ==
                expected_motif_scope_strong
            )
            else "FAIL"
    },

    {
        "Metric":
            "Used_motifs_missing_Step08D_metadata",
        "Value":
            len(
                missing_motif_dictionary
            ),
        "Expected":
            0,
        "Status":
            "PASS"
            if len(
                missing_motif_dictionary
            ) == 0
            else "FAIL"
    },

    {
        "Metric":
            "Tier3plus_pairs",
        "Value":
            len(
                tier3plus
            ),
        "Expected":
            "INFO",
        "Status":
            "INFO"
    },

    {
        "Metric":
            "Tier4_pairs",
        "Value":
            len(
                tier4_only
            ),
        "Expected":
            "INFO",
        "Status":
            "INFO"
    }
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09B7_QC_summary.tsv"
    ),
    qc_rows,
    [
        "Metric",
        "Value",
        "Expected",
        "Status"
    ]
)


# ============================================================
# 23. Column dictionary
# ============================================================

column_dictionary = [

    {
        "Column":
            "Evidence_tier",
        "Meaning":
            "Integrated evidence category Tier0-Tier4; not a statistical significance level."
    },

    {
        "Column":
            "Tier0_subtype",
        "Meaning":
            "Distinguishes motif-evaluated pairs lacking selected sites from tissues not evaluated in Step08 motif analysis."
    },

    {
        "Column":
            "Step07_Pair_specific_class",
        "Meaning":
            "Promoter-first genomic class relative specifically to the UROPA-associated gene."
    },

    {
        "Column":
            "Step07_Peak_to_UROPA_TSS_min_abs_distance_bp",
        "Meaning":
            "Minimum absolute distance from the peak to the representative-transcript TSS of the UROPA-associated gene."
    },

    {
        "Column":
            "Reported_motif_*",
        "Meaning":
            "Selected motif sites reported by FIMO at p<=1e-4."
    },

    {
        "Column":
            "Q005_motif_*",
        "Meaning":
            "Selected motif sites passing FIMO q<=0.05 within the corresponding FIMO run."
    },

    {
        "Column":
            "Primary_q005_*",
        "Meaning":
            "q<=0.05 motif evidence from Step08D primary-prioritized motif families."
    },

    {
        "Column":
            "Candidate_TFs",
        "Meaning":
            "TF candidates from known motif identity or Tomtom family-level annotation; not validated TF occupancy."
    },

    {
        "Column":
            "High_RNA_candidate_TFs",
        "Meaning":
            "Candidate TFs supported by RNA expression in the corresponding tissue."
    },

    {
        "Column":
            "Tissue_preferred_candidate_TFs",
        "Meaning":
            "Candidate TFs with target-tissue expression preference in Step08C."
    },

    {
        "Column":
            "Any_*_family_promoter_enrichment",
        "Meaning":
            "Family-level Step08E3A promoter enrichment FDR<0.05; this is not a pair-specific statistical test."
    },

    {
        "Column":
            "Any_*_family_TSS_proximity",
        "Meaning":
            "Family-level Step08E3A significantly shorter TSS distance FDR<0.05."
    },

    {
        "Column":
            "Primary_q005_has_structural_support",
        "Meaning":
            "At least one primary q005 motif family shows promoter enrichment or significant TSS proximity."
    },

    {
        "Column":
            "Any_reported_AME_formal_E_lt_0.05",
        "Meaning":
            "Whether any reported known motif had formal AME database-wide E<0.05."
    },

    {
        "Column":
            "Any_reported_STREME_holdout_formal_E_lt_0.05",
        "Meaning":
            "Whether any reported de novo motif passed independent STREME holdout E<0.05."
    },

    {
        "Column":
            "Step09A_tissue_*",
        "Meaning":
            "Tissue-level matched-background GO/KEGG context from Step09A; not pair-specific functional evidence."
    },

    {
        "Column":
            "Publication_core_candidate",
        "Meaning":
            "TRUE for Tier3 or Tier4 pairs; convenience flag for downstream figure/network selection."
    },

    {
        "Column":
            "Publication_structural_candidate",
        "Meaning":
            "TRUE for Tier4 pairs; convenience flag only."
    }
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09B8_column_dictionary.tsv"
    ),
    column_dictionary,
    [
        "Column",
        "Meaning"
    ]
)


# ============================================================
# 24. Console output
# ============================================================

print("\n" + "=" * 80)
print("STEP 09B COMPLETED")
print("=" * 80)

print(
    f"Final master rows                  : {len(master_rows)}"
)

print(
    f"Final unique SAF peak IDs          : "
    f"{len(set(x['SAF_peak_id'] for x in master_rows))}"
)

print(
    f"Unique associated genes            : "
    f"{len(set(x['gene_id'] for x in master_rows))}"
)

print(
    f"Reported motif-positive pairs      : "
    f"{observed_reported_unique}"
)

print(
    f"q<=0.05 motif-positive pairs       : "
    f"{observed_q005_unique}"
)

print("\nEvidence tiers:")

for row in tier_summary:

    print(
        f"{row['Evidence_tier']:6s} "
        f"{row['Pair_N']:4d} pairs | "
        f"{row['Unique_gene_N']:4d} genes | "
        f"{row['Percent_of_702']:.2f}%"
    )


print("\nTissue evidence summary:")

for row in tissue_summary:

    print(
        f"{row['Tissue']:12s} "
        f"strong={row['Strong_pair_N']:3d} "
        f"reported={row['Reported_motif_positive_pair_N']:3d} "
        f"q005={row['Q005_motif_positive_pair_N']:3d} "
        f"primary_q005={row['Primary_q005_pair_N']:3d} "
        f"Tier4={row['Tier4_pair_N']:3d}"
    )


print("\nInterpretation boundaries:")

print(
    "1. One row is one strong UROPA-associated peak-gene pair."
)

print(
    "2. Associated genes are candidate genomic associations, "
    "not validated regulatory targets."
)

print(
    "3. Candidate TFs do not establish direct TF occupancy."
)

print(
    "4. Evidence tiers are integration labels, "
    "not significance levels."
)

print(
    "5. Tier0 in non-Step08 tissues means motif not evaluated, "
    "not motif absent."
)

print(
    "6. Family-level structural support is inherited from "
    "Step08E3A and is not a new pair-specific test."
)


print("\nMain outputs:")

for name in [

    "09B1_702_strong_pair_multi_evidence_master.tsv",

    "09B2_evidence_tier_summary.tsv",

    "09B3_tissue_evidence_summary.tsv",

    "09B4_used_motif_family_dictionary.tsv",

    "09B5_Tier3_Tier4_publication_candidate_pairs.tsv",

    "09B6_Tier4_structurally_supported_candidate_pairs.tsv",

    "09B7_QC_summary.tsv",

    "09B8_column_dictionary.tsv"
]:

    print(
        os.path.join(
            OUTDIR,
            name
        )
    )


print("=" * 80)

