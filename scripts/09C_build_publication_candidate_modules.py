#!/usr/bin/env python3

# ============================================================
# STEP 09C
#
# Publication-ready candidate regulatory modules
# and network-ready node/edge tables
#
# SOFTWARE:
#   Python >= 3
#
# No third-party Python package is required.
# No pip/conda installation is required.
#
# ============================================================
#
# SOURCE OF TRUTH:
#
#   09B1_702_strong_pair_multi_evidence_master.tsv
#
# One source row =
#   one strong UROPA-associated peak-gene pair.
#
# ============================================================
#
# MODULE ROLE DEFINITIONS
#
# Primary_structural:
#   Step08D primary motif
#   + FIMO q <= 0.05 sites
#   + family-level promoter enrichment and/or
#     significant TSS proximity in Step08E3A.
#
# Primary_integrated:
#   Step08D primary motif
#   + FIMO q <= 0.05 sites
#   + no significant family-level structural shift.
#
# Supporting_q005_structural:
#   Step08D secondary motif
#   + FIMO q <= 0.05
#   + family-level structural support.
#
# Supporting_q005:
#   Step08D secondary motif
#   + FIMO q <= 0.05
#   + no significant family-level structural support.
#
# Exploratory_reported:
#   no q <= 0.05 selected sites;
#   use reported FIMO p <= 1e-4 sites only.
#
# IMPORTANT:
# These roles are integrated evidence labels,
# NOT statistical significance levels.
#
# ============================================================
#
# REPRESENTATIVE GENE RANKING
#
# There is NO arbitrary weighted score.
#
# Each motif-associated gene is first deduplicated.
# If one gene has multiple motif-positive strong peaks,
# one representative pair is selected deterministically.
#
# For structurally supported modules:
#
#   1. named gene required for publication shortlist
#   2. pair-specific promoter first
#   3. absolute TSS distance <= 2 kb
#   4. higher minimum(P348 rho, P350 rho)
#   5. higher rho_mean
#   6. smaller |rho_P348-rho_P350|
#   7. more motif-positive peaks for the gene
#   8. smaller TSS distance
#
# For non-structural modules:
#
#   1. named gene required for publication shortlist
#   2. higher minimum(P348 rho, P350 rho)
#   3. higher rho_mean
#   4. smaller |rho_P348-rho_P350|
#   5. more motif-positive peaks for the gene
#   6. promoter status only as a late tie-breaker
#   7. smaller TSS distance
#
# This prevents Spleen modules from being artificially
# forced into a promoter-centered interpretation.
#
# ============================================================
#
# NETWORK SEMANTICS
#
# candidate TF -> motif:
#   candidate annotation only;
#   NOT demonstrated TF occupancy.
#
# motif -> peak:
#   FIMO motif-site localization.
#
# peak -> associated gene:
#   strong RNA-ATAC coupled UROPA genomic association;
#   NOT validated causal regulation.
#
# ============================================================


import os
import csv
import re
import math
from collections import defaultdict, Counter


# ============================================================
# 0. Files
# ============================================================

MASTER_FILE = (
    "09B_702_pair_multi_evidence_master/"
    "09B1_702_strong_pair_multi_evidence_master.tsv"
)

MOTIF_DICT_FILE = (
    "09B_702_pair_multi_evidence_master/"
    "09B4_used_motif_family_dictionary.tsv"
)

OUTDIR = "09C_publication_candidate_modules"

MODULE_DIR = os.path.join(
    OUTDIR,
    "09C_modules"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)

os.makedirs(
    MODULE_DIR,
    exist_ok=True
)


# ============================================================
# 1. Helpers
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
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def write_tsv(
    path,
    rows,
    fields
):

    with open(
        path,
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:

            writer.writerow({
                field: row.get(
                    field,
                    ""
                )
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

    x = as_float(x)

    if x is None:
        return 0

    return int(x)


def as_bool(x):

    if x is None:
        return False

    return str(x).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y"
    }


def split_semicolon(x):

    if x is None:
        return []

    x = str(x).strip()

    if x == "":
        return []

    return [
        z.strip()
        for z in x.split(";")
        if z.strip()
    ]


def unique_sorted(values):

    return sorted({
        str(x).strip()
        for x in values
        if x is not None
        and str(x).strip() != ""
    })


def join_values(values):

    return ";".join(
        unique_sorted(values)
    )


def safe_name(x):

    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(x)
    )


def finite_or_inf(x):

    x = as_float(x)

    if x is None:
        return math.inf

    return x


def finite_or_minus_inf(x):

    x = as_float(x)

    if x is None:
        return -math.inf

    return x


def gene_has_name(row):

    return (
        str(
            row.get(
                "gene_name",
                ""
            )
        ).strip()
        != ""
    )


def get_min_animal_rho(row):

    x = as_float(
        row.get(
            "rho_P348"
        )
    )

    y = as_float(
        row.get(
            "rho_P350"
        )
    )

    if x is None or y is None:
        return -math.inf

    return min(
        x,
        y
    )


def get_animal_difference(row):

    x = as_float(
        row.get(
            "rho_P348"
        )
    )

    y = as_float(
        row.get(
            "rho_P350"
        )
    )

    if x is None or y is None:
        return math.inf

    return abs(
        x - y
    )


def exact_motif_present(
    field_value,
    motif_id
):

    return (
        motif_id
        in split_semicolon(
            field_value
        )
    )


# ============================================================
# 2. Read data
# ============================================================

print("=" * 80)
print("STEP 09C")
print("Publication-ready candidate regulatory modules")
print("=" * 80)

print(
    "No third-party Python package is required."
)

print("\nReading Step09B 702-pair master...")

master = read_tsv(
    MASTER_FILE
)

if len(master) != 702:

    raise RuntimeError(
        "Expected 702 rows in Step09B master, "
        f"observed {len(master)}."
    )


if len({
    row["SAF_peak_id"]
    for row in master
}) != 702:

    raise RuntimeError(
        "Step09B master SAF_peak_id is not unique."
    )


print(
    "Master rows:",
    len(master)
)

print(
    "Unique associated genes:",
    len({
        row["gene_id"]
        for row in master
    })
)


print("\nReading Step09B motif-family dictionary...")

motifs = read_tsv(
    MOTIF_DICT_FILE
)

if len(motifs) != 11:

    raise RuntimeError(
        "Expected 11 used motif families, "
        f"observed {len(motifs)}."
    )


print(
    "Selected motif families:",
    len(motifs)
)


# ============================================================
# 3. Define module roles automatically
# ============================================================

module_rows = []


for row in motifs:

    tissue = row["Tissue"]

    motif_id = row[
        "Representative_motif_ID"
    ]

    fimo_tier = str(
        row.get(
            "FIMO_tier",
            ""
        )
    ).strip().lower()

    reported_n = as_int(
        row.get(
            "Reported_motif_peak_link_N"
        )
    )

    q005_n = as_int(
        row.get(
            "Q005_motif_peak_link_N"
        )
    )

    q005_promoter = as_bool(
        row.get(
            "Q005_promoter_enrichment"
        )
    )

    q005_tss = as_bool(
        row.get(
            "Q005_TSS_proximity"
        )
    )

    q005_structural = (
        q005_promoter
        or
        q005_tss
    )


    if (
        q005_n > 0
        and
        fimo_tier == "primary"
        and
        q005_structural
    ):

        module_role = (
            "Primary_structural"
        )

        threshold = (
            "Q005"
        )

        threshold_definition = (
            "FIMO q<=0.05"
        )

        role_priority = 1


    elif (
        q005_n > 0
        and
        fimo_tier == "primary"
    ):

        module_role = (
            "Primary_integrated"
        )

        threshold = (
            "Q005"
        )

        threshold_definition = (
            "FIMO q<=0.05"
        )

        role_priority = 2


    elif (
        q005_n > 0
        and
        q005_structural
    ):

        module_role = (
            "Supporting_q005_structural"
        )

        threshold = (
            "Q005"
        )

        threshold_definition = (
            "FIMO q<=0.05"
        )

        role_priority = 3


    elif q005_n > 0:

        module_role = (
            "Supporting_q005"
        )

        threshold = (
            "Q005"
        )

        threshold_definition = (
            "FIMO q<=0.05"
        )

        role_priority = 4


    else:

        module_role = (
            "Exploratory_reported"
        )

        threshold = (
            "Reported"
        )

        threshold_definition = (
            "FIMO p<=1e-4"
        )

        role_priority = 5


    candidate_tfs = split_semicolon(
        row.get(
            "Candidate_TFs",
            ""
        )
    )

    high_rna_tfs = split_semicolon(
        row.get(
            "High_RNA_candidate_TFs",
            ""
        )
    )

    preferred_tfs = split_semicolon(
        row.get(
            "Tissue_preferred_candidate_TFs",
            ""
        )
    )


    if preferred_tfs:

        representative_tfs = (
            preferred_tfs
        )

        tf_selection_basis = (
            "Tissue_preferred_candidate_TFs"
        )


    elif high_rna_tfs:

        representative_tfs = (
            high_rna_tfs
        )

        tf_selection_basis = (
            "High_RNA_candidate_TFs"
        )


    else:

        representative_tfs = (
            candidate_tfs
        )

        tf_selection_basis = (
            "Motif_annotation_candidates"
        )


    family_id = row.get(
        "Motif_family_ID",
        ""
    )


    if ".DENOVO." in family_id.upper():

        tf_assignment_type = (
            "Tomtom_family_level_candidate"
        )

    else:

        tf_assignment_type = (
            "Known_motif_identity_candidate"
        )


    module_id = (
        f"{tissue}__{motif_id}"
    )


    module = dict(
        row
    )

    module.update({

        "Module_ID":
            module_id,

        "Module_role":
            module_role,

        "Module_role_priority":
            role_priority,

        "Candidate_threshold":
            threshold,

        "Candidate_threshold_definition":
            threshold_definition,

        "Q005_structural_support":
            str(
                q005_structural
            ),

        "Representative_candidate_TFs":
            ";".join(
                representative_tfs
            ),

        "TF_selection_basis":
            tf_selection_basis,

        "TF_assignment_type":
            tf_assignment_type
    })


    module_rows.append(
        module
    )


module_rows.sort(
    key=lambda x: (
        as_int(
            x[
                "Module_role_priority"
            ]
        ),
        x["Tissue"],
        x[
            "Representative_motif_ID"
        ]
    )
)


module_fields = list(
    motifs[0].keys()
)

for field in [
    "Module_ID",
    "Module_role",
    "Module_role_priority",
    "Candidate_threshold",
    "Candidate_threshold_definition",
    "Q005_structural_support",
    "Representative_candidate_TFs",
    "TF_selection_basis",
    "TF_assignment_type"
]:

    if field not in module_fields:
        module_fields.append(
            field
        )


write_tsv(
    os.path.join(
        OUTDIR,
        "09C1_module_definition.tsv"
    ),
    module_rows,
    module_fields
)


# ============================================================
# 4. Extract motif-specific candidate pairs
# ============================================================

all_module_pairs = []

module_pair_map = defaultdict(
    list
)


for module in module_rows:

    tissue = module[
        "Tissue"
    ]

    motif_id = module[
        "Representative_motif_ID"
    ]

    module_id = module[
        "Module_ID"
    ]

    threshold = module[
        "Candidate_threshold"
    ]


    if threshold == "Q005":

        motif_field = (
            "Q005_motif_IDs"
        )

    else:

        motif_field = (
            "Reported_motif_IDs"
        )


    selected = []

    for row in master:

        if row["Tissue"] != tissue:
            continue

        if exact_motif_present(
            row.get(
                motif_field,
                ""
            ),
            motif_id
        ):

            selected.append(
                row
            )


    expected_n = (
        as_int(
            module.get(
                "Q005_motif_peak_link_N"
            )
        )
        if threshold == "Q005"
        else
        as_int(
            module.get(
                "Reported_motif_peak_link_N"
            )
        )
    )


    if len(selected) != expected_n:

        raise RuntimeError(
            f"{module_id}: expected {expected_n} "
            f"candidate peak-gene pairs, "
            f"observed {len(selected)}."
        )


    for row in selected:

        out = dict(
            row
        )

        out.update({

            "Module_ID":
                module_id,

            "Module_role":
                module[
                    "Module_role"
                ],

            "Module_role_priority":
                module[
                    "Module_role_priority"
                ],

            "Module_motif_ID":
                motif_id,

            "Module_motif_name":
                module.get(
                    "Representative_motif_name",
                    ""
                ),

            "Module_motif_family_label":
                module.get(
                    "Motif_family_label",
                    ""
                ),

            "Module_candidate_threshold":
                threshold,

            "Module_candidate_threshold_definition":
                module[
                    "Candidate_threshold_definition"
                ],

            "Module_structural_support":
                module[
                    "Q005_structural_support"
                ],

            "Module_representative_candidate_TFs":
                module[
                    "Representative_candidate_TFs"
                ],

            "Module_TF_selection_basis":
                module[
                    "TF_selection_basis"
                ],

            "Module_TF_assignment_type":
                module[
                    "TF_assignment_type"
                ]
        })


        all_module_pairs.append(
            out
        )

        module_pair_map[
            module_id
        ].append(
            out
        )


# ============================================================
# 5. Module-specific gene multiplicity
# ============================================================

for module in module_rows:

    module_id = module[
        "Module_ID"
    ]

    pairs = module_pair_map[
        module_id
    ]

    gene_counts = Counter(
        row[
            "gene_id"
        ]
        for row in pairs
    )

    for row in pairs:

        row[
            "Module_motif_positive_peak_N_for_gene"
        ] = gene_counts[
            row[
                "gene_id"
            ]
        ]


pair_fields = list(
    master[0].keys()
)

extra_pair_fields = [
    "Module_ID",
    "Module_role",
    "Module_role_priority",
    "Module_motif_ID",
    "Module_motif_name",
    "Module_motif_family_label",
    "Module_candidate_threshold",
    "Module_candidate_threshold_definition",
    "Module_structural_support",
    "Module_representative_candidate_TFs",
    "Module_TF_selection_basis",
    "Module_TF_assignment_type",
    "Module_motif_positive_peak_N_for_gene"
]

for field in extra_pair_fields:

    if field not in pair_fields:
        pair_fields.append(
            field
        )


write_tsv(
    os.path.join(
        OUTDIR,
        "09C2_all_module_peak_gene_candidates.tsv"
    ),
    all_module_pairs,
    pair_fields
)


# ============================================================
# 6. Representative ranking functions
# ============================================================

def pair_sort_key(
    row,
    structural_module
):

    named = gene_has_name(
        row
    )

    promoter = as_bool(
        row.get(
            "Pair_is_Promoter_TSSproximal"
        )
    )

    tss2kb = as_bool(
        row.get(
            "Pair_abs_TSS_distance_le_2kb"
        )
    )

    min_rho = get_min_animal_rho(
        row
    )

    rho_mean = finite_or_minus_inf(
        row.get(
            "rho_mean"
        )
    )

    animal_diff = get_animal_difference(
        row
    )

    peak_n = as_int(
        row.get(
            "Module_motif_positive_peak_N_for_gene"
        )
    )

    tss_dist = finite_or_inf(
        row.get(
            "Step07_Peak_to_UROPA_TSS_min_abs_distance_bp"
        )
    )


    if structural_module:

        return (

            0 if named else 1,

            0 if promoter else 1,

            0 if tss2kb else 1,

            -min_rho,

            -rho_mean,

            animal_diff,

            -peak_n,

            tss_dist,

            row.get(
                "gene_id",
                ""
            ),

            row.get(
                "SAF_peak_id",
                ""
            )
        )


    return (

        0 if named else 1,

        -min_rho,

        -rho_mean,

        animal_diff,

        -peak_n,

        0 if promoter else 1,

        tss_dist,

        row.get(
            "gene_id",
            ""
        ),

        row.get(
            "SAF_peak_id",
            ""
        )
    )


# ============================================================
# 7. Deduplicate by gene and rank representatives
# ============================================================

gene_ranking_rows = []

top10_rows = []

figure_top5_rows = []


for module in module_rows:

    module_id = module[
        "Module_ID"
    ]

    role = module[
        "Module_role"
    ]

    structural_module = as_bool(
        module[
            "Q005_structural_support"
        ]
    )

    pairs = module_pair_map[
        module_id
    ]

    by_gene = defaultdict(
        list
    )


    for row in pairs:

        by_gene[
            row["gene_id"]
        ].append(
            row
        )


    best_gene_rows = []


    for gid, gene_pairs in by_gene.items():

        gene_pairs_sorted = sorted(
            gene_pairs,
            key=lambda x: pair_sort_key(
                x,
                structural_module
            )
        )

        best = dict(
            gene_pairs_sorted[0]
        )


        best[
            "Representative_eligible_named_gene"
        ] = str(
            gene_has_name(
                best
            )
        )


        best[
            "Representative_min_animal_rho"
        ] = get_min_animal_rho(
            best
        )


        best[
            "Representative_animal_rho_abs_difference"
        ] = get_animal_difference(
            best
        )


        best[
            "Representative_pair_selection_rule"
        ] = (
            "Structural_module_hierarchical_rule"
            if structural_module
            else
            "Nonstructural_module_hierarchical_rule"
        )


        best_gene_rows.append(
            best
        )


    best_gene_rows.sort(
        key=lambda x: pair_sort_key(
            x,
            structural_module
        )
    )


    rank = 0

    publication_rank = 0


    for row in best_gene_rows:

        rank += 1

        row[
            "Gene_rank_all_within_module"
        ] = rank


        if gene_has_name(
            row
        ):

            publication_rank += 1

            row[
                "Publication_named_gene_rank_within_module"
            ] = publication_rank

        else:

            row[
                "Publication_named_gene_rank_within_module"
            ] = ""


        gene_ranking_rows.append(
            row
        )


        if (
            gene_has_name(
                row
            )
            and
            publication_rank <= 10
        ):

            top10_rows.append(
                dict(
                    row
                )
            )


        if (
            role
            !=
            "Exploratory_reported"
            and
            gene_has_name(
                row
            )
            and
            publication_rank <= 5
        ):

            figure_top5_rows.append(
                dict(
                    row
                )
            )


ranking_extra_fields = [
    "Representative_eligible_named_gene",
    "Representative_min_animal_rho",
    "Representative_animal_rho_abs_difference",
    "Representative_pair_selection_rule",
    "Gene_rank_all_within_module",
    "Publication_named_gene_rank_within_module"
]

ranking_fields = list(
    pair_fields
)

for field in ranking_extra_fields:

    if field not in ranking_fields:
        ranking_fields.append(
            field
        )


write_tsv(
    os.path.join(
        OUTDIR,
        "09C3_gene_level_candidate_ranking.tsv"
    ),
    gene_ranking_rows,
    ranking_fields
)


write_tsv(
    os.path.join(
        OUTDIR,
        "09C4_top10_representative_genes_per_module.tsv"
    ),
    top10_rows,
    ranking_fields
)


write_tsv(
    os.path.join(
        OUTDIR,
        "09C5_top5_primary_supporting_figure_candidates.tsv"
    ),
    figure_top5_rows,
    ranking_fields
)


# ============================================================
# 8. Module summary
# ============================================================

module_summary = []


for module in module_rows:

    module_id = module[
        "Module_ID"
    ]

    pairs = module_pair_map[
        module_id
    ]

    genes = {
        row["gene_id"]
        for row in pairs
    }

    named_genes = {
        row["gene_id"]
        for row in pairs
        if gene_has_name(
            row
        )
    }

    promoter_n = sum(
        as_bool(
            row.get(
                "Pair_is_Promoter_TSSproximal"
            )
        )
        for row in pairs
    )

    tss2kb_n = sum(
        as_bool(
            row.get(
                "Pair_abs_TSS_distance_le_2kb"
            )
        )
        for row in pairs
    )

    tss_values = [
        as_float(
            row.get(
                "Step07_Peak_to_UROPA_TSS_min_abs_distance_bp"
            )
        )
        for row in pairs
    ]

    tss_values = [
        x
        for x in tss_values
        if x is not None
    ]

    rho_values = [
        as_float(
            row.get(
                "rho_mean"
            )
        )
        for row in pairs
    ]

    rho_values = [
        x
        for x in rho_values
        if x is not None
    ]


    def median(values):

        if not values:
            return ""

        z = sorted(
            values
        )

        n = len(z)

        if n % 2 == 1:
            return z[n // 2]

        return (
            z[n // 2 - 1]
            +
            z[n // 2]
        ) / 2


    module_summary.append({

        "Module_ID":
            module_id,

        "Tissue":
            module[
                "Tissue"
            ],

        "Motif_ID":
            module[
                "Representative_motif_ID"
            ],

        "Motif_name":
            module.get(
                "Representative_motif_name",
                ""
            ),

        "Motif_family_label":
            module.get(
                "Motif_family_label",
                ""
            ),

        "Module_role":
            module[
                "Module_role"
            ],

        "Candidate_threshold":
            module[
                "Candidate_threshold_definition"
            ],

        "Structural_support":
            module[
                "Q005_structural_support"
            ],

        "Representative_candidate_TFs":
            module[
                "Representative_candidate_TFs"
            ],

        "TF_selection_basis":
            module[
                "TF_selection_basis"
            ],

        "Candidate_pair_N":
            len(
                pairs
            ),

        "Unique_gene_N":
            len(
                genes
            ),

        "Named_gene_N":
            len(
                named_genes
            ),

        "Promoter_pair_N":
            promoter_n,

        "Promoter_pair_fraction":
            (
                promoter_n /
                len(pairs)
                if pairs
                else ""
            ),

        "TSS_le_2kb_pair_N":
            tss2kb_n,

        "TSS_le_2kb_pair_fraction":
            (
                tss2kb_n /
                len(pairs)
                if pairs
                else ""
            ),

        "Median_abs_TSS_distance_bp":
            median(
                tss_values
            ),

        "Median_rho_mean":
            median(
                rho_values
            )
    })


module_summary.sort(
    key=lambda x: (
        next(
            m[
                "Module_role_priority"
            ]
            for m in module_rows
            if m[
                "Module_ID"
            ] == x[
                "Module_ID"
            ]
        ),
        x[
            "Tissue"
        ],
        x[
            "Motif_ID"
        ]
    )
)


module_summary_fields = [
    "Module_ID",
    "Tissue",
    "Motif_ID",
    "Motif_name",
    "Motif_family_label",
    "Module_role",
    "Candidate_threshold",
    "Structural_support",
    "Representative_candidate_TFs",
    "TF_selection_basis",
    "Candidate_pair_N",
    "Unique_gene_N",
    "Named_gene_N",
    "Promoter_pair_N",
    "Promoter_pair_fraction",
    "TSS_le_2kb_pair_N",
    "TSS_le_2kb_pair_fraction",
    "Median_abs_TSS_distance_bp",
    "Median_rho_mean"
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09C6_module_candidate_summary.tsv"
    ),
    module_summary,
    module_summary_fields
)


# ============================================================
# 9. Within-tissue module overlap
# ============================================================

overlap_rows = []


for i in range(
    len(module_rows)
):

    a = module_rows[i]

    for j in range(
        i + 1,
        len(module_rows)
    ):

        b = module_rows[j]

        if (
            a["Tissue"]
            !=
            b["Tissue"]
        ):
            continue


        a_pairs = module_pair_map[
            a["Module_ID"]
        ]

        b_pairs = module_pair_map[
            b["Module_ID"]
        ]


        a_peaks = {
            x["SAF_peak_id"]
            for x in a_pairs
        }

        b_peaks = {
            x["SAF_peak_id"]
            for x in b_pairs
        }


        a_genes = {
            x["gene_id"]
            for x in a_pairs
        }

        b_genes = {
            x["gene_id"]
            for x in b_pairs
        }


        peak_intersection = (
            a_peaks
            &
            b_peaks
        )

        gene_intersection = (
            a_genes
            &
            b_genes
        )


        peak_union = (
            a_peaks
            |
            b_peaks
        )

        gene_union = (
            a_genes
            |
            b_genes
        )


        overlap_rows.append({

            "Tissue":
                a["Tissue"],

            "Module_A":
                a["Module_ID"],

            "Module_B":
                b["Module_ID"],

            "Module_A_peak_N":
                len(
                    a_peaks
                ),

            "Module_B_peak_N":
                len(
                    b_peaks
                ),

            "Shared_peak_N":
                len(
                    peak_intersection
                ),

            "Peak_Jaccard":
                (
                    len(
                        peak_intersection
                    ) /
                    len(
                        peak_union
                    )
                    if peak_union
                    else 0
                ),

            "Module_A_gene_N":
                len(
                    a_genes
                ),

            "Module_B_gene_N":
                len(
                    b_genes
                ),

            "Shared_gene_N":
                len(
                    gene_intersection
                ),

            "Gene_Jaccard":
                (
                    len(
                        gene_intersection
                    ) /
                    len(
                        gene_union
                    )
                    if gene_union
                    else 0
                )
        })


overlap_fields = [
    "Tissue",
    "Module_A",
    "Module_B",
    "Module_A_peak_N",
    "Module_B_peak_N",
    "Shared_peak_N",
    "Peak_Jaccard",
    "Module_A_gene_N",
    "Module_B_gene_N",
    "Shared_gene_N",
    "Gene_Jaccard"
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09C7_within_tissue_module_overlap.tsv"
    ),
    overlap_rows,
    overlap_fields
)


# ============================================================
# 10. Write per-module files
# ============================================================

for module in module_rows:

    module_id = module[
        "Module_ID"
    ]

    d = os.path.join(
        MODULE_DIR,
        safe_name(
            module_id
        )
    )

    os.makedirs(
        d,
        exist_ok=True
    )


    module_pairs = module_pair_map[
        module_id
    ]


    module_ranked = [
        x
        for x in gene_ranking_rows
        if x[
            "Module_ID"
        ] == module_id
    ]


    module_top10 = [
        x
        for x in top10_rows
        if x[
            "Module_ID"
        ] == module_id
    ]


    write_tsv(
        os.path.join(
            d,
            "all_peak_gene_candidates.tsv"
        ),
        module_pairs,
        pair_fields
    )


    write_tsv(
        os.path.join(
            d,
            "gene_level_ranking.tsv"
        ),
        module_ranked,
        ranking_fields
    )


    write_tsv(
        os.path.join(
            d,
            "top10_representative_genes.tsv"
        ),
        module_top10,
        ranking_fields
    )


# ============================================================
# 11. Network builder
# ============================================================

module_lookup = {
    x["Module_ID"]:
        x
    for x in module_rows
}


def build_network(
    candidate_pairs,
    prefix
):

    nodes = {}

    edges = []

    edge_seen = set()


    def add_node(
        node_id,
        node_type,
        tissue,
        label,
        **attrs
    ):

        if node_id in nodes:
            return

        row = {

            "node_id":
                node_id,

            "node_type":
                node_type,

            "Tissue":
                tissue,

            "label":
                label
        }

        row.update(
            attrs
        )

        nodes[
            node_id
        ] = row


    def add_edge(
        source,
        target,
        edge_type,
        tissue,
        module_id,
        **attrs
    ):

        key = (
            source,
            target,
            edge_type,
            module_id
        )

        if key in edge_seen:
            return

        edge_seen.add(
            key
        )

        edge_id = (
            f"E{len(edges)+1:06d}"
        )

        row = {

            "edge_id":
                edge_id,

            "source":
                source,

            "target":
                target,

            "edge_type":
                edge_type,

            "Tissue":
                tissue,

            "Module_ID":
                module_id
        }

        row.update(
            attrs
        )

        edges.append(
            row
        )


    # --------------------------------------------------------
    # Process each module separately
    # --------------------------------------------------------

    pairs_by_module = defaultdict(
        list
    )

    for row in candidate_pairs:

        pairs_by_module[
            row[
                "Module_ID"
            ]
        ].append(
            row
        )


    for module_id, pairs in pairs_by_module.items():

        module = module_lookup[
            module_id
        ]

        tissue = module[
            "Tissue"
        ]

        motif_id = module[
            "Representative_motif_ID"
        ]

        motif_label = (
            module.get(
                "Representative_motif_name",
                ""
            )
            or
            motif_id
        )


        motif_node = (
            f"Motif:{tissue}:{motif_id}"
        )


        add_node(

            motif_node,

            "Motif",

            tissue,

            motif_label,

            Motif_ID =
                motif_id,

            Motif_family_label =
                module.get(
                    "Motif_family_label",
                    ""
                ),

            Module_role =
                module[
                    "Module_role"
                ],

            Candidate_threshold =
                module[
                    "Candidate_threshold_definition"
                ],

            Structural_support =
                module[
                    "Q005_structural_support"
                ]
        )


        # ----------------------------------------------------
        # Candidate TF -> motif
        # ----------------------------------------------------

        tf_candidates = split_semicolon(
            module[
                "Representative_candidate_TFs"
            ]
        )


        for tf in tf_candidates:

            tf_node = (
                f"TF:{tissue}:{tf}"
            )


            add_node(

                tf_node,

                "Candidate_TF",

                tissue,

                tf,

                TF_selection_basis =
                    module[
                        "TF_selection_basis"
                    ],

                TF_assignment_type =
                    module[
                        "TF_assignment_type"
                    ]
            )


            add_edge(

                tf_node,

                motif_node,

                "candidate_TF_to_motif",

                tissue,

                module_id,

                Evidence =
                    module[
                        "TF_selection_basis"
                    ],

                Interpretation = (
                    "Candidate TF annotation; "
                    "does not establish TF occupancy."
                )
            )


        # ----------------------------------------------------
        # Motif -> peak -> gene
        # ----------------------------------------------------

        for row in pairs:

            peak = row[
                "SAF_peak_id"
            ]

            gid = row[
                "gene_id"
            ]

            gname = (
                row.get(
                    "gene_name",
                    ""
                )
                or
                gid
            )


            peak_node = (
                f"Peak:{tissue}:{peak}"
            )


            gene_node = (
                f"Gene:{tissue}:{gid}"
            )


            add_node(

                peak_node,

                "ATAC_peak",

                tissue,

                peak,

                SAF_peak_id =
                    peak,

                Pair_specific_class =
                    row.get(
                        "Step07_Pair_specific_class",
                        ""
                    ),

                Abs_TSS_distance_bp =
                    row.get(
                        "Step07_Peak_to_UROPA_TSS_min_abs_distance_bp",
                        ""
                    ),

                Evidence_tier =
                    row.get(
                        "Evidence_tier",
                        ""
                    )
            )


            add_node(

                gene_node,

                "Associated_gene",

                tissue,

                gname,

                gene_id =
                    gid,

                gene_name =
                    row.get(
                        "gene_name",
                        ""
                    ),

                Strong_peak_N_for_gene_within_tissue =
                    row.get(
                        "Strong_peak_N_for_gene_within_tissue",
                        ""
                    )
            )


            add_edge(

                motif_node,

                peak_node,

                "motif_to_peak",

                tissue,

                module_id,

                Motif_ID =
                    motif_id,

                FIMO_evidence =
                    module[
                        "Candidate_threshold_definition"
                    ],

                Interpretation = (
                    "FIMO motif-site localization "
                    "within a strong accessible region."
                )
            )


            add_edge(

                peak_node,

                gene_node,

                "peak_to_associated_gene",

                tissue,

                module_id,

                rho_P348 =
                    row.get(
                        "rho_P348",
                        ""
                    ),

                rho_P350 =
                    row.get(
                        "rho_P350",
                        ""
                    ),

                rho_mean =
                    row.get(
                        "rho_mean",
                        ""
                    ),

                Pair_specific_class =
                    row.get(
                        "Step07_Pair_specific_class",
                        ""
                    ),

                Abs_TSS_distance_bp =
                    row.get(
                        "Step07_Peak_to_UROPA_TSS_min_abs_distance_bp",
                        ""
                    ),

                Evidence_tier =
                    row.get(
                        "Evidence_tier",
                        ""
                    ),

                Interpretation = (
                    "Strong RNA-ATAC coupled UROPA-associated "
                    "peak-gene candidate; not validated causal regulation."
                )
            )


    node_rows = list(
        nodes.values()
    )


    node_fields = [
        "node_id",
        "node_type",
        "Tissue",
        "label",
        "Motif_ID",
        "Motif_family_label",
        "Module_role",
        "Candidate_threshold",
        "Structural_support",
        "TF_selection_basis",
        "TF_assignment_type",
        "SAF_peak_id",
        "Pair_specific_class",
        "Abs_TSS_distance_bp",
        "Evidence_tier",
        "gene_id",
        "gene_name",
        "Strong_peak_N_for_gene_within_tissue"
    ]


    edge_fields = [
        "edge_id",
        "source",
        "target",
        "edge_type",
        "Tissue",
        "Module_ID",
        "Evidence",
        "Motif_ID",
        "FIMO_evidence",
        "rho_P348",
        "rho_P350",
        "rho_mean",
        "Pair_specific_class",
        "Abs_TSS_distance_bp",
        "Evidence_tier",
        "Interpretation"
    ]


    write_tsv(
        os.path.join(
            OUTDIR,
            f"{prefix}_nodes.tsv"
        ),
        node_rows,
        node_fields
    )


    write_tsv(
        os.path.join(
            OUTDIR,
            f"{prefix}_edges.tsv"
        ),
        edges,
        edge_fields
    )


    return (
        node_rows,
        edges
    )


# ============================================================
# 12. Full module network
# ============================================================

full_nodes, full_edges = build_network(

    all_module_pairs,

    "09C8_full_candidate_network"
)


# ============================================================
# 13. Publication figure network
#
# Top 5 named genes from all non-exploratory modules.
# ============================================================

figure_gene_keys = {
    (
        row[
            "Module_ID"
        ],
        row[
            "gene_id"
        ]
    )
    for row in figure_top5_rows
}


figure_network_pairs = [

    row
    for row in all_module_pairs

    if (
        row[
            "Module_role"
        ]
        !=
        "Exploratory_reported"
        and
        (
            row[
                "Module_ID"
            ],
            row[
                "gene_id"
            ]
        )
        in
        figure_gene_keys
    )
]


figure_nodes, figure_edges = build_network(

    figure_network_pairs,

    "09C9_representative_figure_network"
)


# ============================================================
# 14. QC
# ============================================================

expected_module_pair_rows = sum(

    (
        as_int(
            row.get(
                "Q005_motif_peak_link_N"
            )
        )
        if row[
            "Candidate_threshold"
        ] == "Q005"

        else

        as_int(
            row.get(
                "Reported_motif_peak_link_N"
            )
        )
    )

    for row in module_rows
)


module_observed_pass = True


for module in module_rows:

    module_id = module[
        "Module_ID"
    ]

    observed = len(
        module_pair_map[
            module_id
        ]
    )

    expected = (
        as_int(
            module.get(
                "Q005_motif_peak_link_N"
            )
        )

        if module[
            "Candidate_threshold"
        ] == "Q005"

        else

        as_int(
            module.get(
                "Reported_motif_peak_link_N"
            )
        )
    )

    if observed != expected:

        module_observed_pass = False


role_counts = Counter(
    x[
        "Module_role"
    ]
    for x in module_rows
)


qc = [

    {
        "Metric":
            "Step09B_master_rows",

        "Value":
            len(
                master
            ),

        "Expected":
            702,

        "Status":
            "PASS"
            if len(master) == 702
            else "FAIL"
    },

    {
        "Metric":
            "Selected_motif_modules",

        "Value":
            len(
                module_rows
            ),

        "Expected":
            11,

        "Status":
            "PASS"
            if len(module_rows) == 11
            else "FAIL"
    },

    {
        "Metric":
            "Module_pair_rows",

        "Value":
            len(
                all_module_pairs
            ),

        "Expected":
            expected_module_pair_rows,

        "Status":
            "PASS"
            if (
                len(
                    all_module_pairs
                )
                ==
                expected_module_pair_rows
            )
            else "FAIL"
    },

    {
        "Metric":
            "Per_module_pair_reconstruction",

        "Value":
            "PASS"
            if module_observed_pass
            else "FAIL",

        "Expected":
            "PASS",

        "Status":
            "PASS"
            if module_observed_pass
            else "FAIL"
    },

    {
        "Metric":
            "Primary_structural_modules",

        "Value":
            role_counts[
                "Primary_structural"
            ],

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Primary_integrated_modules",

        "Value":
            role_counts[
                "Primary_integrated"
            ],

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Supporting_q005_structural_modules",

        "Value":
            role_counts[
                "Supporting_q005_structural"
            ],

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Supporting_q005_modules",

        "Value":
            role_counts[
                "Supporting_q005"
            ],

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Exploratory_reported_modules",

        "Value":
            role_counts[
                "Exploratory_reported"
            ],

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Gene_level_module_rows",

        "Value":
            len(
                gene_ranking_rows
            ),

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Top10_representative_rows",

        "Value":
            len(
                top10_rows
            ),

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Top5_primary_supporting_rows",

        "Value":
            len(
                figure_top5_rows
            ),

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Full_network_nodes",

        "Value":
            len(
                full_nodes
            ),

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Full_network_edges",

        "Value":
            len(
                full_edges
            ),

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Representative_network_nodes",

        "Value":
            len(
                figure_nodes
            ),

        "Expected":
            "INFO",

        "Status":
            "INFO"
    },

    {
        "Metric":
            "Representative_network_edges",

        "Value":
            len(
                figure_edges
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
        "09C10_QC_summary.tsv"
    ),
    qc,
    [
        "Metric",
        "Value",
        "Expected",
        "Status"
    ]
)


# ============================================================
# 15. Column / interpretation dictionary
# ============================================================

dictionary = [

    {
        "Item":
            "Primary_structural",

        "Meaning":
            "Primary prioritized motif with q<=0.05 FIMO localization and family-level promoter/TSS structural support."
    },

    {
        "Item":
            "Primary_integrated",

        "Meaning":
            "Primary prioritized motif with q<=0.05 FIMO localization but no significant family-level promoter/TSS shift."
    },

    {
        "Item":
            "Supporting_q005_structural",

        "Meaning":
            "Secondary motif with q<=0.05 FIMO localization and family-level structural support."
    },

    {
        "Item":
            "Supporting_q005",

        "Meaning":
            "Secondary motif with q<=0.05 FIMO localization and no significant family-level structural shift."
    },

    {
        "Item":
            "Exploratory_reported",

        "Meaning":
            "Exploratory module based on FIMO p<=1e-4 reported sites because no selected q<=0.05 sites were available."
    },

    {
        "Item":
            "Representative_candidate_TFs",

        "Meaning":
            "Uses tissue-preferred candidates when available, otherwise high-RNA candidates, otherwise motif annotation candidates."
    },

    {
        "Item":
            "Tomtom_family_level_candidate",

        "Meaning":
            "TF assignment for a de novo motif is family-level sequence similarity and does not demonstrate TF occupancy."
    },

    {
        "Item":
            "Known_motif_identity_candidate",

        "Meaning":
            "Candidate TF identity derives from a known JASPAR motif, but FIMO localization still does not prove TF binding."
    },

    {
        "Item":
            "Module_motif_positive_peak_N_for_gene",

        "Meaning":
            "Number of strong motif-positive peaks associated with the same gene within the same motif module."
    },

    {
        "Item":
            "Representative_pair_selection_rule",

        "Meaning":
            "Deterministic hierarchical ranking; no arbitrary weighted evidence score was used."
    },

    {
        "Item":
            "candidate_TF_to_motif",

        "Meaning":
            "Candidate annotation edge only; not experimentally validated TF occupancy."
    },

    {
        "Item":
            "motif_to_peak",

        "Meaning":
            "FIMO motif-site localization edge."
    },

    {
        "Item":
            "peak_to_associated_gene",

        "Meaning":
            "Strong RNA-ATAC coupled UROPA genomic association; not validated causal regulation."
    }
]


write_tsv(
    os.path.join(
        OUTDIR,
        "09C11_interpretation_dictionary.tsv"
    ),
    dictionary,
    [
        "Item",
        "Meaning"
    ]
)


# ============================================================
# 16. Console summary
# ============================================================

print("\n" + "=" * 80)
print("STEP 09C COMPLETED")
print("=" * 80)

print(
    f"Selected modules                  : {len(module_rows)}"
)

print(
    f"Module-specific peak-gene rows    : {len(all_module_pairs)}"
)

print(
    f"Gene-level module rows            : {len(gene_ranking_rows)}"
)

print(
    f"Top10 representative rows         : {len(top10_rows)}"
)

print(
    f"Top5 primary/supporting rows       : {len(figure_top5_rows)}"
)


print("\nModule roles:")

for role in [
    "Primary_structural",
    "Primary_integrated",
    "Supporting_q005_structural",
    "Supporting_q005",
    "Exploratory_reported"
]:

    print(
        f"{role:32s}: "
        f"{role_counts[role]}"
    )


print("\nModule summary:")

for row in module_summary:

    print(

        f"{row['Tissue']:12s} | "
        f"{row['Motif_ID']:22s} | "
        f"{row['Module_role']:28s} | "
        f"pairs={row['Candidate_pair_N']:3d} | "
        f"genes={row['Unique_gene_N']:3d} | "
        f"promoter={float(row['Promoter_pair_fraction']):.3f}"
    )


print("\nNetwork tables:")

print(
    f"Full network nodes                : {len(full_nodes)}"
)

print(
    f"Full network edges                : {len(full_edges)}"
)

print(
    f"Representative network nodes      : {len(figure_nodes)}"
)

print(
    f"Representative network edges      : {len(figure_edges)}"
)


print("\nInterpretation boundaries:")

print(
    "1. Module roles are integrated evidence labels, not significance levels."
)

print(
    "2. Representative gene ranking uses deterministic hierarchical rules, not a weighted score."
)

print(
    "3. De novo motif TF assignments remain family-level candidates."
)

print(
    "4. Motif->peak is sequence localization, not TF occupancy."
)

print(
    "5. Peak->gene remains a strong coupled UROPA association, not validated causal regulation."
)


print("\nMain outputs:")

for x in [

    "09C1_module_definition.tsv",

    "09C2_all_module_peak_gene_candidates.tsv",

    "09C3_gene_level_candidate_ranking.tsv",

    "09C4_top10_representative_genes_per_module.tsv",

    "09C5_top5_primary_supporting_figure_candidates.tsv",

    "09C6_module_candidate_summary.tsv",

    "09C7_within_tissue_module_overlap.tsv",

    "09C8_full_candidate_network_nodes.tsv",

    "09C8_full_candidate_network_edges.tsv",

    "09C9_representative_figure_network_nodes.tsv",

    "09C9_representative_figure_network_edges.tsv",

    "09C10_QC_summary.tsv",

    "09C11_interpretation_dictionary.tsv"
]:

    print(
        os.path.join(
            OUTDIR,
            x
        )
    )


print("=" * 80)

