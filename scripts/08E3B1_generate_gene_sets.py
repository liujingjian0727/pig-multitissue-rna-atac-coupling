#!/usr/bin/env python3

import os
import re
import csv
from collections import defaultdict

# ============================================================
# STEP 08E3B1
#
# Generate:
#   1. Deduplicated motif-associated gene sets
#   2. Same-tissue strong-gene universes
#
# Main sets:
#   FIMO q <= 0.05
#
# Exploratory sets:
#   FIMO p <= 1e-4
#
# Universe:
#   all unique genes associated with strong peak-gene pairs
#   from the SAME tissue.
# ============================================================

EPS = 1e-12

COUPLING = (
    "05_RNA_ATAC_quantitative_coupling/"
    "05_concordant_HC_peak_gene_correlations.tsv"
)

ALL_LINKS = (
    "08E2_motif_peak_gene_TF_network/"
    "08E2A_all_reported_motif_peak_gene_links.tsv"
)

Q005_LINKS = (
    "08E2_motif_peak_gene_TF_network/"
    "08E2B_high_confidence_q005_motif_peak_gene_links.tsv"
)

GTF = "Sus_longest.gtf"

OUTDIR = "08E3B_gene_sets"

os.makedirs(OUTDIR, exist_ok=True)

os.makedirs(
    os.path.join(OUTDIR, "universes"),
    exist_ok=True
)

os.makedirs(
    os.path.join(OUTDIR, "main_q005"),
    exist_ok=True
)

os.makedirs(
    os.path.join(OUTDIR, "exploratory_reported"),
    exist_ok=True
)


# ============================================================
# Frozen analysis sets
# ============================================================

SELECTED_SETS = [

    # --------------------------------------------------------
    # MAIN q<=0.05 gene sets
    # --------------------------------------------------------

    {
        "set_id":
            "Muscle_KLF9_MA1107.3_q005",

        "tier":
            "Main",

        "evidence":
            "High_confidence_q005",

        "tissue":
            "Muscle",

        "motif_id":
            "MA1107.3",

        "motif_name":
            "KLF9",

        "expected_gene_n":
            96
    },

    {
        "set_id":
            "Muscle_STREME3_CCTCCGCCCCTG_q005",

        "tier":
            "Main",

        "evidence":
            "High_confidence_q005",

        "tissue":
            "Muscle",

        "motif_id":
            "3-CCTCCGCCCCTG",

        "motif_name":
            "STREME-3_GC-rich_SP-KLF-like",

        "expected_gene_n":
            108
    },

    {
        "set_id":
            "Spleen_STREME1_MGGGGSAGGAGCMG_q005",

        "tier":
            "Main",

        "evidence":
            "High_confidence_q005",

        "tissue":
            "Spleen",

        "motif_id":
            "1-MGGGGSAGGAGCMG",

        "motif_name":
            "STREME-1_GC-rich_ZNF-SP-like",

        "expected_gene_n":
            93
    },

    {
        "set_id":
            "Spleen_STREME2_GCCCAGCCCMGCCCM_q005",

        "tier":
            "Main",

        "evidence":
            "High_confidence_q005",

        "tissue":
            "Spleen",

        "motif_id":
            "2-GCCCAGCCCMGCCCM",

        "motif_name":
            "STREME-2_GC-rich_KLF-SP-like",

        "expected_gene_n":
            83
    },

    {
        "set_id":
            "Spleen_KLF12_MA0742.2_q005",

        "tier":
            "Main",

        "evidence":
            "High_confidence_q005",

        "tissue":
            "Spleen",

        "motif_id":
            "MA0742.2",

        "motif_name":
            "KLF12",

        "expected_gene_n":
            48
    },

    {
        "set_id":
            "Liver_STREME3_AGGGGGAGGGGRGGR_q005",

        "tier":
            "Main",

        "evidence":
            "High_confidence_q005",

        "tissue":
            "Liver",

        "motif_id":
            "3-AGGGGGAGGGGRGGR",

        "motif_name":
            "STREME-3_GC-rich_KLF-SP-ZNF-like",

        "expected_gene_n":
            23
    },


    # --------------------------------------------------------
    # EXPLORATORY reported p<=1e-4 sets
    # --------------------------------------------------------

    {
        "set_id":
            "Spleen_ELF4_MA0641.1_reported",

        "tier":
            "Exploratory",

        "evidence":
            "Reported_p1e-4",

        "tissue":
            "Spleen",

        "motif_id":
            "MA0641.1",

        "motif_name":
            "ELF4",

        "expected_gene_n":
            49
    },

    {
        "set_id":
            "Muscle_ZNF684_MA1600.2_reported",

        "tier":
            "Exploratory",

        "evidence":
            "Reported_p1e-4",

        "tissue":
            "Muscle",

        "motif_id":
            "MA1600.2",

        "motif_name":
            "ZNF684",

        "expected_gene_n":
            86
    },
]


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

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
            lineterminator="\n"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow({
                x: row.get(x, "")
                for x in fields
            })


def as_float(x):

    try:
        return float(x)
    except Exception:
        return None


# ============================================================
# Gene names from GTF
# ============================================================

print("=" * 80)
print("STEP 08E3B1")
print("=" * 80)

print("\nReading representative GTF...")

gene_names = {}

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

        attr = fields[8]

        gid_match = gene_id_re.search(
            attr
        )

        if not gid_match:
            continue

        gid = gid_match.group(1)

        gname_match = gene_name_re.search(
            attr
        )

        gname = (
            gname_match.group(1)
            if gname_match
            else ""
        )

        if gid not in gene_names:
            gene_names[gid] = gname


print(
    "GTF gene IDs:",
    len(gene_names)
)


# ============================================================
# Reconstruct official 702 strong pairs
# ============================================================

print("\nReading Step05 coupling table...")

coupling = read_tsv(
    COUPLING
)

strong = []

for row in coupling:

    r348 = as_float(
        row["rho_P348"]
    )

    r350 = as_float(
        row["rho_P350"]
    )

    rmean = as_float(
        row["rho_mean"]
    )

    if None in (
        r348,
        r350,
        rmean
    ):
        continue

    if (
        r348 >= 0.5 - EPS and
        r350 >= 0.5 - EPS and
        rmean >= 0.7 - EPS
    ):
        strong.append(row)


print(
    "Strong pairs reconstructed:",
    len(strong)
)

if len(strong) != 702:

    raise RuntimeError(
        "Expected 702 strong pairs, "
        f"observed {len(strong)}"
    )


# ============================================================
# Same-tissue strong-gene universes
# ============================================================

universe_by_tissue = defaultdict(set)

for row in strong:

    tissue = row[
        "ATAC_Max_tissue"
    ]

    gid = row[
        "gene_id"
    ]

    universe_by_tissue[
        tissue
    ].add(gid)


print("\nSame-tissue strong-gene universes:")

universe_manifest = []

for tissue in sorted(
    universe_by_tissue
):

    genes = sorted(
        universe_by_tissue[tissue]
    )

    txt_file = os.path.join(
        OUTDIR,
        "universes",
        f"{tissue}.same_tissue_strong_gene_universe.txt"
    )

    tsv_file = os.path.join(
        OUTDIR,
        "universes",
        f"{tissue}.same_tissue_strong_gene_universe.tsv"
    )

    with open(
        txt_file,
        "w"
    ) as f:

        for gid in genes:
            print(gid, file=f)

    write_tsv(
        tsv_file,
        [
            {
                "gene_id": gid,
                "gene_name":
                    gene_names.get(
                        gid,
                        ""
                    )
            }
            for gid in genes
        ],
        [
            "gene_id",
            "gene_name"
        ]
    )

    universe_manifest.append({
        "Tissue":
            tissue,

        "Universe_gene_N":
            len(genes),

        "Universe_txt":
            txt_file,

        "Universe_tsv":
            tsv_file
    })

    print(
        f"{tissue:12s}: "
        f"{len(genes)} unique genes"
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "08E3B1_same_tissue_universe_manifest.tsv"
    ),
    universe_manifest,
    [
        "Tissue",
        "Universe_gene_N",
        "Universe_txt",
        "Universe_tsv"
    ]
)


# ============================================================
# Read motif links
# ============================================================

all_links = read_tsv(
    ALL_LINKS
)

q_links = read_tsv(
    Q005_LINKS
)


# ============================================================
# Generate selected gene sets
# ============================================================

print("\nGenerating frozen motif gene sets...")

manifest = []

for config in SELECTED_SETS:

    tissue = config[
        "tissue"
    ]

    motif = config[
        "motif_id"
    ]

    evidence = config[
        "evidence"
    ]

    if evidence == "High_confidence_q005":

        source_rows = q_links

        subdir = "main_q005"

    else:

        source_rows = all_links

        subdir = "exploratory_reported"


    selected_rows = [

        row
        for row in source_rows

        if (
            row["Tissue"] == tissue
            and
            row["Motif_ID"] == motif
        )
    ]


    genes = sorted(set(
        row["target_gene_id"]
        for row in selected_rows
        if row["target_gene_id"]
    ))


    expected = config[
        "expected_gene_n"
    ]


    if len(genes) != expected:

        raise RuntimeError(
            f"{config['set_id']}: "
            f"expected {expected} unique genes, "
            f"observed {len(genes)}"
        )


    universe = universe_by_tissue[
        tissue
    ]


    not_in_universe = [
        gid
        for gid in genes
        if gid not in universe
    ]


    if not_in_universe:

        raise RuntimeError(
            f"{config['set_id']}: "
            "target genes are not a subset "
            "of the same-tissue universe."
        )


    txt_file = os.path.join(
        OUTDIR,
        subdir,
        config["set_id"] +
        ".genes.txt"
    )


    tsv_file = os.path.join(
        OUTDIR,
        subdir,
        config["set_id"] +
        ".genes.tsv"
    )


    with open(
        txt_file,
        "w"
    ) as f:

        for gid in genes:
            print(gid, file=f)


    write_tsv(
        tsv_file,
        [
            {
                "gene_id":
                    gid,

                "gene_name":
                    gene_names.get(
                        gid,
                        ""
                    )
            }
            for gid in genes
        ],
        [
            "gene_id",
            "gene_name"
        ]
    )


    universe_txt = os.path.join(
        OUTDIR,
        "universes",
        f"{tissue}.same_tissue_strong_gene_universe.txt"
    )


    manifest.append({

        "Set_ID":
            config["set_id"],

        "Analysis_tier":
            config["tier"],

        "Evidence_level":
            evidence,

        "Tissue":
            tissue,

        "Motif_ID":
            motif,

        "Motif_name":
            config["motif_name"],

        "Gene_N":
            len(genes),

        "Expected_gene_N":
            expected,

        "Universe_gene_N":
            len(universe),

        "Gene_fraction_of_universe":
            len(genes) /
            len(universe),

        "Gene_list_txt":
            txt_file,

        "Gene_list_tsv":
            tsv_file,

        "Universe_txt":
            universe_txt,

        "QC":
            "PASS"
    })


    print(
        f"{config['set_id']:45s} "
        f"genes={len(genes):3d} "
        f"universe={len(universe):3d} "
        "PASS"
    )


# ============================================================
# Write manifest
# ============================================================

manifest_file = os.path.join(
    OUTDIR,
    "08E3B1_gene_set_manifest.tsv"
)

write_tsv(
    manifest_file,
    manifest,
    [
        "Set_ID",
        "Analysis_tier",
        "Evidence_level",
        "Tissue",
        "Motif_ID",
        "Motif_name",
        "Gene_N",
        "Expected_gene_N",
        "Universe_gene_N",
        "Gene_fraction_of_universe",
        "Gene_list_txt",
        "Gene_list_tsv",
        "Universe_txt",
        "QC"
    ]
)


# ============================================================
# QC
# ============================================================

main_n = sum(
    x["Analysis_tier"] ==
    "Main"
    for x in manifest
)

explore_n = sum(
    x["Analysis_tier"] ==
    "Exploratory"
    for x in manifest
)

qc = [

    {
        "Metric":
            "Official_strong_pairs",

        "Value":
            len(strong),

        "Expected":
            702,

        "Status":
            "PASS"
    },

    {
        "Metric":
            "Selected_gene_sets",

        "Value":
            len(manifest),

        "Expected":
            8,

        "Status":
            (
                "PASS"
                if len(manifest) == 8
                else "FAIL"
            )
    },

    {
        "Metric":
            "Main_q005_gene_sets",

        "Value":
            main_n,

        "Expected":
            6,

        "Status":
            (
                "PASS"
                if main_n == 6
                else "FAIL"
            )
    },

    {
        "Metric":
            "Exploratory_reported_gene_sets",

        "Value":
            explore_n,

        "Expected":
            2,

        "Status":
            (
                "PASS"
                if explore_n == 2
                else "FAIL"
            )
    }
]


write_tsv(
    os.path.join(
        OUTDIR,
        "08E3B1_QC_summary.tsv"
    ),
    qc,
    [
        "Metric",
        "Value",
        "Expected",
        "Status"
    ]
)


print("\n" + "=" * 80)
print("STEP 08E3B1 COMPLETED")
print("=" * 80)

print(
    "Gene-set manifest:",
    manifest_file
)

print(
    "All gene sets are deduplicated "
    "by Ensembl gene ID."
)

print(
    "All target sets are verified subsets "
    "of the matched same-tissue strong-gene universe."
)

print("=" * 80)

