#!/usr/bin/env python3

import os
import re
import csv
import glob
import math
import statistics
from collections import defaultdict

# ============================================================
# Configuration
# ============================================================

MAP_FILE = (
    "08E1b_unique_peak_FASTA/"
    "08E1b_sequence_ID_mapping.tsv"
)

COUPLING_FILE = (
    "05_RNA_ATAC_quantitative_coupling/"
    "05_concordant_HC_peak_gene_correlations.tsv"
)

MOTIF_FAMILY_FILE = (
    "08D_motif_family_TF_master/"
    "08D1_publication_motif_family_master.tsv"
)

TF_MASTER_FILE = (
    "08D_motif_family_TF_master/"
    "08D2_publication_TF_candidate_master.tsv"
)

FIMO_ROOT = "08E1c_FIMO_unique_peak_ids"

FIMO_SUMMARY_FILE = (
    "08E1c_FIMO_unique_peak_ids/"
    "08E1c_FIMO_per_motif_summary.tsv"
)

GTF_FILE = "Sus_longest.gtf"

STEP07_FILE = (
    "07_promoter_first_classification/"
    "07A_all_UROPA_assigned_promoter_first.tsv"
)

OUTDIR = "08E2_motif_peak_gene_TF_network"

TISSUES = [
    "Muscle",
    "Spleen",
    "Liver",
    "Cerebellum"
]

EPS = 1e-12

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# Utility functions
# ============================================================

def read_tsv(path, skip_comments=False):
    rows = []

    with open(path, "r", encoding="utf-8") as f:
        lines = []

        for line in f:
            if skip_comments and line.startswith("#"):
                continue

            if not line.strip():
                continue

            lines.append(line)

    if not lines:
        return [], []

    reader = csv.DictReader(lines, delimiter="\t")

    fields = reader.fieldnames

    for row in reader:
        rows.append(row)

    return rows, fields


def write_tsv(path, rows, fields):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:
            out = {}

            for col in fields:
                value = row.get(col, "")

                if value is None:
                    value = ""

                out[col] = value

            writer.writerow(out)


def as_float(x):
    if x is None:
        return None

    x = str(x).strip()

    if x in ("", "NA", "NaN", "nan", "."):
        return None

    try:
        v = float(x)

        if math.isnan(v):
            return None

        return v
    except Exception:
        return None


def as_int(x):
    try:
        return int(float(x))
    except Exception:
        return None


def is_true(x):
    return str(x).strip().lower() in (
        "true", "t", "1", "yes", "y"
    )


def median_numeric(values):
    vals = [
        x for x in values
        if x is not None
    ]

    if not vals:
        return ""

    return statistics.median(vals)


def split_semicolon(x):
    if x is None:
        return []

    return [
        y.strip()
        for y in str(x).split(";")
        if y.strip()
    ]


# ============================================================
# 1. Read sequence-ID mapping
# ============================================================

print("=" * 80)
print("STEP 08E2")
print("=" * 80)

print("\nReading unique FASTA sequence mapping...")

mapping_rows, mapping_fields = read_tsv(MAP_FILE)

required_map = {
    "Tissue",
    "Set",
    "Sequence_ID",
    "BED_chr",
    "BED_start",
    "BED_end",
    "BED_name"
}

missing = required_map - set(mapping_fields)

if missing:
    raise RuntimeError(
        "Missing mapping columns: " + ",".join(sorted(missing))
    )

fg_mapping = [
    x for x in mapping_rows
    if x["Set"] == "FG" and x["Tissue"] in TISSUES
]

seq_map = {}

for row in fg_mapping:
    sid = row["Sequence_ID"]

    if sid in seq_map:
        raise RuntimeError(
            f"Duplicated Sequence_ID: {sid}"
        )

    seq_map[sid] = row

print(f"Foreground sequence mappings: {len(fg_mapping)}")

for tissue in TISSUES:
    n = sum(
        1 for x in fg_mapping
        if x["Tissue"] == tissue
    )

    print(f"{tissue:12s}: {n}")


# ============================================================
# 2. Read concordant HC peak-gene pairs and reconstruct
#    the official 702 strong set
# ============================================================

print("\nReading Step05 coupling table...")

coupling_rows, coupling_fields = read_tsv(
    COUPLING_FILE
)

required_coupling = {
    "SAF_peak_id",
    "BED_peak_id",
    "gene_id",
    "ATAC_Max_tissue",
    "rho_P348",
    "rho_P350",
    "rho_mean"
}

missing = required_coupling - set(coupling_fields)

if missing:
    raise RuntimeError(
        "Missing coupling columns: " +
        ",".join(sorted(missing))
    )

strong_rows = []

for row in coupling_rows:

    r348 = as_float(row["rho_P348"])
    r350 = as_float(row["rho_P350"])
    rmean = as_float(row["rho_mean"])

    if None in (r348, r350, rmean):
        continue

    if (
        r348 >= 0.5 - EPS and
        r350 >= 0.5 - EPS and
        rmean >= 0.7 - EPS
    ):
        strong_rows.append(row)

print(f"Strong peak-gene pairs reconstructed: {len(strong_rows)}")

if len(strong_rows) != 702:
    raise RuntimeError(
        "Strong-pair reconstruction did not produce 702 rows. "
        f"Observed={len(strong_rows)}"
    )

strong_saf = {}
strong_bed = {}

for row in strong_rows:

    saf = row["SAF_peak_id"]
    bed = row["BED_peak_id"]

    if saf in strong_saf:
        raise RuntimeError(
            f"Duplicated strong SAF peak: {saf}"
        )

    if bed in strong_bed:
        raise RuntimeError(
            f"Duplicated strong BED peak: {bed}"
        )

    strong_saf[saf] = row
    strong_bed[bed] = row


# ============================================================
# 3. Determine whether BED_name corresponds to SAF or BED ID
# ============================================================

fg_names = set(
    x["BED_name"] for x in fg_mapping
)

saf_match = len(
    fg_names.intersection(strong_saf.keys())
)

bed_match = len(
    fg_names.intersection(strong_bed.keys())
)

print("\nPeak ID bridge check:")
print(f"BED_name matching SAF_peak_id: {saf_match}")
print(f"BED_name matching BED_peak_id: {bed_match}")

if saf_match == len(fg_mapping):
    mapping_peak_mode = "SAF"

elif bed_match == len(fg_mapping):
    mapping_peak_mode = "BED"

else:
    raise RuntimeError(
        "Sequence mapping cannot be fully linked to "
        "the reconstructed strong peak set."
    )

print(
    "Mapping BED_name corresponds to:",
    mapping_peak_mode + "_peak_id"
)


# ============================================================
# 4. Verify all 671 FIMO foreground peaks correspond to
#    the strong set and correct tissue
# ============================================================

for row in fg_mapping:

    peak_name = row["BED_name"]

    if mapping_peak_mode == "SAF":
        strong = strong_saf[peak_name]
    else:
        strong = strong_bed[peak_name]

    if strong["ATAC_Max_tissue"] != row["Tissue"]:
        raise RuntimeError(
            "Tissue mismatch for peak "
            f"{peak_name}: mapping={row['Tissue']} "
            f"strong={strong['ATAC_Max_tissue']}"
        )

print("PASS: all foreground sequences map to strong peak-gene pairs.")


# ============================================================
# 5. Gene ID -> gene name from representative GTF
# ============================================================

print("\nReading target-gene names from GTF...")

gene_name_map = {}

gene_id_re = re.compile(
    r'gene_id "([^"]+)"'
)

gene_name_re = re.compile(
    r'gene_name "([^"]+)"'
)

with open(
    GTF_FILE,
    "r",
    encoding="utf-8"
) as f:

    for line in f:

        if line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")

        if len(cols) < 9:
            continue

        attr = cols[8]

        m_id = gene_id_re.search(attr)

        if not m_id:
            continue

        gid = m_id.group(1)

        m_name = gene_name_re.search(attr)

        if m_name:
            gname = m_name.group(1)
        else:
            gname = ""

        if gid not in gene_name_map:
            gene_name_map[gid] = gname

print(
    "Gene IDs with GTF mapping:",
    len(gene_name_map)
)


# ============================================================
# 6. Read Step08D motif-family master
# ============================================================

print("\nReading Step08D motif-family master...")

family_rows, family_fields = read_tsv(
    MOTIF_FAMILY_FILE
)

family_index = {}

for row in family_rows:

    key = (
        row["Tissue"],
        row["Representative_motif_ID"]
    )

    if key in family_index:
        raise RuntimeError(
            f"Duplicate motif family key: {key}"
        )

    family_index[key] = row

print(
    "Motif families:",
    len(family_rows)
)


# ============================================================
# 7. Read TF candidate master and index by motif family
# ============================================================

print("\nReading Step08D TF candidate master...")

tf_rows, tf_fields = read_tsv(
    TF_MASTER_FILE
)

tf_by_family = defaultdict(list)

for row in tf_rows:

    tissue = row["Tissue"]

    families = split_semicolon(
        row["Motif_family_IDs"]
    )

    for family_id in families:
        tf_by_family[
            (tissue, family_id)
        ].append(row)

print(
    "Mapped TF candidate rows:",
    len(tf_rows)
)


# ============================================================
# 8. Optional Step07 promoter-first annotation
# ============================================================

step07_index = {}
step07_selected_fields = []
step07_key_mode = None

if os.path.exists(STEP07_FILE):

    print("\nReading Step07 promoter-first annotation...")

    s7_rows, s7_fields = read_tsv(
        STEP07_FILE
    )

    if (
        "SAF_peak_id" in s7_fields and
        "gene_id" in s7_fields
    ):
        step07_key_mode = "SAF"

        for row in s7_rows:
            step07_index[
                (
                    row["SAF_peak_id"],
                    row["gene_id"]
                )
            ] = row

    elif (
        "BED_peak_id" in s7_fields and
        "gene_id" in s7_fields
    ):
        step07_key_mode = "BED"

        for row in s7_rows:
            step07_index[
                (
                    row["BED_peak_id"],
                    row["gene_id"]
                )
            ] = row

    else:
        print(
            "WARNING: Step07 file found but no usable "
            "peak-ID + gene-ID key detected."
        )

    # Carry biologically useful classification / TSS / distance fields
    for col in s7_fields:

        low = col.lower()

        if (
            "class" in low or
            "tss" in low or
            "distance" in low
        ):
            step07_selected_fields.append(col)

    print(
        "Step07 join mode:",
        step07_key_mode
    )

    print(
        "Step07 fields carried forward:",
        ", ".join(step07_selected_fields)
    )

else:
    print(
        "\nWARNING: Step07 file not found. "
        "Step08E2 will continue without promoter-first fields."
    )


# ============================================================
# 9. Read all selected FIMO results
# ============================================================

print("\nReading FIMO results...")

fimo_files = sorted(
    glob.glob(
        os.path.join(
            FIMO_ROOT,
            "*",
            "*",
            "*",
            "fimo.tsv"
        )
    )
)

if not fimo_files:
    raise RuntimeError(
        "No FIMO result files found."
    )

site_rows = []

for path in fimo_files:

    rel = os.path.relpath(
        path,
        FIMO_ROOT
    )

    parts = rel.split(os.sep)

    if len(parts) != 4:
        raise RuntimeError(
            f"Unexpected FIMO path: {path}"
        )

    tissue, tier, source, _ = parts

    rows, fields = read_tsv(
        path,
        skip_comments=True
    )

    required_fimo = {
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p-value",
        "q-value",
        "matched_sequence"
    }

    missing = required_fimo - set(fields)

    if missing:
        raise RuntimeError(
            f"{path}: missing FIMO columns "
            + ",".join(sorted(missing))
        )

    for row in rows:

        p = as_float(row["p-value"])

        if p is None:
            continue

        # Formal reported-site threshold used in Step08E1c
        if p > 1e-4 + EPS:
            continue

        row["_Tissue"] = tissue
        row["_Tier"] = tier
        row["_Source"] = source

        site_rows.append(row)

print(
    "FIMO reported sites p<=1e-4:",
    len(site_rows)
)


# ============================================================
# 10. Collapse FIMO sites to Tissue × Motif × Peak
# ============================================================

grouped_sites = defaultdict(list)

for row in site_rows:

    key = (
        row["_Tissue"],
        row["_Tier"],
        row["_Source"],
        row["motif_id"],
        row["sequence_name"]
    )

    grouped_sites[key].append(row)

motif_peak_rows = []

missing_family = 0
missing_sequence = 0
missing_strong = 0

for key, sites in grouped_sites.items():

    tissue, tier, source, motif_id, sequence_id = key

    if sequence_id not in seq_map:
        missing_sequence += 1
        continue

    maprow = seq_map[sequence_id]

    peak_name = maprow["BED_name"]

    if mapping_peak_mode == "SAF":

        if peak_name not in strong_saf:
            missing_strong += 1
            continue

        strong = strong_saf[peak_name]

    else:

        if peak_name not in strong_bed:
            missing_strong += 1
            continue

        strong = strong_bed[peak_name]

    family = family_index.get(
        (tissue, motif_id)
    )

    if family is None:
        missing_family += 1
        continue

    def site_sort_key(x):

        p = as_float(x["p-value"])
        q = as_float(x["q-value"])
        score = as_float(x["score"])

        if q is None:
            q = float("inf")

        if score is None:
            score = float("-inf")

        return (
            p if p is not None else float("inf"),
            q,
            -score
        )

    best = sorted(
        sites,
        key=site_sort_key
    )[0]

    q_values = [
        as_float(x["q-value"])
        for x in sites
    ]

    q_values = [
        x for x in q_values
        if x is not None
    ]

    best_q = (
        min(q_values)
        if q_values
        else None
    )

    best_p = min(
        as_float(x["p-value"])
        for x in sites
        if as_float(x["p-value"]) is not None
    )

    scores = [
        as_float(x["score"])
        for x in sites
        if as_float(x["score"]) is not None
    ]

    best_score = (
        max(scores)
        if scores
        else None
    )

    q005 = (
        best_q is not None and
        best_q <= 0.05 + EPS
    )

    bed_start = as_int(
        maprow["BED_start"]
    )

    site_start = as_int(
        best["start"]
    )

    site_stop = as_int(
        best["stop"]
    )

    genomic_start_1based = ""
    genomic_end_1based = ""

    if (
        bed_start is not None and
        site_start is not None and
        site_stop is not None
    ):

        rel_start = min(
            site_start,
            site_stop
        )

        rel_end = max(
            site_start,
            site_stop
        )

        genomic_start_1based = (
            bed_start + rel_start
        )

        genomic_end_1based = (
            bed_start + rel_end
        )

    target_gene_id = strong["gene_id"]

    out = {
        "Tissue": tissue,
        "FIMO_tier": tier,
        "Motif_source": source,

        "Motif_ID": motif_id,
        "Motif_alt_ID": best["motif_alt_id"],

        "Motif_family_ID":
            family["Motif_family_ID"],

        "Motif_family_label":
            family["Motif_family_label"],

        "Representative_motif_name":
            family["Representative_motif_name"],

        "Representative_source":
            family["Representative_source"],

        "Consensus":
            family["Consensus"],

        "Candidate_TFs":
            family["Candidate_TFs"],

        "High_RNA_candidate_TFs":
            family["High_RNA_candidate_TFs"],

        "Tissue_preferred_candidate_TFs":
            family["Tissue_preferred_candidate_TFs"],

        "AME_exploratory_support":
            family["AME_exploratory_support"],

        "AME_best_adj_p":
            family["AME_best_adj_p"],

        "AME_best_E":
            family["AME_best_E"],

        "Default_STREME_holdout_E":
            family["Default_STREME_holdout_E"],

        "CrossRun_reproduced":
            family["CrossRun_reproduced"],

        "Best_default_Tomtom_q":
            family["Best_default_Tomtom_q"],

        "Sequence_ID":
            sequence_id,

        "SAF_peak_id":
            strong["SAF_peak_id"],

        "BED_peak_id":
            strong["BED_peak_id"],

        "peak_chr":
            strong.get(
                "peak_chr",
                maprow["BED_chr"]
            ),

        "BED_start":
            strong.get(
                "BED_start",
                maprow["BED_start"]
            ),

        "BED_end":
            strong.get(
                "BED_end",
                maprow["BED_end"]
            ),

        "target_gene_id":
            target_gene_id,

        "target_gene_name":
            gene_name_map.get(
                target_gene_id,
                ""
            ),

        "UROPA_feature":
            strong.get(
                "UROPA_feature",
                ""
            ),

        "UROPA_relative_location":
            strong.get(
                "UROPA_relative_location",
                ""
            ),

        "UROPA_distance":
            strong.get(
                "UROPA_distance",
                ""
            ),

        "rho_P348":
            strong.get(
                "rho_P348",
                ""
            ),

        "rho_P350":
            strong.get(
                "rho_P350",
                ""
            ),

        "rho_mean":
            strong.get(
                "rho_mean",
                ""
            ),

        "rho_all16":
            strong.get(
                "rho_all16",
                ""
            ),

        "Mean_animal_rho":
            strong.get(
                "Mean_animal_rho",
                ""
            ),

        "Min_animal_rho":
            strong.get(
                "Min_animal_rho",
                ""
            ),

        "ATAC_Tau":
            strong.get(
                "ATAC_Tau",
                ""
            ),

        "ATAC_Max_SPM":
            strong.get(
                "ATAC_Max_SPM",
                ""
            ),

        "RNA_Tau":
            strong.get(
                "RNA_Tau",
                ""
            ),

        "RNA_Max_SPM":
            strong.get(
                "RNA_Max_SPM",
                ""
            ),

        "FIMO_site_count_p_le_1e-4":
            len(sites),

        "FIMO_site_count_q_le_0.05":
            sum(
                1
                for x in sites
                if (
                    as_float(x["q-value"]) is not None
                    and
                    as_float(x["q-value"])
                    <= 0.05 + EPS
                )
            ),

        "FIMO_best_score":
            best_score,

        "FIMO_best_p":
            best_p,

        "FIMO_best_q":
            (
                best_q
                if best_q is not None
                else ""
            ),

        "FIMO_q_le_0.05":
            q005,

        "FIMO_evidence_level":
            (
                "High_confidence_q_le_0.05"
                if q005
                else
                "Reported_p_le_1e-4_only"
            ),

        "FIMO_best_start":
            best["start"],

        "FIMO_best_stop":
            best["stop"],

        "FIMO_best_strand":
            best["strand"],

        "FIMO_best_matched_sequence":
            best["matched_sequence"],

        "FIMO_best_site_chr":
            maprow["BED_chr"],

        "FIMO_best_site_start_1based":
            genomic_start_1based,

        "FIMO_best_site_end_1based":
            genomic_end_1based,
    }

    # --------------------------------------------------------
    # Attach Step07 structural/TSS classification where possible
    # --------------------------------------------------------

    if step07_key_mode is not None:

        if step07_key_mode == "SAF":
            s7key = (
                strong["SAF_peak_id"],
                target_gene_id
            )
        else:
            s7key = (
                strong["BED_peak_id"],
                target_gene_id
            )

        s7 = step07_index.get(s7key)

        if s7 is not None:

            for col in step07_selected_fields:
                out[
                    "Step07_" + col
                ] = s7.get(col, "")

        else:

            for col in step07_selected_fields:
                out[
                    "Step07_" + col
                ] = ""

    motif_peak_rows.append(out)


if missing_sequence:
    raise RuntimeError(
        f"FIMO sequence IDs without mapping: {missing_sequence}"
    )

if missing_strong:
    raise RuntimeError(
        f"FIMO peaks without strong pair: {missing_strong}"
    )

if missing_family:
    raise RuntimeError(
        f"FIMO motifs without motif-family mapping: {missing_family}"
    )


# ============================================================
# 11. Sort publication master
# ============================================================

def tier_order(x):
    return {
        "primary": 1,
        "secondary": 2
    }.get(
        str(x).lower(),
        9
    )


motif_peak_rows.sort(
    key=lambda x: (
        TISSUES.index(x["Tissue"]),
        tier_order(x["FIMO_tier"]),
        x["Motif_source"],
        x["Motif_ID"],
        x["SAF_peak_id"]
    )
)


# ============================================================
# 12. Core output fields
# ============================================================

core_fields = [
    "Tissue",
    "FIMO_tier",
    "Motif_source",

    "Motif_family_ID",
    "Motif_family_label",
    "Motif_ID",
    "Motif_alt_ID",
    "Representative_motif_name",
    "Representative_source",
    "Consensus",

    "Candidate_TFs",
    "High_RNA_candidate_TFs",
    "Tissue_preferred_candidate_TFs",

    "AME_exploratory_support",
    "AME_best_adj_p",
    "AME_best_E",
    "Default_STREME_holdout_E",
    "CrossRun_reproduced",
    "Best_default_Tomtom_q",

    "Sequence_ID",

    "SAF_peak_id",
    "BED_peak_id",
    "peak_chr",
    "BED_start",
    "BED_end",

    "target_gene_id",
    "target_gene_name",

    "UROPA_feature",
    "UROPA_relative_location",
    "UROPA_distance",

    "rho_P348",
    "rho_P350",
    "rho_mean",
    "rho_all16",
    "Mean_animal_rho",
    "Min_animal_rho",

    "ATAC_Tau",
    "ATAC_Max_SPM",
    "RNA_Tau",
    "RNA_Max_SPM",

    "FIMO_site_count_p_le_1e-4",
    "FIMO_site_count_q_le_0.05",

    "FIMO_best_score",
    "FIMO_best_p",
    "FIMO_best_q",

    "FIMO_q_le_0.05",
    "FIMO_evidence_level",

    "FIMO_best_start",
    "FIMO_best_stop",
    "FIMO_best_strand",
    "FIMO_best_matched_sequence",

    "FIMO_best_site_chr",
    "FIMO_best_site_start_1based",
    "FIMO_best_site_end_1based",
]

step07_output_fields = [
    "Step07_" + x
    for x in step07_selected_fields
]

master_fields = (
    core_fields +
    step07_output_fields
)


# ============================================================
# 13. Write all reported motif–peak–gene links
# ============================================================

out_all = os.path.join(
    OUTDIR,
    "08E2A_all_reported_motif_peak_gene_links.tsv"
)

write_tsv(
    out_all,
    motif_peak_rows,
    master_fields
)


# ============================================================
# 14. High-confidence q<=0.05 motif–peak–gene links
# ============================================================

q005_rows = [
    x for x in motif_peak_rows
    if x["FIMO_q_le_0.05"]
]

out_q = os.path.join(
    OUTDIR,
    "08E2B_high_confidence_q005_motif_peak_gene_links.tsv"
)

write_tsv(
    out_q,
    q005_rows,
    master_fields
)


# ============================================================
# 15. Motif-level target summary
# ============================================================

motif_groups = defaultdict(list)

for row in motif_peak_rows:

    key = (
        row["Tissue"],
        row["FIMO_tier"],
        row["Motif_source"],
        row["Motif_ID"]
    )

    motif_groups[key].append(row)


summary_rows = []

for key, rows in motif_groups.items():

    tissue, tier, source, motif_id = key

    qrows = [
        x for x in rows
        if x["FIMO_q_le_0.05"]
    ]

    target_genes = set(
        x["target_gene_id"]
        for x in rows
    )

    q_target_genes = set(
        x["target_gene_id"]
        for x in qrows
    )

    family_id = rows[0]["Motif_family_ID"]

    summary_rows.append({
        "Tissue":
            tissue,

        "FIMO_tier":
            tier,

        "Motif_source":
            source,

        "Motif_family_ID":
            family_id,

        "Motif_family_label":
            rows[0]["Motif_family_label"],

        "Motif_ID":
            motif_id,

        "Motif_alt_ID":
            rows[0]["Motif_alt_ID"],

        "Candidate_TFs":
            rows[0]["Candidate_TFs"],

        "High_RNA_candidate_TFs":
            rows[0]["High_RNA_candidate_TFs"],

        "Tissue_preferred_candidate_TFs":
            rows[0]["Tissue_preferred_candidate_TFs"],

        "Reported_peak_N_p_le_1e-4":
            len(rows),

        "High_confidence_peak_N_q_le_0.05":
            len(qrows),

        "Unique_target_gene_N":
            len(target_genes),

        "High_confidence_unique_target_gene_N":
            len(q_target_genes),

        "Median_rho_mean_all_reported":
            median_numeric([
                as_float(x["rho_mean"])
                for x in rows
            ]),

        "Median_rho_mean_q_le_0.05":
            median_numeric([
                as_float(x["rho_mean"])
                for x in qrows
            ]),
    })


summary_rows.sort(
    key=lambda x: (
        TISSUES.index(x["Tissue"]),
        tier_order(x["FIMO_tier"]),
        x["Motif_ID"]
    )
)

summary_fields = [
    "Tissue",
    "FIMO_tier",
    "Motif_source",
    "Motif_family_ID",
    "Motif_family_label",
    "Motif_ID",
    "Motif_alt_ID",
    "Candidate_TFs",
    "High_RNA_candidate_TFs",
    "Tissue_preferred_candidate_TFs",
    "Reported_peak_N_p_le_1e-4",
    "High_confidence_peak_N_q_le_0.05",
    "Unique_target_gene_N",
    "High_confidence_unique_target_gene_N",
    "Median_rho_mean_all_reported",
    "Median_rho_mean_q_le_0.05",
]

write_tsv(
    os.path.join(
        OUTDIR,
        "08E2C_motif_level_target_summary.tsv"
    ),
    summary_rows,
    summary_fields
)


# ============================================================
# 16. Expand motif–peak–gene links to TF candidates
# ============================================================

tf_expanded_rows = []

for row in motif_peak_rows:

    tissue = row["Tissue"]
    family_id = row["Motif_family_ID"]

    candidates = tf_by_family.get(
        (tissue, family_id),
        []
    )

    if not candidates:

        new = dict(row)

        new.update({
            "TF_symbol": "",
            "TF_gene_id": "",
            "TF_GTF_gene_name": "",
            "TF_Evidence_sources": "",
            "TF_target_tissue_mean_TPM": "",
            "TF_target_tissue_P348_TPM": "",
            "TF_target_tissue_P350_TPM": "",
            "TF_High_RNA_expression_support": "",
            "TF_Tau": "",
            "TF_target_tissue_SPM": "",
            "TF_target_tissue_expression_rank": "",
            "TF_max_expression_tissue": "",
            "TF_target_is_max_expression_tissue": "",
            "TF_SPM_ge_0.5": "",
            "TF_RNA_support_role": "",
            "TF_Interpretation": "",
        })

        tf_expanded_rows.append(new)

        continue

    for tf in candidates:

        new = dict(row)

        new.update({
            "TF_symbol":
                tf.get(
                    "TF_symbol",
                    ""
                ),

            "TF_gene_id":
                tf.get(
                    "gene_id",
                    ""
                ),

            "TF_GTF_gene_name":
                tf.get(
                    "GTF_gene_name",
                    ""
                ),

            "TF_Evidence_sources":
                tf.get(
                    "Evidence_sources",
                    ""
                ),

            "TF_target_tissue_mean_TPM":
                tf.get(
                    "Target_tissue_mean_TPM",
                    ""
                ),

            "TF_target_tissue_P348_TPM":
                tf.get(
                    "Target_tissue_P348_TPM",
                    ""
                ),

            "TF_target_tissue_P350_TPM":
                tf.get(
                    "Target_tissue_P350_TPM",
                    ""
                ),

            "TF_High_RNA_expression_support":
                tf.get(
                    "High_RNA_expression_support",
                    ""
                ),

            "TF_Tau":
                tf.get(
                    "Tau",
                    ""
                ),

            "TF_target_tissue_SPM":
                tf.get(
                    "Target_tissue_SPM",
                    ""
                ),

            "TF_target_tissue_expression_rank":
                tf.get(
                    "Target_tissue_expression_rank",
                    ""
                ),

            "TF_max_expression_tissue":
                tf.get(
                    "Max_expression_tissue",
                    ""
                ),

            "TF_target_is_max_expression_tissue":
                tf.get(
                    "Target_is_max_expression_tissue",
                    ""
                ),

            "TF_SPM_ge_0.5":
                tf.get(
                    "SPM_ge_0.5",
                    ""
                ),

            "TF_RNA_support_role":
                tf.get(
                    "RNA_support_role",
                    ""
                ),

            "TF_Interpretation":
                tf.get(
                    "Interpretation",
                    ""
                ),
        })

        tf_expanded_rows.append(new)


tf_fields_out = master_fields + [
    "TF_symbol",
    "TF_gene_id",
    "TF_GTF_gene_name",
    "TF_Evidence_sources",

    "TF_target_tissue_mean_TPM",
    "TF_target_tissue_P348_TPM",
    "TF_target_tissue_P350_TPM",

    "TF_High_RNA_expression_support",

    "TF_Tau",
    "TF_target_tissue_SPM",
    "TF_target_tissue_expression_rank",

    "TF_max_expression_tissue",
    "TF_target_is_max_expression_tissue",
    "TF_SPM_ge_0.5",

    "TF_RNA_support_role",
    "TF_Interpretation",
]


write_tsv(
    os.path.join(
        OUTDIR,
        "08E2D_TF_motif_peak_gene_master.tsv"
    ),
    tf_expanded_rows,
    tf_fields_out
)


# ============================================================
# 17. High-RNA TF subset
# ============================================================

high_rna_tf_rows = [
    x for x in tf_expanded_rows
    if is_true(
        x.get(
            "TF_High_RNA_expression_support",
            ""
        )
    )
]

write_tsv(
    os.path.join(
        OUTDIR,
        "08E2E_high_RNA_TF_motif_peak_gene_master.tsv"
    ),
    high_rna_tf_rows,
    tf_fields_out
)


# ============================================================
# 18. Tissue-level network summary
# ============================================================

tissue_summary = []

for tissue in TISSUES:

    rows = [
        x for x in motif_peak_rows
        if x["Tissue"] == tissue
    ]

    qrows = [
        x for x in rows
        if x["FIMO_q_le_0.05"]
    ]

    tfrows = [
        x for x in tf_expanded_rows
        if x["Tissue"] == tissue
        and x["TF_symbol"] != ""
    ]

    high_tf = [
        x for x in high_rna_tf_rows
        if x["Tissue"] == tissue
    ]

    tissue_summary.append({
        "Tissue":
            tissue,

        "Selected_motif_N":
            len(set(
                x["Motif_ID"]
                for x in rows
            )),

        "Reported_motif_peak_links":
            len(rows),

        "Q005_motif_peak_links":
            len(qrows),

        "Unique_FIMO_positive_peaks":
            len(set(
                x["SAF_peak_id"]
                for x in rows
            )),

        "Unique_Q005_peaks":
            len(set(
                x["SAF_peak_id"]
                for x in qrows
            )),

        "Unique_target_genes":
            len(set(
                x["target_gene_id"]
                for x in rows
            )),

        "Unique_Q005_target_genes":
            len(set(
                x["target_gene_id"]
                for x in qrows
            )),

        "TF_motif_peak_gene_rows":
            len(tfrows),

        "High_RNA_TF_motif_peak_gene_rows":
            len(high_tf),

        "Unique_candidate_TFs":
            len(set(
                x["TF_symbol"]
                for x in tfrows
                if x["TF_symbol"]
            )),

        "Unique_high_RNA_candidate_TFs":
            len(set(
                x["TF_symbol"]
                for x in high_tf
                if x["TF_symbol"]
            )),
    })


tissue_summary_fields = [
    "Tissue",
    "Selected_motif_N",
    "Reported_motif_peak_links",
    "Q005_motif_peak_links",
    "Unique_FIMO_positive_peaks",
    "Unique_Q005_peaks",
    "Unique_target_genes",
    "Unique_Q005_target_genes",
    "TF_motif_peak_gene_rows",
    "High_RNA_TF_motif_peak_gene_rows",
    "Unique_candidate_TFs",
    "Unique_high_RNA_candidate_TFs",
]

write_tsv(
    os.path.join(
        OUTDIR,
        "08E2F_tissue_network_summary.tsv"
    ),
    tissue_summary,
    tissue_summary_fields
)


# ============================================================
# 19. Validate against Step08E1c per-motif summary
# ============================================================

expected_rows, expected_fields = read_tsv(
    FIMO_SUMMARY_FILE
)

observed_by_motif = {}

for row in motif_peak_rows:

    key = (
        row["Tissue"],
        row["FIMO_tier"],
        row["Motif_source"],
        row["Motif_ID"]
    )

    if key not in observed_by_motif:
        observed_by_motif[key] = {
            "sites_p": 0,
            "sites_q": 0,
            "peaks_p": 0,
            "peaks_q": 0
        }

    observed_by_motif[key]["sites_p"] += int(
        row["FIMO_site_count_p_le_1e-4"]
    )

    observed_by_motif[key]["sites_q"] += int(
        row["FIMO_site_count_q_le_0.05"]
    )

    observed_by_motif[key]["peaks_p"] += 1

    if row["FIMO_q_le_0.05"]:
        observed_by_motif[key]["peaks_q"] += 1


validation_status = "PASS"

validation_messages = []

for exp in expected_rows:

    key = (
        exp["Tissue"],
        exp["Tier"],
        exp["Motif_source"],
        exp["Motif_ID"]
    )

    obs = observed_by_motif.get(key)

    if obs is None:
        validation_status = "FAIL"
        validation_messages.append(
            f"Missing observed motif: {key}"
        )
        continue

    expected_values = {
        "sites_p":
            as_int(exp["Sites_p_le_1e-4"]),

        "sites_q":
            as_int(exp["Sites_q_le_0.05"]),

        "peaks_p":
            as_int(exp["Peaks_p_le_1e-4"]),

        "peaks_q":
            as_int(exp["Peaks_q_le_0.05"]),
    }

    for metric in expected_values:

        if obs[metric] != expected_values[metric]:
            validation_status = "FAIL"

            validation_messages.append(
                f"{key} {metric}: "
                f"observed={obs[metric]} "
                f"expected={expected_values[metric]}"
            )


# ============================================================
# 20. QC summary
# ============================================================

expected_reported_peaks = sum(
    as_int(x["Peaks_p_le_1e-4"])
    for x in expected_rows
)

expected_q005_peaks = sum(
    as_int(x["Peaks_q_le_0.05"])
    for x in expected_rows
)

four_tissue_strong = [
    x for x in strong_rows
    if x["ATAC_Max_tissue"] in TISSUES
]

qc_rows = [
    {
        "Metric":
            "Official_strong_peak_gene_pairs",
        "Value":
            len(strong_rows),
        "Expected":
            702,
        "Status":
            "PASS" if len(strong_rows) == 702 else "FAIL"
    },

    {
        "Metric":
            "Strong_pairs_in_FIMO_tissues",
        "Value":
            len(four_tissue_strong),
        "Expected":
            len(fg_mapping),
        "Status":
            (
                "PASS"
                if len(four_tissue_strong) == len(fg_mapping)
                else "FAIL"
            )
    },

    {
        "Metric":
            "Foreground_unique_sequences",
        "Value":
            len(fg_mapping),
        "Expected":
            671,
        "Status":
            (
                "PASS"
                if len(fg_mapping) == 671
                else "FAIL"
            )
    },

    {
        "Metric":
            "Reported_motif_peak_gene_links",
        "Value":
            len(motif_peak_rows),
        "Expected":
            expected_reported_peaks,
        "Status":
            (
                "PASS"
                if len(motif_peak_rows)
                == expected_reported_peaks
                else "FAIL"
            )
    },

    {
        "Metric":
            "High_confidence_q005_motif_peak_gene_links",
        "Value":
            len(q005_rows),
        "Expected":
            expected_q005_peaks,
        "Status":
            (
                "PASS"
                if len(q005_rows)
                == expected_q005_peaks
                else "FAIL"
            )
    },

    {
        "Metric":
            "FIMO_per_motif_reconstruction",
        "Value":
            validation_status,
        "Expected":
            "PASS",
        "Status":
            validation_status
    },

    {
        "Metric":
            "Sequence_mapping_peak_ID_mode",
        "Value":
            mapping_peak_mode,
        "Expected":
            "SAF",
        "Status":
            (
                "PASS"
                if mapping_peak_mode == "SAF"
                else "CHECK"
            )
    },

    {
        "Metric":
            "TF_expanded_rows",
        "Value":
            len(tf_expanded_rows),
        "Expected":
            "",
        "Status":
            "INFO"
    },

    {
        "Metric":
            "High_RNA_TF_expanded_rows",
        "Value":
            len(high_rna_tf_rows),
        "Expected":
            "",
        "Status":
            "INFO"
    },
]

write_tsv(
    os.path.join(
        OUTDIR,
        "08E2G_QC_summary.tsv"
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
# 21. Final console report
# ============================================================

print("\n" + "=" * 80)
print("STEP 08E2 COMPLETED")
print("=" * 80)

print(
    f"Official strong pairs                : {len(strong_rows)}"
)

print(
    f"Strong pairs in four FIMO tissues   : {len(four_tissue_strong)}"
)

print(
    f"Foreground peak sequences           : {len(fg_mapping)}"
)

print(
    f"Reported motif-peak-gene links      : {len(motif_peak_rows)}"
)

print(
    f"High-confidence q<=0.05 links       : {len(q005_rows)}"
)

print(
    f"TF-expanded links                   : {len(tf_expanded_rows)}"
)

print(
    f"High-RNA TF-expanded links          : {len(high_rna_tf_rows)}"
)

print(
    f"FIMO reconstruction QC              : {validation_status}"
)

if validation_messages:

    print("\nValidation problems:")

    for msg in validation_messages:
        print(" -", msg)

print("\nTissue network summary:")

print(
    "\t".join(
        tissue_summary_fields
    )
)

for row in tissue_summary:

    print(
        "\t".join(
            str(row[x])
            for x in tissue_summary_fields
        )
    )

print("\nMain outputs:")

for f in [
    "08E2A_all_reported_motif_peak_gene_links.tsv",
    "08E2B_high_confidence_q005_motif_peak_gene_links.tsv",
    "08E2C_motif_level_target_summary.tsv",
    "08E2D_TF_motif_peak_gene_master.tsv",
    "08E2E_high_RNA_TF_motif_peak_gene_master.tsv",
    "08E2F_tissue_network_summary.tsv",
    "08E2G_QC_summary.tsv",
]:

    print(
        os.path.join(
            OUTDIR,
            f
        )
    )

print("=" * 80)
