#!/usr/bin/env python3

import os
import re
import sys
import math
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np
import pandas as pd


# ============================================================
# STEP 08C
# TF RNA EXPRESSION + TAU/SPM INTEGRATION
#
# Candidate evidence:
#   1. AME exploratory motifs:
#         adj_p < 0.05
#      NOTE:
#         this is NOT database-wide significance.
#
#   2. Default STREME -> JASPAR Tomtom q < 0.05
#
#   3. No-holdout STREME -> JASPAR Tomtom q < 0.05
#
#   4. Cross-run STREME reproducibility
#         default vs no-holdout Tomtom q < 0.05
#
# RNA:
#   linear TMM.TPM
#
# Tau / SPM:
#   calculated using 8 tissue mean TPM values
#
# IMPORTANT:
#   Tissue specificity is SUPPORTIVE, not mandatory.
# ============================================================


# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------

TPM_FILE = "sus.TMM.TPM.matrix"

META_FILE = "sample_pair_metadata.tsv"

JASPAR_MAP_FILE = "JASPAR2024_motif_ID_to_TF.tsv"

AME_ROOT = "08B1c_AME_final_summary"

DEFAULT_TOMTOM_ROOT = "08B3_Tomtom_JASPAR2024"

NOHOLD_TOMTOM_ROOT = "08B3b_Tomtom_noholdout_JASPAR2024"

REPRO_ROOT = "08B4B_STREME_reproducibility"

DEFAULT_STREME_ROOT = "08B2_STREME_de_novo"

NOHOLD_STREME_ROOT = "08B2b_STREME_noholdout"

OUTDIR = "08C_TF_RNA_integration"


TISSUES = [
    "Adipose",
    "Cerebellum",
    "Cortex",
    "Hypothalamus",
    "Liver",
    "Lung",
    "Muscle",
    "Spleen"
]

MOTIF_TISSUES = [
    "Muscle",
    "Spleen",
    "Liver",
    "Cerebellum"
]


# ------------------------------------------------------------
# RNA support thresholds
#
# These are prioritization thresholds, NOT statistical tests.
# ------------------------------------------------------------

MEAN_TPM_HIGH = 1.0
ANIMAL_TPM_DETECTED = 0.5

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def canonical_gene_id(x):
    x = str(x)
    return re.sub(r"\.\d+$", "", x)


def normalize_col(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def find_col(df, aliases, required=True):
    lookup = {
        normalize_col(c): c
        for c in df.columns
    }

    for a in aliases:
        key = normalize_col(a)
        if key in lookup:
            return lookup[key]

    if required:
        raise ValueError(
            "Cannot find column. Tried: " +
            ", ".join(aliases) +
            "\nAvailable columns:\n" +
            "\n".join(df.columns.astype(str))
        )

    return None


def split_tf_label(x):
    """
    Split composite JASPAR labels such as:
        ETV2::DRGX
        RXRA::VDR
    Preserve ordinary gene symbols.
    """

    if pd.isna(x):
        return []

    x = str(x).strip()

    if x == "" or x.upper() in {"NA", "NONE"}:
        return []

    parts = re.split(r"::|;", x)

    cleaned = []

    for p in parts:

        p = p.strip()

        if p:
            cleaned.append(p)

    return cleaned


def read_table(path):
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        dtype=str
    )


def numeric(x):
    try:
        return float(x)
    except:
        return np.nan


def parse_streme_xml(path):

    stats = {}

    if not os.path.exists(path):
        return stats

    root = ET.parse(path).getroot()

    for elem in root.iter():

        tag = elem.tag.split("}")[-1]

        if tag != "motif":
            continue

        a = elem.attrib

        motif_id = a.get("id")

        stats[motif_id] = {
            "train_pvalue":
                numeric(a.get("train_pvalue")),
            "test_pvalue":
                numeric(a.get("test_pvalue")),
            "test_evalue":
                numeric(a.get("test_evalue")),
            "train_pos_count":
                numeric(a.get("train_pos_count")),
            "train_neg_count":
                numeric(a.get("train_neg_count")),
            "test_pos_count":
                numeric(a.get("test_pos_count")),
            "test_neg_count":
                numeric(a.get("test_neg_count")),
            "width":
                numeric(a.get("width")),
            "seed":
                a.get("seed", "")
        }

    return stats


# ============================================================
# INPUT QC
# ============================================================

required_files = [
    TPM_FILE,
    META_FILE,
    JASPAR_MAP_FILE
]

for f in required_files:

    if not os.path.exists(f):
        sys.exit(
            f"ERROR: required file missing: {f}"
        )


print("=" * 72)
print("STEP 08C: TF RNA EXPRESSION + TAU/SPM INTEGRATION")
print("=" * 72)
print()


# ============================================================
# JASPAR MOTIF -> TF MAP
# ============================================================

jaspar = pd.read_csv(
    JASPAR_MAP_FILE,
    sep="\t",
    header=None,
    names=["motif_ID", "TF_label"],
    dtype=str
)

jaspar_map = dict(
    zip(
        jaspar["motif_ID"],
        jaspar["TF_label"]
    )
)

print(
    "JASPAR motif mappings:",
    len(jaspar_map)
)


# ============================================================
# STREME STATISTICS
# ============================================================

default_streme_stats = {}
nohold_streme_stats = {}

for tissue in MOTIF_TISSUES:

    default_streme_stats[tissue] = parse_streme_xml(
        os.path.join(
            DEFAULT_STREME_ROOT,
            tissue,
            "streme.xml"
        )
    )

    nohold_streme_stats[tissue] = parse_streme_xml(
        os.path.join(
            NOHOLD_STREME_ROOT,
            tissue,
            "streme.xml"
        )
    )


# ============================================================
# CROSS-RUN REPRODUCIBLE MOTIFS
# ============================================================

default_repro = defaultdict(set)
nohold_repro = defaultdict(set)

repro_pairs = []

for tissue in MOTIF_TISSUES:

    f = os.path.join(
        REPRO_ROOT,
        tissue,
        "tomtom.tsv"
    )

    if not os.path.exists(f):
        continue

    x = read_table(f)

    if x.empty:
        continue

    query_col = find_col(
        x,
        ["Query_ID"]
    )

    target_col = find_col(
        x,
        ["Target_ID"]
    )

    q_col = find_col(
        x,
        ["q-value", "q_value"]
    )

    x[q_col] = pd.to_numeric(
        x[q_col],
        errors="coerce"
    )

    x = x[
        x[q_col] < 0.05
    ].copy()

    for _, row in x.iterrows():

        q = row[query_col]
        t = row[target_col]

        default_repro[tissue].add(q)
        nohold_repro[tissue].add(t)

        repro_pairs.append({
            "Tissue": tissue,
            "Default_STREME_motif": q,
            "NoHoldout_STREME_motif": t,
            "CrossRun_Tomtom_q": row[q_col]
        })


pd.DataFrame(
    repro_pairs
).to_csv(
    os.path.join(
        OUTDIR,
        "08C1_STREME_crossrun_reproducible_pairs.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# BUILD LONG-FORM MOTIF -> TF EVIDENCE
# ============================================================

evidence = []


# ------------------------------------------------------------
# AME exploratory candidates
#
# adj_p < 0.05:
# threshold-optimization-adjusted AME statistic.
#
# NOT interpreted as database-wide significance.
# ------------------------------------------------------------

for tissue in MOTIF_TISSUES:

    f = os.path.join(
        AME_ROOT,
        f"{tissue}.AME_all_results.tsv"
    )

    if not os.path.exists(f):
        print(
            "WARNING: AME file missing:",
            f
        )
        continue

    x = read_table(f)

    motif_col = find_col(
        x,
        ["motif_ID", "motif id"]
    )

    alt_col = find_col(
        x,
        [
            "motif_alt_ID",
            "motif_alt_id",
            "motif_alt_ID_DB",
            "motif_alt"
        ],
        required=False
    )

    p_col = find_col(
        x,
        ["p_value", "p-value"],
        required=False
    )

    adj_col = find_col(
        x,
        [
            "adj_p_value",
            "adj_p-value",
            "adjP"
        ]
    )

    e_col = find_col(
        x,
        ["E_value", "E-value"],
        required=False
    )

    or_col = find_col(
        x,
        [
            "OR_HA",
            "odds_ratio_HA",
            "odds_ratio"
        ],
        required=False
    )

    x[adj_col] = pd.to_numeric(
        x[adj_col],
        errors="coerce"
    )

    x = x[
        x[adj_col] < 0.05
    ].copy()

    for _, row in x.iterrows():

        motif_id = row[motif_col]

        if alt_col is not None:
            tf_label = row[alt_col]
        else:
            tf_label = jaspar_map.get(
                motif_id,
                ""
            )

        if (
            pd.isna(tf_label)
            or str(tf_label).strip() == ""
        ):
            tf_label = jaspar_map.get(
                motif_id,
                ""
            )

        components = split_tf_label(
            tf_label
        )

        for tf in components:

            evidence.append({
                "Tissue":
                    tissue,

                "TF_component":
                    tf,

                "Original_TF_label":
                    tf_label,

                "Evidence_source":
                    "AME_exploratory_adjP_lt_0.05",

                "Motif_ID":
                    motif_id,

                "STREME_motif_ID":
                    "",

                "CrossRun_reproducible_motif":
                    False,

                "AME_p":
                    numeric(
                        row[p_col]
                        if p_col else np.nan
                    ),

                "AME_adj_p":
                    numeric(
                        row[adj_col]
                    ),

                "AME_E":
                    numeric(
                        row[e_col]
                        if e_col else np.nan
                    ),

                "AME_OR_HA":
                    numeric(
                        row[or_col]
                        if or_col else np.nan
                    ),

                "Tomtom_q":
                    np.nan,

                "Tomtom_E":
                    np.nan,

                "STREME_train_p":
                    np.nan,

                "STREME_holdout_p":
                    np.nan,

                "STREME_holdout_E":
                    np.nan
            })


# ------------------------------------------------------------
# TOMTOM DEFAULT STREME -> JASPAR
# ------------------------------------------------------------

for tissue in MOTIF_TISSUES:

    f_annot = os.path.join(
        DEFAULT_TOMTOM_ROOT,
        tissue,
        f"{tissue}.tomtom.annotated.tsv"
    )

    f_raw = os.path.join(
        DEFAULT_TOMTOM_ROOT,
        tissue,
        "tomtom.tsv"
    )

    f = (
        f_annot
        if os.path.exists(f_annot)
        else f_raw
    )

    if not os.path.exists(f):
        continue

    x = read_table(f)

    if x.empty:
        continue

    query_col = find_col(
        x,
        ["Query_ID"]
    )

    target_col = find_col(
        x,
        ["Target_ID"]
    )

    q_col = find_col(
        x,
        ["q-value", "q_value"]
    )

    e_col = find_col(
        x,
        ["E-value", "E_value"],
        required=False
    )

    tf_col = find_col(
        x,
        ["Target_TF"],
        required=False
    )

    x[q_col] = pd.to_numeric(
        x[q_col],
        errors="coerce"
    )

    x = x[
        x[q_col] < 0.05
    ].copy()

    for _, row in x.iterrows():

        query_id = row[query_col]
        target_id = row[target_col]

        if tf_col:
            tf_label = row[tf_col]
        else:
            tf_label = jaspar_map.get(
                target_id,
                ""
            )

        stat = (
            default_streme_stats
            .get(tissue, {})
            .get(query_id, {})
        )

        crossrun = (
            query_id
            in default_repro[tissue]
        )

        for tf in split_tf_label(
            tf_label
        ):

            evidence.append({
                "Tissue":
                    tissue,

                "TF_component":
                    tf,

                "Original_TF_label":
                    tf_label,

                "Evidence_source":
                    "Tomtom_default_STREME",

                "Motif_ID":
                    target_id,

                "STREME_motif_ID":
                    query_id,

                "CrossRun_reproducible_motif":
                    crossrun,

                "AME_p":
                    np.nan,

                "AME_adj_p":
                    np.nan,

                "AME_E":
                    np.nan,

                "AME_OR_HA":
                    np.nan,

                "Tomtom_q":
                    numeric(row[q_col]),

                "Tomtom_E":
                    numeric(
                        row[e_col]
                        if e_col else np.nan
                    ),

                "STREME_train_p":
                    stat.get(
                        "train_pvalue",
                        np.nan
                    ),

                "STREME_holdout_p":
                    stat.get(
                        "test_pvalue",
                        np.nan
                    ),

                "STREME_holdout_E":
                    stat.get(
                        "test_evalue",
                        np.nan
                    )
            })


# ------------------------------------------------------------
# TOMTOM NO-HOLDOUT STREME -> JASPAR
# ------------------------------------------------------------

for tissue in MOTIF_TISSUES:

    f = os.path.join(
        NOHOLD_TOMTOM_ROOT,
        tissue,
        "tomtom.tsv"
    )

    if not os.path.exists(f):
        continue

    x = read_table(f)

    if x.empty:
        continue

    query_col = find_col(
        x,
        ["Query_ID"]
    )

    target_col = find_col(
        x,
        ["Target_ID"]
    )

    q_col = find_col(
        x,
        ["q-value", "q_value"]
    )

    e_col = find_col(
        x,
        ["E-value", "E_value"],
        required=False
    )

    x[q_col] = pd.to_numeric(
        x[q_col],
        errors="coerce"
    )

    x = x[
        x[q_col] < 0.05
    ].copy()

    for _, row in x.iterrows():

        query_id = row[query_col]
        target_id = row[target_col]

        tf_label = jaspar_map.get(
            target_id,
            ""
        )

        stat = (
            nohold_streme_stats
            .get(tissue, {})
            .get(query_id, {})
        )

        crossrun = (
            query_id
            in nohold_repro[tissue]
        )

        for tf in split_tf_label(
            tf_label
        ):

            evidence.append({
                "Tissue":
                    tissue,

                "TF_component":
                    tf,

                "Original_TF_label":
                    tf_label,

                "Evidence_source":
                    "Tomtom_noholdout_STREME",

                "Motif_ID":
                    target_id,

                "STREME_motif_ID":
                    query_id,

                "CrossRun_reproducible_motif":
                    crossrun,

                "AME_p":
                    np.nan,

                "AME_adj_p":
                    np.nan,

                "AME_E":
                    np.nan,

                "AME_OR_HA":
                    np.nan,

                "Tomtom_q":
                    numeric(row[q_col]),

                "Tomtom_E":
                    numeric(
                        row[e_col]
                        if e_col else np.nan
                    ),

                # exploratory only
                "STREME_train_p":
                    stat.get(
                        "train_pvalue",
                        np.nan
                    ),

                "STREME_holdout_p":
                    np.nan,

                "STREME_holdout_E":
                    np.nan
            })


evidence_df = pd.DataFrame(
    evidence
)

if evidence_df.empty:
    sys.exit(
        "ERROR: no candidate TF evidence found."
    )

evidence_df.to_csv(
    os.path.join(
        OUTDIR,
        "08C2_candidate_motif_TF_evidence_long.tsv"
    ),
    sep="\t",
    index=False
)

print(
    "Motif-TF evidence rows:",
    len(evidence_df)
)


# ============================================================
# COLLAPSE MOTIF EVIDENCE TO TISSUE x TF
# ============================================================

collapsed = []

for (tissue, tf), g in evidence_df.groupby(
    ["Tissue", "TF_component"],
    sort=False
):

    sources = sorted(
        set(g["Evidence_source"])
    )

    motif_ids = sorted(
        set(
            x
            for x in g["Motif_ID"]
            if str(x).strip()
        )
    )

    streme_ids = sorted(
        set(
            x
            for x in g["STREME_motif_ID"]
            if str(x).strip()
        )
    )

    has_ame = any(
        g["Evidence_source"].str.startswith(
            "AME"
        )
    )

    has_default_tomtom = any(
        g["Evidence_source"]
        == "Tomtom_default_STREME"
    )

    has_nohold_tomtom = any(
        g["Evidence_source"]
        == "Tomtom_noholdout_STREME"
    )

    crossrun_tomtom = bool(
        (
            g["CrossRun_reproducible_motif"]
            == True
        ).any()
        and
        (
            g["Evidence_source"].str.startswith(
                "Tomtom"
            )
        ).any()
    )

    collapsed.append({

        "Tissue":
            tissue,

        "TF_symbol":
            tf,

        "Evidence_sources":
            ";".join(sources),

        "Motif_IDs":
            ";".join(motif_ids),

        "STREME_motif_IDs":
            ";".join(streme_ids),

        "AME_exploratory_support":
            has_ame,

        "Default_Tomtom_support":
            has_default_tomtom,

        "NoHoldout_Tomtom_support":
            has_nohold_tomtom,

        "CrossRun_reproducible_Tomtom_support":
            crossrun_tomtom,

        "Best_AME_p":
            pd.to_numeric(
                g["AME_p"],
                errors="coerce"
            ).min(),

        "Best_AME_adj_p":
            pd.to_numeric(
                g["AME_adj_p"],
                errors="coerce"
            ).min(),

        "Best_AME_E":
            pd.to_numeric(
                g["AME_E"],
                errors="coerce"
            ).min(),

        "Max_AME_OR_HA":
            pd.to_numeric(
                g["AME_OR_HA"],
                errors="coerce"
            ).max(),

        "Best_Tomtom_q":
            pd.to_numeric(
                g["Tomtom_q"],
                errors="coerce"
            ).min(),

        "Best_Tomtom_E":
            pd.to_numeric(
                g["Tomtom_E"],
                errors="coerce"
            ).min(),

        "Best_default_STREME_holdout_E":
            pd.to_numeric(
                g["STREME_holdout_E"],
                errors="coerce"
            ).min(),

        "Best_STREME_train_p":
            pd.to_numeric(
                g["STREME_train_p"],
                errors="coerce"
            ).min()
    })


candidate_df = pd.DataFrame(
    collapsed
)


# ============================================================
# GENE SYMBOL -> ENSEMBL GENE ID
# ============================================================

symbol_to_gene = defaultdict(set)
gene_to_symbol = {}


gtf_candidates = [
    "Sus_longest.gtf",
    "Sus_longest.nochr.gtf"
]

gtf_file = None

for f in gtf_candidates:
    if os.path.exists(f):
        gtf_file = f
        break


if gtf_file:

    print(
        "Gene-name mapping GTF:",
        gtf_file
    )

    with open(
        gtf_file,
        "r",
        errors="ignore"
    ) as fh:

        for line in fh:

            if (
                line.startswith("#")
                or "\t" not in line
            ):
                continue

            parts = line.rstrip(
                "\n"
            ).split("\t")

            if len(parts) < 9:
                continue

            attrs = parts[8]

            gid_m = re.search(
                r'gene_id "([^"]+)"',
                attrs
            )

            name_m = re.search(
                r'gene_name "([^"]+)"',
                attrs
            )

            if not name_m:
                name_m = re.search(
                    r'gene_symbol "([^"]+)"',
                    attrs
                )

            if (
                gid_m
                and name_m
            ):

                gid = canonical_gene_id(
                    gid_m.group(1)
                )

                name = name_m.group(1)

                symbol_to_gene[
                    name.upper()
                ].add(gid)

                gene_to_symbol[
                    gid
                ] = name


# ------------------------------------------------------------
# Optional UROPA fallback
# ------------------------------------------------------------

UROPA_FILE = (
    "pig_8tissues_master_peaks_finalhits.txt"
)

if os.path.exists(UROPA_FILE):

    try:

        uropa = pd.read_csv(
            UROPA_FILE,
            sep="\t",
            usecols=[
                "gene_id",
                "name"
            ],
            dtype=str
        )

        uropa = uropa.dropna()

        for _, row in uropa.iterrows():

            gid = canonical_gene_id(
                row["gene_id"]
            )

            name = str(
                row["name"]
            ).strip()

            if (
                name
                and name.upper() != "NA"
            ):

                symbol_to_gene[
                    name.upper()
                ].add(gid)

                gene_to_symbol.setdefault(
                    gid,
                    name
                )

    except Exception as e:

        print(
            "WARNING: UROPA fallback not used:",
            e
        )


print(
    "Unique mapped gene symbols:",
    len(symbol_to_gene)
)


# ============================================================
# READ RNA TPM
# ============================================================

print()
print("Reading RNA TMM.TPM...")

rna = pd.read_csv(
    TPM_FILE,
    sep="\t",
    index_col=0
)

rna.index = [
    canonical_gene_id(x)
    for x in rna.index
]

rna = rna.apply(
    pd.to_numeric,
    errors="coerce"
).fillna(0)


if rna.index.duplicated().any():

    print(
        "WARNING: duplicated canonical gene IDs; summing."
    )

    rna = rna.groupby(
        level=0
    ).sum()


print(
    "RNA dimensions:",
    rna.shape[0],
    "genes x",
    rna.shape[1],
    "samples"
)


# ============================================================
# METADATA
# ============================================================

meta = pd.read_csv(
    META_FILE,
    sep="\t",
    dtype=str
)

required_meta = [
    "Tissue",
    "Animal",
    "RNA_SRR"
]

for c in required_meta:

    if c not in meta.columns:
        sys.exit(
            f"ERROR: metadata missing column {c}"
        )


missing_samples = (
    set(meta["RNA_SRR"])
    - set(rna.columns)
)

if missing_samples:

    sys.exit(
        "ERROR: RNA samples missing from TPM matrix:\n"
        + "\n".join(
            sorted(missing_samples)
        )
    )


# ============================================================
# BUILD 8-TISSUE MEAN MATRIX
# ============================================================

tissue_mean = pd.DataFrame(
    index=rna.index
)

animal_tpm = {}

for tissue in TISSUES:

    sub = meta[
        meta["Tissue"] == tissue
    ]

    samples = sub["RNA_SRR"].tolist()

    if len(samples) == 0:
        sys.exit(
            f"ERROR: no RNA samples for {tissue}"
        )

    tissue_mean[tissue] = (
        rna[samples]
        .mean(axis=1)
    )

    for _, row in sub.iterrows():

        animal_tpm[
            (
                tissue,
                row["Animal"]
            )
        ] = rna[
            row["RNA_SRR"]
        ]


# ============================================================
# TAU
# ============================================================

max_expr = tissue_mean.max(
    axis=1
)

relative = tissue_mean.div(
    max_expr.replace(
        0,
        np.nan
    ),
    axis=0
)

tau = (
    (1 - relative)
    .sum(axis=1)
    / (len(TISSUES) - 1)
)

tau[max_expr == 0] = 0


# ============================================================
# SPM
# ============================================================

l2_norm = np.sqrt(
    (tissue_mean ** 2)
    .sum(axis=1)
)

spm = tissue_mean.div(
    l2_norm.replace(
        0,
        np.nan
    ),
    axis=0
)


# ============================================================
# EXPRESSION RANK
# ============================================================

tissue_rank = tissue_mean.rank(
    axis=1,
    ascending=False,
    method="min"
)

max_tissue = tissue_mean.idxmax(
    axis=1
)


# ============================================================
# MAP CANDIDATE TF SYMBOLS TO RNA GENES
# ============================================================

integrated_rows = []
unmapped_rows = []

for _, row in candidate_df.iterrows():

    tissue = row["Tissue"]
    tf = row["TF_symbol"]

    gene_ids = sorted(
        symbol_to_gene.get(
            str(tf).upper(),
            []
        )
    )

    # retain only IDs present in RNA
    gene_ids_in_rna = [
        gid
        for gid in gene_ids
        if gid in rna.index
    ]

    if not gene_ids_in_rna:

        unmapped = row.to_dict()

        unmapped[
            "Reason"
        ] = (
            "No_gene_symbol_to_RNA_gene_mapping"
        )

        unmapped_rows.append(
            unmapped
        )

        continue


    for gid in gene_ids_in_rna:

        target_mean = tissue_mean.loc[
            gid,
            tissue
        ]

        p348 = (
            animal_tpm[
                (
                    tissue,
                    "P348"
                )
            ].loc[gid]
            if (
                tissue,
                "P348"
            ) in animal_tpm
            else np.nan
        )

        p350 = (
            animal_tpm[
                (
                    tissue,
                    "P350"
                )
            ].loc[gid]
            if (
                tissue,
                "P350"
            ) in animal_tpm
            else np.nan
        )

        both_detected = (
            p348 >= ANIMAL_TPM_DETECTED
            and
            p350 >= ANIMAL_TPM_DETECTED
        )

        mean_ge1 = (
            target_mean >= MEAN_TPM_HIGH
        )

        high_rna_support = (
            both_detected
            and mean_ge1
        )

        target_is_max = (
            max_tissue.loc[gid]
            == tissue
        )

        target_spm = spm.loc[
            gid,
            tissue
        ]

        target_rank = tissue_rank.loc[
            gid,
            tissue
        ]

        tau_value = tau.loc[
            gid
        ]


        # ----------------------------------------------------
        # Exploratory evidence priority
        #
        # IMPORTANT:
        # These are prioritization categories,
        # NOT formal statistical significance tiers.
        # ----------------------------------------------------

        if (
            bool(
                row[
                    "CrossRun_reproducible_Tomtom_support"
                ]
            )
            and high_rna_support
        ):

            priority = (
                "P1_reproducible_motif_Tomtom_plus_RNA"
            )

        elif (
            (
                bool(
                    row[
                        "Default_Tomtom_support"
                    ]
                )
                or
                bool(
                    row[
                        "NoHoldout_Tomtom_support"
                    ]
                )
            )
            and high_rna_support
        ):

            priority = (
                "P2_Tomtom_plus_RNA"
            )

        elif (
            bool(
                row[
                    "AME_exploratory_support"
                ]
            )
            and high_rna_support
        ):

            priority = (
                "P3_AME_exploratory_plus_RNA"
            )

        else:

            priority = (
                "P4_motif_candidate_low_RNA_support"
            )


        out = row.to_dict()

        out.update({

            "gene_id":
                gid,

            "GTF_gene_name":
                gene_to_symbol.get(
                    gid,
                    tf
                ),

            "Gene_symbol_mapping_count":
                len(gene_ids_in_rna),

            "Target_tissue_mean_TPM":
                target_mean,

            "Target_tissue_P348_TPM":
                p348,

            "Target_tissue_P350_TPM":
                p350,

            "Both_animals_TPM_ge_0.5":
                both_detected,

            "Target_mean_TPM_ge_1":
                mean_ge1,

            "High_RNA_expression_support":
                high_rna_support,

            "Tau":
                tau_value,

            "Target_tissue_SPM":
                target_spm,

            "Target_tissue_expression_rank":
                target_rank,

            "Max_expression_tissue":
                max_tissue.loc[gid],

            "Target_is_max_expression_tissue":
                target_is_max,

            "SPM_ge_0.5":
                (
                    target_spm >= 0.5
                    if not pd.isna(
                        target_spm
                    )
                    else False
                ),

            "Tissue_preference_support":
                (
                    target_is_max
                    and
                    not pd.isna(
                        target_spm
                    )
                    and
                    target_spm >= 0.5
                ),

            "Exploratory_priority":
                priority
        })


        # add all eight tissue means
        for t in TISSUES:

            out[
                f"Mean_TPM_{t}"
            ] = tissue_mean.loc[
                gid,
                t
            ]

            out[
                f"SPM_{t}"
            ] = spm.loc[
                gid,
                t
            ]


        integrated_rows.append(
            out
        )


integrated = pd.DataFrame(
    integrated_rows
)

unmapped = pd.DataFrame(
    unmapped_rows
)


# ============================================================
# PRIORITY ORDER
# ============================================================

priority_order = {
    "P1_reproducible_motif_Tomtom_plus_RNA": 1,
    "P2_Tomtom_plus_RNA": 2,
    "P3_AME_exploratory_plus_RNA": 3,
    "P4_motif_candidate_low_RNA_support": 4
}

if not integrated.empty:

    integrated[
        "_priority_order"
    ] = integrated[
        "Exploratory_priority"
    ].map(priority_order)

    integrated = integrated.sort_values(
        [
            "Tissue",
            "_priority_order",
            "Target_tissue_mean_TPM",
            "Target_tissue_SPM"
        ],
        ascending=[
            True,
            True,
            False,
            False
        ]
    )

    integrated = integrated.drop(
        columns=[
            "_priority_order"
        ]
    )


# ============================================================
# WRITE OUTPUTS
# ============================================================

integrated_file = os.path.join(
    OUTDIR,
    "08C3_TF_RNA_Tau_SPM_integration_all.tsv"
)

integrated.to_csv(
    integrated_file,
    sep="\t",
    index=False
)


unmapped_file = os.path.join(
    OUTDIR,
    "08C4_unmapped_TF_candidates.tsv"
)

unmapped.to_csv(
    unmapped_file,
    sep="\t",
    index=False
)


if not integrated.empty:

    priority = integrated[
        integrated[
            "Exploratory_priority"
        ].isin(
            [
                "P1_reproducible_motif_Tomtom_plus_RNA",
                "P2_Tomtom_plus_RNA",
                "P3_AME_exploratory_plus_RNA"
            ]
        )
    ].copy()

else:

    priority = integrated.copy()


priority_file = os.path.join(
    OUTDIR,
    "08C5_priority_TF_candidates.tsv"
)

priority.to_csv(
    priority_file,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Candidate TF tissue-mean TPM matrix
# ------------------------------------------------------------

if not integrated.empty:

    candidate_genes = (
        integrated[
            [
                "gene_id",
                "GTF_gene_name"
            ]
        ]
        .drop_duplicates()
    )

    matrix_rows = []

    for _, r in candidate_genes.iterrows():

        gid = r["gene_id"]

        x = {
            "gene_id":
                gid,
            "gene_name":
                r["GTF_gene_name"],
            "Tau":
                tau.loc[gid],
            "Max_expression_tissue":
                max_tissue.loc[gid]
        }

        for t in TISSUES:
            x[t] = tissue_mean.loc[
                gid,
                t
            ]

        matrix_rows.append(x)

    pd.DataFrame(
        matrix_rows
    ).to_csv(
        os.path.join(
            OUTDIR,
            "08C6_candidate_TF_tissue_mean_TPM.tsv"
        ),
        sep="\t",
        index=False
    )


# ============================================================
# QC SUMMARY
# ============================================================

qc_rows = []

for tissue in MOTIF_TISSUES:

    sub_e = candidate_df[
        candidate_df[
            "Tissue"
        ] == tissue
    ]

    sub_i = integrated[
        integrated[
            "Tissue"
        ] == tissue
    ] if not integrated.empty else integrated

    sub_u = unmapped[
        unmapped[
            "Tissue"
        ] == tissue
    ] if not unmapped.empty else unmapped

    qc_rows.append({

        "Tissue":
            tissue,

        "Candidate_TF_symbols":
            sub_e[
                "TF_symbol"
            ].nunique(),

        "Mapped_TF_gene_rows":
            len(sub_i),

        "Mapped_unique_TF_symbols":
            sub_i[
                "TF_symbol"
            ].nunique()
            if not sub_i.empty
            else 0,

        "Unmapped_TF_symbols":
            sub_u[
                "TF_symbol"
            ].nunique()
            if not sub_u.empty
            else 0,

        "High_RNA_support":
            int(
                sub_i[
                    "High_RNA_expression_support"
                ].sum()
            )
            if not sub_i.empty
            else 0,

        "P1_reproducible_Tomtom_RNA":
            int(
                (
                    sub_i[
                        "Exploratory_priority"
                    ]
                    ==
                    "P1_reproducible_motif_Tomtom_plus_RNA"
                ).sum()
            )
            if not sub_i.empty
            else 0,

        "P2_Tomtom_RNA":
            int(
                (
                    sub_i[
                        "Exploratory_priority"
                    ]
                    ==
                    "P2_Tomtom_plus_RNA"
                ).sum()
            )
            if not sub_i.empty
            else 0,

        "P3_AME_RNA":
            int(
                (
                    sub_i[
                        "Exploratory_priority"
                    ]
                    ==
                    "P3_AME_exploratory_plus_RNA"
                ).sum()
            )
            if not sub_i.empty
            else 0,

        "Target_is_max_expression_tissue":
            int(
                sub_i[
                    "Target_is_max_expression_tissue"
                ].sum()
            )
            if not sub_i.empty
            else 0,

        "SPM_ge_0.5":
            int(
                sub_i[
                    "SPM_ge_0.5"
                ].sum()
            )
            if not sub_i.empty
            else 0
    })


qc = pd.DataFrame(
    qc_rows
)

qc_file = os.path.join(
    OUTDIR,
    "08C7_QC_summary.tsv"
)

qc.to_csv(
    qc_file,
    sep="\t",
    index=False
)


# ============================================================
# PER-TISSUE PRIORITY FILES
# ============================================================

if not integrated.empty:

    for tissue in MOTIF_TISSUES:

        z = integrated[
            integrated[
                "Tissue"
            ] == tissue
        ].copy()

        z.to_csv(
            os.path.join(
                OUTDIR,
                f"08C_{tissue}_TF_candidates.tsv"
            ),
            sep="\t",
            index=False
        )


# ============================================================
# LOG SUMMARY
# ============================================================

print()
print("=" * 72)
print("STEP 08C COMPLETED")
print("=" * 72)

print()
print(qc.to_string(index=False))

print()
print("IMPORTANT INTERPRETATION:")
print(
    "P1-P4 are exploratory prioritization categories, "
    "NOT statistical significance tiers."
)
print(
    "Tau/SPM describe tissue preference; "
    "they are NOT mandatory filters for TF plausibility."
)
print(
    "AME adj_p<0.05 is treated as exploratory because "
    "AME E<0.05 was not achieved."
)
print(
    "Default STREME hold-out E-values remain the formal "
    "de novo validation statistic."
)

print()
print("Main outputs:")
print(integrated_file)
print(priority_file)
print(unmapped_file)
print(qc_file)
print("=" * 72)

