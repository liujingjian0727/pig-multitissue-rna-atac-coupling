#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import math
import xml.etree.ElementTree as ET
from collections import defaultdict

# ============================================================
# STEP 08D
# Publication-ready motif-family / TF-candidate integration
#
# Inputs:
#   Step 08B1c AME
#   Step 08B2 default STREME
#   Step 08B3 default Tomtom
#   Step 08B4B STREME cross-run reproducibility
#   Step 08C TF RNA/Tau/SPM integration
#
# Main principles:
#   1. AME E<0.05 remains the formal known-motif criterion.
#   2. Default STREME hold-out E<0.05 remains the formal
#      de novo motif criterion.
#   3. AME adj_p<0.05 is exploratory only.
#   4. no-holdout STREME is sensitivity/reproducibility only.
#   5. FIMO priority is NOT a statistical-significance tier.
#   6. FIMO scans motif occurrences; it does not establish
#      TF binding or TF->target causality.
# ============================================================

TISSUES = ["Muscle", "Spleen", "Liver", "Cerebellum"]

ROOT = Path(".")
OUT = ROOT / "08D_motif_family_TF_master"
FIMO_ROOT = OUT / "FIMO_motif_sets"

OUT.mkdir(exist_ok=True)
FIMO_ROOT.mkdir(exist_ok=True)

RNA_ALL = Path(
    "08C_TF_RNA_integration/"
    "08C3_TF_RNA_Tau_SPM_integration_all.tsv"
)

RNA_UNMAPPED = Path(
    "08C_TF_RNA_integration/"
    "08C4_unmapped_TF_candidates.tsv"
)

JASPAR_MAP = Path("JASPAR2024_motif_ID_to_TF.tsv")

JASPAR_DB = Path(
    "/datadisk2/liujingjian_data/ATAC-seq/Sus/"
    "footprints/footprinting/motif_databases/JASPAR/"
    "JASPAR2024_CORE_vertebrates_non-redundant_v2.meme"
)

# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n"
        )
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in columns})


def as_bool(x):
    if x is None:
        return False
    return str(x).strip().lower() in {
        "true", "t", "1", "yes", "y"
    }


def as_float(x):
    if x is None:
        return math.nan
    x = str(x).strip()
    if x == "":
        return math.nan
    try:
        return float(x)
    except Exception:
        return math.nan


def finite(x):
    return not math.isnan(x)


def split_multi(x):
    if x is None:
        return []
    x = str(x).strip()
    if not x:
        return []
    return [
        z.strip()
        for z in re.split(r"[;,]", x)
        if z.strip()
    ]


def uniq_sorted(x):
    return sorted(set(z for z in x if z))


def join_values(x):
    return ";".join(uniq_sorted(x))


def fmt_num(x):
    if isinstance(x, float):
        if math.isnan(x):
            return ""
        return f"{x:.8g}"
    return str(x)


def minimum_numeric(values):
    vals = [
        as_float(x)
        for x in values
        if finite(as_float(x))
    ]
    return min(vals) if vals else math.nan


def maximum_numeric(values):
    vals = [
        as_float(x)
        for x in values
        if finite(as_float(x))
    ]
    return max(vals) if vals else math.nan


def tf_preference_support(row):
    return (
        as_bool(row.get("Target_is_max_expression_tissue"))
        or
        as_bool(row.get("SPM_ge_0.5"))
    )


def high_rna(row):
    return as_bool(row.get("High_RNA_expression_support"))


def sanitize_id(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


# ------------------------------------------------------------
# JASPAR ID -> TF/motif name
# ------------------------------------------------------------

jaspar_name = {}

with JASPAR_MAP.open() as fh:
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            jaspar_name[parts[0]] = parts[1]

print(f"JASPAR mappings loaded: {len(jaspar_name)}")

# ------------------------------------------------------------
# Step 08C evidence
# ------------------------------------------------------------

mapped_rows = read_tsv(RNA_ALL)
unmapped_rows = read_tsv(RNA_UNMAPPED)

# Normalize unmapped rows so later code can treat all evidence similarly.
for r in unmapped_rows:
    r.setdefault("gene_id", "")
    r.setdefault("GTF_gene_name", "")
    r.setdefault("High_RNA_expression_support", "False")
    r.setdefault("Target_is_max_expression_tissue", "False")
    r.setdefault("SPM_ge_0.5", "False")
    r.setdefault("Target_tissue_mean_TPM", "")
    r.setdefault("Target_tissue_SPM", "")
    r.setdefault("Tau", "")
    r.setdefault("Exploratory_priority", "Unmapped")

evidence_rows = mapped_rows + unmapped_rows

print(f"Mapped Step08C rows  : {len(mapped_rows)}")
print(f"Unmapped Step08C rows: {len(unmapped_rows)}")

# ------------------------------------------------------------
# Parse default STREME XML
# ------------------------------------------------------------

def parse_streme_xml(path):
    root = ET.parse(path).getroot()

    model = root.find("model")
    test_name = ""
    if model is not None:
        test_el = model.find("test")
        if test_el is not None and test_el.text:
            test_name = test_el.text.strip()

    result = []

    motifs_el = root.find("motifs")
    if motifs_el is None:
        return result

    for m in motifs_el.findall("motif"):
        a = dict(m.attrib)

        motif_id = a.get("id", "")
        consensus = motif_id.split("-", 1)[1] if "-" in motif_id else motif_id

        result.append({
            "STREME_motif_ID": motif_id,
            "STREME_alt": a.get("alt", ""),
            "Consensus": consensus,
            "Width": a.get("width", ""),
            "Train_p": as_float(a.get("train_pvalue")),
            "Holdout_p": as_float(a.get("test_pvalue")),
            "Holdout_E": as_float(a.get("test_evalue")),
            "Train_pos": a.get("train_pos_count", ""),
            "Train_neg": a.get("train_neg_count", ""),
            "Test_pos": a.get("test_pos_count", ""),
            "Test_neg": a.get("test_neg_count", ""),
            "Test_method": test_name,
            "Npassing": a.get("npassing", ""),
        })

    return result


streme_stats = {}

for tissue in TISSUES:
    xml = Path(
        f"08B2_STREME_de_novo/{tissue}/streme.xml"
    )
    motifs = parse_streme_xml(xml)

    streme_stats[tissue] = {
        m["STREME_motif_ID"]: m
        for m in motifs
    }

    print(
        f"{tissue:12s} default STREME motifs: {len(motifs)}"
    )

# ------------------------------------------------------------
# Parse default Tomtom JASPAR annotation
# ------------------------------------------------------------

default_tomtom = {
    t: defaultdict(list)
    for t in TISSUES
}

for tissue in TISSUES:

    f = Path(
        f"08B3_Tomtom_JASPAR2024/{tissue}/"
        f"{tissue}.tomtom.annotated.tsv"
    )

    rows = read_tsv(f)

    for r in rows:
        qid = r.get("Query_ID", "").strip()

        if not qid:
            continue

        q = as_float(r.get("q-value"))

        if not finite(q) or q >= 0.05:
            continue

        default_tomtom[tissue][qid].append({
            "Target_ID": r.get("Target_ID", ""),
            "Target_TF": r.get("Target_TF", ""),
            "q": q,
            "E": as_float(r.get("E-value")),
            "p": as_float(r.get("p-value")),
            "Overlap": r.get("Overlap", ""),
        })

# ------------------------------------------------------------
# Parse default-vs-noholdout STREME reproducibility
# ------------------------------------------------------------

crossrun_reproduced = {
    t: set()
    for t in TISSUES
}

crossrun_best_q = {
    t: {}
    for t in TISSUES
}

for tissue in TISSUES:

    f = Path(
        f"08B4B_STREME_reproducibility/"
        f"{tissue}/tomtom.tsv"
    )

    rows = read_tsv(f)

    for r in rows:
        qid = r.get("Query_ID", "").strip()
        if not qid:
            continue

        q = as_float(r.get("q-value"))

        if finite(q) and q < 0.05:
            crossrun_reproduced[tissue].add(qid)

            old = crossrun_best_q[tissue].get(
                qid,
                math.inf
            )
            crossrun_best_q[tissue][qid] = min(old, q)

# ------------------------------------------------------------
# Family label heuristic
#
# Important:
# This is used only to collapse redundant TF identities into
# biologically interpretable motif-family labels.
# It does NOT assign exact TF occupancy.
# ------------------------------------------------------------

def infer_denovo_family(consensus, tf_symbols):

    tfs = [x.upper() for x in tf_symbols]

    gc = 0
    seq = consensus.upper()

    if seq:
        gc = sum(x in {"G", "C", "S"} for x in seq) / len(seq)

    has_klfsp = any(
        re.match(r"^(KLF\d+|SP\d+)$", x)
        for x in tfs
    )

    has_znf_gc = any(
        x.startswith("ZNF")
        or x in {
            "MAZ", "PATZ1", "ZBED4", "VEZF1",
            "ZKSCAN3", "RREB1"
        }
        for x in tfs
    )

    has_ets = any(
        re.match(
            r"^(ELF|ETS|ETV|ERG|FLI|ELK|GABP)",
            x
        )
        for x in tfs
    )

    has_ap1 = any(
        x.startswith("JUN")
        or x.startswith("FOS")
        or x in {"BATF", "ATF3"}
        for x in tfs
    )

    if has_ets:
        return "ETS-like motif"

    if has_klfsp:
        if gc >= 0.60:
            return "GC-rich KLF/SP-like motif"
        return "KLF/SP-like motif"

    if has_znf_gc:
        if gc >= 0.60:
            return "GC-rich zinc-finger-like motif"
        return "zinc-finger-like motif"

    if has_ap1:
        return "AP-1-like motif"

    if tf_symbols:
        return "/".join(tf_symbols[:3]) + "-like motif"

    return "Unannotated de novo motif"


# ------------------------------------------------------------
# Index Step08C TF rows
# ------------------------------------------------------------

by_tissue_tf = {}

for r in evidence_rows:
    tissue = r.get("Tissue", "")
    tf = r.get("TF_symbol", "")

    if tissue and tf:
        by_tissue_tf[(tissue, tf)] = r


# ============================================================
# PART A. Known JASPAR motif families
# ============================================================

known_groups = defaultdict(list)

for r in evidence_rows:

    if not as_bool(r.get("AME_exploratory_support")):
        continue

    motif_ids = split_multi(r.get("Motif_IDs", ""))

    for motif_id in motif_ids:

        if not motif_id.startswith("MA"):
            continue

        known_groups[
            (r.get("Tissue", ""), motif_id)
        ].append(r)


motif_family_rows = []

for (tissue, motif_id), rows in sorted(known_groups.items()):

    all_tfs = uniq_sorted([
        r.get("TF_symbol", "")
        for r in rows
    ])

    mapped = [
        r for r in rows
        if r.get("gene_id", "").strip()
    ]

    high_rows = [
        r for r in mapped
        if high_rna(r)
    ]

    preferred_rows = [
        r for r in high_rows
        if tf_preference_support(r)
    ]

    high_tfs = uniq_sorted([
        r.get("TF_symbol", "")
        for r in high_rows
    ])

    preferred_tfs = uniq_sorted([
        r.get("TF_symbol", "")
        for r in preferred_rows
    ])

    ame_p = minimum_numeric(
        r.get("Best_AME_p") for r in rows
    )

    ame_adj = minimum_numeric(
        r.get("Best_AME_adj_p") for r in rows
    )

    ame_e = minimum_numeric(
        r.get("Best_AME_E") for r in rows
    )

    ame_or = maximum_numeric(
        r.get("Max_AME_OR_HA") for r in rows
    )

    formal_sig = finite(ame_e) and ame_e < 0.05

    if preferred_rows:
        tier = "Primary"
        include = True
        reason = (
            "AME exploratory motif + high RNA + "
            "target-tissue preference"
        )

    elif high_rows:
        tier = "Secondary"
        include = True
        reason = (
            "AME exploratory motif + high RNA; "
            "no strong target-tissue preference"
        )

    else:
        tier = "Not_selected"
        include = False
        reason = (
            "AME exploratory motif without sufficient "
            "RNA-supported TF evidence"
        )

    name = jaspar_name.get(motif_id, motif_id)

    motif_family_rows.append({
        "Tissue": tissue,
        "Motif_family_ID":
            f"{tissue}.KNOWN.{motif_id}",
        "Evidence_track": "Known_JASPAR_AME",
        "Motif_family_label":
            f"{name} known motif",
        "Representative_motif_ID": motif_id,
        "Representative_motif_name": name,
        "Representative_source": "JASPAR2024",
        "Consensus": "",
        "Candidate_TFs": join_values(all_tfs),
        "High_RNA_candidate_TFs":
            join_values(high_tfs),
        "Tissue_preferred_candidate_TFs":
            join_values(preferred_tfs),
        "Candidate_TF_N": len(all_tfs),
        "High_RNA_candidate_TF_N": len(high_tfs),
        "Tissue_preferred_candidate_TF_N":
            len(preferred_tfs),
        "AME_exploratory_support": True,
        "AME_best_p": fmt_num(ame_p),
        "AME_best_adj_p": fmt_num(ame_adj),
        "AME_best_E": fmt_num(ame_e),
        "AME_max_OR_HA": fmt_num(ame_or),
        "AME_formal_E_lt_0.05": formal_sig,
        "Default_STREME_motif_ID": "",
        "Default_STREME_train_p": "",
        "Default_STREME_holdout_p": "",
        "Default_STREME_holdout_E": "",
        "Default_STREME_formal_E_lt_0.05": "",
        "Default_Tomtom_q_lt_0.05": "",
        "Best_default_Tomtom_q": "",
        "CrossRun_reproduced": "",
        "CrossRun_best_q": "",
        "FIMO_tier": tier,
        "FIMO_include": include,
        "FIMO_rationale": reason,
        "Interpretation":
            "Known motif candidate; exact TF binding "
            "is not established"
    })


# ============================================================
# PART B. Default STREME de novo motif families
# ============================================================

for tissue in TISSUES:

    for motif_id, stat in streme_stats[tissue].items():

        hits = default_tomtom[tissue].get(
            motif_id,
            []
        )

        hit_tfs = uniq_sorted([
            h["Target_TF"]
            for h in hits
            if h["Target_TF"]
        ])

        candidate_rows = []

        for tf in hit_tfs:
            r = by_tissue_tf.get((tissue, tf))
            if r is not None:
                candidate_rows.append(r)

        high_rows = [
            r for r in candidate_rows
            if high_rna(r)
        ]

        preferred_rows = [
            r for r in high_rows
            if tf_preference_support(r)
        ]

        high_tfs = uniq_sorted([
            r.get("TF_symbol", "")
            for r in high_rows
        ])

        preferred_tfs = uniq_sorted([
            r.get("TF_symbol", "")
            for r in preferred_rows
        ])

        best_q = minimum_numeric(
            h["q"] for h in hits
        )

        reproduced = (
            motif_id in crossrun_reproduced[tissue]
        )

        cr_q = crossrun_best_q[
            tissue
        ].get(
            motif_id,
            math.nan
        )

        has_tomtom = len(hits) > 0

        holdout_e = stat["Holdout_E"]

        formal_sig = (
            finite(holdout_e)
            and holdout_e < 0.05
        )

        if (
            has_tomtom
            and reproduced
            and len(high_rows) > 0
        ):
            tier = "Primary"
            include = True
            reason = (
                "Default STREME motif with JASPAR "
                "Tomtom support, cross-run motif "
                "reproducibility, and high-RNA TF candidate"
            )

        elif (
            has_tomtom
            and len(high_rows) > 0
        ):
            tier = "Secondary"
            include = True
            reason = (
                "Default STREME motif with JASPAR "
                "Tomtom + high-RNA TF candidate, "
                "but no cross-run motif reproducibility"
            )

        else:
            tier = "Not_selected"
            include = False

            if not has_tomtom and reproduced:
                reason = (
                    "Reproducible de novo motif but "
                    "no significant JASPAR Tomtom annotation"
                )
            elif not has_tomtom:
                reason = (
                    "No significant JASPAR Tomtom annotation"
                )
            elif not high_rows:
                reason = (
                    "Tomtom annotation present but "
                    "no high-RNA TF candidate"
                )
            else:
                reason = "Insufficient integrated support"

        label = infer_denovo_family(
            stat["Consensus"],
            hit_tfs
        )

        motif_family_rows.append({
            "Tissue": tissue,
            "Motif_family_ID":
                f"{tissue}.DENOVO.{sanitize_id(motif_id)}",
            "Evidence_track": "De_novo_STREME",
            "Motif_family_label": label,
            "Representative_motif_ID": motif_id,
            "Representative_motif_name":
                stat["STREME_alt"],
            "Representative_source":
                "STREME_default",
            "Consensus": stat["Consensus"],
            "Candidate_TFs": join_values(hit_tfs),
            "High_RNA_candidate_TFs":
                join_values(high_tfs),
            "Tissue_preferred_candidate_TFs":
                join_values(preferred_tfs),
            "Candidate_TF_N": len(hit_tfs),
            "High_RNA_candidate_TF_N":
                len(high_tfs),
            "Tissue_preferred_candidate_TF_N":
                len(preferred_tfs),
            "AME_exploratory_support": "",
            "AME_best_p": "",
            "AME_best_adj_p": "",
            "AME_best_E": "",
            "AME_max_OR_HA": "",
            "AME_formal_E_lt_0.05": "",
            "Default_STREME_motif_ID":
                motif_id,
            "Default_STREME_train_p":
                fmt_num(stat["Train_p"]),
            "Default_STREME_holdout_p":
                fmt_num(stat["Holdout_p"]),
            "Default_STREME_holdout_E":
                fmt_num(stat["Holdout_E"]),
            "Default_STREME_formal_E_lt_0.05":
                formal_sig,
            "Default_Tomtom_q_lt_0.05":
                has_tomtom,
            "Best_default_Tomtom_q":
                fmt_num(best_q),
            "CrossRun_reproduced":
                reproduced,
            "CrossRun_best_q":
                fmt_num(cr_q),
            "FIMO_tier": tier,
            "FIMO_include": include,
            "FIMO_rationale": reason,
            "Interpretation":
                "De novo motif family; Tomtom TF "
                "matches indicate motif similarity, "
                "not exact TF occupancy"
        })


# ------------------------------------------------------------
# Sort motif-family master
# ------------------------------------------------------------

tier_order = {
    "Primary": 1,
    "Secondary": 2,
    "Not_selected": 3
}

track_order = {
    "Known_JASPAR_AME": 1,
    "De_novo_STREME": 2
}

tissue_order = {
    x: i
    for i, x in enumerate(TISSUES, start=1)
}

motif_family_rows.sort(
    key=lambda r: (
        tissue_order.get(r["Tissue"], 99),
        tier_order.get(r["FIMO_tier"], 99),
        track_order.get(r["Evidence_track"], 99),
        r["Representative_motif_ID"]
    )
)

family_columns = [
    "Tissue",
    "Motif_family_ID",
    "Evidence_track",
    "Motif_family_label",
    "Representative_motif_ID",
    "Representative_motif_name",
    "Representative_source",
    "Consensus",
    "Candidate_TFs",
    "High_RNA_candidate_TFs",
    "Tissue_preferred_candidate_TFs",
    "Candidate_TF_N",
    "High_RNA_candidate_TF_N",
    "Tissue_preferred_candidate_TF_N",
    "AME_exploratory_support",
    "AME_best_p",
    "AME_best_adj_p",
    "AME_best_E",
    "AME_max_OR_HA",
    "AME_formal_E_lt_0.05",
    "Default_STREME_motif_ID",
    "Default_STREME_train_p",
    "Default_STREME_holdout_p",
    "Default_STREME_holdout_E",
    "Default_STREME_formal_E_lt_0.05",
    "Default_Tomtom_q_lt_0.05",
    "Best_default_Tomtom_q",
    "CrossRun_reproduced",
    "CrossRun_best_q",
    "FIMO_tier",
    "FIMO_include",
    "FIMO_rationale",
    "Interpretation"
]

write_tsv(
    OUT / "08D1_publication_motif_family_master.tsv",
    motif_family_rows,
    family_columns
)


# ============================================================
# PART C. TF candidate master
# ============================================================

family_by_tissue_tf = defaultdict(list)

for fam in motif_family_rows:

    for tf in split_multi(fam["Candidate_TFs"]):

        family_by_tissue_tf[
            (fam["Tissue"], tf)
        ].append(fam)


tf_master = []

for r in evidence_rows:

    tissue = r.get("Tissue", "")
    tf = r.get("TF_symbol", "")

    if not tissue or not tf:
        continue

    fams = family_by_tissue_tf.get(
        (tissue, tf),
        []
    )

    if not fams:
        continue

    tiers = [
        f["FIMO_tier"]
        for f in fams
    ]

    if "Primary" in tiers:
        best_tier = "Primary"
    elif "Secondary" in tiers:
        best_tier = "Secondary"
    else:
        best_tier = "Not_selected"

    family_ids = [
        f["Motif_family_ID"]
        for f in fams
    ]

    family_labels = [
        f["Motif_family_label"]
        for f in fams
    ]

    rep_motifs = [
        f["Representative_motif_ID"]
        for f in fams
    ]

    if high_rna(r) and tf_preference_support(r):
        rna_role = "High_RNA_plus_tissue_preference"
    elif high_rna(r):
        rna_role = "High_RNA"
    elif r.get("gene_id", "").strip():
        rna_role = "Low_RNA_support"
    else:
        rna_role = "Unmapped_to_RNA_gene"

    tf_master.append({
        "Tissue": tissue,
        "TF_symbol": tf,
        "gene_id": r.get("gene_id", ""),
        "GTF_gene_name": r.get("GTF_gene_name", ""),
        "Motif_family_IDs": join_values(family_ids),
        "Motif_family_labels":
            join_values(family_labels),
        "Representative_motif_IDs":
            join_values(rep_motifs),
        "Evidence_sources":
            r.get("Evidence_sources", ""),
        "AME_exploratory_support":
            r.get("AME_exploratory_support", ""),
        "Default_Tomtom_support":
            r.get("Default_Tomtom_support", ""),
        "NoHoldout_Tomtom_support":
            r.get("NoHoldout_Tomtom_support", ""),
        "CrossRun_reproducible_Tomtom_support":
            r.get(
                "CrossRun_reproducible_Tomtom_support",
                ""
            ),
        "Best_AME_p":
            r.get("Best_AME_p", ""),
        "Best_AME_adj_p":
            r.get("Best_AME_adj_p", ""),
        "Best_AME_E":
            r.get("Best_AME_E", ""),
        "Best_Tomtom_q":
            r.get("Best_Tomtom_q", ""),
        "Target_tissue_mean_TPM":
            r.get("Target_tissue_mean_TPM", ""),
        "Target_tissue_P348_TPM":
            r.get("Target_tissue_P348_TPM", ""),
        "Target_tissue_P350_TPM":
            r.get("Target_tissue_P350_TPM", ""),
        "High_RNA_expression_support":
            r.get("High_RNA_expression_support", ""),
        "Tau":
            r.get("Tau", ""),
        "Target_tissue_SPM":
            r.get("Target_tissue_SPM", ""),
        "Target_tissue_expression_rank":
            r.get(
                "Target_tissue_expression_rank",
                ""
            ),
        "Max_expression_tissue":
            r.get("Max_expression_tissue", ""),
        "Target_is_max_expression_tissue":
            r.get(
                "Target_is_max_expression_tissue",
                ""
            ),
        "SPM_ge_0.5":
            r.get("SPM_ge_0.5", ""),
        "RNA_support_role": rna_role,
        "FIMO_best_tier": best_tier,
        "Interpretation":
            "RNA-supported candidate TF for motif family; "
            "not proof of direct DNA binding"
    })


tf_master.sort(
    key=lambda r: (
        tissue_order.get(r["Tissue"], 99),
        tier_order.get(r["FIMO_best_tier"], 99),
        -as_float(r["Target_tissue_mean_TPM"])
        if finite(as_float(r["Target_tissue_mean_TPM"]))
        else 0,
        r["TF_symbol"]
    )
)

tf_columns = [
    "Tissue",
    "TF_symbol",
    "gene_id",
    "GTF_gene_name",
    "Motif_family_IDs",
    "Motif_family_labels",
    "Representative_motif_IDs",
    "Evidence_sources",
    "AME_exploratory_support",
    "Default_Tomtom_support",
    "NoHoldout_Tomtom_support",
    "CrossRun_reproducible_Tomtom_support",
    "Best_AME_p",
    "Best_AME_adj_p",
    "Best_AME_E",
    "Best_Tomtom_q",
    "Target_tissue_mean_TPM",
    "Target_tissue_P348_TPM",
    "Target_tissue_P350_TPM",
    "High_RNA_expression_support",
    "Tau",
    "Target_tissue_SPM",
    "Target_tissue_expression_rank",
    "Max_expression_tissue",
    "Target_is_max_expression_tissue",
    "SPM_ge_0.5",
    "RNA_support_role",
    "FIMO_best_tier",
    "Interpretation"
]

write_tsv(
    OUT / "08D2_publication_TF_candidate_master.tsv",
    tf_master,
    tf_columns
)


# ============================================================
# PART D. FIMO manifests
# ============================================================

primary_manifest = [
    r for r in motif_family_rows
    if r["FIMO_tier"] == "Primary"
]

secondary_manifest = [
    r for r in motif_family_rows
    if r["FIMO_tier"] == "Secondary"
]

selected_manifest = (
    primary_manifest +
    secondary_manifest
)

manifest_columns = [
    "Tissue",
    "FIMO_tier",
    "Evidence_track",
    "Motif_family_ID",
    "Motif_family_label",
    "Representative_motif_ID",
    "Representative_motif_name",
    "Representative_source",
    "Consensus",
    "High_RNA_candidate_TFs",
    "Tissue_preferred_candidate_TFs",
    "AME_best_adj_p",
    "AME_best_E",
    "Default_STREME_holdout_E",
    "CrossRun_reproduced",
    "Best_default_Tomtom_q",
    "FIMO_rationale"
]

write_tsv(
    OUT / "08D3_FIMO_primary_motif_manifest.tsv",
    primary_manifest,
    manifest_columns
)

write_tsv(
    OUT / "08D4_FIMO_secondary_motif_manifest.tsv",
    secondary_manifest,
    manifest_columns
)

write_tsv(
    OUT / "08D5_FIMO_all_selected_motif_manifest.tsv",
    selected_manifest,
    manifest_columns
)


# ============================================================
# PART E. Extract motif blocks into per-tissue MEME files
# ============================================================

def parse_meme_blocks(path):

    with path.open() as fh:
        lines = fh.readlines()

    header = []
    blocks = {}
    current_id = None
    current = []

    seen_motif = False

    for line in lines:

        if line.startswith("MOTIF "):

            if current_id is not None:
                blocks[current_id] = current

            seen_motif = True

            parts = line.strip().split()
            current_id = parts[1]

            current = [line]

        elif seen_motif:
            current.append(line)

        else:
            header.append(line)

    if current_id is not None:
        blocks[current_id] = current

    return header, blocks


def write_selected_meme(
    source,
    selected_ids,
    output
):

    selected_ids = list(dict.fromkeys(selected_ids))

    if not selected_ids:
        return 0

    header, blocks = parse_meme_blocks(source)

    missing = [
        x for x in selected_ids
        if x not in blocks
    ]

    if missing:
        raise RuntimeError(
            f"Motifs missing from {source}: "
            + ",".join(missing)
        )

    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with output.open("w") as fh:

        for line in header:
            fh.write(line)

        if header and not header[-1].endswith("\n"):
            fh.write("\n")

        for motif_id in selected_ids:
            for line in blocks[motif_id]:
                fh.write(line)

    return len(selected_ids)


fimo_file_records = []

for tissue in TISSUES:

    tdir = FIMO_ROOT / tissue
    tdir.mkdir(
        parents=True,
        exist_ok=True
    )

    for tier in ["Primary", "Secondary"]:

        known_ids = [
            r["Representative_motif_ID"]
            for r in motif_family_rows
            if (
                r["Tissue"] == tissue
                and
                r["FIMO_tier"] == tier
                and
                r["Evidence_track"]
                == "Known_JASPAR_AME"
            )
        ]

        denovo_ids = [
            r["Representative_motif_ID"]
            for r in motif_family_rows
            if (
                r["Tissue"] == tissue
                and
                r["FIMO_tier"] == tier
                and
                r["Evidence_track"]
                == "De_novo_STREME"
            )
        ]

        known_out = (
            tdir /
            f"{tissue}.known.{tier.lower()}.meme"
        )

        denovo_out = (
            tdir /
            f"{tissue}.denovo.{tier.lower()}.meme"
        )

        nk = 0
        nd = 0

        if known_ids:
            nk = write_selected_meme(
                JASPAR_DB,
                known_ids,
                known_out
            )

        if denovo_ids:
            nd = write_selected_meme(
                Path(
                    f"08B2_STREME_de_novo/"
                    f"{tissue}/streme.txt"
                ),
                denovo_ids,
                denovo_out
            )

        fimo_file_records.append({
            "Tissue": tissue,
            "FIMO_tier": tier,
            "Known_motif_N": nk,
            "Known_MEME_file":
                str(known_out) if nk else "",
            "De_novo_motif_N": nd,
            "De_novo_MEME_file":
                str(denovo_out) if nd else ""
        })


write_tsv(
    OUT / "08D6_FIMO_MEME_file_manifest.tsv",
    fimo_file_records,
    [
        "Tissue",
        "FIMO_tier",
        "Known_motif_N",
        "Known_MEME_file",
        "De_novo_motif_N",
        "De_novo_MEME_file"
    ]
)


# ============================================================
# PART F. QC summary
# ============================================================

qc = []

for tissue in TISSUES:

    fam = [
        r for r in motif_family_rows
        if r["Tissue"] == tissue
    ]

    known = [
        r for r in fam
        if r["Evidence_track"] ==
        "Known_JASPAR_AME"
    ]

    denovo = [
        r for r in fam
        if r["Evidence_track"] ==
        "De_novo_STREME"
    ]

    qc.append({
        "Tissue": tissue,
        "Known_motif_families": len(known),
        "De_novo_motif_families": len(denovo),
        "Known_primary_FIMO": sum(
            r["FIMO_tier"] == "Primary"
            for r in known
        ),
        "Known_secondary_FIMO": sum(
            r["FIMO_tier"] == "Secondary"
            for r in known
        ),
        "De_novo_primary_FIMO": sum(
            r["FIMO_tier"] == "Primary"
            for r in denovo
        ),
        "De_novo_secondary_FIMO": sum(
            r["FIMO_tier"] == "Secondary"
            for r in denovo
        ),
        "Primary_FIMO_total": sum(
            r["FIMO_tier"] == "Primary"
            for r in fam
        ),
        "Secondary_FIMO_total": sum(
            r["FIMO_tier"] == "Secondary"
            for r in fam
        ),
        "Mapped_TF_candidates": len({
            r["TF_symbol"]
            for r in tf_master
            if (
                r["Tissue"] == tissue
                and
                r["gene_id"]
            )
        }),
        "High_RNA_TF_candidates": len({
            r["TF_symbol"]
            for r in tf_master
            if (
                r["Tissue"] == tissue
                and
                as_bool(
                    r[
                        "High_RNA_expression_support"
                    ]
                )
            )
        })
    })


write_tsv(
    OUT / "08D7_QC_summary.tsv",
    qc,
    [
        "Tissue",
        "Known_motif_families",
        "De_novo_motif_families",
        "Known_primary_FIMO",
        "Known_secondary_FIMO",
        "De_novo_primary_FIMO",
        "De_novo_secondary_FIMO",
        "Primary_FIMO_total",
        "Secondary_FIMO_total",
        "Mapped_TF_candidates",
        "High_RNA_TF_candidates"
    ]
)


# ============================================================
# PART G. Publication / methods note
# ============================================================

methods = """STEP 08D MOTIF-FAMILY AND TF-CANDIDATE INTEGRATION

Analysis role
-------------
This step prioritizes motif families for downstream FIMO site
mapping. It does not introduce a new statistical significance
test.

Formal motif statistics
-----------------------
Known motifs:
AME E-value < 0.05 remains the predefined database-wide formal
criterion. AME adjusted p-values < 0.05 are treated as
exploratory motif evidence because no tested tissue achieved
AME E < 0.05.

De novo motifs:
The default STREME run with an internal hold-out set remains the
primary de novo motif analysis. Its hold-out E-value is the
formal discovery statistic. The no-holdout STREME run is used
only as a sensitivity analysis for motif reproducibility.

FIMO Primary tier
-----------------
Known motif:
AME exploratory support + high RNA expression of at least one
candidate TF + target-tissue expression preference
(max-expression tissue or SPM >= 0.5).

De novo motif:
default STREME motif + Tomtom q < 0.05 JASPAR similarity +
cross-run motif reproducibility + at least one high-RNA TF
candidate.

FIMO Secondary tier
-------------------
Known motif:
AME exploratory support + high RNA expression, without strong
target-tissue preference.

De novo motif:
default STREME motif + Tomtom q < 0.05 + high-RNA TF candidate,
but without cross-run motif reproducibility.

Interpretation
--------------
Tomtom similarity does not establish exact TF identity.
Multiple KLF/SP/ZNF-family TFs matching the same GC-rich motif
are treated as alternative candidate TFs for a motif family,
rather than as independent regulatory events.

FIMO motif occurrence does not prove in vivo TF occupancy or a
causal TF-to-target-gene interaction.
"""

with (
    OUT / "08D8_methods_and_interpretation.txt"
).open("w") as fh:
    fh.write(methods)


# ============================================================
# Console summary
# ============================================================

print()
print("=" * 78)
print("STEP 08D COMPLETED")
print("=" * 78)

header = [
    "Tissue",
    "KnownPri",
    "KnownSec",
    "DeNovoPri",
    "DeNovoSec",
    "PrimaryTotal",
    "SecondaryTotal"
]

print(
    "{:<12s} {:>8s} {:>8s} {:>10s} {:>10s} {:>12s} {:>14s}".format(
        *header
    )
)

for q in qc:
    print(
        "{:<12s} {:>8d} {:>8d} {:>10d} {:>10d} {:>12d} {:>14d}".format(
            q["Tissue"],
            q["Known_primary_FIMO"],
            q["Known_secondary_FIMO"],
            q["De_novo_primary_FIMO"],
            q["De_novo_secondary_FIMO"],
            q["Primary_FIMO_total"],
            q["Secondary_FIMO_total"]
        )
    )

print()
print("Primary FIMO motifs:")
for r in primary_manifest:
    print(
        f"{r['Tissue']:12s} | "
        f"{r['Evidence_track']:18s} | "
        f"{r['Representative_motif_ID']:20s} | "
        f"{r['Motif_family_label']} | "
        f"TF={r['High_RNA_candidate_TFs']}"
    )

print()
print("Secondary FIMO motifs:")
for r in secondary_manifest:
    print(
        f"{r['Tissue']:12s} | "
        f"{r['Evidence_track']:18s} | "
        f"{r['Representative_motif_ID']:20s} | "
        f"{r['Motif_family_label']} | "
        f"TF={r['High_RNA_candidate_TFs']}"
    )

print()
print("Main outputs:")
print(
    OUT / "08D1_publication_motif_family_master.tsv"
)
print(
    OUT / "08D2_publication_TF_candidate_master.tsv"
)
print(
    OUT / "08D3_FIMO_primary_motif_manifest.tsv"
)
print(
    OUT / "08D4_FIMO_secondary_motif_manifest.tsv"
)
print(
    OUT / "08D6_FIMO_MEME_file_manifest.tsv"
)
print(
    OUT / "08D7_QC_summary.tsv"
)
print("=" * 78)

