#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import re
import sys

print("=" * 78)
print("STEP10F4A — SUPPLEMENTARY TABLE SOURCE LOCK")
print("Frozen-result source inspection only")
print("No biological analysis rerun")
print("No statistic recomputed")
print("No Supplementary Table generated")
print("=" * 78)

ROOT = Path(".")

OUTDIR = Path(
    "10_publication_figures/"
    "Step10F_final_consistency/"
    "Step10F4_Supplementary_Tables/"
    "Step10F4A_source_lock"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================================
# 1. Helpers
# ======================================================================

def unique_find(basename, preferred_parent=None):
    """
    Locate exact basename.
    Prefer exact path under preferred_parent if supplied.
    Otherwise require unique project-wide match.
    """

    if preferred_parent is not None:
        p = ROOT / preferred_parent / basename

        if p.exists():
            return p

    hits = [
        x for x in ROOT.rglob(basename)
        if x.is_file()
        and "Step10F4A_source_lock" not in str(x)
    ]

    # Prefer files outside publication-figure derived outputs.
    primary = [
        x for x in hits
        if "10_publication_figures" not in str(x)
    ]

    if len(primary) == 1:
        return primary[0]

    if len(hits) == 1:
        return hits[0]

    return None


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace"
        )

    return open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    )


def table_info(path):
    """
    Count non-comment result rows and return header preview.
    Works for TSV/CSV and simple gz tables.
    """

    if path is None or not path.exists():
        return {
            "rows": None,
            "cols": None,
            "header": "",
            "delimiter": ""
        }

    try:
        with open_text(path) as f:

            lines = []

            for line in f:
                line = line.rstrip("\n\r")

                if not line:
                    continue

                if line.startswith("#"):
                    continue

                lines.append(line)

                if len(lines) >= 2:
                    break

        if not lines:
            return {
                "rows": 0,
                "cols": 0,
                "header": "",
                "delimiter": ""
            }

        header = lines[0]

        if "\t" in header:
            delim = "\t"
        elif "," in header:
            delim = ","
        else:
            delim = None

        if delim:
            cols = len(header.split(delim))
        else:
            cols = 1

        n = 0

        with open_text(path) as f:

            first_data_seen = False

            for line in f:

                line = line.rstrip("\n\r")

                if not line:
                    continue

                if line.startswith("#"):
                    continue

                if not first_data_seen:
                    first_data_seen = True
                    continue

                n += 1

        return {
            "rows": n,
            "cols": cols,
            "header": header,
            "delimiter": (
                "TAB" if delim == "\t"
                else "COMMA" if delim == ","
                else "OTHER"
            )
        }

    except Exception as e:

        return {
            "rows": None,
            "cols": None,
            "header": f"ERROR: {e}",
            "delimiter": ""
        }


def status_exact(path, expected_rows=None):

    if path is None or not path.exists():
        return "FAIL"

    info = table_info(path)

    if expected_rows is not None:
        if info["rows"] != expected_rows:
            return "FAIL"

    return "PASS"


def candidate_files(root_prefixes, keywords=None):

    roots = []

    for prefix in root_prefixes:

        exact = ROOT / prefix

        if exact.exists():
            roots.append(exact)

        else:
            roots.extend(
                [
                    x for x in ROOT.glob(prefix)
                    if x.is_dir()
                ]
            )

    files = []

    for r in roots:

        for ext in (
            "*.tsv",
            "*.csv",
            "*.txt",
            "*.tsv.gz",
            "*.csv.gz"
        ):

            for f in r.rglob(ext):

                if not f.is_file():
                    continue

                if keywords:

                    name = f.name.lower()

                    if not any(
                        k.lower() in name
                        for k in keywords
                    ):
                        continue

                files.append(f)

    return sorted(
        set(files),
        key=lambda x: str(x)
    )


# ======================================================================
# 2. Exact canonical sources already known
# ======================================================================

exact_sources = []


# ----------------------------------------------------------------------
# Table S1 — paired sample metadata
# ----------------------------------------------------------------------

S1 = unique_find(
    "sample_pair_metadata.tsv"
)

exact_sources.append({
    "Table": "S1",
    "Role": "Primary",
    "Description":
        "Paired multi-tissue RNA-seq/ATAC-seq sample metadata",
    "Path":
        str(S1) if S1 else "NOT_FOUND",
    "Expected_rows": 16,
    "Status":
        status_exact(S1, 16)
})


# ----------------------------------------------------------------------
# Table S3 — 702 strong peak-gene associations
# ----------------------------------------------------------------------

S3 = unique_find(
    "09B1_702_strong_pair_multi_evidence_master.tsv",
    "09B_702_pair_multi_evidence_master"
)

exact_sources.append({
    "Table": "S3",
    "Role": "Primary",
    "Description":
        "702 strongly coupled tissue-specific peak-gene associations",
    "Path":
        str(S3) if S3 else "NOT_FOUND",
    "Expected_rows": 702,
    "Status":
        status_exact(S3, 702)
})


# ----------------------------------------------------------------------
# Table S5 — matched-background QC component
# ----------------------------------------------------------------------

S5_MATCH = unique_find(
    "08A_v2_matching_QC_summary.tsv",
    "08A_matched_motif_background_v2"
)

exact_sources.append({
    "Table": "S5",
    "Role": "Matched_background_QC",
    "Description":
        "GC/length/class-matched motif-background QC summary",
    "Path":
        str(S5_MATCH) if S5_MATCH else "NOT_FOUND",
    "Expected_rows": 4,
    "Status":
        status_exact(S5_MATCH, 4)
})


# ----------------------------------------------------------------------
# Table S6 — integrated motif-peak-gene-TF evidence candidate source
# ----------------------------------------------------------------------

S6_MASTER = unique_find(
    "08E2E_high_RNA_TF_motif_peak_gene_master.tsv",
    "08E2_motif_peak_gene_TF_network"
)

exact_sources.append({
    "Table": "S6",
    "Role": "Candidate_primary",
    "Description":
        "Frozen high-RNA motif-peak-gene-TF integrated evidence master",
    "Path":
        str(S6_MASTER) if S6_MASTER else "NOT_FOUND",
    "Expected_rows": 3178,
    "Status":
        status_exact(S6_MASTER, 3178)
})


# ----------------------------------------------------------------------
# Table S7 — publication TF candidate master
# ----------------------------------------------------------------------

S7 = unique_find(
    "08D2_publication_TF_candidate_master.tsv",
    "08D_motif_family_TF_master"
)

exact_sources.append({
    "Table": "S7",
    "Role": "Primary",
    "Description":
        "Publication candidate TF RNA-support master",
    "Path":
        str(S7) if S7 else "NOT_FOUND",
    "Expected_rows": 48,
    "Status":
        status_exact(S7, 48)
})


# ----------------------------------------------------------------------
# Table S8 — module summary
# ----------------------------------------------------------------------

S8_MODULE = unique_find(
    "09C6_module_candidate_summary.tsv",
    "09C_publication_candidate_modules"
)

exact_sources.append({
    "Table": "S8",
    "Role": "Module_summary",
    "Description":
        "Frozen publication candidate-module summary",
    "Path":
        str(S8_MODULE) if S8_MODULE else "NOT_FOUND",
    "Expected_rows": 11,
    "Status":
        status_exact(S8_MODULE, 11)
})


# ======================================================================
# 3. Exact known structural sources for Table S4
# ======================================================================

s4_known = [
    (
        "06_UROPA_relative_location_summary.tsv",
        "Step06 UROPA relative-location summary"
    ),
    (
        "06_structural_class_summary.tsv",
        "Step06 pair-specific structural-class summary"
    ),
    (
        "06B_structural_enrichment_tests.tsv",
        "Step06B structural enrichment tests"
    ),
    (
        "06_coupling_by_structural_class.tsv",
        "Coupling by structural class"
    ),
    (
        "06_flanking_distance_summary.tsv",
        "Flanking-distance summary"
    ),
    (
        "06_702_strong_pairs_with_structure.tsv",
        "702 strong pairs with frozen structural annotation"
    ),
]

for basename, description in s4_known:

    p = unique_find(
        basename
    )

    exact_sources.append({
        "Table": "S4",
        "Role": "Component",
        "Description": description,
        "Path":
            str(p) if p else "NOT_FOUND",
        "Expected_rows": "",
        "Status":
            "PASS" if p and p.exists()
            else "CHECK"
    })


# ======================================================================
# 4. Restricted inventories for composite tables
#
# IMPORTANT:
# Classification is based on exact analysis-root prefixes.
# We do NOT classify AME by substring "ame", avoiding the old "same" bug.
# ======================================================================

inventory_specs = {

    "S2": {
        "roots": [
            "01*",
            "02*",
            "03*",
            "03B_SPM_analysis"
        ],

        "keywords": [
            "summary",
            "qc",
            "lrt",
            "tau",
            "spm",
            "specific"
        ]
    },

    "S4": {
        "roots": [
            "06*",
            "07*"
        ],

        "keywords": [
            "structure",
            "structural",
            "promoter",
            "tss",
            "annotation",
            "enrichment",
            "flanking",
            "702"
        ]
    },

    "S5": {
        "roots": [
            "08A_matched_motif_background_v2",
            "08B1_AME_known_motif",
            "08B1b_AME_all_motifs",
            "08B1c_AME_final_summary",
            "08B2_STREME_de_novo",
            "08B2b_STREME_noholdout",
            "08B4B_STREME_reproducibility"
        ],

        "keywords": None
    },

    "S6": {
        "roots": [
            "08E2_motif_peak_gene_TF_network"
        ],

        "keywords": None
    },

    "S8": {
        "roots": [
            "09C_publication_candidate_modules"
        ],

        "keywords": None
    }
}


candidate_inventory = []

for table, spec in inventory_specs.items():

    files = candidate_files(
        spec["roots"],
        spec["keywords"]
    )

    for f in files:

        info = table_info(f)

        candidate_inventory.append({
            "Table": table,
            "Path": str(f),
            "Basename": f.name,
            "Rows":
                info["rows"]
                if info["rows"] is not None
                else "NA",
            "Columns":
                info["cols"]
                if info["cols"] is not None
                else "NA",
            "Delimiter":
                info["delimiter"],
            "Header_preview":
                info["header"][:500]
        })


# ======================================================================
# 5. Exact-source information
# ======================================================================

exact_info = []

for x in exact_sources:

    p = Path(x["Path"]) \
        if x["Path"] != "NOT_FOUND" \
        else None

    info = table_info(p) \
        if p else {
            "rows": None,
            "cols": None,
            "header": "",
            "delimiter": ""
        }

    exact_info.append({
        **x,
        "Observed_rows":
            info["rows"]
            if info["rows"] is not None
            else "NA",
        "Observed_columns":
            info["cols"]
            if info["cols"] is not None
            else "NA",
        "Delimiter":
            info["delimiter"],
        "Header_preview":
            info["header"][:700]
    })


# ======================================================================
# 6. Core hard QC
# ======================================================================

hard_required = [
    x for x in exact_info
    if (
        x["Table"] in {
            "S1",
            "S3",
            "S7"
        }
        or
        (
            x["Table"] == "S5"
            and x["Role"] ==
                "Matched_background_QC"
        )
        or
        (
            x["Table"] == "S6"
            and x["Role"] ==
                "Candidate_primary"
        )
        or
        (
            x["Table"] == "S8"
            and x["Role"] ==
                "Module_summary"
        )
    )
]

hard_pass = all(
    x["Status"] == "PASS"
    for x in hard_required
)


# Candidate inventory availability
inventory_counts = {}

for table in [
    "S2",
    "S4",
    "S5",
    "S6",
    "S8"
]:

    inventory_counts[table] = sum(
        x["Table"] == table
        for x in candidate_inventory
    )


inventory_pass = all(
    inventory_counts[x] > 0
    for x in inventory_counts
)


# ======================================================================
# 7. Write outputs
# ======================================================================

exact_file = (
    OUTDIR /
    "10F4A_exact_source_lock.tsv"
)

with open(
    exact_file,
    "w",
    newline=""
) as f:

    fields = [
        "Table",
        "Role",
        "Description",
        "Path",
        "Expected_rows",
        "Observed_rows",
        "Observed_columns",
        "Delimiter",
        "Status",
        "Header_preview"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for x in exact_info:
        writer.writerow(
            {
                k: x.get(k, "")
                for k in fields
            }
        )


inventory_file = (
    OUTDIR /
    "10F4A_component_candidate_inventory.tsv"
)

with open(
    inventory_file,
    "w",
    newline=""
) as f:

    fields = [
        "Table",
        "Path",
        "Basename",
        "Rows",
        "Columns",
        "Delimiter",
        "Header_preview"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        candidate_inventory
    )


count_file = (
    OUTDIR /
    "10F4A_component_inventory_counts.tsv"
)

with open(
    count_file,
    "w"
) as f:

    f.write(
        "Table\tCandidate_files\n"
    )

    for table in [
        "S2",
        "S4",
        "S5",
        "S6",
        "S8"
    ]:

        f.write(
            f"{table}\t"
            f"{inventory_counts[table]}\n"
        )


# ======================================================================
# 8. Human-readable preview
# ======================================================================

preview_file = (
    OUTDIR /
    "10F4A_source_lock_preview.txt"
)

with open(
    preview_file,
    "w"
) as f:

    f.write(
        "=" * 78 + "\n"
    )
    f.write(
        "STEP10F4A — SUPPLEMENTARY TABLE SOURCE LOCK PREVIEW\n"
    )
    f.write(
        "=" * 78 + "\n\n"
    )

    f.write(
        "================ EXACT / CANONICAL SOURCES ================\n\n"
    )

    for x in exact_info:

        f.write(
            f"TABLE       : {x['Table']}\n"
        )
        f.write(
            f"ROLE        : {x['Role']}\n"
        )
        f.write(
            f"DESCRIPTION : {x['Description']}\n"
        )
        f.write(
            f"PATH        : {x['Path']}\n"
        )
        f.write(
            f"ROWS        : {x['Observed_rows']}\n"
        )
        f.write(
            f"COLUMNS     : {x['Observed_columns']}\n"
        )
        f.write(
            f"STATUS      : {x['Status']}\n"
        )
        f.write(
            f"HEADER      : {x['Header_preview']}\n\n"
        )

    f.write(
        "================ COMPONENT INVENTORY COUNTS ================\n\n"
    )

    for table in [
        "S2",
        "S4",
        "S5",
        "S6",
        "S8"
    ]:

        f.write(
            f"{table}: "
            f"{inventory_counts[table]} candidate files\n"
        )

    f.write("\n")

    for table in [
        "S2",
        "S4",
        "S5",
        "S6",
        "S8"
    ]:

        f.write(
            f"================ {table} CANDIDATES ================\n"
        )

        subset = [
            x for x in candidate_inventory
            if x["Table"] == table
        ]

        for x in subset:

            f.write(
                f"FILE : {x['Path']}\n"
            )
            f.write(
                f"ROWS : {x['Rows']}\n"
            )
            f.write(
                f"COLS : {x['Columns']}\n"
            )
            f.write(
                f"HEAD : {x['Header_preview']}\n\n"
            )


# ======================================================================
# 9. Overall status
# ======================================================================

overall = (
    "PASS"
    if hard_pass and inventory_pass
    else "CHECK"
)


status_file = (
    OUTDIR /
    "10F4A_overall_status.tsv"
)

with open(
    status_file,
    "w"
) as f:

    f.write(
        "Metric\tValue\n"
    )

    f.write(
        f"Hard_exact_sources_PASS\t"
        f"{int(hard_pass)}\n"
    )

    f.write(
        f"Composite_inventory_available\t"
        f"{int(inventory_pass)}\n"
    )

    f.write(
        f"Step10F4A_source_lock\t"
        f"{overall}\n"
    )

    f.write(
        "Biological_results_modified\t0\n"
    )


# ======================================================================
# 10. Console report
# ======================================================================

print()
print("=" * 78)
print("EXACT / CANONICAL SOURCE QC")
print("=" * 78)

for x in exact_info:

    print(
        f"{x['Table']:>2} "
        f"{x['Role']:<24} "
        f"rows={str(x['Observed_rows']):>6} "
        f"cols={str(x['Observed_columns']):>4} "
        f"{x['Status']}"
    )


print()
print("=" * 78)
print("COMPOSITE SOURCE INVENTORY")
print("=" * 78)

for table in [
    "S2",
    "S4",
    "S5",
    "S6",
    "S8"
]:

    print(
        f"{table}: "
        f"{inventory_counts[table]} candidate files"
    )


print()
print("=" * 78)

if overall == "PASS":

    print(
        "STEP10F4A STATUS: PASS"
    )

    print(
        "Canonical sources are locked and composite-table "
        "components are ready for final selection."
    )

else:

    print(
        "STEP10F4A STATUS: CHECK"
    )

    print(
        "Review missing/ambiguous source files before building tables."
    )

print(
    "No Supplementary Table was generated."
)

print(
    "Biological results modified: 0"
)

print("=" * 78)

print("\nOutputs:")
print(exact_file)
print(inventory_file)
print(count_file)
print(preview_file)
print(status_file)
