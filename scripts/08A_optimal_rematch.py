#!/usr/bin/env python3

import csv
import math
import os
import subprocess
import sys
from collections import defaultdict, Counter

try:
    import numpy as np
    from scipy.optimize import linear_sum_assignment
except ImportError:
    sys.exit(
        "ERROR: numpy/scipy unavailable.\n"
        "This script requires scipy.optimize.linear_sum_assignment."
    )

# ============================================================
# STEP 08A-v2
# OPTIMAL GC/LENGTH/GENOMIC-CLASS MATCHING
#
# Foreground:
#   tissue-specific strong peaks
#
# Background:
#   same-tissue HC non-strong peaks
#
# Matching:
#   exact genomic class
#   + optimal GC/log2(length) assignment
#   + local distribution-balance optimization
#
# No replacement.
# ============================================================

MASTER = (
    "08A_matched_motif_background/"
    "08A_all_candidates_GC_length_class.tsv"
)

FASTA = "Sus_genome.fasta"

OUTDIR = "08A_matched_motif_background_v2"

TISSUES = [
    "Muscle",
    "Spleen",
    "Liver",
    "Cerebellum"
]

CLASSES = [
    "Promoter_TSSproximal",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Distal_to_associated_gene"
]

os.makedirs(
    OUTDIR,
    exist_ok=True
)

if not os.path.exists(MASTER):
    sys.exit(
        f"ERROR: missing {MASTER}"
    )

if not os.path.exists(FASTA):
    sys.exit(
        f"ERROR: missing {FASTA}"
    )

if not os.path.exists(FASTA + ".fai"):
    subprocess.run(
        ["samtools", "faidx", FASTA],
        check=True
    )


# ============================================================
# READ FASTA INDEX
# ============================================================

fasta_names = set()

with open(FASTA + ".fai") as f:
    for line in f:
        if line.strip():
            fasta_names.add(
                line.split("\t")[0]
            )


def resolve_chr(chrom):

    candidates = [
        chrom,
        "chr" + chrom
    ]

    if chrom in ("M", "MT"):
        candidates += [
            "chrM",
            "chrMT",
            "M",
            "MT"
        ]

    for x in candidates:
        if x in fasta_names:
            return x

    raise RuntimeError(
        f"Cannot resolve chromosome {chrom}"
    )


# ============================================================
# READ MASTER CANDIDATES
# ============================================================

rows = []

with open(MASTER) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for r in reader:

        r["Peak_length_bp"] = int(
            r["Peak_length_bp"]
        )

        r["GC_fraction"] = float(
            r["GC_fraction"]
        )

        r["N_fraction"] = float(
            r["N_fraction"]
        )

        r["log2_length"] = math.log2(
            r["Peak_length_bp"]
        )

        chrom, coords = r[
            "BED_peak_id"
        ].rsplit(":", 1)

        start, end = coords.split("-")

        r["chrom"] = chrom
        r["start"] = int(start)
        r["end"] = int(end)

        rows.append(r)


print("=" * 72)
print("STEP 08A-v2: OPTIMAL MATCHING")
print("=" * 72)

print(
    "Candidate rows:",
    len(rows)
)


# ============================================================
# BASIC STATISTICS
# ============================================================

def sample_variance(x):

    x = np.asarray(
        x,
        dtype=float
    )

    if len(x) < 2:
        return 0.0

    return float(
        np.var(
            x,
            ddof=1
        )
    )


def smd(x, y):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    pooled = math.sqrt(
        (
            sample_variance(x) +
            sample_variance(y)
        ) / 2
    )

    if pooled == 0:
        return 0.0

    return float(
        (
            np.mean(x) -
            np.mean(y)
        ) / pooled
    )


# ============================================================
# HUNGARIAN MATCHING FOR ONE CLASS
# ============================================================

def optimal_class_match(
    fg,
    bg_pool
):

    if len(bg_pool) < len(fg):

        raise RuntimeError(
            "Insufficient background pool."
        )

    # Standardize using the combined class-specific population.
    combined_gc = np.array(
        [
            x["GC_fraction"]
            for x in fg + bg_pool
        ],
        dtype=float
    )

    combined_len = np.array(
        [
            x["log2_length"]
            for x in fg + bg_pool
        ],
        dtype=float
    )

    gc_sd = np.std(
        combined_gc,
        ddof=1
    )

    len_sd = np.std(
        combined_len,
        ddof=1
    )

    if gc_sd == 0:
        gc_sd = 1.0

    if len_sd == 0:
        len_sd = 1.0


    nfg = len(fg)
    nbg = len(bg_pool)

    cost = np.zeros(
        (nfg, nbg),
        dtype=float
    )


    for i, a in enumerate(fg):

        for j, b in enumerate(bg_pool):

            dg = abs(
                a["GC_fraction"] -
                b["GC_fraction"]
            )

            dl = abs(
                a["log2_length"] -
                b["log2_length"]
            )

            base = (
                (
                    dg / gc_sd
                ) ** 2
                +
                (
                    dl / len_sd
                ) ** 2
            )


            ratio = (
                b["Peak_length_bp"] /
                a["Peak_length_bp"]
            )


            # Strongly favour the previous Strict definition.
            if (
                dg <= 0.02
                and
                0.80 <= ratio <= 1.25
            ):

                penalty = 0.0


            elif (
                dg <= 0.04
                and
                0.67 <= ratio <= 1.50
            ):

                penalty = 4.0


            else:

                penalty = 20.0


            cost[i, j] = (
                base +
                penalty
            )


    row_ind, col_ind = (
        linear_sum_assignment(
            cost
        )
    )


    matches = []

    for i, j in zip(
        row_ind,
        col_ind
    ):

        matches.append(
            {
                "fg": fg[i],
                "bg": bg_pool[j],
                "cost": float(
                    cost[i, j]
                )
            }
        )


    return matches


# ============================================================
# GROUP OBJECTIVE FOR LOCAL SWAPS
# ============================================================

def balance_objective(
    fg_rows,
    bg_rows
):

    gc_smd = smd(
        [
            x["GC_fraction"]
            for x in fg_rows
        ],
        [
            x["GC_fraction"]
            for x in bg_rows
        ]
    )

    len_smd = smd(
        [
            x["log2_length"]
            for x in fg_rows
        ],
        [
            x["log2_length"]
            for x in bg_rows
        ]
    )

    return (
        abs(gc_smd) +
        abs(len_smd),
        gc_smd,
        len_smd
    )


# ============================================================
# LOCAL BACKGROUND SWAP OPTIMIZATION
#
# Preserve exact genomic class while improving tissue-wide
# GC + log2(length) balance.
# ============================================================

def optimize_selected_background(
    tissue_fg,
    class_matches,
    class_pools,
    max_rounds=100
):

    selected = {}

    for cls in CLASSES:

        selected[cls] = [
            m["bg"]
            for m in class_matches.get(
                cls,
                []
            )
        ]


    def flatten():

        out = []

        for cls in CLASSES:
            out.extend(
                selected.get(
                    cls,
                    []
                )
            )

        return out


    current_bg = flatten()

    current_obj, _, _ = (
        balance_objective(
            tissue_fg,
            current_bg
        )
    )


    for round_no in range(
        max_rounds
    ):

        best = None
        best_obj = current_obj


        for cls in CLASSES:

            chosen = selected.get(
                cls,
                []
            )

            if not chosen:
                continue


            chosen_ids = {
                x["SAF_peak_id"]
                for x in chosen
            }


            unused = [
                x
                for x in class_pools[
                    cls
                ]
                if x["SAF_peak_id"]
                not in chosen_ids
            ]


            if not unused:
                continue


            # Evaluate one-for-one swaps.
            for i, old in enumerate(
                chosen
            ):

                for new in unused:

                    trial_bg = []

                    for cls2 in CLASSES:

                        if cls2 != cls:

                            trial_bg.extend(
                                selected.get(
                                    cls2,
                                    []
                                )
                            )

                        else:

                            z = list(
                                selected[
                                    cls
                                ]
                            )

                            z[i] = new

                            trial_bg.extend(z)


                    obj, g_smd, l_smd = (
                        balance_objective(
                            tissue_fg,
                            trial_bg
                        )
                    )


                    if obj < (
                        best_obj -
                        1e-10
                    ):

                        best_obj = obj

                        best = (
                            cls,
                            i,
                            new,
                            g_smd,
                            l_smd
                        )


        if best is None:
            break


        cls, i, new, _, _ = best

        selected[
            cls
        ][i] = new

        current_obj = best_obj


        bg_now = flatten()

        _, g_now, l_now = (
            balance_objective(
                tissue_fg,
                bg_now
            )
        )


        # More than sufficient for our QC criterion.
        if (
            abs(g_now) < 0.05
            and
            abs(l_now) < 0.05
        ):
            break


    return selected


# ============================================================
# SEQUENCE EXTRACTION
# ============================================================

def region_for(row):

    fasta_chr = resolve_chr(
        row["chrom"]
    )

    # BED start is 0-based; faidx is 1-based closed
    return (
        f"{fasta_chr}:"
        f"{row['start'] + 1}-"
        f"{row['end']}"
    )


def extract_sequences(data):

    if not data:
        return {}


    regions = [
        region_for(x)
        for x in data
    ]


    region_file = os.path.join(
        OUTDIR,
        ".tmp_regions.txt"
    )


    with open(
        region_file,
        "w"
    ) as out:

        for r in regions:
            out.write(
                r + "\n"
            )


    result = subprocess.run(
        [
            "samtools",
            "faidx",
            FASTA,
            "-r",
            region_file
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )


    seqs = {}

    header = None
    parts = []


    for line in result.stdout.splitlines():

        if line.startswith(">"):

            if header is not None:

                seqs[
                    header
                ] = "".join(
                    parts
                ).upper()


            header = line[
                1:
            ].split()[0]

            parts = []

        else:

            parts.append(
                line.strip()
            )


    if header is not None:

        seqs[
            header
        ] = "".join(
            parts
        ).upper()


    try:
        os.remove(
            region_file
        )
    except OSError:
        pass


    return seqs


def write_fasta(
    data,
    path,
    tissue
):

    seqs = extract_sequences(
        data
    )


    with open(
        path,
        "w"
    ) as out:

        for x in data:

            region = region_for(x)

            seq = seqs.get(
                region
            )

            if seq is None:

                raise RuntimeError(
                    f"Sequence missing: {region}"
                )


            if len(seq) != x[
                "Peak_length_bp"
            ]:

                raise RuntimeError(
                    f"Length mismatch: "
                    f"{x['SAF_peak_id']}"
                )


            out.write(
                f">{x['SAF_peak_id']}"
                f"|Tissue={tissue}"
                f"|Class={x['Pair_specific_class']}\n"
            )


            for i in range(
                0,
                len(seq),
                60
            ):

                out.write(
                    seq[
                        i:
                        i + 60
                    ] +
                    "\n"
                )


def write_bed(
    data,
    path
):

    # IMPORTANT:
    # use FASTA-compatible chromosome names (chr1 etc.)
    with open(
        path,
        "w"
    ) as out:

        for x in data:

            chrom = resolve_chr(
                x["chrom"]
            )

            out.write(
                "\t".join(
                    [
                        chrom,
                        str(
                            x["start"]
                        ),
                        str(
                            x["end"]
                        ),
                        x["SAF_peak_id"],
                        "0",
                        "."
                    ]
                ) +
                "\n"
            )


# ============================================================
# RUN FOUR TISSUES
# ============================================================

global_qc = []


for tissue in TISSUES:

    print(
        "\n" +
        "=" * 72
    )

    print(
        "Tissue:",
        tissue
    )


    fg_tissue = [
        x
        for x in rows
        if (
            x["Tissue"] ==
            tissue
            and
            x["Set"] ==
            "Strong_foreground"
        )
    ]


    bg_tissue = [
        x
        for x in rows
        if (
            x["Tissue"] ==
            tissue
            and
            x["Set"] ==
            "HC_nonstrong_background_pool"
        )
    ]


    class_matches = {}
    class_pools = {}


    for cls in CLASSES:

        fg = [
            x
            for x in fg_tissue
            if x[
                "Pair_specific_class"
            ] == cls
        ]


        bg = [
            x
            for x in bg_tissue
            if x[
                "Pair_specific_class"
            ] == cls
        ]


        class_pools[
            cls
        ] = bg


        if not fg:

            class_matches[
                cls
            ] = []

            continue


        print(
            f"{cls:28s}",
            f"FG={len(fg):3d}",
            f"BGpool={len(bg):4d}"
        )


        class_matches[
            cls
        ] = optimal_class_match(
            fg,
            bg
        )


    # --------------------------------------------------------
    # Optimize overall GC / length means using same-class swaps
    # --------------------------------------------------------

    selected = optimize_selected_background(
        fg_tissue,
        class_matches,
        class_pools
    )


    selected_bg = []

    for cls in CLASSES:
        selected_bg.extend(
            selected[
                cls
            ]
        )


    # --------------------------------------------------------
    # Re-pair foreground with FINAL selected BG within class,
    # solely for match-level QC/output.
    # --------------------------------------------------------

    final_matches = []


    for cls in CLASSES:

        fg = [
            x
            for x in fg_tissue
            if x[
                "Pair_specific_class"
            ] == cls
        ]


        bg = selected[
            cls
        ]


        if not fg:
            continue


        m = optimal_class_match(
            fg,
            bg
        )

        final_matches.extend(
            m
        )


    if len(
        final_matches
    ) != len(
        fg_tissue
    ):

        raise RuntimeError(
            f"{tissue}: matching count error."
        )


    # --------------------------------------------------------
    # Match tiers
    # --------------------------------------------------------

    strict = 0
    relaxed = 0
    nearest = 0

    match_rows = []


    for m in final_matches:

        fg = m["fg"]
        bg = m["bg"]

        dg = abs(
            fg["GC_fraction"] -
            bg["GC_fraction"]
        )

        ratio = (
            bg["Peak_length_bp"] /
            fg["Peak_length_bp"]
        )


        if (
            dg <= 0.02
            and
            0.80 <= ratio <= 1.25
        ):

            tier = "Strict"
            strict += 1


        elif (
            dg <= 0.04
            and
            0.67 <= ratio <= 1.50
        ):

            tier = "Relaxed"
            relaxed += 1


        else:

            tier = "Nearest_same_class"
            nearest += 1


        match_rows.append(
            {
                "Tissue":
                    tissue,

                "Genomic_class":
                    fg[
                        "Pair_specific_class"
                    ],

                "Foreground_peak_id":
                    fg[
                        "SAF_peak_id"
                    ],

                "Foreground_gene_id":
                    fg[
                        "gene_id"
                    ],

                "Foreground_length":
                    fg[
                        "Peak_length_bp"
                    ],

                "Foreground_GC":
                    fg[
                        "GC_fraction"
                    ],

                "Background_peak_id":
                    bg[
                        "SAF_peak_id"
                    ],

                "Background_gene_id":
                    bg[
                        "gene_id"
                    ],

                "Background_length":
                    bg[
                        "Peak_length_bp"
                    ],

                "Background_GC":
                    bg[
                        "GC_fraction"
                    ],

                "GC_absolute_difference":
                    dg,

                "Length_ratio_BG_over_FG":
                    ratio,

                "Match_tier":
                    tier,

                "Optimal_cost":
                    m[
                        "cost"
                    ]
            }
        )


    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    fg_gc = [
        x["GC_fraction"]
        for x in fg_tissue
    ]

    bg_gc = [
        x["GC_fraction"]
        for x in selected_bg
    ]

    fg_len = [
        x["log2_length"]
        for x in fg_tissue
    ]

    bg_len = [
        x["log2_length"]
        for x in selected_bg
    ]


    gc_smd = smd(
        fg_gc,
        bg_gc
    )

    len_smd = smd(
        fg_len,
        bg_len
    )


    qc = {
        "Tissue":
            tissue,

        "Foreground_N":
            len(
                fg_tissue
            ),

        "Background_pool_N":
            len(
                bg_tissue
            ),

        "Matched_background_N":
            len(
                selected_bg
            ),

        "Strict_matches":
            strict,

        "Relaxed_matches":
            relaxed,

        "Nearest_same_class_matches":
            nearest,

        "Strict_match_percent":
            100 *
            strict /
            len(
                fg_tissue
            ),

        "Foreground_mean_GC":
            np.mean(
                fg_gc
            ),

        "Background_mean_GC":
            np.mean(
                bg_gc
            ),

        "GC_SMD":
            gc_smd,

        "Foreground_mean_log2_length":
            np.mean(
                fg_len
            ),

        "Background_mean_log2_length":
            np.mean(
                bg_len
            ),

        "Log2_length_SMD":
            len_smd,

        "GC_balance_PASS_absSMD_lt_0.10":
            abs(
                gc_smd
            ) < 0.10,

        "Length_balance_PASS_absSMD_lt_0.10":
            abs(
                len_smd
            ) < 0.10
    }


    global_qc.append(
        qc
    )


    print(
        "Final GC SMD:",
        f"{gc_smd:.6f}"
    )

    print(
        "Final log2(length) SMD:",
        f"{len_smd:.6f}"
    )

    print(
        "Strict:",
        strict,
        "/",
        len(
            fg_tissue
        )
    )


    # --------------------------------------------------------
    # OUTPUT DIRECTORY
    # --------------------------------------------------------

    tissue_dir = os.path.join(
        OUTDIR,
        tissue
    )

    os.makedirs(
        tissue_dir,
        exist_ok=True
    )


    # Class QC
    class_qc_file = os.path.join(
        tissue_dir,
        f"{tissue}.genomic_class_matching_QC.tsv"
    )


    with open(
        class_qc_file,
        "w",
        newline=""
    ) as out:

        writer = csv.writer(
            out,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(
            [
                "Genomic_class",
                "Foreground_N",
                "Background_N",
                "Difference"
            ]
        )


        for cls in CLASSES:

            fgn = sum(
                x[
                    "Pair_specific_class"
                ] == cls
                for x in fg_tissue
            )

            bgn = sum(
                x[
                    "Pair_specific_class"
                ] == cls
                for x in selected_bg
            )

            writer.writerow(
                [
                    cls,
                    fgn,
                    bgn,
                    bgn - fgn
                ]
            )


    # Match table
    match_file = os.path.join(
        tissue_dir,
        f"{tissue}.foreground_background_matches.tsv"
    )


    with open(
        match_file,
        "w",
        newline=""
    ) as out:

        writer = csv.DictWriter(
            out,
            fieldnames=list(
                match_rows[
                    0
                ].keys()
            ),
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()

        writer.writerows(
            match_rows
        )


    # BED uses chr names that DIRECTLY match Sus_genome.fasta.
    write_bed(
        fg_tissue,
        os.path.join(
            tissue_dir,
            f"{tissue}.strong.foreground.fasta_compatible.bed"
        )
    )

    write_bed(
        selected_bg,
        os.path.join(
            tissue_dir,
            f"{tissue}.matched.background.fasta_compatible.bed"
        )
    )


    write_fasta(
        fg_tissue,
        os.path.join(
            tissue_dir,
            f"{tissue}.strong.foreground.fa"
        ),
        tissue
    )

    write_fasta(
        selected_bg,
        os.path.join(
            tissue_dir,
            f"{tissue}.matched.background.fa"
        ),
        tissue
    )


# ============================================================
# GLOBAL QC
# ============================================================

qc_file = os.path.join(
    OUTDIR,
    "08A_v2_matching_QC_summary.tsv"
)


with open(
    qc_file,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(
        out,
        fieldnames=list(
            global_qc[
                0
            ].keys()
        ),
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    writer.writerows(
        global_qc
    )


print(
    "\n" +
    "=" * 72
)

print(
    "STEP 08A-v2 COMPLETED"
)

print(
    "=" * 72
)

print(
    "QC:"
)

print(
    qc_file
)

print(
    "\nRequired final balance:"
)

print(
    "abs(GC_SMD) < 0.10"
)

print(
    "abs(Log2_length_SMD) < 0.10"
)

print(
    "=" * 72
)

