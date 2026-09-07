#!/usr/bin/env python3

import argparse
import csv
import math
import os
import random
import subprocess
import sys
from collections import defaultdict, Counter


# ============================================================
# STEP 08A
# TISSUE-SPECIFIC MATCHED MOTIF FOREGROUND / BACKGROUND
#
# FOREGROUND:
#   Strong RNA-ATAC associations
#
# BACKGROUND POOL:
#   Concordant HC associations excluding strong associations
#
# Matching:
#   1. same tissue
#   2. EXACT Pair_specific_class
#   3. similar GC
#   4. similar peak length
#
# 1:1 matching without replacement
#
# Primary tissues:
#   Muscle
#   Spleen
#   Liver
#   Cerebellum
#
# Input genomic class comes from Step07:
# promoter-first pair-specific classification.
# ============================================================


parser = argparse.ArgumentParser()

parser.add_argument(
    "--fasta",
    required=True,
    help="Sscrofa11.1 reference genome FASTA"
)

parser.add_argument(
    "--hc",
    default=(
        "07_promoter_first_enrichment/"
        "07B_4674_HC_promoter_first.tsv"
    )
)

parser.add_argument(
    "--strong",
    default=(
        "07_promoter_first_enrichment/"
        "07B_702_strong_promoter_first.tsv"
    )
)

parser.add_argument(
    "--outdir",
    default="08A_matched_motif_background"
)

parser.add_argument(
    "--seed",
    type=int,
    default=20260830
)

args = parser.parse_args()


TISSUES = [
    "Muscle",
    "Spleen",
    "Liver",
    "Cerebellum"
]


CLASS_LEVELS = [
    "Promoter_TSSproximal",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Distal_to_associated_gene"
]


random.seed(
    args.seed
)


# ============================================================
# CHECK
# ============================================================

print("=" * 70)
print("STEP 08A: MATCHED MOTIF FOREGROUND / BACKGROUND")
print("=" * 70)

print("\nReference FASTA:")
print(args.fasta)

print("\nHC table:")
print(args.hc)

print("\nStrong table:")
print(args.strong)


for f in [
    args.fasta,
    args.hc,
    args.strong
]:
    if not os.path.exists(f):
        sys.exit(
            f"\nERROR: file not found: {f}"
        )


if subprocess.call(
    ["which", "samtools"],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL
) != 0:

    sys.exit(
        "\nERROR: samtools was not found in PATH."
    )


os.makedirs(
    args.outdir,
    exist_ok=True
)


# ============================================================
# INDEX FASTA
# ============================================================

fai = (
    args.fasta + ".fai"
)


if not os.path.exists(fai):

    print(
        "\nFASTA index not found. Running samtools faidx..."
    )

    subprocess.run(
        [
            "samtools",
            "faidx",
            args.fasta
        ],
        check=True
    )


# ============================================================
# FASTA SEQUENCE NAMES
# ============================================================

fasta_names = set()


with open(fai) as fh:

    for line in fh:

        fields = line.rstrip(
            "\n"
        ).split(
            "\t"
        )

        fasta_names.add(
            fields[0]
        )


print(
    "\nFASTA sequences:",
    len(
        fasta_names
    )
)


def fasta_chromosome_name(chrom):

    candidates = [
        chrom,
        "chr" + chrom
    ]


    if chrom in (
        "MT",
        "M"
    ):

        candidates.extend(
            [
                "MT",
                "M",
                "chrM",
                "chrMT"
            ]
        )


    for x in candidates:

        if x in fasta_names:
            return x


    raise RuntimeError(
        "Cannot find chromosome/contig "
        f"'{chrom}' in FASTA index."
    )


# ============================================================
# BED ID PARSER
# ============================================================

def parse_bed_peak_id(x):

    chrom, coords = x.rsplit(
        ":",
        1
    )

    start, end = coords.split(
        "-"
    )

    start = int(start)
    end = int(end)

    if end <= start:

        raise ValueError(
            f"Invalid BED peak: {x}"
        )


    return (
        chrom,
        start,
        end
    )


# ============================================================
# READ TABLE
# ============================================================

def read_table(path):

    data = []


    with open(path) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )


        required = {
            "SAF_peak_id",
            "BED_peak_id",
            "gene_id",
            "ATAC_Max_tissue",
            "Pair_specific_class"
        }


        missing = (
            required -
            set(
                reader.fieldnames
            )
        )


        if missing:

            raise RuntimeError(
                f"{path} missing columns: "
                + ",".join(
                    sorted(missing)
                )
            )


        for row in reader:

            chrom, start, end = (
                parse_bed_peak_id(
                    row[
                        "BED_peak_id"
                    ]
                )
            )


            row["_chrom"] = chrom
            row["_start"] = start
            row["_end"] = end
            row["_length"] = (
                end - start
            )


            data.append(
                row
            )


    return data


hc = read_table(
    args.hc
)

strong = read_table(
    args.strong
)


print(
    "\nHC associations:",
    len(hc)
)

print(
    "Strong associations:",
    len(strong)
)


if len(hc) != 4674:

    print(
        "WARNING: expected 4674 HC associations."
    )


if len(strong) != 702:

    print(
        "WARNING: expected 702 strong associations."
    )


# ============================================================
# UNIQUE PEAK QC
# ============================================================

def duplicated_ids(data):

    ids = [
        x[
            "SAF_peak_id"
        ]
        for x in data
    ]

    return (
        len(ids) -
        len(set(ids))
    )


if duplicated_ids(hc) > 0:

    sys.exit(
        "ERROR: duplicated HC SAF_peak_id."
    )


if duplicated_ids(strong) > 0:

    sys.exit(
        "ERROR: duplicated strong SAF_peak_id."
    )


strong_ids = {
    x[
        "SAF_peak_id"
    ]
    for x in strong
}


hc_ids = {
    x[
        "SAF_peak_id"
    ]
    for x in hc
}


if not strong_ids.issubset(
    hc_ids
):

    sys.exit(
        "ERROR: strong peaks are not a subset of HC peaks."
    )


# ============================================================
# RESTRICT TO FOUR TISSUES
# ============================================================

foreground = [
    x
    for x in strong
    if x[
        "ATAC_Max_tissue"
    ] in TISSUES
]


background_pool = [
    x
    for x in hc
    if (
        x[
            "ATAC_Max_tissue"
        ] in TISSUES
        and
        x[
            "SAF_peak_id"
        ] not in strong_ids
    )
]


# ============================================================
# PRE-MATCHING COUNTS
# ============================================================

print(
    "\nForeground / background pool:"
)


for tissue in TISSUES:

    fg = [
        x
        for x in foreground
        if x[
            "ATAC_Max_tissue"
        ] == tissue
    ]


    bg = [
        x
        for x in background_pool
        if x[
            "ATAC_Max_tissue"
        ] == tissue
    ]


    print(
        f"{tissue:12s}",
        f"foreground={len(fg):4d}",
        f"background_pool={len(bg):4d}"
    )


# ============================================================
# EXACT GENOMIC-CLASS FEASIBILITY
# ============================================================

feasibility_rows = []

matching_impossible = False


for tissue in TISSUES:

    fg = [
        x
        for x in foreground
        if x[
            "ATAC_Max_tissue"
        ] == tissue
    ]


    bg = [
        x
        for x in background_pool
        if x[
            "ATAC_Max_tissue"
        ] == tissue
    ]


    fg_counts = Counter(
        x[
            "Pair_specific_class"
        ]
        for x in fg
    )


    bg_counts = Counter(
        x[
            "Pair_specific_class"
        ]
        for x in bg
    )


    for cls in CLASS_LEVELS:

        fgn = fg_counts[
            cls
        ]

        bgn = bg_counts[
            cls
        ]


        feasible = (
            bgn >= fgn
        )


        if (
            fgn > 0 and
            not feasible
        ):

            matching_impossible = True


        feasibility_rows.append(
            {
                "Tissue":
                    tissue,

                "Genomic_class":
                    cls,

                "Foreground_N":
                    fgn,

                "Background_pool_N":
                    bgn,

                "Exact_1to1_feasible":
                    feasible
            }
        )


feasibility_file = os.path.join(
    args.outdir,
    "08A_exact_class_matching_feasibility.tsv"
)


with open(
    feasibility_file,
    "w",
    newline=""
) as out:

    fields = [
        "Tissue",
        "Genomic_class",
        "Foreground_N",
        "Background_pool_N",
        "Exact_1to1_feasible"
    ]


    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writeheader()

    writer.writerows(
        feasibility_rows
    )


print(
    "\nExact genomic-class matching feasibility:"
)


for row in feasibility_rows:

    if (
        row[
            "Foreground_N"
        ] > 0
    ):

        print(
            row[
                "Tissue"
            ],
            row[
                "Genomic_class"
            ],
            "FG=",
            row[
                "Foreground_N"
            ],
            "BG=",
            row[
                "Background_pool_N"
            ],
            "PASS" if row[
                "Exact_1to1_feasible"
            ] else "FAIL"
        )


if matching_impossible:

    sys.exit(
        "\nERROR: exact 1:1 genomic-class matching "
        "is impossible for at least one tissue/class.\n"
        "Inspect:\n"
        + feasibility_file
    )


# ============================================================
# PREPARE REGIONS FOR FASTA EXTRACTION
# ============================================================

all_candidates = (
    foreground +
    background_pool
)


region_to_records = defaultdict(
    list
)


regions = []


for x in all_candidates:

    fasta_chr = fasta_chromosome_name(
        x[
            "_chrom"
        ]
    )


    # BED 0-based half-open:
    # start/end
    #
    # samtools faidx:
    # 1-based inclusive

    region = (
        f"{fasta_chr}:"
        f"{x['_start'] + 1}-"
        f"{x['_end']}"
    )


    x[
        "_fasta_chr"
    ] = fasta_chr

    x[
        "_region"
    ] = region


    region_to_records[
        region
    ].append(
        x
    )


    regions.append(
        region
    )


regions = list(
    dict.fromkeys(
        regions
    )
)


region_file = os.path.join(
    args.outdir,
    "08A_regions_for_GC.faidx.txt"
)


with open(
    region_file,
    "w"
) as out:

    for region in regions:

        out.write(
            region +
            "\n"
        )


# ============================================================
# EXTRACT SEQUENCES
# ============================================================

sequence_file = os.path.join(
    args.outdir,
    "08A_all_candidate_sequences.tmp.fa"
)


print(
    "\nExtracting sequences for GC calculation..."
)


with open(
    sequence_file,
    "w"
) as out:

    subprocess.run(
        [
            "samtools",
            "faidx",
            args.fasta,
            "-r",
            region_file
        ],
        stdout=out,
        check=True
    )


# ============================================================
# READ EXTRACTED FASTA
# ============================================================

def read_fasta(path):

    result = {}

    header = None
    seq = []


    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue


            if line.startswith(
                ">"
            ):

                if header is not None:

                    result[
                        header
                    ] = "".join(
                        seq
                    ).upper()


                header = (
                    line[
                        1:
                    ].split()[0]
                )

                seq = []


            else:

                seq.append(
                    line
                )


    if header is not None:

        result[
            header
        ] = "".join(
            seq
        ).upper()


    return result


sequences = read_fasta(
    sequence_file
)


missing_regions = [
    x
    for x in regions
    if x not in sequences
]


if missing_regions:

    sys.exit(
        "\nERROR: some regions were not recovered "
        "from the FASTA.\nExamples:\n"
        + "\n".join(
            missing_regions[:20]
        )
    )


# ============================================================
# GC / N FRACTION
# ============================================================

for x in all_candidates:

    seq = sequences[
        x[
            "_region"
        ]
    ]


    expected_length = x[
        "_length"
    ]


    if len(seq) != expected_length:

        raise RuntimeError(
            "Sequence-length mismatch for "
            + x[
                "SAF_peak_id"
            ]
            + f": BED={expected_length}, FASTA={len(seq)}"
        )


    a = seq.count(
        "A"
    )

    c = seq.count(
        "C"
    )

    g = seq.count(
        "G"
    )

    t = seq.count(
        "T"
    )

    valid = (
        a + c + g + t
    )


    if valid > 0:

        gc = (
            g + c
        ) / valid

    else:

        gc = float(
            "nan"
        )


    n_fraction = (
        len(seq) -
        valid
    ) / len(seq)


    x[
        "_GC"
    ] = gc

    x[
        "_N_fraction"
    ] = n_fraction

    x[
        "_sequence"
    ] = seq


# ============================================================
# WRITE MASTER GC/LENGTH TABLE
# ============================================================

master_file = os.path.join(
    args.outdir,
    "08A_all_candidates_GC_length_class.tsv"
)


with open(
    master_file,
    "w",
    newline=""
) as out:

    fields = [
        "Set",
        "Tissue",
        "SAF_peak_id",
        "BED_peak_id",
        "gene_id",
        "Pair_specific_class",
        "Peak_length_bp",
        "GC_fraction",
        "N_fraction"
    ]


    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writeheader()


    for x in all_candidates:

        set_name = (
            "Strong_foreground"
            if x[
                "SAF_peak_id"
            ] in strong_ids
            else
            "HC_nonstrong_background_pool"
        )


        writer.writerow(
            {
                "Set":
                    set_name,

                "Tissue":
                    x[
                        "ATAC_Max_tissue"
                    ],

                "SAF_peak_id":
                    x[
                        "SAF_peak_id"
                    ],

                "BED_peak_id":
                    x[
                        "BED_peak_id"
                    ],

                "gene_id":
                    x[
                        "gene_id"
                    ],

                "Pair_specific_class":
                    x[
                        "Pair_specific_class"
                    ],

                "Peak_length_bp":
                    x[
                        "_length"
                    ],

                "GC_fraction":
                    f"{x['_GC']:.8f}",

                "N_fraction":
                    f"{x['_N_fraction']:.8f}"
            }
        )


# ============================================================
# MATCHING FUNCTIONS
# ============================================================

def match_cost(
    fg,
    bg
):

    gc_diff = abs(
        fg[
            "_GC"
        ] -
        bg[
            "_GC"
        ]
    )


    length_log2_ratio = abs(
        math.log2(
            bg[
                "_length"
            ] /
            fg[
                "_length"
            ]
        )
    )


    # Standardized matching cost.
    #
    # 0.02 GC difference and
    # 25% length ratio are approximately
    # one cost unit.

    cost = math.sqrt(
        (
            gc_diff /
            0.02
        ) ** 2
        +
        (
            length_log2_ratio /
            math.log2(
                1.25
            )
        ) ** 2
    )


    return (
        cost,
        gc_diff,
        length_log2_ratio
    )


def match_tier(
    fg,
    bg
):

    gc_diff = abs(
        fg[
            "_GC"
        ] -
        bg[
            "_GC"
        ]
    )


    ratio = (
        bg[
            "_length"
        ] /
        fg[
            "_length"
        ]
    )


    # Strict
    if (
        gc_diff <= 0.02
        and
        0.80 <= ratio <= 1.25
    ):
        return 1


    # Relaxed
    if (
        gc_diff <= 0.04
        and
        0.67 <= ratio <= 1.50
    ):
        return 2


    # Exact class but nearest available
    return 3


def match_label(tier):

    return {
        1:
            "Strict",

        2:
            "Relaxed",

        3:
            "Nearest_same_class"
    }[
        tier
    ]


# ============================================================
# STANDARDIZED MEAN DIFFERENCE
# ============================================================

def mean(values):

    if not values:
        return float(
            "nan"
        )

    return sum(
        values
    ) / len(
        values
    )


def variance(values):

    if len(
        values
    ) < 2:

        return 0.0


    m = mean(
        values
    )


    return sum(
        (
            x - m
        ) ** 2
        for x in values
    ) / (
        len(
            values
        ) - 1
    )


def smd(
    x,
    y
):

    mx = mean(
        x
    )

    my = mean(
        y
    )


    pooled = math.sqrt(
        (
            variance(x) +
            variance(y)
        ) /
        2.0
    )


    if pooled == 0:

        return 0.0


    return (
        mx - my
    ) / pooled


def quantile(
    values,
    q
):

    if not values:
        return float(
            "nan"
        )


    values = sorted(
        values
    )


    pos = (
        len(values) - 1
    ) * q


    low = math.floor(
        pos
    )

    high = math.ceil(
        pos
    )


    if low == high:

        return values[
            int(pos)
        ]


    return (
        values[low] *
        (
            high - pos
        )
        +
        values[high] *
        (
            pos - low
        )
    )


# ============================================================
# MATCH EACH TISSUE
# ============================================================

global_qc = []


for tissue in TISSUES:

    print(
        "\n" +
        "=" * 70
    )

    print(
        "Matching tissue:",
        tissue
    )


    tissue_dir = os.path.join(
        args.outdir,
        tissue
    )


    os.makedirs(
        tissue_dir,
        exist_ok=True
    )


    fg_tissue = [
        x
        for x in foreground
        if x[
            "ATAC_Max_tissue"
        ] == tissue
    ]


    bg_tissue = [
        x
        for x in background_pool
        if x[
            "ATAC_Max_tissue"
        ] == tissue
    ]


    matches = []


    for cls in CLASS_LEVELS:

        fg_class = [
            x
            for x in fg_tissue
            if x[
                "Pair_specific_class"
            ] == cls
        ]


        bg_class = [
            x
            for x in bg_tissue
            if x[
                "Pair_specific_class"
            ] == cls
        ]


        if not fg_class:
            continue


        if len(
            bg_class
        ) < len(
            fg_class
        ):

            raise RuntimeError(
                f"{tissue}/{cls}: "
                "insufficient same-class background."
            )


        # ----------------------------------------------------
        # Determine how difficult each foreground peak is.
        #
        # Peaks with fewer strict candidates are matched first.
        # ----------------------------------------------------

        difficulty = []


        for fg in fg_class:

            strict_n = 0
            relaxed_n = 0


            for bg in bg_class:

                tier = match_tier(
                    fg,
                    bg
                )


                if tier == 1:
                    strict_n += 1


                if tier <= 2:
                    relaxed_n += 1


            difficulty.append(
                (
                    strict_n,
                    relaxed_n,
                    random.random(),
                    fg
                )
            )


        difficulty.sort(
            key=lambda x: (
                x[0],
                x[1],
                x[2]
            )
        )


        available = {
            x[
                "SAF_peak_id"
            ]:
                x
            for x in bg_class
        }


        for (
            strict_n,
            relaxed_n,
            _,
            fg
        ) in difficulty:


            candidates = []


            for bg in available.values():

                cost, gc_diff, log2_ratio = (
                    match_cost(
                        fg,
                        bg
                    )
                )


                tier = match_tier(
                    fg,
                    bg
                )


                candidates.append(
                    (
                        tier,
                        cost,
                        random.random(),
                        gc_diff,
                        log2_ratio,
                        bg
                    )
                )


            if not candidates:

                raise RuntimeError(
                    f"No available background for "
                    f"{tissue}/{cls}/"
                    f"{fg['SAF_peak_id']}"
                )


            candidates.sort(
                key=lambda x: (
                    x[0],
                    x[1],
                    x[2]
                )
            )


            (
                tier,
                cost,
                _,
                gc_diff,
                log2_ratio,
                bg
            ) = candidates[0]


            del available[
                bg[
                    "SAF_peak_id"
                ]
            ]


            matches.append(
                {
                    "foreground":
                        fg,

                    "background":
                        bg,

                    "tier":
                        tier,

                    "cost":
                        cost,

                    "gc_diff":
                        gc_diff,

                    "abs_log2_length_ratio":
                        log2_ratio
                }
            )


    if len(
        matches
    ) != len(
        fg_tissue
    ):

        raise RuntimeError(
            f"{tissue}: expected "
            f"{len(fg_tissue)} matches, "
            f"obtained {len(matches)}."
        )


    # ========================================================
    # OUTPUT MATCH PAIRS
    # ========================================================

    match_file = os.path.join(
        tissue_dir,
        f"{tissue}.foreground_background_matches.tsv"
    )


    with open(
        match_file,
        "w",
        newline=""
    ) as out:

        fields = [
            "Tissue",
            "Genomic_class",

            "Foreground_peak_id",
            "Foreground_gene_id",
            "Foreground_length",
            "Foreground_GC",

            "Background_peak_id",
            "Background_gene_id",
            "Background_length",
            "Background_GC",

            "GC_absolute_difference",
            "Length_ratio_BG_over_FG",
            "Abs_log2_length_ratio",
            "Match_tier",
            "Match_cost"
        ]


        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )


        writer.writeheader()


        for m in matches:

            fg = m[
                "foreground"
            ]

            bg = m[
                "background"
            ]


            writer.writerow(
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
                            "_length"
                        ],

                    "Foreground_GC":
                        f"{fg['_GC']:.8f}",

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
                            "_length"
                        ],

                    "Background_GC":
                        f"{bg['_GC']:.8f}",

                    "GC_absolute_difference":
                        f"{m['gc_diff']:.8f}",

                    "Length_ratio_BG_over_FG":
                        f"{bg['_length']/fg['_length']:.8f}",

                    "Abs_log2_length_ratio":
                        f"{m['abs_log2_length_ratio']:.8f}",

                    "Match_tier":
                        match_label(
                            m[
                                "tier"
                            ]
                        ),

                    "Match_cost":
                        f"{m['cost']:.8f}"
                }
            )


    # ========================================================
    # BED + FASTA
    # ========================================================

    selected_bg = [
        m[
            "background"
        ]
        for m in matches
    ]


    def write_bed(
        data,
        path
    ):

        with open(
            path,
            "w"
        ) as out:

            for x in data:

                out.write(
                    "\t".join(
                        [
                            x[
                                "_chrom"
                            ],
                            str(
                                x[
                                    "_start"
                                ]
                            ),
                            str(
                                x[
                                    "_end"
                                ]
                            ),
                            x[
                                "SAF_peak_id"
                            ],
                            "0",
                            "."
                        ]
                    ) +
                    "\n"
                )


    def write_fasta(
        data,
        path
    ):

        with open(
            path,
            "w"
        ) as out:

            for x in data:

                out.write(
                    ">" +
                    x[
                        "SAF_peak_id"
                    ] +
                    "|Tissue=" +
                    tissue +
                    "|Class=" +
                    x[
                        "Pair_specific_class"
                    ] +
                    "\n"
                )


                seq = x[
                    "_sequence"
                ]


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


    fg_bed = os.path.join(
        tissue_dir,
        f"{tissue}.strong.foreground.bed"
    )


    bg_bed = os.path.join(
        tissue_dir,
        f"{tissue}.matched.background.bed"
    )


    fg_fa = os.path.join(
        tissue_dir,
        f"{tissue}.strong.foreground.fa"
    )


    bg_fa = os.path.join(
        tissue_dir,
        f"{tissue}.matched.background.fa"
    )


    write_bed(
        fg_tissue,
        fg_bed
    )


    write_bed(
        selected_bg,
        bg_bed
    )


    write_fasta(
        fg_tissue,
        fg_fa
    )


    write_fasta(
        selected_bg,
        bg_fa
    )


    # ========================================================
    # CLASS COUNTS
    # ========================================================

    fg_class_counts = Counter(
        x[
            "Pair_specific_class"
        ]
        for x in fg_tissue
    )


    bg_class_counts = Counter(
        x[
            "Pair_specific_class"
        ]
        for x in selected_bg
    )


    class_file = os.path.join(
        tissue_dir,
        f"{tissue}.genomic_class_matching_QC.tsv"
    )


    with open(
        class_file,
        "w",
        newline=""
    ) as out:

        fields = [
            "Genomic_class",
            "Foreground_N",
            "Background_N",
            "Difference"
        ]


        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )


        writer.writeheader()


        for cls in CLASS_LEVELS:

            writer.writerow(
                {
                    "Genomic_class":
                        cls,

                    "Foreground_N":
                        fg_class_counts[
                            cls
                        ],

                    "Background_N":
                        bg_class_counts[
                            cls
                        ],

                    "Difference":
                        (
                            bg_class_counts[
                                cls
                            ]
                            -
                            fg_class_counts[
                                cls
                            ]
                        )
                }
            )


    # ========================================================
    # MATCHING QC
    # ========================================================

    fg_gc = [
        x[
            "_GC"
        ]
        for x in fg_tissue
    ]


    bg_gc = [
        x[
            "_GC"
        ]
        for x in selected_bg
    ]


    fg_loglen = [
        math.log2(
            x[
                "_length"
            ]
        )
        for x in fg_tissue
    ]


    bg_loglen = [
        math.log2(
            x[
                "_length"
            ]
        )
        for x in selected_bg
    ]


    tiers = Counter(
        match_label(
            x[
                "tier"
            ]
        )
        for x in matches
    )


    gc_diffs = [
        x[
            "gc_diff"
        ]
        for x in matches
    ]


    loglen_diffs = [
        x[
            "abs_log2_length_ratio"
        ]
        for x in matches
    ]


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
            tiers[
                "Strict"
            ],

        "Relaxed_matches":
            tiers[
                "Relaxed"
            ],

        "Nearest_same_class_matches":
            tiers[
                "Nearest_same_class"
            ],

        "Strict_match_percent":
            (
                100 *
                tiers[
                    "Strict"
                ] /
                len(
                    matches
                )
            ),

        "Foreground_mean_GC":
            mean(
                fg_gc
            ),

        "Background_mean_GC":
            mean(
                bg_gc
            ),

        "GC_SMD":
            smd(
                fg_gc,
                bg_gc
            ),

        "Median_abs_GC_difference":
            quantile(
                gc_diffs,
                0.50
            ),

        "Q95_abs_GC_difference":
            quantile(
                gc_diffs,
                0.95
            ),

        "Foreground_mean_log2_length":
            mean(
                fg_loglen
            ),

        "Background_mean_log2_length":
            mean(
                bg_loglen
            ),

        "Log2_length_SMD":
            smd(
                fg_loglen,
                bg_loglen
            ),

        "Median_abs_log2_length_ratio":
            quantile(
                loglen_diffs,
                0.50
            ),

        "Q95_abs_log2_length_ratio":
            quantile(
                loglen_diffs,
                0.95
            ),

        "Foreground_N_fraction_gt_0.05":
            sum(
                x[
                    "_N_fraction"
                ] > 0.05
                for x in fg_tissue
            ),

        "Background_N_fraction_gt_0.05":
            sum(
                x[
                    "_N_fraction"
                ] > 0.05
                for x in selected_bg
            )
    }


    qc_file = os.path.join(
        tissue_dir,
        f"{tissue}.matching_QC.tsv"
    )


    with open(
        qc_file,
        "w",
        newline=""
    ) as out:

        fields = list(
            qc.keys()
        )


        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )


        writer.writeheader()

        writer.writerow(
            qc
        )


    global_qc.append(
        qc
    )


    print(
        "Foreground:",
        len(
            fg_tissue
        )
    )

    print(
        "Matched background:",
        len(
            selected_bg
        )
    )

    print(
        "Strict matches:",
        tiers[
            "Strict"
        ],
        f"({qc['Strict_match_percent']:.2f}%)"
    )

    print(
        "GC SMD:",
        f"{qc['GC_SMD']:.4f}"
    )

    print(
        "log2(length) SMD:",
        f"{qc['Log2_length_SMD']:.4f}"
    )


# ============================================================
# GLOBAL QC TABLE
# ============================================================

global_qc_file = os.path.join(
    args.outdir,
    "08A_matching_QC_summary.tsv"
)


with open(
    global_qc_file,
    "w",
    newline=""
) as out:

    fields = list(
        global_qc[
            0
        ].keys()
    )


    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writeheader()

    writer.writerows(
        global_qc
    )


# ============================================================
# FINAL
# ============================================================

print(
    "\n" +
    "=" * 70
)

print(
    "STEP 08A COMPLETED"
)

print(
    "=" * 70
)

print(
    "\nPrimary QC:"
)

print(
    global_qc_file
)

print(
    "\nExact class feasibility:"
)

print(
    feasibility_file
)

print(
    "\nForeground and matched background BED/FASTA files "
    "were written under each tissue directory."
)

print(
    "\nInterpretation:"
)

print(
    "GC_SMD and Log2_length_SMD should ideally have "
    "absolute values < 0.10."
)

print(
    "Exact genomic-class counts must be identical "
    "between foreground and matched background."
)

print(
    "=" * 70
)

