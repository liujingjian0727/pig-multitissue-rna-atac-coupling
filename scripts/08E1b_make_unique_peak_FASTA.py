#!/usr/bin/env python3

from pathlib import Path
import sys

ROOT = Path("08A_matched_motif_background_v2")
OUT_ROOT = Path("08E1b_unique_peak_FASTA")

TISSUES = ["Muscle", "Spleen", "Liver", "Cerebellum"]

SETS = {
    "FG": {
        "fa_suffix": "strong.foreground.fa",
        "bed_suffix": "strong.foreground.fasta_compatible.bed",
        "out_suffix": "strong.foreground.unique.fa"
    },
    "BG": {
        "fa_suffix": "matched.background.fa",
        "bed_suffix": "matched.background.fasta_compatible.bed",
        "out_suffix": "matched.background.unique.fa"
    }
}

OUT_ROOT.mkdir(parents=True, exist_ok=True)

map_file = OUT_ROOT / "08E1b_sequence_ID_mapping.tsv"


def read_fasta(path):
    records = []

    header = None
    seq_parts = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))

                header = line[1:]
                seq_parts = []

            else:
                seq_parts.append(line.strip())

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def read_bed(path):
    rows = []

    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                raise RuntimeError(
                    f"BED has fewer than 3 columns: {path}"
                )

            rows.append(fields)

    return rows


with open(map_file, "w") as map_out:

    print(
        "\t".join([
            "Tissue",
            "Set",
            "Sequence_ID",
            "BED_chr",
            "BED_start",
            "BED_end",
            "BED_name",
            "Original_FASTA_header",
            "Sequence_length"
        ]),
        file=map_out
    )

    print("=" * 72)
    print("STEP 08E1b: rebuilding FASTA with unique sequence IDs")
    print("=" * 72)

    total = 0

    for tissue in TISSUES:

        tissue_out = OUT_ROOT / tissue
        tissue_out.mkdir(parents=True, exist_ok=True)

        for set_name, cfg in SETS.items():

            fa = (
                ROOT /
                tissue /
                f"{tissue}.{cfg['fa_suffix']}"
            )

            bed = (
                ROOT /
                tissue /
                f"{tissue}.{cfg['bed_suffix']}"
            )

            out_fa = (
                tissue_out /
                f"{tissue}.{cfg['out_suffix']}"
            )

            if not fa.exists():
                raise FileNotFoundError(
                    f"Missing FASTA: {fa}"
                )

            if not bed.exists():
                raise FileNotFoundError(
                    f"Missing BED: {bed}"
                )

            fasta_records = read_fasta(fa)
            bed_rows = read_bed(bed)

            if len(fasta_records) != len(bed_rows):
                raise RuntimeError(
                    f"{tissue} {set_name}: "
                    f"FASTA N={len(fasta_records)} "
                    f"but BED N={len(bed_rows)}"
                )

            sequence_ids = []

            with open(out_fa, "w") as out:

                for i, ((old_header, seq), bedrow) in enumerate(
                    zip(fasta_records, bed_rows),
                    start=1
                ):

                    chrom = bedrow[0]
                    start = int(bedrow[1])
                    end = int(bedrow[2])

                    bed_name = (
                        bedrow[3]
                        if len(bedrow) >= 4
                        else f"{chrom}:{start}-{end}"
                    )

                    expected_len = end - start

                    if len(seq) != expected_len:
                        raise RuntimeError(
                            f"{tissue} {set_name} row {i}: "
                            f"sequence length={len(seq)}, "
                            f"BED length={expected_len}, "
                            f"{chrom}:{start}-{end}"
                        )

                    sequence_id = (
                        f"{tissue}__{set_name}__{i:06d}"
                    )

                    sequence_ids.append(sequence_id)

                    print(
                        f">{sequence_id}",
                        file=out
                    )

                    for pos in range(0, len(seq), 80):
                        print(
                            seq[pos:pos+80],
                            file=out
                        )

                    print(
                        "\t".join([
                            tissue,
                            set_name,
                            sequence_id,
                            chrom,
                            str(start),
                            str(end),
                            bed_name,
                            old_header,
                            str(len(seq))
                        ]),
                        file=map_out
                    )

            if len(sequence_ids) != len(set(sequence_ids)):
                raise RuntimeError(
                    f"{tissue} {set_name}: duplicated Sequence_ID"
                )

            print(
                f"{tissue:12s} "
                f"{set_name:2s} "
                f"N={len(sequence_ids):4d} "
                f"PASS"
            )

            total += len(sequence_ids)

    print("-" * 72)
    print(f"Total sequences rebuilt: {total}")
    print(f"Mapping table: {map_file}")
    print("=" * 72)
    print("STEP 08E1b COMPLETED")
    print("=" * 72)

