#!/usr/bin/env python3

import os
import csv
import xml.etree.ElementTree as ET

ROOT = "08B2_STREME_de_novo"
OUTDIR = "08B4_integrated_motif_evidence"

TISSUES = [
    "Muscle",
    "Spleen",
    "Liver",
    "Cerebellum"
]

os.makedirs(OUTDIR, exist_ok=True)

all_rows = []


def local_tag(tag):
    return tag.split("}")[-1]


for tissue in TISSUES:

    xml_file = os.path.join(
        ROOT,
        tissue,
        "streme.xml"
    )

    if not os.path.exists(xml_file):
        raise SystemExit(
            f"ERROR: missing {xml_file}"
        )

    tree = ET.parse(xml_file)
    root = tree.getroot()

    motifs = []

    for elem in root.iter():

        if local_tag(elem.tag) == "motif":

            row = {
                "Tissue": tissue
            }

            for k, v in elem.attrib.items():
                row[k] = v

            motifs.append(row)
            all_rows.append(row)

    print(
        tissue,
        "motifs parsed =",
        len(motifs)
    )


# union of all XML attributes
fields = ["Tissue"]

for row in all_rows:
    for key in row:
        if key not in fields:
            fields.append(key)


outfile = os.path.join(
    OUTDIR,
    "08B4A_STREME_motif_statistics_raw.tsv"
)

with open(
    outfile,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore"
    )

    writer.writeheader()
    writer.writerows(all_rows)


print()
print("Output:")
print(outfile)
print()
print("Columns:")
print("\t".join(fields))
