# Paired multi-tissue RNA-seq and ATAC-seq integration in pig

This repository contains custom scripts, metadata templates, workflow notes, and audit records associated with the manuscript:

"Paired multi-tissue RNA-seq and ATAC-seq integration identifies tissue-specific chromatin-transcription coupling in the pig"

The study reanalyzed publicly available paired RNA-seq and ATAC-seq datasets from two pigs across eight tissues: adipose, cerebellum, cerebral cortex, hypothalamus, liver, lung, skeletal muscle, and spleen.

## Scope of this repository

This repository is intended to document the computational workflow and support reproducibility of the manuscript-level analyses. It includes:

- custom R, Python, and shell scripts used for quality control, tissue-specificity analysis, RNA-ATAC integration, permutation testing, promoter-first annotation, motif-support summarization, figure generation, supplementary-table construction, and manuscript-level numeric auditing;
- sample-pair metadata and public accession information;
- software-version records;
- locked numeric-audit records used to verify manuscript consistency;
- figure and supplementary-table source manifests.

Large raw sequencing files, BAM files, bigWig files, and other large intermediate files are not stored in this repository. Raw sequencing data are publicly available from GEO/SRA under the accessions listed in `metadata/public_accessions.tsv`.

## Main data sources

RNA-seq:

- GEO: GSE158412
- BioProject: PRJNA665193

ATAC-seq:

- GEO: GSE158414
- BioProject: PRJNA665194

SuperSeries:

- GSE158430

Reference assembly:

- Sus scrofa Sscrofa11.1
- NCBI assembly: GCF_000003025.6
- Ensembl release 115

## Interpretation note

Peak-gene associations, motif matches, candidate transcription factors, and network modules are computationally prioritized hypotheses. They should not be interpreted as experimentally validated regulatory interactions, transcription-factor occupancy, or causal enhancer-gene relationships.

## Citation

A formal citation and archived DOI will be added after repository release.
