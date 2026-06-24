# PinPis

PinPis is an internal reproducible pipeline for estimating population-level coding-sequence dN/dS or Pin/Pis summaries from VCF data.

This pipeline was developed for *Plasmodium falciparum* population genomic analyses using Pf7 data. It reconstructs coding sequences from population VCFs, estimates pairwise dN/dS using PAML `yn00`, counts population-level variant sites, and produces a clean gene-by-population result table.

## Current status

This is a wrapped legacy analysis pipeline intended for reproducible internal use. It is not yet intended as a general public software package.

The public-facing methodological package is the simulation/null-calibration component, not this full VCF-to-CDS pipeline.

## Main workflow

```text
sample IDs
→ cut VCF by gene and population group
→ reconstruct CDS FASTA files
→ estimate dN/dS with yn00
→ count population-level variant sites
→ merge final result table
