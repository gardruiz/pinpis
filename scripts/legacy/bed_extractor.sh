#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gene_ID> <gff_file> <output_dir>"
    exit 1
fi

# Assign command-line arguments to variables
gene_ID="$1"
gtf_file="$2"
output_dir="$3"

# Generate BED file name based on the gene ID
bed_file="$output_dir/${gene_ID}.bed"

# Extract coding regions of the specific gene from GTF
awk -v gene="$gene_ID" '$3 == "CDS" && $9 ~ gene|| $3 == "gene" && $9 ~ gene {print $1, $3, $4, $5, $7}' "$gtf_file" > "$bed_file"

echo "BED file '$bed_file' created."



