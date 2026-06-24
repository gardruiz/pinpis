#!/bin/bash

#THIS SCRIPT TAKES A VCF FILE, A SAMPLE LIST FILE AND AN OUTPUT FILE NAME AS INPUTS, IT WILL RETURN A VCF FILE CONTAINING ONLY THE INFORMATION FROM THE SPECIFIED SAMPLES#

# Check for correct number of arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_vcf_file> <sample_list_file> <output_file>"
    exit 1
fi

# Get input arguments
#input as vcf
input_file=$1
#input as compressed vcf.gz uncoment line below
#input_file=`zcat $1 | awk '$0="$1"'`
sample_list=$2
output_file=$3
nthreads=$(nproc)

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"

# Check if the output file already exists and remove it
if [ -e "$output_file" ]; then
    rm "$output_file"
    echo "Removed existing $output_file"
fi


# Extract information for specified samples using BCFtools and .gz compressed files
"$BCFTOOLS_BIN" view --threads $nthreads -S $sample_list -Oz -o $output_file  $input_file

# Check if vcftools completed successfully
if [ $? -eq 0 ]; then
    echo "Information for specified samples was successfully extracted and saved to $output_file."
else
    echo "Error: bcftools failed to extract information for specified samples."
    exit 1
fi

