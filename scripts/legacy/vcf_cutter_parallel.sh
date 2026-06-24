#!/bin/bash

if [ $# -ne 5 ]; then
  echo "Usage: $0 <gene_names_file> <sample_list_directory> <vcf_file> <GFF_file> <use_gene_name>" 
  exit 1
fi

gene_locations_file="$1"
sample_list_parent_dir="$2"
input_file="$3"
gff_file="$4"
use_gene_name="$5"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
TABIX_BIN="${TABIX_BIN:-tabix}"

gene_locations_output_directory=$(dirname "$gene_locations_file")
gene_locations_output_file="${gene_locations_file%.*}_locations.txt"
num_threads=$(( $(nproc) / 2 ))

# Create output directory
mkdir -p "$gene_locations_output_directory"

# Check if the output file exists and handle it
if [ -e "$gene_locations_output_file" ]; then
    echo "File $gene_locations_output_file already exists. Replacing it..."
    rm "$gene_locations_output_file"
fi


process_gene() {
    word="$1"
    sample_list_parent_dir="$2"
    input_file="$3"
    gff_file="$4"
    use_gene_name="$5"
    gene_locations_output_file="$6"
    num_threads="$7"

    if [ "$use_gene_name" = "True" ]; then
        result=$(grep "Name=$word" "$gff_file" | grep gene | cut -f1,4,5)
    else
        result=$(grep "ID=$word" "$gff_file" | grep gene | cut -f1,4,5)
    fi

    # Output gene locations
    while read line; do
        echo -e "$word\t$line" >> "$gene_locations_output_file"
    done <<< "$result"

    while IFS=$'\t' read -r field1 field2 field3; do
        region="$field1:$field2-$field3"

	# Create the output file directory
	vcf_output_directory="variants/$word"
        mkdir -p "$vcf_output_directory"

	# Construct the outout file name
        vcf_output_file="$vcf_output_directory/$word.vcf.gz"

	if [ ! -f "$vcf_output_file" ]; then
        echo "Creating gene-level VCF for $word in region $region"
        
	    if ! "$BCFTOOLS_BIN" view --threads "$num_threads" -r "$region" "$input_file" -Oz -o "$vcf_output_file"; then
                echo "bcftools failed for gene: $word"
                return 1
            fi
	    
	    "$TABIX_BIN" -p vcf "$vcf_output_file"
        else
            echo "Gene-level VCF for $word already exists. Skipping extraction."
        fi


        # For each sample subdirectory
        for sample_list_dir in "$sample_list_parent_dir"/*/; do
            [ -d "$sample_list_dir" ] || continue
            samples_dir=$(basename "$sample_list_dir")

            for sample_list_file in "$sample_list_dir"/*.txt; do
                sample_name=$(basename "$sample_list_file" .txt)
                output_dir="$vcf_output_directory/$samples_dir"
                mkdir -p "$output_dir"
                output_file="$output_dir/$sample_name.vcf.gz"
		
                "$script_dir/vcf_by_sample.sh" "$vcf_output_file" "$sample_list_file" "$output_file"
                echo "Processed $sample_list_file for gene $word"
            done
        done
    done <<< "$result"
}

export -f process_gene
export script_dir
export BCFTOOLS_BIN
export TABIX_BIN

# Use GNU Parallel if available
if command -v parallel &> /dev/null; then
    parallel process_gene {} "$sample_list_parent_dir" "$input_file" "$gff_file" "$use_gene_name" "$gene_locations_output_file" "$num_threads" :::: "$gene_locations_file"
    wait
else
    echo "GNU Parallel not found. Falling back to xargs..."
    cat "$gene_locations_file" | xargs -n 1 -P "$num_threads" -I {} bash -c \
    'process_gene "$@"' _ {} "$sample_list_parent_dir" "$input_file" "$gff_file" "$use_gene_name" "$gene_locations_output_file" "$num_threads"
fi

