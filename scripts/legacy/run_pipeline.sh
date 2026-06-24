#!/bin/bash


# Function to read the control file and export variables
parse_control_file() {
    local control_file="Config/control_file.txt"
    while IFS='=' read -r key value; do
        # Strip leading/trailing whitespace and export variable
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs)
        export "$key=$value"
    done < "$control_file"
}


run_command() {
    echo "Running command: $*"
    "$@"
    if [ $? -ne 0 ]; then
        echo "Error: Command failed: $*"
        return 1
    fi
}



control_file="Config/control_file.txt"

# Parse the control file
parse_control_file "$control_file"

# Check if use_gene_name is set, if not default to using IDs
if [ -z "$use_gene_name" ]; then
    echo "use_gene_name is not set, defaulting to using IDs."
    use_gene_name="False"
fi

# initialize variant sites file with header
echo -e "gene_name\tpopulation\tvariant_sites" >> "results/variant_sites.txt"

# Check if there is a previos pMK test file when performin a pMK test
if [ "$pMK" = "True" ]; then
    echo "pMK test enabled."
                    
    # Ensure the results/pmk.txt file exists
    pmk_results_file="results/pmk.txt"
    if [ -f "$pmk_results_file" ]; then
            rm "$pmk_results_file"
            touch "$pmk_results_file" 
	    echo -e "Gene\tA\tB\tC1\tC2\tD1\tD2\tpMKContinet_1\tpMKContinent_2" > "$pmk_results_file"  # Add header
    else 
	    echo -e "Gene\tA\tB\tC1\tC2\tD1\tD2\tpMKConrinet_1\tpMKContinent_2" > "$pmk_results_file"  # Add header
    fi
fi

echo "Extracting samples ids from: $Populations_1_file"
# Run the  script and check for success
if ! run_command python3.8 Scripts/multysampler.py $Populations_1_file $Number_of_samples $Samples_ID_file $Years_file $Samples_1_output_dir; then
    echo "Error: Failed to extract sample ids"
fi

echo "Extracting samples ids from: $Populations_2_file"
# Run the  script and check for success
if ! run_command python3.8 Scripts/multysampler.py $Populations_2_file $Number_of_samples $Samples_ID_file $Years_file $Samples_2_output_dir; then
    echo "Error: Failed to extract sample ids"
fi


echo "Cutting VCF into individual gene files"
# Run the vcf cutter  script and check for success
if ! Scripts/vcf_cutter_parallel.sh $Gene_names_file $Samples_file_dir $VCF_file $GFF_file $use_gene_name; then
    echo "Error: Failed to extract genes VCF files"
fi



# Check which pipelines to run
IFS=' ' read -r -a pipelines <<< "$Pipelines"

# Parse gene locations
gene_locations_file="${Gene_names_file%.*}_locations.txt"

for pipeline in "${pipelines[@]}"; do
    case "$pipeline" in
	dnds)
	
	    echo "Starting dN/dS pipeline"

            # Function to calculate dN/dS for each gene

	    extract_cds() {
    	        gene_name="$1"
    		echo "Extracting CDS for gene: $gene_name"

    		gene_dir="variants/$gene_name"
    		subdirs=($(find "$gene_dir" -mindepth 1 -maxdepth 1 -type d))
    		Pop1_dir="${subdirs[0]}"
    		Pop2_dir="${subdirs[1]}"

 	        if [ "$use_gene_name" = "True" ]; then
        	    gene_info=$(grep "Name=$gene_name" "$GFF_file")
                    gene_ID=$(echo "$gene_info" | awk -F 'ID=' '{print $2}' | awk -F ';' '{print $1}')
   	 	else
        	    gene_info=$(grep "ID=$gene_name" "$GFF_file"|grep gene)
        	    gene_ID=$gene_name
    		fi

    		Genomic_region=$(echo "$gene_info" | awk -F '\t' '{print $1":"$4"-"$5}' | head -n1)
    		genome_dir="Genome/Genes"
    		BED_file="$genome_dir/$gene_ID.bed"

    		run_command Scripts/bed_extractor.sh "$gene_ID" "$GFF_file" "$genome_dir"
		
	       	# Extract CDS from populations dir
                if ! run_command Scripts/vcf_to_cds.sh "$Pop1_dir" "$Genomic_region" "$BED_file"; then
                    echo "Error: Failed to extract CDS for population 1 of $gene_name! $Pop1_dir $Genomic_region Skipping this gene."
                    return  # Skip to the next gene
                fi

                if ! run_command Scripts/vcf_to_cds.sh "$Pop2_dir" "$Genomic_region" "$BED_file"; then
                    echo "Error: Failed to extract CDS for population 2 of $gene_name! $Pop2_dir $Genomic_region  Skipping this gene."
                    return  # Skip to the next gene
                fi


    		#run_command "Scripts/vcf_to_cds.sh $Pop1_dir $Genomic_region $BED_file"
    		#run_command "Scripts/vcf_to_cds.sh $Pop2_dir $Genomic_region $BED_file"
	    }


	    calculate_dnds() {
		gene_name="$1"
                echo "Processing gene: $gene_name"
		
		# Identify subdirectories inside the gene directory
                gene_dir="variants/$gene_name"
                subdirs=($(find "$gene_dir" -mindepth 1 -maxdepth 1 -type d))

                if [ "${#subdirs[@]}" -lt 2 ]; then
                    echo "Error: Less than 2 subdirectories found in $gene_dir!"
                    exit 1
                fi

                # Dynamically assign Pop1_files and Pop2_files from the first two subdirectories
                Pop1_dir="${subdirs[0]}"
                Pop2_dir="${subdirs[1]}"

		# Extract the ID from the GFF file
		if [ "$use_gene_name" = "True" ]; then
		    echo "Using gene names."
		    gene_info=$(grep "Name=$gene_name" "$GFF_file")
		    # Extract the ID from the gene_info using awk
                    gene_ID=$(echo "$gene_info" | awk -F 'ID=' '{print $2}' | awk -F ';' '{print $1}')
                    echo "Gene ID: $gene_ID"

	        else 
		    echo "Using IDs names."
                    gene_info=$(grep "ID=$gene_name" "$GFF_file"|grep gene)
		    gene_ID=$gene_name
		fi

		# Extract the genomic region (chromosome, start, and end) using awk
		Genomic_region=$(echo "$gene_info" | awk -F '\t' '{print $1":"$4"-"$5}' | head -n1)
		genome_dir="Genome/Genes"
		BED_file="$genome_dir/$gene_ID.bed"

		# Extract samples ids parent directory
		sample_list_dir=$(dirname "$Samples_1_output_dir")
    		
		# Calculate dnds
    		echo "Calculating dN/dS"
		

		# Create gene-specific output file
                output_dnds_file="results/${gene_name}_dnds.txt"

       		if [ -f "$output_dnds_file" ]; then
                    echo "File $output_dnds_file exists, replacing it."
                    rm "$output_dnds_file"  # Remove the existing file
                fi

    		if ! run_command Scripts/dnds2.sh "$Pop1_dir" "$Pop2_dir" "$output_dnds_file" "$gene_name"; then
        	    echo "Error: Failed to calculate dN/dS for $gene_name! Skipping this gene."
                    return  # Skip to the next gene
    		fi
	 
            }
	
	    pmk() {
	       gene_name="$1"

               # Clean YN00 output
    	       if ! Scripts/clean_yn00_output.sh "results/$gene_name"; then
                  echo "Error cleaning YN00 output for $gene_name"
                  return
               fi

               wait
	       if [ "$pMK" = "True" ]; then
                   echo "Running pMK test for $gene_name"

                   pmk_output=$(python3.8 Scripts/pmk.py "results/$gene_name/clean" "$Populations_1_file" "$Populations_2_file")
                   if [ $? -ne 0 ]; then
                      echo "Error running pMK test for $gene_name"
                      return
	   	   else
		      # Append the gene name and pMK output to the results file
                      echo -e "$gene_name\t$pmk_output" >> "$pmk_results_file"
		      return 1
                   fi
		fi
            } 
	    count_variant_sites() {
    	        local gene_name=$1
    		local folder="variants/$gene_name"
    		local output_file="results/variant_sites.txt"

    		mkdir -p tmp results

    		# Add header if file doesn't exist
    		if [ ! -f "$output_file" ]; then
        	    echo -e "gene_name\tpopulation\tvariant_sites\tlevel" > "$output_file"
    		fi

    		# === 1. Per-file logic ===
    		for fasta_file in "$folder"/*/sequences/*_CDS.fasta; do
        	    echo "Processing individual file: $fasta_file"
        	    filename=$(basename "$fasta_file")
        	    population="${filename%%_*}"

        	variable_sites=$(python3.8 - <<EOF
from Bio import AlignIO
alignment = AlignIO.read("$fasta_file", "fasta")
count = 0
for i in range(alignment.get_alignment_length()):
    col = alignment[:, i].upper()
    bases = [b for b in col if b in "ACGT"]
    if len(set(bases)) > 1:
        count += 1
print(count)
EOF
)
        	    echo -e "$gene_name\t$population\t$variable_sites\tfile" >> "$output_file"
    	        done

    	        # === 2. Per-population logic ===
    		for pop_dir in "$folder"/*/; do
        	    if [ -d "$pop_dir" ]; then
            		pop_name=$(basename "$pop_dir")
            		echo "Processing population: $pop_name"
            		seq_dir="$pop_dir/sequences"
            		concat_fasta="tmp/${gene_name}_${pop_name}.fasta"

            		cat "$seq_dir"/*.fasta > "$concat_fasta"

            		variable_sites=$(python3.8 - <<EOF
from Bio import AlignIO
alignment = AlignIO.read("$concat_fasta", "fasta")
count = 0
for i in range(alignment.get_alignment_length()):
    col = alignment[:, i].upper()
    bases = [b for b in col if b in "ACGT"]
    if len(set(bases)) > 1:
        count += 1
print(count)
EOF
)
            		echo -e "$gene_name\t$pop_name\t$variable_sites\tpopulation" >> "$output_file"
        	    fi
                done
	    }


            
	    has_variants_fast() {
	        local gene_name=$1
   		local vcf_file="variants/$gene_name/$gene_name.vcf.gz"
    
                if [ ! -s "$vcf_file" ]; then
                    echo "Skipping $gene_name: Missing or empty VCF" >&2
                    return 1
                fi

                  if bcftools view -H "$vcf_file" | head -n 1 | grep -q .; then
                      return 0
                  else
                      echo "Skipping $gene_name: No variants" >&2
                      return 1
                  fi
            }
	    
	    # Export necessary functions and variables
	    export -f extract_cds
	    export -f calculate_dnds
            export -f run_command
            export gene_ID GFF_file pMK $pmk_output gene_locations_file pmk_results_file
            export -f has_variants_fast
            export -f pmk
	    export -f count_variant_sites

            # Phase 1: Extract CDS
	    echo "Phase 1: Extracting CDS"
	    cat "$gene_locations_file" | awk '{print $1}' | \
                parallel --line-buffer -j "$(nproc)" 'extract_cds {}'
	    wait

	    # Phase 2: dN/dS calculation
            echo "Phase 2: Calculating dN/dS"
            cat "$gene_locations_file" | awk '{print $1}' | \
                parallel --line-buffer -j $(nproc) '
                    if has_variants_fast {}; then
                        calculate_dnds {}
                    fi
              '
	    wait
	    #  Phase 3: pMK calculation
            echo "Phase 3: Calculating pMK"
            cat "$gene_locations_file" | awk '{print $1}' | \
               parallel --line-buffer -j $(nproc) '
                    if has_variants_fast {}; then
                        pmk {}
                    fi
                '
	    wait
            #  Phase 4: variant site calculation
            echo "Phase 4: Calculating variant sites"
            cat "$gene_locations_file" | awk '{print $1}' | \
               parallel --line-buffer -j $(nproc) '
                    if has_variants_fast {}; then
                        count_variant_sites {}
                    fi
                '

        ;;
    esac
done


echo "Pipeline executed successfully."
