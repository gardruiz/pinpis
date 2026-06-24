#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 pop1_fasta_files_dir pop2_fasta_files_dir output_file gene_name" >&2
    exit 1
fi

echo "Estimating dN/dS"

pop1_directory="$1"
pop2_directory="$2"
output_file="$3"
gene_name="$4"

# Directory where this script lives.
# This lets the script find fasta_to_phylip_sequential.py inside scripts/legacy/
# instead of relying on the old hard-coded Scripts/ directory.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Allow the Python package to pass specific tools through the environment.
PYTHON_BIN="${PYTHON_BIN:-python3}"
YN00_BIN="${YN00_BIN:-yn00}"

INCLUDE_BETWEEN_GROUPS="${INCLUDE_BETWEEN_GROUPS:-false}"

phylip_script="$script_dir/fasta_to_phylip_sequential.py"

if [ ! -f "$phylip_script" ]; then
    echo "Error: FASTA-to-PHYLIP script not found: $phylip_script" >&2
    exit 1
fi

# Find CDS FASTA files.
mapfile -t pop1_files < <(find "$pop1_directory" -type f -name "*_CDS.fasta" | sort)
mapfile -t pop2_files < <(find "$pop2_directory" -type f -name "*_CDS.fasta" | sort)

if [ "${#pop1_files[@]}" -eq 0 ]; then
    echo "Error: no CDS FASTA files found in $pop1_directory" >&2
    exit 1
fi

if [ "${#pop2_files[@]}" -eq 0 ]; then
    echo "Error: no CDS FASTA files found in $pop2_directory" >&2
    exit 1
fi

# Main result folders.
results_folder="results/$gene_name"
tmp_folder="$results_folder/tmp"

mkdir -p "$results_folder" "$tmp_folder"
mkdir -p "$(dirname "$output_file")"

# Optional: add a header if the output file does not exist.
if [ ! -f "$output_file" ]; then
    echo -e "population_1\tpopulation_2\tdN\tdS\tdN_dS" > "$output_file"
fi

run_dnds() {
    file1="$1"
    file2="$2"

    uuid=$(uuidgen)

    concatenated_base="align_${uuid}.fasta"
    yn00_input_base="align_${uuid}.phy"
    yn00_output_base="yn00_output_${uuid}.txt"
    control_base="yn00_control_file_${uuid}.ctl"

    concatenated_file="$tmp_folder/$concatenated_base"
    yn00_input_file="$tmp_folder/$yn00_input_base"
    yn00_output_file="$tmp_folder/$yn00_output_base"
    control_file="$tmp_folder/$control_base"

    country1=$(basename "$file1" | sed 's/_CDS\.fasta$//')
    country2=$(basename "$file2" | sed 's/_CDS\.fasta$//')

    echo "[dnds] Running pair: $country1 vs $country2"

    # If the same file is compared to itself, do not duplicate the sequences.
    if [ "$file1" = "$file2" ]; then
        cp "$file1" "$concatenated_file"
    else
        cat "$file1" "$file2" > "$concatenated_file"
    fi

    # Convert FASTA to PHYLIP.
    "$PYTHON_BIN" "$phylip_script" "$concatenated_file" "$yn00_input_file"

    if [ ! -s "$yn00_input_file" ]; then
        echo "Error: PHYLIP file was not created or is empty: $yn00_input_file" >&2
        return 1
    fi

    # IMPORTANT:
    # The control file uses basenames, because yn00 will be run from inside tmp_folder.
    # This keeps PAML auxiliary files out of the project root.
    cat > "$control_file" <<EOL
seqfile = $yn00_input_base
outfile = $yn00_output_base

verbose = 1
icode = 0
weighting = 0
commonf3x4 = 0
EOL

    (
        cd "$tmp_folder"
        "$YN00_BIN" "$control_base"
    )

    if [ ! -s "$yn00_output_file" ]; then
        echo "Error: yn00 output was not created or is empty: $yn00_output_file" >&2
        return 1
    fi

    tmp_output="${results_folder}/tmp_${uuid}_dnds.txt"
    summary_output="${tmp_output}_summary.txt"

    # Parse YN00 output.
    awk '/seq\. seq./,/LWL85/ {
        if ($0 ~ /LWL85|LPB93|LWLm/) next;

        if ($8 != "" && $11 != "") {
            sumdN += $8;
            sumdS += $11;
            count++;
        }
    }
    END {
        if (count > 0) {
            avg_dN = sumdN / count;
            avg_dS = sumdS / count;
            dnds = (avg_dS > 0 ? avg_dN / avg_dS : "NA");

            print "Average dN:", avg_dN > "'"$tmp_output"'";
            print "Average dS:", avg_dS > "'"$tmp_output"'";
            print "dN/dS ratio:", dnds > "'"$tmp_output"'";

            print avg_dN, avg_dS, dnds;
        } else {
            print "NA\tNA\tNA";
        }
    }' "$yn00_output_file" | tail -n1 > "$summary_output"

    awk -v c1="$country1" -v c2="$country2" '
        {print c1"\t"c2"\t"$1"\t"$2"\t"$3}
    ' "$summary_output" >> "$output_file"

    # Save full yn00 output.
    final_out="$results_folder/${country1}-${country2}.out"
    mv "$yn00_output_file" "$final_out"

    # Clean PAML auxiliary files written inside tmp_folder.
    rm -f "$tmp_folder"/rst* \
          "$tmp_folder"/2YN.* \
          "$tmp_folder"/rub \
          "$tmp_folder"/lnf \
          "$tmp_folder"/yn00.ctl \
          2>/dev/null || true

    # Clean pair-specific temporary files.
    rm -f "$concatenated_file" \
          "$yn00_input_file" \
          "$control_file" \
          "$summary_output"
}

# Within population 1.
for ((i = 0; i < ${#pop1_files[@]}; i++)); do
    run_dnds "${pop1_files[i]}" "${pop1_files[i]}"

    for ((j = i + 1; j < ${#pop1_files[@]}; j++)); do
        run_dnds "${pop1_files[i]}" "${pop1_files[j]}"
    done
done

# Within population 2.
for ((i = 0; i < ${#pop2_files[@]}; i++)); do
    run_dnds "${pop2_files[i]}" "${pop2_files[i]}"

    for ((j = i + 1; j < ${#pop2_files[@]}; j++)); do
        run_dnds "${pop2_files[i]}" "${pop2_files[j]}"
    done
done

# Between population 1 and population 2.
if [[ "$INCLUDE_BETWEEN_GROUPS" == "true" || "$INCLUDE_BETWEEN_GROUPS" == "True" || "$INCLUDE_BETWEEN_GROUPS" == "1" ]]; then
    echo "Running between-group dN/dS comparisons."

    for file1 in "${pop1_files[@]}"; do
        for file2 in "${pop2_files[@]}"; do
            run_dnds "$file1" "$file2"
        done
    done
else
    echo "Skipping between-group dN/dS comparisons."
fi

echo "dN/dS estimation completed for $gene_name"
