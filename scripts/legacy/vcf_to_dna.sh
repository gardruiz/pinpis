#!/bin/bash
# Scripts/vcf_to_dna.sh
# Builds per-sample DNA consensus for a region and concatenates to sequences/<country>.fasta
# sample.fa  = per-sample region slice
# sample.fas = per-sample after remove_indels.py
# country.fasta = concatenation of sample.fas (in VCF sample order)

set -euo pipefail

# ===== ARGS =====
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input.vcf[.gz]> <genomic_region>" >&2
  exit 1
fi

reference="$REFERENCE_GENOME"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

vcf="$1"
genomic_region="$(printf '%s' "$2" | tr -d '\r\n[:space:]')"

# ===== VCF INDEX =====
if [ ! -e "${vcf}.tbi" ] && [ ! -e "${vcf}.csi" ]; then
  "$BCFTOOLS_BIN" index "$vcf"
fi

# ===== PATHS =====
vcf_filename=$(basename -- "$vcf")
vcf_filename_noext="${vcf_filename%%.*}"

vcf_filename_clean="$vcf_filename_noext"
vcf_filename_clean="${vcf_filename_clean/_normalized_snps/}"
vcf_filename_clean="${vcf_filename_clean/_normalized/}"
vcf_dir=$(dirname -- "$vcf")
output_dir="$vcf_dir"

if [[ "$vcf_dir" == */sequences ]]; then
  sequences_dir="$vcf_dir"
else
  sequences_dir="$vcf_dir/sequences"
fi
mkdir -p "$sequences_dir"

python_script="$script_dir/remove_indels.py"

# ===== PREP REFERENCE SLICE =====
mkdir -p tmp
region_fa="$(mktemp -p tmp region_XXXXXX.fa)"
"$SAMTOOLS_BIN" faidx "$reference" "$genomic_region" > "$region_fa" || exit 3

# ===== CAPTURE SAMPLE LIST ONCE =====
mapfile -t samples < <("$BCFTOOLS_BIN" query -l "$vcf")

# ===== PER-SAMPLE LOOP =====
for sample in "${samples[@]}"; do
  [ -n "$sample" ] || continue
  fa_out="$output_dir/${sample}.fa"
  fas_out="$output_dir/${sample}.fas"

  
  tmp_consensus_full="$(mktemp -p tmp cons_full_XXXXXX.fa)"

  # Generate consensus sequence
  if ! "$BCFTOOLS_BIN" consensus -H A --mark-ins lc --mark-del '*' \
      -s "$sample" -f "$reference" "$vcf" > "$tmp_consensus_full"; then

    echo "[ERROR] Consensus reconstruction failed: sample=$sample region=$genomic_region" >&2
    rm -f "$tmp_consensus_full" "$fa_out" "$fas_out"
    exit 1
  fi

  # Extract genomic region
  if ! "$SAMTOOLS_BIN" faidx "$tmp_consensus_full" "$genomic_region" > "$fa_out"; then

    echo "[ERROR] Region extraction failed: sample=$sample region=$genomic_region" >&2
    rm -f "$tmp_consensus_full" "$tmp_consensus_full.fai" "$fa_out" "$fas_out"
    exit 1
  fi

  # Clean up temporary files
  rm -f "$tmp_consensus_full" "$tmp_consensus_full.fai"


  if head -n1 "$fa_out" | grep -q '^>'; then
    sed -i "1s|^>.*$|>$sample|" "$fa_out"
  else
    tmpfa="$(mktemp)"
    { printf ">%s\n" "$sample"; cat "$fa_out"; } > "$tmpfa"
    mv "$tmpfa" "$fa_out"
  fi

  if [ -f "$python_script" ]; then
    ref_content="$(cat "$region_fa")"
    "$PYTHON_BIN" "$python_script" "$fa_out" "$ref_content" "$fas_out"
  else
    cp "$fa_out" "$fas_out"
  fi
done

# ===== CONCAT: sequences/<country>.fasta using ONLY sample.fas =====
output_combined="$sequences_dir/${vcf_filename_clean}.fasta"
tmp_combined="$(mktemp)"
: > "$tmp_combined"

for sample in "${samples[@]}"; do
  f="$output_dir/${sample}.fas"
  [ -s "$f" ] && cat "$f" >> "$tmp_combined"
done

mv -f "$tmp_combined" "$output_combined"

