#!/usr/bin/env bash
# vcf_to_cds.sh — per-country normalization + DNA consensus + CDS extraction (no filtering)

set -euo pipefail
IFS=$'\n\t'

# --- Args ---
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <vcf_dir> <genomic_region> <bed_file>" >&2
  exit 1
fi

vcf_dir="$1"
genomic_region_raw="$2"
bed_file="$3"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REFERENCE_GENOME="$REFERENCE_GENOME"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
BGZIP_BIN="${BGZIP_BIN:-bgzip}"
TABIX_BIN="${TABIX_BIN:-tabix}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

# --- Paths ---
sequences_dir="$vcf_dir/sequences"
mkdir -p "$sequences_dir" tmp

dna_script="$script_dir/vcf_to_dna.sh"
cds_script="$script_dir/CDS.py"

# --- Sanitize region (avoid stray spaces/newlines/CR) ---
genomic_region="$(printf "%s" "$genomic_region_raw" | tr -d '\r' | sed 's/[[:space:]]//g')"

# --- Loop all country VCFs ---
shopt -s nullglob
cds_created=0
for vcf_file in "$vcf_dir"/*.vcf.gz; do
  [[ -f "$vcf_file" ]] || continue

  # Country name (e.g. benin from benin.vcf.gz)
  country="$(basename "${vcf_file%%.vcf.gz}")"

  echo "[INFO] Processing $country (VCF: $(basename "$vcf_file"))"
  variants_total=$("$BCFTOOLS_BIN" view -H "$vcf_file" | wc -l || true)
  echo "[INFO] Raw variant rows in $country VCF: $variants_total"

  # --- Skip filtering: use original VCF directly ---
  unfiltered_vcf="$vcf_file"
  copied_vcf="$sequences_dir/${country}_filtered.vcf.gz"
  cp "$unfiltered_vcf" "$copied_vcf"

  
  normalised_vcf_file="$sequences_dir/${country}_normalized_snps.vcf.gz"

  echo "[INFO] Normalizing and keeping SNP-only records for $country"

  "$BCFTOOLS_BIN" norm \
    -m -any \
    -f "$REFERENCE_GENOME" \
    "$copied_vcf" \
    -Ou \
  | "$BCFTOOLS_BIN" view \
    -v snps \
    -Oz \
    -o "$normalised_vcf_file"

  "$TABIX_BIN" -f -p vcf "$normalised_vcf_file"

  compressed_vcf_file="$normalised_vcf_file"

  variants_clean=$("$BCFTOOLS_BIN" view -H "$compressed_vcf_file" | wc -l || true)
  echo "[INFO] SNP rows after normalization/filtering ($country): $variants_clean"
  
  

  # === Consensus DNA per sample ===
  echo "[INFO] Generating per-sample DNA consensus for $country in region $genomic_region"
  DEBUG=0 REFERENCE_GENOME="$REFERENCE_GENOME" BCFTOOLS_BIN="$BCFTOOLS_BIN" SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}" PYTHON_BIN="$PYTHON_BIN" bash "$dna_script" "$compressed_vcf_file" "$genomic_region"
  
  # Combined DNA file the DNA script writes (country.fasta)
  dna_output="$sequences_dir/${country}.fasta"
  if [[ ! -s "$dna_output" ]]; then
    echo "[ERROR] DNA combined fasta missing for $country: $dna_output"
    rm -f "$copied_vcf" "$compressed_vcf_file" "${compressed_vcf_file}.tbi" || true
    continue
  fi

  # === CDS from DNA ===
  cds_file="$sequences_dir/${country}_CDS.fasta"
  echo "[INFO] Translating DNA to CDS for $country → $(basename "$cds_file")"
  "$PYTHON_BIN" "$cds_script" "$dna_output" "$bed_file" "$cds_file"
  if [[ -s "$cds_file" ]]; then
    cds_created=$((cds_created + 1))
  else
    echo "[ERROR] CDS file missing or empty for $country: $cds_file" >&2
  fi
  # --- Optional cleanup (currently disabled) ---
  # rm -f "$copied_vcf" "$compressed_vcf_file" "${compressed_vcf_file}.tbi" "$dna_output" || true

done

if [[ "$cds_created" -eq 0 ]]; then
  echo "[ERROR] No CDS FASTA files were created for directory: $vcf_dir" >&2
  exit 1
fi

shopt -u nullglob

