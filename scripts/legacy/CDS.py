#!usr/bin/env python3

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
import pandas as pd




### This script takes a DNA sequence, a gene name and a BED file as input to produce the CDS. The input sequence should be in fasta format ###



def reverse_complement(sequence):
    complement = {'*': '*', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'r': 'y', 'y': 'r', 'w': 'w', 's': 's', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b'}
    reversed_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement[base] for base in reversed_sequence)
    return reverse_complement_sequence


def read_fasta(file):
   sequences=[]
   with open(file, "r") as fasta_file:
       for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record)
   return sequences


def read_bed(file):
    cds_positions =[]
    gene_start=''
    sence=''
    bed_df=pd.read_csv(file, sep=" ", header=None, names=['Chr', "element", 'start', 'end', 'sense'])        
    # Extract gene start and sense direction
    gene_row = bed_df[bed_df['element'] == 'gene']
    if not gene_row.empty:
        gene_start = int(gene_row['start'].values[0])
        sense = gene_row['sense'].values[0]
    else:
        raise ValueError("No 'gene' row found in the BED file.")

    for _, row in bed_df.iterrows():
        if row['element'] != 'gene':
            start=int(row['start'])-gene_start
            end=int(row['end'])-gene_start
            cds_positions.append([start,end])
    return cds_positions, sense


def get_cds(sequence, positions, sence):
    cds = ""
    for start, end in positions:
        cds += sequence[start:end + 1]
    if sence == "-":
        cds = reverse_complement(cds)

    stop_codons = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}  # Set of stop codons

    # Check if the last codon is a stop codon
    last_codon = cds[-3:]  
    if last_codon in stop_codons:
        cds = cds[:-3]  # Remove the last codon

    # Check for premature stop codons directly
    for i in range(0, len(cds) - 2, 3):  # Iterate over the CDS in steps of 3
        codon = cds[i:i + 3]  # Extract codon
        if codon in stop_codons:
            print("Stop codon found, excluding sequence")
            return None  # Found a premature stop codon

    return cds if cds else None

if len(sys.argv) != 4:
    print("Usage: python CDS.py fasta_file bed_file output_file")
    exit(1)

file=sys.argv[1]
sequences=read_fasta(file)
cds_positions, sence = read_bed(sys.argv[2])

# Initialize an empty list to store SeqRecord objecs
cds_records = []

# Generate the sequence ID

for sequence in sequences:
    header = sequence.description  # Get the header from the original FASTA record
    seq_id = header.split()[0]  # Use the first part of the header as the seq_id
    print(f"Obtaining CDS from sample: {header}")
    cds = get_cds(sequence.seq, cds_positions, sence)
   
    file_name_without_extension = os.path.basename(file).split('.')[0]
    if cds is None:
        print(f"Premature stop codon detected in sequence ID {seq_id}; skipping this sequence.")
        continue
    # Extract the file name without extension
    cds_record   = SeqRecord(Seq(cds), id=seq_id)
    
    # Add the SeqRecord to the list
    cds_records.append(cds_record)
    
output_file=sys.argv[3]
SeqIO.write(cds_records, output_file, "fasta")
print(f"CDS sequences saved to {output_file}")
