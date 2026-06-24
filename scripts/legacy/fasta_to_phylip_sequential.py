from Bio import SeqIO
import sys

def truncate_and_split_id(record_id):
    # Split the record ID by "_" and keep the first part
    first_part = record_id.split("_")[0]
    # Truncate to maximum of 10 characters
    truncated_id = first_part[:10]
     # Pad with spaces if the length is less than 10
    truncated_id = truncated_id.ljust(10, "X")
    return truncated_id

def fasta_to_phylip(fasta_file, phylip_file):
    # Open input FASTA file and output Phylip file
    with open(fasta_file, "r") as fasta_handle, open(phylip_file, "w") as phylip_handle:
        # Parse sequences from FASTA file
        records = list(SeqIO.parse(fasta_handle, "fasta"))
        
        # Write the number of sequences and sequence length to Phylip file
        phylip_handle.write(f"{len(records)} {len(records[0].seq)}\n")
        
        # Write sequences in Phylip format
        for record in records:
            # Truncate and split the ID
            truncated_id = truncate_and_split_id(record.id)
            # Write the truncated ID and sequence, padding ID to 20 characters
            phylip_handle.write(f"{truncated_id:<20}{record.seq}\n")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.phy")
        sys.exit(1)
    print("Converting CDS file to Phylip format") 
    fasta_file = sys.argv[1]
    phylip_file = sys.argv[2]
    fasta_to_phylip(fasta_file, phylip_file)

