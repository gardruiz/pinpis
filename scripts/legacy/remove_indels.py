import os
import sys

def remove_insertions(consensus_file, temp_file):
    print("Removing insertions")
    with open(consensus_file, 'r') as infile:
        lines = infile.readlines()

    processed_lines = []
    for line in lines:
        if line.startswith('>'):
            processed_lines.append(line)
        else:
            processed_line = ''.join([char for char in line if not char.islower()])
            processed_lines.append(processed_line)

    with open(temp_file, 'w') as outfile:
        outfile.writelines(processed_lines)


def add_deletions(temp_file, reference, header, output_file):
    print("Ading deletions using reference file")
    ref_lines = reference.split('\n')
    
    with open(temp_file, 'r') as cons_file:
        cons_lines = cons_file.readlines()

    header
    ref_seq = ''
    cons_seq = ''
    
    # Extract headers and sequences
    for line in ref_lines:
        if not line.startswith('>'):
            ref_seq += line.strip()

    for line in cons_lines:
        if not line.startswith('>'):
            cons_seq += line.strip()
    
    print("Consensus sequence length: ", len(cons_seq))
    print("Reference sequence length: ", len(ref_seq))
    # Process sequences
    final_seq = ''.join([ref_seq[i] if char == '*' else char for i, char in enumerate(cons_seq)])
    print("Final sequencelength: ", len(final_seq))
    # Format as FASTA string
    final_fasta = ">"+ header + '\n' + final_seq + '\n'

    with open(output_file, 'w') as outfile:

        outfile.write(final_fasta)


if __name__ == "__main__":
    
    print(os.getcwd())  # Check script's working directory
    print("Removing indels")
    consensus_file = sys.argv[1] # Input consensus sequence file
    reference = sys.argv[2]   # Reference sequence file
    output_file = sys.argv[3]      # Output final consensus sequence file
    pid = os.getpid()  # Get process ID
    temp_file = f'./tmp/temp_no_insertions{pid}.fa'  # Temporary file to store intermediate result
    filename = consensus_file.split('/')[-1]
    header = filename.split('.')[0]
    remove_insertions(consensus_file, temp_file)
    add_deletions(temp_file, reference, header, output_file)
  
    # Clean up temporary file
    os.remove(temp_file)


