from Bio import SeqIO

def extract_proteins_with_eight_cysteines(input_fasta_file, output_fasta_file):
    """
    Extracts protein sequences from a FASTA file that contain exactly 8 cysteines.

    Args:
        input_fasta_file (str): Path to the input FASTA file.
        output_fasta_file (str): Path to the output FASTA file for filtered sequences.
    """
    eight_cysteine_proteins = []

    try:
        for record in SeqIO.parse(input_fasta_file, "fasta"):
            cysteine_count = record.seq.count('C')
            if cysteine_count == 8:
                eight_cysteine_proteins.append(record)
    except FileNotFoundError:
        print(f"Error: Input file '{input_fasta_file}' not found.")
        return

    if eight_cysteine_proteins:
        with open(output_fasta_file, "w") as output_handle:
            SeqIO.write(eight_cysteine_proteins, output_handle, "fasta")
        print(f"Successfully extracted {len(eight_cysteine_proteins)} proteins with 8 cysteines.")
        print(f"Filtered sequences saved to '{output_fasta_file}'.")
    else:
        print("No proteins with exactly 8 cysteines were found in the input file.")

if __name__ == "__main__":
    # --- User Configuration ---
    # Replace 'your_input_file.fasta' with the actual name of your FASTA file
    input_file = "your_input_file.fasta"
    # Replace 'eight_cysteine_proteins.fasta' with your desired output file name
    output_file = "eight_cysteine_proteins.fasta"
    # ------------------------

    extract_proteins_with_eight_cysteines(input_file, output_file)