import re
import os

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary of protein sequences.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where keys are protein IDs and values are sequences.
              Returns an empty dictionary if the file cannot be read or is empty.
    """
    proteins = {}
    current_protein_id = None
    current_sequence = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines

                if line.startswith('>'):
                    # New protein header
                    if current_protein_id and current_sequence:
                        proteins[current_protein_id] = "".join(current_sequence)
                    current_protein_id = line[1:].split(' ')[0] # Get ID, ignore description
                    current_sequence = []
                else:
                    # Sequence line
                    current_sequence.append(line)

            # Add the last protein after the loop finishes
            if current_protein_id and current_sequence:
                proteins[current_protein_id] = "".join(current_sequence)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return {}
    except Exception as e:
        print(f"An error occurred while reading the FASTA file: {e}")
        return {}

    return proteins

def detect_two_or_more_8cm_motifs(proteins_data):
    """
    Detects proteins that contain two or more 8CM motifs.

    The 8CM motif pattern used is C-Xn-C-Xn-CC-Xn-CXC-Xn-C-Xn-C,
    where 'n' ranges from 2 to 70 for all Xn regions, including the X in CXC.

    Args:
        proteins_data (dict): A dictionary of protein IDs and sequences.

    Returns:
        dict: A dictionary where keys are protein IDs and values are the count
              of 8CM motifs found for proteins with 2 or more motifs.
    """
    # Define the regular expression for the 8CM motif.
    # We use non-capturing groups (?:...) where we don't need to extract
    # the content of the Xn regions, only count the matches.
    # The pattern is (C)(.{2,70})(C)(.{2,70})(C)(C)(.{2,70})(C)(.{0,70})(C)(.{2,70})(C)(.{2,70})(C)
    # The problem asks for "two 8CM motifs", implying non-overlapping.
    # re.findall by default finds non-overlapping matches.
    eight_cysteine_motif_regex = re.compile(
        r'C.{2,70}C.{2,70}CC.{2,70}C.{0,70}C.{2,70}C.{2,70}C'
    )

    proteins_with_multiple_motifs = {}
    print(f"Searching for proteins with two or more 8CM motifs using pattern: {eight_cysteine_motif_regex.pattern}")

    for protein_id, sequence in proteins_data.items():
        # Find all non-overlapping occurrences of the motif
        matches = eight_cysteine_motif_regex.findall(sequence)
        motif_count = len(matches)

        if motif_count >= 2:
            proteins_with_multiple_motifs[protein_id] = motif_count
            print(f"  Found {motif_count} 8CM motifs in protein: {protein_id}")

    return proteins_with_multiple_motifs

def main():
    """
    Main function to run the protein detection script.
    Prompts the user for a FASTA file path and displays results.
    """
    print("--- Protein Detector for Two or More 8CM Motifs ---")
    fasta_file = input("Please enter the path to your FASTA file: ").strip()

    if not os.path.exists(fasta_file):
        print(f"Error: The specified file '{fasta_file}' does not exist.")
        return

    print(f"\nParsing FASTA file: {fasta_file}...")
    proteins = parse_fasta(fasta_file)

    if not proteins:
        print("No proteins found or file could not be parsed. Exiting.")
        return

    print(f"Successfully parsed {len(proteins)} proteins.")

    print("\nDetecting proteins with two or more 8CM motifs...")
    multi_motif_proteins = detect_two_or_more_8cm_motifs(proteins)

    if multi_motif_proteins:
        print(f"\n--- Found {len(multi_motif_proteins)} proteins with two or more 8CM motifs ---")
        for protein_id, count in multi_motif_proteins.items():
            print(f"Protein ID: {protein_id}, Number of 8CM Motifs: {count}")
    else:
        print("\nNo proteins with two or more 8CM motifs were found.")

if __name__ == "__main__":
    main()
