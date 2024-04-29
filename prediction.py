import pandas as pd
import os
from glob import glob

# Get list of all CSV files in the ../output_csvs directory
csv_files = glob("../output_csvs/*.csv")

# Populate similar_proteins with data from CSV files
similar_proteins = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file, header=None)
    # Extract the amino acid sequences from the second column and store them as lists
    protein_sequence = df.iloc[:, 1].tolist()
    similar_proteins.append(protein_sequence)

print("Number of protein sequences loaded:", len(similar_proteins))

stored_sasa_val=[]
for csv_file in csv_files:
    df = pd.read_csv(csv_file, header=None)
    # Extract the amino acid sequences from the third column and store them as lists
    sasa_val = df.iloc[:, 2].tolist()
    stored_sasa_val.append(sasa_val)
# print(similar_proteins[1])
# Mapping between one-letter and three-letter amino acid representations
one_letter_to_three_letter = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

three_letter_to_one_letter = {v: k for k, v in one_letter_to_three_letter.items()}

# Load query protein sequence from FASTA file
with open("rcsb_pdb_1A3N.fasta", "r") as f:
    lines = f.readlines()
    query_sequence = lines[1].strip()  # Get sequence from the second line

# Define function to calculate similarity score for residue types
def calculate_similarity(residue1, residue2):
    # Define similarity scores based on domain knowledge
    similarity_scores = {
        "hydrophobic": {"A", "I", "L", "M", "F", "W", "V"},
        "hydrophilic": {"R", "N", "D", "E", "Q", "H", "K", "S", "T", "Y"},
        # Add more similarity scores as needed
    }
    # Convert three-letter residue to one-letter residue
    residue2 = three_letter_to_one_letter.get(residue2, residue2)
    if residue1 == residue2:
        return 1https://freesasa.github.io/
    elif residue1 in similarity_scores["hydrophobic"] and residue2 in similarity_scores["hydrophobic"]:
        return 0.8
    elif residue1 in similarity_scores["hydrophilic"] and residue2 in similarity_scores["hydrophilic"]:
        return 0.8
    elif residue1 in similarity_scores["hydrophilic"] and residue2 in similarity_scores["hydrophobic"]:
        return 2
    elif residue1 in similarity_scores["hydrophobic"] and residue2 in similarity_scores["hydrophilic"]:
        return 0.2
    else:
        return 0.5

def local_alignment(seq1, seq2, gap_open_penalty):
    """
    Perform local sequence alignment between two sequences.

    Arguments:
    seq1 (str): The query sequence.
    seq2 (str): The protein sequence.
    gap_open_penalty (int): Penalty for opening a gap.

    Returns:
    aligned_seq1 (str): Aligned query sequence.
    aligned_seq2 (str): Aligned protein sequence.
    start_index (int): Index in seq2 where the alignment starts.
    """

    # Initialize variables
    m = len(seq1)
    n = len(seq2)
    max_score = 0
    max_i = 0
    max_j = 0

    # Initialize score matrix and traceback matrix
    score = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[0] * (n + 1) for _ in range(m + 1)]

    # Fill score and traceback matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + (1 if seq1[i - 1] == seq2[j - 1] else -1)
            delete = score[i - 1][j] - gap_open_penalty
            insert = score[i][j - 1] - gap_open_penalty
            score[i][j] = max(0, match, delete, insert)
            if score[i][j] > max_score:
                max_score = score[i][j]
                max_i = i
                max_j = j
            if score[i][j] == match:
                traceback[i][j] = 1  # Match or mismatch
            elif score[i][j] == delete:
                traceback[i][j] = 2  # Gap in seq1
            elif score[i][j] == insert:
                traceback[i][j] = 3  # Gap in seq2

    # Traceback to find alignment start point
    i = max_i
    j = max_j
    while i > 0 and j > 0 and score[i][j] > 0:
        if traceback[i][j] == 1:  # Match or mismatch
            break
        elif traceback[i][j] == 2:  # Gap in seq1
            i -= 1
        elif traceback[i][j] == 3:  # Gap in seq2
            j -= 1

    # Perform alignment from the alignment start point
    aligned_seq1 = seq1[:i]
    aligned_seq2 = seq2[j:max_j]

    # Convert aligned_seq2 to string
    aligned_seq2 = ''.join(aligned_seq2)

    # Extend aligned_seq2 with gaps if necessary
    if i == 0:  # Alignment starts at the beginning of seq1
        aligned_seq2 = "-" * (max_j - len(aligned_seq2)) + aligned_seq2
        start_index = 0
    else:  # Alignment starts within seq1
        start_index = j

    return aligned_seq1, aligned_seq2, start_index


# Perform local alignment with each similar protein sequence
predicted_sasa_values = {}
i=0
for protein_sequence in similar_proteins:
    print("Query sequence:", query_sequence)
    print("Protein sequence:", protein_sequence)
    aligned_query, aligned_protein, start_index = local_alignment(query_sequence, protein_sequence, gap_open_penalty=5)  # Adjust the gap_open_penalty as needed
    print("Aligned query sequence:", aligned_query)
    print("Aligned protein sequence:", aligned_protein)
    print("Start index of alignment in protein sequence:", start_index)
    aligned_query_pos = 0
    for i in range(len(aligned_query)):
        query_residue = aligned_query[i]
        protein_residue = aligned_protein[i]
        if query_residue != "-":  # Ignore gaps
            
            aligned_query_pos += 1
            if protein_residue != "-":  # Only update position in similar sequence if it's not a gap
                similarity_score = calculate_similarity(query_residue, protein_residue)
                residue_position = aligned_query_pos
                residue_position1 = start_index + aligned_query_pos - 1  # Adjust the residue position using the start index
                residue_key = (residue_position, query_residue)
                if residue_key not in predicted_sasa_values:
                    predicted_sasa_values[residue_key] = 0
                predicted_sasa_values[residue_key] += stored_sasa_val[i][residue_position1] * similarity_score / len(similar_proteins)  # Assuming the SASA value is in the third column
    i=i+1
    
print("Number of predicted SASA values:", len(predicted_sasa_values))
# Output predicted SASA values for the query sequence
for residue_key, sasa in predicted_sasa_values.items():
    residue_position, residue_type = residue_key
    print(f"Residue {residue_type} at position {residue_position}: Predicted SASA = {sasa}")
