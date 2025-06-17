import os
import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
import numpy as np

# Directories
fingerprint_and_activity_dir = "../1-data/2-fingerprints+activity/"
activity_similarity_dir = "../1-data/3-similarity+activity/"

def list_to_bitvector(bitlist):
    """Convert a list of 0s and 1s to an RDKit ExplicitBitVect."""
    bitvect = ExplicitBitVect(len(bitlist))
    for i, bit in enumerate(bitlist):
        if bit == 1:
            bitvect.SetBit(i)
    return bitvect

def calculate_similarity(fp1, fp2, measure):
    """Calculates various similarities based on the selected measure."""
    if measure == "tanimoto":
        return f"{DataStructs.TanimotoSimilarity(fp1, fp2):.6f}"
    elif measure == "dice":
        return f"{DataStructs.DiceSimilarity(fp1, fp2):.6f}"
    elif measure == "cosine":
        return f"{DataStructs.CosineSimilarity(fp1, fp2):.6f}"
    elif measure == "sokal":
        return f"{DataStructs.SokalSimilarity(fp1, fp2):.6f}"
    elif measure == "tversky_lig_hit":
        return f"{DataStructs.TverskySimilarity(fp1, fp2, 0.7, 0.3):.6f}"  # Ligand-focused  # α=0.7, β=0.3 - emphasizes query features ==== Tversky (Ligand-weighted)
    elif measure == "tversky_hit_lig":
        return f"{DataStructs.TverskySimilarity(fp1, fp2, 0.3, 0.7):.6f}"  # Hit-focused # α=0.3, β=0.7 - emphasizes reference/hit features ===== Tversky (Reference-weighted)
    elif measure == "kulczynski":
        return f"{DataStructs.KulczynskiSimilarity(fp1, fp2):.6f}"
    elif measure == "russell_rao":
        return f"{DataStructs.RusselSimilarity(fp1, fp2):.6f}"
    elif measure == "mcconnaughey":
        return f"{DataStructs.McConnaugheySimilarity(fp1, fp2):.6f}"
    elif measure == "braun_blanquet":
        return f"{DataStructs.BraunBlanquetSimilarity(fp1, fp2):.6f}"
    elif measure == "euclidean":
        return f"{calculate_euclidean_similarity(fp1, fp2):.6f}"
    elif measure == "manhattan":
        return f"{calculate_manhattan_similarity(fp1, fp2):.6f}"
    else:
        raise ValueError("Unsupported similarity measure.")


def calculate_euclidean_similarity(fp1, fp2):
    """Calculates the Euclidean similarity between two bit vectors."""
    diff = np.array(fp1) - np.array(fp2)
    distance = np.sqrt(np.sum(diff ** 2))
    max_distance = np.sqrt(len(fp1))  # Maximum possible Euclidean distance
    similarity = 1 - (distance / max_distance)
    return similarity


def calculate_manhattan_similarity(fp1, fp2):
    """Calculates the Manhattan similarity between two bit vectors."""
    distance = np.sum(np.abs(np.array(fp1) - np.array(fp2)))
    max_distance = len(fp1)  # Maximum possible Manhattan distance
    similarity = 1 - (distance / max_distance)
    return similarity


def process_csv_files():
    # Similarity measures to calculate
    measures = ["tanimoto", "dice", "cosine", "sokal", "tversky_lig_hit", "tversky_hit_lig", 
                "kulczynski", "russell_rao", "mcconnaughey", "braun_blanquet"]

    # "euclidean", "manhattan", 

    # Loop through all subdirectories and CSV.bz2 files in fingerprint_and_activity_dir
    for target_dir in os.listdir(fingerprint_and_activity_dir):
        target_path = os.path.join(fingerprint_and_activity_dir, target_dir)

        if os.path.isdir(target_path):
            for fingerprint_file in os.listdir(target_path):
                if fingerprint_file.endswith('.csv.bz2'):
                    # Read the compressed CSV file
                    file_path = os.path.join(target_path, fingerprint_file)
                    df = pd.read_csv(file_path, compression='bz2')

                    # Select hits where hit == 1
                    hits_df = df[df['hit'] == 1]
                    print(f"Processing {fingerprint_file} - Number of hits: {len(hits_df)}")

                    # Convert the dataframe rows to RDKit bit vectors
                    hits_fingerprints = [list_to_bitvector(fplist) for fplist in hits_df.drop(columns=['hit']).values.tolist()]
                    all_fingerprints = [list_to_bitvector(fplist) for fplist in df.drop(columns=['hit']).values.tolist()]

                    # Process for each similarity measure
                    for measure in measures:
                        print(f"Calculating {measure} similarity")

                        # Create list to store similarity results
                        similarity_results = []

                        # For each row, calculate the similarity to all hits
                        for index, fingerprint in enumerate(all_fingerprints):
                            row_similarities = []

                            for hit_fp in hits_fingerprints:
                                similarity = calculate_similarity(fingerprint, hit_fp, measure)
                                row_similarities.append(similarity)

                            # Store the row with the corresponding hit value and similarities
                            similarity_results.append([df.loc[index, 'hit']] + row_similarities)

                        # Convert results to DataFrame and save
                        similarity_df = pd.DataFrame(similarity_results, columns=['hit'] + [f"similarity{i+1}" for i in range(len(hits_fingerprints))])

                        # Define output directory and file path
                        output_dir = os.path.join(activity_similarity_dir, measure, target_dir)
                        os.makedirs(output_dir, exist_ok=True)
                        output_file = os.path.join(output_dir, fingerprint_file)

                        # Save the result to the output file in compressed format
                        print(f"Saving {measure} similarity data to: {output_file}")
                        similarity_df.to_csv(output_file, index=False, compression='bz2')


if __name__ == "__main__":
    process_csv_files()
