import pandas as pd
from map4 import MAP4Calculator
from rdkit import Chem

def build_dataframe(input_file_name):
    # Read SMILES and Name from input file
    df = pd.read_csv(input_file_name, header=None, sep="\t")
    df.columns = ["SMILES", "Name"]

    # Initialize the MAP4Calculator
    m4_calc = MAP4Calculator(is_folded=True, dimensions=fingerprint_size)

    # Generate RDKit Mol objects and MAP4 bit vectors
    df['mol'] = [Chem.MolFromSmiles(x) for x in df.SMILES]
    df['map4'] = [m4_calc.calculate(x) for x in df.mol]

    return df

def save_bitvector_to_csv(df, output_csv):
    # Prepare a bit vector for each molecule
    bitvector_df = pd.DataFrame(df['map4'].apply(lambda x: list(x)).tolist())

    # Add molecule names to the dataframe
    bitvector_df.insert(0, "Name", df["Name"])

    header = ["ID"] + [f"bit{i}" for i in range(fingerprint_size)]

    # Save the bit vectors to CSV
    bitvector_df.to_csv(output_csv, index=False, header=header)

if __name__ == "__main__":
    # Input SMILES file
    smiles_file = "../0-input/robin-normalized_noempty.smi"
    fingerprint_size = 256  # You can change this size if needed (e.g., 2048)
    output_csv = f"../1-data/1-fingerprints/map4_fingerprint_{fingerprint_size}.csv"

    # Build dataframe and calculate MAP4 bit vectors
    df = build_dataframe(smiles_file)

    # Save the bitvector to CSV
    save_bitvector_to_csv(df, output_csv)

