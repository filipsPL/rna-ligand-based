import csv
from openbabel import pybel

def calculate_fingerprint(smiles_file, fingerprint_type, output_csv):
    # Define the fingerprint sizes based on the type
    fingerprint_sizes = {
        "fp2": 1024,
        "fp3": 55,
        "fp4": 307
    }

    # Read the SMILES strings and IDs from the input file
    smiles_list = []
    ids = []
    with open(smiles_file, 'r') as f:
        for line in f:
            if line.strip():
                # Split the line into SMILES and ID by whitespace
                parts = line.strip().split()
                if len(parts) == 2:
                    smiles, mol_id = parts
                    smiles_list.append(smiles)
                    ids.append(mol_id)

    # Initialize molecule objects and calculate fingerprints
    mols = [pybel.readstring("smi", smi) for smi in smiles_list]
    fps = [mol.calcfp(fptype=fingerprint_type) for mol in mols]

    # Get the correct size for the selected fingerprint type
    num_bits = fingerprint_sizes[fingerprint_type]

    # Write the fingerprints to the output CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write the header (ID, bit0, bit1, ...)
        header = ["ID"] + [f"bit{i}" for i in range(num_bits)]
        csvwriter.writerow(header)

        # Write each molecule's ID and full fingerprint
        for i, fp in enumerate(fps):
            # Initialize the full bit vector with 0s
            bits = [0] * num_bits
            # Set the bits that are active to 1
            for bit_index in fp.bits:
                bits[bit_index - 1] = 1  # OpenBabel indices start from 1, so we subtract 1 for Python lists
            csvwriter.writerow([ids[i]] + bits)


if __name__ == "__main__":
    # Specify the input SMILES file, fingerprint type, and output CSV file
    smiles_file = "../0-input/robin-normalized_noempty.smi"

    for fingerprint_type in ['fp2', 'fp3', 'fp4']:

        print(fingerprint_type)

        output_csv = f"../1-data/1-fingerprints/openbabel_{fingerprint_type}.csv"

        # Calculate fingerprints and save to the CSV file
        calculate_fingerprint(smiles_file, fingerprint_type, output_csv)
