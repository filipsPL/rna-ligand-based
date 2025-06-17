import csv
import random

def generate_random_bits(bit_count, density_percentage):
    """Generates a random bitset of length bit_count with a given density of active bits."""
    num_active_bits = int(bit_count * density_percentage / 100)
    bits = [0] * bit_count
    active_bit_indices = random.sample(range(bit_count), num_active_bits)
    for idx in active_bit_indices:
        bits[idx] = 1
    return bits

def process_smiles_file(smiles_file):
    """Reads SMILES strings and IDs from the input file."""
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
    return smiles_list, ids

def write_random_fingerprints(smiles_file):
    bit_lengths = [64, 128, 256, 512, 1024]
    densities = [1, 5, 10, 40]

    # Read SMILES strings and IDs
    smiles_list, ids = process_smiles_file(smiles_file)

    for bit_count in bit_lengths:
        for density in densities:
            output_csv = f"../1-data/1-fingerprints/random_{bit_count}_{density}.csv"
            with open(output_csv, 'w', newline='') as csvfile:
                print(f"{bit_count} - {density}")
                csvwriter = csv.writer(csvfile)
                # Write header
                header = ["ID"] + [f"bit{i}" for i in range(num_bits)]
                csvwriter.writerow(header)

                # Write random fingerprints for each molecule
                for i, mol_id in enumerate(ids):
                    random_bits = generate_random_bits(bit_count, density)
                    csvwriter.writerow([mol_id] + random_bits)

if __name__ == "__main__":
    # Specify the input SMILES file
    smiles_file = "../0-input/robin-normalized_noempty.smi"

    # Generate and save random fingerprints
    write_random_fingerprints(smiles_file)
