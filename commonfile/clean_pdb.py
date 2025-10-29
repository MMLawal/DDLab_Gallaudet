import os

def process_pdb_file(input_file, output_file):
    atom_counter = 1
    in_ligand = False
    ligand_start_atom = 3398

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if any(line.startswith(prefix) for prefix in ['REMARK', 'COMPND', 'AUTHOR', 'CONECT', 'MASTER']):
                continue  # skip unwanted lines

            if line.startswith('TER'):
                outfile.write(line)
                in_ligand = True
                atom_counter = ligand_start_atom  # Reset counter for ligand
                continue

            if line.startswith(('ATOM', 'HETATM')):
                # Update atom number
                line = line[:6] + f"{atom_counter:5d}" + line[11:]

                # Update HETATM to ATOM
                if line.startswith('HETATM'):
                    line = 'ATOM  ' + line[6:]

                # Replace UNL or UNK residue identifiers
                if 'UNL     1' in line:
                    line = line.replace('UNL     1', 'LIG   214')
                if 'UNK     0' in line:
                    line = line.replace('UNK     0', 'LIG   214')

                atom_counter += 1

            outfile.write(line)

# Process files: 4hco_l1_m1.pdb to 4hco_l46_m1.pdb
for i in range(1, 47):
    pdb_file = f"4xb0_bs1_tq.pdb"
    output_file = f"plk2_bs1_tq.pdb"

    if os.path.exists(pdb_file):
        print(f"Processing {pdb_file}...")
        process_pdb_file(pdb_file, output_file)
    else:
        print(f"File not found: {pdb_file}")

