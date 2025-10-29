import os

def insert_chain_id_to_pdb(input_path, output_path, chain_id='A'):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM')):
                line = line[:21] + chain_id + line[22:]
            outfile.write(line)

# Loop over directories and files
for i in range(1, 47):  # Adjust the range as needed
    dir_name = "." #f"lig{i}"
    input_file = os.path.join(f"4hco_bs2_allop-a.pdb")
    output_file = os.path.join(f"plk1_bs2_allop-a.pdb")

    if os.path.exists(input_file):
        print(f"Processing {input_file}...")
        insert_chain_id_to_pdb(input_file, output_file, chain_id='A')
    else:
        print(f"File not found: {input_file}")


