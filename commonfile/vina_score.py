import re
import glob

def extract_affinity_mode1(filename):
    """
    Extracts the affinity value for mode 1 from the specified Vina output file.

    Args:
        filename: The path to the Vina output file.

    Returns:
        The affinity value for mode 1 as a float, or None if not found.
    """
    with open(filename, 'r') as f:
        for line in f:
            match = re.search(r"1\s+(-?\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+", line)
            if match:
                return float(match.group(1))
    return None

# Get all .txt files in the current directory
txt_files = glob.glob("*.txt")

# Open the output file
with open("result.txt", "w") as output_file:
    for file in txt_files:
        affinity_mode1 = extract_affinity_mode1(file)
        if affinity_mode1 is not None:
            output_file.write(f"{file}: {affinity_mode1}\n")
        else:
            output_file.write(f"{file}: NA\n")

