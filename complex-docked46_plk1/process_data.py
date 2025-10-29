import re

def extract_affinity_mode1(filename):
  """Extracts the affinity value for mode 1 from the specified Vina output file.

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

# Process multiple Vina output files and save affinities to a file
with open("result.txt", "w") as output_file:
  for i in range(1, 47):
    vina_result_file = f"c{i}.txt"
    affinity_mode1 = extract_affinity_mode1(vina_result_file)
    if affinity_mode1 is not None:
      output_file.write(f"{affinity_mode1}\n")
    else:
      output_file.write("NA\n")

