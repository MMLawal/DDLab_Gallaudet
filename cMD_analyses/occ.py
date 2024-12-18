import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Load the topology and trajectory files
u = mda.Universe('com.top', 'bs3_c-9.dcd')

# Define hydrophobic residues (common hydrophobic residues)
hydrophobic_residues = ['ARG', 'ASP', 'VAL', 'LEU', 'ILE', 'LYS', 'HIS', 'GLU', 'TYR', 'PHE', 'MET', 'PRO', 'ALA', 'CYS', 'TRP']

# Select hydrophobic atoms in the protein and ligand
protein_hydrophobic = u.select_atoms(f"protein and (resname {' '.join(hydrophobic_residues)}) and name C* and not backbone")
ligand_hydrophobic = u.select_atoms("resname LIG")  # Adjust 'LIG' to match your ligand name

# Define cutoff distance for hydrophobic interactions
hydrophobic_cutoff = 3  # Hydrophobic interaction cutoff (in Angstroms)

# Initialize a dictionary to track interactions for each residue
hydrophobic_counts = defaultdict(int)

# Get total number of frames for percent occupancy calculation
total_frames = len(u.trajectory)

# Iterate over each frame of the trajectory
for ts in u.trajectory:
    # Hydrophobic interaction: Calculate distances between protein and ligand hydrophobic atoms
    hydrophobic_distances = distance_array(protein_hydrophobic.positions, ligand_hydrophobic.positions)
    
    # Find pairs of interacting hydrophobic atoms below the cutoff distance
    hydrophobic_interacting_pairs = np.where(hydrophobic_distances < hydrophobic_cutoff)

    # For each interacting hydrophobic atom, increase the count for the residue it belongs to
    for i in hydrophobic_interacting_pairs[0]:
        # Add 370 to the residue index and format it as 'ResidueName + ResidueIndex'
        residue_name_index = f"{protein_hydrophobic[i].resname}{protein_hydrophobic[i].resid + 370}"
        hydrophobic_counts[residue_name_index] += 1

# Calculate percent occupancy for each residue
hydrophobic_occupancy = {resid: (count / total_frames) * 100 for resid, count in hydrophobic_counts.items()}

# Filter to keep only residues with percent occupancy >= 10%
filtered_hydrophobic_occupancy = {resid: occupancy for resid, occupancy in hydrophobic_occupancy.items() if occupancy >= 10}

# Open a file to write the filtered percent occupancy results
with open('interaction_percent_occupancy_above_10.txt', 'w') as output_file:
    output_file.write('Residue,Hydrophobic Occupancy (%)\n')
    for resid, occupancy in filtered_hydrophobic_occupancy.items():
        output_file.write(f'{resid},{occupancy:.2f}\n')

# Prepare data for plotting
residues = list(filtered_hydrophobic_occupancy.keys())
hydrophobic_percents = list(filtered_hydrophobic_occupancy.values())

# Sort residues by hydrophobic occupancy for better visualization
sorted_residues, sorted_hydro = zip(*sorted(zip(residues, hydrophobic_percents), key=lambda x: x[1], reverse=True))

# Plot the percent occupancy of residues with occupancy >= 10%
plt.figure(figsize=(10, 6))
plt.bar(sorted_residues, sorted_hydro, color='blue')

plt.xlabel('Residue')
plt.ylabel('Percent Occupancy (%)')
plt.title('Percent Occupancy of Hydrophobic Contacts Per Residue (>= 10%)')
plt.xticks(rotation=90)
plt.tight_layout()

# Optionally save the plot
# plt.savefig('percent_occupancy_above_10_plot.png', dpi=300)

# Display the plot
plt.show()

