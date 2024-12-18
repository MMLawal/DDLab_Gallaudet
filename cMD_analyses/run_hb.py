import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from collections import defaultdict

# Load the universe with topology and trajectory files
u = mda.Universe('com.top', 'bs3_c-9.dcd')

# Initialize HydrogenBondAnalysis to analyze hydrogen bonds between ligand and protein
hbonds = HBA(universe=u, between=['resname LIG', 'protein'], d_a_cutoff=4, d_h_a_angle_cutoff=120)

# Guess hydrogen bond donors and acceptors for protein and ligand
protein_hydrogens_sel = hbonds.guess_hydrogens("protein")
protein_acceptors_sel = hbonds.guess_acceptors("protein")

ligand_hydrogens_sel = "resname LIG and name H*"
ligand_acceptors_sel = "resname LIG and name O* N* S* F* Cl*"

# Define hydrogen donors and acceptors
hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({ligand_hydrogens_sel})"
hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({ligand_acceptors_sel})"

# Run the hydrogen bond analysis
hbonds.run()

# Dictionary to store hydrogen bond occurrences per residue (by residue ID)
residue_hbond_counts = defaultdict(int)

# Get the total number of frames in the trajectory
total_frames = len(u.trajectory)

# Process hydrogen bond results and count occurrences per residue
for hbond in hbonds.hbonds:
    frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hbond
    
    # Convert atom indices to integers
    donor_idx = int(donor_idx)
    acceptor_idx = int(acceptor_idx)
    
    # Get the residues of the donor and acceptor atoms
    donor_residue = u.atoms[donor_idx].residue
    acceptor_residue = u.atoms[acceptor_idx].residue
    
    # Increment the count for each residue involved in a hydrogen bond
    residue_hbond_counts[donor_residue.resid] += 1
    residue_hbond_counts[acceptor_residue.resid] += 1

# Calculate percent occupancy per residue
residue_occupancy = {resid: (count / total_frames) * 100 for resid, count in residue_hbond_counts.items()}

# Filter and display the hydrogen bond occupancy per residue (only 10% and above)
print("Hydrogen Bond Percent Occupancy per Residue (≥ 10%):")
for resid, occupancy in residue_occupancy.items():
    if occupancy >= 10:  # Only show occupancies >= 10%
        print(f"Residue {resid}: {occupancy:.2f}%")

# Save the percent occupancy per residue to a CSV file (only 10% and above)
with open('residue_hbond_percent_occupancy.csv', 'w') as f:
    f.write('Residue,Percent Occupancy (%)\n')
    for resid, occupancy in residue_occupancy.items():
        if occupancy >= 10:  # Only write occupancies >= 10%
            f.write(f"{resid},{occupancy:.2f}\n")

print("Results saved to 'residue_hbond_percent_occupancy.csv' with only occupancies ≥ 10%.")

