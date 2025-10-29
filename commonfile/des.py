import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Load Data ===
df = pd.read_csv("10descriptors_50_ligands.csv")

# === Select properties for the radar chart ===
properties = ['logP', 'logD', 'logS', 'Fsp3', 'nHet', 'nHD', 'nHA', 'pka_basic', 'pka_acidic', 'nStereo']

# Create transformed columns for negative logs
df['logP'] = df['logP']
df['logD'] = df['logD']
df['logS'] = df['logS']

# === Choose ligands to visualize ===
# You can modify this list, e.g., ['Ligand 1', 'Ligand 2', 'Ligand 3']
ligands_to_plot = ["Ligand 1", "Ligand 2", "Ligand 3", "Ligand 4", "Ligand 5", "Ligand 6", "Ligand 7", "Ligand 8", "Ligand 9", "Ligand 10", "Ligand 11", "Ligand 12", "Ligand 13", "Ligand 14", "Ligand 15", "Ligand 16", "Ligand 17", "Ligand 18", "Ligand 19", "Ligand 20", "Ligand 21", "Ligand 22", "Ligand 23", "Ligand 24", "Ligand 25", "Ligand 26", "Ligand 27", "Ligand 28", "Ligand 29", "Ligand 30", "Ligand 31", "Ligand 32", "Ligand 33", "Ligand 34", "Ligand 35", "Ligand 36", "Ligand 37", "Ligand 38", "Ligand 39", "Ligand 40", "Ligand 41", "Ligand 42", "Ligand 43", "Ligand 44", "Ligand 45", "Ligand 46", "Ligand 47", "Ligand 48", "Ligand 49", "Ligand 50"]

# === Radar chart setup ===
num_vars = len(properties)
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]  # complete the loop

# === Plot ===
plt.figure(figsize=(8, 8))
ax = plt.subplot(111, polar=True)

for ligand in ligands_to_plot:
    values = df.loc[df['Compounds'] == ligand, properties].values.flatten().tolist()
    values += values[:1]  # close the loop
    ax.plot(angles, values, label=ligand, linewidth=2)
    ax.fill(angles, values, alpha=0.1)

# === Customize ===
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles[:-1]), properties)
plt.title("Radar Chart of Molecular Descriptors", size=14, weight='bold')
plt.legend(loc='upper left', bbox_to_anchor=(1.2, 1.1), fontsize=7)
plt.tight_layout()
plt.show()

