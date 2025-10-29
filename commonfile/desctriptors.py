import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your descriptor data
df = pd.read_csv("10descriptors_50_ligands.csv")

# Define the properties of interest
properties = [
    "logP", "logD", "logS", "Fsp3", "nHet", 
    "nHD", "nHA", "pka_basic", "pka_acidic", "nStereo"
]

# Define lower and upper limits (l, u)
limits = {
    "logP": (0, 3),
    "logD": (1, 3),
    "logS": (-4, 0.5),
    "Fsp3": (0.25, 1),
    "nHet": (1, 15),
    "nHD": (0, 7),
    "nHA": (0, 12),
    "pka_basic": (3, 10),
    "pka_acidic": (3, 10),
    "nStereo": (0, 2)
}

# Choose one or more ligands to visualize
selected_ligands = ["Ligand 1", "Ligand 2", "Ligand 3", "Ligand 4", "Ligand 5", "Ligand 6", "Ligand 7", "Ligand 8", "Ligand 9", "Ligand 10", "Ligand 11", "Ligand 12", "Ligand 13", "Ligand 14", "Ligand 15", "Ligand 16", "Ligand 17", "Ligand 18", "Ligand 19", "Ligand 20", "Ligand 21", "Ligand 22", "Ligand 23", "Ligand 24", "Ligand 25", "Ligand 26", "Ligand 27", "Ligand 28", "Ligand 29", "Ligand 30", "Ligand 31", "Ligand 32", "Ligand 33", "Ligand 34", "Ligand 35", "Ligand 36", "Ligand 37", "Ligand 38", "Ligand 39", "Ligand 40", "Ligand 41", "Ligand 42", "Ligand 43", "Ligand 44", "Ligand 45", "Ligand 46", "Ligand 47", "Ligand 48", "Ligand 49", "Ligand 50"]  # you can modify this list

# Compute number of axes
num_vars = len(properties)

# Compute angle for each axis
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]  # close the radar circle

# Function to normalize data (scale 0-1 between lower and upper bounds)
def normalize(value, lower, upper):
    return (value - lower) / (upper - lower)

# Create Radar Chart
plt.figure(figsize=(8, 8))
ax = plt.subplot(111, polar=True)

# Plot acceptable region polygon
lower_bounds = [0 for _ in properties]
upper_bounds = [1 for _ in properties]
ax.fill(angles, upper_bounds + [upper_bounds[0]], color='lightgreen', alpha=0.2, label='Acceptable Region')

# Plot each ligand’s normalized descriptor profile
for ligand in selected_ligands:
    row = df[df["Compounds"] == ligand]
    if not row.empty:
        values = []
        for prop in properties:
            val = row[prop].values[0]
            l, u = limits[prop]
            norm = normalize(val, l, u)
            norm = np.clip(norm, 0, 1)  # clamp outliers
            values.append(norm)
        values += values[:1]
        ax.plot(angles, values, linewidth=1.5, label=ligand)
        ax.fill(angles, values, alpha=0.1)

# Configure radar labels and appearance
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
plt.xticks(angles[:-1], properties, fontsize=9)
ax.set_rlabel_position(0)
plt.yticks([0, 0.5, 1], ["Low", "Mid", "High"], color="grey", size=8)
plt.ylim(0, 1)

plt.title("Molecular Descriptor Radar Chart", size=14, pad=20)
plt.legend(bbox_to_anchor=(1.3, 1.1))
plt.tight_layout()
plt.show()

