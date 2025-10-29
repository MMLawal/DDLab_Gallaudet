import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Load Data ===
df = pd.read_csv("10descriptors_50_ligands.csv")

# === Define properties and their acceptable ranges (lower, upper) ===
properties = ['logP', 'logD', 'logS', 'Fsp3', 'nHet', 'nHD', 'nHA', 'nStereo', 'pka_basic', 'pka_acidic']
limits = {
    'logP': (0, 3),
    'logD': (1, 3),
    'logS': (-4, 0.5),
    'Fsp3': (0.25, 1),
    'nHet': (1, 15),
    'nHD': (0, 7),
    'nHA': (0, 12),
    'nStereo': (0, 2),
    'pka_basic': (3, 10),
    'pka_acidic': (2, 12)
}

# Note: pka_basic and pka_acidic are excluded as no limits were provided.

# === Radar chart setup ===
num_vars = len(properties)
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]  # Complete the loop

# === Plot ===
plt.figure(figsize=(10, 10))
ax = plt.subplot(111, polar=True)

# --- Plot Lower and Upper Limits ---
# Create arrays for the lower and upper limit polygons
lower_limits = [limits[prop][0] for prop in properties]
upper_limits = [limits[prop][1] for prop in properties]

# Add closing point for the polygons
lower_limits += [lower_limits[0]]
upper_limits += [upper_limits[0]]

# Plot and fill the acceptable region
ax.plot(angles, lower_limits, color='green', linewidth=2, label='Lower Limit')
ax.fill(angles, lower_limits, color='green', alpha=0.2)

ax.plot(angles, upper_limits, color='blue', linewidth=2, label='Upper Limit')
ax.fill(angles, upper_limits, color='blue', alpha=0.2)

# --- Plot Compound Properties ---
for i, ligand in enumerate(df['Compounds']):
    values = df.loc[df['Compounds'] == ligand, properties].values.flatten().tolist()
    values += values[:1]  # Close the loop
    ax.plot(angles, values, label=ligand, linewidth=1, alpha=0.7)
    # Optionally, you can add a filled area for each compound with low alpha
    # ax.fill(angles, values, alpha=0.05)

# === Customize ===
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles[:-1]), properties)

# Set radial limits to encompass all data and limits
all_values = []
for prop in properties:
    all_values.extend([limits[prop][0], limits[prop][1]])
    all_values.extend(df[prop].tolist())
ax.set_ylim(min(all_values), max(all_values))

#plt.title("Radar Chart of Molecular Descriptors with Acceptable Ranges", size=16, weight='bold')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.1), fontsize=7)
plt.tight_layout()
plt.show()
