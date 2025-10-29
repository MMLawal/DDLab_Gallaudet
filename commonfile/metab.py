import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === Load Data ===
df = pd.read_csv("metab_50_supern.csv")

# === Define Metabolic Descriptors for Y-axis ===
metabolic_params = [
    'HLM.Stability',        # HLM stability
    'CYP1A2-inh', 'CYP1A2-sub',
    'CYP2C19-inh', 'CYP2C19-sub',
    'CYP2C9-inh', 'CYP2C9-sub',
    'CYP2D6-inh', 'CYP2D6-sub',
    'CYP3A4-inh', 'CYP3A4-sub',
    'CYP2B6-inh', 'CYP2B6-sub',
    'CYP2C8-inh'
]

# Ensure consistent column names
df.columns = [c.strip() for c in df.columns]

# Select only relevant columns (exclude "Compounds" if present)
if "Compounds" in df.columns:
    df_plot = df[metabolic_params].copy()
else:
    df_plot = df.copy()

# === Convert values into binary classification ===
# Threshold = 0.5 → below = non-inhibitor/substrate (0), above = inhibitor/substrate (1)
threshold = 0.5
binary_df = df_plot.map(lambda x: 1 if x >= threshold else 0)

# === Transpose for heatmap ===
# Rows (y-axis): metabolic descriptors
# Columns (x-axis): Ligands 1–50
binary_df = binary_df.transpose()

# === Plot Heatmap ===
plt.figure(figsize=(8, 7))

sns.heatmap(
    binary_df,
    cmap=["lightgreen", "red"],
    cbar=False,
    linewidths=0.5,
    linecolor='white',
    square=True
)

# === Customize the Plot ===
#plt.title("Metabolic Descriptors for Ligands 1–50", fontsize=12, fontweight='bold', pad=15)
plt.xlabel("Ligands", fontsize=8, fontweight='bold')
plt.ylabel("Metabolic Descriptors", fontsize=8, fontweight='bold')

# Label x-axis as ligand numbers 1–50
plt.xticks(
    ticks=range(binary_df.shape[1]),
    labels=[f"L{i+1}" for i in range(binary_df.shape[1])],
    rotation=90,
    fontsize=8
)

# Label y-axis with the descriptor names (already from DataFrame index)
plt.yticks(fontsize=8)

# === Add Custom Legend ===
legend_elements = [
    Patch(facecolor='lightgreen', edgecolor='black', label='Non-Inhibitor/Non-Substrate'),
    Patch(facecolor='red', edgecolor='black', label='Inhibitor/Substrate')
]
plt.legend(
    handles=legend_elements,
    loc='lower center',
    bbox_to_anchor=(0.5, 1.0),
    ncol=2,
    frameon=False,
    fontsize=8
)

plt.tight_layout()
plt.show()

