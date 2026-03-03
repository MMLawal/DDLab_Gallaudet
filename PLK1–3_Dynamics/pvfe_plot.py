#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PVFE Visualization - Explicit Residue-Based Analysis
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

POCKET_RESIDUES = {
    'PLK1': {
        'BS1': [41, 42, 43, 44, 45, 46, 120, 121, 122, 146, 148, 163, 164, 165, 168, 170, 187],
        'BS2': [47, 51, 100, 101, 102, 103, 104, 105, 107, 108, 111, 115, 119],
        'BS3': [36, 37, 40, 58, 122, 123, 124, 125, 126, 127, 128, 129, 181, 183, 189],
        'BS4': [3, 35, 36, 37, 38, 130, 134, 135, 136, 137, 139, 140, 142, 156, 157, 158, 159, 174, 176, 177],
        'BS5': [20, 24, 25, 27, 28, 31, 182, 183, 184, 185, 188, 189, 190, 191, 192, 202, 203, 204, 205, 206],
        'BS6': [38, 39, 40, 58, 59, 60, 76, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138]
    },
    'PLK2': {
        'BS1': [37, 38, 39, 106, 107, 110, 114, 115, 116, 117, 118, 119, 139, 140, 141, 142, 156, 157, 158, 159, 160, 161, 163, 182],
        'BS2': [42, 45, 46, 73, 103, 104, 105, 106, 109, 110, 113, 114],
        'BS3': [32, 33, 35, 37, 53, 117, 118, 119, 120, 121, 122, 123, 124, 163, 165, 176, 178, 184],
        'BS4': [31, 32, 33, 34, 122, 123, 124, 125, 126, 127, 128, 129, 131, 132, 137, 150, 151, 152, 154, 167, 168, 169, 170, 171, 172, 174],
        'BS5': [15, 19, 20, 22, 23, 24, 26, 28, 31, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187, 198, 199, 200, 201, 203, 204, 207],
        'BS6': [53, 54, 55, 64, 65, 66, 70, 71, 72, 75, 77, 79, 80, 81, 82, 124, 125, 126, 127, 128, 129, 130, 131, 132]
    },
    'PLK3': {
        'BS1': [40, 41, 42, 117, 118, 119, 120, 121, 122, 160, 161, 162, 165, 166, 167, 185],
        'BS2': [43, 44, 45, 48, 49, 109, 106, 113, 116, 117],
        'BS3': [29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 123, 124, 125, 126, 127, 128, 130, 132, 134, 136, 154, 156, 171, 172, 177, 179, 187, 188, 189],
        'BS4': [33, 34, 35, 36, 37, 127, 128, 130, 132, 134, 135, 136, 140, 154, 155, 156, 157, 171, 172, 177],  # Removed duplicates
        'BS5': [24, 25, 26, 29, 166, 180, 182, 183, 184, 186, 188, 190, 201, 202, 203],
        'BS6': [59, 73, 75, 77, 78, 80, 88, 90, 91, 92, 129, 131, 132, 133, 134]
    }
}


# Load results from explicit residue analysis
df = pd.read_csv("pvfe_explicit_residues_results.csv")

# Filter successful calculations only
df_success = df[df["diagnostic"] == "SUCCESS"].copy()

# Create pocket ordering consistent with experimental validation status
pocket_order = ["BS1", "BS2", "BS3", "BS4", "BS5", "BS6"]
isoform_order = ["PLK1", "PLK2", "PLK3"]

# Map pocket categories for visual annotation
df_success["pocket_category"] = df_success["pocket"].apply(
    lambda x: "Experimentally\nverified" if x in ["BS1", "BS2", "BS3"] else "Novel\npocket"
)

# Create figure
fig, axes = plt.subplots(1, 3, figsize=(14, 6), sharey=True, dpi=300)
fig.suptitle("Pocket Volume Fluctuation Entropy (Explicit Residues)\nAcross PLK Isoforms", 
             fontsize=14, fontweight='bold', y=1.02)

colors = {"BS1": "#E64B35", "BS2": "#4DBBD5", "BS3": "#00A087", 
          "BS4": "#F39B7F", "BS5": "#8491B4", "BS6": "#91D1C2"}

for idx, isoform in enumerate(isoform_order):
    ax = axes[idx]
    isoform_data = df_success[df_success["isoform"] == isoform]
    
    # Violin plot with embedded boxplot
    sns.violinplot(data=isoform_data, x="pocket", y="pvfe_entropy", 
                   order=pocket_order, ax=ax, palette=[colors[p] for p in pocket_order],
                   inner=None, linewidth=0.8, alpha=0.7)
    
    # Overlay boxplot for quartiles
    sns.boxplot(data=isoform_data, x="pocket", y="pvfe_entropy",
                order=pocket_order, ax=ax, color="white", 
                width=0.3, showfliers=False, linewidth=1.2)
    
    # Overlay individual ligand points
    sns.stripplot(data=isoform_data, x="pocket", y="pvfe_entropy",
                  order=pocket_order, ax=ax, color="black", 
                  size=3, alpha=0.6, jitter=0.15)
    
    # Highlight extreme outliers (>2 std from mean per pocket)
    for pocket in pocket_order:
        pocket_data = isoform_data[isoform_data["pocket"] == pocket]
        if len(pocket_data) > 0:
            mean_val = pocket_data["pvfe_entropy"].mean()
            std_val = pocket_data["pvfe_entropy"].std()
            outliers = pocket_data[
                (pocket_data["pvfe_entropy"] < mean_val - 2*std_val) | 
                (pocket_data["pvfe_entropy"] > mean_val + 2*std_val)
            ]
            for _, row in outliers.iterrows():
                ax.scatter(pocket_order.index(pocket), row["pvfe_entropy"],
                          s=80, facecolors='none', edgecolors='red', 
                          linewidths=1.5, zorder=10)
    
    ax.set_title(isoform, fontsize=12, fontweight='bold')
    ax.set_xlabel("Binding Pocket", fontsize=11)
    if idx == 0:
        ax.set_ylabel("PVFE Entropy (nats)", fontsize=11)
    else:
        ax.set_ylabel("")
    
    ax.set_ylim(2.2, 3.3)
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    ax.tick_params(axis='both', labelsize=9)

# Add annotation for experimental vs novel pockets
for ax in axes:
    ax.axvspan(-0.5, 2.5, alpha=0.1, color='gray', zorder=0, label='_nolegend_')
    ax.axvspan(2.5, 5.5, alpha=0.05, color='gray', zorder=0, label='_nolegend_')

# Custom legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='gray', alpha=0.1, label='Experimentally verified (BS1-BS3)'),
    Patch(facecolor='gray', alpha=0.05, label='Novel pockets (BS4-BS6)'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='none', 
               markeredgecolor='red', markersize=8, label='Extreme outlier (>2σ)')
]
axes[2].legend(handles=legend_elements, loc='lower right', fontsize=8)

plt.tight_layout()
plt.savefig("pvfe_explicit_distribution.png", bbox_inches='tight', dpi=300)
plt.savefig("pvfe_explicit_distribution.pdf", bbox_inches='tight')
plt.show()

# Supplementary statistical table
summary_stats = df_success.groupby(["isoform", "pocket"])["pvfe_entropy"].agg(
    ['count', 'mean', 'std', 'min', 'max']
).round(3).reset_index()
summary_stats.to_csv("pvfe_explicit_summary_statistics.csv", index=False)
print("Summary statistics saved to pvfe_explicit_summary_statistics.csv")

# Print pocket residue counts for methods section
print("\nPocket residue counts per isoform-pocket:")
for isoform in isoform_order:
    print(f"\n{isoform}:")
    for pocket in pocket_order:
        n_res = len(POCKET_RESIDUES.get(isoform, {}).get(pocket, []))
        print(f"  {pocket}: {n_res} residues")
