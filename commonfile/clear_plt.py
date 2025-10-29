import pandas as pd
import matplotlib.pyplot as plt

# === Load data ===
df = pd.read_csv("clearance_50_supern.csv")

# Ensure column names are correct (you can adjust if needed)
# Expected columns: "Ligands", "cl-plasma", "t0.5"
df.columns = [col.strip().lower().replace(" ", "_") for col in df.columns]

# === Create figure and twin y-axis ===
fig, ax1 = plt.subplots(figsize=(10, 6))

# --- Plot CL-plasma (Left Y-axis) ---
ax1.plot(df.index + 1, df["cl-plasma"], color='tab:blue', marker='o', label='CL-plasma')
ax1.set_xlabel("Ligands", fontsize=8, weight='bold')
ax1.set_ylabel("CL-plasma (mL/min/kg)", color='tab:blue', fontsize=8, weight='bold')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Add red dotted line at CL-plasma = 5
ax1.axhline(y=5, color='red', linestyle='--', linewidth=2, label='Optimal Clearance (<5)')

# --- Plot t0.5 (Right Y-axis) ---
ax2 = ax1.twinx()
ax2.plot(df.index + 1, df["t0.5"], color='tab:green', marker='s', label='t₀.₅ (Half-life)')
ax2.set_ylabel("t₀.₅ (hours)", color='tab:green', fontsize=8, weight='bold')
ax2.tick_params(axis='y', labelcolor='tab:green')

# === Styling ===
#plt.title("CL-plasma and t₀.₅ Profile for 50 Ligands", fontsize=14, weight='bold')
ax1.set_xticks(range(1, len(df) + 1))
ax1.set_xticklabels(df["ligands"] if "ligands" in df.columns else [f"Lig {i}" for i in range(1, len(df) + 1)], rotation=90, fontsize=8)

# Combine legends from both axes
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=9)

plt.tight_layout()
plt.show()

