import pandas as pd
import numpy as np

# === Load dataset ===
df = pd.read_csv("admetlab_met_ime.csv")

# === Convert numeric columns safely ===
for col in df.columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# === Define desirable descriptor limits ===
criteria = {
    'MW': (100, 600),
    'nHA': (0, 12),
    'nHD': (0, 7),
    'nRot': (0, 11),
    'nRing': (0, 6),
    'MaxRing': (0, 18),
    'nHet': (1, 15),
    'fChar': (-4, 4),
    'nRig': (0, 30),
    'nStereo': (None, 2),
    'TPSA': (20, 140),
    'logP': (-0.7, 5),
    'logS': (-6, 0),
    'logD': (1, 3),
    'pka_basic': (2, 10),
    'pka_acidic': (2, 12),
    'mp': (50, 250),
    'Flex': (0, 9),
    'QED': ('>', 0.67),
    'Synth': ('<', 6),
    'gasa': (0, 0),
    'Fsp3': ('>=', 0.42),
    'MCE-18': ('>=', 45),
    'Natural Product-likeness': (-5, 5),
    'Pfizer': (0, 0),
    'GSK': (0, 0),
    'GoldenTriangle': (0, 0),
    'PAINS': (0, 0),
    'Alarm_NMR': ('<=', 3),
    'BMS': (0, 0),
    'Aggregators': ('<', 6),
    'Fluc': ('<', 6),
    'Blue_fluorescence': ('<', 6),
    'Green_fluorescence': ('<', 6),
    'Reactive': ('<', 6),
    'Promiscuous': ('<', 6),
    'Ghose': (0, 0),
    'Veber': (0, 0),
    'Egan': (0, 0),
    'Muegge': (0, 0),
    'hia': ('<', 0.45),
    'caco2': ('>=', -5.15),
    'PAMPA': ('<', 0.6),
    'pgp_inh': ('<', 0.6),
    'pgp_sub': ('<', 0.6),
    'f20': ('<', 0.6),
    'f30': ('<', 0.6),
    'f50': ('<', 0.6),
    'OATP1B1': ('<', 0.6),
    'OATP1B3': ('<', 0.6),
    'BSEP': ('<', 0.6),
    'MRP1': ('<', 0.6),
    'PPB': ('<', 90),
    'logVDss': (0.04, 20),
    'Fu': ('>', 20),
    'log Kp': (-13.89, 0),
    'CYP1A2-inh': ('<', 0.6),
    'CYP1A2-sub': ('<', 0.6),
    'CYP2C19-inh': ('<', 0.6),
    'CYP2C19-sub': ('<', 0.6),
    'CYP2C9-inh': ('<', 0.6),
    'CYP2C9-sub': ('<', 0.6),
    'CYP2D6-inh': ('<', 0.6),
    'CYP2D6-sub': ('<', 0.6),
    'CYP3A4-inh': ('<', 0.6),
    'CYP3A4-sub': ('<', 0.6),
    'CYP2B6-inh': ('<', 0.6),
    'CYP2B6-sub': ('<', 0.6),
    'CYP2C8-inh': ('<', 0.6),
    'LM-Human': ('>', 0.3),
    'cl-plasma': ('<', 10),
    'hERG': ('<', 0.7),
    'hERG-10um': ('<', 0.7),
    'DILI': ('<', 0.7),
    'Ames': ('<', 0.7),
    'ROA': ('<', 0.7),
    'FDAMDD': ('<', 0.7),
    'SkinSen': ('<', 0.7),
    'Carcinogenicity': ('<', 0.7),
    'EC': ('<', 0.7),
    'EI': ('<', 0.7),
    'Respiratory': ('<', 0.7),
    'H-HT': ('<', 0.7),
    'Neurotoxicity-DI': ('<', 0.7),
    'Ototoxicity': ('<', 0.7),
    'Hematotoxicity': ('<', 0.7),
    'Nephrotoxicity-DI': ('<', 0.7),
    'Genotoxicity': ('<', 0.7),
    'RPMI-8226': ('<', 0.7),
    'A549': ('<', 0.7),
    'HEK293': ('<', 0.7),
    'NR-AhR': ('<', 0.7),
    'NR-AR': ('<', 0.7),
    'NR-AR-LBD': ('<', 0.7),
    'NR-Aromatase': ('<', 0.7),
    'NR-ER': ('<', 0.7),
    'NR-ER-LBD': ('<', 0.7),
    'NR-PPAR-gamma': ('<', 0.7),
    'SR-ARE': ('<', 0.7),
    'SR-ATAD5': ('<', 0.7),
    'SR-HSE': ('<', 0.7),
    'SR-MMP': ('<', 0.7),
    'SR-p53': ('<', 0.7),
    'NonBiodegradable': ('<=', 5),
    'NonGenotoxic_Carcinogenicity': ('<=', 5),
    'SureChEMBL': ('<=', 5),
    'LD50_oral': ('<=', 5),
    'Skin_Sensitization': ('<=', 5),
    'Acute_Aquatic_Toxicity': ('<=', 5),
    'FAF-Drugs4 Rule': ('<=', 5),
    'Brenk': ('<=', 1),
    'Leadlikeness': ('<=', 1),
    'Genotoxic_Carcinogenicity_Mutagenicity': ('<=', 5)
}

# === Evaluate criteria ===
results = pd.DataFrame(index=df.index)

for desc, rule in criteria.items():
    if desc not in df.columns:
        continue  # skip if not found
    val = df[desc]
    if val.isnull().all():
        continue  # skip columns entirely missing or non-numeric
    
    # Two-sided range
    if isinstance(rule, tuple) and len(rule) == 2 and isinstance(rule[0], (int, float, type(None))):
        low, high = rule
        if low is None:
            results[desc] = val <= high
        elif high is None:
            results[desc] = val >= low
        else:
            results[desc] = (val >= low) & (val <= high)
    else:
        op, thr = rule
        if op == '<': results[desc] = val < thr
        elif op == '<=': results[desc] = val <= thr
        elif op == '>': results[desc] = val > thr
        elif op == '>=': results[desc] = val >= thr
        elif op in ['=', '==']: results[desc] = val == thr
        else: results[desc] = False

# === Count how many desirable descriptors each compound satisfies ===
df["n_passed"] = results.sum(axis=1)

# === Filter compounds with ≥70 satisfied descriptors ===
filtered = df[df["n_passed"] >= 70]

# === Save filtered compounds ===
filtered.to_csv("filtered_contr_70.csv", index=False)

print(f"✅ Total compounds: {len(df)}")
print(f"✅ Compounds meeting ≥70 desirable descriptors: {len(filtered)}")
print("💾 Saved as 'filtered_contr_70.csv'")

