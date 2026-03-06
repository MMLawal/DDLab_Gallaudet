#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from multiprocessing import Pool, cpu_count

# ============================================================
# CONSTANTS
# ============================================================

BETA = 0.92  # 1/(kT) at 300 K in kcal/mol
WINDOW = 50  # window size for interaction entropy
NPROC = max(cpu_count()-2,1)

# ============================================================
# SYSTEM DISCOVERY
# ============================================================

def discover_systems(base_dir="."):

    systems = []

    for iso_pocket in sorted(glob.glob(os.path.join(base_dir,"plk*-bs*"))):

        iso_raw, pocket_raw = os.path.basename(iso_pocket).split("-")

        isoform = iso_raw.upper()
        pocket = pocket_raw.upper()

        for lig_dir in glob.glob(os.path.join(iso_pocket,"*")):

            ligand = os.path.basename(lig_dir)

            prmtop = glob.glob(os.path.join(lig_dir,"plk*_bs*.parm7"))
            dcd = glob.glob(os.path.join(lig_dir,"plk*_bs*.dcd"))

            if not prmtop or not dcd:
                continue

            systems.append({
                "isoform": isoform,
                "pocket": pocket,
                "ligand": ligand,
                "prmtop": prmtop[0],
                "traj": dcd[0]
            })

    print(f"\n✓ Found {len(systems)} simulation systems\n")

    return systems

# ============================================================
# LOAD LIE ENERGY FILE
# ============================================================

def load_lie_energy(file):

    data = np.loadtxt(file)

    vdw = data[:,1]
    ele = data[:,2]

    return vdw + ele


# ============================================================
# INTERACTION ENTROPY
# ============================================================

def interaction_entropy(energy_series):

    IE_values = []

    for i in range(WINDOW,len(energy_series)):

        window = energy_series[i-WINDOW:i]

        dE = window - np.mean(window)

        val = np.log(np.mean(np.exp(BETA*dE)))

        IE_values.append(val)

    IE = np.mean(IE_values)

    return IE


# ============================================================
# PROCESS SINGLE SYSTEM
# ============================================================

def process_system(system):

    iso = system["isoform"]
    pocket = system["pocket"]
    ligand = system["ligand"]

    lie_file = f"lie_outputs/{iso}_{pocket}_{ligand}_lie.dat"

    if not os.path.exists(lie_file):
        return None

    energies = load_lie_energy(lie_file)

    LIE = np.mean(energies)

    IE = interaction_entropy(energies)

    dG_pred = LIE + IE

    return {
        "isoform":iso,
        "pocket":pocket,
        "ligand":ligand,
        "LIE_kcal_mol":LIE,
        "IE_kcal_mol":IE,
        "dG_pred":dG_pred
    }


# ============================================================
# PARALLEL EXECUTION
# ============================================================

def run_pipeline():

    systems = discover_systems()

    with Pool(NPROC) as pool:
        results = pool.map(process_system,systems)

    results = [r for r in results if r is not None]

    df = pd.DataFrame(results)

    df.to_csv("plk_energy_results.csv",index=False)

    print("✓ Energy analysis complete")

    return df


# ============================================================
# MERGE WITH EXPERIMENTAL DATA
# ============================================================

def merge_experiment(df):

    exp = pd.read_csv("binding_data1.csv")

    merged = pd.merge(
        df,
        exp,
        on=["ligand","isoform","pocket"],
        how="inner"
    )

    merged["exp_dG"] = merged["experimental_value"]

    merged.to_csv("plk_energy_experiment_merged.csv",index=False)

    print(f"✓ {len(merged)} systems matched with experiment")

    return merged


# ============================================================
# CORRELATION ANALYSIS
# ============================================================

def correlation_plots(df):

    os.makedirs("plots",exist_ok=True)

    r,p = pearsonr(df["dG_pred"],df["exp_dG"])

    plt.figure(figsize=(6,6))

    sns.regplot(
        x="exp_dG",
        y="dG_pred",
        data=df,
        scatter_kws={"s":60}
    )

    plt.xlabel("Experimental Binding (kcal/mol)")
    plt.ylabel("Predicted ΔG (LIE + IE)")
    plt.title(f"Free Energy Correlation\nr = {r:.2f}")

    plt.tight_layout()

    plt.savefig("plots/dG_correlation.png",dpi=300)

    print(f"✓ Correlation r = {r:.3f}")


# ============================================================
# ISOFORM ANALYSIS
# ============================================================

def isoform_statistics(df):

    plt.figure(figsize=(8,6))

    sns.boxplot(
        x="pocket",
        y="IE_kcal_mol",
        hue="isoform",
        data=df
    )

    plt.ylabel("Interaction Entropy (kcal/mol)")

    plt.tight_layout()

    plt.savefig("plots/IE_by_isoform.png",dpi=300)


# ============================================================
# MAIN
# ============================================================

def main():

    energy_df = run_pipeline()

    merged = merge_experiment(energy_df)

    correlation_plots(merged)

    isoform_statistics(energy_df)

    print("\n✓ PLK Analysis Pipeline Complete\n")


if __name__ == "__main__":
    main()
