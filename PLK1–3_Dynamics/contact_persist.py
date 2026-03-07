import os
import glob
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from multiprocessing import Pool, cpu_count

def compute_contact_persistence(system):

    isoform = system["isoform"]
    pocket = system["pocket"]
    ligand_name = system["ligand"]
    parm = system["parm"]
    dcd = system["dcd"]

    try:
        u = mda.Universe(parm, dcd)
    except:
        return None

    n_frames = len(u.trajectory)
    if n_frames == 0:
        return None

    # Ligand selection
    ligand = u.select_atoms("resname UNK UNL VIH LIG")
    if len(ligand) == 0:
        return None

    # Pocket residues
    residues = POCKET_RESIDUES[isoform][pocket]

    results = []

    for resid in residues:

        res_atoms = u.select_atoms(f"resid {resid}")
        contact_frames = 0

        for ts in u.trajectory:

            d = distance_array(res_atoms.positions,
                               ligand.positions)

            if np.min(d) <= 5.0:
                contact_frames += 1

        persistence = contact_frames / n_frames

        results.append({
            "isoform": isoform,
            "pocket": pocket,
            "ligand": ligand_name,
            "resid": resid,
            "persistence": persistence
        })

    return results

def build_systems():

    systems = []

    for iso in ["PLK1","PLK2","PLK3"]:
        for pocket in ["BS1","BS2","BS3","BS4","BS5","BS6"]:

            pattern = f"./plk{iso[-1]}-{pocket.lower()}/"

            for ligand_dir in glob.glob(pattern + "*"):

                ligand = os.path.basename(ligand_dir)

                dcd_files = glob.glob(ligand_dir + "/plk*.dcd")
                parm_files = glob.glob(ligand_dir + "/plk*.parm7")

                if len(dcd_files) == 0 or len(parm_files) == 0:
                    continue

                systems.append({
                    "isoform": iso,
                    "pocket": pocket,
                    "ligand": ligand,
                    "parm": parm_files[0],
                    "dcd": dcd_files[0]
                })

    return systems

def run_parallel():

    systems = build_systems()
    print("Total systems:", len(systems))

    with Pool(cpu_count()) as pool:
        results = pool.map(compute_contact_persistence, systems)

    flat = []
    for r in results:
        if r is not None:
            flat.extend(r)

    df = pd.DataFrame(flat)
    df.to_csv("contact_persistence_all.csv", index=False)

    return df

def plot_heatmaps(df):

    os.makedirs("contact_heatmaps", exist_ok=True)

    for iso in df["isoform"].unique():

        sub = df[df["isoform"] == iso]

        pivot = sub.pivot_table(index="resid",
                                columns="ligand",
                                values="persistence",
                                aggfunc="mean")

        plt.figure(figsize=(12, 8))
        sns.heatmap(pivot,
                    cmap="viridis",
                    vmin=0,
                    vmax=1)

        plt.title(f"{iso} Contact Persistence Landscape")
        plt.xlabel("Ligand")
        plt.ylabel("Residue")

        plt.tight_layout()
        plt.savefig(f"contact_heatmaps/{iso}_heatmap.png",
                    dpi=300)
        plt.close()

if __name__ == "__main__":

    df = run_parallel()

    plot_heatmaps(df)

    anchors = identify_anchor_residues(df)

    print("Pipeline Complete.")
