import os
import glob
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import MDAnalysis as mda
from tqdm import tqdm

DCCM_STRIDE = 1
MIN_FRAMES_DCCM = 50
CORRELATION_THRESHOLD = 0.3
OUTPUT_DIR = "dccm_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
        'BS4': [33, 34, 35, 36, 37, 127, 128, 130, 132, 134, 135, 136, 37, 140, 154, 155, 156, 157, 171, 172, 134, 177],
        'BS5': [24, 25, 26, 29, 166, 180, 182, 183, 184, 186, 188, 190, 201, 202, 203],
        'BS6': [59, 73, 75, 77, 78, 80, 88, 90, 91, 92, 129, 131, 132, 133, 134]
    }
}

def discover_systems(base_dir="."):

    systems = []
    iso_dirs = glob.glob(os.path.join(base_dir, "plk[1-3]-bs[1-6]"))

    for iso_dir in iso_dirs:

        iso_pocket = os.path.basename(iso_dir)
        isoform_raw, pocket = iso_pocket.split("-")

        ligand_dirs = glob.glob(os.path.join(iso_dir, "*"))

        for lig_dir in ligand_dirs:

            prmtop = glob.glob(os.path.join(lig_dir, "plk*_bs*.parm7"))
            dcd = glob.glob(os.path.join(lig_dir, "plk*_bs*.dcd"))

            if prmtop and dcd:

                systems.append({
                    "prmtop": prmtop[0],
                    "dcd": dcd[0],
                    "isoform": f"PLK{isoform_raw[3:]}",
                    "pocket": pocket.upper(),
                    "ligand": os.path.basename(lig_dir)
                })

    print(f"Discovered {len(systems)} systems.")
    return systems

def compute_dccm(u, selection="protein and name CA"):

    atoms = u.select_atoms(selection)
    n_atoms = len(atoms)

    coords = []
    for ts in u.trajectory[::DCCM_STRIDE]:
        coords.append(atoms.positions.copy())

    coords = np.array(coords)

    if coords.shape[0] < MIN_FRAMES_DCCM:
        raise ValueError("Insufficient frames for DCCM")

    mean_pos = coords.mean(axis=0)
    fluctuations = coords - mean_pos

    norms = np.zeros(n_atoms)
    for i in range(n_atoms):
        norms[i] = np.sum(np.sum(fluctuations[:, i]**2, axis=1))

    C = np.zeros((n_atoms, n_atoms))

    for i in range(n_atoms):
        for j in range(i, n_atoms):

            num = np.sum(
                np.sum(fluctuations[:, i] * fluctuations[:, j], axis=1)
            )

            denom = np.sqrt(norms[i] * norms[j])
            val = num / denom if denom != 0 else 0

            C[i, j] = val
            C[j, i] = val

    return C, atoms.resids.copy()

def build_network(C):

    G = nx.Graph()
    n = C.shape[0]

    for i in range(n):
        G.add_node(i)

    for i in range(n):
        for j in range(i+1, n):
            if abs(C[i, j]) >= CORRELATION_THRESHOLD:
                G.add_edge(i, j, weight=abs(C[i, j]))

    return G

def detect_communities(G):

    if G.number_of_edges() == 0:
        return []

    comp = nx.algorithms.community.girvan_newman(G)
    communities = next(comp)
    return list(communities)

def compute_inter_pocket_coupling(C, resids, pocket_dict):

    coupling = {}

    for p1, p2 in itertools.combinations(pocket_dict.keys(), 2):

        idx1 = [i for i, r in enumerate(resids) if r in pocket_dict[p1]]
        idx2 = [i for i, r in enumerate(resids) if r in pocket_dict[p2]]

        if not idx1 or not idx2:
            coupling[f"{p1}_{p2}"] = np.nan
            continue

        sub = np.abs(C[np.ix_(idx1, idx2)])
        coupling[f"{p1}_{p2}"] = sub.mean()

    return coupling

def save_isoform_dccm(C, isoform):

    plt.figure(figsize=(8,7))

    sns.heatmap(
        C,
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        cbar_kws={"label": "Correlation"}
    )

    plt.title(f"{isoform} Average DCCM")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/{isoform}_Average_DCCM.png",
        dpi=300
    )

    plt.close()

def run_pipeline_isoform_level():

    systems = discover_systems()

    results = []

    # Store DCCMs by isoform
    isoform_dccm_storage = {
        "PLK1": [],
        "PLK2": [],
        "PLK3": []
    }

    isoform_resids = {
        "PLK1": None,
        "PLK2": None,
        "PLK3": None
    }

    for sys in tqdm(systems):

        print("Processing:", sys["isoform"], sys["pocket"])

        u = mda.Universe(sys["prmtop"], sys["dcd"])

        C, resids = compute_dccm(u)

        # Store matrix for averaging
        isoform_dccm_storage[sys["isoform"]].append(C)

        # Save residue indexing (same across systems of isoform)
        isoform_resids[sys["isoform"]] = resids

        # --- Network + coupling (unchanged) ---
        G = build_network(C)
        communities = detect_communities(G)

        coupling = compute_inter_pocket_coupling(
            C,
            resids,
            POCKET_RESIDUES[sys["isoform"]]
        )

        row = {
            "isoform": sys["isoform"],
            "pocket": sys["pocket"],
            "n_communities": len(communities)
        }

        row.update(coupling)
        results.append(row)

    df = pd.DataFrame(results)
    df.to_csv(f"{OUTPUT_DIR}/inter_pocket_coupling_all.csv", index=False)

    # =====================================================
    # ISOFORM-LEVEL AVERAGED DCCM
    # =====================================================

    for iso in ["PLK1", "PLK2", "PLK3"]:

        matrices = isoform_dccm_storage[iso]

        if len(matrices) == 0:
            continue

        # Stack and average
        avg_dccm = np.mean(np.stack(matrices), axis=0)
        avg_plk1 = np.mean(np.stack(isoform_dccm_storage["PLK1"]), axis=0)
        avg_plk2 = np.mean(np.stack(isoform_dccm_storage["PLK2"]), axis=0)
        avg_plk3 = np.mean(np.stack(isoform_dccm_storage["PLK3"]), axis=0)

        save_isoform_dccm(avg_dccm, iso)

    # Plot coupling grid
    plot_all_couplings(df)

    return df

def plot_all_couplings(df):

    coupling_cols = [c for c in df.columns if "BS" in c]

    fig, axes = plt.subplots(5, 3, figsize=(15, 18))
    axes = axes.flatten()

    for i, col in enumerate(coupling_cols):

        sns.boxplot(data=df,
                    x="isoform",
                    y=col,
                    ax=axes[i])

        axes[i].set_title(col)

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/All_15_InterPocket_Couplings.png",
                dpi=300)
    plt.close()

df_results = run_pipeline_isoform_level()

def save_difference_map(C1, C2, label1, label2):

    diff = C1 - C2

    plt.figure(figsize=(8,7))
    sns.heatmap(
        diff,
        cmap="bwr",
        center=0,
        vmin=-0.5,
        vmax=0.5,
        cbar_kws={"label":"Δ Correlation"}
    )

    save_difference_map(avg_plk1, avg_plk3, "PLK1", "PLK3")
    save_difference_map(avg_plk1, avg_plk2, "PLK1", "PLK2")
    save_difference_map(avg_plk2, avg_plk3, "PLK2", "PLK3")

    
    plt.title(f"{label1} – {label2} DCCM Difference")
    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/{label1}_minus_{label2}_DCCM.png",
        dpi=300
    )
    plt.close()

def bootstrap_significance(matrices, n_boot=200):

    stack = np.stack(matrices)
    n = stack.shape[1]

    mean = stack.mean(axis=0)

    boot_means = []

    for _ in range(n_boot):
        sample = stack[np.random.choice(len(stack), len(stack), replace=True)]
        boot_means.append(sample.mean(axis=0))

    boot_means = np.stack(boot_means)

    lower = np.percentile(boot_means, 2.5, axis=0)
    upper = np.percentile(boot_means, 97.5, axis=0)
    mean1, mask1 = bootstrap_significance(isoform_dccm_storage["PLK1"])
    mean2, mask2 = bootstrap_significance(isoform_dccm_storage["PLK2"])
    mean3, mask3 = bootstrap_significance(isoform_dccm_storage["PLK3"])

    mask = np.logical_or(mean < lower, mean > upper)

    return mean, mask

def save_masked_dccm(mean, mask, isoform):

    masked = np.where(mask, mean, 0)

    plt.figure(figsize=(8,7))
    sns.heatmap(
        masked,
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        cbar_kws={"label":"Significant Correlation"}
    )

    plt.title(f"{isoform} Significant DCCM")
    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/{isoform}_Significant_DCCM.png",
        dpi=300
    )
    plt.close()
    save_masked_dccm(mean1, mask1, "PLK1")
    save_masked_dccm(mean2, mask2, "PLK2")
    save_masked_dccm(mean3, mask3, "PLK3")

def save_dccm_with_communities(C, isoform):

    G = build_network(C)
    communities = detect_communities(G)

    plt.figure(figsize=(8,7))
    sns.heatmap(C, cmap="coolwarm", center=0,
                vmin=-1, vmax=1)

    for comm in communities:
        indices = sorted(list(comm))
        start = min(indices)
        end = max(indices)

        plt.gca().add_patch(
            plt.Rectangle(
                (start, start),
                end-start,
                end-start,
                fill=False,
                edgecolor="black",
                linewidth=2
            )
        )
    save_dccm_with_communities(avg_plk1, "PLK1")
    save_dccm_with_communities(avg_plk2, "PLK2")
    save_dccm_with_communities(avg_plk3, "PLK3") 

    plt.title(f"{isoform} DCCM with Communities")
    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/{isoform}_DCCM_Communities.png",
        dpi=300
    )
    plt.close() 


def compute_correlation_decay(C):

    n = C.shape[0]
    decay = []

    for d in range(1, n):

        vals = []
        for i in range(n-d):
            vals.append(abs(C[i, i+d]))

        if len(vals) > 0:
            decay.append(np.mean(vals))

    return np.array(decay) 

def save_decay_plot(decay_dict):

    plt.figure(figsize=(7,6))

    for iso, decay in decay_dict.items():
        plt.plot(decay, label=iso)

    decay_dict = {
        "PLK1": compute_correlation_decay(avg_plk1),
        "PLK2": compute_correlation_decay(avg_plk2),
        "PLK3": compute_correlation_decay(avg_plk3)
    }

    save_decay_plot(decay_dict)
    plt.xlabel("Residue Separation |i - j|")
    plt.ylabel("Mean |Correlation|")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/Correlation_Decay_Per_Isoform.png",
        dpi=300
    )
    plt.close()
