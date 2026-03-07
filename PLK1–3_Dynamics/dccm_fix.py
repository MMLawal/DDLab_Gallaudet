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

# ==============================
# PARAMETERS
# ==============================

DCCM_STRIDE = 1
MIN_FRAMES_DCCM = 50
CORRELATION_THRESHOLD = 0.3
OUTPUT_DIR = "dccm_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

ISO_LIST = ["PLK1", "PLK2", "PLK3"]

# ==============================
# UTILS
# ==============================

def file_exists(fname):
    return os.path.exists(os.path.join(OUTPUT_DIR, fname))


# ==============================
# SYSTEM DISCOVERY
# ==============================

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
                })

    print(f"Discovered {len(systems)} systems.")
    return systems


# ==============================
# DCCM
# ==============================

def compute_dccm(u, selection="protein and name CA"):

    atoms = u.select_atoms(selection)
    coords = np.array([atoms.positions.copy()
                       for ts in u.trajectory[::DCCM_STRIDE]])

    if coords.shape[0] < MIN_FRAMES_DCCM:
        raise ValueError("Insufficient frames for DCCM")

    fluctuations = coords - coords.mean(axis=0)
    n = fluctuations.shape[1]

    C = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):

            num = np.sum(np.sum(fluctuations[:, i] *
                                fluctuations[:, j], axis=1))

            denom = np.sqrt(
                np.sum(np.sum(fluctuations[:, i]**2, axis=1)) *
                np.sum(np.sum(fluctuations[:, j]**2, axis=1))
            )

            val = num / denom if denom != 0 else 0
            C[i, j] = C[j, i] = val

    return C, atoms.resids.copy()


# ==============================
# NETWORK + COMMUNITIES
# ==============================

def build_network(C):

    G = nx.Graph()
    n = C.shape[0]

    for i in range(n):
        for j in range(i+1, n):
            if abs(C[i, j]) >= CORRELATION_THRESHOLD:
                G.add_edge(i, j, weight=abs(C[i, j]))

    return G


def detect_communities(G):

    if G.number_of_edges() == 0:
        return []

    comp = nx.algorithms.community.girvan_newman(G)
    return list(next(comp))


# ==============================
# PLOTTING
# ==============================

def save_heatmap(matrix, title, fname, vmin=-1, vmax=1):

    plt.figure(figsize=(8,7))
    sns.heatmap(matrix, cmap="coolwarm", center=0,
                vmin=vmin, vmax=vmax,
                cbar_kws={"label": "Correlation"})
    plt.title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300)
    plt.close()


def save_difference_map(C1, res1, C2, res2, label1, label2):

    common = sorted(set(res1).intersection(set(res2)))
    if len(common) == 0:
        print(f"No overlap between {label1}-{label2}")
        return

    idx1 = [np.where(res1 == r)[0][0] for r in common]
    idx2 = [np.where(res2 == r)[0][0] for r in common]

    diff = C1[np.ix_(idx1, idx1)] - C2[np.ix_(idx2, idx2)]

    save_heatmap(diff,
                 f"{label1} – {label2} DCCM Difference",
                 f"{label1}_minus_{label2}_DCCM.png",
                 vmin=-0.5, vmax=0.5)


def save_masked_dccm(mean, mask, isoform):

    masked = np.where(mask, mean, 0)
    save_heatmap(masked,
                 f"{isoform} Significant DCCM",
                 f"{isoform}_Significant_DCCM.png")


def save_dccm_with_communities(C, isoform):

    G = build_network(C)
    communities = detect_communities(G)

    plt.figure(figsize=(8,7))
    sns.heatmap(C, cmap="coolwarm", center=0,
                vmin=-1, vmax=1)

    for comm in communities:
        idx = sorted(comm)
        start, end = min(idx), max(idx)
        plt.gca().add_patch(
            plt.Rectangle((start, start),
                          end-start,
                          end-start,
                          fill=False,
                          edgecolor="black",
                          linewidth=2)
        )

    plt.title(f"{isoform} DCCM with Communities")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR,
                             f"{isoform}_DCCM_Communities.png"),
                dpi=300)
    plt.close()


def compute_correlation_decay(C):

    n = C.shape[0]
    decay = []

    for d in range(1, n):
        vals = [abs(C[i, i+d]) for i in range(n-d)]
        decay.append(np.mean(vals))

    return np.array(decay)


def save_decay_plot(decay_dict):

    plt.figure(figsize=(7,6))
    for iso, decay in decay_dict.items():
        plt.plot(decay, label=iso)

    plt.xlabel("Residue Separation |i-j|")
    plt.ylabel("Mean |Correlation|")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR,
                             "Correlation_Decay_Per_Isoform.png"),
                dpi=300)
    plt.close()


# ==============================
# BOOTSTRAP
# ==============================

def bootstrap_significance(matrices, n_boot=200):

    stack = np.stack(matrices)
    mean = stack.mean(axis=0)

    boot_means = []
    for _ in range(n_boot):
        sample = stack[np.random.choice(len(stack),
                                        len(stack),
                                        replace=True)]
        boot_means.append(sample.mean(axis=0))

    boot_means = np.stack(boot_means)
    lower = np.percentile(boot_means, 2.5, axis=0)
    upper = np.percentile(boot_means, 97.5, axis=0)

    mask = (mean < lower) | (mean > upper)
    return mean, mask


# ==============================
# MAIN PIPELINE
# ==============================

def run_pipeline():

    systems = discover_systems()

    iso_storage = {iso: [] for iso in ISO_LIST}
    iso_resids = {}

    for sys in tqdm(systems):

        u = mda.Universe(sys["prmtop"], sys["dcd"])
        C, resids = compute_dccm(u)

        iso_storage[sys["isoform"]].append(C)
        iso_resids[sys["isoform"]] = resids

    # ---- Average matrices
    avg = {}
    for iso in ISO_LIST:
        if iso_storage[iso]:
            avg[iso] = np.mean(np.stack(iso_storage[iso]), axis=0)

    # ---- Save average DCCMs
    for iso in avg:
        fname = f"{iso}_Average_DCCM.png"
        if not file_exists(fname):
            save_heatmap(avg[iso],
                         f"{iso} Average DCCM",
                         fname)

    # ---- Difference maps
    for a, b in itertools.combinations(ISO_LIST, 2):
        fname = f"{a}_minus_{b}_DCCM.png"
        if not file_exists(fname):
            save_difference_map(avg[a], iso_resids[a],
                                avg[b], iso_resids[b],
                                a, b)

    # ---- Masked DCCM
    for iso in avg:
        fname = f"{iso}_Significant_DCCM.png"
        if not file_exists(fname):
            mean, mask = bootstrap_significance(iso_storage[iso])
            save_masked_dccm(mean, mask, iso)

    # ---- Communities
    for iso in avg:
        fname = f"{iso}_DCCM_Communities.png"
        if not file_exists(fname):
            save_dccm_with_communities(avg[iso], iso)

    # ---- Correlation decay
    if not file_exists("Correlation_Decay_Per_Isoform.png"):
        decay_dict = {
            iso: compute_correlation_decay(avg[iso])
            for iso in avg
        }
        save_decay_plot(decay_dict)


# RUN
run_pipeline()
