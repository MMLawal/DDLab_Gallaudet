# water_cons.py
import os
import glob
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import f_oneway, kruskal
from multiprocessing import Pool, cpu_count

# =============================================================================
# ISOFORM AND POCKET DEFINITIONS
# =============================================================================

ISOFORM_RESIDUE_MAP = {
    "PLK1": {"protein": "1-224", "ligand": "225"},
    "PLK2": {"protein": "1-213", "ligand": "214"},
    "PLK3": {"protein": "1-221", "ligand": "222"},
}

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
        'BS4': [33, 34, 35, 36, 37, 127, 128, 130, 132, 134, 135, 136, 140, 154, 155, 156, 157, 171, 172, 177],
        'BS5': [24, 25, 26, 29, 166, 180, 182, 183, 184, 186, 188, 190, 201, 202, 203],
        'BS6': [59, 73, 75, 77, 78, 80, 88, 90, 91, 92, 129, 131, 132, 133, 134]
    }
}


def survival_lifetime(binary_series):
    """Calculate mean lifetime from binary presence series"""
    lifetimes = []
    count = 0

    for v in binary_series:
        if v == 1:
            count += 1
        else:
            if count > 0:
                lifetimes.append(count)
                count = 0

    if count > 0:
        lifetimes.append(count)

    if len(lifetimes) == 0:
        return 0

    return np.mean(lifetimes)


def bootstrap_ci(data, n_boot=2000):
    """Calculate bootstrap confidence intervals"""
    boot = []

    for _ in range(n_boot):
        sample = np.random.choice(data, size=len(data), replace=True)
        boot.append(np.mean(sample))

    mean = np.mean(boot)
    low = np.percentile(boot, 2.5)
    high = np.percentile(boot, 97.5)

    return mean, low, high


def discover_systems(base_dir="."):
    """Discover all PLK-ligand simulation systems"""
    systems = []

    for iso in [1, 2, 3]:
        for bs in range(1, 7):

            iso_dir = f"{base_dir}/plk{iso}-bs{bs}"

            if not os.path.isdir(iso_dir):
                continue

            for lig_dir in glob.glob(f"{iso_dir}/*"):

                dcd = glob.glob(f"{lig_dir}/plk*_bs*.dcd")
                parm = glob.glob(f"{lig_dir}/plk*_bs*.parm7")

                if dcd and parm:
                    systems.append({
                        "isoform": f"PLK{iso}",
                        "pocket": f"BS{bs}",
                        "ligand": os.path.basename(lig_dir),
                        "dcd": dcd[0],
                        "parm": parm[0]
                    })

    print(f"Discovered {len(systems)} systems")
    return systems


def analyze_system(sys):
    """Analyze water residence networks for a single system"""
    print("Processing:", sys["isoform"], sys["pocket"], sys["ligand"])

    u = mda.Universe(sys["parm"], sys["dcd"])

    # Protein / ligand selections
    config = ISOFORM_RESIDUE_MAP[sys["isoform"]]
    ligand_resid = int(config["ligand"])

    protein = u.select_atoms("protein")
    ligand = u.select_atoms(f"resid {ligand_resid}")

    pocket_res = POCKET_RESIDUES[sys["isoform"]][sys["pocket"]]
    pocket_atoms = u.select_atoms(f"resid {' '.join(map(str, pocket_res))}")

    waters = u.select_atoms("resname WAT TIP3 HOH and name O")

    n_frames = len(u.trajectory)

    if n_frames == 0:
        return None

    water_ids = waters.resids
    presence = {wid: np.zeros(n_frames) for wid in water_ids}

    for frame_i, ts in enumerate(u.trajectory):
        d = distance_array(waters.positions, pocket_atoms.positions)
        close = np.any(d < 5.0, axis=1)

        for wid in water_ids[close]:
            presence[wid][frame_i] = 1

    conserved_waters = []
    water_lifetimes = {}

    for wid, series in presence.items():
        occupancy = series.mean()

        if occupancy > 0.5:
            lifetime = survival_lifetime(series)
            conserved_waters.append(wid)
            water_lifetimes[wid] = lifetime

    if len(conserved_waters) == 0:
        return None

    # Build water-mediated interaction network
    G = nx.Graph()
    G.add_node("Ligand", type="ligand")

    for resid in pocket_res:
        G.add_node(f"P{resid}", type="protein")

    for wid in conserved_waters:
        water_node = f"W{wid}"
        G.add_node(water_node, type="water")

        h = HydrogenBondAnalysis(
            universe=u,
            donors_sel=f"resid {wid}",
            acceptors_sel="protein or resid {}".format(ligand_resid),
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150
        )

        h.run()
        hbonds = h.results.hbonds

        if hbonds is None:
            continue

        hbonds = np.array(hbonds)
        binary_lig = np.zeros(n_frames)
        protein_binary = {}

        for row in hbonds:
            frame = int(row[0])
            acceptor = u.atoms[int(row[3])]
            resid = acceptor.resid

            if resid == ligand_resid:
                binary_lig[frame] = 1
            elif resid in pocket_res:
                if resid not in protein_binary:
                    protein_binary[resid] = np.zeros(n_frames)
                protein_binary[resid][frame] = 1

        lig_lifetime = survival_lifetime(binary_lig)

        if lig_lifetime > 0:
            G.add_edge(water_node, "Ligand", weight=lig_lifetime)

        for resid, binary in protein_binary.items():
            lifetime = survival_lifetime(binary)

            if lifetime > 0:
                G.add_edge(water_node, f"P{resid}", weight=lifetime)

    if len(G.edges) == 0:
        return None

    # Calculate network metrics
    degree = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G, weight="weight")
    eigen = nx.eigenvector_centrality(G, max_iter=2000, weight="weight")

    mean_degree = np.mean(list(degree.values()))
    mean_bet = np.mean(list(betweenness.values()))
    mean_eigen = np.mean(list(eigen.values()))

    communities = nx.algorithms.community.greedy_modularity_communities(G)
    modularity = nx.algorithms.community.modularity(G, communities)
    clustering = nx.average_clustering(G, weight="weight")

    edge_bet = nx.edge_betweenness_centrality(G, weight="weight")

    water_bridges = []
    for (n1, n2), val in edge_bet.items():
        if n1.startswith("W") or n2.startswith("W"):
            water_bridges.append(val)

    mean_bridge = np.mean(water_bridges)
    max_bridge = np.max(water_bridges)

    return {
        "isoform": sys["isoform"],
        "pocket": sys["pocket"],
        "ligand": sys["ligand"],
        "n_conserved_waters": len(conserved_waters),
        "mean_water_lifetime": np.mean(list(water_lifetimes.values())),
        "network_nodes": len(G.nodes),
        "network_edges": len(G.edges),
        "mean_degree": mean_degree,
        "mean_betweenness": mean_bet,
        "mean_eigenvector": mean_eigen,
        "modularity": modularity,
        "weighted_clustering": clustering,
        "mean_bridge_betweenness": mean_bridge,
        "max_bridge_betweenness": max_bridge
    }


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    systems = discover_systems()

    with Pool(cpu_count()) as pool:
        results = pool.map(analyze_system, systems)

    results = [r for r in results if r is not None]

    df = pd.DataFrame(results)
    df.to_csv("water_residence_network_results.csv", index=False)

    for metric in [
        "mean_water_lifetime",
        "weighted_clustering",
        "mean_bridge_betweenness"
    ]:
        groups = [
            df[df.isoform == iso][metric]
            for iso in df.isoform.unique()
        ]

        F, p = f_oneway(*groups)
        print(metric, "ANOVA p =", p)


plt.figure(figsize=(6,5))

data = [df[df.isoform==iso]["weighted_clustering"]
        for iso in sorted(df.isoform.unique())]

plt.boxplot(data, labels=sorted(df.isoform.unique()))

plt.ylabel("Weighted Clustering")
plt.title("Network Compactness Across PLK Isoforms")

plt.savefig("network_compactness_isoforms.png", dpi=300)        


plt.figure(figsize=(6,5))

data = [df[df.isoform==iso]["mean_bridge_betweenness"]
        for iso in sorted(df.isoform.unique())]

plt.boxplot(data, labels=sorted(df.isoform.unique()))

plt.ylabel("Bridge Betweenness")
plt.title("Critical Water Bridge Centrality")

plt.savefig("water_bridge_centrality.png", dpi=300)

plt.figure(figsize=(6,5))

data = [df[df.isoform==iso]["mean_water_lifetime"]
        for iso in sorted(df.isoform.unique())]

plt.boxplot(data, labels=sorted(df.isoform.unique()))

plt.ylabel("Residence Lifetime (frames)")
plt.title("Isoform Comparison of Water Residence Time")

plt.savefig("water_residence_lifetime.png", dpi=300)
