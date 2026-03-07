import os
import glob
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

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

# ============================================================
# PARAMETERS
# ============================================================

GRID_SPACING = 1.0      # Å grid resolution
POCKET_CUTOFF = 5.0     # Å pocket definition
KT = 0.593              # kcal/mol at 300 K
BULK_DENSITY = 0.0334   # waters per Å^3

# ============================================================
# SYSTEM DISCOVERY
# ============================================================

def discover_systems(base_dir="."):

    systems = []

    for iso in [1,2,3]:
        for bs in range(1,7):

            iso_dir = f"{base_dir}/plk{iso}-bs{bs}"

            if not os.path.isdir(iso_dir):
                continue

            for lig_dir in glob.glob(f"{iso_dir}/*"):

                parm = glob.glob(f"{lig_dir}/plk*_bs*.parm7")
                dcd = glob.glob(f"{lig_dir}/plk*_bs*.dcd")

                if parm and dcd:

                    systems.append({
                        "isoform": f"PLK{iso}",
                        "pocket": f"BS{bs}",
                        "ligand": os.path.basename(lig_dir),
                        "parm": parm[0],
                        "dcd": dcd[0]
                    })

    print("Systems discovered:", len(systems))
    return systems

# ============================================================
# POCKET SELECTION
# ============================================================

def select_pocket_atoms(u, isoform, pocket):

    residues = POCKET_RESIDUES[isoform][pocket]

    res_string = " ".join([str(r) for r in residues])

    return u.select_atoms(f"resid {res_string}")

# ============================================================
# GRID GENERATION
# ============================================================

def generate_grid(coords):

    min_xyz = coords.min(axis=0) - 4
    max_xyz = coords.max(axis=0) + 4

    xs = np.arange(min_xyz[0], max_xyz[0], GRID_SPACING)
    ys = np.arange(min_xyz[1], max_xyz[1], GRID_SPACING)
    zs = np.arange(min_xyz[2], max_xyz[2], GRID_SPACING)

    grid = np.zeros((len(xs), len(ys), len(zs)))

    return grid, xs, ys, zs

# ============================================================
# WATER DENSITY ACCUMULATION
# ============================================================

def accumulate_water_density(u, pocket_atoms):

    water = u.select_atoms("resname WAT TIP3 HOH and name O")

    grid, xs, ys, zs = generate_grid(pocket_atoms.positions)

    n_frames = len(u.trajectory)

    for ts in u.trajectory:

        pocket_coords = pocket_atoms.positions

        d = distance_array(water.positions, pocket_coords)

        close = np.any(d < POCKET_CUTOFF, axis=1)

        waters = water.positions[close]

        for w in waters:

            ix = int((w[0] - xs[0]) / GRID_SPACING)
            iy = int((w[1] - ys[0]) / GRID_SPACING)
            iz = int((w[2] - zs[0]) / GRID_SPACING)

            if (
                0 <= ix < grid.shape[0]
                and 0 <= iy < grid.shape[1]
                and 0 <= iz < grid.shape[2]
            ):
                grid[ix,iy,iz] += 1

    grid /= n_frames

    return grid, xs, ys, zs

# ============================================================
# CONVERT DENSITY → FREE ENERGY
# ============================================================

def density_to_free_energy(grid):

    rho = grid / (GRID_SPACING**3)

    with np.errstate(divide='ignore'):

        deltaG = -KT * np.log(rho / BULK_DENSITY)

    deltaG[np.isinf(deltaG)] = np.nan

    return deltaG

# ============================================================
# VISUALIZATION
# ============================================================

def plot_pocket_map(deltaG, isoform, pocket):

    slice_z = deltaG.shape[2] // 2

    plt.figure(figsize=(6,5))

    plt.imshow(
        gaussian_filter(deltaG[:,:,slice_z], sigma=1),
        cmap="coolwarm",
        origin="lower"
    )

    plt.colorbar(label="ΔG_water (kcal/mol)")

    plt.title(f"{isoform} {pocket} Water Thermodynamic Map")

    plt.tight_layout()

    plt.savefig(
        f"GIST_map_{isoform}_{pocket}.png",
        dpi=300
    )

    plt.close()

# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_system(sys):

    print("Analyzing", sys["isoform"], sys["pocket"], sys["ligand"])

    u = mda.Universe(sys["parm"], sys["dcd"])

    pocket_atoms = select_pocket_atoms(
        u,
        sys["isoform"],
        sys["pocket"]
    )

    grid, xs, ys, zs = accumulate_water_density(
        u,
        pocket_atoms
    )

    deltaG = density_to_free_energy(grid)

    plot_pocket_map(
        deltaG,
        sys["isoform"],
        sys["pocket"]
    )

    np.save(
        f"GIST_grid_{sys['isoform']}_{sys['pocket']}_{sys['ligand']}.npy",
        deltaG
    )

# ============================================================
# RUN PIPELINE
# ============================================================

systems = discover_systems()

for sys in systems:

    analyze_system(sys)

