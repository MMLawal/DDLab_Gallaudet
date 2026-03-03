#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pocket Volume Fluctuation Entropy (PVFE) Analysis
For PLK1/2/3 PBD-ligand complexes (18 μs aggregate simulation time)

"""

import os
import glob
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances
from scipy.spatial import ConvexHull
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================================================
# POCKET RESIDUE DEFINITIONS (Explicit residue lists per isoform-pocket)
# =============================================================================
# Residue IDs adjusted: PLK1 (-370), PLK2 (-468), PLK3 (-425)
# These represent the POLO-BOX DOMAIN only (not full protein)

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

# Analysis parameters
LIGAND_RESNAMES = ["UNK", "UNL", "VIH", "LIG"]
VOLUME_BINS = 25     # For entropy calculation
MIN_POCKET_ATOMS = 5  # Minimum atoms for stable convex hull (reduced for explicit residues)


# =============================================================================
# SYSTEM LOADING
# =============================================================================

def load_system(prmtop, dcd):
    """Load stripped trajectory with robust ligand selection"""
    u = mda.Universe(prmtop, dcd, in_memory=False)
    protein = u.select_atoms("protein")
    # Select ligand by resname (water/ions already stripped per user spec)
    ligand = u.select_atoms(f"(not protein) and resname {' '.join(LIGAND_RESNAMES)}")
    
    if len(ligand) == 0:
        raise ValueError(f"No ligand found with resnames {LIGAND_RESNAMES} in {prmtop}")
    return u, protein, ligand


# =============================================================================
# POCKET VOLUME CALCULATION (REVISED: Explicit residue-based)
# =============================================================================

def get_pocket_atom_indices(protein, isoform, pocket):
    """
    Get atom indices for explicit pocket residues
    
    Returns:
        pocket_atom_indices: list of atom indices belonging to pocket residues
        n_pocket_residues: number of unique residues in pocket
    """
    pocket_residues = POCKET_RESIDUES.get(isoform, {}).get(pocket, [])
    
    if not pocket_residues:
        return [], 0
    
    # Find atoms belonging to pocket residues
    pocket_atom_indices = []
    pocket_resids_found = set()
    
    for atom in protein:
        if atom.resid in pocket_residues:
            pocket_atom_indices.append(atom.index)
            pocket_resids_found.add(atom.resid)
    
    return pocket_atom_indices, len(pocket_resids_found)


def compute_frame_volume_explicit(protein, isoform, pocket, box):
    """
    Compute pocket volume via convex hull of EXPLICIT pocket residue atoms
    (not distance-based selection)
    
    Returns volume in Å³ or np.nan if pocket undefined
    """
    # Get atom indices for explicit pocket residues
    pocket_atom_indices, n_residues = get_pocket_atom_indices(protein, isoform, pocket)
    
    # Safety check: need sufficient residues/atoms for meaningful hull
    if n_residues < 3 or len(pocket_atom_indices) < MIN_POCKET_ATOMS:
        return np.nan
    
    # Extract positions of pocket atoms
    pocket_positions = protein.positions[pocket_atom_indices]
    
    # Compute convex hull volume
    try:
        # Use QJ option to handle near-coplanar points robustly
        hull = ConvexHull(pocket_positions, qhull_options="QJ Pp")
        return hull.volume
    except Exception:
        return np.nan


def calculate_pvfe(system):
    """
    Compute Pocket Volume Fluctuation Entropy for a single system
    Using EXPLICIT pocket residue definitions (not ligand-radius)
    
    Returns dict with entropy and diagnostic metrics
    """
    try:
        u, protein, ligand = load_system(system["prmtop"], system["dcd"])
        isoform = system["isoform"]
        pocket = system["pocket"]
        
        volumes = []
        n_frames_processed = 0
        
        # Frame-by-frame pocket volume calculation using explicit residues
        for ts in u.trajectory:
            vol = compute_frame_volume_explicit(protein, isoform, pocket, ts.dimensions)
            volumes.append(vol)
            n_frames_processed += 1
        
        volumes = np.array(volumes)
        valid_volumes = volumes[~np.isnan(volumes)]
        
        # Entropy calculation from volume distribution
        if len(valid_volumes) < 0.8 * len(volumes):  # >20% frames failed
            entropy = np.nan
            diagnostic = "HIGH_NAN_RATE"
        elif len(valid_volumes) < 50:  # Insufficient sampling
            entropy = np.nan
            diagnostic = "INSUFFICIENT_SAMPLES"
        else:
            # Shannon entropy with natural log (units: nats)
            hist, _ = np.histogram(valid_volumes, bins=VOLUME_BINS, density=False)
            probs = hist[hist > 0] / hist.sum()
            entropy = -np.sum(probs * np.log(probs))
            diagnostic = "SUCCESS"
        
        # Additional diagnostics for J. Phys. Chem. B rigor
        return {
            "isoform": isoform,
            "pocket": pocket,
            "ligand_id": system["ligand"],
            "pvfe_entropy": entropy,
            "mean_volume": np.mean(valid_volumes) if len(valid_volumes) > 0 else np.nan,
            "volume_std": np.std(valid_volumes) if len(valid_volumes) > 0 else np.nan,
            "frames_total": n_frames_processed,
            "frames_valid": len(valid_volumes),
            "n_pocket_residues": len(POCKET_RESIDUES.get(isoform, {}).get(pocket, [])),
            "diagnostic": diagnostic
        }
    
    except Exception as e:
        return {
            "isoform": system["isoform"],
            "pocket": system["pocket"],
            "ligand_id": system["ligand"],
            "pvfe_entropy": np.nan,
            "mean_volume": np.nan,
            "volume_std": np.nan,
            "frames_total": 0,
            "frames_valid": 0,
            "n_pocket_residues": 0,
            "diagnostic": f"ERROR:{str(e)[:50]}"
        }


# =============================================================================
# SYSTEM DISCOVERY (Unchanged)
# =============================================================================

def discover_systems(base_dir="."):
    """Discover all PLK-ligand simulation systems with validation"""
    systems = []
    isoform_dirs = glob.glob(os.path.join(base_dir, "plk[1-3]-bs[1-6]"))
    
    for iso_pocket_dir in isoform_dirs:
        isoform_pocket = os.path.basename(iso_pocket_dir)
        isoform, pocket = isoform_pocket.split("-")
        
        # Validate pocket naming
        if pocket.upper() not in ["BS1", "BS2", "BS3", "BS4", "BS5", "BS6"]:
            continue
            
        ligand_dirs = glob.glob(os.path.join(iso_pocket_dir, "*"))
        for lig_dir in ligand_dirs:
            ligand_id = os.path.basename(lig_dir)
            
            # Find topology and trajectory files (Amber format)
            prmtop = glob.glob(os.path.join(lig_dir, "plk*_bs*.parm7"))
            dcd = glob.glob(os.path.join(lig_dir, "plk*_bs*.dcd"))
            
            if not prmtop or not dcd:
                continue
                
            systems.append({
                "prmtop": prmtop[0],
                "dcd": dcd[0],
                "ligand": ligand_id,
                "isoform": f"PLK{isoform[3:]}",  # Convert 'plk1' -> 'PLK1'
                "pocket": pocket.upper()
            })
    
    print(f"Discovered {len(systems)} simulation systems")
    return systems


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main(base_dir=".", n_processes=None):
    """Orchestrate PVFE analysis across all systems"""
    if n_processes is None:
        n_processes = min(cpu_count(), 32)  # Cap for memory safety
    
    systems = discover_systems(base_dir)
    
    if not systems:
        raise RuntimeError("No simulation systems found. Check directory structure.")
    
    # Parallel processing with progress tracking
    print(f"Computing PVFE (explicit residues) for {len(systems)} systems using {n_processes} processes...")
    with Pool(processes=n_processes) as pool:
        results = pool.map(calculate_pvfe, systems)
    
    # Compile results and save
    df = pd.DataFrame(results)
    output_file = "pvfe_explicit_residues_results.csv"
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    
    # Summary statistics for J. Phys. Chem. B reporting
    print("\nPVFE Summary Statistics (explicit residues, nats):")
    print(df.groupby(["isoform", "pocket"])["pvfe_entropy"].describe().round(3))
    
    # Flag systems requiring manual inspection
    failed = df[df["diagnostic"] != "SUCCESS"]
    if not failed.empty:
        print(f"\n⚠️  {len(failed)} systems require inspection (see 'diagnostic' column)")
        print(failed[["isoform", "pocket", "ligand_id", "diagnostic"]].head(10))
    
    return df


if __name__ == "__main__":
    # Example usage: python pvfe_explicit.py /path/to/simulations
    import sys
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    results_df = main(base_dir)
