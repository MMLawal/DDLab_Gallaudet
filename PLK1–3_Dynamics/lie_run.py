#run within AmberTools25 environment

import os
import pandas as pd
from collections import defaultdict
import glob

ISOFORM_RESIDUE_MAP = {
    "PLK1": {"protein": "1-224", "ligand": "225"},
    "PLK2": {"protein": "1-213", "ligand": "214"},
    "PLK3": {"protein": "1-221", "ligand": "222"},
}
def load_system(prmtop, dcd, ligand_resname):
    u = mda.Universe(prmtop, dcd)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms("not protein")
    return u, protein, ligand

def discover_systems(base_dir="."):
    systems = []

    # Example: ./plk1-bs1/3v/
    for iso_pocket in glob.glob(os.path.join(base_dir, "plk*-bs*")):
        isoform, pocket = os.path.basename(iso_pocket).split("-")

        for lig_dir in glob.glob(os.path.join(iso_pocket, "*")):
            ligand = os.path.basename(lig_dir)

            prmtop_files = glob.glob(os.path.join(lig_dir, "plk*_bs*.parm7"))
            dcd_files = glob.glob(os.path.join(lig_dir, "plk*_bs*.dcd"))

            if not prmtop_files or not dcd_files:
                continue  # skip incomplete systems

            systems.append({
                "prmtop": prmtop_files[0],
                "dcd": dcd_files[0],
                "ligand": ligand,
                "isoform": isoform.upper(),   # PLK1 / PLK2 / PLK3
                "pocket": pocket.upper()      # BS1–BS6
            })

    return systems

systems = discover_systems(base_dir=".")

print(f"Discovered {len(systems)} systems")
print(systems[0])

import subprocess
import numpy as np

def compute_interaction_energy(sys, workdir="lie_output"):

    os.makedirs(workdir, exist_ok=True)

    iso = sys["isoform"]
    prmtop = sys["prmtop"]
    dcd = sys["dcd"]

    res_map = ISOFORM_RESIDUE_MAP[iso]
    protein_sel = res_map["protein"]
    ligand_sel = res_map["ligand"]

    # Unique output name
    out_prefix = f"{iso}_{sys['pocket']}_{sys['ligand']}"
    lie_out = os.path.join(workdir, f"{out_prefix}_lie.dat")
    lie_in = os.path.join(workdir, f"{out_prefix}_lie.in")

    # Build cpptraj input
    with open(lie_in, "w") as f:
        f.write(f"parm {prmtop}\n")
        f.write(f"trajin {dcd}\n")
        f.write(
            f"lie :{ligand_sel} :{protein_sel} "
            f"out {lie_out} cutvdw 12.0 cutelec 12.0\n"
        )
        f.write("run\n")

    # Run cpptraj
    try:
        subprocess.run(
            ["cpptraj", "-i", lie_in],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        print(f"cpptraj failed for {out_prefix}")
        print(e.stderr.decode())
        return None

    # Parse lie.dat
    if not os.path.exists(lie_out):
        print(f"Missing LIE output for {out_prefix}")
        return None

for sys in systems:
    print(f"Running LIE for {sys['isoform']} {sys['pocket']} {sys['ligand']}")
    compute_interaction_energy(sys)
