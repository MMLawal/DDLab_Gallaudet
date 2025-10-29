# Requires RDKit. Example: conda install -c conda-forge rdkit
from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd
import numpy as np
import os

# Paths (change if needed)
hits_sdf = 'rep_struct.csv'          # your uploaded hits
known_sdf = 'all_struct.smi'  # provide/upload this file

def read_mols(path):
    if path.lower().endswith('.sdf'):
        suppl = Chem.SDMolSupplier(path)
        mols = [m for m in suppl if m is not None]
    elif path.lower().endswith('.smi') or path.lower().endswith('.csv'):
        mols = []
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    smi = parts[0]
                    m = Chem.MolFromSmiles(smi)
                    if m:
                        mols.append(m)
    else:
        raise ValueError("Unsupported file type; use .sdf or .smi")
    return mols

def canonical_smiles(m):
    try:
        return Chem.MolToSmiles(m, isomericSmiles=True)
    except:
        return None

def make_fp(m, radius=2, nBits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nBits)

# Load molecules
hits = read_mols(hits_sdf)
knowns = read_mols(known_sdf)

# Basic checks
print(f"Loaded {len(hits)} hits and {len(knowns)} known inhibitors")

# Generate canonical SMILES and fingerprints
hits_data = []
for i,m in enumerate(hits):
    smi = canonical_smiles(m)
    fp = make_fp(m)
    hits_data.append({'idx': i, 'smiles': smi, 'mol': m, 'fp': fp})

known_data = []
for j,m in enumerate(knowns):
    smi = canonical_smiles(m)
    fp = make_fp(m)
    # try to get name from SDF properties
    name = m.GetProp('_Name') if m.HasProp('_Name') else f'known_{j}'
    known_data.append({'idx': j, 'name': name, 'smiles': smi, 'mol': m, 'fp': fp})

# Compute pairwise similarities, for each hit keep top3 known matches
rows = []
for h in hits_data:
    scores = []
    for k in known_data:
        s = TanimotoSimilarity(h['fp'], k['fp'])
        scores.append((k['idx'], k['name'], k['smiles'], s))
    scores_sorted = sorted(scores, key=lambda x: x[3], reverse=True)
    top3 = scores_sorted[:3]
    row = {
        'hit_idx': h['idx'],
        'hit_smiles': h['smiles'],
        'top1_known_idx': top3[0][0],
        'top1_known_name': top3[0][1],
        'top1_known_smiles': top3[0][2],
        'top1_score': top3[0][3],
        'top2_known_idx': top3[1][0],
        'top2_known_name': top3[1][1],
        'top2_known_smiles': top3[1][2],
        'top2_score': top3[1][3],
        'top3_known_idx': top3[2][0],
        'top3_known_name': top3[2][1],
        'top3_known_smiles': top3[2][2],
        'top3_score': top3[2][3],
    }
    rows.append(row)

df = pd.DataFrame(rows)
out_csv = 'known_similarity.csv'
df.to_csv(out_csv, index=False)
print("Wrote similarity table to", out_csv)

# Summary stats
all_top_scores = df['top1_score']
print("Top-1 similarity: mean {:.3f}, median {:.3f}, max {:.3f}".format(all_top_scores.mean(), all_top_scores.median(), all_top_scores.max()))

