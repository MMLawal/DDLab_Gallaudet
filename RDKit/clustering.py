#Call your SMILES list - make it readable for RDKit analysis
from rdkit import Chem

def read_smiles_list(filename):
 """
 This function reads a list of SMILES strings from a text file.

 Args:
   filename: The path to the text file containing SMILES strings.

 Returns:
   A list of SMILES strings read from the file.
 """
 # Open the file for reading
 with open(filename, 'r') as f:
  # Read all lines from the file
  lines = f.readlines()

 # Remove leading/trailing whitespaces from each line (optional)
 smiles_list = [line.strip() for line in lines]

 return smiles_list

# Example usage:
filename = "all.csv" # Replace with your actual filename
smiles_list = read_smiles_list(filename)

# Now you can use the smiles_list for further processing
print(f"Read {len(smiles_list)} SMILES strings from the file.")


#RDKit and scikit-learn to cluster small molecules using the Tanimoto coefficient and select representative structures
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import SparseBitVect
from sklearn.cluster import KMeans

def cluster_molecules(smiles_list, num_clusters):
  """
  Clusters molecules based on Tanimoto similarity and selects representative structures.

  Args:
    smiles_list: A list of SMILES strings.
    num_clusters: The desired number of clusters.

  Returns:
    A dictionary containing cluster assignments, centroids, and representative structures.
  """

  # Convert SMILES to molecules and fingerprints
  mols = [Chem.MolFromSmiles(smile) for smile in smiles_list]
  fps = [AllChem.RDKFingerprint(mol) for mol in mols]

  # Calculate Tanimoto similarity matrix efficiently
  tanimoto_matrix = []
  for i in range(len(fps)):
    row = [DataStructs.TanimotoSimilarity(fps[i], fp) for fp in fps]
    tanimoto_matrix.append(row)

  # Cluster molecules using KMeans
  kmeans = KMeans(n_clusters=num_clusters, random_state=0)
  kmeans.fit(tanimoto_matrix)
  cluster_labels = kmeans.labels_
  centroids = kmeans.cluster_centers_ # Note: Centroids are not fingerprints here

  # Select representative structures (closest to centroid)
  representatives = []
  for i in range(num_clusters):
    cluster_mols = [mols[j] for j, label in enumerate(cluster_labels) if label == i]

    # Find molecule with minimum Euclidean distance to centroid (since centroids aren't fingerprints)
    min_dist = np.inf
    rep_index = None
    for idx, mol in enumerate(cluster_mols):
      dist = np.linalg.norm(tanimoto_matrix[idx] - centroids[i])
      if dist < min_dist:
        min_dist = dist
        rep_index = idx
    representatives.append(Chem.MolToSmiles(cluster_mols[rep_index]))

  return {
    "cluster_labels": cluster_labels,
    "centroids": centroids,
    "representatives": representatives
  }

# Example usage:
smiles_list = smiles_list # Replace with your 183 SMILES strings
num_clusters = 9 #11 or whatever size

# Cluster molecules and get results
results = cluster_molecules(smiles_list, num_clusters)

# Access cluster assignments, centroids, and representatives
cluster_labels = results["cluster_labels"]
centroids = results["centroids"]
representatives = results["representatives"]

# Print or use the results for further analysis
print("Cluster assignments:")
for i, label in enumerate(cluster_labels):
 print(f"Molecule {i}: Cluster {label}")

print("\nRepresentative structures (SMILES):")
for rep in representatives:
 print(rep)
