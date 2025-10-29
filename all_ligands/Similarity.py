import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def read_smiles_from_csv(file_path):
    """Reads SMILES strings from a CSV file."""
    try:
        df = pd.read_csv(file_path, header=None)
        # Assuming SMILES are in the first column
        return df[0].tolist()
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return []

def read_smiles_from_smi(file_path):
    """Reads SMILES strings from a .smi file."""
    try:
        with open(file_path, 'r') as f:
            # SMILES strings are typically the first part of each line, separated by a space or tab
            smiles_list = [line.split()[0] for line in f if line.strip()]
        return smiles_list
    except Exception as e:
        print(f"Error reading SMI file: {e}")
        return []

def calculate_similarity(smiles_list1, smiles_list2):
    """Calculates Tanimoto similarity between two lists of SMILES strings."""
    mols1 = [Chem.MolFromSmiles(s) for s in smiles_list1 if s]
    mols2 = [Chem.MolFromSmiles(s) for s in smiles_list2 if s]

    # Handle failed conversions
    if not all(mols1) or not all(mols2):
        print("Warning: Some SMILES strings could not be converted to molecules.")
        mols1 = [m for m in mols1 if m]
        mols2 = [m for m in mols2 if m]

    # Handle cases where one or both lists are empty after filtering
    if not mols1 or not mols2:
        return [], [], []

    # Calculate Morgan fingerprints
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols1]
    fps2 = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols2]

    similarity_matrix = []
    for fp1 in fps1:
        row = []
        for fp2 in fps2:
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            row.append(similarity)
        similarity_matrix.append(row)

    # Return the matrix and labels
    return similarity_matrix, [f'Hit_{i+1}' for i in range(len(mols1))], [f'Inhibitor_{i+1}' for i in range(len(mols2))]

def main():
    """Main function to run the analysis, calculate similarity, and save similar pairs."""
    hits_file = 'rep_struct.csv'
    inhibitors_file = 'all_20.smi'
    output_file = 'similar_pairs_summary.csv' # New output file name
    similarity_threshold = 0.25 # Define a threshold for what you consider "similar"

    hits_smiles = read_smiles_from_csv(hits_file)
    inhibitors_smiles = read_smiles_from_smi(inhibitors_file)

    if not hits_smiles or not inhibitors_smiles:
        print("Could not proceed with analysis due to file reading errors.")
        return

    # Calculate the similarity matrix and labels
    similarity_matrix, hit_labels, inhibitor_labels = calculate_similarity(hits_smiles, inhibitors_smiles)

    if not similarity_matrix:
        print("No valid molecules found for similarity calculation.")
        return

    # Create a DataFrame for better visualization and processing
    similarity_df = pd.DataFrame(similarity_matrix, index=hit_labels, columns=inhibitor_labels)

    print("Tanimoto Similarity Matrix:")
    print(similarity_df)

    # --- New Logic: Identify and Store Similar Pairs ---
    
    similar_pairs_data = []

    print(f"\nIdentifying pairs with Tanimoto similarity >= {similarity_threshold:.2f}...")

    # Iterate through the similarity matrix rows (Hits)
    for i, hit_label in enumerate(similarity_df.index):
        # Iterate through the columns (Inhibitors)
        for j, inhibitor_label in enumerate(similarity_df.columns):
            score = similarity_df.iloc[i, j]
            
            if score >= similarity_threshold:
                # Store the data for the similar pair
                hit_smiles = hits_smiles[i]
                inhibitor_smiles = inhibitors_smiles[j]
                
                similar_pairs_data.append({
                    'Hit_Label': hit_label,
                    'Inhibitor_Label': inhibitor_label,
                    'Tanimoto_Similarity': score,
                    'Hit_SMILES': hit_smiles,
                    'Inhibitor_SMILES': inhibitor_smiles
                })

    # Convert the list of dictionaries into a DataFrame
    similar_pairs_df = pd.DataFrame(similar_pairs_data)

    if not similar_pairs_df.empty:
        # Sort by similarity score in descending order
        similar_pairs_df = similar_pairs_df.sort_values(by='Tanimoto_Similarity', ascending=False)
        
        # Save the DataFrame to a CSV file
        similar_pairs_df.to_csv(output_file, index=False)
        print(f"\n✅ Successfully saved {len(similar_pairs_df)} similar pairs (Tanimoto >= {similarity_threshold:.2f}) to {output_file}")
        print("\nSummary of the most similar pairs:")
        print(similar_pairs_df.head())
    else:
        print(f"\nNo pairs found with Tanimoto similarity >= {similarity_threshold:.2f}.")

if __name__ == "__main__":
    main()
