import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

def calculate_similarity(target_smiles, library_smiles_list):
    """
    Computes Tanimoto similarity for a list of molecules against a target.
    """
    target_mol = Chem.MolFromSmiles(target_smiles)
    if not target_mol:
        raise ValueError("Invalid target SMILES")
    
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    target_fp = mfpgen.GetFingerprint(target_mol)
    
    results = []
    
    for entry in library_smiles_list:
        smiles = entry['smiles']
        mol_id = entry['id']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:
            mol_fp = mfpgen.GetFingerprint(mol)
            similarity = DataStructs.TanimotoSimilarity(target_fp, mol_fp)
            
            if similarity >= 0.40:
                results.append({
                    'ChEMBL_ID': mol_id,
                    'SMILES': smiles,
                    'Similarity': similarity
                })
                
    return pd.DataFrame(results).sort_values(by='Similarity', ascending=False)

# Example Usage
if __name__ == "__main__":
    # Target: Imatinib
    target = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
    
    # Mock library of FDA approved molecules (representative subset)
    # In a production environment, this would be loaded from a CSV or fetched from ChEMBL API
    library = [
        {'id': 'CHEMBL941', 'smiles': 'CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5'}, # Imatinib (itself)
        {'id': 'CHEMBL1201138', 'smiles': 'CC1=CN=C(C=N1)C2=C(C=C(C=C2)C)NC3=NC=CC(=N3)C4=CN=CC=C4'}, # Nilotinib
        {'id': 'CHEMBL123', 'smiles': 'CN1CCN(CC1)C2=CC=C(C=C2)NC(=O)C3=CC=CC=C3'}, # Random molecule
        {'id': 'CHEMBL192', 'smiles': 'CCN(CC)CCNC(=O)C1=CC=C(C=C1)N'} # Another random molecule
    ]
    
    df_results = calculate_similarity(target, library)
    print("Virtual Screening Results:")
    print(df_results)
