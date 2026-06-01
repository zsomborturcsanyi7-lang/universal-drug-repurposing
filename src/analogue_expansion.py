import pandas as pd
import numpy as np
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import logging
import os
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from prep_ligands import prepare_ligands

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def find_analogues_offline(hits_csv, library_csv, n_analogues=200, similarity_threshold=0.7):
    """
    Finds analogues in the local library based on structural similarity to current hits.
    Used when ChEMBL API is unstable or for local refinement.
    """
    if not os.path.exists(hits_csv) or not os.path.exists(library_csv):
        logger.error("Required CSV files (hits or library) missing.")
        return pd.DataFrame()

    hits_df = pd.read_csv(hits_csv)
    lib_df = pd.read_csv(library_csv)
    
    logger.info(f"Searching for analogues of {len(hits_df)} hits in library of {len(lib_df)} molecules...")

    # Generate fingerprints for hits
    hit_mols = [Chem.MolFromSmiles(s) for s in hits_df['SMILES'].tolist()]
    hit_mols = [m for m in hit_mols if m]
    hit_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in hit_mols]

    # Generate fingerprints for library
    lib_mols_data = []
    for _, row in lib_df.iterrows():
        m = Chem.MolFromSmiles(row['SMILES'])
        if m:
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
            lib_mols_data.append({'Name': row['Name'], 'SMILES': row['SMILES'], 'FP': fp})

    analogues = []
    seen_names = set(hits_df['Name'].tolist())

    for lib_item in lib_mols_data:
        if lib_item['Name'] in seen_names:
            continue
            
        # Check similarity against all hits
        max_sim = 0
        for h_fp in hit_fps:
            sim = DataStructs.TanimotoSimilarity(lib_item['FP'], h_fp)
            if sim > max_sim:
                max_sim = sim
        
        if max_sim >= similarity_threshold:
            lib_item['Max_Similarity'] = max_sim
            analogues.append(lib_item)

    if not analogues:
        logger.warning("No analogues found meeting the threshold.")
        return pd.DataFrame()

    analogues_df = pd.DataFrame(analogues).sort_values(by='Max_Similarity', ascending=False).head(n_analogues)
    
    # Remove FP column for saving
    analogues_df = analogues_df.drop(columns=['FP'])
        
    return analogues_df

def run_expansion():
    hits_file = 'ai_hits_for_docking.csv'
    lib_file = 'full_fda_library.csv'
    output_file = 'analogue_expansion_set.csv'

    # 1. Find analogues - Lower threshold to 0.3 to ensure we get data from a small library
    analogues_df = find_analogues_offline(hits_file, lib_file, n_analogues=200, similarity_threshold=0.3)
    
    if analogues_df.empty:
        logger.warning("No analogues found with current settings.")
        return

    logger.info(f"Found {len(analogues_df)} analogues. Saving to {output_file}")
    analogues_df.to_csv(output_file, index=False)

    # 2. Prepare for docking
    print(f"\nPreparing {len(analogues_df)} analogues for docking...")
    prepare_ligands(output_file)
    
    print("\nAnalogue expansion complete. Run 'python run_parallel_docking.py' to generate 'Expert' training data.")

if __name__ == "__main__":
    run_expansion()
