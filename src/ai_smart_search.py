import pandas as pd
import numpy as np
import joblib
import os
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import logging
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, Descriptors
from prep_ligands import prepare_ligands

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_physchem_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return [0] * 5
    return [Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol), Descriptors.TPSA(mol)]

def generate_boosted_features(smiles_list):
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    features = []
    valid_indices = []
    for i, smiles in enumerate(smiles_list):
        if not smiles or not isinstance(smiles, str): continue
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = fp_gen.GetFingerprintAsNumPy(mol)
            phys = get_physchem_descriptors(smiles)
            features.append(np.concatenate([fp, phys]))
            valid_indices.append(i)
    return np.array(features), valid_indices

def run_expert_ai_search(n_candidates=50):
    model_path = 'vina_affinity_proxy.pkl'
    if not os.path.exists(model_path):
        logger.error("Expert model not found. Run 'python train_boosted_proxy.py' first.")
        return

    logger.info("Loading Expert AI Model (HistGradientBoosting)...")
    model = joblib.load(model_path)
    
    # Load candidate pool
    lib_file = 'full_fda_library.csv'
    if not os.path.exists(lib_file):
        logger.error(f"{lib_file} missing.")
        return
        
    candidates_df = pd.read_csv(lib_file)
    logger.info(f"Screening {len(candidates_df)} drugs with Expert AI...")
    
    # Predict
    X, valid_idx = generate_boosted_features(candidates_df['SMILES'].tolist())
    predictions = model.predict(X)
    
    results = candidates_df.iloc[valid_idx].copy()
    results['Predicted_Affinity'] = predictions
    
    # Get top 50
    hits = results.sort_values(by='Predicted_Affinity').head(n_candidates)
    
    print("\n" + "*"*50)
    print(f"EXPERT AI TOP {n_candidates} CANDIDATES")
    print("*"*50)
    print(hits[['Name', 'Predicted_Affinity']].to_string(index=False))
    print("*"*50)
    
    # Save and Prepare
    output_file = 'expert_ai_hits.csv'
    hits.to_csv(output_file, index=False)
    logger.info(f"Saved expert hits to {output_file}")
    
    print(f"\nTriggering preparation for {len(hits)} high-probability ligands...")
    prepare_ligands(output_file)
    
    print("\nExpert screening complete. Next: Run 'python run_parallel_docking.py'")

if __name__ == "__main__":
    run_expert_ai_search(50)
