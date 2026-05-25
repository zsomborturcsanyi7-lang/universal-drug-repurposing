import logging
import pandas as pd
import random
from rdkit import Chem
from chembl_webresource_client.new_client import new_client
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def fetch_diverse_training_set(n_compounds=50):
    """
    Fetches a diverse set of FDA-approved compounds for training data collection.
    """
    logger.info("Connecting to ChEMBL API...")
    molecule = new_client.molecule
    
    logger.info(f"Fetching FDA-approved small-molecule drugs...")
    # Fetching compounds with max_phase=4 (FDA Approved)
    res = molecule.filter(max_phase=4).only('molecule_chembl_id', 'pref_name', 'molecule_structures')
    
    all_drugs = []
    for mol in res:
        if mol['molecule_structures'] and mol['molecule_structures'].get('canonical_smiles'):
            all_drugs.append({
                'Name': mol['pref_name'] if mol['pref_name'] else mol['molecule_chembl_id'],
                'SMILES': mol['molecule_structures'].get('canonical_smiles'),
                'ChEMBL_ID': mol['molecule_chembl_id']
            })
    
    logger.info(f"Retrieved {len(all_drugs)} approved drugs. Selecting {n_compounds} random compounds for diversity...")
    
    # Random sample for diverse training set
    if len(all_drugs) > n_compounds:
        training_set = random.sample(all_drugs, n_compounds)
    else:
        training_set = all_drugs
        
    return pd.DataFrame(training_set)

if __name__ == "__main__":
    output_csv = "training_set_50.csv"
    
    try:
        # Step 1: Fetch 50 diverse drugs
        df_training = fetch_diverse_training_set(50)
        df_training.to_csv(output_csv, index=False)
        logger.info(f"Saved 50 diverse compounds to {output_csv}")
        
    except Exception as e:
        logger.error(f"Failed to fetch training set: {e}")
