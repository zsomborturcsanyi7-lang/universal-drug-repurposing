import joblib
import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

def run_prediction():
    print("AI Modell betöltése...")
    if not os.path.exists('vina_affinity_proxy.pkl'):
        print("Hiba: 'vina_affinity_proxy.pkl' nem található.")
        return

    model = joblib.load('vina_affinity_proxy.pkl')
    
    print("Offline könyvtár betöltése...")
    if not os.path.exists('full_fda_library.csv'):
        print("Hiba: 'full_fda_library.csv' nem található.")
        return
        
    lib = pd.read_csv('full_fda_library.csv')
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    fps = []
    valid_indices = []
    
    print("Szerkezeti ujjlenyomatok generálása...")
    for i, smiles in enumerate(lib['SMILES']):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = fp_gen.GetFingerprintAsNumPy(mol)
            fps.append(fp)
            valid_indices.append(i)
        else:
            print(f"Figyelem: Érvénytelen SMILES kihagyva: {lib.iloc[i]['Name']}")
            
    if not fps:
        print("Hiba: Egyetlen érvényes molekula sem maradt.")
        return
        
    print(f"Predikció futtatása {len(fps)} molekulára...")
    preds = model.predict(np.array(fps))
    
    # Eredmények mentése
    results = lib.iloc[valid_indices].copy()
    results['Predicted_Affinity'] = preds
    
    hits = results.sort_values('Predicted_Affinity').head(15)
    hits.to_csv('ai_hits_for_docking.csv', index=False)
    
    print("\n=== Top 15 AI Jelölt ===")
    print(hits[['Name', 'Predicted_Affinity']])
    print("\nEredmények elmentve: ai_hits_for_docking.csv")

if __name__ == "__main__":
    run_prediction()
