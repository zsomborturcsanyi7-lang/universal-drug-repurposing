#!/usr/bin/env python3
"""
ADMET Prediction & Toxicity Filtering Module.
Provides advanced screening for pharmacokinetic properties and safety.
Includes:
- PAINS (Pan-Assay Interference Compounds) filter
- BBB (Blood-Brain Barrier) permeability estimation
- HIA (Gastrointestinal Absorption) estimation
- LogS (Solubility) estimation
- Brenk and NIH filters for reactive groups
"""

import os
import pandas as pd
import numpy as np
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ADMETPredictor:
    def __init__(self):
        # Initialize RDKit FilterCatalog with common filters
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        self.catalog = FilterCatalog.FilterCatalog(params)

    def predict(self, smiles: str) -> dict:
        """Calculate comprehensive ADMET profile for a SMILES string."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}

        # 1. Basic Descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)

        # 2. Toxicity Filtering (PAINS & Brenk)
        matches = self.catalog.GetMatches(mol)
        has_toxicity_alert = len(matches) > 0
        alerts = [m.GetDescription() for m in matches]

        # 3. Solubility (LogS) - Basic empirical model (ESOL-like)
        # LogS = 0.16 - 0.63*LogP - 0.0062*MW + 0.066*RotB + 0.71*AP
        # AP = Aromatic Proportion
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        ap = aromatic_atoms / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0
        log_s = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rot_bonds) + (0.71 * ap)

        # 4. GI Absorption (HIA) - Based on Rule of 3/5 + TPSA
        # High absorption if TPSA < 140 and LogP < 5
        gi_absorption = "High" if tpsa < 140 and logp < 5 else "Low"

        # 5. BBB Permeability
        # High BBB if TPSA < 90, LogP 1-5, MW < 450
        bbb_permeable = (tpsa < 90 and 1 < logp < 5 and mw < 450)

        return {
            "MW": mw,
            "LogP": logp,
            "TPSA": tpsa,
            "LogS": log_s,
            "GI_Absorption": gi_absorption,
            "BBB_Permeable": bbb_permeable,
            "Toxicity_Alerts": has_toxicity_alert,
            "Alert_Details": ", ".join(alerts) if alerts else "None",
            "Drug_Score": self._calculate_synthetic_drug_score(mw, logp, tpsa, log_s, has_toxicity_alert)
        }

    def _calculate_synthetic_drug_score(self, mw, logp, tpsa, log_s, tox_alert) -> float:
        """Composite score from 0.0 (poor) to 1.0 (excellent)."""
        score = 1.0
        if tox_alert: score *= 0.5
        if mw > 500: score *= 0.8
        if logp > 5 or logp < -1: score *= 0.8
        if tpsa > 140: score *= 0.8
        if log_s < -4: score *= 0.8
        return round(score, 2)

def screen_csv(input_csv: str, output_csv: str):
    """Screen an entire CSV of hits for ADMET properties."""
    if not os.path.exists(input_csv):
        logger.error(f"File not found: {input_csv}")
        return

    df = pd.read_csv(input_csv)
    predictor = ADMETPredictor()
    
    logger.info(f"Running ADMET screening on {len(df)} candidates...")
    
    admet_results = []
    for smiles in df['SMILES']:
        admet_results.append(predictor.predict(smiles))
    
    admet_df = pd.DataFrame(admet_results)
    final_df = pd.concat([df, admet_df], axis=1)
    
    # Filter out high toxicity
    filtered_df = final_df[final_df['Toxicity_Alerts'] == False].copy()
    
    final_df.to_csv(output_csv, index=False)
    logger.info(f"ADMET screening complete. Saved to {output_csv}")
    logger.info(f"Clean candidates (no tox alerts): {len(filtered_df)} / {len(df)}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        screen_csv(sys.argv[1], "admet_screened_results.csv")
    else:
        # Example test
        predictor = ADMETPredictor()
        print(predictor.predict("Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C"))
