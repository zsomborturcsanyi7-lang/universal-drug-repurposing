#!/usr/bin/env python3
"""
Docked Pose Stability & Strain Analysis Module.
Provides a computational proxy for Molecular Dynamics (MD) stability.
Calculates Ligand Strain Energy:
- Compares the docked pose geometry vs the relaxed, natural geometry.
- High strain energy (> 5.0 kcal/mol) suggests a "forced" pose that may be unstable in MD.
"""

import os
import logging
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_pose_strain(docked_pdbqt: str, smiles: str) -> float:
    """
    Calculate the internal strain energy of a docked pose.
    Returns strain in kcal/mol.
    """
    try:
        # 1. Load docked pose (convert PDBQT to Mol)
        # Note: We use the docked PDBQT as a reference for coordinates
        raw_mol = Chem.MolFromPDBBlock(open(docked_pdbqt).read())
        if not raw_mol:
            return 999.9

        # 2. Assign bond orders from SMILES (PDB/PDBQT lacks them)
        template = Chem.MolFromSmiles(smiles)
        docked_mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
        docked_mol = Chem.AddHs(docked_mol, addCoords=True)

        # 3. Calculate energy of docked pose (Force Field: UFF or MMFF)
        ff_docked = AllChem.MMFFGetMoleculeForceField(docked_mol, AllChem.MMFFGetMoleculeProperties(docked_mol))
        if not ff_docked:
            return 999.9
        e_docked = ff_docked.CalcEnergy()

        # 4. Relax the molecule (find local minimum)
        relaxed_mol = Chem.Mol(docked_mol)
        AllChem.MMFFOptimizeMolecule(relaxed_mol)
        ff_relaxed = AllChem.MMFFGetMoleculeForceField(relaxed_mol, AllChem.MMFFGetMoleculeProperties(relaxed_mol))
        e_relaxed = ff_relaxed.CalcEnergy()

        # 5. Strain = Docked Energy - Minimum Energy
        strain = e_docked - e_relaxed
        return round(strain, 2)
    except Exception as e:
        logger.debug(f"Strain calculation failed: {e}")
        return 999.9

def evaluate_stability(docking_results_csv: str, output_csv: str):
    """Analyze stability for all docked results."""
    import pandas as pd
    
    if not os.path.exists(docking_results_csv):
        logger.error(f"Results not found: {docking_results_csv}")
        return

    df = pd.read_csv(docking_results_csv)
    if 'SMILES' not in df.columns:
        logger.error("CSV missing SMILES column for strain analysis.")
        return

    strains = []
    stability_labels = []
    
    logger.info(f"Analyzing pose stability for {len(df)} compounds...")
    
    for _, row in df.iterrows():
        name = row['Ligand_Name']
        smiles = row['SMILES']
        
        # Look for the docked PDBQT file
        pdbqt_path = f"docking_results/{name}_out.pdbqt"
        if not os.path.exists(pdbqt_path):
            # Try alternative path
            pdbqt_path = f"targets/{row.get('PDB_ID', '6LU7')}/docking_results/{name}_out.pdbqt"
            
        if os.path.exists(pdbqt_path):
            strain = calculate_pose_strain(pdbqt_path, smiles)
            strains.append(strain)
            
            if strain < 3.0:
                stability_labels.append("Highly Stable")
            elif strain < 6.0:
                stability_labels.append("Stable")
            else:
                stability_labels.append("Strained / Unstable")
        else:
            strains.append(999.9)
            stability_labels.append("Missing Data")

    df['Pose_Strain_kcal_mol'] = strains
    df['MD_Stability_Proxy'] = stability_labels
    
    df.to_csv(output_csv, index=False)
    logger.info(f"Stability analysis complete. Saved to {output_csv}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        evaluate_stability(sys.argv[1], "final_validated_results.csv")
