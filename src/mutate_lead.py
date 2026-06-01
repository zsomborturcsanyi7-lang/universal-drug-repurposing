import pandas as pd
import numpy as np
import os
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from prep_ligands import prepare_ligands

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def mutate_smiles(smiles, n_variants=50):
    """
    Generates variants of a lead SMILES by substituting H atoms with functional groups.
    """
    parent = Chem.MolFromSmiles(smiles)
    if not parent:
        logger.error(f"Invalid SMILES: {smiles}")
        return []
    
    # Target groups for mutation - Use valid SMILES or atom symbols
    mutations = [
        (Chem.MolFromSmiles('F'), 'F'),
        (Chem.MolFromSmiles('Cl'), 'Cl'),
        (Chem.MolFromSmiles('C'), 'CH3'),
        (Chem.MolFromSmiles('O'), 'OH'),
        (Chem.MolFromSmiles('N'), 'NH2'),
        (Chem.MolFromSmiles('CC(F)(F)F'), 'CF3') # Valid SMILES for CF3 derivative fragment
    ]
    
    # Filter out any None values if SMILES failed (shouldn't happen with these)
    mutations = [(m, n) for m, n in mutations if m is not None]
    
    # Get indices of H atoms (implicit or explicit)
    # RDKit usually works better with explicit H for this type of substitution
    parent_with_h = Chem.AddHs(parent)
    h_indices = [atom.GetIdx() for atom in parent_with_h.GetAtoms() if atom.GetSymbol() == 'H']
    
    variants = []
    seen_smiles = {Chem.MolToSmiles(parent)}
    
    # We want a variety of single and double mutations
    logger.info(f"Generating {n_variants} variants of lead scaffold...")
    
    attempt = 0
    while len(variants) < n_variants and attempt < 500:
        attempt += 1
        # Copy parent
        mol = Chem.RWMol(parent_with_h)
        
        # Decide how many mutations (1 to 2)
        num_muts = np.random.choice([1, 2], p=[0.7, 0.3])
        target_h_indices = np.random.choice(h_indices, num_muts, replace=False)
        
        valid_variant = True
        for h_idx in target_h_indices:
            # Pick a random mutation group
            mut_mol, mut_name = mutations[np.random.randint(len(mutations))]
            
            # Find the atom the H is attached to
            h_atom = mol.GetAtomWithIdx(int(h_idx))
            neighbors = h_atom.GetNeighbors()
            if not neighbors:
                valid_variant = False
                break
            
            attached_atom_idx = neighbors[0].GetIdx()
            
            # Remove the H
            mol.RemoveAtom(int(h_idx))
            
            # Add the new group
            # Note: This is a simplified substitution logic
            new_atom_idx = mol.AddAtom(mut_mol.GetAtomWithIdx(0))
            mol.AddBond(attached_atom_idx, new_atom_idx, Chem.rdchem.BondType.SINGLE)
            
        if valid_variant:
            try:
                # Cleanup and convert back to SMILES
                final_mol = mol.GetMol()
                Chem.SanitizeMol(final_mol)
                v_smiles = Chem.MolToSmiles(Chem.RemoveHs(final_mol))
                
                if v_smiles not in seen_smiles:
                    seen_smiles.add(v_smiles)
                    variants.append({
                        'Name': f"Nilotinib_Var_{len(variants)+1}",
                        'SMILES': v_smiles
                    })
            except:
                continue
                
    return variants

def run_lead_optimization():
    # Our current winner: Nilotinib
    lead_smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C"
    lead_name = "Nilotinib"
    
    # 1. Generate Variants
    variants = mutate_smiles(lead_smiles, n_variants=50)
    
    if not variants:
        logger.error("Failed to generate variants.")
        return
        
    df_variants = pd.DataFrame(variants)
    output_csv = "nilotinib_variants.csv"
    df_variants.to_csv(output_csv, index=False)
    logger.info(f"Saved {len(variants)} variants to {output_csv}")
    
    # 2. Prepare for Docking
    print(f"\nPreparing {len(variants)} Nilotinib derivatives for docking...")
    prepare_ligands(output_csv)
    
    print("\nMutation phase complete. Next step: Run 'python src/run_parallel_docking.py'")

if __name__ == "__main__":
    run_lead_optimization()
