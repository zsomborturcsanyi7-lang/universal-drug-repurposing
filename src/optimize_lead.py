#!/usr/bin/env python3
"""
Lead Optimizer — Generate structural variants of top docked drugs.
Takes best docking hits and creates variants to improve binding affinity.

Strategies:
1. Bioisostere replacement (ring variants, linker swaps)
2. Methyl/fluoro scans
3. Fragment growing/shrinking
4. Random mutations via SMILES manipulation

Usage:
    python optimize_lead.py --smiles "CC(=O)..." --name Nilotinib
    python optimize_lead.py --csv expert_ai_hits.csv --top 3
    python optimize_lead.py --pdb 6LU7  # Auto-reads top docked hits
"""

import os
import sys
import random
import logging
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import Recap, BRICS

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# Variant generation
# ──────────────────────────────────────────────

def generate_methyl_scan(smiles: str, n_variants: int = 20) -> list:
    """
    Add/remove methyl groups at various positions.
    Simple approach: replace H with CH3 at random positions.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    variants = []
    
    # Add methyl at each possible position
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() < 4:  # Carbon with free bonds
            new_mol = Chem.RWMol(mol)
            idx = new_mol.AddAtom(Chem.Atom(6))  # Carbon
            new_mol.AddBond(atom.GetIdx(), idx, Chem.BondType.SINGLE)
            
            try:
                new_mol = new_mol.GetMol()
                Chem.SanitizeMol(new_mol)
                variants.append(Chem.MolToSmiles(new_mol))
            except:
                pass
    
    # De-duplicate and limit
    variants = list(set(variants))[:n_variants]
    
    return variants


def generate_bioisostere_replacements(smiles: str) -> list:
    """
    Replace common functional groups with bioisosteres.
    """
    replacements = [
        # Carboxylic acid isosteres
        ("C(=O)O", "C(=O)NS(=O)(=O)"),  # → sulfonamide
        ("C(=O)O", "C(=O)N"),            # → amide
        ("C(=O)O", "c1nn[nH]n1"),       # → tetrazole
        ("C(=O)O", "P(=O)(O)O"),        # → phosphonate
        
        # Amide isosteres
        ("C(=O)N", "C(=S)N"),           # → thioamide
        ("C(=O)N", "C(=O)O"),           # → acid
        ("C(=O)N", "c1noc1"),           # → oxadiazole
        
        # Phenyl isosteres
        ("c1ccccc1", "c1ccccn1"),       # → pyridine
        ("c1ccccc1", "c1ccncc1"),       # → pyrazine
        ("c1ccccc1", "c1cnccn1"),       # → pyrimidine
        ("c1ccccc1", "C1CCCCC1"),       # → cyclohexane
        
        # Ester isosteres
        ("C(=O)OC", "C(=O)NC"),         # → amide
        ("C(=O)OC", "C(=O)C"),          # → ketone
        
        # Methyl isosteres
        ("C", "Cl"),                    # Methyl → Chloro
        ("C", "F"),                     # Methyl → Fluoro
        ("C", "Br"),                    # Methyl → Bromo
    ]
    
    variants = []
    for old, new in replacements:
        if old in smiles:
            variants.append(smiles.replace(old, new, 1))
    
    return variants


def generate_fragment_variants(smiles: str, n_variants: int = 20) -> list:
    """
    Use BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures)
    to fragment the molecule and recombine.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    variants = []
    
    try:
        # BRICS decomposition
        fragments = BRICS.BRICSDecompose(mol, returnMols=True)
        frag_smiles = [Chem.MolToSmiles(f) for f in fragments if f]
        
        if len(frag_smiles) >= 2:
            # Recombine different fragments
            for i in range(min(len(frag_smiles), 5)):
                for j in range(i + 1, min(len(frag_smiles), 5)):
                    # Simple concatenation (BRICS BUILD would be better)
                    combined = frag_smiles[i] + "." + frag_smiles[j]
                    variants.append(combined)
    except:
        pass
    
    return variants[:n_variants]


def generate_random_structure_mutations(smiles: str, n_variants: int = 30) -> list:
    """
    Generate random structural mutations:
    - Ring expansion/contraction
    - Atom type changes (C→N, O→S, etc.)
    - Bond order changes
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    variants = []
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    
    for _ in range(n_variants):
        try:
            new_mol = Chem.RWMol(mol)
            mutation_type = random.choice(['atom', 'bond', 'add', 'remove'])
            
            if mutation_type == 'atom' and atoms:
                atom = random.choice(atoms)
                # Only mutate carbon atoms
                if atom.GetAtomicNum() == 6:
                    new_atom = random.choice([7, 8, 16])  # N, O, S
                    new_mol.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(new_atom)
            
            elif mutation_type == 'bond' and bonds:
                bond = random.choice(bonds)
                # Toggle between single and double
                current = bond.GetBondType()
                new_type = Chem.BondType.DOUBLE if current == Chem.BondType.SINGLE else Chem.BondType.SINGLE
                new_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                new_mol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), new_type)
            
            new_mol = new_mol.GetMol()
            Chem.SanitizeMol(new_mol)
            variants.append(Chem.MolToSmiles(new_mol))
        except:
            continue
    
    return list(set(variants))[:n_variants]


def generate_fluorine_scan(smiles: str) -> list:
    """
    Replace H with F at various positions (common med-chem optimization).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    variants = []
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            # Check if this carbon has hydrogen neighbors
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 1:  # Hydrogen
                    try:
                        new_mol = Chem.RWMol(mol)
                        new_mol.GetAtomWithIdx(neighbor.GetIdx()).SetAtomicNum(9)  # H → F
                        new_mol = new_mol.GetMol()
                        Chem.SanitizeMol(new_mol)
                        variants.append(Chem.MolToSmiles(new_mol))
                    except:
                        continue
    
    return list(set(variants))


# ──────────────────────────────────────────────
# Multi-target optimization
# ──────────────────────────────────────────────

def generate_variants(smiles: str, name: str, n_total: int = 50) -> pd.DataFrame:
    """
    Generate a diverse set of structural variants.
    
    Args:
        smiles: SMILES string of the lead molecule
        name: Drug name
        n_total: Target number of variants (~50 total)
    
    Returns:
        DataFrame with Name, SMILES, Variant_Type columns
    """
    logger.info(f"Generating variants for {name}...")
    
    all_variants = []
    variant_types = []
    
    # Strategy 1: Methyl scan (10 variants)
    methyl_variants = generate_methyl_scan(smiles, 10)
    all_variants.extend(methyl_variants)
    variant_types.extend(['methyl_scan'] * len(methyl_variants))
    logger.info(f"  Methyl scan: {len(methyl_variants)} variants")
    
    # Strategy 2: Bioisostere replacements (as many as found)
    bio_variants = generate_bioisostere_replacements(smiles)
    all_variants.extend(bio_variants)
    variant_types.extend(['bioisostere'] * len(bio_variants))
    logger.info(f"  Bioisosteres: {len(bio_variants)} variants")
    
    # Strategy 3: Fluoro scan (as many as found)
    fluoro_variants = generate_fluorine_scan(smiles)
    all_variants.extend(fluoro_variants)
    variant_types.extend(['fluoro_scan'] * len(fluoro_variants))
    logger.info(f"  Fluoro scan: {len(fluoro_variants)} variants")
    
    # Strategy 4: Random mutations (fill remaining)
    remaining = max(0, n_total - len(all_variants))
    if remaining > 0:
        random_variants = generate_random_structure_mutations(smiles, remaining)
        all_variants.extend(random_variants)
        variant_types.extend(['random_mutation'] * len(random_variants))
        logger.info(f"  Random mutations: {len(random_variants)} variants")
    
    # Strategy 5: Additional fragment-based variants
    remaining2 = max(0, n_total - len(all_variants))
    if remaining2 > 0:
        frag_variants = generate_fragment_variants(smiles, remaining2)
        all_variants.extend(frag_variants)
        variant_types.extend(['fragment'] * len(frag_variants))
    
    # Deduplicate
    seen = set()
    unique_variants = []
    unique_types = []
    for v, t in zip(all_variants, variant_types):
        if v not in seen and v != smiles:
            seen.add(v)
            unique_variants.append(v)
            unique_types.append(t)
    
    # Build DataFrame
    df = pd.DataFrame({
        'Name': [f"{name}_Var_{i+1}" for i in range(len(unique_variants))],
        'SMILES': unique_variants,
        'Variant_Type': unique_types,
        'Parent_Drug': [name] * len(unique_variants),
    })
    
    logger.info(f"  Total unique variants: {len(df)}")
    
    return df


# ──────────────────────────────────────────────
# Optimize from docking results
# ──────────────────────────────────────────────

def optimize_from_results(results_csv: str = "final_3d_docking_insights.csv",
                           top_n: int = 3, 
                           variants_per_lead: int = 50) -> pd.DataFrame:
    """
    Read docking results and generate variants for the top hits.
    """
    if not os.path.exists(results_csv):
        logger.error(f"No results file: {results_csv}")
        return pd.DataFrame()
    
    df = pd.read_csv(results_csv)
    
    # Take top hits
    top_hits = df.head(top_n)
    logger.info(f"Optimizing top {top_n} docking hits:")
    
    all_variants = []
    
    for _, row in top_hits.iterrows():
        name = row.iloc[0]
        # Try to find SMILES
        smiles = None
        
        # Check expert_ai_hits or full_fda_library
        for lib in ["expert_ai_hits.csv", "full_fda_library.csv", "ai_hits_for_docking.csv"]:
            if os.path.exists(lib):
                lib_df = pd.read_csv(lib)
                match = lib_df[lib_df['Name'].str.contains(name.split('_')[0], case=False, na=False)]
                if not match.empty and 'SMILES' in match.columns:
                    smiles = match.iloc[0]['SMILES']
                    break
        
        if not smiles:
            logger.warning(f"  No SMILES found for {name}")
            continue
        
        logger.info(f"  {name}: generating variants...")
        variants = generate_variants(smiles, name, variants_per_lead)
        all_variants.append(variants)
    
    if not all_variants:
        logger.warning("No variants generated")
        return pd.DataFrame()
    
    result = pd.concat(all_variants, ignore_index=True)
    result.to_csv("optimized_variants.csv", index=False)
    logger.info(f"\nSaved {len(result)} variants to optimized_variants.csv")
    
    return result


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lead Optimizer")
    parser.add_argument('--smiles', help='SMILES string of lead molecule')
    parser.add_argument('--name', default='Lead', help='Molecule name')
    parser.add_argument('--csv', help='Docking results CSV')
    parser.add_argument('--top', type=int, default=3, help='Top N hits to optimize')
    parser.add_argument('--variants', type=int, default=50, help='Variants per lead')
    parser.add_argument('--pdb', help='PDB ID (auto-reads results)')
    parser.add_argument('--output', default='optimized_variants.csv')
    
    args = parser.parse_args()
    
    if args.smiles:
        # Single molecule mode
        variants = generate_variants(args.smiles, args.name, args.variants)
        variants.to_csv(args.output, index=False)
        print(f"\nGenerated {len(variants)} variants → {args.output}")
        
        # Show top 5
        print("\n  Top 5 variants:")
        for i, row in variants.head(5).iterrows():
            print(f"  {row['Name']}: {row['SMILES'][:60]}...")
    
    elif args.csv or args.pdb:
        # Batch mode
        csv_path = args.csv or "final_3d_docking_insights.csv"
        variants = optimize_from_results(csv_path, args.top, args.variants)
        
        if len(variants) > 0:
            print(f"\nTotal variants: {len(variants)}")
            print(f"Variant types: {variants['Variant_Type'].value_counts().to_dict()}")
    
    else:
        parser.print_help()
        print("\n  Examples:")
        print("  python optimize_lead.py --smiles 'CC(=O)OC1=CC=CC=C1C(=O)O' --name Aspirin")
        print("  python optimize_lead.py --csv final_3d_docking_insights.csv --top 3")
