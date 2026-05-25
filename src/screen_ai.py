#!/usr/bin/env python3
"""
AI-Powered Drug Screening.
Uses ML proxy to rapidly screen large drug libraries and predict
binding affinity without running actual docking.

Three screening modes:
1. AI Proxy — Predict affinity using trained RandomForest/Boosting model
2. Similarity Search — Tanimoto similarity to known actives
3. Combined — Merge both, rank by consensus

Usage:
    python screen_ai.py                        # Screen full_fda_library.csv
    python screen_ai.py --library custom.csv   # Custom library
    python screen_ai.py --top 20               # Top 20 candidates
    python screen_ai.py --model vina_affinity_proxy_boosted.pkl
"""

import os
import sys
import logging
import argparse
import numpy as np
import pandas as pd
import joblib
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, Descriptors, DataStructs
from rdkit import DataStructs as DS

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# Feature extraction (shared with train_ scripts)
# ──────────────────────────────────────────────

def get_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048):
    """Morgan fingerprint (ECFP4)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return fp_gen.GetFingerprintAsNumPy(mol)


def get_physchem_descriptors(smiles: str):
    """Physicochemical descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [0.0] * 5
    
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
    ]


def extract_features(smiles_list: list):
    """Extract fingerprint + physchem features."""
    features = []
    valid_idx = []
    
    for i, smi in enumerate(smiles_list):
        if not isinstance(smi, str) or not smi.strip():
            continue
        fp = get_fingerprint(smi)
        if fp is None:
            continue
        phys = get_physchem_descriptors(smi)
        features.append(np.concatenate([fp, phys]).astype(np.float32))
        valid_idx.append(i)
    
    if not features:
        return np.array([]).reshape(0, 0), []
    return np.array(features), valid_idx


# ──────────────────────────────────────────────
# AI Proxy screening
# ──────────────────────────────────────────────

def screen_ai_proxy(library_csv: str, model_path: str, top_n: int = 50) -> pd.DataFrame:
    """
    Screen a drug library using the ML affinity predictor.
    """
    if not os.path.exists(model_path):
        logger.error(f"Model not found: {model_path}")
        logger.info("Train it with: python train_ml_proxy.py")
        return pd.DataFrame()
    
    if not os.path.exists(library_csv):
        logger.error(f"Library not found: {library_csv}")
        return pd.DataFrame()
    
    # Load model
    package = joblib.load(model_path)
    if isinstance(package, dict) and 'model' in package:
        model = package['model']
        model_info = package.get('training_info', {})
        logger.info(f"Loaded {package.get('model_type', 'model')}")
        logger.info(f"  Training MAE: {model_info.get('training_mae', '?')} kcal/mol")
    else:
        model = package
    
    # Load library
    df = pd.read_csv(library_csv)
    logger.info(f"Screening {len(df)} compounds from {library_csv}...")
    
    if 'SMILES' not in df.columns:
        logger.error("CSV must have a 'SMILES' column")
        return pd.DataFrame()
    
    name_col = 'Name' if 'Name' in df.columns else df.columns[0]
    
    # Extract features and predict
    X, valid_idx = extract_features(df['SMILES'].tolist())
    
    if X.shape[0] == 0:
        logger.error("No valid SMILES in library")
        return pd.DataFrame()
    
    predictions = model.predict(X)
    
    # Build results
    results = df.iloc[valid_idx].copy()
    results['Predicted_Affinity'] = predictions
    
    # Also compute drug-likeness for ranking
    results['MolWt'] = results['SMILES'].apply(
        lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else 0
    )
    results['LogP'] = results['SMILES'].apply(
        lambda s: Descriptors.MolLogP(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else 0
    )
    
    # Filter drug-like compounds (Lipinski's Rule of 5 heuristic)
    drug_like = (
        (results['MolWt'] < 500) & 
        (results['LogP'] < 5) &
        (results['MolWt'] > 150)
    )
    
    # Rank: best predicted affinity, prefer drug-like
    results['Score'] = results['Predicted_Affinity'] + (drug_like * 0.3)  # Bonus for drug-likeness
    results = results.sort_values('Score')
    
    # Top hits
    hits = results.head(top_n)
    
    logger.info(f"\n  TOP {top_n} AI-PREDICTED HITS:")
    logger.info(f"  {'Rank':5s} {'Name':35s} {'Predicted Affinity'}")
    logger.info(f"  {'-'*5} {'-'*35} {'-'*18}")
    
    for i, (_, row) in enumerate(hits.iterrows()):
        name = str(row.get(name_col, '?'))[:35]
        aff = row['Predicted_Affinity']
        logger.info(f"  {i+1:<5d} {name:35s} {aff:>10.3f} kcal/mol")
    
    # Save
    output = f"ai_screened_hits.csv"
    hits.to_csv(output, index=False)
    logger.info(f"\n  Saved {len(hits)} hits to {output}")
    
    return hits


# ──────────────────────────────────────────────
# Similarity screening
# ──────────────────────────────────────────────

def screen_similarity(library_csv: str, reference_smiles: list = None, 
                       threshold: float = 0.4, top_n: int = 50) -> pd.DataFrame:
    """
    Find structurally similar drugs to known actives.
    """
    if not os.path.exists(library_csv):
        logger.error(f"Library not found: {library_csv}")
        return pd.DataFrame()
    
    df = pd.read_csv(library_csv)
    
    if 'SMILES' not in df.columns:
        logger.error("CSV must have a 'SMILES' column")
        return pd.DataFrame()
    
    # Default reference: Nilotinib (best-known binder for many targets)
    if reference_smiles is None:
        reference_smiles = [
            "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C",  # Nilotinib
            "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(C(=C4)OCC5=CC(=CC=C5)F)Cl",  # Lapatinib
        ]
    
    # Compute reference fingerprints
    ref_fps = []
    for smi in reference_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            ref_fps.append(fp_gen.GetFingerprint(mol))
    
    if not ref_fps:
        logger.error("No valid reference molecules")
        return pd.DataFrame()
    
    # Screen library
    results = []
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    for _, row in df.iterrows():
        smiles = str(row.get('SMILES', ''))
        name = str(row.get('Name', '?'))
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        
        fp = fp_gen.GetFingerprint(mol)
        
        # Max similarity to any reference
        max_sim = max(DataStructs.TanimotoSimilarity(fp, ref) for ref in ref_fps)
        
        if max_sim >= threshold:
            results.append({
                'Name': name,
                'SMILES': smiles,
                'Similarity': max_sim,
            })
    
    if not results:
        logger.warning(f"No molecules above similarity threshold {threshold}")
        return pd.DataFrame()
    
    sim_df = pd.DataFrame(results)
    sim_df = sim_df.sort_values('Similarity', ascending=False)
    sim_df = sim_df.head(top_n)
    
    logger.info(f"\n  TOP {top_n} SIMILARITY HITS (threshold > {threshold}):")
    logger.info(f"  {'Rank':5s} {'Name':35s} {'Tanimoto'}")
    logger.info(f"  {'-'*5} {'-'*35} {'-'*10}")
    
    for i, (_, row) in enumerate(sim_df.iterrows()):
        logger.info(f"  {i+1:<5d} {str(row['Name'])[:35]:35s} {row['Similarity']:>10.3f}")
    
    sim_df.to_csv("similarity_hits.csv", index=False)
    
    return sim_df


# ──────────────────────────────────────────────
# Combined screening
# ──────────────────────────────────────────────

def screen_combined(library_csv: str, model_path: str = None, top_n: int = 50) -> pd.DataFrame:
    """
    Combined screening: AI proxy + similarity.
    Merges results and ranks by consensus.
    """
    all_results = []
    
    # AI proxy
    if model_path and os.path.exists(model_path):
        ai_hits = screen_ai_proxy(library_csv, model_path, top_n=top_n * 2)
        if not ai_hits.empty:
            ai_hits['Source'] = 'AI_Proxy'
            all_results.append(ai_hits)
    else:
        logger.info("No ML model available — using auto-generated predictions")
        # Generate approximate predictions based on drug properties
        df = pd.read_csv(library_csv)
        if 'SMILES' in df.columns:
            name_col = 'Name' if 'Name' in df.columns else df.columns[0]
            
            results = []
            for _, row in df.iterrows():
                smi = str(row.get('SMILES', ''))
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    hbd = Descriptors.NumHDonors(mol)
                    hba = Descriptors.NumHAcceptors(mol)
                    tpsa = Descriptors.TPSA(mol)
                    
                    # Heuristic: drug-like molecules with moderate size bind better
                    # This is a VERY rough approximation
                    score = -5.0 - (0.5 * (300 < mw < 500)) - (0.3 * (logp < 5))
                    
                    results.append({
                        'Name': str(row.get(name_col, '?')),
                        'SMILES': smi,
                        'Predicted_Affinity': score,
                        'Source': 'Heuristic',
                    })
            
            if results:
                heur_df = pd.DataFrame(results)
                heur_df = heur_df.sort_values('Predicted_Affinity')
                heur_df = heur_df.head(top_n * 2)
                all_results.append(heur_df)
    
    # Similarity
    sim_hits = screen_similarity(library_csv, top_n=top_n * 2)
    if not sim_hits.empty:
        sim_hits['Source'] = 'Similarity'
        all_results.append(sim_hits)
    
    if not all_results:
        logger.error("No screening results obtained")
        return pd.DataFrame()
    
    # Merge and deduplicate
    merged = pd.concat(all_results, ignore_index=True)
    
    # Save
    merged.to_csv("screening_results.csv", index=False)
    logger.info(f"\n  Total candidates (before dedup): {len(merged)}")
    
    # Generate docking-ready hits file
    # Prefer AI hits but include similarity discoveries
    docking_candidates = merged.drop_duplicates(subset=['Name'])
    docking_candidates = docking_candidates.head(top_n)
    
    output_cols = ['Name', 'SMILES']
    for col in ['Predicted_Affinity', 'Similarity']:
        if col in docking_candidates.columns:
            output_cols.append(col)
    
    docking_hits = docking_candidates[output_cols]
    docking_hits.to_csv("ai_hits_for_docking.csv", index=False)
    logger.info(f"  Docking candidates: {len(docking_hits)} → ai_hits_for_docking.csv")
    
    return docking_hits


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AI Drug Screening")
    parser.add_argument('--library', default='full_fda_library.csv', 
                       help='Drug library CSV (Name, SMILES columns)')
    parser.add_argument('--model', default='vina_affinity_proxy.pkl',
                       help='ML proxy model')
    parser.add_argument('--top', type=int, default=50, help='Top N candidates')
    parser.add_argument('--mode', default='combined',
                       choices=['ai', 'similarity', 'combined'],
                       help='Screening mode')
    parser.add_argument('--threshold', type=float, default=0.4,
                       help='Tanimoto similarity threshold')
    
    args = parser.parse_args()
    
    if args.mode == 'ai':
        screen_ai_proxy(args.library, args.model, args.top)
    elif args.mode == 'similarity':
        screen_similarity(args.library, threshold=args.threshold, top_n=args.top)
    else:
        screen_combined(args.library, args.model, args.top)
