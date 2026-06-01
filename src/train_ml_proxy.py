#!/usr/bin/env python3
"""
Train ML Proxy Model: RandomForestRegressor.
Predicts AutoDock Vina binding affinity from SMILES strings.

Features: Morgan fingerprints (ECFP4, radius=2, 2048 bits) 
          + 5 physicochemical descriptors

Training data: training_set_50.csv (SMILES + known Vina affinities)
Output: vina_affinity_proxy.pkl

Usage:
    python train_ml_proxy.py [training_set_50.csv]
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import joblib
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# Feature extraction
# ──────────────────────────────────────────────

def get_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048):
    """Generate Morgan fingerprint (ECFP4) for a molecule."""
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return fp_gen.GetFingerprintAsNumPy(mol)


def get_physchem_descriptors(smiles: str):
    """Extract 5 key physicochemical descriptors."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [0.0] * 5
    
    return [
        Descriptors.MolWt(mol),           # Molecular weight
        Descriptors.MolLogP(mol),         # Lipophilicity
        Descriptors.NumHDonors(mol),      # H-bond donors
        Descriptors.NumHAcceptors(mol),   # H-bond acceptors
        Descriptors.TPSA(mol),            # Polar surface area
    ]


def extract_features(smiles_list: list):
    """
    Convert list of SMILES to feature matrix.
    Returns (X, valid_indices) — invalid SMILES are skipped.
    """
    features = []
    valid_idx = []
    
    for i, smi in enumerate(smiles_list):
        if not isinstance(smi, str) or not smi.strip():
            continue
        
        fp = get_fingerprint(smi)
        if fp is None:
            continue
        
        phys = get_physchem_descriptors(smi)
        feat = np.concatenate([fp, phys]).astype(np.float32)
        features.append(feat)
        valid_idx.append(i)
    
    if not features:
        return np.array([]).reshape(0, 0), []
    
    return np.array(features), valid_idx


# ──────────────────────────────────────────────
# Training
# ──────────────────────────────────────────────

def train_random_forest(train_csv: str, model_path: str = "vina_affinity_proxy.pkl"):
    """
    Train a RandomForestRegressor to predict Vina binding affinity.
    """
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import mean_absolute_error, r2_score
    
    # ── Load data ──
    if not os.path.exists(train_csv):
        logger.error(f"Training data not found: {train_csv}")
        logger.info("Generate it with: python generate_offline_lib.py")
        return None
    
    df = pd.read_csv(train_csv)
    logger.info(f"Loaded {len(df)} training compounds")
    
    # ── Extract features ──
    X, valid_idx = extract_features(df['SMILES'].tolist())
    df_valid = df.iloc[valid_idx].copy()
    y = df_valid['Binding_Affinity_kcal_mol'].values  # Column name from training set
    
    logger.info(f"Valid compounds: {len(X)}")
    logger.info(f"Feature dimension: {X.shape[1]} (2048 fingerprints + 5 physchem)")
    
    # ── Train model ──
    model = RandomForestRegressor(
        n_estimators=200,
        max_depth=20,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features='sqrt',
        random_state=42,
        n_jobs=-1,
    )
    
    # Cross-validation
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='neg_mean_absolute_error')
    logger.info(f"CV MAE: {-cv_scores.mean():.3f} ± {cv_scores.std():.3f} kcal/mol")
    
    # Final fit
    model.fit(X, y)
    
    # Training metrics
    y_pred = model.predict(X)
    mae = mean_absolute_error(y, y_pred)
    r2 = r2_score(y, y_pred)
    
    logger.info(f"Training MAE: {mae:.3f} kcal/mol")
    logger.info(f"Training R²: {r2:.3f}")
    
    # ── Feature importance ──
    importances = model.feature_importances_
    fp_importance = np.sum(importances[:2048])
    phys_importance = np.sum(importances[2048:])
    logger.info(f"Feature importance: Fingerprint={fp_importance:.1%}, PhysChem={phys_importance:.1%}")
    
    # ── Save ──
    # Wrap model with metadata
    model_package = {
        'model': model,
        'model_type': 'RandomForestRegressor',
        'feature_config': {
            'fingerprint_radius': 2,
            'fingerprint_bits': 2048,
            'physchem_descriptors': ['MolWt', 'MolLogP', 'NumHDonors', 'NumHAcceptors', 'TPSA'],
        },
        'training_info': {
            'n_samples': len(X),
            'n_features': X.shape[1],
            'cv_mae_mean': float(-cv_scores.mean()),
            'cv_mae_std': float(cv_scores.std()),
            'training_mae': float(mae),
            'training_r2': float(r2),
        },
        'predict': None,  # Will be replaced by a cleaner API
    }
    
    joblib.dump(model_package, model_path)
    logger.info(f"Model saved to: {model_path}")
    
    return model_package


def predict_affinity(smiles_list: list, model_path: str = "vina_affinity_proxy.pkl"):
    """
    Predict binding affinities for a list of SMILES strings.
    Returns: list of (affinity, None) tuples — None for invalid SMILES.
    """
    if not os.path.exists(model_path):
        logger.error(f"Model not found: {model_path}")
        return []
    
    package = joblib.load(model_path)
    
    # Handle both new format (dict) and old format (raw model)
    if isinstance(package, dict) and 'model' in package:
        model = package['model']
    else:
        model = package
    
    # Handle both package types
    if hasattr(model, 'predict'):
        X, valid_idx = extract_features(smiles_list)
        if X.shape[0] == 0:
            return [None] * len(smiles_list)
        
        preds = model.predict(X)
        
        results = [None] * len(smiles_list)
        for vi, pi in zip(valid_idx, preds):
            results[vi] = float(pi)
        return results
    else:
        logger.error("Model doesn't have predict method")
        return [None] * len(smiles_list)


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) > 1:
        train_csv = sys.argv[1]
    else:
        train_csv = "training_set_50.csv"
    
    train_random_forest(train_csv)
    
    # Quick test
    print("\n" + "=" * 50)
    print("  Quick predictions:")
    
    test_smiles = [
        ("Nilotinib", "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C"),
        ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ]
    
    for name, smi in test_smiles:
        pred = predict_affinity([smi], "vina_affinity_proxy.pkl")
        if pred and pred[0] is not None:
            print(f"  {name}: {pred[0]:.3f} kcal/mol (predicted)")
        else:
            print(f"  {name}: prediction failed")
    
    print("=" * 50)
