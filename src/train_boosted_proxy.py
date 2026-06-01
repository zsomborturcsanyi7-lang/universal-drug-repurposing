#!/usr/bin/env python3
"""
Train Boosted Proxy Model: HistGradientBoostingRegressor.
More accurate than RandomForest — uses tree-based gradient boosting
with the same Morgan fingerprint + physicochemical features.

Training data: training_set_50.csv (SMILES + known Vina affinities)
Output: vina_affinity_proxy_boosted.pkl

Usage:
    python train_boosted_proxy.py [training_set_50.csv]
"""

import os
import sys
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import logging
import numpy as np
import pandas as pd
import joblib
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Reuse feature extraction from train_ml_proxy
from train_ml_proxy import extract_features, get_fingerprint, get_physchem_descriptors


def train_boosted(train_csv: str, model_path: str = "vina_affinity_proxy_boosted.pkl"):
    """
    Train a HistGradientBoostingRegressor for more accurate affinity prediction.
    """
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import mean_absolute_error, r2_score
    
    # ── Load data ──
    if not os.path.exists(train_csv):
        logger.error(f"Training data not found: {train_csv}")
        return None
    
    df = pd.read_csv(train_csv)
    logger.info(f"Loaded {len(df)} training compounds")
    
    # ── Extract features ──
    X, valid_idx = extract_features(df['SMILES'].tolist())
    df_valid = df.iloc[valid_idx].copy()
    y = df_valid['Binding_Affinity_kcal_mol'].values
    
    logger.info(f"Valid compounds: {len(X)}")
    logger.info(f"Feature dimension: {X.shape[1]}")
    
    # ── Train model ──
    model = HistGradientBoostingRegressor(
        loss='squared_error',
        learning_rate=0.05,
        max_iter=500,
        max_depth=8,
        min_samples_leaf=3,
        max_leaf_nodes=31,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=20,
        random_state=42,
        verbose=0,
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
    
    # ── Save ──
    model_package = {
        'model': model,
        'model_type': 'HistGradientBoostingRegressor',
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
    }
    
    joblib.dump(model_package, model_path)
    logger.info(f"Model saved to: {model_path}")
    
    # ── Compare to RandomForest ──
    rf_path = "vina_affinity_proxy.pkl"
    if os.path.exists(rf_path):
        rf_pkg = joblib.load(rf_path)
        rf_info = rf_pkg.get('training_info', {}) if isinstance(rf_pkg, dict) else None
        
        logger.info("\n" + "=" * 50)
        logger.info("  Model Comparison:")
        if rf_info:
            logger.info(f"  RandomForest        MAE: {rf_info.get('training_mae', '?')}  R²: {rf_info.get('training_r2', '?')}")
        logger.info(f"  GradientBoosting    MAE: {mae:.3f} kcal/mol  R²: {r2:.3f}")
        
        if rf_info and mae < rf_info.get('training_mae', float('inf')):
            logger.info("  → GradientBoosting is MORE accurate")
        else:
            logger.info("  → RandomForest is comparable or better (boosting needs more data)")
        logger.info("=" * 50)
    
    return model_package


if __name__ == "__main__":
    if len(sys.argv) > 1:
        train_csv = sys.argv[1]
    else:
        train_csv = "training_set_50.csv"
    
    train_boosted(train_csv)
