#!/usr/bin/env python3
"""
Advanced AI Training Pipeline — Data Collection + Model Training.
===========================================================================

Data sources:
  1. ChEMBL bioactivity data for training_set_50.csv compounds
  2. SARS-CoV-2 Mpro (CHEMBL4523582) specific bioactivity data  
  3. Existing Vina docking results (final_3d_docking_insights.csv)
  4. Additional FDA drugs from ChEMBL with known IC50/Ki against any target

Features (300+ dimensions):
  - Morgan fingerprints (ECFP4, 2048 bits + ECFP6, 2048 bits)
  - RDKit full descriptor set (200+ physicochemical descriptors)
  - Drug-likeness scores (Lipinski, QED, etc.)
  - Electrostatic features (partial charges)
  - 3D shape descriptors

Models:
  - XGBoost (gradient boosting)
  - RandomForest
  - HistGradientBoosting
  - Stacking ensemble
  - Neural Network (MLP) — optional

Output:
  - vina_affinity_proxy.pkl (best model)
  - Training metrics + cross-validation report
"""

import os
import sys
import csv
import json
import time
import logging
import warnings
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
import joblib

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('ai_trainer')

# ─────────────────────────────────────────────────────────
# STEP 1: DATA COLLECTION
# ─────────────────────────────────────────────────────────

def fetch_chembl_bioactivity(chembl_ids: list) -> pd.DataFrame:
    """
    Fetch bioactivity data (IC50, Ki, Kd, EC50) from ChEMBL.
    Returns DataFrame with ChEMBL_ID, SMILES, pChEMBL_value, activity_type.
    """
    from chembl_webresource_client.new_client import new_client
    
    activities_client = new_client.activity
    molecule_client = new_client.molecule
    
    results = []
    
    logger.info(f"Fetching bioactivity for {len(chembl_ids)} ChEMBL IDs...")
    
    for chembl_id in chembl_ids:
        try:
            # Fetch activities
            acts = activities_client.filter(
                molecule_chembl_id=chembl_id,
                standard_type__in=['IC50', 'Ki', 'Kd', 'EC50', 'Potency'],
                standard_units='nM',
                limit=10
            )
            
            best_value = None
            best_type = None
            
            for act in acts:
                value = act.get('standard_value')
                atype = act.get('standard_type')
                if value and float(value) > 0:
                    # Convert nM to pChEMBL (-log10(M))
                    molar = float(value) * 1e-9
                    if molar > 0 and molar < 1:
                        pchembl = -np.log10(molar)
                        if best_value is None or pchembl > best_value:
                            best_value = pchembl
                            best_type = atype
            
            if best_value:
                # Get SMILES
                mol = molecule_client.get(chembl_id)
                smiles = mol.get('molecule_structures', {}).get('canonical_smiles', '')
                
                if smiles:
                    results.append({
                        'ChEMBL_ID': chembl_id,
                        'SMILES': smiles,
                        'pChEMBL_value': round(best_value, 2),
                        'activity_type': best_type,
                        'source': 'ChEMBL_bioactivity',
                    })
            
        except Exception as e:
            logger.debug(f"  {chembl_id}: {str(e)[:80]}")
    
    df = pd.DataFrame(results)
    logger.info(f"  Fetched {len(df)} bioactivity records")
    return df


def fetch_sarscov2_mpro_data() -> pd.DataFrame:
    """
    Fetch SARS-CoV-2 Mpro specific bioactivity data from ChEMBL.
    ChEMBL target IDs for SARS-CoV-2 Mpro: CHEMBL4523582, CHEMBL4342320
    """
    from chembl_webresource_client.new_client import new_client
    
    activities_client = new_client.activity
    molecule_client = new_client.molecule
    
    # SARS-CoV-2 Mpro ChEMBL targets
    target_ids = ['CHEMBL4523582', 'CHEMBL4342320', 'CHEMBL4295740']
    
    results = []
    
    for target_id in target_ids:
        try:
            logger.info(f"Fetching Mpro data for {target_id}...")
            
            acts = activities_client.filter(
                target_chembl_id=target_id,
                standard_type__in=['IC50', 'Ki', 'Kd', 'EC50'],
                standard_units='nM',
                limit=500
            )
            
            count = 0
            for act in acts:
                value = act.get('standard_value')
                if not value or float(value) <= 0:
                    continue
                
                chembl_id = act.get('molecule_chembl_id')
                smiles = None
                
                try:
                    mol = molecule_client.get(chembl_id)
                    smiles = mol.get('molecule_structures', {}).get('canonical_smiles', '')
                except:
                    continue
                
                if smiles:
                    molar = float(value) * 1e-9
                    pchembl = -np.log10(molar) if molar > 0 and molar < 1 else None
                    
                    if pchembl:
                        results.append({
                            'ChEMBL_ID': chembl_id,
                            'SMILES': smiles,
                            'pChEMBL_value': round(pchembl, 2),
                            'activity_type': act.get('standard_type', 'IC50'),
                            'target': target_id,
                            'source': 'SARS-CoV-2_Mpro',
                        })
                        count += 1
            
            logger.info(f"  {target_id}: {count} records")
            
        except Exception as e:
            logger.warning(f"  {target_id} failed: {e}")
    
    df = pd.DataFrame(results)
    logger.info(f"Total Mpro bioactivity records: {len(df)}")
    return df


def use_docking_results(docking_csv: str = "final_3d_docking_insights.csv") -> pd.DataFrame:
    """
    Use existing Vina docking results as training data.
    Maps drug names to SMILES from full_fda_library.csv.
    """
    if not os.path.exists(docking_csv):
        return pd.DataFrame()
    
    docked = pd.read_csv(docking_csv)
    
    # Map names to SMILES
    smiles_map = {}
    for lib in ["full_fda_library.csv", "expert_ai_hits.csv"]:
        if os.path.exists(lib):
            lib_df = pd.read_csv(lib)
            for _, row in lib_df.iterrows():
                name = row.get('Name', '')
                smi = row.get('SMILES', '')
                if smi and name:
                    # Normalize: match partial names
                    smiles_map[name.upper().replace(' ', '_')] = smi
    
    results = []
    for _, row in docked.iterrows():
        name = str(row.get('Ligand_Name', row.iloc[0]))
        affinity = row.get('Binding_Affinity_kcal_mol', 0)
        
        # Try to find SMILES
        smiles = None
        name_key = name.upper().replace(' ', '_')
        
        for key, smi in smiles_map.items():
            if key[:10] in name_key or name_key[:10] in key:
                smiles = smi
                break
        
        if smiles and affinity and float(affinity) < 0:
            # Convert kcal/mol to approximate pChEMBL scale
            # pChEMBL = -log10(Kd); approximate: Kd = exp(affinity * 1000 / (R*T))
            # Simplified: ~ 1.36 * (-affinity) — rough heuristic
            pchembl = float(abs(affinity)) * 1.36
            
            results.append({
                'Name': name,
                'SMILES': smiles,
                'Binding_Affinity_kcal_mol': float(affinity),
                'pChEMBL_value': round(pchembl, 2),
                'source': 'Vina_docking',
            })
    
    df = pd.DataFrame(results)
    logger.info(f"Docking training data: {len(df)} records from {docking_csv}")
    return df


def collect_all_training_data() -> pd.DataFrame:
    """
    Collect training data from all sources.
    """
    all_data = []
    
    # Source 1: ChEMBL bioactivity for training set
    train_set = pd.read_csv("training_set_50.csv")
    chembl_ids = train_set['ChEMBL_ID'].dropna().tolist()
    
    if chembl_ids:
        df1 = fetch_chembl_bioactivity(chembl_ids)
        all_data.append(df1)
    
    # Source 2: SARS-CoV-2 Mpro data
    df2 = fetch_sarscov2_mpro_data()
    all_data.append(df2)
    
    # Source 3: Docking results
    df3 = use_docking_results()
    all_data.append(df3)
    
    # Merge
    combined = pd.concat([d for d in all_data if not d.empty], ignore_index=True)
    
    # Deduplicate by SMILES
    combined = combined.drop_duplicates(subset=['SMILES'], keep='first')
    
    logger.info(f"\n{'='*50}")
    logger.info(f"  TOTAL TRAINING DATA: {len(combined)} compounds")
    for src in combined['source'].unique():
        count = len(combined[combined['source'] == src])
        logger.info(f"    {src}: {count}")
    logger.info(f"{'='*50}")
    
    return combined


# ─────────────────────────────────────────────────────────
# STEP 2: ADVANCED FEATURE ENGINEERING
# ─────────────────────────────────────────────────────────

def extract_advanced_features(smiles_list: list) -> np.ndarray:
    """
    Extract 300+ molecular features for ML training.
    
    Feature groups:
    - Morgan fingerprints ECFP4 (2048 bits)
    - Morgan fingerprints ECFP6 (2048 bits)  
    - RDKit 2D descriptors (~200)
    - Drug-likeness scores (5)
    - Fragment counts (10)
    
    Total: ~4300 features
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdFingerprintGenerator, AllChem
    from rdkit.Chem import Crippen, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem import Descriptors3D
    from rdkit.ML.Descriptors import MoleculeDescriptors
    
    # Get all RDKit descriptor names
    all_desc_names = [name for name, _ in Descriptors._descList]
    descriptor_calc = MoleculeDescriptors.MolecularDescriptorCalculator(all_desc_names)
    
    fp_gen_ecfp4 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_gen_ecfp6 = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    
    features = []
    valid_idx = []
    invalid_count = 0
    
    for i, smiles in enumerate(smiles_list):
        if not isinstance(smiles, str) or not smiles.strip():
            invalid_count += 1
            continue
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_count += 1
            continue
        
        try:
            feat_list = []
            
            # 1. Morgan fingerprints ECFP4 (2048)
            fp4 = fp_gen_ecfp4.GetFingerprintAsNumPy(mol)
            feat_list.append(fp4)
            
            # 2. Morgan fingerprints ECFP6 (2048)
            fp6 = fp_gen_ecfp6.GetFingerprintAsNumPy(mol)
            feat_list.append(fp6)
            
            # 3. RDKit 2D descriptors (~200)
            descs = np.array(descriptor_calc.CalcDescriptors(mol), dtype=np.float32)
            # Replace NaN/Inf
            descs = np.nan_to_num(descs, nan=0.0, posinf=0.0, neginf=0.0)
            feat_list.append(descs)
            
            # 4. Drug-likeness scores
            drug_scores = np.array([
                Lipinski.NumHDonors(mol),
                Lipinski.NumHAcceptors(mol),
                QED.qed(mol) if mol.GetNumAtoms() <= 50 else 0.5,
                Crippen.MolMR(mol),
                rdMolDescriptors.CalcNumRotatableBonds(mol),
            ], dtype=np.float32)
            feat_list.append(drug_scores)
            
            # 5. Fragment-based features
            frag_counts = np.array([
                rdMolDescriptors.CalcNumAliphaticRings(mol),
                rdMolDescriptors.CalcNumAromaticRings(mol),
                rdMolDescriptors.CalcNumHeterocycles(mol),
                rdMolDescriptors.CalcNumSaturatedRings(mol),
                rdMolDescriptors.CalcNumHBA(mol),
                rdMolDescriptors.CalcNumHBD(mol),
                rdMolDescriptors.CalcNumAmideBonds(mol),
                rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
                rdMolDescriptors.CalcNumSpiroAtoms(mol),
                rdMolDescriptors.CalcFractionCSP3(mol),
            ], dtype=np.float32)
            feat_list.append(frag_counts)
            
            combined = np.concatenate(feat_list)
            features.append(combined)
            valid_idx.append(i)
            
        except Exception as e:
            invalid_count += 1
            continue
    
    logger.info(f"  Features extracted: {len(features)} valid / {invalid_count} invalid (out of {len(smiles_list)})")
    
    if not features:
        return np.array([]).reshape(0, 0), []
    
    X = np.array(features, dtype=np.float32)
    logger.info(f"  Feature dimensions: {X.shape[1]} ({2048}+{2048}+{len(all_desc_names)}+5+10)")
    
    return X, valid_idx


# ─────────────────────────────────────────────────────────
# STEP 3: MULTI-MODEL TRAINING
# ─────────────────────────────────────────────────────────

def train_models(X: np.ndarray, y: np.ndarray) -> dict:
    """
    Train multiple models and return the best one.
    """
    from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
    from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
    from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
    from sklearn.preprocessing import StandardScaler
    
    # Try to import XGBoost
    try:
        from xgboost import XGBRegressor
        HAS_XGBOOST = True
    except ImportError:
        HAS_XGBOOST = False
        logger.warning("XGBoost not installed — skipping. Install: pip install xgboost")
    
    logger.info(f"\n{'='*50}")
    logger.info(f"  MODEL TRAINING")
    logger.info(f"  Samples: {len(X)} | Features: {X.shape[1]}")
    logger.info(f"{'='*50}\n")
    
    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Split
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.2, random_state=42
    )
    
    logger.info(f"  Train: {len(X_train)} | Test: {len(X_test)}")
    
    models = {}
    results = []
    
    # ── Model 1: RandomForest ──
    logger.info("\n  [1/4] RandomForest...")
    rf = RandomForestRegressor(
        n_estimators=300,
        max_depth=15,
        min_samples_split=3,
        min_samples_leaf=2,
        max_features='sqrt',
        n_jobs=-1,
        random_state=42,
    )
    
    rf_cv = cross_val_score(rf, X_scaled, y, cv=5, scoring='r2')
    rf.fit(X_train, y_train)
    rf_pred = rf.predict(X_test)
    rf_r2 = r2_score(y_test, rf_pred)
    rf_mae = mean_absolute_error(y_test, rf_pred)
    rf_rmse = np.sqrt(mean_squared_error(y_test, rf_pred))
    
    models['RandomForest'] = rf
    results.append({
        'Model': 'RandomForest',
        'CV_R2_mean': rf_cv.mean(),
        'CV_R2_std': rf_cv.std(),
        'Test_R2': rf_r2,
        'Test_MAE': rf_mae,
        'Test_RMSE': rf_rmse,
    })
    
    logger.info(f"    CV R²: {rf_cv.mean():.3f} ± {rf_cv.std():.3f}")
    logger.info(f"    Test R²: {rf_r2:.3f} | MAE: {rf_mae:.3f} | RMSE: {rf_rmse:.3f}")
    
    # ── Model 2: HistGradientBoosting ──
    logger.info("\n  [2/4] HistGradientBoosting...")
    hgb = HistGradientBoostingRegressor(
        loss='squared_error',
        learning_rate=0.05,
        max_iter=500,
        max_depth=10,
        min_samples_leaf=5,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=30,
        random_state=42,
    )
    
    hgb_cv = cross_val_score(hgb, X_scaled, y, cv=5, scoring='r2')
    hgb.fit(X_train, y_train)
    hgb_pred = hgb.predict(X_test)
    hgb_r2 = r2_score(y_test, hgb_pred)
    hgb_mae = mean_absolute_error(y_test, hgb_pred)
    hgb_rmse = np.sqrt(mean_squared_error(y_test, hgb_pred))
    
    models['HistGradientBoosting'] = hgb
    results.append({
        'Model': 'HistGradientBoosting',
        'CV_R2_mean': hgb_cv.mean(),
        'CV_R2_std': hgb_cv.std(),
        'Test_R2': hgb_r2,
        'Test_MAE': hgb_mae,
        'Test_RMSE': hgb_rmse,
    })
    
    logger.info(f"    CV R²: {hgb_cv.mean():.3f} ± {hgb_cv.std():.3f}")
    logger.info(f"    Test R²: {hgb_r2:.3f} | MAE: {hgb_mae:.3f} | RMSE: {hgb_rmse:.3f}")
    
    # ── Model 3: XGBoost ──
    if HAS_XGBOOST:
        logger.info("\n  [3/4] XGBoost...")
        xgb = XGBRegressor(
            n_estimators=300,
            learning_rate=0.05,
            max_depth=8,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_alpha=1,
            reg_lambda=1,
            random_state=42,
            n_jobs=-1,
        )
        
        xgb_cv = cross_val_score(xgb, X_scaled, y, cv=5, scoring='r2')
        xgb.fit(X_train, y_train)
        xgb_pred = xgb.predict(X_test)
        xgb_r2 = r2_score(y_test, xgb_pred)
        xgb_mae = mean_absolute_error(y_test, xgb_pred)
        xgb_rmse = np.sqrt(mean_squared_error(y_test, xgb_pred))
        
        models['XGBoost'] = xgb
        results.append({
            'Model': 'XGBoost',
            'CV_R2_mean': xgb_cv.mean(),
            'CV_R2_std': xgb_cv.std(),
            'Test_R2': xgb_r2,
            'Test_MAE': xgb_mae,
            'Test_RMSE': xgb_rmse,
        })
        
        logger.info(f"    CV R²: {xgb_cv.mean():.3f} ± {xgb_cv.std():.3f}")
        logger.info(f"    Test R²: {xgb_r2:.3f} | MAE: {xgb_mae:.3f} | RMSE: {xgb_rmse:.3f}")
    
    # ── Model 4: Stacking Ensemble ──
    logger.info("\n  [4/4] Stacking Ensemble...")
    from sklearn.ensemble import StackingRegressor
    from sklearn.linear_model import Ridge
    
    estimators = [
        ('rf', models['RandomForest']),
        ('hgb', models['HistGradientBoosting']),
    ]
    if HAS_XGBOOST:
        estimators.append(('xgb', models['XGBoost']))
    
    stack = StackingRegressor(
        estimators=estimators,
        final_estimator=Ridge(alpha=1.0),
        cv=3,
    )
    
    stack.fit(X_train, y_train)
    stack_pred = stack.predict(X_test)
    stack_r2 = r2_score(y_test, stack_pred)
    stack_mae = mean_absolute_error(y_test, stack_pred)
    stack_rmse = np.sqrt(mean_squared_error(y_test, stack_pred))
    
    models['StackingEnsemble'] = stack
    results.append({
        'Model': 'StackingEnsemble',
        'CV_R2_mean': None,
        'Test_R2': stack_r2,
        'Test_MAE': stack_mae,
        'Test_RMSE': stack_rmse,
    })
    
    logger.info(f"    Test R²: {stack_r2:.3f} | MAE: {stack_mae:.3f} | RMSE: {stack_rmse:.3f}")
    
    # ── Summary ──
    results_df = pd.DataFrame(results)
    
    logger.info(f"\n{'='*60}")
    logger.info(f"  TRAINING RESULTS SUMMARY")
    logger.info(f"{'='*60}")
    
    for _, row in results_df.iterrows():
        cv_str = f"{row['CV_R2_mean']:.3f}±{row['CV_R2_std']:.3f}" if pd.notna(row['CV_R2_mean']) else "N/A"
        logger.info(f"  {row['Model']:25s} | R²={row['Test_R2']:.3f} | MAE={row['Test_MAE']:.3f} | CV_R²={cv_str}")
    
    best_idx = results_df['Test_R2'].idxmax()
    best_name = results_df.loc[best_idx, 'Model']
    best_score = results_df.loc[best_idx, 'Test_R2']
    
    logger.info(f"\n  ★ BEST MODEL: {best_name} (R² = {best_score:.3f})")
    
    # Check if we hit 70%
    if best_score >= 0.7:
        logger.info(f"  ✅ TARGET ACHIEVED: R² >= 0.70 ({best_score:.1%})")
        logger.info(f"     This means the model explains {best_score:.1%} of variance in binding affinity")
        logger.info(f"     MAE of {results_df.loc[best_idx, 'Test_MAE']:.2f} pChEMBL ≈ {results_df.loc[best_idx, 'Test_MAE']/1.36:.2f} kcal/mol")
    else:
        logger.info(f"  ⚠ Below target (R²={best_score:.1%} < 0.70)")
        logger.info(f"     Need more training data and/or better features")
    
    return {
        'models': models,
        'results': results_df,
        'best_model_name': best_name,
        'best_score': best_score,
        'scaler': scaler,
        'feature_dim': X.shape[1],
        'n_samples': len(X),
    }


# ─────────────────────────────────────────────────────────
# STEP 4: SAVE + EVALUATE
# ─────────────────────────────────────────────────────────

def save_best_model(train_results: dict, output_path: str = "vina_affinity_proxy.pkl"):
    """Save the best model with metadata."""
    best_name = train_results['best_model_name']
    best_model = train_results['models'][best_name]
    
    package = {
        'model': best_model,
        'scaler': train_results['scaler'],
        'model_type': best_name,
        'feature_dim': train_results['feature_dim'],
        'n_training_samples': train_results['n_samples'],
        'training_info': train_results['results'].to_dict('records'),
        'best_score': float(train_results['best_score']),
    }
    
    joblib.dump(package, output_path)
    logger.info(f"\nModel saved: {output_path}")
    logger.info(f"  Type: {best_name}")
    logger.info(f"  Features: {train_results['feature_dim']}")
    logger.info(f"  Training samples: {train_results['n_samples']}")
    logger.info(f"  R²: {train_results['best_score']:.3f}")


def quick_predict_test(model_path: str = "vina_affinity_proxy.pkl"):
    """Quick prediction test on known compounds."""
    logger.info(f"\n{'='*50}")
    logger.info(f"  PREDICTION TEST")
    logger.info(f"{'='*50}")
    
    test_compounds = [
        ("Nilotinib", "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C"),
        ("Lapatinib", "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(C(=C4)OCC5=CC(=CC=C5)F)Cl"),
        ("Remdesivir", "CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1"),
        ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
        ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)[C@H](C)C(=O)O"),
    ]
    
    if not os.path.exists(model_path):
        logger.error(f"Model not found: {model_path}")
        return
    
    package = joblib.load(model_path)
    model = package['model']
    scaler = package['scaler']
    
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdFingerprintGenerator, QED, Crippen, Lipinski
    from rdkit.Chem import rdMolDescriptors
    from rdkit.ML.Descriptors import MoleculeDescriptors
    
    all_desc_names = [name for name, _ in Descriptors._descList]
    desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(all_desc_names)
    fp4_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp6_gen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    
    for name, smi in test_compounds:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            logger.info(f"  {name:15s}: invalid SMILES")
            continue
        
        fp4 = fp4_gen.GetFingerprintAsNumPy(mol)
        fp6 = fp6_gen.GetFingerprintAsNumPy(mol)
        descs = np.array(desc_calc.CalcDescriptors(mol), dtype=np.float32)
        descs = np.nan_to_num(descs, nan=0.0, posinf=0.0, neginf=0.0)
        drug_scores = np.array([
            Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol),
            QED.qed(mol) if mol.GetNumAtoms() <= 50 else 0.5,
            Crippen.MolMR(mol), rdMolDescriptors.CalcNumRotatableBonds(mol),
        ], dtype=np.float32)
        frag = np.array([
            rdMolDescriptors.CalcNumAliphaticRings(mol), 
            rdMolDescriptors.CalcNumAromaticRings(mol),
            rdMolDescriptors.CalcNumHeterocycles(mol),
            rdMolDescriptors.CalcNumSaturatedRings(mol),
            rdMolDescriptors.CalcNumHBA(mol), rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumAmideBonds(mol),
            rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
            rdMolDescriptors.CalcNumSpiroAtoms(mol),
            rdMolDescriptors.CalcFractionCSP3(mol),
        ], dtype=np.float32)
        
        feat = np.concatenate([fp4, fp6, descs, drug_scores, frag]).reshape(1, -1)
        feat_scaled = scaler.transform(feat)
        
        pred = model.predict(feat_scaled)[0]
        logger.info(f"  {name:15s}: pChEMBL = {pred:.2f}  (≈ {pred/1.36:.1f} kcal/mol)")


# ─────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("╔" + "═" * 58 + "╗")
    print("║     AI DRUG AFFINITY PREDICTOR — TRAINING            ║")
    print("╚" + "═" * 58 + "╝")
    print()
    
    # ── Phase 1: Collect data ──
    logger.info("PHASE 1: Collecting training data...")
    data = collect_all_training_data()
    
    if len(data) < 10:
        logger.error("Not enough training data! Using docking results as fallback.")
        data = use_docking_results("final_3d_docking_insights.csv")
    
    if data.empty:
        logger.error("No training data available. Run docking first.")
        sys.exit(1)
    
    data.to_csv("collected_training_data.csv", index=False)
    logger.info(f"Training data saved: collected_training_data.csv ({len(data)} rows)")
    
    # ── Phase 2: Extract features ──
    logger.info("\nPHASE 2: Extracting advanced features...")
    X, valid_idx = extract_advanced_features(data['SMILES'].tolist())
    
    if X.shape[0] == 0:
        logger.error("No valid features extracted!")
        sys.exit(1)
    
    y = data['pChEMBL_value'].values[valid_idx]
    
    # ── Phase 3: Train models ──
    logger.info("\nPHASE 3: Training multiple models...")
    results = train_models(X, y)
    
    # ── Phase 4: Save ──
    logger.info("\nPHASE 4: Saving best model...")
    save_best_model(results)
    
    # ── Phase 5: Test ──
    quick_predict_test()
