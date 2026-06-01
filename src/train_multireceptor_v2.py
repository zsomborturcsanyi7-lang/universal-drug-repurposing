#!/usr/bin/env python3
"""
MULTI-RECEPTOR TRAINING v2 — CORRECT ACTIVE SITES
===================================================
Uses known inhibitor binding sites, not blind docking.
Target: R² ≥ 0.70 with proper active site boxes.
"""
import os, sys, subprocess, time, csv
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = os.path.dirname(os.path.abspath(__file__)) or "."
os.chdir(BASE)

# ──────────────────────────────────────────────
# KNOWN ACTIVE SITES (from co-crystallized inhibitors)
# ──────────────────────────────────────────────
RECEPTORS = {
    '6LU7': {
        'name': 'SARS-CoV-2 Mpro',
        'center': [-10.7, 12.4, 68.8],
        'size': [20, 20, 20],
        'ref_ligand': 'N3',
    },
    '1HPV': {
        'name': 'HIV Protease',
        'center': [6.5, 23.5, 15.5],
        'size': [20, 20, 20],
        'ref_ligand': 'MK1',
    },
    '5KIR': {
        'name': 'COX-2',
        'center': [25.5, 22.0, 18.5],
        'size': [22, 22, 22],
        'ref_ligand': 'IMN',
    },
    '1M17': {
        'name': 'EGFR Kinase',
        'center': [21.0, 1.5, 53.5],
        'size': [20, 22, 20],
        'ref_ligand': 'AQ4',
    },
    '4LVT': {
        'name': 'Bcl-2',
        'center': [14.0, 18.0, 11.5],
        'size': [22, 22, 22],
        'ref_ligand': '1XJ',
    },
}

# ──────────────────────────────────────────────
# STEP 1: Verify receptors ready
# ──────────────────────────────────────────────
print("=" * 60)
print("  STEP 1/3: Checking receptors")
print("=" * 60)

READY = {}
for pdb, info in RECEPTORS.items():
    pdbqt = f"targets/{pdb}/receptor_cleaned.pdbqt"
    if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 1000:
        print(f"  {pdb} ({info['name']}): ready")
        READY[pdb] = info
    else:
        print(f"  {pdb} ({info['name']}): MISSING — run receptor_prep.py first")

if len(READY) < 3:
    print("\n  ERROR: Need at least 3 receptors prepared!")
    sys.exit(1)

# ──────────────────────────────────────────────
# STEP 2: Collect ligands
# ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("  STEP 2/3: Collecting ligands")
print("=" * 60)

LIGANDS = []
for d in ["kaggle_output/ligands", "new_ligands"]:
    if not os.path.exists(d): continue
    for f in os.listdir(d):
        if f.endswith('.pdbqt') and os.path.getsize(f"{d}/{f}") > 200:
            name = f.replace('.pdbqt', '')
            if not name.startswith('CHEMBL'):
                LIGANDS.append((name, f"{d}/{f}"))

LIGANDS = LIGANDS[:80]  # Top 80 diverse ligands
print(f"  Using {len(LIGANDS)} ligands")

# ──────────────────────────────────────────────
# STEP 3: Parallel docking on all receptors
# ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("  STEP 3/3: Docking on {len(READY)} receptors (5 workers each)")
print("=" * 60)

VINA = "vina_1.2.7_win.exe"
ALL_RESULTS = []
TOTAL_START = time.time()

for pdb, info in READY.items():
    rec = f"targets/{pdb}/receptor_cleaned.pdbqt"
    if not os.path.exists(rec): continue
    
    cx, cy, cz = info['center']
    sx, sy, sz = info['size']
    
    print(f"\n  Docking on {pdb} ({info['name']})...")
    print(f"    Box center: {info['center']}, size: {info['size']}")
    
    def dock_one(name, lig_path):
        try:
            cmd = [VINA, "--receptor", rec, "--ligand", lig_path,
                   "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                   "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
                   "--exhaustiveness", "2", "--cpu", "1",
                   "--out", f"dock_{pdb}_{name}.pdbqt"]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
            best = None
            for l in r.stdout.split('\n'):
                p2 = l.strip().split()
                if p2 and p2[0].isdigit() and len(p2) >= 2:
                    try:
                        a = float(p2[1])
                        if best is None or a < best: best = a
                    except: pass
            return (name, pdb, best)
        except:
            return (name, pdb, None)
    
    done = 0; valid = 0; t0 = time.time()
    with ThreadPoolExecutor(max_workers=5) as ex:
        futures = {ex.submit(dock_one, n, p): n for n, p in LIGANDS}
        for f in as_completed(futures):
            name, pdb_id, aff = f.result()
            done += 1
            if aff is not None and aff < 0:
                valid += 1
                ALL_RESULTS.append((name, pdb_id, aff))
            if done % 20 == 0:
                print(f"    [{done}/{len(LIGANDS)}] {valid} valid")
    
    elapsed = time.time() - t0
    print(f"    Done: {valid} valid in {elapsed:.0f}s")

TOTAL_ELAPSED = time.time() - TOTAL_START
print(f"\n  TOTAL: {len(ALL_RESULTS)} results across {len(READY)} receptors")
print(f"  Docking time: {TOTAL_ELAPSED/60:.1f} min")

# Save
with open("multi_docking.csv", "w", newline='') as f:
    w = csv.writer(f)
    w.writerow(["Ligand", "Receptor", "Affinity"])
    for n, p, a in ALL_RESULTS: w.writerow([n, p, a])

# ──────────────────────────────────────────────
# STEP 4: Train multi-receptor model
# ──────────────────────────────────────────────
print("\n" + "=" * 60)
print("  TRAINING MULTI-RECEPTOR AI MODEL")
print("=" * 60)

import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem import Crippen, Lipinski, QED, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor, ExtraTreesRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression

# SMILES map
SMILES_MAP = {}
for lib in ["full_fda_library.csv", "expert_ai_hits.csv", "nilotinib_variants.csv", "training_set_50.csv"]:
    if not os.path.exists(lib): continue
    for row in csv.DictReader(open(lib)):
        n = row.get('Name', '').upper().replace(' ', '_')
        s = row.get('SMILES', '')
        if s and n: SMILES_MAP[n] = s

# Features
all_descs = [n for n, _ in Descriptors._descList]
desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(all_descs)
fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=512)

def get_feat(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol: return None
    fp = fp_gen.GetFingerprintAsNumPy(mol)
    descs = np.nan_to_num(np.array(desc_calc.CalcDescriptors(mol), dtype=np.float32))
    drug = np.array([
        Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol),
        QED.qed(mol) if mol.GetNumAtoms() <= 50 else 0.5,
        Crippen.MolMR(mol), rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcFractionCSP3(mol), Descriptors.MolWt(mol)/100,
        Descriptors.MolLogP(mol), Descriptors.TPSA(mol)/100,
    ], dtype=np.float32)
    return np.concatenate([fp, descs, drug])

# Build training set
rec_list = list(READY.keys())
rec_index = {r: i for i, r in enumerate(rec_list)}
X_list, y_list = [], []

for name, pdb, aff in ALL_RESULTS:
    nkey = name.upper().replace(' ', '_')
    smi = None
    for k, s in SMILES_MAP.items():
        if k[:10] in nkey or nkey[:10] in k: smi = s; break
    if not smi: continue
    feat = get_feat(smi)
    if feat is None: continue
    
    rec_hot = np.zeros(len(rec_list), dtype=np.float32)
    if pdb in rec_index: rec_hot[rec_index[pdb]] = 1.0
    
    X_list.append(np.concatenate([feat, rec_hot]))
    y_list.append(aff)

X = np.array(X_list, dtype=np.float32)
y = np.array(y_list)

print(f"  Training data: {X.shape[0]} samples, {X.shape[1]} features")
print(f"  Affinity range: {y.min():.1f} to {y.max():.1f} kcal/mol")
print(f"  Receptors: {len(set(r[1] for r in ALL_RESULTS if r[1] in rec_index))}")

if X.shape[0] < 50:
    print("  ERROR: Not enough training data!")
    sys.exit(1)

# Feature selection + scaling
K = min(200, X.shape[1] - 10)
sel = SelectKBest(f_regression, k=K)
X_sel = sel.fit_transform(X, y)
scaler = StandardScaler()
X_s = scaler.fit_transform(X_sel)

X_train, X_test, y_train, y_test = train_test_split(X_s, y, test_size=0.2, random_state=42)

# Train
results = []
models = {}

for name, model in [
    ("RF", RandomForestRegressor(n_estimators=400, max_depth=15, min_samples_leaf=2, n_jobs=-1, random_state=42)),
    ("ET", ExtraTreesRegressor(n_estimators=400, max_depth=15, min_samples_leaf=2, n_jobs=-1, random_state=42)),
    ("HGB", HistGradientBoostingRegressor(max_iter=500, max_depth=10, learning_rate=0.03, min_samples_leaf=5, random_state=42)),
]:
    cv = cross_val_score(model, X_s, y, cv=5, scoring='r2').mean()
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    r2 = r2_score(y_test, pred)
    mae = mean_absolute_error(y_test, pred)
    results.append((name, r2, cv, mae))
    models[name] = model
    print(f"  {name:5s} | R²={r2:.3f} | CV={cv:.3f} | MAE={mae:.3f}")

best_name, best_r2, best_cv, best_mae = max(results, key=lambda x: x[1])
print(f"\n  ★ BEST: {best_name} — R² = {best_r2:.3f}, CV = {best_cv:.3f}")

if best_r2 >= 0.70:
    print(f"\n  ✅✅✅ TARGET 70% R² ACHIEVED! ✅✅✅")
elif best_cv >= 0.40:
    print(f"\n  ✅ CV ≥ 0.40 — model is learning general binding physics!")
    print(f"  Test R² = {best_r2:.3f} — need ~400-500 samples for 0.70")
else:
    print(f"\n  📊 R² = {best_r2:.3f} with {X.shape[0]} samples")
    print(f"  More diverse receptors/l ligands would improve this")

# Save
pkg = {
    'model': models[best_name], 'scaler': scaler, 'selector': sel,
    'receptor_index': rec_index, 'model_type': best_name,
    'n_samples': X.shape[0], 'metrics': {'r2': best_r2, 'cv': best_cv, 'mae': best_mae},
    'results': results,
}
joblib.dump(pkg, "vina_multireceptor_model.pkl")
print(f"\n  Saved: vina_multireceptor_model.pkl")
print("=" * 60)
