#!/usr/bin/env python3
"""
Ligand Preparation: SMILES → 3D conformer → PDBQT for AutoDock Vina.
Uses RDKit for 3D generation and OpenBabel for PDBQT conversion.
Cross-platform: auto-detects OpenBabel on Windows, Linux, and macOS.
"""

import os
import sys
import platform
import shutil
import subprocess
import logging
import csv
import tempfile
from io import StringIO

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# OpenBabel path detection
# ──────────────────────────────────────────────

def _find_obabel():
    """Find OpenBabel executable. Returns path or raises FileNotFoundError."""
    system = platform.system()
    
    if system == "Windows":
        # Check known install locations
        candidates = [
            r"C:\Program Files\OpenBabel-3.1.1\obabel.exe",
            r"C:\Program Files (x86)\OpenBabel-3.1.1\obabel.exe",
            r"C:\Program Files\OpenBabel\obabel.exe",
            r"C:\Program Files (x86)\OpenBabel\obabel.exe",
            r"C:\OpenBabel\obabel.exe",
        ]
        for p in candidates:
            if os.path.exists(p):
                return p
        # Check PATH
        found = shutil.which("obabel.exe") or shutil.which("obabel")
        if found:
            return found
        raise FileNotFoundError(
            "OpenBabel not found on Windows. "
            "Download from https://github.com/openbabel/openbabel/releases "
            "and install to C:\\Program Files\\OpenBabel-3.1.1\\"
        )
    else:
        # Linux/macOS: use system obabel
        found = shutil.which("obabel")
        if found:
            return found
        raise FileNotFoundError(
            "OpenBabel not found. Install with: "
            "sudo apt-get install openbabel  (Linux)  or  brew install open-babel  (macOS)"
        )

OBABEL_PATH = None  # Lazily detected on first use

def get_obabel():
    global OBABEL_PATH
    if OBABEL_PATH is None:
        OBABEL_PATH = _find_obabel()
        logger.info(f"OpenBabel found at: {OBABEL_PATH}")
    return OBABEL_PATH

# ──────────────────────────────────────────────
# Core preparation functions
# ──────────────────────────────────────────────

def smiles_to_3d_pdb(smiles: str, name: str, output_dir: str = "ligands") -> str:
    """
    Convert a SMILES string to a 3D PDB file using RDKit.
    Returns path to the generated PDB file.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for '{name}': {smiles}")

    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)
    
    if status != 0:
        # Fallback: try with different params
        params = AllChem.ETKDG()
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise ValueError(f"Failed to generate 3D conformer for '{name}'")

    # Optimize geometry
    AllChem.MMFFOptimizeMolecule(mol)

    os.makedirs(output_dir, exist_ok=True)
    pdb_path = os.path.join(output_dir, f"{name}.pdb")
    Chem.MolToPDBFile(mol, pdb_path)
    
    return pdb_path


def pdb_to_pdbqt(pdb_path: str, output_path: str) -> bool:
    """
    Convert PDB to PDBQT using OpenBabel.
    Adds Gasteiger charges (required by AutoDock Vina).
    """
    obabel = get_obabel()
    cmd = [
        obabel, pdb_path,
        "-O", output_path,
        "-xr",                          # Preserve residue/atom names
        "--partialcharge", "gasteiger", # Required for Vina
        "-h",                           # Add hydrogens
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            logger.error(f"OpenBabel failed for {pdb_path}: {result.stderr.strip()}")
            return False
        return True
    except subprocess.TimeoutExpired:
        logger.error(f"OpenBabel timeout for {pdb_path}")
        return False
    except Exception as e:
        logger.error(f"OpenBabel error for {pdb_path}: {e}")
        return False


def prepare_ligand(smiles: str, name: str, output_dir: str = "ligands") -> str:
    """
    Full pipeline: SMILES → 3D → PDBQT.
    Returns path to the PDBQT file, or None on failure.
    """
    # Sanitize name for filesystem
    safe_name = name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    pdbqt_path = os.path.join(output_dir, f"{safe_name}.pdbqt")
    
    # Skip if already prepared
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 0:
        logger.debug(f"Skipping {name} — already prepared")
        return pdbqt_path
    
    try:
        # Step 1: SMILES → 3D PDB
        pdb_path = smiles_to_3d_pdb(smiles, safe_name, output_dir)
        
        # Step 2: PDB → PDBQT
        if not pdb_to_pdbqt(pdb_path, pdbqt_path):
            return None
        
        # Cleanup intermediate PDB
        if os.path.exists(pdb_path):
            os.remove(pdb_path)
            
        logger.info(f"Prepared: {name} → {pdbqt_path}")
        return pdbqt_path
        
    except Exception as e:
        logger.error(f"Failed to prepare {name}: {e}")
        return None


def prepare_ligands(csv_path: str, output_dir: str = "ligands") -> dict:
    """
    Batch prepare all ligands from a CSV file.
    CSV must have columns: Name, SMILES
    
    Returns: {name: pdbqt_path or None for failures}
    """
    if not os.path.exists(csv_path):
        logger.error(f"CSV not found: {csv_path}")
        return {}
    
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        logger.warning(f"Empty CSV: {csv_path}")
        return {}
    
    required_cols = {'Name', 'SMILES'}
    if not required_cols.issubset(reader.fieldnames or []):
        logger.error(f"CSV missing required columns. Need: {required_cols}")
        return {}
    
    logger.info(f"Starting preparation of {len(rows)} ligands from {os.path.basename(csv_path)}...")
    
    results = {}
    success = 0
    for row in rows:
        name = row.get('Name', '').strip()
        smiles = row.get('SMILES', '').strip()
        
        if not name or not smiles:
            continue
            
        pdbqt_path = prepare_ligand(smiles, name, output_dir)
        results[name] = pdbqt_path
        if pdbqt_path:
            success += 1
    
    logger.info(f"Successfully prepared {success} ligands.")
    return results


# ──────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        # Default: prepare the AI hits
        csv_file = "ai_hits_for_docking.csv"
        if not os.path.exists(csv_file):
            csv_file = "expert_ai_hits.csv"
        if not os.path.exists(csv_file):
            print("Usage: python prep_ligands.py <csv_file>")
            print("CSV must have columns: Name, SMILES")
            sys.exit(1)
    
    results = prepare_ligands(csv_file)
    
    if results:
        prepared = sum(1 for v in results.values() if v is not None)
        print(f"\nDone. Prepared {prepared}/{len(results)} ligands.")
    else:
        print(f"\nNo ligands prepared. Check that {csv_file} has Name and SMILES columns.")
