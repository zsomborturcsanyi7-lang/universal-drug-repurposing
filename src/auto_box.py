#!/usr/bin/env python3
"""
Automatic Docking Box Detection.
Replaces hardcoded coordinates with automatic binding pocket detection.

Strategies (in priority order):
1. Reference ligand — if co-crystallized ligand exists, box around it + 8A
2. FPocket — external cavity detection tool (if installed)
3. P2Rank — ML-based pocket prediction (if installed)  
4. Grid Blind Docking — cover entire protein surface (fallback)

Usage:
    python auto_box.py --pdb 6LU7
    python auto_box.py --pdb 6LU7 --method reference --ligand N3
    python auto_box.py --receptor targets/6LU7/receptor_cleaned.pdb
"""

import os
import sys
import argparse
import logging
import yaml
import subprocess
import shutil
import numpy as np
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────
# Strategy 1: Reference ligand box
# ──────────────────────────────────────────────

def box_from_reference_ligand(pdb_path: str, ligand_resname: str = None, 
                                padding: float = 8.0) -> dict:
    """
    Extract docking box from a co-crystallized ligand.
    
    Args:
        pdb_path: Path to PDB file
        ligand_resname: 3-letter residue name of the ligand (e.g., 'N3')
                       If None, auto-detect from HETATM records
        padding: Extra padding around ligand in Angstroms
    
    Returns:
        {'center': [x, y, z], 'size': [x, y, z], 'method': 'reference_ligand'}
        or None if no ligand found
    """
    coords = []
    found_resname = None
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('HETATM'):
                resname = line[17:20].strip()
                
                # Skip water and ions
                if resname in ('HOH', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'SO4'):
                    continue
                
                if ligand_resname and resname != ligand_resname:
                    continue
                
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                    found_resname = resname
                except ValueError:
                    continue
    
    if not coords:
        logger.warning(f"No reference ligand found in {pdb_path}")
        return None
    
    coords = np.array(coords)
    center = coords.mean(axis=0)
    
    # Box size: span of ligand + padding
    ligand_size = coords.max(axis=0) - coords.min(axis=0)
    box_size = ligand_size + padding
    
    # Ensure minimum box size
    box_size = np.maximum(box_size, 15.0)
    # Cap maximum
    box_size = np.minimum(box_size, 40.0)
    
    result = {
        'center': center.round(1).tolist(),
        'size': box_size.round(1).tolist(),
        'method': 'reference_ligand',
        'ligand': found_resname,
        'ligand_atoms': len(coords),
        'padding': padding,
    }
    
    logger.info(f"Box from reference ligand '{found_resname}':")
    logger.info(f"  Center: {result['center']}")
    logger.info(f"  Size:   {result['size']}")
    logger.info(f"  Atoms:  {result['ligand_atoms']}")
    
    return result


# ──────────────────────────────────────────────
# Strategy 2: FPocket
# ──────────────────────────────────────────────

def box_from_fpocket(pdb_path: str, pocket_index: int = 0) -> dict:
    """
    Run fpocket to detect binding pockets.
    Requires: fpocket installed (apt-get install fpocket)
    """
    if not shutil.which('fpocket'):
        logger.warning("FPocket not installed. Install with: apt-get install fpocket")
        logger.warning("Install from: https://github.com/Discngine/fpocket")
        return None
    
    # Run fpocket
    cmd = ['fpocket', '-f', pdb_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        # fpocket output goes to <pdb_name>_out/
        base = os.path.splitext(pdb_path)[0]
        out_dir = f"{base}_out"
        pockets_dir = os.path.join(out_dir, "pockets")
        
        if not os.path.exists(pockets_dir):
            logger.warning("FPocket produced no pockets")
            return None
        
        # Read pocket info
        info_file = os.path.join(out_dir, f"{os.path.basename(base)}_info.txt")
        if not os.path.exists(info_file):
            return None
        
        # Parse pockets — pick the largest one (usually the active site)
        pockets = []
        with open(info_file) as f:
            pocket = None
            for line in f:
                if 'Pocket' in line and ':' in line:
                    if pocket:
                        pockets.append(pocket)
                    pocket = {'id': len(pockets)}
                elif 'Score' in line and pocket:
                    try:
                        pocket['score'] = float(line.split(':')[1])
                    except:
                        pass
                elif 'Druggability Score' in line and pocket:
                    try:
                        pocket['druggability'] = float(line.split(':')[1])
                    except:
                        pass
            if pocket:
                pockets.append(pocket)
        
        if not pockets:
            return None
        
        # Use the best pocket (highest druggability or score)
        best_pocket = max(pockets, key=lambda p: p.get('druggability', p.get('score', 0)))
        pocket_id = best_pocket['id']
        
        # Read pocket PDB to get coordinates
        pocket_pdb = os.path.join(pockets_dir, f"pocket{pocket_id}_atm.pdb")
        if not os.path.exists(pocket_pdb):
            return None
        
        coords = []
        with open(pocket_pdb) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError:
                        continue
        
        if not coords:
            return None
        
        coords = np.array(coords)
        center = coords.mean(axis=0)
        size = (coords.max(axis=0) - coords.min(axis=0)) + 6.0
        size = np.maximum(size, 18.0)
        size = np.minimum(size, 35.0)
        
        result = {
            'center': center.round(1).tolist(),
            'size': size.round(1).tolist(),
            'method': 'fpocket',
            'pocket_score': round(best_pocket.get('score', 0), 2),
            'druggability': round(best_pocket.get('druggability', 0), 2),
        }
        
        logger.info(f"Box from FPocket (pocket {pocket_id}):")
        logger.info(f"  Center: {result['center']}")
        logger.info(f"  Size:   {result['size']}")
        
        # Cleanup
        shutil.rmtree(out_dir, ignore_errors=True)
        
        return result
        
    except Exception as e:
        logger.warning(f"FPocket error: {e}")
        return None


# ──────────────────────────────────────────────
# Strategy 3: P2Rank (ML-based)
# ──────────────────────────────────────────────

def box_from_p2rank(pdb_path: str) -> dict:
    """
    Use P2Rank for ML-based binding site prediction.
    Requires: p2rank installed (pip install p2rank)
    """
    try:
        import p2rank
    except ImportError:
        logger.debug("P2Rank not installed")
        return None
    
    # P2Rank typically runs as CLI tool
    # This is a simplified wrapper
    logger.info("P2Rank integration not yet implemented — install from https://github.com/rdk/p2rank")
    return None


# ──────────────────────────────────────────────
# Strategy 4: Grid blind docking
# ──────────────────────────────────────────────

def box_from_blind(pdb_path: str) -> dict:
    """
    Cover the entire protein surface.
    Used as last-resort fallback.
    """
    coords = []
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                # Only use CA atoms for approximate surface
                if line[12:16].strip() == 'CA':
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError:
                        continue
    
    if not coords:
        # Fallback: use all atoms
        with open(pdb_path) as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError:
                        continue
    
    if not coords:
        logger.error(f"No atoms found in {pdb_path}")
        return None
    
    coords = np.array(coords)
    center = coords.mean(axis=0)
    size = (coords.max(axis=0) - coords.min(axis=0)) + 4.0
    
    # Cap at reasonable size for Vina
    size = np.minimum(size, 50.0)
    
    result = {
        'center': center.round(1).tolist(),
        'size': size.round(1).tolist(),
        'method': 'blind_docking',
        'note': 'Covers entire protein — slower but guaranteed to find binding site',
    }
    
    logger.info(f"Box from blind docking (CA atoms):")
    logger.info(f"  Center: {result['center']}")
    logger.info(f"  Size:   {result['size']}")
    
    return result


# ──────────────────────────────────────────────
# Auto-detect (try strategies in order)
# ──────────────────────────────────────────────

def auto_detect_box(pdb_path: str, method: str = 'auto', 
                     ligand_resname: str = None) -> dict:
    """
    Auto-detect the optimal docking box using available strategies.
    
    Args:
        pdb_path: Path to cleaned PDB file
        method: 'auto', 'reference', 'fpocket', 'p2rank', 'blind'
        ligand_resname: Ligand residue name for reference method
    
    Returns:
        Box dict with center, size, and method used
    """
    if method == 'reference' or (method == 'auto' and ligand_resname):
        box = box_from_reference_ligand(pdb_path, ligand_resname)
        if box:
            return box
        if method == 'reference':
            logger.error("Reference ligand not found")
            return None
    
    if method in ('fpocket', 'auto'):
        box = box_from_fpocket(pdb_path)
        if box:
            return box
    
    if method in ('p2rank', 'auto'):
        box = box_from_p2rank(pdb_path)
        if box:
            return box
    
    if method in ('blind', 'auto'):
        box = box_from_blind(pdb_path)
        if box:
            return box
    
    return None


def box_to_vina_config(box: dict, receptor_pdbqt: str) -> str:
    """Convert box dict to AutoDock Vina config file content."""
    center = box['center']
    size = box['size']
    
    return f"""# Auto-generated by auto_box.py (method: {box.get('method', 'unknown')})
receptor = {receptor_pdbqt}
center_x = {center[0]}
center_y = {center[1]}
center_z = {center[2]}
size_x = {size[0]}
size_y = {size[1]}
size_z = {size[2]}
exhaustiveness = 8
num_modes = 9
energy_range = 3
cpu = 0
"""


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto-detect docking box")
    parser.add_argument('--pdb', help='PDB ID (looks in targets/<id>/)')
    parser.add_argument('--receptor', help='Direct path to PDB file')
    parser.add_argument('--method', default='auto', 
                       choices=['auto', 'reference', 'fpocket', 'blind'],
                       help='Detection method')
    parser.add_argument('--ligand', help='Reference ligand residue name (e.g., N3)')
    parser.add_argument('--output', help='Output config file path')
    
    args = parser.parse_args()
    
    # Determine PDB path
    if args.receptor:
        pdb_path = args.receptor
    elif args.pdb:
        pdb_id = args.pdb.upper()
        pdb_path = f"targets/{pdb_id}/receptor_cleaned.pdb"
        if not os.path.exists(pdb_path):
            logger.error(f"Receptor not prepared. Run: python receptor_prep.py --pdb {pdb_id}")
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)
    
    if not os.path.exists(pdb_path):
        logger.error(f"File not found: {pdb_path}")
        sys.exit(1)
    
    # Detect box
    box = auto_detect_box(pdb_path, args.method, args.ligand)
    
    if box is None:
        logger.error("Could not detect docking box with any method")
        sys.exit(1)
    
    # Print
    print(f"\n  Method: {box['method']}")
    print(f"  Center: {box['center']}")
    print(f"  Size:   {box['size']}")
    
    # Save config
    if args.output or args.pdb:
        if args.output:
            config_path = args.output
        else:
            config_path = f"targets/{args.pdb.upper()}/vina_config.txt"
        
        receptor_pdbqt = f"targets/{args.pdb.upper()}/receptor_cleaned.pdbqt" if args.pdb else "receptor/receptor_cleaned.pdbqt"
        config_content = box_to_vina_config(box, receptor_pdbqt)
        
        with open(config_path, 'w') as f:
            f.write(config_content)
        
        print(f"\n  Config saved: {config_path}")
    
    # Also print Vina command
    print(f"\n  Vina command:")
    print(f"  vina --receptor {receptor_pdbqt if args.pdb else 'RECEPTOR.pdbqt'} \\")
    print(f"       --center_x {box['center'][0]} --center_y {box['center'][1]} --center_z {box['center'][2]} \\")
    print(f"       --size_x {box['size'][0]} --size_y {box['size'][1]} --size_z {box['size'][2]}")
