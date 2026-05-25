#!/usr/bin/env python3
"""
Universal Receptor Preparation.
Fetches any PDB structure from RCSB, cleans it, and converts to PDBQT.

Usage:
    python receptor_prep.py --pdb 6LU7
    python receptor_prep.py --pdb 6LU7 --chain A
    python receptor_prep.py --pdb-file protein.pdb
    python receptor_prep.py --list              # List cached receptors

Output:
    targets/<pdb_id>/
        receptor_raw.pdb        # Original download
        receptor_cleaned.pdb    # Cleaned: no water, no heteroatoms, hydrogens added
        receptor_cleaned.pdbqt  # For AutoDock Vina
        target_info.yaml        # Metadata
"""

import os
import sys
import argparse
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import logging
import shutil
import subprocess
import platform
import yaml
import urllib.request

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────
# PDB fetching
# ──────────────────────────────────────────────

def fetch_pdb(pdb_id: str, output_path: str) -> bool:
    """
    Download PDB structure from RCSB.
    Returns True on success.
    """
    pdb_id = pdb_id.lower().strip()
    
    # RCSB download URLs
    urls = [
        f"https://files.rcsb.org/download/{pdb_id}.pdb",
        f"https://files.rcsb.org/download/{pdb_id}.pdb1",
    ]
    
    for url in urls:
        try:
            logger.info(f"Downloading {url}...")
            req = urllib.request.Request(url, headers={'User-Agent': 'DrugRepurposing/1.0'})
            with urllib.request.urlopen(req, timeout=30) as response:
                data = response.read().decode('utf-8')
            
            if len(data) > 100:
                with open(output_path, 'w') as f:
                    f.write(data)
                logger.info(f"Downloaded {pdb_id}.pdb ({len(data)} bytes)")
                return True
        except Exception as e:
            logger.debug(f"Failed {url}: {e}")
            continue
    
    logger.error(f"Could not download {pdb_id} from RCSB")
    return False


# ──────────────────────────────────────────────
# PDB cleaning with BioPython
# ──────────────────────────────────────────────

def clean_pdb(input_path: str, output_path: str, keep_chain: str = None) -> bool:
    """
    Clean a PDB file:
    1. Remove water molecules (HOH)
    2. Remove heteroatoms (keep only standard ATOM records)
    3. Keep only specified chain (if given)
    4. Keep only first model (NMR structures)
    5. Remove hydrogens (OpenBabel will re-add them properly)
    
    Uses BioPython if available, falls back to simple text parsing.
    """
    try:
        from Bio.PDB import PDBParser, PDBIO, Select
        return _clean_pdb_biopython(input_path, output_path, keep_chain)
    except ImportError:
        logger.warning("BioPython not available, using basic text-based cleaning")
        return _clean_pdb_basic(input_path, output_path, keep_chain)


def _clean_pdb_biopython(input_path: str, output_path: str, keep_chain: str = None) -> bool:
    """Clean using BioPython (preferred, handles edge cases)."""
    from Bio.PDB import PDBParser, PDBIO, Select
    
    class CleanSelect(Select):
        def accept_residue(self, residue):
            # Remove water
            if residue.get_resname() == 'HOH':
                return False
            # Remove heteroatoms
            if residue.id[0] != ' ':
                return False
            return True
        
        def accept_chain(self, chain):
            if keep_chain and chain.id != keep_chain:
                return False
            return True
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', input_path)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_path, CleanSelect())
        
        # Count residues
        residue_count = sum(1 for _ in structure.get_residues())
        logger.info(f"Cleaned: {residue_count} residues → {output_path}")
        return True
    except Exception as e:
        logger.error(f"BioPython cleaning failed: {e}")
        return _clean_pdb_basic(input_path, output_path, keep_chain)


def _clean_pdb_basic(input_path: str, output_path: str, keep_chain: str = None) -> bool:
    """Basic text-based PDB cleaning. Falls back when BioPython is unavailable."""
    with open(input_path) as f:
        lines = f.readlines()
    
    cleaned = []
    atom_count = 0
    seen_model = False
    
    for line in lines:
        # Skip ENDMDL / multiple models (only keep MODEL 1)
        if line.startswith('MODEL'):
            if seen_model:
                break
            seen_model = True
            continue
        if line.startswith('ENDMDL'):
            continue
        
        # Keep only ATOM records (skip HETATM = heteroatoms, water, ions)
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Skip water
            if 'HOH' in line[17:20]:
                continue
            
            # Chain filter
            if keep_chain and len(line) >= 22:
                chain_id = line[21]
                if chain_id != keep_chain:
                    continue
            
            # For HETATM, only keep if it looks like a co-factor (not water/ions)
            if line.startswith('HETATM'):
                resname = line[17:20].strip()
                if resname in ('HOH', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'SO4', 'PO4', 'GOL', 'EDO'):
                    continue
            
            atom_count += 1
            cleaned.append(line)
    
    # Add END
    cleaned.append('END\n')
    
    with open(output_path, 'w') as f:
        f.writelines(cleaned)
    
    logger.info(f"Cleaned: {atom_count} atoms → {output_path}")
    return True


# ──────────────────────────────────────────────
# PDB → PDBQT conversion
# ──────────────────────────────────────────────

def pdb_to_pdbqt(input_pdb: str, output_pdbqt: str) -> bool:
    """
    Convert cleaned PDB to PDBQT using OpenBabel.
    Adds hydrogens and Gasteiger charges.
    """
    from prep_ligands import get_obabel
    
    obabel = get_obabel()
    cmd = [
        obabel, input_pdb,
        "-O", output_pdbqt,
        "-xr",                          # Preserve residue names
        "-h",                           # Add hydrogens
        "--partialcharge", "gasteiger", # Required by Vina
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            logger.error(f"PDB→PDBQT failed: {result.stderr.strip()}")
            return False
        
        if os.path.exists(output_pdbqt) and os.path.getsize(output_pdbqt) > 100:
            logger.info(f"Converted: {output_pdbqt}")
            return True
        else:
            logger.error(f"PDBQT output is empty or missing")
            return False
    except Exception as e:
        logger.error(f"OpenBabel error: {e}")
        return False


# ──────────────────────────────────────────────
# Main preparation pipeline
# ──────────────────────────────────────────────

def prepare_receptor(pdb_id: str = None, pdb_file: str = None, 
                     chain: str = None, base_dir: str = "targets") -> dict:
    """
    Prepare a receptor for docking.
    
    Args:
        pdb_id: PDB ID to fetch (e.g., '6LU7')
        pdb_file: Local PDB file (alternative to pdb_id)
        chain: Specific chain to keep (e.g., 'A')
        base_dir: Output directory for targets
    
    Returns:
        Dictionary with receptor info (paths, metadata)
    """
    # ── Determine source ──
    if pdb_id:
        pdb_id = pdb_id.upper().strip()
        target_dir = os.path.join(base_dir, pdb_id)
        raw_pdb = os.path.join(target_dir, f"{pdb_id}_raw.pdb")
    elif pdb_file:
        pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
        target_dir = os.path.join(base_dir, pdb_id)
        raw_pdb = pdb_file  # Use directly
    else:
        raise ValueError("Must provide --pdb or --pdb-file")
    
    os.makedirs(target_dir, exist_ok=True)
    
    cleaned_pdb = os.path.join(target_dir, "receptor_cleaned.pdb")
    pdbqt_path = os.path.join(target_dir, "receptor_cleaned.pdbqt")
    info_path = os.path.join(target_dir, "target_info.yaml")
    
    info = {
        'pdb_id': pdb_id,
        'chain': chain,
        'prepared': None,
        'status': 'preparing',
    }
    
    # ── Step 1: Fetch (if needed) ──
    if not os.path.exists(raw_pdb) and pdb_file is None:
        if not fetch_pdb(pdb_id, raw_pdb):
            info['status'] = 'error'
            info['error'] = 'PDB download failed'
            return info
    
    logger.info(f"Raw PDB: {raw_pdb}")
    
    # ── Step 2: Clean ──
    if not clean_pdb(raw_pdb, cleaned_pdb, keep_chain=chain):
        info['status'] = 'error'
        info['error'] = 'PDB cleaning failed'
        return info
    
    # ── Step 3: Convert to PDBQT ──
    if not pdb_to_pdbqt(cleaned_pdb, pdbqt_path):
        # Fallback: if OpenBabel fails, try without hydrogen addition
        logger.warning("Retrying conversion without -h flag...")
        from prep_ligands import get_obabel
        obabel = get_obabel()
        cmd = [obabel, cleaned_pdb, "-O", pdbqt_path, "-xr", "--partialcharge", "gasteiger"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            info['status'] = 'error'
            info['error'] = 'PDB→PDBQT conversion failed'
            return info
    
    # ── Step 4: Save metadata ──
    from datetime import datetime
    info['status'] = 'ready'
    info['prepared'] = datetime.now().isoformat()
    info['cleaned_pdb'] = cleaned_pdb
    info['pdbqt'] = pdbqt_path
    info['residues'] = _count_residues(cleaned_pdb)
    
    with open(info_path, 'w') as f:
        yaml.dump(info, f, default_flow_style=False)
    
    logger.info(f"Receptor ready: {pdb_id}")
    logger.info(f"  Cleaned PDB: {cleaned_pdb}")
    logger.info(f"  PDBQT: {pdbqt_path}")
    logger.info(f"  Residues: {info.get('residues', '?')}")
    
    return info


def _count_residues(pdb_path: str) -> int:
    """Count unique residues in a PDB file."""
    residues = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21] if len(line) > 21 else ' '
                resi = line[22:26].strip()
                resn = line[17:20].strip()
                residues.add(f"{chain}:{resn}{resi}")
    return len(residues)


def list_receptors(base_dir: str = "targets"):
    """List all prepared receptors."""
    if not os.path.exists(base_dir):
        print("No targets prepared yet.")
        return
    
    print(f"\n  {'PDB ID':8s}  {'Status':10s}  {'Residues':8s}  {'Prepared'}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*20}")
    
    for pdb_id in sorted(os.listdir(base_dir)):
        info_path = os.path.join(base_dir, pdb_id, "target_info.yaml")
        if os.path.exists(info_path):
            with open(info_path) as f:
                info = yaml.safe_load(f)
            status = info.get('status', '?')
            residues = info.get('residues', '?')
            prepared = info.get('prepared', '?')[:10] if info.get('prepared') else '?'
            print(f"  {pdb_id:8s}  {status:10s}  {str(residues):8s}  {prepared}")
        else:
            print(f"  {pdb_id:8s}  no_info")


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Universal Receptor Preparation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python receptor_prep.py --pdb 6LU7              # SARS-CoV-2 Mpro
  python receptor_prep.py --pdb 6LU7 --chain A     # Specific chain
  python receptor_prep.py --pdb-file protein.pdb    # Local file
  python receptor_prep.py --list                    # Show cached
        """
    )
    parser.add_argument('--pdb', help='PDB ID to fetch (e.g., 6LU7)')
    parser.add_argument('--pdb-file', help='Local PDB file path')
    parser.add_argument('--chain', help='Specific chain to keep (e.g., A)')
    parser.add_argument('--output-dir', default='targets', help='Output directory')
    parser.add_argument('--list', action='store_true', help='List prepared receptors')
    
    args = parser.parse_args()
    
    if args.list:
        list_receptors(args.output_dir)
    elif args.pdb or args.pdb_file:
        info = prepare_receptor(
            pdb_id=args.pdb,
            pdb_file=args.pdb_file,
            chain=args.chain,
            base_dir=args.output_dir,
        )
        if info.get('status') == 'ready':
            print(f"\nReceptor ready for docking!")
            print(f"  Target directory: {os.path.dirname(info['pdbqt'])}")
            print(f"  PDBQT: {info['pdbqt']}")
        else:
            print(f"\nERROR: {info.get('error', 'Unknown error')}")
            sys.exit(1)
    else:
        parser.print_help()
