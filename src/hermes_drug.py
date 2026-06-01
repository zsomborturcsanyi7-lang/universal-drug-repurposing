#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════╗
║       UNIVERSAL DRUG REPURPOSING PLATFORM                  ║
║       One command to screen any protein against any drug   ║
║       (v2.0: Scientifically Corrected & Validated)         ║
╚══════════════════════════════════════════════════════════════╝


Usage:
    python hermes_drug.py --pdb 6LU7
    python hermes_drug.py --pdb 6LU7 --deep
    python hermes_drug.py --pdb 6LU7 --cloud
    python hermes_drug.py --interactive
    python hermes_drug.py --list-targets
    python hermes_drug.py --pdb 6LU7 --drug-library custom_drugs.csv

Pipeline:
    1. receptor_prep   → Fetch & clean PDB structure
    2. auto_box        → Auto-detect binding pocket
    3. screen_ai       → ML proxy screening (fast)
    4. prep_ligands    → SMILES → 3D → PDBQT
    5. run_docking     → Parallel AutoDock Vina
    6. parse_results   → Rank & report
    7. optimize_leads  → Generate better variants (--deep mode)
"""

import os
import sys
import yaml
import argparse
import logging
import subprocess
import time
import json
from datetime import datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('hermes_drug')

# ──────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────

def load_config():
    """Load user config or use defaults."""
    config_path = "config.yaml"
    defaults = {
        'defaults': {
            'exhaustiveness': 8,
            'cpu': 0,
            'grid_spacing': 0.375,
            'num_modes': 9,
            'energy_range': 3,
        },
        'screening': {
            'quick_max_ligands': 10,
            'standard_max_ligands': 50,
            'deep_max_ligands': 200,
            'ml_proxy_model': 'vina_affinity_proxy.pkl',
            'similarity_threshold': 0.4,
            'variant_count': 50,
        },
    }
    
    if os.path.exists(config_path):
        with open(config_path, encoding='utf-8') as f:
            user_config = yaml.safe_load(f) or {}
        # Deep merge
        for section in defaults:
            if section in user_config:
                defaults[section].update(user_config[section])
    
    return defaults


# ──────────────────────────────────────────────
# Pipeline steps
# ──────────────────────────────────────────────

def step_receptor(pdb_id: str, config: dict):
    """Step 1: Prepare receptor."""
    logger.info("=" * 60)
    logger.info(f"  STEP 1: Preparing receptor {pdb_id}")
    logger.info("=" * 60)
    
    target_dir = f"targets/{pdb_id}"
    pdbqt_path = f"{target_dir}/receptor_cleaned.pdbqt"
    
    # Check cache
    if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 100:
        logger.info(f"  Receptor already prepared: {pdbqt_path}")
        return pdbqt_path
    
    # Prepare
    cmd = [sys.executable, "src/receptor_prep.py", "--pdb", pdb_id]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Receptor preparation failed:\n{result.stderr}")
        return None
    
    logger.info(result.stdout.strip())
    return pdbqt_path


def step_box(pdb_id: str, config: dict):
    """Step 2: Auto-detect docking box."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 2: Detecting binding pocket for {pdb_id}")
    logger.info("=" * 60)
    
    # Check if we already have a config
    config_path = f"targets/{pdb_id}/vina_config.txt"
    
    # Read existing box from config.yaml
    targets_config = config.get('targets', {})
    pdb_key = pdb_id.lower()
    
    target_info = targets_config.get(pdb_key, {})
    if target_info.get('center') and target_info.get('size'):
        box = {
            'center': target_info['center'],
            'size': target_info['size'],
            'method': 'cached',
        }
        logger.info(f"  Using cached box from config.yaml")
        logger.info(f"  Center: {box['center']}  Size: {box['size']}")
        return box, config_path
    
    # Auto-detect
    cmd = [
        sys.executable, "src/auto_box.py",
        "--pdb", pdb_id,
        "--output", config_path,
        "--method", "auto",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Box detection failed:\n{result.stderr}")
        return None, None
    
    # Parse output
    output = result.stdout
    center = None
    size = None
    method = 'auto'
    
    for line in output.split('\n'):
        if 'Center:' in line:
            try:
                center = eval(line.split(':', 1)[1].strip())
            except:
                pass
        elif 'Size:' in line:
            try:
                size = eval(line.split(':', 1)[1].strip())
            except:
                pass
        elif 'Method:' in line:
            method = line.split(':', 1)[1].strip()
    
    if not center or not size:
        logger.error("Could not parse box coordinates from auto_box output")
        return None, None
    
    box = {'center': center, 'size': size, 'method': method}
    logger.info(f"  Method: {method}")
    logger.info(f"  Center: {center}  Size: {size}")
    
    # Save to config
    config['targets'][pdb_key] = {
        'pdb_id': pdb_id,
        'center': center,
        'size': size,
        'prepared': datetime.now().isoformat(),
    }
    with open("config.yaml", 'w') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    
    return box, config_path


def step_screen(pdb_id: str, config: dict, depth: str = 'standard'):
    """Step 3: Screen drug library with AI proxy."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 3: AI screening ({depth} depth)")
    logger.info("=" * 60)
    
    screening = config['screening']
    max_ligands = screening.get(f'{depth}_max_ligands', 50)
    model_path = screening.get('ml_proxy_model', 'vina_affinity_proxy.pkl')
    
    if not os.path.exists(model_path):
        logger.warning(f"ML model not found: {model_path}")
        logger.info("Using similarity-based screening instead...")
        return step_screen_similarity(pdb_id, config, max_ligands)
    
    # Run AI screening
    if os.path.exists("src/ai_smart_search.py"):
        cmd = [sys.executable, "src/ai_smart_search.py"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        logger.info(result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
        
        # Return the AI hits file
        if os.path.exists("expert_ai_hits.csv"):
            return "expert_ai_hits.csv"
    
    # Fallback: use pre-existing hits
    for f in ["expert_ai_hits.csv", "ai_hits_for_docking.csv"]:
        if os.path.exists(f):
            logger.info(f"Using existing {f}")
            return f
    
    logger.error("No screening results available")
    return None


def step_screen_similarity(pdb_id: str, config: dict, max_ligands: int):
    """Fallback: Tanimoto similarity screening."""
    if os.path.exists("src/production_screening.py"):
        cmd = [sys.executable, "src/production_screening.py"]
        subprocess.run(cmd, capture_output=True, text=True)
    
    for f in ["expert_ai_hits.csv", "ai_hits_for_docking.csv", "global_repurposing_hits.csv"]:
        if os.path.exists(f):
            return f
    return None


def step_prepare_ligands(hits_csv: str, config: dict):
    """Step 4: Prepare ligand PDBQT files."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 4: Preparing ligands from {hits_csv}")
    logger.info("=" * 60)
    
    if not hits_csv or not os.path.exists(hits_csv):
        logger.warning("No hits to prepare")
        return
    
    cmd = [sys.executable, "src/prep_ligands.py", hits_csv]
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.info(result.stdout.strip())
    if result.stderr:
        logger.warning(result.stderr.strip())


def step_dock(pdb_id: str, config: dict):
    """Step 5: Run AutoDock Vina docking."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 5: Molecular docking for {pdb_id}")
    logger.info("=" * 60)
    
    if os.path.exists("src/run_parallel_docking.py"):
        cmd = [sys.executable, "src/run_parallel_docking.py", "--pdb", pdb_id]
    elif os.path.exists("src/run_docking_parallel.py"):
        cmd = [sys.executable, "src/run_docking_parallel.py", "--pdb", pdb_id]
    else:
        # Run a single test docking
        logger.warning("No parallel docking script found. Running single test dock...")
        return _run_single_dock(pdb_id, config)
    
    result = subprocess.run(cmd)
    
    if result.returncode != 0:
        logger.error("Docking failed")
        return
    
    logger.info("Docking complete!")


def _run_single_dock(pdb_id: str, config: dict):
    """Run a single docking as fallback."""
    cfg = config['defaults']
    targets = config.get('targets', {}).get(pdb_id.lower(), {})
    
    if not targets.get('center'):
        logger.error("No docking box coordinates found")
        return
    
    # Find a ligand
    ligand = None
    if os.path.exists("ligands"):
        for f in sorted(os.listdir("ligands")):
            if f.endswith('.pdbqt'):
                ligand = f"ligands/{f}"
                break
    
    if not ligand:
        logger.warning("No ligands found for docking")
        return
    
    center = targets['center']
    size = targets['size']
    receptor = f"targets/{pdb_id}/receptor_cleaned.pdbqt"
    
    vina_bin = os.environ.get('VINA_BINARY', 'vina')
    
    cmd = [
        vina_bin,
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(cfg.get('exhaustiveness', 8)),
        "--cpu", str(cfg.get('cpu', 0)),
        "--out", f"docking_results/{pdb_id}_docked.pdbqt",
    ]
    
    logger.info(f"  Ligand: {ligand}")
    logger.info(f"  Command: {' '.join(cmd[:6])}...")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.info(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)


def step_parse_results(pdb_id: str):
    """Step 6: Parse docking results."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 6: Parsing results")
    logger.info("=" * 60)
    
    if os.path.exists("src/parse_docking_results.py"):
        cmd = [sys.executable, "src/parse_docking_results.py"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        logger.info(result.stdout.strip())
    
    # Show top results
    results_file = "final_3d_docking_insights.csv"
    if os.path.exists(results_file):
        with open(results_file) as f:
            lines = f.readlines()
        
        if len(lines) > 1:
            logger.info(f"\n  TOP DOCKING RESULTS ({pdb_id}):")
            logger.info(f"  {'Rank':5s} {'Drug':30s} {'Affinity'}")
            logger.info(f"  {'-'*5} {'-'*30} {'-'*10}")
            
            for i, line in enumerate(lines[1:11]):  # Top 10
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    name = parts[0][:30]
                    aff = parts[1] if len(parts) > 1 else '?'
                    logger.info(f"  {i+1:<5d} {name:30s} {aff:>10s} kcal/mol")


def step_optimize(pdb_id: str, config: dict):
    """Step 7 (optional): Generate optimized variants."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 7: Lead optimization")
    logger.info("=" * 60)
    
    if os.path.exists("src/optimize_lead.py"):
        cmd = [sys.executable, "src/optimize_lead.py", "--pdb", pdb_id]
        subprocess.run(cmd)
    elif os.path.exists("src/mutate_lead.py"):
        cmd = [sys.executable, "src/mutate_lead.py"]
        subprocess.run(cmd)
    else:
        logger.info("  Lead optimization not available (no optimize_lead.py)")


# ──────────────────────────────────────────────
# Main orchestrator
# ──────────────────────────────────────────────

def run_pipeline(pdb_id: str, depth: str = 'standard', 
                  use_cloud: bool = False, 
                  drug_library: str = None,
                  skip_receptor: bool = False):
    """
    Run the full drug repurposing pipeline.
    
    Args:
        pdb_id: PDB identifier (e.g., '6LU7')
        depth: 'quick', 'standard', or 'deep'
        use_cloud: Push heavy workloads to Kaggle
        drug_library: Custom drug CSV (optional)
        skip_receptor: Skip receptor prep if already done
    """
    config = load_config()
    
    start_time = time.time()
    
    print()
    print("╔" + "═" * 58 + "╗")
    print(f"║  UNIVERSAL DRUG REPURPOSING PLATFORM" + " " * 22 + "║")
    print(f"║  Target: {pdb_id}" + " " * (47 - len(pdb_id)) + "║")
    print(f"║  Depth: {depth}" + " " * (48 - len(depth)) + "║")
    print(f"║  Cloud: {'ON' if use_cloud else 'OFF (local)'}" + " " * 30 + "║")
    print("╚" + "═" * 58 + "╝")
    print()
    
    # ── Step 1: Receptor ──
    if not skip_receptor:
        pdbqt_path = step_receptor(pdb_id, config)
        if not pdbqt_path:
            logger.error("Receptor preparation failed. Aborting.")
            return
    else:
        logger.info("STEP 1: SKIPPED (using cached receptor)")
    
    # ── Step 2: Box ──
    if 'targets' not in config:
        config['targets'] = {}
    box, config_path = step_box(pdb_id, config)
    
    # ── Adjust exhaustiveness based on depth ──
    if depth == 'quick':
        config['defaults']['exhaustiveness'] = 4
    elif depth == 'standard':
        config['defaults']['exhaustiveness'] = 16  # Increased from 8
    elif depth == 'deep':
        config['defaults']['exhaustiveness'] = 32  # High-quality validation

    # ── Step 3: AI Screening ──
    hits_csv = step_screen(pdb_id, config, depth)

    
    # ── Step 4: Prepare ligands ──
    step_prepare_ligands(hits_csv, config)
    
    # ── Step 5: Dock ──
    if use_cloud:
        step_cloud_dock(pdb_id, config)
    else:
        step_dock(pdb_id, config)
    
    # ── Step 6: Results ──
    step_parse_results(pdb_id)
    
    # ── Step 7: Optimize (deep mode only) ──
    if depth == 'deep':
        step_optimize(pdb_id, config)
    
    elapsed = time.time() - start_time
    
    print()
    print("╔" + "═" * 58 + "╗")
    print(f"║  PIPELINE COMPLETE" + " " * 39 + "║")
    print(f"║  Time: {elapsed/60:.1f} minutes" + " " * 35 + "║")
    target_dir = f"targets/{pdb_id}"
    print(f"║  Results: {target_dir}/" + " " * (35 - len(target_dir)) + "║")
    print("╚" + "═" * 58 + "╝")
    print()


def step_cloud_dock(pdb_id: str, config: dict):
    """Push docking to Kaggle cloud."""
    logger.info("\n" + "=" * 60)
    logger.info(f"  STEP 5 (CLOUD): Pushing to Kaggle")
    logger.info("=" * 60)
    
    if os.path.exists("cloud_runner.py"):
        cmd = [sys.executable, "cloud_runner.py", "--pdb", pdb_id]
        subprocess.run(cmd)
    else:
        logger.warning("cloud_runner.py not available. Running locally instead.")
        step_dock(pdb_id, config)


# ──────────────────────────────────────────────
# Interactive mode
# ──────────────────────────────────────────────

def interactive_mode():
    """Guided mode for non-programmers."""
    print()
    print("╔" + "═" * 58 + "╗")
    print("║     UNIVERSAL DRUG REPURPOSING PLATFORM         ║")
    print("║     Interactive Mode — No coding needed         ║")
    print("╚" + "═" * 58 + "╝")
    print()
    print("  This tool screens existing drugs against")
    print("  any protein target to find potential treatments.")
    print()
    
    # ── Mode selection ──
    print("  What would you like to do?")
    print()
    print("  1. Screen drugs against a new virus protein")
    print("  2. Optimize a known drug for better binding")
    print("  3. View previously prepared targets")
    print()
    
    choice = input("  Choice [1]: ").strip() or "1"
    
    if choice == "3":
        if os.path.exists("targets"):
            os.system(f"{sys.executable} receptor_prep.py --list")
        return
    
    # ── PDB input ──
    print()
    print("  Enter the PDB ID of the target protein.")
    print("  Examples: 6LU7 (SARS-CoV-2 Mpro), 6M0J (SARS-CoV-2 spike)")
    print()
    pdb_id = input("  PDB ID: ").strip().upper()
    
    if not pdb_id:
        print("  No PDB ID provided. Exiting.")
        return
    
    # ── Depth ──
    print()
    print("  Screening depth:")
    print("  1. Quick  — FDA-approved drugs only (~2 min)")
    print("  2. Normal — FDA + experimental drugs (~10 min)")
    print("  3. Deep   — Full library + AI optimization (~30 min)")
    print()
    depth_choice = input("  Depth [2]: ").strip() or "2"
    
    depth_map = {"1": "quick", "2": "standard", "3": "deep"}
    depth = depth_map.get(depth_choice, "standard")
    
    # ── Cloud ──
    print()
    if depth == "deep":
        print("  Deep screening is recommended to run on cloud.")
        cloud_choice = input("  Use Kaggle cloud? [Y/n]: ").strip().lower()
        use_cloud = cloud_choice != 'n'
    else:
        use_cloud = False
    
    # ── Run ──
    print()
    print(f"  Starting pipeline for {pdb_id} ({depth} mode)...")
    print()
    
    run_pipeline(pdb_id, depth=depth, use_cloud=use_cloud)


# ──────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Universal Drug Repurposing Platform",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python hermes_drug.py --pdb 6LU7              # Screen SARS-CoV-2 Mpro
  python hermes_drug.py --pdb 6LU7 --deep       # Full optimization
  python hermes_drug.py --pdb 6LU7 --cloud      # Use Kaggle acceleration
  python hermes_drug.py --interactive           # Guided mode for non-programmers
  python hermes_drug.py --list-targets          # Show prepared receptors
        """
    )
    
    parser.add_argument('--pdb', help='PDB ID of target protein (e.g., 6LU7)')
    parser.add_argument('--depth', default='standard',
                       choices=['quick', 'standard', 'deep'],
                       help='Screening depth (default: standard)')
    parser.add_argument('--cloud', action='store_true',
                       help='Use Kaggle cloud for docking')
    parser.add_argument('--drug-library',
                       help='Custom CSV file with drugs to screen')
    parser.add_argument('--skip-receptor', action='store_true',
                       help='Skip receptor preparation (use cached)')
    parser.add_argument('--interactive', '-i', action='store_true',
                       help='Guided interactive mode')
    parser.add_argument('--list-targets', action='store_true',
                       help='List prepared receptor targets')
    
    args = parser.parse_args()
    
    if args.list_targets:
        os.system(f"{sys.executable} receptor_prep.py --list")
    elif args.interactive:
        interactive_mode()
    elif args.pdb:
        run_pipeline(
            pdb_id=args.pdb,
            depth=args.depth,
            use_cloud=args.cloud,
            drug_library=args.drug_library,
            skip_receptor=args.skip_receptor,
        )
    else:
        parser.print_help()
        print("\n  Try interactive mode: python hermes_drug.py --interactive")
