#!/usr/bin/env python3
"""
Parallel Docking Runner.
Executes AutoDock Vina in parallel across multiple ligands.
"""

import os
import sys
import subprocess
import concurrent.futures
import logging
import argparse
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_vina(ligand_path, receptor_path, config_path, output_dir):
    ligand_name = Path(ligand_path).stem
    output_pdbqt = os.path.join(output_dir, f"{ligand_name}_out.pdbqt")
    log_file = os.path.join(output_dir, f"{ligand_name}_log.txt")
    
    # Check if already done
    if os.path.exists(output_pdbqt) and os.path.getsize(output_pdbqt) > 500:
        return f"Skipped {ligand_name} (already exists)"
    
    # Build command
    # Try to find vina in path, fallback to common windows names
    vina_bin = "vina"
    if os.name == 'nt':
        # Check for common windows names
        for v in ["vina.exe", "vina_1.2.5_win.exe", "vina_1.2.3_win.exe", "vina_1.2.7_win.exe"]:
            if subprocess.run(["where", v], capture_output=True).returncode == 0:
                vina_bin = v
                break

    cmd = [
        vina_bin,
        "--config", config_path,
        "--ligand", ligand_path,
        "--out", output_pdbqt
    ]
    
    try:
        with open(log_file, "w") as f:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            f.write(result.stdout)
            f.write(result.stderr)
            
        if result.returncode == 0:
            return f"Completed {ligand_name}"
        else:
            return f"Failed {ligand_name} (check log)"
    except Exception as e:
        return f"Error docking {ligand_name}: {e}"

def main():
    parser = argparse.ArgumentParser(description="Parallel Vina Runner")
    parser.add_argument('--pdb', help='PDB ID for context')
    parser.add_argument('--ligands', default='ligands', help='Directory with PDBQT ligands')
    parser.add_argument('--results', default='docking_results', help='Output directory')
    parser.add_argument('--cpus', type=int, default=os.cpu_count(), help='Number of parallel tasks')
    
    args = parser.parse_args()
    
    # Paths
    if args.pdb:
        pdb_id = args.pdb.upper()
        config_path = f"targets/{pdb_id}/vina_config.txt"
        receptor_path = f"targets/{pdb_id}/receptor_cleaned.pdbqt"
    else:
        # Fallback to current dir if no PDB specified
        config_path = "vina_config.txt"
        receptor_path = "receptor_cleaned.pdbqt"
        
    if not os.path.exists(config_path):
        logger.error(f"Config file not found: {config_path}")
        # Try to find ANY vina_config.txt in targets/
        configs = list(Path("targets").glob("**/vina_config.txt"))
        if configs:
            config_path = str(configs[0])
            logger.info(f"Using fallback config: {config_path}")
        else:
            sys.exit(1)
            
    if not os.path.exists(args.results):
        os.makedirs(args.results)
        
    # Gather ligands
    ligands = [os.path.join(args.ligands, f) for f in os.listdir(args.ligands) if f.endswith(".pdbqt")]
    
    if not ligands:
        logger.error(f"No ligands found in {args.ligands}")
        sys.exit(1)
        
    logger.info(f"Starting parallel docking of {len(ligands)} ligands using {args.cpus} CPUs...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.cpus // 2)) as executor:
        # Vina itself is multi-threaded, so we don't want too many parallel instances
        # unless we set --cpu 1 in Vina. Assuming default Vina behavior.
        futures = [executor.submit(run_vina, l, receptor_path, config_path, args.results) for l in ligands]
        
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            logger.info(f"[{i}/{len(ligands)}] {future.result()}")

    logger.info("All docking tasks finished.")

if __name__ == "__main__":
    main()
