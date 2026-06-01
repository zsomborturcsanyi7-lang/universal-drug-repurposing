#!/usr/bin/env python3
"""
Docking Result Parser.
Reads AutoDock Vina output logs and extracts binding affinities.

Usage:
    python parse_docking_results.py [docking_results_dir]
    
Output:
    final_3d_docking_insights.csv — ranked list of all docked ligands
"""

import os
import re
import sys
import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_vina_log(log_path: str) -> dict:
    """
    Parse a single AutoDock Vina output log.
    
    Returns:
        {'ligand': name, 'affinity': kcal/mol, 'modes': count, 'status': 'ok'/'failed'}
    """
    if not os.path.exists(log_path):
        return {'status': 'missing', 'ligand': Path(log_path).stem.replace('_log', '')}
    
    try:
        with open(log_path) as f:
            content = f.read()
    except:
        return {'status': 'error', 'ligand': Path(log_path).stem.replace('_log', '')}
    
    # Check for errors
    if 'Parse error' in content or 'ERROR' in content.upper():
        return {'status': 'error', 'ligand': Path(log_path).stem.replace('_log', '')}
    
    # Extract affinity values
    # Vina output format:
    #   1       -8.295      0.000      0.000
    #   2       -7.809      2.164      2.619
    affinities = []
    
    for line in content.split('\n'):
        line = line.strip()
        
        # Match the docking results table
        # Pattern: <number> <affinity> <rmsd_lb> <rmsd_ub>
        if line and line[0].isdigit():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    affinities.append(affinity)
                except ValueError:
                    continue
    
    if not affinities:
        # Try broader pattern matching
        affinity_pattern = re.findall(r'-\d+\.\d+', content)
        affinities = [float(a) for a in affinity_pattern if float(a) < 0]
    
    if not affinities:
        return {'status': 'no_results', 'ligand': Path(log_path).stem.replace('_log', '')}
    
    best_affinity = min(affinities)  # Most negative = best
    
    return {
        'status': 'ok',
        'ligand': Path(log_path).stem.replace('_log', ''),
        'affinity': best_affinity,
        'modes': len(affinities),
        'all_affinities': sorted(affinities),
    }


def parse_all_logs(log_dir: str = "docking_results") -> pd.DataFrame:
    """
    Parse all Vina output logs in a directory.
    
    Returns:
        DataFrame ranked by binding affinity (best first)
    """
    if not os.path.exists(log_dir):
        logger.error(f"Directory not found: {log_dir}")
        return pd.DataFrame()
    
    log_files = [f for f in os.listdir(log_dir) if f.endswith('_log.txt') or f.endswith('.log')]
    
    if not log_files:
        # Also look for _log.txt files in subdirectories
        for root, dirs, files in os.walk(log_dir):
            for f in files:
                if f.endswith('_log.txt'):
                    log_files.append(os.path.join(root, f))
            if log_files:
                break
    
    if not log_files:
        logger.warning(f"No docking log files found in {log_dir}")
        return pd.DataFrame()
    
    logger.info(f"Parsing {len(log_files)} docking logs...")
    
    results = []
    for log_file in log_files:
        log_path = os.path.join(log_dir, log_file) if not os.path.dirname(log_file) else log_file
        result = parse_vina_log(log_path)
        
        if result['status'] == 'ok':
            results.append({
                'Ligand_Name': result['ligand'],
                'Binding_Affinity_kcal_mol': result['affinity'],
                'Num_Modes': result['modes'],
                'Best_Affinity': result['all_affinities'][0] if result['all_affinities'] else None,
            })
        elif result['status'] == 'error':
            logger.warning(f"  Failed: {result['ligand']}")
        elif result['status'] == 'no_results':
            logger.debug(f"  No affinity found: {result['ligand']}")
    
    if not results:
        logger.warning("No valid docking results found")
        return pd.DataFrame()
    
    # Build DataFrame and rank
    df = pd.DataFrame(results)
    df = df.sort_values('Binding_Affinity_kcal_mol', ascending=True)  # Best first
    
    # ── Merge with Molecular Metadata (Scientific Scores) ──
    metadata_files = ["optimized_variants.csv", "expert_ai_hits.csv", "full_fda_library.csv", "ai_hits_for_docking.csv"]
    for meta_file in metadata_files:
        if os.path.exists(meta_file):
            try:
                meta_df = pd.read_csv(meta_file)
                # Find the name column
                name_col = None
                for col in ['Name', 'Ligand_Name', 'drug_name']:
                    if col in meta_df.columns:
                        name_col = col
                        break
                
                if name_col:
                    # Clean names for matching (remove .pdbqt extension if present in metadata)
                    meta_df[name_col] = meta_df[name_col].astype(str).str.replace('.pdbqt', '', regex=False)
                    
                    df = df.merge(meta_df, left_on='Ligand_Name', right_on=name_col, how='left', suffixes=('', '_meta'))
                    
                    # Drop redundant columns
                    if name_col != 'Ligand_Name':
                        df = df.drop(columns=[name_col])
                    
                    logger.info(f"  Enriched results with metadata from {meta_file}")
                    break
            except Exception as e:
                logger.warning(f"  Could not merge metadata from {meta_file}: {e}")

    df = df.reset_index(drop=True)
    
    # Save
    output_file = "final_3d_docking_insights.csv"
    df.to_csv(output_file, index=False)
    
    logger.info(f"\n{'='*50}")
    logger.info(f"  DOCKING RESULTS SUMMARY")
    logger.info(f"{'='*50}")
    logger.info(f"  Total docked:  {len(df)}")
    logger.info(f"  Best affinity: {df.iloc[0]['Binding_Affinity_kcal_mol']:.3f} kcal/mol ({df.iloc[0]['Ligand_Name']})")
    logger.info(f"  Saved to:      {output_file}")
    logger.info(f"{'='*50}")
    
    # Show top 10
    print(f"\n  TOP 10 CANDIDATES:")
    print(f"  {'Rank':5s} {'Drug':35s} {'Affinity (kcal/mol)'}")
    print(f"  {'-'*5} {'-'*35} {'-'*18}")
    
    for i, row in df.head(10).iterrows():
        print(f"  {i+1:<5d} {row['Ligand_Name']:35s} {row['Binding_Affinity_kcal_mol']:>10.3f}")
    
    return df


if __name__ == "__main__":
    if len(sys.argv) > 1:
        log_dir = sys.argv[1]
    else:
        log_dir = "docking_results"
    
    parse_all_logs(log_dir)
