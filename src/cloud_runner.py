#!/usr/bin/env python3
"""
Kaggle Cloud Runner.
Automates dataset push → kernel execution → result pull.

Usage:
    python cloud_runner.py --pdb 6LU7
    python cloud_runner.py --pdb 6LU7 --deep
    python cloud_runner.py --status         # Check current kernel status
    python cloud_runner.py --pull           # Pull latest results
"""

import os
import sys
import json
import time
import shutil
import logging
import argparse
import subprocess
import tempfile

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('cloud_runner')

# ──────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────

DEFAULT_USER = "horvatjanos"
DEFAULT_DATASET = "ai-drug-workspace"
DEFAULT_KERNEL = "ai-drug-docking-execution"

KAGGLE_CLI = "kaggle"


# ──────────────────────────────────────────────
# Kaggle helpers
# ──────────────────────────────────────────────

def kaggle_status(kernel_id: str) -> str:
    """Get kernel status."""
    result = subprocess.run(
        [KAGGLE_CLI, "kernels", "status", kernel_id],
        capture_output=True, text=True, timeout=30
    )
    return result.stdout.strip()


def kaggle_push(kernel_dir: str) -> bool:
    """Push kernel to Kaggle."""
    result = subprocess.run(
        [KAGGLE_CLI, "kernels", "push", "-p", kernel_dir],
        capture_output=True, text=True, timeout=30
    )
    return result.returncode == 0


def kaggle_pull(kernel_id: str, output_dir: str) -> bool:
    """Pull kernel output."""
    result = subprocess.run(
        [KAGGLE_CLI, "kernels", "output", kernel_id, "-p", output_dir, "--force"],
        capture_output=True, text=True, timeout=300
    )
    return result.returncode == 0


# ──────────────────────────────────────────────
# Remote entrypoint builder
# ──────────────────────────────────────────────

def build_remote_entry(pdb_id: str, depth: str = "standard") -> str:
    """
    Create a remote_entry.py that runs the full pipeline on Kaggle.
    This is the file that actually executes on Kaggle's servers.
    """
    
    script = f'''#!/usr/bin/env python3
"""Remote docking pipeline — runs on Kaggle."""

import os, sys, subprocess, shutil, zipfile

print("=" * 60)
print(f"  DRUG REPURPOSING — Target: {pdb_id}")
print(f"  Depth: {depth}")
print("=" * 60)

# ── 1. Find & copy dataset ──
INPUT = "/kaggle/input/ai-drug-workspace"
if not os.path.exists(INPUT):
    # Fallback search
    for r, d, f in os.walk("/kaggle/input"):
        if "config.yaml" in f or "targets.zip" in f: 
            INPUT = r
            break

print(f"Using input from: {{INPUT}}")

# Copy scripts and csv files
for item in os.listdir(INPUT):
    src = os.path.join(INPUT, item)
    if item.endswith(".py") or item.endswith(".csv") or item.endswith(".yaml") or item.endswith(".txt"):
        if os.path.isfile(src): shutil.copy2(src, '.')
        print(f"Copied {{item}}")

# Extract zipped folders if present
for folder in ["ligands", "targets"]:
    zip_path = os.path.join(INPUT, f"{{folder}}.zip")
    if os.path.exists(zip_path):
        print(f"Extracting {{zip_path}}...")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall('.')
    elif os.path.isdir(os.path.join(INPUT, folder)):
        print(f"Copying {{folder}} directory...")
        shutil.copytree(os.path.join(INPUT, folder), folder, dirs_exist_ok=True)

# ── 2. Install deps ──
print("Installing dependencies...")
subprocess.run(["pip", "install", "-q", "rdkit", "biopython", "tqdm", "pandas", "pyyaml"], check=True)
subprocess.run(["apt-get", "update", "-qq"], check=True)
subprocess.run(["apt-get", "install", "-y", "-qq", "openbabel"], check=True)

# ── 3. Setup Vina ──
print("Setting up AutoDock Vina...")
subprocess.run(["wget", "-q",
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64",
    "-O", "vina"], check=True)
subprocess.run(["chmod", "+x", "vina"], check=True)
os.environ["VINA_BINARY"] = os.path.abspath("vina")
# Create a symlink in path for convenience
subprocess.run(["ln", "-s", os.path.abspath("vina"), "/usr/local/bin/vina"], check=True)

# ── 4. Run docking ──
print("Starting docking process...")
# We already have the ligands prepared in the dataset.
# Run the parallel docking script directly for the target.
subprocess.run([sys.executable, "run_parallel_docking.py", "--pdb", "{pdb_id}"], check=True)

# ── 5. Parse results ──
print("Parsing results...")
subprocess.run([sys.executable, "parse_docking_results.py"], check=True)

# ── 6. Summary ──
results = "final_3d_docking_insights.csv"
if os.path.exists(results):
    with open(results) as f: lines = f.readlines()
    with open("RESULT.txt", "w") as f:
        f.write(f"DOCKING_COMPLETE: target={pdb_id}\\n")
        f.write(f"LIGANDS_DOCKED: {{len(lines)-1 if lines else 0}}\\n")
        if len(lines) > 1:
            top = lines[1].strip().split(',')
            f.write(f"BEST_DRUG: {{top[0]}}\\n")
            f.write(f"BEST_AFFINITY: {{top[1]}} kcal/mol\\n")
    print(f"\\nDone. {{len(lines)-1 if lines else 0}} ligands docked.")
else:
    print("Error: results file not found!")
    with open("RESULT.txt", "w") as f:
        f.write("DOCKING_FAILED\\n")

# Clean up Vina to keep download small
if os.path.exists("vina"): os.remove("vina")
'''

    return script


# ──────────────────────────────────────────────
# Main runner
# ──────────────────────────────────────────────

def run_cloud(pdb_id: str, depth: str = "standard", 
               username: str = DEFAULT_USER,
               kernel_slug: str = DEFAULT_KERNEL):
    """Run the full pipeline on Kaggle cloud."""
    
    kernel_id = f"{username}/{kernel_slug}"
    
    # ── Create upload directory ──
    upload_dir = tempfile.mkdtemp(prefix="kaggle_push_")
    logger.info(f"Upload directory: {upload_dir}")
    
    # Write remote_entry.py
    remote_entry = build_remote_entry(pdb_id, depth)
    with open(os.path.join(upload_dir, "remote_entry.py"), "w", encoding='utf-8') as f:
        f.write(remote_entry)
    
    # Write kernel-metadata.json
    kernel_meta = {
        "id": kernel_id,
        "title": f"Drug Repurposing: {pdb_id}",
        "code_file": "remote_entry.py",
        "language": "python",
        "kernel_type": "script",
        "is_private": "true",
        "enable_gpu": "false",
        "enable_internet": "true",
        "dataset_sources": [f"{username}/{DEFAULT_DATASET}"],
        "competition_sources": [],
        "kernel_sources": [],
    }
    with open(os.path.join(upload_dir, "kernel-metadata.json"), "w", encoding='utf-8') as f:
        json.dump(kernel_meta, f, indent=2)
    
    # ── Push kernel ──
    logger.info(f"Pushing kernel to {kernel_id}...")
    
    result = subprocess.run(
        [KAGGLE_CLI, "kernels", "push", "-p", upload_dir],
        capture_output=True, text=True, timeout=30
    )
    
    # Cleanup
    shutil.rmtree(upload_dir, ignore_errors=True)
    
    if result.returncode != 0:
        logger.error(f"Kernel push failed: {result.stderr}")
        return
    
    logger.info(result.stdout.strip())
    
    # ── Poll for completion ──
    logger.info("Waiting for kernel to complete...")
    
    max_wait = 1800  # 30 minutes
    poll_interval = 30
    waited = 0
    
    while waited < max_wait:
        time.sleep(poll_interval)
        waited += poll_interval
        
        status = kaggle_status(kernel_id)
        logger.info(f"  Status ({waited}s): {status}")
        
        if "COMPLETE" in status:
            logger.info("Kernel completed!")
            break
        elif "ERROR" in status:
            logger.error("Kernel failed!")
            break
    
    # ── Pull results ──
    output_dir = f"cloud_results/{pdb_id}"
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"Downloading results to {output_dir}...")
    
    if kaggle_pull(kernel_id, output_dir):
        logger.info(f"Results downloaded to {output_dir}")
        
        # Check result
        result_file = os.path.join(output_dir, "RESULT.txt")
        if os.path.exists(result_file):
            with open(result_file) as f:
                logger.info(f.read())
    else:
        logger.error("Failed to download results")


def check_status(username: str = DEFAULT_USER, kernel_slug: str = DEFAULT_KERNEL):
    """Check current kernel status."""
    kernel_id = f"{username}/{kernel_slug}"
    status = kaggle_status(kernel_id)
    print(f"Kernel: {kernel_id}")
    print(f"Status: {status}")


def pull_results(pdb_id: str = None, username: str = DEFAULT_USER, 
                  kernel_slug: str = DEFAULT_KERNEL):
    """Pull latest results from Kaggle."""
    kernel_id = f"{username}/{kernel_slug}"
    output_dir = f"cloud_results/{pdb_id or 'latest'}"
    os.makedirs(output_dir, exist_ok=True)
    
    if kaggle_pull(kernel_id, output_dir):
        print(f"Results downloaded to {output_dir}")
    else:
        print("Failed to download results")


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Kaggle Cloud Runner")
    parser.add_argument('--pdb', help='PDB ID')
    parser.add_argument('--depth', default='standard', 
                       choices=['quick', 'standard', 'deep'])
    parser.add_argument('--status', action='store_true', help='Check kernel status')
    parser.add_argument('--pull', action='store_true', help='Pull latest results')
    
    args = parser.parse_args()
    
    if args.status:
        check_status()
    elif args.pull:
        pull_results(args.pdb)
    elif args.pdb:
        run_cloud(args.pdb, args.depth)
    else:
        parser.print_help()
