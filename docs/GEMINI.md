# AI Drug Repurposing Pipeline

## Project Overview
This project provides an automated, Python-based pipeline for virtual drug repurposing. It leverages open-source chemoinformatics tools to perform two levels of screening:
1. **Structural Similarity Screening:** Using RDKit to compute Morgan Fingerprints (ECFP4) and Tanimoto similarity against a reference drug.
2. **Molecular Docking:** Utilizing AutoDock Vina to predict binding affinity (kcal/mol) of drugs against specific protein targets (e.g., SARS-CoV-2 Mpro).

## Architecture
- **Data Acquisition:** Integrates with the ChEMBL API to fetch FDA-approved small-molecule drugs.
- **Structural Screening:** Uses `production_screening.py` to identify structural analogs.
- **Docking Workflow:** 
    - `receptor_prep.py`: Fetches and cleans PDB structures (e.g., 6LU7).
    - `prep_ligands.py`: Prepares ligand libraries (SMILES -> 3D -> PDBQT) for docking.
    - `generate_docking_commands.py`: Orchestrates the AutoDock Vina execution via batch scripting.
- **Analysis:** `parse_docking_results.py` parses Vina output logs to extract binding affinities and generate a ranked report.

## Building and Running
### Dependencies
Ensure the following are installed:
- Python 3.x
- RDKit (`pip install rdkit`)
- Biopython (`pip install biopython`)
- pandas (`pip install pandas`)
- chembl_webresource_client (`pip install chembl_webresource_client`)
- Open Babel (Must be in system PATH)
- AutoDock Vina (Binary provided as `vina_1.2.7_win.exe`)

### Workflow
1. **Prepare Receptor:** `python receptor_prep.py`
2. **Screen by Similarity:** `python production_screening.py`
3. **Prepare Ligands for Docking:** `python prep_ligands.py`
4. **Docking:**
   - Generate configuration: `python generate_docking_commands.py`
   - Run docking (Parallel): `python run_parallel_docking.py`
   - *Legacy (Sequential):* `.\run_batch_docking.bat`
5. **Analyze Results:** `python parse_docking_results.py`

## Development Conventions
- **Naming:** Uses clear, descriptive filenames for pipeline stages.
- **Logging:** Employs Python `logging` for real-time tracking of pipeline progress.
- **Robustness:** Includes regex-based parsing and strict PDBQT formatting for Vina compliance.
