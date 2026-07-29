# universal-drug-repurposing

Automated molecular docking and machine learning pipeline for drug screening.

## Overview & Purpose
universal-drug-repurposing provides an automated bioinformatic computational pipeline that uses AutoDock Vina alongside trained ML proxy models to screen candidate compound libraries against designated target receptor proteins.

## Key Features
- Automated receptor PDB structure fetching and grid box alignment.
- Fast ligand affinity pre-screening via machine learning surrogate models.
- Parallel AutoDock Vina docking execution.
- Automated score aggregation and candidate ranking reports.

## Tech Stack & Dependencies
- **Language**: Python 3.9+
- **Bioinformatics**: AutoDock Vina, RDKit, Biopython
- **Machine Learning**: Scikit-Learn, NumPy, Pandas

## Project Structure
```text
universal-drug-repurposing/
├── run_pipeline.py
├── models/
├── scripts/
├── requirements.txt
└── README.md
```

## Installation & Setup

### Prerequisites
- Python 3.9+
- AutoDock Vina 1.2+ installed

### Steps
```bash
git clone https://github.com/zsomborturcsanyi7-lang/universal-drug-repurposing.git
cd universal-drug-repurposing

pip install -r requirements.txt
```

## Usage Examples
```bash
python run_pipeline.py --receptor receptor.pdbqt --ligands ligands.sdf --out results/
```

## Status & License
Status: In Silico Research Pipeline / Experimental.
License: MIT
