# universal-drug-repurposing

Automated molecular docking and machine learning pipeline for drug screening.

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Experimental](https://img.shields.io/badge/Status-Experimental%20%2F%20In%20Development-orange)](https://github.com/zsomborturcsanyi7-lang/universal-drug-repurposing)

---

>  **IMPORTANT: Security & Architectural Warning**
>
> The current legacy architecture of this repository contains known **security and structural vulnerabilities** (including unvalidated CLI input execution, hardcoded file path dependencies, lack of sandboxed subprocess execution, and insecure model/data serialization).
> 
> **Do NOT use the current version in untrusted production environments.** 
>
>  **Upcoming Refactors & Architectural Improvements:**
> - **Input Sanitization & Security Hardening:** Strict validation for PDB/SDF inputs and safe execution bounds for external CLI tools (OpenBabel, Vina).
> - **Modular Architecture:** Full decoupling of data preprocessing, surrogate ML inference, and docking execution layers.
> - **Containerization:** Docker container support for isolated, reproducible execution environments.
> - **Enhanced Error Isolation:** Process-level error handling to prevent pipeline halts on corrupt ligand conformers.

---

## Overview & Purpose
`universal-drug-repurposing` provides an automated bioinformatic computational pipeline that uses AutoDock Vina alongside trained ML proxy models to screen candidate compound libraries against designated target receptor proteins.

## Key Features
- **Automated Receptor Prep:** PDB structure fetching and dynamic active-site grid box alignment.
- **ML Pre-Screening:** Fast ligand binding affinity pre-screening via machine learning surrogate models.
- **Parallel Docking:** Multi-core parallel AutoDock Vina execution.
- **Structured Reporting:** Automated score aggregation, Gibbs free energy ($\Delta G$) calculation, and candidate ranking reports.

## Tech Stack & Dependencies
- **Language**: Python 3.9+
- **Bioinformatics**: AutoDock Vina, RDKit, Biopython, OpenBabel
- **Machine Learning**: Scikit-Learn, PyTorch / NumPy, Pandas

## Project Structure
```text
universal-drug-repurposing/
├── run_pipeline.py      # Pipeline entry point
├── models/             # Trained surrogate ML models
├── scripts/            # Parsing, docking, and preprocessing utilities
├── requirements.txt    # Python dependencies
└── README.md
