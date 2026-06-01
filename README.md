# AI Drug Repurposing: Multi-Target Virtual Screening Platform

> **⚠️ PROJECT OVERHAUL NOTICE (June 2026):** 
> This is a **new, significantly modified and updated version** of the platform. Previous versions contained methodological flaws (low docking exhaustiveness and lack of chemical property filtering) that led to over-optimistic or scientifically invalid results. 
> 
> **Key improvements in this version:**
> *   **Corrected Results:** Previous top hits (like Nilotinib_Var_17) have been re-evaluated and refuted using high-precision docking.
> *   **Scientific Filtering:** Integrated Lipinski's Rule of Five and QED scoring to ensure variants are drug-like and synthesizable.
> *   **Increased Precision:** Default AutoDock Vina exhaustiveness increased by 400% (from 8 to 32) to eliminate false positives.
> *   **Technical Stability:** Fixed all pathing, encoding, and execution bugs.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![RDKit](https://img.shields.io/badge/RDKit-2022.03+-green.svg)](https://www.rdkit.org/)

An automated, scientifically-validated pipeline for computational drug repurposing. This platform integrates Machine Learning pre-screening with high-precision Molecular Docking and Lead Optimization to identify potential therapeutic candidates for any protein target.

*   **ADMET Profiling:** Automated estimation of Absorption, Distribution, Metabolism, Excretion, and Toxicity.
*   **MD-Proxy Stability Analysis:** Evaluates docked pose strain energy to predict if a binding mode is physically realistic and stable.
*   **Toxicity Filtering:** Integrated PAINS and Brenk filters to identify and remove reactive or interfering compounds.

##  Scientific Validation

This platform enforces professional-grade pharmaceutical standards:
1.  **Lipinski's Rule of Five:** Filters for bioavailability.
2.  **QED (Quantitative Estimate of Drug-likeness):** Scores drug-likeness.
3.  **PAINS Filtering:** Removes Pan-Assay Interference Compounds.
4.  **Precision Docking:** (Exhaustiveness 32) minimizes false positives.
5.  **Pose Strain Analysis (MD Proxy):** Ensures the docked geometry is physically stable and not energetically "forced".


##  Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/ai-drug-repurposing.git
    cd ai-drug-repurposing
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **External Requirements:**
    *   [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)
    *   [OpenBabel](http://openbabel.org/) (for PDBQT conversion)

## 📖 Usage

### 1. Simple Screening (Standard)
```bash
python src/hermes_drug.py --pdb 6LU7
```

### 2. Deep Optimization (Scientific Mode)
```bash
python src/hermes_drug.py --pdb 6LU7 --depth deep
```

### 3. Running in the Cloud (Kaggle)
```bash
python src/cloud_runner.py --pdb 6LU7 --depth deep
```

##  Repository Structure

*   `src/`: Core Python modules for screening, docking, and optimization.
*   `data/`: Drug libraries, training sets, and ChEMBL exports.
*   `targets/`: Processed protein structures and docking configuration.
*   `docs/`: Detailed technical white paper and methodology guides.

##  License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

##  Citation

If you use this platform in your research, please cite:
```bibtex
@software{horvat_hermes_2026,
  author = {Turcsányi Zsombor Leonard and Hermes AI},
  title = {AI Drug Repurposing: Multi-Target Virtual Screening Platform},
  year = {2026},
  url = {https://github.com/horvatjanos/universal-drug-repurposing}
}
```
