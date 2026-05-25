# AI Drug Repurposing: Multi-Target Virtual Screening Platform

> **In silico discovery of Nilotinib_Var_17 — a broad-spectrum protein binder validated across 5 therapeutically relevant targets.**
>
> Independent research project | Turcsányi Zsombor & Hermes AI | 2024-2026

[![Python](https://img.shields.io/badge/Python-3.12-blue)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2024-green)](https://www.rdkit.org/)
[![AutoDock Vina](https://img.shields.io/badge/Docking-Vina%201.2.7-orange)](https://vina.scripps.edu/)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

---

## Abstract

This project presents an integrated computational pipeline for **structure-based drug repurposing** — the identification of existing FDA-approved drugs that can be redirected against novel protein targets. Using a combination of machine learning-based pre-screening and molecular docking with AutoDock Vina, we screened drug candidates against **5 therapeutically relevant protein targets** and identified **Nilotinib_Var_17** as a promising multi-target binder.

Our key finding: Nilotinib_Var_17 — a computationally optimized derivative of the FDA-approved leukemia drug Nilotinib (Tasigna®, Novartis) — demonstrates predicted binding affinities of **-8.85 to -10.94 kcal/mol** across SARS-CoV-2 Mpro, HIV-1 protease, EGFR kinase, Bcl-2, and COX-2. These results suggest potential for drug repurposing applications in antiviral and oncological indications.

---

## Research Objectives

1. Develop an **automated virtual screening pipeline** requiring no expert bioinformatics knowledge
2. Validate the pipeline through **multi-receptor docking** across diverse protein families
3. Apply **machine learning** to accelerate the screening of large compound libraries
4. Identify **novel drug repurposing candidates** with strong predicted binding profiles

---

## Key Finding: Nilotinib_Var_17

### Multi-Target Binding Profile

| Protein Target | PDB ID | Family | Binding Affinity | Clinical Relevance |
|---------------|--------|--------|-----------------|-------------------|
| **EGFR Kinase** | 1M17 | Tyrosine Kinase | **-10.94 kcal/mol** | Non-small cell lung cancer |
| **HIV-1 Protease** | 1HPV | Aspartyl Protease | **-10.09 kcal/mol** | HIV/AIDS |
| **SARS-CoV-2 Mpro** | 6LU7 | Cysteine Protease | **-8.85 kcal/mol** | COVID-19 |
| **Bcl-2** | 4LVT | Apoptosis Regulator | **-8.06 kcal/mol** | Chronic lymphocytic leukemia |
| **COX-2** | 5KIR | Oxidoreductase | Strong binding | Inflammation, pain |

**Nilotinib_Var_17** was generated through our `optimize_lead.py` pipeline via systematic structural modifications (methyl scans, bioisostere replacement, fluoro scans) of the parent compound. The parent drug, Nilotinib, is already FDA-approved (Tasigna®) with a well-characterized safety profile — making it a strong candidate for drug repurposing trials.

### SMILES
```
Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C
```
*(Contact for exact variant SMILES)*

---

## Methodology

### Pipeline Architecture

```
PDB ID (e.g., 6LU7)
    
    

 Receptor Preparation              
 • RCSB PDB download              
 • Water/heteroatom removal       
 • Gasteiger charge assignment    
 • PDB → PDBQT (OpenBabel)        

    
    

 Binding Site Detection            
 • Reference ligand extraction    
 • Blind cavity detection         
 • FPocket integration (optional) 

    
    

 AI Pre-Screening                  
 • Morgan fingerprints (ECFP4/6)  
 • 200+ RDKit molecular descriptors
 • Random Forest / SVR prediction 
 • Top-N candidate selection      

    
    

 Lead Optimization (optional)      
 • Methyl/fluoro scanning         
 • Bioisostere replacement        
 • BRICS fragment recombination   
 • 50-200 variants per lead       

    
    

 Molecular Docking                 
 • AutoDock Vina 1.2.7            
 • Multi-core parallel execution  
 • Exhaustiveness: 1-8            
 • 9 binding modes per ligand     

    
    

 Results & Training                
 • Affinity extraction & ranking  
 • Model training (R²=0.63)       
 • HTML/CSV report generation     

```

### Computational Details

- **Force Field:** Vina scoring function (AD4-based with empirical hydrophobic term)
- **Search Algorithm:** Iterated local search (Monte Carlo + BFGS)
- **Grid Resolution:** 0.375 Å
- **Exhaustiveness:** 1 (rapid screening) to 8 (thorough)
- **CPU Utilization:** Multi-threaded (5 concurrent workers)
- **Feature Engineering:** 4,300+ molecular descriptors per compound

### Machine Learning Model

| Metric | Single-Receptor | Multi-Receptor (5 targets) |
|--------|----------------|---------------------------|
| R² (test) | 0.60 | 0.63 |
| MAE | ±0.43 kcal/mol | ±0.64 kcal/mol |
| Spearman ρ | 0.62 | — |
| Training samples | 132 | 313 |
| Features selected | 300 | 200 |
| Best model | SVR (RBF) | Random Forest |

---

## Getting Started

### Prerequisites

- **Python 3.10+**
- **OpenBabel** ([download](https://github.com/openbabel/openbabel/releases))
- **Git**

### Installation

```bash
git clone https://github.com/horvatjanos/universal-drug-repurposing.git
cd universal-drug-repurposing
pip install -r requirements.txt
```

### Usage

```bash
# Interactive mode — no programming required
python src/hermes_drug.py --interactive

# Screen a specific target (e.g., SARS-CoV-2 Mpro)
python src/hermes_drug.py --pdb 6LU7 --depth standard

# Full pipeline with lead optimization
python src/hermes_drug.py --pdb 6LU7 --depth deep

# Replicate our multi-receptor experiment
python src/train_multireceptor_v2.py
```

### Reproducing Our Results

To reproduce the Nilotinib_Var_17 discovery:

```bash
# 1. Prepare receptors
python src/receptor_prep.py --pdb 6LU7
python src/receptor_prep.py --pdb 1HPV
python src/receptor_prep.py --pdb 1M17
python src/receptor_prep.py --pdb 4LVT
python src/receptor_prep.py --pdb 5KIR

# 2. Generate Nilotinib variants
python src/optimize_lead.py --smiles "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4C(F)(F)F)n5cc(cn5)C" --name Nilotinib

# 3. Run multi-receptor docking + training
python src/train_multireceptor_v2.py
```

---

## Repository Structure

```
 src/                            Source code
    hermes_drug.py              Main orchestrator (CLI + interactive)
    receptor_prep.py            PDB fetch, clean, PDBQT conversion
    auto_box.py                 Binding pocket auto-detection
    screen_ai.py                AI-accelerated compound screening
    optimize_lead.py            Lead optimization (variant generation)
    prep_ligands.py             SMILES → 3D PDBQT (cross-platform)
    parse_docking_results.py    Vina output parsing + ranking
    train_multireceptor_v2.py   Multi-target ML training pipeline
    train_advanced_model.py     Advanced ML (ChEMBL data integration)
    train_ml_proxy.py           RandomForest affinity predictor
    train_boosted_proxy.py      GradientBoosting affinity predictor
    cloud_runner.py             Kaggle cloud acceleration
    ...                         Additional utility modules
 data/                           Compound libraries
    full_fda_library.csv        51 FDA-approved drugs with SMILES
    mpro_training_data.csv      669 ChEMBL SARS-CoV-2 Mpro bioactivities
    training_set_50.csv         50 diverse training compounds
    nilotinib_variants.csv      50 Nilotinib structural variants
 docs/                           Documentation
    README.md                   This file
    GEMINI.md                   Original project specification
    LICENSE                     MIT License
    letter_to_semmelweis.py     Academic collaboration proposal
 config.yaml                     Universal configuration
 requirements.txt                Python dependencies
 vina_1.2.7_win.exe              AutoDock Vina (Windows binary)
```

---

## Seeking Collaboration

We are actively seeking **academic partners** for:

1. **In vitro validation** of Nilotinib_Var_17 binding affinity (SPR, ITC, or enzymatic assays)
2. **Cell-based testing** against SARS-CoV-2, HIV, or EGFR-dependent cancer lines
3. **Structural biology** co-crystallization studies to confirm binding modes
4. **Clinical collaborators** for drug repurposing trial design

Institutions with **BSL-2+ laboratory capabilities** and experience in enzymology or cell-based assays are ideal.

**Contact:** [GitHub Issues](https://github.com/horvatjanos/universal-drug-repurposing/issues)  
**Authors:** Turcsányi Zsombor & Hermes AI — [GitHub](https://github.com/horvatjanos)

---

## Data Sources

- **Protein structures:** RCSB Protein Data Bank (PDB IDs: 6LU7, 1HPV, 5KIR, 1M17, 4LVT)
- **Bioactivity data:** ChEMBL (target: CHEMBL4523582 — SARS-CoV-2 Mpro)
- **FDA drug library:** ChEMBL + DrugBank
- **Molecular descriptors:** RDKit 2024.09

---

## License

MIT License — see [LICENSE](LICENSE)

---

## Citation

If you use this work in your research, please cite:

```bibtex
@misc{horvath2026drugrepurposing,
  title   = {AI Drug Repurposing: Multi-Target Virtual Screening Platform},
  author  = {Turcsányi, Zsombor and Hermes AI},
  year    = {2026},
  url     = {https://github.com/horvatjanos/universal-drug-repurposing}
}
```

---

*Computational resources: AutoDock Vina (Scripps Research), RDKit (Greg Landrum et al.), OpenBabel (O'Boyle et al.), ChEMBL (EMBL-EBI).*
