# universal-drug-repurposing

Automated molecular docking and machine learning screening pipeline for drug repurposing (AutoDock Vina).

## 📌 Overview & Purpose
Performs automated screening of existing drug databases against target receptor proteins to identify potential novel therapeutic applications.

## ⚙️ Tech Stack & Architecture
- **Language**: Python
- **Bioinformatics**: AutoDock Vina, RDKit, Biopython
- **Machine Learning**: Scikit-Learn proxy model for rapid affinity estimation

## 🚀 Installation & Quickstart

### Prerequisites
- Python 3.9+
- AutoDock Vina binary installed

### Steps
```bash
git clone https://github.com/zsomborturcsanyi7-lang/universal-drug-repurposing.git
cd universal-drug-repurposing

pip install -r requirements.txt
python run_pipeline.py --receptor protein.pdbqt --ligands dataset.sdf
```

## 📊 Project Status
⚠️ **In Silico Experimental Pipeline**.
