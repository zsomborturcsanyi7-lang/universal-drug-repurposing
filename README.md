# universal-drug-repurposing

Automatizált molekuláris dokkolási és gépi tanulási szűrési pipeline gél- és gyógyszer-újrahasznosításhoz (AutoDock Vina).

## 📌 A projekt célja
A meglévő hatóanyag-adatbázisok automatizált szűrése célfehérjék ellen a potenciális új terápiás alkalmazások gyors azonosítására.

## ⚙️ Technológiai stakk & Működés
- **Nyelv**: Python
- **Bioinformatika**: AutoDock Vina, RDKit, Biopython
- **Gépi tanulás**: Scikit-Learn proxy modell szelektivitáshoz

## 🚀 Telepítés és Használat

### Előfeltételek
- Python 3.9+
- AutoDock Vina telepítve

### Lépések
```bash
git clone https://github.com/zsomborturcsanyi7-lang/universal-drug-repurposing.git
cd universal-drug-repurposing

pip install -r requirements.txt
python run_pipeline.py --receptor protein.pdbqt --ligands dataset.sdf
```

## 📊 Status
⚠️ **In silico kísérleti pipeline**.
