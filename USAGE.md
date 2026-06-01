# 📖 Használati Útmutató / User Guide (v3.0)

Ez a dokumentum részletes útmutatást nyújt az **AI Drug Repurposing Platform** telepítéséhez és használatához.
This document provides detailed instructions for installing and using the **AI Drug Repurposing Platform**.

---

## 🇭🇺 MAGYAR NYELVŰ ÚTMUTATÓ

### 1. Rendszerkövetelmények
A futtatáshoz Python 3.10 vagy újabb verzió szükséges.
*   **Operációs rendszer:** Windows 10/11, Linux (Ubuntu ajánlott) vagy macOS.
*   **Hardver:** Minimum 8GB RAM, többmagos CPU (a párhuzamos dokkoláshoz).

### 2. Telepítés lépésről lépésre

1.  **Klónozd a tárolót:**
    ```bash
    git clone https://github.com/zsomborturcsanyi7-lang/universal-drug-repurposing.git
    cd universal-drug-repurposing
    ```

2.  **Telepítsd a Python függőségeket:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Külső tudományos szoftverek telepítése (KÖTELEZŐ):**
    *   **AutoDock Vina:** Töltsd le a [hivatalos oldalról](https://github.com/ccsb-scripps/AutoDock-Vina/releases). Windows-on nevezd át a fájlt `vina.exe`-re és add hozzá a rendszer PATH-hoz, vagy másold a `src/` mappába.
    *   **OpenBabel:** Töltsd le az [OpenBabel-t](http://openbabel.org/). Szükséges a molekulák 3D formátumra alakításához.

### 3. Használat

#### A. Egyszerű szűrés (Alap üzemmód)
Ha egy adott fehérje (PDB ID) ellen szeretnél gyorsan szűrni:
```bash
python src/hermes_drug.py --pdb 6LU7
```

#### B. Mély Tudományos Validáció (Ajánlott)
Ez az üzemmód lefuttatja a v3.0 összes újdonságát: megnövelt precizitású dokkolás (Exhaustiveness: 32), ADMET jóslás és stabilitási vizsgálat.
```bash
python src/hermes_drug.py --pdb 6LU7 --depth deep
```

#### C. Interaktív mód (Kezdőknek)
Vezetett folyamat, ahol a program kérdéseket tesz fel:
```bash
python src/hermes_drug.py --interactive
```

#### D. Felhő alapú futtatás (Kaggle)
Nagyobb szűrésekhez használd a Kaggle ingyenes erőforrásait:
```bash
python src/cloud_runner.py --pdb 6LU7 --depth deep
```

### 4. Eredmények értelmezése
A futtatás végén a rendszer egy `final_validated_results.csv` fájlt generál:
*   **Binding_Affinity:** Minél negatívabb (pl. -9.0), annál erősebb a kötődés.
*   **Drug_Score:** 0 és 1 közötti érték. 0.5 felett már jó gyógyszerjelölt.
*   **Toxicity_Alerts:** Ha "True", a molekula potenciálisan mérgező.
*   **Pose_Strain:** Ha alacsony (< 3.0), a kötődési mód fizikailag stabil.

---

## 🇬🇧 ENGLISH USER GUIDE

### 1. System Requirements
Python 3.10 or higher is required.
*   **OS:** Windows 10/11, Linux (Ubuntu recommended), or macOS.
*   **Hardware:** Min. 8GB RAM, Multi-core CPU (for parallel docking).

### 2. Installation Step-by-Step

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/zsomborturcsanyi7-lang/universal-drug-repurposing.git
    cd universal-drug-repurposing
    ```

2.  **Install Python dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **External Scientific Software (MANDATORY):**
    *   **AutoDock Vina:** Download from the [official releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases). On Windows, rename to `vina.exe` and add to PATH or copy to `src/`.
    *   **OpenBabel:** Download from [OpenBabel](http://openbabel.org/). Required for 3D ligand preparation.

### 3. Usage

#### A. Basic Screening
Quickly screen a protein target using default settings:
```bash
python src/hermes_drug.py --pdb 6LU7
```

#### B. Deep Scientific Validation (Recommended)
This mode executes the full v3.0 pipeline: High-precision docking (Exhaustiveness: 32), ADMET profiling, and Pose Stability analysis.
```bash
python src/hermes_drug.py --pdb 6LU7 --depth deep
```

#### C. Interactive Mode
Guided wizard for non-programmers:
```bash
python src/hermes_drug.py --interactive
```

#### D. Cloud Execution (Kaggle)
Push heavy workloads to Kaggle's infrastructure:
```bash
python src/cloud_runner.py --pdb 6LU7 --depth deep
```

### 4. Interpreting Results
The system generates a `final_validated_results.csv` file:
*   **Binding_Affinity:** More negative values (e.g., -9.0) indicate stronger binding.
*   **Drug_Score:** Value between 0 and 1. Scores above 0.5 are promising leads.
*   **Toxicity_Alerts:** If "True", the molecule contains potentially toxic fragments.
*   **Pose_Strain:** If low (< 3.0), the binding mode is physically stable (MD Proxy).

---

## 📜 Licenc / License
Distributed under the MIT License. See `LICENSE` for more information.
