# Universal Drug Repurposing — AI gyógyszer-újrahasznosítás: AutoDock Vina screening + ML

**Status:** ⚠️ Prototype — teljes pipeline validálva, Kaggle GPU kell a teljes DrugBank futtatáshoz

Teljes, automatizált gyógyszer-újrahasznosító platform: AutoDock Vina docking + ML proxy gyors előszűrés. Receptor előkészítés → docking → ML kiértékelés → riport generálás, egyetlen paranccsal.

**Megjegyzés:** Ez a repo a fő/átfogó verzió. A `drug-repurposing-pipeline` egy kisebb al-repo ugyanezzel a koncepcióval.

## ⚠️ THIS PROJECT IS UNFINISHED — FEEL FREE TO CONTINUE IT ⚠️

**Ez a projekt NINCS KÉSZEN. Bárki folytathatja, aki akarja!**
Ezt a projektet Zsombi & Hermes Agent (Nous Research) közösen fejlesztette, de egyik projekt sincs 100%-osan befejezve. Ha tetszik az ötlet és tovább fejlesztenéd, nyugodtan fork-old, folytasd, és csinálj belőle valami nagyszerűt!

---

## Pipeline lépések

1. **Receptor Preparation** — PDB struktúrák letöltése és tisztítása
2. **Auto Box** — Automatikus kötőhely detektálás
3. **AI Screening** — ML proxy alapú gyors előszűrés
4. **Ligand Preparation** — SMILES → 3D konformáció → PDBQT
5. **Parallel Docking** — Párhuzamosított AutoDock Vina futtatás
6. **Result Evaluation** — Rangsorolás és riport generálás

## Tech stack
- AutoDock Vina 1.2.7
- ML proxy (tanított modell)
- Python pipeline
- YAML konfiguráció

## Fájlok
| Fájl | Leírás |
|------|--------|
| `src/` | Forráskód |
| `data/` | Adatok |
| `docs/` | Dokumentáció |
| `config.yaml` | Konfiguráció |

## Kapcsolódó repók
- `drug-repurposing-pipeline` — Kisebb, átfedő verzió

## Fejlesztő
Zsombi & Hermes Agent (Nous Research)
