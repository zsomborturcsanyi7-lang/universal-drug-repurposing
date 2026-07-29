# universal-drug-repurposing

**Status:** ⚠️ Prototype — full pipeline validated, Kaggle GPU needed for full DrugBank run

Complete automated drug repurposing platform: AutoDock Vina docking + ML proxy pre-screening. Receptor prep → docking → ML evaluation → report generation, all in one command.

**Note:** This is the main/comprehensive version. `drug-repurposing-pipeline` is a smaller overlapping repo with the same concept.

## ⚠️ THIS PROJECT IS UNFINISHED — FEEL FREE TO CONTINUE IT ⚠️

This project was developed by Zsombi & Hermes Agent (Nous Research).

---

## Pipeline
1. **Receptor Preparation** — Download and clean PDB structures
2. **Auto Box** — Automatic binding pocket detection
3. **AI Screening** — ML proxy-based fast pre-screening
4. **Ligand Preparation** — SMILES → 3D conformation → PDBQT
5. **Parallel Docking** — Parallelized AutoDock Vina execution
6. **Result Evaluation** — Ranking and report generation

## Tech stack
- AutoDock Vina 1.2.7
- ML proxy (trained model)
- Python pipeline, YAML config

## Related
- `drug-repurposing-pipeline` — Smaller overlapping version

## Developer
Zsombi & Hermes Agent (Nous Research)
