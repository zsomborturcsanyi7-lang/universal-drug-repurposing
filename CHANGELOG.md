# Changelog

All notable changes to this project will be documented in this file.

## [2.0.0] - 2026-06-01
### ⚠️ SCIENTIFIC OVERHAUL & CRITICAL BUG FIXES

This version represents a major restructuring of the project to address methodological flaws in previous releases.

### Added
- **ADMET Prediction Module:** New `src/admet_prediction.py` to estimate drug safety and pharmacokinetics (GI absorption, BBB permeability, Solubility).
- **Pose Stability Analysis (MD Proxy):** New `src/pose_stability.py` that calculates ligand strain energy to ensure realistic binding modes.
- **Scientific Lead Filtering:** All generated variants now pass through Lipinski's Rule of Five, QED, and PAINS/Brenk toxicophore filters.
- **Parallel Docking Engine:** New `src/run_parallel_docking.py` for efficient, multi-core screening of variants.
- **Enhanced Result Parsing:** `final_3d_docking_insights.csv` now includes molecular descriptors (MW, LogP, QED) alongside binding energy.
- **Auto-Correction for Kinase Inhibitors:** Lenient filtering mode for large drugs like Nilotinib (allowing up to 2 Lipinski violations).

### Changed
- **Increased Docking Precision:** Default `exhaustiveness` for Vina increased from 8 to **32** in deep mode, drastically reducing false positive affinity scores.
- **Improved PDBQT Generation:** RDKit connectivity records (CONECT) are now preserved to ensure proper torsion tree generation by OpenBabel.
- **Unified Pipeline Paths:** All internal script calls now correctly reference the `src/` directory.

### Fixed
- **Refuted Invalid Results:** High-precision docking proved that previous "winner" candidates (e.g., Nilotinib_Var_17) were artifacts of low-precision simulation.
- **Unicode Errors:** Fixed crashes when reading `config.yaml` on Windows systems.
- **PDBQT Parsing Errors:** Resolved "Unknown tag" errors in Vina 1.2.x by fixing ligand preparation headers.
- **Import Errors:** Fixed `UnboundLocalError` in `auto_box.py`.

## [1.0.0] - 2026-05-24
- Initial release (Legacy / Flawed version)
