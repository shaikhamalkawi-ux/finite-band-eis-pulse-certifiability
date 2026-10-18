# Finite-Band EIS Pulse Certifiability

Reviewer-facing reproducibility repository for the manuscript:

**Certifiability Before Prediction: Finite-Band EIS Bounds for Lithium-Ion Battery Pulse Diagnostics**

## Current reviewer archive (R37)

The all-case directed-rounding verifier and systematic-correction sensitivity are supplied in the accompanying **`Batteries_R37_Reviewer_Replay.zip`**, SHA-256 `277c73ff96923904d9eff9980dece66b52d4484add76b1913b9944962a06e1a1`. See [the R37 reviewer guide and exact manuscript/archive pairing](reviewer_tools/R37/README.md) and [the archive checker](reviewer_tools/R37/verify_review_archive.py).

**This repository does not yet mirror the complete R36/R37 archive.** The files below are the earlier public analysis layer; in particular, the old representative-case verifier is not the new all-42 verifier. For full review, use the computational archive accompanying the manuscript, not the manuscript PDFs alone. The new result is paired-reference model-set exclusion, not an EIS-only stopping policy, calibrated instrument validation or a multiple-cell population claim.

This repository contains the analysis code and derived tables used for the finite-band EIS-to-pulse certifiability study. It intentionally does **not** redistribute the original third-party RWTH Aachen raw dataset.

## Source data

The source dataset is publicly available from RWTH Aachen University:

- DOI: https://doi.org/10.18154/RWTH-2025-04811

For scripts that replay raw-trace or EIS calculations, download the official source archive and place the original `MatFiles.zip` in the repository root.

Expected `MatFiles.zip` SHA-256:

`875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f`

Run:

```bash
python verify_raw_data.py
```

before any raw-data-dependent replay.

## Contents

- **12 Python scripts** in the repository root for analysis, validation, and source-integrity checks.
- **28 derived CSV tables** in `derived_data/`, with the principal script inputs also mirrored at the repository root for direct replay.
- The repository contains derived/reproducibility materials only; it does not contain the original third-party RWTH raw archive.
- `RWTH_Data_Integrity_Manifest.txt` - source archive integrity record.
- `R10_DATA_DERIVED_PRECISION_METHOD.txt` - precision-audit method note.
- `R19_SELECTION_AND_EQ9_README.txt` - selection/Eq. (9) technical note.

## Python dependencies

Install the required packages with:

```bash
python -m pip install -r requirements.txt
```

The scripts use Python with NumPy, pandas, SciPy, matplotlib, and mpmath. Reviewer-facing files are directly browsable because the repository is public.

## Reproduction notes

Run scripts from the repository root. Some scripts use only the included derived tables; others require the official `MatFiles.zip` source archive.

Examples that can be run directly from the included derived tables:

```bash
python R19_build_historical_selection_matrix.py
python R23_compute_vehicle_design_table.py
python R25_build_adaptive_fast_eis_map.py
```

Raw-data-dependent scripts will stop if `MatFiles.zip` is absent. The original raw data are not committed to this repository.

## Scientific scope

The repository reproduces the finite-band information-radius calculations, tolerance-sensitivity analyses, representative continuum-validation checks, and retrospective/adaptive fast-EIS screening analyses reported in the manuscript. The code and tables are provided for verification of the reported computational results; they do not change the scientific baseline of the manuscript.
