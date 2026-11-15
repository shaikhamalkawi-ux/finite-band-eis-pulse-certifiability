# Finite-Band EIS Pulse Certifiability

## Corrected numerical release: R42 baseline / R43 manuscript interpretation

The corrected-current reproducibility materials are in **[reviewer_tools/R43](reviewer_tools/R43/README.md)**. Start there rather than with the legacy root-level tables. This release preserves the R42 calculation lock and supplies an executable worked engineering example used in manuscript revision R43.

The study tests **EIS--pulse model-set consistency**, not battery failure, causal mechanism, a calibrated confidence-level rejection, or an onboard BMS policy. Its primary 42 observations are operating conditions on one cell. The sparse-frequency result is retrospective certificate compression at a fixed tolerance, not measured acquisition-time savings.

### Reproduce the corrected results

```bash
cd reviewer_tools/R43
python -B -m pip install -r requirements.txt
python -B code/verify_public.py
python -B code/worked_engineering_example.py
```

The release guide distinguishes read-only evidence checking, numerical regeneration from the original source archive, and fresh MPFR certificate replay. The public mirror includes processed cases, source-window extractions, certificates, code, complete saved verification records and hashes. No draft manuscript, full third-party raw archive, private correspondence or new public software-license declaration is included.

### Source data

- RWTH: https://doi.org/10.18154/RWTH-2025-04811
- Bayreuth version 1: https://doi.org/10.5281/zenodo.15755725

The required RWTH `MatFiles.zip` SHA-256 is `875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f`. Download original raw data from their providers and preserve their attribution.

### Historical material

[reviewer_tools/R38](reviewer_tools/R38/README.md), [reviewer_tools/R37](reviewer_tools/R37/README.md), and the older root scripts/derived tables remain for provenance. They are not the current corrected-current results and must not be mixed with the R42/R43 companion release. The previous published baseline was commit `e27d428bc0ada47b5bfb23c4b5bd749b5cac39bf`.
