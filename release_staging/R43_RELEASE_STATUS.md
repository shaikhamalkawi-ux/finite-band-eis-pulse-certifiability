# R43 corrected-current reproducibility synchronization — STAGING ONLY

Date: 2026-11-16.

**Not a completed numerical release. Do not cite this branch as the source of the R43 calculations.**

The manuscript is **Certifiability Before Prediction: Testing EIS–Pulse Consistency in Lithium-Ion Batteries** (R43 author-review). The corrected numerical baseline is R42, derived from measured-current normalization and fresh raw-trace precision checks. Existing public `main` contains historical R38 materials and must not be presented as the corrected release.

The complete locally prepared R43 public overlay has **2,456 files**, ~18.84 MB uncompressed, and is held outside GitHub in `Project52_Batteries_R43_GitHub_Sync_READY_NOT_PUSHED.zip`. It comprises reproducibility code, processed source windows, fixed computational instances, signed dual certificates, MPFR verifier source, stored verification records, checksums, and a README. The overlay intentionally excludes the complete third-party raw archives, draft manuscript PDF, and private correspondence.

**Publication acceptance gate:** verify the complete payload's SHA-256 manifest, upload all 2,456 source-approved files without modifications, run `python -B code/verify_public.py` inside `reviewer_tools/R43`, run the worked example and fresh verifier as documented, verify content readback and the final immutable commit. Only then promote to `main` and update manuscript Data and Code Availability and Supplement S9 with that exact commit URL.

Current public source DOIs: RWTH 10.18154/RWTH-2025-04811; Bayreuth 10.5281/zenodo.15755725. No claim of instrument calibration, independent 42-cell testing, onboard BMS deployment, or measured acquisition time saving.

The staging branch is a transparent audit waypoint, not the final journal-facing reproducibility mirror.
