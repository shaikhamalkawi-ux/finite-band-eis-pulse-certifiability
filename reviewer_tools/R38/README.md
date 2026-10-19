# R38 — Verified computational supplement and external-cell check

This directory supports the R38 EP1 submission edition of
**Certifiability Before Prediction: Finite-Band EIS Bounds for Lithium-Ion Battery Pulse Diagnostics**,
prepared for IEEE Open Journal of Vehicular Technology (OJVT). No claim is made that a journal submission has occurred.

## What has been tested

- Primary RWTH evidence: 42 paired operating conditions measured on **one** Samsung cell, with 84 saved full-band case–tolerance witnesses.
- Additional Bayreuth evidence: 12 paired pulse conditions from **three** Varta CP1654 cells, with 24 saved signed dual witnesses and 48 saved full-domain verification checks at two tolerance factors. These are 12 conditions, *not* 12 independent cells.
- A conservative 10-case timing sensitivity omits the F0001 end-of-life pair (two pulse cases), which had a long EIS–pulse interval.
- Producer correspondence clarified the nominal 25 °C experiment and outer-archive cell identifiers and identified the later F0002 beginning-of-life EIS as the intended file. Continuous state matching during long rest periods and absolute instrument calibration are **not established**.
- The certificates establish conditional paired-reference **model-set exclusion** under the stated positive-RC/Stieltjes class and computational tolerances, not an EIS-only early-stopping rule, exact model identification, causal mechanism, calibrated error rate, or cell-population result.

## Download/reconstruct the exact reviewer supplement

GitHub transfer limits required splitting the **exact 3,747,969-byte final computational ZIP** into four files in this directory, without altering its bytes.

```bash
git clone https://github.com/shaikhamalkawi-ux/finite-band-eis-pulse-certifiability.git
cd finite-band-eis-pulse-certifiability/reviewer_tools/R38
python assemble_archive.py
```

The script verifies SHA-256 for every part and for the assembled ZIP, then tests ZIP CRC. Required combined archive SHA-256:

`e8b47aad8f1698a8af614a93403c522a6c735639638cde9b0786e8091a3da347`

The combined ZIP contains:
- `01_RWTH_R37_Frozen_Replay.zip` (retained original R37/RWTH 42-case witness archive)
- `02_External_E3_Replay.zip` (the R38 new external-cell witness archive)
- `README.txt`

The external ZIP alone is also mirrored as `E3_External_Replay.zip` here (SHA-256 `5966ddd0b9232edd6c4266d6f5876ea57bf380c39471b92c0e951006a0af1826`).

To replay the external saved-input evidence, unpack `E3_External_Replay.zip` and run from its `Project52_External_Verification` directory:

```bash
python code/replay_external.py --workers 4
```

Dependencies: Python 3.11+, C++17 compiler, MPFR and GMP development headers, and Python dependencies listed in the archive's `requirements.txt`. On Ubuntu/Debian, `g++ libmpfr-dev libgmp-dev` are needed. Do not use `python -O` or `PYTHONOPTIMIZE`. Raw-source rebuilding is a separate workflow, not part of the saved-certificate replay.

## Evidence and source constraints

`results/` holds concise 12-case and 24-certificate tables. `code/` holds the exact E3 computational code; `qa/` contains selected saved audits. The ZIP contains the full signed witnesses, computational instances, interval partitions and test records.

Original third-party RWTH and Bayreuth large raw archives and private email correspondence are **not redistributed**. See RWTH source DOI https://doi.org/10.18154/RWTH-2025-04811 and Bayreuth data DOI https://doi.org/10.5281/zenodo.15755725, associated with Ramasubramanian et al. (2025), https://doi.org/10.1016/j.est.2025.119113 . The producer's source dataset is CC BY 4.0; this repository does not assign an additional software license.

These are archived, AI-assisted computations and replay tests, not an independent human mathematical or instrument-calibration audit.