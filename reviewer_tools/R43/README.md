# Corrected numerical release: R42 baseline / R43 engineering interpretation

This is the numerical companion to manuscript revision R43. The underlying corrected-current results are locked to R42; the only added numerical illustration is the worked laboratory decision. The manuscript, draft PDFs, private correspondence and full third-party raw ZIP archives are **not** deposited here.

## Start here

Use Python 3.11+ and the dependencies in `requirements.txt`. Set `PYTHONDONTWRITEBYTECODE=1` (or `python -B`) to keep the release immutable.

```bash
python -B code/verify_public.py
python -B code/worked_engineering_example.py
```

These read-only commands check the complete public manifest, arithmetic and source-window consistency, stored proof records and exact-zero sparse-support identities. They do **not** claim a fresh MPFR run.

## Fresh continuous-domain verification

Install a C++17 compiler and MPFR/GMP development libraries. The unchanged verifier source is compiled locally; no precompiled executable is distributed. Output directories must be new and outside this release.

```bash
python -B code/replay_saved.py --family primary --workers 4 --output ../fresh_primary
python -B code/replay_saved.py --family per_case --workers 4 --output ../fresh_sparse_case
python -B code/replay_saved.py --family soc_union --workers 4 --output ../fresh_sparse_soc
python -B code/replay_saved.py --family truncated --workers 4 --output ../fresh_truncated
python -B code/replay_external_fresh.py --workers 4 --output ../fresh_external
```

Expected positive-witness runs are 168 full-band, 168 per-certificate sparse, 168 SOC-common sparse, 8 truncated-band and 48 external. The external runner additionally checks six feasible models (852 EIS coordinates), 24 arithmetic witnesses and 12 invalid controls. Counts describe computations on operating conditions, not independent cell experiments. `replay_saved.py --case SOC30_discharge_0.7C` provides a clearly delimited test of the engineering example's witnesses. No cache reuse is allowed by that runner.

## Raw-to-instance regeneration

Work in a scratch copy because generation writes result files. Obtain the source archive from the original RWTH record, extract `MatFiles.zip`, and check its exact SHA-256:

`875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f`

```bash
python -B code/extract_rwth_r42.py --raw /path/to/MatFiles.zip
python -B code/recompute_temperature.py
python -B code/independent_raw_audit.py --raw /path/to/MatFiles.zip
```

The last command can also audit the supplied selected windows without `--raw`, but that narrower mode does not reopen the original archive. Numerical solver and library versions may affect last-place reproduction; no fresh numerical run is claimed from a manifest-only check. The original R42 reference environment is in `requirements-lock.txt`.

## Results and boundaries

The primary 42 conditions belong to one RWTH cell. All 84 full-band saved witnesses remain positive at the two stated EIS tolerance factors. Endpoint-window local 3-sigma expansion leaves minima 1.2359252691862501 and 0.009321772570726054 mOhm. It is recorded scatter, not a 99.7% confidence interval or calibrated instrument uncertainty. The temperature slice is an exploratory grid calculation. Sparse supports are retrospective and use the original full-band tolerance, not a prospective frequency-acquisition policy. The external result remains conditional on source-reconciled pairing and state persistence.

`evidence/R43_Worked_Engineering_Example.json` brackets the weakest doubled-tolerance local-scatter margin with **hypothetical** additional pulse-only bounds. Certificate loss is not proof of compatibility or cell failure.

## Source attribution and provenance

RWTH DOI: https://doi.org/10.18154/RWTH-2025-04811
Official record: https://publications.rwth-aachen.de/record/1012023

Bayreuth version DOI: https://doi.org/10.5281/zenodo.15755725
Official record: https://zenodo.org/records/15755725

See `SOURCE_ATTRIBUTION.md` for provider attribution. Source-derived windows retain the original data attribution. This publication does not assign a new software license or override any provider's rights.

`SOURCE_FILE_HASHES.json` maps unchanged scientific files to the R42 source package. `baselines/` preserves numerical provenance; old development scripts are not the canonical entry points. Use only the commands above for this corrected release. The complete external archive remains byte-identical, including its original internal manifest. Legacy R38 root files elsewhere in the repository must not be mixed with this corrected numerical release.
