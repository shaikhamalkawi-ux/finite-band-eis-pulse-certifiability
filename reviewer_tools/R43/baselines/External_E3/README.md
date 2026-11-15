# Project 52: Bayreuth external saved-certificate verification (E3)

This archive verifies 12 paired pulse cases on 3 source cells, at 2 EIS tolerances.
It contains processed computational instances, all 24 signed dual witnesses,
6 explicit EIS-feasible atomic models, complete interval partitions, candidate
history, verifier code, invalid controls, and results. It is not 24 independent
battery experiments, a calibration certificate, or a measured state-drift bound.

## Saved-input verification
Install Python 3.11+, a C++17 compiler, MPFR/GMP development libraries, and the
Python dependencies in `requirements.txt`. On Debian/Ubuntu the system packages
are `g++ libmpfr-dev libgmp-dev`.

Run from any directory:

```bash
python code/replay_external.py --workers 4
```

The runner checks every SHA-256 entry, compiles the verifier, checks all 48
full-domain enclosures, independently checks objectives/finite coordinates/tails,
checks 852 constraints on six nonempty models, and rejects 12 invalid witnesses.
It writes to a new temporary output directory and never reoptimizes a witness.
Do not use Python -O or PYTHONOPTIMIZE. Expected outcome: PASS.

## Scope and source preparation
The producer clarified the nominal 25 C condition and use of outer archive cell
labels. The later `F0002/f0002_bol_2.txt` is the selected BOL EIS file. Nominal
50% SOC and state continuity during OCV rest remain conditional; the longest
EIS-pulse delay is approximately 149 hours. The two F0001 EOL cases are also
reported separately through the 10-case timing sensitivity. All outcomes are retained.

All stored binary64 input numbers define the computational instance exactly.
The observed recording envelope is not an absolute calibration or a confidence
interval. It does not propagate every possible waveform/time/state uncertainty.
No optimizer success flag enters acceptance of the saved continuum certificates.

The code supports signed initial current samples without clipping. The slow-tail
bound uses the absolute first time moment; six prior positive-waveform witnesses
and two signed analytic toy checks were separately regressed. Their records are
included; rerunning those six controls requires the original R36 archive.

## Original data
Source dataset: DOI 10.5281/zenodo.15755725, record version 15755725.
Paper: Ramasubramanian et al., Journal of Energy Storage 140 (2025) 119113,
DOI 10.1016/j.est.2025.119113. Original dataset license: CC BY 4.0.
The original large cell ZIPs and private correspondence are not redistributed
here. To recreate raw parsing, obtain F0001.zip, F0002.zip and F0003.zip, place
in `raw/`, and consult `code/prepare_external_e3.py`. Its historical-delta ledger uses the twelve preserved E1 pulse records included
under history/. Raw preparation should be run in a separate scratch copy, since
it regenerates pre-tolerance instances. The saved-input verification above is
the frozen, solver-free acceptance path. New code has no new
public license assigned by this package. No repository push or DOI is claimed.
