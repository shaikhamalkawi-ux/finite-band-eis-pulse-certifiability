# R37 reviewer replay: exact artifact and scope

This directory publishes a small, SHA-pinned checker for the complete reviewer archive supplied with the manuscript. **The complete R36/R37 witness archive is not yet mirrored in this GitHub repository.** The older scripts at the repository root are not a substitute for the all-case verifier.

## Required accompanying files

Request or use the single `Batteries_R37_Reviewer_Kit.zip` supplied for review. It contains the main manuscript, Supporting Information, the complete computational archive, and audit instructions. The following pairing identifies the R37 edition reviewed in this conversation; do not match files by title alone.

| File | SHA-256 |
|---|---|
| `Batteries_R37_Main.pdf` | `17e569a5443e8d93f8c46a188c66cce85b277d60364a49e53aaf3422ac3b75d8` |
| `Batteries_R37_Supplement.pdf` | `eb76e41947aa353552b816e1feaaf278290407acfcfe87d024015e1f0fe4e9b2` |
| `Batteries_R37_Reviewer_Replay.zip` | `277c73ff96923904d9eff9980dece66b52d4484add76b1913b9944962a06e1a1` |

The paper has 11 pages, the supplement 18 pages and the bibliography 37 references. The frozen replay archive contains 746 entries, including 745 manifest-checked files.

## Run

On Linux, install Python 3.11+, a C++17 compiler, GMP and MPFR development packages, and `mpmath==1.3.0`. For Debian/Ubuntu the system packages are `g++ libmpfr-dev libgmp-dev`.

From this directory, after obtaining the accompanying archive:

```bash
python verify_review_archive.py /path/to/Batteries_R37_Reviewer_Replay.zip
python verify_review_archive.py /path/to/Batteries_R37_Reviewer_Replay.zip --run --workers 4
```

The first command checks the exact archive and every manifest entry. The second extracts into a new directory and runs the unchanged frozen verifier and sensitivity scripts. It rejects Python optimization mode because inherited scientific checks use assertions.

Expected core outcomes: 84 saved case-tolerance witnesses; 168 full-domain checks; 42 invalid controls rejected; seven nonempty EIS models checked; seven sensitivity outputs reproduced byte-for-byte. The archive includes input instances, witnesses, continuum partitions, development history, programs and outputs. No optimizer or battery-data download is needed for this saved-input replay. Raw-to-instance preparation is a separate documented procedure requiring the original RWTH archive.

## What a PASS means

A fresh-copy replay by the same AI-assisted preparation workflow reproduced those outcomes on 19 October 2026. Four deliberate file corruptions and one missing-file test were rejected before calculations. All frozen files were restored to their original hashes afterwards.

These checks are **not an independent human review, new battery experiments, instrument calibration, or external-cell validation**. There are 42 within-cell operating conditions on one physical cell, not 84 independent experiments. A lost fixed-witness certificate under hypothetical corrections is inconclusive; it does not construct a compatible model. Paired-reference exclusion needs measured pulse data and is not EIS-only prediction accuracy or an acquisition-stop policy.

No scientific data, equation, witness, threshold or result in the frozen archive has been modified by this checker. No new license is assigned. No raw third-party archive or manuscript PDF is committed here.
