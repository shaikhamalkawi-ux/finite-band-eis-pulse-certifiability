# Project 52 — E3 external verification protocol, 19 October 2026

R37 EP1-T1 is preserved, not replaced. E1 grid-only results have already been seen; this protocol is frozen before new E3 optimizations, not externally preregistered. No outcome-based case selection.

## Panel and source rule
F0001, F0002, F0003; BOL and EOL; charge and discharge. All 12 originally selected candidate pairs must be reported. Producer clarification in user-supplied correspondence supports archive identities and nominal 25 C; the earlier F0002 BOL file is superseded by f0002_bol_2.txt. f0056_eol inside F0001 remains excluded. Raw headers are not rewritten. Retain the full 12 and the predesignated 10-case timing sensitivity omitting F0001 EOL. The latter is post-intake, not an admissibility time threshold. Absolute calibration and state persistence remain unverified.

## Raw processing
Parse all rows without skipping malformed records. EIS: final descending branch from the maximum frequency (71 points), ascending sort only for computation; retain both signed components in mOhm. Pulse command and approximately 10-second duration identify the unique charge/discharge pulse per checkup. Normalize measured current by magnitude of its median in [2,8] seconds; retain small signed initial samples. Baseline is the median of the last 10 near-zero-current records before the pulse (absolute current <0.005 A). For the primary certified instance, T is the exact binary64 elapsed time of the final active pulse record; u is left ZOH through that T; iT is its left value. Compare to the recorded voltage at that same T. This corrects the E1 pilot's use of a nominal 10-second kernel with a last-at-or-before-10s endpoint. It needs no extrapolation and no outcome choice. The nominal duration remains 10 seconds, but actual horizons are reported. Preserve an E1-to-E3 preprocessing ledger. Verify the right-endpoint coefficient as a separate robustness check.

Local recorded envelope only: baseline last10 minima/maxima, endpoint last10 active samples minima/maxima, plateau [2,8] current magnitude minima/maxima. Outward round the resulting response interval. This is not a confidence interval, calibration interval, or uncertainty on the entire waveform; the waveform/time computational input remains fixed.

## EIS-only tolerance and non-vacuity
For each of six spectra use the original E1 1201-atom logspace [-10,10] finite-grid EIS-only minimum common componentwise tolerance plus the pre-existing 0.001 mOhm numerical buffer. Retain explicit nonnegative primal atoms and interval-check all 142 constraints. The fixed literal tolerance is a verified feasible grid-derived level, not an exact continuum admission minimum or physical calibration. Report factors 1 and 2 separately. Never increase physical tolerances to get a favorable outcome.

## Certificates
Reuse R36 candidate LP with initial 1001 time constants, 0.005 common numerical surplus and finite-coordinate surplus 1e-7. Two full-domain MPFR enclosures: 160-bit Taylor and 256-bit mean-value, with analytic fast/slow tails. Sign-changing ZOH extension: remove u>0 input restriction; replace the signed first delay moment in the slow-tail error by integral |u(s)|(T-s) ds. Keep the exact signed Q. Fast tail and compact derivatives are unchanged. Record source diff and verify old positive-waveform replay plus signed-waveform unit checks.

A finite grid generates candidates only. Inconclusive candidates may be refined by uniform numerical margin ladder [0.005, 0.001, 0.0002] and add violated relaxation times (no data/tolerance changes). Every attempt is logged; reported outcome is proved exclusion, constructive paired compatibility, or unresolved. For factors with no positive dual gap, seek a nonnegative finite atomic model matching the pulse endpoint; verify EIS and the pulse against a recorded envelope using directed intervals. Never turn a failed certificate into a compatibility assertion.

Domain ladder shifts [0,1,2,4,8], Taylor then mean-value with shift+1. Failure/depth/runtime limits do not mean acceptance. Independently check objectives, finite constraints and both tails using 75-decimal-digit interval arithmetic. Check exact dyadic partition coverage. Fresh solver-free replay. Retain negative controls, malformed-input controls, witnesses, processed instances, code, source hashes and all non-exclusions.

## Scope
Three independent physical source cells; 12 operating-condition pairs, not 12 independent cells. Comparison is nominal-state, delayed, producer-reconciled, model/processed-input conditional. No actual calibration, causal mechanism, EV pack validity or EIS-only stopping. No author approval, submission, public correspondence release or GitHub upload inferred.
