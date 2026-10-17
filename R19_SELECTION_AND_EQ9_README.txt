52 Batteries R19 — selection-screen and Eq. (9) additions

1. Historical exploratory selection screen
Files beginning R19_Historical_Exploratory_* are derived from the preserved R3 full-state exploratory gate that was used before the final continuum-searched primary analysis was fixed. That historical screen used an earlier finite-relaxation-grid / ideal-step analysis and is included ONLY to document why 5 °C was selected. It is not used as final continuum evidence and does not replace the R4D/R9/R10 analysis.

Across interior SOC 20–80%, the historical screen's median extra componentwise tolerance was largest at 5 °C at every SOC. The aggregate temperature medians were approximately 2.320, 0.709, 0.207, 0.010 and 0 mΩ at 5, 15, 25, 35 and 45 °C, respectively.

2. Eq. (9) numerical connection
R19_RWTH_5C_Eq9_SlowMode_Radius_42cases.csv evaluates 2*pi*epsilon*f_min*|Q_T| for the 42 primary 5 °C measured-current waveforms using f_min=0.01 Hz and the primary optimization-derived EIS admission floor delta_EIS as the stated common componentwise tolerance. This is not an instrument-calibrated uncertainty claim.

Across the 42 cases the resulting asymptotic slow-mode lower-bound scale is 0.367–0.428 mΩ, median 0.376 mΩ. Eq. (9) is a necessary lower bound on unavoidable ambiguity, not an upper bound on the complete EIS-compatible response range; continuum optimization remains necessary.

The replay script requires the preserved official MatFiles.zip and does not redistribute that raw archive in this compact package.
