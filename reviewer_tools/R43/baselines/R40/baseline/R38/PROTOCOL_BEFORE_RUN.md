# All-42 proof-oriented validation: execution protocol
Frozen before candidate generation (19 October 2026 project date).
Baseline R35; no manuscript changes before verified computational evidence.
All 42 primary 5 deg C rows, SOC 20:10:80, both directions, rates 0.7,1,1.4 C.
Source raw MatFiles.zip SHA-256: 875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f.
Use frozen R20 measured-waveform normalization (raw voltage endpoint and locked response to recover current scale); audit all mappings.
Scientific tolerances: locked primary eps and previously reported 2*eps; never optimize tolerance for the new certificate.
Candidate search: floating point LP with strictly positive witness slack. Slack is only a conservative numerical device and is not an altered measurement tolerance.
Search slack ladder 0.02,0.01,0.005,0.002,0.001,0.0005; use the first certificate that proves separation. Store attempts, including failures.
Verification must be independent of the LP status: directed-rounding interval checks for finite R0,L,c and continuum residual over ALL relaxation times; explicit analytic fast and slow tails.
Try interval natural extension, mean-value/Taylor enclosure and increased precision/subdivision as needed; random/dense meshes alone cannot confer PASS.
Recompute lower bound and margin with outward rounding from saved coefficients.
Do not claim interval endpoint optimality, instrument calibration, online stopping, 42 independent cells, or strengthened temperature coverage.
Negative controls: sign-reversed or perturbed invalid certificates must be rejected.
Record all 42 statuses, conservative lower bounds, precision-expanded gaps, covered domains and numerical witness hashes.
