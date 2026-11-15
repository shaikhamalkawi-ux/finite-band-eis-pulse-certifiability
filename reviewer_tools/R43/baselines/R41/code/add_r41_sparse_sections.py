#!/usr/bin/env python3
"""Bounded R41 research-manuscript update: compressed continuous dual witnesses.
Operates on R41 finite-atom candidate manuscript source. Not a public release.
"""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
main=ROOT/'overleaf'/'main.tex'
sup=ROOT/'overleaf'/'supplementary.tex'
s=main.read_text()

def change(old,new,label):
 global s
 assert s.count(old)==1,(label,s.count(old));s=s.replace(old,new)

# A coherent abstract within IEEE conventions. This preserves previous quantitative conclusions.
abst=s.split('\\begin{abstract}',1)[1].split('\\end{abstract}',1)[0]
new_abst=r'''
An electrochemical impedance spectroscopy (EIS) fit does not establish that finite-band data support a specified battery pulse-response claim. We test certifiability over positive resistor--capacitor (RC)/Stieltjes models without selecting a unique circuit or predictor. A worst-case information test uses both an asymptotically slow relaxation and a finite-time atom; continuum dual witnesses then test recorded-EIS/pulse compatibility. The finite-time witness strengthens the median observation-uniform lower bound from 0.375 to 0.427~m$\Omega$ over 42 corrected-current waveform configurations. In the RWTH Aachen 5~$^\circ$C panel, all 42 measured pulse responses fall below the EIS-compatible interval (median nominal gap 5.984~m$\Omega$). Correct measured-current normalization resolves three tapered charging cases. Directed-rounding certificates verify exclusion at the EIS model-admission tolerance and its double, with minimum fixed-current recorded-sample gaps 1.2874 and 0.0608~m$\Omega$. Retrospective zero-dual-weight pruning retains each 48-frequency exclusion certificate using only 15--21 frequencies, or a shared 24--30 per state of charge (SOC) for all six pulse conditions, under the original full-band tolerance; 336 reduced-input MPFR verification runs pass. A twelve-condition check on three Bayreuth cells retains exclusion at both source-specific tolerances after recorded-sample expansion, conditional on source-reconciled pairing and state persistence. Historical local $3\sigma$ bounds remain conditional and do not establish calibration. The results show model-set inconsistency and retrospective certificate compressibility, not population inference, EIS-only stopping, or measured acquisition-time savings.
'''
change(abst,new_abst,'abstract')

# Following introduction and generalization of the method, include sparse support appropriately.
change(r'The study contributes a finite-band necessary information criterion with a closed-form finite-time single-atom refinement, continuum response bounds, and a paired-data model-set diagnosis;',
       r'The study contributes a finite-band necessary information criterion with a closed-form finite-time single-atom refinement, continuum response bounds, paired-data model-set diagnosis, and lossless dual-support projection of those certificates;', 'contribution enumeration')

change(r'''The certified upper endpoint is obtained by minimizing the right-hand side of Eq.~\eqref{eq:dualbound} over nonnegative multipliers satisfying Eq.~\eqref{eq:dualfeas}; applying the same construction to $-F_T$ gives the lower endpoint.''',
       r'''The certified upper endpoint is obtained by minimizing the right-hand side of Eq.~\eqref{eq:dualbound} over nonnegative multipliers satisfying Eq.~\eqref{eq:dualfeas}; applying the same construction to $-F_T$ gives the lower endpoint. A further consequence is \emph{lossless witness-support projection}: if both real and imaginary dual coefficients at a given frequency are exactly zero, that frequency contributes neither to the objective nor to the continuum inequality. Removing it while holding the original EIS tolerance fixed therefore preserves the same sufficient dual certificate. This does not prove that an a priori measurement schedule can identify the retained frequencies.''','dual support math')

# One informative main-table and results paragraph, separate from contiguous truncation.
anchor=r'''\subsection{Historical temperature slice at 50\% SOC (exploratory)}'''
addition=r'''\subsection{Certificate-support compression at fixed EIS tolerance}
The contiguous low-frequency truncations above address a different question from how many EIS coordinates an \emph{already verified} exclusion certificate actually uses. Each 48-frequency lower dual has separate real and imaginary multipliers at each frequency. We removed a frequency only when \emph{both} multipliers were exactly zero in the archived witness, without re-optimizing the model, altering the waveform, or changing the full-band EIS-admission tolerance. All 84 certificates retained their exact mathematical objective and continuum residual using only 15--21 frequency points each (median 18). To permit a common per-SOC EIS subset for all six charge/discharge pulse conditions and both tolerance factors, we also took the union of their twelve supports, producing seven 24--30-frequency subsets (Table~\ref{tab:sparse}). These are sufficient subsets for the \emph{retrospective paired-reference negative certificates}; they are neither proven minimal nor universal.

\begin{table}[!t]
\centering
\caption{Sufficient SOC-specific common frequency support for all six pulse cases at both frozen tolerance factors. Original spectra have 48 frequency points.}
\label{tab:sparse}
\small
\begin{tabular}{rrrrrrrr}
\toprule
SOC (\%) & 20 & 30 & 40 & 50 & 60 & 70 & 80 \\
\midrule
Retained points & 24 & 30 & 29 & 25 & 27 & 29 & 28 \\
\bottomrule
\end{tabular}
\end{table}

Both the per-certificate projection and the per-SOC common-subset projection were replayed with independent 160-bit Taylor and 256-bit mean-value full-domain interval enclosures: 84 certificates $\times$ two projections $\times$ two verification methods $=336$ successful directed-rounding runs. All recorded conservative lower bounds agree numerically with their original witnesses, and the minimum doubled-tolerance point gap remains $0.076987\,\mohm$. The original seven feasible full-band EIS models also remain feasible after removing constraints. As all seven SOC-common sets still retain the lowest-frequency measurement, and the retained sets were chosen \emph{after} viewing full EIS and pulse evidence, this does not demonstrate prospective selection, acquisition-time savings, optimal sparsity, or an EIS-only early-stopping policy. Full index-level records and clean verifier inputs are in Supplementary Section~S15 and the accompanying computational archive.

'''
change(anchor,addition+anchor,'sparse result')

# One sentence in discussion to integrate without repeating limitations many times.
old=r'''The truncation study adds an engineering qualification to this result.'''
new=r'''The truncation and dual-support studies add different engineering qualifications to this result. Truncation removes a contiguous low-frequency band and requires a new optimization and continuum check; lossless support projection instead removes exact-zero-weight coordinates from \emph{already certified} witnesses at a frozen EIS tolerance. The former probes a possible low-frequency acquisition tradeoff; the latter establishes retrospective certificate compression, not a prospective reduction in measurement duration. The truncation study adds a further qualification.'''
change(old,new,'discussion differentiation')

# Concise conclusion without inflating evidence.
old=r'''The practical consequence is a retrospective, reference-conditional method for identifying which truncated EIS bands still support a negative cross-domain consistency certificate.'''
new=r'''The practical consequence is a retrospective, reference-conditional method for identifying which truncated EIS bands still support a negative cross-domain consistency certificate. In addition, exact-zero dual-coordinate pruning preserves all 84 verified lower certificates with 15--21 frequency points individually or 24--30 common points per SOC, at the full-band tolerance. This is a statement about certificate support, not validated acquisition acceleration.'''
change(old,new,'conclusion compression')

main.write_text(s)

u=sup.read_text()
assert u.endswith('\\end{document}\n')
soc=[(20,24),(30,30),(40,29),(50,25),(60,27),(70,29),(80,28)]
rows='\n'.join(f'{a} & {b} & {48-b} \\\\' for a,b in soc)
extra=r'''\section*{S15. Lossless support projection of the 84 corrected-current continuum certificates}
The original lower certificates use a 48-frequency EIS spectrum at each SOC, giving 96 real and imaginary coordinate multipliers per dual witness. Let $y^R_k$ and $y^I_k$ denote the stored signed multipliers. Let $J=\{k:(y^R_k,y^I_k)\ne(0,0)\}$ be the set of EIS frequencies appearing with a nonzero dual coefficient. For any omitted $k\notin J$, both multipliers are exactly zero in the serialized witness. All terms involving frequency $k$ in the dual objective, finite-coordinate inequalities, and continuous relaxation-kernel inequality therefore vanish identically. Deleting those rows while preserving all other coefficients and the full-band analysis tolerance leaves exactly the same \emph{mathematical} certificate for the relaxed EIS data subset. No new optimizer, reduced-band model-admission tolerance, or re-estimated current waveform is required. Because dropping EIS constraints cannot destroy an already feasible original full-band primal model, nonvacuity carries over as well.

For each of the 84 original case--tolerance witnesses, exact-zero-weight removal retains 15--21 of 48 frequency indices (median 18). Each such per-certificate support can differ across pulse conditions. Since all six pulse conditions at one SOC share an EIS spectrum, a practical common-support question is to take the union of the 12 supports at that SOC (six pulse conditions at two tolerance factors). Table~\ref{tab:sparse_support} lists the resulting seven common sets. These sets do not require measurements at every full-band frequency to \emph{verify the already selected} duals, but their selection used full-data optimization and independent pulse results and hence does not establish a prospective EIS-only subset-selection policy.

\begin{table}[!ht]
\centering
\caption{Retrospective frozen-tolerance common supports for all twelve lower certificates at each SOC.}
\label{tab:sparse_support}
\begin{tabular}{rrr}
\toprule
SOC (\%) & Frequencies retained & Frequencies removed \\
\midrule
''' + rows + r'''
\bottomrule
\end{tabular}
\end{table}

The accompanying archive contains: (i) all 84 original full-band witness text files; (ii) the 84 per-certificate frequency-projected witnesses; (iii) 84 further witnesses using the SOC-common union sets; (iv) the exact retained frequency indices and frequency values; (v) the seven original feasible EIS primal witnesses; and (vi) results of a fresh MPFR replay of every compressed witness using both the Taylor160 and the independent mean-value256 domain partition method. All $84\times2\times2=336$ verifier invocations passed. The maximum recorded deviation of the certified lower bound from the original R39 value is zero at exported precision, and the minimum positive double-tolerance point gap is $0.0769875208573$~m$\Omega$. These figures were generated under the fixed original full-spectrum EIS tolerance. They are not an optimization of minimum required frequencies; some additional frequency rows may be removable by finding \emph{new} dual witnesses, but that is not claimed here.

The same clean verifier automatically confirms the finite $R_0$, $L$, and $c$ coordinates, positive residual lower bounds over every positive relaxation time, analytic fast and slow tails, and positivity of the final measured-response separation. Separately verifying the stored full-band feasible primal witnesses establishes that the original data model class is nonempty at each SOC. These statements concern exact processed numerical inputs rather than metrological calibration. The lowest measured frequency remains in every SOC-common union, so this support compression cannot be translated directly into shorter low-frequency acquisition or tested online stopping. Only a separately designed prospective measurement protocol could evaluate those applications.

'''
u=u.replace('\\end{document}\n',extra+'\\end{document}\n')
sup.write_text(u)
print('abstract words',len(new_abst.split()),'main size',len(s),'supp size',len(u))
