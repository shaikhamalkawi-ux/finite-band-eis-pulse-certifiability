#!/usr/bin/env python3
"""R41: Exact single-frequency atom witness for common EIS componentwise tolerance.

The Proposition is WORST CASE over synthetic EIS observation centers, NOT a
lower bound on the width of the particular observed RWTH spectral box. The
41/42 etc MUST NOT be read as rejection counts; only R39 MPFR certificates do so.

Inputs: immutable 42 serialized, measured-current-normalized RWTH cases from R39.
No third-party raw archive or empirical pulse resampling is performed here.
"""
import csv, json, math, pathlib, hashlib, statistics
import numpy as np
import mpmath as mp
from decimal import Decimal, ROUND_FLOOR
mp.mp.dps = 80
ROOT = pathlib.Path(__file__).resolve().parents[1]
IN = ROOT / 'evidence' / 'inputs'; OUT = ROOT / 'evidence'
C = lambda x: mp.mpf(str(x))

def zoh_response(q, tau):
    # Compute exactly as given by archived decimal input coefficients using high precision.
    T = C(q['T_s']); g = mp.mpf(0)
    for aa, bb, uu in zip(q['a_s'], q['b_s'], q['u']):
        a,b,u=C(aa),C(bb),C(uu)
        g += u*mp.exp(-(T-b)/tau)*(-mp.expm1(-(b-a)/tau))
    return g

def zoh_response_independent(q,tau):
    T=C(q['T_s']); g=mp.mpf(0)
    for aa,bb,uu in zip(q['a_s'],q['b_s'],q['u']):
        a,b,u=C(aa),C(bb),C(uu)
        g+=u*(mp.exp(-(T-b)/tau)-mp.exp(-(T-a)/tau))
    return g

def hidden_atom_cap(eps, fmin, tau):
    w = 2*mp.pi*fmin; x=w*tau
    return eps*(1+x*x) if x<=1 else eps*(x+1/x)

def original_atom_cap(eps, freqs, tau):
    w=[2*mp.pi*C(f) for f in freqs]
    return min(min(eps*(1+(v*tau)**2), eps*(1+(v*tau)**2)/(v*tau)) for v in w)

def run():
    paths=sorted(IN.glob('*.json')); assert len(paths)==42, len(paths)
    rows=[]; ids=set(); tests={'inputs':len(paths),'zoh_independent_max_abs':mp.mpf(0),'closed_rhide_max_abs':mp.mpf(0),
                           'tau_star_atom_check_max_abs':mp.mpf(0),'charge_decomposition_max_abs':mp.mpf(0)}
    for file in paths:
        q=json.loads(file.read_text(),parse_float=str)
        case=q['case_id']; assert case not in ids; ids.add(case)
        assert q['current_reference'].startswith('absolute median measured pulse current')
        freqs=list(q['frequency_Hz']); fmin=min(map(C,freqs)); eps=C(q['eps_mOhm']); w=2*mp.pi*fmin; tau=1/w
        assert C(q['T_s'])==20
        # At tau_* all EIS real/imag weights are <=1/2; min-frequency has both exactly 1/2.
        point_mass=original_atom_cap(eps,freqs,tau)
        tests['tau_star_atom_check_max_abs']=max(tests['tau_star_atom_check_max_abs'],abs(point_mass-2*eps))
        g=zoh_response(q,tau); gd=zoh_response_independent(q,tau)
        tests['zoh_independent_max_abs']=max(tests['zoh_independent_max_abs'],abs(g-gd))
        Q=sum((C(uu)*(C(bb)-C(aa)) for aa,bb,uu in zip(q['a_s'],q['b_s'],q['u'])),mp.mpf(0))
        tests['charge_decomposition_max_abs']=max(tests['charge_decomposition_max_abs'],abs(Q-C(q['Q_T_float_s'])))
        slow=eps*w*abs(Q); finite=2*eps*abs(g)
        assert slow<mp.mpf('.4') and finite>mp.mpf('.4'),case
        assert finite>slow,case
        # Relative to synthetic-center worst case only, not to this case's EIS observation center.
        rows.append(dict(case_id=case,SOC=q['SOC'],direction=q['direction'],C_rate=q['C_rate'],
                         eps_mOhm=float(eps),fmin_Hz=float(fmin),tau_star_s=float(tau),
                         measured_Q_T_s=float(Q),unit_atom_g_tau=float(g),
                         R41_slow_limit_mOhm=float(slow),R41_finite_atom_onepoint_mOhm=float(finite),
                         R41_finite_over_slow=float(finite/slow),
                         worst_case_uniform_target_0p4_fails_slow=bool(slow>mp.mpf('.4')),
                         worst_case_uniform_target_0p4_fails_finite_atom=bool(finite>mp.mpf('.4')),
                         finite_atom_minus_0p4_mOhm=float(finite-mp.mpf('.4'))))
        # Independent randomized but deterministic checks on the complete 48-frequency direct min.
        rng=np.random.default_rng(20261115+len(ids))
        exponent_samples=[-8,-6,-4,-2,-1,0,1,2,3,5,7,9]+rng.uniform(-8,9,70).tolist()
        for expo in exponent_samples:
            test_tau=C(10)**C(expo)
            direct=original_atom_cap(eps,freqs,test_tau)
            closed=hidden_atom_cap(eps,fmin,test_tau)
            tests['closed_rhide_max_abs']=max(tests['closed_rhide_max_abs'],abs(direct-closed))
            assert mp.almosteq(direct,closed,rel_eps=mp.mpf('1e-60'),abs_eps=mp.mpf('1e-60')),(case,expo)
    rows.sort(key=lambda x:x['case_id'])
    csv_path=OUT/'R41_FiniteAtom_42Cases.csv'
    with csv_path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={
      'version':'R41 working, exploratory mathematical strengthening, not journal-release clearance',
      'source':'42 immutable R39 measured-current 5--15-second normalized JSON pulse/EIS cases',
      'model':'positive-RC/Stieltjes, identical box half-tolerance for all real/imag EIS coordinates',
      'theorem_scope':'worst-case over allowed EIS observation centers, not width of a particular observed center',
      'analytic_cap':'r_hide(eps,tau) = eps*(1+x*x) when x<=1 else eps*(x+1/x), x=2*pi*f_min*tau',
      'analytic_witness':'tau_star = 1/(2*pi*f_min), r_hide(tau_star)=2*eps',
      'N':42,'target_mOhm':0.4,
      'slow_screen_fail_count':sum(x['worst_case_uniform_target_0p4_fails_slow'] for x in rows),
      'finite_screen_fail_count':sum(x['worst_case_uniform_target_0p4_fails_finite_atom'] for x in rows),
      'slow_min_mOhm':min(x['R41_slow_limit_mOhm'] for x in rows),
      'slow_median_mOhm':statistics.median(x['R41_slow_limit_mOhm'] for x in rows),
      'slow_max_mOhm':max(x['R41_slow_limit_mOhm'] for x in rows),
      'finite_min_mOhm':min(x['R41_finite_atom_onepoint_mOhm'] for x in rows),
      'finite_median_mOhm':statistics.median(x['R41_finite_atom_onepoint_mOhm'] for x in rows),
      'finite_max_mOhm':max(x['R41_finite_atom_onepoint_mOhm'] for x in rows),
      'median_strengthening_pct':100*(statistics.median(x['R41_finite_atom_onepoint_mOhm'] for x in rows)/statistics.median(x['R41_slow_limit_mOhm'] for x in rows)-1),
      'tests':{k:str(v) for k,v in tests.items()},
      'tests_common_uniform_48_frequency':42*(12+70),
      'not_claimed':['specific RWTH measured EIS-box half-width > 0.4 mOhm', 'global B_atom maximization',
                     'instrument calibrated tolerance', 'new data','causal mechanism','prospective pulse design']}
    (OUT/'R41_Math_Reproduction_Summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # explicit non-uniform tolerance control: with a very precise higher frequency, result fails.
    x0=1;omega=[1,10];t=1;eps_real=[1,0.01];eps_imag=[1,0.01]
    cap=min([eps_real[i]*(1+(w*t)**2) for i,w in enumerate(omega)]+[eps_imag[i]*(1+(w*t)**2)/(w*t) for i,w in enumerate(omega)])
    assert cap<0.11 and 2*eps_real[0]==2
    summary['nonuniform_counterexample']={'omega':[1,10],'epsilon_each':[1,0.01], 'tau':t,
                                               'actual_rhide':cap,'false_fmin_only_prediction':2.0}
    (OUT/'R41_Math_Reproduction_Summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2)); return rows

if __name__=='__main__':run()
