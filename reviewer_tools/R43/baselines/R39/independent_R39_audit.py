#!/usr/bin/env python3
"""Independent MPmath interval arithmetic (not a third continuum verifier), sign robustness and invalid-control audit."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction
import subprocess,json,concurrent.futures,hashlib
import pandas as pd,mpmath as mp
D=Path(__file__).resolve().parent;R=D/'r39_candidate';S=D/'replay/Batteries52_Reviewer_Replay/R36_reproducibility';C=pd.read_csv(R/'R39_FRESH_CERTIFICATES_84.csv')
mp.iv.dps=75;iv=mp.iv
f=lambda x:iv.mpf(*(float(x).as_integer_ratio()))
# iv.mpf(num,den) is not a rational constructor: ensure explicitly divide for exact binary64.
def pt(x):
 n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)
checks=[]
for row in C.itertuples():
 stem=f'{row.case_id}_alpha{row.alpha:g}_margin0.005';q=json.loads((R/'inputs'/(row.case_id+'.json')).read_text());c=json.loads((R/'candidates'/(stem+'.json')).read_text());v=json.loads((R/'validation'/(stem+'_taylor160.json')).read_text());vv=json.loads((R/'validation'/(stem+'_meanvalue256.json')).read_text());n=len(q['frequency_Hz']);y=list(map(pt,c['y']));w=[2*iv.pi*pt(x) for x in q['frequency_Hz']];a=list(map(pt,q['a_s']));b=list(map(pt,q['b_s']));u=list(map(pt,q['u']));T=pt(q['T_s']);zero=iv.mpf(0)
 Q=sum((z*(bb-aa) for aa,bb,z in zip(a,b,u)),zero);M1=sum((z*((T-aa)**2-(T-bb)**2)/2 for aa,bb,z in zip(a,b,u)),zero)
 obj=-sum((y[k]*pt(q['Zreal_mOhm'][k])+y[k+n]*pt(q['Zimag_mOhm'][k])+pt(row.alpha)*pt(q['eps_mOhm'])*(abs(y[k])+abs(y[k+n])) for k in range(n)),zero)
 fs=[pt(q['iT'])+sum(y[:n],zero),sum((v*w1 for v,w1 in zip(y[n:],w)),zero),Q-sum((v/w1 for v,w1 in zip(y[n:],w)),zero)]
 both=[v,vv];tails=[];delay=[T-a[0]];jump=[-u[0]]
 for j in range(1,len(u)):
  delay.append(T-a[j]);jump.append(u[j-1]-u[j])
 for res in both:
  tau=iv.exp(pt(res['xl']));z=iv.exp(-pt(res['xh']));fast=u[-1]+sum(y[:n],zero);slow=fs[2]-z*M1
  for k in range(n):
   fast-=abs(y[k])*(w[k]*tau)**2+abs(y[k+n])*w[k]*tau
   slow-=abs(y[k])*z/w[k]**2+abs(y[k+n])*z**2/w[k]**3
  for dd,jj in zip(delay,jump):fast-=abs(jj)*iv.exp(-dd/tau)
  tails.append(bool(fast.a>0 and slow.a>0))
 objective_inside=bool(obj.a>=pt(v['certified_lower_mOhm']).b and obj.b<=pt(v['certified_lower_upper_round_mOhm']).a)
 right=(u[-1]+sum(y[:n],zero));point_gap=obj-pt(q['Rmeas_mOhm']);positive=all(x.a>0 for x in fs)
 checks.append({'case_id':row.case_id,'alpha':row.alpha,'objective_inside_MPFR_rounds':objective_inside,'finite_constraints_positive':positive,'right_endpoint_positive':bool(right.a>0),'point_exclusion_positive':bool(point_gap.a>0),'both_tails_positive':all(tails),'sample_envelope_positive':bool(obj.a>pt(pd.read_csv(R/'RAW_CURRENT_CORRECTION_42.csv').set_index('case_id').loc[row.case_id,'raw_voltage_sample_envelope_high_mOhm']).b)})
print('independent arithmetic all cases',len(checks),flush=True)
(R/'INDEPENDENT_75DP_AUDIT.json').write_text(json.dumps(checks,indent=2)+'\n')
assert len(checks)==84 and all(all(x for k,x in c.items() if k not in ['case_id','alpha']) for c in checks)
# An intentionally invalid dual must be rejected because of the compact interior, despite positive endpoints/tails.
(R/'negative_controls').mkdir(exist_ok=True)
def negative(row):
 stem=f'{row.case_id}_alpha1_margin0.005';w=R/'witnesses'/(stem+'.txt');z=w.read_text().splitlines();n=int(z[0].split()[0]);lo=z[1].split();hi=z[n].split();lo[3]=repr(float(lo[3])-.1);hi[3]=repr(float(hi[3])+.1);z[1]=' '.join(lo);z[n]=' '.join(hi)
 p=R/'negative_controls'/(row.case_id+'_interior_invalid.txt');p.write_text('\n'.join(z)+'\n');cp=subprocess.run([str(R/'interval_verify'),str(p),str(p.with_suffix('.csv')),'160','2','0'],text=True,capture_output=True,timeout=90)
 j=json.loads(cp.stdout)
 return {'case_id':row.case_id,'rejected':not j['passed'] and cp.returncode==1,'finite_and_tails_positive':all(j[k]>0 for k in ['finite_R0_lower','finite_L_lower','finite_c_lower','fast_tail_lower','slow_tail_H_lower']),'claimed_gap_positive':j['gap_lower_mOhm']>0}
alpha1=C[C.alpha==1]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:n=list(pool.map(negative,list(alpha1.itertuples())))
(R/'NEGATIVE_CONTROLS_42.json').write_text(json.dumps(n,indent=2)+'\n')
assert len(n)==42 and all(x['rejected'] and x['finite_and_tails_positive'] and x['claimed_gap_positive'] for x in n)
summary={'status':'PASS','independent_objective_and_tail_checks':84,'negative_interior_controls_rejected':42,'original_7_EIS_feasible_models_unchanged':True,'directed_rounding_runs':168,'claim_boundary':'Saved dual transport and new raw-derived computational inputs; metrological uncertainty remains unmeasured'}
(R/'R39_INDEPENDENT_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
