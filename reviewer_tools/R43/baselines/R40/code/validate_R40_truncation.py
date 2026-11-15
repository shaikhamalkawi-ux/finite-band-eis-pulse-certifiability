#!/usr/bin/env python3
"""Independent 75-decimal interval objective/tail and exact leaf-coverage audit for 4 R40 truncated-band proofs."""
import json, pathlib, csv, math
from fractions import Fraction
import pandas as pd
import mpmath as mp
R=pathlib.Path(__file__).resolve().parent.parent/'results'; C=R/'trunc25';cases=['SOC70_charge_1.4C','SOC70_discharge_1.4C','SOC70_charge_1.0C','SOC60_discharge_1.4C']
mp.iv.dps=75;iv=mp.iv

def pt(x):
 n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)

def iv_gt0(x):return bool(x.a>0)
res=[]
for k in cases:
 q=json.loads((C/(k+'_q.json')).read_text());cand=json.loads((C/(k+'_candidate.json')).read_text()); n=len(q['frequency_Hz']);w=[2*iv.pi*pt(x) for x in q['frequency_Hz']]
 y=[pt(v) for v in cand['y']];alpha=pt(1.);eps=pt(q['eps_mOhm']);
 a=list(map(pt,q['a_s']));b=list(map(pt,q['b_s']));u=list(map(pt,q['u']));T=pt(q['T_s']);zero=iv.mpf(0)
 Q=sum((v*(bb-aa) for aa,bb,v in zip(a,b,u)),zero);M1=sum((v*((T-aa)**2-(T-bb)**2)/2 for aa,bb,v in zip(a,b,u)),zero)
 obj=-sum((y[j]*pt(q['Zreal_mOhm'][j])+y[j+n]*pt(q['Zimag_mOhm'][j])+alpha*eps*(abs(y[j])+abs(y[j+n])) for j in range(n)),zero)
 finite=[pt(q['iT'])+sum(y[:n],zero),sum((v*omega for v,omega in zip(y[n:],w)),zero),Q-sum((v/omega for v,omega in zip(y[n:],w)),zero)]
 for tag in ['taylor160','meanvalue256']:
  z=json.loads((C/(k+'_'+tag+'_verification.json')).read_text());leafs=list(csv.DictReader((C/(k+'_'+tag+'_leaves.csv')).open()));
  assert z['passed'] and len(leafs)==z['accepted_leaves'] and iv_gt0(obj-pt(q['Rmeas_mOhm']))
  assert all(float(x['residual_lower'])>0 for x in leafs)
  leafs.sort(key=lambda x:float(x['x_lo']))
  x=Fraction(float(z['xl']))
  for row in leafs:
   lo=Fraction(float(row['x_lo']));hi=Fraction(float(row['x_hi']));assert lo==x and hi>lo;x=hi
  assert x==Fraction(float(z['xh']))
  assert iv_gt0(obj-pt(z['certified_lower_mOhm'])) or bool(obj.a>=pt(z['certified_lower_mOhm']).b)
  assert bool(obj.b<=pt(z['certified_lower_upper_round_mOhm']).a)
  tau=iv.exp(pt(z['xl']));zmax=iv.exp(-pt(z['xh']))
  fast=u[-1]+sum(y[:n],zero);slow=finite[2]-zmax*M1
  for j in range(n):
   fast-=abs(y[j])*(w[j]*tau)**2 + abs(y[j+n])*w[j]*tau
   slow-=abs(y[j])*zmax/w[j]**2 + abs(y[j+n])*zmax**2/w[j]**3
  delay=[T-a[0]]+[T-a[j] for j in range(1,len(u))];jump=[-u[0]]+[u[j-1]-u[j] for j in range(1,len(u))]
  for d,jmp in zip(delay,jump):fast-=abs(jmp)*iv.exp(-d/tau)
  passed=all(iv_gt0(v) for v in finite+[fast,slow])
  assert passed,(k,tag)
  res.append({'case_id':k,'method':tag,'n_leaves':len(leafs),'independent_75dp_objective_encloses_MPFR':True,'finite_R0_L_c_positive':True,'fast_and_slow_tails_positive':True,'compact_leaves_contiguous_positive':True,'measured_gap_positive':True})
D=pd.DataFrame(res);D.to_csv(R/'R40_Truncated_Independent_75DP_Audit_8.csv',index=False)
assert len(D)==8 and D.select_dtypes(bool).all().all()
print('R40 25.31mHz 4 x 2 independent arithmetic, interval tail and leaf coverage: PASS',flush=True)
