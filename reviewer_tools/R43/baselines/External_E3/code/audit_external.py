#!/usr/bin/env python3
"""Independent endpoint/objective/tail arithmetic, partitions, and negative controls.
This is a separate arithmetic implementation, not an independent human review.
"""
from pathlib import Path
from fractions import Fraction
import json,subprocess,concurrent.futures,copy
import numpy as np,pandas as pd,mpmath as mp,sympy as sp
R=Path(__file__).resolve().parents[1];mp.iv.dps=75;iv=mp.iv

def pt(x):n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)
def pos(x):return bool(x.a>0)
def save(p,x):p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def independent(q,c,v,v2):
 n=len(q['frequency_Hz']);y=list(map(pt,c['y']));w=[2*iv.pi*pt(f) for f in q['frequency_Hz']];a=list(map(pt,q['a_s']));b=list(map(pt,q['b_s']));u=list(map(pt,q['u']));T=pt(q['T_s']);zero=iv.mpf(0)
 Q=sum((vv*(bb-aa) for aa,bb,vv in zip(a,b,u)),zero);M1=sum((abs(vv)*((T-aa)**2-(T-bb)**2)/2 for aa,bb,vv in zip(a,b,u)),zero)
 lb=-sum((y[k]*pt(q['Zreal_mOhm'][k])+y[k+n]*pt(q['Zimag_mOhm'][k])+pt(c['alpha'])*pt(q['eps_mOhm'])*(abs(y[k])+abs(y[k+n])) for k in range(n)),zero)
 fs=[pt(q['iT'])+sum(y[:n],zero),sum((yy*ww for yy,ww in zip(y[n:],w)),zero),Q-sum((yy/ww for yy,ww in zip(y[n:],w)),zero)]
 right=pt(q['iT_right_recorded'])+sum(y[:n],zero)
 delay=[T-a[0]];jump=[-u[0]]
 for j in range(1,len(u)):
  if q['u'][j-1]!=q['u'][j]:delay.append(T-a[j]);jump.append(u[j-1]-u[j])
 tail=[]
 for j in [v,v2]:
  tau=iv.exp(pt(j['xl']));z=iv.exp(-pt(j['xh']));fast=u[-1]+sum(y[:n],zero);slow=fs[2]-z*M1
  for k in range(n):
   fast-=abs(y[k])*(w[k]*tau)**2+abs(y[k+n])*w[k]*tau
   slow-=abs(y[k])*z/(w[k]**2)+abs(y[k+n])*z**2/(w[k]**3)
  for dd,jj in zip(delay,jump):fast-=abs(jj)*iv.exp(-dd/tau)
  tail.append({'xl':j['xl'],'xh':j['xh'],'fast_interval':str(fast),'slow_H_interval':str(slow),'pass':pos(fast) and pos(slow)})
 return {'case_id':q['case_id'],'alpha':c['alpha'],'objective_interval':str(lb),'objective_inside_MPFR':bool(lb.a>=pt(v['certified_lower_mOhm']).b and lb.b<=pt(v['certified_lower_upper_round_mOhm']).a),'finite_constraints_positive':all(map(pos,fs)),'right_endpoint_positive':pos(right),'both_tails_positive':all(t['pass'] for t in tail),'tails':tail,'recorded_envelope_gap_interval':str(lb-pt(q['recorded_envelope_upper_mOhm'])),'recorded_envelope_excluded':pos(lb-pt(q['recorded_envelope_upper_mOhm'])),'negative_samples':sum(x<0 for x in q['u'])}

def negative_one(row):
 q=json.loads((R/'inputs'/(row.case_id+'.json')).read_text());c=json.loads((R/row.candidate_file).read_text());n=len(q['frequency_Hz']);c['y'][0]-=.1;c['y'][n-1]+=.1
 from run_external_certificates import export
 p=R/'qa/negative_controls'/(row.case_id+'.txt');p.parent.mkdir(exist_ok=True,parents=True);export(q,c,p)
 cp=subprocess.run([str(R/'code/interval_verify_signed'),str(p),str(p.with_suffix('.csv')),'160','2',str(int(row.tail_shift))],text=True,capture_output=True,timeout=120)
 vv=json.loads(cp.stdout);save(p.with_suffix('.json'),vv)
 # Independent interval evaluation at a point: certify actual invalidity, not merely numerical inconclusiveness.
 om=2*np.pi*np.array(q['frequency_Hz']);xx=np.linspace(-8,8,2001);tau=np.exp(xx);wt=om[:,None]*tau;yn=np.array(c['y']);res=yn[:n]@(1/(1+wt*wt))+yn[n:]@(-wt/(1+wt*wt))
 for aa,bb,u in zip(q['a_s'],q['b_s'],q['u']):res+=u*np.exp(-(q['T_s']-bb)/tau)*(-np.expm1(-(bb-aa)/tau))
 xstar=float(xx[int(res.argmin())]);tv=iv.exp(pt(xstar));r=iv.mpf(0)
 for k in range(n):
  wt=2*iv.pi*pt(q['frequency_Hz'][k])*tv;r+=(pt(yn[k])-pt(yn[k+n])*wt)/(1+wt*wt)
 for a,b,u in zip(q['a_s'],q['b_s'],q['u']):r+=pt(u)*iv.exp(-(pt(q['T_s'])-pt(b))/tv)*(-iv.expm1(-(pt(b)-pt(a))/tv))
 return {'case_id':row.case_id,'verifier_rejected':not vv['passed'],'finite_and_tails_pass':all(vv[k]>0 for k in ['finite_R0_lower','finite_L_lower','finite_c_lower','fast_tail_lower','slow_tail_H_lower']),'claimed_gap_positive':vv['gap_lower_mOhm']>0,'negative_point_logtau':xstar,'negative_point_interval':str(r),'negative_point_upper_lt_zero':bool(r.b<0)}

if __name__=='__main__':
 d=pd.read_csv(R/'results/EXTERNAL_CERTIFICATES_24.csv');assert len(d)==24
 audits=[];covers=[]
 for _,r in d.iterrows():
  if r.status!='EXCLUDED_INTERVAL_VERIFIED':continue
  q=json.loads((R/'inputs'/(r.case_id+'.json')).read_text());c=json.loads((R/r.candidate_file).read_text());stem=Path(r.witness_file).stem
  v=json.loads((R/'results/validation'/(stem+'_taylor160.json')).read_text());v2=json.loads((R/'results/validation'/(stem+'_meanvalue256.json')).read_text())
  audits.append(independent(q,c,v,v2))
  for suffix,j in [('_taylor160',v),('_meanvalue256',v2)]:
   f=pd.read_csv(R/'results/leaves'/(stem+suffix+'.csv')).sort_values('x_lo');lo=f.x_lo.values;hi=f.x_hi.values
   ok=(lo[0]==j['xl'] and hi[-1]==j['xh'] and np.array_equal(hi[:-1],lo[1:]) and np.all(hi>lo) and np.all(f.residual_lower>0))
   covers.append({'case_id':r.case_id,'alpha':r.alpha,'method':suffix,'leaves':len(f),'complete_dyadic_cover':bool(ok)})
  print('ARITHMETIC',r.case_id,r.alpha,all(audits[-1][k] for k in ['objective_inside_MPFR','finite_constraints_positive','right_endpoint_positive','both_tails_positive','recorded_envelope_excluded']),flush=True)
 save(R/'qa/INDEPENDENT_ARITHMETIC.json',audits);save(R/'qa/DOMAIN_COVERAGE.json',covers)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:neg=list(pool.map(negative_one,[r for _,r in d[d.alpha==1].iterrows() if r.status=='EXCLUDED_INTERVAL_VERIFIED']))
 save(R/'qa/NEGATIVE_CONTROLS.json',neg)
 sym=sp.symbols('z',positive=True);Dr=lambda e:sp.simplify(sym*sp.diff(e,sym));kr=1/(1+sym**2);ki=-sym/(1+sym**2)
 syms=[sp.simplify(Dr(kr)+2*sym**2/(1+sym**2)**2)==0,sp.simplify(Dr(Dr(kr))-4*sym**2*(sym**2-1)/(1+sym**2)**3)==0,sp.simplify(Dr(ki)-sym*(sym**2-1)/(1+sym**2)**2)==0,sp.simplify(Dr(Dr(ki))-sym*(-sym**4+6*sym**2-1)/(1+sym**2)**3)==0]
 checks={'24_unique_case_factor_pairs':len(d)==24 and not d.duplicated(['case_id','alpha']).any(),'24_double_method_certificates':bool((d.status=='EXCLUDED_INTERVAL_VERIFIED').all()),'objective_independent_interval':all(a['objective_inside_MPFR'] for a in audits),'finite_and_both_endpoints':all(a['finite_constraints_positive'] and a['right_endpoint_positive'] for a in audits),'independent_both_tail_bounds':all(a['both_tails_positive'] for a in audits),'recorded_envelopes_excluded':all(a['recorded_envelope_excluded'] for a in audits),'48_complete_interval_partitions':len(covers)==48 and all(c['complete_dyadic_cover'] for c in covers),'12_invalid_controls_rejected':len(neg)==12 and all(n['verifier_rejected'] for n in neg),'12_invalid_points_certified_negative':len(neg)==12 and all(n['negative_point_upper_lt_zero'] for n in neg),'negative_controls_isolate_compact_interior':all(n['finite_and_tails_pass'] and n['claimed_gap_positive'] for n in neg),'symbolic_kernel_derivatives':all(syms)}
 save(R/'qa/AUDIT_SUMMARY.json',{'checks':checks,'passed':sum(checks.values()),'total':len(checks),'overall_pass':all(checks.values())});print(json.dumps(checks,indent=2),flush=True)
