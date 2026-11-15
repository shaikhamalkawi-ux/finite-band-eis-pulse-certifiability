#!/usr/bin/env python3
"""All external cases and both factors, same surplus; retain all outcomes."""
from pathlib import Path
import json,time,subprocess,concurrent.futures
import numpy as np,pandas as pd
from fractions import Fraction
from generate_candidates import solve
R=Path(__file__).resolve().parents[1]
for name in ['candidates','witnesses','results/validation','results/leaves']: (R/name).mkdir(exist_ok=True,parents=True)
def low_difference(a,b):
 exact=Fraction.from_float(float(a))-Fraction.from_float(float(b));v=float(exact)
 if Fraction.from_float(v)>exact:v=np.nextafter(v,-np.inf)
 return float(v)
def export(q,c,p):
 n=len(q['frequency_Hz']);y=c['y']
 lines=[' '.join(map(str,[n,len(q['u']),q['T_s'],q['iT'],q['eps_mOhm'],q['Rmeas_mOhm'],c['alpha']]))]
 for k in range(n):lines.append(' '.join(repr(float(v)) for v in [q['frequency_Hz'][k],q['Zreal_mOhm'][k],q['Zimag_mOhm'][k],y[k],y[k+n]]))
 for a,b,v in zip(q['a_s'],q['b_s'],q['u']):lines.append(' '.join(map(repr,[a,b,v])))
 p.write_text('\n'.join(lines)+'\n')
def one(path,alpha):
 q=json.loads(path.read_text());tt=time.time();attempts=[]
 result={'case_id':q['case_id'],'alpha':alpha,'eps_mOhm':q['eps_mOhm'],'measured_mOhm':q['Rmeas_mOhm'],'recorded_envelope_upper_mOhm':q['recorded_envelope_upper_mOhm'],'status':'UNRESOLVED'}
 for margin in [.005,.001,.0002]:
  c=solve(q,alpha,margin);stem=f"{q['case_id']}_alpha{alpha:g}_margin{margin:g}"
  (R/'candidates'/(stem+'.json')).write_text(json.dumps(c,indent=2))
  if not c['solver_success']:attempts.append({'margin':margin,'candidate':c});continue
  if c['floating_gap_mOhm']<=0:
   attempts.append({'margin':margin,'floating_gap_mOhm':c['floating_gap_mOhm'],'stage':'nonpositive_candidate_bound'});continue
  witness=R/'witnesses'/(stem+'.txt');export(q,c,witness)
  for shift in [0,1,2,4,8]:
   cp=subprocess.run([str(R/'code/interval_verify_signed'),str(witness),str(R/'results/leaves'/(stem+'_taylor160.csv')),'160','2',str(shift)],text=True,capture_output=True,timeout=180)
   if not cp.stdout.strip():attempts.append({'margin':margin,'shift':shift,'stderr':cp.stderr,'returncode':cp.returncode});continue
   v=json.loads(cp.stdout);attempts.append({'margin':margin,'shift':shift,'validation':v})
   if v['passed']:break
  else:continue
  cp2=subprocess.run([str(R/'code/interval_verify_signed'),str(witness),str(R/'results/leaves'/(stem+'_meanvalue256.csv')),'256','1',str(shift+1)],text=True,capture_output=True,timeout=180)
  if not cp2.stdout.strip():attempts.append({'margin':margin,'stage':'meanvalue','stderr':cp2.stderr});continue
  v2=json.loads(cp2.stdout)
  (R/'results/validation'/(stem+'_taylor160.json')).write_text(json.dumps(v,indent=2));(R/'results/validation'/(stem+'_meanvalue256.json')).write_text(json.dumps(v2,indent=2))
  if not v2['passed']:continue
  envgap=low_difference(v['certified_lower_mOhm'],q['recorded_envelope_upper_mOhm'])
  result.update({'status':'EXCLUDED_INTERVAL_VERIFIED','numerical_surplus':margin,'certified_lower_mOhm':v['certified_lower_mOhm'],'point_gap_lower_mOhm':v['gap_lower_mOhm'],'recorded_envelope_gap_lower_mOhm':envgap,'recorded_envelope_excluded':envgap>0,'pass_taylor160':True,'pass_meanvalue256':True,'tail_shift':shift,'taylor_leaves':v['accepted_leaves'],'meanvalue_leaves':v2['accepted_leaves'],'witness_file':str(witness.relative_to(R)),'candidate_file':str((R/'candidates'/(stem+'.json')).relative_to(R))})
  break
 (R/'results/validation'/(q['case_id']+f'_alpha{alpha:g}_attempts.json')).write_text(json.dumps(attempts,indent=2))
 result['elapsed_seconds']=time.time()-tt
 print(json.dumps(result),flush=True);return result
if __name__=='__main__':
 tasks=[(p,a) for p in sorted((R/'inputs').glob('*.json')) if not p.name.endswith('_EIS.json') for a in [1.,2.]]
 assert len(tasks)==24;out=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  fs=[pool.submit(one,p,a) for p,a in tasks]
  for f in concurrent.futures.as_completed(fs):
   out.append(f.result());pd.DataFrame(out).sort_values(['case_id','alpha']).to_csv(R/'results/EXTERNAL_CERTIFICATES_24.csv',index=False)
 print('ALL COMPLETED',len(out),flush=True)
