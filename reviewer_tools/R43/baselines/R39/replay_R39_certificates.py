#!/usr/bin/env python3
"""Recheck transported dual witnesses with fresh MPFR Taylor160 and mean-value256 enclosures."""
from pathlib import Path
from fractions import Fraction
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import subprocess,csv,json,math
import pandas as pd
D=Path(__file__).resolve().parent;R=D/'r39_candidate';S=D/'replay/Batteries52_Reviewer_Replay/R36_reproducibility';C=pd.read_csv(S/'FINAL_ALL42_CERTIFICATES.csv');P=pd.read_csv(R/'RAW_CURRENT_CORRECTION_42.csv').set_index('case_id')
assert len(C)==84 and len(P)==42
for t in ['leaves','validation']:(R/t).mkdir(exist_ok=True)
def run(i):
 old=C.iloc[i];q=json.loads((R/'inputs'/(old.case_id+'.json')).read_text());pr=P.loc[old.case_id];stem=Path(old.witness_file).stem;w=R/'witnesses'/(stem+'.txt');rows=[line.split() for line in w.read_text().splitlines()];m=len(q['frequency_Hz']);n=len(q['u']);c=json.loads((R/'candidates'/(stem+'.json')).read_text())
 assert len(rows)==1+m+n and [float(x) for x in rows[0]]==[m,n,q['T_s'],q['iT'],q['eps_mOhm'],q['Rmeas_mOhm'],float(old.alpha)]
 for k in range(m):assert [float(x) for x in rows[1+k]]==[q['frequency_Hz'][k],q['Zreal_mOhm'][k],q['Zimag_mOhm'][k],c['y'][k],c['y'][k+m]]
 for k in range(n):assert [float(x) for x in rows[1+m+k]]==[q['a_s'][k],q['b_s'][k],q['u'][k]]
 v=[]
 for bits,mode,shift,tag in [(160,2,0,'taylor160'),(256,1,1,'meanvalue256')]:
  leaf=R/'leaves'/(stem+'_'+tag+'.csv');p=subprocess.run([str(R/'interval_verify'),str(w),str(leaf),str(bits),str(mode),str(shift)],capture_output=True,text=True,timeout=150)
  assert p.stdout.strip(),(old.case_id,old.alpha,p.stderr)
  z=json.loads(p.stdout);assert z['passed'] and p.returncode==0,(old.case_id,old.alpha,bits,z)
  t=sorted(csv.DictReader(leaf.open()),key=lambda a:float(a['x_lo']));x=Fraction(float(z['xl']));assert len(t)==z['accepted_leaves']
  for row in t:
   a,b=Fraction(float(row['x_lo'])),Fraction(float(row['x_hi']));assert a==x and b>a and float(row['residual_lower'])>0;x=b
  assert x==Fraction(float(z['xh']))
  (R/'validation'/(stem+'_'+tag+'.json')).write_text(json.dumps(z,indent=2)+'\n');v.append(z)
 low=min(z['certified_lower_mOhm'] for z in v)
 g=Fraction(float(low))-Fraction(float(q['Rmeas_mOhm']))
 e=Fraction(float(low))-Fraction(float(pr['raw_voltage_sample_envelope_high_mOhm']))
 assert g>0 and e>0
 return {'case_id':old.case_id,'alpha':old.alpha,'SOC':q['SOC'],'direction':q['direction'],'C_rate':q['C_rate'],'measured_mOhm':q['Rmeas_mOhm'],'certified_lower_mOhm':low,'point_gap_lower_mOhm':math.nextafter(float(g),-math.inf),'recorded_voltage_envelope_gap_lower_mOhm':math.nextafter(float(e),-math.inf),'taylor160_leaves':v[0]['accepted_leaves'],'meanvalue256_leaves':v[1]['accepted_leaves'],'minimum_compact_residual':min(z['compact_residual_lower'] for z in v),'minimum_finite_R0_residual':min(z['finite_R0_lower'] for z in v),'minimum_finite_L_residual':min(z['finite_L_lower'] for z in v),'minimum_finite_c_residual':min(z['finite_c_lower'] for z in v),'minimum_fast_tail_residual':min(z['fast_tail_lower'] for z in v),'minimum_slow_tail_residual':min(z['slow_tail_H_lower'] for z in v),'witness_sha256':sha256(w.read_bytes()).hexdigest(),'scale':pr['rescale']}
if __name__=='__main__':
 with ThreadPoolExecutor(max_workers=8) as pool:
  d=[]
  for i,row in enumerate(pool.map(run,range(len(C))),1):
   d.append(row)
   if i%12==0:print('MPFR VERIFIED',i,'/84',flush=True)
 d=pd.DataFrame(d).sort_values(['case_id','alpha']);d.to_csv(R/'R39_FRESH_CERTIFICATES_84.csv',index=False,float_format='%.17g')
 result={'status':'PASS','cases':42,'dual_witnesses':84,'fresh_full_domain_MPFR_runs':168,'source':'RWTH exact published raw archive','normalization':'median measured current 5-15s in each pulse','solver':'No new optimization; original dual witnesses scaled and numerically reverified','recorded_range':'fixed-normalizer voltage-sample envelope; no calibration or confidence-level interpretation','all_point_and_sample_envelope_gaps_positive':True}
 for a in [1.,2.]:
  t=d[d.alpha==a];result[f'factor_{a:g}']={'min_point_gap_mOhm':float(t.point_gap_lower_mOhm.min()),'median_point_gap_mOhm':float(t.point_gap_lower_mOhm.median()),'min_sample_envelope_gap_mOhm':float(t.recorded_voltage_envelope_gap_lower_mOhm.min()),'min_case':t.loc[t.point_gap_lower_mOhm.idxmin(),'case_id'],'min_sample_case':t.loc[t.recorded_voltage_envelope_gap_lower_mOhm.idxmin(),'case_id']}
 (R/'R39_FRESH_VERIFY_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL',json.dumps(result,indent=2),flush=True)
