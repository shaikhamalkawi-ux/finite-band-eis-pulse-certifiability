#!/usr/bin/env python3
"""Exact-rational sensitivity of R39 measured-normalizer dual portfolios. Not instrument calibration."""
from pathlib import Path
from fractions import Fraction as F
import json,math,itertools,statistics,csv
import pandas as pd
R=Path(__file__).resolve().parent/'r39_candidate';C=pd.read_csv(R/'R39_FRESH_CERTIFICATES_84.csv');P=pd.read_csv(R/'RAW_CURRENT_CORRECTION_42.csv').set_index('case_id')
def f(x):return F.from_float(float(x))
def down(q):return math.nextafter(float(q),-math.inf) if F.from_float(float(q))>q else float(q)
def write(name,arr):
 with (R/name).open('w',newline='') as g:
  w=csv.DictWriter(g,fieldnames=list(arr[0]));w.writeheader();w.writerows(arr)
portfolios={}
for row in C.itertuples():
 q=json.loads((R/'inputs'/(row.case_id+'.json')).read_text());stem=f'{row.case_id}_alpha{row.alpha:g}_margin0.005';c=json.loads((R/'candidates'/(stem+'.json')).read_text());n=len(q['frequency_Hz']);v=list(map(f,c['y']));vr=v[:n];vi=v[n:];b=list(map(f,q['Zreal_mOhm']+q['Zimag_mOhm']));pr=P.loc[row.case_id]
 dot=sum((z*x for z,x in zip(v,b)),F(0));rot=sum((-x*y for x,y in zip(vr,b[n:])),F(0))+sum((x*y for x,y in zip(vi,b[:n])),F(0))
 o={'case_id':row.case_id,'from_alpha':row.alpha,'dot':dot,'l1':sum(map(abs,v),F(0)),'real':abs(sum(vr)),'imag':abs(sum(vi)),'gain':abs(dot),'quadrature':abs(rot),'eps':f(q['eps_mOhm']),'point':f(q['Rmeas_mOhm']),'recorded':f(pr.raw_voltage_sample_envelope_high_mOhm),'current':f(q['normalization_current_A'])}
 portfolios.setdefault(row.case_id,[]).append(o)
assert len(portfolios)==42 and all(len(x)==2 for x in portfolios.values())
grids={'pulse_bias_mOhm':['0','.005','.01','.025','.05','.1','.2','.5','1','2','5'],
 'EIS_box_mOhm':['0','.0001','.0005','.001','.002','.005','.01','.02','.05','.1','.2','.5','1'],
 'EIS_common_real_mOhm':['0','.005','.01','.025','.05','.1','.2','.5','1','2','5'],
 'EIS_common_imag_mOhm':['0','.005','.01','.025','.05','.1','.2','.5','1','2','5'],
 'EIS_gain_percent':['0','.005','.01','.025','.05','.1','.2','.5','1','2','5'],
 'pulse_each_gain_percent':['0','.005','.01','.025','.05','.1','.2','.5','1','2','5'],
 'pulse_differential_voltage_mV':['0','.01','.025','.05','.1','.2','.5','1','2','5']}
def margin(c,alpha,env,kind,val,second=F(0)):
 u=c[env];base=-c['dot']-alpha*c['eps']*c['l1']-u
 if kind=='pulse_bias_mOhm':cost=val
 elif kind=='EIS_box_mOhm':cost=val*c['l1']
 elif kind=='EIS_common_real_mOhm':cost=val*c['real']
 elif kind=='EIS_common_imag_mOhm':cost=val*c['imag']
 elif kind=='EIS_gain_percent':cost=val/100*c['gain']
 elif kind=='pulse_each_gain_percent':
  x=val/100;cost=u*((1+x)/(1-x)-1)
 elif kind=='pulse_differential_voltage_mV':cost=val/c['current']
 elif kind=='joint_EISbox_pulsebias':cost=val*c['l1']+second
 else:raise ValueError(kind)
 return base-cost
rows=[];summary=[];budgets=[]
for alpha in [F(1),F(2)]:
 for env in ['point','recorded']:
  for cid,bank in sorted(portfolios.items()):
   d={'case_id':cid,'alpha':float(alpha),'expansion':env}
   for kind in grids:
    vals=[]
    for c in bank:
     gap=margin(c,alpha,env,kind,F(0));
     if kind=='pulse_bias_mOhm':v=gap
     elif kind=='EIS_box_mOhm':v=gap/c['l1']
     elif kind=='EIS_common_real_mOhm':v=gap/c['real'] if c['real'] else F(10**30)
     elif kind=='EIS_common_imag_mOhm':v=gap/c['imag'] if c['imag'] else F(10**30)
     elif kind=='EIS_gain_percent':v=100*gap/c['gain']
     elif kind=='pulse_each_gain_percent':v=100*gap/(2*c[env]+gap)
     elif kind=='pulse_differential_voltage_mV':v=gap*c['current']
     vals.append(v)
    d[kind]=down(max(vals))
   budgets.append(d)
  for kind,grid in grids.items():
   counts=[]
   for s in grid:
    v=F(s);ms=[]
    for cid,bank in sorted(portfolios.items()):
     mm=max(margin(c,alpha,env,kind,v) for c in bank);ms.append(mm)
     rows.append({'case_id':cid,'alpha':float(alpha),'expansion':env,'scenario':kind,'magnitude':s,'margin_lower_mOhm':down(mm),'certified':mm>0})
    count=sum(m>0 for m in ms);counts.append(count)
    summary.append({'alpha':float(alpha),'expansion':env,'scenario':kind,'magnitude':s,'certified_cases':count,'min_margin_lower_mOhm':down(min(ms)),'median_margin_lower_mOhm':down(statistics.median(ms))})
   assert all(a>=b for a,b in zip(counts,counts[1:]))
write('R39_SENSITIVITY_CASES.csv',rows);write('R39_SENSITIVITY_SUMMARY.csv',summary);write('R39_SENSITIVITY_BUDGETS.csv',budgets)
summ=pd.DataFrame(summary);budget=pd.DataFrame(budgets)
for alpha in [1.,2.]:
 t=budget[(budget.alpha==alpha)&(budget.expansion=='recorded')]
 print('\nalpha',alpha,'min budgets recorded envelope')
 print(t[list(grids)].min().to_string())
 for scenario,val in [('pulse_bias_mOhm','0'),('pulse_bias_mOhm','.01'),('pulse_bias_mOhm','.1'),('pulse_bias_mOhm','.5'),('EIS_box_mOhm','.1'),('pulse_each_gain_percent','.1')]:
  t1=summ[(summ.alpha==alpha)&(summ.expansion=='recorded')&(summ.scenario==scenario)&(summ.magnitude==val)]
  print(scenario,val,'certified',t1.certified_cases.iloc[0])
result={'status':'PASS','cases':42,'witnesses':84,'scenarios':len(rows),'summaries':len(summary),'budget_entries':len(budgets),'recorded_range':'Fixed 5–15 s median current reference with recorded local voltage-sample envelope. No metrological calibration.'}
(R/'R39_SENSITIVITY_QA.json').write_text(json.dumps(result,indent=2)+'\n')
