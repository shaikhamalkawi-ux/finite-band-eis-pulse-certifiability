#!/usr/bin/env python3
"""R41 portable nonvacuity check, seven original atomic EIS feasible witnesses.
Original models feasible at full 48 frequencies imply feasibility for all support subsets.
MPMath directed interval evaluation at 75 decimal digits; no optimizer or raw source.
"""
from pathlib import Path
import json, mpmath as mp
R=Path(__file__).resolve().parents[1]/'evidence';mp.iv.dps=75;iv=mp.iv

def pt(x):
 n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)

summary=[]
for soc in range(20,81,10):
 q=json.loads((R/'inputs'/f'SOC{soc}_charge_0.7C.json').read_text())
 m=json.loads((R/'primal'/f'SOC{soc}_primal.json').read_text())
 assert all(x>=0 for x in [m['R0'],m['L'],m['c']]+m['mass']) and all(t>0 for t in m['tau'])
 count=0; min_slack=1e100
 for f,zr,zi in zip(q['frequency_Hz'],q['Zreal_mOhm'],q['Zimag_mOhm']):
  w=2*iv.pi*pt(f);rr=pt(m['R0']);ii=w*pt(m['L'])-pt(m['c'])/w
  for tau,mass in zip(m['tau'],m['mass']):
   wt=w*pt(tau);rr+=pt(mass)/(1+wt*wt);ii-=pt(mass)*wt/(1+wt*wt)
  for value,recorded in [(rr,zr),(ii,zi)]:
   s=pt(q['eps_mOhm'])-abs(value-pt(recorded));assert s.a>0,(soc,count)
   min_slack=min(min_slack,float(s.a));count+=1
 assert count==96
 summary.append({'SOC':soc,'eis_coordinates':count,'positive_margin_verified':True,'minimum_margin_mOhm':min_slack})
print(json.dumps({'status':'PASS','SOC_models':7,'coordinates_verified':672,'rows':summary},indent=2))
(R/'R41_Seven_FullBand_Primal_Nonvacuity_Replay.json').write_text(json.dumps(summary,indent=2)+'\n')
