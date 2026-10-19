#!/usr/bin/env python3
"""Construct EIS-only feasible atomic models; fix tolerance before response optimization."""
from pathlib import Path
import json,numpy as np,pandas as pd,mpmath as mp
from scipy.optimize import linprog
from scipy import sparse
from generate_candidates import model
R=Path(__file__).resolve().parents[1];mp.iv.dps=75;iv=mp.iv

def pt(x):n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)
def save(p,d):p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
rows=[]
for cell in ['F0001','F0002','F0003']:
 for stage in ['BOL','EOL']:
  q=json.loads((R/'inputs'/f'{cell}_{stage}_charge.json').read_text())
  tau=np.logspace(-10,10,1201);M,_=model(q,tau);scale=np.maximum(np.max(abs(M),axis=0),1e-200);MS=M/scale
  obs=np.r_[q['Zreal_mOhm'],q['Zimag_mOhm']];n=MS.shape[1]
  A=sparse.vstack([sparse.hstack([sparse.csr_matrix(MS),-np.ones((len(obs),1))]),sparse.hstack([sparse.csr_matrix(-MS),-np.ones((len(obs),1))])],format='csr')
  sol=linprog(np.r_[np.zeros(n),1.],A_ub=A,b_ub=np.r_[obs,-obs],bounds=(0,None),method='highs-ds',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
  if not sol.success: raise RuntimeError((cell,stage,sol.message))
  x=sol.x[:-1]/scale; eps=float(sol.fun)+.001
  assert x.min()>=0 and eps>0
  inds=np.flatnonzero(x[3:]>0)
  pr={'cell':cell,'stage':stage,'R0':float(x[0]),'L':float(x[1]),'c':float(x[2]),'tau':tau[inds].tolist(),'mass':x[3:][inds].tolist(),'grid_minimum_mOhm':float(sol.fun),'numerical_buffer_mOhm':.001,'eps_mOhm':eps}
  checks=[];smins=[]
  for freq,zr,zi in zip(q['frequency_Hz'],q['Zreal_mOhm'],q['Zimag_mOhm']):
   w=2*iv.pi*pt(freq);re=pt(x[0]);im=w*pt(x[1])-pt(x[2])/w
   for t,mass in zip(pr['tau'],pr['mass']):
    u=w*pt(t);re+=pt(mass)/(1+u*u);im-=pt(mass)*u/(1+u*u)
   for k,val,record in [('real',re,zr),('imag',im,zi)]:
    slack=pt(eps)-abs(val-pt(record));ok=bool(slack.a>0);assert ok,(cell,stage,freq,slack)
    smins.append(float(slack.a));checks.append({'frequency_Hz':freq,'component':k,'slack_interval':str(slack),'pass':ok})
  save(R/'results/primal'/f'{cell}_{stage}_model.json',pr);save(R/'results/primal'/f'{cell}_{stage}_interval.json',checks)
  rows.append({'cell':cell,'stage':stage,'grid_minimum_mOhm':float(sol.fun),'eps_mOhm':eps,'support_atoms':len(inds),'constraints_verified':len(checks),'minimum_slack_mOhm':min(smins),'nonempty_verified':True})
  for direction in ['charge','discharge']:
   p=R/'inputs'/f'{cell}_{stage}_{direction}.json';qq=json.loads(p.read_text());qq['eps_mOhm']=eps;qq['grid_eis_minimum_mOhm']=float(sol.fun);save(p,qq)
  print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(R/'results/EIS_NONVACUITY_6.csv',index=False)
