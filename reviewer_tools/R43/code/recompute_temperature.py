#!/usr/bin/env python3
"""Explicit finite-grid temperature sensitivity; never a continuum certificate."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import csr_matrix,vstack
ROOT=Path(__file__).resolve().parents[1]
def matrix_and_functional(q,n):
 tau=np.logspace(-10,10,n);w=2*np.pi*np.asarray(q['frequency_Hz']);x=w[:,None]*tau
 re=np.zeros((len(w),n+3));im=np.zeros_like(re);re[:,0]=1;im[:,1]=w;im[:,2]=-1/w;re[:,3:]=1/(1+x*x);im[:,3:]=-x/(1+x*x)
 g=np.zeros(n+3);g[0]=q['iT'];g[2]=sum(u*(b-a) for a,b,u in zip(q['a_s'],q['b_s'],q['u']))
 for a,b,u in zip(q['a_s'],q['b_s'],q['u']):g[3:]+=u*(np.exp(-(q['T_s']-b)/tau)-np.exp(-(q['T_s']-a)/tau))
 return np.vstack([re,im]),g,tau

def main():
 out=ROOT/'evidence';rows=[];certs={}
 scatter=pd.read_csv(out/'R42_Raw_Scatter_66_Conditions.csv').set_index('trace_id')
 for path in sorted((out/'fresh_inputs').glob('T*_SOC50_*.json')):
  q=json.loads(path.read_text());trace=path.stem;vals={}
  for n in [1201,2401]:
   M,g,tau=matrix_and_functional(q,n);b=np.r_[q['Zreal_mOhm'],q['Zimag_mOhm']];eps=q['eps_mOhm'];column_scale=np.max(abs(M),axis=0);Ms=M/column_scale;gs=g/column_scale;A=vstack([csr_matrix(Ms),csr_matrix(-Ms)],format='csr');ub=np.r_[b+eps,-b+eps]
   endpoints=[]
   for side,sg in [('lower',1),('upper',-1)]:
    r=linprog(sg*gs,A_ub=A,b_ub=ub,bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-8,'dual_feasibility_tolerance':1e-8})
    if not r.success:raise RuntimeError((trace,n,side,r.message))
    xx=r.x/column_scale
    error=float(np.max(abs(M@xx-b))-eps)
    if error>2e-7:raise RuntimeError((trace,n,side,error))
    value=float(g@xx);endpoints.append(value)
    nz=np.flatnonzero(xx>0)
    certs[f'{trace}_n{n}_{side}']={'trace_id':trace,'grid_n':n,'side':side,'objective_mOhm':value,'maximum_primal_excess_mOhm':error,'nonzero_indices':nz.tolist(),'nonzero_values':xx[nz].tolist(),'maximum_dual_infeasibility':float(max(0.,-(sg*gs-A.T@r.ineqlin.marginals).min())),'solver_message':r.message}
   lo,hi=endpoints;rv=q['Rmeas_mOhm'];gap=max(0.,lo-rv,rv-hi)
   vals[n]=(lo,hi,gap)
  lo,hi,gap=vals[1201];h=float(scatter.loc[trace,'halfwidth_3sigma_endpoint_mOhm'])
  rows.append({'trace_id':trace,'temperature_C':q['temperature_C'],'SOC':50,'direction':q['direction'],'C_rate':q['C_rate'],'eps_mOhm':q['eps_mOhm'],'Rmeas_mOhm':q['Rmeas_mOhm'],'grid1201_lower_mOhm':lo,'grid1201_upper_mOhm':hi,'grid1201_gap_mOhm':gap,'grid2401_lower_mOhm':vals[2401][0],'grid2401_upper_mOhm':vals[2401][1],'grid2401_gap_mOhm':vals[2401][2],'maximum_endpoint_change_mOhm':max(abs(lo-vals[2401][0]),abs(hi-vals[2401][1])),'local_halfwidth_3sigma_mOhm':h,'grid_gap_after_3sigma_mOhm':max(0.,gap-h)})
  print(trace,'gap',gap,'grid_change',rows[-1]['maximum_endpoint_change_mOhm'],flush=True)
 df=pd.DataFrame(rows).sort_values(['temperature_C','C_rate','direction']);df.to_csv(out/'R42_Temperature_Fresh_Grid_30.csv',index=False,float_format='%.17g')
 (out/'R42_Temperature_Grid_Primal_120.json').write_text(json.dumps(certs,indent=2)+'\n')
 summary={'status':'PASS','conditions':30,'linear_programs':120,'relaxation_grids':[1201,2401],'tau_domain_s':[1e-10,1e10],'maximum_endpoint_change_mOhm':float(df.maximum_endpoint_change_mOhm.max()),'maximum_primal_excess_mOhm':max(x['maximum_primal_excess_mOhm'] for x in certs.values()),'scope':'Exploratory finite-grid sensitivity only. No warm-case continuum exclusion or optimality claim.'}
 (out/'R42_Temperature_Grid_Summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
 print(df.groupby(['temperature_C','direction']).grid1201_gap_mOhm.agg(['median','min','max']).to_string())
if __name__=='__main__':main()
