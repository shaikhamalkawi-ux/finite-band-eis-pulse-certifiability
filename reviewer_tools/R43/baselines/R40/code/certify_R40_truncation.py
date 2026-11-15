#!/usr/bin/env python3
"""Attempt fresh strict full-domain certificates for four R24-selected 25.31mHz cases using R39 inputs.
Not confirmation until every independent method and tail passes.
"""
import os, pathlib, subprocess, json, csv, sys, math, importlib.util
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import pandas as pd
R=pathlib.Path(__file__).resolve().parent.parent/'results'
SRC=R.parent/'baseline'/'R38'
from compute_R40_lockstep import inputs
spec=importlib.util.spec_from_file_location('candidate_generation',SRC/'generate_candidates.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
CASES=['SOC70_charge_1.4C','SOC70_discharge_1.4C','SOC70_charge_1.0C','SOC60_discharge_1.4C']
S=pd.read_csv(R/'R40_FiniteBand_Grid_Screen_All42.csv')
C=R/'trunc25';C.mkdir(exist_ok=True)

exe=C/'interval_verify'
if not exe.exists():
    p=subprocess.run(['g++','-O2','-std=c++17',str(SRC/'interval_verify.cpp'),'-Wl,-l:libmpfr.so.6','-lgmp','-o',str(exe)],capture_output=True,text=True)
    if p.returncode:raise RuntimeError(p.stderr)
print('verifier compiled',flush=True)
rows=[]
for k in CASES:
    q=inputs()[k]
    s=S[(S.case_id==k)&(S.requested_cutoff_mHz==20.0)].iloc[0]
    f=q['frequency_Hz'];mask=[v>=.02*.999 for v in f]
    qt={**q,'frequency_Hz':[v for v,b in zip(q['frequency_Hz'],mask) if b],
          'Zreal_mOhm':[v for v,b in zip(q['Zreal_mOhm'],mask) if b],
          'Zimag_mOhm':[v for v,b in zip(q['Zimag_mOhm'],mask) if b],
          'eps_mOhm':float(s.analysis_eps_mOhm)}
    assert len(qt['frequency_Hz'])==44
    (C/(k+'_q.json')).write_text(json.dumps(qt,indent=2)+'\n')
    solved=None
    for margin in [.02,.01,.005,.002,.001]:
        r=mod.solve(qt,1.,margin)
        print(k,'candidate',margin, r.get('solver_success'),r.get('floating_gap_mOhm'),r.get('finite_slacks_float'),flush=True)
        if r.get('solver_success') and r['floating_gap_mOhm']>0:
            solved=r;break
    if solved is None:
        rows.append({'case_id':k,'status':'CANDIDATE_SOLVER_FAIL'});continue
    (C/(k+'_candidate.json')).write_text(json.dumps(solved,indent=2)+'\n')
    y=solved['y'];n=len(qt['frequency_Hz']);assert len(y)==2*n
    lines=[' '.join(map(str,[n,len(qt['u']),qt['T_s'],qt['iT'],qt['eps_mOhm'],qt['Rmeas_mOhm'],1.]))]
    for j in range(n): lines.append(' '.join(repr(float(v)) for v in [qt['frequency_Hz'][j],qt['Zreal_mOhm'][j],qt['Zimag_mOhm'][j],y[j],y[j+n]]))
    for a,b,v in zip(qt['a_s'],qt['b_s'],qt['u']):lines.append(' '.join(repr(float(x)) for x in [a,b,v]))
    wf=C/(k+'_witness.txt');wf.write_text('\n'.join(lines)+'\n')
    v=[]
    for bits,mode,shift,tag in [(160,2,0,'taylor160'),(256,1,1,'meanvalue256')]:
        leaf=C/(k+'_'+tag+'_leaves.csv')
        p=subprocess.run([str(exe),str(wf),str(leaf),str(bits),str(mode),str(shift)],capture_output=True,text=True,timeout=120)
        try:j=json.loads(p.stdout)
        except Exception as e:
            j={'passed':False,'stderr':p.stderr,'stdout':p.stdout[:1000],'exception':str(e)}
        (C/(k+'_'+tag+'_verification.json')).write_text(json.dumps(j,indent=2)+'\n')
        v.append(j)
        print(k,tag,'passed',j.get('passed'),'gap',j.get('gap_lower_mOhm'),'leaves',j.get('accepted_leaves'),'stderr',p.stderr[:100],flush=True)
    if all(x.get('passed',False) for x in v):
        lb=min(z['certified_lower_mOhm'] for z in v);gap=lb-qt['Rmeas_mOhm'];print(k,'CERTIFIED g=',gap,flush=True)
        rows.append({'case_id':k,'status':'PASS','measured_mOhm':qt['Rmeas_mOhm'],'certified_lower_mOhm':lb,'gap_mOhm':gap,'eps_mOhm':qt['eps_mOhm'],'taylor_leaves':v[0]['accepted_leaves'],'meanvalue_leaves':v[1]['accepted_leaves']})
    else:rows.append({'case_id':k,'status':'HOLD'})
D=pd.DataFrame(rows);D.to_csv(R/'R40_Truncated_FreshMPFR_4Case.csv',index=False)
print(D.to_string(index=False),flush=True)
