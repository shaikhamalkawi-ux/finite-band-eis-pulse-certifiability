#!/usr/bin/env python3
"""Read-only verification of the corrected public numerical mirror.
No optimizer or MPFR subprocess is run. Use replay_saved.py for a fresh proof run.
"""
from pathlib import Path
import hashlib,json,subprocess,sys
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];checks=[]
def check(ok,name):
 if not ok:raise AssertionError(name)
 checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
manifest=load(R/'PUBLIC_MANIFEST_SHA256.json')
for f,h in manifest['files'].items():
 check((R/f).is_file() and sha(R/f)==h,'public hash '+f)
expected=set(manifest['files'])|{'PUBLIC_MANIFEST_SHA256.json'}
actual={str(p.relative_to(R)) for p in R.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
check(actual==expected,'complete public manifest coverage')
e=R/'evidence'
raw=pd.read_csv(e/'R42_Raw_Scatter_66_Conditions.csv');rec=pd.read_csv(e/'R42_R39_Raw_Input_Reconciliation_42.csv');g=pd.read_csv(e/'R42_Raw3sigma_Certificate_Gaps_84.csv')
check(len(raw)==66 and raw.trace_id.nunique()==66,'66 unique raw source conditions')
check(len(rec)==42 and rec.max_abs_difference.max()==0,'42 exact raw-to-frozen input reconstructions')
check(raw[raw.temperature_C==5].backward_timestamp_steps.sum()==0,'primary timestamps no reversal')
check((raw[raw.temperature_C!=5].backward_timestamp_steps>0).sum()==9,'nine warmer endpoint timestamp reversals disclosed')
check(raw.canonical_minus_legacy_R_mOhm.abs().max()<.067064,'timestamp sensitivity bounded in manuscript')
for x in raw.itertuples():
 w=load(e/'source_windows'/(x.trace_id+'.json'));q=load(e/'fresh_inputs'/(x.trace_id+'.json'))
 check(len(w['baseline_V'])==10 and all(v==0 for v in w['baseline_I_A']),x.trace_id+' baseline')
 check(abs(q['Rmeas_mOhm']-x.Rmeas_mOhm)<1e-11,x.trace_id+' response')
 check(all(b>a for a,b in zip(q['a_s'],q['b_s'])) and all(a==b for a,b in zip(q['a_s'][1:],q['b_s'][:-1])),x.trace_id+' contiguous waveform intervals')
 check(q['a_s'][0]==0 and q['b_s'][-1]==20,x.trace_id+' exact horizon')
 h=3*np.hypot(1000*np.hypot(x.sigma_baseline_V,x.sigma_endpoint_V)/x.I_p_A,x.Rmeas_mOhm*x.sigma_endpoint_I_A/x.I_p_A)
 check(abs(h-x.halfwidth_3sigma_endpoint_mOhm)<1e-10,x.trace_id+' local propagation')
check(len(g)==84 and len(g[['case_id','alpha']].drop_duplicates())==84,'84 primary case-factor combinations')
for kind in ['endpoint','plateau']:
 check((g[f'gap_after_3sigma_{kind}_mOhm']>0).all(),'84 positive '+kind+' expanded gaps')
 check(np.max(abs(g[f'gap_after_3sigma_{kind}_mOhm']-(g.point_gap_lower_mOhm-g[f'halfwidth_3sigma_{kind}_mOhm'])))<1e-12,'gap arithmetic '+kind)
check(abs(g[g.alpha==2].gap_after_3sigma_endpoint_mOhm.min()-.009321772570726054)<1e-11,'critical doubled local gap')
t=pd.read_csv(e/'R42_Temperature_Fresh_Grid_30.csv');check(len(t)==30,'30 temperature rows')
check(t.maximum_endpoint_change_mOhm.max()<.000290,'grid sensitivity maximum')
check(len(load(e/'R42_Temperature_Grid_Primal_120.json'))==120,'120 stored R42 LP outputs')
check(load(R/'qa/R42_Independent_Arithmetic_Summary.json')['status']=='PASS','independent 75-decimal arithmetic PASS')
tr=pd.read_csv(e/'R42_Truncated25mHz_Local3sigma_4.csv');check(len(tr)==4 and (tr.gap_after_endpoint_3sigma_mOhm>0).all() and (tr.gap_after_plateau_3sigma_mOhm>0).all(),'four truncated certificates survive both local windows')
ext=e/'External_Fresh_Replay';summary=load(ext/'REPLAY_SUMMARY.json')
check(summary['status']=='PASS' and summary['full_domain_runs']==48 and summary['all_outputs_identical'],'stored R42 external replay summary')
check(summary['primal_constraints']==852 and summary['invalid_controls_rejected']==12,'external primal and negative controls')
proofs=[p for p in ext.glob('*.json') if '_taylor160' in p.name or '_meanvalue256' in p.name]
check(len(proofs)==48,'48 stored R42 external proof outputs')
for p in proofs:
 check(load(p)==load(R/'baselines/External_E3/results/validation'/p.name),'stored R42 external exact result '+p.name)

for script in ['check_sparse_saved.py','worked_engineering_example.py']:
 result=subprocess.run([sys.executable,str(R/'code'/script)],capture_output=True,text=True,check=True)
 check('"status": "PASS"' in result.stdout,'independent check '+script)
print(json.dumps({'status':'PASS','checks_passed':len(checks),'public_files_hashed':len(manifest['files']),
 'numerical_baseline':'R42 unchanged; R43 worked example only',
 'scope':'Read-only hashes, raw-window arithmetic, stored proof records and sparse-support identities. Not fresh MPFR or calibration.'},indent=2))
