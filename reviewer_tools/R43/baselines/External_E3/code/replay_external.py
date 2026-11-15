#!/usr/bin/env python3
"""Replay frozen external certificates without optimization or raw battery downloads.
A valid run checks SHA-256 before compiling and writes only to a new output directory.
It verifies the declared saved computational instances, not instrument calibration.
"""
from __future__ import annotations
import argparse,concurrent.futures,csv,hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]

def require(ok:bool,message:str)->None:
    if not ok: raise RuntimeError(message)
def read(path:Path): return json.loads(path.read_text())
def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()
def pt(x): n,d=float(x).as_integer_ratio();return mp.iv.mpf(n)/mp.iv.mpf(d)
def nonvac(q,p):
    require(all(float(p[k])>=0 for k in ['R0','L','c']) and all(v>=0 for v in p['mass']) and all(v>0 for v in p['tau']),'Nonpositive primal coordinates')
    slacks=[]
    for f,zr,zi in zip(q['frequency_Hz'],q['Zreal_mOhm'],q['Zimag_mOhm']):
        w=2*mp.iv.pi*pt(f);re=pt(p['R0']);im=w*pt(p['L'])-pt(p['c'])/w
        for t,m in zip(p['tau'],p['mass']):
            x=w*pt(t);re+=pt(m)/(1+x*x);im-=pt(m)*x/(1+x*x)
        for v,b in [(re,zr),(im,zi)]:
            s=pt(q['eps_mOhm'])-abs(v-pt(b));require(s.a>0,'EIS model is not interval feasible');slacks.append(float(s.a))
    return {'cell':q['cell'],'stage':q['stage'],'constraints':len(slacks),'minimum_slack_mOhm':min(slacks),'passed':True}
def main()->None:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workers',type=int,default=4);p.add_argument('--output',type=Path);a=p.parse_args()
    require(1<=a.workers<=16,'workers must be 1..16')
    require(__debug__ and os.getenv('PYTHONOPTIMIZE','') in ('','0'),'Run without python -O or PYTHONOPTIMIZE')
    manifest=read(ROOT/'MANIFEST.json');require(len(manifest)>100,'Missing complete manifest')
    for name,digest in manifest.items():
        path=ROOT/name;require(path.is_file() and sha(path)==digest,'Missing or changed artifact: '+name)
    print('MANIFEST PASS',len(manifest),flush=True)
    out=a.output or Path(tempfile.mkdtemp(prefix='project52_external_replay_'))
    if a.output:out.mkdir(parents=True,exist_ok=False)
    exe=out/'interval_verify_signed'
    cmd=['g++','-O2','-std=c++17',str(ROOT/'code/interval_verify_signed.cpp'),'-lmpfr','-lgmp','-o',str(exe)]
    c=subprocess.run(cmd,capture_output=True,text=True)
    if c.returncode:
        cmd[4]='-Wl,-l:libmpfr.so.6';c=subprocess.run(cmd,capture_output=True,text=True)
    require(c.returncode==0,'Compiler/MPFR dependencies missing: '+c.stderr)
    rows=list(csv.DictReader((ROOT/'results/EXTERNAL_CERTIFICATES_24.csv').open()))
    require(len(rows)==24 and len({(r['case_id'],r['alpha']) for r in rows})==24,'Case coverage mismatch')
    tasks=[]
    for r in rows:
        q=read(ROOT/'inputs'/(r['case_id']+'.json'));c=read(ROOT/r['candidate_file']);w=(ROOT/r['witness_file']).read_text().splitlines();n=len(q['frequency_Hz']);npulse=len(q['u']);h=list(map(float,w[0].split()))
        require(h==[n,npulse,q['T_s'],q['iT'],q['eps_mOhm'],q['Rmeas_mOhm'],c['alpha']],'Witness input header mismatch')
        require(len(w)==1+n+npulse,'Unexpected witness row count')
        for k in range(n):require(list(map(float,w[1+k].split()))==[q['frequency_Hz'][k],q['Zreal_mOhm'][k],q['Zimag_mOhm'][k],c['y'][k],c['y'][k+n]],'Witness EIS/dual mismatch')
        for k in range(npulse):require(list(map(float,w[1+n+k].split()))==[q['a_s'][k],q['b_s'][k],q['u'][k]],'Witness waveform mismatch')
        for bits,mode,suffix,add in [(160,2,'taylor160',0),(256,1,'meanvalue256',1)]:tasks.append((r,bits,mode,suffix,int(r['tail_shift'])+add))
    def run(t):
        r,bits,mode,suffix,shift=t;stem=Path(r['witness_file']).stem;csvpath=out/(stem+'_'+suffix+'.csv')
        cp=subprocess.run([str(exe),str(ROOT/r['witness_file']),str(csvpath),str(bits),str(mode),str(shift)],text=True,capture_output=True,timeout=300)
        require(cp.returncode==0 and cp.stdout.strip(),'Verification failed: '+stem+' '+cp.stderr)
        v=json.loads(cp.stdout);expected=read(ROOT/'results/validation'/(stem+'_'+suffix+'.json'))
        require(v==expected,'Saved verification output differs: '+stem+' '+suffix)
        leaves=sorted(csv.DictReader(csvpath.open()),key=lambda k:float(k['x_lo']))
        require(float(leaves[0]['x_lo'])==v['xl'] and float(leaves[-1]['x_hi'])==v['xh'],'Incomplete compact endpoints')
        require(all(float(k['residual_lower'])>0 and float(k['x_hi'])>float(k['x_lo']) for k in leaves),'Nonpositive compact leaf')
        require(all(leaves[k]['x_hi']==leaves[k+1]['x_lo'] for k in range(len(leaves)-1)),'Gap or overlap in compact cover')
        (out/(stem+'_'+suffix+'.json')).write_text(json.dumps(v,indent=2));print('PASS',stem,suffix,flush=True);return {'case_id':r['case_id'],'alpha':float(r['alpha']),'method':suffix,'passed':True}
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:ver=list(pool.map(run,tasks))
    mp.iv.dps=75
    import audit_external
    arithmetic=[];primals=[]
    for r in rows:
        q=read(ROOT/'inputs'/(r['case_id']+'.json'));c=read(ROOT/r['candidate_file']);stem=Path(r['witness_file']).stem
        v=read(out/(stem+'_taylor160.json'));v2=read(out/(stem+'_meanvalue256.json'));j=audit_external.independent(q,c,v,v2)
        require(all(j[k] for k in ['objective_inside_MPFR','finite_constraints_positive','right_endpoint_positive','both_tails_positive','recorded_envelope_excluded']),'Independent arithmetic failed')
        arithmetic.append(j)
        if q['direction']=='charge' and float(r['alpha'])==1:primals.append(nonvac(q,read(ROOT/'results/primal'/(q['cell']+'_'+q['stage']+'_model.json'))))
    def neg(r):
        cp=subprocess.run([str(exe),str(ROOT/'qa/negative_controls'/(r['case_id']+'.txt')),str(out/(r['case_id']+'_negative.csv')),'160','2',str(int(r['tail_shift']))],capture_output=True,text=True,timeout=180)
        v=json.loads(cp.stdout);require(cp.returncode==1 and not v['passed'],'Invalid control was not rejected');return {'case_id':r['case_id'],'rejected':True}
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:bad=list(pool.map(neg,[r for r in rows if float(r['alpha'])==1]))
    summary={'status':'PASS','manifest_files':len(manifest),'input_witness_pairs':24,'full_domain_runs':len(ver),'independent_arithmetic_witnesses':len(arithmetic),'feasible_models':len(primals),'primal_constraints':sum(x['constraints'] for x in primals),'invalid_controls_rejected':len(bad),'all_outputs_identical':True,'scope':'Saved external computational inputs; no optimization, new experiment, calibration or independent-human review.'}
    (out/'REPLAY_SUMMARY.json').write_text(json.dumps(summary,indent=2));(out/'ARITHMETIC_REPLAY.json').write_text(json.dumps(arithmetic,indent=2));(out/'PRIMAL_REPLAY.json').write_text(json.dumps(primals,indent=2));print(json.dumps(summary,indent=2));print('RESULTS',out)
if __name__=='__main__':
    try:main()
    except Exception as e:print('FAIL',e,file=sys.stderr);sys.exit(1)
