#!/usr/bin/env python3
"""Fresh MPFR replay of the unchanged RWTH full-band, sparse, or truncated witnesses.

All results are newly computed in a NEW output directory and compared with saved
outputs. This is distinct from read-only verification of stored proof records.
External Bayreuth replay uses code/replay_external_fresh.py instead.
"""
from __future__ import annotations
import argparse, concurrent.futures, csv, hashlib, json, os, subprocess, sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
METHODS=[(160,2,0,'taylor160'),(256,1,1,'meanvalue256')]

def require(ok:bool,msg:str)->None:
    if not ok: raise RuntimeError(msg)

def jobs(family:str, selected:str|None):
    out=[]
    kinds=['primary','per_case','soc_union','truncated'] if family=='all-rwth' else [family]
    for kind in kinds:
        if kind=='primary':folder=ROOT/'baselines/R39/witnesses'
        elif kind in ['per_case','soc_union']:folder=ROOT/'baselines/R41/evidence/sparse'/kind
        else:folder=ROOT/'baselines/R40/results/trunc25'
        for p in sorted(folder.glob('*_witness.txt' if kind=='truncated' else '*.txt')):
            cid=p.stem.removesuffix('_witness') if kind=='truncated' else p.stem.split('_alpha')[0]
            if selected and cid!=selected:continue
            for bits,mode,shift,tag in METHODS:
                if kind=='primary':expected=ROOT/'baselines/R39/validation'/(p.stem+'_'+tag+'.json')
                elif kind=='truncated':expected=folder/(cid+'_'+tag+'_verification.json')
                else:expected=ROOT/'baselines/R41/evidence/sparse/mpfr_replay'/kind/(p.stem+f'_{bits}.json')
                require(expected.is_file(),'Missing expected proof '+str(expected))
                out.append((kind,cid,p,expected,bits,mode,shift,tag))
    require(bool(out),'No matching saved witness')
    if not selected:
        require(len(out)=={'primary':168,'per_case':168,'soc_union':168,'truncated':8,'all-rwth':512}[family], 'Unexpected witness count')
    return out

def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--family', choices=['primary','per_case','soc_union','truncated','all-rwth'],default='primary')
    ap.add_argument('--case',help='Exact case ID; omitted means all cases in the selected family')
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--output',type=Path,required=True,help='Must not already exist and must be outside this release')
    args=ap.parse_args();require(1<=args.workers<=16,'workers must be 1..16')
    require(__debug__,'Do not run with Python -O')
    dest=args.output.resolve();require(not dest.is_relative_to(ROOT.resolve()),'Output must be outside immutable release')
    require(not dest.exists(),'Refusing to reuse a directory or cached proof output')
    manifest=json.loads((ROOT/'PUBLIC_MANIFEST_SHA256.json').read_text())['files']
    tasks=jobs(args.family,args.case)
    needed={Path(t[2]).relative_to(ROOT).as_posix() for t in tasks}|{Path(t[3]).relative_to(ROOT).as_posix() for t in tasks}|{'baselines/R41/code/interval_verify.cpp'}
    for rel in needed:require(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==manifest[rel],'Input/source hash mismatch '+rel)
    dest.mkdir(parents=True);binary=dest/'interval_verify'
    src=ROOT/'baselines/R41/code/interval_verify.cpp'
    cmd=['g++','-O2','-std=c++17',str(src),'-lmpfr','-lgmp','-o',str(binary)]
    cp=subprocess.run(cmd,capture_output=True,text=True)
    if cp.returncode:
        cmd[4]='-Wl,-l:libmpfr.so.6';cp=subprocess.run(cmd,capture_output=True,text=True)
    require(cp.returncode==0,'C++/MPFR compilation failed: '+cp.stderr)
    def run(t):
        kind,cid,p,expected,bits,mode,shift,tag=t
        sub=dest/kind;sub.mkdir(exist_ok=True)
        leaf=sub/(p.stem+'_'+tag+'.csv')
        cp=subprocess.run([str(binary),str(p),str(leaf),str(bits),str(mode),str(shift)],capture_output=True,text=True,timeout=1200)
        require(cp.returncode==0 and bool(cp.stdout.strip()),'Verifier failed '+p.name+' '+tag+' '+cp.stderr)
        result=json.loads(cp.stdout);old=json.loads(expected.read_text())
        require(result==old,'Fresh output does not match saved proof '+p.name+' '+tag)
        require(result['passed'] and result['gap_lower_mOhm']>0,'Nonpositive certificate')
        rows=sorted(csv.DictReader(leaf.open()),key=lambda r:float(r['x_lo']))
        require(len(rows)==result['accepted_leaves'],'Leaf count mismatch')
        edge=Fraction(float(result['xl']))
        for row in rows:
            lo,hi=Fraction(float(row['x_lo'])),Fraction(float(row['x_hi']))
            require(lo==edge and hi>lo and float(row['residual_lower'])>0,'Incomplete or nonpositive compact cover')
            edge=hi
        require(edge==Fraction(float(result['xh'])),'Incomplete compact endpoint')
        leaf.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
        return {'family':kind,'case_id':cid,'witness':str(p.relative_to(ROOT)),'method':tag,'passed':True,'exact_saved_output_match':True}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results=list(pool.map(run,tasks))
    summary={'status':'PASS','fresh_MPFR_runs':len(results),'family':args.family,'case_filter':args.case,
      'all_outputs_identical_to_saved':True,'scope':'Fresh computation on unchanged saved inputs; not raw-data extraction, calibration or a new experiment.', 'runs':results}
    (dest/'REPLAY_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
