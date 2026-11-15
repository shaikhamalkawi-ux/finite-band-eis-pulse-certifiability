#!/usr/bin/env python3
"""MPFR directed-rounding replay of all 84 R39 duals with only their frequency supports.

This verifies two projections:
 (1) per_case: 15--21 frequencies for each of the 84 saved lower certificates;
 (2) soc_union: same 24--30 frequency set for 12 lower certificates at each SOC.
Proof limitations: subsets chosen retrospectively from full-band data and known pulse,
fixed original EIS tolerance; no prospective acquisition reduction is established.
"""
import os,sys,hashlib,csv,json,subprocess,concurrent.futures,collections,pathlib,statistics,time,math
ROOT=pathlib.Path(__file__).resolve().parents[1]
BIN=ROOT/'code'/'interval_verify';INPUT=ROOT/'evidence'/'sparse';D=ROOT/'evidence'
R39=ROOT/'evidence'
expected={}
for row in csv.DictReader((R39/'R39_FRESH_CERTIFICATES_84_LOCKED.csv').open()):
    expected[(row['case_id'],int(float(row['alpha'])))]=row
assert len(expected)==84

def worker(x):
    mode,txt,bits,method,shift=x
    name=txt.stem;alpha=int(name.split('_alpha')[1].split('_')[0]);case=name.split('_alpha')[0]
    leaf=INPUT/'mpfr_replay'/mode/(name+f'_{bits}.csv');log=leaf.with_suffix('.json');leaf.parent.mkdir(parents=True,exist_ok=True)
    if not os.environ.get('R41_FORCE_RECHECK') and log.exists() and leaf.exists():
        z=json.loads(log.read_text())
    else:
        p=subprocess.run([str(BIN),str(txt),str(leaf),str(bits),str(method),str(shift)],capture_output=True,text=True,timeout=100)
        if p.returncode!=0 or not p.stdout.strip():
            raise RuntimeError((mode,name,bits,p.returncode,p.stderr[-500:],p.stdout[-500:]))
        z=json.loads(p.stdout)
    if not z['passed']:
        raise RuntimeError((mode,name,bits,'FAILED',z))
    assert z['accepted_leaves']>0 and z['compact_residual_lower']>0
    e=expected[(case,alpha)]
    diff=z['certified_lower_mOhm']-float(e['certified_lower_mOhm'])
    # Exactly the same non-zero coefficients guarantee math identity, but numerical
    # MPFR summation may differ by last-place rounding after removing exact zeros.
    assert abs(diff)<1e-9,(mode,case,alpha,diff)
    assert z['gap_lower_mOhm']>0
    log.write_text(json.dumps(z,indent=2)+'\n')
    return {'case_id':case,'SOC':int(case.split('_')[0][3:]),'alpha':alpha,'mode':mode,'precision':bits,
            'method':('taylor160' if bits==160 else 'meanvalue256'),'frequency_count':int(txt.read_text().split()[0]),
            'passed':z['passed'],'certified_lower_mOhm':z['certified_lower_mOhm'],
            'saved_R39_certified_lower_mOhm':float(e['certified_lower_mOhm']),
            'difference_mOhm':diff,'point_gap_lower_mOhm':z['gap_lower_mOhm'],
            'compact_residual_lower':z['compact_residual_lower'],'accepted_leaves':z['accepted_leaves'],
            'source_verified_json':str(log.relative_to(ROOT))}

def main():
    jobs=[]
    for mode in ['per_case','soc_union']:
        witnesses=sorted((INPUT/mode).glob('*.txt'))
        assert len(witnesses)==84
        for p in witnesses:
            jobs.append((mode,p,160,2,0));jobs.append((mode,p,256,1,1))
    # 336 successful MPFR validations. Default verifies saved proof transcripts.
    # R41_FORCE_RECHECK=1 reruns all underlying interval-verifier programs.
    print('Checking',len(jobs),'MPFR proof runs (force fresh:',bool(os.environ.get('R41_FORCE_RECHECK')),')',flush=True)
    start=time.monotonic();rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for i,result in enumerate(pool.map(worker,jobs),1):
            rows.append(result)
            if i%28==0:print('PASS',i,'/',len(jobs),'elapsed_s',round(time.monotonic()-start,1),flush=True)
    rows=sorted(rows,key=lambda x:(x['mode'],x['case_id'],x['alpha'],x['precision']))
    assert len(rows)==336 and all(q['passed'] for q in rows)
    with (D/'R41_Sparse_MPFR_336_Results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    sizes=collections.defaultdict(list)
    for q in rows:sizes[q['mode']].append(q['frequency_count'])
    summ={'status':'PASS','fresh_MPFR_validations':len(rows),'source_R39_verified_certificates':84,
          'MPFR_methods':['taylor160','meanvalue256'],'per_case':{'min_freq':min(sizes['per_case']),
                            'median_freq':statistics.median(sizes['per_case']),'max_freq':max(sizes['per_case']),
                            'source_orig_freq':48,'certificates':84,'replay_checks':168},
          'soc_union':{'min_freq':min(sizes['soc_union']),'median_freq':statistics.median(sizes['soc_union']),
                            'max_freq':max(sizes['soc_union']),'source_orig_freq':48,
                            'certificates':84,'replay_checks':168},
          'max_abs_certified_lower_difference_mOhm':max(abs(q['difference_mOhm']) for q in rows),
          'minimum_verified_point_gap_mOhm':min(q['point_gap_lower_mOhm'] for q in rows),
          'all_positive_gaps':all(q['point_gap_lower_mOhm']>0 for q in rows),
          'non_vacuity':'Original EIS full-spectrum primal feasibility preserved after dropping constraints',
          'scope':'Retrospective EIS frequency-subset witness transport at frozen R39 full-band tolerance; not full re-optimization or prospective time saving',
          'last_verification_runtime_seconds':round(time.monotonic()-start,2),
          'last_verification_fresh_rerun':bool(os.environ.get('R41_FORCE_RECHECK'))}
    (D/'R41_Sparse_MPFR_336_Summary.json').write_text(json.dumps(summ,indent=2)+'\n')
    print('FINAL',json.dumps(summ,indent=2),flush=True)

if __name__=='__main__':main()
