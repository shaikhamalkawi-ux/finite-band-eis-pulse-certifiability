#!/usr/bin/env python3
"""R41: Lossless dual-support projection of 84 pre-verified signed lower witnesses.

For each original full-band dual witness, the frequency rows whose BOTH real and
imaginary multipliers are serialized exactly `0.0` are discarded. The remaining
rows, original tolerances, pulse waveform, endpoint and coefficients are kept
byte-exact (except first header token n -> reduced n). This is not a new solver
nor a proof of measurement-time savings, sparse optimality, or an a priori plan.

The SOC-union mode retains a common frequency-index set across 6 waveform
conditions x 2 tolerances for a given SOC, preserving all 12 certificates.
"""
from pathlib import Path
import csv,json,hashlib,collections
ROOT=Path(__file__).resolve().parents[1]
R39=ROOT/'evidence'/'sparse'
OUT=ROOT/'evidence'/'sparse'

def text_hash(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def parse_witness(p):
    lines=p.read_text().splitlines();h=lines[0].split();n=int(h[0]);pulses=int(h[1]);assert len(lines)==1+n+pulses
    freq=lines[1:1+n];pulse=lines[1+n:];idx=[]
    for j,line in enumerate(freq):
        words=line.split();assert len(words)==5
        if words[3]!='0.0' or words[4]!='0.0':
            # Some signed zeros can be literal -0.0; equivalent but require exact numeric equality.
            if float(words[3])!=0 or float(words[4])!=0:idx.append(j)
    assert 1<=len(idx)<=n
    return h,freq,pulse,idx

def main():
    sources=sorted((R39/'original_fullband').glob('*.txt'));assert len(sources)==84
    parsed={p:parse_witness(p) for p in sources};group=collections.defaultdict(list)
    for p,(h,fr,pulse,idx) in parsed.items():
        stem=p.stem; soc=stem.split('_')[0];group[soc].append((p,idx))
    assert sorted(group)==['SOC20','SOC30','SOC40','SOC50','SOC60','SOC70','SOC80']
    OUT.mkdir(exist_ok=True)
    for mode in ['per_case','soc_union']:(OUT/mode).mkdir(exist_ok=True)
    soc_idx={soc:sorted(set.union(*(set(ids) for _,ids in group[soc]))) for soc in group}
    rows=[]
    for path,(head,fr,pulse,ids) in parsed.items():
        stem=path.stem; soc=stem.split('_')[0];name=stem.rsplit('_alpha',1)[0]
        alpha=int(stem.rsplit('_alpha',1)[1].split('_')[0]);assert alpha in (1,2)
        for mode,kept in [('per_case',ids),('soc_union',soc_idx[soc])]:
            h=[str(len(kept))]+head[1:];lines=[' '.join(h)]+[fr[i] for i in kept]+pulse
            dest=OUT/mode/(stem+'.txt'); dest.write_text('\n'.join(lines)+'\n')
            assert len(lines)==1+len(kept)+int(h[1]);assert all((i in kept) for i in ids)
            rows.append({'case_id':name,'SOC':int(soc[3:]),'alpha':alpha,'mode':mode,'original_freq_n':len(fr),
                         'retained_freq_n':len(kept),'removed_freq_n':len(fr)-len(kept),
                         'removed_pct':100*(len(fr)-len(kept))/len(fr),
                         'support_original_nonzero_n':len(ids), 'retained_zero_coord_n':len(kept)-len(ids),
                         'retained_frequency_indices_0based':';'.join(map(str,kept)),
                         'retained_frequency_Hz':';'.join(fr[i].split()[0] for i in kept),
                         'original_witness_sha256':text_hash(path),'projected_witness_sha256':text_hash(dest)})
    with (OUT/'R41_Sparse_Witness_168Rows.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(sorted(rows,key=lambda d:(d['mode'],d['case_id'],d['alpha'])))
    for soc,idx in soc_idx.items():
        qlist=group[soc];p0=parsed[qlist[0][0]][1]
        record={'SOC':int(soc[3:]),'case_count':len(qlist)//2,'certificate_count':len(qlist),
                'original_frequency_count':len(p0),'union_frequency_count':len(idx),
                'zero_based_indices':idx,'frequency_Hz':[float(p0[i].split()[0]) for i in idx]}
        (OUT/f'R41_{soc}_Common_Frequency_Set.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'case_certificate_count':len(sources), 'case_support_min':min(len(v[3]) for v in parsed.values()),
          'case_support_max':max(len(v[3]) for v in parsed.values()),
          'case_support_median':sorted(len(v[3]) for v in parsed.values())[len(parsed)//2],
          'soc_common_subsets':{k:len(v) for k,v in soc_idx.items()},
          'retained_all_original_nonzero_coefficients':True},indent=2))

if __name__=='__main__':main()
