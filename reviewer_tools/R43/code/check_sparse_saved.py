#!/usr/bin/env python3
"""R41 scientific / production integrity QA; strict and non-mutating.
Test only frozen numerical inputs, stored MPFR proofs, source and generated PDFs.
No claim to verify raw-source extraction, hardware calibration or prospectively
measured frequency-subset performance.
"""
from pathlib import Path
import csv,json,hashlib,statistics,collections,re,math,subprocess,sys
R=Path(__file__).resolve().parents[1]/'baselines/R41'; E=R/'evidence'; S=E/'sparse'; SRC=R/'overleaf'
log=[]
def check(flag,name):
 if not flag:raise AssertionError('FAILED: '+name)
 log.append({'test':name,'status':'PASS'})
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
def rows(p):return list(csv.DictReader(p.open()))
# immutable R39 corrected-current processed instances
casefiles=sorted((E/'inputs').glob('*.json'));check(len(casefiles)==42,'exactly 42 corrected-current R39 processed input cases')
cases=[json.loads(p.read_text()) for p in casefiles]
check(len({q['case_id'] for q in cases})==42,'all 42 pulse case IDs are distinct')
check(all('5 <= t <= 15 seconds' in q['current_reference'] for q in cases),'same measured 5-15-second current-reference definition')
check(all(len(q['frequency_Hz'])==48 for q in cases),'all spectra have 48 original frequencies')
# finite-time atom values (recapitulation, NOT raw-source repro)
atom=rows(E/'R41_FiniteAtom_42Cases.csv');check(len(atom)==42,'42 finite-atom case result rows')
summ=json.loads((E/'R41_Math_Reproduction_Summary.json').read_text());check(summ['N']==42,'finite-atom summary N=42')
check(summ['tests_common_uniform_48_frequency']==3444,'3,444 checked direct frequency-cap formula configurations')
check(float(summ['tests']['zoh_independent_max_abs'])<1e-60,'two high-precision ZOH formulas agree')
check(float(summ['tests']['closed_rhide_max_abs'])<1e-60,'piecewise finite atom cap matches all-frequency direct min at 3,444 points')
check(float(summ['tests']['tau_star_atom_check_max_abs'])<1e-60,'one-point cap 2 epsilon at tau-star')
slow=[float(q['R41_slow_limit_mOhm']) for q in atom]
finite=[float(q['R41_finite_atom_onepoint_mOhm']) for q in atom]
check(len({q['case_id'] for q in atom})==42,'finite-atom case IDs distinct')
check(max(slow)<.4 and min(finite)>.4,'0/42 slow and 42/42 finite screens at 0.4 mOhm')
check(math.isclose(statistics.median(slow),summ['slow_median_mOhm'],abs_tol=1e-11),'slow-median matches summary')
check(math.isclose(statistics.median(finite),summ['finite_median_mOhm'],abs_tol=1e-11),'finite-median matches summary')
check(summ['theorem_scope'].startswith('worst-case over allowed'),'critical synthetic-observation-centers scope retained')
check(summ['nonuniform_counterexample']['actual_rhide']<summ['nonuniform_counterexample']['false_fmin_only_prediction'],'nonuniform-tolerance counterexample passes')
# original and compressed signed-dual witness structural exactness
orig=sorted((S/'original_fullband').glob('*.txt'))
check(len(orig)==84,'84 original R39 signed dual witnesses included')
source=rows(E/'R39_FRESH_CERTIFICATES_84_LOCKED.csv')
check(len(source)==84,'84 source R39 certified output records')
source_by={ (v['case_id'],int(v['alpha'])):v for v in source}
check(len(source_by)==84,'one source certificate per 42 case x 2 tolerances')
mpfr=rows(E/'R41_Sparse_MPFR_336_Results.csv')
check(len(mpfr)==336,'336 MPFR proof-result rows recorded')
mpfr_by={(q['case_id'],int(q['alpha']),q['mode'],int(q['precision'])):q for q in mpfr}
check(len(mpfr_by)==336,'336 MPFR combinations unique')
check(all(q['passed']=='True' for q in mpfr),'all 336 stored MPFR results passed')
check(all(int(q['accepted_leaves'])>0 and float(q['compact_residual_lower'])>0 for q in mpfr),'all MPFR proofs show full accepted domain with positive minimum')
check(all(float(q['point_gap_lower_mOhm'])>0 for q in mpfr),'positive gaps at all 336 proofs')
check(max(abs(float(q['difference_mOhm'])) for q in mpfr)<1e-9,'all projected bound evaluations agree with R39')
proj=rows(S/'R41_Sparse_Witness_168Rows.csv');check(len(proj)==168,'168 frequency-projected witness records')
proj_by={(q['case_id'],int(q['alpha']),q['mode']):q for q in proj}
check(len(proj_by)==168,'168 unique witness projections')
by_soc=collections.defaultdict(list)
per_n=[];union_n=[]
for p in orig:
 stem=p.stem;case=stem.split('_alpha')[0];alpha=int(stem.split('_alpha')[1].split('_')[0]);soc=int(case.split('_')[0][3:]);
 src=source_by[(case,alpha)]
 check(sha(p)==src['witness_sha256'],'R39 source witness SHA256 '+stem)
 l=p.read_text().splitlines();h=l[0].split();n=int(h[0]);pulses=int(h[1]);check(n==48 and len(l)==1+n+pulses,'original witness structure '+stem)
 freq=l[1:1+n];wave=l[1+n:]
 support=[i for i,x in enumerate(freq) if float(x.split()[3])!=0 or float(x.split()[4])!=0]
 check(15<=len(support)<=21,'original dual sparse support in 15-21 range '+stem)
 by_soc[soc].append((stem,support))
 for mode in ['per_case','soc_union']:
  newp=S/mode/p.name;v=proj_by[(case,alpha,mode)]
  new=newp.read_text().splitlines();kept=[int(z) for z in v['retained_frequency_indices_0based'].split(';')]
  check(len(new)==1+len(kept)+pulses,'projected rows '+mode+' '+stem)
  check(new[0].split()==[str(len(kept))]+h[1:],'projected header preserves physics '+mode+' '+stem)
  check(new[1:1+len(kept)]==[freq[i] for i in kept],'retained EIS measurements/coefficients identical '+mode+' '+stem)
  check(new[1+len(kept):]==wave,'waveform rows identical '+mode+' '+stem)
  check(all(i in kept for i in support),'every nonzero weight retained '+mode+' '+stem)
  check(all(float(freq[i].split()[3])==0 and float(freq[i].split()[4])==0 for i in range(n) if i not in kept),'all removed coordinates are exact-zero '+mode+' '+stem)
  if mode=='soc_union': check(kept[0]==0,'lowest frequency retained by common SOC union '+stem)
  check(sha(newp)==v['projected_witness_sha256'],'projected witness SHA256 '+mode+' '+stem)
  for bits in [160,256]:
   rec=mpfr_by[(case,alpha,mode,bits)];j= S/'mpfr_replay'/mode/(stem+'_'+str(bits)+'.json');leaf=j.with_suffix('.csv')
   check(j.is_file() and leaf.is_file(),'saved output plus domain enclosure leaves '+mode+' '+str(bits)+' '+stem)
   result=json.loads(j.read_text())
   check(result['passed'] and result['accepted_leaves']>0 and result['compact_residual_lower']>0,'saved MPFR certified complete domain '+mode+' '+str(bits)+' '+stem)
   check(abs(result['certified_lower_mOhm']-float(src['certified_lower_mOhm']))<1e-9,'saved MPFR lower equals source R39 '+mode+' '+str(bits)+' '+stem)
  (per_n if mode=='per_case' else union_n).append(len(kept))
check(min(per_n)==15 and statistics.median(per_n)==18 and max(per_n)==21,'individual witness support sizes exactly 15-21; median 18')
check(min(union_n)==24 and statistics.median(union_n)==28 and max(union_n)==30,'SOC union support sizes exactly 24-30; median 28')
for soc,ss in sorted(by_soc.items()):
 intended=sorted(set(i for _,indices in ss for i in indices))
 j=json.loads((S/f'R41_SOC{soc}_Common_Frequency_Set.json').read_text())
 check(j['certificate_count']==12 and len(ss)==12,'twelve certificates per SOC '+str(soc))
 check(j['zero_based_indices']==intended and j['union_frequency_count']==len(intended),'exact union of 12 supports '+str(soc))
check(len(by_soc)==7,'seven SOC groups with valid union records')
mpfrsum=json.loads((E/'R41_Sparse_MPFR_336_Summary.json').read_text())
check(mpfrsum['status']=='PASS' and mpfrsum['last_verification_fresh_rerun'],'336 MPFR proofs all rerun freshly')
check(mpfrsum['max_abs_certified_lower_difference_mOhm']==0,'MPFR lower bound difference recorded as zero')
check((E/'R41_336_replay_fullfresh.log').read_text().count('PASS 336 / 336')==1,'full-fresh MPFR terminal execution log recorded')
# Nonvacuuity proof from previously independently replayed 7 full-band feasible primal witnesses.
primal=json.loads((E/'R41_Seven_FullBand_Primal_Nonvacuity_Replay.json').read_text());check(len(primal)==7,'seven SOC primal certificates in archive')
check(sum(v['eis_coordinates'] for v in primal)==672,'672 full-band EIS coordinates in primal summary')
check(all(v['positive_margin_verified'] and v['minimum_margin_mOhm']>0 for v in primal),'all primal witness constraints have positive margin')
check(len(list((E/'primal').glob('SOC*_primal.json')))==7,'all seven source primal JSON witnesses supplied')

print(json.dumps({'status':'PASS','stored_scientific_checks_passed':len(log),'scope':'Read-only exact-zero support identity and stored-proof audit; no new MPFR execution.'},indent=2))
