#!/usr/bin/env python3
"""Raw-source reconstruction and positive-homogeneous transport of R36 dual witnesses.
Scientific R39 candidate: no optimizer; MPFR revalidation is separate. Historical data remain untouched.
"""
from __future__ import annotations
import io,hashlib,json,zipfile
from pathlib import Path
import numpy as np,pandas as pd,scipy.io as sio
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'replay/Batteries52_Reviewer_Replay/R36_reproducibility'
RAW=Path('/mnt/data/MatFiles.zip')
assert hashlib.sha256(RAW.read_bytes()).hexdigest()=='875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f'
OUT=ROOT/'r39_candidate';(OUT/'inputs').mkdir(parents=True,exist_ok=True);(OUT/'witnesses').mkdir(parents=True,exist_ok=True);(OUT/'candidates').mkdir(parents=True,exist_ok=True)
oldtab=pd.read_csv(SRC/'baseline/R4D_RWTH_5C_ContinuumSafe_42cases.csv');oldcert=pd.read_csv(SRC/'FINAL_ALL42_CERTIFICATES.csv')
with zipfile.ZipFile(RAW) as z:
 assert z.testzip() is None
 ov=sio.loadmat(io.BytesIO(z.read('5deg/5deg_Overview.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
t=np.asarray(ov.Zeit).squeeze();I=np.asarray(ov.Strom).squeeze();V=np.asarray(ov.Spannung).squeeze();steps=np.asarray(ov.Schritt).squeeze()
on=np.abs(I)>.5;d=np.diff(on.astype(int));starts=np.where(d==1)[0]+1;ends=np.where(d==-1)[0]+1
if on[0]:starts=np.r_[0,starts]
if on[-1]:ends=np.r_[ends,len(on)]
segs=[]
for s,e in zip(starts,ends):
 if 15<=t[e-1]-t[s]<=25 and abs(np.median(I[s:e]))>1.5:segs.append((s,e,int(steps[s])))
groups=[];cur=[];prev=None
for x in segs:
 if prev is not None and x[0]-prev>10000:groups.append(cur);cur=[]
 cur.append(x);prev=x[0]
groups.append(cur)
SOCG={90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10}
STEP={('charge',0.7):3,('discharge',0.7):8,('charge',1.0):13,('discharge',1.0):18,('charge',1.4):23,('discharge',1.4):28}
summary=[];worst=0
for row in oldtab.itertuples():
 soc=int(row.SOC);rate=float(row.C_rate);direction=row.direction;case=f'SOC{soc}_{direction}_{rate:.1f}C'
 s,e=[x[:2] for x in groups[SOCG[soc]] if x[2]==STEP[(direction,rate)]][0]
 tt=t[s:e]-t[s];ii=I[s:e];vv=V[s:e];pre=np.arange(max(0,s-200),s);pre=pre[np.abs(I[pre])<=.5]
 base=V[pre[-10:]];v0=float(np.median(base));VT=float(np.interp(20.,tt,vv));direction_sign=np.sign(np.median(ii[(tt>=5)&(tt<=15)]));Ip=float(abs(np.median(ii[(tt>=5)&(tt<=15)])))
 assert Ip>0 and direction_sign*(VT-v0)>0
 measured=float(1000*direction_sign*(VT-v0)/Ip)
 oldq=json.loads((SRC/'inputs'/(case+'.json')).read_text());oldip=float(1000*abs(VT-v0)/row.Rmeas_mOhm);scale=oldip/Ip
 norm=ii/(direction_sign*Ip)
 aa=[];bb=[];u=[]
 for j in range(len(tt)-1):
  a=float(tt[j]);b=float(tt[j+1]);
  if a>=20.:break
  b=min(b,20.)
  if b<=a:continue
  aa.append(a);bb.append(b);u.append(float(norm[j]))
  if b>=20.:break
 jt=np.searchsorted(tt,20.,side='right')-1;iT=float(norm[max(0,min(jt,len(norm)-1))]);Q=float(np.sum(np.array(u)*(np.array(bb)-np.array(aa))))
 assert aa==oldq['a_s'] and bb==oldq['b_s']
 assert len(u)==len(oldq['u']) and all(np.isfinite(u)) and min(u)>0
 err=np.max(abs(np.array(u)-np.array(oldq['u'])*scale))
 assert err<5e-14,(case,err)
 assert abs(measured-oldq['Rmeas_mOhm']*scale)<3e-11
 assert abs(iT-oldq['iT']*scale)<5e-14
 q={**oldq,'Rmeas_mOhm':measured,'iT':iT,'Q_T_float_s':Q,'u':u,
 'normalization_current_A':Ip,'current_reference':'absolute median measured pulse current on 5 <= t <= 15 seconds',
 'current_reference_interval_s':[5.,15.], 'voltage_baseline_last10_V':v0,
 'voltage_endpoint_interpolated_T20_V':VT,'old_frozen_Rmeas_mOhm':float(row.Rmeas_mOhm),
 'old_inferred_current_A':oldip,'old_to_new_rescaling_factor':scale}
 (OUT/'inputs'/(case+'.json')).write_text(json.dumps(q,indent=2)+'\n')
 # Fixed endpoint recorded-sample envelope: ten active records preceding T,
 # plus first active record at or after T. No outcome-based window selection.
 before=np.flatnonzero(tt<=20.)[-10:]
 after=np.flatnonzero(tt>=20.)[:1]
 en=vv[np.unique(np.r_[before,after])]
 assert len(en)>=2 and en.min()-1e-12<=VT<=en.max()+1e-12,(case,len(en),VT,en.min(),en.max())
 if direction_sign>0:
  dmin=float(en.min()-base.max());dmax=float(en.max()-base.min())
 else:
  dmin=float(base.min()-en.max());dmax=float(base.max()-en.min())
 lower=1000*dmin/Ip;upper=1000*dmax/Ip
 assert lower<=measured<=upper,(case,lower,measured,upper)
 rowdict={'case_id':case,'SOC':soc,'direction':direction,'C_rate':rate,'Ip_measured_median_5to15_A':Ip,'Ip_old_inferred_A':oldip,
 'rescale':scale,'R_old_mOhm':row.Rmeas_mOhm,'R_new_measured_normalizer_mOhm':measured,
 'raw_voltage_sample_envelope_low_mOhm':lower,'raw_voltage_sample_envelope_high_mOhm':upper,
 'number_endpoint_samples':len(en),'old_lower_mOhm':row.lower_mOhm,'old_upper_mOhm':row.upper_mOhm,
 'source_mat_sha256':hashlib.sha256(zipfile.ZipFile(RAW).read('5deg/5deg_Overview.mat')).hexdigest()}
 summary.append(rowdict)
 for r in oldcert[oldcert.case_id==case].itertuples():
  stem=Path(r.witness_file).stem
  c=json.loads((SRC/'candidates'/(stem+'.json')).read_text());assert c['case_id']==case
  y=[float(float(z)*scale) for z in c['y']]
  cnew={**c,'y':y,'source_witness_sha256':hashlib.sha256((SRC/r.witness_file).read_bytes()).hexdigest(),'scaling_factor_old_to_new':scale,
        'candidate_method':'Exact linear-cone rescaling of previously feasible signed dual, followed by fresh directed-rounding verification; no solver optimization'}
  (OUT/'candidates'/(stem+'.json')).write_text(json.dumps(cnew,indent=2)+'\n')
  n=len(q['frequency_Hz']);lines=[' '.join(map(str,[n,len(u),q['T_s'],q['iT'],q['eps_mOhm'],q['Rmeas_mOhm'],c['alpha']]))]
  for k in range(n):lines.append(' '.join(repr(float(z)) for z in [q['frequency_Hz'][k],q['Zreal_mOhm'][k],q['Zimag_mOhm'][k],y[k],y[k+n]]))
  for a,b,z in zip(aa,bb,u):lines.append(' '.join(repr(float(v)) for v in [a,b,z]))
  (OUT/'witnesses'/(stem+'.txt')).write_text('\n'.join(lines)+'\n')

df=pd.DataFrame(summary).sort_values(['SOC','direction','C_rate']);df.to_csv(OUT/'RAW_CURRENT_CORRECTION_42.csv',index=False,float_format='%.17g')
print('reconstructed all',len(df),'source cases and',len(list((OUT/'witnesses').glob('*.txt'))),'witnesses')
print('old/new R',df[['R_old_mOhm','R_new_measured_normalizer_mOhm']].describe());print('3 affected',df.assign(diff=df.R_new_measured_normalizer_mOhm-df.R_old_mOhm).sort_values('diff').head(5)[['case_id','Ip_measured_median_5to15_A','Ip_old_inferred_A','R_old_mOhm','R_new_measured_normalizer_mOhm','diff']].to_string(index=False));print('minimum raw voltage envelope check', (df.raw_voltage_sample_envelope_high_mOhm-df.R_new_measured_normalizer_mOhm).min())
