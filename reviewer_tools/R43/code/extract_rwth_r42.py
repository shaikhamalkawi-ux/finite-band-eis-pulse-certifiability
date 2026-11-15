#!/usr/bin/env python3
"""Fresh raw-waveform reconstruction for Project 52 R42.

No historical local-scatter value is an input. The 42 primary cases, 5--15 s
normalizer, 19--20.05 s endpoint window, and last-ten-sample baseline are fixed.
The 30-condition SOC50 slice uses the SAME definitions at all five temperatures.
Local scatter is a recorded-data sensitivity diagnostic, not a confidence interval.
"""
from __future__ import annotations
import argparse, hashlib, io, json, math, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]
EXPECTED='875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f'
SG={90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10}
ST={('charge',.7):3,('discharge',.7):8,('charge',1.):13,('discharge',1.):18,('charge',1.4):23,('discharge',1.4):28}
def require(ok:bool,msg:str)->None:
    if not ok: raise RuntimeError(msg)
def dump(p:Path,obj)->None:
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def robust_scale(y):
    y=np.asarray(y,float)
    return float(1.4826*np.median(np.abs(y-np.median(y))))
def detrended_scale(t,y):
    t=np.asarray(t,float);y=np.asarray(y,float)
    require(len(t)>=3 and np.ptp(t)>0,'Insufficient detrending samples')
    x=t-np.mean(t);v=y-np.mean(y)
    slope=float(np.dot(x,v)/np.dot(x,x))
    residual=v-slope*x
    return robust_scale(residual),slope

def load_overview(z,temp):
    name=f'{temp}deg/{temp}deg_Overview.mat';raw=z.read(name)
    d=loadmat(io.BytesIO(raw),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
    a={k:np.asarray(getattr(d,k)).ravel() for k in ['Zeit','Strom','Spannung','Schritt']}
    t,I,V,step=[a[k] for k in ['Zeit','Strom','Spannung','Schritt']]
    require(all(np.isfinite(a[k]).all() for k in a),'Nonfinite source')
    on=np.abs(I)>.5;d=np.diff(on.astype(int));ss=np.flatnonzero(d==1)+1;ee=np.flatnonzero(d==-1)+1
    if on[0]:ss=np.r_[0,ss]
    if on[-1]:ee=np.r_[ee,len(on)]
    require(len(ss)==len(ee),'Segmentation boundary mismatch')
    segments=[(int(s),int(e),int(step[s])) for s,e in zip(ss,ee) if 15<=t[e-1]-t[s]<=25 and abs(np.median(I[s:e]))>1.5]
    groups=[];cur=[];prev=None
    for q in segments:
        if prev is not None and q[0]-prev>10000:groups.append(cur);cur=[]
        cur.append(q);prev=q[0]
    groups.append(cur)
    return a,groups,{'file':name,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'group_sizes':list(map(len,groups))}

def load_eis(z,temp,soc):
    name=f'{temp}deg/{temp}deg_SOC{soc}.mat';raw=z.read(name)
    d=loadmat(io.BytesIO(raw),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
    f,zr,zi=[np.asarray(getattr(d,k)).ravel() for k in ['ActFreq','Zreal1','Zimg1']]
    m=(f>0)&((zr!=0)|(zi!=0));F=np.unique(f[m]);ix=np.array([np.flatnonzero(m&(f==x))[0] for x in F])
    require(len(ix)==48,'Unexpected EIS frequency count')
    return {'frequency_Hz':f[ix].tolist(),'Zreal_mOhm':zr[ix].tolist(),'Zimag_mOhm':zi[ix].tolist(),'source_EIS_file':name,'source_EIS_sha256':hashlib.sha256(raw).hexdigest()}

def extract(a,groups,info,eis,temp,soc,direction,rate,eps):
    t,I,V=a['Zeit'],a['Strom'],a['Spannung']
    matches=[x[:2] for x in groups[SG[soc]] if x[2]==ST[(direction,rate)]]
    require(len(matches)==1,'Ambiguous pulse identifier')
    s,e=matches[0];tt=t[s:e]-t[s];ii=I[s:e];vv=V[s:e]
    require(tt[0]==0 and np.max(tt)>=20,'Active pulse does not bracket 20 s')
    # Record-order duplicates are real source records. Preserve them for scatter;
    # use stable timestamp order and the last record at each repeated time for
    # interpolation and the ZOH current kernel. At 5 C this is exactly the
    # archived convention of skipping zero-duration intervals.
    order=np.argsort(tt,kind='stable');stt=tt[order]
    retain=np.r_[stt[1:]!=stt[:-1],True];ix=order[retain]
    tc,ic,vc=tt[ix],ii[ix],vv[ix]
    require(np.all(np.diff(tc)>0),'Canonical timestamps not strictly increasing')
    pre=np.arange(max(0,s-200),s);pre=pre[np.abs(I[pre])<=.5][-10:]
    require(len(pre)==10 and np.all(I[pre]==0),'Nonzero or insufficient baseline')
    pl=(tt>=5)&(tt<=15);end=(tt>=19)&(tt<=20.05)
    sg=1 if direction=='charge' else -1
    Ip=abs(float(np.median(ii[pl])));v0=float(np.median(V[pre]));vt=float(np.interp(20.,tc,vc));dv=float(sg*(vt-v0));legacy_vt=float(np.interp(20.,tt,vv))
    require(Ip>0 and sg*np.median(ii[pl])>0 and dv>0,'Direction/normalizer mismatch')
    norm=sg*ic/Ip
    aa=[];bb=[];u=[]
    for j in range(len(tc)-1):
        x=float(tc[j]);y=min(float(tc[j+1]),20.)
        if x>=20:break
        if y<=x:continue
        aa.append(x);bb.append(y);u.append(float(norm[j]))
        if y==20:break
    k=int(np.searchsorted(tc,20.,side='right')-1);iT=float(norm[k])
    s0=robust_scale(V[pre]);sv,slopeV=detrended_scale(tt[end],vv[end]);si,slopeI=detrended_scale(tt[end],ii[end]);sip,slopeIp=detrended_scale(tt[pl],ii[pl])
    vterm=1000*np.hypot(s0,sv)/Ip
    iterm=1000*dv*si/(Ip*Ip);pterm=1000*dv*sip/(Ip*Ip)
    h=float(3*np.hypot(vterm,iterm));hp=float(3*np.hypot(vterm,pterm))
    cid=f'SOC{soc}_{direction}_{rate:.1f}C';tid=f'T{temp}_{cid}'
    q={'case_id':cid,'temperature_C':temp,'SOC':soc,'direction':direction,'C_rate':rate,**eis,'T_s':20.,'a_s':aa,'b_s':bb,'u':u,'iT':iT,'eps_mOhm':eps,'Rmeas_mOhm':1000*dv/Ip,'normalization_current_A':Ip,'Q_T_float_s':float(np.sum(np.asarray(u)*(np.asarray(bb)-aa))),'source_overview_file':info['file'],'source_overview_sha256':info['sha256'],'raw_start_index':s,'raw_stop_index':e}
    windows={'trace_id':tid,'baseline_t_s':(t[pre]-t[s]).tolist(),'baseline_V':V[pre].tolist(),'baseline_I_A':I[pre].tolist(),'pulse_t_s':tt.tolist(),'pulse_V':vv.tolist(),'pulse_I_A':ii.tolist(),'normalizer_window_s':[5.,15.],'endpoint_scatter_window_s':[19.,20.05],'baseline_indices':pre.tolist(),'source_overview_sha256':info['sha256'],'raw_start_index':s,'raw_stop_index':e,'canonical_source_indices_within_pulse':ix.tolist()}
    stats={'trace_id':tid,'case_id':cid,'temperature_C':temp,'SOC':soc,'direction':direction,'C_rate':rate,'I_p_A':Ip,'Rmeas_mOhm':q['Rmeas_mOhm'],'sigma_baseline_V':s0,'sigma_endpoint_V':sv,'sigma_endpoint_I_A':si,'sigma_plateau_I_A':sip,'endpoint_voltage_slope_V_s':slopeV,'endpoint_current_slope_A_s':slopeI,'plateau_current_slope_A_s':slopeIp,'voltage_contribution_mOhm':vterm,'endpoint_current_contribution_mOhm':iterm,'plateau_current_contribution_mOhm':pterm,'halfwidth_3sigma_endpoint_mOhm':h,'halfwidth_3sigma_plateau_mOhm':hp,'baseline_n':len(pre),'endpoint_n':int(end.sum()),'plateau_n':int(pl.sum()),'baseline_V':v0,'endpoint_interpolated_V':vt,'last_active_sample_s':float(tt[-1]),'duplicate_timestamp_steps':int((np.diff(tt)==0).sum()),'backward_timestamp_steps':int((np.diff(tt)<0).sum()),'legacy_order_endpoint_V':legacy_vt,'canonical_minus_legacy_R_mOhm':float(1000*sg*(vt-legacy_vt)/Ip)}
    return q,windows,stats

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,required=True);args=ap.parse_args()
    require(hashlib.sha256(args.raw.read_bytes()).hexdigest()==EXPECTED,'Original archive SHA256 mismatch')
    old=pd.read_csv(ROOT/'inputs/R4D_SOC50_Temperature_Historical.csv')
    oldcases={p.stem:json.loads(p.read_text()) for p in (ROOT/'inputs/R39_cases').glob('*.json')}
    oldcor=pd.read_csv(ROOT/'inputs/RAW_CURRENT_CORRECTION_42.csv').set_index('case_id')
    cert=pd.read_csv(ROOT/'inputs/R39_FRESH_CERTIFICATES_84.csv')
    out=ROOT/'evidence';source=[];stats=[];reconcile=[];tempout=[]
    with zipfile.ZipFile(args.raw) as z:
        require(z.testzip() is None,'Raw ZIP CRC failed')
        for temp in [5,15,25,35,45]:
            a,groups,info=load_overview(z,temp);source.append(info)
            socs=range(20,81,10) if temp==5 else [50]
            for soc in socs:
                ei=load_eis(z,temp,soc);source.append({'file':ei['source_EIS_file'],'sha256':ei['source_EIS_sha256']})
                for direction,rate in ST:
                    cid=f'SOC{soc}_{direction}_{rate:.1f}C'
                    eps=oldcases[cid]['eps_mOhm'] if temp==5 else float(old[old.temp==temp].eps_mOhm.iloc[0])
                    q,windows,st=extract(a,groups,info,ei,temp,soc,direction,rate,eps)
                    dump(out/'source_windows'/(st['trace_id']+'.json'),windows);dump(out/'fresh_inputs'/(st['trace_id']+'.json'),q);stats.append(st)
                    if temp==5:
                        oldq=oldcases[cid];errs={}
                        for k in ['frequency_Hz','Zreal_mOhm','Zimag_mOhm','a_s','b_s','u','iT','T_s','Rmeas_mOhm','eps_mOhm']:
                            v=np.asarray(q[k]);w=np.asarray(oldq[k]);require(v.shape==w.shape,'Input shape mismatch '+cid+' '+k);err=float(np.max(abs(v-w)));errs[k]=err
                            require(err<3e-11,'R39 raw reconstruction mismatch '+cid+' '+k)
                        reconcile.append({'case_id':cid,**errs,'max_abs_difference':max(errs.values())})
                    if soc==50:
                        r=old[(old.temp==temp)&(old.direction==direction)&np.isclose(old.C_rate,rate)].iloc[0]
                        scale=1000*abs(st['legacy_order_endpoint_V']-st['baseline_V'])/(float(r.Rmeas_mOhm)*q['normalization_current_A'])
                        tempout.append({'trace_id':st['trace_id'],'temp':temp,'SOC':soc,'direction':direction,'C_rate':rate,'eps_mOhm':eps,'old_Rmeas_mOhm':r.Rmeas_mOhm,'measured_current_A':q['normalization_current_A'],'new_Rmeas_mOhm':q['Rmeas_mOhm'],'old_to_new_scale':scale,'nominal_lower_mOhm':float(r.lower_mOhm)*scale,'nominal_upper_mOhm':float(r.upper_mOhm)*scale,'nominal_gap_mOhm':max(0.,float(r.lower_mOhm)*scale-q['Rmeas_mOhm'],q['Rmeas_mOhm']-float(r.upper_mOhm)*scale),'halfwidth_3sigma_mOhm':st['halfwidth_3sigma_endpoint_mOhm'],'canonical_minus_legacy_R_mOhm':st['canonical_minus_legacy_R_mOhm']})
    require(len(stats)==66 and len(reconcile)==42 and len(tempout)==30,'Unexpected case coverage')
    df=pd.DataFrame(stats);df.to_csv(out/'R42_Raw_Scatter_66_Conditions.csv',index=False,float_format='%.17g')
    primary=df[df.temperature_C==5].copy();primary.to_csv(out/'R42_Raw3sigma_42.csv',index=False,float_format='%.17g')
    d=cert.merge(primary,on='case_id',validate='many_to_one')
    for kind in ['endpoint','plateau']:
        d[f'gap_after_3sigma_{kind}_mOhm']=np.nextafter(d.point_gap_lower_mOhm-d[f'halfwidth_3sigma_{kind}_mOhm'],-np.inf)
    d.to_csv(out/'R42_Raw3sigma_Certificate_Gaps_84.csv',index=False,float_format='%.17g')
    pd.DataFrame(reconcile).to_csv(out/'R42_R39_Raw_Input_Reconciliation_42.csv',index=False,float_format='%.17g')
    td=pd.DataFrame(tempout).sort_values(['temp','C_rate','direction']);td.to_csv(out/'R42_Temperature_SOC50_Corrected_30.csv',index=False,float_format='%.17g')
    s={'status':'PASS','raw_sha256':EXPECTED,'raw_source_crc':'PASS','primary_cases':42,'temperature_slice_cases':30,'unique_source_conditions':66,'source_members':source,'primary_source_max_input_difference':max(r['max_abs_difference'] for r in reconcile),'local_3sigma':{},'scope':'Raw recorded scatter; diagonal first-order sensitivity, not absolute calibration, confidence coverage, or full-waveform uncertainty.'}
    for kind in ['endpoint','plateau']:
        s['local_3sigma'][kind]={'median_halfwidth_mOhm':float(primary[f'halfwidth_3sigma_{kind}_mOhm'].median()),'max_halfwidth_mOhm':float(primary[f'halfwidth_3sigma_{kind}_mOhm'].max()),'by_alpha':{str(alpha):{'positive_cases':int((g[f'gap_after_3sigma_{kind}_mOhm']>0).sum()),'minimum_gap_mOhm':float(g[f'gap_after_3sigma_{kind}_mOhm'].min())} for alpha,g in d.groupby('alpha')}}
    dump(out/'R42_Raw_Reconstruction_Summary.json',s)
    print(json.dumps({k:v for k,v in s.items() if k!='source_members'},indent=2))
    print(td.groupby(['temp','direction']).nominal_gap_mOhm.agg(['median','min','max']).to_string())
if __name__=='__main__':main()
