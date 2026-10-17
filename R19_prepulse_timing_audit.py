#!/usr/bin/env python3
"""R19 pre-pulse timing audit for the 42 primary RWTH 5 C pulse cases.

Requires the official RWTH MatFiles.zip archive (not redistributed here).
The active-pulse threshold |I| > 0.5 A matches the segmentation used by the
representative continuum-validation script.
"""
import argparse, io, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io as sio

PULSE_THRESHOLD_A = 0.5
SOCS = [20,30,40,50,60,70,80]
SOC_TO_GROUP = {90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10}
STEP_MAP = {
    ('charge',0.7):3, ('discharge',0.7):8,
    ('charge',1.0):13, ('discharge',1.0):18,
    ('charge',1.4):23, ('discharge',1.4):28,
}

def load_overview(path):
    with zipfile.ZipFile(path) as z:
        ov=sio.loadmat(io.BytesIO(z.read('5deg/5deg_Overview.mat')),
                       squeeze_me=False, struct_as_record=False)['diga'][0,0].daten[0,0]
    return (np.asarray(ov.Zeit).squeeze().astype(float),
            np.asarray(ov.Strom).squeeze().astype(float),
            np.asarray(ov.Spannung).squeeze().astype(float),
            np.asarray(ov.Schritt).squeeze())

def pulse_groups(t,I,V,steps):
    on=np.abs(I)>PULSE_THRESHOLD_A
    d=np.diff(on.astype(int))
    starts=np.where(d==1)[0]+1; ends=np.where(d==-1)[0]+1
    if on[0]: starts=np.r_[0,starts]
    if on[-1]: ends=np.r_[ends,len(on)]
    segs=[]
    for s,e in zip(starts,ends):
        dur=t[e-1]-t[s]; med=np.median(I[s:e])
        if 15<=dur<=25 and abs(med)>1.5:
            segs.append((s,e,dur,med,int(steps[s])))
    groups=[]; cur=[]; prev=None
    for sg in segs:
        if prev is not None and sg[0]-prev>10000:
            groups.append(cur); cur=[]
        cur.append(sg); prev=sg[0]
    groups.append(cur)
    return groups

def run(matfiles):
    t,I,V,steps=load_overview(matfiles)
    groups=pulse_groups(t,I,V,steps)
    rows=[]
    for soc in SOCS:
        g=groups[SOC_TO_GROUP[soc]]
        for direction in ['charge','discharge']:
            for rate in [0.7,1.0,1.4]:
                st=STEP_MAP[(direction,rate)]
                s,e,dur,med,_=[x for x in g if x[4]==st][0]
                j=s-1
                while j>=0 and abs(I[j])<=PULSE_THRESHOLD_A:
                    j-=1
                rest_start=j+1
                prezero=t[s]-t[rest_start]
                idx=np.arange(s-10,s)
                slope=np.polyfit(t[idx]-t[idx][0],V[idx],1)[0]*1000.0
                vspan=(V[idx].max()-V[idx].min())*1000.0
                tt=t[s:e]-t[s]
                Iend=np.median(I[e-10:e])
                I20=float(np.interp(20.0,tt,I[s:e]))
                rows.append({
                    'SOC_pct':soc,'direction':direction,'C_rate':rate,
                    'pre_zero_current_duration_s':prezero,
                    'baseline10_span_s':t[s-1]-t[s-10],
                    'baseline10_abs_current_max_A':float(np.max(np.abs(I[idx]))),
                    'baseline10_voltage_slope_mV_per_s':slope,
                    'baseline10_voltage_range_mV':vspan,
                    'active_last_sample_rel_s':t[e-1]-t[s],
                    'first_postpulse_sample_rel_s':(t[e]-t[s]) if e<len(t) else np.nan,
                    'I_plateau_last10_A':Iend,
                    'I_at_20s_active_interp_A':I20,
                    'rel_I20_vs_plateau':abs(I20-Iend)/max(abs(Iend),1e-12),
                    '20s_within_active_sample_span':bool(20.0<=tt[-1]+1e-12),
                })
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('matfiles_zip')
    ap.add_argument('--out',default='R19_RWTH_5C_PrePulse_Timing_Audit.csv')
    args=ap.parse_args()
    df=run(args.matfiles_zip)
    df.to_csv(args.out,index=False)
    print(f'n={len(df)}')
    print(f'min pre-zero duration [s]={df.pre_zero_current_duration_s.min():.6f}')
    print(f'median pre-zero duration [s]={df.pre_zero_current_duration_s.median():.6f}')
    print(f'median |baseline slope| [mV/s]={df.baseline10_voltage_slope_mV_per_s.abs().median():.6g}')
    print(f'max |baseline slope| [mV/s]={df.baseline10_voltage_slope_mV_per_s.abs().max():.6f}')
    print(f'max baseline voltage range [mV]={df.baseline10_voltage_range_mV.max():.6f}')
    print(f'20 s inside active span={int(df["20s_within_active_sample_span"].sum())}/{len(df)}')
    print(f'active last sample rel range [s]={df.active_last_sample_rel_s.min():.6f}..{df.active_last_sample_rel_s.max():.6f}')
    print(f'first post-pulse sample rel range [s]={df.first_postpulse_sample_rel_s.min():.6f}..{df.first_postpulse_sample_rel_s.max():.6f}')
    print(f'max relative I20-vs-plateau difference={df.rel_I20_vs_plateau.max():.6g}')
if __name__=='__main__': main()
