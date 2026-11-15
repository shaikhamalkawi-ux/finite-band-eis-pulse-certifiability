#!/usr/bin/env python3
"""Read original publisher ZIPs; preserve provenance; export matched E3 pulse inputs."""
from pathlib import Path
import io, re, json, hashlib, zipfile, datetime as dt
import numpy as np
import pandas as pd
from fractions import Fraction as F
R=Path(__file__).resolve().parents[1]
CELLS=['F0001','F0002','F0003']
MD5={'F0001':'8181a89aee2e525ca52b60491ee3209b','F0002':'3eaca9d8eff8ccf2d1caf55575a3af76','F0003':'4f3c7587d697b2b863b36a13a7ab8f25'}
COLUMNS=['t_h','dataset','t_step','t_set','line','command','U_V','I_A','Ah','Ah_step','Wh','cycle','state']
def save(path,data): path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def ratio(x):return F.from_float(float(x))
def lower(x):
 v=float(x);return float(np.nextafter(v,-np.inf)) if F.from_float(v)>x else v
def upper(x):
 v=float(x);return float(np.nextafter(v,np.inf)) if F.from_float(v)<x else v

def parse_eis(raw,path):
 lines=raw.decode('utf-8-sig').splitlines();rows=[];header={}
 for line in lines:
  for k in ['File','System','Temperature','Potential','Current','Time']:
   m=re.match('^'+k+r'\.*:\s*(.*)$',line)
   if m:header[k]=m[1]
  parts=line.split()
  if len(parts)==8 and parts[0].isdigit():rows.append([float(x) for x in parts])
 a=np.array(rows);assert len(a)==91 and np.isfinite(a).all(),path
 # first/only maximum is the turning point shared by the rising and falling branches
 idx=int(a[:,1].argmax());d=a[idx:];assert len(d)==71 and np.all(np.diff(d[:,1])<0)
 e={'source_file':path,'source_sha256':hashlib.sha256(raw).hexdigest(),'header':header,'recorded_rows':len(a),'selected_rows':len(d),'branch_rule':'descending branch beginning at maximum frequency','frequency_Hz':d[::-1,1].tolist(),'Zreal_mOhm':(1000*d[::-1,2]).tolist(),'Zimag_mOhm':(1000*d[::-1,3]).tolist(),'EIS_recorded_DC_voltage_range_V':[float(min(d[:,6])),float(max(d[:,6]))]}
 return e

rows=[];manifest=[];amend=[]
for cell in CELLS:
 raw=(R/'raw'/f'{cell}.zip').read_bytes();assert hashlib.md5(raw).hexdigest()==MD5[cell]
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  assert z.testzip() is None
  manifest.append({'cell':cell,'bytes':len(raw),'md5':hashlib.md5(raw).hexdigest(),'sha256':hashlib.sha256(raw).hexdigest(),'producer_MD5_matches':True,'CRC_pass':True})
  for stage,checkup in [('BOL','M0001_checkup'),('EOL','M0003_checkup')]:
   eisname=f'{cell}/{cell.lower()}_{stage.lower()}'+('_2' if (cell,stage)==('F0002','BOL') else '')+'.txt'
   e=parse_eis(z.read(eisname),eisname)
   cr=z.read(f'{cell}/{checkup}.txt');text=cr.decode('utf-8-sig');sep=',' if cell!='F0003' else r'\s+'
   d=pd.read_csv(io.StringIO(text),comment='~',sep=sep,header=None,names=COLUMNS)
   assert len(d)>1000 and d[['t_h','U_V','I_A']].notna().all().all()
   assert np.isfinite(d[['t_h','U_V','I_A']].values).all()
   assert np.all(np.diff(d.t_h.values)>=0)
   match=re.search(r'~Start of Test:\s*([^\r\n]+)',text);teststart=dt.datetime.strptime(match[1],'%d.%m.%Y %H:%M:%S')
   # Use recorded measurement header, NOT ZIP modification/export dates. Timezone is unspecified.
   dm=re.search(r'([A-Z]{3}),([0-9]{2})[.]([0-9]{4})',e['header']['File'].upper())
   assert dm, eisname
   eisdate=dt.datetime.strptime(dm[2]+' '+dm[1]+' '+dm[3]+' '+e['header']['Time'].split(' - ')[0],'%d %b %Y %H:%M:%S')
   for direction,cmd,sign in [('charge','Charge',1),('discharge','Discharge',-1)]:
    g=d[d.command==cmd].groupby('line',sort=False).agg(t0=('t_h','min'),t1=('t_h','max'),imax=('I_A','max'),imin=('I_A','min'))
    durations=(g.t1-g.t0)*3600;eligible=g[durations.between(9.95,10.05)&(g[['imax','imin']].abs().max(axis=1)>0.06)]
    assert len(eligible)==1
    p=d[d.line==eligible.index[0]].copy();tt=(p.t_h.values-p.t_h.iloc[0])*3600;T=float(tt[-1]);assert len(p)>900 and np.all(np.diff(tt)>0)
    plateau=p[(tt>=2)&(tt<=8)];ip=abs(float(plateau.I_A.median()));assert ip>0
    assert sign*p.I_A.iloc[-1]>0.8*ip
    base=d.loc[:p.index[0]-1];base=base[base.I_A.abs()<.005].tail(10);assert len(base)==10
    v0=float(base.U_V.median());vT=float(p.U_V.iloc[-1]);meas=1000*sign*(vT-v0)/ip
    u=sign*p.I_A.values/ip
    # Every interval uses the actually recorded signed current, with a matched left endpoint at T.
    q={'case_id':cell+'_'+stage+'_'+direction,'cell':cell,'stage':stage,'direction':direction,'frequency_Hz':e['frequency_Hz'],'Zreal_mOhm':e['Zreal_mOhm'],'Zimag_mOhm':e['Zimag_mOhm'],'T_s':T,'a_s':tt[:-1].tolist(),'b_s':tt[1:].tolist(),'u':u[:-1].tolist(),'iT':float(u[-2]),'iT_right_recorded':float(u[-1]),'Rmeas_mOhm':float(meas),'normalizing_current_A':ip,'pulse_raw_sha256':hashlib.sha256(cr).hexdigest(),'eis_raw_sha256':e['source_sha256'],'source_EIS_file':eisname,'source_pulse_file':f'{cell}/{checkup}.txt','nominal_SOC_pct':50,'nominal_temperature_C':25,'temperature_basis':'producer-supplied correspondence and published nominal protocol, not sensor calibration','state_continuity':'assumed under producer-reported OCV chamber rest; not independently measured'}
    tail=p.tail(10);Imin=min(abs(plateau.I_A));Imax=max(abs(plateau.I_A))
    if sign==1:Dmin=ratio(tail.U_V.min())-ratio(base.U_V.max());Dmax=ratio(tail.U_V.max())-ratio(base.U_V.min())
    else:Dmin=ratio(base.U_V.min())-ratio(tail.U_V.max());Dmax=ratio(base.U_V.max())-ratio(tail.U_V.min())
    assert Dmin>0 and Imin>0
    q['recorded_envelope_lower_mOhm']=lower(1000*Dmin/ratio(Imax));q['recorded_envelope_upper_mOhm']=upper(1000*Dmax/ratio(Imin));assert q['recorded_envelope_lower_mOhm']<=meas<=q['recorded_envelope_upper_mOhm']
    save(R/'inputs'/(q['case_id']+'.json'),q)
    # Two timestamps for each EIS/pulse pairing, not charge-only delay reused for discharge.
    pulse_dt=teststart+dt.timedelta(hours=float(p.t_h.iloc[0]));delay=(eisdate-pulse_dt).total_seconds()/3600
    r={k:q[k] for k in ['case_id','cell','stage','direction','T_s','Rmeas_mOhm','normalizing_current_A','recorded_envelope_lower_mOhm','recorded_envelope_upper_mOhm','source_EIS_file']}
    r.update({'pulse_datetime':pulse_dt.isoformat(),'EIS_datetime':eisdate.isoformat(),'delay_hours':delay,'EIS_points':len(e['frequency_Hz']),'waveform_intervals':len(q['u']),'negative_current_intervals':int(sum(u[:-1]<0)),'initial_normalized_current':float(u[0]),'longest_delay_sensitivity_excluded':cell=='F0001' and stage=='EOL','voltage_baseline_V':v0,'voltage_endpoint_V':vT,'baseline_range_V':float(base.U_V.max()-base.U_V.min()),'endpoint_tail_range_V':float(tail.U_V.max()-tail.U_V.min()),'EIS_DC_min_V':e['EIS_recorded_DC_voltage_range_V'][0],'EIS_DC_max_V':e['EIS_recorded_DC_voltage_range_V'][1]})
    rows.append(r)
    oldp=R/'history/Project52_ExternalIntake_E1/processed'/(q['case_id']+'_pulse.json')
    old=json.loads(oldp.read_text())
    amend.append({'case_id':q['case_id'],'E1_kernel_T_s':10.0,'E1_voltage_sample_T_s':old['T_minus_endpoint_s'],'E3_matched_kernel_and_voltage_T_s':T,'measured_response_delta_mOhm':meas-old['endpoint_mV_per_A'],'E1_Q_s':old['Q_norm_s'],'E3_Q_s':float(np.sum(u[:-1]*np.diff(tt))),'EIS_reselected':cell=='F0002' and stage=='BOL'})
   save(R/'inputs'/(cell+'_'+stage+'_EIS.json'),e)
pd.DataFrame(rows).to_csv(R/'results/INPUT_CASES_12.csv',index=False)
pd.DataFrame(amend).to_csv(R/'qa/E1_E3_PREPROCESSING_LEDGER.csv',index=False)
save(R/'qa/SOURCE_INTEGRITY.json',manifest)
print(pd.DataFrame(rows)[['case_id','T_s','Rmeas_mOhm','recorded_envelope_upper_mOhm','delay_hours','negative_current_intervals']].to_string(index=False))
