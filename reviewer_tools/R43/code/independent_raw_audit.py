#!/usr/bin/env python3
"""Independent 75-decimal raw-window arithmetic and unscaled LP re-evaluation.
No calls to the production scatter, segmentation or LP-building functions.
"""
from pathlib import Path
from decimal import Decimal,localcontext
import csv,hashlib,io,json,math,zipfile,argparse
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]
def D(x):return Decimal.from_float(float(x))
def median(x):
 y=sorted(x);n=len(y)
 return y[n//2] if n%2 else (y[n//2-1]+y[n//2])/2
def mad(x):
 m=median(x);return Decimal('1.4826')*median([abs(v-m) for v in x])
def residual_scale(t,y):
 n=Decimal(len(t));xm=sum(t)/n;ym=sum(y)/n
 slope=sum((x-xm)*(v-ym) for x,v in zip(t,y))/sum((x-xm)**2 for x in t)
 return mad([v-ym-slope*(x-xm) for x,v in zip(t,y)])
def check(flag,name):
 if not flag:raise AssertionError(name)
 LOG.append(name)
LOG=[]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path);args=ap.parse_args()
 sc={x['trace_id']:x for x in csv.DictReader((ROOT/'evidence/R42_Raw_Scatter_66_Conditions.csv').open())}
 results=[];source={}
 if args.raw:
  z=zipfile.ZipFile(args.raw)
  check(hashlib.sha256(args.raw.read_bytes()).hexdigest()=='875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f','original archive SHA256')
  for temp in [5,15,25,35,45]:
   d=loadmat(io.BytesIO(z.read(f'{temp}deg/{temp}deg_Overview.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
   source[temp]={k:np.asarray(getattr(d,k)).ravel() for k in ['Zeit','Strom','Spannung','Schritt']}
 for p in sorted((ROOT/'evidence/source_windows').glob('*.json')):
  w=json.loads(p.read_text());row=sc[p.stem];temp=int(row['temperature_C'])
  if args.raw:
   a=source[temp];s,e=w['raw_start_index'],w['raw_stop_index'];ids=w['baseline_indices']
   for key,src in [('pulse_I_A','Strom'),('pulse_V','Spannung')]:check(np.array_equal(w[key],a[src][s:e]),p.stem+' original '+key)
   check(np.array_equal(w['pulse_t_s'],a['Zeit'][s:e]-a['Zeit'][s]),p.stem+' original relative timestamps')
   check(np.array_equal(w['baseline_V'],a['Spannung'][ids]),p.stem+' original baseline')
  with localcontext() as ctx:
   ctx.prec=75
   t=list(map(D,w['pulse_t_s']));v=list(map(D,w['pulse_V']));i=list(map(D,w['pulse_I_A']));b=list(map(D,w['baseline_V']))
   # Independently canonicalize by a map: last record at each exact time wins.
   bytime={x:(y,cur) for x,y,cur in zip(t,v,i)};tt=sorted(bytime);T=Decimal(20)
   left=max(x for x in tt if x<=T);right=min(x for x in tt if x>=T)
   vt=bytime[left][0] if left==right else bytime[left][0]+(bytime[right][0]-bytime[left][0])*(T-left)/(right-left)
   v0=median(b);dv=abs(vt-v0)
   aidx=[k for k,x in enumerate(t) if Decimal(5)<=x<=Decimal(15)];eidx=[k for k,x in enumerate(t) if Decimal(19)<=x<=Decimal('20.05')]
   ip=abs(median([i[k] for k in aidx]));s0=mad(b);sv=residual_scale([t[k] for k in eidx],[v[k] for k in eidx]);si=residual_scale([t[k] for k in eidx],[i[k] for k in eidx]);sip=residual_scale([t[k] for k in aidx],[i[k] for k in aidx])
   h=Decimal(3000)*((s0*s0+sv*sv)/(ip*ip)+(dv*si/(ip*ip))**2).sqrt()
   hp=Decimal(3000)*((s0*s0+sv*sv)/(ip*ip)+(dv*sip/(ip*ip))**2).sqrt()
   vals={'sigma_baseline_V':s0,'sigma_endpoint_V':sv,'sigma_endpoint_I_A':si,'sigma_plateau_I_A':sip,'halfwidth_3sigma_endpoint_mOhm':h,'halfwidth_3sigma_plateau_mOhm':hp,'Rmeas_mOhm':1000*dv/ip,'I_p_A':ip}
   errors={k:abs(float(x)-float(row[k])) for k,x in vals.items()}
   for k,err in errors.items():check(err<5e-10,p.stem+' decimal '+k)
   results.append({'trace_id':p.stem,'h_endpoint_75dp':str(h),'h_plateau_75dp':str(hp),'maximum_difference':max(errors.values()),'h_endpoint_upper_float':math.nextafter(float(h),math.inf),'h_plateau_upper_float':math.nextafter(float(hp),math.inf)})
 check(len(results)==66,'exactly 66 independent source-window audits')
 # Re-evaluate stored finite-grid LP solutions independently, without the solver.
 sols=json.loads((ROOT/'evidence/R42_Temperature_Grid_Primal_120.json').read_text());pe=[]
 for name,c in sols.items():
  q=json.loads((ROOT/'evidence/fresh_inputs'/(c['trace_id']+'.json')).read_text());idx=c['nonzero_indices'];val=c['nonzero_values'];n=c['grid_n'];tau=np.logspace(-10,10,n)
  R,L,C=[next((v for j,v in zip(idx,val) if j==k),0.) for k in [0,1,2]];atoms=[(tau[j-3],v) for j,v in zip(idx,val) if j>=3]
  w=2*np.pi*np.asarray(q['frequency_Hz']);pred=R+1j*w*L+C/(1j*w)
  for t0,mass in atoms:pred+=mass/(1+1j*w*t0)
  excess=max(np.max(abs(pred.real-q['Zreal_mOhm'])),np.max(abs(pred.imag-q['Zimag_mOhm'])))-q['eps_mOhm']
  obj=R*q['iT']+C*sum(u*(b-a) for a,b,u in zip(q['a_s'],q['b_s'],q['u']))
  for t0,mass in atoms:
   obj+=mass*sum(u*math.exp(-(q['T_s']-b)/t0)*(-math.expm1(-(b-a)/t0)) for a,b,u in zip(q['a_s'],q['b_s'],q['u']))
  check(min(val)>=0,name+' positive parameters');check(excess<2e-7,name+' original unscaled EIS constraints');check(abs(obj-c['objective_mOhm'])<2e-6,name+' independent stable waveform objective')
  pe.append({'case':name,'excess_mOhm':float(excess),'objective_difference_mOhm':float(abs(obj-c['objective_mOhm']))})
 check(len(pe)==120,'120 finite-grid solutions independently evaluated')
 (ROOT/'qa/R42_Independent_75DP_66.json').write_text(json.dumps(results,indent=2)+'\n')
 (ROOT/'qa/R42_Temperature_Independent_120.json').write_text(json.dumps(pe,indent=2)+'\n')
 summary={'status':'PASS','checks':len(LOG),'raw_source_reopened':bool(args.raw),'source_conditions':66,'decimal_precision':75,'maximum_scatter_or_response_difference':max(x['maximum_difference'] for x in results),'grid_solutions':120,'maximum_grid_primal_excess_mOhm':max(x['excess_mOhm'] for x in pe),'maximum_grid_objective_difference_mOhm':max(x['objective_difference_mOhm'] for x in pe),'scope':'Arithmetic and raw-window provenance checks, not measurement calibration or an independent human audit.'}
 (ROOT/'qa/R42_Independent_Arithmetic_Summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
