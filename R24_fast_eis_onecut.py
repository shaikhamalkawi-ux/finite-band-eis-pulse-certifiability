import sys, os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np, pandas as pd, scipy.io as sio, zipfile, io
from scipy.optimize import linprog
from scipy import sparse as sp
req=float(sys.argv[1]); tag=sys.argv[2] if len(sys.argv)>2 else str(req)
ZP='MatFiles.zip'; TAB='R4D_RWTH_5C_ContinuumSafe_42cases.csv'; T=20.; NTAU=801; MARGIN=0.0002
TAU=np.logspace(-10,10,NTAU); tab=pd.read_csv(TAB)
with zipfile.ZipFile(ZP) as z:
 ov=sio.loadmat(io.BytesIO(z.read('5deg/5deg_Overview.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
t=np.asarray(ov.Zeit).squeeze();I=np.asarray(ov.Strom).squeeze();V=np.asarray(ov.Spannung).squeeze();steps=np.asarray(ov.Schritt).squeeze()
on=np.abs(I)>0.5; d=np.diff(on.astype(int)); starts=np.where(d==1)[0]+1; ends=np.where(d==-1)[0]+1
if on[0]: starts=np.r_[0,starts]
if on[-1]: ends=np.r_[ends,len(on)]
segs=[]
for s,e in zip(starts,ends):
 dur=t[e-1]-t[s]; med=np.median(I[s:e])
 if 15<=dur<=25 and abs(med)>1.5: segs.append((s,e,int(steps[s])))
groups=[];cur=[];prev=None
for sg in segs:
 if prev is not None and sg[0]-prev>10000: groups.append(cur);cur=[]
 cur.append(sg);prev=sg[0]
groups.append(cur)
SOCG={90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10};STEP={('charge',0.7):3,('discharge',0.7):8,('charge',1.0):13,('discharge',1.0):18,('charge',1.4):23,('discharge',1.4):28}
EIS={}
with zipfile.ZipFile(ZP) as z:
 for soc in sorted(tab.SOC.unique()):
  e=sio.loadmat(io.BytesIO(z.read(f'5deg/5deg_SOC{int(soc)}.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
  f=np.asarray(e.ActFreq).squeeze();zr=np.asarray(e.Zreal1).squeeze();zi=np.asarray(e.Zimg1).squeeze();mask=(f>0)&((zr!=0)|(zi!=0));uf=np.unique(f[mask]);F=[];ZR=[];ZI=[]
  for ff in uf:
   j=np.where(mask&(f==ff))[0][0];F.append(ff);ZR.append(zr[j]);ZI.append(zi[j])
  EIS[int(soc)]=(np.array(F),np.array(ZR),np.array(ZI))

def pulse_f(row):
 soc=int(row.SOC);st=STEP[(row.direction,float(row.C_rate))];g=groups[SOCG[soc]];s,e=[x[:2] for x in g if x[2]==st][0]
 tt=t[s:e]-t[s];ii=I[s:e];vv=V[s:e];pre=np.arange(max(0,s-200),s);zz=pre[np.abs(I[pre])<=0.5];V0=float(np.median(V[zz[-10:]]));VT=float(np.interp(T,tt,vv));Ip_abs=abs(VT-V0)/(float(row.Rmeas_mOhm)/1000.);sgn=np.sign(np.median(ii[(tt>=5)&(tt<=15)]));norm=ii/(sgn*Ip_abs)
 aa=[];bb=[];vals=[]
 for j in range(len(tt)-1):
  a=float(tt[j]);b=float(tt[j+1]);
  if a>=T:break
  b=min(b,T)
  if b<=a:continue
  aa.append(a);bb.append(b);vals.append(float(norm[j]))
  if b>=T:break
 aa=np.array(aa);bb=np.array(bb);vals=np.array(vals);jt=np.searchsorted(tt,T,side='right')-1;iT=float(norm[max(0,min(jt,len(norm)-1))]);Q=float(np.sum(vals*(bb-aa)))
 kw=np.zeros(NTAU)
 for a,bb_,v in zip(aa,bb,vals):kw += v*(np.exp(-(T-bb_)/TAU)-np.exp(-(T-a)/TAU))
 f=np.zeros(3+NTAU);f[0]=iT;f[2]=Q;f[3:]=kw
 return f

def model_matrix(F):
 om=2*np.pi*F;x=om[:,None]*TAU[None,:];Ar=np.zeros((len(F),3+NTAU));Ar[:,0]=1;Ar[:,3:]=1/(1+x*x);Ai=np.zeros_like(Ar);Ai[:,1]=om;Ai[:,2]=-1/om;Ai[:,3:]=-x/(1+x*x);return np.vstack([Ar,Ai])

def min_eps(M,b):
 n=M.shape[1];A1=sp.hstack([sp.csr_matrix(M),-np.ones((M.shape[0],1))],format='csr');A2=sp.hstack([sp.csr_matrix(-M),-np.ones((M.shape[0],1))],format='csr');A=sp.vstack([A1,A2],format='csr');r=linprog(np.r_[np.zeros(n),1.],A_ub=A,b_ub=np.r_[b,-b],bounds=[(0,None)]*n+[(0,None)],method='highs',options={'presolve':True});
 if not r.success: raise RuntimeError(r.message)
 return float(r.fun)

def min_response(M,b,eps,f):
 A=sp.vstack([sp.csr_matrix(M),sp.csr_matrix(-M)],format='csr');r=linprog(f,A_ub=A,b_ub=np.r_[b+eps,-(b-eps)],bounds=[(0,None)]*M.shape[1],method='highs',options={'presolve':True});
 if not r.success: raise RuntimeError(r.message)
 return float(r.fun)

fobjs={idx:pulse_f(row) for idx,row in tab.iterrows()}
rows=[]
for soc in sorted(tab.SOC.unique()):
 F,ZR,ZI=EIS[int(soc)];mask=F>=req*0.999;Fc=F[mask];M=model_matrix(Fc);b=np.r_[ZR[mask],ZI[mask]];eps0=min_eps(M,b);eps=eps0+MARGIN
 for idx,row in tab[tab.SOC==soc].iterrows():
  lower=min_response(M,b,eps,fobjs[idx]);meas=float(row.Rmeas_mOhm);rows.append({'requested_cutoff_Hz':req,'actual_fmin_Hz':float(Fc.min()),'n_freq':len(Fc),'SOC':int(soc),'direction':row.direction,'C_rate':float(row.C_rate),'grid_min_eps_mOhm':eps0,'analysis_eps_mOhm':eps,'Rmeas_mOhm':meas,'grid_lower_mOhm':lower,'grid_gap_mOhm':max(0.,lower-meas),'excluded_grid':lower>meas+1e-9,'locked10_lower_mOhm':float(row.lower_mOhm),'locked10_gap_mOhm':float(row.gap_mOhm)})
df=pd.DataFrame(rows);df.to_csv(f'R24_FastEIS_{tag}_Cases.csv',index=False)
print(df.groupby(['direction','C_rate']).agg(n=('excluded_grid','size'),excluded=('excluded_grid','sum'),median_gap=('grid_gap_mOhm','median')).to_string())
print('SUMMARY', {'req_mHz':req*1000,'actual_mHz':df.actual_fmin_Hz.iloc[0]*1000,'n_freq':int(df.n_freq.iloc[0]),'excluded':int(df.excluded_grid.sum()),'median_gap':float(df.grid_gap_mOhm.median()),'min_gap':float(df.grid_gap_mOhm.min()),'max_gap':float(df.grid_gap_mOhm.max()),'median_eps':float(df.analysis_eps_mOhm.median())})
if req<=0.0101:
 print('REPLAY maxabs',float(np.max(np.abs(df.grid_lower_mOhm-df.locked10_lower_mOhm))),'medianabs',float(np.median(np.abs(df.grid_lower_mOhm-df.locked10_lower_mOhm))))
sys.stdout.flush(); os._exit(0)
