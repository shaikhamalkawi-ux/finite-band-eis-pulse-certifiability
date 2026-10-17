import os, math, io, zipfile, json
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np, pandas as pd, scipy.io as sio
from scipy.optimize import linprog
from scipy import sparse as sp

ZP='MatFiles.zip'; TAB='R4D_RWTH_5C_ContinuumSafe_42cases.csv'; tab=pd.read_csv(TAB); T=20.; MARGIN=0.0002
with zipfile.ZipFile(ZP) as z:
    ov=sio.loadmat(io.BytesIO(z.read('5deg/5deg_Overview.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
t=np.asarray(ov.Zeit).squeeze();I=np.asarray(ov.Strom).squeeze();V=np.asarray(ov.Spannung).squeeze();steps=np.asarray(ov.Schritt).squeeze()
on=np.abs(I)>0.5;d=np.diff(on.astype(int));starts=np.where(d==1)[0]+1;ends=np.where(d==-1)[0]+1
if on[0]:starts=np.r_[0,starts]
if on[-1]:ends=np.r_[ends,len(on)]
segs=[]
for s,e in zip(starts,ends):
 dur=t[e-1]-t[s];med=np.median(I[s:e])
 if 15<=dur<=25 and abs(med)>1.5:segs.append((s,e,int(steps[s])))
groups=[];cur=[];prev=None
for sg in segs:
 if prev is not None and sg[0]-prev>10000:groups.append(cur);cur=[]
 cur.append(sg);prev=sg[0]
groups.append(cur)
SOCG={90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10};STEP={('charge',0.7):3,('discharge',0.7):8,('charge',1.0):13,('discharge',1.0):18,('charge',1.4):23,('discharge',1.4):28}

def load_eis(soc,req):
 with zipfile.ZipFile(ZP) as z:
  e=sio.loadmat(io.BytesIO(z.read(f'5deg/5deg_SOC{soc}.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
 f=np.asarray(e.ActFreq).squeeze();zr=np.asarray(e.Zreal1).squeeze();zi=np.asarray(e.Zimg1).squeeze();mask=(f>0)&((zr!=0)|(zi!=0));uf=np.unique(f[mask]);F=[];ZR=[];ZI=[]
 for ff in uf:
  j=np.where(mask&(f==ff))[0][0];F.append(ff);ZR.append(zr[j]);ZI.append(zi[j])
 F=np.array(F);ZR=np.array(ZR);ZI=np.array(ZI);keep=F>=req*0.999
 return F[keep],ZR[keep],ZI[keep]

def pulse(soc,direction,rate,row):
 st=STEP[(direction,rate)];g=groups[SOCG[soc]];s,e=[x[:2] for x in g if x[2]==st][0]
 tt=t[s:e]-t[s];ii=I[s:e];vvraw=V[s:e];pre=np.arange(max(0,s-200),s);zz=pre[np.abs(I[pre])<=0.5];V0=float(np.median(V[zz[-10:]]));VT=float(np.interp(T,tt,vvraw));Ip_abs=abs(VT-V0)/(float(row.Rmeas_mOhm)/1000.);sgn=np.sign(np.median(ii[(tt>=5)&(tt<=15)]));norm=ii/(sgn*Ip_abs)
 aa=[];bb=[];vv=[]
 for j in range(len(tt)-1):
  a=float(tt[j]);b=float(tt[j+1]);
  if a>=T:break
  b=min(b,T)
  if b<=a:continue
  aa.append(a);bb.append(b);vv.append(float(norm[j]))
  if b>=T:break
 jt=np.searchsorted(tt,T,side='right')-1;iT=float(norm[max(0,min(jt,len(norm)-1))])
 return np.array(aa),np.array(bb),np.array(vv),iT

def min_eps_grid(F,ZR,ZI,nt=801):
 tau=np.logspace(-10,10,nt);om=2*np.pi*F;x=om[:,None]*tau[None,:];Ar=np.zeros((len(F),3+nt));Ar[:,0]=1;Ar[:,3:]=1/(1+x*x);Ai=np.zeros_like(Ar);Ai[:,1]=om;Ai[:,2]=-1/om;Ai[:,3:]=-x/(1+x*x);M=np.vstack([Ar,Ai]);b=np.r_[ZR,ZI];n=M.shape[1];A1=sp.hstack([sp.csr_matrix(M),-np.ones((M.shape[0],1))],format='csr');A2=sp.hstack([sp.csr_matrix(-M),-np.ones((M.shape[0],1))],format='csr');r=linprog(np.r_[np.zeros(n),1.],A_ub=sp.vstack([A1,A2]),b_ub=np.r_[b,-b],bounds=[(0,None)]*n+[(0,None)],method='highs');return float(r.fun)

def check(req,soc,direction,rate,margin=0.03):
 row=tab[(tab.SOC==soc)&(tab.direction==direction)&(np.isclose(tab.C_rate,rate))].iloc[0]; measured=float(row.Rmeas_mOhm)
 F,ZR,ZI=load_eis(soc,req);eps=min_eps_grid(F,ZR,ZI)+MARGIN;om=2*np.pi*F;aa,bb,vv,iT=pulse(soc,direction,rate,row);area=float(np.sum(vv*(bb-aa)))
 def kw(tau):
  tau=np.asarray(tau,float);out=np.zeros_like(tau)
  for a,b,v in zip(aa,bb,vv):out += v*(np.exp(-(T-b)/tau)-np.exp(-(T-a)/tau))
  return out
 bdat=np.r_[ZR,ZI];lower=bdat-eps;upper=bdat+eps;cdual=np.r_[upper,-lower];m=len(bdat)
 def buildM(active):
  x=om[:,None]*active[None,:];Ar=np.zeros((len(F),3+len(active)));Ar[:,0]=1;Ar[:,3:]=1/(1+x*x);Ai=np.zeros_like(Ar);Ai[:,1]=om;Ai[:,2]=-1/om;Ai[:,3:]=-x/(1+x*x);return np.vstack([Ar,Ai])
 def fobj(active):
  f=np.zeros(3+len(active));f[0]=iT;f[2]=area;f[3:]=kw(active);return f
 active=np.logspace(-10,10,801);target=-fobj(active);target[3:]+=margin/(1+active);M=buildM(active)
 r=linprog(cdual,A_ub=np.hstack([-M.T,M.T]),b_ub=-target,bounds=[(0,None)]*(2*m),method='highs')
 if not r.success:return {'success':False,'message':r.message}
 y=r.x[:m]-r.x[m:]
 def resid(tau):
  tau=np.asarray(tau,float);x=om[:,None]*tau[None,:];meas=np.vstack([1/(1+x*x),-x/(1+x*x)]);return y@meas+kw(tau)
 xs=np.linspace(-30,30,30000);rv=resid(np.exp(xs));j=int(np.argmin(rv));tv=float(abs(vv[0])+np.sum(np.abs(np.diff(vv))));Mlog=0.5*float(np.sum(np.abs(y)))+tv/math.e
 # compact coarse proof-like h=1e-4 to keep runtime manageable
 h=5e-6;minv=1e99
 for start in np.arange(-20,6,h*30000):
  stop=min(6,start+h*30000);xx=np.arange(start,stop+h/2,h);rr=resid(np.exp(xx));minv=min(minv,float(np.min(rr)))
 compact_lb=minv-Mlog*h/2-1e-8
 yr=y[:len(om)];yi=y[len(om):];tau0=math.exp(-20);c0=float(np.sum(yr)+iT);small_tail_lb=c0-float(np.sum(np.abs(yr)*(om*tau0)**2))-float(np.sum(np.abs(yi)*om*tau0))
 # large-tail interval check as in R9 validator
 import mpmath as mp
 mp.iv.dps=40; zmax=math.exp(-6)
 def Hmeas_iv(zl,zh):
  z=mp.iv.mpf([zl,zh]);tot=mp.iv.mpf([0,0])
  for cc,w in zip(yr,om):tot += cc*z/(w*w+z*z)
  for cc,w in zip(yi,om):tot += cc*(-w/(w*w+z*z))
  return float(tot.a)
 def J(z):
  if z==0:return area
  return sum(v*(math.exp(-(T-b)*z)-math.exp(-(T-a)*z))/z for a,b,v in zip(aa,bb,vv))
 large_lb=1e99
 for q in range(200):
  zl=zmax*q/200;zh=zmax*(q+1)/200;large_lb=min(large_lb,Hmeas_iv(zl,zh)+J(zh))
 finite=(buildM(active).T@y-(-fobj(active)))[:3]
 passes=bool(compact_lb>0 and small_tail_lb>0 and large_lb>0 and np.min(finite)>=-1e-8 and -r.fun>measured)
 return {'success':True,'requested_cutoff_mHz':req*1000,'actual_fmin_mHz':float(F.min()*1000),'n_freq':len(F),'SOC':soc,'direction':direction,'C_rate':rate,'analysis_eps_mOhm':eps,'measured_mOhm':measured,'dual_certified_lower_mOhm':float(-r.fun),'dual_gap_mOhm':float(-r.fun-measured),'numeric_residual_min':float(rv[j]),'compact_lower_bound':float(compact_lb),'small_tail_lower_bound':float(small_tail_lb),'large_tail_H_lower_bound':float(large_lb),'finite_slack_min':float(np.min(finite)),'passes_validation':passes}

cases=[(0.03,70,'charge',1.4,0.03)]
out=[]
for c in cases:
 res=check(*c);out.append(res);print(json.dumps(res),flush=True)
pd.DataFrame(out).to_csv('R25_33mHz_RefinedProof.csv',index=False)
os._exit(0)
