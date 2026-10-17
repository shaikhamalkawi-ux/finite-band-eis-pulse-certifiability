import numpy as np, pandas as pd, scipy.io as sio, zipfile, io, scipy.optimize as opt, scipy.sparse as sp, time
ZP='MatFiles.zip'; TAB='R4D_RWTH_5C_ContinuumSafe_42cases.csv'; T=20.;NTAU=1201
TAU=np.logspace(-10,10,NTAU);tab=pd.read_csv(TAB)
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
EIS={}
with zipfile.ZipFile(ZP) as z:
 for soc in sorted(tab.SOC.unique()):
  e=sio.loadmat(io.BytesIO(z.read(f'5deg/5deg_SOC{int(soc)}.mat')),squeeze_me=False,struct_as_record=False)['diga'][0,0].daten[0,0]
  f=np.asarray(e.ActFreq).squeeze();zr=np.asarray(e.Zreal1).squeeze();zi=np.asarray(e.Zimg1).squeeze();mask=(f>0)&((zr!=0)|(zi!=0));uf=np.unique(f[mask]);F=[];ZR=[];ZI=[]
  for ff in uf:
   j=np.where(mask&(f==ff))[0][0];F.append(ff);ZR.append(zr[j]);ZI.append(zi[j])
  EIS[int(soc)]=(np.array(F),np.array(ZR),np.array(ZI))

def case_mats(row):
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
 F,ZR,ZI=EIS[soc];om=2*np.pi*F;x=om[:,None]*TAU[None,:];kr=1/(1+x*x);ki=-x/(1+x*x)
 Ar=np.zeros((len(F),3+NTAU));Ar[:,0]=1;Ar[:,3:]=kr;Ai=np.zeros_like(Ar);Ai[:,1]=om;Ai[:,2]=-1/om;Ai[:,3:]=ki;M=np.vstack([Ar,Ai]);b=np.r_[ZR,ZI]
 kw=np.zeros(NTAU)
 for a,bb_,v in zip(aa,bb,vals):kw += v*(np.exp(-(T-bb_)/TAU)-np.exp(-(T-a)/TAU))
 f=np.zeros(3+NTAU);f[0]=iT;f[1]=0.;f[2]=Q;f[3:]=kw
 return M,b,f,Ip_abs,iT,Q

out=[];t0=time.time()
for idx,row in tab.iterrows():
 M,b,f,Ip,iT,Q=case_mats(row);eps=float(row.eps_mOhm);meas=float(row.Rmeas_mOhm);n=M.shape[1]
 # alpha-star: vars [x(n), alpha]. Mx - eps alpha <= b; -Mx - eps alpha <= -b; f'x <= meas; -alpha <= -1
 A1=sp.hstack([sp.csr_matrix(M), -eps*np.ones((M.shape[0],1))],format='csr');A2=sp.hstack([sp.csr_matrix(-M), -eps*np.ones((M.shape[0],1))],format='csr');A3=sp.csr_matrix(np.r_[f,0.][None,:]);A4=sp.csr_matrix(np.r_[np.zeros(n),-1.][None,:]);A=sp.vstack([A1,A2,A3,A4],format='csr');ub=np.r_[b,-b,meas,-1.]
 c=np.r_[np.zeros(n),1.];bounds=[(0,None)]*n+[(1,None)]
 r=opt.linprog(c,A_ub=A,b_ub=ub,bounds=bounds,method='highs')
 if not r.success: raise RuntimeError((idx,r.message))
 alpha=float(r.x[-1])
 # factor1 replay lower
 Afix=sp.vstack([sp.csr_matrix(M),sp.csr_matrix(-M)],format='csr');ubfix=np.r_[b+eps,-(b-eps)];rr=opt.linprog(f,A_ub=Afix,b_ub=ubfix,bounds=[(0,None)]*n,method='highs');lb=float(rr.fun) if rr.success else np.nan
 out.append({'SOC':int(row.SOC),'direction':row.direction,'C_rate':float(row.C_rate),'eps0_mOhm':eps,'Rmeas_mOhm':meas,'locked_lower_mOhm':float(row.lower_mOhm),'grid1201_lower_mOhm':lb,'replay_diff_mOhm':lb-float(row.lower_mOhm),'alpha_star_grid1201':alpha,'eps_star_mOhm':alpha*eps,'Ip_recovered_A':Ip,'iT_norm':iT,'Q_T_s':Q})
 print(f"{idx+1:02d}/42 SOC{int(row.SOC)} {row.direction} {row.C_rate}: a*={alpha:.4f}, replay={lb-float(row.lower_mOhm):+.4f}",flush=True)
df=pd.DataFrame(out);df.to_csv('R20_RWTH_5C_CriticalToleranceMultiplier_42cases.csv',index=False)
# survival
fac=[1,1.1,1.25,1.5,2,3,5,8,10,15,20]
s=[]
for a in fac:
 n=int((df.alpha_star_grid1201>a+1e-9).sum());s.append({'factor':a,'cases_still_excluded':n,'pct_still_excluded':100*n/42})
pd.DataFrame(s).to_csv('R20_RWTH_5C_CriticalToleranceMultiplier_Survival.csv',index=False)
print('\nSURVIVAL\n',pd.DataFrame(s).to_string(index=False));print('\nALPHA\n',df.alpha_star_grid1201.describe().to_string());print('\nREPLAY maxabs',df.replay_diff_mOhm.abs().max(),'medianabs',df.replay_diff_mOhm.abs().median());print('elapsed',time.time()-t0)
