import numpy as np, scipy.io as sio, zipfile, io, scipy.optimize as opt, math, pandas as pd, json, os
ZP='MatFiles.zip'
TAB='R4D_RWTH_5C_ContinuumSafe_42cases.csv'
assert os.path.exists(ZP) and os.path.exists(TAB)
tab=pd.read_csv(TAB)
with zipfile.ZipFile(ZP) as z:
    ov=sio.loadmat(io.BytesIO(z.read('5deg/5deg_Overview.mat')), squeeze_me=False, struct_as_record=False)['diga'][0,0].daten[0,0]
t=np.asarray(ov.Zeit).squeeze(); I=np.asarray(ov.Strom).squeeze(); V=np.asarray(ov.Spannung).squeeze(); steps=np.asarray(ov.Schritt).squeeze()
on=np.abs(I)>0.5; d=np.diff(on.astype(int)); starts=np.where(d==1)[0]+1; ends=np.where(d==-1)[0]+1
if on[0]: starts=np.r_[0,starts]
if on[-1]: ends=np.r_[ends,len(on)]
segs=[]
for s,e in zip(starts,ends):
    dur=t[e-1]-t[s]; med=np.median(I[s:e])
    if 15<=dur<=25 and abs(med)>1.5:
        segs.append((s,e,dur,med,float(np.median(V[max(0,s-20):s])),int(steps[s])))
groups=[]; cur=[]; prev=None
for sg in segs:
    if prev is not None and sg[0]-prev>10000:
        groups.append(cur); cur=[]
    cur.append(sg); prev=sg[0]
groups.append(cur)
soc_to_group={90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10}
step_map={('charge',0.7):3,('discharge',0.7):8,('charge',1.0):13,('discharge',1.0):18,('charge',1.4):23,('discharge',1.4):28}

def load_eis(soc):
    with zipfile.ZipFile(ZP) as z:
        e=sio.loadmat(io.BytesIO(z.read(f'5deg/5deg_SOC{soc}.mat')), squeeze_me=False, struct_as_record=False)['diga'][0,0].daten[0,0]
    f=np.asarray(e.ActFreq).squeeze(); zr=np.asarray(e.Zreal1).squeeze(); zi=np.asarray(e.Zimg1).squeeze()
    mask=(f>0)&((zr!=0)|(zi!=0)); uf=np.unique(f[mask])
    F=[];ZR=[];ZI=[]
    for ff in uf:
        ii=np.where(mask&(f==ff))[0][0]; F.append(ff);ZR.append(zr[ii]);ZI.append(zi[ii])
    return np.array(F),np.array(ZR),np.array(ZI)

def pulse(soc,direction,rate,T=20.0):
    st=step_map[(direction,rate)]; g=groups[soc_to_group[soc]]
    s,e=[x[:2] for x in g if x[5]==st][0]
    tt=t[s:e]-t[s]; ii=I[s:e]; Iref=np.median(ii[-10:]); norm=ii/Iref
    aa=[];bb=[];vv=[]
    for j in range(len(tt)-1):
        a=float(tt[j]);b=float(tt[j+1])
        if a>=T: break
        b=min(b,T)
        if b<=a: continue
        aa.append(a);bb.append(b);vv.append(float(norm[j]))
        if b>=T: break
    return np.array(tt),np.array(norm),np.array(aa),np.array(bb),np.array(vv),float(np.interp(T,tt,norm))

def solve_case(soc,direction,rate,margin,T=20.0):
    row=tab[(tab.SOC==soc)&(tab.direction==direction)&(np.isclose(tab.C_rate,rate))].iloc[0]
    eps=float(row.eps_mOhm); measured=float(row.Rmeas_mOhm)
    F,ZR,ZI=load_eis(soc); om=2*np.pi*F
    tt,norm,aa,bb,vv,iT=pulse(soc,direction,rate,T)
    area=float(np.sum(vv*(bb-aa)))
    def kw(tau):
        tau=np.asarray(tau,float); out=np.zeros_like(tau)
        for a,b,v in zip(aa,bb,vv): out += v*(np.exp(-(T-b)/tau)-np.exp(-(T-a)/tau))
        return out
    bdat=np.r_[ZR,ZI]; lower=bdat-eps; upper=bdat+eps; cdual=np.r_[upper,-lower]; m=len(bdat)
    def buildM(active):
        x=om[:,None]*active[None,:]; kr=1/(1+x*x);ki=-x/(1+x*x)
        Ar=np.zeros((len(F),3+len(active)));Ar[:,0]=1;Ar[:,3:]=kr
        Ai=np.zeros_like(Ar);Ai[:,1]=om;Ai[:,2]=-1/om;Ai[:,3:]=ki
        return np.vstack([Ar,Ai])
    def fobj(active):
        f=np.zeros(3+len(active));f[0]=iT;f[2]=area;f[3:]=kw(active);return f
    active=np.logspace(-10,10,801)
    target=-fobj(active); target[3:]+=margin/(1+active)
    M=buildM(active)
    r=opt.linprog(cdual,A_ub=np.hstack([-M.T,M.T]),b_ub=-target,bounds=[(0,None)]*(2*m),method='highs')
    if not r.success: raise RuntimeError(r.message)
    y=r.x[:m]-r.x[m:]
    def resid(tau):
        tau=np.asarray(tau,float); x=om[:,None]*tau[None,:]
        meas=np.vstack([1/(1+x*x),-x/(1+x*x)])
        return y@meas + kw(tau) # upper certificate for -F
    xs=np.linspace(-30,30,30000); rv=resid(np.exp(xs)); j=int(np.argmin(rv))
    # analytic global log-Lipschitz bound: each RC real/imag log-derivative <= 1/2; waveform <= TV/e
    tv=float(abs(vv[0])+np.sum(np.abs(np.diff(vv))))
    Mlog=0.5*float(np.sum(np.abs(y)))+tv/math.e
    # compact proof-like mesh over [-20,6]
    h=5e-5; minv=1e99; minx=None
    for start in np.arange(-20,6,h*50000):
        stop=min(6,start+h*50000)
        xx=np.arange(start,stop+h/2,h); rr=resid(np.exp(xx)); jj=int(np.argmin(rr))
        if rr[jj]<minv: minv=float(rr[jj]);minx=float(xx[jj])
    compact_lb=minv-Mlog*h/2-1e-8
    # small tail x<-20: limit at tau=0 and crude deviations
    yr=y[:len(om)]; yi=y[len(om):]; tau0=math.exp(-20)
    c0=float(np.sum(yr)+iT)
    Breal=float(np.sum(np.abs(yr)*(om*tau0)**2))
    Bimag=float(np.sum(np.abs(yi)*om*tau0))
    # waveform final-current convergence; safe zero here due exponential underflow at tau0
    small_tail_lb=c0-Breal-Bimag
    # large tail x>=6: factor residual = z H(z), z=1/tau. Validate H over 200 bins using interval arithmetic for rational measurement and monotonic positive waveform integral endpoint.
    import mpmath as mp
    mp.iv.dps=40; zmax=math.exp(-6)
    def Hmeas_iv(zl,zh):
        z=mp.iv.mpf([zl,zh]); tot=mp.iv.mpf([0,0])
        for c,w in zip(yr,om): tot += c*z/(w*w+z*z)
        for c,w in zip(yi,om): tot += c*(-w/(w*w+z*z))
        return float(tot.a)
    def J(z):
        if z==0: return area
        return sum(v*(math.exp(-(T-b)*z)-math.exp(-(T-a)*z))/z for a,b,v in zip(aa,bb,vv))
    large_lb=1e99
    for q in range(200):
        zl=zmax*q/200; zh=zmax*(q+1)/200
        lb=Hmeas_iv(zl,zh)+J(zh) # J decreases since normalized pulse current is positive
        large_lb=min(large_lb,float(lb))
    # finite-coordinate dual margins for R0,L,c
    finite_slack=(buildM(active).T@y - (-fobj(active)))[:3]
    return {
        'SOC':soc,'direction':direction,'C_rate':rate,'margin':margin,'eps_mOhm':eps,'measured_mOhm':measured,
        'certified_lower_mOhm':float(-r.fun),'certified_gap_mOhm':float(-r.fun-measured),
        'global_numeric_residual_min':float(rv[j]),'global_numeric_tau_at_min_s':float(math.exp(xs[j])),
        'compact_lower_bound':float(compact_lb),'small_tail_lower_bound':float(small_tail_lb),'large_tail_H_lower_bound':float(large_lb),
        'finite_R0_L_c_slack_min':float(np.min(finite_slack)),'Mlog':float(Mlog),'mesh_h':h,
        'passes_validation':bool(compact_lb>0 and small_tail_lb>0 and large_lb>0 and np.min(finite_slack)>=-1e-8 and -r.fun>measured)
    }

cases=[(30,'discharge',0.7,[0.02]),(50,'charge',1.0,[0.02]),(70,'charge',1.4,[0.1])]
out=[]
for soc,direction,rate,margins in cases:
    chosen=None
    for margin in margins:
        res=solve_case(soc,direction,rate,margin)
        print(json.dumps(res))
        if res['passes_validation']:
            chosen=res
            break
    if chosen is None:
        chosen=res
    out.append(chosen)
pd.DataFrame(out).to_csv('R9_representative_validated_certificates.csv',index=False)
with open('R9_representative_validated_certificates.json','w') as f: json.dump(out,f,indent=2)
print('\nCHOSEN')
print(pd.DataFrame(out).to_string(index=False))
