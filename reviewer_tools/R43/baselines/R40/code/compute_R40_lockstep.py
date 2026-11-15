#!/usr/bin/env python3
"""Project52 R40 corrected-current propagation (research calculations; separate release gate).

Cite R39 source archives and R38 R10/R20/R24 protocols; uses no new third-party raw extraction.
Zero inference of metrological calibration. A truncated grid-screen is NOT a continuum certificate.
"""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import json, pathlib, argparse, math
import numpy as np, pandas as pd
from scipy.optimize import linprog
from scipy import sparse as sp
ROOT=pathlib.Path(__file__).resolve().parent.parent/'results'
R39=ROOT.parent/'baseline'/'R39'
R38=ROOT.parent/'baseline'/'R38'

def inputs():
    return {p.stem:json.loads(p.read_text()) for p in sorted((R39/'inputs').glob('*.json'))}

def build_kernel(q,tau):
    a=np.array(q['a_s'],float);b=np.array(q['b_s'],float);u=np.array(q['u'],float);T=float(q['T_s'])
    g=np.zeros_like(tau)
    # Stable exp difference formula (reproduces generator)
    for aa,bb,v in zip(a,b,u):g+=v*np.exp(-(T-bb)/tau)*(-np.expm1(-(bb-aa)/tau))
    return np.r_[q['iT'],0.,np.sum(u*(b-a)),g]

def spectral_matrix(f,tau):
    w=2*np.pi*np.asarray(f,float);x=w[:,None]*tau[None,:];nr=len(f);nt=len(tau)
    M=np.zeros((2*nr,3+nt));M[:nr,0]=1.;M[nr:,1]=w;M[nr:,2]=-1/w
    M[:nr,3:]=1/(1+x*x);M[nr:,3:]=-x/(1+x*x)
    return M

def critical_multiplier(q,M,f):
    b=np.r_[q['Zreal_mOhm'],q['Zimag_mOhm']];eps=float(q['eps_mOhm']);n=M.shape[1]
    A=sp.vstack([sp.hstack([sp.csr_matrix(M),-eps*np.ones((len(b),1))]),
                 sp.hstack([sp.csr_matrix(-M),-eps*np.ones((len(b),1))]),
                 sp.csr_matrix(np.r_[f,0.][None,:])],format='csr')
    UB=np.r_[b,-b,q['Rmeas_mOhm']]
    c=np.r_[np.zeros(n),1.]
    r=linprog(c,A_ub=A,b_ub=UB,bounds=[(0,None)]*n+[(1,None)],method='highs')
    if not r.success:raise RuntimeError((q['case_id'],r.message))
    return float(r.fun)

def min_eps(M,b):
    n=M.shape[1];H=sp.csr_matrix(M);A=sp.vstack([sp.hstack([H,-np.ones((len(b),1))]),sp.hstack([-H,-np.ones((len(b),1))])]);
    r=linprog(np.r_[np.zeros(n),1.],A_ub=A,b_ub=np.r_[b,-b],bounds=[(0,None)]*n+[(0,None)],method='highs')
    if not r.success:raise RuntimeError(r.message)
    return float(r.fun)

def lower_on_grid(M,b,eps,f):
    A=sp.vstack([sp.csr_matrix(M),sp.csr_matrix(-M)],format='csr')
    r=linprog(f,A_ub=A,b_ub=np.r_[b+eps,-b+eps],bounds=[(0,None)]*M.shape[1],method='highs')
    if not r.success:raise RuntimeError(r.message)
    return float(r.fun)

def run_alpha(case_map, old_alpha=None):
    tau=np.logspace(-10,10,1201);rows=[]
    for i,case in enumerate(sorted(case_map),1):
        q=case_map[case];M=spectral_matrix(q['frequency_Hz'],tau);f=build_kernel(q,tau)
        new=critical_multiplier(q,M,f)
        # Invariance audit at the *same coefficients*: old f and R exactly /scale.
        factor=float(q['old_to_new_rescaling_factor']);old_f=f/factor
        qold={**q,'Rmeas_mOhm':q['Rmeas_mOhm']/factor}
        aold=critical_multiplier(qold,M,old_f)
        rows.append({'case_id':case,'SOC':q['SOC'],'direction':q['direction'],'C_rate':q['C_rate'],
                     'alpha_star_R39_grid1201':new,'alpha_star_backtransformed_R38_grid1201':aold,
                     'absolute_alpha_star_solver_difference':abs(new-aold),'rescale':factor})
        if i%7==0:print('critical-multiplier',i,'/',len(case_map),flush=True)
    d=pd.DataFrame(rows);d.to_csv(ROOT/'R40_Critical_Multiplier_42.csv',index=False,float_format='%.17g')
    fac=[1,1.1,1.25,1.5,2,3,5,8,10,15,20]
    out=pd.DataFrame([{'factor':x,'cases_excluded_grid':int((d.alpha_star_R39_grid1201>x+1e-9).sum())} for x in fac]);out.to_csv(ROOT/'R40_Critical_Multiplier_Survival.csv',index=False)
    print('ALPHA summary',d.alpha_star_R39_grid1201.describe().to_string(), 'max numerical difference',d.absolute_alpha_star_solver_difference.max(),flush=True)
    print(out.to_string(index=False),flush=True)
    return d

def run_screen(case_map):
    # Replicate original R24 cutoffs: 10,20,30,40,50,100 mHz and min grid EIS tolerance+0.0002.
    tau=np.logspace(-10,10,801);rows=[]
    for req in [0.01,0.02,0.03,0.04,0.05,0.10]:
        for soc in sorted({q['SOC'] for q in case_map.values()}):
            cases=[q for q in case_map.values() if q['SOC']==soc];q0=cases[0];F=np.array(q0['frequency_Hz']);
            mask=F>=req*.999;fm=F[mask];M=spectral_matrix(fm,tau);
            allb=np.r_[q0['Zreal_mOhm'],q0['Zimag_mOhm']];keep=np.r_[mask,mask];b=allb[keep]
            eps0=min_eps(M,b);eps=eps0+.0002
            for q in cases:
                assert np.max(abs(np.array(q['frequency_Hz'])-F))<1e-11
                assert np.max(abs(np.r_[q['Zreal_mOhm'],q['Zimag_mOhm']][keep]-b))<1e-10
                f=build_kernel(q,tau);lo=lower_on_grid(M,b,eps,f);g=lo-q['Rmeas_mOhm']
                rows.append({'requested_cutoff_mHz':req*1000,'actual_fmin_mHz':float(fm.min()*1000),'frequencies_retained':len(fm),
                             'case_id':q['case_id'],'SOC':soc,'direction':q['direction'],'C_rate':q['C_rate'],
                             'analysis_eps_mOhm':eps,'Rmeas_R39_mOhm':q['Rmeas_mOhm'],'grid_lower_R39_mOhm':lo,
                             'signed_grid_gap_R39_mOhm':g,'grid_gap_R39_mOhm':max(0.,g),'excluded_grid':bool(g>1e-9),
                             'rescale':q['old_to_new_rescaling_factor'],
                             'old_backtransformed_grid_gap_mOhm':g/q['old_to_new_rescaling_factor']})
            print('cutoff',req*1000,'SOC',soc,'cases',len(cases),flush=True)
    d=pd.DataFrame(rows);d.to_csv(ROOT/'R40_FiniteBand_Grid_Screen_All42.csv',index=False,float_format='%.17g')
    summary=d.groupby('requested_cutoff_mHz').agg(actual_fmin_mHz=('actual_fmin_mHz','first'),freq_n=('frequencies_retained','first'),
             screen_excluded=('excluded_grid','sum'),median_grid_gap_mOhm=('grid_gap_R39_mOhm','median')).reset_index()
    summary.to_csv(ROOT/'R40_FiniteBand_Grid_Summary.csv',index=False,float_format='%.17g')
    print('SCREEN',summary.to_string(index=False),flush=True)
    return d

def run_precision_bound(case_map):
    old=pd.read_csv(R38/'baseline/R10_RWTH_5C_DataDerived_Precision_All42.csv')
    cert=pd.read_csv(R39/'R39_FRESH_CERTIFICATES_84.csv')
    old['case_id']=old.apply(lambda q:f"SOC{int(q.SOC)}_{q.direction}_{q.C_rate:.1f}C",axis=1)
    out=[]
    for k,q in sorted(case_map.items()):
        p=old.set_index('case_id').loc[k]
        s=q['old_to_new_rescaling_factor']
        # Conditional conservative bound under unchanged first-order voltage/current propagation.
        # sigma_old^2=(sigma_dV/I_old)^2 + (dV*sigma_I/I_old^2)^2.
        # sigma_new^2=s^2*(voltage contribution)^2+s^4*(current contribution)^2.
        # A split-free upper bound is max(s,s^2)*sigma_old.
        h=float(p.halfwidth_3sigma_mOhm);hb=math.nextafter(h*max(s,s*s),math.inf)
        for a in [1.,2.]:
            c=cert[(cert.case_id==k)&(cert.alpha==a)].iloc[0];gap=math.nextafter(float(c.point_gap_lower_mOhm)-hb,-math.inf)
            out.append({'case_id':k,'SOC':q['SOC'],'direction':q['direction'],'C_rate':q['C_rate'],'alpha':a,
                        'rescale':s,'R38_3sigma_halfwidth_mOhm':h,'R40_3sigma_CONDITIONAL_conservative_halfwidth_ub_mOhm':hb,
                        'R39_certified_point_gap_lower_mOhm':c.point_gap_lower_mOhm,
                        'R40_CONDITIONAL_3sigma_gap_lower_mOhm':gap,
                        'positive_under_same_first_order_propagation':bool(gap>0)})
    d=pd.DataFrame(out);d.to_csv(ROOT/'R40_Conditional_3sigma_Bound_84.csv',index=False,float_format='%.17g')
    print('CONDITIONAL 3sigma\n',d.groupby('alpha').agg(n=('positive_under_same_first_order_propagation','sum'),min=('R40_CONDITIONAL_3sigma_gap_lower_mOhm','min'),median_halfwidth=('R40_3sigma_CONDITIONAL_conservative_halfwidth_ub_mOhm','median'),max_halfwidth=('R40_3sigma_CONDITIONAL_conservative_halfwidth_ub_mOhm','max')).to_string(),flush=True)
    return d

if __name__=='__main__':
    args=argparse.ArgumentParser();args.add_argument('--stage',choices=['alpha','screen','precision','all'],default='all');c=args.parse_args()
    cases=inputs();assert len(cases)==42
    if c.stage in ('alpha','all'):run_alpha(cases)
    if c.stage in ('screen','all'):run_screen(cases)
    if c.stage in ('precision','all'):run_precision_bound(cases)
