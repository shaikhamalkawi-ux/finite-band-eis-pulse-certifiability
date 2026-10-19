from pathlib import Path
import json, numpy as np, scipy.optimize as opt, argparse
ROOT=Path(__file__).resolve().parent

def model(q,tau):
 om=2*np.pi*np.array(q['frequency_Hz']);x=om[:,None]*tau[None,:]
 M=np.zeros((2*len(om),3+len(tau)))
 M[:len(om),0]=1;M[len(om):,1]=om;M[len(om):,2]=-1/om
 M[:len(om),3:]=1/(1+x*x);M[len(om):,3:]=-x/(1+x*x)
 a=np.array(q['a_s']);b=np.array(q['b_s']);u=np.array(q['u']);T=q['T_s']
 g=np.sum(u[:,None]*np.exp(-(T-b[:,None])/tau)*(-np.expm1(-(b-a)[:,None]/tau)),axis=0)
 f=np.r_[q['iT'],0.,np.sum(u*(b-a)),g]
 return M,f

def solve(q,alpha,margin):
 tau=np.logspace(-10,10,1001);M,f=model(q,tau);n=M.shape[0]
 b=np.r_[q['Zreal_mOhm'],q['Zimag_mOhm']];eps=q['eps_mOhm']*alpha;c=np.r_[b+eps,-b+eps]
 scale=np.max(np.abs(M),axis=0);scale=np.maximum(scale,1e-200)
 # Normalized finite constraints; conservative positive margins do not affect admissible data.
 tighten=np.r_[np.array([1e-7,1e-7,1e-7])*scale[:3],margin/(1+tau)]
 A=np.hstack([-M.T,M.T])/scale[:,None];B=(f-tighten)/scale
 r=opt.linprog(c,A_ub=A,b_ub=B,bounds=(0,None),method='highs-ds',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
 if not r.success:return {'case_id':q['case_id'],'alpha':alpha,'margin':margin,'solver_success':False,'message':r.message}
 y=r.x[:n]-r.x[n:]
 fn={'case_id':q['case_id'],'alpha':alpha,'margin':margin,'solver_success':True,'solver':'scipy highs-ds','y':y.tolist(),'floating_lower_mOhm':float(-(b@y+eps*np.sum(np.abs(y)))),'finite_slacks_float':(M[:,:3].T@y+f[:3]).tolist(),'tau_candidate_count':len(tau)}
 fn['floating_gap_mOhm']=fn['floating_lower_mOhm']-q['Rmeas_mOhm']
 return fn

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--alpha',type=float,default=1.);p.add_argument('--margin',type=float,default=.02);p.add_argument('--limit',type=int,default=42);a=p.parse_args()
 out=ROOT/'candidates';out.mkdir(exist_ok=True)
 for i,path in enumerate(sorted((ROOT/'inputs').glob('*.json'))[:a.limit]):
  q=json.loads(path.read_text());r=solve(q,a.alpha,a.margin);name=q['case_id']+f'_alpha{a.alpha:g}_margin{a.margin:g}.json';(out/name).write_text(json.dumps(r,indent=2));print(i+1,name,r.get('floating_gap_mOhm'),r.get('finite_slacks_float'),flush=True)
