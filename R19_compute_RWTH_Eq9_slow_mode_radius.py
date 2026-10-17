import io, math, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io as sio

# Inputs expected beside the script or supplied by path:
ROOT = Path(__file__).resolve().parent
TAB = ROOT / 'R4D_RWTH_5C_ContinuumSafe_42cases.csv'
# The official RWTH MatFiles.zip is not redistributed in this compact package.
# Place it at ROOT/MatFiles.zip or edit this path to the preserved raw archive.
ZP = ROOT / 'MatFiles.zip'
OUT = ROOT / 'R19_RWTH_5C_Eq9_SlowMode_Radius_42cases_recomputed.csv'

if not ZP.exists():
    raise FileNotFoundError('Place the preserved official MatFiles.zip beside this script before replay.')

tab = pd.read_csv(TAB)
with zipfile.ZipFile(ZP) as z:
    ov = sio.loadmat(io.BytesIO(z.read('5deg/5deg_Overview.mat')), squeeze_me=False, struct_as_record=False)['diga'][0,0].daten[0,0]

t = np.asarray(ov.Zeit).squeeze()
I = np.asarray(ov.Strom).squeeze()
steps = np.asarray(ov.Schritt).squeeze()
on = np.abs(I) > 0.5
d = np.diff(on.astype(int))
starts = np.where(d == 1)[0] + 1
ends = np.where(d == -1)[0] + 1
if on[0]: starts = np.r_[0, starts]
if on[-1]: ends = np.r_[ends, len(on)]
segs=[]
for s,e in zip(starts,ends):
    dur=t[e-1]-t[s]; med=np.median(I[s:e])
    if 15 <= dur <= 25 and abs(med) > 1.5:
        segs.append((s,e,dur,med,int(steps[s])))
groups=[]; cur=[]; prev=None
for sg in segs:
    if prev is not None and sg[0]-prev > 10000:
        groups.append(cur); cur=[]
    cur.append(sg); prev=sg[0]
groups.append(cur)

soc_to_group={90:1,80:2,70:3,60:4,50:5,40:6,30:7,20:8,10:9,0:10}
step_map={('charge',0.7):3,('discharge',0.7):8,('charge',1.0):13,('discharge',1.0):18,('charge',1.4):23,('discharge',1.4):28}

rows=[]
for _,r in tab.iterrows():
    soc=int(r.SOC); direction=str(r.direction); rate=float(r.C_rate); T=20.0
    st=step_map[(direction,rate)]
    g=groups[soc_to_group[soc]]
    s,e=[x[:2] for x in g if x[4] == st][0]
    tt=t[s:e]-t[s]; ii=I[s:e]
    Iref=float(np.median(ii[-10:]))
    norm=ii/Iref
    q=0.0
    for j in range(len(tt)-1):
        a=float(tt[j]); b=float(tt[j+1])
        if a >= T: break
        b=min(b,T)
        if b <= a: continue
        q += float(norm[j])*(b-a)
        if b >= T: break
    eps=float(r.eps_mOhm)
    fmin=0.01
    radius=2*math.pi*eps*fmin*abs(q)
    rows.append({**r.to_dict(), 'Q_T_s':q, 'f_min_Hz':fmin, 'slow_radius_mOhm':radius, 'slow_to_gap_ratio':radius/float(r.gap_mOhm)})

out=pd.DataFrame(rows)
out.to_csv(OUT,index=False)
print(out['slow_radius_mOhm'].describe().to_string())
