import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
base=Path('.')
cutoffs=[10,20,30,40,50,100]
frames=[]
for c in cutoffs:
    d=pd.read_csv(base/f'R24_FastEIS_{c}mHz_Cases.csv'); d['requested_cutoff_label_mHz']=c; frames.append(d)
allc=pd.concat(frames,ignore_index=True)
records=[]
for (soc,direction,c_rate),g in allc.groupby(['SOC','direction','C_rate']):
    g=g.sort_values('actual_fmin_Hz'); kept=g[g['excluded_grid'].astype(bool)]; last=kept.iloc[-1]
    records.append({'SOC':int(soc),'direction':direction,'C_rate':float(c_rate),'max_tested_retained_fmin_mHz':float(last['actual_fmin_Hz']*1000),'frequencies_retained_at_stop':int(last['n_freq']),'grid_gap_at_stop_mOhm':float(last['grid_gap_mOhm'])})
adaptive=pd.DataFrame(records)
proof=pd.read_csv(base/'R24_FastEIS_25mHz_ProofOriented_ValidatedCases.csv')
keys={(int(r.SOC),r.direction,float(r.C_rate)) for _,r in proof[proof['passes_validation'].astype(bool)].iterrows()}
adaptive['proof_oriented_validated_at_25p31mHz']=[(int(r.SOC),r.direction,float(r.C_rate)) in keys for _,r in adaptive.iterrows()]
eng=pd.read_csv(base/'R24_FastEIS_Engineering_Summary.csv')
proxy=dict(zip(eng['actual_fmin_mHz'].round(6),eng['one_cycle_time_proxy_reduction_pct']))
adaptive['one_cycle_proxy_reduction_pct_at_screen_stop']=adaptive['max_tested_retained_fmin_mHz'].round(6).map(proxy)
adaptive.to_csv('R25_AdaptiveFastEIS_CaseMap.csv',index=False)
rate=adaptive.groupby('C_rate').agg(median_max_tested_fmin_mHz=('max_tested_retained_fmin_mHz','median'),minimum_max_tested_fmin_mHz=('max_tested_retained_fmin_mHz','min'),maximum_max_tested_fmin_mHz=('max_tested_retained_fmin_mHz','max'),median_one_cycle_proxy_reduction_pct=('one_cycle_proxy_reduction_pct_at_screen_stop','median'),cases=('SOC','size')).reset_index()
rate.to_csv('R25_AdaptiveFastEIS_CRateSummary.csv',index=False)
conditions=[('charge',0.7),('discharge',0.7),('charge',1.0),('discharge',1.0),('charge',1.4),('discharge',1.4)]; socs=[20,30,40,50,60,70,80]
mat=np.full((6,7),np.nan); pm=np.zeros((6,7),bool)
for i,(direction,ratev) in enumerate(conditions):
    for j,soc in enumerate(socs):
        row=adaptive[(adaptive.SOC==soc)&(adaptive.direction==direction)&(np.isclose(adaptive.C_rate,ratev))].iloc[0]
        mat[i,j]=row.max_tested_retained_fmin_mHz; pm[i,j]=bool(row.proof_oriented_validated_at_25p31mHz)
fig,ax=plt.subplots(figsize=(9,4.8)); im=ax.imshow(mat,aspect='auto'); ax.set_xticks(range(7),labels=[str(s) for s in socs]); ax.set_yticks(range(6),labels=[f'{r:.1f}C {"chg" if d=="charge" else "dis"}' for d,r in conditions]); ax.set_xlabel('State of charge (%)'); ax.set_ylabel('Pulse condition'); ax.set_title('Adaptive fast-EIS screening map at 5 °C')
for i in range(6):
    for j in range(7): ax.text(j,i,f'{mat[i,j]:.1f}'+(' *' if pm[i,j] else ''),ha='center',va='center',fontsize=8)
cb=fig.colorbar(im,ax=ax); cb.set_label('Highest tested retained fmin (mHz)'); fig.tight_layout(); fig.savefig('R25_Figure_6_AdaptiveFastEIS_Map.png',dpi=300,bbox_inches='tight')
