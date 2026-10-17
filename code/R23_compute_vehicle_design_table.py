import pandas as pd
from pathlib import Path
inp=Path('R19_RWTH_5C_Eq9_SlowMode_Radius_42cases.csv')
df=pd.read_csv(inp)
f_current_hz=0.01
rows=[]
for eta in [0.25,0.50,1.00,2.00]:
    fmax=eta*f_current_hz/df['slow_radius_mOhm']
    rows.append({
        'target_half_width_mOhm':eta,
        'conservative_fmin_ceiling_all42_mHz':1000*fmax.min(),
        'median_case_fmin_ceiling_mHz':1000*fmax.median(),
        'conservative_characteristic_period_s':1/fmax.min(),
        'median_case_characteristic_period_s':1/fmax.median(),
    })
pd.DataFrame(rows).to_csv('R23_Vehicular_Diagnostic_Design_Table.csv',index=False)
