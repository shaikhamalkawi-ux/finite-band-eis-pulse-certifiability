from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parent
src=pd.read_csv(ROOT/'R19_Historical_Exploratory_FullState_Screen.csv')
interior=src[src['SOC'].between(20,80)].copy()
matrix=(interior.groupby(['temp','SOC'])['extra20'].median().unstack('SOC').sort_index().sort_index(axis=1))
matrix.to_csv(ROOT/'R19_Historical_Exploratory_Temperature_by_SOC_recomputed.csv')
summary=(interior.groupby('temp').agg(n=('extra20','size'),extra20_median=('extra20','median'),extra20_p90=('extra20',lambda x:x.quantile(0.9)),extra20_max=('extra20','max')))
print(matrix.round(6).to_string())
print('\nTemperature summary:')
print(summary.round(6).to_string())
