#!/usr/bin/env python3
"""Reproduce the R43 laboratory example from unchanged R42 stored inputs.

This is arithmetic on one saved certificate, not a new measurement or proof of
compatibility when its sufficient exclusion condition fails. No optimizer is used.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from decimal import Decimal, localcontext, ROUND_FLOOR
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASE = 'SOC30_discharge_0.7C'

def build(baseline: Path) -> dict:
    source = baseline / 'evidence/R42_Raw3sigma_Certificate_Gaps_84.csv'
    rows = list(csv.DictReader(source.open(encoding='utf-8', newline='')))
    selected = [r for r in rows if r['case_id'] == CASE]
    if len(selected) != 2 or {r['alpha'] for r in selected} != {'1', '2'}:
        raise ValueError('Expected the same case at both frozen EIS tolerance factors.')
    records = []
    with localcontext() as ctx:
        ctx.prec = 75
        for row in sorted(selected, key=lambda r: int(r['alpha'])):
            D = lambda k: Decimal(row[k])
            h = Decimal(3) * ((Decimal(1000) / D('I_p_A')) ** 2 *
                  (D('sigma_baseline_V')**2 + D('sigma_endpoint_V')**2) +
                  (D('Rmeas_mOhm') * D('sigma_endpoint_I_A') / D('I_p_A'))**2).sqrt()
            if abs(h - D('halfwidth_3sigma_endpoint_mOhm')) > Decimal('1e-12'):
                raise ValueError('Independent local-scatter arithmetic disagreement.')
            if D('point_gap_lower_mOhm') > D('certified_lower_mOhm')-D('measured_mOhm')+Decimal('1e-12'):
                raise ValueError('Stored point-gap bound is not conservative.')
            gap = D('gap_after_3sigma_endpoint_mOhm')
            if abs(gap - (D('point_gap_lower_mOhm')-D('halfwidth_3sigma_endpoint_mOhm'))) > Decimal('1e-14'):
                raise ValueError('Stored expanded-gap arithmetic disagreement.')
            scenarios = []
            for rho in map(Decimal, ['0', '0.005', '0.010']):
                residual = gap-rho
                scenarios.append({'additional_pulse_only_bound_mOhm': str(rho),
                    'remaining_saved_witness_margin_mOhm': str(residual),
                    'remaining_margin_downward_6dp_mOhm': str(residual.quantize(Decimal('.000001'), rounding=ROUND_FLOOR)),
                    'saved_witness_excludes': residual > 0,
                    'interpretation': 'Exclusion retained for this declared scenario' if residual > 0 else
                      'This saved witness does not certify exclusion; compatibility is not established'})
            records.append({'case_id': CASE, 'alpha': int(row['alpha']),
              'measured_response_mOhm': row['measured_mOhm'],
              'saved_lower_bound_mOhm': row['certified_lower_mOhm'],
              'local_halfwidth_mOhm': row['halfwidth_3sigma_endpoint_mOhm'],
              'local_expanded_gap_mOhm': str(gap),
              'local_expanded_gap_microOhm': str(gap * 1000),
              'current_normalizer_A': row['I_p_A'],
              'scenarios': scenarios})
    if [s['saved_witness_excludes'] for s in records[1]['scenarios']] != [True,True,False]:
        raise AssertionError('Doubled-tolerance worked-case decisions changed.')
    return {'status':'PASS', 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'case_selection':'Weakest endpoint-window local-scatter margin at doubled EIS tolerance in R42',
            'scope':'Illustration from existing data; hypothetical extra pulse-only bounds. EIS center, tolerance, waveform and state matching held fixed. No calibrated coverage, new experiment, or tested decision-policy performance.',
            'records':records}

def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    default=ROOT/'baseline_R42' if (ROOT/'baseline_R42').is_dir() else ROOT
    ap.add_argument('--baseline', type=Path, default=default)
    ap.add_argument('--write', action='store_true', help='Write the example evidence files; otherwise read-only.')
    args=ap.parse_args(); data=build(args.baseline)
    output=ROOT/'evidence/R43_Worked_Engineering_Example.json'
    if args.write:
        output.parent.mkdir(exist_ok=True,parents=True)
        output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
        csvpath=output.with_suffix('.csv')
        flat=[]
        for row in data['records']:
            for sc in row['scenarios']:
                flat.append({k:v for k,v in row.items() if k!='scenarios'}|sc)
        with csvpath.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    elif output.exists() and json.loads(output.read_text()) != data:
        raise AssertionError('Saved worked-example evidence differs from current reproduction.')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
