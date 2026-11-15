#!/usr/bin/env python3
"""Fresh replay of all 48 archived Bayreuth MPFR jobs (no raw download).
Uses the unchanged, manifest-checked original program with only portability/
resource-limit substitutions. Specify a new --output directory and --workers 4.
"""
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[1]/'baselines/External_E3'
p=root/'code/replay_external.py'
s=p.read_text()
assert 'timeout=300' in s and 'ROOT=Path(__file__).resolve().parents[1]' in s
s=s.replace('ROOT=Path(__file__).resolve().parents[1]',f'ROOT=Path({str(root)!r})')
s=s.replace('timeout=300','timeout=1200')
sys.path.insert(0,str(root/'code'))
exec(compile(s,str(p),'exec'))
