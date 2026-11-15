#!/usr/bin/env python3
"""Check the frozen R41 package SHA256 manifest (not including itself)."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1]
D=json.loads((R/'report'/'R41_SHA256SUMS.json').read_text());bad=[]
for name,s in D['files'].items():
 p=R/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=s:bad.append(name)
assert not bad,('SHA256 verification FAILED',bad[:15]);print('R41 SHA256 manifest: PASS',len(D['files']),'files')
