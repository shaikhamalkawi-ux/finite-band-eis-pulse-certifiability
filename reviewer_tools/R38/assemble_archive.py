#!/usr/bin/env python3
"""Assemble and hash-check the R38 reviewer computational supplement from Git-stored parts.

This utility neither downloads source data nor runs a scientific calculation.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys
import zipfile


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        while True:
            block = fh.read(1024 * 1024)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def main() -> None:
    here = Path(__file__).resolve().parent
    m = json.loads((here / 'ARCHIVE_PARTS_MANIFEST.json').read_text())
    out = here / m['archive_name']
    for part in m['parts']:
        fp = here / part['filename']
        if not fp.is_file() or fp.stat().st_size != part['bytes']:
            raise RuntimeError(f"Missing or wrong-size archive part: {fp}")
        if digest(fp) != part['sha256']:
            raise RuntimeError(f"Part hash mismatch: {fp}")
    with out.open('wb') as fout:
        for part in m['parts']:
            with (here / part['filename']).open('rb') as fin:
                while True:
                    block = fin.read(1024 * 1024)
                    if not block:
                        break
                    fout.write(block)
    if out.stat().st_size != m['archive_bytes'] or digest(out) != m['archive_sha256']:
        out.unlink(missing_ok=True)
        raise RuntimeError('Combined archive checksum mismatch')
    with zipfile.ZipFile(out) as z:
        bad = z.testzip()
        if bad:
            raise RuntimeError(f'Combined archive ZIP CRC failed: {bad}')
    print('PASS:', str(out), m['archive_sha256'])


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        sys.exit(f'FAIL: {exc}')
