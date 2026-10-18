#!/usr/bin/env python3
"""Check or replay the exact R37 archive supplied with the reviewer kit.
The manuscript, mathematical code, witnesses and results are not modified.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import zipfile

EXPECTED_SHA256 = "277c73ff96923904d9eff9980dece66b52d4484add76b1913b9944962a06e1a1"
EXPECTED_ROOT = "Batteries52_Reviewer_Replay"


def verify(archive: Path) -> dict[str, object]:
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError("Wrong or changed archive. Expected SHA-256 " + EXPECTED_SHA256)
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
        if len(names) != 746 or len(set(names)) != len(names):
            raise ValueError("Unexpected or duplicate archive entries")
        for info in z.infolist():
            p = Path(info.filename)
            if p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] != EXPECTED_ROOT:
                raise ValueError("Unsafe archive path")
            if stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError("Symlink archive member is not allowed")
        bad = z.testzip()
        if bad:
            raise ValueError("ZIP CRC failure: " + bad)
        manifest = json.loads(z.read(EXPECTED_ROOT + "/MANIFEST.json"))
        if len(manifest) != 745:
            raise ValueError("Unexpected manifest size")
        for name, expected in manifest.items():
            actual = hashlib.sha256(z.read(EXPECTED_ROOT + "/" + name)).hexdigest()
            if actual != expected:
                raise ValueError("Internal hash mismatch: " + name)
    return {"archive_sha256": digest, "archive_entries": len(names),
            "verified_manifest_files": len(manifest), "integrity": "PASS"}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("--run", action="store_true", help="Compile and replay after integrity checks")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--output-parent", type=Path, default=Path.cwd())
    args = p.parse_args()
    if not 1 <= args.workers <= 32:
        p.error("--workers must be 1..32")
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE", "") not in ("", "0"):
        p.error("Do not use Python -O or PYTHONOPTIMIZE: inherited checks use assertions")
    print(json.dumps(verify(args.archive), indent=2), flush=True)
    if not args.run:
        return
    args.output_parent.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix="R37_review_", dir=args.output_parent.resolve()))
    with zipfile.ZipFile(args.archive) as z:
        z.extractall(out)
    root = out / EXPECTED_ROOT
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONOPTIMIZE="0")
    subprocess.run([sys.executable, "run_review.py", "--workers", str(args.workers)],
                   cwd=root, env=env, check=True)
    report = root / "review_replay" / "REVIEW_REPLAY_SUMMARY.json"
    print("REPLAY REPORT:", report)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, subprocess.CalledProcessError) as exc:
        print("FAIL:", exc, file=sys.stderr)
        sys.exit(1)
