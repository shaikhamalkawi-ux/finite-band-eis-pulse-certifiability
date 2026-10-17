from pathlib import Path
import hashlib
import sys

EXPECTED = "875c266c9ac163d02803d3354453509bcd06e9255c1d063dc8ee97c9e0879e0f"
path = Path(__file__).resolve().parent / "MatFiles.zip"

if not path.exists():
    print("MatFiles.zip not found. Download the official RWTH source archive from DOI 10.18154/RWTH-2025-04811 and place MatFiles.zip in the repository root.")
    sys.exit(2)

h = hashlib.sha256()
with path.open("rb") as f:
    for chunk in iter(lambda: f.read(1024 * 1024), b""):
        h.update(chunk)
actual = h.hexdigest()
print(f"MatFiles.zip SHA-256: {actual}")
if actual != EXPECTED:
    print("FAIL: SHA-256 does not match the preserved source-integrity manifest.")
    sys.exit(1)
print("PASS: MatFiles.zip matches the preserved source-integrity manifest.")
