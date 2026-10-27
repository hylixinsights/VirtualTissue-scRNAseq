"""Resume the public object in bounded HTTP ranges and verify published MD5.

No object is promoted to seurat_full.qs2 until every range and its complete
source checksum match. Retry-After is respected; existing partials are retained.
"""
import concurrent.futures
import hashlib
import json
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT.parent / "data/Lin2026/raw"
META = json.loads((ROOT / "reports/evidence/cohort-expansion/zenodo-19586163.json").read_text())
ENTRY = next(x for x in META["files"] if x["key"] == "seurat_full.qs2")
URL = "https://zenodo.org/records/19586163/files/seurat_full.qs2?download=1"
SIZE = ENTRY["size"]
CHUNK = 64 * 1024 * 1024
PARTS = RAW / "download_parts"
PARTS.mkdir(exist_ok=True)


def fetch(index):
    first, last = index * CHUNK, min((index + 1) * CHUNK, SIZE) - 1
    target = PARTS / f"{index:04d}.bin"
    header = PARTS / f"{index:04d}.headers"
    expected = last - first + 1
    if target.exists() and header.exists() and target.stat().st_size == expected:
        text = header.read_text(errors="replace")
        if re.search(rf"content-range:\s*bytes {first}-{last}/{SIZE}\s", text, re.I):
            return expected
    temporary = PARTS / f"{index:04d}.partial"
    result = subprocess.run([
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--retry", "4", "--retry-delay", "5", "--connect-timeout", "30",
        "--max-time", "900", "--range", f"{first}-{last}",
        "--dump-header", str(header), "--output", str(temporary), URL,
    ], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"Range {index}: {result.stderr.strip()}")
    text = header.read_text(errors="replace")
    if not re.search(rf"content-range:\s*bytes {first}-{last}/{SIZE}\s", text, re.I):
        raise RuntimeError(f"Wrong Content-Range for part {index}")
    if temporary.stat().st_size != expected:
        raise RuntimeError(f"Wrong length for part {index}")
    temporary.rename(target)
    return expected


if __name__ == "__main__":
    completed = 0
    started = time.monotonic()
    indices = list(range((SIZE + CHUNK - 1) // CHUNK))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(fetch, i): i for i in indices}
        for future in concurrent.futures.as_completed(futures):
            completed += future.result()
            print(json.dumps({"completed_bytes": completed, "expected_bytes": SIZE,
                              "elapsed_s": round(time.monotonic() - started)}), flush=True)
    staged = RAW / "seurat_full.qs2.assembled"
    md5, sha = hashlib.md5(), hashlib.sha256()
    with staged.open("wb") as out:
        for i in indices:
            with (PARTS / f"{i:04d}.bin").open("rb") as source:
                while block := source.read(8 * 1024 * 1024):
                    md5.update(block)
                    sha.update(block)
                    out.write(block)
    assert staged.stat().st_size == SIZE
    assert f"md5:{md5.hexdigest()}" == ENTRY["checksum"], "Published checksum mismatch"
    staged.rename(RAW / "seurat_full.qs2")
    result = {"source_url": URL, "bytes": SIZE, "md5": md5.hexdigest(),
              "sha256": sha.hexdigest(), "status": "verified", "http_ranges": len(indices)}
    (ROOT / "reports/evidence/cohort-expansion/lin-object-integrity.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result), flush=True)
