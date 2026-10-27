"""Hash saved outputs and sources; deliberately exclude this index itself."""
from pathlib import Path
from vtsc.provenance import file_hashes, utc_now, write_json

ROOT = Path(__file__).resolve().parents[1]
records = []
for folder in ["reports", "tables", "figures", "data", "prereg", "config", "docs"]:
    for path in sorted((ROOT / folder).rglob("*")):
        if not path.is_file() or path.name in {"EVIDENCE_SHA256.json", ".gitkeep"} or "__pycache__" in path.parts:
            continue
        records.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": file_hashes(path)["sha256"]})
for filename in ["requirements.lock.txt", "environment.yml", "renv.lock", "sessionInfo.txt", "DATA_LICENSE.md", ".python-version"]:
    path = ROOT / filename
    records.append({"path": filename, "bytes": path.stat().st_size, "sha256": file_hashes(path)["sha256"]})
write_json(ROOT / "reports/EVIDENCE_SHA256.json", {"created_at_utc": utc_now(), "files": records})
print(f"Recorded hashes for {len(records)} evidence files.")
