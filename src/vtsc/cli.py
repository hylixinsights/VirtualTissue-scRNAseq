"""Project status, integrity registration and immutable pre-analysis freeze."""

import argparse
import json
from pathlib import Path

from .provenance import file_hashes, inspect_download, utc_now, verify_freeze, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    sub.add_parser("verify-freeze")
    register = sub.add_parser("register")
    register.add_argument("filename")
    register.add_argument("--data-root", type=Path, required=True)
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--file", action="append", required=True)
    freeze.add_argument("--commit", required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.command == "status":
        print((root / "reports/STATUS.md").read_text())
        return
    if args.command == "verify-freeze":
        payload = json.loads((root / "prereg/FREEZE.json").read_text())
        failures = verify_freeze(root, payload)
        print(json.dumps({"valid": not failures, "mismatches": failures}))
        raise SystemExit(bool(failures))
    if args.command == "freeze":
        target = root / "prereg/FREEZE.json"
        if target.exists():
            parser.error("Freeze already exists. Record an amendment; never overwrite it.")
        records = []
        for name in args.file:
            path = (root / name).resolve()
            if not path.is_relative_to(root):
                parser.error("Freeze inputs must be inside the repository")
            records.append({"path": str(path.relative_to(root)), "sha256": file_hashes(path)["sha256"]})
        write_json(target, {"schema_version": 1, "stage": "T0.4", "created_at_utc": utc_now(),
                            "source_commit": args.commit, "public_registration": False, "files": records})
        print("Created prereg/FREEZE.json (local freeze, not public registration).")
        return
    manifest_path = root / "data/MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    records = [r for r in manifest["files"] if r["filename"] == args.filename]
    if len(records) != 1:
        parser.error("Filename must identify one manifest record")
    record = records[0]
    data_root = args.data_root.resolve()
    path = (data_root / record["relative_path"]).resolve()
    if not path.is_relative_to(data_root):
        parser.error("Data path escapes data root")
    result = inspect_download(path, record.get("expected_bytes"), record.get("source_checksum"))
    # Never silently replace a previously recorded local hash after modification.
    previous_hash = record.get("local", {}).get("sha256")
    if previous_hash and previous_hash != result.get("sha256"):
        parser.error("Registered file has changed; investigate and log a deviation")
    record["local"] = result
    write_json(manifest_path, manifest)
    print(json.dumps({"filename": args.filename, **result}))
    raise SystemExit(result["status"] not in {"verified", "hashed_only"})


if __name__ == "__main__":
    main()
