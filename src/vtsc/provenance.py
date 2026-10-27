"""Streaming integrity checks; incomplete downloads cannot become verified data."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def file_hashes(path):
    digests = {name: hashlib.new(name) for name in ("sha256", "md5")}
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            for digest in digests.values():
                digest.update(chunk)
    return {name: digest.hexdigest() for name, digest in digests.items()}


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temporary.replace(path)


def inspect_download(path, expected_bytes=None, expected_checksum=None):
    path = Path(path)
    if not path.is_file():
        return {"status": "missing", "bytes": None, "sha256": None}
    if path.name.endswith(".part"):
        return {"status": "partial", "bytes": path.stat().st_size, "sha256": None}
    size = path.stat().st_size
    hashes = file_hashes(path)
    problems = []
    if expected_bytes is not None and size != expected_bytes:
        problems.append("size_mismatch")
    if expected_checksum:
        algorithm = expected_checksum["type"].lower().replace("-", "")
        if algorithm not in hashes:
            raise ValueError(f"Unsupported checksum: {algorithm}")
        if hashes[algorithm] != expected_checksum["value"].lower():
            problems.append("checksum_mismatch")
    return {
        "status": "invalid" if problems else ("verified" if expected_checksum else "hashed_only"),
        "bytes": size,
        **hashes,
        "problems": problems,
        "checked_at_utc": utc_now(),
    }


def verify_freeze(root, freeze):
    root = Path(root).resolve()
    failures = []
    for record in freeze["files"]:
        path = (root / record["path"]).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Freeze paths must be inside the repository")
        if not path.is_file() or file_hashes(path)["sha256"] != record["sha256"]:
            failures.append(record["path"])
    if not freeze["files"]:
        failures.append("empty_freeze")
    return failures
