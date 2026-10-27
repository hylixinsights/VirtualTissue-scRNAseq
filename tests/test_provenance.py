import hashlib

from vtsc.provenance import file_hashes, inspect_download, verify_freeze


def test_truncated_download_is_rejected(tmp_path):
    path = tmp_path / "immune.rds"
    path.write_bytes(b"truncated")
    result = inspect_download(path, 1000, {"type": "MD5", "value": "0" * 32})
    assert result["status"] == "invalid"
    assert result["problems"] == ["size_mismatch", "checksum_mismatch"]


def test_same_size_corruption_is_rejected(tmp_path):
    path = tmp_path / "immune.rds"
    path.write_bytes(b"abc")
    result = inspect_download(path, 3, {"type": "MD5", "value": hashlib.md5(b"xyz").hexdigest()})
    assert result["status"] == "invalid"
    assert result["problems"] == ["checksum_mismatch"]


def test_verified_requires_source_checksum(tmp_path):
    path = tmp_path / "immune.rds"
    path.write_bytes(b"abc")
    assert inspect_download(path, 3)["status"] == "hashed_only"
    verified = inspect_download(path, 3, {"type": "MD5", "value": hashlib.md5(b"abc").hexdigest()})
    assert verified["status"] == "verified"
    assert verified["sha256"] == hashlib.sha256(b"abc").hexdigest()


def test_missing_and_in_progress_are_not_verified(tmp_path):
    assert inspect_download(tmp_path / "missing")["status"] == "missing"
    path = tmp_path / "immune.rds.part"
    path.write_bytes(b"abc")
    assert inspect_download(path)["status"] == "partial"


def test_freeze_detects_config_changes(tmp_path):
    path = tmp_path / "study.yaml"
    path.write_text("lockbox_n: 4\n")
    freeze = {"files": [{"path": "study.yaml", "sha256": file_hashes(path)["sha256"]}]}
    assert verify_freeze(tmp_path, freeze) == []
    path.write_text("lockbox_n: 0\n")
    assert verify_freeze(tmp_path, freeze) == ["study.yaml"]
    assert verify_freeze(tmp_path, {"files": []}) == ["empty_freeze"]
