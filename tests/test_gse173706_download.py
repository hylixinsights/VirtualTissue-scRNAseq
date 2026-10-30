"""The GSE173706 downloader resumes partial files and retries after a dropped connection."""
import http.server
import importlib.util
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("download_gse173706", ROOT / "scripts/download_gse173706.py")
dl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dl)

PAYLOAD = bytes(range(256)) * 4096  # 1 MiB


class RangeHandler(http.server.BaseHTTPRequestHandler):
    fail_next = 0
    ranges = []

    def log_message(self, *args):
        pass

    def do_HEAD(self):
        self.send_response(200)
        self.send_header("Content-Length", str(len(PAYLOAD)))
        self.end_headers()

    def do_GET(self):
        start = int(self.headers.get("Range", "bytes=0-")[6:].rstrip("-") or 0)
        RangeHandler.ranges.append(start)
        body = PAYLOAD[start:]
        self.send_response(206 if start else 200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if RangeHandler.fail_next:  # send half, then drop the connection
            RangeHandler.fail_next -= 1
            self.wfile.write(body[: len(body) // 2])
            self.wfile.flush()
            self.connection.shutdown(2)
            return
        self.wfile.write(body)


def serve():
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_address[1]}/file.tar"


def test_fetch_resumes_after_drop(tmp_path, monkeypatch):
    monkeypatch.setattr(dl.time, "sleep", lambda s: None)
    srv, url = serve()
    try:
        RangeHandler.fail_next, RangeHandler.ranges = 1, []
        dest = tmp_path / "file.tar"
        dl.fetch(url, dest)
        assert dest.read_bytes() == PAYLOAD
        assert RangeHandler.ranges[0] == 0 and RangeHandler.ranges[1] > 0  # second request resumed
        assert not (tmp_path / "file.tar.part").exists()
        RangeHandler.ranges = []
        dl.fetch(url, dest)  # complete file is not downloaded again
        assert RangeHandler.ranges == []
    finally:
        srv.shutdown()


def test_fetch_replaces_truncated_file(tmp_path, monkeypatch):
    monkeypatch.setattr(dl.time, "sleep", lambda s: None)
    srv, url = serve()
    try:
        RangeHandler.fail_next, RangeHandler.ranges = 0, []
        dest = tmp_path / "file.tar"
        dest.write_bytes(PAYLOAD[:100])
        dl.fetch(url, dest)
        assert dest.read_bytes() == PAYLOAD
    finally:
        srv.shutdown()


def test_classify_pairs_donors_from_geo_tag():
    meta = lambda t: {"title": t, "source": "", "characteristics": "", "library": ""}
    pn = dl.classify("GSM1", meta("Peripheral normal skin from donor 3696 [PN-30696]"))  # GEO text typo
    pp = dl.classify("GSM2", meta("Psoriatic skin from donor 30696 [PP-30696]"))
    ns = dl.classify("GSM3", meta("Normal skin from donor AR001 [NS-AR001]"))
    assert (pn["condition"], pp["condition"], ns["condition"]) == ("non_lesional", "lesional", "healthy")
    assert pn["donor_id"] == pp["donor_id"] == "donor_30696"
    assert ns["donor_id"] == "donor_AR001"
