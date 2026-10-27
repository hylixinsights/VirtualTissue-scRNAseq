"""Extract only the clinical table on PDF page 11; never inspect state figures.

Run with the bundled document Python (pdfplumber), or an environment containing
the versions written into the extraction record. The source PDF is not edited.
"""
import hashlib
import json
from pathlib import Path

import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
source = ROOT.parent / "data/Lin2026/raw/mmc1.pdf"
with pdfplumber.open(source) as doc:
    page = doc.pages[10]
    tables = page.extract_tables()
    assert len(tables) == 1
    assert tables[0][0][:2] == ["Patient", "SampleID"]
    payload = {
        "source": "../data/Lin2026/raw/mmc1.pdf",
        "source_url": "https://ars.els-cdn.com/content/image/1-s2.0-S1535610826003843-mmc1.pdf",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "pdf_page": 11,
        "pdfplumber_version": pdfplumber.__version__,
        "visual_review": "Page rendered and manually checked on 2026-10-27",
        "table": tables[0],
    }
(ROOT / "reports/evidence/cohort-expansion/lin-s1-extraction.json").write_text(
    json.dumps(payload, indent=2, allow_nan=False) + "\n"
)
