# Software validation — v0.2.0, 2026-10-27

This report concerns implementation and user workflows, not clinical validation.
Source/result hashes are indexed in `reports/EVIDENCE_SHA256.json`; the compact
execution record is `reports/software/validation.json`.

| Check | Result / scope |
|---|---|
| Fresh isolated installation | Python 3.12, macOS ARM; `pip install -e '.[test,notebooks]'` succeeded and `pip check` passed |
| Automated tests | 25 tests passed, including original integrity tests and numerical/input/result tests; custom own-data H5AD counts layer/state model covered; `pytest.xml` |
| Mathematical checks | Exact two-state transition expectation; linear birth/death/influx expectation; Monte Carlo mean within a predeclared five-standard-error tolerance |
| Seeded synthetic workflow | 300 initial cells → 302 final cells, 394 events over 48 simulated hours, seed 21 |
| No-events control | 300 initial/final cells, no events, full horizon completed |
| PBMC3k source | 7,621,991-byte archive matched pinned SHA-256; downloaded from 10x |
| PBMC3k import/prepare/run | 2,700 × 32,738 source matrix; fixed demo recipe selected 427 cells; 417 final cells and 490 events at 48 h, seed 21 |
| Notebooks | All three executed successfully through real local kernels, synthetic route, with package already installed; `notebook-validation.json` |
| Browser check | Offline-capable explorer served locally: initial 300 cells; final slider shows 48 h / 302 cells; clicking synthetic_cell_00000 shows its GZMK-like → Cytotoxic-like event at 5.865 h |
| Original research freeze | Still valid; no clinical gate changed |

The browser screenshot is `docs/images/explorer_demo.png`. The independent
static scientific figure is `docs/images/synthetic_demo.png`. Both use generated
synthetic data. Large real input matrices and run folders are excluded from Git.

Observed nonfatal warnings concern upstream AnnData/Pandas deprecations and the
deliberate duplicate-ID rejection fixture. Jupyter warned about local TCP kernel
transport; notebook execution used temporary local kernels and synthetic data.
No model API was called and no private source data were uploaded.

Not yet verified: hosted Google Colab sessions, actual GitHub Actions jobs,
Linux/Windows runtime behavior, Python 3.13 runtime behavior, PBMC1k import,
every Seurat version/export variant, large-atlas performance or biological
validity of state/rate choices. CI is configured, not claimed as already green.

Reproduce locally:

```bash
python -m pip install -e ".[test,notebooks]"
python -m pytest -q
python scripts/validate_notebooks.py
vtsc verify-freeze
vtsc demo outputs/retest --cells 300 --hours 48 --seed 21
```

For the real-data path, use the exact PBMC commands in README. The source
downloader rejects checksum/size changes; do not update its checksum merely
to accept an unexplained new file. Test on new output paths. Reopening existing
outputs verifies checksums instead of overwriting them.
