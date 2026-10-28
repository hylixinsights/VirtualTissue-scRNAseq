# Melanoma software validation — v0.3.0

Date: 2026-10-27 (local). Scope: implementation tests and example execution;
not biological/clinical validation. Supersedes the v0.2 PBMC demonstration as
the primary workflow. The older validation report remains historical evidence.

## Executed checks

- **34 tests passed** in the local Python 3.12 environment. JUnit and console
  output: `melanoma-pytest.xml`, `melanoma-test-console.txt`.
- Tests include contact/recognition/fatigue gates, allowed action selection,
  no resurrection, birth/death accounting, real division events, position
  exclusion and movement limits, frozen/recognition-off controls, seeded replay,
  Cell Ranger H5 orientation and modality filtering, H5AD raw-count layers,
  sample separation, bundle integrity and offline HTML escaping.
- Three notebooks executed in real local kernels using the offline synthetic
  melanoma route. See `notebook-validation.json`. The default public-melanoma
  notebook branch uses the same functions separately exercised on the real matrix.
  Hosted Colab and remote GitHub Actions were not executed.
- Built the v0.3 wheel with its normal isolated Hatch backend. Verified the HTML
  template is packaged and ran both scenarios from an extracted wheel in an
  isolated temporary path. The first no-build-isolation attempt lacked Hatch;
  the normal isolated build succeeded.
- Verified the frozen SPECIAL plan/config hashes; no mismatches.
- Reopened and verified the real melanoma results ZIP. Replayed the final engine
  against the real input; frames, decisions, events and summary match exactly.

## Real dissociated melanoma run

Source: official 10x DTC melanoma matrix; 10,645 cells × 38,606 features;
56,774,638 bytes; locally computed SHA-256
`a2472fd15fbc708c50b95aff7f762872929a0365b9d06454cc558d29d4400091`.
Downloaded source and derivatives remain ignored local artifacts. Source licence
and the site's anatomical metadata discrepancy are recorded in the dataset guide.

Uniform 400-cell subset, seed 21; provisional markers, 15 melanoma-like cells.
No immune or tumor enrichment. Two assumed nests; 0.25 h steps; 24 h total.

| Scenario | Initial | Final alive | Attacks | Deaths | Decisions |
|---|---:|---:|---:|---:|---:|
| Declared recognition | 400 | 394 | 142 | 6 | 38,084 |
| Recognition off | 400 | 400 | 0 | 0 | 38,400 |

The 142 entries are individual damage applications, not 142 killed cells.
Both runs start from the same input and positions. No division occurred under
these default settings in this 24 h real-data example. A faster-cycle regression
test independently exercises successful division. Phagocytosis was exercised
in the synthetic example (40 events); the real run recorded timeout clearance.
These values describe the implemented assumptions, not patient outcomes.

## Browser inspection

Opened the self-contained HTML through a localhost-only server in the in-app
browser. Verified initial 400/15 counts, Find interaction jumping to 0.75 h,
contact distance 14.24 µm and an allowed attack decision; focused zoom, signal
overlay and playback; then paused/scrubbed to 24 h and verified 394 live cells,
9 melanoma-like cells, 142 attacks and 6 deaths. The screenshot in
`docs/images/melanoma_demo.png` shows the recorded CD8–melanoma contact and inspector.
Full file opening without a server is supported by the self-contained structure;
the manual UI check used localhost. No external API calls occur in the viewer.

## Remaining limits

Marker annotation/malignancy require expert review; positions and cell shapes
are assumed; rules/rates are uncalibrated; no drug or antigen inference, Jev,
ligand–receptor model, measured spatial recovery or clinical prediction is claimed.
CI is configured for Linux/Windows and Python 3.12/3.13 but only local macOS
Python 3.12 execution is confirmed here. Existing dependency deprecation warnings
are retained in the test log rather than suppressed.
