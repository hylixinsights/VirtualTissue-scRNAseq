# Lin deposited-object audit: reproduction and limits

This is a lineage-count feasibility audit. The original `renv.lock`,
`config/study.yaml` and `prereg/FREEZE.json` remain unchanged. The separately
approved scope is in `config/cohort-expansion-audit.yaml`.

## Source and memory constraint

The source is Zenodo v2 record 19586163, `seurat_full.qs2` (CC BY 4.0 at deposit
level; upstream cohort restrictions still apply). Its published MD5 and computed
SHA-256 are recorded in `data/COHORT_EXPANSION_MANIFEST.json`. Range downloading
is implemented in `scripts/download_lin_object.py`; incomplete files are never
accepted as the source. No spatial archive is needed for this audit.

Direct R loading and the official `qs_to_rds` conversion exceeded configured
memory limits on this host; their failure logs are preserved. The successful
path decodes bounded qs2 blocks to disk, then parses R serialization while
discarding numeric arrays too large for metadata. It converts only `meta.data`,
checks complete input consumption and unique cell IDs, and refuses to convert
metadata if any of its numeric arrays were discarded. Actual matrix values and
matrix-to-cell alignment have not been verified by this route.

## Isolated reader extensions

Restore the base Python/R environment first, as described in the README. Run
from the repository root:

```sh
.venv/bin/python -m pip install --target .pydeps-expansion --no-deps -r requirements-expansion.lock.txt
Rscript r/restore_expansion_environment.R
.venv/bin/python scripts/validate_qs2_blocks.py
.venv/bin/python scripts/validate_metadata_stream.py
```

The decoder fixtures must be byte-identical to the official R converter with
and without shuffle and reject a corrupted checksum. Metadata fixtures use a
real Seurat object serialized in native and XDR formats, including factors and
missing logical values; every retained metadata value must match the R export.
A deliberately unsafe truncation limit must fail. Saved validation JSON/logs
are in `reports/evidence/cohort-expansion/`. These tests do not substitute for
a full actual-source matrix round trip. The Python overlay imports pinned rdata,
xarray, zstandard and xxhash; its other dependencies come from the base lock.

## Metadata extraction

The following destinations must not already contain an accepted output. Keep
the originals and the successful audit outputs; use a fresh directory when
reproducing. The uncompressed intermediate requires about 29 GB of disk space.

```sh
.venv/bin/python scripts/decode_qs2_blocks.py ../data/Lin2026/raw/seurat_full.qs2 ../data/Lin2026/derived/audit-quarantine/source_stream.rds
.venv/bin/python scripts/extract_seurat_metadata_stream.py ../data/Lin2026/derived/audit-quarantine/source_stream.rds ../data/Lin2026/derived/audit-quarantine/source_metadata.parquet
.venv/bin/python scripts/audit_lin_object_counts.py
```

The source's internal XXH3 and decoded stream SHA-256 are verified/recorded by
the decoder. Metadata extraction records discarded-array offsets and dimensions
and checks the metadata subtree. The provenance record
`reports/evidence/cohort-expansion/lin-reader-provenance.json` hashes the chain,
scripts, dependency locks and metadata. Reader scope and failure logs are retained.

The parquet is quarantined and retains source fields, including response and
author state labels. The count script whitelists only cohort, patient, sample,
timepoint, treatment and annotation columns. It inspects label vocabulary to
define CD8 membership, never per-state compositions or response. The source
annotations are not a fold-fitted state reference or model input.

## Interpretation and next dependency

Clinical patient-number/timepoint joins remain provisional because library IDs
are absent. P2 has two possible On dates; missing patient/timepoint labels do
not imply biological zeros. SPECIAL reused within the atlas is excluded from
new counts. Patient lists, sensitivity thresholds and input/output hashes are
saved in `lin-object-count-audit.json` and the corresponding tables.

The primary count prerequisite fails. No lockbox assignment, state-change
analysis, model fitting or spatial/Jev step is enabled. The next dependency is
new source evidence or independently acquired qualifying data, followed by a
new eligibility/design assessment. Use only public sources, as directed by the user on 2026-10-28; no author
contact or reply is a dependency. See the review and the new-conversation
prompt for the independent software implementation track.
