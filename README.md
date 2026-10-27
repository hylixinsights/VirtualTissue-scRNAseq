# VirtualTissue-SC

Patient-level SPECIAL analysis following the execution plan in the parent
`scRNAseq` folder. Current outcome: the initial CD8 count audit cannot satisfy
the planned feasibility threshold. The user selected cohort expansion and
separate Mono/Combo strata. The deposited-object audit finds only 12 count
candidates at ≥30 CD8 per biopsy (4 Mono, 8 Combo), below the existing 10 + 4
allocation. The [cohort review](reports/cohort_expansion_review.md) and
[draft amendment](prereg/amendment-A02-cohort-expansion-draft-2026-10-27.md)
document the evidence and unresolved data dependencies. Read
[the original feasibility report](reports/T1.7_preliminary_feasibility.md) and
[task status](reports/STATUS.md) before resuming.

This is a working acquisition, provenance and cohort-audit project. Dynamics
models, state calibration, lockbox selection and spatial/Jev integration have
not been implemented or run. No state changes have been computed. Published
Sade-Feldman state summaries were exposed during an initial workbook inventory;
the review records why that cohort cannot be claimed as an untouched final test.

## Layout and data

- `config/study.yaml`: pre-specified settings, including the user-approved
  inclusion of nivolumab and pembrolizumab monotherapy.
- `prereg/`: original plan, amended frozen plan and local freeze record.
- `data/MANIFEST.json`: source URLs, versions, sizes and verification hashes.
- `data/COHORT_EXPANSION_MANIFEST.json`: candidate sources, verified Lin object and metadata derivatives.
- `../data/SPECIAL/raw/`: unmodified downloads, outside this repository.
- `../data/SPECIAL/derived/immune/`: exact per-sample RNA count and metadata exports.
- `tables/`: clinical, candidate-pair and provisional eligibility tables.
- `reports/`: task reports and saved evidence; `EVIDENCE_SHA256.json` indexes hashes.
- `r/` and `scripts/`: reproducible technical export and count-only audit.

Source-derived tables/figures carry the SPECIAL data-use policy. See
[DATA_LICENSE.md](DATA_LICENSE.md). Response fields are for stratification only.
Author embeddings are archival; never reuse them as a fold-fitted reference.

## Reproduce the environment

Python interpreter and all package versions are pinned. On this host:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock.txt
.venv/bin/python -m pip install --no-deps --no-build-isolation -e .
.venv/bin/python -m pytest -q
```

Use the R version in `renv.lock`. On a fresh checkout bootstrap renv into a local
library, then restore the pinned packages:

```sh
Rscript -e 'dir.create(".Rlib", showWarnings=FALSE); install.packages("renv", lib=".Rlib", repos="https://cloud.r-project.org")'
Rscript r/restore_environment.R
```

The exact environments were restored into separate empty environments during
T0.1. The R object export also passed a real-data round trip. Installation logs
preserve nonfatal patch-version warnings. Linux/Colab compatibility has not yet
been tested. `r/setup_environment.R` documents initial discovery; use the restore
script, not unpinned installation, for reproduction.

## Commands

Run from the repository root:

```sh
.venv/bin/vtsc status
.venv/bin/vtsc verify-freeze
.venv/bin/vtsc register GC_all_immune_sharing.rds --data-root ../data
Rscript r/inventory_seurat.R ../data/SPECIAL/raw/GC_all_immune_sharing.rds reports/evidence/immune-structure.json
```

The export requires an empty destination and checks the registered input hash:

```sh
Rscript r/export_seurat.R ../data/SPECIAL/raw/GC_all_immune_sharing.rds ../data/SPECIAL/derived/immune
.venv/bin/python scripts/extract_clinical_sources.py
.venv/bin/python scripts/audit_cohort.py
Rscript r/validate_cohort_counts.R
MPLCONFIGDIR=/tmp/vtsc-matplotlib .venv/bin/python figures/plot_feasibility.py
.venv/bin/python scripts/record_evidence.py
```

The first workbook extraction used the bundled document Python runtime in
accordance with the spreadsheet skill; `openpyxl` is also pinned locally for
standalone reproduction. Workbooks are never edited.

For cohort expansion, first run `scripts/extract_lin_clinical.py` with the
bundled document Python (`pdfplumber`), then
`.venv/bin/python scripts/audit_cohort_expansion.py`. This audits clinical
candidates and total CD8 counts only. It does not certify new patients as
eligible or freeze a new model. The separate-stratum audit is user approved.
For the Lin object, follow [metadata reproduction](docs/lin-metadata-reproduction.md),
then run `.venv/bin/python scripts/audit_lin_object_counts.py`. See the review for source versions,
PDF page/worksheet provenance and remaining checks.

The missing-data request in `docs/data-request-draft.txt` is a draft, not sent.

Do not run state analyses while the feasibility decision is unresolved. A future
resumed pipeline must implement and test the T1.4 access guard before using any
lockbox data. It is not provided by these technical-audit scripts.
