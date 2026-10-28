# A05 data audit and unresolved sources

The baseline inventory is historical (`baseline.json`). It records both verified
sources and partial/quarantined intermediates. This addendum records new work;
it does not rewrite accepted immune exports or the original clinical inventory.

| Source | Current product status | Evidence |
|---|---|---|
| 10x melanoma | Verified raw H5; 10645 cells; provisional lineage/state annotation; 400-cell stratified example | melanoma-input-audit.json |
| SPECIAL immune | Previously verified 25178-cell export; neutral adapter reads a sample preserving order; clinical/state metadata not loaded | neutral-adapter.json; original immune-export-validation.json |
| Entire TME | Verified 2718597537 bytes and published MD5; 59170 cells, 22897 genes, 46 samples exported into a new directory | tme-acquisition.json; tme-export.json |
| Malignant cells | Verified 864438613 bytes and published MD5; 14208 cells, 22897 genes, 35 samples exported into a new directory | tme-acquisition.json; malignant-export.json |
| Lin full qs2 | Previously verified 4 GB source; metadata-only extraction; actual expression alignment remains unverified | ../cohort_expansion_review.md; ../evidence/cohort-expansion/lin-object-integrity.json |
| Lin partial download | Invalid historical intermediate; never used as a matrix | baseline.json |
| Lin 28.6 GB expanded stream | Quarantined audit asset; not reconverted; prior full-memory failure remains concrete | ../cohort_expansion_review.md; ../../docs/lin-metadata-reproduction.md |
| Sade-Feldman / other complementary cohorts | Prior public audit retained; no newly established eligible early paired cohort | ../cohort_expansion_review.md |

`r/export_tme_blind.R` fixes the assay and export columns before reading the
objects, requires a verified manifest source, recomputes technical QC, and
round-trips sparse values, dimensions, sums and barcode order. A coarse T/CD8
marker screen is fixed and labelled provisional. It does not assign the four
states, export author embeddings, inspect pre/on changes or read response
columns. Exports are CC BY-NC-SA 4.0, outside Git. Per-file hashes are in the
export reports and the product export manifest. R used R_MAX_VSIZE=12Gb.

The Python neutral adapter reads only identity/sample/QC parquet columns before
constructing AnnData; author states/outcomes/embeddings are not exposed through
this route. Three real sample exports were read to verify alignment and totals.
The strict own-data path still requires reviewed lineages and CD8 states or an
explicit provisional annotation recipe; a neutral export alone is not validated
annotation.

No additional clinical pair was manufactured from metadata. The 19 reused
SPECIAL patients in the Lin atlas are still excluded from new-patient totals.
Missing library IDs, P2 on-biopsy ambiguity, exact drugs, lesion matching and
chemistry conflict remain unresolved in the earlier public audit. No author
contact was made. No reprocessing of the 29 GB stream was attempted. P1.3 stays
partial until expression alignment is actually demonstrated. The frozen
scientific feasibility/gates and original T-task status are unchanged.
