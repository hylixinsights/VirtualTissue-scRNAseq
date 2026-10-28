# Data licence and permitted project use

## Melanoma software examples (v0.4; v0.3 retained)

Original code and generated synthetic melanoma fixtures are MIT-licensed. The
primary real example is the official 10x Genomics dissociated melanoma matrix,
explicitly [CC BY 4.0](https://www.10xgenomics.com/datasets/10k-human-dtc-melanoma-GEM-X).
Source identifiers, attribution, license URL, a locally computed checksum and
modifications are recorded in `examples/datasets.json`, downloaded `SOURCE.json`
and derived input/run provenance. Retain these notices with shared derivatives.
See `docs/05_example_datasets.md` for attribution text and the source discrepancy.

The v0.2 PBMC example remains a legacy regression fixture with its original
CC BY 4.0 terms; it is not a melanoma substitute. Own-data inputs retain their
source permissions. A public software repository does not grant redistribution
rights to private data. No SPECIAL/Lin input is required for the tissue demo.

## Historical SPECIAL research assets

The SPECIAL processed-data deposit (Joanna Poźniak, KU Leuven RDR,
[doi:10.48804/GSAXBN](https://doi.org/10.48804/GSAXBN), version 8) declares
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
The deposited metadata snapshot is `reports/evidence/rdr-v8-metadata.json`.

For SPECIAL-derived matrices, metadata exports, avatars and simulation-input
bundles, this project retains attribution, the licence link, and a description
of modifications; distributes adaptations under CC BY-NC-SA 4.0; and prohibits
commercial use. Keep these assets separate from Hylix commercial services.
Do not label them MIT. The licence does not require new, independent analysis
code to use the data licence or change the existing VirtualTissue code licence.

Suggested attribution: “Processed single-cell data from the SPECIAL cohort
(Poźniak et al., Cell 2024; Poźniak, Roda et al., Cancer Discovery 2025),
Joanna Poźniak, KU Leuven RDR, doi:10.48804/GSAXBN, version 8,
CC BY-NC-SA 4.0. Changes: [describe the export, filtering or transformation].”

Clinical article supplements retain their own source terms, recorded separately
in the manifest; the dataset licence is not automatically assigned to them.
Only public processed data are acquired here. No EGA controlled-access raw
sequencing files or direct patient identifiers are requested. Study pseudonyms
are retained. Original downloads and derived data live outside the Git repository.

This is the project data-use policy. Any release must include the applicable
source notices and transformation provenance for each asset.

## Cohort-expansion source screening

Lin2026 Zenodo v2 (`doi:10.5281/zenodo.19586163`) declares CC BY 4.0 in its
downloaded metadata. Its integrated object also incorporates upstream studies,
including SPECIAL according to Table S2A. Preserve source-specific terms;
do not treat the deposit-level licence as proof that SPECIAL-derived content
has lost its CC BY-NC-SA restrictions. Article supplements, GEO files and
preprints keep their own source terms. Unverified terms remain unresolved in
`data/COHORT_EXPANSION_MANIFEST.json`; this screening is not a redistribution
licence audit. No candidate data have been released or assigned the code's MIT
licence.

## A05 outcome-blind exports (2026-10-28)

Entire_TME.rds and Malignant_cells.rds were verified against published size/MD5
and locally hashed. New exports are in `../data/SPECIAL/derived/tme-A05-v1` and
`malignant-A05-v1`. They retain CC BY-NC-SA 4.0 and source attribution and are
not included in the MIT software wheel. Metadata export is restricted to IDs,
recomputed technical QC and a fixed coarse CD8 marker screen. No author state,
embedding or clinical outcome was used to tune the product. See the two
`reports/product-v04/*-export.json` files for per-file checksums and round trips.
