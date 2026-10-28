# Current v0.4 sampling note

The current `prepare-tissue` defaults to proportional stratified sampling by
lineage (largest-remainder allocation), with explicit seed and no enrichment.
The 400-cell example has 27 CD8 and 8 provisionally melanoma-labelled cells.
The full 10x matrix has 10645 cells; marker annotation retains Unknown and does
not confirm malignancy by CNV. Rank-state coverage and held-out-marker audits
are in `reports/product-v04/melanoma-input-audit.json`. The exact v0.4 workflow
is in the README; the following source attribution and v0.3 history are retained.

# Melanoma example and provenance

The primary example is [10k Human DTC Melanoma, Chromium GEM-X Single Cell 3′](https://www.10xgenomics.com/datasets/10k-human-dtc-melanoma-GEM-X),
from 10x Genomics. Its description explicitly identifies dissociated melanoma
tumor cells. The page declares CC BY 4.0. Download and licence verified on
2026-10-27; pinned source metadata are in `examples/datasets.json`.

The filtered Cell Ranger 8.0.0 matrix was downloaded and inspected: 10,645 barcodes,
38,606 Gene Expression features, 56,774,638 bytes. SHA-256:
`a2472fd15fbc708c50b95aff7f762872929a0365b9d06454cc558d29d4400091`.
The hash was computed locally from the official download; it is not a checksum
published by the vendor. The downloader rejects a differing file.

The source's structured Anatomical Entity field says “Blood”, while the title,
description and selected cell type say dissociated melanoma. This inconsistency
is recorded in SOURCE.json. No biopsy site, tumor purity, clinical response,
measured geometry or paired treatment design is invented from this page.

The fixed demonstration recipe produces the following **provisional marker
labels**, not vendor annotations or confirmed malignancy:

| Label | Full matrix | Uniform 400-cell subset, seed 21 |
|---|---:|---:|
| Melanoma-like | 219 | 15 |
| CD8 T | 723 | 31 |
| Other T | 3,109 | 100 |
| B cell | 5,402 | 203 |
| Myeloid | 167 | 10 |
| NK | 59 | 2 |
| Endothelial | 14 | 2 |
| Fibroblast | 12 | 1 |
| Unknown | 940 | 36 |

This is a lymphocyte-rich example. The software does not enrich the malignant
compartment to make a more impressive image. The full and sampled compositions
are both exported, and users should replace provisional labels with reviewed
annotations for their own experiments. Uniform subsampling changes exact proportions.

Attribution for distributed examples/results: “10x Genomics (2024), 10k Human
DTC Melanoma, Chromium GEM-X Single Cell 3′, Cell Ranger 8.0.0,
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: uniform
subsampling, illustrative marker annotation/program scoring and assumed spatial
agent simulation by VirtualTissue-scRNAseq.”

The offline fixture is generated melanoma/microenvironment counts, licensed MIT,
with no patient measurements. It is used for CI and regression tests.

PBMC3k is a legacy v0.2 input-format fixture only, reachable with
`download-example --dataset pbmc3k-legacy`. It is not a melanoma example and
is absent from the primary workflow and notebooks. SPECIAL/Lin are a separate
clinical research audit with their own source terms and feasibility gates.
