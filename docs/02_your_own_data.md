# Own-data preparation (v0.4)

Prepare exactly one biopsy. Do upstream QC/doublet review; the software validates
counts and identifiers, but does not certify cell viability, tumor malignancy or
state annotation. The example recipe is provisional, not a trained classifier.

## Formats and metadata

- H5AD: cells × genes. `--layer counts` selects a layer; `--layer raw` selects
  `.raw.X` with its own gene metadata. If omitted, X must already be raw counts.
- Cell Ranger v3+ H5: genes × barcodes on disk, transposed by the reader. Only
  Gene Expression modality is retained. Antibody counts do not become genes.
- MTX: `matrix.mtx[.gz]`, `features.tsv[.gz]` (or `genes.tsv`) and
  `barcodes.tsv[.gz]`. Matrix is genes × cells, with exact row/column alignment.
- SPECIAL neutral technical export: `counts.mtx`, headered `var.tsv` (`gene_id`),
  `obs.parquet` (`cell_id`). Install `pip install -e '.[neutral]'`. Only barcode,
  sample and technical QC columns are read, never outcomes, states or embeddings.
  This adapter is not permission to fit SPECIAL or compute paired changes.

`--metadata cells.tsv` replaces obs with an exact barcode join; missing, extra or
duplicate IDs are errors. Column order is irrelevant. Example:

```text
cell_id\tsample_id\tpatient_id\tcell_type\tcd8_state\tclone_id
AAAC-1\tP01_pre\tP01\tCD8_T\tMemory\t
AAAG-1\tP01_pre\tP01\tMelanoma\t\t
```

Use actual tabs. Lineages: `Melanoma`, `CD8_T`, `T_cell`, `NK`, `Myeloid`,
`Fibroblast`, `Endothelial`, `B_cell`, `Unknown`. CD8 requires `Memory`,
`Activated`, `Cytotoxic`, `Dysfunctional` or `Unknown`; other lineages leave state
empty. Clone is optional and never implies antigen specificity or a mother–daughter
relationship. Patient identity is nullable; one supplied patient/timepoint per
selected sample is enforced. Use `--sample-id` to explicitly subset mixed inputs.

Without reviewed states, choose `--state-recipe rank-markers-v1`. Without reviewed
lineages, choose `--recipe melanoma-markers`, which also requests state inference.
These explicit recipe options supersede existing annotation columns. The bundle
reports the choice, coverage, uncertainty and held-out-marker diagnostics.
No states or lineages are manufactured to complete a desired composition.

Gene symbols go in `var_names` or `var['gene_symbols']`. Otherwise supply
`--gene-map mapping.tsv` with `gene_id` and `gene_symbol`, one unique row per
feature ID and full coverage. Duplicate feature IDs are rejected; duplicate gene
symbols use the best positive feature rank in the rank scorer and their total
counts in the observed marker display. Species other than human require a new,
reviewed program library; this version does not silently translate orthologs.

Sampling defaults to proportional largest-remainder quotas by lineage, followed
by seeded draws without replacement. `--sampling uniform` is available. The
full and sampled composition, scheme and seed are exported. An exceptionally
rare tumor lineage can be absent after integer allocation; increase sample size.
No tumor/immune enrichment occurs silently.

## Seurat export

Export the RNA **counts**, not SCT, integrated, scale.data or log-normalized data.
One simple route (run in your existing R/Seurat environment):

```r
library(SeuratObject)
library(Matrix)
x <- LayerData(object, assay="RNA", layer="counts") # join RNA count layers explicitly if needed
stopifnot(!anyDuplicated(colnames(x)), !anyDuplicated(rownames(x)))
dir.create("my10x") # choose a new empty directory
writeMM(x, "my10x/matrix.mtx")
write.table(data.frame(rownames(x), rownames(x), "Gene Expression"),
            "my10x/features.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(colnames(x), "my10x/barcodes.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE)
metadata <- data.frame(cell_id=colnames(x), object[[]][colnames(x),,drop=FALSE])
# Add reviewed cell_type and cd8_state; remove private/outcome fields before sharing.
write.table(metadata, "cells.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

```bash
vtsc prepare-tissue my10x outputs/biopsy.zip --metadata cells.tsv --sample-id P01_pre
vtsc run-melanoma outputs/biopsy.zip outputs/run --arm anti-pd1
```

## Errors and limits

Fractional, negative, NaN and zero-total cell counts are rejected. Integer-valued
normalization cannot be proven raw from numbers alone: the user must identify
the genuine count layer. Missing feature coverage is reported, not interpreted
as absent biological activity. Unknown state has no invented effector capacity.

Bundles contain observed **selected marker** counts and source hashes, not the
entire count matrix. Retain the original source for full RNA inspection. Generated
daughters point to that immutable molecular source. Influx has no observed RNA.
Public software does not grant rights to redistribute your data; provenance
retains the supplied licence. There are no uploads or provider calls.
