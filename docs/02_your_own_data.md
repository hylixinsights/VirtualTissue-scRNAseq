# Bring your own scRNA-seq data

The simulator starts from **reviewed, discrete cell-state annotations** and a
declared rate model. Prepare one sample/timepoint at a time. Upstream QC,
doublet removal, cell-type selection and biological annotation remain your
responsibility; the software does not automatically solve those analyses.

## Supported inputs

| Format | Counts | Annotations |
|---|---|---|
| H5AD | `X` or an explicitly selected layer (`--layer counts`), cells × genes | `obs`, or a replacement metadata TSV |
| 10x Matrix Market | `matrix.mtx[.gz]`, genes × cells; `barcodes.tsv[.gz]`; `features.tsv[.gz]` or `genes.tsv[.gz]` | Metadata TSV with one row for every input cell |

Only Gene Expression features are retained when a 10x feature-type column is
present. Do not provide log-normalized expression, scaled residuals or a PCA
embedding as counts. Noninteger/negative/nonfinite values, zero-count cells,
duplicate IDs and mismatched dimensions fail validation. Integer-like values
alone cannot prove a matrix's processing history: you must identify the raw
count source correctly. Direct `.rds`, 10x HDF5 and CSV count matrices are not
implemented in this version.

## Annotation contract

Minimal TSV (tabs between columns):

```text
cell_id	state	sample_id
AAAC-1	Resting-like	sample_01
AAAG-1	GZMK-like	sample_01
AAAT-1	Cytotoxic-like	sample_01
```

`cell_id` must match the exact barcode/observation identifier. External metadata
replaces H5AD `obs` rather than silently mixing conflicting columns. All matrix
cells must have a metadata row and vice versa; row order is irrelevant.

`state` must match a state in your model YAML. Use `--state-column` for another
column name. The three default names are demonstration labels; use your own
scientifically reviewed vocabulary and model where appropriate. Select the
desired lineage upstream: B, T and myeloid cell types are not interchangeable
states of one cell.

`sample_id` is strongly recommended. Multiple values require `--sample-id`;
`--sample-column` changes the field name. If no sample column exists, the input
is explicitly treated as a user-declared single sample. Conventional `patient_id`
or `timepoint` columns containing multiple values after selection are rejected.
Use one sample ID per patient/biopsy and do not pool baseline and follow-up.

## H5AD example

```python
import anndata as ad

data = ad.read_h5ad('my_annotated_cells.h5ad')
# Choose the reviewed lineage, one sample and one timepoint upstream if needed.
# Keep raw counts in data.layers['counts'] or data.X.
data.obs['state'] = data.obs['reviewed_state']
data.obs['sample_id'] = data.obs['biopsy_id']
data.write_h5ad('my_input.h5ad')
```

```bash
vtsc model-template my_model.yaml
# Edit states, transitions, division/death/influx rates and rate_provenance.
vtsc validate my_input.h5ad --layer counts
vtsc prepare my_input.h5ad outputs/my_input.zip --layer counts \
  --state-column state --sample-id sample_01 --model my_model.yaml
vtsc check outputs/my_input.zip --hours 24
vtsc run outputs/my_input.zip outputs/my_run --hours 24 --seed 42
```

If `X` already contains raw counts, omit `--layer`. The simulator consumes the
prepared labels, not raw expression trajectories. No RNA values are sent to an API.

## Seurat export without Python conversion dependencies

For a Seurat v5 object with a single raw-count layer, the following exports a
10x-style directory. Review the object/assay first. A split-layer object requires
an intentional layer selection/join before this example; do not silently choose
the first batch. The original object is not modified.

```r
library(SeuratObject)
library(Matrix)
obj <- readRDS('my_object.rds')
counts <- LayerData(obj, assay='RNA', layer='counts')
meta <- obj[[]][colnames(counts), , drop=FALSE]
stopifnot(identical(rownames(meta), colnames(counts)))
stopifnot('reviewed_state' %in% names(meta), 'sample_id' %in% names(meta))
dir.create('my_matrix', showWarnings=FALSE)
writeMM(counts, 'my_matrix/matrix.mtx')
write.table(data.frame(rownames(counts), rownames(counts), 'Gene Expression'),
            'my_matrix/features.tsv', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(colnames(counts), 'my_matrix/barcodes.tsv', sep='\t',
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(data.frame(cell_id=colnames(counts), state=meta$reviewed_state,
                       sample_id=meta$sample_id), 'my_cells.tsv', sep='\t',
            quote=FALSE, row.names=FALSE)
```

```bash
vtsc prepare my_matrix outputs/my_input.zip --metadata my_cells.tsv \
  --sample-id sample_01 --model my_model.yaml
```

The export template must be adapted/tested against your Seurat object; the new
software tests validate the resulting Matrix Market contract, not every Seurat
version/layer configuration.

## Size and sharing

Preparation loads the input matrix into memory. This is not an out-of-core atlas
pipeline. Subset a very large atlas upstream. `--max-cells 1000 --seed 21` applies
an explicit uniform subsample **after** sample/lineage selection; it does not
reduce initial input-loading memory. Source/selected/exported counts are recorded.
This subsampling is for software tests, not the frozen clinical eligibility rule.

Bundles omit expression but retain cell IDs, sample IDs and state annotations.
Use pseudonymous identifiers and share only data you are entitled to share.
Use `--source-license` to document your source terms; it records your statement,
not an automatic legal clearance. Your private data do not need a public licence
for local use.
