# Your own dissociated melanoma biopsy

The tissue pipeline accepts one sample/timepoint at a time. Prepare upstream QC,
doublet removal and reviewed cell annotations before running it. Cell Ranger's
filtered matrix alone is not a substitute for biological annotation.

| Input | Counts | Annotation |
|---|---|---|
| Cell Ranger v3+ H5 | `matrix`, Gene Expression features | Separate metadata TSV |
| 10x Matrix Market directory | `matrix.mtx[.gz]`, `features.tsv[.gz]`, `barcodes.tsv[.gz]` | Separate metadata TSV |
| H5AD | Integer `X` or `--layer counts` | `obs` or explicit replacement metadata TSV |

Matrix Market files have genes on rows and cells on columns. H5AD has cells
on rows. Gene symbols belong in `var_names` or `var['gene_symbols']`.
Gene Expression features are selected from multi-assay 10x inputs. Normalized
noninteger values, zero-count cells, duplicate identifiers and mismatched
metadata are rejected. External metadata replaces `obs` after an exact barcode
join; its order may differ but its set of cell IDs must match the whole input.

Example `metadata.tsv` (tab-delimited):

```text
cell_id	cell_type	sample_id	patient_id	timepoint
AAAC...-1	Melanoma	P01_pre	P01	pre
AAAG...-1	CD8_T	P01_pre	P01	pre
```

Canonical types are `Melanoma`, `CD8_T`, `T_cell`, `NK`, `Myeloid`, `Fibroblast`,
`Endothelial`, `B_cell`, `Unknown`. Map your reviewed labels explicitly.
The software does not prove malignancy. Retain uncertain identities as Unknown;
no absent lineage is invented. At least one melanoma-labelled cell must survive
subsampling. A metadata `sample_id` with multiple values requires `--sample-id`;
multiple patients/timepoints within that selected sample are rejected. Without
sample metadata, the input is treated as your declaration of one biopsy.

```bash
vtsc prepare-tissue my_10x_directory outputs/my-input.zip --metadata metadata.tsv --sample-id P01_pre --max-cells 400
vtsc run-tissue outputs/my-input.zip outputs/my-control --scenario recognition-off
vtsc run-tissue outputs/my-input.zip outputs/my-recognition --scenario declared-recognition
```

`--recipe melanoma-markers` is an explicit demonstration option for unannotated
human melanoma counts. It uses the fixed panels in `src/vtsc/tissue_inputs.py`,
keeps unresolved cells, and reports composition. It is not a reference-based
classifier, CNV analysis, doublet detector or substitute for reviewed annotation.
Missing human program genes cause a clear error; do not relabel Ensembl IDs as
symbols without a real identifier mapping.

The bundle contains selected barcodes, labels, initial program scores and
source checksums. It does not contain the full expression matrix. IDs are not
anonymized automatically. Keep your own data and derived results in ignored
local directories and apply your institution's sharing rules.

For Seurat, export the RNA raw counts as a sparse Matrix Market matrix with
matching features and barcodes, plus a metadata TSV. Alternatively create an
H5AD that retains raw counts and annotations. `.rds` is not directly accepted
by this independent pipeline. The study-specific RDS audit scripts belong to
the historical SPECIAL workflow.
