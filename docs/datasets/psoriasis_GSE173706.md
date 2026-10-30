# Psoriasis skin dataset: GEO GSE173706

Ma F, Plazyo O, Billi AC, et al. *Single cell and spatial sequencing define
processes by which keratinocytes and fibroblasts amplify inflammatory responses in
psoriasis.* Nat Commun 14, 3455 (2023).
[doi:10.1038/s41467-023-39020-4](https://doi.org/10.1038/s41467-023-39020-4).
GEO: [GSE173706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173706).

## Why this dataset

- Whole-skin biopsies (not enriched for one lineage), so cell-type frequencies can
  seed a virtual tissue (context: `library: whole_tissue`).
- Lesional (`PP`) and non-lesional (`PN`) biopsies from the **same patients**, plus
  healthy skin (`NS`). Paired donors let the lesion be tested as a "prompt"
  perturbation on the same starting tissue.
- Published as 14 PP, 11 PN and 8 NS samples. The series also contains spatial
  transcriptomics samples, which the loader drops by default.

## Download and build

**On your own computer, one command:** `bash scripts/local_psoriasis.sh` sets up a
virtual environment, downloads, builds and runs every step below (add `download` to
stop after the raw `.h5ad`). Downloads resume after an interruption: rerun the same
command. Expected files:

```
data/raw/GSE173706/GSE173706_RAW.tar              GEO supplementary archive
data/raw/GSE173706/GSE173706_series_matrix.txt.gz sample metadata
data/raw/GSE173706/RAW/                           one count matrix per GSM
data/raw/GSE173706/SHA256SUMS, SOURCE.json        checksums and provenance of this download
data/psoriasis/GSE173706_raw.h5ad                 raw counts + sample/donor/condition
data/psoriasis/GSE173706_raw.samples.tsv          sample table used (check it)
data/psoriasis/annotated/ trajectory/ communication/
data/psoriasis/psoriasis_results_small.zip        tables and reports only, easy to share
```

GEO publishes no checksums for these files; `SHA256SUMS` records the ones of your
download so that a second download can be compared with it. If you already have the
files from elsewhere, put them at the paths above and run only `build` and later steps.

Manual equivalent:

Needs `ftp.ncbi.nlm.nih.gov` reachable over HTTPS (about 250 MB download).
Data goes under `data/raw/` and `data/psoriasis/`, which are git-ignored.

```bash
python scripts/download_gse173706.py download data/raw/GSE173706
python scripts/download_gse173706.py build data/raw/GSE173706 data/psoriasis/GSE173706_raw.h5ad
```

`download` writes `SHA256SUMS` and `SOURCE.json` beside the files; commit those
two records if you want the exact download pinned, never the data.

`build` produces one AnnData:

| Field | Content |
|---|---|
| `X`, `layers["counts"]` | raw counts, cells x genes (float32 sparse) |
| `obs["gsm"]` | GEO sample accession |
| `obs["sample_id"]` | GEO sample title |
| `obs["donor_id"]` | title with the PP/PN/NS tag removed, so paired biopsies share a donor |
| `obs["condition"]` | `lesional`, `non_lesional` or `healthy` |
| `obs["assay"]` | `scRNA` (spatial samples are excluded unless `--include-spatial`) |
| `uns["source"]` | accession and citation |

It also writes `GSE173706_raw.samples.tsv`, the sample table it used. Check it
before going further: sample classification is parsed from GEO titles and
characteristics, and **the loader has only been tested on synthetic files that
mimic GEO's layout, because NCBI was not reachable from the environment where it
was written.** If a sample is misclassified, or the donor pairing is wrong, copy
the TSV, fix the rows and rerun with `--samples fixed.tsv`. The loader stops with
a message if a condition cannot be parsed, or if a GSM has more than one matrix.
It warns if values are not integers (normalised rather than raw counts).
Accepted per-sample formats: CSV/TSV (either orientation, gzipped or not, read in
chunks), 10x `.h5`, and 10x `matrix.mtx(.gz)` with `barcodes.tsv` and
`features.tsv`/`genes.tsv` sharing the same file prefix. Files that look like
spatial assets (positions, scale factors, images) are ignored.

Expected checks after the build:

- Sample counts per condition match 14 / 11 / 8 (or the subset you selected).
- Every PN sample has a PP sample with the same `donor_id`.
- Cell numbers per sample are in the range the paper reports.

## Cell-type annotation, frequencies, trajectory, communication

GEO has **no per-cell labels** for this series, so the pipeline annotates them.
Install the extras once with `pip install -e ".[scrna]"`, then run everything in
order with `bash scripts/run_psoriasis_inputs.sh`, or step by step:

| Step | Plan task | Command | Depends on | Output under `data/psoriasis/` |
|---|---|---|---|---|
| 1 | A1 | `download_gse173706.py download` then `build` | network | `GSE173706_raw.h5ad`, `GSE173706_raw.samples.tsv` |
| 2 | A2 | `sc_annotate.py` | 1 | `annotated/` |
| 3 | A3 | `sc_trajectory.py --lineage keratinocyte:keratinocyte_basal --lineage fibroblast:fibroblast_SFRP2` | 2 | `trajectory/` |
| 4 | A4 | `sc_communication.py` | 2 (parallel with 3) | `communication/` |

**Step 2, `scripts/sc_annotate.py`** (generic; works for any tissue given a marker file):

1. Per-sample QC: fewer than 200 genes, MAD outliers (5 MADs) of log counts and log
   genes, more than 20 % mitochondrial counts; Scrublet doublets per sample.
2. Normalise (10^4, log1p), 3000 HVGs (Seurat v3, per batch), PCA, Harmony on
   `donor_id`, Leiden clustering, UMAP.
3. Level 1: each Leiden cluster gets the type with the highest median marker score
   from `resources/markers/skin_v1.yaml`. Level 2: cells of each level-1 type with
   subtypes are re-clustered and assigned subtypes (keratinocyte basal / spinous /
   granular / cycling / inflammatory; fibroblast SFRP2+ / SFRP1+COL11A1+ /
   inflammatory; T CD4 / CD8 / Treg / IL17 / NK; macrophage, cDC1, cDC2, mature
   LAMP3+ DC, Langerhans). Clusters with a weak or ambiguous best score stay
   `unassigned_<cluster>` instead of being forced into a type.
4. Labels: `obs.cell_type_l1`, `obs.cell_type_l2`, `obs.cell_type` (= level 2),
   `obs.cl_id` (Cell Ontology; closest parent term where no exact term exists,
   marked `cl_exact: false` in the marker file).
5. Writes `annotation_report.md`, cluster score tables, detected marker genes, and
   frequency tables per sample and per condition (pooled fraction, mean and SD of
   per-sample fractions). The per-sample fractions are what the tissue initializer
   samples from.

Read the report before using the labels: check unassigned clusters, small margins,
and whether the detected genes per label make sense. Then compare the per-condition
fractions with the paper's figures and record the agreement.

**Step 3, `scripts/sc_trajectory.py`**: diffusion pseudotime per lineage rooted in
the named level-2 state, and PAGA connectivities between level-2 states (all cells and
per condition). `transitions_<lineage>.tsv` lists edges directed from lower to higher
median pseudotime; these become the allowed state transitions in the guidebook.

**Step 4, `scripts/sc_communication.py`**: ligand-receptor interactions per
condition. It uses liana-py (`rank_aggregate`, consensus resource) when installed;
otherwise a built-in CellPhoneDB-style score (ligand mean in source x receptor mean
in target, label-permutation p-values) over the curated skin/psoriasis table
`resources/ligand_receptor/skin_core_v1.tsv` (45 pairs: IL17, IL23, IL22, IL36, TNF,
IFN-II, chemokines, growth factors, contact pairs). The built-in table is a starting
set, not a replacement for CellChatDB; a CellChat run in R can be exported to the same
columns and used instead.

Two optional cross-checks: ask the authors for their annotated object, or transfer
labels from the annotated Reynolds et al. 2021 skin atlas (ArrayExpress E-MTAB-8142),
which also includes psoriasis lesional and non-lesional skin (plan task A5).

## Where the data lives

Nothing in `data/raw/` or `data/psoriasis/` is committed (see `.gitignore`). To
regenerate everything, rerun `scripts/run_psoriasis_inputs.sh`; `SHA256SUMS` and
`SOURCE.json` in `data/raw/GSE173706/` pin the download, and each output folder
records its parameters (`uns["annotation"]` in the h5ad, `summary.json`,
`SOURCE.json`).

## Status

The code is tested on synthetic files only (`tests/test_sc_pipeline.py`): loader
formats, QC, annotation, frequencies, trajectory and communication run end to end and
recover the planted labels. It has **not yet run on the real GSE173706 files**,
because `ftp.ncbi.nlm.nih.gov` was not reachable from the environment where it was
written.
