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

Expected checks after the build:

- Sample counts per condition match 14 / 11 / 8 (or the subset you selected).
- Every PN sample has a PP sample with the same `donor_id`.
- Cell numbers per sample are in the range the paper reports.

## Cell-type annotation is required

GEO has **no per-cell labels** for this series. Before a guidebook can be built
(plan task A2):

1. QC per sample (genes, UMIs, mitochondrial fraction), doublet removal.
2. Integrate across donors (for example scVI or Harmony with `donor_id` as batch).
3. Cluster and annotate with the paper's marker sets, at the resolution the
   guidebook needs: keratinocyte layers (basal KRT5/KRT14, spinous KRT1/KRT10,
   supraspinous/granular, proliferating MKI67), psoriasis-activated keratinocytes
   (KRT16, KRT6A, S100A7/8/9, SERPINB4), fibroblasts including SFRP2+ and
   WISP2+ subsets, endothelial cells, pericytes/smooth muscle, melanocytes,
   T-cell subsets (CD4, CD8, Treg, IL17A+ Tc17/Th17), NK/ILC, dendritic cells
   (including LAMP3+ mature DC), macrophages, mast cells.
4. Map labels to the canonical vocabulary with Cell Ontology IDs (plan task B2).
5. Compare per-condition type frequencies with the paper's figures; record the
   agreement in an annotation report.

Two shortcuts, both optional: ask the authors for their annotated object, or
transfer labels from the annotated Reynolds et al. 2021 skin atlas (ArrayExpress
E-MTAB-8142), which also includes psoriasis lesional and non-lesional skin.
