"""Download GEO GSE173706 (psoriasis skin, Ma et al. 2023) and build one raw-count AnnData.

Two steps, both idempotent:

    python scripts/download_gse173706.py download data/raw/GSE173706
    python scripts/download_gse173706.py build data/raw/GSE173706 data/psoriasis/GSE173706_raw.h5ad

`download` fetches the supplementary RAW.tar and the series matrix from the NCBI
FTP server over HTTPS, extracts the tar and writes SHA256SUMS and SOURCE.json.
`build` reads the series matrix to map each GSM to its title and characteristics,
classifies samples as lesional (PP), non-lesional (PN) or healthy (NS), drops
spatial samples, and concatenates every scRNA-seq count matrix into one AnnData
with raw integer counts in X and obs columns sample_id, gsm, donor_id, condition.
No cell-type labels exist on GEO; annotate afterwards (see
docs/datasets/psoriasis_GSE173706.md).

If automatic classification is wrong for a sample, pass --samples a TSV with
columns gsm, sample_id, donor_id, condition, assay (assay = scRNA or spatial);
rows in that file override the parsed values. The build writes the sample table
it used next to the output so the choice is recorded.

Network: needs ftp.ncbi.nlm.nih.gov (HTTPS). Data never goes into Git.
"""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import numpy as np

ACCESSION = "GSE173706"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173706"
FILES = {
    "GSE173706_RAW.tar": f"{BASE}/suppl/GSE173706_RAW.tar",
    "GSE173706_series_matrix.txt.gz": f"{BASE}/matrix/GSE173706_series_matrix.txt.gz",
}
CITATION = ("Ma F, Plazyo O, Billi AC, et al. Single cell and spatial sequencing define processes by which "
            "keratinocytes and fibroblasts amplify inflammatory responses in psoriasis. Nat Commun 14, 3455 (2023). "
            "doi:10.1038/s41467-023-39020-4")
CONDITIONS = {"PP": "lesional", "PN": "non_lesional", "NS": "healthy"}
SPATIAL_WORDS = re.compile(r"spatial|visium|\bST\b", re.I)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists: {dest.name}")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"downloading {url}")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 20)
    tmp.rename(dest)


def download(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        fetch(url, outdir / name)
    extracted = outdir / "RAW"
    if not extracted.exists():
        with tarfile.open(outdir / "GSE173706_RAW.tar") as tar:
            members = [m for m in tar.getmembers() if m.isfile() and "/" not in m.name and not m.name.startswith(".")]
            tar.extractall(extracted, members=members, filter="data")
    sums = {p.relative_to(outdir).as_posix(): sha256(p)
            for p in sorted(outdir.rglob("*")) if p.is_file() and p.name not in ("SHA256SUMS", "SOURCE.json")}
    (outdir / "SHA256SUMS").write_text("".join(f"{v}  {k}\n" for k, v in sums.items()))
    (outdir / "SOURCE.json").write_text(json.dumps({
        "accession": ACCESSION, "urls": FILES, "citation": CITATION,
        "downloaded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "files": len(sums)}, indent=1) + "\n")
    print(f"{len(sums)} files recorded in {outdir / 'SHA256SUMS'}")


def parse_series_matrix(path):
    """Return {gsm: {title, source, characteristics, library}} from the series matrix header."""
    rows = {}
    with gzip.open(path, "rt", encoding="utf8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *values = line.rstrip("\n").split("\t")
            rows.setdefault(key[len("!Sample_"):], []).append([v.strip('"') for v in values])
    gsms = rows["geo_accession"][0]
    get = lambda k, i: [r[i] for r in rows.get(k, [])]
    samples = {}
    for i, gsm in enumerate(gsms):
        samples[gsm] = {
            "title": get("title", i)[0] if "title" in rows else "",
            "source": " ".join(get("source_name_ch1", i)),
            "characteristics": " | ".join(v for v in get("characteristics_ch1", i) if v),
            "library": " ".join(get("library_strategy", i) + get("description", i) + get("data_processing", i)[:1]),
        }
    return samples


def classify(gsm, meta):
    text = " ".join([meta["title"], meta["source"], meta["characteristics"]])
    condition = next((v for k, v in CONDITIONS.items() if re.search(rf"(^|[^A-Za-z]){k}([^A-Za-z]|$)", text)), "")
    if not condition:
        low = text.lower()
        condition = ("non_lesional" if re.search(r"non[- ]?lesion|uninvolved", low)
                     else "lesional" if "lesion" in low or "psoria" in low
                     else "healthy" if re.search(r"healthy|normal|control", low) else "")
    assay = "spatial" if SPATIAL_WORDS.search(text + " " + meta["library"]) else "scRNA"
    donor = re.sub(r"[^A-Za-z0-9]+", "_", re.sub(rf"(^|[^A-Za-z])(PP|PN|NS)([^A-Za-z]|$)", r"\1\3", meta["title"])).strip("_")
    return {"gsm": gsm, "sample_id": meta["title"] or gsm, "donor_id": f"donor_{donor}" if donor else gsm,
            "condition": condition, "assay": assay}


HOUSEKEEPING = {"ACTB", "GAPDH", "B2M", "MALAT1", "RPLP0", "TMSB4X", "KRT14", "KRT10", "COL1A1", "PTPRC"}


def genes_are_columns(index, columns):
    """True when the matrix is cells x genes. Uses known genes, then barcode-like names."""
    known = lambda idx: len(HOUSEKEEPING & {str(x).split("|")[-1].upper() for x in idx})
    barcode = lambda idx: sum(bool(re.search(r"[ACGT]{10,}", str(x))) for x in list(idx)[:200])
    if known(columns) != known(index):
        return known(columns) > known(index)
    return barcode(index) > barcode(columns)


def read_counts(path):
    """Read one sample's counts; return (sparse cells x genes, cells, genes).

    Accepts a CSV/TSV matrix (optionally gzipped, either orientation, read in chunks so a
    dense genes x cells file never sits in memory as float64), a 10x .h5, or a 10x
    matrix.mtx(.gz) whose barcodes/features files sit beside it with the same prefix.
    """
    import pandas as pd
    import scipy.sparse as sp
    name = path.name
    if name.endswith(".h5"):
        import scanpy as sc
        a = sc.read_10x_h5(path)
        a.var_names_make_unique()
        return sp.csr_matrix(a.X), list(a.obs_names), list(a.var_names)
    if re.search(r"matrix\.mtx(\.gz)?$", name):
        import scipy.io
        prefix = name[:name.index("matrix.mtx")]
        side = lambda stems: next((path.parent / f"{prefix}{s}{z}" for s in stems for z in (".gz", "")
                                   if (path.parent / f"{prefix}{s}{z}").exists()), None)
        bc, ft = side(["barcodes.tsv"]), side(["features.tsv", "genes.tsv"])
        if bc is None or ft is None:
            sys.exit(f"{name}: barcodes/features files with prefix '{prefix}' not found")
        m = scipy.io.mmread(path).tocsr()  # genes x cells
        cells = pd.read_csv(bc, sep="\t", header=None)[0].astype(str).tolist()
        feats = pd.read_csv(ft, sep="\t", header=None)
        genes = feats[1 if feats.shape[1] > 1 else 0].astype(str).tolist()
        return m.T.tocsr(), cells, genes
    sep = "\t" if re.search(r"\.(tsv|txt)(\.gz)?$", name) else ","
    header = pd.read_csv(path, sep=sep, index_col=0, nrows=5)
    blocks, rows = [], []
    for chunk in pd.read_csv(path, sep=sep, index_col=0, chunksize=2000):
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))
        rows += [str(x) for x in chunk.index]
    m = sp.vstack(blocks).tocsr()
    cols = [str(c) for c in header.columns]
    if genes_are_columns(rows, cols):  # rows are cells
        return m, rows, cols
    return m.T.tocsr(), cols, rows


MATRIX = re.compile(r"\.(csv|tsv|txt)(\.gz)?$|\.h5$|matrix\.mtx(\.gz)?$")


def build(indir, out, samples_tsv=None, include_spatial=False):
    import anndata as ad
    import pandas as pd
    indir, out = Path(indir), Path(out)
    meta = parse_series_matrix(indir / "GSE173706_series_matrix.txt.gz")
    table = pd.DataFrame([classify(g, m) for g, m in meta.items()]).set_index("gsm")
    if samples_tsv:
        override = pd.read_csv(samples_tsv, sep="\t", dtype=str).set_index("gsm")
        table.update(override)
    files = {}
    for p in sorted((indir / "RAW").iterdir()):
        m = re.match(r"(GSM\d+)", p.name)
        if m and MATRIX.search(p.name) and not re.search(r"(barcodes|features|genes)\.tsv", p.name) \
                and not re.search(r"spatial|position|scalefactor|image", p.name, re.I):
            files.setdefault(m.group(1), []).append(p)
    table["files"] = [";".join(f.name for f in files.get(g, [])) for g in table.index]
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out.with_suffix(".samples.tsv"), sep="\t")
    keep = table[(table["files"] != "") & (include_spatial | (table["assay"] == "scRNA"))]
    unknown = keep[keep["condition"] == ""]
    if len(unknown):
        sys.exit(f"Condition not recognised for {list(unknown.index)}; fix them in a --samples TSV "
                 f"(template written to {out.with_suffix('.samples.tsv')}).")
    adatas = []
    for gsm, row in keep.iterrows():
        paths = files[gsm]
        if len(paths) != 1:
            sys.exit(f"{gsm} has {len(paths)} matrix files ({row['files']}); list one per GSM via --samples.")
        X, cells, genes = read_counts(paths[0])
        a = ad.AnnData(X=X.astype(np.float32), obs=pd.DataFrame(index=[f"{gsm}_{c}" for c in cells]),
                       var=pd.DataFrame(index=genes))
        a.var_names_make_unique()
        for col in ("sample_id", "donor_id", "condition", "assay"):
            a.obs[col] = row[col]
        a.obs["gsm"] = gsm
        print(f"{gsm} {row['sample_id']} {row['condition']}: {a.n_obs} cells x {a.n_vars} genes")
        adatas.append(a)
    if not adatas:
        sys.exit("No scRNA-seq matrices found; inspect the samples table.")
    adata = ad.concat(adatas, join="outer", fill_value=0)
    if not np.allclose(adata.X.data, np.round(adata.X.data)):
        print("warning: matrix values are not integers; these may be normalised, not raw counts")
    adata.layers["counts"] = adata.X.copy()
    adata.uns["source"] = {"accession": ACCESSION, "citation": CITATION,
                           "note": "Raw counts from GEO; no cell-type annotation. Re-annotate before use."}
    for col in ("sample_id", "donor_id", "condition", "assay", "gsm"):
        adata.obs[col] = adata.obs[col].astype("category")
    adata.write_h5ad(out, compression="gzip")
    print(f"wrote {out}: {adata.n_obs} cells, {adata.n_vars} genes")
    print(adata.obs.groupby("condition", observed=True)["sample_id"].nunique().to_string())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("download")
    d.add_argument("outdir")
    b = sub.add_parser("build")
    b.add_argument("indir")
    b.add_argument("out")
    b.add_argument("--samples", help="TSV overriding gsm, sample_id, donor_id, condition, assay")
    b.add_argument("--include-spatial", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "download":
        download(args.outdir)
    else:
        build(args.indir, args.out, args.samples, args.include_spatial)


if __name__ == "__main__":
    main()
