"""End-to-end checks of the single-cell input pipeline on synthetic skin-like counts.

The synthetic data only exercise the code paths (loader formats, QC, two-level marker
annotation, frequencies, trajectory, communication); they say nothing about real biology.
"""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sc = pytest.importorskip("scanpy")
ad = pytest.importorskip("anndata")
import scipy.io
import scipy.sparse as sp
import yaml

ROOT = Path(__file__).resolve().parents[1]
MARKERS = ROOT / "resources/markers/skin_v1.yaml"
PAIRS = ROOT / "resources/ligand_receptor/skin_core_v1.tsv"


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# truth label -> (level-1 parent, extra genes); level-2 markers come from the YAML itself
TYPES = {
    "keratinocyte_basal": 300, "keratinocyte_spinous": 300, "keratinocyte_granular": 150,
    "keratinocyte_inflammatory": 200, "fibroblast_SFRP2": 200, "fibroblast_inflammatory": 120,
    "vascular_endothelial_cell": 150, "T_CD4": 150, "T_CD8": 120, "T_IL17": 100,
    "macrophage": 120, "dendritic_cell_cDC2": 100, "mast_cell": 80, "melanocyte": 80,
}


def synthetic(seed=0):
    spec = yaml.safe_load(MARKERS.read_text())["level1"]
    parent, markers = {}, {}
    for p, s in spec.items():
        for st, ss in (s.get("subtypes") or {p: s}).items():
            parent[st] = p
            markers[st] = list(s["markers"]) + (list(ss["markers"]) if st != p else [])
    genes = sorted({g for m in markers.values() for g in m} | {"IL17RA", "IL17RC", "MT-CO1", "MT-ND1"}
                   | {f"GENE{i}" for i in range(400)})
    gi = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(seed)
    blocks, obs = [], []
    samples = [(f"D{d}_{c}", f"donor_{d}", cond) for d in range(3) for c, cond in (("PP", "lesional"), ("PN", "non_lesional"))]
    for sid, donor, cond in samples:
        for t, n in TYPES.items():
            n = int(n * (1.6 if cond == "lesional" and t in ("keratinocyte_inflammatory", "T_IL17") else
                         0.3 if cond == "non_lesional" and t in ("keratinocyte_inflammatory", "T_IL17") else 1) / 3)
            mu = np.full(len(genes), 0.3)
            mu[[gi[g] for g in markers[t]]] = 8.0
            mu[gi["MT-CO1"]] = mu[gi["MT-ND1"]] = 2.0
            if parent[t] == "keratinocyte":
                mu[[gi["IL17RA"], gi["IL17RC"]]] = 4.0
            blocks.append(rng.poisson(mu * rng.gamma(4, 0.25, (n, 1)) * rng.uniform(0.7, 1.3, len(genes))))
            obs += [{"sample_id": sid, "donor_id": donor, "condition": cond, "truth": t, "truth_l1": parent[t]}] * n
    X = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    o = pd.DataFrame(obs)
    o.index = [f"c{i}" for i in range(len(o))]
    return ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=genes))


@pytest.fixture(scope="module")
def annotated(tmp_path_factory):
    out = tmp_path_factory.mktemp("annot")
    raw = synthetic()
    raw.write_h5ad(out / "raw.h5ad")
    load("sc_annotate").run(out / "raw.h5ad", out, MARKERS, "sample_id", "donor_id", "condition",
                            max_mito=20, min_genes=50, nmads=5, doublets=False, n_hvg=500, n_pcs=20,
                            res1=1.0, res2=0.8, min_score=0.3, min_margin=0.1, seed=0)
    return out


def test_annotation_recovers_types(annotated):
    a = ad.read_h5ad(annotated / "annotated.h5ad")
    acc1 = (a.obs["cell_type_l1"].astype(str) == a.obs["truth_l1"].astype(str)).mean()
    acc2 = (a.obs["cell_type_l2"].astype(str) == a.obs["truth"].astype(str)).mean()
    assert acc1 > 0.95, acc1
    assert acc2 > 0.85, acc2
    assert (a.obs["cl_id"].astype(str).str.startswith("CL:") | a.obs["cell_type"].astype(str).str.contains("unassigned")).all()
    assert "counts" in a.layers and a.X.max() < 20


def test_frequency_tables(annotated):
    ps = pd.read_csv(annotated / "frequencies_l2_per_sample.tsv", sep="\t")
    assert np.allclose(ps.groupby("sample_id")["fraction"].sum(), 1)
    pc = pd.read_csv(annotated / "frequencies_l2_per_condition.tsv", sep="\t").set_index(["condition", "cell_type_l2"])
    infl = pc["mean_sample_fraction"]
    assert infl[("lesional", "keratinocyte_inflammatory")] > 2 * infl[("non_lesional", "keratinocyte_inflammatory")]
    assert (annotated / "annotation_report.md").read_text().startswith("# Annotation report")


def test_trajectory(annotated, tmp_path):
    load("sc_trajectory").main([str(annotated / "annotated.h5ad"), str(tmp_path),
                                "--lineage", "keratinocyte:keratinocyte_basal"])
    pt = pd.read_csv(tmp_path / "pseudotime_keratinocyte.tsv", sep="\t")
    assert pt["dpt_pseudotime"].between(0, 1).all()
    med = pt.groupby("cell_type_l2")["dpt_pseudotime"].median()
    assert med["keratinocyte_basal"] == med.min()
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["lineages"]["keratinocyte"]["root"] == "keratinocyte_basal"


def test_communication_builtin(annotated, tmp_path):
    load("sc_communication").main([str(annotated / "annotated.h5ad"), str(tmp_path), "--method", "builtin",
                                   "--pairs", str(PAIRS), "--n-perm", "50", "--max-pvalue", "1"])
    t = pd.read_csv(tmp_path / "communication_lesional.tsv", sep="\t")
    hit = t[(t.ligand == "IL17A") & (t.source == "T_IL17") & t.target.str.startswith("keratinocyte")]
    assert len(hit) and hit["lr_score"].gt(0).all()
    assert set(pd.read_csv(tmp_path / "communication_all.tsv", sep="\t")["condition"]) == {"lesional", "non_lesional"}


def test_loader_formats(tmp_path):
    dl = load("download_gse173706")
    rng = np.random.default_rng(1)
    genes = ["ACTB", "GAPDH", "KRT14", "KRT10", "COL1A1", "PTPRC"] + [f"G{i}" for i in range(20)]
    cells = [f"AAACCTGAGAAGGCCT-{i}" for i in range(15)]
    m = rng.poisson(2, (len(genes), len(cells)))
    csv = tmp_path / "GSM1_PP1.csv.gz"
    pd.DataFrame(m, index=genes, columns=cells).to_csv(csv)
    X, c, g = dl.read_counts(csv)
    assert X.shape == (15, len(genes)) and c == cells and g == genes and X.sum() == m.sum()
    tsv = tmp_path / "GSM2_PN1.tsv"
    pd.DataFrame(m.T, index=cells, columns=genes).to_csv(tsv, sep="\t")
    X2, c2, g2 = dl.read_counts(tsv)
    assert (X2 != X).nnz == 0 and g2 == genes
    scipy.io.mmwrite(tmp_path / "GSM3_matrix.mtx", sp.coo_matrix(m))
    with gzip.open(tmp_path / "GSM3_matrix.mtx.gz", "wb") as f:
        f.write((tmp_path / "GSM3_matrix.mtx").read_bytes())
    (tmp_path / "GSM3_matrix.mtx").unlink()
    pd.Series(cells).to_csv(tmp_path / "GSM3_barcodes.tsv.gz", sep="\t", header=False, index=False)
    pd.DataFrame({"id": genes, "name": genes, "t": "Gene Expression"}).to_csv(
        tmp_path / "GSM3_features.tsv.gz", sep="\t", header=False, index=False)
    X3, c3, g3 = dl.read_counts(tmp_path / "GSM3_matrix.mtx.gz")
    assert (X3 != X).nnz == 0 and c3 == cells and g3 == genes
