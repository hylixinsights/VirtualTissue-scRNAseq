"""QC, integrate, cluster and marker-annotate a raw-count AnnData; write cell-type frequencies.

    python scripts/sc_annotate.py data/psoriasis/GSE173706_raw.h5ad data/psoriasis/annotated \
        --markers resources/markers/skin_v1.yaml --batch-key donor_id --condition-key condition

Input: one AnnData with raw counts in X (or layers["counts"]) and obs columns for sample,
batch and condition (scripts/download_gse173706.py writes such a file). Output directory:

    annotated.h5ad            counts in layers["counts"], log-normalised X, obs labels
                              cell_type_l1, cell_type_l2, cell_type (= l2), cl_id, qc columns
    qc_per_sample.tsv         cells before/after each filter, per sample
    clusters_l1.tsv           Leiden cluster -> level-1 label with every marker score
    clusters_l2.tsv           sub-cluster -> level-2 label, per level-1 type
    frequencies_<level>_per_sample.tsv     counts and fractions per sample
    frequencies_<level>_per_condition.tsv  pooled fraction and mean/sd of per-sample fractions
    markers_detected.tsv      top genes per final label (Wilcoxon), for manual review
    annotation_report.md      the above in one readable page

Labels are assigned per cluster from the marker YAML (see resources/markers/skin_v1.yaml),
not per cell, so a single cell's noise cannot flip its label. A cluster whose best score is
below --min-score, or whose best and second-best scores are within --min-margin, is kept as
"unassigned_<cluster>" for manual review instead of being forced into a type.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def load_markers(path):
    spec = yaml.safe_load(Path(path).read_text())
    return spec


def mad_outlier(x, nmads):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    return (x < med - nmads * mad) | (x > med + nmads * mad)


def qc(adata, sample_key, max_mito, min_genes, nmads, doublets, seed):
    import scanpy as sc
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True)
    rows = []
    keep = np.ones(adata.n_obs, bool)
    for s, idx in adata.obs.groupby(sample_key, observed=True).indices.items():
        o = adata.obs.iloc[idx]
        low_genes = o["n_genes_by_counts"].to_numpy() < min_genes
        outlier = (mad_outlier(o["log1p_total_counts"].to_numpy(), nmads)
                   | mad_outlier(o["log1p_n_genes_by_counts"].to_numpy(), nmads))
        mito = o["pct_counts_mt"].to_numpy() > max_mito
        k = ~(low_genes | outlier | mito)
        keep[idx] = k
        rows.append({"sample": s, "cells_in": len(idx), "low_genes": int(low_genes.sum()),
                     "count_outlier": int((outlier & ~low_genes).sum()), "high_mito": int((mito & ~low_genes & ~outlier).sum()),
                     "median_genes": float(np.median(o["n_genes_by_counts"])),
                     "median_umis": float(np.median(o["total_counts"])),
                     "median_pct_mito": float(np.median(o["pct_counts_mt"]))})
    adata = adata[keep].copy()
    table = pd.DataFrame(rows).set_index("sample")
    table["doublets"] = 0
    if doublets:
        sc.pp.scrublet(adata, batch_key=sample_key, random_state=seed)
        d = adata.obs["predicted_doublet"].astype(bool)
        table["doublets"] = d.groupby(adata.obs[sample_key], observed=True).sum().reindex(table.index).fillna(0).astype(int)
        adata = adata[~d.to_numpy()].copy()
    table["cells_out"] = adata.obs[sample_key].value_counts().reindex(table.index).fillna(0).astype(int)
    sc.pp.filter_genes(adata, min_cells=3)
    return adata, table


def embed(adata, batch_key, n_hvg, n_pcs, seed):
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    batch = batch_key if batch_key and adata.obs[batch_key].nunique() > 1 else None
    sc.pp.highly_variable_genes(adata, n_top_genes=min(n_hvg, adata.n_vars - 1), flavor="seurat_v3",
                                layer="counts", batch_key=batch, subset=False)
    n_pcs = min(n_pcs, adata.n_obs - 1, int(adata.var["highly_variable"].sum()) - 1)
    sc.pp.pca(adata, n_comps=n_pcs, mask_var="highly_variable", random_state=seed)
    rep = "X_pca"
    if batch:
        try:
            import harmonypy
            z = np.asarray(harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [batch],
                                                 random_state=seed, verbose=False).Z_corr)
            # harmonypy < 0.1 returns PCs x cells, newer versions cells x PCs
            adata.obsm["X_pca_harmony"] = z if z.shape[0] == adata.n_obs else z.T
            rep = "X_pca_harmony"
        except ImportError:
            print("warning: harmonypy not installed; clustering on uncorrected PCA", file=sys.stderr)
    return rep


def neighbours_and_cluster(adata, rep, resolution, seed, key):
    import scanpy as sc
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=min(15, adata.n_obs - 1), random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, key_added=key, random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)


def score_types(adata, types, prefix):
    """Add one obs column per type: mean z-scored expression of markers minus half of negatives."""
    import scipy.sparse as sp
    genes = {g.upper(): i for i, g in enumerate(adata.var_names)}
    used = {}
    X = adata.X
    for name, spec in types.items():
        pos = [genes[g.upper()] for g in spec["markers"] if g.upper() in genes]
        neg = [genes[g.upper()] for g in spec.get("negative", []) if g.upper() in genes]
        used[name] = [adata.var_names[i] for i in pos]
        score = np.zeros(adata.n_obs)
        for idx, w in ((pos, 1.0), (neg, -0.5)):
            if not idx:
                continue
            m = X[:, idx]
            m = m.toarray() if sp.issparse(m) else np.asarray(m)
            sd = m.std(0)
            sd[sd == 0] = 1
            score += w * ((m - m.mean(0)) / sd).mean(1)
        adata.obs[f"{prefix}{name}"] = score if pos else np.nan
    return used


def label_clusters(adata, cluster_key, types, prefix, min_score, min_margin):
    cols = [f"{prefix}{t}" for t in types]
    med = adata.obs.groupby(cluster_key, observed=True)[cols].median()
    med.columns = list(types)
    rows = []
    for cl, r in med.iterrows():
        r = r.dropna().sort_values(ascending=False)
        best, second = (r.index[0], r.iloc[0]), (r.index[1], r.iloc[1]) if len(r) > 1 else (None, -np.inf)
        ok = best[1] >= min_score and best[1] - second[1] >= min_margin
        rows.append({"cluster": cl, "label": best[0] if ok else f"unassigned_{cl}", "best": best[0],
                     "best_score": round(float(best[1]), 3), "second": second[0],
                     "margin": round(float(best[1] - second[1]), 3),
                     "n_cells": int((adata.obs[cluster_key] == cl).sum()), **med.loc[cl].round(3).to_dict()})
    table = pd.DataFrame(rows).set_index("cluster")
    return adata.obs[cluster_key].map(table["label"]).astype(str), table


def frequencies(obs, label, sample_key, condition_key):
    counts = obs.groupby([sample_key, label], observed=True).size().unstack(fill_value=0)
    frac = counts.div(counts.sum(axis=1), axis=0)
    cond = obs.groupby(sample_key, observed=True)[condition_key].first()
    per_sample = (counts.stack().rename("n_cells").to_frame()
                  .join(frac.stack().rename("fraction")).reset_index())
    per_sample.insert(1, condition_key, per_sample[sample_key].map(cond))
    pooled = obs.groupby([condition_key, label], observed=True).size().unstack(fill_value=0)
    pooled_frac = pooled.div(pooled.sum(axis=1), axis=0)
    f = frac.join(cond)
    mean, sd, n = (f.groupby(condition_key, observed=True).agg(a) for a in ("mean", "std", "count"))
    per_condition = pd.concat({"n_cells": pooled.stack(), "pooled_fraction": pooled_frac.stack(),
                               "mean_sample_fraction": mean.stack(), "sd_sample_fraction": sd.stack(),
                               "n_samples": n.stack()}, axis=1)
    per_condition.index.names = [condition_key, label]
    per_condition = per_condition.reset_index()
    return per_sample, per_condition


def md_table(df, floatfmt=3):
    df = df.copy()
    for c in df.columns:
        if pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].round(floatfmt)
    head = "| " + " | ".join(map(str, [df.index.name or ""] + list(df.columns))) + " |"
    sep = "|" + "---|" * (len(df.columns) + 1)
    body = ["| " + " | ".join(map(str, [i] + list(r))) + " |" for i, r in zip(df.index, df.itertuples(index=False))]
    return "\n".join([head, sep] + body)


def run(inp, outdir, markers, sample_key, batch_key, condition_key, max_mito, min_genes, nmads,
        doublets, n_hvg, n_pcs, res1, res2, min_score, min_margin, seed):
    import anndata as ad
    import scanpy as sc
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    spec = load_markers(markers)
    adata = ad.read_h5ad(inp)
    for k in (sample_key, condition_key) + ((batch_key,) if batch_key else ()):
        if k not in adata.obs:
            sys.exit(f"obs column '{k}' missing; available: {list(adata.obs.columns)}")
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    adata.layers["counts"] = adata.X.copy()

    adata, qc_table = qc(adata, sample_key, max_mito, min_genes, nmads, doublets, seed)
    qc_table.to_csv(outdir / "qc_per_sample.tsv", sep="\t")
    print(f"QC kept {adata.n_obs} cells x {adata.n_vars} genes")

    rep = embed(adata, batch_key, n_hvg, n_pcs, seed)
    neighbours_and_cluster(adata, rep, res1, seed, "leiden_l1")
    sc.tl.umap(adata, random_state=seed)
    l1 = spec["level1"]
    used = {"level1": score_types(adata, l1, "score_l1_")}
    adata.obs["cell_type_l1"], t1 = label_clusters(adata, "leiden_l1", l1, "score_l1_", min_score, min_margin)
    t1.to_csv(outdir / "clusters_l1.tsv", sep="\t")

    adata.obs["cell_type_l2"] = adata.obs["cell_type_l1"].astype(str)
    adata.obs["leiden_l2"] = adata.obs["leiden_l1"].astype(str)
    t2_all = []
    for parent, pspec in l1.items():
        subtypes = pspec.get("subtypes")
        mask = (adata.obs["cell_type_l1"] == parent).to_numpy()
        if not subtypes or mask.sum() < 30:
            continue
        sub = adata[mask].copy()
        neighbours_and_cluster(sub, rep, res2, seed, "leiden_sub")
        used[parent] = score_types(sub, subtypes, "score_l2_")
        labels, t2 = label_clusters(sub, "leiden_sub", subtypes, "score_l2_", min_score, min_margin)
        labels = labels.where(~labels.str.startswith("unassigned_"), parent + "_" + labels)
        adata.obs.loc[mask, "cell_type_l2"] = labels.to_numpy()
        adata.obs.loc[mask, "leiden_l2"] = (parent + ":" + sub.obs["leiden_sub"].astype(str)).to_numpy()
        t2.insert(0, "parent", parent)
        t2_all.append(t2.reset_index())
    t2 = pd.concat(t2_all) if t2_all else pd.DataFrame()
    t2.to_csv(outdir / "clusters_l2.tsv", sep="\t", index=False)

    cl = {t: s["cl_id"] for t, s in l1.items()}
    cl.update({st: ss["cl_id"] for s in l1.values() for st, ss in (s.get("subtypes") or {}).items()})
    for col in ("cell_type_l1", "cell_type_l2"):
        adata.obs[col] = adata.obs[col].astype("category")
    adata.obs["cell_type"] = adata.obs["cell_type_l2"]
    adata.obs["cl_id"] = adata.obs["cell_type_l2"].astype(str).map(cl).fillna("").astype("category")
    adata.obs = adata.obs.drop(columns=[c for c in adata.obs if c.startswith("score_l2_")])

    for level in ("cell_type_l1", "cell_type_l2"):
        ps, pc = frequencies(adata.obs, level, sample_key, condition_key)
        ps.to_csv(outdir / f"frequencies_{level[-2:]}_per_sample.tsv", sep="\t", index=False)
        pc.to_csv(outdir / f"frequencies_{level[-2:]}_per_condition.tsv", sep="\t", index=False)

    groups = adata.obs["cell_type"].value_counts()
    valid = groups[groups >= 10].index.tolist()
    top = pd.DataFrame()
    if len(valid) > 1:
        sc.tl.rank_genes_groups(adata, "cell_type", groups=valid, method="wilcoxon", pts=True)
        top = (sc.get.rank_genes_groups_df(adata, None).groupby("group", observed=True).head(15))
        top.to_csv(outdir / "markers_detected.tsv", sep="\t", index=False)
        del adata.uns["rank_genes_groups"]

    adata.uns["annotation"] = json.loads(json.dumps({
        "markers": spec.get("version", str(markers)), "markers_used": used,
        "params": dict(sample_key=sample_key, batch_key=batch_key, condition_key=condition_key,
                       max_mito=max_mito, min_genes=min_genes, nmads=nmads, doublets=doublets,
                       n_hvg=n_hvg, n_pcs=n_pcs, res_l1=res1, res_l2=res2, min_score=min_score,
                       min_margin=min_margin, seed=seed, representation=rep),
        "method": "cluster-level marker scoring (z-scored mean), two levels"}))
    adata.write_h5ad(outdir / "annotated.h5ad", compression="gzip")
    write_report(outdir, adata, qc_table, t1, t2, condition_key, sample_key, top)
    print(f"wrote {outdir / 'annotated.h5ad'}")
    return adata


def write_report(outdir, adata, qc_table, t1, t2, condition_key, sample_key, top):
    pc1 = pd.read_csv(outdir / "frequencies_l1_per_condition.tsv", sep="\t")
    pc2 = pd.read_csv(outdir / "frequencies_l2_per_condition.tsv", sep="\t")
    wide = lambda pc, lab: pc.pivot(index=lab, columns=condition_key, values="mean_sample_fraction").fillna(0)
    unassigned = adata.obs["cell_type"].astype(str).str.contains("unassigned").mean()
    p = adata.uns["annotation"]["params"]
    lines = [
        "# Annotation report", "",
        f"Cells after QC: **{adata.n_obs}** from {adata.obs[sample_key].nunique()} samples; "
        f"genes: {adata.n_vars}. Unassigned after annotation: {unassigned:.1%}.", "",
        f"Parameters: `{json.dumps(p)}`", "",
        "Labels come from marker scores per cluster. Treat this as a first pass: check the "
        "unassigned clusters, the margin column and the detected markers before using labels "
        "in a guidebook.", "",
        "## QC per sample", "", md_table(qc_table), "",
        "## Level-1 clusters", "", md_table(t1[["label", "n_cells", "best_score", "second", "margin"]]), "",
    ]
    if len(t2):
        lines += ["## Level-2 sub-clusters", "",
                  md_table(t2.set_index("cluster")[["parent", "label", "n_cells", "best_score", "second", "margin"]]), ""]
    lines += ["## Mean per-sample fraction by condition, level 1", "", md_table(wide(pc1, "cell_type_l1")), "",
              "## Mean per-sample fraction by condition, level 2", "", md_table(wide(pc2, "cell_type_l2")), ""]
    if len(top):
        best = top.groupby("group", observed=True)["names"].apply(lambda s: ", ".join(s.head(8)))
        lines += ["## Top detected genes per label (Wilcoxon)", "", md_table(best.to_frame("top genes")), ""]
    (outdir / "annotation_report.md").write_text("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input")
    ap.add_argument("outdir")
    ap.add_argument("--markers", default="resources/markers/skin_v1.yaml")
    ap.add_argument("--sample-key", default="sample_id")
    ap.add_argument("--batch-key", default="donor_id", help="'' to skip integration")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--max-mito", type=float, default=20.0, help="max %% mitochondrial counts")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--nmads", type=float, default=5.0, help="MAD cut-off for counts/genes outliers per sample")
    ap.add_argument("--no-doublets", action="store_true", help="skip Scrublet")
    ap.add_argument("--n-hvg", type=int, default=3000)
    ap.add_argument("--n-pcs", type=int, default=40)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--resolution-sub", type=float, default=0.8)
    ap.add_argument("--min-score", type=float, default=0.3)
    ap.add_argument("--min-margin", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    run(a.input, a.outdir, a.markers, a.sample_key, a.batch_key or None, a.condition_key, a.max_mito,
        a.min_genes, a.nmads, not a.no_doublets, a.n_hvg, a.n_pcs, a.resolution, a.resolution_sub,
        a.min_score, a.min_margin, a.seed)


if __name__ == "__main__":
    main()
