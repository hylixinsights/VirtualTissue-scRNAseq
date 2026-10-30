"""Pseudotime and state transitions within lineages of an annotated AnnData (plan task A3).

    python scripts/sc_trajectory.py data/psoriasis/annotated/annotated.h5ad data/psoriasis/trajectory \
        --lineage keratinocyte:keratinocyte_basal --lineage fibroblast:fibroblast_SFRP2

Each --lineage is PARENT:ROOT, where PARENT is a cell_type_l1 value and ROOT a cell_type_l2
value inside it. For each lineage the script recomputes neighbours on the integrated
representation, runs diffusion pseudotime (DPT) rooted at the root cell with the highest
root-marker score, and PAGA on the level-2 labels.

Outputs per lineage:
    pseudotime_<parent>.tsv    cell_id, cell_type_l2, condition, dpt_pseudotime
    transitions_<parent>.tsv   from, to, paga_connectivity, median pseudotime of each end,
                               plus per-condition connectivities; direction is from lower
                               to higher median pseudotime
and summary.json. The guidebook builder reads the transitions table as allowed state edges.
DPT assumes a continuous process; edges with low connectivity (< --min-connectivity) are
dropped and should be read as "not supported by these data", not as impossible.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def paga_edges(a, label, min_conn):
    import scanpy as sc
    if a.obs[label].nunique() < 2:
        return pd.DataFrame(columns=["a", "b", "paga_connectivity"])
    try:
        sc.tl.paga(a, groups=label)
    except ValueError:  # scanpy/igraph fail when no two states share a neighbour edge
        print(f"warning: no neighbour edges between {label} states; no transitions", file=sys.stderr)
        return pd.DataFrame(columns=["a", "b", "paga_connectivity"])
    cats = list(a.obs[label].cat.categories)
    conn = a.uns["paga"]["connectivities"].toarray()
    rows = [{"a": cats[i], "b": cats[j], "paga_connectivity": float(conn[i, j])}
            for i in range(len(cats)) for j in range(i + 1, len(cats)) if conn[i, j] >= min_conn]
    return pd.DataFrame(rows, columns=["a", "b", "paga_connectivity"])


def lineage(adata, parent, root, rep, condition_key, min_conn, seed):
    import scanpy as sc
    a = adata[adata.obs["cell_type_l1"].astype(str) == parent].copy()
    if a.n_obs < 50:
        raise ValueError(f"lineage {parent}: only {a.n_obs} cells")
    a.obs["cell_type_l2"] = a.obs["cell_type_l2"].astype(str).astype("category")
    roots = np.flatnonzero((a.obs["cell_type_l2"] == root).to_numpy())
    if not len(roots):
        raise ValueError(f"lineage {parent}: no cells labelled {root}; labels are {list(a.obs['cell_type_l2'].cat.categories)}")
    sc.pp.neighbors(a, use_rep=rep, n_neighbors=min(30, a.n_obs - 1), random_state=seed)
    sc.tl.diffmap(a, n_comps=min(15, a.n_obs - 2), random_state=seed)
    # root: the root-labelled cell at the extreme of the first diffusion component on the root side
    dc1 = a.obsm["X_diffmap"][:, 1]
    sign = 1 if np.median(dc1[roots]) >= np.median(dc1) else -1
    a.uns["iroot"] = int(roots[np.argmax(sign * dc1[roots])])
    sc.tl.dpt(a)
    pt = a.obs[["cell_type_l2", condition_key, "dpt_pseudotime"]].copy()
    med = pt.groupby("cell_type_l2", observed=True)["dpt_pseudotime"].median()

    edges = paga_edges(a, "cell_type_l2", min_conn)
    for cond in a.obs[condition_key].astype(str).unique():
        sub = a[a.obs[condition_key].astype(str) == cond].copy()
        sub.obs["cell_type_l2"] = sub.obs["cell_type_l2"].astype(str).astype("category")
        if sub.n_obs < 30:
            continue
        sc.pp.neighbors(sub, use_rep=rep, n_neighbors=min(30, sub.n_obs - 1), random_state=seed)
        e = paga_edges(sub, "cell_type_l2", 0.0).rename(columns={"paga_connectivity": f"connectivity_{cond}"})
        edges = edges.merge(e, on=["a", "b"], how="left")
    rows = []
    for r in edges.to_dict("records"):
        fwd = med[r["a"]] <= med[r["b"]]
        src, dst = (r["a"], r["b"]) if fwd else (r["b"], r["a"])
        rows.append({"from": src, "to": dst, "median_pseudotime_from": round(float(med[src]), 4),
                     "median_pseudotime_to": round(float(med[dst]), 4),
                     **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k not in ("a", "b")}})
    trans = pd.DataFrame(rows)
    by_cond = pt.groupby([condition_key, "cell_type_l2"], observed=True)["dpt_pseudotime"].median().unstack(0)
    return pt, trans, med, by_cond


def main(argv=None):
    import anndata as ad
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input")
    ap.add_argument("outdir")
    ap.add_argument("--lineage", action="append", required=True, help="PARENT_L1:ROOT_L2")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--rep", default=None, help="obsm key; default X_pca_harmony if present else X_pca")
    ap.add_argument("--min-connectivity", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    adata = ad.read_h5ad(a.input)
    rep = a.rep or ("X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    summary = {"input": str(a.input), "representation": rep, "method": "scanpy DPT + PAGA", "lineages": {}}
    for spec in a.lineage:
        parent, root = spec.split(":", 1)
        try:
            pt, trans, med, by_cond = lineage(adata, parent, root, rep, a.condition_key, a.min_connectivity, a.seed)
        except ValueError as e:  # one missing lineage should not cancel the others
            print(f"skipped: {e}", file=sys.stderr)
            summary["lineages"][parent] = {"root": root, "skipped": str(e)}
            continue
        pt.rename_axis("cell_id").to_csv(out / f"pseudotime_{parent}.tsv", sep="\t")
        trans.to_csv(out / f"transitions_{parent}.tsv", sep="\t", index=False)
        summary["lineages"][parent] = {
            "root": root, "n_cells": int(len(pt)),
            "median_pseudotime": {k: round(float(v), 4) for k, v in med.sort_values().items()},
            "median_pseudotime_by_condition": json.loads(by_cond.round(4).to_json()),
            "n_edges": int(len(trans))}
        print(f"{parent}: {len(pt)} cells, order {' -> '.join(med.sort_values().index)}")
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
