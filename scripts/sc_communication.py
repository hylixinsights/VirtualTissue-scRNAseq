"""Cell-cell communication per condition from an annotated AnnData (plan task A4).

    python scripts/sc_communication.py data/psoriasis/annotated/annotated.h5ad data/psoriasis/communication

Writes one TSV per condition, communication_<condition>.tsv, with columns
    source, target, ligand, receptor, pathway, lr_score, pvalue, ligand_frac, receptor_frac, method
plus communication_all.tsv (all conditions, with a condition column) and SOURCE.json.

Method (--method):
  liana     liana-py rank_aggregate with its consensus resource (needs `pip install liana`);
            the preferred method, and the one closest to CellChat in coverage.
  builtin   no extra dependency. For each ligand-receptor pair of a curated table
            (default resources/ligand_receptor/skin_core_v1.tsv; complexes written A_B use
            the minimum subunit) the score is mean ligand expression in the source type times
            mean receptor expression in the target type (log-normalised X). P-values come from
            --n-perm label permutations within the condition, as in CellPhoneDB. Pairs are
            kept when both sides are detected in at least --min-frac of the cells.
  auto      liana if installed, else builtin (default).

A CellChat run in R gives the same table through `subsetCommunication()`; save it as
communication_<condition>.tsv with these column names and the guidebook builder treats it
the same way.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def type_means(X, labels, genes_idx=None):
    """Return (mean expression, fraction expressing) per label for the given gene indices (all if None)."""
    import scipy.sparse as sp
    cats = pd.Categorical(labels)
    onehot = sp.csr_matrix((np.ones(len(cats)), (cats.codes, np.arange(len(cats)))),
                           shape=(len(cats.categories), len(cats)))
    sizes = np.asarray(onehot.sum(1)).ravel()[:, None]
    m = X if genes_idx is None else X[:, genes_idx]
    m = sp.csr_matrix(m) if not sp.issparse(m) else m.tocsr()
    mean = np.asarray((onehot @ m).todense()) / sizes
    frac = np.asarray((onehot @ (m > 0).astype(np.float32)).todense()) / sizes
    return list(cats.categories), mean, frac


def builtin(adata, label, pairs, min_frac, n_perm, min_cells, seed):
    genes = {g.upper(): i for i, g in enumerate(adata.var_names)}
    rows = []
    for _, p in pairs.iterrows():
        lig, rec = p["ligand"].split("_"), p["receptor"].split("_")
        if not all(g.upper() in genes for g in lig + rec):
            continue
        rows.append((p, [genes[g.upper()] for g in lig], [genes[g.upper()] for g in rec]))
    if not rows:
        return pd.DataFrame()
    counts = adata.obs[label].value_counts()
    keep = adata.obs[label].isin(counts[counts >= min_cells].index).to_numpy()
    a = adata[keep]
    labels = a.obs[label].astype(str).to_numpy()
    all_idx = sorted({i for _, l, r in rows for i in l + r})
    pos = {g: k for k, g in enumerate(all_idx)}
    sub = a.X[:, all_idx]  # slice once; permutations reuse it
    cats, mean, frac = type_means(sub, labels)

    def complex_stat(arr, idx):
        return arr[:, [pos[i] for i in idx]].min(1)

    obs = np.stack([np.outer(complex_stat(mean, l), complex_stat(mean, r)) for _, l, r in rows])
    rng = np.random.default_rng(seed)
    exceed = np.zeros_like(obs)
    for _ in range(n_perm):
        _, pm, _ = type_means(sub, rng.permutation(labels))
        perm = np.stack([np.outer(complex_stat(pm, l), complex_stat(pm, r)) for _, l, r in rows])
        exceed += perm >= obs
    pval = (exceed + 1) / (n_perm + 1)
    out = []
    for k, (p, l, r) in enumerate(rows):
        lf, rf = complex_stat(frac, l), complex_stat(frac, r)
        for i, src in enumerate(cats):
            for j, dst in enumerate(cats):
                if lf[i] < min_frac or rf[j] < min_frac or obs[k, i, j] <= 0:
                    continue
                out.append({"source": src, "target": dst, "ligand": p["ligand"], "receptor": p["receptor"],
                            "pathway": p["pathway"], "lr_score": float(obs[k, i, j]), "pvalue": float(pval[k, i, j]),
                            "ligand_frac": float(lf[i]), "receptor_frac": float(rf[j]), "method": "builtin"})
    return pd.DataFrame(out)


def with_liana(adata, label, min_frac, seed):
    import liana as li
    li.mt.rank_aggregate(adata, groupby=label, resource_name="consensus", expr_prop=min_frac,
                         use_raw=False, verbose=False, seed=seed)
    r = adata.uns["liana_res"].rename(columns={"ligand_complex": "ligand", "receptor_complex": "receptor"})
    return pd.DataFrame({"source": r["source"], "target": r["target"], "ligand": r["ligand"],
                         "receptor": r["receptor"], "pathway": "", "lr_score": r["magnitude_rank"].rsub(1),
                         "pvalue": r["specificity_rank"], "ligand_frac": r.get("ligand_props", np.nan),
                         "receptor_frac": r.get("receptor_props", np.nan), "method": "liana_rank_aggregate"})


def main(argv=None):
    import anndata as ad
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input")
    ap.add_argument("outdir")
    ap.add_argument("--label", default="cell_type")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--method", choices=["auto", "liana", "builtin"], default="auto")
    ap.add_argument("--pairs", default="resources/ligand_receptor/skin_core_v1.tsv")
    ap.add_argument("--min-frac", type=float, default=0.1)
    ap.add_argument("--min-cells", type=int, default=20, help="drop types with fewer cells in a condition")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--max-pvalue", type=float, default=0.05, help="filter for the written tables (1 keeps all)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    method = a.method
    if method == "auto":
        try:
            import liana  # noqa: F401
            method = "liana"
        except ImportError:
            method = "builtin"
    adata = ad.read_h5ad(a.input)
    if adata.X.max() > 50:
        raise SystemExit("X looks like raw counts; pass the annotated.h5ad (log-normalised X)")
    pairs = pd.read_csv(a.pairs, sep="\t", dtype=str).fillna("")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = []
    for cond in sorted(adata.obs[a.condition_key].astype(str).unique()):
        sub = adata[adata.obs[a.condition_key].astype(str) == cond].copy()
        if method == "liana":
            counts = sub.obs[a.label].value_counts()
            sub = sub[sub.obs[a.label].isin(counts[counts >= a.min_cells].index)].copy()
            t = with_liana(sub, a.label, a.min_frac, a.seed)
        else:
            t = builtin(sub, a.label, pairs, a.min_frac, a.n_perm, a.min_cells, a.seed)
        if len(t):
            t = t[t["pvalue"] <= a.max_pvalue].sort_values(["pvalue", "lr_score"], ascending=[True, False])
        t.to_csv(out / f"communication_{cond}.tsv", sep="\t", index=False)
        print(f"{cond}: {len(t)} interactions ({method})")
        tables.append(t.assign(condition=cond))
    allt = pd.concat(tables) if tables else pd.DataFrame()
    allt.to_csv(out / "communication_all.tsv", sep="\t", index=False)
    (out / "SOURCE.json").write_text(json.dumps({
        "input": str(a.input), "method": method, "label": a.label, "pairs": a.pairs if method == "builtin" else "liana consensus",
        "min_frac": a.min_frac, "n_perm": a.n_perm, "max_pvalue": a.max_pvalue, "seed": a.seed}, indent=1) + "\n")


if __name__ == "__main__":
    main()
