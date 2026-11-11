"""Step 3b: intervals for the recall@k comparison.

Repeats the phylogeny-grouped cross-validation with permuted group labels
(GroupKFold is deterministic, so permuting the labels changes the partition)
and reports mean, spread, and paired tests against the prevalence baseline.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

DATA = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/data"
THR = int(sys.argv[2]) if len(sys.argv) > 2 else 1
N_REP = int(sys.argv[3]) if len(sys.argv) > 3 else 5
KS = [1, 2, 3, 5, 10]

mat = pd.read_csv(f"{DATA}/interaction_matrix.csv", sep=";").set_index("bacteria")
host = pd.read_csv(f"{DATA}/picard_collection.csv", sep=";").set_index("bacteria")
dist = pd.read_csv(f"{DATA}/phylo_dist.tsv", sep="\t").set_index("bacteria")
ds = pd.read_csv(f"{DATA}/defense_subtypes.csv", sep=";").set_index("bacteria")

B = (mat >= THR).astype(float).where(mat.notna())
strains = [s for s in B.index if s in dist.index and s in dist.columns]
B = B.loc[strains]
H = host.reindex(strains)
D = dist.loc[strains, strains].astype(float).to_numpy()
np.fill_diagonal(D, 0.0)
D = (D + D.T) / 2
Y = B.to_numpy()

PROXY = ["O-type", "ST_Warwick", "H-type", "Clermont_Phylo"]
MECH = ["ABC_serotype", "LPS_type"]
DEF = (ds.reindex(strains) > 0).astype(float)
DEF = DEF.loc[:, DEF.sum() >= 10]

Z = linkage(squareform(D, checks=False), method="single")
base_groups = fcluster(Z, t=1e-4, criterion="distance")
n_g = len(set(base_groups))
print(f"strains {len(strains)}  phages {B.shape[1]}  groups {n_g}  "
      f"threshold score>={THR}  repeats {N_REP}")


def design(cols, idx_fit):
    parts = []
    for c in cols:
        v = H[c].astype("object").astype(str).fillna("__missing__")
        keep = v.iloc[idx_fit].value_counts()
        keep = set(keep[keep >= 3].index)
        v = v.where(v.isin(keep), "__other__")
        parts.append(pd.get_dummies(v, prefix=c))
    return pd.concat(parts, axis=1).astype(float).to_numpy()


def scores_for(kind, cols, extra, groups):
    S = np.full(B.shape, np.nan)
    gkf = GroupKFold(n_splits=10)
    idx = np.arange(len(strains))
    for tr, te in gkf.split(idx, groups=groups):
        if kind == "prev":
            for j in range(B.shape[1]):
                ytr = Y[tr, j][~np.isnan(Y[tr, j])]
                S[te, j] = ytr.mean() if len(ytr) else 0.0
            continue
        if kind == "nn":
            prev = np.nan_to_num(np.array(
                [np.nanmean(Y[tr, j]) for j in range(B.shape[1])]))
            for i in te:
                nn = tr[np.argsort(D[i, tr])[:5]]
                prof = np.nanmean(Y[nn], axis=0)
                S[i, :] = np.where(np.isnan(prof), prev, prof) + 1e-6 * prev
            continue
        X = design(cols, tr)
        if extra is not None:
            X = np.hstack([X, extra.to_numpy()])
        for j in range(B.shape[1]):
            y = Y[:, j]
            m_tr = tr[~np.isnan(y[tr])]
            ytr = y[m_tr]
            if len(ytr) == 0 or ytr.sum() < 2 or ytr.sum() == len(ytr):
                S[te, j] = ytr.mean() if len(ytr) else 0.0
                continue
            clf = LogisticRegression(max_iter=2000, C=1.0)
            clf.fit(X[m_tr], ytr)
            S[te, j] = clf.predict_proba(X[te])[:, 1]
    return S


PHAGE_RANK = np.argsort(np.argsort(np.array(list(B.columns), dtype=object)))


def _order(scores, mask):
    """score descending, ties broken by phage identifier ascending (declared)."""
    return np.lexsort((PHAGE_RANK[mask], -scores[mask]))


def per_strain_recall(S):
    """returns a per-strain recall@k array for each k (for paired tests)"""
    out = {k: [] for k in KS}
    for i in range(Y.shape[0]):
        y, s = Y[i], S[i]
        m = ~np.isnan(y) & ~np.isnan(s)
        tot = y[m].sum()
        if tot < 1:
            continue
        yy = y[m][_order(s, m)]
        for k in KS:
            out[k].append(yy[:k].sum() / tot)
    return {k: np.array(v) for k, v in out.items()}


def oracle():
    out = {k: [] for k in KS}
    for i in range(Y.shape[0]):
        y = Y[i]
        tot = y[~np.isnan(y)].sum()
        if tot < 1:
            continue
        for k in KS:
            out[k].append(min(k, tot) / tot)
    return {k: float(np.mean(v)) for k, v in out.items()}


ORC = oracle()
CONFIGS = {
    "prevalence":   ("prev", None, None),
    "nn_phylo":     ("nn", None, None),
    "proxy":        ("lr", PROXY, None),
    "mechanism":    ("lr", MECH, None),
    "all_annot":    ("lr", PROXY + MECH, None),
    "mech+defense": ("lr", MECH, DEF),
    "all+defense":  ("lr", PROXY + MECH, DEF),
}

rng = np.random.default_rng(0)
acc = {n: {k: [] for k in KS} for n in CONFIGS}
paired = {n: {k: [] for k in KS} for n in CONFIGS}

for rep in range(N_REP):
    if rep == 0:
        g = base_groups
    else:
        perm = rng.permutation(n_g) + 1
        remap = {old: new for old, new in zip(sorted(set(base_groups)), perm)}
        g = np.array([remap[x] for x in base_groups])
    per_strain = {}
    for name, (kind, cols, extra) in CONFIGS.items():
        S = scores_for(kind, cols, extra, g)
        r = per_strain_recall(S)
        per_strain[name] = r
        for k in KS:
            acc[name][k].append(float(r[k].mean()))
    for name in CONFIGS:
        for k in KS:
            paired[name][k].append(per_strain[name][k] - per_strain["prevalence"][k])
    print(f"  repeat {rep+1}/{N_REP} done", flush=True)

print(f"\noracle recall@k: " + "  ".join(f"k={k}:{ORC[k]:.3f}" for k in KS))

print("\nmean recall@k over repeats, with min-max across the 5 partitions")
hdr = f"{'feature set':<15}" + "".join(f"{'k='+str(k):>18}" for k in KS)
print(hdr)
print("-" * len(hdr))
for name in CONFIGS:
    cells = []
    for k in KS:
        v = np.array(acc[name][k])
        cells.append(f"{v.mean():.3f} [{v.min():.3f},{v.max():.3f}]")
    print(f"{name:<15}" + "".join(f"{c:>18}" for c in cells))

print("\nfraction of oracle (mean over repeats)")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>9}" for k in KS))
for name in CONFIGS:
    print(f"{name:<15}" + "".join(
        f"{np.mean(acc[name][k])/ORC[k]:>9.3f}" for k in KS))

print("\ngain over prevalence in percentage points, paired across strains")
print("  (Wilcoxon on per-strain differences, pooled over repeats)")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>22}" for k in KS))
for name in CONFIGS:
    if name == "prevalence":
        continue
    cells = []
    for k in KS:
        d = np.concatenate(paired[name][k])
        nz = d[d != 0]
        p = stats.wilcoxon(nz).pvalue if len(nz) > 10 else np.nan
        cells.append(f"{100*d.mean():+.1f}  p={p:.1g}")
    print(f"{name:<15}" + "".join(f"{c:>22}" for c in cells))

print("\ndoes adding defense systems help or hurt? paired per-strain difference")
print("  (positive = adding defense systems improves recall)")
for base, withdef in [("mechanism", "mech+defense"), ("all_annot", "all+defense")]:
    for k in [2, 5]:
        d = np.concatenate(paired[withdef][k]) - np.concatenate(paired[base][k])
        nz = d[d != 0]
        p = stats.wilcoxon(nz).pvalue if len(nz) > 10 else np.nan
        print(f"  {withdef} - {base:<12} k={k}: {100*d.mean():+.2f} pp   p={p:.2g}")
