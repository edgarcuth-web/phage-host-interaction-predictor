"""Step 2c: recall at a fixed experimental budget.

Clinical framing: a patient isolate arrives, you can afford to plaque-test k
phages. How many of the phages that would actually have worked do you test?

Feature sets compared under identical phylogeny-grouped cross-validation:
  prevalence   training positive rate per phage, no host information
  nn_phylo     susceptibility profile of the nearest training strain
  proxy        annotations that are largely a relabelling of phylogeny
  mechanism    annotations that are not (ABC capsule, LPS outer core)
  all          both annotation groups
  mech+defense mechanism plus the DefenseFinder repertoire
"""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

DATA = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/data"
THR = int(sys.argv[2]) if len(sys.argv) > 2 else 1
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

PROXY = ["O-type", "ST_Warwick", "H-type", "Clermont_Phylo"]
MECH = ["ABC_serotype", "LPS_type"]
DEF = (ds.reindex(strains) > 0).astype(float)
DEF = DEF.loc[:, DEF.sum() >= 10]          # systems present in >= 10 strains

# phylogeny-grouped folds, single linkage at 1e-4 as in the paper
Z = linkage(squareform(D, checks=False), method="single")
groups = fcluster(Z, t=1e-4, criterion="distance")
print(f"strains {len(strains)}   phages {B.shape[1]}   "
      f"phylogenetic groups at 1e-4: {len(set(groups))}")


def design(cols, idx_fit, idx_all):
    """one-hot the categorical columns, rare levels (fit set) collapsed"""
    parts = []
    for c in cols:
        v = H[c].astype("object").astype(str).fillna("__missing__")
        keep = v.iloc[idx_fit].value_counts()
        keep = set(keep[keep >= 3].index)
        v = v.where(v.isin(keep), "__other__")
        parts.append(pd.get_dummies(v, prefix=c))
    X = pd.concat(parts, axis=1).astype(float)
    return X.to_numpy()


def fit_predict(cols, extra=None):
    """per-phage logistic regression; returns a strain x phage score matrix"""
    S = np.full(B.shape, np.nan)
    gkf = GroupKFold(n_splits=10)
    for tr, te in gkf.split(np.arange(len(strains)), groups=groups):
        X = design(cols, tr, None) if cols else np.zeros((len(strains), 0))
        if extra is not None:
            X = np.hstack([X, extra.to_numpy()])
        for j, p in enumerate(B.columns):
            y = B.iloc[:, j].to_numpy()
            m_tr = tr[~np.isnan(y[tr])]
            m_te = te
            ytr = y[m_tr]
            if ytr.sum() < 2 or ytr.sum() == len(ytr):
                S[m_te, j] = np.nanmean(ytr) if len(ytr) else 0.0
                continue
            if X.shape[1] == 0:
                S[m_te, j] = ytr.mean()
                continue
            clf = LogisticRegression(max_iter=2000, C=1.0)
            clf.fit(X[m_tr], ytr)
            S[m_te, j] = clf.predict_proba(X[m_te])[:, 1]
    return S


def prevalence_scores():
    S = np.full(B.shape, np.nan)
    gkf = GroupKFold(n_splits=10)
    for tr, te in gkf.split(np.arange(len(strains)), groups=groups):
        for j in range(B.shape[1]):
            y = B.iloc[:, j].to_numpy()
            ytr = y[tr][~np.isnan(y[tr])]
            S[te, j] = ytr.mean() if len(ytr) else 0.0
    return S


def nn_phylo_scores():
    """nearest training strain's profile, ties broken by training prevalence"""
    S = np.full(B.shape, np.nan)
    gkf = GroupKFold(n_splits=10)
    for tr, te in gkf.split(np.arange(len(strains)), groups=groups):
        prev = np.array([np.nanmean(B.iloc[tr, j].to_numpy()) for j in range(B.shape[1])])
        prev = np.nan_to_num(prev)
        for i in te:
            d = D[i, tr]
            nn = tr[np.argsort(d)[:5]]
            prof = np.nanmean(B.iloc[nn].to_numpy(), axis=0)
            prof = np.where(np.isnan(prof), prev, prof)
            S[i, :] = prof + 1e-6 * prev
    return S


def recall_at_k(S):
    """per-strain: rank the 96 phages, recall of that strain's true positives"""
    Y = B.to_numpy()
    out = {k: [] for k in KS}
    n_used = 0
    for i in range(Y.shape[0]):
        y, s = Y[i], S[i]
        m = ~np.isnan(y) & ~np.isnan(s)
        if y[m].sum() < 1:
            continue
        n_used += 1
        order = np.argsort(-s[m])
        yy = y[m][order]
        tot = yy.sum()
        for k in KS:
            out[k].append(yy[:k].sum() / tot)
    return {k: float(np.mean(v)) for k, v in out.items()}, n_used


configs = {
    "prevalence":   ("prev", None),
    "nn_phylo":     ("nn", None),
    "proxy":        (PROXY, None),
    "mechanism":    (MECH, None),
    "all_annot":    (PROXY + MECH, None),
    "mech+defense": (MECH, DEF),
    "all+defense":  (PROXY + MECH, DEF),
}

print(f"\nrecall at k, averaged over strains with >= 1 susceptible phage (score >= {THR})")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>9}" for k in KS))
print("-" * (15 + 9 * len(KS)))
results = {}
for name, (cols, extra) in configs.items():
    if cols == "prev":
        S = prevalence_scores()
    elif cols == "nn":
        S = nn_phylo_scores()
    else:
        S = fit_predict(cols, extra)
    r, n_used = recall_at_k(S)
    results[name] = r
    print(f"{name:<15}" + "".join(f"{r[k]:>9.3f}" for k in KS))
print(f"\nstrains evaluated: {n_used} of {len(strains)}")

print("\ngain over the prevalence baseline (percentage points of recall):")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>9}" for k in KS))
for name in configs:
    if name == "prevalence":
        continue
    print(f"{name:<15}" + "".join(
        f"{100*(results[name][k]-results['prevalence'][k]):>9.1f}" for k in KS))

print("\nnote: per-phage models give per-phage probabilities, so ranking phages")
print("for one strain compares uncalibrated scores. The flaw is identical across")
print("feature sets, so the comparison between rows is fair; absolute values are not.")


# ---------------------------------------------------------------- ceiling
def oracle_and_more(S):
    """oracle recall@k (perfect ranking), plus precision@k and hit@k"""
    Y = B.to_numpy()
    orc = {k: [] for k in KS}
    prec = {k: [] for k in KS}
    hit = {k: [] for k in KS}
    for i in range(Y.shape[0]):
        y, s = Y[i], S[i]
        m = ~np.isnan(y) & ~np.isnan(s)
        tot = y[m].sum()
        if tot < 1:
            continue
        order = np.argsort(-s[m])
        yy = y[m][order]
        for k in KS:
            orc[k].append(min(k, tot) / tot)
            prec[k].append(yy[:k].sum() / k)
            hit[k].append(1.0 if yy[:k].sum() > 0 else 0.0)
    return ({k: float(np.mean(v)) for k, v in orc.items()},
            {k: float(np.mean(v)) for k, v in prec.items()},
            {k: float(np.mean(v)) for k, v in hit.items()})


print("\n" + "=" * 60)
print("ceiling, precision@k and hit@k")
print("=" * 60)
Sprev = prevalence_scores()
orc, _, _ = oracle_and_more(Sprev)
print(f"{'oracle recall@k':<15}" + "".join(f"{orc[k]:>9.3f}" for k in KS))
print("  (perfect ranking; bounded by each strain's own number of positives)")

print(f"\nfraction of the oracle achieved:")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>9}" for k in KS))
for name in configs:
    print(f"{name:<15}" + "".join(f"{results[name][k]/orc[k]:>9.3f}" for k in KS))

print(f"\nprecision@k   (of the k phages you test, how many work)")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>9}" for k in KS))
store = {}
for name, (cols, extra) in configs.items():
    if cols == "prev":
        S = Sprev
    elif cols == "nn":
        S = nn_phylo_scores()
    else:
        S = fit_predict(cols, extra)
    store[name] = S
    _, pr, ht = oracle_and_more(S)
    print(f"{name:<15}" + "".join(f"{pr[k]:>9.3f}" for k in KS))

print(f"\nhit@k   (at least one working phage among the k tested)")
print(f"{'feature set':<15}" + "".join(f"{'k='+str(k):>9}" for k in KS))
for name in configs:
    _, _, ht = oracle_and_more(store[name])
    print(f"{name:<15}" + "".join(f"{ht[k]:>9.3f}" for k in KS))
