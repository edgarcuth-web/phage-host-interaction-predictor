"""Step 2a: generalise the paper's K1 experiment.

For every phage, derive the best single-feature annotation rule
("lytic if and only if the strain's <feature> equals <level>"), then measure
how that rule transfers to other phages of the same genus / family.

Step 2b: decompose the ABC capsule signal into "has a capsule at all"
versus "which serotype".
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

DATA = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/data"
THR = int(sys.argv[2]) if len(sys.argv) > 2 else 1

mat = pd.read_csv(f"{DATA}/interaction_matrix.csv", sep=";").set_index("bacteria")
host = pd.read_csv(f"{DATA}/picard_collection.csv", sep=";").set_index("bacteria")
phage = pd.read_csv(f"{DATA}/guelin_collection.csv", sep=";").set_index("phage")

B = (mat >= THR).astype(float).where(mat.notna())
strains = list(B.index)
H = host.reindex(strains)

FEATS = ["ABC_serotype", "LPS_type", "O-type", "H-type", "ST_Warwick", "Clermont_Phylo"]
MIN_SUPPORT = 10          # a level must cover at least this many strains
MIN_SOURCE_PREC = 0.50    # the derived rule must be at least this precise on its own phage


def prf(pred, y):
    """precision, recall, n_pred, n_pos on co-observed cells"""
    m = ~np.isnan(y)
    pred, y = pred[m].astype(bool), y[m].astype(bool)
    tp = int((pred & y).sum())
    return (tp / pred.sum() if pred.sum() else np.nan,
            tp / y.sum() if y.sum() else np.nan,
            int(pred.sum()), int(y.sum()))


print("=" * 72)
print(f"2a. single-feature annotation rules and how they transfer (score >= {THR})")
print("=" * 72)

rules = []
for p in B.columns:
    y = B[p].to_numpy()
    if np.nansum(y) < 5:
        continue
    best = None
    for f in FEATS:
        if f not in H.columns:
            continue
        col = H[f].astype("object")
        for lev, cnt in col.value_counts().items():
            if cnt < MIN_SUPPORT:
                continue
            pred = (col == lev).to_numpy()
            pr, rc, npred, npos = prf(pred, y)
            if np.isnan(pr) or pr < MIN_SOURCE_PREC:
                continue
            f1 = 0 if (pr + rc) == 0 else 2 * pr * rc / (pr + rc)
            if best is None or f1 > best["f1"]:
                best = dict(phage=p, feature=f, level=str(lev), f1=f1,
                            prec=pr, rec=rc, n_strains=int(cnt), n_pos=npos)
    if best:
        rules.append(best)

R = pd.DataFrame(rules)
print(f"\nphages for which a rule with precision >= {MIN_SOURCE_PREC} exists: "
      f"{len(R)} of {B.shape[1]}")
print(f"  on its own phage: median precision {R['prec'].median():.3f}, "
      f"median recall {R['rec'].median():.3f}")
print("\n  which feature supplies the rule:")
for f, c in R["feature"].value_counts().items():
    print(f"    {f:<18} {c:>3} phages")

print("\n  10 strongest rules by F1:")
top = R.sort_values("f1", ascending=False).head(10)
for _, r in top.iterrows():
    print(f"    {r['phage']:<14} {r['feature']}=={r['level']:<10} "
          f"prec {r['prec']:.2f}  rec {r['rec']:.2f}  ({r['n_pos']:>3} positives)")

# ---- transfer to related phages -------------------------------------------
print("\n" + "-" * 72)
print("transfer of each rule to OTHER phages, grouped by phage relatedness")
print("-" * 72)

ph = phage.reindex(B.columns)
rows = []
for _, r in R.iterrows():
    src = r["phage"]
    f, lev = r["feature"], r["level"]
    pred = (H[f].astype("object").astype(str) == lev).to_numpy()
    g_src = ph.loc[src, "Genus"] if src in ph.index else None
    fam_src = ph.loc[src, "Family"] if src in ph.index else None
    for q in B.columns:
        if q == src:
            continue
        y = B[q].to_numpy()
        if np.nansum(y) < 5:
            continue
        pr, rc, _, npos = prf(pred, y)
        same_genus = (g_src is not None and ph.loc[q, "Genus"] == g_src)
        same_family = (fam_src is not None and ph.loc[q, "Family"] == fam_src)
        rows.append(dict(src=src, tgt=q, prec=pr, rec=rc, n_pos=npos,
                         same_genus=bool(same_genus), same_family=bool(same_family)))

T = pd.DataFrame(rows)


def summarise(sub, label):
    if len(sub) == 0:
        print(f"  {label:<34} (no pairs)")
        return
    print(f"  {label:<34} n={len(sub):>5}   median precision {sub['prec'].median():.3f}   "
          f"median recall {sub['rec'].median():.3f}")


print(f"\nsource phage itself:               n={len(R):>5}   "
      f"median precision {R['prec'].median():.3f}   median recall {R['rec'].median():.3f}")
summarise(T[T.same_genus], "transferred, same genus")
summarise(T[~T.same_genus & T.same_family], "transferred, same family only")
summarise(T[~T.same_family], "transferred, different family")

print("\n  interpretation: a mechanism rule derived on one phage, applied to a")
print("  relative of that phage. The paper's K1 example is one row of this table.")

# the paper's own K1 example, reproduced
k1 = R[(R.feature == "ABC_serotype") & (R.level == "K1")]
if len(k1):
    print(f"\n  ABC_serotype==K1 rules were derived for {len(k1)} phages: "
          f"{', '.join(k1['phage'])}")
    k1t = T[T.src.isin(k1["phage"])]
    print(f"    on their own phages : median precision {k1['prec'].median():.3f}  "
          f"median recall {k1['rec'].median():.3f}")
    sg = k1t[k1t.same_genus]
    if len(sg):
        print(f"    same genus elsewhere: median precision {sg['prec'].median():.3f}  "
              f"median recall {sg['rec'].median():.3f}  (n={len(sg)})")

# ---- 2b. ABC capsule decomposition ---------------------------------------
print("\n" + "=" * 72)
print("2b. is the ABC capsule signal presence, or serotype identity?")
print("=" * 72)

abc = H["ABC_serotype"]
has_abc = abc.notna().to_numpy().astype(float)
print(f"strains with an ABC serotype: {int(has_abc.sum())} / {len(H)}")
print(f"serotype levels among them  : {abc.nunique()}")
print("  level counts:", ", ".join(f"{k}={v}" for k, v in abc.value_counts().head(8).items()))


def auroc(score, y):
    m = ~np.isnan(y)
    y2, s2 = y[m], score[m]
    if y2.sum() == 0 or y2.sum() == len(y2):
        return np.nan
    pos, neg = s2[y2 == 1], s2[y2 == 0]
    # Mann-Whitney U / (n1*n2)
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


pres_auc, sero_auc, npos_list = [], [], []
for p in B.columns:
    y = B[p].to_numpy()
    if np.nansum(y) < 20:
        continue
    npos_list.append(np.nansum(y))
    # (a) presence only
    pres_auc.append(auroc(has_abc, y))
    # (b) serotype identity, restricted to capsulated strains: use per-level
    #     positive rate as the score (in-sample, so this is an upper bound)
    sub = has_abc == 1
    ysub = y[sub]
    lev = abc[sub].astype(str).to_numpy()
    rate = pd.Series(ysub).groupby(lev).transform("mean").to_numpy()
    sero_auc.append(auroc(rate, ysub))

pa = np.array(pres_auc, dtype=float)
sa = np.array(sero_auc, dtype=float)
ok = ~np.isnan(pa) & ~np.isnan(sa)
print(f"\nphages with >= 20 positives: {ok.sum()}")
print(f"  AUROC of 'has an ABC capsule' alone        : median {np.median(pa[ok]):.3f}")
print(f"  AUROC of serotype identity within capsulated: median {np.median(sa[ok]):.3f}"
      "   (in-sample, an upper bound)")
w = stats.wilcoxon(pa[ok], sa[ok])
print(f"  Wilcoxon on the pair: p = {w.pvalue:.2g}")
print("\n  presence alone above 0.5 means capsulated strains differ in susceptibility")
print("  from uncapsulated ones; the serotype column is the authors' load-bearing")
print("  feature, so which half carries it matters for any mechanism claim.")
