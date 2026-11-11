"""Step 5c: 80% of Klebsiella positives are out of reach of host-locus
similarity. Is the missing signal on the phage side?

Same leave-one-out coverage question, but over phages: for a positive
(strain, phage), is there another phage that infects this same strain and
whose receptor-binding proteins resemble this phage's? RBP similarity is a
3-mer profile cosine over the protein sequences, a crude stand-in for
alignment, taking each phage's best-matching RBP pair.
"""
import sys
import numpy as np
import pandas as pd
from itertools import product
from scipy import stats

D = sys.argv[1]
M = pd.read_csv(f"{D}/phage_host_interactions.csv", index_col=0)
M = M.dropna(axis=1, how="all").dropna(axis=0, how="all")
R = pd.read_csv(f"{D}/RBPbase.csv")

phages = [p for p in M.columns if p in set(R["phage_ID"])]
M = M[phages]
Y = M.to_numpy(float)
print(f"phages with RBP data: {len(phages)} of {len(R['phage_ID'].unique())} in RBPbase")
print(f"matrix {Y.shape[0]} strains x {Y.shape[1]} phages   "
      f"observed {int(np.sum(~np.isnan(Y)))}  positives {int(np.nansum(Y))}")

AAS = "ACDEFGHIKLMNPQRSTVWY"
kmers = {"".join(t): i for i, t in enumerate(product(AAS, repeat=3))}


def profile(seq):
    v = np.zeros(len(kmers), dtype=np.float32)
    for i in range(len(seq) - 2):
        j = kmers.get(seq[i:i + 3])
        if j is not None:
            v[j] += 1
    n = np.linalg.norm(v)
    return v / n if n else v


prot = {p: [profile(s) for s in R.loc[R["phage_ID"] == p, "protein_sequence"]]
        for p in phages}
print(f"RBPs per phage: median "
      f"{int(np.median([len(prot[p]) for p in phages]))}")

n_p = len(phages)
S = np.full((n_p, n_p), np.nan)
for a in range(n_p):
    for b in range(a + 1, n_p):
        best = max((float(u @ v) for u in prot[phages[a]] for v in prot[phages[b]]),
                   default=np.nan)
        S[a, b] = S[b, a] = best
print(f"phage RBP similarity (best matching pair): median "
      f"{np.nanmedian(S):.3f}  min {np.nanmin(S):.3f}  max {np.nanmax(S):.3f}")

obs = ~np.isnan(Y)
ypos = obs & (Y == 1)

# leave-one-out: for each observed cell, best RBP similarity to another phage
# that infects the same strain
best = np.full(Y.shape, np.nan)
for i in range(Y.shape[0]):
    y = Y[i]
    o = ~np.isnan(y)
    pos = np.where(o & (y == 1))[0]
    if len(pos) == 0:
        continue
    for j in np.where(o)[0]:
        v = S[j, pos]
        v = v[~np.isnan(v)]
        if len(v):
            best[i, j] = v.max()

m = obs & ~np.isnan(best)
s, yy = best[m], Y[m]
p_, n_ = s[yy == 1], s[yy == 0]
u = stats.mannwhitneyu(p_, n_, alternative="two-sided").statistic
print(f"\npooled AUROC of nearest-co-infecting-phage RBP similarity: "
      f"{u/(len(p_)*len(n_)):.3f}")
print(f"  positives {len(p_)}  negatives {len(n_)}")

print("\n  threshold   recall  precision  n predicted")
for t in [0.99, 0.95, 0.90, 0.80, 0.70, 0.60]:
    pred = m & (best >= t)
    tp = int(np.sum(pred & ypos))
    print(f"  {t:>9.2f}{tp/len(p_):>9.3f}"
          f"{tp/pred.sum() if pred.sum() else np.nan:>11.3f}{int(pred.sum()):>13}")

print("\n  how many positives have no co-infecting phage with similar RBPs?")
v = best[ypos & ~np.isnan(best)]
for t in [0.95, 0.90, 0.80, 0.70]:
    print(f"    none at RBP similarity >= {t:.2f}: "
          f"{int(np.sum(v < t)):>4} / {len(v)}  ({100*np.mean(v < t):.1f}%)")
