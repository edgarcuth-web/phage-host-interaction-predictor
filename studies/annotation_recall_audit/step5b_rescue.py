"""Step 5b: fix the embedding comparison and score it properly.

Raw ESM-2 cosine between K-locus embeddings is degenerate (median 0.999), so
centre the embeddings first. Compare the two similarity routes by how well
each ranks positives above negatives for the same phage, which avoids
threshold artefacts entirely.
"""
import json, sys
import numpy as np
import pandas as pd
from scipy import stats

D = sys.argv[1]
M = pd.read_csv(f"{D}/phage_host_interactions.csv", index_col=0)
M = M.dropna(axis=1, how="all").dropna(axis=0, how="all")
with open(f"{D}/Locibase.json") as fh:
    lk = list(json.load(fh).keys())
A = pd.read_csv(f"{D}/all_loci_score_matrix.txt", sep="\t", header=None).to_numpy(float)
E = pd.read_csv(f"{D}/esm2_embeddings_loci.csv", index_col=0)

strains = [s for s in M.index if s in lk]
M = M.loc[strains]
pi = {s: i for i, s in enumerate(lk)}
K = A[np.ix_([pi[s] for s in strains], [pi[s] for s in strains])]

X = E.reindex(strains).to_numpy(float)
Xc = X - X.mean(axis=0, keepdims=True)
Xc = Xc / np.linalg.norm(Xc, axis=1, keepdims=True)
C = Xc @ Xc.T
np.fill_diagonal(K, np.nan); np.fill_diagonal(C, np.nan)
print(f"centred ESM-2 cosine: median {np.nanmedian(C):.3f}  "
      f"min {np.nanmin(C):.3f}  max {np.nanmax(C):.3f}")
ku, cu = K[~np.isnan(K)], C[~np.isnan(C)]
print(f"Spearman(K-locus identity, centred ESM-2 cosine) = "
      f"{stats.spearmanr(ku, cu).statistic:.3f}")

Y = M.to_numpy(float)
n_s, n_p = Y.shape

def nps(S):
    best = np.full(Y.shape, np.nan)
    for j in range(n_p):
        y = Y[:, j]; obs = ~np.isnan(y)
        pos = np.where(obs & (y == 1))[0]
        if len(pos) == 0: continue
        sub = S[:, pos]
        for i in np.where(obs)[0]:
            v = sub[i]; v = v[~np.isnan(v)]
            if len(v): best[i, j] = v.max()
    return best

bK, bC = nps(K), nps(C)
obs = ~np.isnan(Y); ypos = obs & (Y == 1)

print("\nper-phage AUROC of nearest-infected-host similarity (phages with >=3 positives")
print("and >=3 negatives among observed cells), leave-one-out:")
rows = []
for j in range(n_p):
    y = Y[:, j]; o = ~np.isnan(y)
    for best, nm in [(bK, "K"), (bC, "C")]:
        pass
    m = o & ~np.isnan(bK[:, j]) & ~np.isnan(bC[:, j])
    yy = y[m]
    if yy.sum() < 3 or (len(yy) - yy.sum()) < 3: continue
    r = {}
    for best, nm in [(bK, "K-locus"), (bC, "ESM-2")]:
        s = best[m, j]
        p_, n_ = s[yy == 1], s[yy == 0]
        u = stats.mannwhitneyu(p_, n_, alternative="two-sided").statistic
        r[nm] = u / (len(p_) * len(n_))
    r["phage"] = M.columns[j]; r["n_pos"] = int(yy.sum())
    rows.append(r)
R = pd.DataFrame(rows)
print(f"  phages evaluated: {len(R)}  (total positives in them: {R['n_pos'].sum()})")
for nm in ["K-locus", "ESM-2"]:
    print(f"  {nm:<9} median AUROC {R[nm].median():.3f}   mean {R[nm].mean():.3f}")
w = stats.wilcoxon(R["K-locus"], R["ESM-2"])
print(f"  paired Wilcoxon: p = {w.pvalue:.3g}")
print(f"  phages where ESM-2 beats K-locus: {int((R['ESM-2']>R['K-locus']).sum())}/{len(R)}")

print("\npooled over all phages, ranking every observed cell:")
for best, nm in [(bK, "K-locus"), (bC, "ESM-2")]:
    m = obs & ~np.isnan(best)
    s, yy = best[m], Y[m]
    p_, n_ = s[yy == 1], s[yy == 0]
    u = stats.mannwhitneyu(p_, n_, alternative="two-sided").statistic
    print(f"  {nm:<9} AUROC {u/(len(p_)*len(n_)):.3f}   "
          f"positives {len(p_)}  negatives {len(n_)}")

print("\nrescue at matched budget: take the top N cells by each route,")
print("N = number of cells the K-locus rule selects at identity >= 0.90")
Npick = int(np.sum(obs & ~np.isnan(bK) & (bK >= 0.90)))
both = obs & ~np.isnan(bK) & ~np.isnan(bC)
def topN(best, N):
    m = both.copy()
    s = np.where(m, best, -np.inf)
    flat = np.argsort(s.ravel())[::-1][:N]
    out = np.zeros_like(m)
    out.ravel()[flat] = True
    return out & m
pk, pc = topN(bK, Npick), topN(bC, Npick)
tot = int(np.sum(ypos & both))
print(f"  N = {Npick}   total positives {tot}")
print(f"  found by both                {int(np.sum(ypos & pk & pc))}")
print(f"  only by K-locus identity     {int(np.sum(ypos & pk & ~pc))}")
print(f"  only by ESM-2 (rescued)      {int(np.sum(ypos & ~pk & pc))}")
print(f"  by neither                   {int(np.sum(ypos & ~pk & ~pc))}")
print(f"  union recall {int(np.sum(ypos & (pk|pc)))/tot:.3f} vs "
      f"K-locus alone {int(np.sum(ypos & pk))/tot:.3f}")
