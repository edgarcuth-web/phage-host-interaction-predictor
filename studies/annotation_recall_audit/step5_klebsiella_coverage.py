"""Step 5: the Klebsiella question the data can actually answer.

Per-phage power is too low for the rule analysis used on Escherichia
(median 2 positives per phage). So pool across phages and ask a coverage
question instead, which needs no fitting at all:

  For an experimentally positive (strain, phage) pair, is there another
  strain that this same phage also infects and whose capsule locus is
  similar enough that a "match the K-locus" rule would have found it?

The annotation route here is K-locus identity, given directly as the
all_loci_score_matrix. The annotation-free alternative is cosine similarity
of the ESM-2 embeddings of the same loci -- PhageHostLearn's own
representation. Neither is trained here, so the comparison is between two
similarity measures over identical sequence data, not between a model and a
rule. Leave-one-out throughout: a pair never sees itself.
"""
import json
import sys
import numpy as np
import pandas as pd

D = sys.argv[1]
THRESHOLDS = [0.99, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]

M = pd.read_csv(f"{D}/phage_host_interactions.csv", index_col=0)
M = M.dropna(axis=1, how="all").dropna(axis=0, how="all")
with open(f"{D}/Locibase.json") as fh:
    loci_keys = list(json.load(fh).keys())
A = pd.read_csv(f"{D}/all_loci_score_matrix.txt", sep="\t", header=None).to_numpy(float)
E = pd.read_csv(f"{D}/esm2_embeddings_loci.csv", index_col=0)

strains = [s for s in M.index if s in loci_keys]
M = M.loc[strains]
pos_in_loci = {s: i for i, s in enumerate(loci_keys)}
idx = [pos_in_loci[s] for s in strains]
K = A[np.ix_(idx, idx)]                       # K-locus identity, annotation route

Ee = E.reindex(strains)
print(f"strains with locus data and interactions: {len(strains)}")
print(f"ESM-2 embeddings available for {int(Ee.notna().all(axis=1).sum())} of them")
X = Ee.to_numpy(float)
X = X / np.linalg.norm(X, axis=1, keepdims=True)
C = X @ X.T                                   # cosine similarity, learned route
np.fill_diagonal(K, np.nan)
np.fill_diagonal(C, np.nan)

Y = M.to_numpy(float)                          # NaN = untested
n_s, n_p = Y.shape
print(f"matrix {n_s} x {n_p}   observed {int(np.sum(~np.isnan(Y)))}   "
      f"positives {int(np.nansum(Y))}")

off = K[~np.isnan(K)]
print(f"\nK-locus identity between strains: median {np.median(off):.3f}, "
      f"{100*np.mean(off>=0.9):.2f}% of pairs >= 0.90")
offc = C[~np.isnan(C)]
print(f"ESM-2 cosine between strains:      median {np.median(offc):.3f}")


def nearest_positive_similarity(S):
    """for each observed cell, the highest similarity to another strain that
    this same phage is observed to infect (leave-one-out)"""
    best = np.full(Y.shape, np.nan)
    for j in range(n_p):
        y = Y[:, j]
        obs = ~np.isnan(y)
        pos = np.where(obs & (y == 1))[0]
        if len(pos) == 0:
            continue
        sub = S[:, pos]                        # every strain vs this phage's hosts
        for i in np.where(obs)[0]:
            v = sub[i]
            v = v[~np.isnan(v)]
            if len(v):
                best[i, j] = v.max()
    return best


print("\ncomputing leave-one-out nearest-infected-host similarity ...", flush=True)
bestK = nearest_positive_similarity(K)
bestC = nearest_positive_similarity(C)

obs = ~np.isnan(Y)
ypos = obs & (Y == 1)
yneg = obs & (Y == 0)
print(f"cells usable: {int(np.sum(obs & ~np.isnan(bestK)))}  "
      f"of which positive {int(np.sum(ypos & ~np.isnan(bestK)))}")


def rule_table(best, label):
    print(f"\n{label}")
    print(f"  {'threshold':>10}{'recall':>9}{'precision':>11}{'n predicted':>13}")
    for t in THRESHOLDS:
        pred = obs & ~np.isnan(best) & (best >= t)
        tp = int(np.sum(pred & ypos))
        np_ = int(pred.sum())
        npos = int(np.sum(ypos & ~np.isnan(best)))
        print(f"  {t:>10.2f}{tp/npos if npos else np.nan:>9.3f}"
              f"{tp/np_ if np_ else np.nan:>11.3f}{np_:>13}")


rule_table(bestK, "K-locus identity rule: predict lytic if some infected host is >= t similar")
rule_table(bestC, "ESM-2 cosine rule: same, with learned locus embeddings")

print("\n" + "=" * 72)
print("how many positives lie outside the reach of K-locus identity?")
print("=" * 72)
m = ypos & ~np.isnan(bestK)
v = bestK[m]
print(f"positives considered: {len(v)}")
for t in THRESHOLDS:
    out = int(np.sum(v < t))
    print(f"  no infected host at K-locus identity >= {t:.2f}: "
          f"{out:>4} / {len(v)}  ({100*out/len(v):.1f}%)")
print(f"\n  distribution of each positive's nearest infected host (K-locus identity):")
for q in [5, 25, 50, 75, 95]:
    print(f"    {q:>2}th percentile {np.percentile(v, q):.3f}")

print("\n" + "=" * 72)
print("does the learned embedding reach positives that K-locus identity misses?")
print("=" * 72)
both = ypos & ~np.isnan(bestK) & ~np.isnan(bestC)
# compare at matched precision rather than matched threshold
for target in [0.90, 0.80, 0.70, 0.60, 0.50]:
    row = []
    for best, nm in [(bestK, "K-locus"), (bestC, "ESM-2")]:
        cand = np.sort(best[obs & ~np.isnan(best)])[::-1]
        chosen, rec = None, np.nan
        for t in cand[::max(1, len(cand)//400)]:
            pred = obs & ~np.isnan(best) & (best >= t)
            tp = int(np.sum(pred & ypos))
            if pred.sum() == 0:
                continue
            if tp / pred.sum() >= target:
                chosen, rec = t, tp / int(np.sum(ypos & ~np.isnan(best)))
        row.append((nm, chosen, rec))
    s = "   ".join(f"{nm}: t={t if t is None else round(t,3)} recall="
                   f"{'na' if np.isnan(r) else round(r,3)}" for nm, t, r in row)
    print(f"  at precision >= {target:.2f}   {s}")

print("\n  rescued positives: reached by ESM-2 but not by K-locus identity,")
print("  at each route's own threshold for precision >= 0.70")


def thresh_for_precision(best, target):
    cand = np.sort(best[obs & ~np.isnan(best)])[::-1]
    out = None
    for t in cand[::max(1, len(cand)//400)]:
        pred = obs & ~np.isnan(best) & (best >= t)
        if pred.sum() and np.sum(pred & ypos) / pred.sum() >= target:
            out = t
    return out


tk, tc = thresh_for_precision(bestK, 0.70), thresh_for_precision(bestC, 0.70)
if tk is not None and tc is not None:
    pk = obs & ~np.isnan(bestK) & (bestK >= tk)
    pc = obs & ~np.isnan(bestC) & (bestC >= tc)
    npos = int(np.sum(ypos & both))
    print(f"    K-locus threshold {tk:.3f}, ESM-2 threshold {tc:.3f}")
    print(f"    positives found by both          {int(np.sum(ypos & pk & pc))}")
    print(f"    found only by K-locus identity   {int(np.sum(ypos & pk & ~pc))}")
    print(f"    found only by ESM-2 (rescued)    {int(np.sum(ypos & ~pk & pc))}")
    print(f"    found by neither                 {int(np.sum(ypos & ~pk & ~pc))}")
    print(f"    total positives                  {npos}")
