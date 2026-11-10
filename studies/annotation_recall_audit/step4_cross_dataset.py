"""Step 4: what can each GenoPHI dataset actually support?

For every shipped interaction matrix: size, sparsity, per-phage and per-strain
positive counts, the oracle recall@k ceiling, the prevalence baseline, and
how many phages clear the threshold needed to derive a single-feature rule.
"""
import sys
import glob
import os
import numpy as np
import pandas as pd

MATDIR = sys.argv[1]
KS = [1, 2, 3, 5, 10]
MIN_POS_FOR_RULE = 5

files = sorted(glob.glob(os.path.join(MATDIR, "*.csv")))
rows = []
for f in files:
    name = os.path.basename(f).replace("_interaction_matrix.csv", "").replace(".csv", "")
    d = pd.read_csv(f)
    if not {"strain", "phage", "interaction"} <= set(d.columns):
        print(f"skip {name}: columns {list(d.columns)}")
        continue
    M = d.pivot_table(index="strain", columns="phage", values="interaction",
                      aggfunc="max")
    Y = M.to_numpy(dtype=float)
    per_phage = np.nansum(Y, axis=0)
    per_strain = np.nansum(Y, axis=1)
    obs = np.sum(~np.isnan(Y))
    pos = int(np.nansum(Y))

    # oracle and prevalence recall@k, per strain
    prev = np.nan_to_num(np.nanmean(Y, axis=0))
    orc, pre, hit = {k: [] for k in KS}, {k: [] for k in KS}, {k: [] for k in KS}
    for i in range(Y.shape[0]):
        y = Y[i]
        m = ~np.isnan(y)
        tot = y[m].sum()
        if tot < 1:
            continue
        yy_prev = y[m][np.argsort(-prev[m])]
        for k in KS:
            orc[k].append(min(k, tot) / tot)
            pre[k].append(yy_prev[:k].sum() / tot)
            hit[k].append(1.0 if yy_prev[:k].sum() > 0 else 0.0)

    rows.append(dict(
        dataset=name,
        strains=M.shape[0], phages=M.shape[1],
        pairs=int(obs), positives=pos,
        pos_rate=100 * pos / obs,
        med_pos_per_phage=float(np.median(per_phage)),
        phages_ge5=int((per_phage >= MIN_POS_FOR_RULE).sum()),
        phages_ge20=int((per_phage >= 20).sum()),
        med_pos_per_strain=float(np.median(per_strain)),
        strains_zero=int((per_strain == 0).sum()),
        strains_eval=len(orc[1]),
        oracle_r2=float(np.mean(orc[2])), oracle_r5=float(np.mean(orc[5])),
        prev_r2=float(np.mean(pre[2])), prev_r5=float(np.mean(pre[5])),
        prev_hit2=float(np.mean(hit[2])), prev_hit5=float(np.mean(hit[5])),
    ))

S = pd.DataFrame(rows)

print("=" * 100)
print("size and sparsity")
print("=" * 100)
print(S[["dataset", "strains", "phages", "pairs", "positives", "pos_rate"]]
      .to_string(index=False, float_format=lambda x: f"{x:.2f}"))

print("\n" + "=" * 100)
print("power for the per-phage rule analysis")
print("=" * 100)
print(S[["dataset", "med_pos_per_phage", "phages_ge5", "phages_ge20", "phages"]]
      .to_string(index=False, float_format=lambda x: f"{x:.1f}"))
print("\n  phages_ge5 is the pool in which a single-feature rule can be derived at all")
print("  (the Escherichia analysis used >= 5 positives; 75 of 96 phages qualified).")

print("\n" + "=" * 100)
print("ranking ceiling and the prevalence baseline, per strain")
print("=" * 100)
print(S[["dataset", "strains_eval", "strains_zero", "med_pos_per_strain",
         "oracle_r2", "prev_r2", "oracle_r5", "prev_r5", "prev_hit2", "prev_hit5"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\n  oracle_rK: perfect ranking recall@K. prev_rK: prevalence baseline recall@K.")
print("  prev_hitK: prevalence baseline, at least one working phage in the top K.")

print("\n" + "=" * 100)
print("prevalence baseline as a fraction of the oracle")
print("=" * 100)
for _, r in S.iterrows():
    print(f"  {r['dataset']:<30} k=2 {r['prev_r2']/r['oracle_r2']:.3f}   "
          f"k=5 {r['prev_r5']/r['oracle_r5']:.3f}   hit@2 {r['prev_hit2']:.3f}")
print("\n  In Escherichia the full annotation set reached 0.62 of the oracle at k=2-5")
print("  against prevalence's 0.43-0.50, and hit@2 0.910 against 0.759.")
