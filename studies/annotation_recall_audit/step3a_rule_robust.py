"""Step 3a: make the rule-transfer result solid.

Three additions over step 2a:
  1. all four binarisation thresholds
  2. a paired test on the same labels -- the transferred rule versus the
     target phage's OWN best rule, both evaluated on the target phage
  3. an analytic chance level: a rule covering m of n strains has expected
     recall m/n under independence, so "lift" = recall - m/n
  4. cluster bootstrap over source phages for confidence intervals
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

DATA = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/data"
RNG = np.random.default_rng(0)
N_BOOT = 2000

mat = pd.read_csv(f"{DATA}/interaction_matrix.csv", sep=";").set_index("bacteria")
host = pd.read_csv(f"{DATA}/picard_collection.csv", sep=";").set_index("bacteria")
phage = pd.read_csv(f"{DATA}/guelin_collection.csv", sep=";").set_index("phage")

FEATS = ["ABC_serotype", "LPS_type", "O-type", "H-type", "ST_Warwick", "Clermont_Phylo"]
MIN_SUPPORT = 10
MIN_PREC = 0.50


def derive_rules(B, H):
    """best single-feature rule per phage, by F1, in-sample"""
    out = {}
    cols = {f: H[f].astype("object").astype(str).fillna("__missing__").to_numpy()
            for f in FEATS if f in H.columns}
    for p in B.columns:
        y = B[p].to_numpy()
        obs = ~np.isnan(y)
        if np.nansum(y) < 5:
            continue
        best = None
        for f, col in cols.items():
            for lev in pd.unique(col):
                pred = col == lev
                if pred.sum() < MIN_SUPPORT:
                    continue
                m = obs & pred
                tp = np.nansum(y[m])
                if m.sum() == 0:
                    continue
                pr = tp / m.sum()
                rc = tp / np.nansum(y[obs])
                if pr < MIN_PREC:
                    continue
                f1 = 0 if (pr + rc) == 0 else 2 * pr * rc / (pr + rc)
                if best is None or f1 > best["f1"]:
                    best = dict(feature=f, level=lev, f1=f1, prec=pr, rec=rc,
                                support=int(pred.sum()))
        if best:
            out[p] = best
    return out


def eval_rule(pred, y):
    obs = ~np.isnan(y)
    m = obs & pred
    npos = np.nansum(y[obs])
    if npos == 0 or m.sum() == 0:
        return np.nan, np.nan
    tp = np.nansum(y[m])
    return tp / m.sum(), tp / npos


def boot_ci(vals, groups, stat=np.median):
    """cluster bootstrap over source phages"""
    vals, groups = np.asarray(vals), np.asarray(groups)
    uniq = np.unique(groups)
    out = []
    for _ in range(N_BOOT):
        pick = RNG.choice(uniq, size=len(uniq), replace=True)
        v = np.concatenate([vals[groups == g] for g in pick])
        v = v[~np.isnan(v)]
        if len(v):
            out.append(stat(v))
    return (np.percentile(out, 2.5), np.percentile(out, 97.5)) if out else (np.nan, np.nan)


print("=" * 78)
print("rule transfer across binarisation thresholds")
print("=" * 78)

summary = []
for THR in [1, 2, 3, 4]:
    B = (mat >= THR).astype(float).where(mat.notna())
    H = host.reindex(B.index)
    n_strains = len(B)
    rules = derive_rules(B, H)
    ph = phage.reindex(B.columns)
    cols = {f: H[f].astype("object").astype(str).fillna("__missing__").to_numpy()
            for f in FEATS if f in H.columns}

    rows = []
    for src, r in rules.items():
        pred = cols[r["feature"]] == r["level"]
        chance_src = r["support"] / n_strains
        g_src = ph.loc[src, "Genus"] if src in ph.index else None
        for tgt in B.columns:
            if tgt == src or tgt not in rules:
                continue
            same_genus = (g_src is not None and ph.loc[tgt, "Genus"] == g_src)
            if not same_genus:
                continue
            y = B[tgt].to_numpy()
            pr_t, rc_t = eval_rule(pred, y)
            own = rules[tgt]
            pred_own = cols[own["feature"]] == own["level"]
            pr_o, rc_o = eval_rule(pred_own, y)
            rows.append(dict(src=src, tgt=tgt,
                             rec_transfer=rc_t, prec_transfer=pr_t,
                             rec_own=rc_o, prec_own=pr_o,
                             lift_transfer=rc_t - pred.sum() / n_strains,
                             lift_own=rc_o - pred_own.sum() / n_strains))
    T = pd.DataFrame(rows)
    src_rec = np.array([r["rec"] for r in rules.values()])
    src_prec = np.array([r["prec"] for r in rules.values()])

    if len(T) < 10:
        print(f"\nscore >= {THR}   positives {int(np.nansum(B.to_numpy())):>5}   "
              f"phages with a rule {len(rules)}/{B.shape[1]}   "
              f"same-genus pairs {len(T)} -- too few to test")
        continue

    ok = T["rec_transfer"].notna() & T["rec_own"].notna()
    Tk = T[ok]
    w = stats.wilcoxon(Tk["rec_own"], Tk["rec_transfer"])
    lo_t, hi_t = boot_ci(Tk["rec_transfer"], Tk["src"])
    lo_o, hi_o = boot_ci(Tk["rec_own"], Tk["src"])

    print(f"\nscore >= {THR}   positives {int(np.nansum(B.to_numpy())):>5}   "
          f"phages with a rule {len(rules)}/{B.shape[1]}   same-genus pairs {len(Tk)}")
    print(f"  rule on its own phage (in-sample upper bound)")
    print(f"     median recall    {np.median(src_rec):.3f}   median precision {np.median(src_prec):.3f}")
    print(f"  SAME-GENUS TARGET, evaluated on identical labels:")
    print(f"     own rule        median recall {np.median(Tk['rec_own']):.3f}  "
          f"95% CI [{lo_o:.3f}, {hi_o:.3f}]")
    print(f"     transferred rule median recall {np.median(Tk['rec_transfer']):.3f}  "
          f"95% CI [{lo_t:.3f}, {hi_t:.3f}]")
    print(f"     paired Wilcoxon (own vs transferred): p = {w.pvalue:.2e}")
    print(f"     relative recall retained: "
          f"{np.median(Tk['rec_transfer'])/np.median(Tk['rec_own']):.2f}")
    print(f"  lift over chance (recall minus the rule's strain coverage):")
    print(f"     own rule        median {np.median(Tk['lift_own']):+.3f}")
    print(f"     transferred     median {np.median(Tk['lift_transfer']):+.3f}")
    lo_l, hi_l = boot_ci(Tk["lift_transfer"], Tk["src"])
    print(f"     transferred lift 95% CI [{lo_l:+.3f}, {hi_l:+.3f}]"
          f"   {'(includes 0: no better than chance)' if lo_l <= 0 <= hi_l else '(excludes 0)'}")

    summary.append(dict(thr=THR, n_rules=len(rules), n_pairs=len(Tk),
                        src_rec=np.median(src_rec),
                        own_rec=np.median(Tk["rec_own"]),
                        trans_rec=np.median(Tk["rec_transfer"]),
                        trans_lo=lo_t, trans_hi=hi_t,
                        lift=np.median(Tk["lift_transfer"]),
                        lift_lo=lo_l, lift_hi=hi_l,
                        p=w.pvalue))

print("\n" + "=" * 78)
print("summary across thresholds")
print("=" * 78)
S = pd.DataFrame(summary)
print(S.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
