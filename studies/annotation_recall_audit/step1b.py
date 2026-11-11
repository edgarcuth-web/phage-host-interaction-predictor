"""Step 1b: annotation/phylogeny redundancy, and whether the published model's
low recall is a small-sample artefact or a real screening failure."""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_mutual_info_score

DATA = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/data"


def sec(t):
    print("\n" + "=" * 70)
    print(t)
    print("=" * 70)


mat = pd.read_csv(f"{DATA}/interaction_matrix.csv", sep=";").set_index("bacteria")
host = pd.read_csv(f"{DATA}/picard_collection.csv", sep=";").set_index("bacteria")
perf = pd.read_csv(f"{DATA}/per_phage_perf.csv").set_index("phage")
B = (mat >= 1).astype(int)

sec("D. are the curated annotations redundant with phylogeny?")
pg = host["Clermont_Phylo"].astype(str).fillna("NA")
print("adjusted mutual information with Clermont phylogroup")
print("  (0 = independent of phylogeny, 1 = a relabelling of phylogeny)\n")
for c in ["O-type", "H-type", "LPS_type", "ST_Warwick", "ABC_serotype", "Klebs_capsule_type"]:
    if c not in host.columns:
        continue
    v = host[c].astype(str).fillna("__missing__")
    ami = adjusted_mutual_info_score(pg, v)
    print(f"  {c:<20} AMI {ami:5.3f}   levels {v.nunique():>4}   "
          f"missing {int(host[c].isna().sum()):>3}/{len(host)}")

sec("D2. defense systems")
ds = pd.read_csv(f"{DATA}/defense_subtypes.csv", sep=";").set_index("bacteria")
pres = (ds > 0).sum(axis=0).sort_values(ascending=False)
cnt = (ds > 0).sum(axis=1)
print(f"  {ds.shape[0]} strains x {ds.shape[1]} systems")
print(f"  systems never present: {int((pres == 0).sum())}")
print(f"  per-strain count: median {int(cnt.median())}  min {int(cnt.min())}  max {int(cnt.max())}")
print("  most common:", ", ".join(f"{k}={v}" for k, v in pres.head(8).items()))
# does defense repertoire size relate to how many phages infect the strain?
shared = [s for s in B.index if s in cnt.index]
rho, p = stats.spearmanr(cnt.loc[shared], B.loc[shared].sum(axis=1))
print(f"\n  Spearman(n defense systems, n susceptible phages) = {rho:.3f} (p={p:.2g})")
print("  the paper's claim: defense systems play only a marginal role")

sec("E. is the published model's low recall a small-sample artefact?")
pp = B.sum(axis=0)
perf = perf.join(pp.rename("n_pos_public"))
print(f"phages in performance table: {len(perf)}")
print(f"\nrecall: median {perf['recall'].median():.3f}  mean {perf['recall'].mean():.3f}")
print(f"precision: median {perf['precision'].median():.3f}")
print(f"auc: median {perf['auc'].median():.3f}")

rho, p = stats.spearmanr(perf["n_infections"], perf["recall"])
print(f"\nSpearman(n_infections, recall) = {rho:.3f} (p={p:.2g})")
rho2, p2 = stats.spearmanr(perf["n_infections"], perf["auc"])
print(f"Spearman(n_infections, auc)    = {rho2:.3f} (p={p2:.2g})")

print("\nrecall by how many strains the phage actually infects:")
edges = [0, 10, 20, 50, 100, 200, 1e9]
labs = ["<10", "10-19", "20-49", "50-99", "100-199", ">=200"]
g = pd.cut(perf["n_infections"], bins=edges, labels=labs, right=False)
tab = perf.groupby(g, observed=True).agg(
    n_phages=("recall", "size"),
    median_recall=("recall", "median"),
    median_precision=("precision", "median"),
    median_auc=("auc", "median"),
)
print(tab.to_string())

# the decisive question: among phages with PLENTY of positives, is recall still bad?
big = perf[perf["n_infections"] >= 50]
print(f"\nphages with >= 50 infections: {len(big)}")
print(f"  median recall    {big['recall'].median():.3f}")
print(f"  median precision {big['precision'].median():.3f}")
print(f"  median auc       {big['auc'].median():.3f}")
print(f"  phages with recall < 0.5: {int((big['recall'] < 0.5).sum())} / {len(big)}")

# total missed positives, weighted by phage size
perf["missed"] = perf["n_infections"] * (1 - perf["recall"])
print(f"\nimplied missed positives summed over phages: {perf['missed'].sum():.0f}"
      f" of {perf['n_infections'].sum():.0f} ({100*perf['missed'].sum()/perf['n_infections'].sum():.1f}%)")
print(f"  of which from phages with >= 50 infections: "
      f"{big['n_infections'].mul(1-big['recall']).sum():.0f}")
