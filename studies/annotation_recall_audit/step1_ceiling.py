"""Step 1: characterise the ceiling of the Gaborieau Escherichia interaction data.

Reads staged copies of the public tables and reports:
  A. score distribution, binarisation sensitivity, missingness
  B. per-phage and per-strain positive counts (statistical power)
  C. how much of susceptibility is explained by core-genome phylogeny
  D. how much of the curated annotations is redundant with phylogeny
"""
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


# ---------------------------------------------------------------- load
mat = pd.read_csv(f"{DATA}/interaction_matrix.csv", sep=";").set_index("bacteria")
host = pd.read_csv(f"{DATA}/picard_collection.csv", sep=";").set_index("bacteria")
phage = pd.read_csv(f"{DATA}/guelin_collection.csv", sep=";").set_index("phage")
dist = pd.read_csv(f"{DATA}/phylo_dist.tsv", sep="\t").set_index("bacteria")

sec("A. shape, scores, missingness")
print(f"matrix: {mat.shape[0]} strains x {mat.shape[1]} phages = {mat.size} cells")
print(f"picard_collection: {host.shape[0]} strains   guelin_collection: {phage.shape[0]} phages")
print(f"phylo distance matrix: {dist.shape}")

n_nan = int(mat.isna().sum().sum())
print(f"\nNaN cells: {n_nan} ({100 * n_nan / mat.size:.2f}%)")
vals = mat.stack().value_counts().sort_index()
print("\nscore distribution (observed cells):")
for v, c in vals.items():
    print(f"  score {v:>3}: {c:>7}  ({100 * c / vals.sum():5.2f}%)")

print("\nbinarisation sensitivity:")
obs = mat.notna().sum().sum()
for thr in [1, 2, 3, 4]:
    pos = int((mat >= thr).sum().sum())
    print(f"  score >= {thr}: {pos:>6} positives  rate {100 * pos / obs:5.2f}%")

# strains present in matrix but not in collection, and vice versa
only_mat = sorted(set(mat.index) - set(host.index))
only_host = sorted(set(host.index) - set(mat.index))
print(f"\nstrains in matrix not in picard_collection: {len(only_mat)} {only_mat[:5]}")
print(f"strains in picard_collection not in matrix: {len(only_host)} {only_host[:5]}")

# ---------------------------------------------------------------- power
sec("B. per-phage / per-strain positive counts (statistical power)")
B = (mat >= 1).astype(int).where(mat.notna())

per_phage = B.sum(axis=0).sort_values()
per_strain = B.sum(axis=1).sort_values()

print("per-phage positive counts (out of ~402 strains):")
print(f"  median {per_phage.median():.0f}   mean {per_phage.mean():.1f}   "
      f"min {per_phage.min():.0f}   max {per_phage.max():.0f}")
for k in [5, 10, 20, 30, 50]:
    n = int((per_phage < k).sum())
    print(f"  phages with < {k:>3} positives: {n:>2} / {len(per_phage)}  ({100*n/len(per_phage):.0f}%)")
print("\n  10 narrowest phages:", ", ".join(f"{p}={int(v)}" for p, v in per_phage.head(10).items()))
print("  10 broadest  phages:", ", ".join(f"{p}={int(v)}" for p, v in per_phage.tail(10).items()))

print("\nper-strain positive counts (out of 96 phages):")
print(f"  median {per_strain.median():.0f}   mean {per_strain.mean():.1f}   "
      f"min {per_strain.min():.0f}   max {per_strain.max():.0f}")
for k in [1, 3, 5, 10]:
    n = int((per_strain < k).sum())
    print(f"  strains with < {k:>2} susceptible phages: {n:>3} / {len(per_strain)}  ({100*n/len(per_strain):.0f}%)")

# fully resistant strains are untreatable targets; they also carry no rescue signal
print(f"\n  strains with ZERO susceptible phage: {int((per_strain == 0).sum())}")

# ---------------------------------------------------------------- phylogeny
sec("C. how much of susceptibility is explained by core-genome phylogeny")
shared = [s for s in B.index if s in dist.index and s in dist.columns]
print(f"strains with both susceptibility profile and phylogenetic distance: {len(shared)}")

Bs = B.loc[shared]
D = dist.loc[shared, shared].astype(float)

# pairwise Jaccard similarity of susceptibility profiles, on co-observed phages only
A = Bs.to_numpy(dtype=float)
M = Bs.notna().to_numpy()
A = np.nan_to_num(A)

n = len(shared)
iu = np.triu_indices(n, k=1)
both = (M.astype(float) @ M.astype(float).T)          # co-observed phage count
inter = (A @ A.T)                                      # both positive
rowpos = (A * M).sum(axis=1)
union = rowpos[:, None] + rowpos[None, :] - inter
with np.errstate(divide="ignore", invalid="ignore"):
    jac = np.where(union > 0, inter / union, np.nan)

d_flat = D.to_numpy()[iu]
j_flat = jac[iu]
co_flat = both[iu]

ok = (~np.isnan(j_flat)) & (~np.isnan(d_flat)) & (co_flat >= 50)
print(f"usable strain pairs (co-observed >= 50 phages): {ok.sum()} of {len(j_flat)}")

rho, p = stats.spearmanr(d_flat[ok], j_flat[ok])
print(f"\nSpearman(phylogenetic distance, susceptibility-profile Jaccard) = {rho:.3f}  (p={p:.2e})")
print("  negative = closer relatives share susceptibility profiles")

# binned view
bins = [0, 0.001, 0.005, 0.01, 0.02, 0.03, 1.0]
lab = ["<0.001", "0.001-0.005", "0.005-0.01", "0.01-0.02", "0.02-0.03", ">0.03"]
idx = np.digitize(d_flat[ok], bins) - 1
print("\n  mean profile Jaccard by phylogenetic distance bin:")
for i, l in enumerate(lab):
    m = idx == i
    if m.sum() > 0:
        print(f"    {l:>12}: Jaccard {j_flat[ok][m].mean():.3f}  (n={m.sum():>6})")

# same vs different phylogroup
pg = host["Clermont_Phylo"].reindex(shared)
same_pg = (pg.to_numpy()[:, None] == pg.to_numpy()[None, :])[iu][ok]
print(f"\n  same phylogroup   : mean Jaccard {j_flat[ok][same_pg].mean():.3f}  (n={same_pg.sum()})")
print(f"  different phylogroup: mean Jaccard {j_flat[ok][~same_pg].mean():.3f}  (n={(~same_pg).sum()})")

# ---------------------------------------------------------------- redundancy
sec("D. are the curated annotations redundant with phylogeny?")
cols = ["O-type", "H-type", "LPS_type", "ST_Warwick", "ABC_serotype", "Klebs_capsule_type"]
pg_all = host["Clermont_Phylo"].astype(str)
print("adjusted mutual information with Clermont phylogroup (0 = independent, 1 = identical):")
for c in cols:
    if c not in host.columns:
        continue
    v = host[c].astype(str).fillna("NA")
    ami = adjusted_mutual_info_score(pg_all, v)
    nlev = v.nunique()
    nmiss = int((host[c].isna()).sum())
    print(f"  {c:<20} AMI {ami:5.3f}   levels {nlev:>4}   missing {nmiss:>3}/{len(host)}")

print("\ndefense systems:")
ds = pd.read_csv(f"{DATA}/defense_subtypes.csv", sep=";").set_index("bacteria")
print(f"  matrix {ds.shape[0]} strains x {ds.shape[1]} systems")
pres = (ds > 0).sum(axis=0).sort_values(ascending=False)
print(f"  systems present in 0 strains: {int((pres == 0).sum())}")
print(f"  per-strain system count: median {int((ds>0).sum(axis=1).median())}  "
      f"min {int((ds>0).sum(axis=1).min())}  max {int((ds>0).sum(axis=1).max())}")
print("  10 most common:", ", ".join(f"{k}={v}" for k, v in pres.head(10).items()))
