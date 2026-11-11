"""Inspect the PhageHostLearn archive: shapes, sparsity, and whether the
strain identifiers line up across the interaction matrix, Locibase and the
K-locus similarity matrix."""
import json
import sys
import numpy as np
import pandas as pd

D = sys.argv[1]

M = pd.read_csv(f"{D}/phage_host_interactions.csv", index_col=0)
print("=" * 72)
print("interaction matrix")
print("=" * 72)
print(f"raw shape {M.shape}")
# drop all-empty columns and rows
M = M.dropna(axis=1, how="all").dropna(axis=0, how="all")
print(f"after dropping all-empty rows/cols: {M.shape[0]} strains x {M.shape[1]} phages")
cells = M.size
obs = int(M.notna().sum().sum())
pos = int((M == 1).sum().sum())
print(f"cells {cells}   observed {obs} ({100*obs/cells:.1f}%)   "
      f"positives {pos} ({100*pos/obs:.2f}% of observed)")
print(f"values present: {sorted(pd.unique(M.values.ravel()[~pd.isna(M.values.ravel())]))}")

pp = (M == 1).sum(axis=0)
ps = (M == 1).sum(axis=1)
obs_p = M.notna().sum(axis=0)
obs_s = M.notna().sum(axis=1)
print(f"\nper-phage positives: median {pp.median():.0f}  min {pp.min()}  max {pp.max()}")
for k in [5, 10, 20]:
    print(f"  phages with >= {k:>2} positives: {int((pp>=k).sum()):>3} / {len(pp)}")
print(f"per-phage observed cells: median {obs_p.median():.0f}  min {obs_p.min()}  max {obs_p.max()}")
print(f"\nper-strain positives: median {ps.median():.0f}  min {ps.min()}  max {ps.max()}")
print(f"  strains with 0 positives: {int((ps==0).sum())} / {len(ps)}")
print(f"per-strain observed cells: median {obs_s.median():.0f}  min {obs_s.min()}  max {obs_s.max()}")

print("\n" + "=" * 72)
print("Locibase and the K-locus similarity matrix")
print("=" * 72)
for fn in ["Locibase.json", "Locibase_invitro.json"]:
    with open(f"{D}/{fn}") as fh:
        L = json.load(fh)
    ks = list(L.keys())
    print(f"\n{fn}: {len(ks)} strains")
    print(f"  first keys: {ks[:5]}")
    v = L[ks[0]]
    print(f"  value type {type(v).__name__}, "
          f"{len(v) if hasattr(v,'__len__') else '?'} entries for {ks[0]}")
    if isinstance(v, list) and v:
        print(f"  first protein, {len(v[0])} aa: {str(v[0])[:70]}...")
    nprot = [len(L[k]) for k in ks]
    print(f"  proteins per strain: median {np.median(nprot):.0f} "
          f"min {min(nprot)} max {max(nprot)}")

S = pd.read_csv(f"{D}/all_loci_score_matrix.txt", sep="\t", header=None)
print(f"\nall_loci_score_matrix: {S.shape}")
A = S.to_numpy(dtype=float)
print(f"  symmetric: {np.allclose(A, A.T, atol=1e-9)}   "
      f"diagonal all 1: {np.allclose(np.diag(A), 1.0)}")
off = A[np.triu_indices(len(A), 1)]
print(f"  off-diagonal: min {off.min():.3f}  median {np.median(off):.3f}  max {off.max():.3f}")
for t in [0.99, 0.95, 0.90, 0.80]:
    print(f"  pairs with similarity >= {t}: {int((off>=t).sum())} of {len(off)}")

with open(f"{D}/Locibase.json") as fh:
    Lkeys = list(json.load(fh).keys())
print(f"\nLocibase keys {len(Lkeys)} vs score matrix dim {len(A)}: "
      f"{'MATCH' if len(Lkeys)==len(A) else 'MISMATCH'}")

print("\n" + "=" * 72)
print("identifier overlap")
print("=" * 72)
mi = set(M.index)
print(f"matrix strains {len(mi)}   Locibase {len(Lkeys)}   "
      f"overlap {len(mi & set(Lkeys))}")
missing = sorted(mi - set(Lkeys))
print(f"matrix strains absent from Locibase: {len(missing)} {missing[:8]}")

R = pd.read_csv(f"{D}/RBPbase.csv")
rp = set(R["phage_ID"])
mc = set(M.columns)
print(f"\nmatrix phages {len(mc)}   RBPbase phages {len(rp)}   "
      f"overlap {len(mc & rp)}")
print(f"matrix phages with no RBP: {len(mc - rp)} {sorted(mc - rp)[:8]}")
print(f"RBPs per phage: median "
      f"{R.groupby('phage_ID').size().median():.0f} max {R.groupby('phage_ID').size().max()}")

E = pd.read_csv(f"{D}/esm2_embeddings_loci.csv", index_col=0, nrows=3)
print(f"\nesm2_embeddings_loci: {E.shape[1]} columns, index e.g. {list(E.index)[:3]}")
