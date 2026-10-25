"""Train-label-only neighbour transfer and descriptive domain-shift diagnostics."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .cohort import host_aliases, host_key
from .data import validate_interactions


def align_core_distances(source, hosts):
    """Align the pinned source matrix; never fuzzy-match host identities."""
    hosts = sorted(hosts)
    if len(set(hosts)) != len(hosts) or not hosts:
        raise ValueError("distance hosts must be non-empty and unique")
    rows, cols = host_aliases(source.index), host_aliases(source.columns)
    try:
        values = source.loc[[rows[host_key(h)] for h in hosts],
                            [cols[host_key(h)] for h in hosts]].to_numpy(dtype=float)
    except KeyError as exc:
        raise ValueError("missing core distance host") from exc
    if (not np.isfinite(values).all() or (values < 0).any()
            or not np.allclose(values, values.T, rtol=0, atol=1e-10)
            or not np.allclose(np.diag(values), 0, rtol=0, atol=1e-10)):
        raise ValueError("core distances must be finite, non-negative, symmetric with zero diagonal")
    return pd.DataFrame(values, index=hosts, columns=hosts)


def neighbour_scores(train, candidates, distances, *, k=5):
    """Fixed-k fractional boundary weights; missing labels never become negatives.

    Choose neighbours once per target, without consulting any candidate label.
    For ties at the kth distance, distribute the remaining neighbour mass evenly.
    Each phage uses Beta(1,1) smoothing over observed weighted labels, falling back
    to its train-only prevalence when all selected labels are unknown.
    """
    validate_interactions(train)
    if type(k) is not int or k < 1 or train.empty:
        raise ValueError("neighbour transfer needs positive k and non-empty training")
    train_hosts = sorted(train.bacterium.unique())
    targets = sorted(candidates.bacterium.unique())
    if set(train_hosts) & set(targets):
        raise ValueError("training and target hosts must be disjoint")
    if set(candidates.phage) - set(train.phage):
        raise ValueError("candidate has no training phage support")
    if set(targets) - set(distances.index) or set(train_hosts) - set(distances.columns):
        raise ValueError("missing core distance host")
    panel = sorted(train.phage.unique())
    labels = train.pivot(index="bacterium", columns="phage", values="infection").reindex(
        index=train_hosts, columns=panel).to_numpy(dtype=float)
    known = np.isfinite(labels)
    stats = train.groupby("phage").infection.agg(["sum", "count"])
    prevalence = ((stats["sum"] + 1) / (stats["count"] + 2)).reindex(panel).to_numpy()
    scores, diagnostics = {}, []
    for target in targets:
        d = distances.loc[target, train_hosts].to_numpy(dtype=float)
        if not np.isfinite(d).all() or (d < 0).any():
            raise ValueError("invalid target-to-train distances")
        effective_k = min(k, len(d))
        boundary = np.sort(d)[effective_k - 1]
        closer, tied = d < boundary, d == boundary
        weights = closer.astype(float) + tied * ((effective_k - closer.sum()) / tied.sum())
        support = weights @ known
        positive = weights @ np.nan_to_num(labels, nan=0.)
        values = np.where(support > 0, (positive + 1) / (support + 2), prevalence)
        scores[target] = dict(zip(panel, values))
        diagnostics.append({
            "bacterium": target, "nearest_train_distance": float(d.min()),
            "nearest_train_hosts": "|".join(h for h, v in zip(train_hosts, d) if v == d.min()),
            "neighbour_boundary_distance": float(boundary), "effective_k": effective_k,
            "selected_train_hosts": "|".join(h for h, w in zip(train_hosts, weights) if w > 0),
            "boundary_tied_hosts": int(tied.sum()),
            "fallback_phages": int((support == 0).sum()),
        })
    return (np.array([scores[r.bacterium][r.phage] for r in candidates.itertuples(index=False)]),
            pd.DataFrame(diagnostics))


def prevalence_shift(train, held, panel):
    """Post-scoring descriptive label drift. Never feed this into a fitted model."""
    result = pd.DataFrame(index=pd.Index(panel, name="phage"))
    for name, pairs in (("train", train), ("held", held)):
        stats = pairs.groupby("phage").infection.agg(["sum", "count"])
        result[f"{name}_observed"] = stats["count"].reindex(panel, fill_value=0)
        result[f"{name}_positive"] = stats["sum"].reindex(panel, fill_value=0)
        result[f"{name}_positive_fraction"] = result[f"{name}_positive"] / result[f"{name}_observed"].replace(0, np.nan)
    result["held_minus_train_positive_fraction"] = result.held_positive_fraction - result.train_positive_fraction
    return result.reset_index()
