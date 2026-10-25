"""Fixed-panel, host-macro retrieval metrics for observed binary endpoints.

Exact score ties are averaged over all within-tie permutations. Unknown labels
exclude a whole host from primary metrics, never just the unknown candidates.
"""
from __future__ import annotations

from math import comb

import numpy as np
import pandas as pd

from .data import validate_interactions

DEFAULT_K = (1, 2, 3, 5)


def _identifiers(values, label):
    values = list(values)
    if not values or any(not isinstance(v, str) or not v.strip() for v in values):
        raise ValueError(f"{label} must be non-empty string identifiers")
    if len(set(values)) != len(values):
        raise ValueError(f"{label} must be unique")
    return sorted(values)


def _expected_metrics(labels, scores, ks):
    order = np.argsort(-scores, kind="stable")
    y, s = labels[order], scores[order]
    boundaries = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1, len(s)]
    positives = int(y.sum())
    result = {}
    for k in ks:
        hits, retrieved = 0., 0.
        for start, end in zip(boundaries[:-1], boundaries[1:]):
            n, p = int(end - start), int(y[start:end].sum())
            slots = min(max(k - int(start), 0), n)
            retrieved += slots * p / n
            if p and slots:
                hits = 1. if y[:start].sum() else 1. - comb(n - p, slots) / comb(n, slots)
            if end >= k:
                break
        result[f"hit@{k}"] = hits
        result[f"recall@{k}"] = retrieved / positives if positives else np.nan
    result["mrr"] = 0.
    if positives:
        for start, end in zip(boundaries[:-1], boundaries[1:]):
            n, p = int(end - start), int(y[start:end].sum())
            if p:
                denominator = comb(n, p)
                result["mrr"] = sum(
                    comb(n - j, p - 1) / denominator / (int(start) + j)
                    for j in range(1, n - p + 2)
                )
                break
    return result


def evaluate_retrieval(observed, scores, *, hosts, candidate_phages, ks=DEFAULT_K):
    """Return a support summary and one metric row per *complete-panel* host.

    ``scores`` must contain exactly the host x candidate library Cartesian
    product with finite ``score`` values. It must be made without test labels.
    ``observed`` contains only known labels for that product. The library must
    be fixed independently of evaluation outcomes; no observed-only ranking.
    K larger than the panel uses the whole panel and is reported as such.
    """
    hosts = _identifiers(hosts, "hosts")
    panel = _identifiers(candidate_phages, "candidate_phages")
    ks = tuple(ks)
    if not ks or any(type(k) is not int or k < 1 for k in ks) or len(set(ks)) != len(ks):
        raise ValueError("K values must be unique positive integers")
    ks = tuple(sorted(ks))
    validate_interactions(observed)
    for column in ("bacterium", "phage", "score"):
        if column not in scores:
            raise ValueError(f"missing score column: {column}")
    expected = pd.MultiIndex.from_product([hosts, panel], names=["bacterium", "phage"])
    actual = pd.MultiIndex.from_frame(scores[["bacterium", "phage"]])
    if actual.has_duplicates or set(actual) != set(expected):
        raise ValueError("scores must cover exactly the full fixed host x candidate panel")
    values = scores.set_index(["bacterium", "phage"]).score.reindex(expected).to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("every candidate score must be finite; partial model support is not comparable")
    known = pd.MultiIndex.from_frame(observed[["bacterium", "phage"]])
    if not set(known).issubset(set(expected)):
        raise ValueError("observed labels must belong to the evaluation host x candidate panel")
    labels = observed.set_index(["bacterium", "phage"]).infection.reindex(expected).to_numpy(dtype=float)
    labels, values = labels.reshape(len(hosts), len(panel)), values.reshape(len(hosts), len(panel))
    rows, incomplete = [], []
    for host, y, s in zip(hosts, labels, values):
        if np.isnan(y).any():
            incomplete.append(host)
            continue
        rows.append({"bacterium": host, "positive_candidates": int(y.sum()),
                     "eligible": bool(y.sum()), **_expected_metrics(y, s, ks)})
    columns = ["bacterium", "positive_candidates", "eligible",
               *[name for k in ks for name in (f"hit@{k}", f"recall@{k}")], "mrr"]
    per_host = pd.DataFrame(rows, columns=columns)
    eligible = per_host[per_host.eligible.astype(bool)]
    summary = {
        "held_out_hosts": len(hosts), "candidate_phages": len(panel),
        "observed_pairs": len(observed), "unknown_pairs": int(np.isnan(labels).sum()),
        "complete_panel_hosts": len(per_host), "incomplete_hosts_excluded": incomplete,
        "eligible_hosts": len(eligible), "all_negative_hosts": len(per_host) - len(eligible),
        "complete_panel_host_fraction": len(per_host) / len(hosts),
        "eligible_host_fraction": len(eligible) / len(per_host) if len(per_host) else None,
        "tie_policy": "exact expectation over uniform within-score-tie permutations",
        "effective_k": {str(k): min(k, len(panel)) for k in ks},
        "all_hosts": {}, "eligible_hosts_only": {},
    }
    for name in [*[f"hit@{k}" for k in ks], "mrr"]:
        summary["all_hosts"][name] = float(per_host[name].mean()) if len(per_host) else None
        summary["eligible_hosts_only"][name] = float(eligible[name].mean()) if len(eligible) else None
    for k in ks:
        summary["eligible_hosts_only"][f"recall@{k}"] = (
            float(eligible[f"recall@{k}"].mean()) if len(eligible) else None
        )
    return summary, per_host


def paired_retrieval_bootstrap(left, right, groups, *, metric="hit@2", eligible_only=False,
                               resamples=2000, seed=20261024, confidence=0.95):
    """Paired host-group percentile interval; retain every host in a drawn group.

    Conditions on the fitted models, fixed library and held-out split. Missing
    eligible support remains undefined and its skipped resamples are counted.
    """
    if type(resamples) is not int or resamples < 2 or not 0 < confidence < 1:
        raise ValueError("invalid bootstrap configuration")
    if left.bacterium.duplicated().any() or right.bacterium.duplicated().any():
        raise ValueError("bootstrap requires unique host metric rows")
    if set(left.bacterium) != set(right.bacterium):
        raise ValueError("paired metrics must cover identical hosts")
    a, b = left.set_index("bacterium").sort_index(), right.set_index("bacterium").sort_index()
    if not a.eligible.equals(b.eligible) or not a.positive_candidates.equals(b.positive_candidates):
        raise ValueError("paired metrics must use identical label support")
    if set(a.index) - set(groups):
        raise ValueError("every metric host needs a bootstrap group")
    group_ids = [groups[h] for h in a.index]
    if group_ids:
        _identifiers(list(set(group_ids)), "bootstrap groups")
    names = sorted(set(group_ids))
    if not np.isfinite(a[metric].to_numpy(dtype=float)).all() or not np.isfinite(b[metric].to_numpy(dtype=float)).all():
        raise ValueError("paired metric must be finite")
    delta = (a[metric] - b[metric]).to_numpy()
    keep = a.eligible.to_numpy(dtype=bool) if eligible_only else np.ones(len(a), dtype=bool)
    result = {"metric": metric, "eligible_only": eligible_only, "groups": len(names),
              "hosts": int(keep.sum()), "requested_resamples": resamples, "seed": seed,
              "confidence_level": confidence, "delta": float(delta[keep].mean()) if keep.any() else None}
    if len(names) < 2 or not keep.any():
        return {**result, "status": "insufficient_support", "interval": None,
                "valid_resamples": 0, "undefined_resamples": resamples}
    positions = [np.flatnonzero(np.asarray(group_ids) == name) for name in names]
    rng, draws = np.random.default_rng(seed), []
    for _ in range(resamples):
        indices = np.concatenate([positions[i] for i in rng.integers(len(names), size=len(names))])
        indices = indices[keep[indices]]
        if len(indices):
            draws.append(float(delta[indices].mean()))
    alpha = (1 - confidence) / 2
    return {**result, "status": "conditional_development_interval",
            "interval": np.quantile(draws, [alpha, 1 - alpha]).tolist() if draws else None,
            "valid_resamples": len(draws), "undefined_resamples": resamples - len(draws)}
