"""Exact paired host-group sign randomization, retaining fractional retrieval metrics."""
from bisect import bisect_left, bisect_right
from functools import reduce
from math import gcd

import numpy as np


def _signed_sums(weights):
    sums = [0]
    for weight in weights:
        sums = [value + weight for value in sums] + [value - weight for value in sums]
    return sums


def paired_group_signflip(left, right, groups, *, metric="hit@2", eligible_only=False,
                          max_exact_groups=30):
    """Exact two-sided sign-flip p-value for the host-macro paired difference.

    Each frozen group's SUM of host differences receives one sign. Converting
    binary64 inputs to a shared integer scale retains every fractional host and
    gives exact tail comparisons, including near-zero and equal-statistic ties.
    Meet-in-the-middle counting needs O(2**(groups/2)) storage, not 2**groups.
    """
    if (left.bacterium.duplicated().any() or right.bacterium.duplicated().any()
            or set(left.bacterium) != set(right.bacterium)):
        raise ValueError("sign-flip requires unique, identically supported hosts")
    a, b = left.set_index("bacterium").sort_index(), right.set_index("bacterium").sort_index()
    if not a.eligible.equals(b.eligible) or not a.positive_candidates.equals(b.positive_candidates):
        raise ValueError("paired metrics must use identical label support")
    if set(a.index) - set(groups) or any(not isinstance(groups[h], str) or not groups[h] for h in a.index):
        raise ValueError("every host requires a non-empty frozen group")
    if type(max_exact_groups) is not int or max_exact_groups < 1:
        raise ValueError("invalid exact-enumeration bound")
    x, y = a[metric].to_numpy(dtype=float), b[metric].to_numpy(dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("paired sign-flip metrics must be finite")
    keep = a.eligible.to_numpy(dtype=bool) if eligible_only else np.ones(len(a), dtype=bool)
    names = [groups[h] for h in a.index[keep]]
    x, y = x[keep], y[keep]
    binary = np.isin(x, [0., 1.]) & np.isin(y, [0., 1.])
    result = {
        "metric": metric, "eligible_only": eligible_only, "hosts": len(x),
        "groups": len(set(names)), "binary_hosts": int(binary.sum()),
        "non_binary_hosts": int((~binary).sum()),
        "left_only": int((binary & (x == 1) & (y == 0)).sum()),
        "right_only": int((binary & (x == 0) & (y == 1)).sum()),
        "both_positive": int((binary & (x == 1) & (y == 1)).sum()),
        "both_negative": int((binary & (x == 0) & (y == 0)).sum()),
        "discordance_scope": "Both-binary host pairs only; all selected hosts enter sign-flip and bootstrap",
        "delta": float((x - y).mean()) if len(x) else None,
        "alternative": "two-sided", "unit": "frozen host group",
        "arithmetic": "exact common-denominator integer arithmetic on stored binary64 values",
    }
    if not len(x):
        return {**result, "status": "insufficient_support", "p_value": None,
                "nonzero_group_totals": 0, "sign_assignments": 0}
    pairs = [(float(v).as_integer_ratio(), float(w).as_integer_ratio()) for v, w in zip(x, y)]
    denominator = max(d for pair in pairs for _, d in pair)
    totals = {name: 0 for name in sorted(set(names))}
    discordant_groups = set()
    for name, ((n1, d1), (n2, d2)), is_binary, v, w in zip(names, pairs, binary, x, y):
        totals[name] += n1 * (denominator // d1) - n2 * (denominator // d2)
        if is_binary and v != w:
            discordant_groups.add(name)
    weights = [value for value in totals.values() if value]
    if len(weights) > max_exact_groups:
        raise ValueError("exact group sign-flip bound exceeded; do not silently approximate")
    if weights:
        divisor = reduce(gcd, map(abs, weights))
        weights = [value // divisor for value in weights]
    observed = abs(sum(weights))
    assignments = 1 << len(weights)
    if observed == 0:
        tail = assignments
    else:
        middle = len(weights) // 2
        first, second = _signed_sums(weights[:middle]), sorted(_signed_sums(weights[middle:]))
        tail = sum(len(second) - bisect_left(second, observed - value)
                   + bisect_right(second, -observed - value) for value in first)
    return {
        **result, "status": "exact_conditional_group_signflip", "p_value": tail / assignments,
        "tail_assignments": tail, "sign_assignments": assignments,
        "nonzero_group_totals": len(weights), "zero_total_groups": len(totals) - len(weights),
        "discordant_host_groups": len(discordant_groups),
    }


def holm_adjust(p_values):
    values = list(map(float, p_values))
    if not values or any(not np.isfinite(p) or not 0 <= p <= 1 for p in values):
        raise ValueError("Holm adjustment requires finite p-values in [0,1]")
    order = sorted(range(len(values)), key=lambda i: (values[i], i))
    adjusted, previous = [None] * len(values), 0.
    for position, index in enumerate(order):
        previous = min(1., max(previous, (len(values) - position) * values[index]))
        adjusted[index] = previous
    return adjusted


def classify_primary_contrasts(records, registration):
    """Apply the registered conjunction; no successful-seed subset is eligible."""
    seeds = registration["sensitivity_seeds"]
    if any(s not in records or records[s].get("status") != "verified" for s in seeds):
        return {"status": "blocked", "seed_status": {str(s): records.get(s, {}).get("status", "missing") for s in seeds}}
    adjusted = holm_adjust([records[s]["p_value"] for s in seeds])
    alpha = registration["decision_rule"]["alpha"]
    evidence = {str(s): {**records[s], "holm_p_value": p} for s, p in zip(seeds, adjusted)}
    if any(records[s]["delta"] <= 0 for s in seeds):
        status = "does_not_hold_directionally"
    elif all(records[s]["interval"][0] > 0 and p <= alpha for s, p in zip(seeds, adjusted)):
        status = "holds_on_prespecified_seeds"
    else:
        status = "not_established"
    return {"status": status, "primary_contrasts": evidence, "multiplicity": "Holm across the two new test contrasts"}


def classify_prevalence_swings(swings, registration):
    seeds, rule = registration["sensitivity_seeds"], registration["prevalence_swing_rule"]
    if any(s not in swings or swings[s] is None for s in seeds):
        return {"status": "blocked"}
    if all(swings[s] <= rule["large_drop_threshold"] for s in seeds):
        status = "reproduced"
    elif all(abs(swings[s]) <= rule["stable_absolute_threshold"] for s in seeds):
        status = "seed42_drop_not_reproduced_with_stable_sensitivities"
    else:
        status = "mixed_or_different"
    return {"status": status, "sensitivity_swings": {str(s): swings[s] for s in seeds},
            "reference_seed42": rule["reference_seed42"],
            "large_drop_threshold": rule["large_drop_threshold"],
            "stable_absolute_threshold": rule["stable_absolute_threshold"]}
