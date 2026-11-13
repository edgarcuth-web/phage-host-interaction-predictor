"""Predeclared Hit@K identification bounds with incomplete source labels.

This is a synthetic-testable estimand helper, NOT a model or scientific result.
Lower/upper bounds are logical identification bounds, not confidence intervals.
"""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np

from .retrieval import _expected_metrics


def bounded_hit_at_k(
    candidate_phages: list[str], observed: Mapping[str, int],
    scores: Mapping[str, float], ks: tuple[int, ...] = (2, 1, 3, 5),
) -> dict:
    """Treat absent observations as unknown, not as measured negatives.

    The candidate panel and score support must remain complete, independent of
    the observed-label subset. Exact ties use the retrieval evaluator's uniform
    within-tie permutation expectation.
    """
    panel = list(candidate_phages)
    if not panel or len(panel) != len(set(panel)) or any(not isinstance(x, str) or not x for x in panel):
        raise ValueError("candidate_phages must be unique nonempty strings")
    if set(scores) != set(panel):
        raise ValueError("every frozen candidate must have a score; no observed-only filtering")
    if not set(observed).issubset(panel):
        raise ValueError("observations contain IDs outside the frozen candidate panel")
    if not ks or len(set(ks)) != len(ks) or any(type(k) is not int or k < 1 or k > len(panel) for k in ks):
        raise ValueError("K must be unique positive integers within the fixed panel")
    for v in observed.values():
        if type(v) is not int or v not in (0, 1):
            raise ValueError("observed status must be 0 or 1; omit unknown pairs")

    ordered = sorted(panel)
    values = []
    for p in ordered:
        v = scores[p]
        if type(v) not in (int, float) or not math.isfinite(v):
            raise ValueError("scores must be finite numbers, not probabilities by implication")
        values.append(float(v))
    score_vector = np.asarray(values, dtype=float)
    lower_labels = np.asarray([observed.get(p, 0) for p in ordered], dtype=int)
    upper_labels = np.asarray([observed.get(p, 1) for p in ordered], dtype=int)
    k_values = tuple(sorted(ks))
    lower = _expected_metrics(lower_labels, score_vector, k_values)
    upper = _expected_metrics(upper_labels, score_vector, k_values)
    bounds = {}
    for k in k_values:
        a, b = float(lower[f"hit@{k}"]), float(upper[f"hit@{k}"])
        if a > b + 1e-12:
            raise AssertionError("monotonicity of partial-label bounds violated")
        bounds[str(k)] = {"lower": a, "upper": b}
    return {
        "status": "partial_label_identification_bounds_not_confidence_intervals",
        "candidate_count": len(ordered),
        "observed_positive": sum(v == 1 for v in observed.values()),
        "observed_negative": sum(v == 0 for v in observed.values()),
        "unknown_count": len(ordered) - len(observed),
        "hit_at_k_bounds": bounds,
    }
