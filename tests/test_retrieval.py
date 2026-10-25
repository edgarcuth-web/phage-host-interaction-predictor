"""Ranking semantics: exact ties, fixed panels, denominators and uncertainty."""
from itertools import permutations

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.retrieval import evaluate_retrieval, paired_retrieval_bootstrap


def single(labels, scores, ks=(1, 2, 3, 5)):
    panel = [f"p{i}" for i in range(len(labels))]
    ids = pd.DataFrame({"bacterium": "host", "phage": panel})
    return evaluate_retrieval(ids.assign(infection=labels), ids.assign(score=scores),
                              hosts=["host"], candidate_phages=panel, ks=ks)


@pytest.mark.parametrize("labels,scores", [
    ([0, 1, 0, 1], [1, 1, 1, 1]),
    ([0, 1, 0, 1], [2, 1, 1, 1]),
    ([1, 0, 1, 0], [2, 2, 1, 1]),
    ([0, 0, 1, 0], [3, 2, 1, 1]),
    ([0, 0, 0, 0], [1, 1, 1, 1]),
    ([1, 1, 1, 1], [1, 1, 1, 1]),
])
def test_tie_expectations_match_exhaustive_rank_permutations(labels, scores):
    summary, _ = single(labels, scores)
    valid = [order for order in permutations(range(4))
             if all(scores[order[i]] >= scores[order[i + 1]] for i in range(3))]
    for k in (1, 2, 3, 5):
        hit = np.mean([any(labels[i] for i in order[:k]) for order in valid])
        assert summary["all_hosts"][f"hit@{k}"] == pytest.approx(hit)
        if sum(labels):
            recall = np.mean([sum(labels[i] for i in order[:k]) / sum(labels) for order in valid])
            assert summary["eligible_hosts_only"][f"recall@{k}"] == pytest.approx(recall)
    rr = [next((1 / (i + 1) for i, j in enumerate(order) if labels[j]), 0.) for order in valid]
    assert summary["all_hosts"]["mrr"] == pytest.approx(np.mean(rr))
    assert summary["effective_k"]["5"] == 4


def test_unknown_host_is_excluded_whole_and_negative_host_remains_in_denominator():
    ids = pd.MultiIndex.from_product([list("ABC"), list("pq")],
                                     names=["bacterium", "phage"]).to_frame(index=False)
    known = ids.assign(infection=[0, 1, 0, 0, 0, 0]).iloc[:-1]
    summary, rows = evaluate_retrieval(known, ids.assign(score=[2, 1] * 3),
                                      hosts=list("ABC"), candidate_phages=list("pq"))
    assert summary["incomplete_hosts_excluded"] == ["C"]
    assert summary["unknown_pairs"] == 1
    assert summary["complete_panel_hosts"] == 2
    assert summary["all_negative_hosts"] == 1
    assert summary["eligible_hosts"] == 1
    assert summary["all_hosts"]["hit@2"] == 0.5
    assert summary["eligible_hosts_only"]["hit@2"] == 1.
    assert np.isnan(rows.set_index("bacterium").loc["B", "recall@1"])


def test_no_complete_or_no_eligible_support_is_null_not_success():
    summary, _ = single([0, 0], [1, 1])
    assert summary["all_hosts"]["mrr"] == 0
    assert summary["eligible_hosts_only"]["hit@1"] is None
    scores = pd.DataFrame({"bacterium": ["H"], "phage": ["P"], "score": [1.]})
    observed = scores.iloc[:0].drop(columns="score").assign(infection=[])
    summary, rows = evaluate_retrieval(observed, scores, hosts=["H"], candidate_phages=["P"])
    assert summary["all_hosts"]["hit@1"] is None
    assert rows.empty
    interval = paired_retrieval_bootstrap(rows, rows, {})
    assert interval["status"] == "insufficient_support"
    assert interval["delta"] is None


@pytest.mark.parametrize("change", ["missing", "duplicate", "extra", "nan", "inf"])
def test_partial_or_invalid_score_support_is_rejected(change):
    ids = pd.DataFrame({"bacterium": ["H", "H"], "phage": ["P", "Q"]})
    scores = ids.assign(score=[1., 0.])
    if change == "missing": scores = scores.iloc[:1]
    if change == "duplicate": scores = pd.concat([scores, scores.iloc[:1]])
    if change == "extra": scores = pd.concat([scores, scores.iloc[:1].assign(phage="EXTRA")])
    if change == "nan": scores.loc[0, "score"] = np.nan
    if change == "inf": scores.loc[0, "score"] = np.inf
    with pytest.raises(ValueError, match="full fixed|finite"):
        evaluate_retrieval(ids.assign(infection=[1, 0]), scores, hosts=["H"], candidate_phages=["P", "Q"])


def test_row_and_identifier_order_cannot_break_score_ties():
    ids = pd.MultiIndex.from_product([list("ABC"), list("pqr")],
                                     names=["bacterium", "phage"]).to_frame(index=False)
    labels, scores = ids.assign(infection=[0, 1, 0] * 3), ids.assign(score=1.)
    a, ar = evaluate_retrieval(labels, scores, hosts=list("ABC"), candidate_phages=list("pqr"))
    b, br = evaluate_retrieval(labels.sample(frac=1, random_state=1), scores.sample(frac=1, random_state=2),
                              hosts=list("CBA"), candidate_phages=list("rqp"))
    assert a == b
    pd.testing.assert_frame_equal(ar, br)
    assert a["all_hosts"]["hit@2"] == pytest.approx(2 / 3)


def test_group_bootstrap_keeps_hosts_together_and_pairs_models():
    # Two hosts in group X, one in Y. Group draws XX, XY, YY give 1, 2/3, 0.
    a = pd.DataFrame({"bacterium": list("ABC"), "positive_candidates": 1,
                      "eligible": True, "hit@2": [1., 1., 0.]})
    b = a.assign(**{"hit@2": 0.})
    groups = dict(zip("ABC", ["X", "X", "Y"]))
    result = paired_retrieval_bootstrap(a, b, groups, resamples=100, seed=7)
    assert result["delta"] == pytest.approx(2 / 3)
    assert result["interval"] == [0., 1.]
    assert result["valid_resamples"] == 100
    same = paired_retrieval_bootstrap(a, a.sample(frac=1, random_state=4), groups, resamples=100)
    assert same["interval"] == [0., 0.]
    with pytest.raises(ValueError, match="identical hosts"):
        paired_retrieval_bootstrap(a, b.iloc[:-1], groups)


def test_eligible_bootstrap_reports_undefined_resamples():
    a = pd.DataFrame({"bacterium": list("AB"), "positive_candidates": [1, 0],
                      "eligible": [True, False], "hit@2": [1., 0.]})
    b = a.assign(**{"hit@2": 0.})
    result = paired_retrieval_bootstrap(a, b, dict(zip("AB", "XY")),
                                        eligible_only=True, resamples=100, seed=7)
    assert result["delta"] == 1.
    assert result["undefined_resamples"] > 0
    assert result["undefined_resamples"] + result["valid_resamples"] == 100
    result = paired_retrieval_bootstrap(a, b, {"A": "X", "B": "X"})
    assert result["status"] == "insufficient_support"
    assert result["interval"] is None
