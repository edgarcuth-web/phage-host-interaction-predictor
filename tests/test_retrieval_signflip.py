from fractions import Fraction
from itertools import product
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import binomtest

from phage_cocktail.retrieval_signflip import (
    paired_group_signflip, holm_adjust, classify_primary_contrasts, classify_prevalence_swings,
)


def paired(left, right):
    def frame(values):
        return pd.DataFrame({"bacterium": [f"H{i}" for i in range(len(values))], "eligible": True,
                             "positive_candidates": 1, "hit@2": values})
    return frame(left), frame(right)


@pytest.mark.parametrize("wins,losses", [(6, 0), (0, 7), (1, 0), (5, 3), (0, 0)])
def test_distinct_discordant_groups_reduce_to_exact_mcnemar(wins, losses):
    a, b = paired([1] * wins + [0] * losses + [1], [0] * wins + [1] * losses + [1])
    groups = {h: h for h in a.bacterium}
    value = paired_group_signflip(a, b, groups)
    expected = binomtest(wins, wins + losses, .5).pvalue if wins + losses else 1.
    assert value["p_value"] == expected
    assert value["non_binary_hosts"] == 0


def test_grouped_binary_hosts_are_not_treated_as_independent_mcnemar_trials():
    a, b = paired([1, 1, 1], [0, 0, 0])
    result = paired_group_signflip(a, b, {"H0": "A", "H1": "A", "H2": "B"})
    assert result["p_value"] == .5
    assert binomtest(3, 3, .5).pvalue == .25
    assert result["discordant_host_groups"] == 2


def test_fractional_hosts_match_brute_force_without_rounding_or_exclusion():
    left, right = [1., 1/3, .25, .75, 0.], [0., .5, .75, .25, 1.]
    a, b = paired(left, right)
    groups = {"H0": "A", "H1": "A", "H2": "B", "H3": "C", "H4": "D"}
    totals = {}
    for h, x, y in zip(a.bacterium, left, right):
        totals[groups[h]] = totals.get(groups[h], Fraction(0)) + Fraction(x) - Fraction(y)
    observed = abs(sum(totals.values()))
    sums = [sum(s * d for s, d in zip(signs, totals.values())) for signs in product([-1, 1], repeat=len(totals))]
    expected = sum(abs(s) >= observed for s in sums) / len(sums)
    result = paired_group_signflip(a, b, groups)
    assert result["p_value"] == expected
    assert result["hosts"] == 5 and result["non_binary_hosts"] == 3
    assert result["left_only"] == result["right_only"] == 1


def test_one_ulp_differences_are_not_thresholded_to_zero():
    a, b = paired([np.nextafter(.5, 1.)] * 8, [.5] * 8)
    result = paired_group_signflip(a, b, {h: h for h in a.bacterium})
    assert result["non_binary_hosts"] == 8
    assert result["nonzero_group_totals"] == 8
    assert result["p_value"] == 2 / 256


def test_eligibility_and_support_are_preserved():
    a, b = paired([1, .5, 0], [0, 0, 0])
    a.loc[2, "eligible"] = b.loc[2, "eligible"] = False
    groups = {h: h for h in a.bacterium}
    assert paired_group_signflip(a, b, groups, eligible_only=True)["hosts"] == 2
    with pytest.raises(ValueError, match="label support"):
        paired_group_signflip(a, b.assign(positive_candidates=2), groups)


def test_registered_decision_is_not_a_majority_vote_or_an_unadjusted_p_value():
    registration = json.loads((Path(__file__).resolve().parents[1] / "data/studies/genophi_seed_sensitivity.v1.json").read_text())
    record = {"status": "verified", "delta": .1, "interval": [.01, .2], "p_value": .03125}
    assert classify_primary_contrasts({43: record, 44: record}, registration)["status"] == "not_established"
    good = {**record, "p_value": .01}
    assert classify_primary_contrasts({43: good, 44: good}, registration)["status"] == "holds_on_prespecified_seeds"
    assert classify_primary_contrasts({43: good, 44: {**good, "delta": -.01}}, registration)["status"] == "does_not_hold_directionally"
    assert classify_primary_contrasts({43: good}, registration)["status"] == "blocked"
    assert holm_adjust([.04, .01]) == [.04, .02]
    assert classify_prevalence_swings({43: -.11, 44: -.12}, registration)["status"] == "reproduced"
    assert classify_prevalence_swings({43: .02, 44: -.03}, registration)["status"] == "seed42_drop_not_reproduced_with_stable_sensitivities"
