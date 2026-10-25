"""Neighbour tie handling, missing-label semantics and leakage boundaries."""
import numpy as np
import pandas as pd
import pytest

from phage_cocktail.host_similarity import align_core_distances, neighbour_scores, prevalence_shift


def test_boundary_ties_are_fractional_and_order_invariant():
    train = pd.DataFrame({"bacterium": list("ABC") * 2, "phage": ["P"] * 3 + ["Q"] * 3,
                          "infection": [1, 0, 1, 0, 0, 1]})
    candidates = pd.DataFrame({"bacterium": ["D", "D"], "phage": ["P", "Q"]})
    d = pd.DataFrame([[1., 2., 2.]], index=["D"], columns=list("ABC"))
    scores, audit = neighbour_scores(train, candidates, d, k=2)
    np.testing.assert_allclose(scores, [2.5 / 4, 1.5 / 4])
    assert audit.iloc[0].boundary_tied_hosts == 2
    assert audit.iloc[0].selected_train_hosts == "A|B|C"
    other, _ = neighbour_scores(train.sample(frac=1, random_state=4),
                                 candidates.assign(infection=[0, 1]), d, k=2)
    np.testing.assert_array_equal(scores, other)


def test_unknown_neighbour_labels_fall_back_without_choosing_distant_hosts():
    train = pd.DataFrame({"bacterium": ["A", "B", "B"], "phage": ["P", "P", "Q"],
                          "infection": [1, 0, 1]})
    candidates = pd.DataFrame({"bacterium": ["C", "C"], "phage": ["P", "Q"]})
    d = pd.DataFrame([[0., 10.]], index=["C"], columns=["A", "B"])
    scores, audit = neighbour_scores(train, candidates, d, k=1)
    np.testing.assert_allclose(scores, [2/3, 2/3])
    assert audit.iloc[0].fallback_phages == 1
    assert audit.iloc[0].selected_train_hosts == "A"
    with pytest.raises(ValueError, match="disjoint"):
        neighbour_scores(train, candidates.assign(bacterium="A"), d, k=1)
    with pytest.raises(ValueError, match="positive k"):
        neighbour_scores(train, candidates, d, k=0)


def test_distance_alignment_rejects_collisions_nonfinite_and_asymmetry():
    source = pd.DataFrame([[0., 0.2], [0.2, 0.]], index=["a-1", "B"], columns=["a-1", "B"])
    aligned = align_core_distances(source, ["B", "A1"])
    assert aligned.loc["A1", "B"] == 0.2
    source.iloc[0, 1] = np.nan
    with pytest.raises(ValueError, match="symmetric"):
        align_core_distances(source, ["B", "A1"])
    source.iloc[0, 1] = 0.5
    with pytest.raises(ValueError, match="symmetric"):
        align_core_distances(source, ["B", "A1"])
    with pytest.raises(ValueError, match="ambiguous"):
        align_core_distances(pd.DataFrame(np.eye(2), index=["a-1", "A1"], columns=["a-1", "A1"]), ["A1"])


def test_prevalence_shift_retains_unobserved_candidates():
    train = pd.DataFrame({"bacterium": ["A", "B"], "phage": ["P", "P"], "infection": [0, 1]})
    held = pd.DataFrame({"bacterium": ["C"], "phage": ["P"], "infection": [1]})
    result = prevalence_shift(train, held, ["P", "Q"]).set_index("phage")
    assert result.loc["P", "held_minus_train_positive_fraction"] == 0.5
    assert result.loc["Q", "held_observed"] == 0
    assert pd.isna(result.loc["Q", "held_positive_fraction"])
