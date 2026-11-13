"""Synthetic-only checks for frozen partial-label Hit@K semantics."""
import json
from pathlib import Path

import pytest

from phage_cocktail.partial_hit_bounds import bounded_hit_at_k


def test_missing_labels_form_bounds_not_false_negatives():
    r = bounded_hit_at_k(["P1", "P2"], {"P2": 0},
                         {"P1": 0.9, "P2": 0.1}, (1, 2))
    assert r["unknown_count"] == 1
    assert r["observed_negative"] == 1
    assert r["hit_at_k_bounds"]["1"] == {"lower": 0.0, "upper": 1.0}


def test_tied_positive_and_unknown_have_exact_expected_bounds():
    r = bounded_hit_at_k(["P1", "P2"], {"P1": 1},
                         {"P1": 0.4, "P2": 0.4}, (1, 2))
    assert r["hit_at_k_bounds"]["1"] == {"lower": 0.5, "upper": 1.0}
    assert r["hit_at_k_bounds"]["2"] == {"lower": 1.0, "upper": 1.0}


def test_complete_panel_has_point_identified_hit():
    r = bounded_hit_at_k(["P1", "P2"], {"P1": 0, "P2": 1},
                         {"P1": 0.9, "P2": 0.1}, (1, 2))
    assert r["unknown_count"] == 0
    assert r["hit_at_k_bounds"]["1"] == {"lower": 0.0, "upper": 0.0}
    assert r["hit_at_k_bounds"]["2"] == {"lower": 1.0, "upper": 1.0}


def test_all_unknown_is_unidentifiable():
    r = bounded_hit_at_k(["P1", "P2"], {},
                         {"P1": 0.9, "P2": 0.1}, (2,))
    assert r["hit_at_k_bounds"]["2"] == {"lower": 0.0, "upper": 1.0}


def test_incomplete_score_support_fails_closed():
    with pytest.raises(ValueError, match="every frozen candidate"):
        bounded_hit_at_k(["P1", "P2"], {"P1": 1}, {"P1": 0.8}, (1,))


def test_unexpected_label_outside_panel_fails():
    with pytest.raises(ValueError, match="outside"):
        bounded_hit_at_k(["P1"], {"P2": 1}, {"P1": 0.7}, (1,))


@pytest.mark.parametrize("bad", [None, -1, 2, True])
def test_unknown_or_invalid_label_must_be_omitted(bad):
    with pytest.raises(ValueError, match="0 or 1"):
        bounded_hit_at_k(["P1"], {"P1": bad}, {"P1": 0.7}, (1,))


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), "0.4", True])
def test_nonfinite_or_nonnumeric_scores_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        bounded_hit_at_k(["P1"], {}, {"P1": bad}, (1,))


def test_frozen_protocol_contains_no_imputed_negative_rule():
    protocol = json.loads((Path(__file__).resolve().parents[1] /
                           "data/studies/basel_ecor_evaluation_protocol.v1.json").read_text())
    assert protocol["status"] == "protocol_frozen_execution_blocked_on_sequence_identity"
    assert protocol["target"]["primary_K"] == 2
    assert protocol["target"]["denominator"] == "all_25_hosts_including_hosts_with_no_observed_positive_and_hosts_with_unknown_labels"
    assert protocol["target"]["missing"].startswith("unreported_or_ambiguous=unknown")
    assert protocol["target"]["bound_lower"].startswith("replace_unknown_labels_with_negative_only_for_mathematical_lower_bound")
    assert protocol["generalization"]["excluded_claims"] == [
        "cold_host", "double_cold", "independent_untouched_release_validation",
        "productive_infection", "cocktail_efficacy"
    ]
