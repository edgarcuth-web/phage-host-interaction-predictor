"""Sequence-cache identity, host-conditioned ranking and missing-label semantics."""
from hashlib import sha256
from itertools import permutations

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.basel_controls import (
    cosine_panel, evaluate_control, expected_slot_fraction, sequence_feature, validate_feature,
)
from phage_cocktail.basel_sequences import assembly_digest
from phage_cocktail.cohort import parse_dna_fasta

SPEC = sha256(b"frozen ordered 4-mer representation").hexdigest()


def feature(payload=b">strain\nACGTTACGTTNNNNACGT\n"):
    assembly = assembly_digest(parse_dna_fasta(payload.decode()))
    vector, receipt = sequence_feature(payload, assembly, SPEC)
    return vector, receipt, assembly


def test_same_host_name_does_not_admit_different_assembly():
    vector, receipt, assembly = feature()
    validate_feature(vector, receipt, assembly, SPEC)
    other = assembly_digest(["TTTTTTTT"])
    with pytest.raises(ValueError, match="assembly/representation"):
        validate_feature(vector, receipt, other, SPEC)
    with pytest.raises(ValueError, match="host name"):
        sequence_feature(b">strain\nTTTTTTTT\n", assembly, SPEC)


def test_different_representation_or_changed_vector_cannot_use_cache():
    vector, receipt, assembly = feature()
    with pytest.raises(ValueError, match="assembly/representation"):
        validate_feature(vector, receipt, assembly, sha256(b"different column order").hexdigest())
    changed = vector.copy()
    positions = np.flatnonzero(changed)
    changed[positions[0]] += .001
    changed[positions[1]] -= .001
    with pytest.raises(ValueError, match="receipt"):
        validate_feature(changed, receipt, assembly, SPEC)


def test_feature_content_identity_ignores_headers_and_strand_but_keeps_ambiguity():
    vector, receipt, assembly = feature(b">old\nACGTTNNACGT\n")
    other, other_receipt = sequence_feature(b">renamed\nACGTNNAACGT\n", assembly, SPEC)
    assert np.array_equal(vector, other)
    assert receipt["feature_key"] == other_receipt["feature_key"]
    assert receipt["ambiguous_bases"] == 2
    with pytest.raises(ValueError, match="no valid"):
        feature(b">strain\nNNNNNNNN\n")


def test_composition_control_can_change_phage_order_with_the_host():
    scores = cosine_panel({"h1": [1., 0.], "h2": [0., 1.]}, {"p1": [1., 0.], "p2": [0., 1.]})
    top = scores.sort_values("score", ascending=False).groupby("bacterium").first().phage.to_dict()
    assert top == {"h1": "p1", "h2": "p2"}
    assert len(scores) == 4
    # Identical phage features retain exact score ties; no ID-dependent perturbation.
    tied = cosine_panel({"h": [.4, .6]}, {"p1": [.2, .8], "p2": [.2, .8]})
    assert tied.score.nunique() == 1


@pytest.mark.parametrize("bad", [[0., 0.], [np.nan, 1.], [-1., 1.], [1.]])
def test_unsupported_sequence_vector_blocks_complete_scoring(bad):
    with pytest.raises(ValueError, match="cosine"):
        cosine_panel({"h": [1., 0.]}, {"p": bad})


def tied_source():
    observed = pd.DataFrame([dict(bacterium="h", phage="p1", infection=1),
                             dict(bacterium="h", phage="p2", infection=0)])
    scores = pd.DataFrame([dict(bacterium="h", phage=p, score=0.) for p in ("p1", "p2", "p3")])
    return observed, scores


def test_tie_bounds_and_capped_costs_match_exhaustive_unknown_assignments():
    observed, scores = tied_source()
    summary, rows = evaluate_control(observed, scores, ["h"], ["p1", "p2", "p3"], ks=(1, 2, 3))
    for k in (1, 2, 3):
        hit, cost = [], []
        for unknown in (0, 1):
            hvalues, cvalues = [], []
            labels = {"p1": 1, "p2": 0, "p3": unknown}
            for order in permutations(labels):
                first = next((j for j, p in enumerate(order[:k], 1) if labels[p]), None)
                hvalues.append(float(first is not None))
                cvalues.append(first if first is not None else k)
            hit.append(np.mean(hvalues)); cost.append(np.mean(cvalues))
        m = summary["macro"][str(k)]
        assert m["hit_lower"] == pytest.approx(min(hit))
        assert m["hit_upper"] == pytest.approx(max(hit))
        assert m["tests_lower"] == pytest.approx(min(cost))
        assert m["tests_upper"] == pytest.approx(max(cost))
        assert m["unknown_slot_fraction"] == pytest.approx(1 / 3)


def test_all_negative_host_stays_in_denominator_and_unknown_is_not_negative():
    observed = pd.DataFrame([dict(bacterium="h1", phage="p2", infection=1),
                             dict(bacterium="h2", phage="p1", infection=0),
                             dict(bacterium="h2", phage="p2", infection=0)])
    scores = pd.DataFrame([dict(bacterium=h, phage=p, score=s)
                           for h in ("h1", "h2") for p, s in (("p1", 1.), ("p2", 0.))])
    summary, rows = evaluate_control(observed, scores, ["h1", "h2"], ["p1", "p2"], ks=(1, 2))
    assert summary["hosts"] == 2 and summary["unknown_pairs"] == 1
    assert summary["macro"]["1"]["hit_lower"] == 0.
    assert summary["macro"]["1"]["hit_upper"] == .5
    assert summary["macro"]["2"]["hit_lower"] == .5
    assert rows.query("bacterium == 'h1'").observed_negative.eq(0).all()


def test_observed_only_candidate_filtering_and_duplicate_scores_are_rejected():
    observed, scores = tied_source()
    with pytest.raises(ValueError, match="complete"):
        evaluate_control(observed, scores.iloc[:2], ["h"], ["p1", "p2", "p3"], ks=(1,))
    with pytest.raises(ValueError, match="complete"):
        evaluate_control(observed, pd.concat([scores, scores.iloc[:1]]), ["h"], ["p1", "p2", "p3"], ks=(1,))


@pytest.mark.parametrize("bad", [None, pd.NA, True, 2])
def test_unknown_or_invalid_label_cannot_become_measured_negative(bad):
    observed, scores = tied_source()
    observed["infection"] = observed.infection.astype(object)
    observed.loc[0, "infection"] = bad
    with pytest.raises(ValueError, match="omit unknown"):
        evaluate_control(observed, scores, ["h"], ["p1", "p2", "p3"], ks=(1,))


def test_marked_candidate_exposure_uses_slots_not_probability_of_any_hit():
    assert expected_slot_fraction({"a": 1., "b": 1., "c": 1.}, {"a"}, 2) == pytest.approx(1 / 3)
    assert expected_slot_fraction({"a": 2., "b": 1., "c": 0.}, {"a"}, 2) == .5
