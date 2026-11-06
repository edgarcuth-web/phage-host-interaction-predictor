"""Offline retrieval rescue audit: distinct endpoints, ties, unknowns and budgets."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import sha256_file
from phage_cocktail.retrieval_rescue import (
    compare_discovery_routes, run_discovery_audit,
)


HOSTS = ["H1", "H2", "H3"]
PHAGES = ["p1", "p2", "p3", "p4"]
GROUPS = {"H1": "group-A", "H2": "group-B", "H3": "group-C"}


def _panels(*, missing_h3=False, identical=False):
    ids = pd.MultiIndex.from_product(
        [HOSTS, PHAGES], names=["bacterium", "phage"],
    ).to_frame(index=False)
    labels = {
        ("H1", "p3"): 1,
        ("H2", "p1"): 0,  # H2 is an observed all-negative host.
    }
    observed = ids.copy()
    observed["infection"] = [
        labels.get((h, p), 0)
        for h, p in zip(ids.bacterium, ids.phage)
    ]
    if missing_h3:
        observed = observed.loc[
            ~((observed.bacterium == "H3") & (observed.phage == "p4"))
        ]
    evidence = ids.assign(score=[4., 3., 2., 1.] * len(HOSTS))
    learned = ids.assign(
        score=([4., 3., 2., 1.] if identical else [2., 1., 4., 3.]) * len(HOSTS)
    )
    return observed.reset_index(drop=True), {"evidence": evidence, "learned": learned}


def _compare(observed, tables, **kwargs):
    return compare_discovery_routes(
        observed, tables,
        hosts=HOSTS, phages=PHAGES, host_groups=GROUPS,
        ks=(1, 2, 3), union_schedule=("evidence", "learned"),
        bootstrap_resamples=30, bootstrap_seed=3, **kwargs,
    )


def test_learned_route_recovers_positive_discarded_by_evidence_with_equal_k():
    observed, tables = _panels()
    report, rankings, rescue, metrics = _compare(observed, tables)
    row = rescue.set_index(["bacterium", "k"]).loc["H1", 2]
    assert row.learned_only_positive_phages == ["p3"]
    assert row.union_only_positive_phages == ["p3"]
    assert row.union_selected == ["p1", "p3"]
    assert len(row.union_selected) == 2
    assert row.learned_only_positive_count == 1
    assert row.union_only_positive_count == 1
    assert report["primary_k_rescue"]["observed_positive_phages"] == 1
    assert report["primary_k_rescue"]["learned_only_positive_instances"] == 1
    assert report["primary_k_rescue"]["union_only_positive_instances"] == 1
    assert report["method_summaries"]["evidence"]["all_hosts"]["hit@2"] == 0
    assert report["method_summaries"]["union"]["all_hosts"]["hit@2"] == pytest.approx(1 / 3)
    assert set(rankings.method) == {"learned", "evidence", "union"}
    assert len(rankings) == len(HOSTS) * len(PHAGES) * 3
    assert set(metrics.method) == {"learned", "evidence", "union"}


def test_identical_input_orders_have_zero_rescue_and_union_uses_no_extra_slots():
    observed, tables = _panels(identical=True)
    report, _, rescue, _ = _compare(observed, tables)
    assert report["primary_k_rescue"]["learned_only_positive_instances"] == 0
    assert report["primary_k_rescue"]["union_only_positive_instances"] == 0
    for _, row in rescue.iterrows():
        assert len(row.union_selected) == row.k
        assert len(set(row.union_selected)) == row.k
        assert set(row.union_selected) == set(row.evidence_selected)


def test_missing_label_excludes_entire_host_but_keeps_scores_and_audit_rows():
    observed, tables = _panels(missing_h3=True)
    report, ranking, rescue, _ = _compare(observed, tables)
    assert report["unknown_pairs"] == 1
    assert report["excluded_incomplete_hosts"] == ["H3"]
    assert report["method_summaries"]["evidence"]["complete_panel_hosts"] == 2
    assert report["method_summaries"]["evidence"]["all_negative_hosts"] == 1
    assert report["method_summaries"]["evidence"]["eligible_hosts"] == 1
    assert report["primary_k_rescue"]["complete_hosts"] == 2
    assert report["primary_k_rescue"]["observed_positive_phages"] == 1
    assert rescue[(rescue.bacterium == "H3") & (rescue.k == 2)].iloc[0].learned_only_positive_phages is None
    assert len(ranking[ranking.bacterium == "H3"]) == len(PHAGES) * 3
    assert set(ranking.loc[(ranking.bacterium == "H3") & (ranking.phage == "p4"), "label_status"]) == {"unknown"}


def test_ties_are_exact_for_aggregate_but_lexical_for_rescue_membership():
    observed, tables = _panels()
    tables["evidence"]["score"] = 1.
    tables["learned"]["score"] = 1.
    report, _, rescue, _ = _compare(observed, tables)
    assert report["method_summaries"]["evidence"]["all_hosts"]["hit@2"] == pytest.approx(1 / 6)
    assert report["membership_tie_policy"] == "score_descending_then_phage_id_ascending"
    assert rescue.loc[(rescue.bacterium == "H1") & (rescue.k == 2), "evidence_selected"].iloc[0] == ["p1", "p2"]
    assert report["primary_k_rescue"]["learned_only_positive_instances"] == 0


@pytest.mark.parametrize("problem", ["missing", "duplicate", "nonfinite", "extra"])
def test_rejects_invalid_full_panel_scores(problem):
    observed, tables = _panels()
    if problem == "missing":
        tables["learned"] = tables["learned"].iloc[:-1]
    elif problem == "duplicate":
        tables["learned"] = pd.concat([tables["learned"], tables["learned"].iloc[:1]])
    elif problem == "nonfinite":
        tables["learned"].loc[0, "score"] = np.inf
    else:
        tables["learned"].loc[0, "phage"] = "p-other"
    with pytest.raises(ValueError, match="complete fixed panel|finite"):
        _compare(observed, tables)


def test_rejects_invalid_schedule_or_observations():
    observed, tables = _panels()
    with pytest.raises(ValueError, match="union_schedule"):
        compare_discovery_routes(observed, tables, hosts=HOSTS, phages=PHAGES, host_groups=GROUPS,
                                 ks=(2,), union_schedule=("learned",))
    with pytest.raises(ValueError, match="duplicate"):
        _compare(pd.concat([observed, observed.iloc[:1]]), tables)
    extra = observed.copy()
    extra.loc[0, "bacterium"] = "other-host"
    with pytest.raises(ValueError, match="outside the fixed panel"):
        _compare(extra, tables)


def _write_fixture(tmp_path):
    observed, tables = _panels(missing_h3=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    observed.to_csv(inputs / "observed.csv", index=False)
    for name, frame in tables.items():
        frame.to_csv(inputs / f"{name}.csv", index=False)
    endpoint = "test:observed-source-mlc:v1"

    def entry(name, *, producer=None):
        path = inputs / f"{name}.csv"
        value = {"path": f"inputs/{name}.csv", "sha256": sha256_file(path),
                 "endpoint": endpoint}
        if producer is not None:
            value.update(producer=producer, training_hosts=["OTHER-TRAIN-HOST"])
        return value

    manifest = {
        "schema_version": 1, "research_status": "development_only",
        "study_id": "synthetic-comparison-for-software-tests",
        "endpoint": endpoint, "hosts": HOSTS, "phages": PHAGES,
        "host_groups": GROUPS,
        "observations": entry("observed"),
        "score_sources": {
            "evidence": entry("evidence", producer="synthetic-evidence-v1"),
            "learned": entry("learned", producer="synthetic-model-v1"),
        },
        "ks": [1, 2, 3], "union_schedule": ["evidence", "learned"],
        "bootstrap": {"resamples": 30, "seed": 3},
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path, manifest


def test_cli_manifest_publishes_immutable_hash_bound_reports(tmp_path):
    path, manifest = _write_fixture(tmp_path)
    out = tmp_path / "output"
    report = run_discovery_audit(path, out)
    assert report["endpoint"] == manifest["endpoint"]
    assert report["primary_k_rescue"]["learned_only_positive_instances"] == 1
    assert (out / "report.json").exists()
    assert json.loads((out / "report.json").read_text()) == report
    for name, digest in report["output_sha256"].items():
        assert sha256_file(out / name) == digest
    assert report["manifest_sha256"] == sha256_file(path)
    with pytest.raises(ValueError, match="must be new"):
        run_discovery_audit(path, out)


@pytest.mark.parametrize("mutation,match", [
    ("changed_checksum", "changed frozen source bytes"),
    ("endpoint", "different source endpoint"),
    ("path_traversal", "confined"),
    ("leak", "overlap evaluation hosts"),
    ("bad_k", "K values"),
    ("bad_source", "producer"),
])
def test_cli_fails_closed_before_output(tmp_path, mutation, match):
    path, manifest = _write_fixture(tmp_path)
    if mutation == "changed_checksum":
        (tmp_path / "inputs" / "learned.csv").write_text("tampered!")
    elif mutation == "endpoint":
        manifest["score_sources"]["learned"]["endpoint"] = "wrong:endpoint"
    elif mutation == "path_traversal":
        manifest["score_sources"]["learned"]["path"] = "../outside.csv"
    elif mutation == "leak":
        manifest["score_sources"]["learned"]["training_hosts"] = ["H2"]
    elif mutation == "bad_k":
        manifest["ks"] = [5]
    elif mutation == "bad_source":
        manifest["score_sources"]["learned"].pop("producer")
    path.write_text(json.dumps(manifest), encoding="utf-8")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match=match):
        run_discovery_audit(path, output)
    assert not output.exists()
