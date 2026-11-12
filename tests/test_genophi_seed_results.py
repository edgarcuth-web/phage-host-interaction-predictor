"""Committed sensitivity bundles retain fractional support and preregistered decisions."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from phage_cocktail.retrieval_signflip import classify_primary_contrasts, classify_prevalence_swings

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "data/studies/results/genophi_seed_sensitivity.v1"


def test_published_seed_payloads_and_unchanged_seed42_are_hash_bound():
    report = json.loads((ROOT / "report.json").read_text())
    for name, digest in report["output_sha256"].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest
    for seed in (43, 44):
        for stem in ("genophi_comparator", "genome_host_retrieval"):
            root = REPO / f"data/studies/results/{stem}.seed{seed}.v1"
            qa = json.loads((root / "NUMERIC_QA.json").read_text())
            assert qa["status"] == "verified"
            for name, digest in qa["published_output_sha256"].items():
                assert sha256((root / name).read_bytes()).hexdigest() == digest
    registration = json.loads((REPO / "data/studies/genophi_seed_sensitivity.v1.json").read_text())
    source = registration["frozen_comparator"]
    assert sha256((REPO / source["path"]).read_bytes()).hexdigest() == source["sha256"]
    assert sha256((REPO / source["report"]).read_bytes()).hexdigest() == source["report_sha256"]


def test_aggregate_numeric_qa_binds_the_report_and_independent_verifier():
    qa = json.loads((ROOT / "NUMERIC_QA.json").read_text())
    assert qa["status"] == "verified" and qa["additional_model_fits"] == 0
    assert qa["primary_Hit2_p_values_independently_enumerated"] == 6
    assert qa["pairwise_contrasts"] == 360
    assert sha256((ROOT / "report.json").read_bytes()).hexdigest() == qa["report_sha256"]
    assert sha256((REPO / "scripts/verify_genophi_seed_sensitivity.py").read_bytes()).hexdigest() == qa["verification_script_sha256"]


def test_all_requested_contrasts_exist_and_fractional_Hit2_is_retained():
    values = pd.read_csv(ROOT / "paired_contrasts.csv", float_precision="round_trip")
    assert len(values) == 3 * 2 * 6 * 5 * 2
    assert values.non_binary_hosts.notna().all()
    selected = values[values.seed.eq(44) & values.partition.eq("test") & values.metric.eq("hit@2")
                      & values.left.eq("genophi_median") & values.right.eq("phage_prevalence") & ~values.eligible_only]
    row = selected.iloc[0]
    assert row.hosts == 22 and row.non_binary_hosts == 5 and row.binary_hosts == 17
    assert row.left_only == 2 and row.right_only == 0
    assert row.p_value == .1875  # Binary-subset McNemar would incorrectly be .5.


def test_independent_fraction_enumeration_reproduces_real_primary_p_values():
    spec = importlib.util.spec_from_file_location("verify_seed_results", REPO / "scripts/verify_genophi_seed_sensitivity.py")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    per_host = pd.read_csv(ROOT / "per_host_metrics.csv", float_precision="round_trip")
    manifest = pd.read_csv(REPO / "data/studies/results/genophi_host_split.seed44.v1.csv", dtype=str)
    selected = per_host[per_host.seed.eq(44) & per_host.partition.eq("test")]
    assert verifier.independent_group_p(selected[selected.model.eq("genophi_median")],
        selected[selected.model.eq("phage_prevalence")], manifest.set_index("bacterium").host_group.to_dict(),
        metric="hit@2", eligible_only=False) == .1875


def test_decisions_reproduce_registered_rules_without_moving_the_swing_threshold():
    report = json.loads((ROOT / "report.json").read_text())
    registration = json.loads((REPO / "data/studies/genophi_seed_sensitivity.v1.json").read_text())
    evidence = report["primary_decision"]["primary_contrasts"]
    records = {int(seed): {k: v for k, v in value.items() if k != "holm_p_value"} for seed, value in evidence.items()}
    assert classify_primary_contrasts(records, registration) == report["primary_decision"]
    swings = {int(k): v for k, v in report["prevalence_swing_decision"]["sensitivity_swings"].items()}
    assert classify_prevalence_swings(swings, registration) == report["prevalence_swing_decision"]
    assert report["primary_decision"]["status"] == "not_established"
    assert swings[44] > -.10 and swings[44] == pytest.approx(-.09848484848484851, rel=0, abs=1e-15)
    assert report["prevalence_swing_decision"]["status"] == "mixed_or_different"
