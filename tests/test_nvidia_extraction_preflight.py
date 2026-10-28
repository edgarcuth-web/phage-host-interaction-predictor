"""Offline preflight of real-plan metadata; no provider or credentials."""
from pathlib import Path
import importlib.util

import pandas as pd
import pytest

from phage_cocktail.artifacts import atomic_write_json
from phage_cocktail import nvidia_host_extraction as ex
# The standalone CLI is not an installed package; load it by its checked-in path.
SCRIPT = Path(__file__).resolve().parents[1] / "scripts/preflight_nvidia_host_extraction.py"
SPEC = importlib.util.spec_from_file_location("preflight_nvidia_host_extraction", SCRIPT)
preflight = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(preflight)
from test_nvidia_host_extraction import extraction


@pytest.fixture
def prepared(extraction, monkeypatch):
    repo, source, store = extraction
    atomic_write_json({"jobs": []}, repo / ex.PILOT)
    plan = ex.make_plan(
        repo, source, selection="original_32_including_short_tails",
        max_new_calls=12, max_new_response_bytes=1200,
    )
    plan["budget"]["max_response_bytes"] = 99

    reference = repo / preflight.REFERENCE
    reference.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(
        {"response": {
            "http_status": 200, "sequence_length": 8192,
            "response_bytes": 1000, "endpoint": plan["endpoint"],
            "layer": plan["layer"],
        }}, reference,
    )
    projected_bytes = round(
        sum(job["length"] for job in plan["jobs"]) * 1000 / 8192
    )
    expected = repo / preflight.PROJECTION
    expected.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "strategy": preflight.SELECTION, "hosts": 3,
        "two_strand_calls": 12,
        "estimated_json_bytes_no_reuse_credit": projected_bytes,
    }]).to_csv(expected, index=False)
    return repo, source, store, plan


def test_complete_plan_cost_preflight_is_offline_and_exact(prepared):
    repo, source, store, plan = prepared
    report = preflight.inspect_plan(repo, plan)
    assert report["status"] == "offline_budget_preflight_only"
    assert report["new_nvidia_requests"] == 0
    assert report["plan_ready_to_run"] is False
    assert report["source_fasta_verification"] == "not_performed"
    assert report["hosts"] == 3
    assert report["selected_windows"] == 6
    assert report["two_strand_jobs"] == 12
    assert report["selected_bases_two_strands"] == 61440
    assert report["estimated_json_bytes_no_reuse_credit"] == 7500
    assert report["conservative_reserved_bytes_no_reuse"] == 1200
    assert report["full_cap_budget_sufficient_no_reuse"] is True
    assert report["potential_pilot_job_matches_unverified"] == 0
    assert len(report["per_host"]) == 3
    assert all(row["windows"] == 2 and row["sampled_fraction"] == 1
               for row in report["per_host"])
    assert store.get("plan.json") is None


def test_insufficient_budget_reported_not_silently_increased(prepared):
    repo, _, _, plan = prepared
    plan["budget"]["max_new_calls"] = 11
    plan["budget"]["max_new_response_bytes"] = 1100
    report = preflight.inspect_plan(repo, plan)
    assert report["full_cap_budget_sufficient_no_reuse"] is False
    assert report["frozen_max_new_calls"] == 11
    assert report["frozen_max_new_response_bytes"] == 1100
    assert report["conservative_reserved_bytes_no_reuse"] == 1200


@pytest.mark.parametrize("change,match", [
    ("partial", "complete frozen host panel"),
    ("selection", "original_32"),
    ("projection", "transfer projection differs"),
    ("jobs", "job count differs"),
])
def test_preflight_rejects_unreviewed_protocol_or_cohort_change(prepared, change, match):
    repo, _, _, plan = prepared
    if change == "partial":
        plan["hosts"] = ["A", "B"]
        plan["jobs"] = [job for job in plan["jobs"] if job["host"] in plan["hosts"]]
        for index, job in enumerate(plan["jobs"]):
            job["id"] = f"{index:05d}"
    elif change == "selection":
        plan["selection"] = "first_full_window"
    else:
        path = repo / preflight.PROJECTION
        table = pd.read_csv(path)
        field = ("estimated_json_bytes_no_reuse_credit"
                 if change == "projection" else "two_strand_calls")
        table.loc[0, field] += 1
        table.to_csv(path, index=False)
    with pytest.raises(ValueError, match=match):
        preflight.inspect_plan(repo, plan)
