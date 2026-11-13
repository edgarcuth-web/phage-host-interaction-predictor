"""Synthetic proof that an identifier join cannot clear phage sequence novelty."""
import copy
import json
from pathlib import Path

import pytest

from phage_cocktail.basel_identity import audit_basel_identity, write_report_new

HEX = "a" * 64


def source(hosts=None, phages=None):
    hosts = hosts or ["ECOR12", "ECOR13"]
    phages = phages or ["Bas01", "Bas02"]
    return {
        "status": "source_intake_verified_development_only",
        "source_dataset_id": "verified-source",
        "upstream_blob_sha1": "f" * 40,
        "independent_release_validation_established": False,
        "host_ids": hosts, "phage_ids": phages,
        "source_summary": {"hosts": len(hosts), "phages": len(phages)},
    }


def resources():
    return {"resources": [
        {"kind": "host_genome", "id": "ECOR12", "sha256": HEX, "url": "https://example.invalid/h1",
         "size_bytes": 30},
        {"kind": "host_genome", "id": "ECOR13", "sha256": "b" * 64,
         "url": "https://example.invalid/h2", "size_bytes": 31},
        {"kind": "phage_genome", "id": "P1", "sha256": "c" * 64},
        {"kind": "table", "id": "source-data"},
    ]}


def test_matching_host_names_never_clears_cold_phage():
    result = audit_basel_identity(source(), resources())
    assert result["host_name_overlap_count"] == 2
    assert result["host_name_unmatched_count"] == 0
    assert result["phage_name_overlap_count"] == 0
    assert result["phage_genome_accessions_missing"] == 2
    assert result["exact_phage_sequence_overlap_result"] == "not_evaluated"
    assert result["cold_phage_claim_allowed"] is False
    assert result["cold_host_claim_allowed"] is False
    assert result["host_ledger"][0]["normalized_sequence_sha256"] is None


def test_shared_phage_name_is_detected_not_promoted():
    result = audit_basel_identity(source(phages=["P1", "Bas02"]), resources())
    assert result["phage_name_overlap_ids"] == ["P1"]
    assert result["phage_genome_accessions_resolved"] == 0
    assert result["cold_phage_claim_allowed"] is False


def test_unmatched_host_is_reported():
    result = audit_basel_identity(source(hosts=["ECOR12", "new-host"]), resources())
    assert result["host_name_unmatched_count"] == 1
    assert result["host_ledger"][-1]["development_source_file_sha256"] is None


@pytest.mark.parametrize("corruption", ["unverified", "yes"])
def test_unverified_source_receipt_blocked(corruption):
    bad = source()
    if corruption == "unverified":
        bad["status"] = "unverified"
    else:
        bad["independent_release_validation_established"] = True
    with pytest.raises(ValueError):
        audit_basel_identity(bad, resources())


def test_incomplete_source_roster_rejected():
    bad = source()
    bad["source_summary"]["hosts"] = 25
    with pytest.raises(ValueError, match="roster"):
        audit_basel_identity(bad, resources())


def test_invalid_resource_hash_rejected():
    bad = copy.deepcopy(resources())
    bad["resources"][0]["sha256"] = "bad"
    with pytest.raises(ValueError, match="hash"):
        audit_basel_identity(source(), bad)


def test_duplicate_dev_resource_ids_rejected():
    bad = copy.deepcopy(resources())
    bad["resources"].append(copy.deepcopy(bad["resources"][0]))
    with pytest.raises(ValueError, match="duplicate"):
        audit_basel_identity(source(), bad)


def test_cannot_overwrite_receipt(tmp_path):
    out = tmp_path / "report.json"
    write_report_new({"status": "blocked"}, out)
    assert json.loads(out.read_text()) == {"status": "blocked"}
    with pytest.raises(FileExistsError):
        write_report_new({"status": "other"}, out)


def test_real_locked_training_host_identifier_coverage():
    lock = json.loads((Path(__file__).resolve().parents[1] /
                       "data/studies/genophi_host_baseline.resources.v1.json").read_text())
    hosts = [f"ECOR{i}" for i in range(12, 36)] + ["ECOR38"]
    audit = audit_basel_identity(source(hosts=hosts), lock)
    assert audit["host_name_overlap_count"] == 25
    assert audit["host_name_unmatched_count"] == 0
    assert all(row["development_source_file_sha256"] for row in audit["host_ledger"])
    assert audit["cold_phage_claim_allowed"] is False
