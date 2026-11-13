"""Identity/lineage release-gate tests use invented sequences and identifiers only."""
import csv

import pytest

from phage_cocktail.external_overlap import COLUMNS, audit_overlap, read_intake


def record(**overrides):
    base = {
        "host_id": "H1", "phage_id": "P1", "host_assembly_sha256": "a" * 64,
        "phage_assembly_sha256": "b" * 64, "source_study_id": "study-A",
        "source_observation_id": "reading-1", "endpoint_id": "source-MLC-v1",
        "observation_status": "observed_positive",
    }
    return dict(base, **overrides)


def load(tmp_path, name, rows):
    file = tmp_path / name
    with file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return read_intake(file)


def test_same_underlying_study_blocks_repackaged_source(tmp_path):
    a = load(tmp_path, "a.csv", [record()])
    b = load(tmp_path, "b.csv", [record(source_observation_id="reading-2", host_id="X")])
    report = audit_overlap(a, b, regime="cold_host")
    assert report["status"] == "blocked"
    assert "shared_primary_source_lineage" in report["blocking_reasons"]
    assert "reused_source_observations" not in report["blocking_reasons"]


def test_exact_sequence_pair_blocks_even_with_different_accessions(tmp_path):
    a = load(tmp_path, "a.csv", [record()])
    b = load(tmp_path, "b.csv", [record(host_id="H2", phage_id="P2", source_study_id="study-B",
                                        source_observation_id="different", endpoint_id="spot-v2")])
    report = audit_overlap(a, b, regime="cold_host")
    assert report["exact_overlap_counts"]["shared_exact_host_phage_pairs"] == 1
    assert "exact_pair_overlap" in report["blocking_reasons"]
    assert report["independent_release_validation_established"] is False


def test_cold_host_allows_known_phage_overlap_but_not_claim_of_independence(tmp_path):
    a = load(tmp_path, "a.csv", [record()])
    b = load(tmp_path, "b.csv", [record(host_id="H2", host_assembly_sha256="c" * 64,
                                        source_study_id="study-B", source_observation_id="reading-2")])
    report = audit_overlap(a, b, regime="cold_host")
    assert report["exact_overlap_counts"]["shared_exact_phage_assemblies"] == 1
    assert report["status"] == "requires_manual_review"
    assert report["independent_release_validation_established"] is False
    assert audit_overlap(a, b, regime="double_cold")["status"] == "blocked"


def test_missing_assembly_blocks_and_unknown_label_not_negative(tmp_path):
    a = load(tmp_path, "a.csv", [record()])
    b = load(tmp_path, "b.csv", [record(host_assembly_sha256="", source_study_id="study-B",
                                        observation_status="unknown", source_observation_id="")])
    report = audit_overlap(a, b, regime="cold_host")
    assert "incomplete_assembly_identity_evidence" in report["blocking_reasons"]
    assert report["missingness"]["candidate"]["unknown_observation_status"] == 1
    assert report["missingness"]["candidate"]["missing_source_observation_id"] == 1


@pytest.mark.parametrize("bad", [
    {"host_id": ""},
    {"host_assembly_sha256": "ABC"},
    {"phage_assembly_sha256": "f" * 63},
    {"observation_status": "negative"},
])
def test_invalid_rows_fail_closed(tmp_path, bad):
    with pytest.raises(ValueError):
        load(tmp_path, "invalid.csv", [record(**bad)])


def test_duplicate_input_row_fails_closed(tmp_path):
    with pytest.raises(ValueError, match="duplicate row"):
        load(tmp_path, "duplicate.csv", [record(), record()])


def test_same_assay_record_reuse_detected(tmp_path):
    a = load(tmp_path, "a.csv", [record()])
    b = load(tmp_path, "b.csv", [record(host_assembly_sha256="c" * 64,
                                        phage_assembly_sha256="d" * 64,
                                        endpoint_id="other", source_study_id="study-A")])
    report = audit_overlap(a, b, regime="double_cold")
    assert report["exact_overlap_counts"]["shared_source_observations"] == 1
    assert "reused_source_observations" in report["blocking_reasons"]
