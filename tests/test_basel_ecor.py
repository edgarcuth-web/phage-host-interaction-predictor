"""Synthetic offline tests of the published BASEL/ECOR source intake."""
import csv
from hashlib import sha1
from io import StringIO
from pathlib import Path

import pytest

from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock


def fixture(rows):
    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("strain", "phage", "interaction"))
    writer.writerows(rows)
    payload = output.getvalue().encode()
    blob = sha1(b"blob " + str(len(payload)).encode() + b"\0" + payload).hexdigest()
    lock = {
        "schema_version": 1, "dataset_id": "synthetic", "status": "published_inspected_reference",
        "endpoint_id": "synthetic-assay", "source": {
            "repository": "example/repo", "commit": "f" * 40, "path": "test.csv",
            "url": "https://example.invalid/test.csv", "git_blob_sha1": blob,
            "byte_size": len(payload)},
        "expected": {"hosts": 2, "phages": 2, "observed_pairs": len(rows),
                     "positive": sum(r[2] == "1" for r in rows),
                     "negative": sum(r[2] == "0" for r in rows),
                     "unreported_pairs": 4 - len(rows)},
    }
    return payload, lock


def test_missing_pairs_are_unreported_not_negative():
    payload, lock = fixture([("h1", "p1", "1"), ("h1", "p2", "0"), ("h2", "p1", "1")])
    result = inspect_basel_ecor(payload, lock)
    assert result["source_summary"]["positive"] == 2
    assert result["source_summary"]["negative"] == 1
    assert result["source_summary"]["unreported_pairs"] == 1
    assert result["unreported_pair_ids"] == [{"host_id": "h2", "phage_id": "p2"}]
    assert not result["independent_release_validation_established"]
    assert result["complete_panel_host_ids"] == ["h1"]
    assert result["incomplete_panel_host_count"] == 1


def test_tampered_bytes_rejected_even_if_labels_look_plausible():
    payload, lock = fixture([("h1", "p1", "1"), ("h1", "p2", "0"), ("h2", "p1", "1")])
    with pytest.raises(ValueError, match="SHA-1 mismatch"):
        inspect_basel_ecor(payload.replace(b"h1", b"x1"), lock)


def test_changed_size_rejected():
    payload, lock = fixture([("h1", "p1", "1"), ("h1", "p2", "0"), ("h2", "p1", "1")])
    with pytest.raises(ValueError, match="byte-size mismatch"):
        inspect_basel_ecor(payload + b"x", lock)


def test_duplicates_fail_even_with_consistent_lock():
    payload, lock = fixture([("h1", "p1", "1"), ("h1", "p2", "0"),
                             ("h2", "p1", "1"), ("h1", "p1", "1")])
    with pytest.raises(ValueError, match="duplicate source pair"):
        inspect_basel_ecor(payload, lock)


def test_unreviewed_category_not_silently_mapped_to_negative():
    payload, lock = fixture([("h1", "p1", "1"), ("h1", "p2", "2"), ("h2", "p1", "0")])
    with pytest.raises(ValueError, match="invalid IDs or source-coded label"):
        inspect_basel_ecor(payload, lock)


def test_cannot_change_locked_expected_counts():
    payload, lock = fixture([("h1", "p1", "1"), ("h1", "p2", "0"), ("h2", "p1", "1")])
    lock["expected"]["negative"] = 2
    with pytest.raises(ValueError, match="summary mismatch"):
        inspect_basel_ecor(payload, lock)


def test_malformed_or_whitespace_ids_are_not_normalized():
    payload, lock = fixture([("h1 ", "p1", "1"), ("h1", "p2", "0"), ("h2", "p1", "1")])
    with pytest.raises(ValueError, match="invalid IDs"):
        inspect_basel_ecor(payload, lock)


def test_pinned_source_lock_is_well_formed():
    path = Path(__file__).resolve().parents[1] / "data/studies/basel_ecor_intake.v1.json"
    lock = load_lock(path)
    assert lock["expected"] == {
        "hosts": 25, "phages": 52, "observed_pairs": 1240,
        "positive": 302, "negative": 938, "unreported_pairs": 60,
    }
    assert lock["source"]["git_blob_sha1"] == "dcbbca4b4e4b1dd0bc543c5e2348189813e55bb0"
    assert lock["scientific_scope"]["release_validation"] == "not_established"
