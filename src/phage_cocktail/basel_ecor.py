"""Source-locked, read-only intake of GenoPHI's published BASEL/ECOR matrix.

Only audits a source-coded assay phenotype. It never infers infection mechanisms,
converts unknown pairs to negatives, or certifies an independent release holdout.
"""
from __future__ import annotations

from collections import Counter
import csv
from hashlib import sha1, sha256
from io import StringIO
import json
from pathlib import Path
import re

_COLS = ("strain", "phage", "interaction")
_HEX40 = re.compile(r"^[0-9a-f]{40}$")


def load_lock(path: str | Path) -> dict:
    lock = json.loads(Path(path).read_text(encoding="utf-8"))
    if lock.get("schema_version") != 1 or lock.get("status") != "published_inspected_reference":
        raise ValueError("unsupported BASEL/ECOR intake lock/status")
    if not isinstance(lock.get("source"), dict) or not isinstance(lock.get("expected"), dict):
        raise ValueError("missing source/expected fields")
    required_source = ("repository", "commit", "path", "url", "git_blob_sha1", "byte_size")
    for field in required_source:
        if field not in lock["source"]:
            raise ValueError(f"source.{field} is required")
    for field in ("commit", "git_blob_sha1"):
        if not isinstance(lock["source"][field], str) or not _HEX40.fullmatch(lock["source"][field]):
            raise ValueError(f"invalid source.{field}")
    if lock["source"]["byte_size"] <= 0:
        raise ValueError("invalid byte_size")
    for field in ("hosts", "phages", "observed_pairs", "positive", "negative", "unreported_pairs"):
        if type(lock["expected"].get(field)) is not int or lock["expected"][field] < 0:
            raise ValueError(f"invalid expected.{field}")
    if lock["expected"]["positive"] + lock["expected"]["negative"] != lock["expected"]["observed_pairs"]:
        raise ValueError("observed class counts do not sum")
    if lock["expected"]["hosts"] * lock["expected"]["phages"] != (lock["expected"]["observed_pairs"] + lock["expected"]["unreported_pairs"]):
        raise ValueError("expected panel size does not reconcile")
    return lock


def inspect_basel_ecor(payload: bytes, lock: dict) -> dict:
    """Fail closed on unexpected upstream bytes, schema, pair IDs, or counts."""
    source, expected = lock["source"], lock["expected"]
    if len(payload) != source["byte_size"]:
        raise ValueError("source byte-size mismatch")
    blob = sha1(b"blob " + str(len(payload)).encode("ascii") + b"\0" + payload).hexdigest()
    if blob != source["git_blob_sha1"]:
        raise ValueError("source Git blob SHA-1 mismatch")
    text = payload.decode("utf-8-sig")
    reader = csv.DictReader(StringIO(text, newline=""))
    if reader.fieldnames != list(_COLS):
        raise ValueError("source must have exact strain,phage,interaction columns")
    seen: set[tuple[str, str]] = set()
    hosts, phages = set(), set()
    values: Counter[str] = Counter()
    per_host: Counter[str] = Counter()
    for line, row in enumerate(reader, 2):
        if None in row or any(v is None for v in row.values()):
            raise ValueError(f"malformed source row {line}")
        h, p, label = (row[field] for field in _COLS)
        if not h or not p or h != h.strip() or p != p.strip() or label not in {"0", "1"}:
            raise ValueError(f"invalid IDs or source-coded label in row {line}")
        if (h, p) in seen:
            raise ValueError(f"duplicate source pair in row {line}")
        seen.add((h, p))
        hosts.add(h)
        phages.add(p)
        values[label] += 1
        per_host[h] += 1
    observed = len(seen)
    total = len(hosts) * len(phages)
    unknown = total - observed
    counts = dict(hosts=len(hosts), phages=len(phages), observed_pairs=observed,
                  positive=values["1"], negative=values["0"], unreported_pairs=unknown)
    if counts != expected:
        raise ValueError(f"source matrix summary mismatch: observed {counts}, expected {expected}")
    # This is the union-defined Cartesian panel, not a claim all missing pairs were tested.
    full_hosts = sorted(h for h in hosts if per_host[h] == len(phages))
    missing = sorted((h, p) for h in hosts for p in phages if (h, p) not in seen)
    if len(missing) != unknown:
        raise AssertionError("inconsistent missing pair accounting")
    return {
        "status": "source_intake_verified_development_only",
        "independent_release_validation_established": False,
        "source_dataset_id": lock["dataset_id"],
        "upstream_revision": source["commit"],
        "upstream_blob_sha1": blob,
        "local_input_sha256": sha256(payload).hexdigest(),
        "endpoint_id": lock["endpoint_id"],
        "source_summary": counts,
        "complete_panel_host_count": len(full_hosts),
        "complete_panel_host_ids": full_hosts,
        "incomplete_panel_host_count": len(hosts) - len(full_hosts),
        "host_ids": sorted(hosts),
        "phage_ids": sorted(phages),
        "unreported_pair_ids": [{"host_id": h, "phage_id": p} for h, p in missing],
        "cautions": [
            "Unreported pairs are unknown, not source-negative; the publication reports 60 ambiguous assays.",
            "This is an already published and inspected reference, not a new independent release holdout.",
            "Host/phage FASTA content, accession joins, relatedness and label conditions require separate verification.",
            "Source binary interaction scores cannot be relabeled productive infection or cocktail efficacy.",
        ],
    }
