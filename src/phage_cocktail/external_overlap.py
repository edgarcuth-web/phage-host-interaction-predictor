"""Conservative identity/lineage checks before claiming an external evaluation.

This checks exact evidence only. A clean result is NOT independent validation.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re

_SHA = re.compile(r"[0-9a-f]{64}\Z")
COLUMNS = (
    "host_id", "phage_id", "host_assembly_sha256", "phage_assembly_sha256",
    "source_study_id", "source_observation_id", "endpoint_id", "observation_status",
)
STATUSES = frozenset({"observed_positive", "observed_negative", "unknown"})
REGIMES = frozenset({"cold_host", "cold_phage", "double_cold"})


@dataclass(frozen=True)
class IntakeRow:
    host_id: str
    phage_id: str
    host_assembly_sha256: str
    phage_assembly_sha256: str
    source_study_id: str
    source_observation_id: str
    endpoint_id: str
    observation_status: str


def read_intake(path: str | Path) -> list[IntakeRow]:
    """Read an explicit source-row ledger without imputing identities or labels."""
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None or set(reader.fieldnames) != set(COLUMNS):
            raise ValueError(f"{path}: expected exactly these columns: {', '.join(COLUMNS)}")
        rows: list[IntakeRow] = []
        seen = set()
        for lineno, source in enumerate(reader, 2):
            if None in source or any(value is None for value in source.values()):
                raise ValueError(f"{path}:{lineno}: malformed CSV row")
            row = IntakeRow(**{field: source[field].strip() for field in COLUMNS})
            for field in ("host_id", "phage_id", "source_study_id", "endpoint_id"):
                if not getattr(row, field):
                    raise ValueError(f"{path}:{lineno}: {field} is required")
            for field in ("host_assembly_sha256", "phage_assembly_sha256"):
                value = getattr(row, field)
                if value and not _SHA.fullmatch(value):
                    raise ValueError(f"{path}:{lineno}: invalid {field}; use lowercase hex SHA-256")
            if row.observation_status not in STATUSES:
                raise ValueError(f"{path}:{lineno}: invalid observation_status")
            if row in seen:
                raise ValueError(f"{path}:{lineno}: duplicate row")
            seen.add(row)
            rows.append(row)
    if not rows:
        raise ValueError(f"{path}: empty source-row ledger")
    return rows


def _nonblank(rows: list[IntakeRow], field: str) -> set[str]:
    return {value for row in rows if (value := getattr(row, field))}


def _pairs(rows: list[IntakeRow]) -> set[tuple[str, str]]:
    return {(r.host_assembly_sha256, r.phage_assembly_sha256)
            for r in rows if r.host_assembly_sha256 and r.phage_assembly_sha256}


def _observations(rows: list[IntakeRow]) -> set[tuple[str, str]]:
    return {(r.source_study_id, r.source_observation_id)
            for r in rows if r.source_observation_id}


def audit_overlap(
    development: list[IntakeRow], candidate: list[IntakeRow], *, regime: str
) -> dict:
    """Block known overlap, never certify independence from exact-hash checks alone."""
    if regime not in REGIMES:
        raise ValueError(f"regime must be one of {', '.join(sorted(REGIMES))}")
    if not development or not candidate:
        raise ValueError("both ledgers must contain at least one row")

    overlaps = {
        "shared_source_studies": len(_nonblank(development, "source_study_id") &
                                     _nonblank(candidate, "source_study_id")),
        "shared_source_observations": len(_observations(development) & _observations(candidate)),
        "shared_exact_host_assemblies": len(_nonblank(development, "host_assembly_sha256") &
                                            _nonblank(candidate, "host_assembly_sha256")),
        "shared_exact_phage_assemblies": len(_nonblank(development, "phage_assembly_sha256") &
                                             _nonblank(candidate, "phage_assembly_sha256")),
        "shared_exact_host_phage_pairs": len(_pairs(development) & _pairs(candidate)),
    }
    missing = {
        role: {
            "missing_host_assembly_sha256": sum(not r.host_assembly_sha256 for r in rows),
            "missing_phage_assembly_sha256": sum(not r.phage_assembly_sha256 for r in rows),
            "missing_source_observation_id": sum(not r.source_observation_id for r in rows),
            "unknown_observation_status": sum(r.observation_status == "unknown" for r in rows),
        }
        for role, rows in (("development", development), ("candidate", candidate))
    }
    reasons = []
    if overlaps["shared_source_studies"]:
        reasons.append("shared_primary_source_lineage")
    if overlaps["shared_source_observations"]:
        reasons.append("reused_source_observations")
    if overlaps["shared_exact_host_phage_pairs"]:
        reasons.append("exact_pair_overlap")
    if regime in {"cold_host", "double_cold"} and overlaps["shared_exact_host_assemblies"]:
        reasons.append("cold_host_identity_overlap")
    if regime in {"cold_phage", "double_cold"} and overlaps["shared_exact_phage_assemblies"]:
        reasons.append("cold_phage_identity_overlap")
    if any(missing[role][field] for role in ("development", "candidate")
           for field in ("missing_host_assembly_sha256", "missing_phage_assembly_sha256")):
        reasons.append("incomplete_assembly_identity_evidence")

    return {
        "schema_version": 1,
        "regime": regime,
        "status": "blocked" if reasons else "requires_manual_review",
        "independent_release_validation_established": False,
        "blocking_reasons": reasons,
        "row_counts": {"development": len(development), "candidate": len(candidate)},
        "source_endpoint_ids": {
            "development": sorted(_nonblank(development, "endpoint_id")),
            "candidate": sorted(_nonblank(candidate, "endpoint_id")),
        },
        "exact_overlap_counts": overlaps,
        "missingness": missing,
        "manual_review_still_required": [
            "source_study_aliases_and_upstream_dataset_derivations",
            "near_duplicate_genomes_and_homology_groups",
            "assay_endpoint_and_label_semantics",
            "training_feature_pretraining_and_model_selection_leakage",
            "license_and_provenance_of_original_files",
            "prospectively_locked_release_cohort_and_protocol",
        ],
        "interpretation": (
            "Known overlap/incomplete identity evidence blocks a release-independence claim."
            if reasons else
            "No disqualifying exact overlap was found; independence remains unverified."
        ),
    }
