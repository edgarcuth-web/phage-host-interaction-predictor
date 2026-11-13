"""Conservative accession/sequence identity gate for the published BASEL–ECOR panel.

A shared strain name is *not* proof of sequence identity. A different phage name
is *not* proof of sequence novelty. No sequence is fetched or aligned here.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re

_HEX64 = re.compile(r"^[0-9a-f]{64}$")


def audit_basel_identity(source_report: dict, development_resources: dict) -> dict:
    """Join verified source IDs to the locked development-resource roster.

    This is intentionally a *blocking* first-stage identity receipt. Only a
    separate accession/assembly review and sequence comparison can clear it.
    """
    if source_report.get("status") != "source_intake_verified_development_only":
        raise ValueError("BASEL/ECOR source intake must be byte-verified first")
    hosts, phages = source_report.get("host_ids"), source_report.get("phage_ids")
    counts = source_report.get("source_summary", {})
    if (not isinstance(hosts, list) or not hosts or not isinstance(phages, list)
            or not phages or len(hosts) != len(set(hosts))
            or len(phages) != len(set(phages))):
        raise ValueError("source host/phage ID roster is invalid")
    if counts.get("hosts") != len(hosts) or counts.get("phages") != len(phages):
        raise ValueError("source roster does not match byte-verified source summary")
    if source_report.get("independent_release_validation_established") is not False:
        raise ValueError("source report cannot establish independent release validation")
    resources = development_resources.get("resources")
    if not isinstance(resources, list):
        raise ValueError("development resource lock must have a resource list")

    training_hosts, training_phages = {}, {}
    for resource in resources:
        kind = resource.get("kind")
        if kind not in ("host_genome", "phage_genome"):
            continue
        name, digest = resource.get("id"), resource.get("sha256")
        if not isinstance(name, str) or not name or not isinstance(digest, str) or not _HEX64.fullmatch(digest):
            raise ValueError("development resource identity/hash missing or invalid")
        store = training_hosts if kind == "host_genome" else training_phages
        if name in store:
            raise ValueError(f"duplicate development {kind} ID: {name}")
        store[name] = resource

    host_ledger = []
    for host in sorted(hosts):
        match = training_hosts.get(host)
        host_ledger.append({
            "host_id": host,
            "found_in_development_resource_lock_by_exact_id": match is not None,
            "development_source_file_sha256": match["sha256"] if match else None,
            "development_source_url": match.get("url") if match else None,
            "development_source_file_size": match.get("size_bytes") if match else None,
            "genbank_assembly_accession": None,
            "normalized_sequence_sha256": None,
            "same_biological_strain_or_assembly_verified": False,
        })
    name_overlap = sorted(set(phages) & set(training_phages))
    candidate_ledger = [
        {
            "phage_id": phage,
            "name_in_development_resource_lock": phage in training_phages,
            "genbank_accession": None,
            "genbank_accession_version": None,
            "normalized_sequence_sha256": None,
            "sequence_novelty_verified": False,
        }
        for phage in sorted(phages)
    ]
    total_known = sum(entry["found_in_development_resource_lock_by_exact_id"] for entry in host_ledger)
    return {
        "schema_version": 1,
        "source_dataset_id": source_report["source_dataset_id"],
        "source_upstream_git_blob_sha1": source_report["upstream_blob_sha1"],
        "development_resource_lock_status": "resource_manifest_join_only_not_raw_sequence_reverification",
        "status": "blocked_pending_phage_accessions_sequence_comparison_and_host_assembly_review",
        "host_name_overlap_count": total_known,
        "host_name_unmatched_count": len(hosts) - total_known,
        "phage_name_overlap_count": len(name_overlap),
        "phage_name_overlap_ids": name_overlap,
        "phage_genome_accessions_resolved": 0,
        "phage_genome_accessions_missing": len(phages),
        "exact_phage_sequence_overlap_result": "not_evaluated",
        "phage_near_duplicate_result": "not_evaluated",
        "cold_host_claim_allowed": False,
        "cold_phage_claim_allowed": False,
        "independent_release_validation_established": False,
        "host_ledger": host_ledger,
        "phage_ledger": candidate_ledger,
        "limitations": [
            "Matching host IDs is name evidence, not a fresh equality check of assembled sequences.",
            "The development host file SHA-256 is a raw-byte identity, not a normalized sequence hash.",
            "Different phage names are not evidence of sequence or evolutionary independence.",
            "The 52 BASEL GenBank accession/version joins and genome sequences remain unresolved.",
            "A joint 94+52 phage similarity audit is required before cold-phage claims.",
            "Published reference data remain inspected development evidence, not an untouched release holdout.",
        ],
    }


def write_report_new(report: dict, out: str | Path) -> None:
    """Write exactly one immutable machine-readable intake receipt."""
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as file:
        json.dump(report, file, indent=2, sort_keys=True)
        file.write("\n")
