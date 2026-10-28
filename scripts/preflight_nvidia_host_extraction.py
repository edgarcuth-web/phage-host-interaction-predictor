"""Offline budget and integrity preflight for a full hosted Evo2 host extraction.

Reads a previously frozen extraction plan. Never starts or authorizes requests.
"""
import argparse
from collections import defaultdict
from hashlib import sha256
from pathlib import Path

import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.nvidia_host_extraction import (
    GENOMES, PILOT, validate_plan, verify_sources,
)

REFERENCE = "data/studies/results/nvidia_evo2_window_feasibility.v1.json"
PROJECTION = "data/studies/results/host_window_sampling.v1/budget_projection.csv"
SELECTION = "original_32_including_short_tails"


def inspect_plan(repo, plan, *, source_verified=False):
    """Compute a reviewable estimate; not an approval to make provider calls."""
    repo = Path(repo)
    validate_plan(plan)
    if plan["selection"] != SELECTION:
        raise ValueError("full hosted study requires original_32_including_short_tails")

    source = read_json_object(repo / GENOMES)
    manifest = pd.read_csv(repo / source["host_manifest"], dtype=str)
    expected_hosts = sorted(manifest.bacterium.tolist())
    if len(expected_hosts) != len(set(expected_hosts)) or plan["hosts"] != expected_hosts:
        raise ValueError("plan does not cover the complete frozen host panel")

    reference_path = repo / REFERENCE
    ref = read_json_object(reference_path)["response"]
    if (ref["http_status"] != 200 or ref["sequence_length"] != 8192
            or ref["layer"] != plan["layer"] or ref["endpoint"] != plan["endpoint"]
            or type(ref["response_bytes"]) is not int or ref["response_bytes"] < 1):
        raise ValueError("real transfer-reference metadata is incompatible with the plan")

    per_host = defaultdict(lambda: {"windows": 0, "selected_bases": 0, "assembly_bases": set()})
    forward_jobs = [job for job in plan["jobs"] if job["strand"] == "forward"]
    for job in forward_jobs:
        row = per_host[job["host"]]
        row["windows"] += 1
        row["selected_bases"] += job["length"]
        row["assembly_bases"].add(job["assembly_bases"])

    host_rows = []
    for host in expected_hosts:
        row = per_host[host]
        if row["windows"] < 1 or row["windows"] > 32 or len(row["assembly_bases"]) != 1:
            raise ValueError("invalid hosted selection or assembly-bases identity")
        assembly_bases = next(iter(row["assembly_bases"]))
        if assembly_bases <= 0 or row["selected_bases"] > assembly_bases:
            raise ValueError("invalid selected-base coverage")
        host_rows.append({
            "bacterium": host,
            "windows": row["windows"],
            "selected_bases": row["selected_bases"],
            "assembly_bases": assembly_bases,
            "sampled_fraction": row["selected_bases"] / assembly_bases,
        })

    count = len(plan["jobs"])
    selected_bases_two_strands = sum(job["length"] for job in plan["jobs"])
    projected_json_bytes = round(
        selected_bases_two_strands * ref["response_bytes"] / ref["sequence_length"]
    )
    cap = plan["budget"]["max_response_bytes"]
    reserved_no_reuse = count * (cap + 1)

    # Compare with the existing independent, FASTA-verified sampling-control
    # projection. A difference is a review blocker, never a silent amendment.
    recorded = pd.read_csv(repo / PROJECTION)
    matches = recorded.loc[recorded.strategy.eq(SELECTION)]
    if len(matches) != 1:
        raise ValueError("missing uniquely recorded original-window projection")
    previous = matches.iloc[0]
    if int(previous.hosts) != len(expected_hosts) or int(previous.two_strand_calls) != count:
        raise ValueError("full-plan job count differs from the recorded verified cohort")
    if int(previous.estimated_json_bytes_no_reuse_credit) != projected_json_bytes:
        raise ValueError("full-plan transfer projection differs from the recorded verified cohort")

    pilot = read_json_object(repo / PILOT)
    reference_inputs = {
        tuple(sorted(record["input"].items())) for record in pilot["jobs"]
    }
    candidate_imports = sum(
        tuple(sorted({k: v for k, v in job.items() if k != "id"}.items())) in reference_inputs
        for job in plan["jobs"]
    )
    budget = plan["budget"]
    return {
        "status": "offline_budget_preflight_only",
        "study": "full-host-nvidia-evo2-40b:v1",
        "new_nvidia_requests": 0,
        "source_fasta_verification": "completed" if source_verified else "not_performed",
        "plan_ready_to_run": False,
        "selection": SELECTION,
        "hosts": len(expected_hosts),
        "selected_windows": len(forward_jobs),
        "two_strand_jobs": count,
        "selected_bases_two_strands": selected_bases_two_strands,
        "median_assembly_fraction": float(pd.DataFrame(host_rows).sampled_fraction.median()),
        "min_assembly_fraction": min(r["sampled_fraction"] for r in host_rows),
        "max_assembly_fraction": max(r["sampled_fraction"] for r in host_rows),
        "reference_response_bytes": ref["response_bytes"],
        "reference_response_bases": ref["sequence_length"],
        "reference_sha256": sha256_file(reference_path),
        "recorded_projection_sha256": sha256_file(repo / PROJECTION),
        "estimated_json_bytes_no_reuse_credit": projected_json_bytes,
        "max_response_cap_bytes": cap,
        "conservative_reserved_bytes_no_reuse": reserved_no_reuse,
        "frozen_max_new_calls": budget["max_new_calls"],
        "frozen_max_new_response_bytes": budget["max_new_response_bytes"],
        "full_cap_budget_sufficient_no_reuse": (
            budget["max_new_calls"] >= count
            and budget["max_new_response_bytes"] >= reserved_no_reuse
        ),
        "potential_pilot_job_matches_unverified": candidate_imports,
        "per_host": host_rows,
        "limitations": [
            "Projection is linear in input bases from one measured 8192-base hosted response; not a bill, quota guarantee or measured whole-cohort run.",
            "Max-response reservation assumes every unimported job can reserve the full configured response cap plus one sentinel byte; actual transfer may be lower.",
            "Potential pilot imports are identified by manifest metadata only. No credit is taken without the verified pilot cache.",
            "Provider model weights and server runtime digests are unavailable and cannot be verified as a local checkpoint.",
            "This preflight never authorizes, starts or resumes a hosted extraction; a human must separately review budget and deployment.",
            "Exact FASTA bytes and executor code are rechecked by run_extraction before any requests, regardless of this optional offline verification.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--verify-sources", action="store_true",
                        help="Also hash/verify the exact FASTA files and executor code locally")
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    args = parser.parse_args()
    if args.out.exists():
        parser.error("output already exists; preflight records are immutable")

    repo = Path(__file__).resolve().parents[1]
    plan = read_json_object(args.plan)
    if args.verify_sources:
        verify_sources(repo, args.data_root, plan)
    result = inspect_plan(repo, plan, source_verified=args.verify_sources)
    result["plan_file_sha256"] = sha256_file(args.plan)
    result["plan_content_sha256"] = sha256(
        args.plan.read_bytes()
    ).hexdigest()
    atomic_write_json(result, args.out)
    print({k: result[k] for k in (
        "status", "hosts", "two_strand_jobs", "estimated_json_bytes_no_reuse_credit",
        "conservative_reserved_bytes_no_reuse", "full_cap_budget_sufficient_no_reuse",
        "source_fasta_verification", "new_nvidia_requests",
    )})


if __name__ == "__main__":
    main()
