"""Verified adapter from hosted Evo2 extraction checkpoints to retrieval features.

Hosted NVIDIA forward extraction has deliberately different numerical pooling
semantics from the local-checkpoint Evo2 v1 cache.  This module preserves that
distinction: it verifies a completed hosted extraction and exposes its host
vectors as a separately identified representation.  It never rewrites the
hosted output into the local-v1 cache format.
"""
from hashlib import sha256
import io
from pathlib import Path

import numpy as np
import pandas as pd

from .artifacts import read_json_object
from .host_genome import host_genome_features
from .host_study import verify_resources
from .nvidia_host_extraction import (
    GENOMES,
    json_bytes,
    read_object,
    validate_plan,
    verify_extraction,
    verify_sources,
)


def load_hosted_cache(repo, data_root, plan_path, store, manifest):
    """Load one complete full-panel hosted extraction for retrieval evaluation.

    The extraction plan remains the representation protocol.  The retrieval
    runner may reuse the same frozen host splits, candidate panel and metrics,
    but must label this representation separately from the local-checkpoint v1
    Evo2 cache.
    """

    repo, plan_path = Path(repo), Path(plan_path)
    plan = read_json_object(plan_path)
    validate_plan(plan)
    if plan["selection"] != "original_32_including_short_tails":
        raise ValueError(
            "hosted retrieval requires the frozen original 32-window selection"
        )

    hosts = sorted(manifest.bacterium.astype(str).tolist())
    if len(hosts) != len(set(hosts)) or plan["hosts"] != hosts:
        raise ValueError("hosted extraction must cover the exact full host manifest")

    # Bind the plan to the exact source FASTAs and extraction code before using
    # any persisted vectors.  This performs no NVIDIA requests.
    verify_sources(repo, data_root, plan)

    status = verify_extraction(plan, store)
    if status["status"] != "complete_engineering_extraction":
        raise ValueError("hosted extraction is not complete and verified")

    plan_payload = store.get("plan.json")
    report_payload = store.get("report.json")
    matrix_payload = store.get("host_vectors.npy")
    if plan_payload != json_bytes(plan):
        raise ValueError("stored hosted extraction plan differs from the supplied plan")
    if report_payload is None or matrix_payload is None:
        raise ValueError("complete hosted extraction is missing final artifacts")

    report = read_object(report_payload)
    expected_shape = [len(hosts), plan["dimension"]]
    if (
        report.get("hosts") != hosts
        or report.get("host_vector_shape") != expected_shape
        or report.get("local_v1_cache_compatible") is not False
        or report.get("prediction_metrics_computed") is not False
        or report.get("host_vectors_sha256") != sha256(matrix_payload).hexdigest()
    ):
        raise ValueError("hosted extraction report differs from evaluation inputs")

    values = np.load(io.BytesIO(matrix_payload), allow_pickle=False)
    if (
        values.dtype != np.float32
        or values.shape != tuple(expected_shape)
        or not np.isfinite(values).all()
    ):
        raise ValueError("invalid hosted host-vector matrix")

    source = read_json_object(repo / GENOMES)
    resources = verify_resources(source, data_root)
    sequence_features, genome_audit = host_genome_features(manifest, resources)
    if sequence_features.index.tolist() != hosts:
        raise ValueError("host genome feature order differs from hosted extraction")

    # Count each sampled genomic window once; the adjacent reverse-complement
    # job is a second orientation of the same sequence, not additional coverage.
    windows = [job for job in plan["jobs"] if job["strand"] == "forward"]
    rows = []
    for host in hosts:
        selected = [job for job in windows if job["host"] == host]
        if not selected:
            raise ValueError("hosted extraction host has no sampled windows")
        assembly_bases = {job["assembly_bases"] for job in selected}
        assembly_sha = {job["assembly_sha256"] for job in selected}
        fasta_sha = {job["fasta_sha256"] for job in selected}
        if len(assembly_bases) != 1 or len(assembly_sha) != 1 or len(fasta_sha) != 1:
            raise ValueError("hosted extraction host identity is inconsistent across windows")
        selected_bases = sum(job["length"] for job in selected)
        total_bases = next(iter(assembly_bases))
        rows.append(
            {
                "bacterium": host,
                "extraction_assembly_sha256": next(iter(assembly_sha)),
                "extraction_fasta_sha256": next(iter(fasta_sha)),
                "selected_windows": len(selected),
                "selected_bases": selected_bases,
                "sampled_fraction": selected_bases / total_bases,
            }
        )
    qc = genome_audit.merge(pd.DataFrame(rows), on="bacterium", validate="one_to_one")
    if not qc.sampled_fraction.between(0, 1).all():
        raise ValueError("invalid hosted sampled-genome coverage")

    features = pd.DataFrame(
        values,
        index=hosts,
        columns=[f"hosted_evo2:{i}" for i in range(plan["dimension"])],
    )
    plan_sha = sha256(plan_payload).hexdigest()
    source_identity = sha256(
        json_bytes(
            {
                "kind": "nvidia_hosted_evo2_evaluation_source",
                "plan_sha256": plan_sha,
                "report_sha256": sha256(report_payload).hexdigest(),
                "host_vectors_sha256": sha256(matrix_payload).hexdigest(),
            }
        )
    ).hexdigest()
    representation = {
        "source": "nvidia_hosted_forward",
        "model_name": "evo2_40b",
        "endpoint": plan["endpoint"],
        "layer": plan["layer"],
        "dimension": plan["dimension"],
        "selection": plan["selection"],
        "pooling": plan["pooling"],
        "weights_sha256": plan["weights_sha256"],
        "server_runtime_version": plan["server_runtime_version"],
        "local_v1_cache_compatible": False,
        "source_plan_sha256": plan_sha,
    }
    return features, sequence_features, qc, source_identity, representation
