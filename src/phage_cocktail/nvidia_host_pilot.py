"""Bounded real-host engineering pilot with a durable no-resubmit journal."""
from datetime import datetime, timezone
from hashlib import sha256
import io
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .artifacts import atomic_write_json, code_provenance, read_json_object, sha256_file
from .cohort import canonical_assembly_digest, parse_dna_fasta, reverse_complement
from .evo2_host import host_windows
from .host_study import verify_resources
from .nvidia_evo2 import NvidiaForwardError, forward

PROTOCOL = "data/studies/nvidia_evo2_host_pilot.v1.json"


def prepare_jobs(repo, data_root, plan):
    for path, digest in plan["repo_files_sha256"].items():
        if sha256_file(repo / path) != digest:
            raise ValueError(f"changed frozen pilot input: {path}")
    genomes = read_json_object(repo / plan["genome_protocol"])
    manifest = pd.read_csv(repo / genomes["host_manifest"], dtype=str).sort_values("bacterium")
    selected = manifest.head(3)
    if selected.bacterium.tolist() != plan["hosts"]:
        raise ValueError("pilot hosts differ from frozen first three")
    entries = [r for r in genomes["resources"] if r["kind"] == "host_genome" and r["id"] in plan["hosts"]]
    resources = verify_resources({"resources": entries}, data_root)
    jobs = []
    for row in selected.itertuples(index=False):
        path = resources[("host_genome", row.bacterium)]
        if canonical_assembly_digest(path) != row.assembly_sha256:
            raise ValueError("pilot host assembly identity changed")
        payload = path.read_bytes()
        total_bases = sum(map(len, parse_dna_fasta(payload.decode("ascii"))))
        windows = host_windows(payload, window_bp=plan["window_bp"], max_windows_per_host=plan["sampling_max_windows"])
        i, contig, start, sequence = next((i, c, start, s) for i, (c, start, s) in enumerate(windows)
                                         if len(s) == plan["window_bp"])
        for strand, oriented in (("forward", sequence), ("reverse_complement", reverse_complement(sequence))):
            jobs.append({"host": row.bacterium, "assembly_sha256": row.assembly_sha256,
                         "fasta_sha256": sha256_file(path), "window_ordinal": i, "contig_ordinal": contig,
                         "start": start, "strand": strand, "length": len(oriented),
                         "assembly_bases": total_bases, "sampled_assembly_fraction": len(sequence) / total_bases,
                         "N_bases": oriented.count("N"),
                         "input_sha256": sha256(oriented.encode("ascii")).hexdigest(), "sequence": oriented})
    return jobs


def _save_vector(path, vector):
    stream = io.BytesIO()
    np.save(stream, vector, allow_pickle=False)
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(stream.getvalue())
    temporary.replace(path)


def run_pilot(repo, data_root, output, api_key, feasibility_cache, *, resume=False, call_forward=forward):
    repo, output, feasibility_cache = Path(repo), Path(output), Path(feasibility_cache)
    if output.exists() and not resume:
        raise ValueError("pilot output must be new, or use explicit resume")
    if resume and not output.exists():
        raise ValueError("cannot resume a missing pilot")
    plan = read_json_object(repo / PROTOCOL)
    jobs = prepare_jobs(repo, data_root, plan)
    protocol_sha = sha256_file(repo / PROTOCOL)
    provenance = code_provenance(repo, repo / "scripts/run_nvidia_host_pilot.py")
    parent = read_json_object(repo / plan["feasibility_report"])
    if parent["status"] != "success" or parent != read_json_object(feasibility_cache / "report.json"):
        raise ValueError("pilot requires the exact successful feasibility run")
    source_vector = feasibility_cache / "pooled_mean.npy"
    if sha256_file(source_vector) != parent["output_sha256"]["pooled_mean.npy"]:
        raise ValueError("reused feasibility vector hash mismatch")
    if jobs[0]["input_sha256"] != parent["provenance"]["input"]["input_sha256"]:
        raise ValueError("reused feasibility window differs")
    output.mkdir(parents=True, exist_ok=True)
    journal_path = output / "journal.json"
    if resume:
        journal = read_json_object(journal_path)
        if journal["protocol_sha256"] != protocol_sha or journal["code_sha256"] != provenance["code_sha256"]:
            raise ValueError("resume protocol or source code changed")
        if (output / "report.json").exists():
            report = read_json_object(output / "report.json")
            if report["protocol_sha256"] != protocol_sha or report["status"] != "complete_engineering_pilot":
                raise ValueError("completed pilot report identity changed")
            for name, digest in report["output_sha256"].items():
                if sha256_file(output / name) != digest:
                    raise ValueError("completed pilot output hash mismatch")
            return report
    else:
        journal = {"protocol_sha256": protocol_sha, "code_sha256": provenance["code_sha256"],
                   "new_post_attempts": 0, "new_response_bytes": 0, "jobs": {}}
        atomic_write_json(plan, output / "protocol.json")
        atomic_write_json(journal, journal_path)
    records, vectors = [], []
    budget = plan["request_budget"]
    for index, job in enumerate(jobs):
        name = f"job_{index:02d}"
        metadata_path, vector_path = output / (name + ".json"), output / (name + ".npy")
        public = {k: v for k, v in job.items() if k != "sequence"}
        if name in journal["jobs"]:
            if journal["jobs"][name] != "complete":
                raise ValueError("unresolved attempted job blocks resume; no automatic resubmission")
            record = read_json_object(metadata_path)
            if record["input"] != public or sha256_file(vector_path) != record["pooled_sha256"]:
                raise ValueError("completed pilot cache input or vector changed")
            vector = np.load(vector_path, allow_pickle=False)
        else:
            if index == 0:
                vector = np.load(source_vector, allow_pickle=False)
                record = {"input": public, "origin": "reused_feasibility_forward", "response": parent["response"]}
            else:
                remaining = budget["max_new_response_bytes"] - journal["new_response_bytes"]
                if journal["new_post_attempts"] >= budget["max_new_inference_posts"] or remaining <= 0:
                    raise ValueError("frozen pilot request or transfer budget exhausted")
                journal["new_post_attempts"] += 1
                journal["jobs"][name] = "attempted_unresolved"
                atomic_write_json(journal, journal_path)  # Before entering the network call.
                started = time.monotonic()
                try:
                    result = call_forward(job["sequence"], api_key, timeout=budget["io_timeout_seconds"],
                                          poll_seconds=budget["server_wait_seconds"],
                                          max_response_bytes=min(budget["max_response_bytes"], remaining),
                                          max_tensor_bytes=budget["max_tensor_bytes"])
                except NvidiaForwardError as exc:
                    atomic_write_json({"input": public, "error": str(exc), "http_status": exc.http_status,
                                       "request_id": exc.request_id, "response_bytes_read": exc.response_bytes}, metadata_path)
                    raise
                mean = result.tokens.mean(axis=0, dtype=np.float64)
                vector = mean.astype(np.float32)
                record = {"input": public, "origin": "new_hosted_forward", "response": result.metadata,
                          "client_elapsed_seconds": time.monotonic() - started,
                          "pool_float32_max_abs_error": float(np.max(np.abs(mean - vector)))}
                journal["new_response_bytes"] += result.metadata["response_bytes"]
                del result
            if vector.shape != (plan["dimension"],) or vector.dtype != np.float32 or not np.isfinite(vector).all() or not np.any(vector):
                raise ValueError("invalid pilot pooled vector")
            _save_vector(vector_path, vector)
            record["pooled_sha256"] = sha256_file(vector_path)
            atomic_write_json(record, metadata_path)
            journal["jobs"][name] = "complete"
            atomic_write_json(journal, journal_path)
            print(f"Completed {job['host']} {job['strand']}: {record['origin']}", flush=True)
        if vector.shape != (plan["dimension"],) or vector.dtype != np.float32 or not np.isfinite(vector).all() or not np.any(vector):
            raise ValueError("invalid cached pilot vector")
        records.append(record)
        vectors.append(vector)
    host_vectors = np.stack([np.mean(vectors[i:i+2], axis=0, dtype=np.float64).astype(np.float32)
                             for i in range(0, len(vectors), 2)])
    _save_vector(output / "host_vectors.npy", host_vectors)
    report = {"status": "complete_engineering_pilot", "observed_at_utc": datetime.now(timezone.utc).isoformat(),
              "protocol_sha256": protocol_sha, "provenance": provenance, "hosts": plan["hosts"],
              "host_vector_shape": list(host_vectors.shape), "new_post_attempts": journal["new_post_attempts"],
              "new_response_bytes": journal["new_response_bytes"], "reused_jobs": 1, "jobs": records,
              "identical_host_vector_pairs": [[plan["hosts"][i], plan["hosts"][j]]
                                               for i in range(len(host_vectors)) for j in range(i)
                                               if np.array_equal(host_vectors[i], host_vectors[j])],
              "local_v1_cache_compatible": False, "prediction_metrics_computed": False,
              "output_sha256": {p.name: sha256_file(p) for p in sorted(output.iterdir()) if p.name != "report.json"}}
    atomic_write_json(report, output / "report.json")
    return report
