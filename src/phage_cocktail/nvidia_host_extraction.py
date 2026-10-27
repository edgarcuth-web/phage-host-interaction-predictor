"""Ordered immutable checkpoints for bounded hosted Evo2 extraction.

No live requests by default. Completed receipts contain their small float32
vector, so a vector and its completion metadata are published atomically.
"""
import base64
from hashlib import sha256
import io
import json
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd

from .artifacts import canonical_json, code_provenance, read_json_object, sha256_file
from .cohort import canonical_assembly_digest, parse_dna_fasta, reverse_complement
from .evo2_host import host_windows
from .extraction_store import put_identical
from .host_study import verify_resources
from .nvidia_evo2 import DIMENSION, FORWARD_URL, LAYER, forward

GENOMES = "data/studies/genome_host_retrieval.v1.json"
PILOT = "data/studies/results/nvidia_evo2_host_pilot.v1/report.json"
SCRIPT = "scripts/run_nvidia_host_extraction.py"
POOLING = "float64 token mean -> float32 orientation; float64 strand mean -> float32 window; length-weighted float64 host mean -> float32"


def json_bytes(value):
    return (canonical_json(value) + "\n").encode("utf-8")


def read_object(payload):
    def reject(value):
        raise ValueError("nonfinite checkpoint JSON")
    value = json.loads(payload, parse_constant=reject)
    if not isinstance(value, dict):
        raise ValueError("checkpoint JSON must be an object")
    return value


def _array_bytes(values):
    stream = io.BytesIO()
    np.save(stream, values, allow_pickle=False)
    return stream.getvalue()


def _valid_vector(vector):
    if (vector.shape != (DIMENSION,) or vector.dtype != np.float32
            or not np.isfinite(vector).all() or not np.any(vector)):
        raise ValueError("invalid pooled orientation vector")


def prepare_inputs(repo, data_root, hosts, selection):
    repo, data_root = Path(repo), Path(data_root)
    source = read_json_object(repo / GENOMES)
    for name, digest in source.get("repo_files_sha256", {}).items():
        if sha256_file(repo / name) != digest:
            raise ValueError("frozen parent study input changed")
    manifest = pd.read_csv(repo / source["host_manifest"], dtype=str).sort_values("bacterium")
    if not hosts or hosts != sorted(set(hosts)) or not set(hosts) <= set(manifest.bacterium):
        raise ValueError("hosts must be a non-empty sorted unique subset of the frozen panel")
    if selection not in ("first_full_window", "original_32_including_short_tails"):
        raise ValueError("unknown fixed window selection")
    entries = [r for r in source["resources"] if r["kind"] == "host_genome" and r["id"] in hosts]
    resources = verify_resources({"resources": entries}, data_root)
    jobs, sequences = [], {}
    for row in manifest[manifest.bacterium.isin(hosts)].itertuples(index=False):
        path = resources[("host_genome", row.bacterium)]
        if canonical_assembly_digest(path) != row.assembly_sha256:
            raise ValueError("assembly differs from frozen panel")
        payload = path.read_bytes()
        fasta_sha = sha256_file(path)
        bases = sum(map(len, parse_dna_fasta(payload.decode("ascii"))))
        windows = list(enumerate(host_windows(payload, window_bp=8192, max_windows_per_host=32)))
        if selection == "first_full_window":
            windows = next(([item] for item in windows if len(item[1][2]) == 8192), [])
            if not windows:
                raise ValueError("host has no full sampled window")
        for ordinal, (contig, start, sequence) in windows:
            for strand, oriented in (("forward", sequence), ("reverse_complement", reverse_complement(sequence))):
                name = f"{len(jobs):05d}"
                jobs.append({"id": name, "host": row.bacterium, "assembly_sha256": row.assembly_sha256,
                             "fasta_sha256": fasta_sha, "window_ordinal": ordinal,
                             "contig_ordinal": contig, "start": start, "strand": strand,
                             "length": len(oriented), "assembly_bases": bases,
                             "sampled_assembly_fraction": len(oriented) / bases,
                             "N_bases": oriented.count("N"),
                             "input_sha256": sha256(oriented.encode("ascii")).hexdigest()})
                sequences[name] = oriented
    return jobs, sequences


def make_plan(repo, data_root, *, hosts=None, selection="first_full_window",
              max_new_calls=0, max_new_response_bytes=0):
    repo = Path(repo)
    for value in (max_new_calls, max_new_response_bytes):
        if type(value) is not int or value < 0:
            raise ValueError("request budgets must be non-negative integers")
    source = read_json_object(repo / GENOMES)
    panel = sorted(pd.read_csv(repo / source["host_manifest"], dtype=str).bacterium)
    hosts = panel[:3] if hosts is None else hosts
    jobs, _ = prepare_inputs(repo, data_root, hosts, selection)
    return {"schema_version": 1, "scope": "hosted engineering extraction; no prediction metrics",
            "endpoint": FORWARD_URL, "layer": LAYER, "dimension": DIMENSION,
            "selection": selection, "hosts": hosts, "jobs": jobs, "pooling": POOLING,
            "repo_files_sha256": {name: sha256_file(repo / name) for name in
                                  (GENOMES, source["host_manifest"], PILOT)},
            "provenance": code_provenance(repo, repo / SCRIPT),
            "budget": {"max_new_calls": max_new_calls, "max_new_response_bytes": max_new_response_bytes,
                       "max_response_bytes": 256 * 1024 * 1024, "max_tensor_bytes": 640 * 1024 * 1024,
                       "io_timeout_seconds": 180, "server_wait_seconds": 60},
            "weights_sha256": None, "server_runtime_version": None,
            "local_v1_cache_compatible": False}


def validate_plan(plan):
    if (plan["schema_version"] != 1 or plan["endpoint"] != FORWARD_URL or plan["layer"] != LAYER
            or plan["dimension"] != DIMENSION or plan["pooling"] != POOLING
            or plan["selection"] not in ("first_full_window", "original_32_including_short_tails")
            or plan["weights_sha256"] is not None or plan["server_runtime_version"] is not None
            or plan["local_v1_cache_compatible"] is not False):
        raise ValueError("unsupported hosted extraction protocol")
    budget = plan["budget"]
    for key, upper in (("max_new_calls", None), ("max_new_response_bytes", None),
                       ("max_response_bytes", 256 * 1024 * 1024), ("max_tensor_bytes", 640 * 1024 * 1024),
                       ("io_timeout_seconds", 180), ("server_wait_seconds", 60)):
        v = budget[key]
        minimum = 0 if key.startswith("max_new_") else 1
        if type(v) is not int or v < minimum or (upper is not None and v > upper):
            raise ValueError("invalid frozen request budget")
    jobs = plan["jobs"]
    if not jobs or len(jobs) % 2 or plan["hosts"] != sorted(set(plan["hosts"])):
        raise ValueError("invalid job/host support")
    if sorted({j["host"] for j in jobs}) != plan["hosts"]:
        raise ValueError("job/host support differs")
    for i, job in enumerate(jobs):
        if job["id"] != f"{i:05d}" or type(job["length"]) is not int or not 1 <= job["length"] <= 8192:
            raise ValueError("invalid ordered job identity")
        if i % 2 == 0:
            other = jobs[i+1]
            if (job["strand"] != "forward" or other["strand"] != "reverse_complement"
                    or {k:v for k,v in job.items() if k not in ("id", "strand", "input_sha256")}
                    != {k:v for k,v in other.items() if k not in ("id", "strand", "input_sha256")}):
                raise ValueError("window strand pairing differs")


def verify_sources(repo, data_root, plan):
    validate_plan(plan)
    repo = Path(repo)
    for name, digest in {**plan["repo_files_sha256"], **plan["provenance"]["code_sha256"]}.items():
        if not (repo / name).resolve().is_relative_to(repo.resolve()) or sha256_file(repo / name) != digest:
            raise ValueError("frozen extraction source or code changed")
    current = code_provenance(repo, repo / SCRIPT)
    if current["code_sha256"] != plan["provenance"]["code_sha256"]:
        raise ValueError("frozen extraction package changed")
    jobs, sequences = prepare_inputs(repo, data_root, plan["hosts"], plan["selection"])
    if jobs != plan["jobs"]:
        raise ValueError("prepared inputs differ from frozen extraction jobs")
    return sequences


def pilot_receipts(repo, cache, plan):
    if cache is None:
        return {}
    repo, cache = Path(repo), Path(cache)
    report = read_json_object(cache / "report.json")
    if report != read_json_object(repo / PILOT) or report["status"] != "complete_engineering_pilot":
        raise ValueError("requires the exact committed real pilot")
    for name, digest in report["output_sha256"].items():
        if Path(name).name != name or sha256_file(cache / name) != digest:
            raise ValueError("pilot payload hash mismatch")
    imported = {}
    for i, record in enumerate(report["jobs"]):
        payload = (cache / f"job_{i:02d}.npy").read_bytes()
        if sha256(payload).hexdigest() != record["pooled_sha256"]:
            raise ValueError("pilot vector hash mismatch")
        vector = np.load(io.BytesIO(payload), allow_pickle=False)
        _valid_vector(vector)
        for job in plan["jobs"]:
            if {k:v for k,v in job.items() if k != "id"} == record["input"]:
                imported[job["id"]] = {"vector": vector, "response": record["response"],
                                       "pilot_report_sha256": sha256_file(cache / "report.json")}
    return imported


def _keys(job):
    prefix = "jobs/" + job["id"]
    return prefix + "/claim.json", prefix + "/result.json"


def _receipt(job, claim_payload, vector, response):
    _valid_vector(vector)
    payload = _array_bytes(vector)
    return json_bytes({"input": job, "claim_sha256": sha256(claim_payload).hexdigest(),
                       "vector_sha256": sha256(payload).hexdigest(),
                       "vector_npy_base64": base64.b64encode(payload).decode("ascii"),
                       "response": response})


def _response_valid(response, job):
    if (response.get("endpoint") != FORWARD_URL or response.get("layer") != LAYER
            or response.get("http_status") != 200 or response.get("finite") is not True
            or response.get("shape") != [job["length"], DIMENSION]
            or response.get("input_sha256") != job["input_sha256"]
            or response.get("sequence_length") != job["length"]
            or type(response.get("response_bytes")) is not int or response["response_bytes"] < 1
            or response.get("weights_sha256") is not None):
        raise ValueError("response metadata differs from extraction input")


def load_state(plan, store):
    """Read every known key, including later jobs, before any new inference."""
    validate_plan(plan)
    plan_payload = json_bytes(plan)
    if store.get("plan.json") != plan_payload:
        raise ValueError("stored extraction plan differs")
    plan_sha = sha256(plan_payload).hexdigest()
    states = []
    calls = reserved = actual = 0
    for job in plan["jobs"]:
        claim_key, result_key = _keys(job)
        claim_payload, payload = store.get(claim_key), store.get(result_key)
        claim = read_object(claim_payload) if claim_payload is not None else None
        vector = None
        if claim is not None:
            if claim["plan_sha256"] != plan_sha or claim["input"] != job:
                raise ValueError("checkpoint claim identity differs")
            if claim["origin"] == "live":
                cap = claim["response_cap"]
                if type(cap) is not int or not 1 <= cap <= plan["budget"]["max_response_bytes"]:
                    raise ValueError("checkpoint response reservation differs")
                calls += 1
                reserved += cap + 1
            elif claim["origin"] != "pilot":
                raise ValueError("unknown checkpoint origin")
        if payload is not None:
            if claim is None:
                raise ValueError("result has no preceding claim")
            record = read_object(payload)
            if record["input"] != job or record["claim_sha256"] != sha256(claim_payload).hexdigest():
                raise ValueError("checkpoint receipt identity differs")
            array_payload = base64.b64decode(record["vector_npy_base64"], validate=True)
            if sha256(array_payload).hexdigest() != record["vector_sha256"]:
                raise ValueError("checkpoint vector hash mismatch")
            vector = np.load(io.BytesIO(array_payload), allow_pickle=False)
            _valid_vector(vector)
            _response_valid(record["response"], job)
            if claim["origin"] == "live":
                if record["response"]["response_bytes"] > claim["response_cap"]:
                    raise ValueError("response exceeds its reservation")
                actual += record["response"]["response_bytes"]
        states.append({"claim": claim, "claim_payload": claim_payload, "payload": payload, "vector": vector})
    if calls > plan["budget"]["max_new_calls"] or reserved > plan["budget"]["max_new_response_bytes"]:
        raise ValueError("persisted extraction exceeds its frozen budget")
    return states, {"new_post_attempts": calls, "reserved_response_bytes": reserved,
                    "successful_new_response_bytes": actual}


def _host_matrix(plan, states):
    totals = {h: 0 for h in plan["hosts"]}
    sums = {h: np.zeros(DIMENSION, dtype=np.float64) for h in plan["hosts"]}
    for i in range(0, len(states), 2):
        job = plan["jobs"][i]
        window = np.mean([states[i]["vector"], states[i+1]["vector"]], axis=0, dtype=np.float64).astype(np.float32)
        sums[job["host"]] += window.astype(np.float64) * job["length"]
        totals[job["host"]] += job["length"]
    return np.stack([(sums[h] / totals[h]).astype(np.float32) for h in plan["hosts"]])


def verify_extraction(plan, store):
    states, counts = load_state(plan, store)
    complete = all(s["vector"] is not None for s in states)
    status = {"status": "ready_to_finalize" if complete else "incomplete", **counts,
              "completed_jobs": sum(s["vector"] is not None for s in states),
              "planned_jobs": len(states),
              "unresolved_live_jobs": [j["id"] for j,s in zip(plan["jobs"], states)
                                       if s["claim"] is not None and s["claim"]["origin"] == "live" and s["vector"] is None]}
    stored_report = store.get("report.json")
    if stored_report is not None:
        if not complete:
            raise ValueError("completed report has incomplete receipts")
        expected, matrix = _final_payloads(plan, states, counts)
        if store.get("host_vectors.npy") != matrix or stored_report != json_bytes(expected):
            raise ValueError("final extraction export differs from verified receipts")
        status["status"] = "complete_engineering_extraction"
    return status


def _final_payloads(plan, states, counts):
    matrix = _array_bytes(_host_matrix(plan, states))
    report = {"status": "complete_engineering_extraction", "hosts": plan["hosts"],
              "host_vector_shape": [len(plan["hosts"]), DIMENSION], **counts,
              "plan_sha256": sha256(json_bytes(plan)).hexdigest(),
              "receipt_sha256": {j["id"]: sha256(s["payload"]).hexdigest() for j,s in zip(plan["jobs"], states)},
              "host_vectors_sha256": sha256(matrix).hexdigest(),
              "prediction_metrics_computed": False, "local_v1_cache_compatible": False,
              "weights_sha256": None, "server_runtime_version": None}
    return report, matrix


def run_extraction(repo, data_root, plan, store, *, pilot_cache=None, allow_live=False,
                   key_provider=None, max_jobs=None, call_forward=forward):
    if max_jobs is not None and (type(max_jobs) is not int or max_jobs < 1):
        raise ValueError("max_jobs must be positive")
    sequences = verify_sources(repo, data_root, plan)
    imported = pilot_receipts(repo, pilot_cache, plan)
    if hasattr(store, "check_conditional_writes"):
        store.check_conditional_writes()
    put_identical(store, "plan.json", json_bytes(plan))
    states, counts = load_state(plan, store)
    # Check the whole run first. Never skip an ambiguous POST and advance.
    for job, state in zip(plan["jobs"], states):
        claim = state["claim"]
        if claim is not None and state["vector"] is None:
            if claim["origin"] == "live":
                raise ValueError("unresolved live attempt blocks resume; never automatically resubmit")
            if job["id"] not in imported:
                raise ValueError("unfinished pilot import requires the same verified pilot cache")
    done, stop_reason = 0, None
    plan_sha = sha256(json_bytes(plan)).hexdigest()
    for job, state in zip(plan["jobs"], states):
        if state["vector"] is not None:
            continue
        if max_jobs is not None and done >= max_jobs:
            stop_reason = "invocation_job_limit"
            break
        claim_key, result_key = _keys(job)
        source = imported.get(job["id"])
        if source is not None:
            claim = {"plan_sha256": plan_sha, "input": job, "origin": "pilot",
                     "pilot_report_sha256": source["pilot_report_sha256"]}
            claim_payload = json_bytes(claim)
            put_identical(store, claim_key, claim_payload)
            vector, response = source["vector"], source["response"]
        else:
            if not allow_live:
                stop_reason = "live_requests_disabled"
                break
            budget = plan["budget"]
            remaining = budget["max_new_response_bytes"] - counts["reserved_response_bytes"]
            if counts["new_post_attempts"] >= budget["max_new_calls"] or remaining <= 1:
                stop_reason = "frozen_budget_exhausted"
                break
            key = key_provider() if key_provider is not None else None
            if not isinstance(key, str) or not key or any(c.isspace() for c in key):
                raise ValueError("live execution requires a non-empty NVIDIA API key")
            cap = min(budget["max_response_bytes"], remaining - 1)
            claim = {"plan_sha256": plan_sha, "input": job, "origin": "live",
                     "response_cap": cap, "attempt_id": str(uuid4())}
            claim_payload = json_bytes(claim)
            if not store.create(claim_key, claim_payload):
                raise ValueError("another worker claimed this job; stop without submitting")
            # Reservation includes the client's one-byte oversize sentinel read.
            result = call_forward(sequences[job["id"]], key, timeout=budget["io_timeout_seconds"],
                                  poll_seconds=budget["server_wait_seconds"], max_response_bytes=cap,
                                  max_tensor_bytes=budget["max_tensor_bytes"])
            key = None
            response = result.metadata
            _response_valid(response, job)
            if result.tokens.shape != (job["length"], DIMENSION) or not np.isfinite(result.tokens).all():
                raise ValueError("returned token tensor differs")
            if response["response_bytes"] > cap:
                raise ValueError("returned response exceeds reserved download")
            vector = result.tokens.mean(axis=0, dtype=np.float64).astype(np.float32)
            del result
            counts["new_post_attempts"] += 1
            counts["reserved_response_bytes"] += cap + 1
        _response_valid(response, job)
        receipt = _receipt(job, claim_payload, vector, response)
        put_identical(store, result_key, receipt)
        done += 1
    # Re-read persisted state: publication, not an in-memory counter, defines completion.
    states, counts = load_state(plan, store)
    if all(s["vector"] is not None for s in states):
        report, matrix = _final_payloads(plan, states, counts)
        put_identical(store, "host_vectors.npy", matrix)
        put_identical(store, "report.json", json_bytes(report))
    status = verify_extraction(plan, store)
    if stop_reason is not None:
        status["stop_reason"] = stop_reason
    return status
