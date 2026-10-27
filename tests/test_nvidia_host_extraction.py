"""CPU fixtures exercise paid-call boundaries; no real NVIDIA or S3 requests."""
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import atomic_write_json, sha256_file
from phage_cocktail.cohort import canonical_assembly_digest
from phage_cocktail.extraction_store import LocalStore, S3Store, StoreError
from phage_cocktail.nvidia_evo2 import ForwardResult, NvidiaForwardError
from phage_cocktail import nvidia_host_extraction as ex


@pytest.fixture
def extraction(tmp_path, monkeypatch):
    # Two feature columns keep the synthetic CPU contract tests small.
    monkeypatch.setattr(ex, "DIMENSION", 2)
    repo, source = tmp_path / "repo", tmp_path / "source"
    (repo / "data/studies/results/nvidia_evo2_host_pilot.v1").mkdir(parents=True)
    source.mkdir()
    manifest, resources = [], []
    for host, base in zip("ABC", "ACG"):
        path = source / f"{host}.fasta"
        path.write_text(">fixture\n" + base * 10240 + "\n")
        manifest.append({"bacterium": host, "assembly_sha256": canonical_assembly_digest(path)})
        resources.append({"kind": "host_genome", "id": host, "path": path.name,
                          "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    pd.DataFrame(manifest).to_csv(repo / "data/studies/hosts.csv", index=False)
    atomic_write_json({"host_manifest": "data/studies/hosts.csv", "resources": resources}, repo / ex.GENOMES)
    atomic_write_json({}, repo / ex.PILOT)
    (repo / "fixture.py").write_text("# frozen fixture\n")
    monkeypatch.setattr(ex, "code_provenance", lambda repo, script: {
        "code_base_revision": "cpu-fixture", "code_sha256": {"fixture.py": sha256_file(repo / "fixture.py")}})
    return repo, source, LocalStore(tmp_path / "store")


def plan_for(fixture, calls=6, total=600, selection="first_full_window"):
    repo, source, _ = fixture
    plan = ex.make_plan(repo, source, max_new_calls=calls, max_new_response_bytes=total, selection=selection)
    plan["budget"]["max_response_bytes"] = 99
    return plan


def response(sequence, values):
    return ForwardResult(np.broadcast_to(np.array(values, dtype=np.float64), (len(sequence), 2)), {
        "endpoint": ex.FORWARD_URL, "layer": ex.LAYER, "http_status": 200, "finite": True,
        "shape": [len(sequence), 2], "input_sha256": sha256(sequence.encode()).hexdigest(),
        "sequence_length": len(sequence), "response_bytes": 10, "weights_sha256": None})


def pilot_cache(fixture):
    repo, source, _ = fixture
    jobs, sequences = ex.prepare_inputs(repo, source, list("ABC"), "first_full_window")
    cache = source.parent / "pilot"
    cache.mkdir()
    records, hashes, vectors = [], {}, []
    for i, job in enumerate(jobs):
        path = cache / f"job_{i:02d}.npy"
        vector = np.array([i+1, i+2], dtype=np.float32)
        np.save(path, vector, allow_pickle=False)
        hashes[path.name] = sha256_file(path)
        records.append({"input": {k:v for k,v in job.items() if k != "id"},
                        "pooled_sha256": hashes[path.name],
                        "response": response(sequences[job["id"]], vector).metadata})
        vectors.append(vector)
    report = {"status": "complete_engineering_pilot", "jobs": records, "output_sha256": hashes}
    atomic_write_json(report, repo / ex.PILOT)
    atomic_write_json(report, cache / "report.json")
    return cache


def forbidden(*args, **kwargs):
    pytest.fail("no new request or key access allowed")


def test_default_plan_and_run_cannot_call_api_or_request_key(extraction):
    repo, source, store = extraction
    plan = ex.make_plan(repo, source)
    assert plan["budget"]["max_new_calls"] == plan["budget"]["max_new_response_bytes"] == 0
    result = ex.run_extraction(repo, source, plan, store, key_provider=forbidden, call_forward=forbidden)
    assert result["completed_jobs"] == result["new_post_attempts"] == 0
    assert store.get("jobs/00000/claim.json") is None
    # Even --allow-live cannot increase a zero frozen budget.
    ex.run_extraction(repo, source, plan, store, allow_live=True, key_provider=forbidden, call_forward=forbidden)


def test_checkpoint_resume_reuses_completed_calls_and_exports_exact_means(extraction):
    repo, source, store = extraction
    plan, calls = plan_for(extraction), []
    def encode(seq, key, **kwargs):
        calls.append(kwargs)
        assert key == "fixture-secret"
        return response(seq, [len(calls), len(calls)+1])
    first = ex.run_extraction(repo, source, plan, store, allow_live=True,
                              key_provider=lambda: "fixture-secret", call_forward=encode, max_jobs=2)
    assert first["completed_jobs"] == first["new_post_attempts"] == 2
    result = ex.run_extraction(repo, source, plan, store, allow_live=True,
                               key_provider=lambda: "fixture-secret", call_forward=encode)
    assert len(calls) == 6 and result["status"] == "complete_engineering_extraction"
    assert result["reserved_response_bytes"] == 600 and result["successful_new_response_bytes"] == 60
    assert all(c["max_response_bytes"] == 99 for c in calls)
    np.testing.assert_array_equal(np.load(io.BytesIO(store.get("host_vectors.npy"))),
                                  [[1.5,2.5], [3.5,4.5], [5.5,6.5]])
    snapshot = {p.relative_to(store.root).as_posix(): p.read_bytes() for p in store.root.rglob("*") if p.is_file()}
    ex.run_extraction(repo, source, plan, store, key_provider=forbidden, call_forward=forbidden)
    assert snapshot == {p.relative_to(store.root).as_posix(): p.read_bytes() for p in store.root.rglob("*") if p.is_file()}
    assert all(b"fixture-secret" not in payload for payload in snapshot.values())


@pytest.mark.parametrize("calls,total,expected,reserved", [(2,600,2,200), (6,250,3,250)])
def test_both_frozen_budgets_survive_resume(extraction, calls, total, expected, reserved):
    repo, source, store = extraction
    plan, submitted = plan_for(extraction, calls, total), []
    def encode(seq, key, **kwargs):
        submitted.append(kwargs["max_response_bytes"])
        return response(seq, [1.,2.])
    for _ in range(2):
        result = ex.run_extraction(repo, source, plan, store, allow_live=True,
                                   key_provider=lambda: "fixture", call_forward=encode)
    assert len(submitted) == expected
    assert result["reserved_response_bytes"] == reserved
    assert result["status"] == "incomplete"


def test_ambiguous_post_blocks_whole_run_without_retry(extraction):
    repo, source, store = extraction
    plan, calls = plan_for(extraction), []
    def fail(*args, **kwargs):
        assert store.get("jobs/00000/claim.json") is not None
        calls.append(1)
        raise NvidiaForwardError("HTTP 202", http_status=202)
    with pytest.raises(NvidiaForwardError):
        ex.run_extraction(repo, source, plan, store, allow_live=True,
                           key_provider=lambda: "fixture", call_forward=fail)
    with pytest.raises(ValueError, match="unresolved live"):
        ex.run_extraction(repo, source, plan, store, allow_live=True,
                           key_provider=forbidden, call_forward=forbidden)
    assert calls == [1]
    state = ex.verify_extraction(plan, store)
    assert state["unresolved_live_jobs"] == ["00000"] and state["new_post_attempts"] == 1


class FailReceipt:
    def __init__(self, store):
        self.store = store
    def get(self, key):
        return self.store.get(key)
    def create(self, key, payload):
        if key.endswith("result.json"):
            raise StoreError("fixture storage outage")
        return self.store.create(key, payload)


def test_storage_failure_after_inference_does_not_resubmit(extraction):
    repo, source, store = extraction
    plan, calls = plan_for(extraction), []
    def encode(seq, key, **kwargs):
        calls.append(1)
        return response(seq, [1.,2.])
    with pytest.raises(StoreError):
        ex.run_extraction(repo, source, plan, FailReceipt(store), allow_live=True,
                           key_provider=lambda: "fixture", call_forward=encode)
    with pytest.raises(ValueError, match="unresolved live"):
        ex.run_extraction(repo, source, plan, store, allow_live=True,
                           key_provider=forbidden, call_forward=forbidden)
    assert calls == [1]


def test_interrupted_pilot_import_can_finish_without_key_or_inference(extraction):
    repo, source, store = extraction
    cache = pilot_cache(extraction)
    plan = plan_for(extraction, 0, 0)
    with pytest.raises(StoreError):
        ex.run_extraction(repo, source, plan, FailReceipt(store), pilot_cache=cache,
                           key_provider=forbidden, call_forward=forbidden)
    result = ex.run_extraction(repo, source, plan, store, pilot_cache=cache,
                               key_provider=forbidden, call_forward=forbidden)
    assert result["status"] == "complete_engineering_extraction" and result["new_post_attempts"] == 0
    np.testing.assert_array_equal(np.load(io.BytesIO(store.get("host_vectors.npy"))),
                                  [[1.5,2.5], [3.5,4.5], [5.5,6.5]])


def test_final_export_failure_is_repaired_from_complete_receipts(extraction):
    repo, source, store = extraction
    cache = pilot_cache(extraction)
    plan = plan_for(extraction, 0, 0)
    class FailReport:
        def get(self, key):
            return store.get(key)
        def create(self, key, payload):
            if key == "report.json":
                raise StoreError("fixture report upload failed")
            return store.create(key, payload)
    with pytest.raises(StoreError):
        ex.run_extraction(repo, source, plan, FailReport(), pilot_cache=cache, call_forward=forbidden)
    assert ex.verify_extraction(plan, store)["status"] == "ready_to_finalize"
    result = ex.run_extraction(repo, source, plan, store, key_provider=forbidden, call_forward=forbidden)
    assert result["status"] == "complete_engineering_extraction"


@pytest.mark.parametrize("target", ["vector", "final", "claim", "plan"])
def test_corruption_is_rejected_before_any_new_call(extraction, target):
    repo, source, store = extraction
    cache = pilot_cache(extraction)
    plan = plan_for(extraction, 0, 0)
    ex.run_extraction(repo, source, plan, store, pilot_cache=cache, call_forward=forbidden)
    if target == "vector":
        path = store.root / "jobs/00000/result.json"
        record = json.loads(path.read_text())
        record["vector_sha256"] = "0" * 64
        path.write_text(json.dumps(record))
    elif target == "final":
        (store.root / "host_vectors.npy").write_bytes(b"corrupted")
    elif target == "claim":
        path = store.root / "jobs/00000/claim.json"
        record = json.loads(path.read_text())
        record["input"]["length"] -= 1
        path.write_text(json.dumps(record))
    else:
        plan["budget"]["max_new_calls"] = 1
    with pytest.raises(ValueError):
        ex.run_extraction(repo, source, plan, store, allow_live=True,
                           key_provider=forbidden, call_forward=forbidden)


@pytest.mark.parametrize("target", ["code", "fasta"])
def test_changed_inputs_rejected_before_checkpoint_or_network(extraction, target):
    repo, source, store = extraction
    plan = plan_for(extraction)
    path = repo / "fixture.py" if target == "code" else source / "A.fasta"
    path.write_text("changed")
    with pytest.raises(ValueError):
        ex.run_extraction(repo, source, plan, store, allow_live=True,
                           key_provider=forbidden, call_forward=forbidden)
    assert not store.root.exists()


def test_losing_a_claim_race_never_submits_inference(extraction):
    repo, source, store = extraction
    plan = plan_for(extraction)
    class OtherWriter:
        def get(self, key):
            return store.get(key)
        def create(self, key, payload):
            if key.endswith("claim.json"):
                store.create(key, payload)
                return False
            return store.create(key, payload)
    with pytest.raises(ValueError, match="another worker"):
        ex.run_extraction(repo, source, plan, OtherWriter(), allow_live=True,
                           key_provider=lambda: "fixture", call_forward=forbidden)


def test_short_tail_is_weighted_by_length_after_strand_mean(extraction):
    repo, source, store = extraction
    plan, calls = plan_for(extraction, 12, 1200, "original_32_including_short_tails"), []
    def encode(seq, key, **kwargs):
        values = [2.,4.] if len(calls) % 4 < 2 else [10.,0.]
        calls.append(len(seq))
        return response(seq, values)
    result = ex.run_extraction(repo, source, plan, store, allow_live=True,
                               key_provider=lambda: "fixture", call_forward=encode)
    assert len(calls) == 12 and result["status"] == "complete_engineering_extraction"
    np.testing.assert_array_equal(np.load(io.BytesIO(store.get("host_vectors.npy"))),
                                  np.array([[3.6,3.2]] * 3, dtype=np.float32))


class S3Failure(Exception):
    def __init__(self, code):
        self.response = {"Error": {"Code": code}}
        super().__init__("fixture-secret must not escape")


def test_s3_conditional_writes_and_errors_fail_closed():
    class Client:
        def __init__(self):
            self.objects, self.calls = {}, []
            self.error = None
        def put_object(self, **kwargs):
            self.calls.append(kwargs)
            if self.error:
                raise S3Failure(self.error)
            assert kwargs["IfNoneMatch"] == "*"
            key = kwargs["Key"]
            if key in self.objects:
                raise S3Failure("PreconditionFailed")
            self.objects[key] = kwargs["Body"]
        def get_object(self, **kwargs):
            if self.error:
                raise S3Failure(self.error)
            if kwargs["Key"] not in self.objects:
                raise S3Failure("NoSuchKey")
            return {"Body": io.BytesIO(self.objects[kwargs["Key"]])}
    client = Client()
    store = S3Store("fixture-bucket", "run-1", client=client)
    assert store.get("plan.json") is None
    assert store.create("plan.json", b"one") and not store.create("plan.json", b"two")
    assert store.get("plan.json") == b"one"
    store.check_conditional_writes()
    client.error = "AccessDenied"
    with pytest.raises(StoreError) as error:
        store.get("plan.json")
    assert "fixture-secret" not in str(error.value)
    client.error = "ConditionalRequestConflict"
    before = len(client.calls)
    with pytest.raises(StoreError):
        store.create("plan.json", b"three")
    assert len(client.calls) == before + 1


def test_s3_backend_that_ignores_conditions_is_rejected():
    class UnsafeClient:
        def __init__(self):
            self.objects = {}
        def put_object(self, **kwargs):
            self.objects[kwargs["Key"]] = kwargs["Body"]
    with pytest.raises(StoreError, match="conditional-create"):
        S3Store("fixture", "run", client=UnsafeClient()).check_conditional_writes()


def test_real_sdk_accepts_conditional_put_contract_without_network():
    boto3 = pytest.importorskip("boto3")
    from botocore.config import Config
    from botocore.stub import Stubber
    client = boto3.client("s3", region_name="us-east-1", aws_access_key_id="fixture",
                          aws_secret_access_key="fixture", config=Config(retries={"total_max_attempts": 1}))
    request = {"Bucket": "fixture-bucket", "Key": "run/plan.json", "Body": b"one",
               "IfNoneMatch": "*", "ContentType": "application/octet-stream"}
    with Stubber(client) as stub:
        stub.add_response("put_object", {}, request)
        stub.add_client_error("put_object", service_error_code="PreconditionFailed", http_status_code=412,
                              expected_params=request)
        store = S3Store("fixture-bucket", "run", client=client)
        assert store.create("plan.json", b"one")
        assert not store.create("plan.json", b"one")
        stub.assert_no_pending_responses()


def test_local_store_does_not_overwrite_or_allow_escape(tmp_path):
    store = LocalStore(tmp_path / "objects")
    assert store.create("jobs/one", b"first")
    assert not store.create("jobs/one", b"second")
    assert store.get("jobs/one") == b"first"
    for key in ("../outside", "/absolute", "jobs//one", "jobs/./one", "jobs\\one"):
        with pytest.raises(ValueError):
            store.create(key, b"bad")
