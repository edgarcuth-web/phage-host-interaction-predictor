"""Budget, durable resume and real-assembly admission boundaries; CPU fixtures only."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import atomic_write_json, sha256_file
from phage_cocktail.cohort import canonical_assembly_digest
from phage_cocktail.nvidia_evo2 import ForwardResult, NvidiaForwardError
from phage_cocktail.nvidia_host_pilot import PROTOCOL, run_pilot


@pytest.fixture
def pilot_fixture(tmp_path, monkeypatch):
    repo, source, cache, output = [tmp_path / n for n in ("repo", "source", "feasibility", "pilot")]
    for p in (repo / "data/studies/results", source, cache):
        p.mkdir(parents=True)
    manifest, resources = [], []
    for i, host in enumerate("ABC"):
        path = source / (host + ".fasta")
        path.write_text(">fixture\n" + "ACGT" * (i + 3) + "NN\n")
        manifest.append({"bacterium": host, "assembly_sha256": canonical_assembly_digest(path)})
        resources.append({"kind": "host_genome", "id": host, "path": path.name,
                          "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    mp = repo / "data/studies/hosts.csv"
    pd.DataFrame(manifest).to_csv(mp, index=False)
    gp = repo / "data/studies/genomes.json"
    atomic_write_json({"host_manifest": "data/studies/hosts.csv", "resources": resources}, gp)
    from phage_cocktail.evo2_host import host_windows
    seq = next(s for _, _, s in host_windows((source / "A.fasta").read_bytes(), window_bp=4, max_windows_per_host=32) if len(s)==4)
    import hashlib
    np.save(cache / "pooled_mean.npy", np.array([1., 2.], dtype=np.float32), allow_pickle=False)
    parent = {"status": "success", "provenance": {"input": {"input_sha256": hashlib.sha256(seq.encode()).hexdigest()}},
              "response": {"response_bytes": 10},
              "output_sha256": {"pooled_mean.npy": sha256_file(cache / "pooled_mean.npy")}}
    pp = repo / "data/studies/results/feasibility.json"
    atomic_write_json(parent, pp)
    atomic_write_json(parent, cache / "report.json")
    plan = {"hosts": list("ABC"), "genome_protocol": "data/studies/genomes.json",
            "feasibility_report": "data/studies/results/feasibility.json", "window_bp": 4,
            "sampling_max_windows": 32, "dimension": 2,
            "repo_files_sha256": {str(p.relative_to(repo)): sha256_file(p) for p in (mp,gp,pp)},
            "request_budget": {"max_new_response_bytes": 100, "max_new_inference_posts": 5,
                               "max_response_bytes": 30, "max_tensor_bytes": 1000,
                               "io_timeout_seconds": 180, "server_wait_seconds": 60}}
    atomic_write_json(plan, repo / PROTOCOL)
    monkeypatch.setattr("phage_cocktail.nvidia_host_pilot.code_provenance",
                        lambda *a: {"code_sha256": {"fixture": "frozen"}, "code_base_revision": "fixture"})
    return repo, source, cache, output, plan


def test_five_requests_one_reuse_pooling_and_zero_request_resume(pilot_fixture):
    repo, source, cache, output, _ = pilot_fixture
    calls = []
    def encode(seq, key, **kwargs):
        calls.append(kwargs)
        return ForwardResult(np.tile([3., 4.], (len(seq),1)), {"response_bytes": 10})
    report = run_pilot(repo, source, output, "fixture-secret", cache, call_forward=encode)
    assert len(calls) == 5 and report["new_response_bytes"] == 50
    assert report["host_vector_shape"] == [3,2]
    np.testing.assert_array_equal(np.load(output / "host_vectors.npy"), [[2.,3.],[3.,4.],[3.,4.]])
    assert "fixture-secret" not in (output / "report.json").read_text()
    repeated = run_pilot(repo, source, output, None, cache, resume=True,
                         call_forward=lambda *a,**k: pytest.fail("completed cache must not call API"))
    assert repeated == report
    (output / "job_02.npy").write_bytes(b"changed")
    with pytest.raises(ValueError, match="output hash"):
        run_pilot(repo, source, output, None, cache, resume=True)


def test_unresolved_attempt_blocks_resume_without_duplicate_inference(pilot_fixture):
    repo, source, cache, output, _ = pilot_fixture
    calls = []
    def fail(*args, **kwargs):
        calls.append(1)
        raise NvidiaForwardError("pending", http_status=202, request_id="fixture-id")
    with pytest.raises(NvidiaForwardError):
        run_pilot(repo, source, output, "fixture-secret", cache, call_forward=fail)
    assert calls == [1]
    assert json.loads((output / "journal.json").read_text())["new_post_attempts"] == 1
    with pytest.raises(ValueError, match="unresolved"):
        run_pilot(repo, source, output, "fixture-secret", cache, resume=True, call_forward=fail)
    assert calls == [1]


def test_transfer_budget_limits_next_download_and_blocks_further_calls(pilot_fixture):
    repo, source, cache, output, plan = pilot_fixture
    plan["request_budget"]["max_new_response_bytes"] = 15
    atomic_write_json(plan, repo / PROTOCOL)
    limits = []
    def encode(seq,key,**kwargs):
        limits.append(kwargs["max_response_bytes"])
        if len(limits)==2:
            raise NvidiaForwardError("download limit")
        return ForwardResult(np.tile([3.,4.], (len(seq),1)), {"response_bytes": 10})
    with pytest.raises(NvidiaForwardError):
        run_pilot(repo, source, output, "fixture-secret", cache, call_forward=encode)
    assert limits == [15,5]


def test_changed_source_is_rejected_before_network(pilot_fixture):
    repo, source, cache, output, _ = pilot_fixture
    (source / "B.fasta").write_text(">changed\nACGT\n")
    with pytest.raises(ValueError, match="locked resource"):
        run_pilot(repo, source, output, "fixture-secret", cache,
                  call_forward=lambda *a,**k: pytest.fail("changed input must not call API"))
    assert not output.exists()
