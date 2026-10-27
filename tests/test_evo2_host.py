"""CPU contract tests only: these fixtures are not Evo2 performance evidence."""
import json

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import sha256_file
from phage_cocktail.cohort import canonical_assembly_digest, reverse_complement
from phage_cocktail.evo2_host import PROTOCOL, extract_cache, host_windows, load_cache
from phage_cocktail.host_baseline import HostFeatureEncoder
from phage_cocktail.retrieval_study import run_retrieval_study
from test_retrieval_study import local_study


def test_window_sampling_preserves_boundaries_ambiguity_and_strand_identity():
    payload = b">one\nACGTRNACGT\n>two\nTTTTCC\n"
    w = host_windows(payload, window_bp=4, max_windows_per_host=32)
    assert len(w) == 5 and sum(len(s) for _, _, s in w) == 16
    assert any("N" in s for _, _, s in w)
    assert all(set(s) <= set("ACGTN") and len(s) <= 4 for _, _, s in w)
    altered = (">two\n" + reverse_complement("TTTTCC") + "\n>one\n" + reverse_complement("ACGTRNACGT") + "\n").encode()
    assert host_windows(altered, window_bp=4, max_windows_per_host=32) == w
    assert host_windows(payload, window_bp=4, max_windows_per_host=3) == [w[0], w[2], w[4]]
    with pytest.raises(ValueError, match="positive integers"):
        host_windows(payload, window_bp=0, max_windows_per_host=32)


def test_signed_embeddings_use_train_only_statistics_without_relaxing_counts():
    train = pd.DataFrame({"evo2:0": [-2., 2.], "defense:RM": [0., 4.]}, index=["A", "B"])
    encoder = HostFeatureEncoder.fit(train, list(train), [])
    assert encoder.means == [0., 2.]
    np.testing.assert_allclose(encoder.transform(train), [[-1, -1], [1, 1]])
    with pytest.raises(ValueError, match="non-negative"):
        encoder.transform(train.assign(**{"defense:RM": -1}))
    with pytest.raises(ValueError, match="finite"):
        encoder.transform(train.assign(**{"evo2:0": np.nan}))


@pytest.fixture
def fixture_cache(local_study, monkeypatch):
    repo, source, output, plan = local_study
    manifest_path = repo / plan["host_manifest"]
    manifest = pd.read_csv(manifest_path)
    resources = []
    for i, host in enumerate(manifest.bacterium):
        path = source / f"{host}.fasta"
        path.write_text(">fixture\n" + "ACGT" * (i + 3) + "NNNNATT\n")
        manifest.loc[manifest.bacterium.eq(host), "assembly_sha256"] = canonical_assembly_digest(path)
        resources.append({"kind": "host_genome", "id": host, "path": path.name,
                          "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    manifest.to_csv(manifest_path, index=False)
    gp = repo / "data/studies/genome_host_retrieval.v1.json"
    gp.write_text(json.dumps({"resources": resources}))
    distance = source / "core.tsv"
    pd.DataFrame(np.abs(np.arange(9)[:, None] - np.arange(9)[None, :]),
                 index=manifest.bacterium, columns=manifest.bacterium).to_csv(distance, sep="\t")
    plan["resources"].append({"kind": "table", "id": "core_distances", "path": distance.name,
                              "size_bytes": distance.stat().st_size, "sha256": sha256_file(distance)})
    plan["repo_files_sha256"].update({str(manifest_path.relative_to(repo)): sha256_file(manifest_path),
                                      str(gp.relative_to(repo)): sha256_file(gp)})
    plan["embedding"] = {"dimension": 2, "window_bp": 4, "max_windows_per_host": 3}
    plan["neighbours"] = {"k": 5}
    (repo / PROTOCOL).write_text(json.dumps(plan))
    identity = {key: "a" * 64 for key in ("checkpoint_sha256", "config_sha256", "evo2_module_sha256", "vortex_module_sha256")}
    identity.update({key: "CPU fixture" for key in ("evo2_version", "vortex_version", "torch_version", "cuda_version", "transformer_engine_version", "hardware", "precision")})
    identity["model_name"] = "evo2_40b"
    monkeypatch.setattr("phage_cocktail.evo2_host.code_provenance", lambda *args: {"code_base_revision": "fixture"})
    # Constant signed feature checks extraction/pooling and fallback plumbing only.
    def fixture_encoder(s):
        return np.tile(np.array([-1., 2.], dtype=np.float32), (len(s), 1))
    return repo, source, output, manifest, identity, fixture_encoder


def test_cache_full_roundtrip_and_all_representation_scoring(fixture_cache):
    repo, source, output, manifest, identity, encode = fixture_cache
    cache = output.with_name("cache")
    extract_cache(repo, source, cache, encode, identity)
    features, fourmers, qc, _ = load_cache(cache, repo, manifest)
    np.testing.assert_array_equal(features.to_numpy(), np.tile([-1., 2.], (9, 1)))
    assert fourmers.shape == (9, 136)
    assert qc.sampled_fraction.between(0, 1).all()
    report = run_retrieval_study(repo, source, output, embedding_cache=cache)
    assert len(report["metrics"]) == 63  # Nine methods, seven scenario/partition panels.
    assert "Evo2" not in report["not_run"]
    scores = pd.read_csv(output / "predictions.csv.gz")
    assert set(scores.groupby(["scenario", "partition", "model", "bacterium"]).size()) == {3}
    (cache / "embeddings.npy").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_cache(cache, repo, manifest)


def test_pilot_cannot_be_used_as_full_evaluation(fixture_cache):
    repo, source, output, manifest, identity, encode = fixture_cache
    extract_cache(repo, source, output, encode, identity, pilot=True)
    assert len(pd.read_csv(output / "qc.csv")) == 3
    with pytest.raises(ValueError, match="not a pilot"):
        load_cache(output, repo, manifest)


def test_pooling_weights_short_windows_and_both_orientations(fixture_cache):
    repo, source, output, manifest, identity, _ = fixture_cache
    def encode(s):
        return np.array([["ACGTN".index(base), position] for position, base in enumerate(s)], dtype=np.float32)
    extract_cache(repo, source, output, encode, identity)
    features, _, _, _ = load_cache(output, repo, manifest)
    host = manifest.bacterium.iloc[0]
    windows = host_windows((source / f"{host}.fasta").read_bytes(), window_bp=4, max_windows_per_host=3)
    assert len({len(s) for _, _, s in windows}) > 1
    all_token_vectors = np.concatenate([encode(strand) for _, _, s in windows
                                        for strand in (s, reverse_complement(s))])
    np.testing.assert_allclose(features.loc[host], all_token_vectors.mean(axis=0, dtype=np.float64), atol=1e-6)


def test_bad_layer_or_assembly_cannot_publish_partial_cache(fixture_cache):
    repo, source, output, manifest, identity, _ = fixture_cache
    with pytest.raises(ValueError, match="dimension"):
        extract_cache(repo, source, output, lambda s: np.ones((len(s), 3)), identity)
    assert not output.exists()
    with pytest.raises(ValueError, match="real checkpoint_sha256"):
        extract_cache(repo, source, output, lambda s: None,
                      {**identity, "checkpoint_sha256": "unknown"})
    assert not output.exists()


def test_hosted_evo2_uses_distinct_methods_and_provenance(fixture_cache, monkeypatch):
    repo, source, output, manifest, identity, encode = fixture_cache
    local_cache = output.with_name("local_cache_for_hosted_fixture")
    extract_cache(repo, source, local_cache, encode, identity)
    features, fourmers, qc, _ = load_cache(local_cache, repo, manifest)
    hosted_features = features.copy()
    hosted_features.columns = [f"evo2:hosted:{i}" for i in range(hosted_features.shape[1])]
    representation = {
        "source": "nvidia_hosted_forward",
        "model_name": "evo2_40b",
        "endpoint": "fixture",
        "layer": "fixture",
        "dimension": hosted_features.shape[1],
        "selection": "original_32_including_short_tails",
        "pooling": "fixture hosted pooling",
        "weights_sha256": None,
        "server_runtime_version": None,
        "local_v1_cache_compatible": False,
        "source_plan_sha256": "a" * 64,
    }

    def fake_hosted_loader(*args, **kwargs):
        return hosted_features, fourmers, qc, "b" * 64, representation

    monkeypatch.setattr(
        "phage_cocktail.retrieval_study.load_hosted_cache", fake_hosted_loader
    )
    hosted_plan = output.with_name("hosted_plan.json")
    hosted_plan.write_text("{}")
    report = run_retrieval_study(
        repo,
        source,
        output,
        hosted_embedding=(hosted_plan, object()),
    )
    models = {row["model"] for row in report["metrics"]}
    assert "host_hosted_evo2_logistic" in models
    assert "hybrid_hosted_evo2_logistic" in models
    assert "host_evo2_logistic" not in models
    assert report["provenance"]["hosted_embedding_sha256"] == "b" * 64
    assert "Evo2" not in report["not_run"]
    protocol = json.loads((output / "protocol.json").read_text())
    assert protocol["study_id"] == "hosted-evo2-host-representation:v1"
    assert protocol["embedding"]["source"] == "nvidia_hosted_forward"
    assert protocol["embedding"]["local_v1_cache_compatible"] is False
