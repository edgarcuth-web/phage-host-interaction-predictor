"""Registered inheritance, immutable feature reuse, and the frozen manifest guard."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import sha256_file
from phage_cocktail import genophi_comparator as comparator
from phage_cocktail.genophi_seed_protocols import (
    BASE_COMPARATOR, BASE_GENOME, BASE_HOST, CACHE_LOCK, INPUT_FREEZE, REGISTRATION,
    expected_variant, import_feature_cache, load_study_protocol, seed_paths,
    assert_baseline_model_split,
)
from phage_cocktail.retrieval_study import retrieval_splits

REPO = Path(__file__).resolve().parents[1]


def test_registered_variants_retain_original_numeric_decisions():
    original = json.loads((REPO / BASE_COMPARATOR).read_text())
    for seed in (43, 44):
        _, plan = load_study_protocol(REPO, seed_paths(seed)["comparator"], "comparator")
        assert plan["scenarios"] == ["grouped_host"]
        for section in ("feature_selection", "modeling", "repeats", "prediction", "clustering", "compute", "evaluation"):
            assert plan[section] == original[section]
        parent = json.loads((REPO / plan["parent_protocol"]).read_text())
        assert parent["split"]["seed"] == seed
        _, baseline = load_study_protocol(REPO, seed_paths(seed)["genome"], "genome")
        old_baseline = json.loads((REPO / BASE_GENOME).read_text())
        for section in ("estimator", "representation", "categorical_features", "bootstrap", "numerical_environment"):
            assert baseline[section] == old_baseline[section]


def test_mutated_cutoff_or_wrong_seed_reference_is_rejected(tmp_path):
    files = [BASE_COMPARATOR, BASE_GENOME, BASE_HOST, REGISTRATION, INPUT_FREEZE,
             *[seed_paths(43)[k] for k in ("manifest", "host", "genome", "comparator")]]
    for relative in files:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((REPO / relative).read_bytes())
    path = tmp_path / seed_paths(43)["comparator"]
    original = json.loads(path.read_text())
    for change in ("cutoff", "reference"):
        modified = deepcopy(original)
        if change == "cutoff":
            modified["feature_selection"]["cutoff"] = 4
        else:
            modified["reference_report"] = "data/studies/results/genome_host_retrieval.v1/report.json"
        path.write_text(json.dumps(modified))
        with pytest.raises(ValueError, match="inherited variant"):
            load_study_protocol(tmp_path, path, "comparator")


def test_B2_is_identical_after_reassigning_the_grouped_manifest():
    table = pd.read_csv(REPO / "data/studies/results/genophi_comparator.v1/scenario_splits.csv", dtype=str)
    metadata = table[table.scenario.eq("grouped_host")].set_index("bacterium")[["phylogroup"]].rename(
        columns={"phylogroup": "Clermont_Phylo"})
    reference = None
    for seed in (42, 43, 44):
        path = ("data/studies/results/genophi_host_split.v1.csv" if seed == 42 else seed_paths(seed)["manifest"])
        manifest = pd.read_csv(REPO / path, dtype=str)
        result = retrieval_splits(manifest, metadata, seed=seed, fractions=[.7, .15, .15])["phylogroup:B2"]
        if reference is None:
            reference = result
        pd.testing.assert_frame_equal(result, reference)


def test_import_reuses_verified_bytes_and_rejects_tampered_cache(tmp_path):
    repo, source, destination = tmp_path / "repo", tmp_path / "source", tmp_path / "target"
    source.mkdir()
    (repo / "data/studies").mkdir(parents=True)
    (source / "checkpoints").mkdir()
    (source / "provenance.json").write_text("{}")
    (source / "checkpoints/features.json").write_text("{}")
    payload = source / "full_feature_table.csv"
    payload.write_text("strain,phage,interaction,sc_0\nA,P,1,1\n")
    lock = {"proteome_receipt_sha256": "proteome",
            "source_provenance_sha256": sha256_file(source / "provenance.json"),
            "source_checkpoint_sha256": {"features": sha256_file(source / "checkpoints/features.json")},
            "artifacts": {payload.name: {"sha256": sha256_file(payload), "size_bytes": payload.stat().st_size}}}
    (repo / CACHE_LOCK).write_text(json.dumps(lock))
    (repo / INPUT_FREEZE).write_text(json.dumps({"files_sha256": {CACHE_LOCK: sha256_file(repo / CACHE_LOCK)}}))
    outputs, record = import_feature_cache(repo, source, destination, "proteome")
    assert outputs[0].read_bytes() == payload.read_bytes()
    assert record["new_clustering_search_assignment_calls"] == 0
    payload.write_text("changed")
    with pytest.raises(ValueError, match="changed source feature"):
        import_feature_cache(repo, source, destination, "proteome")


def test_original_freeze_cohort_still_blocks_a_mismatched_manifest(tmp_path, monkeypatch):
    from phage_cocktail.cohort import canonical_assembly_digest

    resources = {}
    for name, content in [("subset", "strain,phage,interaction\nA,P,1\nB,P,0\nC,P,1\n"),
                           ("mlc", ";P\nA;1\nB;0\nC;1\n"),
                           ("core_distances", "\tA\tB\tC\nA\t0\t.01\t.01\nB\t.01\t0\t.01\nC\t.01\t.01\t0\n")]:
        path = tmp_path / f"{name}.csv"
        path.write_text(content)
        resources[("table", name)] = path
    for h, sequence in zip("ABC", ["ACGT", "AAGT", "AAAA"]):
        path = tmp_path / f"{h}.fna"
        path.write_text(f">contig\n{sequence}\n")
        resources[("host_genome", h)] = path
    parent = {"subset_manifest": "unused", "grouping": {"core_distance_threshold": 1e-4},
              "split": {"seed": 43, "fractions": [.7, .15, .15]}}
    (tmp_path / "parent.json").write_text(json.dumps(parent))
    manifest = pd.DataFrame({"bacterium": list("ABC"), "host_group": list("ABC"), "split": "train",
                             "source_bacterium": list("ABC"),
                             "assembly_sha256": [canonical_assembly_digest(resources[("host_genome", h)]) for h in "ABC"]})
    manifest.to_csv(tmp_path / "wrong.csv", index=False)
    monkeypatch.setattr(comparator, "verify_resources", lambda *a, **k: resources)
    monkeypatch.setattr(comparator, "verify_registered_dataset", lambda *a: None)
    with pytest.raises(comparator.ComparatorBlocked, match="recomputed host split differs from the frozen study"):
        comparator._freeze_cohort(tmp_path, tmp_path, tmp_path / "run",
                                 {"parent_protocol": "parent.json", "host_manifest": "wrong.csv"}, {})


def test_baseline_fitted_under_one_seed_cannot_be_scored_as_another_seed_refit():
    from phage_cocktail.host_baseline import fit_host_model
    from phage_cocktail.retrieval_study import score_known_panel

    features = pd.DataFrame({"f": [0., 1., 2., 3.]}, index=list("ABCD"))
    old_training = pd.DataFrame({"bacterium": ["A", "B"], "phage": "P", "infection": [0, 1]})
    old_model = fit_host_model(old_training, features, ["f"], [], C=.1, max_iter=2000)
    new_manifest = pd.DataFrame({"bacterium": list("ABCD"), "split": ["train", "test", "train", "validation"]})
    with pytest.raises(ValueError, match="this seed's frozen training partition"):
        assert_baseline_model_split(old_model, new_manifest)
    candidates = pd.DataFrame({"bacterium": ["B", "D"], "phage": "P"})
    with pytest.raises(ValueError, match="disjoint"):
        score_known_panel(old_training, candidates, features, ["f"], [], {"C": .1, "max_iter": 2000})


def test_failed_seed_is_reported_and_does_not_cancel_the_other_seed(tmp_path):
    import importlib.util
    from types import SimpleNamespace

    spec = importlib.util.spec_from_file_location("seed_driver", REPO / "scripts/run_genophi_seed_sensitivity.py")
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    assert driver.interpreter_path(tmp_path, Path(".venv/bin/python")) == tmp_path / ".venv/bin/python"
    calls = []

    def run(command, **kwargs):
        calls.append(command[0])
        return SimpleNamespace(returncode=1 if command[0] == "fail43" else 0)

    jobs = {43: [("baseline", ["fail43"]), ("comparator", ["never43"])],
            44: [("baseline", ["baseline44"]), ("comparator", ["comparator44"])]}
    root = tmp_path / "run"
    value = driver.execute_studies(root, jobs, repo=tmp_path, run=run)
    assert calls == ["fail43", "baseline44", "comparator44"]
    assert value["seeds"]["43"]["status"] == "blocked"
    assert value["seeds"]["44"]["status"] == "verified"
    with pytest.raises(ValueError, match="do not rerun"):
        driver.execute_studies(root, jobs, repo=tmp_path, run=run)


def test_sensitivity_stage_imports_cache_and_never_calls_native_feature_workflows(tmp_path, monkeypatch):
    from phage_cocktail import genophi_seed_protocols as protocols

    protocol = tmp_path / "protocol.json"
    plan = {"split_seed": 43, "repo_files_sha256": {}, "resource_lock": "lock.json",
            "versions": {"genophi": "1.0.3"}, "reference_report": "reference.json"}
    protocol.write_text(json.dumps(plan))
    original_reader = comparator.read_json_object

    def read(path):
        if Path(path).name == "proteomes.json":
            return {"resource_lock_sha256": "digest"}
        if Path(path).name == "lock.json":
            return {}
        if Path(path) == tmp_path / CACHE_LOCK:
            return {"versions": plan["versions"]}
        return original_reader(path)

    monkeypatch.setattr(protocols, "load_study_protocol", lambda *a: (protocol, plan))
    monkeypatch.setattr(protocols, "verify_seed_baseline_reference", lambda *a: {"status": "verified"})
    monkeypatch.setattr(comparator, "read_json_object", read)
    monkeypatch.setattr(comparator, "sha256_file", lambda *a: "digest")
    monkeypatch.setattr(comparator, "verify_proteomes", lambda *a: None)
    monkeypatch.setattr(comparator, "runtime_identity", lambda: {"versions": plan["versions"]})
    monkeypatch.setattr(comparator, "code_provenance", lambda *a: {})
    monkeypatch.setattr(comparator, "_freeze_cohort", lambda *a: ([], {"status": "frozen"}))
    for name in ("_cluster_set", "_assign_set", "_merge_features"):
        monkeypatch.setattr(comparator, name, lambda *a: pytest.fail("sensitivity rebuilt feature stages"))
    imported = []

    def import_cache(repo, source, destination, receipt):
        imported.append(source)
        path = destination / "full_feature_table.csv"
        path.write_text("strain,phage,interaction,sc_0\nA,P,1,1\n")
        return [path], {"new_clustering_search_assignment_calls": 0}

    monkeypatch.setattr(protocols, "import_feature_cache", import_cache)
    comparator.run_comparator(tmp_path, tmp_path, tmp_path, tmp_path / "proteomes.json", tmp_path / "run",
                              stage="cluster", protocol=protocol, feature_cache=tmp_path / "source")
    assert imported == [tmp_path / "source"]
    assert (tmp_path / "run/checkpoints/reuse_features.json").is_file()
