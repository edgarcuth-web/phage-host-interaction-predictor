"""Failure boundaries and train-only fitting for the retrieval runner."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.analysis import summarize_interactions
from phage_cocktail.artifacts import sha256_file
from phage_cocktail.retrieval_study import retrieval_splits, score_known_panel, run_retrieval_study

ROOT = Path(__file__).resolve().parents[1]


def test_phylogroup_holdout_purges_touching_near_clonal_groups():
    manifest = pd.DataFrame({"bacterium": list("ABCDEF"), "host_group": list("AACCDE"),
                             "split": ["train"] * 4 + ["validation", "test"]})
    metadata = pd.DataFrame({"Clermont_Phylo": list("XYXXYY")}, index=list("ABCDEF"))
    splits = retrieval_splits(manifest, metadata, seed=42, fractions=[0.7, 0.15, 0.15])
    x = splits["phylogroup:X"].set_index("bacterium")
    assert x.loc["B", "split"] == "purged"
    assert set(x.loc[x.split.eq("test"), "phylogroup"]) == {"X"}
    for name, table in splits.items():
        if name != "random_host":
            train = set(table.loc[table.split.eq("train"), "host_group"])
            test = set(table.loc[table.split.eq("test"), "host_group"])
            assert not train & test
    shuffled = retrieval_splits(manifest.sample(frac=1, random_state=1), metadata,
                                 seed=42, fractions=[0.7, 0.15, 0.15])
    for name in splits:
        pd.testing.assert_frame_equal(splits[name], shuffled[name])


def test_candidate_scoring_uses_train_labels_only_and_no_test_metadata_fit():
    train = pd.DataFrame({"bacterium": list("ABCD") * 2, "phage": ["P"] * 4 + ["Q"] * 4,
                          "infection": [0, 0, 1, 1, 0, 0, 0, 0]})
    features = pd.DataFrame({"defense": [0, 1, 2, 3, 100]}, index=list("ABCDE"))
    candidates = pd.DataFrame({"bacterium": ["E", "E"], "phage": ["P", "Q"]})
    model, scores = score_known_panel(train, candidates, features, ["defense"], [],
                                      {"C": 0.1, "max_iter": 2000})
    altered, other = score_known_panel(train, candidates.assign(infection=[1, 1]), features,
                                       ["defense"], [], {"C": 0.1, "max_iter": 2000})
    assert model == altered
    assert model["preprocessing"]["means"] == [1.5]
    for method in scores:
        np.testing.assert_array_equal(scores[method], other[method])
    np.testing.assert_allclose(scores["phage_prevalence"], [0.5, 1 / 6])
    with pytest.raises(ValueError, match="disjoint"):
        score_known_panel(train, candidates.assign(bacterium="A"), features, ["defense"], [],
                           {"C": 0.1, "max_iter": 2000})


@pytest.fixture
def local_study(tmp_path, monkeypatch):
    repo, source, output = [tmp_path / name for name in ("repo", "source", "output")]
    source.mkdir()
    hosts, phages = list("ABCDEFGHI"), list("PQR")
    pairs = pd.MultiIndex.from_product([hosts, phages], names=["bacterium", "phage"]).to_frame(index=False)
    pairs["infection"] = [0, 1, 0, 1, 0, 0, 0, 0, 0] * 3
    pairs.to_csv(source / "subset.csv", index=False)
    pairs.pivot(index="bacterium", columns="phage", values="infection").to_csv(source / "mlc.csv", sep=";")
    metadata = pd.DataFrame({"Clermont_Phylo": ["X", "Y", "Z"] * 3, "O-type": "O1"}, index=hosts)
    metadata.to_csv(source / "metadata.csv", sep=";")
    pd.DataFrame({"RM": range(9)}, index=hosts).to_csv(source / "defenses.csv", sep=";")
    plan = json.loads((ROOT / "data/studies/host_retrieval.v1.json").read_text())
    plan.update(expected_candidate_phages=3, categorical_features=["O-type"])
    plan["bootstrap"]["resamples"] = 20
    manifest_path = repo / plan["host_manifest"]
    manifest_path.parent.mkdir(parents=True)
    pd.DataFrame({"bacterium": hosts, "host_group": hosts,
                  "split": ["train"] * 5 + ["validation"] * 2 + ["test"] * 2}).to_csv(manifest_path, index=False)
    registry = json.loads((ROOT / plan["subset_manifest"]).read_text())
    registry["sha256"] = sha256_file(source / "subset.csv")
    registry["expected_summary"] = {k: v for k, v in summarize_interactions(pairs).items()
                                     if k in registry["expected_summary"]}
    registry_path = repo / plan["subset_manifest"]
    registry_path.parent.mkdir(parents=True)
    registry_path.write_text(json.dumps(registry))
    plan["repo_files_sha256"] = {str(p.relative_to(repo)): sha256_file(p)
                                  for p in [manifest_path, registry_path]}
    plan["resources"] = [{"kind": "table", "id": name, "path": f"{name}.csv",
                           "sha256": sha256_file(source / f"{name}.csv"),
                           "size_bytes": (source / f"{name}.csv").stat().st_size}
                          for name in ("subset", "mlc", "metadata", "defenses")]
    path = repo / "data/studies/host_retrieval.v1.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setattr("phage_cocktail.retrieval_study.code_provenance", lambda *args: {"code_base_revision": "fixture"})
    return repo, source, output, plan


def test_study_publishes_fixed_panels_metrics_and_hashes(local_study):
    repo, source, output, _ = local_study
    report = run_retrieval_study(repo, source, output)
    assert report["audit"]["observed_pairs"] == 27
    assert len(report["metrics"]) == 28  # 2x2 partitions + 3 phylogroups, four models
    for name, digest in report["output_sha256"].items():
        assert sha256_file(output / name) == digest
    assert json.loads((output / "report.json").read_text()) == report
    predictions = pd.read_csv(output / "predictions.csv.gz")
    assert set(predictions.groupby(["scenario", "partition", "model", "bacterium"]).size()) == {3}
    for row in report["metrics"]:
        if row["model"] == "random":
            global_row = next(r for r in report["metrics"] if r["scenario"] == row["scenario"]
                              and r["partition"] == row["partition"] and r["model"] == "global_prevalence")
            assert row["all_hosts"] == global_row["all_hosts"]
    with pytest.raises(ValueError, match="must be new"):
        run_retrieval_study(repo, source, output)


def test_changed_manifest_is_rejected_before_fitting(local_study):
    repo, source, output, plan = local_study
    path = repo / plan["host_manifest"]
    path.write_text(path.read_text().replace("train", "other"))
    with pytest.raises(ValueError, match="changed frozen"):
        run_retrieval_study(repo, source, output)
    assert not output.exists()


def test_fit_failure_does_not_publish_partial_run(local_study, monkeypatch):
    repo, source, output, _ = local_study
    def fail(*args, **kwargs):
        raise ValueError("fit failed")
    monkeypatch.setattr("phage_cocktail.retrieval_study.score_known_panel", fail)
    with pytest.raises(ValueError, match="fit failed"):
        run_retrieval_study(repo, source, output)
    assert not output.exists()
    assert not list(output.parent.glob(".output.*"))


def test_similarity_study_keeps_diagnostics_and_old_scores(local_study):
    repo, source, output, plan = local_study
    original = run_retrieval_study(repo, source, output)
    hosts = list("ABCDEFGHI")
    pd.DataFrame(np.abs(np.arange(9)[:, None] - np.arange(9)[None, :]),
                 index=hosts, columns=hosts).to_csv(source / "core.tsv", sep="\t")
    plan["resources"].append({"kind": "table", "id": "core_distances", "path": "core.tsv",
                               "size_bytes": (source / "core.tsv").stat().st_size,
                               "sha256": sha256_file(source / "core.tsv")})
    plan["neighbours"] = {"k": 5}
    (repo / "data/studies/host_similarity.v1.json").write_text(json.dumps(plan))
    out = output.with_name("similarity")
    result = run_retrieval_study(repo, source, out, similarity=True)
    assert [r for r in result["metrics"] if r["model"] != "core_neighbour"] == original["metrics"]
    d = pd.read_csv(out / "host_diagnostics.csv")
    assert len(d) == len(pd.read_csv(out / "per_host_metrics.csv"))
    assert np.isfinite(d.nearest_train_distance).all()
    assert d.loc[d.model.eq("phage_prevalence"), "hit@2_delta_vs_prevalence"].eq(0).all()
    shifts = pd.read_csv(out / "phage_prevalence_shift.csv")
    assert set(shifts.groupby(["scenario", "partition"]).size()) == {3}
