"""Leakage, memory-reduction equivalence, fixed panels and resumable receipts."""
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.genophi_comparator import (
    ComparatorBlocked, _fit_repeat, assemble_feature_table, candidate_features, checkpoint,
    evaluate_score_runs, predict_fitted_model, select_best_hits_streaming,
    selected_occurrences, training_table,
)


def test_chunked_best_hits_preserve_float32_keys_and_first_tie(tmp_path):
    # q1 spans chunks, and its final row improves the e-value. q2's identity
    # rounds to a float32 tie, so its first hit (not the later float64 maximum)
    # must survive. Distinct targets map to different families.
    rows = [
        ("q1", "a", 0.8, 1e-20), ("q2", "b", 0.800000001, 1e-10),
        ("q1", "b", 0.9, 1e-10), ("q2", "c", 0.800000002, 1e-10),
        ("q3", "c", 1.0, 1e-100), ("q1", "c", 0.9, 1e-15),
        ("q3", "a", 1.0, 1e-110),
    ]
    assigned, clusters = tmp_path / "assigned.tsv", tmp_path / "clusters.tsv"
    assigned.write_text("".join(f"{q}\t{t}\t100\t{i}\t{e}\t1\t90\t90\t1\t90\t90\n"
                                for q, t, i, e in rows))
    clusters.write_text("F1\ta\nF2\tb\nF3\tc\n")
    expected = b"q1\tF3\nq2\tF2\nq3\tF3\n"
    for chunk in (1, 2, 4, 100):
        out = tmp_path / f"best{chunk}.tsv"
        select_best_hits_streaming(assigned, out, clusters, chunksize=chunk)
        assert out.read_bytes() == expected


def test_chunked_reducer_matches_installed_upstream_on_randomized_alignments(tmp_path):
    upstream = pytest.importorskip("genophi.mmseqs2_clustering")
    rng = np.random.default_rng(42)
    rows = []
    for _ in range(1000):
        rows.append([f"Q{rng.integers(25)}", f"T{rng.integers(10)}", 100,
                     rng.choice([0.4, 0.9, 1]), rng.choice([1e-5, 1e-20, 1e-80]), *([1] * 6)])
    assigned, clusters = tmp_path / "assigned.tsv", tmp_path / "clusters.tsv"
    pd.DataFrame(rows).to_csv(assigned, sep="\t", header=False, index=False)
    clusters.write_text("".join(f"C{i // 2}\tT{i}\n" for i in range(10)))
    upstream.select_best_hits(str(assigned), str(tmp_path / "upstream.tsv"), str(clusters))
    select_best_hits_streaming(assigned, tmp_path / "bounded.tsv", clusters, chunksize=17)
    assert (tmp_path / "upstream.tsv").read_bytes() == (tmp_path / "bounded.tsv").read_bytes()


@pytest.fixture
def feature_inputs():
    phenotype = pd.DataFrame([
        (h, p, int((h == "A") == (p == "P"))) for h in "ABC" for p in ["P", "Q"]
    ], columns=["strain", "phage", "interaction"])
    hosts = pd.DataFrame({"strain": list("ABC"), "sc_0": [0, 1, 1]})
    phages = pd.DataFrame({"phage": ["P", "Q"], "pc_0": [0, 1]})
    split = pd.DataFrame({"bacterium": list("ABC"), "host_group": list("ABC"),
                          "split": ["train", "validation", "test"]})
    return phenotype, hosts, phages, split


def test_outer_labels_cannot_change_training_or_prediction_features(feature_inputs):
    phenotype, hosts, phages, split = feature_inputs
    full = assemble_feature_table(phenotype, hosts, phages)
    train = training_table(full, split)
    assert set(train.strain) == {"A"}
    changed = phenotype.copy()
    changed.loc[changed.strain.ne("A"), "interaction"] ^= 1
    pd.testing.assert_frame_equal(train, training_table(assemble_feature_table(changed, hosts, phages), split))
    candidates = candidate_features(hosts, phages, ["B", "C"], ["P", "Q"])
    assert "interaction" not in candidates
    assert len(candidates) == 4
    with pytest.raises(ValueError, match="groups cross"):
        training_table(full, split.assign(host_group="same"))
    with pytest.raises(ValueError, match="only identifiers"):
        training_table(full.assign(source_mlc=4), split)


def test_feature_join_rejects_silent_host_loss_and_outcome_metadata(feature_inputs):
    phenotype, hosts, phages, _ = feature_inputs
    with pytest.raises(ValueError, match="genome sets"):
        assemble_feature_table(phenotype, hosts.iloc[:2], phages)
    with pytest.raises(ValueError, match="genome sets"):
        assemble_feature_table(phenotype, hosts.assign(n_infections=99), phages)
    with pytest.raises(ValueError, match="binary"):
        assemble_feature_table(phenotype, hosts.assign(sc_0=np.nan), phages)
    with pytest.raises(ValueError, match="binary"):
        candidate_features(hosts, phages, ["MISSING"], ["P", "Q"])


def test_prediction_never_passes_ids_or_labels_and_rejects_partial_output(feature_inputs):
    _, hosts, phages, _ = feature_inputs
    candidates = candidate_features(hosts, phages, ["B", "C"], ["P", "Q"])

    class Model:
        feature_names_ = ["pc_0", "sc_0"]
        classes_ = [0, 1]

        def predict_proba(self, frame, **kwargs):
            assert list(frame) == self.feature_names_
            return np.tile([0.25, 0.75], (len(frame), 1))

    scores = predict_fitted_model(Model(), candidates, threads=1)
    assert list(scores) == ["bacterium", "phage", "score"]
    assert scores.score.eq(0.75).all()
    with pytest.raises(ValueError, match="held-out labels|non-protein|only identifiers"):
        predict_fitted_model(Model(), candidates.assign(interaction=0), threads=1)
    bad = SimpleNamespace(feature_names_=["sc_0"], classes_=[0, 1],
                          predict_proba=lambda *a, **k: np.full((3, 2), 0.5))
    with pytest.raises(ComparatorBlocked, match="incomplete"):
        predict_fitted_model(bad, candidates, threads=1)


def test_occurrence_is_counted_once_per_completed_run_including_resumes(tmp_path):
    for i, features in enumerate([["sc_0", "pc_0"], ["sc_0", "pc_1"], ["sc_0", "pc_0"]]):
        directory = tmp_path / f"run_{i}"
        directory.mkdir()
        pd.DataFrame({"Feature": features, "Importance": 0.5}).to_csv(directory / "feature_importances.csv", index=False)
    result = selected_occurrences(tmp_path, num_runs=3, allowed_features={"sc_0", "pc_0", "pc_1"})
    assert result.set_index("Feature").Occurrence.to_dict() == {"sc_0": 3, "pc_0": 2, "pc_1": 1}
    (tmp_path / "run_1/feature_importances.csv").write_text("Feature\nsc_0\nsc_0\n")
    with pytest.raises(ComparatorBlocked, match="selected feature"):
        selected_occurrences(tmp_path, num_runs=3, allowed_features={"sc_0", "pc_0", "pc_1"})


def test_checkpoint_reuse_reverifies_bytes_and_failed_stage_is_not_published(tmp_path):
    output = tmp_path / "payload.txt"

    def action():
        output.write_text("checkpointed output")
        return [output], {"rows": 1}

    first = checkpoint(tmp_path, "fixture", action)
    assert checkpoint(tmp_path, "fixture", lambda: pytest.fail("completed stage reran")) == first
    output.write_text("tampered")
    with pytest.raises(ValueError, match="changed checkpoint"):
        checkpoint(tmp_path, "fixture", action)

    def fail():
        raise ComparatorBlocked("no supported feature stratum")

    with pytest.raises(ComparatorBlocked):
        checkpoint(tmp_path, "blocked", fail)
    assert not (tmp_path / "checkpoints/blocked.json").exists()
    failure = json.loads((tmp_path / "logs/blocked.00.failure.json").read_text())
    assert failure["status"] == "blocked"


def test_retrieval_uses_full_panel_unknown_host_exclusion_and_exact_ties():
    # T has a positive, N is all negative, U has an unknown second label.
    observed = pd.DataFrame([
        ("A", "P", 1), ("A", "Q", 0), ("T", "P", 1), ("T", "Q", 0),
        ("N", "P", 0), ("N", "Q", 0), ("U", "P", 1),
    ], columns=["bacterium", "phage", "infection"])
    splits = pd.DataFrame({"bacterium": ["A", "T", "N", "U"], "host_group": ["A", "T", "N", "U"],
                           "split": ["train", "test", "test", "test"], "scenario": "s"})
    panel = pd.MultiIndex.from_product([["T", "N", "U"], ["P", "Q"]], names=["bacterium", "phage"]).to_frame(index=False).assign(score=0.5)
    plan = {"scenarios": ["s"], "repeats": 2, "modeling": {"num_runs_modeling": 2},
            "evaluation": {"partitions": ["test"], "k": [1, 2],
                           "bootstrap": {"resamples": 10, "seed": 42, "confidence": 0.95}}}
    runs = {("s", r, i): panel.copy() for r in range(2) for i in range(2)}
    summaries, per_host, individual, _ = evaluate_score_runs(observed, splits, runs, plan)
    result = next(s for s in summaries if s["model"] == "genophi_median")
    assert result["all_hosts"]["hit@1"] == 0.25
    assert result["eligible_hosts_only"]["hit@1"] == 0.5
    assert result["all_hosts"]["hit@2"] == 0.5
    assert result["incomplete_hosts_excluded"] == ["U"]
    assert result["all_negative_hosts"] == 1
    assert len(individual) == 4
    assert set(per_host.bacterium) == {"T", "N"}
    runs[("s", 0, 0)] = panel.iloc[:-1]
    with pytest.raises(ValueError, match="full fixed"):
        evaluate_score_runs(observed, splits, runs, plan)


def test_protocol_cutoff_and_supervised_budget_are_frozen():
    repo = Path(__file__).resolve().parents[1]
    plan = json.loads((repo / "data/studies/genophi_comparator.v1.json").read_text())
    fs = plan["feature_selection"]
    assert fs["cutoff"] == 3 <= fs["num_runs_fs"]
    assert fs["min_features"] == 20 and fs["max_features"] == 501
    assert fs["num_features"] * fs["num_runs_fs"] / fs["cutoff"] < fs["max_features"]


def test_upstream_fit_apis_receive_only_outer_train_and_predictions_have_full_panel(
        tmp_path, monkeypatch, feature_inputs):
    import joblib

    phenotype, hosts, phages, split = feature_inputs
    for folder in ("cohort", "clustering/strain/features", "clustering/phage/features"):
        (tmp_path / folder).mkdir(parents=True)
    split.assign(scenario="s").to_csv(tmp_path / "cohort/scenario_splits.csv", index=False)
    hosts.to_csv(tmp_path / "clustering/strain/features/feature_table.csv", index=False)
    phages.to_csv(tmp_path / "clustering/phage/features/feature_table.csv", index=False)
    assemble_feature_table(phenotype, hosts, phages).to_csv(tmp_path / "full_feature_table.csv", index=False)
    seen = []

    def assert_training(path):
        frame = pd.read_csv(path)
        assert set(frame.strain) == {"A"}
        seen.append(frame.copy())
        return frame

    def select(input_path, output, **kwargs):
        assert_training(input_path)
        directory = Path(output) / "run_0"
        directory.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"Feature": ["sc_0", "pc_0"], "Importance": [0.5, 0.5]}).to_csv(
            directory / "feature_importances.csv", index=False)

    def generate(selection, input_path, output, **kwargs):
        assert kwargs["cut_offs"] == [1]
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        assert_training(input_path).to_csv(directory / "select_feature_table_cutoff_1.csv", index=False)

    def fit(input_path, output, *args, **kwargs):
        frame = assert_training(input_path)
        directory = Path(output)
        (directory / "best_model.pkl").write_bytes(b"fake test model")
        frame[["strain", "phage"]].to_csv(directory / "best_model_predictions.csv", index=False)

    feature_module, model_module = ModuleType("genophi.feature_selection"), ModuleType("genophi.select_feature_modeling")
    feature_module.run_feature_selection_iterations = select
    feature_module.generate_feature_tables = generate
    model_module.model_testing_select_MCC = fit
    monkeypatch.setitem(sys.modules, "genophi", ModuleType("genophi"))
    monkeypatch.setitem(sys.modules, "genophi.feature_selection", feature_module)
    monkeypatch.setitem(sys.modules, "genophi.select_feature_modeling", model_module)

    def predict(frame, **kwargs):
        assert set(frame) == {"sc_0", "pc_0"}
        assert len(frame) == 4
        return np.tile([0.3, 0.7], (len(frame), 1))

    model = SimpleNamespace(feature_names_=["sc_0", "pc_0"], classes_=[0, 1], predict_proba=predict)
    monkeypatch.setattr(joblib, "load", lambda path: model)
    plan = {"feature_selection": {"n_clusters": 2, "num_features": 2, "num_runs_fs": 1,
                                  "method": "rfe", "cutoff": 1, "min_features": 1, "max_features": 3},
            "compute": {"threads": 1, "max_ram_gb": 1}, "modeling": {"num_runs_modeling": 1}}
    output = _fit_repeat(tmp_path, "s", 0, plan)
    predictions = pd.read_csv(output / "models/run_0/candidate_scores.csv")
    assert set(predictions.bacterium) == {"B", "C"}
    assert len(predictions) == 4 and len(seen) == 3
    # A resumed invocation verifies and reuses all supervised stages.
    _fit_repeat(tmp_path, "s", 0, plan)
    assert len(seen) == 3
