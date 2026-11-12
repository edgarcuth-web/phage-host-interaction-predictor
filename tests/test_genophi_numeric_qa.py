"""The independent GenoPHI QA must catch score, support and artifact drift."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/verify_genophi_comparator.py"
SPEC = importlib.util.spec_from_file_location("verify_genophi_comparator", SCRIPT)
qa = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(qa)


def test_qa_reloads_pair_alignment_and_requires_bit_exact_inference():
    candidates = pd.DataFrame({"strain": ["H", "H"], "phage": ["P", "Q"], "sc_0": [1, 1], "pc_0": [0, 1]})

    def predict(frame, **kwargs):
        assert list(frame) == ["pc_0", "sc_0"]
        return np.array([[0.75, 0.25], [0.25, 0.75]])

    model = SimpleNamespace(feature_names_=["pc_0", "sc_0"], classes_=[0, 1], predict_proba=predict)
    saved = pd.DataFrame({"bacterium": ["H", "H"], "phage": ["Q", "P"], "score": [0.75, 0.25]})
    np.testing.assert_array_equal(qa.verify_predicted_scores(model, candidates, saved, threads=1), [0.25, 0.75])
    saved.loc[0, "score"] = np.nextafter(0.75, 1.)
    with pytest.raises(ValueError, match="exactly"):
        qa.verify_predicted_scores(model, candidates, saved, threads=1)
    with pytest.raises(ValueError, match="complete candidate"):
        qa.verify_predicted_scores(model, candidates, saved.iloc[:1], threads=1)


def test_qa_rejects_altered_artifacts_and_changed_support(tmp_path):
    payload = b"fixed run artifact"
    path = tmp_path / "output.txt"
    path.write_bytes(payload)
    receipt = {"status": "complete", "output_sha256": {"output.txt": sha256(payload).hexdigest()}}
    assert qa.verify_artifact_hashes(tmp_path, [receipt]) == 1
    path.write_bytes(b"modified")
    with pytest.raises(ValueError, match="changed"):
        qa.verify_artifact_hashes(tmp_path, [receipt])
    with pytest.raises(ValueError, match="support"):
        qa.verify_summary({"complete_panel_hosts": 34}, {"complete_panel_hosts": 33})
    with pytest.raises(ValueError, match="metric differs"):
        qa.verify_summary({"all_hosts": {"hit@2": 0.5}}, {"all_hosts": {"hit@2": 0.75}})


def test_upstream_singleton_filter_requires_an_explicit_genome_selection(tmp_path):
    upstream = pytest.importorskip("genophi.mmseqs2_clustering")
    source = tmp_path / "presence.csv"
    pd.DataFrame({"Genome": list("ABC"), "family_common": [1, 1, 1],
                  "family_singleton": [0, 0, 1]}).to_csv(source, index=False)
    upstream.run_feature_assignment(str(source), str(tmp_path / "default"))
    default = pd.read_csv(tmp_path / "default/feature_table.csv").drop(columns="strain")
    assert 1 in default.sum().tolist()
    genomes = tmp_path / "genomes.csv"
    pd.DataFrame({"strain": list("ABC")}).to_csv(genomes, index=False)
    upstream.run_feature_assignment(str(source), str(tmp_path / "filtered"), select=str(genomes))
    filtered = pd.read_csv(tmp_path / "filtered/feature_table.csv").drop(columns="strain")
    assert filtered.sum().tolist() == [3]


def test_publication_retains_verified_results_and_excludes_genomes_and_models(tmp_path):
    root, published = tmp_path / "run", tmp_path / "published"
    (root / "evaluation").mkdir(parents=True)
    (root / "checkpoints").mkdir()
    (root / "cohort").mkdir()
    metric = b"model,hit@2\nfixture,0.5\n"
    (root / "evaluation/metrics.csv").write_bytes(metric)
    report = {"output_sha256": {"metrics.csv": sha256(metric).hexdigest()}}
    (root / "evaluation/report.json").write_text(json.dumps(report))
    for name in ("provenance.json", "protocol.json"):
        (root / name).write_text("{}")
    for name in ("host_manifest.csv", "scenario_splits.csv"):
        (root / "cohort" / name).write_text("bacterium\nH\n")
    (root / "evaluation/DO_NOT_PUBLISH.faa").write_text(">H_1\nMAA\n")
    (root / "best_model.pkl").write_bytes(b"model test fixture")
    runtime = tmp_path / "runtime.json"
    runtime.write_text(json.dumps({"status": "completed", "elapsed_wall_seconds": 100}))
    receipt = {"status": "verified", "report_sha256": sha256((root / "evaluation/report.json").read_bytes()).hexdigest()}
    result = qa.publish_verified_results(root, published, receipt, runtime_receipt=runtime)
    assert (published / "metrics.csv").read_bytes() == metric
    assert result["published_output_sha256"]["metrics.csv"] == sha256(metric).hexdigest()
    assert json.loads((published / "runtime.json").read_text())["elapsed_wall_seconds"] == 100
    assert not list(published.rglob("*.faa")) and not list(published.rglob("*.pkl"))
    with pytest.raises(ValueError, match="already exists"):
        qa.publish_verified_results(root, published, receipt)


def test_committed_result_bundle_matches_qa_hashes_and_keeps_full_candidate_support():
    repo = Path(__file__).resolve().parents[1]
    root = repo / "data/studies/results/genophi_comparator.v1"
    receipt = json.loads((root / "NUMERIC_QA.json").read_text())
    report = json.loads((root / "report.json").read_text())
    assert receipt["status"] == "verified" and receipt["models_reloaded"] == 60
    for relative, digest in receipt["published_output_sha256"].items():
        assert sha256((root / relative).read_bytes()).hexdigest() == digest, relative
    assert sha256((repo / "scripts/verify_genophi_comparator.py").read_bytes()).hexdigest() == receipt["verification_script_sha256"]
    per_host = pd.read_csv(root / "per_host_metrics.csv")
    for summary in report["metrics"]:
        part = per_host[per_host.scenario.eq(summary["scenario"]) & per_host.partition.eq(summary["partition"])
                        & per_host.model.eq(summary["model"])]
        assert len(part) == summary["complete_panel_hosts"]
        assert part.eligible.sum() == summary["eligible_hosts"]
        assert not set(part.bacterium) & set(summary["incomplete_hosts_excluded"])
        for metric, expected in summary["all_hosts"].items():
            assert part[metric].mean() == pytest.approx(expected, abs=1e-14, rel=0)
    for scenario in report["protocol"]["scenarios"]:
        scores = pd.read_csv(root / f"{scenario.replace(':', '_')}_candidate_scores.csv")
        hosts = scores.bacterium.nunique()
        assert scores.phage.nunique() == 94
        assert len(scores) == hosts * 94
        assert not scores.duplicated(["bacterium", "phage"]).any()
        assert np.isfinite(scores.score).all()
        for summary in report["metrics"]:
            if summary["scenario"] == scenario:
                assert set(summary["incomplete_hosts_excluded"]).issubset(set(scores.bacterium))
