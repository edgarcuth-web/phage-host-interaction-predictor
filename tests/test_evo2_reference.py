import csv
import json
from pathlib import Path
import subprocess

import pytest

from phage_cocktail import evo2_reference as reference
from phage_cocktail.artifacts import sha256_file


@pytest.fixture
def protocol():
    path = Path(__file__).resolve().parents[1] / "data/studies/evo2_retrieval_reference.v1.json"
    return json.loads(path.read_text())


def write_outputs(output, protocol):
    for method in protocol["matrices"]:
        row = {"evaluated_phages": 4490, **protocol["expected_metrics"][method]}
        (output / f"{method}.json").write_text(json.dumps(row))
    row = {"evaluated_phages": 4490, **protocol["expected_metrics"]["rrf"]}
    with (output / "rrf.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def test_prespecified_commands_do_not_search_fusion_combinations(tmp_path, protocol):
    commands = dict(reference.build_commands(tmp_path / "source", tmp_path / "out", protocol))
    assert set(commands) == {"blastn", "virhostmatcher", "phist", "wish", "evo2", "rrf"}
    fusion = commands["rrf"]
    assert fusion[fusion.index("--sizes") + 1] == "4"
    assert fusion[fusion.index("--fusion-rrf-k") + 1] == "60"
    assert fusion.count("--matrix") == 4
    assert not any(value.startswith("wish=") for value in fusion)
    assert "--matrix-orientation" not in commands["evo2"]


def test_summary_accepts_rounding_and_detects_missing_coverage(tmp_path, protocol):
    write_outputs(tmp_path, protocol)
    path = tmp_path / "evo2.json"
    row = json.loads(path.read_text())
    row["hit@10"] = 0.5540521220512888
    path.write_text(json.dumps(row))
    assert all(r["matches_reference"] for r in reference.summarize_outputs(tmp_path, protocol).values())
    row["evaluated_phages"] -= 1
    path.write_text(json.dumps(row))
    assert not reference.summarize_outputs(tmp_path, protocol)["evo2"]["matches_reference"]


@pytest.mark.parametrize("value", [0.5, -0.1, 1.1, None, "nan", "inf"])
def test_changed_or_invalid_metrics_cannot_pass(tmp_path, protocol, value):
    write_outputs(tmp_path, protocol)
    path = tmp_path / "evo2.json"
    row = json.loads(path.read_text())
    row["mrr"] = value
    path.write_text(json.dumps(row))
    try:
        result = reference.summarize_outputs(tmp_path, protocol)
    except (ValueError, TypeError):
        return
    assert not result["evo2"]["matches_reference"]


def test_multiple_fusion_rows_rejected(tmp_path, protocol):
    write_outputs(tmp_path, protocol)
    path = tmp_path / "rrf.csv"
    lines = path.read_text().splitlines()
    path.write_text("\n".join([*lines, lines[-1]]) + "\n")
    with pytest.raises(ValueError, match="exactly one"):
        reference.summarize_outputs(tmp_path, protocol)


@pytest.fixture
def local_run(tmp_path, protocol, monkeypatch):
    repo, source, output = (tmp_path / name for name in ("repo", "source", "output"))
    source.mkdir()
    resource = source / "locked.txt"
    resource.write_text("original")
    protocol["resources"] = [{"kind": "test", "id": "fixture", "path": "locked.txt",
                              "size_bytes": 8, "sha256": sha256_file(resource)}]
    path = repo / "data/studies/evo2_retrieval_reference.v1.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(protocol))
    monkeypatch.setattr(reference, "code_provenance", lambda *args: {"code_base_revision": "test"})
    return repo, source, output, protocol


def test_same_size_tampering_rejected_before_execution(local_run, monkeypatch):
    repo, source, output, _ = local_run
    (source / "locked.txt").write_text("tampered")
    monkeypatch.setattr(reference.subprocess, "run", lambda *a, **k: pytest.fail("executed changed input"))
    with pytest.raises(ValueError, match="changed locked resource"):
        reference.run_reference(repo, source, output)
    assert not output.exists()


def test_existing_outputs_rejected(local_run):
    repo, source, output, _ = local_run
    output.mkdir()
    with pytest.raises(ValueError, match="must be new"):
        reference.run_reference(repo, source, output)


def test_run_publishes_hashed_outputs(local_run, monkeypatch):
    repo, source, output, protocol = local_run
    def fake_run(command, **kwargs):
        assert kwargs["check"] and kwargs["cwd"] == source
        flag = "--output" if "--output" in command else "--output-json"
        write_outputs(Path(command[command.index(flag) + 1]).parent, protocol)
    monkeypatch.setattr(reference.subprocess, "run", fake_run)
    report = reference.run_reference(repo, source, output)
    assert report["matches_reference"]
    assert not report["generated_embeddings"]
    assert len(report["executions"]) == 6
    assert json.loads((output / "report.json").read_text()) == report
    for name, digest in report["output_sha256"].items():
        assert sha256_file(output / name) == digest


def test_subprocess_failure_does_not_publish_partial_result(local_run, monkeypatch):
    repo, source, output, _ = local_run
    def fail(command, **kwargs):
        raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(reference.subprocess, "run", fail)
    with pytest.raises(subprocess.CalledProcessError):
        reference.run_reference(repo, source, output)
    assert not output.exists()
    assert not list(output.parent.glob(".output.*"))
