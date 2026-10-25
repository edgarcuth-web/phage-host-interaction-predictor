"""Run the pinned upstream species-level retrieval evaluation unchanged."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd

from .artifacts import (
    atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory,
)
from .host_study import verify_resources

METRICS = ("mrr", "hit@1", "hit@5", "hit@10")


def build_commands(source_root, output_dir, protocol):
    """Invoke upstream CLIs; do not reimplement their ranking or weighting."""
    root, output = Path(source_root).resolve(), Path(output_dir).resolve()
    mapping = str(root / protocol["mapping_path"])
    commands = []
    for method, matrix in protocol["matrices"].items():
        commands.append((method, [
            sys.executable, str(root / "code/compute_metrics.py"),
            "--matrix", str(root / matrix), "--mapping", mapping,
            "--output-json", str(output / f"{method}.json"),
        ]))
    fusion = protocol["fusion"]
    command = [sys.executable, str(root / "code/evaluate_rank_fusion_overall.py")]
    for method in fusion["methods"]:
        command.extend(["--matrix", f"{method}={root / protocol['matrices'][method]}"])
    command.extend([
        "--mapping", mapping, "--sizes", str(len(fusion["methods"])),
        "--fusion-method", "rrf", "--fusion-rrf-k", str(fusion["rrf_k"]),
        "--phage-policy", "intersection", "--output", str(output / "rrf.csv"),
    ])
    commands.append(("rrf", command))
    return commands


def summarize_outputs(output_dir, protocol):
    """Compare full-precision outputs with rounded upstream README targets."""
    output = Path(output_dir)
    observed = {method: read_json_object(output / f"{method}.json")
                for method in protocol["matrices"]}
    with (output / "rrf.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError("expected exactly one prespecified fusion result")
    observed["rrf"] = rows[0]
    comparisons = {}
    for method, expected in protocol["expected_metrics"].items():
        row = observed[method]
        count = float(row["evaluated_phages"])
        values = {metric: float(row[metric]) for metric in METRICS}
        if not count.is_integer() or count < 0 or any(
            not math.isfinite(value) or not 0 <= value <= 1 for value in values.values()
        ):
            raise ValueError(f"invalid metric output for {method}")
        deltas = {key: values[key] - expected[key] for key in METRICS}
        comparisons[method] = {
            "evaluated_phages": int(count), "observed": values,
            "expected": expected, "delta": deltas,
            "matches_reference": (
                count == protocol["expected_evaluated_phages"]
                and all(abs(delta) <= protocol["absolute_tolerance"]
                        for delta in deltas.values())
            ),
        }
    return comparisons


def run_reference(repo, source_root, output_dir, *, fetch=False):
    """Verify locked inputs, run six evaluations, then publish a hashed report."""
    repo = Path(repo).resolve()
    root, output = Path(source_root).resolve(), Path(output_dir).resolve()
    if output.exists():
        raise ValueError("reference output directory must be new")
    protocol_path = repo / "data/studies/evo2_retrieval_reference.v1.json"
    protocol = read_json_object(protocol_path)
    verify_resources(protocol, root, fetch=fetch)
    with staged_directory(output) as stage:
        executions = []
        for name, command in build_commands(root, stage, protocol):
            with (stage / f"{name}.log").open("w", encoding="utf-8") as log:
                subprocess.run(command, check=True, cwd=root, stdout=log, stderr=subprocess.STDOUT)
            executions.append({"name": name, "argv": command})
        comparisons = summarize_outputs(stage, protocol)
        report = {
            "study_id": protocol["study_id"],
            "task": protocol["task"],
            "source_revision": protocol["source_revision"],
            "protocol_sha256": sha256_file(protocol_path),
            **code_provenance(repo, repo / "scripts/run_evo2_reference.py"),
            "versions": {"python": platform.python_version(),
                         "numpy": np.__version__, "pandas": pd.__version__},
            "input_sha256": {entry["path"]: entry["sha256"]
                             for entry in protocol["resources"]},
            "executions": executions,
            "comparisons": comparisons,
            "matches_reference": all(row["matches_reference"] for row in comparisons.values()),
            "generated_embeddings": False,
            "trained_model": False,
            "limitations": protocol["limitations"],
            "output_sha256": {path.name: sha256_file(path) for path in sorted(stage.iterdir())},
        }
        atomic_write_json(report, stage / "report.json")
    return report
