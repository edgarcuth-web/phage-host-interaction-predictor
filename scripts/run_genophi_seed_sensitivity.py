"""Execute both preregistered seeds once, with matching baseline refits and QA."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import time

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.genophi_seed_protocols import REGISTRATION, INPUT_FREEZE, seed_paths, load_study_protocol


def interpreter_path(repo, path):
    # Resolving the venv's Python symlink would select the base interpreter and
    # silently lose the installed numerical environment.
    path = Path(path)
    return path if path.is_absolute() else Path(repo) / path


def execute_studies(root, jobs, *, repo, run=subprocess.run):
    """A failed seed is terminal, but never cancels the other mandatory seed."""
    root = Path(root)
    if root.exists():
        raise ValueError("run-accounting directory exists; do not rerun a started seed")
    root.mkdir(parents=True)
    ledger = {"status": "started", "started_utc": datetime.now(timezone.utc).isoformat(),
              "seeds": {str(s): {"status": "planned", "components": {name: {"status": "planned"}
                         for name, _ in commands}} for s, commands in jobs.items()}}
    path = root / "execution.json"
    atomic_write_json(ledger, path)
    started = time.perf_counter()
    for seed, commands in jobs.items():
        state = ledger["seeds"][str(seed)]
        state["status"] = "started"
        seed_root = root / f"seed{seed}"
        seed_root.mkdir()
        for index, (name, command) in enumerate(commands):
            component = state["components"][name]
            component.update(status="started", command=list(map(str, command)))
            atomic_write_json(ledger, path)
            print(f"Seed {seed}: starting {name}", flush=True)
            tick = time.perf_counter()
            with (seed_root / f"{name}.log").open("x", encoding="utf-8", buffering=1) as log:
                try:
                    result = run(list(map(str, command)), cwd=repo, stdout=log, stderr=subprocess.STDOUT)
                    code = result.returncode
                except OSError as exc:
                    code = -1
                    component["error"] = f"{type(exc).__name__}: {exc}"
            component.update(returncode=code, wall_seconds=time.perf_counter() - tick,
                             status="completed" if code == 0 else "blocked")
            atomic_write_json(ledger, path)
            if code:
                state["status"] = "blocked"
                for later, _ in commands[index + 1:]:
                    state["components"][later]["status"] = "not_started_due_to_prior_block"
                print(f"Seed {seed}: blocked at {name}; retained log and all prior outputs", flush=True)
                break
            print(f"Seed {seed}: completed {name} ({component['wall_seconds']:.1f} s)", flush=True)
        else:
            state["status"] = "verified"
        atomic_write_json(ledger, path)
    ledger.update(status="verified" if all(s["status"] == "verified" for s in ledger["seeds"].values()) else "blocked",
                  wall_seconds=time.perf_counter() - started)
    atomic_write_json(ledger, path)
    return ledger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--baseline-python", type=Path, required=True)
    parser.add_argument("--genophi-python", type=Path, required=True)
    parser.add_argument("--feature-cache", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--proteomes", type=Path, default=Path("data/derived/genophi_proteomes"))
    parser.add_argument("--proteome-receipt", type=Path, default=Path("data/receipts/genophi_proteomes.v1.json"))
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    registration, inputs = read_json_object(repo / REGISTRATION), read_json_object(repo / INPUT_FREEZE)
    if sha256_file(repo / REGISTRATION) != inputs["registration_sha256"]:
        raise ValueError("changed registered decision rule")
    jobs = {}
    for seed in registration["sensitivity_seeds"]:
        paths = seed_paths(seed)
        for kind in ("genome", "comparator"):
            load_study_protocol(repo, paths[kind], kind)
        seed_root = args.outdir.resolve() / f"seed{seed}"
        baseline_root, comparator_root = seed_root / "baseline", seed_root / "comparator"
        shared = ["--data-root", args.data_root.resolve()]
        jobs[seed] = [
            ("baseline_fit", [interpreter_path(repo, args.baseline_python), repo / "scripts/run_genome_host_retrieval.py",
                "--protocol", paths["genome"], "--scenario", "grouped_host", "--outdir", baseline_root, *shared]),
            ("baseline_qa", [interpreter_path(repo, args.baseline_python), repo / "scripts/verify_genophi_comparator.py",
                "--baseline-protocol", paths["genome"], "--run", baseline_root,
                "--out", seed_root / "BASELINE_QA.json", "--publish", repo / paths["baseline_output"], *shared]),
            ("comparator_fit", [interpreter_path(repo, args.genophi_python), repo / "scripts/run_genophi_comparator.py",
                "--protocol", paths["comparator"], "--feature-cache", args.feature_cache.resolve(),
                "--proteomes", args.proteomes.resolve(), "--proteome-receipt", args.proteome_receipt.resolve(),
                "--outdir", comparator_root, *shared]),
            ("comparator_qa", [interpreter_path(repo, args.genophi_python), repo / "scripts/verify_genophi_comparator.py",
                "--run", comparator_root, "--proteomes", args.proteomes.resolve(),
                "--proteome-receipt", args.proteome_receipt.resolve(), "--out", seed_root / "COMPARATOR_QA.json",
                "--publish", repo / paths["comparator_output"], *shared]),
        ]
    result = execute_studies(args.outdir.resolve(), jobs, repo=repo)
    print({s: record["status"] for s, record in result["seeds"].items()})
    if result["status"] != "verified":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
