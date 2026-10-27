"""Replay six committed real vectors through the executor with zero API calls."""
import argparse
from pathlib import Path

import numpy as np

from phage_cocktail.artifacts import atomic_write_json, code_provenance, sha256_file
from phage_cocktail.extraction_store import LocalStore
from phage_cocktail.nvidia_host_extraction import make_plan, run_extraction, verify_extraction


def verify(repo, data_root, output):
    repo, output = Path(repo), Path(output)
    if output.exists():
        raise ValueError("replay output must be a new directory")
    output.mkdir(parents=True)
    cache = repo / "data/studies/results/nvidia_evo2_host_pilot.v1"
    plan = make_plan(repo, data_root)
    atomic_write_json(plan, output / "plan.json")
    store = LocalStore(output / "objects")
    calls = []
    def no_network(*args, **kwargs):
        calls.append(1)
        raise AssertionError("offline replay attempted a request/key access")
    first = run_extraction(repo, data_root, plan, store, pilot_cache=cache, max_jobs=2,
                           key_provider=no_network, call_forward=no_network)
    if first["completed_jobs"] != 2 or first["new_post_attempts"]:
        raise ValueError("replay pause boundary differs")
    complete = run_extraction(repo, data_root, plan, store, pilot_cache=cache,
                              key_provider=no_network, call_forward=no_network)
    if complete["completed_jobs"] != 6 or complete["status"] != "complete_engineering_extraction":
        raise ValueError("replay did not complete the exact six jobs")
    matrix_path = store.root / "host_vectors.npy"
    np.testing.assert_array_equal(np.load(matrix_path, allow_pickle=False),
                                  np.load(cache / "host_vectors.npy", allow_pickle=False))
    def hashes():
        return {p.relative_to(store.root).as_posix(): sha256_file(p)
                for p in sorted(store.root.rglob("*")) if p.is_file()}
    before = hashes()
    repeated = run_extraction(repo, data_root, plan, store, key_provider=no_network, call_forward=no_network)
    if before != hashes() or repeated != complete or verify_extraction(plan, store) != complete or calls:
        raise ValueError("completed replay changed artifacts or accessed a provider")
    return {"status": "passed", "scope": "offline executor replay of previously measured real pilot",
            "new_nvidia_calls": len(calls), "live_s3_tested": False,
            "paused_after_completed_jobs": first["completed_jobs"], "resumed_completed_jobs": complete["completed_jobs"],
            "three_host_vectors_match_committed_pilot_exactly": True,
            "completed_resume_preserves_every_object_hash": True,
            "budget": plan["budget"], "final_state": complete,
            "stored_object_bytes": sum(p.stat().st_size for p in store.root.rglob("*") if p.is_file()),
            "object_sha256": before, "plan_file_sha256": sha256_file(output / "plan.json"),
            "source_pilot_report_sha256": sha256_file(cache / "report.json"),
            "provenance": code_provenance(repo, repo / "scripts/verify_nvidia_host_extraction.py"),
            "prediction_metrics_computed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--qa-out", type=Path, required=True)
    args = parser.parse_args()
    report = verify(Path(__file__).resolve().parents[1], args.data_root, args.outdir)
    atomic_write_json(report, args.qa_out)
    print({k: report[k] for k in ("status", "new_nvidia_calls", "live_s3_tested", "stored_object_bytes")})


if __name__ == "__main__":
    main()
