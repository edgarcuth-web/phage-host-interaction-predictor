"""Offline numeric verification of a completed hosted engineering pilot."""
import argparse
from pathlib import Path

import numpy as np

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file


def verify(repo, cache):
    report = read_json_object(cache / "report.json")
    plan = read_json_object(cache / "protocol.json")
    if report["status"] != "complete_engineering_pilot" or report["prediction_metrics_computed"]:
        raise ValueError("requires completed engineering-only pilot")
    if sha256_file(repo / "data/studies/nvidia_evo2_host_pilot.v1.json") != report["protocol_sha256"]:
        raise ValueError("pilot protocol changed")
    for name, digest in report["output_sha256"].items():
        if Path(name).name != name or sha256_file(cache / name) != digest:
            raise ValueError("pilot payload hash mismatch")
    for name, digest in report["provenance"]["code_sha256"].items():
        if sha256_file(repo / name) != digest:
            raise ValueError("pilot source code hash mismatch")
    if report["hosts"] != plan["hosts"] or len(report["jobs"]) != 6:
        raise ValueError("pilot support differs")
    vectors, new_responses, records = [], [], []
    for i, row in enumerate(report["jobs"]):
        vector_path = cache / f"job_{i:02d}.npy"
        vector = np.load(vector_path, allow_pickle=False)
        if sha256_file(vector_path) != row["pooled_sha256"] or vector.shape != (plan["dimension"],):
            raise ValueError("pilot pooled vector identity differs")
        if vector.dtype != np.float32 or not np.isfinite(vector).all() or not np.any(vector):
            raise ValueError("pilot pooled vector invalid")
        response = row["response"]
        if (response["shape"] != [plan["window_bp"], plan["dimension"]] or not response["finite"]
                or response["http_status"] != 200 or response["input_sha256"] != row["input"]["input_sha256"]):
            raise ValueError("pilot response identity differs")
        if row["origin"] == "new_hosted_forward":
            new_responses.append(response["response_bytes"])
        vectors.append(vector)
        records.append({"host": row["input"]["host"], "strand": row["input"]["strand"],
                        "sampled_assembly_fraction": row["input"]["sampled_assembly_fraction"],
                        "pooled_l2": float(np.linalg.norm(vector.astype(np.float64)))})
    reconstructed = np.stack([np.mean(vectors[i:i+2], axis=0, dtype=np.float64).astype(np.float32)
                              for i in (0,2,4)])
    stored = np.load(cache / "host_vectors.npy", allow_pickle=False)
    np.testing.assert_array_equal(stored, reconstructed)
    if len(new_responses) != report["new_post_attempts"] or sum(new_responses) != report["new_response_bytes"]:
        raise ValueError("pilot budget accounting differs")
    if (sum(new_responses) > plan["request_budget"]["max_new_response_bytes"]
            or len(new_responses) > plan["request_budget"]["max_new_inference_posts"]):
        raise ValueError("pilot exceeds frozen budget")
    parent = read_json_object(repo / plan["feasibility_report"])
    if report["jobs"][0]["pooled_sha256"] != parent["output_sha256"]["pooled_mean.npy"]:
        raise ValueError("reused feasibility vector differs")
    return {"status": "passed", "six_oriented_vectors_verified": True,
            "three_host_vectors_reconstructed_exactly": True, "new_post_attempts": len(new_responses),
            "new_response_bytes": sum(new_responses), "within_frozen_budget": True,
            "report_sha256": sha256_file(cache / "report.json"), "vector_qc": records,
            "prediction_performance_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = verify(Path(__file__).resolve().parents[1], args.cache)
    atomic_write_json(report, args.out)
    print(report)


if __name__ == "__main__":
    main()
