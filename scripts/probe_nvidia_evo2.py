"""Make four bounded real forward requests; never store the NVIDIA API key."""
import argparse
from datetime import datetime, timezone
import getpass
import os
from pathlib import Path

import numpy as np

from phage_cocktail.artifacts import atomic_write_json, code_provenance, staged_directory
from phage_cocktail.nvidia_evo2 import DIMENSION, forward


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        raise ValueError("probe output directory must be new")
    # Environment injection or a hidden prompt, never a command-line secret.
    api_key = os.environ.get("NVIDIA_API_KEY") or getpass.getpass("NVIDIA API key (hidden): ")
    sequences = ["ACGTTGCATGTCGCAT", "TGCACGATACGTAGCT", "ACGTTGCATGTCGCAT", "ACGTTGCATGTCGCAT" * 16]
    results = []
    try:
        for sequence in sequences:
            result = forward(sequence, api_key)
            results.append(result)
            print({k: result.metadata[k] for k in ("http_status", "sequence_length", "shape", "returned_dtype", "elapsed_ms", "response_bytes")}, flush=True)
    finally:
        api_key = None
    first, second, repeated, longer = results
    if np.array_equal(first.tokens, second.tokens):
        raise ValueError("different input controls returned identical layer outputs")
    means = [r.tokens.mean(axis=0, dtype=np.float64) for r in results]
    pooled_difference = float(np.linalg.norm(means[0] - means[1]))
    if pooled_difference == 0:
        raise ValueError("different input controls returned identical pooled representations")
    repo = Path(__file__).resolve().parents[1]
    for result, mean in zip(results, means):
        result.metadata.update(token_rms=float(np.sqrt(np.mean(result.tokens ** 2))),
                               pooled_l2=float(np.linalg.norm(mean)),
                               min_value=float(result.tokens.min()), max_value=float(result.tokens.max()))
    with staged_directory(args.outdir) as stage:
        atomic_write_json({
            "purpose": "Hosted API connectivity, layer shape, input sensitivity and transfer feasibility only; not a genome or host-prediction experiment",
            "observed_at_utc": datetime.now(timezone.utc).isoformat(),
            **code_provenance(repo, Path(__file__)),
            "requests": [r.metadata for r in results],
            "distinct_input_pooled_l2": pooled_difference,
            "distinct_input_pooled_relative_l2": pooled_difference / max(np.linalg.norm(means[0]), np.linalg.norm(means[1])),
            "same_input_repeat_pooled_l2": float(np.linalg.norm(means[0] - means[2])),
            "same_input_repeat_max_abs": float(np.max(np.abs(first.tokens - repeated.tokens))),
            "same_final_base_context_l2": float(np.linalg.norm(first.tokens[-1] - second.tokens[-1])),
            "token_axis_varies": [bool(np.any(r.tokens[1:] != r.tokens[:-1])) for r in results],
            "hosted_service_weight_identity_verified": False,
            "local_v1_cache_compatible": False,
            "transfer_extrapolation": {
                "basis": "256-base response bytes linearly scaled by input bases; rough estimate, not a measured 8192-base response",
                "pilot_3hosts_32windows_8192bp_2strands_json_bytes": longer.metadata["response_bytes"] * (3 * 32 * 8192 * 2 / 256),
                "full_177hosts_32windows_8192bp_2strands_json_bytes": longer.metadata["response_bytes"] * (177 * 32 * 8192 * 2 / 256),
                "single_8192bp_float64_tensor_bytes": 8192 * DIMENSION * 8,
            },
        }, stage / "report.json")


if __name__ == "__main__":
    main()
