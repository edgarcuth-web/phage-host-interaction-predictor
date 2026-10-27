"""One real-host 8192-base feasibility request, with a frozen explicit budget."""
import argparse
from datetime import datetime, timezone
import getpass
from hashlib import sha256
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory
from phage_cocktail.cohort import canonical_assembly_digest
from phage_cocktail.evo2_host import host_windows
from phage_cocktail.host_study import verify_resources
from phage_cocktail.nvidia_evo2 import NvidiaForwardError, forward

PROTOCOL = "data/studies/nvidia_evo2_window_feasibility.v1.json"


def selected_window(repo, data_root, plan):
    """Verify one frozen host and select a full window before looking at outcomes."""
    for path, digest in plan["repo_files_sha256"].items():
        if sha256_file(repo / path) != digest:
            raise ValueError(f"changed frozen input: {path}")
    genomes = read_json_object(repo / plan["genome_protocol"])
    manifest = pd.read_csv(repo / genomes["host_manifest"], dtype=str).sort_values("bacterium")
    row = manifest.iloc[0]
    if row.bacterium != plan["host"]:
        raise ValueError("frozen feasibility host no longer matches manifest")
    source = next(r for r in genomes["resources"] if r["kind"] == "host_genome" and r["id"] == row.bacterium)
    files = verify_resources({"resources": [source]}, data_root)
    path = files[("host_genome", row.bacterium)]
    if canonical_assembly_digest(path) != row.assembly_sha256:
        raise ValueError("feasibility host assembly identity changed")
    windows = host_windows(path.read_bytes(), window_bp=plan["window_bp"],
                           max_windows_per_host=plan["sampling_max_windows"])
    selected = next(((i, c, start, s) for i, (c, start, s) in enumerate(windows)
                     if len(s) == plan["window_bp"]), None)
    if selected is None:
        raise ValueError("no full-sized window in frozen host sampling")
    i, contig, start, sequence = selected
    return sequence, {"host": row.bacterium, "assembly_sha256": row.assembly_sha256,
                      "fasta_sha256": source["sha256"], "selected_window_ordinal": i,
                      "contig_ordinal": contig, "start": start, "length": len(sequence),
                      "input_sha256": sha256(sequence.encode("ascii")).hexdigest(),
                      "ambiguous_bases_mapped_to_N": sequence.count("N")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        raise ValueError("feasibility output directory must be new")
    repo = Path(__file__).resolve().parents[1]
    plan = read_json_object(repo / PROTOCOL)
    sequence, source = selected_window(repo, args.data_root, plan)
    identity = {"protocol_sha256": sha256_file(repo / PROTOCOL),
                **code_provenance(repo, Path(__file__)), "input": source}
    print({"input": source, "request_budget": plan["request_budget"]}, flush=True)
    api_key = os.environ.get("NVIDIA_API_KEY") or getpass.getpass("NVIDIA API key (hidden): ")
    with staged_directory(args.outdir) as stage:
        atomic_write_json(plan, stage / "protocol.json")
        started = time.monotonic()
        result = None
        try:
            result = forward(sequence, api_key,
                             timeout=plan["request_budget"]["io_timeout_seconds"],
                             poll_seconds=plan["request_budget"]["server_wait_seconds"],
                             max_response_bytes=plan["request_budget"]["max_response_bytes"],
                             max_tensor_bytes=plan["request_budget"]["max_tensor_bytes"])
            pooled = result.tokens.mean(axis=0, dtype=np.float64)
            stored = pooled.astype(np.float32)
            if not np.isfinite(stored).all() or not np.any(stored):
                raise ValueError("pooled vector is nonfinite or all zero")
            np.save(stage / "pooled_mean.npy", stored, allow_pickle=False)
            report = {"status": "success", "response": result.metadata,
                      "pooled_shape": list(stored.shape), "pooled_dtype": str(stored.dtype),
                      "pooled_l2": float(np.linalg.norm(pooled)),
                      "float32_pool_max_abs_error": float(np.max(np.abs(pooled - stored))),
                      "token_rms": float(np.sqrt(np.mean(result.tokens ** 2))),
                      "transfer_extrapolation": {
                          "basis": "One measured real-host 8192-base response; other windows/hosts can compress differently",
                          "three_hosts_32windows_two_strands_json_bytes": result.metadata["response_bytes"] * 3 * 32 * 2,
                          "full_177hosts_32windows_two_strands_json_bytes": result.metadata["response_bytes"] * 177 * 32 * 2}}
        except NvidiaForwardError as exc:
            report = {"status": "request_not_completed", "error": str(exc),
                      "http_status": exc.http_status, "request_id": exc.request_id,
                      "response_bytes_read": exc.response_bytes,
                      "automatic_retry_submitted": False}
        finally:
            api_key = None
        report.update(observed_at_utc=datetime.now(timezone.utc).isoformat(),
                      client_elapsed_seconds=time.monotonic() - started,
                      provenance=identity, inference_post_attempts=1,
                      scope="Single-window feasibility only; not a host embedding, full pilot or prediction benchmark",
                      output_sha256={p.name: sha256_file(p) for p in sorted(stage.iterdir())})
        atomic_write_json(report, stage / "report.json")
        print(report, flush=True)


if __name__ == "__main__":
    main()
