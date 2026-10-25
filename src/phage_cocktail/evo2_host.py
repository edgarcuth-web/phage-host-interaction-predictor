"""Deterministic sampled host windows and an integrity-bound embedding cache.

The encoder is injected: CPU tests never synthesize claimed Evo2 results.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd

from .artifacts import atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory
from .cohort import canonical_assembly_digest, parse_dna_fasta, reverse_complement
from .host_genome import CANONICAL_WORDS, host_fourmers
from .host_study import verify_resources

PROTOCOL = "data/studies/evo2_host_representation.v1.json"


def host_windows(payload, *, window_bp, max_windows_per_host):
    """Yield deterministic (contig ordinal, start, sequence) without crossing contigs."""
    if any(type(v) is not int or v < 1 for v in (window_bp, max_windows_per_host)):
        raise ValueError("window parameters must be positive integers")
    contigs = sorted(min(c, reverse_complement(c)) for c in parse_dna_fasta(payload.decode("ascii")))
    windows = [(i, start, c[start:start + window_bp]) for i, c in enumerate(contigs)
               for start in range(0, len(c), window_bp)]
    count = min(max_windows_per_host, len(windows))
    selected = ([0] if count == 1 else
                [i * (len(windows) - 1) // (count - 1) for i in range(count)])
    return [(windows[i][0], windows[i][1],
             "".join(b if b in "ACGT" else "N" for b in windows[i][2])) for i in selected]


def validate_deployment(identity):
    for name in ("checkpoint_sha256", "config_sha256", "evo2_module_sha256", "vortex_module_sha256"):
        if not re.fullmatch(r"[0-9a-f]{64}", str(identity.get(name, ""))):
            raise ValueError(f"deployment requires a real {name}")
    for name in ("model_name", "evo2_version", "vortex_version", "torch_version",
                 "cuda_version", "transformer_engine_version", "hardware", "precision"):
        if not isinstance(identity.get(name), str) or not identity[name].strip():
            raise ValueError(f"deployment requires {name}")
    if identity["model_name"] != "evo2_40b":
        raise ValueError("this protocol requires evo2_40b")


def extract_cache(repo, data_root, output, encode, deployment, *, pilot=False, runtime_metrics=None):
    """Extract from verified local FASTAs; no downloads, labels or fitting here.

    ``encode(sequence)`` returns a token-by-hidden numpy array for the fixed layer.
    Every selected orientation must return exactly one vector per input base.
    """
    repo, output = Path(repo), Path(output)
    if output.exists():
        raise ValueError("cache output must be new")
    validate_deployment(deployment)
    plan = read_json_object(repo / PROTOCOL)
    spec = plan["embedding"]
    for path, digest in plan["repo_files_sha256"].items():
        if sha256_file(repo / path) != digest:
            raise ValueError(f"changed frozen embedding input: {path}")
    manifest = pd.read_csv(repo / plan["host_manifest"], dtype=str).sort_values("bacterium")
    if pilot:
        manifest = manifest.head(3)
    genome_plan = read_json_object(repo / "data/studies/genome_host_retrieval.v1.json")
    entries = [r for r in genome_plan["resources"]
               if r["kind"] == "host_genome" and r["id"] in set(manifest.bacterium)]
    resources = verify_resources({"resources": entries}, data_root)
    vectors, fourmers, qc, window_records = [], [], [], []
    with staged_directory(output) as stage:
        # Snapshot all settings before the first model call.
        atomic_write_json(deployment, stage / "deployment.json")
        (stage / "protocol.json").write_bytes((repo / PROTOCOL).read_bytes())
        for row in manifest.itertuples(index=False):
            path = resources[("host_genome", row.bacterium)]
            if canonical_assembly_digest(path) != row.assembly_sha256:
                raise ValueError("assembly differs from frozen host manifest")
            payload = path.read_bytes()
            composition, genome_qc = host_fourmers(payload)
            windows = host_windows(payload, window_bp=spec["window_bp"],
                                   max_windows_per_host=spec["max_windows_per_host"])
            total, vector = 0, np.zeros(spec["dimension"], dtype=np.float64)
            started = time.monotonic()
            for contig, start, sequence in windows:
                orientation_means = []
                for oriented in (sequence, reverse_complement(sequence)):
                    values = np.asarray(encode(oriented), dtype=np.float32)
                    if values.shape != (len(sequence), spec["dimension"]) or not np.isfinite(values).all():
                        raise ValueError("embedding token count, hidden dimension or finiteness differs")
                    orientation_means.append(values.mean(axis=0, dtype=np.float32))
                vector += np.mean(orientation_means, axis=0, dtype=np.float32).astype(np.float64) * len(sequence)
                total += len(sequence)
                window_records.append({"bacterium": row.bacterium, "contig": contig,
                                       "start": start, "length": len(sequence),
                                       "forward_input_sha256": sha256(sequence.encode("ascii")).hexdigest()})
            vectors.append((vector / total).astype(np.float32))
            fourmers.append(composition)
            qc.append({"bacterium": row.bacterium, "assembly_sha256": row.assembly_sha256,
                       "fasta_sha256": sha256_file(path), "selected_windows": len(windows),
                       "selected_bases": total, "sampled_fraction": total / genome_qc["bases"],
                       "encoding_seconds": time.monotonic() - started, **genome_qc})
            print(f"Encoded {row.bacterium}: {len(windows)} windows, {total} bases per strand", flush=True)
        np.save(stage / "embeddings.npy", np.stack(vectors), allow_pickle=False)
        pd.DataFrame(fourmers, index=manifest.bacterium,
                     columns=["host_4mer:" + w for w in CANONICAL_WORDS]).to_csv(
                         stage / "fourmers.csv", index_label="bacterium")
        pd.DataFrame(qc).to_csv(stage / "qc.csv", index=False)
        pd.DataFrame(window_records).to_csv(stage / "windows.csv", index=False)
        atomic_write_json(runtime_metrics() if runtime_metrics else {}, stage / "runtime.json")
        atomic_write_json({
            "schema_version": 1, "pilot": pilot, "hosts": manifest.bacterium.tolist(),
            "protocol_sha256": sha256_file(repo / PROTOCOL),
            **code_provenance(repo, repo / "scripts/extract_evo2_hosts.py"),
            "files_sha256": {p.name: sha256_file(p) for p in sorted(stage.iterdir())},
        }, stage / "cache.json")


def load_cache(cache, repo, manifest):
    """Require a complete aligned cache; never silently impute absent embeddings."""
    cache, repo = Path(cache), Path(repo)
    lock = read_json_object(cache / "cache.json")
    if lock.get("schema_version") != 1 or lock.get("pilot") is not False:
        raise ValueError("evaluation requires a full cache, not a pilot")
    if lock["protocol_sha256"] != sha256_file(repo / PROTOCOL):
        raise ValueError("embedding protocol hash mismatch")
    required = {"protocol.json", "deployment.json", "embeddings.npy", "fourmers.csv", "qc.csv", "windows.csv", "runtime.json"}
    if set(lock["files_sha256"]) != required:
        raise ValueError("cache requires the complete payload")
    for name, digest in lock["files_sha256"].items():
        if sha256_file(cache / name) != digest:
            raise ValueError(f"embedding cache hash mismatch: {name}")
    if sha256_file(cache / "protocol.json") != lock["protocol_sha256"]:
        raise ValueError("cached protocol differs")
    validate_deployment(read_json_object(cache / "deployment.json"))
    hosts = lock["hosts"]
    if hosts != sorted(manifest.bacterium) or len(set(hosts)) != len(hosts):
        raise ValueError("embedding cache host panel differs")
    qc = pd.read_csv(cache / "qc.csv", dtype={"bacterium": str})
    if qc.bacterium.tolist() != hosts:
        raise ValueError("embedding QC host alignment differs")
    expected = manifest.set_index("bacterium").loc[hosts, "assembly_sha256"].tolist()
    if qc.assembly_sha256.tolist() != expected:
        raise ValueError("embedding assembly identities differ")
    source = read_json_object(repo / "data/studies/genome_host_retrieval.v1.json")
    fasta_hashes = {r["id"]: r["sha256"] for r in source["resources"] if r["kind"] == "host_genome"}
    if qc.fasta_sha256.tolist() != [fasta_hashes[h] for h in hosts]:
        raise ValueError("embedding source bytes differ")
    dimension = read_json_object(repo / PROTOCOL)["embedding"]["dimension"]
    values = np.load(cache / "embeddings.npy", allow_pickle=False)
    if values.dtype != np.float32 or values.shape != (len(hosts), dimension) or not np.isfinite(values).all():
        raise ValueError("invalid embedding matrix")
    features = pd.DataFrame(values, index=hosts, columns=[f"evo2:{i}" for i in range(dimension)])
    fourmers = pd.read_csv(cache / "fourmers.csv", index_col="bacterium")
    if (fourmers.index.tolist() != hosts or fourmers.columns.tolist() != ["host_4mer:" + w for w in CANONICAL_WORDS]
            or not np.isfinite(fourmers.to_numpy()).all() or (fourmers.to_numpy() < 0).any()
            or not np.allclose(fourmers.sum(axis=1), 1)):
        raise ValueError("invalid cached genome composition")
    return features, fourmers, qc, sha256_file(cache / "cache.json")
