"""Immutable protocol variants and split-independent feature-cache provenance."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import shutil

import pandas as pd

from .artifacts import atomic_write_json, read_json_object, sha256_file
from .cohort import canonical_assembly_digest, curate_mlc, group_hosts
from .data import load_interactions
from .host_study import verify_resources
from .split import split_by_host_group

REGISTRATION = "data/studies/genophi_seed_sensitivity.v1.json"
CACHE_LOCK = "data/studies/genophi_features_reuse.v1.json"
INPUT_FREEZE = "data/studies/genophi_seed_sensitivity.inputs.v1.json"
BASE_COMPARATOR = "data/studies/genophi_comparator.v1.json"
BASE_GENOME = "data/studies/genome_host_retrieval.v1.json"
BASE_HOST = "data/studies/genophi_host_baseline.v1.json"
CACHE_STAGES = ("cluster_strain", "features_strain", "cluster_phage", "features_phage", "merge_features")


def seed_paths(seed):
    if type(seed) is not int or seed not in (43, 44):
        raise ValueError("only the preregistered sensitivity seeds 43 and 44 are supported")
    return {
        "manifest": f"data/studies/results/genophi_host_split.seed{seed}.v1.csv",
        "host": f"data/studies/genophi_host_baseline.seed{seed}.v1.json",
        "genome": f"data/studies/genome_host_retrieval.seed{seed}.v1.json",
        "comparator": f"data/studies/genophi_comparator.seed{seed}.v1.json",
        "baseline_output": f"data/studies/results/genome_host_retrieval.seed{seed}.v1",
        "comparator_output": f"data/studies/results/genophi_comparator.seed{seed}.v1",
    }


def _bind_manifest(bindings, old_manifest, new_manifest, repo):
    result = dict(bindings)
    result.pop(old_manifest)
    result[new_manifest] = sha256_file(repo / new_manifest)
    result[REGISTRATION] = sha256_file(repo / REGISTRATION)
    return result


def expected_variant(repo, seed, kind):
    """Construct the complete allowed variant, including all inherited decisions."""
    repo = Path(repo)
    paths = seed_paths(seed)
    if kind == "host":
        result = deepcopy(read_json_object(repo / BASE_HOST))
        result["study_id"] += f":seed{seed}"
        result["split"]["seed"] = seed
        return result
    base = BASE_COMPARATOR if kind == "comparator" else BASE_GENOME
    if kind not in {"comparator", "genome"}:
        raise ValueError("unknown sensitivity protocol kind")
    result = deepcopy(read_json_object(repo / base))
    result["study_id"] += f":seed{seed}"
    result["repo_files_sha256"] = _bind_manifest(
        result["repo_files_sha256"], result["host_manifest"], paths["manifest"], repo)
    result["host_manifest"] = paths["manifest"]
    result["split_seed"] = seed
    result["sensitivity_registration"] = REGISTRATION
    result["repo_files_sha256"][base] = sha256_file(repo / base)
    if kind == "comparator":
        result["parent_protocol"] = paths["host"]
        result["repo_files_sha256"][paths["host"]] = sha256_file(repo / paths["host"])
        result["scenarios"] = ["grouped_host"]
        result["reference_protocol"] = paths["genome"]
        result["repo_files_sha256"][paths["genome"]] = sha256_file(repo / paths["genome"])
        result["reference_report"] = paths["baseline_output"] + "/report.json"
        result["reference_per_host"] = paths["baseline_output"] + "/per_host_metrics.csv"
        result["feature_artifact_lock"] = CACHE_LOCK
    else:
        result["random_split"]["seed"] = seed
    return result


def load_study_protocol(repo, path, kind):
    """Allow frozen v1 or one exact, registered derivative; never an unchecked seed flag."""
    repo = Path(repo).resolve()
    path = Path(path)
    path = path.resolve() if path.is_absolute() else (repo / path).resolve()
    base = BASE_COMPARATOR if kind == "comparator" else BASE_GENOME
    plan = read_json_object(path)
    if path == repo / base:
        return path, plan
    seed = plan.get("split_seed")
    paths = seed_paths(seed)
    if path != repo / paths[kind] or plan != expected_variant(repo, seed, kind):
        raise ValueError("sensitivity protocol differs from the registered inherited variant")
    freeze = read_json_object(repo / INPUT_FREEZE)
    if freeze["files_sha256"].get(paths[kind]) != sha256_file(path):
        raise ValueError("changed frozen sensitivity protocol")
    return path, plan


def _write_once_json(document, path):
    if path.exists():
        if read_json_object(path) != document:
            raise ValueError(f"immutable input already differs: {path}")
    else:
        atomic_write_json(document, path)


def _write_once_bytes(payload, path):
    if path.exists():
        if path.read_bytes() != payload:
            raise ValueError(f"immutable input already differs: {path}")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(payload)


def freeze_seed_inputs(repo, data_root, feature_cache, audit_path):
    """Recompute genome groups/manifests, freeze derivatives, and lock old feature bytes."""
    repo, feature_cache = Path(repo).resolve(), Path(feature_cache).resolve()
    registered = read_json_object(repo / REGISTRATION)
    for key in ("frozen_comparator", "frozen_baseline"):
        parent = registered[key]
        if sha256_file(repo / parent["path"]) != parent["sha256"]:
            raise ValueError("changed preregistered seed-42 parent")
    audit = read_json_object(audit_path)
    if audit["status"] != "verified" or audit["errors"]:
        raise ValueError("seed-42 request audit did not reproduce")
    audit_dest = repo / "data/studies/genophi_seed_sensitivity.audit.v1.json"
    _write_once_json(audit, audit_dest)
    parent = read_json_object(repo / BASE_HOST)
    lock = read_json_object(repo / "data/studies/genophi_host_baseline.resources.v1.json")
    resources = verify_resources(lock, data_root)
    observed, excluded, aliases = curate_mlc(load_interactions(resources[("table", "subset")]),
        pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0))
    hosts = sorted(observed.bacterium.unique())
    digests = {h: canonical_assembly_digest(resources[("host_genome", h)]) for h in hosts}
    groups = group_hosts(hosts,
        pd.read_csv(resources[("table", "core_distances")], sep="\t", index_col=0), digests,
        threshold=parent["grouping"]["core_distance_threshold"])
    archived = pd.read_csv(repo / "data/studies/results/genophi_host_split.v1.csv", dtype=str)
    archived = archived.sort_values("bacterium").reset_index(drop=True)
    summaries, files = {}, [audit_dest]
    for seed in [42, *registered["sensitivity_seeds"]]:
        assigned = split_by_host_group(observed, groups, seed=seed, fractions=parent["split"]["fractions"])
        manifest = assigned[["bacterium", "host_group", "split"]].drop_duplicates().sort_values("bacterium").reset_index(drop=True)
        manifest["source_bacterium"] = manifest.bacterium.map(aliases)
        manifest["assembly_sha256"] = manifest.bacterium.map(digests)
        if seed == 42:
            if not manifest.equals(archived):
                raise ValueError("recomputed seed-42 genome groups/manifest differ")
            continue
        paths = seed_paths(seed)
        destination = repo / paths["manifest"]
        _write_once_bytes(manifest.to_csv(index=False, lineterminator="\n").encode(), destination)
        files.append(destination)
        for kind in ("host", "genome", "comparator"):
            destination = repo / paths[kind]
            _write_once_json(expected_variant(repo, seed, kind), destination)
            files.append(destination)
        previous = archived.set_index("host_group").split.to_dict()
        current = manifest.set_index("host_group").split.to_dict()
        summaries[str(seed)] = {
            "partitions": {p: {"hosts": int(part.bacterium.nunique()), "groups": int(part.host_group.nunique()),
                               "observed_pairs": len(part)} for p, part in assigned.groupby("split")},
            "groups_changed_partition_vs42": sum(current[g] != previous[g] for g in previous),
            "changed_groups": [g for g in sorted(previous) if current[g] != previous[g]],
        }
    committed = repo / "data/studies/results/genophi_comparator.v1"
    artifacts, stage_hashes = {}, {}
    for name in CACHE_STAGES:
        checkpoint = committed / "checkpoints" / f"{name}.json"
        stage_hashes[name] = sha256_file(checkpoint)
        record = read_json_object(checkpoint)
        if sha256_file(feature_cache / "checkpoints" / f"{name}.json") != stage_hashes[name]:
            raise ValueError("local source checkpoint differs from committed v1")
        for relative, digest in record["output_sha256"].items():
            path = (feature_cache / relative).resolve()
            if not path.is_relative_to(feature_cache) or sha256_file(path) != digest:
                raise ValueError(f"changed or missing source feature artifact: {relative}")
            artifacts[relative] = {"sha256": digest, "size_bytes": path.stat().st_size}
    source_provenance = read_json_object(committed / "provenance.json")
    if read_json_object(feature_cache / "provenance.json") != source_provenance:
        raise ValueError("source feature-cache provenance differs from v1")
    cache = {
        "schema_version": 1, "source_study_id": "genophi-frozen-host-comparator:v1",
        "source_protocol_sha256": registered["frozen_comparator"]["sha256"],
        "source_provenance_sha256": sha256_file(committed / "provenance.json"),
        "proteome_receipt_sha256": source_provenance["proteome_receipt_sha256"],
        "versions": source_provenance["versions"], "source_checkpoint_sha256": stage_hashes,
        "artifacts": artifacts, "scope": "The same all-genome vocabulary and full observed-row feature table; no seed-specific supervised artifact is admitted",
    }
    _write_once_json(cache, repo / CACHE_LOCK)
    files.append(repo / CACHE_LOCK)
    archive = repo / "data/studies/code_snapshots/verify_genophi_comparator.v1.py"
    validator = repo / "scripts/verify_genophi_comparator.py"
    expected = read_json_object(committed / "NUMERIC_QA.json")["verification_script_sha256"]
    if not archive.exists():
        if sha256_file(validator) != expected:
            raise ValueError("archive the original validator before extending it")
        _write_once_bytes(validator.read_bytes(), archive)
    if sha256_file(archive) != expected:
        raise ValueError("archived v1 validator differs from historical QA")
    files.append(archive)
    result = {
        "status": "inputs_frozen_before_sensitivity_fit", "study_id": registered["study_id"],
        "registration_sha256": sha256_file(repo / REGISTRATION), "verified_resources": len(resources),
        "observed_pairs": len(observed), "unknown_pairs_excluded": len(excluded), "host_groups": len(set(groups.values())),
        "seed_splits": summaries,
        "files_sha256": {p.relative_to(repo).as_posix(): sha256_file(p) for p in files},
    }
    _write_once_json(result, repo / INPUT_FREEZE)
    return result


def import_feature_cache(repo, source, destination, proteome_receipt_sha256):
    """Verify and materialize original feature-stage artifacts without any reclustering."""
    repo, source, destination = Path(repo), Path(source).resolve(), Path(destination).resolve()
    lock_path = repo / CACHE_LOCK
    lock = read_json_object(lock_path)
    frozen = read_json_object(repo / INPUT_FREEZE)
    if sha256_file(lock_path) != frozen["files_sha256"][CACHE_LOCK]:
        raise ValueError("changed frozen feature reuse lock")
    if lock["proteome_receipt_sha256"] != proteome_receipt_sha256:
        raise ValueError("feature cache was built from different proteomes")
    if sha256_file(source / "provenance.json") != lock["source_provenance_sha256"]:
        raise ValueError("feature cache provenance differs from the frozen source")
    for stage, digest in lock["source_checkpoint_sha256"].items():
        if sha256_file(source / "checkpoints" / f"{stage}.json") != digest:
            raise ValueError("changed feature source checkpoint")
    outputs = []
    for relative, entry in lock["artifacts"].items():
        path, target = (source / relative).resolve(), (destination / relative).resolve()
        if not path.is_relative_to(source) or not target.is_relative_to(destination):
            raise ValueError("feature artifact escapes its root")
        if path.stat().st_size != entry["size_bytes"] or sha256_file(path) != entry["sha256"]:
            raise ValueError(f"changed source feature bytes: {relative}")
        if target.exists():
            if sha256_file(target) != entry["sha256"]:
                raise ValueError("existing imported feature artifact differs")
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_name(target.name + ".import-partial")
            with path.open("rb") as original, temporary.open("xb") as copied:
                shutil.copyfileobj(original, copied, 1024 * 1024)
            if sha256_file(temporary) != entry["sha256"]:
                raise ValueError("feature copy digest mismatch")
            temporary.rename(target)
        outputs.append(target)
    return outputs, {"status": "reused_v1_feature_stages", "source_checkpoint_stages": list(CACHE_STAGES),
                     "feature_artifact_lock_sha256": sha256_file(lock_path), "artifacts": len(outputs),
                     "new_clustering_search_assignment_calls": 0}
