"""Research FASTA-to-ranking command for two explicit fixed-library baselines."""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .artifacts import canonical_content_sha256, canonical_json, sha256_file
from .cohort import parse_dna_fasta, reverse_complement
from .host_baseline import predict_host_model
from .host_genome import CANONICAL_WORDS, host_fourmers

DEFAULT_ROOT = Path(__file__).resolve().parents[2]
REGISTRY = "models/research_ranker.v1.json"
REGISTRY_SHA256 = "12e6e8adb928f7d9dcdbb156c3a7ea3fcabaa86827103051bbcadedcba7d20bb"
METHODS = ("host-4mer", "prevalence")
MAX_FASTA_BYTES = 64 * 1024 * 1024


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _inside(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    _require(path.is_relative_to(root), "model artifact path escapes model root")
    return path


def _validate_weights(weights, registry):
    ids = {r["phage_id"] for r in registry["library"]["candidates"]}
    _require(set(weights["phages"]) == ids, "model does not support the exact frozen candidate library")
    expected_training = {r["host_id"] for r in registry["training"]["known_cohort_members"] if r["training_member"]}
    _require(set(weights["training_hosts"]) == expected_training and len(weights["training_hosts"]) == 117,
             "model training-host scope differs")
    enc = weights["preprocessing"]
    columns = ["host_4mer:" + word for word in CANONICAL_WORDS]
    _require(enc["numeric_columns"] == columns and not enc["categorical_columns"], "model feature schema/order differs")
    _require(len(enc["means"]) == len(enc["scales"]) == 136
             and np.isfinite(enc["means"]).all() and np.isfinite(enc["scales"]).all()
             and (np.asarray(enc["scales"]) > 0).all(), "invalid saved preprocessing")
    _require(weights["score_semantics"] == "uncalibrated_source_mlc_positive", "model source endpoint differs")
    counts = registry["methods"]["prevalence"]["training_counts"]
    for name, model in weights["phages"].items():
        _require(model["training_rows"] == counts[name]["observed"], "model training-label support differs")
        if model["kind"] == "logistic":
            _require(len(model["coefficients"]) == 136 and np.isfinite(model["coefficients"]).all()
                     and np.isfinite(model["intercept"]), "invalid logistic coefficients")
        elif model["kind"] == "single_class_prevalence":
            expected = (counts[name]["positive"] + 1) / (counts[name]["observed"] + 2)
            _require(np.isfinite(model["score"]) and model["score"] == expected, "invalid single-class fallback")
        else:
            raise ValueError("unsupported saved estimator kind")


def load_ranker(model_root=DEFAULT_ROOT, *, method):
    """Load only pinned JSON data; never fit a model or execute serialized code."""
    _require(method in METHODS, "unsupported research method")
    payload = _inside(model_root, REGISTRY).read_bytes()
    _require(sha256(payload).hexdigest() == REGISTRY_SHA256, "research registry differs from its pinned version")
    registry = json.loads(payload)
    _require(registry["schema_version"] == 1 and registry["status"] == "research_only_not_promoted", "unsupported registry version/status")
    library = registry["library"]
    candidates = library["candidates"]
    ids = [r["phage_id"] for r in candidates]
    _require(len(ids) == len(set(ids)) == library["candidate_count"] == 94, "candidate roster must contain exactly 94 unique source IDs")
    _require(len({r["assembly_sha256"] for r in candidates}) == library["unique_assemblies"] == 89,
             "candidate assembly identity count differs")
    digest = sha256(canonical_json(dict(library_id=library["library_id"], candidates=candidates)).encode()).hexdigest()
    _require(digest == library["snapshot_sha256"], "candidate snapshot digest differs")
    counts = registry["methods"]["prevalence"]["training_counts"]
    _require(set(counts) == set(ids), "prevalence does not support the exact candidate library")
    _require(all(type(r["positive"]) is int and type(r["observed"]) is int
                 and 0 <= r["positive"] <= r["observed"] and r["observed"] > 0 for r in counts.values()),
             "invalid prevalence training counts")
    _require(sum(r["observed"] for r in counts.values()) == registry["training"]["observed_pairs"] == 10957,
             "prevalence training observation count differs")
    _require(sha256(canonical_json(counts).encode()).hexdigest() == registry["methods"]["prevalence"]["training_counts_sha256"],
             "prevalence training-count digest differs")
    weights = None
    if method == "host-4mer":
        spec = registry["methods"][method]["weights"]
        data = _inside(model_root, spec["path"]).read_bytes()
        _require(len(data) == spec["size_bytes"] and sha256(data).hexdigest() == spec["sha256"], "host 4-mer weights differ from their pinned bytes")
        weights = json.loads(data)
        _validate_weights(weights, registry)
    return dict(registry=registry, weights=weights, method=method)


def score_feature_vector(bundle, vector=None):
    """Score the full known library; this helper makes no FASTA-identity claim."""
    registry, method = bundle["registry"], bundle["method"]
    phages = [r["phage_id"] for r in registry["library"]["candidates"]]
    if method == "prevalence":
        counts = registry["methods"][method]["training_counts"]
        return {p: (counts[p]["positive"] + 1) / (counts[p]["observed"] + 2) for p in phages}
    values = np.asarray(vector, dtype=float)
    _require(values.shape == (136,) and np.isfinite(values).all() and (values >= 0).all()
             and abs(float(values.sum()) - 1) <= 1e-12, "host 4-mer method requires the complete normalized 136-value vector")
    columns = ["host_4mer:" + word for word in CANONICAL_WORDS]
    hosts = pd.DataFrame([values], index=["query"], columns=columns)
    pairs = pd.DataFrame(dict(bacterium=["query"] * len(phages), phage=phages))
    with threadpool_limits(limits=1):
        scores = predict_host_model(bundle["weights"], pairs, hosts)
    _require(scores.shape == (94,) and np.isfinite(scores).all(), "incomplete or nonfinite model score support")
    return dict(zip(phages, map(float, scores)))


def candidate_rows(library, scores, top_k, *, reason=None):
    _require(type(top_k) is int and 1 <= top_k <= 94, "top-k must be an integer from 1 to 94")
    roster = {r["phage_id"]: r for r in library["candidates"]}
    if reason is not None:
        return [dict(**roster[p], status="not_estimated", reason=reason, score=None,
                     rank_min=None, rank_max=None, display_rank=None, expected_top_k_slot=None)
                for p in sorted(roster)]
    _require(set(scores) == set(roster) and all(np.isfinite(v) for v in scores.values()),
             "ranking requires exactly the complete finite candidate panel")
    ordered = sorted(roster, key=lambda p: (-scores[p], p))
    groups = {}
    for i, phage in enumerate(ordered, start=1):
        groups.setdefault(scores[phage], []).append(i)
    rows = []
    for i, phage in enumerate(ordered, start=1):
        ranks = groups[scores[phage]]
        first, last = ranks[0], ranks[-1]
        fraction = min(max(top_k - first + 1, 0), len(ranks)) / len(ranks)
        rows.append(dict(**roster[phage], status="scored", reason=None, score=float(scores[phage]),
                         rank_min=first, rank_max=last, display_rank=i, expected_top_k_slot=fraction))
    return rows


def rank_fasta(fasta, *, method, top_k=5, model_root=DEFAULT_ROOT, host_species=None):
    bundle = load_ranker(model_root, method=method)
    registry = bundle["registry"]
    _require(type(top_k) is int and 1 <= top_k <= 94, "top-k must be an integer from 1 to 94")
    path = Path(fasta)
    with path.open("rb") as stream:
        payload = stream.read(MAX_FASTA_BYTES + 1)
    _require(len(payload) <= MAX_FASTA_BYTES, "FASTA exceeds the 64 MiB input limit")
    contigs = parse_dna_fasta(payload.decode("ascii"), path.name)
    contig_digests = [sha256(min(c, reverse_complement(c)).encode()).hexdigest() for c in contigs]
    assembly = sha256("\n".join(sorted(contig_digests)).encode()).hexdigest()
    raw_digest = sha256(payload).hexdigest()
    qc = dict(contigs=len(contigs), bases=sum(map(len, contigs)),
              ambiguous_bases=sum(sum(base not in "ACGT" for base in c) for c in contigs),
              possible_windows=sum(max(0, len(c)-3) for c in contigs), valid_windows=0)
    vector = None
    if any(re.search("[ACGT]{4}", c) for c in contigs):
        vector, qc = host_fourmers(payload)
    members = [dict(**r, exact_source_bytes=r["sequence_sha256"] == raw_digest)
               for r in registry["training"]["known_cohort_members"] if r["assembly_sha256"] == assembly]
    envelope = registry["applicability"]["observed_training_bases"]
    ambiguity_fraction = qc["ambiguous_bases"] / qc["bases"]
    supported_species = host_species is None or host_species.strip().casefold() in ("escherichia coli", "e. coli")
    reason = None if supported_species else "unsupported_host_species"
    if reason is None and method == "host-4mer" and vector is None:
        reason = "no_valid_fourmer_windows"
    feature = dict(representation_id=registry["representation"]["id"], representation_sha256=registry["representation"]["sha256"],
                   values_sha256=None if vector is None else sha256(np.asarray(vector, dtype="<f8").tobytes()).hexdigest(),
                   status="available" if vector is not None else "not_estimated")
    feature["feature_key"] = None if vector is None else sha256(canonical_json(dict(assembly_sha256=assembly, representation_sha256=feature["representation_sha256"])).encode()).hexdigest()
    scores = None if reason else score_feature_vector(bundle, vector)
    rows = candidate_rows(registry["library"], scores, top_k, reason=reason)
    method_spec = registry["methods"][method]
    dependencies = ("research_ranker.py", "host_genome.py", "host_baseline.py", "cohort.py", "artifacts.py")
    result = dict(schema_version=1, kind="research_candidate_ranking", status="not_estimated" if reason else "scored", reason=reason,
                  query=dict(input_name=path.name, sequence_sha256=raw_digest, assembly_sha256=assembly, qc=qc, features=feature),
                  model=dict(method=method, model_id=method_spec["model_id"], registry_id=registry["registry_id"], registry_sha256=REGISTRY_SHA256,
                             artifact_sha256=method_spec["weights"]["sha256"] if method == "host-4mer" else method_spec["training_counts_sha256"],
                             genome_conditioned=method_spec["genome_conditioned"], status=registry["status"], training_hosts=117, training_observed_pairs=10957),
                  endpoint=registry["target"], library=dict(library_id=registry["library"]["library_id"], snapshot_sha256=registry["library"]["snapshot_sha256"],
                                                          candidate_count=94, unique_assemblies=89),
                  applicability=dict(declared_host_species=host_species, taxon_provenance="caller_declared" if host_species else "source_dataset_label" if members else "not_verified",
                                     taxon_verified_by_this_command=False, known_cohort_matches=members,
                                     within_observed_training_length=envelope["min"] <= qc["bases"] <= envelope["max"], observed_training_bases=envelope,
                                     ambiguity_fraction=ambiguity_fraction, within_observed_training_ambiguity=ambiguity_fraction <= registry["applicability"]["observed_training_max_ambiguity_fraction"],
                                     range_check_semantics=registry["applicability"]["range_checks"], assembly_completeness_verified=False),
                  coverage=dict(considered_candidates=94, scored_candidates=0 if reason else 94, not_estimated_candidates=94 if reason else 0,
                                filtered_candidates=0, filter_reasons=[]),
                  ranking=dict(requested_top_k=top_k, score_ties="exact equality; rank_min/rank_max retain all tied positions",
                               display_ties="lexicographic phage ID; named top-k is an inspection list", uncertainty="not_estimated; no calibrated query uncertainty model"),
                  top_candidates=[] if reason else rows[:top_k], candidates=rows,
                  provenance=dict(source_model_commit=registry["source_model_commit"], runtime={"python":sys.version.split()[0],
                                  **{p:importlib.metadata.version(p) for p in ("numpy", "pandas", "scikit-learn")}},
                                  scorer_files_sha256={p:sha256_file(Path(__file__).with_name(p)) for p in dependencies}),
                  limitations=registry["limitations"])
    result["content_sha256"] = canonical_content_sha256(result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True, help="assembled ASCII IUPAC DNA FASTA")
    parser.add_argument("--method", choices=METHODS, required=True, help="explicit research baseline; prevalence does not condition on the genome")
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--host-species", help="optional caller declaration; an explicitly unsupported species yields not_estimated")
    parser.add_argument("--model-root", type=Path, default=DEFAULT_ROOT, help="repository or relocated bundle root containing the pinned registry/artifact paths")
    parser.add_argument("--out", type=Path, required=True, help="new JSON output file")
    args = parser.parse_args(argv)
    try:
        _require(not args.out.exists(), "output must be new; refusing to overwrite")
        result = rank_fasta(args.fasta, method=args.method, top_k=args.top_k, model_root=args.model_root, host_species=args.host_species)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        with args.out.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(f"research ranking failed: {error}", file=sys.stderr)
        return 2
    print(f"wrote {args.out}; method={args.method}; status={result['status']}; scored={result['coverage']['scored_candidates']}/94")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
