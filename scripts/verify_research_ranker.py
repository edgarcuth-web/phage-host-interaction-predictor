#!/usr/bin/env python3
"""Check research inference against saved study scores and an optional real FASTA."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, canonical_content_sha256, read_json_object, sha256_file
from phage_cocktail.host_genome import CANONICAL_WORDS
from phage_cocktail.research_ranker import REGISTRY_SHA256, load_ranker, rank_fasta, score_feature_vector


def frame(path):
    return pd.read_csv(path, dtype={"bacterium":str, "phage":str, "entity_id":str}, float_precision="round_trip")


def verify(root, *, fasta=None):
    root = Path(root).resolve()
    bundles = {m:load_ranker(root, method=m) for m in ("host-4mer", "prevalence")}
    registry = bundles["host-4mer"]["registry"]
    for item in registry["provenance"].values():
        if sha256_file(root / item["path"]) != item["sha256"]:
            raise ValueError("changed source-study provenance")
    plan = read_json_object(root / registry["provenance"]["source_protocol"]["path"])
    source_report = read_json_object(root / registry["provenance"]["source_report"]["path"])
    vectors = plan["inputs"]["sequence_vectors"]
    if sha256_file(root / vectors["path"]) != vectors["sha256"]:
        raise ValueError("changed source feature vectors")
    scores_path = root / "data/studies/results/new_host_retrieval.v1/scores.csv.gz"
    if sha256_file(scores_path) != source_report["artifacts_sha256"]["scores.csv.gz"]:
        raise ValueError("changed source score matrix")
    features = frame(root / vectors["path"]).set_index(["entity_kind", "entity_id"])
    scores = frame(scores_path)
    scores = scores[scores.panel.eq("grouped_host")]
    observed = frame(root / registry["provenance"]["observations"]["path"])
    train_ids = {r["host_id"] for r in registry["training"]["known_cohort_members"] if r["training_member"]}
    train = observed[observed.bacterium.isin(train_ids)]
    stats = train.groupby("phage").infection.agg(["sum", "count"])
    counts = registry["methods"]["prevalence"]["training_counts"]
    if len(train) != 10957 or set(stats.index) != set(counts):
        raise ValueError("training count scope differs")
    for p, row in stats.iterrows():
        if counts[p] != dict(positive=int(row["sum"]), observed=int(row["count"])):
            raise ValueError("prevalence does not match source training labels")
    checked, maximum_error = 0, 0.
    for method, source_method in (("host-4mer", "host_4mer"), ("prevalence", "phage_prevalence")):
        reference = scores[scores.method.eq(source_method)]
        if len(reference) != 5640 or reference.bacterium.nunique() != 60:
            raise ValueError("source primary score support differs")
        for host, rows in reference.groupby("bacterium", sort=True):
            vector = features.loc[("host_genome",host), CANONICAL_WORDS].to_numpy(float)
            predicted = score_feature_vector(bundles[method], vector)
            errors = [abs(predicted[r.phage]-r.score) for r in rows.itertuples()]
            if max(errors) > 1e-12:
                raise ValueError("research scorer differs from the frozen study")
            maximum_error = max(maximum_error, max(errors)); checked += len(rows)
    fasta_check = None
    if fasta is not None:
        records, maximum_fasta_error = {}, 0.
        for method, source_method in (("host-4mer", "host_4mer"), ("prevalence", "phage_prevalence")):
            result = rank_fasta(fasta, method=method, top_k=5, model_root=root, host_species="Escherichia coli")
            if result["status"] != "scored" or canonical_content_sha256(result) != result["content_sha256"]:
                raise ValueError("FASTA ranking is unavailable or not hash-bound")
            matches = result["applicability"]["known_cohort_matches"]
            if len(matches) != 1 or matches[0]["training_member"] or not matches[0]["exact_source_bytes"]:
                raise ValueError("FASTA QA requires exact source bytes for one development holdout host")
            host = matches[0]["host_id"]
            reference = scores[scores.method.eq(source_method)&scores.bacterium.eq(host)].set_index("phage").score
            if len(reference) != 94:
                raise ValueError("FASTA host does not belong to the primary held-out panel")
            difference = max(abs(r["score"]-reference[r["phage_id"]]) for r in result["candidates"])
            if difference > 1e-12:
                raise ValueError("full-FASTA ranking differs from study scores")
            maximum_fasta_error = max(maximum_fasta_error, difference)
            records[method] = result["content_sha256"]
        fasta_check = dict(host_id=host, source_sequence_sha256=sha256_file(fasta),
                           assembly_sha256=result["query"]["assembly_sha256"], scored_candidate_checks=188,
                           maximum_absolute_error=maximum_fasta_error, result_content_sha256=records,
                           scope="already inspected development holdout; software conformance, not new biological validation")
    return dict(status="verified_software_conformance", registry_sha256=REGISTRY_SHA256,
                verified_training_observations=10957, methods=["host-4mer","prevalence"],
                held_out_hosts=60, cached_feature_score_checks=checked, maximum_cached_score_error=maximum_error,
                full_fasta=fasta_check, verifier_sha256=sha256_file(Path(__file__)),
                scorer_sha256=sha256_file(root / "src/phage_cocktail/research_ranker.py"), new_model_fits=0,
                independent_release_validation=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("output must be new")
    result = verify(args.root, fasta=args.fasta)
    atomic_write_json(result,args.out)
    print(result)


if __name__ == "__main__":
    main()
