#!/usr/bin/env python3
"""Independently check saved BASEL features, cosine scores and label bounds."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import canonical_json, read_json_object, sha256_file
from phage_cocktail.basel_controls import validate_feature
from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock


def close(actual, expected):
    if not math.isfinite(float(actual)) or not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"numerical mismatch: {actual} versus {expected}")


def independent_metrics(scores, labels, marked, k, unknown_value):
    """Multiply exact hypergeometric no-hit probabilities over tied blocks.

    Separate from the production helper's cumulative-hit implementation.
    """
    order = sorted(scores, key=lambda p: (-scores[p], p))
    failure, marked_slots, start = 1., 0., 0
    while start < k:
        end = start + 1
        while end < len(order) and scores[order[end]] == scores[order[start]]:
            end += 1
        block, slots = order[start:end], min(k - start, end - start)
        positive = sum(labels.get(p, unknown_value) for p in block)
        failure *= math.comb(len(block) - positive, slots) / math.comb(len(block), slots)
        marked_slots += slots * sum(p in marked for p in block) / len(block)
        start = end
    return 1. - failure, marked_slots / k


def verify(bundle, data_root):
    bundle, data_root = Path(bundle), Path(data_root)
    report = read_json_object(bundle / "report.json")
    protocol_path = ROOT / "data/studies/basel_ecor_evaluation_protocol.v2.json"
    protocol = read_json_object(protocol_path)
    if sha256_file(protocol_path) != report["provenance"]["protocol_sha256"]:
        raise ValueError("protocol provenance differs")
    for name, digest in protocol["files_sha256"].items():
        if sha256_file(ROOT / name) != digest:
            raise ValueError(f"frozen input changed: {name}")
    for name, digest in report["artifacts_sha256"].items():
        path = (bundle / name).resolve()
        if not path.is_relative_to(bundle.resolve()) or sha256_file(path) != digest:
            raise ValueError(f"saved output changed: {name}")
    for name, digest in report["provenance"]["code_sha256"].items():
        if sha256_file(ROOT / name) != digest:
            raise ValueError(f"execution code changed: {name}")
    source_path = data_root / "source.csv"
    source = inspect_basel_ecor(source_path.read_bytes(), load_lock(ROOT / "data/studies/basel_ecor_intake.v1.json"))
    if report["source_counts"] != source["source_summary"]:
        raise ValueError("reported source support differs")
    with source_path.open(newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    labels = {h: {} for h in source["host_ids"]}
    for row in source_rows:
        labels[row["strain"]][row["phage"]] = int(row["interaction"])
    spec = sha256(canonical_json(protocol["representation"]).encode()).hexdigest()
    if spec != report["representation_sha256"]:
        raise ValueError("reported representation differs")
    columns = protocol["representation"]["ordered_features"]
    key = ["entity_kind", "source_build", "entity_id"]
    manifest = pd.read_csv(bundle / "feature_manifest.csv").set_index(key, verify_integrity=True)
    features = pd.read_csv(bundle / "sequence_features.csv", float_precision="round_trip").set_index(key, verify_integrity=True)
    if set(features.index) != set(manifest.index) or len(features) != 102 or features.columns.tolist() != columns:
        raise ValueError("feature identities or ordered columns differ")
    audited_hosts = pd.read_csv(ROOT / "data/studies/results/basel_ecor_sequence_audit.v1/host_assemblies.csv").set_index("host")
    audited_phages = pd.read_csv(ROOT / "data/studies/results/basel_ecor_sequence_audit.v1/basel_accessions.csv").set_index("phage")
    for identity, row in manifest.iterrows():
        kind, build, name = identity
        if kind == "phage":
            expected, raw_hash = audited_phages.loc[name, ["assembly_sha256", "raw_fasta_sha256"]]
        elif build == "published_wgs_primary":
            expected, raw_hash = audited_hosts.loc[name, ["published_assembly_sha256", "published_raw_gzip_sha256"]]
        else:
            expected, raw_hash = audited_hosts.loc[name, ["training_assembly_sha256", "training_raw_sha256"]]
        validate_feature(features.loc[identity].to_numpy(float), row.to_dict(), expected, spec)
        if row.source_file_sha256 != raw_hash:
            raise ValueError("feature source file does not match the audited genome")
    per_host = pd.read_csv(bundle / "per_host_bounds.csv", float_precision="round_trip")
    if len(per_host) != 300 or per_host.duplicated(["method", "bacterium", "k"]).any():
        raise ValueError("host metric table must contain exactly three complete 25-host panels at four K values")
    rankings = pd.read_csv(bundle / "deterministic_top2.csv")
    near, checked, independent_top = set(report["training_neighbours_i90_retained"]), 0, {}
    audit = read_json_object(ROOT / "data/studies/results/basel_ecor_sequence_audit.v1/report.json")
    if near != set(next(s["candidates_in_training_components"] for s in audit["phage_overlap"] if s["identity"] == .9)):
        raise ValueError("near-training flags differ from frozen genomic groups")
    for method, summary in report["methods"].items():
        frame = pd.read_csv(bundle / (method + "_scores.csv"), float_precision="round_trip")
        pair_index = pd.MultiIndex.from_frame(frame[["bacterium", "phage"]])
        expected_pairs = set(pd.MultiIndex.from_product([source["host_ids"], source["phage_ids"]]))
        if pair_index.has_duplicates or set(pair_index) != expected_pairs:
            raise ValueError("score panel membership differs")
        if method != "uniform_random":
            build = method.removeprefix("cosine_")
            for row in frame.itertuples():
                a = features.loc[("host", build, row.bacterium)].tolist()
                b = features.loc[("phage", "basel_accession", row.phage)].tolist()
                expected_score = math.fsum(x*y for x, y in zip(a, b)) / math.sqrt(
                    math.fsum(x*x for x in a) * math.fsum(x*x for x in b))
                close(row.score, expected_score)
                checked += 1
        elif not frame.score.eq(0.).all():
            raise ValueError("random comparator must have equal zero scores")
        metric_rows = []
        for host, group in frame.groupby("bacterium", sort=True):
            scores = dict(zip(group.phage, group.score))
            unknown = set(source["phage_ids"]) - set(labels[host])
            values = {j: [independent_metrics(scores, labels[host], unknown, j, v)[0] for v in (0, 1)]
                      for j in range(1, 6)}
            top = sorted(scores, key=lambda p: (-scores[p], p))[:2]
            independent_top[method, host] = set(top)
            saved_top = rankings[rankings.method.eq(method) & rankings.bacterium.eq(host)]
            if len(saved_top) != 1 or saved_top.iloc[0][["first", "second"]].tolist() != top:
                raise ValueError("named top-2 set differs from the fixed score order")
            for k in (1, 2, 3, 5):
                expected = dict(hit_lower=values[k][0], hit_upper=values[k][1],
                                tests_lower=k - sum(values[j][1] for j in range(1, k)),
                                tests_upper=k - sum(values[j][0] for j in range(1, k)),
                                unknown_slot_fraction=independent_metrics(scores, labels[host], unknown, k, 0)[1],
                                near_training_i90_slot_fraction=independent_metrics(scores, labels[host], near, k, 0)[1])
                saved = per_host[per_host.method.eq(method) & per_host.bacterium.eq(host) & per_host.k.eq(k)]
                if len(saved) != 1:
                    raise ValueError("missing or duplicated host metric")
                for metric, value in expected.items():
                    close(saved.iloc[0][metric], value)
                metric_rows.append(dict(k=k, **expected))
        for k, rows in pd.DataFrame(metric_rows).groupby("k"):
            for metric, value in rows.drop(columns="k").mean().items():
                close(summary["macro"][str(k)][metric], value)
        if method != "uniform_random":
            for k, values in summary["macro"].items():
                random = report["methods"]["uniform_random"]["macro"][k]
                saved = summary["hit_minus_random_conservative_outer_bounds"][k]
                close(saved["lower"], values["hit_lower"] - random["hit_upper"])
                close(saved["upper"], values["hit_upper"] - random["hit_lower"])
    changed = sorted(h for h in source["host_ids"] if independent_top["cosine_published_wgs_primary", h]
                     != independent_top["cosine_training_assembly_sensitivity", h])
    if changed != sorted(report["assembly_sensitivity"]["changed_deterministic_top2_hosts"]):
        raise ValueError("assembly sensitivity membership differs")
    return dict(status="verified", feature_vectors=102, independently_checked_cosine_scores=checked,
                independently_checked_host_K_rows=len(per_host), retained_unknown_pairs=source["source_summary"]["unreported_pairs"],
                checks=["published input and execution code hashes", "102 serialized feature vectors and audited sequence identities",
                        "2600 scalar-fsum cosine recomputations", "complete panels and all 25 host denominators",
                        "independent hypergeometric hit/cost bounds and slot exposures", "deterministic ranking and assembly sensitivity"],
                published_output_sha256={p.name: sha256_file(p) for p in sorted(bundle.iterdir()) if p.is_file()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.bundle, args.data_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print("Verified", result["feature_vectors"], "feature vectors and", result["independently_checked_cosine_scores"], "scores")


if __name__ == "__main__":
    main()
