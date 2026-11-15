#!/usr/bin/env python3
"""Independent fixed-panel metrics, evidence and paired group-bootstrap checks."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import statistics
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from verify_admission import ROOT, STUDY, read, require, sha


def close(actual, expected):
    require(abs(float(actual) - float(expected)) < 1e-12, f"numeric mismatch: {actual} != {expected}")


def hit(scores, labels, k):
    cutoff = sorted(scores.values(), reverse=True)[k - 1]
    certain = [p for p in scores if scores[p] > cutoff]
    tied = [p for p in scores if scores[p] == cutoff]
    if any(labels[p] for p in certain):
        return 1.
    n, positives, take = len(tied), sum(labels[p] for p in tied), k - len(certain)
    return 1. if n - positives < take else 1 - math.comb(n - positives, take) / math.comb(n, take)


def interval(difference, group):
    members = defaultdict(list)
    for h in sorted(difference):
        members[group[h]].append(h)
    names = sorted(members)
    if len(names) < 2:
        return None, None
    rng, values = random.Random(20261115), []
    for _ in range(2000):
        expanded = [h for _ in names for h in members[names[rng.randrange(len(names))]]]
        values.append(statistics.mean(difference[h] for h in expanded))
    values.sort()
    def percentile(q):
        index = q * (len(values) - 1)
        low, high = math.floor(index), math.ceil(index)
        return values[low] + (values[high] - values[low]) * (index - low)
    return percentile(.025), percentile(.975)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, default=ROOT / "data/raw/picard_extension_sources")
    p.add_argument("--admission", type=Path, default=ROOT / "data/studies/results/picard_extension.v1/admission")
    p.add_argument("--evaluation", type=Path, required=True)
    args = p.parse_args()
    plan = json.loads((STUDY / "protocol.json").read_text())
    lock = json.loads((STUDY / "evaluation_implementation.json").read_text())
    for path, digest in lock["files_sha256"].items():
        require(sha(ROOT / path) == digest, "frozen evaluation code/input changed")
    freeze = json.loads((STUDY / "cohort_freeze.json").read_text())
    for path, digest in freeze["admission_output_sha256"].items():
        require(sha(args.admission / path) == digest, "admission changed after freezing")
    for name, spec in plan["resources"].items():
        if "sha256" in spec:
            require(sha(args.sources / (name + (".tsv" if name == "core_distances" else ".csv"))) == spec["sha256"], "source input changed")
    matrix = read(args.sources / "source_mlc.csv", ";")
    source = {r["bacteria"]: r for r in matrix}
    coordinates = {r["bacteria"]: i for i, r in enumerate(matrix, 2)}
    columns = {p: i for i, p in enumerate(matrix[0], 1)}
    phages = sorted(set(matrix[0]) - {"bacteria", "T4LD", "T7_Portugal"})
    roster = read(args.admission / "admitted_hosts.csv")
    groups = {r["source_host"]: r["component"] for r in roster}
    require(len(roster) == freeze["eligible_hosts"] and len(set(groups.values())) == freeze["eligible_groups"], "target support differs")
    with gzip.open(ROOT / "data/studies/results/paired_interaction.v1/observed_pairs.csv.gz", "rt") as f:
        observed = list(csv.DictReader(f))
    require(len(observed) == 16597, "old observation support differs")
    for r in observed:
        close(r["source_mlc"], source[r["source_bacterium"]][r["phage"]])
        require(int(r["infection"]) == int(float(r["source_mlc"]) > 0), "old binary label differs")
    splits = read(ROOT / "data/studies/results/genophi_comparator.v1/scenario_splits.csv")
    policies, training_statistics = {}, {}
    for scope, scenario, pair_key in (("primary", "grouped_host", "primary_pair"), ("secondary", "phylogroup:B2", "secondary_pair")):
        train = {r["bacterium"] for r in splits if r["scenario"] == scenario and r["split"] == "train"}
        training = [r for r in observed if r["bacterium"] in train]
        by_phage = {p: [r for r in training if r["phage"] == p] for p in phages}
        policies[scope + "_prevalence"] = {p: float(Fraction(sum(int(r["infection"]) for r in by_phage[p]) + 1, len(by_phage[p]) + 2)) for p in phages}
        pair = plan["evaluation"][pair_key]
        policies[scope + "_complement"] = {p: float(-i) for i, p in enumerate(pair)}
        for p in phages:
            training_statistics[scope, p] = len(train), sum(int(r["infection"]) for r in by_phage[p]), len(by_phage[p])
    expected, expected_evidence = {}, {}
    for host in roster:
        h = host["source_host"]
        require(all(source[h][p].strip().lower() not in ("", "nan", "na") for p in phages), "unknown target entered evaluation")
        y = {p: int(float(source[h][p]) > 0) for p in phages}
        for method, scores in policies.items():
            expected[h, method] = {"hit_at_1": hit(scores, y, 1), "hit_at_2": hit(scores, y, 2),
                "tests_consumed_capped_2": 2 - hit(scores, y, 1), "positive_candidates": sum(y.values())}
            threshold = sorted(scores.values(), reverse=True)[1]
            for p in scores:
                if scores[p] >= threshold:
                    expected_evidence[h, method, p] = 1 + sum(v > scores[p] for v in scores.values()), sum(v >= scores[p] for v in scores.values())
    records = read(args.evaluation / "per_host_metrics.csv")
    require(len(records) == len(expected) and {(r["source_host"], r["method"]) for r in records} == set(expected), "host metric support differs")
    for r in records:
        h, method = r["source_host"], r["method"]
        require(r["component"] == groups[h], "target group differs")
        for field, value in expected[h, method].items():
            close(r[field], value)
    training_rows = read(args.evaluation / "training_prevalence.csv")
    require(len(training_rows) == len(training_statistics), "training prior support differs")
    for r in training_rows:
        h, positives, known = training_statistics[r["scope"], r["phage"]]
        require((int(r["training_hosts"]), int(r["training_positive"]), int(r["training_observed"])) == (h, positives, known), "training prior label scope differs")
        close(r["score"], policies[r["scope"] + "_prevalence"][r["phage"]])
    stratum_hosts = {s: [r["source_host"] for r in roster if s == "all" or (r["source_phylogroup"] == "B2") == (s == "B2")] for s in ("all", "B2", "non_B2")}
    summaries = read(args.evaluation / "summary.csv")
    require(len(summaries) == 12, "summary support differs")
    for r in summaries:
        hosts, method = stratum_hosts[r["stratum"]], r["method"]
        require(int(r["hosts"]) == len(hosts) and int(r["groups"]) == len({groups[h] for h in hosts}), "summary denominator differs")
        require(int(r["all_source_negative_hosts"]) == sum(expected[h, method]["positive_candidates"] == 0 for h in hosts), "all-negative support differs")
        for field in ("hit_at_1", "hit_at_2", "tests_consumed_capped_2"):
            close(r[field], statistics.mean(expected[h, method][field] for h in hosts))
    transitions = read(args.evaluation / "host_transitions.csv")
    require(len(transitions) == len(roster) * 2 and len({(r["source_host"], r["scope"]) for r in transitions}) == len(transitions), "transition support differs")
    for r in transitions:
        h, scope = r["source_host"], r["scope"]
        g, ref = expected[h, scope + "_complement"]["hit_at_2"], expected[h, scope + "_prevalence"]["hit_at_2"]
        close(r["expected_rescued"], g * (1 - ref))
        close(r["expected_newly_missed"], (1 - g) * ref)
        classification = "rescue" if g > ref else "new_miss" if g < ref else "both_hit" if g else "both_miss"
        require(r["transition"] == classification, "transition classification differs")
    comparisons = read(args.evaluation / "paired_comparisons.csv")
    require(len(comparisons) == 6, "paired comparison support differs")
    for r in comparisons:
        scope, hosts = r["scope"], stratum_hosts[r["stratum"]]
        rescued = [h for h in hosts if expected[h, scope + "_complement"]["hit_at_2"] > expected[h, scope + "_prevalence"]["hit_at_2"]]
        lost = [h for h in hosts if expected[h, scope + "_complement"]["hit_at_2"] < expected[h, scope + "_prevalence"]["hit_at_2"]]
        require(r["rescued_host_ids"] == ";".join(rescued) and r["newly_missed_host_ids"] == ";".join(lost), "named paired cases differ")
        require(int(r["rescued_hosts"]) == len(rescued) and int(r["newly_missed_hosts"]) == len(lost), "paired counts differ")
        for field, label in (("hit_at_2", "hit"), ("tests_consumed_capped_2", "cost")):
            differences = {h: expected[h, scope + "_complement"][field] - expected[h, scope + "_prevalence"][field] for h in hosts}
            low, high = interval(differences, groups)
            close(r[label + "_delta_ci_low"], low)
            close(r[label + "_delta_ci_high"], high)
            close(r["hit_at_2_delta" if label == "hit" else "tests_consumed_delta"], statistics.mean(differences.values()))
    evidence = read(args.evaluation / "top2_source_evidence.csv")
    require(len(evidence) == len(expected_evidence) and {(r["source_host"], r["method"], r["phage"]) for r in evidence} == set(expected_evidence), "source evidence support differs")
    for r in evidence:
        h, method, p = r["source_host"], r["method"], r["phage"]
        low, high = expected_evidence[h, method, p]
        require(int(r["rank_min"]) == low and int(r["rank_max"]) == high, "evidence ranks differ")
        close(r["top2_inclusion_probability"], min(3 - low, high - low + 1) / (high - low + 1))
        close(r["source_mlc"], source[h][p])
        require(int(r["source_positive"]) == int(float(source[h][p]) > 0), "source endpoint differs")
        require(int(r["source_matrix_row"]) == coordinates[h] and int(r["source_matrix_column"]) == columns[p], "source coordinates differ")
        require(r["matrix_url"].endswith(f"/interaction_matrix.csv#L{coordinates[h]}") and r["raw_annotations_independently_checked"] == "False", "source evidence scope differs")
    changed = {r["source_host"] for r in transitions if r["transition"] in ("rescue", "new_miss")}
    require(read(args.evaluation / "changed_host_evidence.csv") == [r for r in evidence if r["source_host"] in changed], "changed-host evidence differs")
    marginal = read(args.evaluation / "candidate_marginal_coverage.csv")
    require(len(marginal) == 6, "marginal contribution support differs")
    for r in marginal:
        anchor, candidate = plan["evaluation"][r["scope"] + "_pair"]
        extra = [h for h in stratum_hosts[r["stratum"]] if float(source[h][anchor]) == 0 and float(source[h][candidate]) > 0]
        require(r["extra_host_ids"] == ";".join(extra) and int(r["extra_hosts_beyond_anchor"]) == len(extra), "candidate contribution differs")
    receipt = json.loads((args.evaluation / "QA.json").read_text())
    for name, digest in receipt["output_sha256"].items():
        require(sha(args.evaluation / name) == digest, "evaluation output checksum differs")
    result = {"status": "verified", "imports_evaluation_runner": False,
        "target_source_pairs_checked": len(roster) * 94, "old_observed_pairs_checked": 16597,
        "host_metric_rows_checked": len(records), "paired_transitions_checked": len(transitions),
        "evidence_rows_checked": len(evidence), "paired_interval_draws_checked": 24000,
        "verifier_sha256": sha(Path(__file__)), "checked_output_sha256": {p.name: sha(p) for p in sorted(args.evaluation.iterdir()) if p.suffix == ".csv" or p.name == "QA.json"}}
    (args.evaluation / "INDEPENDENT_QA.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "checked_output_sha256"}, indent=2))


if __name__ == "__main__":
    main()
