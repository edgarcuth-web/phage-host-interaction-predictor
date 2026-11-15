#!/usr/bin/env python3
"""Independently verify the fixed-panel study without importing its runner."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path, delimiter=","):
    with (gzip.open if str(path).endswith(".gz") else open)(path, "rt", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def require(value, message):
    if not value:
        raise ValueError(message)


def close(actual, expected):
    require(math.isfinite(float(actual)) and abs(float(actual) - expected) < 1e-12,
            f"numeric mismatch: {actual} != {expected}")


def cutoff_hit(scores, labels, k, unknown):
    # Only the cutoff block can be partially selected; all blocks above it are certain.
    threshold = sorted(scores.values(), reverse=True)[k - 1]
    above = [p for p in scores if scores[p] > threshold]
    tied = [p for p in scores if scores[p] == threshold]
    positive = lambda p: labels[p] == 1 or (unknown and labels[p] is None)
    if any(positive(p) for p in above):
        return 1.
    n, count, take = len(tied), sum(positive(p) for p in tied), k - len(above)
    return 1. if n - count < take else 1 - math.comb(n - count, take) / math.comb(n, take)


def percentile(values, q):
    values.sort()
    lower, fraction = divmod((len(values) - 1) * q, 1)
    i = int(lower)
    return values[i] * (1 - fraction) + values[min(i + 1, len(values) - 1)] * fraction


def interval(differences, host_groups, draws, seed):
    members = defaultdict(list)
    for h in sorted(differences):
        members[host_groups[h]].append(h)
    names = sorted(members)
    rng, means = random.Random(seed), []
    for _ in range(draws):
        expanded = [h for _ in names for h in members[names[rng.randrange(len(names))]]]
        means.append(statistics.mean(differences[h] for h in expanded))
    return percentile(means.copy(), .025), percentile(means, .975)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--source-matrix", type=Path, default=ROOT / "data/raw/complementary_panel/mlc.csv")
    args = parser.parse_args()
    bundle = args.bundle
    plan = json.loads((ROOT / "studies/complementary_panel/protocol.json").read_text())
    lock = json.loads((ROOT / "studies/complementary_panel/implementation.json").read_text())
    for path, digest in {**plan["inputs_sha256"], **lock["files_sha256"]}.items():
        require(sha(ROOT / path) == digest, f"input/code changed: {path}")
    require(sha(args.source_matrix) == plan["source_matrix"]["sha256"], "source matrix changed")
    base = ROOT / "data/studies/results"
    observed = rows(base / "paired_interaction.v1/observed_pairs.csv.gz")
    missing = rows(base / "paired_interaction.v1/excluded_unknown_pairs.csv")
    labels = {(r["bacterium"], r["phage"]): int(r["infection"]) for r in observed}
    labels.update({(r["bacterium"], r["phage"]): None for r in missing})
    require(len(labels) == 177 * 94 and len(observed) == 16597 and len(missing) == 41, "cohort differs")
    phages = sorted({p for _, p in labels})
    original = rows(args.source_matrix, ";")
    matrix = {r["bacteria"]: (line, r) for line, r in enumerate(original, 2)}
    column = {name: index for index, name in enumerate(original[0], 1)}
    aliases = {r["bacterium"]: r["source_bacterium"] for r in observed + missing}
    for r in observed:
        value = matrix[r["source_bacterium"]][1][r["phage"]]
        close(value, float(r["source_mlc"]))
        require(int(float(value) > 0) == int(r["infection"]), "source label differs")
    for r in missing:
        require(matrix[r["source_bacterium"]][1][r["phage"]].lower() in ("", "nan", "na"), "unknown was observed")
    splits = rows(base / "genophi_comparator.v1/scenario_splits.csv")
    exports = defaultdict(dict)
    for r in rows(base / "new_host_retrieval.v1/scores.csv.gz"):
        if r["method"] in plan["comparators"]:
            key = r["panel"], r["partition"], r["method"], r["bacterium"]
            require(r["phage"] not in exports[key], "duplicate score")
            exports[key][r["phage"]] = float(r["score"])
    panel_rows = rows(bundle / "selected_panels.csv")
    selected = {panel: [r["phage"] for r in sorted(panel_rows, key=lambda r: int(r["step"])) if r["panel"] == panel]
                for panel in plan["panels"]}
    stats, groups = {}, {}
    for panel, config in plan["panels"].items():
        local = [r for r in splits if r["scenario"] == panel]
        train = sorted(r["bacterium"] for r in local if r["split"] == "train")
        groups[panel] = {r["bacterium"]: r["host_group"] for r in local}
        require(len(train) == config["train_hosts"] and len(selected[panel]) == len(set(selected[panel])) == 2, "panel support differs")
        positives = {p: {h for h in train if labels[h, p] == 1} for p in phages}
        known = {p: {h for h in train if labels[h, p] is not None} for p in phages}
        rates = {p: Fraction(len(positives[p]) + 1, len(known[p]) + 2) for p in phages}
        stats[panel] = train, positives, known, rates
        covered = set()
        for step, actual in enumerate(selected[panel], 1):
            allowed = set(phages) - set(selected[panel][:step - 1])
            values = {p: len(positives[p] - covered) for p in allowed}
            optimum = max(values.values())
            ties = sorted(p for p in allowed if values[p] == optimum)
            require(actual == sorted(ties, key=lambda p: (-rates[p], p))[0], "train-only optimum differs")
            row = next(r for r in panel_rows if r["panel"] == panel and int(r["step"]) == step)
            require(row["max_gain_ties"].split(";") == ties, "selection ties differ")
            require(row["new_training_hosts"].split(";") == sorted(positives[actual] - covered), "training coverage differs")
            require(int(row["marginal_training_hosts"]) == optimum, "training gain differs")
            close(row["training_smoothed_prevalence"], float(rates[actual]))
            require((row["below_training_median"] == "True") == (rates[actual] < statistics.median(rates.values())), "prevalence category differs")
            covered |= positives[actual]
        # Independently expand the group bootstrap into host multiplicities.
        names = sorted({groups[panel][h] for h in train})
        members = {g: [h for h in train if groups[panel][h] == g] for g in names}
        rng = random.Random(plan["stability"]["seed"])
        counts, first_counts = Counter(), Counter()
        for _ in range(plan["stability"]["training_group_bootstrap_draws"]):
            sample = [h for _ in names for h in members[names[rng.randrange(len(names))]]]
            multiplicity = Counter(sample)
            weighted_rate = {p: Fraction(sum(multiplicity[h] for h in positives[p]) + 1,
                                         sum(multiplicity[h] for h in known[p]) + 2) for p in phages}
            used, coverage = [], set()
            for step in range(2):
                chosen = sorted((p for p in phages if p not in used), key=lambda p:
                    (-sum(multiplicity[h] for h in positives[p] - coverage), -weighted_rate[p], p))[0]
                counts[chosen] += 1
                first_counts[chosen] += int(step == 0)
                coverage |= positives[chosen]
                used.append(chosen)
        retained = [r for r in rows(bundle / "training_bootstrap_selection.csv") if r["panel"] == panel]
        require(len(retained) == 94 and {r["phage"] for r in retained} == set(phages), "bootstrap candidate support differs")
        for row in retained:
            p = row["phage"]
            require(int(row["selected_draws"]) == counts[p] and int(row["first_draws"]) == first_counts[p], "training bootstrap differs")
            require(int(row["second_draws"]) == counts[p] - first_counts[p], "bootstrap second-position count differs")
            close(row["selection_frequency"], counts[p] / plan["stability"]["training_group_bootstrap_draws"])
    metrics = rows(bundle / "per_host_metrics.csv")
    require(len(metrics) == 488, "target metric support differs")
    by_key, expected_top2 = {}, {}
    for r in metrics:
        panel, partition, method, host = r["panel"], r["partition"], r["method"], r["bacterium"]
        key = panel, partition, method, host
        require(key not in by_key, "duplicated metric row")
        by_key[key] = r
        actual_partition = next(x["split"] for x in splits if x["scenario"] == panel and x["bacterium"] == host)
        require(actual_partition == partition and partition != "train", "target membership differs")
        y = {p: labels[host, p] for p in phages}
        s = exports[key] if method != "greedy_coverage" else {p: float(-i) for i, p in enumerate(selected[panel])}
        if method != "greedy_coverage":
            require(set(s) == set(phages), "score support differs")
        threshold = sorted(s.values(), reverse=True)[1]
        for p in s:
            if s[p] >= threshold:
                rank_min = 1 + sum(v > s[p] for v in s.values())
                rank_max = sum(v >= s[p] for v in s.values())
                expected_top2[panel, partition, host, method, p] = rank_min, rank_max
        if method == "phage_prevalence":
            for p in phages:
                close(s[p], float(stats[panel][3][p]))
        complete = all(v is not None for v in y.values())
        require((r["complete_panel"] == "True") == complete, "unknown host entered primary metrics")
        low1, low2 = cutoff_hit(s, y, 1, False), cutoff_hit(s, y, 2, False)
        high1, high2 = cutoff_hit(s, y, 1, True), cutoff_hit(s, y, 2, True)
        for field, v in {"hit_at_1_lower": low1, "hit_at_1_upper": high1, "hit_at_2_lower": low2,
                         "hit_at_2_upper": high2, "tests_consumed_lower": 2 - high1, "tests_consumed_upper": 2 - low1}.items():
            close(r[field], v)
        for field, v in {"hit_at_1": low1, "hit_at_2": low2, "tests_consumed_capped_2": 2 - low1}.items():
            if complete:
                close(r[field], v)
            else:
                require(r[field] == "", "incomplete host entered primary result")
    transitions = rows(bundle / "host_transitions.csv")
    require(len(transitions) == 363, "paired host support differs")
    for r in transitions:
        key = r["panel"], r["partition"], r["bacterium"]
        a, b = by_key[key[0], key[1], "greedy_coverage", key[2]], by_key[key[0], key[1], r["comparator"], key[2]]
        g, ref = float(a["hit_at_2"]), float(b["hit_at_2"])
        close(r["expected_rescued"], g * (1 - ref))
        close(r["expected_newly_missed"], (1 - g) * ref)
        expected = "rescue" if g == 1 and ref == 0 else "new_miss" if g == 0 and ref == 1 else "tie_dependent" if g != ref else "both_hit" if g else "both_miss"
        require(r["transition"] == expected, "paired transition differs")
    for r in rows(bundle / "summary.csv"):
        local = [x for x in metrics if all(x[k] == r[k] for k in ("panel", "partition", "method")) and x["complete_panel"] == "True"]
        require(len(local) == int(r["complete_hosts"]), "summary denominator differs")
        for field in ("hit_at_1", "hit_at_2", "tests_consumed_capped_2"):
            close(r[field], statistics.mean(float(x[field]) for x in local))
    for r in rows(bundle / "paired_comparisons.csv"):
        panel, partition, method = r["panel"], r["partition"], r["comparator"]
        local = [x for x in transitions if (x["panel"], x["partition"], x["comparator"]) == (panel, partition, method)]
        close(r["expected_rescued_hosts"], sum(float(x["expected_rescued"]) for x in local))
        close(r["expected_newly_missed_hosts"], sum(float(x["expected_newly_missed"]) for x in local))
        for metric, stem in (("hit_at_2", "hit"), ("tests_consumed_capped_2", "cost")):
            difference = {x["bacterium"]: float(by_key[panel, partition, "greedy_coverage", x["bacterium"]][metric]) -
                          float(by_key[panel, partition, method, x["bacterium"]][metric]) for x in local}
            low, high = interval(difference, groups[panel], 2000, plan["evaluation"]["seed"])
            close(r[stem + "_delta_ci_low"], low)
            close(r[stem + "_delta_ci_high"], high)
            close(r["hit_at_2_delta" if stem == "hit" else "tests_consumed_delta"], statistics.mean(difference.values()))
    candidate_rows = rows(bundle / "candidate_marginal_coverage.csv")
    require(len(candidate_rows) == len({(r["panel"], r["partition"], r["phage"]) for r in candidate_rows}) == 282,
            "candidate diagnostic support differs")
    for r in candidate_rows:
        panel, partition, p = r["panel"], r["partition"], r["phage"]
        anchor = selected[panel][0]
        train, positive, known, rates = stats[panel]
        close(r["training_smoothed_prevalence"], float(rates[p]))
        require(int(r["training_observed"]) == len(known[p]) and int(r["training_positive"]) == len(positive[p]), "candidate training statistics differ")
        require(int(r["extra_known_positive_training_hosts"]) == len(positive[p] - positive[anchor]), "candidate training contribution differs")
        require((r["below_training_median"] == "True") == (rates[p] < statistics.median(rates.values())), "candidate prevalence class differs")
        target = [x["bacterium"] for x in splits if x["scenario"] == panel and x["split"] == partition]
        extra = sorted(h for h in target if labels[h, anchor] == 0 and labels[h, p] == 1)
        require(r["target_extra_host_ids"] == ";".join(extra), "candidate target contribution differs")
        require(int(r["target_extra_hosts_anchor_known_negative"]) == len(extra), "candidate target count differs")
        require(int(r["target_extra_groups_anchor_known_negative"]) == len({groups[panel][h] for h in extra}), "candidate target groups differ")
        unresolved = sorted(h for h in target if labels[h, anchor] is None and labels[h, p] == 1)
        require(r["target_anchor_unknown_positive_host_ids"] == ";".join(unresolved), "candidate unknown-anchor support differs")
    evidence = rows(bundle / "top2_source_evidence.csv")
    evidence_key = lambda r: (r["panel"], r["partition"], r["bacterium"], r["method"], r["phage"])
    require(len(evidence) == len({evidence_key(r) for r in evidence}) == len(expected_top2)
            and {evidence_key(r) for r in evidence} == set(expected_top2), "top-two source evidence support differs")
    for r in evidence:
        host, p = r["bacterium"], r["phage"]
        line, source = matrix[aliases[host]]
        require(int(r["source_matrix_row"]) == line and int(r["source_matrix_column"]) == column[p], "source coordinate differs")
        require(r["matrix_url"].endswith(f"/interaction_matrix.csv#L{line}"), "source URL differs")
        low, high = expected_top2[evidence_key(r)]
        require(int(r["rank_min"]) == low and int(r["rank_max"]) == high, "source evidence ranks differ")
        close(r["top2_inclusion_probability"], min(3 - low, high - low + 1) / (high - low + 1))
        if labels[host, p] is None:
            require(r["source_label"] == "" and r["source_mlc"] == "", "unknown evidence relabelled")
        else:
            close(r["source_mlc"], float(source[p]))
            require(int(r["source_label"]) == labels[host, p], "source evidence label differs")
    require(rows(bundle / "focus_source_evidence.csv") == [r for r in evidence if r["bacterium"] in ("NILS53", "NILS78", "AL505")],
            "focus cases differ")
    changed = {(r["panel"], r["partition"], r["bacterium"]) for r in transitions
               if r["transition"] in ("rescue", "new_miss", "tie_dependent")}
    require(rows(bundle / "changed_host_evidence.csv") == [r for r in evidence
            if (r["panel"], r["partition"], r["bacterium"]) in changed], "changed-host evidence differs")
    result = {"status": "verified", "imports_runner": False, "source_pairs_checked": len(labels),
        "host_metric_rows_checked": len(metrics), "paired_transitions_checked": len(transitions),
        "training_group_bootstrap_draws_checked": 4000, "paired_interval_draws_checked": 36000,
        "source_evidence_rows_checked": len(evidence), "verifier_sha256": sha(Path(__file__)),
        "checked_output_sha256": {p.name: sha(p) for p in sorted(bundle.iterdir()) if p.suffix == ".csv" or p.name == "QA.json"}}
    (bundle / "INDEPENDENT_QA.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "checked_output_sha256"}, indent=2))


if __name__ == "__main__":
    main()
