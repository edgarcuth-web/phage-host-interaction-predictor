#!/usr/bin/env python3
"""Select fixed two-phage screening panels using training labels only.

Python 3.11+ standard library; no model fitting or target-driven selection.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import platform
import random
import shutil
import statistics
import tempfile
import urllib.request
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROTOCOL = "studies/complementary_panel/protocol.json"
OBSERVED = "data/studies/results/paired_interaction.v1/observed_pairs.csv.gz"
UNKNOWN = "data/studies/results/paired_interaction.v1/excluded_unknown_pairs.csv"
SPLITS = "data/studies/results/genophi_comparator.v1/scenario_splits.csv"
SCORES = "data/studies/results/new_host_retrieval.v1/scores.csv.gz"
HISTORICAL = "data/studies/results/new_host_retrieval.v1/per_host_metrics.csv"
EVIDENCE = "data/studies/results/missed_candidate_audit.v1/source_evidence.csv"
FOCUS = ("NILS53", "NILS78", "AL505")
SELECTOR = "greedy_coverage"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path, delimiter=","):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"empty output without a schema: {path}")
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def check(condition, message):
    if not condition:
        raise ValueError(message)


def candidate_stats(training, phages, labels, weights=None):
    """Only known observations contribute to prevalence or coverage."""
    weights = {h: 1 for h in training} if weights is None else weights
    check(set(weights) <= set(training), "selection weights include a nontraining host")
    result = {}
    for phage in phages:
        positive = {h for h in training if labels[h, phage] == 1}
        n = sum(weights.get(h, 0) for h in training if labels[h, phage] is not None)
        p = sum(weights.get(h, 0) for h in positive)
        result[phage] = {"positive_hosts": positive, "positive": p, "observed": n,
                         "rate": Fraction(p + 1, n + 2)}
    return result


def select_panel(training, phages, labels, budget=2, weights=None):
    """Return an ordered panel and trace; target labels are never inspected."""
    check(0 < budget <= len(phages), "invalid panel budget")
    stats = candidate_stats(training, phages, labels, weights)
    weights = {h: 1 for h in training} if weights is None else weights
    selected, covered, trace = [], set(), []
    for step in range(1, budget + 1):
        gains = {p: sum(weights.get(h, 0) for h in s["positive_hosts"] - covered)
                 for p, s in stats.items() if p not in selected}
        maximum = max(gains.values())
        tied = sorted(p for p, gain in gains.items() if gain == maximum)
        chosen = min(tied, key=lambda p: (-stats[p]["rate"], p))
        added = stats[chosen]["positive_hosts"] - covered
        covered |= stats[chosen]["positive_hosts"]
        selected.append(chosen)
        trace.append({"step": step, "phage": chosen, "marginal_training_hosts": maximum,
                      "cumulative_training_hosts": sum(weights.get(h, 0) for h in covered),
                      "new_training_hosts": ";".join(sorted(h for h in added if weights.get(h, 0))),
                      "max_gain_ties": ";".join(tied)})
    return selected, trace, stats


def hit_probability(scores, labels, k, *, unknown_positive=False):
    """Exact rational Hit@K under uniform within-score-tie permutations."""
    blocks = defaultdict(list)
    for phage, score in scores.items():
        check(math.isfinite(score), "nonfinite score")
        blocks[score].append(phage)
    remaining, no_hit = k, Fraction(1)
    for score in sorted(blocks, reverse=True):
        members = blocks[score]
        take = min(remaining, len(members))
        positive = sum(labels[p] == 1 or (unknown_positive and labels[p] is None) for p in members)
        negative = len(members) - positive
        if negative < take:
            return Fraction(1)
        no_hit *= Fraction(math.comb(negative, take), math.comb(len(members), take))
        remaining -= take
        if remaining == 0:
            return 1 - no_hit
    raise ValueError("budget exceeds candidate support")


def quantile(values, q):
    values = sorted(values)
    index = (len(values) - 1) * q
    left = math.floor(index)
    return values[left] + (values[math.ceil(index)] - values[left]) * (index - left)


def paired_bootstrap(differences, groups, draws, seed):
    """Resample whole held-out groups; keep their hosts and multiplicity."""
    grouped = defaultdict(list)
    for host in sorted(differences):
        grouped[groups[host]].append(differences[host])
    names = sorted(grouped)
    check(len(names) >= 2, "bootstrap needs at least two held-out groups")
    rng, estimates = random.Random(seed), []
    for _ in range(draws):
        sample = [grouped[rng.choice(names)] for _ in names]
        estimates.append(sum(sum(v) for v in sample) / sum(len(v) for v in sample))
    return [quantile(estimates, .025), quantile(estimates, .975)]


def training_stability(training, phages, labels, groups, draws, seed):
    members = defaultdict(list)
    for host in training:
        members[groups[host]].append(host)
    names = sorted(members)
    rng, selected_counts, first_counts = random.Random(seed), Counter(), Counter()
    for _ in range(draws):
        group_weights = Counter(rng.choice(names) for _ in names)
        weights = {host: n for group, n in group_weights.items() for host in members[group]}
        selected, _, _ = select_panel(training, phages, labels, weights=weights)
        selected_counts.update(selected)
        first_counts.update(selected[:1])
    return [{"phage": p, "draws": draws, "selected_draws": selected_counts[p],
             "selection_frequency": selected_counts[p] / draws,
             "first_draws": first_counts[p], "second_draws": selected_counts[p] - first_counts[p]}
            for p in phages]


def load_inputs(repo, plan, source_path, fetch):
    for name, sha in plan["inputs_sha256"].items():
        check((repo / name).is_file() and digest(repo / name) == sha, f"missing/changed input: {name}")
    source = plan["source_matrix"]
    if not source_path.exists() and fetch:
        source_path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(source["url"], headers={"User-Agent": "complementary-panel-study-v1"})
        with urllib.request.urlopen(request, timeout=60) as response:
            data = response.read(source["size_bytes"] + 1)
        check(len(data) == source["size_bytes"] and hashlib.sha256(data).hexdigest() == source["sha256"],
              "downloaded source matrix differs")
        source_path.write_bytes(data)
    check(source_path.is_file() and digest(source_path) == source["sha256"], "missing/changed source matrix")
    observed, unknown = read_rows(repo / OBSERVED), read_rows(repo / UNKNOWN)
    check(len(observed) == 16597 and len(unknown) == 41, "source support differs")
    labels, mlc, aliases = {}, {}, {}
    for row in observed + unknown:
        host, phage = row["bacterium"], row["phage"]
        check((host, phage) not in labels, "duplicated source pair")
        check(aliases.get(host, row["source_bacterium"]) == row["source_bacterium"], "ambiguous host alias")
        aliases[host] = row["source_bacterium"]
        value = row.get("source_mlc", "")
        mlc[host, phage] = None if value.strip().lower() in ("", "nan") else float(value)
        labels[host, phage] = None if mlc[host, phage] is None else int(mlc[host, phage] > 0)
        if "infection" in row:
            check(labels[host, phage] == int(row["infection"]), "binary label differs from MLC")
    hosts, phages = sorted(aliases), sorted({p for _, p in labels})
    check(len(hosts) == 177 and len(phages) == 94 and len(labels) == len(hosts) * len(phages), "cohort differs")
    source_rows = read_rows(source_path, ";")
    source_index = {row["bacteria"]: (line, row) for line, row in enumerate(source_rows, 2)}
    source_columns = {name: col for col, name in enumerate(source_rows[0], 1)}
    evidence = {}
    audit_evidence = {(r["bacterium"], r["phage"]): r for r in read_rows(repo / EVIDENCE)}
    source_url = source["url"].replace("raw.githubusercontent.com/", "github.com/")
    source_url = source_url.replace("/99656d8c020fcd1dbbc30a02ee10274a07655562/", "/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/")
    for host in hosts:
        line, raw = source_index[aliases[host]]
        for phage in phages:
            value = raw[phage]
            source_mlc = None if value.strip().lower() in ("", "nan", "na") else float(value)
            check(source_mlc == mlc[host, phage], f"source MLC differs: {host}/{phage}")
            previous = audit_evidence.get((host, phage), {})
            evidence[host, phage] = {"evidence_id": host + "::" + phage,
                "source_bacterium": aliases[host], "source_mlc": "" if source_mlc is None else source_mlc,
                "source_label": "" if labels[host, phage] is None else labels[host, phage],
                "source_matrix_row": line, "source_matrix_column": source_columns[phage],
                "matrix_url": source_url + f"#L{line}",
                "raw_annotation_count": previous.get("raw_annotation_count", ""),
                "raw_evidence_status": previous.get("raw_evidence_status", "aggregate matrix; raw annotations not re-audited")}
    splits = read_rows(repo / SPLITS)
    saved_scores = defaultdict(dict)
    for row in read_rows(repo / SCORES):
        if row["method"] in plan["comparators"]:
            key = (row["panel"], row["partition"], row["method"], row["bacterium"])
            check(row["phage"] not in saved_scores[key], "duplicate score")
            saved_scores[key][row["phage"]] = float(row["score"])
    return hosts, phages, labels, splits, saved_scores, evidence


def run(repo, plan, source_path, out, fetch=False):
    hosts, phages, labels, splits, scores, evidence = load_inputs(repo, plan, source_path, fetch)
    historical = {(r["panel"], r["partition"], r["method"], r["bacterium"]): r
                  for r in read_rows(repo / HISTORICAL)}
    selected_rows, stability_rows, candidate_rows, metric_rows, transitions, evidence_rows = [], [], [], [], [], []
    summaries, contrasts, split_receipts = [], [], []
    focus_evidence_keys, checked_historical = set(), 0
    for panel, expected in plan["panels"].items():
        local = [r for r in splits if r["scenario"] == panel]
        check(len(local) == len(hosts) and {r["bacterium"] for r in local} == set(hosts), "split support differs")
        membership = defaultdict(set)
        for row in local:
            membership[row["host_group"]].add(row["split"])
        check(all(len(s) == 1 for s in membership.values()), "host group crosses split")
        groups = {r["bacterium"]: r["host_group"] for r in local}
        phylogroups = {r["bacterium"]: r["phylogroup"] for r in local}
        train = sorted(r["bacterium"] for r in local if r["split"] == "train")
        check(len(train) == expected["train_hosts"], "training host count differs")
        selected, trace, stats = select_panel(train, phages, labels)
        median = statistics.median(s["rate"] for s in stats.values())
        for row in trace:
            info = stats[row["phage"]]
            selected_rows.append({"panel": panel, **row, "training_positive": info["positive"],
                "training_observed": info["observed"], "training_unknown": len(train) - info["observed"],
                "training_smoothed_prevalence": float(info["rate"]),
                "below_training_median": info["rate"] < median})
        stability = training_stability(train, phages, labels, groups,
            plan["stability"]["training_group_bootstrap_draws"], plan["stability"]["seed"])
        stability_rows.extend({"panel": panel, **row} for row in stability)
        split_receipts.append({"panel": panel, "training_hosts": len(train),
            "training_groups": len({groups[h] for h in train}), "selected": selected,
            "training_host_ids_sha256": hashlib.sha256(("\n".join(train) + "\n").encode()).hexdigest()})
        for partition in ("test", "validation"):
            targets = sorted(r["bacterium"] for r in local if r["split"] == partition)
            check(len(targets) == expected[partition + "_hosts"], "target host count differs")
            if not targets:
                continue
            check(not set(targets) & set(train), "target leaked into training")
            # A group is the leakage boundary; the B2 holdout is checked explicitly.
            check(not {groups[h] for h in targets} & {groups[h] for h in train}, "target group leaked into training")
            if panel == "phylogroup:B2":
                check(all(phylogroups[h] == "B2" for h in targets) and all(phylogroups[h] != "B2" for h in train),
                      "B2 holdout membership differs")
            first = selected[0]
            for phage in phages:
                s = stats[phage]
                marginal = sorted(h for h in targets if labels[h, first] == 0 and labels[h, phage] == 1)
                unresolved = sorted(h for h in targets if labels[h, first] is None and labels[h, phage] == 1)
                candidate_rows.append({"panel": panel, "partition": partition, "phage": phage,
                    "selected_step": selected.index(phage) + 1 if phage in selected else "",
                    "training_positive": s["positive"], "training_observed": s["observed"],
                    "training_unknown": len(train) - s["observed"], "training_smoothed_prevalence": float(s["rate"]),
                    "below_training_median": s["rate"] < median, "anchor_phage": first,
                    "extra_known_positive_training_hosts": len(s["positive_hosts"] - stats[first]["positive_hosts"]),
                    "training_anchor_unknown_extra_hosts": ";".join(sorted(h for h in s["positive_hosts"] if labels[h, first] is None)),
                    "target_extra_hosts_anchor_known_negative": len(marginal),
                    "target_extra_groups_anchor_known_negative": len({groups[h] for h in marginal}),
                    "target_extra_host_ids": ";".join(marginal),
                    "target_anchor_unknown_positive_host_ids": ";".join(unresolved)})
            local_metrics = defaultdict(dict)
            for host in targets:
                y = {p: labels[host, p] for p in phages}
                complete = all(v is not None for v in y.values())
                for method in [SELECTOR, *plan["comparators"]]:
                    if method == SELECTOR:
                        s = {p: float(2 - selected.index(p)) if p in selected else -1. for p in phages}
                    else:
                        key = (panel, partition, method, host)
                        check(key in scores and set(scores[key]) == set(phages), "incomplete score support")
                        s = scores[key]
                        if method == "phage_prevalence":
                            check(all(s[p] == float(stats[p]["rate"]) for p in phages), "saved prevalence differs from training observations")
                    low = {k: hit_probability(s, y, k) for k in (1, 2)}
                    high = {k: hit_probability(s, y, k, unknown_positive=True) for k in (1, 2)}
                    row = {"panel": panel, "partition": partition, "bacterium": host,
                        "host_group": groups[host], "phylogroup": phylogroups[host], "method": method,
                        "complete_panel": complete, "source_positive_candidates": sum(v == 1 for v in y.values()),
                        "unknown_pairs": sum(v is None for v in y.values()),
                        "hit_at_1": float(low[1]) if complete else "", "hit_at_2": float(low[2]) if complete else "",
                        "tests_consumed_capped_2": float(2 - low[1]) if complete else "",
                        "hit_at_1_lower": float(low[1]), "hit_at_1_upper": float(high[1]),
                        "hit_at_2_lower": float(low[2]), "hit_at_2_upper": float(high[2]),
                        "tests_consumed_lower": float(2 - high[1]), "tests_consumed_upper": float(2 - low[1])}
                    local_metrics[method][host] = row
                    metric_rows.append(row)
                    if method != SELECTOR and complete:
                        previous = historical[panel, partition, method, host]
                        for k in (1, 2):
                            check(abs(float(previous[f"hit@{k}"]) - float(low[k])) < 1e-12, "saved metric does not reproduce")
                            checked_historical += 1
                    # Keep all candidates in score blocks intersecting the top-two budget.
                    ordered = sorted(s, key=lambda p: (-s[p], p))
                    cutoff = s[ordered[1]]
                    for p in ordered:
                        if s[p] < cutoff:
                            break
                        rank_min = 1 + sum(v > s[p] for v in s.values())
                        rank_max = sum(v >= s[p] for v in s.values())
                        key = (panel, partition, host, method, p)
                        evidence_rows.append({"panel": panel, "partition": partition, "bacterium": host,
                            "method": method, "phage": p, "rank_min": rank_min, "rank_max": rank_max,
                            "top2_inclusion_probability": min(2 - rank_min + 1, rank_max - rank_min + 1) / (rank_max - rank_min + 1),
                            **evidence[host, p]})
                        if host in FOCUS:
                            focus_evidence_keys.add(key)
            complete_hosts = [h for h in targets if local_metrics[SELECTOR][h]["complete_panel"]]
            eligible = [h for h in complete_hosts if local_metrics[SELECTOR][h]["source_positive_candidates"] > 0]
            for method, values in local_metrics.items():
                summaries.append({"panel": panel, "partition": partition, "method": method,
                    "target_hosts": len(targets), "complete_hosts": len(complete_hosts),
                    "complete_groups": len({groups[h] for h in complete_hosts}), "eligible_hosts": len(eligible),
                    "excluded_incomplete_hosts": ";".join(h for h in targets if h not in complete_hosts),
                    "hit_at_1": statistics.mean(values[h]["hit_at_1"] for h in complete_hosts),
                    "hit_at_2": statistics.mean(values[h]["hit_at_2"] for h in complete_hosts),
                    "tests_consumed_capped_2": statistics.mean(values[h]["tests_consumed_capped_2"] for h in complete_hosts),
                    "eligible_hit_at_2": statistics.mean(values[h]["hit_at_2"] for h in eligible)})
            for method in plan["comparators"]:
                deltas, cost_deltas, rescues, losses = {}, {}, 0., 0.
                for host in complete_hosts:
                    a, b = local_metrics[SELECTOR][host], local_metrics[method][host]
                    g, reference = a["hit_at_2"], b["hit_at_2"]
                    rescued, lost = g * (1 - reference), (1 - g) * reference
                    category = ("rescue" if rescued == 1 else "new_miss" if lost == 1 else
                                "tie_dependent" if rescued > 0 or lost > 0 else "both_hit" if g else "both_miss")
                    transitions.append({"panel": panel, "partition": partition, "comparator": method,
                        "bacterium": host, "host_group": groups[host], "phylogroup": phylogroups[host],
                        "transition": category, "greedy_hit_at_2": g, "comparator_hit_at_2": reference,
                        "expected_rescued": rescued, "expected_newly_missed": lost,
                        "source_positive_candidates": a["source_positive_candidates"]})
                    deltas[host] = g - reference
                    cost_deltas[host] = a["tests_consumed_capped_2"] - b["tests_consumed_capped_2"]
                    rescues += rescued
                    losses += lost
                interval = paired_bootstrap(deltas, groups, 2000, plan["evaluation"]["seed"])
                cost_interval = paired_bootstrap(cost_deltas, groups, 2000, plan["evaluation"]["seed"])
                contrasts.append({"panel": panel, "partition": partition, "comparator": method,
                    "complete_hosts": len(complete_hosts), "expected_rescued_hosts": rescues,
                    "expected_newly_missed_hosts": losses, "net_hosts": rescues - losses,
                    "hit_at_2_delta": statistics.mean(deltas.values()), "hit_delta_ci_low": interval[0],
                    "hit_delta_ci_high": interval[1], "tests_consumed_delta": statistics.mean(cost_deltas.values()),
                    "cost_delta_ci_low": cost_interval[0], "cost_delta_ci_high": cost_interval[1]})
    write_csv(out / "selected_panels.csv", selected_rows)
    write_csv(out / "training_bootstrap_selection.csv", stability_rows)
    write_csv(out / "candidate_marginal_coverage.csv", candidate_rows)
    write_csv(out / "per_host_metrics.csv", metric_rows)
    write_csv(out / "summary.csv", summaries)
    write_csv(out / "paired_comparisons.csv", contrasts)
    write_csv(out / "host_transitions.csv", transitions)
    write_csv(out / "top2_source_evidence.csv", evidence_rows)
    write_csv(out / "focus_source_evidence.csv", [r for r in evidence_rows if r["bacterium"] in FOCUS])
    # Evidence for each affected host, including successful controls and both policies.
    affected = {(r["panel"], r["partition"], r["bacterium"])
                for r in transitions if r["transition"] in ("rescue", "new_miss", "tie_dependent")}
    write_csv(out / "changed_host_evidence.csv", [r for r in evidence_rows
              if (r["panel"], r["partition"], r["bacterium"]) in affected])
    write_json(out / "QA.json", {"status": "verified", "protocol_sha256": digest(repo / PROTOCOL),
        "source_matrix_sha256": digest(source_path), "source_observed_pairs": 16597, "source_unknown_pairs": 41,
        "historical_hit_values_reproduced": checked_historical, "maximum_historical_error_tolerance": 1e-12,
        "panel_receipts": split_receipts, "metric_rows": len(metric_rows),
        "independent_validation": False, "python": platform.python_version()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-matrix", type=Path, default=REPO / "data/raw/complementary_panel/mlc.csv")
    parser.add_argument("--fetch", action="store_true", help="fetch only the pinned 158 KB source matrix")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    check(not args.out.exists(), "output must be a new directory")
    plan = json.loads((REPO / PROTOCOL).read_text())
    check(plan["budget"] == 2, "this implementation is locked to budget two")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".complementary-panel-", dir=args.out.parent))
    try:
        run(REPO, plan, args.source_matrix, stage, args.fetch)
        stage.rename(args.out)
    except BaseException:
        shutil.rmtree(stage)
        raise
    print(f"Verified comparison: {args.out}")


if __name__ == "__main__":
    main()
