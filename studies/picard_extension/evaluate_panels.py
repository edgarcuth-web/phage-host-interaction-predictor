#!/usr/bin/env python3
"""Evaluate the two already fixed pairs after genomic admission is sealed."""
import argparse
import gzip
import importlib.util
import json
import shutil
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path

import audit_cohort as audit

ROOT, STUDY = audit.ROOT, audit.STUDY
HELPER = "studies/complementary_panel/run_comparison.py"


def helper_module():
    spec = importlib.util.spec_from_file_location("frozen_complementary_metrics", ROOT / HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(sources, admission, out):
    lock = json.loads((STUDY / "evaluation_implementation.json").read_text())
    for path, digest in lock["files_sha256"].items():
        audit.check(audit.sha(ROOT / path) == digest, f"evaluation implementation changed: {path}")
    freeze = json.loads((STUDY / "cohort_freeze.json").read_text())
    for path, digest in freeze["admission_output_sha256"].items():
        audit.check(audit.sha(admission / path) == digest, f"admission changed after freezing: {path}")
    plan = audit.load_sources(sources)
    roster = audit.read(admission / "admitted_hosts.csv")
    audit.check(len(roster) == freeze["eligible_hosts"] and len({r["component"] for r in roster}) == freeze["eligible_groups"], "admitted support differs")
    metrics = helper_module()
    matrix = audit.read(sources / "source_mlc.csv", ";")
    source = {r["bacteria"]: (line, r) for line, r in enumerate(matrix, 2)}
    phages = sorted(p for p in matrix[0] if p not in ("bacteria", "T4LD", "T7_Portugal"))
    columns = {p: index for index, p in enumerate(matrix[0], 1)}
    labels = {(h, p): None if r[p].strip().lower() in ("", "nan", "na") else int(float(r[p]) > 0)
              for h, (_, r) in source.items() for p in phages}
    splits = audit.read(ROOT / audit.SPLITS)
    old_alias = {r["bacterium"]: r["source_bacterium"] for r in splits}
    # Reconcile every reused training source reading with the original observation ledger.
    with gzip.open(ROOT / "data/studies/results/paired_interaction.v1/observed_pairs.csv.gz", "rt") as f:
        import csv
        observed = list(csv.DictReader(f))
    for r in observed:
        host = old_alias[r["bacterium"]]
        audit.check(float(source[host][1][r["phage"]]) == float(r["source_mlc"]) and labels[host, r["phage"]] == int(r["infection"]), "original MLC labels differ")
    policies, prior_rows = {}, []
    trace = audit.read(ROOT / "data/studies/results/complementary_panel.v1/selected_panels.csv")
    for name, scenario, pair_key in (("primary", "grouped_host", "primary_pair"), ("secondary", "phylogroup:B2", "secondary_pair")):
        pair = plan["evaluation"][pair_key]
        previous = [r["phage"] for r in trace if r["panel"] == scenario]
        audit.check(previous == pair, "fixed pair differs from the original study")
        train = sorted(old_alias[r["bacterium"]] for r in splits if r["scenario"] == scenario and r["split"] == "train")
        stats = metrics.candidate_stats(train, phages, labels)
        policies[name + "_complement"] = {p: float(2 - pair.index(p)) if p in pair else -1. for p in phages}
        policies[name + "_prevalence"] = {p: float(stats[p]["rate"]) for p in phages}
        prior_rows.extend({"scope": name, "phage": p, "training_hosts": len(train),
            "training_positive": stats[p]["positive"], "training_observed": stats[p]["observed"],
            "score": float(stats[p]["rate"])} for p in phages)
    host_rows, evidence = [], []
    by_key = {}
    matrix_url = plan["resources"]["source_mlc"]["url"].replace("raw.githubusercontent.com/", "github.com/")
    matrix_url = matrix_url.replace("/" + plan["source_revision"] + "/", "/blob/" + plan["source_revision"] + "/")
    for host_record in roster:
        h = host_record["source_host"]
        y = {p: labels[h, p] for p in phages}
        audit.check(all(value is not None for value in y.values()), "incomplete host entered evaluation")
        for method, scores in policies.items():
            hit1, hit2 = [float(metrics.hit_probability(scores, y, k)) for k in (1, 2)]
            row = {"source_host": h, "component": host_record["component"],
                "source_phylogroup": host_record["source_phylogroup"], "method": method,
                "positive_candidates": sum(y.values()), "hit_at_1": hit1, "hit_at_2": hit2,
                "tests_consumed_capped_2": 2 - hit1}
            host_rows.append(row)
            by_key[h, method] = row
            cutoff = sorted(scores.values(), reverse=True)[1]
            for p in sorted(scores, key=lambda p: (-scores[p], p)):
                if scores[p] < cutoff:
                    break
                low, high = 1 + sum(v > scores[p] for v in scores.values()), sum(v >= scores[p] for v in scores.values())
                evidence.append({"source_host": h, "component": host_record["component"],
                    "source_phylogroup": host_record["source_phylogroup"], "method": method, "phage": p,
                    "rank_min": low, "rank_max": high, "top2_inclusion_probability": min(3 - low, high - low + 1) / (high - low + 1),
                    "source_mlc": float(source[h][1][p]), "source_positive": labels[h, p],
                    "source_matrix_row": source[h][0], "source_matrix_column": columns[p],
                    "matrix_url": matrix_url + f"#L{source[h][0]}", "raw_annotations_independently_checked": False})
    summaries, contrasts, transitions, marginal = [], [], [], []
    groups = {r["source_host"]: r["component"] for r in roster}
    for stratum in ("all", "B2", "non_B2"):
        hosts = [r["source_host"] for r in roster if stratum == "all" or (r["source_phylogroup"] == "B2") == (stratum == "B2")]
        if not hosts:
            continue
        for method in policies:
            local = [by_key[h, method] for h in hosts]
            summaries.append({"stratum": stratum, "method": method, "hosts": len(hosts),
                "groups": len({groups[h] for h in hosts}), "all_source_negative_hosts": sum(r["positive_candidates"] == 0 for r in local),
                **{field: statistics.mean(r[field] for r in local) for field in ("hit_at_1", "hit_at_2", "tests_consumed_capped_2")}})
        for name, pair_key in (("primary", "primary_pair"), ("secondary", "secondary_pair")):
            difference, cost_difference, rescued, lost = {}, {}, [], []
            for h in hosts:
                a, b = by_key[h, name + "_complement"], by_key[h, name + "_prevalence"]
                g, reference = a["hit_at_2"], b["hit_at_2"]
                difference[h], cost_difference[h] = g - reference, a["tests_consumed_capped_2"] - b["tests_consumed_capped_2"]
                if g > reference:
                    rescued.append(h)
                if g < reference:
                    lost.append(h)
                if stratum == "all":
                    transitions.append({"source_host": h, "component": groups[h], "scope": name,
                        "complement_hit_at_2": g, "prevalence_hit_at_2": reference,
                        "expected_rescued": g * (1 - reference), "expected_newly_missed": (1 - g) * reference,
                        "transition": "rescue" if g > reference else "new_miss" if g < reference else "both_hit" if g else "both_miss"})
            enough = len({groups[h] for h in hosts}) >= 2
            hit_ci = metrics.paired_bootstrap(difference, groups, 2000, 20261115) if enough else (None, None)
            cost_ci = metrics.paired_bootstrap(cost_difference, groups, 2000, 20261115) if enough else (None, None)
            contrasts.append({"stratum": stratum, "scope": name, "hosts": len(hosts), "groups": len({groups[h] for h in hosts}),
                "rescued_hosts": len(rescued), "newly_missed_hosts": len(lost), "net_hosts": len(rescued) - len(lost),
                "hit_at_2_delta": statistics.mean(difference.values()), "hit_delta_ci_low": hit_ci[0], "hit_delta_ci_high": hit_ci[1],
                "tests_consumed_delta": statistics.mean(cost_difference.values()), "cost_delta_ci_low": cost_ci[0], "cost_delta_ci_high": cost_ci[1],
                "rescued_host_ids": ";".join(rescued), "newly_missed_host_ids": ";".join(lost)})
            anchor, second = plan["evaluation"][pair_key]
            extra = [h for h in hosts if labels[h, anchor] == 0 and labels[h, second] == 1]
            marginal.append({"stratum": stratum, "scope": name, "anchor_phage": anchor, "second_phage": second,
                "extra_hosts_beyond_anchor": len(extra), "extra_groups_beyond_anchor": len({groups[h] for h in extra}),
                "extra_host_ids": ";".join(extra)})
    for filename, records in (("per_host_metrics.csv", host_rows), ("summary.csv", summaries),
            ("paired_comparisons.csv", contrasts), ("host_transitions.csv", transitions),
            ("candidate_marginal_coverage.csv", marginal), ("top2_source_evidence.csv", evidence),
            ("training_prevalence.csv", prior_rows)):
        audit.write_csv(out / filename, records)
    changed = {r["source_host"] for r in transitions if r["transition"] in ("rescue", "new_miss")}
    audit.write_csv(out / "changed_host_evidence.csv", [r for r in evidence if r["source_host"] in changed])
    audit.write_json(out / "QA.json", {"status": "verified_by_runner", "admitted_hosts": len(roster),
        "admitted_groups": len(set(groups.values())), "source_known_pairs_reconciled": len(observed),
        "new_models_fitted": 0, "new_candidate_selections": 0, "genome_models_scored": 0,
        "cohort_freeze_sha256": audit.sha(STUDY / "cohort_freeze.json"), "protocol_sha256": audit.sha(STUDY / "protocol.json"),
        "code_sha256": audit.sha(Path(__file__)), "output_sha256": {p.name: audit.sha(p) for p in sorted(out.iterdir())}})
    print(json.dumps(contrasts, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, default=ROOT / "data/raw/picard_extension_sources")
    p.add_argument("--admission", type=Path, default=ROOT / "data/studies/results/picard_extension.v1/admission")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    audit.check(not args.out.exists(), "output must be a new directory")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".panel-evaluation-", dir=args.out.parent))
    try:
        evaluate(args.sources, args.admission, stage)
        stage.rename(args.out)
    except BaseException:
        shutil.rmtree(stage)
        raise


if __name__ == "__main__":
    main()
