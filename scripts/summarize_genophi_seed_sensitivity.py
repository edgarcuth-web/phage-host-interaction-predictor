"""Apply preregistered statistics to immutable matched-seed results; never fit models."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file, staged_directory
from phage_cocktail.genophi_seed_protocols import REGISTRATION, INPUT_FREEZE, seed_paths, assert_baseline_model_split
from phage_cocktail.retrieval import paired_retrieval_bootstrap
from phage_cocktail.retrieval_signflip import (
    paired_group_signflip, classify_primary_contrasts, classify_prevalence_swings,
)


def read_frame(path):
    return pd.read_csv(path, dtype={"bacterium": str}, float_precision="round_trip")


def verify_bundle(root):
    root = Path(root)
    qa = read_json_object(root / "NUMERIC_QA.json")
    if qa["status"] != "verified":
        raise ValueError("an unverified seed bundle cannot enter the sensitivity report")
    for relative, digest in qa["published_output_sha256"].items():
        if sha256_file(root / relative) != digest:
            raise ValueError(f"changed verified published result: {root}/{relative}")
    return qa


def panel_summary(frame, split, method, seed, partition):
    hosts = set(split.loc[split.split.eq(partition), "bacterium"])
    selected = frame[frame.partition.eq(partition) & frame.model.eq(method)].copy()
    if selected.bacterium.duplicated().any() or not set(selected.bacterium).issubset(hosts):
        raise ValueError("invalid per-host panel support")
    eligible = selected[selected.eligible.astype(bool)]
    metrics = ["hit@1", "hit@2", "hit@3", "hit@5", "mrr"]
    return {
        "seed": seed, "scenario": "grouped_host", "partition": partition, "model": method,
        "train_hosts": int(split.split.eq("train").sum()), "held_out_hosts": len(hosts),
        "complete_panel_hosts": len(selected), "eligible_hosts": len(eligible),
        "all_negative_hosts": len(selected) - len(eligible),
        "incomplete_hosts_excluded": sorted(hosts - set(selected.bacterium)), "candidate_phages": 94,
        "positive_candidates_median": float(selected.positive_candidates.median()) if len(selected) else None,
        "all_hosts": {m: float(selected[m].mean()) if len(selected) else None for m in metrics},
        "eligible_hosts_only": {m: float(eligible[m].mean()) if len(eligible) else None
                                for m in [*metrics, "recall@1", "recall@2", "recall@3", "recall@5"]},
        "fractional_hosts": {m: int((~selected[m].isin([0., 1.])).sum()) for m in metrics},
    }, selected


def contrast_record(left, right, groups, *, seed, partition, left_name, right_name, metric, eligible, registration):
    bootstrap = paired_retrieval_bootstrap(left, right, groups, metric=metric, eligible_only=eligible,
                                           **registration["statistics"]["bootstrap"])
    permutation = paired_group_signflip(left, right, groups, metric=metric, eligible_only=eligible,
        max_exact_groups=registration["statistics"]["sign_flip"]["max_exact_nonzero_groups"])
    return {"seed": seed, "scenario": "grouped_host", "partition": partition,
            "left": left_name, "right": right_name, "metric": metric, "eligible_only": eligible,
            **permutation, "bootstrap_interval": bootstrap["interval"],
            "bootstrap_resamples": bootstrap["requested_resamples"], "bootstrap_seed": bootstrap["seed"],
            "bootstrap_valid_resamples": bootstrap["valid_resamples"],
            "p_value_scope": "Descriptive unadjusted, except the registered two new-seed primary contrasts"}


def summarize(repo, execution_root, destination):
    repo, execution_root, destination = Path(repo).resolve(), Path(execution_root).resolve(), Path(destination)
    if destination.exists():
        raise ValueError("sensitivity summary directory exists; keep completed analyses immutable")
    registration = read_json_object(repo / REGISTRATION)
    frozen = read_json_object(repo / INPUT_FREEZE)
    if sha256_file(repo / REGISTRATION) != frozen["registration_sha256"]:
        raise ValueError("registered decision rule changed")
    for relative, digest in frozen["files_sha256"].items():
        if sha256_file(repo / relative) != digest:
            raise ValueError(f"changed frozen sensitivity input: {relative}")
    execution = read_json_object(execution_root / "execution.json")
    summaries, contrasts, all_hosts, provenance, accounts = [], [], [], {}, {}
    reference_comp = repo / "data/studies/results/genophi_comparator.v1"
    reference_base = repo / "data/studies/results/genome_host_retrieval.v1"
    for seed in [42, *registration["sensitivity_seeds"]]:
        if seed == 42:
            comp, baseline = reference_comp, reference_base
            manifest_path = repo / "data/studies/results/genophi_host_split.v1.csv"
            accounts[str(seed)] = {"status": "existing_primary_reference_no_refit"}
        else:
            accounts[str(seed)] = execution["seeds"][str(seed)]
            if accounts[str(seed)]["status"] != "verified":
                continue
            paths = seed_paths(seed)
            comp, baseline = repo / paths["comparator_output"], repo / paths["baseline_output"]
            manifest_path = repo / paths["manifest"]
            for root in (comp, baseline):
                verify_bundle(root)
            manifest = pd.read_csv(manifest_path, dtype=str)
            for name in ("model_00.json", "model_00_host_4mer_logistic.json", "model_00_hybrid_4mer_logistic.json"):
                assert_baseline_model_split(read_json_object(baseline / name), manifest)
        comp_report, baseline_report = read_json_object(comp / "report.json"), read_json_object(baseline / "report.json")
        for root, report in [(comp, comp_report), (baseline, baseline_report)]:
            if sha256_file(root / "per_host_metrics.csv") != report["output_sha256"]["per_host_metrics.csv"]:
                raise ValueError("changed source per-host metrics")
        provenance[str(seed)] = {
            "comparator_report_sha256": sha256_file(comp / "report.json"),
            "baseline_report_sha256": sha256_file(baseline / "report.json"),
            "comparator_per_host_sha256": sha256_file(comp / "per_host_metrics.csv"),
            "baseline_per_host_sha256": sha256_file(baseline / "per_host_metrics.csv"),
            "host_manifest_sha256": sha256_file(manifest_path),
            "comparator_versions": comp_report["provenance"]["versions"],
            "baseline_versions": baseline_report["provenance"]["versions"],
        }
        split = pd.read_csv(manifest_path, dtype=str)
        groups = split.set_index("bacterium").host_group.to_dict()
        frames = pd.concat([
            read_frame(comp / "per_host_metrics.csv").query("scenario == 'grouped_host' and model in ['genophi_median', 'phage_prevalence']"),
            read_frame(baseline / "per_host_metrics.csv").query("scenario == 'grouped_host' and model in ['host_4mer_logistic', 'hybrid_4mer_logistic']"),
        ], ignore_index=True)
        all_hosts.append(frames.assign(seed=seed))
        for partition in registration["statistics"]["partitions"]:
            tables = {}
            for method in registration["reported_methods"]:
                summary, table = panel_summary(frames, split, method, seed, partition)
                summaries.append(summary)
                tables[method] = table
                source = comp_report if method in ["genophi_median", "phage_prevalence"] else baseline_report
                previous = next(m for m in source["metrics"] if m["scenario"] == "grouped_host"
                                and m["partition"] == partition and m["model"] == method)
                for denominator in ["all_hosts", "eligible_hosts_only"]:
                    for metric, value in summary[denominator].items():
                        if not np.isclose(value, previous[denominator][metric], rtol=0, atol=1e-14):
                            raise ValueError("source report and independent per-host aggregation disagree")
            for left, right in combinations(registration["reported_methods"], 2):
                for metric in registration["statistics"]["metrics"]:
                    for eligible in (False, True):
                        contrasts.append(contrast_record(tables[left], tables[right], groups,
                            seed=seed, partition=partition, left_name=left, right_name=right,
                            metric=metric, eligible=eligible, registration=registration))
    primary = {}
    swings = {}
    for seed in registration["sensitivity_seeds"]:
        if accounts[str(seed)]["status"] != "verified":
            primary[seed] = {"status": "blocked"}
            continue
        value = next(c for c in contrasts if c["seed"] == seed and c["partition"] == "test"
                     and c["left"] == "genophi_median" and c["right"] == "phage_prevalence"
                     and c["metric"] == "hit@2" and not c["eligible_only"])
        primary[seed] = {"status": "verified", "delta": value["delta"], "p_value": value["p_value"],
                         "interval": value["bootstrap_interval"]}
        prev = {m["partition"]: m["all_hosts"]["hit@2"] for m in summaries
                if m["seed"] == seed and m["model"] == "phage_prevalence"}
        swings[seed] = prev["test"] - prev["validation"]
    decision = classify_primary_contrasts(primary, registration)
    swing_decision = classify_prevalence_swings(swings, registration)
    b2 = []
    comp = read_frame(reference_comp / "per_host_metrics.csv")
    base = read_frame(reference_base / "per_host_metrics.csv")
    original_groups = pd.read_csv(repo / "data/studies/results/genophi_host_split.v1.csv", dtype=str).set_index("bacterium").host_group.to_dict()
    left = comp[comp.scenario.eq("phylogroup:B2") & comp.model.eq("genophi_median")]
    for right_name, source in [("phage_prevalence", comp), ("host_4mer_logistic", base)]:
        right = source[source.scenario.eq("phylogroup:B2") & source.model.eq(right_name)]
        value = contrast_record(left, right, original_groups, seed=42, partition="test",
            left_name="genophi_median", right_name=right_name, metric="hit@2", eligible=False, registration=registration)
        value["scenario"] = "phylogroup:B2"
        b2.append(value)
    report = {
        "study_id": registration["study_id"], "status": "verified" if execution["status"] == "verified" else "blocked",
        "registration_sha256": sha256_file(repo / REGISTRATION), "input_freeze_sha256": sha256_file(repo / INPUT_FREEZE),
        "reporting_script_sha256": sha256_file(__file__),
        "statistical_code_sha256": sha256_file(repo / "src/phage_cocktail/retrieval_signflip.py"),
        "source_provenance": provenance, "seed_accounting": accounts, "metrics": summaries,
        "primary_decision": decision, "prevalence_swing_decision": swing_decision,
        "seed_splits": frozen["seed_splits"], "B2_seed42_descriptive_no_refit": b2,
        "limits": registration["limits"], "execution_wall_seconds": execution["wall_seconds"],
        "execution_clock_scope": execution.get("wall_time_scope", "monotonic execution clock"),
    }
    with staged_directory(destination) as stage:
        per_host = pd.concat(all_hosts, ignore_index=True)
        per_host.to_csv(stage / "per_host_metrics.csv", index=False)
        pd.DataFrame(contrasts).to_csv(stage / "paired_contrasts.csv", index=False)
        flat = []
        for summary in summaries:
            common = {k: v for k, v in summary.items() if k not in ["all_hosts", "eligible_hosts_only", "fractional_hosts"]}
            flat.append({**common, **{f"all_{k}": v for k, v in summary["all_hosts"].items()},
                         **{f"eligible_{k}": v for k, v in summary["eligible_hosts_only"].items()}})
        pd.DataFrame(flat).to_csv(stage / "metrics.csv", index=False)
        (stage / "execution.json").write_bytes((execution_root / "execution.json").read_bytes())
        interruption = execution_root / "interruption.before_comparator_start.json"
        if interruption.is_file():
            (stage / interruption.name).write_bytes(interruption.read_bytes())
        report["output_sha256"] = {p.name: sha256_file(p) for p in sorted(stage.iterdir())}
        atomic_write_json(report, stage / "report.json")
    print(json.dumps({"primary_decision": decision, "prevalence_swing_decision": swing_decision}, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    summarize(Path(__file__).resolve().parents[1], args.execution, args.outdir)


if __name__ == "__main__":
    main()
