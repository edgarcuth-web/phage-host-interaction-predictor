"""Offline verification of the registered sensitivity analysis, without fitting."""
from fractions import Fraction
from itertools import product
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.genophi_seed_protocols import REGISTRATION, INPUT_FREEZE, seed_paths
from phage_cocktail.retrieval_signflip import classify_primary_contrasts, classify_prevalence_swings


def independent_group_p(left, right, groups, *, metric, eligible_only):
    """Small real primary contrasts: brute-force Fraction arithmetic, no shared engine."""
    a, b = left.set_index("bacterium").sort_index(), right.set_index("bacterium").sort_index()
    if list(a.index) != list(b.index):
        raise ValueError("independent QA requires identical host support")
    totals = {}
    for host in a.index:
        if eligible_only and not a.loc[host, "eligible"]:
            continue
        value = Fraction(float(a.loc[host, metric])) - Fraction(float(b.loc[host, metric]))
        group = groups[host]
        totals[group] = totals.get(group, Fraction(0)) + value
    weights = [v for v in totals.values() if v]
    if len(weights) > 20:
        raise ValueError("independent primary enumeration unexpectedly exceeds its QA budget")
    observed, tails = abs(sum(weights)), 0
    for signs in product((-1, 1), repeat=len(weights)):
        tails += abs(sum(sign * value for sign, value in zip(signs, weights))) >= observed
    return tails / (1 << len(weights))


def verify(repo, root):
    repo, root = Path(repo).resolve(), Path(root).resolve()
    report = read_json_object(root / "report.json")
    registration, frozen = read_json_object(repo / REGISTRATION), read_json_object(repo / INPUT_FREEZE)
    if sha256_file(repo / REGISTRATION) != report["registration_sha256"]:
        raise ValueError("registered decision rule changed")
    if sha256_file(repo / INPUT_FREEZE) != report["input_freeze_sha256"]:
        raise ValueError("input-freeze record changed")
    for relative, digest in frozen["files_sha256"].items():
        if sha256_file(repo / relative) != digest:
            raise ValueError("frozen input changed")
    for name, digest in report["output_sha256"].items():
        if sha256_file(root / name) != digest:
            raise ValueError("changed sensitivity output")
    if sha256_file(repo / "scripts/summarize_genophi_seed_sensitivity.py") != report["reporting_script_sha256"]:
        raise ValueError("reporting script changed")
    if sha256_file(repo / "src/phage_cocktail/retrieval_signflip.py") != report["statistical_code_sha256"]:
        raise ValueError("registered statistical implementation changed")
    per_host = pd.read_csv(root / "per_host_metrics.csv", dtype={"bacterium": str}, float_precision="round_trip")
    contrasts = pd.read_csv(root / "paired_contrasts.csv", float_precision="round_trip")
    expected = 3 * 2 * 6 * 5 * 2
    if len(contrasts) != expected or contrasts.non_binary_hosts.isna().any():
        raise ValueError("missing pairwise contrasts or non-binary support counts")
    primary, swings, checked = {}, {}, 0
    for seed in [42, 43, 44]:
        manifest_path = (repo / "data/studies/results/genophi_host_split.v1.csv" if seed == 42
                         else repo / seed_paths(seed)["manifest"])
        manifest = pd.read_csv(manifest_path, dtype=str)
        groups = manifest.set_index("bacterium").host_group.to_dict()
        if seed != 42:
            for kind in ("baseline_output", "comparator_output"):
                source = repo / seed_paths(seed)[kind]
                qa = read_json_object(source / "NUMERIC_QA.json")
                if qa["status"] != "verified":
                    raise ValueError("unverified seed source")
                for relative, digest in qa["published_output_sha256"].items():
                    if sha256_file(source / relative) != digest:
                        raise ValueError("changed verified seed-source payload")
            imported = read_json_object(repo / seed_paths(seed)["comparator_output"] / "checkpoints/reuse_features.json")
            if imported["result"]["new_clustering_search_assignment_calls"] != 0:
                raise ValueError("feature workflow unexpectedly repeated")
        for partition in ("validation", "test"):
            for method in registration["reported_methods"]:
                rows = per_host[per_host.seed.eq(seed) & per_host.partition.eq(partition) & per_host.model.eq(method)]
                metric = next(m for m in report["metrics"] if (m["seed"], m["partition"], m["model"]) == (seed, partition, method))
                if len(rows) != metric["complete_panel_hosts"] or int(rows.eligible.sum()) != metric["eligible_hosts"]:
                    raise ValueError("aggregate support changed")
                for name in registration["statistics"]["metrics"]:
                    if not np.isclose(rows[name].mean(), metric["all_hosts"][name], rtol=0, atol=1e-14):
                        raise ValueError("aggregate retrieval mean changed")
        for eligible in (False, True):
            selected = per_host[per_host.seed.eq(seed) & per_host.partition.eq("test")]
            a, b = selected[selected.model.eq("genophi_median")], selected[selected.model.eq("phage_prevalence")]
            result = contrasts[contrasts.seed.eq(seed) & contrasts.partition.eq("test")
                               & contrasts.left.eq("genophi_median") & contrasts.right.eq("phage_prevalence")
                               & contrasts.metric.eq("hit@2") & contrasts.eligible_only.eq(eligible)].iloc[0]
            p = independent_group_p(a, b, groups, metric="hit@2", eligible_only=eligible)
            if p != result.p_value:
                raise ValueError("independent real-data exact enumeration disagrees")
            checked += 1
            if seed != 42 and not eligible:
                primary[seed] = {"status": "verified", "delta": float(result.delta), "p_value": p,
                                 "interval": json.loads(result.bootstrap_interval)}
        if seed != 42:
            values = {m["partition"]: m["all_hosts"]["hit@2"] for m in report["metrics"]
                      if m["seed"] == seed and m["model"] == "phage_prevalence"}
            swings[seed] = values["test"] - values["validation"]
    if classify_primary_contrasts(primary, registration) != report["primary_decision"]:
        raise ValueError("published primary decision differs from preregistration")
    if classify_prevalence_swings(swings, registration) != report["prevalence_swing_decision"]:
        raise ValueError("published prevalence-swing classification differs from preregistration")
    return {"status": "verified", "report_sha256": sha256_file(root / "report.json"),
            "verification_script_sha256": sha256_file(__file__), "pairwise_contrasts": len(contrasts),
            "primary_Hit2_p_values_independently_enumerated": checked,
            "fractional_hosts_preserved": True, "registered_decisions_reproduced": True,
            "verified_sensitivity_seeds": [43, 44], "additional_model_fits": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("verification receipt exists; choose a new immutable path")
    result = verify(Path(__file__).resolve().parents[1], args.run)
    atomic_write_json(result, args.out)
    print(result)


if __name__ == "__main__":
    main()
