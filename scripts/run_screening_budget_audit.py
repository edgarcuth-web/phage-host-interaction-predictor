"""Audit immutable rankings at fixed screening budgets; never fit a model."""
import argparse
import importlib.metadata
from pathlib import Path

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import (
    atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory,
)
from phage_cocktail.cohort import curate_mlc
from phage_cocktail.data import load_interactions
from phage_cocktail.host_study import verify_resources
from phage_cocktail.registry import verify_registered_dataset
from phage_cocktail.retrieval import evaluate_retrieval

PROTOCOL = "data/studies/screening_budget.v1.json"
IMPLEMENTATION = "data/studies/screening_budget.implementation.v1.json"
SCRIPT = "scripts/run_screening_budget_audit.py"
HISTORICAL = ("hit@1", "hit@2", "hit@3", "hit@5", "mrr")


def read_frame(path):
    return pd.read_csv(path, dtype={"bacterium": str, "phage": str, "host_group": str},
                       float_precision="round_trip")


def verify_hashes(root, hashes):
    root = Path(root).resolve()
    for relative, digest in hashes.items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file() or sha256_file(path) != digest:
            raise ValueError(f"missing or changed frozen input: {relative}")


def verify_implementation(repo):
    freeze = read_json_object(Path(repo) / IMPLEMENTATION)
    if sha256_file(Path(repo) / PROTOCOL) != freeze["protocol_sha256"]:
        raise ValueError("changed screening-budget protocol")
    verify_hashes(repo, freeze["files_sha256"])
    return freeze


def budget_metrics(frame, budgets):
    """Derive capped costs only when every preceding hit probability exists."""
    result = frame.copy()
    for k in budgets:
        preceding = [f"hit@{j}" for j in range(1, k)]
        if all(name in result and np.isfinite(result[name]).all() for name in preceding):
            result[f"tests@{k}"] = k - result[preceding].sum(axis=1) if preceding else 1.
        else:
            result[f"tests@{k}"] = np.nan
        if f"hit@{k}" not in result:
            result[f"hit@{k}"] = np.nan
    return result


def paired_intervals(left, right, groups, metrics, *, eligible_only=False,
                     resamples=2000, seed=20261024, confidence=.95):
    """Batch the existing paired group-percentile procedure over metrics.

    Draw groups uniformly, then retain all their hosts. The denominator is the
    number of sampled hosts (eligible hosts for the secondary denominator).
    Shared draws across columns yield pointwise, not simultaneous, intervals.
    """
    if resamples < 2 or not 0 < confidence < 1:
        raise ValueError("invalid bootstrap configuration")
    if (left.bacterium.duplicated().any() or right.bacterium.duplicated().any()
            or set(left.bacterium) != set(right.bacterium)):
        raise ValueError("paired metrics require unique, identical host support")
    a, b = left.set_index("bacterium").sort_index(), right.set_index("bacterium").sort_index()
    if not a[["eligible", "positive_candidates"]].equals(b[["eligible", "positive_candidates"]]):
        raise ValueError("paired label support differs")
    if set(a.index) - set(groups) or any(not isinstance(groups[h], str) or not groups[h] for h in a.index):
        raise ValueError("every host requires a frozen group")
    difference = a[metrics].to_numpy(float) - b[metrics].to_numpy(float)
    if not np.isfinite(difference).all():
        raise ValueError("bootstrap metrics must be finite")
    keep = a.eligible.to_numpy(bool) if eligible_only else np.ones(len(a), bool)
    names = sorted({groups[h] for h in a.index})
    group_ids = np.array([groups[h] for h in a.index])
    sums = np.array([difference[(group_ids == g) & keep].sum(axis=0) for g in names])
    sizes = np.array([int(((group_ids == g) & keep).sum()) for g in names])
    means = difference[keep].mean(axis=0) if keep.any() else np.full(len(metrics), np.nan)
    if len(names) < 2 or not keep.any():
        return {m: {"delta": float(v) if np.isfinite(v) else None, "interval": None,
                    "valid_resamples": 0, "status": "insufficient_support"}
                for m, v in zip(metrics, means)}
    draws = np.random.default_rng(seed).integers(len(names), size=(resamples, len(names)))
    multiplicity = np.zeros((resamples, len(names)), dtype=int)
    np.add.at(multiplicity, (np.arange(resamples)[:, None], draws), 1)
    denominators = multiplicity @ sizes
    valid = denominators > 0
    estimates = (multiplicity[valid] @ sums) / denominators[valid, None]
    alpha = (1 - confidence) / 2
    intervals = np.quantile(estimates, [alpha, 1-alpha], axis=0)
    return {m: {"delta": float(v), "interval": intervals[:, i].tolist(),
                "valid_resamples": int(valid.sum()), "status": "conditional_development_interval"}
            for i, (m, v) in enumerate(zip(metrics, means))}


def check_scores(frame, hosts, phages):
    expected = pd.MultiIndex.from_product([sorted(hosts), sorted(phages)])
    actual = pd.MultiIndex.from_frame(frame[["bacterium", "phage"]])
    if actual.has_duplicates or set(actual) != set(expected) or not np.isfinite(frame.score).all():
        raise ValueError("scores must cover the exact finite host x phage panel")


def check_reference(calculated, reference):
    if reference.bacterium.duplicated().any() or set(calculated.bacterium) != set(reference.bacterium):
        raise ValueError("historical complete-host support differs")
    a, b = calculated.set_index("bacterium").sort_index(), reference.set_index("bacterium").sort_index()
    if not np.array_equal(a.positive_candidates, b.positive_candidates) or not np.array_equal(a.eligible, b.eligible):
        raise ValueError("historical label support differs")
    for name in HISTORICAL:
        if name in a and (not np.isfinite(b[name]).all() or
                          not np.allclose(a[name], b[name], rtol=0, atol=1e-12)):
            raise ValueError(f"historical retrieval metric differs: {name}")


def select_reference(frame, panel, partition, method):
    return frame[frame.scenario.eq(panel["scenario"]) & frame.partition.eq(partition)
                 & frame.model.eq(method)][["bacterium", "positive_candidates", "eligible", *HISTORICAL]].copy()


def load_study(repo, data_root, *, fetch=False):
    """Verify inputs and expose source panels without computing new metrics."""
    repo = Path(repo)
    plan = read_json_object(repo / PROTOCOL)
    verify_hashes(repo, plan["files_sha256"])
    for bundle in plan["verified_bundles"]:
        qa = read_json_object(repo / bundle / "NUMERIC_QA.json")
        if qa["status"] != "verified":
            raise ValueError("source bundle has not passed its numerical QA")
        verify_hashes(repo / bundle, qa["published_output_sha256"])
    resources = verify_resources({"resources": plan["resources"]}, data_root, fetch=fetch)
    subset = resources[("table", "subset")]
    verify_registered_dataset(repo / "data/registry/genophi_ecoli_subset.v1.json", subset)
    pairs = load_interactions(subset)
    observed, excluded, _ = curate_mlc(pairs, pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0))
    phages = sorted(pairs.phage.unique())
    if (len(phages) != plan["candidate_phages"] or len(observed) != plan["expected_observed_rows"]
            or len(excluded) != plan["expected_unknown_rows"]):
        raise ValueError("audited source cohort differs")
    panels = []
    for config in plan["panels"]:
        splits = read_frame(repo / config["splits"])
        splits = splits[splits.scenario.eq(config["scenario"])].copy()
        if (splits.bacterium.duplicated().any() or set(splits.bacterium) != set(pairs.bacterium)
                or splits[["bacterium", "host_group", "split"]].isna().any().any()
                or splits.groupby("host_group").split.nunique().max() != 1):
            raise ValueError("invalid frozen split/group manifest")
        training = set(splits.loc[splits.split.eq("train"), "bacterium"])
        targets = set(splits.loc[splits.split.isin(config["partitions"]), "bacterium"])
        if not training or training & targets:
            raise ValueError("training and target hosts must be disjoint")
        candidates = pd.MultiIndex.from_product([sorted(targets), phages], names=["bacterium", "phage"]).to_frame(index=False)
        train = observed[observed.bacterium.isin(training)]
        stats = train.groupby("phage").infection.agg(["sum", "count"])
        if set(stats.index) != set(phages):
            raise ValueError("prevalence lacks training support")
        rates = (stats["sum"] + 1) / (stats["count"] + 2)
        prevalence = candidates.assign(score=candidates.phage.map(rates))
        full = {"phage_prevalence": prevalence, "random": candidates.assign(score=0.),
                "genophi_median": read_frame(repo / config["genophi_scores"])}
        for scores in full.values():
            check_scores(scores, targets, phages)
        if config["saved_prevalence_scores"]:
            saved = read_frame(repo / config["saved_prevalence_scores"])
            check_scores(saved, targets, phages)
            key = ["bacterium", "phage"]
            if not np.array_equal(saved.sort_values(key).score, prevalence.sort_values(key).score):
                raise ValueError("train-only prevalence differs from saved scores")
        baselines = read_frame(repo / config["baseline_metrics"])
        if config["baseline_scores"]:
            exported = read_frame(repo / config["baseline_scores"])
            if set(exported.scenario) != {config["scenario"]}:
                raise ValueError("baseline score scenarios differ")
            expected_partition = splits.set_index("bacterium").split
            if not exported.bacterium.map(expected_partition).equals(exported.partition):
                raise ValueError("baseline score partitions differ")
            for method in ["phage_prevalence", *plan["full_baseline_methods"]]:
                scores = exported.loc[exported.model.eq(method), ["bacterium", "phage", "score"]]
                check_scores(scores, targets, phages)
                if method == "phage_prevalence":
                    if not np.array_equal(scores.sort_values(["bacterium", "phage"]).score,
                                          prevalence.sort_values(["bacterium", "phage"]).score):
                        raise ValueError("baseline export prevalence differs")
                else:
                    full[method] = scores
        geno_reference = read_frame(repo / config["genophi_metrics"])
        neighbour = read_frame(repo / config["neighbour_metrics"]) if config["neighbour_metrics"] else None
        for partition in config["partitions"]:
            hosts = sorted(splits.loc[splits.split.eq(partition), "bacterium"])
            references = {m: select_reference(baselines, config, partition, m)
                          for m in ["random", *plan["full_baseline_methods"]]}
            references.update({m: select_reference(geno_reference, config, partition, m)
                               for m in ("genophi_median", "phage_prevalence")})
            if neighbour is not None:
                references["core_neighbour"] = select_reference(neighbour, config, partition, "core_neighbour")
            panels.append({"id": config["id"], "partition": partition, "hosts": hosts,
                           "groups": splits.set_index("bacterium").host_group.to_dict(),
                           "observed": observed[observed.bacterium.isin(hosts)].copy(),
                           "scores": {m: s[s.bacterium.isin(hosts)].copy() for m, s in full.items()},
                           "references": references})
    return plan, phages, panels


def evaluate_panel(panel, phages, budgets):
    frames, support = {}, None
    for method, scores in panel["scores"].items():
        summary, frame = evaluate_retrieval(panel["observed"], scores, hosts=panel["hosts"],
                                           candidate_phages=phages, ks=tuple(range(1, max(budgets)+1)))
        check_reference(frame, panel["references"][method])
        frames[method] = budget_metrics(frame, budgets).assign(source_kind="full_scores")
        if support is None:
            support = {k: summary[k] for k in ("held_out_hosts", "complete_panel_hosts", "eligible_hosts",
                "all_negative_hosts", "incomplete_hosts_excluded", "unknown_pairs", "candidate_phages")}
    canonical = frames["phage_prevalence"]
    for method, reference in panel["references"].items():
        if method in frames:
            continue
        check_reference(canonical[["bacterium", "positive_candidates", "eligible"]], reference)
        frames[method] = budget_metrics(reference, budgets).assign(source_kind="published_metrics")
    return frames, support


def summarize_panel(panel, frames, plan):
    curves, contrasts, availability = [], [], []
    budgets = plan["budgets"]
    for method, frame in sorted(frames.items()):
        for k in budgets:
            available = {kind: bool(np.isfinite(frame[f"{kind}@{k}"]).all()) for kind in ("hit", "tests")}
            availability.append({"panel": panel["id"], "partition": panel["partition"], "model": method,
                                 "budget": k, "source_kind": frame.source_kind.iloc[0], **available})
        for eligible_only in (False, True):
            denominator = "eligible_complete_hosts" if eligible_only else "all_complete_hosts"
            selected = frame[frame.eligible] if eligible_only else frame
            meta = {"panel": panel["id"], "partition": panel["partition"], "model": method,
                    "denominator": denominator, "hosts": len(selected),
                    "groups": len({panel["groups"][h] for h in selected.bacterium}),
                    "source_kind": frame.source_kind.iloc[0]}
            for k in budgets:
                curves.append({**meta, "budget": k,
                               "hit_probability": float(selected[f"hit@{k}"].mean()) if len(selected) else np.nan,
                               "expected_tests_capped": float(selected[f"tests@{k}"].mean()) if len(selected) else np.nan,
                               "mrr": float(selected.mrr.mean()) if len(selected) else np.nan})
            if method == plan["baseline"]:
                continue
            metrics = [*[f"hit@{k}" for k in budgets], *[f"tests@{k}" for k in budgets], "mrr"]
            metrics = [m for m in metrics if np.isfinite(frame[m]).all()]
            intervals = paired_intervals(frame, frames[plan["baseline"]], panel["groups"], metrics,
                                          eligible_only=eligible_only, **plan["statistics"]["bootstrap"])
            for metric, value in intervals.items():
                interval = value["interval"]
                contrasts.append({**meta, "reference": plan["baseline"], "metric": metric,
                                  "delta_model_minus_reference": value["delta"],
                                  "ci_low": interval[0] if interval else np.nan,
                                  "ci_high": interval[1] if interval else np.nan,
                                  "better_direction": "negative" if metric.startswith("tests@") else "positive",
                                  "valid_resamples": value["valid_resamples"], "status": value["status"]})
    return curves, contrasts, availability


def run(repo, data_root, destination, *, fetch=False):
    repo, destination = Path(repo).resolve(), Path(destination).resolve()
    if destination.exists():
        raise ValueError("output directory already exists")
    verify_implementation(repo)
    plan, phages, panels = load_study(repo, data_root, fetch=fetch)
    per_hosts, curves, contrasts, availability, support = [], [], [], [], []
    for panel in panels:
        frames, counts = evaluate_panel(panel, phages, plan["budgets"])
        support.append({"panel": panel["id"], "partition": panel["partition"], **counts})
        for method, frame in sorted(frames.items()):
            per_hosts.append(frame.assign(panel=panel["id"], partition=panel["partition"], model=method,
                                         host_group=frame.bacterium.map(panel["groups"])))
        c, d, a = summarize_panel(panel, frames, plan)
        curves.extend(c); contrasts.extend(d); availability.extend(a)
    with staged_directory(destination) as stage:
        pd.concat(per_hosts, ignore_index=True).to_csv(stage / "per_host_metrics.csv", index=False)
        for name, rows in (("curves.csv", curves), ("contrasts.csv", contrasts), ("availability.csv", availability)):
            pd.DataFrame(rows).to_csv(stage / name, index=False)
        atomic_write_json({"panels": support}, stage / "support.json")
        report = {"study_id": plan["study_id"], "status": "descriptive_existing_score_audit_complete",
                  "protocol_sha256": sha256_file(repo / PROTOCOL),
                  "implementation_sha256": sha256_file(repo / IMPLEMENTATION),
                  "new_learned_model_fits": 0, "new_provider_requests": 0,
                  "cohort_panels": len(plan["panels"]), "partition_panels": len(panels),
                  "per_host_rows": sum(map(len, per_hosts)), "curve_rows": len(curves),
                  "contrast_rows": len(contrasts), "incomplete_support": support,
                  "interval_scope": plan["statistics"]["interval_scope"],
                  "limitations": plan["limitations"], "source_files_sha256": plan["files_sha256"],
                  "output_sha256": {p.name: sha256_file(p) for p in sorted(stage.iterdir())},
                  "provenance": code_provenance(repo, repo / SCRIPT),
                  "versions": {n: importlib.metadata.version(n) for n in ("numpy", "pandas")}}
        atomic_write_json(report, stage / "report.json")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/screening_budget"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Fetch only the two locked source-label tables")
    args = parser.parse_args()
    result = run(Path(__file__).resolve().parents[1], args.data_root, args.outdir, fetch=args.fetch)
    print({k: result[k] for k in ("status", "partition_panels", "curve_rows", "contrast_rows", "new_learned_model_fits")})


if __name__ == "__main__":
    main()
