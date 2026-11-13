#!/usr/bin/env python3
"""Independently check the matched study's weights, scores, ranks and decisions."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]

import numpy as np
import pandas as pd
from scipy.special import expit
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.retrieval_protocol import PROTOCOL, PROTOCOL_SHA256, check_protocol
from scripts.run_new_host_retrieval import IMPLEMENTATION, METHODS, check_environment, check_implementation, frame


def close(actual, expected, atol=1e-10):
    if not np.isfinite(np.asarray(actual)).all() or not np.allclose(actual, expected, rtol=0, atol=atol, equal_nan=False):
        raise ValueError("independent numerical check failed")


def independent_metrics(labels, scores, ks=(1, 2, 3, 4, 5)):
    """No-hit hypergeometric products and exact fractions, independent of retrieval.py."""
    groups = {}
    for phage, score in scores.items():
        groups.setdefault(score, []).append(int(labels[phage]))
    ordered = [groups[s] for s in sorted(groups, reverse=True)]
    positives = sum(labels.values())
    result = {}
    for k in ks:
        remaining, no_hit, found = k, Fraction(1), Fraction(0)
        for values in ordered:
            n, p = len(values), sum(values)
            take = min(remaining, n)
            if take:
                no_hit *= Fraction(comb(n - p, take), comb(n, take)) if take <= n - p else 0
                found += Fraction(take * p, n)
            remaining -= take
            if not remaining:
                break
        result[f"hit@{k}"] = float(1 - no_hit)
        result[f"recall@{k}"] = float(found / positives) if positives else np.nan
    reciprocal, offset = Fraction(0), 0
    for values in ordered:
        n, p = len(values), sum(values)
        if p:
            reciprocal = sum((Fraction(comb(n - j, p - 1), comb(n, p) * (offset + j))
                              for j in range(1, n - p + 2)), Fraction(0))
            break
        offset += n
    result["mrr"] = float(reciprocal)
    for k in ks:
        result[f"tests@{k}"] = k - sum(result[f"hit@{j}"] for j in range(1, k))
    return result


def matrix_scores(theta, kind, H, P):
    dh, dp, position = H.shape[1], P.shape[1], 1
    a = np.zeros(dh)
    if kind != "phage_only":
        a, position = theta[position:position + dh], position + dh
    b, position = theta[position:position + dp], position + dp
    expected = position + (dh * dp if kind == "bilinear" else 0)
    if len(theta) != expected:
        raise ValueError("pooled coefficient shape differs")
    W = theta[position:].reshape(dh, dp) if kind == "bilinear" else np.zeros((dh, dp))
    return theta[0] + (H @ a)[:, None] + (P @ b)[None, :] + H @ W @ P.T


def bootstrap_expanded(left, right, groups, metrics, *, eligible_only, config):
    """Expand each sampled group's hosts rather than multiplying grouped sums."""
    a, b = left.set_index("bacterium").sort_index(), right.set_index("bacterium").sort_index()
    if list(a.index) != list(b.index) or not a.eligible.equals(b.eligible):
        raise ValueError("paired bootstrap support differs")
    names = sorted({groups[h] for h in a.index})
    ids = np.array([groups[h] for h in a.index])
    positions = [np.flatnonzero(ids == name) for name in names]
    difference = a[metrics].to_numpy(float) - b[metrics].to_numpy(float)
    keep = a.eligible.to_numpy(bool) if eligible_only else np.ones(len(a), bool)
    if not np.isfinite(difference[keep]).all():
        raise ValueError("undefined metric on an included host")
    rng, estimates = np.random.default_rng(config["seed"]), []
    for _ in range(config["resamples"]):
        indices = np.concatenate([positions[i] for i in rng.integers(len(names), size=len(names))])
        indices = indices[keep[indices]]
        if len(indices):
            estimates.append(difference[indices].mean(axis=0))
    alpha = (1 - config["confidence"]) / 2
    return difference[keep].mean(axis=0), np.quantile(estimates, [alpha, 1 - alpha], axis=0), len(estimates)


def annotation_features(hosts, metadata, defenses, preprocessing):
    """Join published source names independently of the production feature builder."""
    def normalized(value):
        return str(value).upper().replace("-", "")
    lookups = []
    for table in (metadata, defenses):
        lookup = {normalized(name): name for name in table.index}
        if len(lookup) != len(table):
            raise ValueError("ambiguous annotation alias")
        lookups.append(lookup)
    rows = []
    for host in hosts:
        m = metadata.loc[lookups[0][normalized(host)]]
        d = defenses.loc[lookups[1][normalized(host)]]
        rows.append({**{c: d[c.removeprefix("defense:")] for c in preprocessing["numeric_columns"]},
                     **{c: m[c] for c in preprocessing["categorical_columns"]}})
    return pd.DataFrame(rows, index=hosts)


def verify(bundle, annotation_data_root, *, root=ROOT):
    root, bundle = Path(root).resolve(), Path(bundle).resolve()
    plan = read_json_object(root / PROTOCOL)
    check_environment(plan); implementation = check_implementation(root)
    check_protocol(root, annotation_data_root=annotation_data_root)
    report = read_json_object(bundle / "report.json")
    if report["provenance"]["protocol_sha256"] != PROTOCOL_SHA256 or report["provenance"]["implementation_sha256"] != sha256_file(root / IMPLEMENTATION):
        raise ValueError("run protocol/implementation identity differs")
    if report["provenance"]["files_sha256"] != implementation["files_sha256"]:
        raise ValueError("run implementation file list differs")
    for name, digest in report["artifacts_sha256"].items():
        path = (bundle / name).resolve()
        if not path.is_relative_to(bundle) or sha256_file(path) != digest:
            raise ValueError(f"changed output: {name}")
    paths = {role: root / value["path"] for role, value in plan["inputs"].items()}
    observed, splits, inner = frame(paths["observations"]), frame(paths["scenario_splits"]), frame(paths["b2_selection"])
    features = frame(paths["sequence_vectors"]).set_index(["entity_kind", "entity_id"])
    receipts = frame(paths["sequence_manifest"]).set_index(["entity_kind", "entity_id"])
    model_info = read_json_object(bundle / "models.json")
    if model_info["feature_columns"] != features.columns.tolist():
        raise ValueError("fitted feature schema differs")
    raw = {"host": features.loc["host_genome"].sort_index(), "phage": features.loc["phage_genome"].sort_index()}
    phages = sorted(raw["phage"].index)
    scope_hosts = {
        "primary_selection": set(splits.loc[splits.scenario.eq("grouped_host") & splits.split.eq("train"), "bacterium"]),
        "b2_selection": set(inner.loc[inner.selection_split.eq("fit"), "bacterium"]),
        "b2_outer_refit": set(splits.loc[splits.scenario.eq("phylogroup:B2") & splits.split.eq("train"), "bacterium"])}
    transformed = {}
    for scope, enc in model_info["encoders"].items():
        transformed[scope] = {}
        for role, kind in (("host", "host_genome"), ("phage", "phage_genome")):
            ids = sorted(scope_hosts[scope] if role == "host" else phages)
            unique = {}
            for name in ids:
                unique.setdefault(receipts.loc[(kind, name), "assembly_sha256"], name)
            values = raw[role].loc[list(unique.values())].to_numpy(float)
            means, scales = values.mean(axis=0), values.std(axis=0)
            scales[scales == 0] = 1.
            close(enc[role]["means"], means, 1e-14); close(enc[role]["scales"], scales, 1e-14)
            if enc[role]["unique_fit_assemblies"] != len(unique):
                raise ValueError("unique-fit scaler scope differs")
            transformed[scope][role] = pd.DataFrame((raw[role].to_numpy(float)-means)/scales,
                                                    index=raw[role].index, columns=raw[role].columns)
    if len(model_info["pooled_models"]) != 21 or sha256_file(bundle / "weights.npz") != model_info["weights_sha256"]:
        raise ValueError("pooled fit archive differs")
    archive = np.load(bundle / "weights.npz", allow_pickle=False)
    if set(archive.files) != set(model_info["pooled_models"]):
        raise ValueError("pooled archive keys differ from model metadata")
    pooled_matrices = {}
    maximum_gradient = 0.
    for key, model in model_info["pooled_models"].items():
        scope = model["encoder_key"]
        if set(model["fit_hosts"]) != scope_hosts[scope] or set(model["fit_phages"]) != set(phages):
            raise ValueError("model training host/phage scope differs")
        fit = observed[observed.bacterium.isin(scope_hosts[scope])]
        if len(fit) != model["fit_rows"]:
            raise ValueError("model training row count differs")
        H, P = transformed[scope]["host"], transformed[scope]["phage"]
        theta = archive[key]
        z = matrix_scores(theta, model["kind"], H.to_numpy(), P.to_numpy())
        pooled_matrices[key] = pd.DataFrame(z, index=H.index, columns=P.index)
        h, p = H.loc[fit.bacterium].to_numpy(), P.loc[fit.phage].to_numpy()
        logits = z[H.index.get_indexer(fit.bacterium), P.index.get_indexer(fit.phage)]
        y = fit.infection.to_numpy()
        residual = (expit(logits) - y) / len(fit)
        blocks = [np.array([residual.sum()])]
        if model["kind"] != "phage_only":
            blocks.append(h.T @ residual)
        blocks.append(p.T @ residual)
        if model["kind"] == "bilinear":
            blocks.append((h.T @ (residual[:, None] * p)).ravel())
        gradient = np.concatenate(blocks); gradient[1:] += model["l2"] * theta[1:]
        norm = float(np.max(np.abs(gradient)))
        maximum_gradient = max(maximum_gradient, norm)
        if norm > plan["methods"]["pooled_sequence"]["optimizer"]["maximum_admitted_gradient_infinity_norm"]:
            raise ValueError("pooled weights fail independent stationarity check")
        objective = np.mean(np.logaddexp(0., logits)-y*logits) + .5*model["l2"]*np.sum(theta[1:]**2)
        close(model["objective"], objective)

    selection = frame(bundle / "selection_scores.csv.gz")
    trials = frame(bundle / "validation_trials.csv")
    if len(trials) != 18 or len(selection) != 39762:
        raise ValueError("validation trial/score budget differs")
    checked_scores = 0
    for row in trials.itertuples():
        targets = set(splits.loc[splits.scenario.eq("grouped_host") & splits.split.eq("validation"), "bacterium"]) if row.panel == "grouped_host" else set(inner.loc[inner.selection_split.eq("validation"), "bacterium"])
        scores = selection[selection.weights_key.eq(row.weights_key)]
        expected_cells = set(product(targets, phages))
        if scores.duplicated(["bacterium", "phage"]).any() or set(map(tuple, scores[["bacterium", "phage"]].to_numpy())) != expected_cells:
            raise ValueError("selection score support differs")
        expected = pooled_matrices[row.weights_key]
        close(scores.score, expected.to_numpy()[expected.index.get_indexer(scores.bacterium), expected.columns.get_indexer(scores.phage)])
        values = []
        for host, panel in scores.groupby("bacterium", sort=True):
            truth = observed[observed.bacterium.eq(host)].set_index("phage").infection.to_dict()
            values.append(independent_metrics(truth, dict(zip(panel.phage, panel.score))))
        close(row.hit2, np.mean([v["hit@2"] for v in values]), 1e-12)
        close(row.tests2, np.mean([v["tests@2"] for v in values]), 1e-12)
        checked_scores += len(scores)
    for (scenario, kind), table in trials.groupby(["panel", "kind"]):
        best = None
        for candidate in table.sort_values("l2", ascending=False).itertuples():
            if (best is None or candidate.hit2 > best.hit2 + 1e-12 or
                    abs(candidate.hit2-best.hit2) <= 1e-12 and
                    (candidate.tests2 < best.tests2-1e-12 or abs(candidate.tests2-best.tests2) <= 1e-12 and candidate.l2 > best.l2)):
                best = candidate
        selected = model_info["pooled_models"][model_info["selected"][scenario][kind]]
        close(selected["l2"], best.l2, 0)

    scores = frame(bundle / "scores.csv.gz")
    if len(scores) != 91744 or scores.duplicated(["panel", "method", "bacterium", "phage"]).any():
        raise ValueError("final score support differs")
    metadata = pd.read_csv(Path(annotation_data_root) / "metadata.csv", sep=";", index_col=0)
    defenses = pd.read_csv(Path(annotation_data_root) / "defenses.csv", sep=";", index_col=0)
    independent_rows = []
    for (scenario, method), table in scores.groupby(["panel", "method"], sort=True):
        split = splits[splits.scenario.eq(scenario)]
        targets = set(split.loc[split.split.isin(["validation", "test"]), "bacterium"])
        if set(map(tuple, table[["bacterium", "phage"]].to_numpy())) != set(product(targets, phages)):
            raise ValueError("final target panel differs")
        partition_map = split.set_index("bacterium").split
        if not table.bacterium.map(partition_map).equals(table.partition):
            raise ValueError("final score partition differs")
        train_ids = set(split.loc[split.split.eq("train"), "bacterium"])
        if method in ("phage_only", "additive", "bilinear"):
            values = pooled_matrices[model_info["selected"][scenario][method]]
            close(table.score, values.to_numpy()[values.index.get_indexer(table.bacterium), values.columns.get_indexer(table.phage)])
        elif method in ("annotation", "host_4mer"):
            name = ("primary" if scenario == "grouped_host" else "b2") + "_" + method
            model = read_json_object(paths["annotation_weights"] if name == "primary_annotation" else bundle / model_info["host_models"][name])
            if set(model["training_hosts"]) != train_ids or set(model["phages"]) != set(phages):
                raise ValueError("per-phage weight training scope differs")
            enc = model["preprocessing"]
            raw_host = raw["host"].rename(columns=lambda c:"host_4mer:"+c) if method == "host_4mer" else annotation_features(sorted(raw["host"].index), metadata, defenses, enc)
            train_values = raw_host.loc[sorted(train_ids), enc["numeric_columns"]].to_numpy(float)
            means, scales = train_values.mean(axis=0), train_values.std(axis=0); scales[scales == 0] = 1.
            close(enc["means"], means, 1e-14); close(enc["scales"], scales, 1e-14)
            pieces = [(raw_host[enc["numeric_columns"]].to_numpy(float)-means)/scales]
            for column in enc["categorical_columns"]:
                values = raw_host[column].fillna("__NOT_REPORTED__").astype(str)
                if sorted(values.loc[sorted(train_ids)].unique()) != enc["categories"][column]:
                    raise ValueError("per-phage category vocabulary differs")
                pieces.append(np.column_stack([(values == c).to_numpy(float) for c in enc["categories"][column]]))
            X = pd.DataFrame(np.column_stack(pieces), index=raw_host.index)
            train = observed[observed.bacterium.isin(train_ids)]
            for phage, group in table.groupby("phage", sort=True):
                weights = model["phages"][phage]
                truth = train[train.phage.eq(phage)]
                if weights["training_rows"] != len(truth):
                    raise ValueError("per-phage training row support differs")
                if weights["kind"] == "single_class_prevalence":
                    expected_score = (truth.infection.sum()+1)/(len(truth)+2)
                    close(group.score, expected_score)
                else:
                    logits = X.loc[group.bacterium].to_numpy() @ np.array(weights["coefficients"]) + weights["intercept"]
                    close(group.score, expit(logits))
        elif method == "genophi":
            source = frame(paths["genophi_primary_scores" if scenario == "grouped_host" else "genophi_b2_scores"]).set_index(["bacterium", "phage"])
            close(table.score, source.loc[pd.MultiIndex.from_frame(table[["bacterium", "phage"]]), "score"], 0)
        elif method == "phage_prevalence":
            train = observed[observed.bacterium.isin(train_ids)]
            counts = train.groupby("phage").infection.agg(["sum", "count"])
            close(table.score, table.phage.map((counts["sum"]+1)/(counts["count"]+2)), 0)
        elif method == "random":
            close(table.score, 0., 0)
        checked_scores += len(table)
        for host, group in table.groupby("bacterium", sort=True):
            truth = observed[observed.bacterium.eq(host)].set_index("phage").infection.to_dict()
            if set(truth) != set(phages):
                continue
            metrics = independent_metrics(truth, dict(zip(group.phage, group.score)))
            independent_rows.append(dict(panel=scenario, partition=group.partition.iloc[0], method=method, bacterium=host,
                                         eligible=bool(sum(truth.values())), positive_candidates=sum(truth.values()), **metrics))
    independent = pd.DataFrame(independent_rows)
    saved = frame(bundle / "per_host_metrics.csv")
    keys = ["panel", "partition", "method", "bacterium"]
    if len(independent) != 968 or len(saved) != 968 or saved.duplicated(keys).any():
        raise ValueError("complete-host metric support differs")
    merged = independent.merge(saved, on=keys, suffixes=("_expected", "_saved"), validate="one_to_one")
    if len(merged) != len(independent):
        raise ValueError("complete-host metric identities differ")
    for metric in [*[f"hit@{k}" for k in range(1,6)], "mrr", *[f"tests@{k}" for k in range(1,6)], *[f"recall@{k}" for k in range(1,6)]]:
        mask = merged.eligible_expected if metric.startswith("recall") else np.ones(len(merged), bool)
        close(merged.loc[mask, metric+"_saved"], merged.loc[mask, metric+"_expected"], 1e-12)
    macro = frame(bundle / "macro_metrics.csv")
    if len(macro) != 48 or macro.duplicated(["panel", "partition", "method", "denominator"]).any():
        raise ValueError("macro metric support differs")
    for row in [*macro.to_dict("records"), *report["macro_metrics"]]:
        group = independent[independent.panel.eq(row["panel"]) & independent.partition.eq(row["partition"]) & independent.method.eq(row["method"])]
        if row["denominator"] == "eligible":
            group = group[group.eligible]
        if len(group) != row["hosts"]:
            raise ValueError("macro denominator differs")
        for metric in [c for c in row if c.startswith(("hit@", "tests@", "recall@")) or c == "mrr"]:
            if metric.startswith("recall") and row["denominator"] != "eligible":
                continue
            close(row[metric], group[metric].mean(), 1e-12)

    contrasts = frame(bundle / "contrasts.csv")
    checked_contrasts = 0
    for (scenario, partition, comparator, denominator), group in contrasts.groupby(["panel", "partition", "right", "denominator"], sort=True):
        base = independent[independent.panel.eq(scenario) & independent.partition.eq(partition)]
        hosts = splits[splits.scenario.eq(scenario)].set_index("bacterium").host_group.to_dict()
        metrics = group.metric.tolist()
        delta, interval, valid = bootstrap_expanded(base[base.method.eq("bilinear")], base[base.method.eq(comparator)], hosts,
                                                    metrics, eligible_only=denominator == "eligible", config=plan["uncertainty"])
        close(group.delta, delta, 1e-12); close(group.lower, interval[0], 1e-12); close(group.upper, interval[1], 1e-12)
        if not group.valid_resamples.eq(valid).all():
            raise ValueError("bootstrap undefined-draw count differs")
        checked_contrasts += len(group)
    primary = contrasts[contrasts.panel.eq("grouped_host") & contrasts.partition.eq("test") & contrasts.right.eq("phage_prevalence") & contrasts.denominator.eq("all_complete")].set_index("metric")
    hit, cost = bool(primary.loc["hit@2", "lower"] > 0), bool(primary.loc["tests@2", "upper"] < 0)
    if report["decisions"] != dict(primary_hit_gain_supported=hit, primary_cost_saving_supported=cost, joint_hit_and_cost_benefit_supported=hit and cost):
        raise ValueError("primary decision differs")
    if not report["no_model_promotion"] or report["fit_counts"]["pooled"] != 21 or report["fit_counts"]["per_phage_slots"] != 282:
        raise ValueError("fit/promotion boundary differs")
    fitted_logistic = sum(m["kind"] == "logistic" for name in model_info["host_models"].values()
                          for m in read_json_object(bundle / name)["phages"].values())
    if fitted_logistic != report["fit_counts"]["logistic"]:
        raise ValueError("per-phage logistic fit count differs")
    archive.close()
    return dict(status="verified", pooled_models=21, maximum_independent_pooled_gradient=maximum_gradient,
                independently_checked_scores=checked_scores, complete_host_metric_rows=len(independent), paired_contrasts=checked_contrasts,
                checks=["frozen protocol/input/code/output identities", "fit-only scalers and every pooled objective/gradient",
                        "independent matrix/JSON model predictions and validation selection", "complete fixed-panel support and unknown-host exclusions",
                        "fraction-based exact-tie retrieval and capped costs", "expanded-host paired group bootstrap and primary decisions"],
                verifier_sha256=sha256_file(Path(__file__)), published_output_sha256={**report["artifacts_sha256"], "report.json":sha256_file(bundle / "report.json")})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--annotation-data-root", type=Path, default=ROOT / "artifacts/annotation_inputs")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("output must be new")
    with threadpool_limits(limits=1):
        report = verify(args.bundle, args.annotation_data_root)
    atomic_write_json(report, args.out)
    print({k:v for k,v in report.items() if k != "published_output_sha256"})


if __name__ == "__main__":
    main()
