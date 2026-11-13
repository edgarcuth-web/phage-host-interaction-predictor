#!/usr/bin/env python3
"""Run the frozen new-host, 94-phage development comparison on one CPU thread."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import importlib.metadata
import json
from pathlib import Path
import re
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file, staged_directory
from phage_cocktail.host_baseline import fit_host_model, predict_host_model
from phage_cocktail.host_genome import CANONICAL_WORDS
from phage_cocktail.host_study import host_features, verify_resources
from phage_cocktail.paired_interaction import fit_entity_scaler, fit_interaction, predict_logits, transform_entities
from phage_cocktail.retrieval import evaluate_retrieval
from phage_cocktail.retrieval_protocol import PROTOCOL, PROTOCOL_SHA256, check_protocol
from scripts.run_screening_budget_audit import budget_metrics, check_scores, paired_intervals

IMPLEMENTATION = "data/studies/new_host_retrieval.implementation.v1.json"
METHODS = ("random", "phage_prevalence", "annotation", "host_4mer", "genophi", "phage_only", "additive", "bilinear")
KS = (1, 2, 3, 4, 5)
TOLERANCE = 1e-12


def frame(path):
    return pd.read_csv(path, dtype={"bacterium": str, "phage": str, "entity_id": str, "host_group": str},
                       float_precision="round_trip")


def save_frame(table, path):
    payload = table.to_csv(index=False, float_format="%.17g").encode()
    path.write_bytes(gzip.compress(payload, mtime=0) if path.suffix == ".gz" else payload)


def check_environment(plan):
    actual = {name: sys.version.split()[0] if name == "python" else importlib.metadata.version(name)
              for name in plan["numerical_environment"]}
    if actual != plan["numerical_environment"]:
        raise ValueError(f"numerical environment differs from the protocol: {actual}")
    return actual


def check_implementation(root=ROOT):
    lock = read_json_object(root / IMPLEMENTATION)
    if lock["protocol_sha256"] != PROTOCOL_SHA256 or sha256_file(root / PROTOCOL) != PROTOCOL_SHA256:
        raise ValueError("implementation targets a different protocol")
    for name, digest in lock["files_sha256"].items():
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()) or sha256_file(path) != digest:
            raise ValueError(f"changed pre-fit implementation: {name}")
    return lock


def assert_selection_scope(fit_hosts, validation_hosts, forbidden_hosts):
    fit_hosts, validation_hosts, forbidden_hosts = map(set, (fit_hosts, validation_hosts, forbidden_hosts))
    if (not fit_hosts or not validation_hosts or fit_hosts & validation_hosts
            or (fit_hosts | validation_hosts) & forbidden_hosts):
        raise ValueError("fit/selection hosts overlap or include outer-test hosts")


def choose_trial(trials):
    """Frozen hierarchy: higher Hit@2, lower two-test cost, stronger L2."""
    if not trials:
        raise ValueError("no validation trials")
    best = None
    for trial in trials:
        if not all(np.isfinite(trial[k]) for k in ("hit2", "tests2", "l2")):
            raise ValueError("nonfinite validation selection metric")
        if best is None:
            best = trial
            continue
        hit = trial["hit2"] - best["hit2"]
        cost = trial["tests2"] - best["tests2"]
        if (hit > TOLERANCE or abs(hit) <= TOLERANCE and
                (cost < -TOLERANCE or abs(cost) <= TOLERANCE and trial["l2"] > best["l2"])):
            best = trial
    return best


def evaluate_scores(observed, scores, hosts, phages):
    held = observed[observed.bacterium.isin(hosts)]
    summary, per_host = evaluate_retrieval(held, scores, hosts=sorted(hosts), candidate_phages=phages, ks=KS)
    return summary, budget_metrics(per_host, KS)


def pooled_encoder(rows, raw, digests):
    return {role: fit_entity_scaler(raw[role], rows[column], digests[role])
            for role, column in (("host", "bacterium"), ("phage", "phage"))}


def matrices(raw, encoder):
    return {role: transform_entities(raw[role], encoder[role]) for role in ("host", "phage")}


def fit_pooled(rows, transformed, kind, l2, optimizer):
    hids, pids = sorted(set(rows.bacterium)), sorted(set(rows.phage))
    H, P = transformed["host"].loc[hids], transformed["phage"].loc[pids]
    return fit_interaction(H.to_numpy(), P.to_numpy(), H.index.get_indexer(rows.bacterium),
                           P.index.get_indexer(rows.phage), rows.infection.to_numpy(), np.ones(len(rows)),
                           kind=kind, l2=l2, optimizer=optimizer)


def score_pooled(model, transformed, hosts, phages):
    hosts = sorted(hosts)
    z = predict_logits(model, transformed["host"].loc[hosts].to_numpy(), transformed["phage"].loc[phages].to_numpy())
    panel = pd.MultiIndex.from_product([hosts, phages], names=["bacterium", "phage"]).to_frame(index=False)
    return panel.assign(score=z.ravel())


def macro_rows(per_host, context):
    rows = []
    for denominator, eligible in (("all_complete", False), ("eligible", True)):
        chosen = per_host[per_host.eligible] if eligible else per_host
        metrics = [*[f"hit@{k}" for k in KS], "mrr", *[f"tests@{k}" for k in KS]]
        if eligible:
            metrics += [f"recall@{k}" for k in KS]
        rows.append(dict(**context, denominator=denominator, hosts=len(chosen),
                         **{m: float(chosen[m].mean()) if len(chosen) else None for m in metrics}))
    return rows


def contrast_rows(per_host, groups, context, plan):
    results = []
    config = plan["uncertainty"]
    comparators = ["phage_prevalence", "host_4mer", "genophi", "annotation", "additive", "phage_only"]
    left = per_host[per_host.method.eq("bilinear")].copy()
    for comparator in comparators:
        right = per_host[per_host.method.eq(comparator)].copy()
        for denominator, eligible in (("all_complete", False), ("eligible", True)):
            metrics = [*[f"hit@{k}" for k in KS], "mrr", *[f"tests@{k}" for k in KS]]
            a, b = left.copy(), right.copy()
            if eligible:
                metrics += [f"recall@{k}" for k in KS]
                # All-negative hosts remain in the bootstrap universe, then are masked.
                # Their undefined Recall is never reported as zero or included in means.
                for table in (a, b):
                    table.loc[~table.eligible, [f"recall@{k}" for k in KS]] = 0.
            contrasts = paired_intervals(a, b, groups, metrics, eligible_only=eligible,
                                         resamples=config["resamples"], seed=config["seed"], confidence=config["confidence"])
            for metric, value in contrasts.items():
                interval = value["interval"]
                results.append(dict(**context, left="bilinear", right=comparator, denominator=denominator,
                                    metric=metric, delta=value["delta"], lower=None if interval is None else interval[0],
                                    upper=None if interval is None else interval[1], valid_resamples=value["valid_resamples"],
                                    status=value["status"]))
    return results


def per_phage_rows(scores, observed, context):
    """Observed-row diagnostics retain unknown slots and every scored host."""
    joined = scores.merge(observed[["bacterium", "phage", "infection"]], on=["bacterium", "phage"], how="left", validate="one_to_one")
    pieces = []
    for _, group in joined.groupby("bacterium", sort=True):
        group = group.copy()
        values = group.score.to_numpy()
        above = (values[:, None] > values[None, :]).sum(axis=0)
        ties = (values[:, None] == values[None, :]).sum(axis=0)
        group["top2_slot"] = np.clip(2 - above, 0, ties) / ties
        pieces.append(group)
    table = pd.concat(pieces, ignore_index=True)
    result = []
    for phage, group in table.groupby("phage", sort=True):
        result.append(dict(**context, phage=phage, scored_hosts=len(group),
                           observed_positives=int(group.infection.eq(1).sum()), observed_negatives=int(group.infection.eq(0).sum()),
                           unknown=int(group.infection.isna().sum()), expected_top2_slots=float(group.top2_slot.sum()),
                           positive_top2_slots=float(group.loc[group.infection.eq(1), "top2_slot"].sum()),
                           negative_top2_slots=float(group.loc[group.infection.eq(0), "top2_slot"].sum()),
                           unknown_top2_slots=float(group.loc[group.infection.isna(), "top2_slot"].sum())))
    return result


def run(out, annotation_data_root, *, fetch=False, code_ref, root=ROOT):
    root, out = Path(root).resolve(), Path(out).resolve()
    if out.exists():
        raise ValueError("output directory must be new")
    if not re.fullmatch(r"[0-9a-f]{40}", code_ref):
        raise ValueError("code-ref must identify the public pre-fit implementation commit")
    started = time.monotonic()
    implementation = check_implementation(root)
    plan = read_json_object(root / PROTOCOL)
    versions = check_environment(plan)
    resources = verify_resources({"resources": plan["annotation_resources"]}, annotation_data_root, fetch=fetch)
    preflight = check_protocol(root, annotation_data_root=annotation_data_root)
    paths = {role: root / value["path"] for role, value in plan["inputs"].items()}
    observed, splits = frame(paths["observations"]), frame(paths["scenario_splits"])
    phages = sorted(frame(paths["library"]).phage)
    features = frame(paths["sequence_vectors"]).set_index(["entity_kind", "entity_id"], verify_integrity=True)
    receipts = frame(paths["sequence_manifest"]).set_index(["entity_kind", "entity_id"], verify_integrity=True)
    raw = {"host": features.loc["host_genome", CANONICAL_WORDS].sort_index(),
           "phage": features.loc["phage_genome", CANONICAL_WORDS].sort_index()}
    digests = {"host": receipts.loc["host_genome"].assembly_sha256.to_dict(),
               "phage": receipts.loc["phage_genome"].assembly_sha256.to_dict()}
    metadata = pd.read_csv(resources["table", "metadata"], sep=";", index_col=0)
    defenses = pd.read_csv(resources["table", "defenses"], sep=";", index_col=0)
    categorical = ["O-type", "H-type", "LPS_type", "ABC_serotype"]
    annotations, annotation_columns = host_features(sorted(raw["host"].index), metadata, defenses, categorical)
    fourmers = raw["host"].rename(columns=lambda word: "host_4mer:" + word)
    config = plan["methods"]["pooled_sequence"]
    archive, all_models, encoders, selected, trials, selection_predictions = {}, {}, {}, {}, [], []
    per_phage_models, fitting_counts, final_scores = {}, Counter(), {}

    def store_model(model, key, encoder_key, train):
        archive[key] = model["theta"]
        all_models[key] = {**{k: v for k, v in model.items() if k != "theta"}, "encoder_key": encoder_key,
                           "fit_hosts": sorted(set(train.bacterium)), "fit_phages": sorted(set(train.phage)), "fit_rows": len(train)}
        fitting_counts["pooled"] += 1

    with staged_directory(out) as stage, threadpool_limits(limits=plan["compute"]["cpu_threads"]):
        for scenario, prefix in (("grouped_host", "primary"), ("phylogroup:B2", "b2")):
            split = splits[splits.scenario.eq(scenario)]
            train_ids = set(split.loc[split.split.eq("train"), "bacterium"])
            test_ids = set(split.loc[split.split.eq("test"), "bacterium"])
            outer_train = observed[observed.bacterium.isin(train_ids)].copy()
            if scenario == "grouped_host":
                fit_ids = train_ids
                validation_ids = set(split.loc[split.split.eq("validation"), "bacterium"])
                targets = validation_ids | test_ids
            else:
                inner = frame(paths["b2_selection"])
                fit_ids = set(inner.loc[inner.selection_split.eq("fit"), "bacterium"])
                validation_ids = set(inner.loc[inner.selection_split.eq("validation"), "bacterium"])
                if fit_ids | validation_ids != train_ids:
                    raise ValueError("B2 inner scope differs from outer training hosts")
                targets = test_ids
            assert_selection_scope(fit_ids, validation_ids, test_ids)
            fit = observed[observed.bacterium.isin(fit_ids)].copy()
            validation = observed[observed.bacterium.isin(validation_ids)].copy()
            encoder_key = prefix + "_selection"
            encoders[encoder_key] = pooled_encoder(fit, raw, digests)
            transformed = matrices(raw, encoders[encoder_key])
            chosen = {}
            for kind in config["kinds"]:
                options = []
                for l2 in config["l2_grid_in_tie_preference_order"]:
                    model = fit_pooled(fit, transformed, kind, l2, config["optimizer"])
                    key = f"{prefix}_{kind}_{str(l2).replace('.', 'p')}"
                    store_model(model, key, encoder_key, fit)
                    scores = score_pooled(model, transformed, validation_ids, phages)
                    summary, metrics = evaluate_scores(validation, scores, validation_ids, phages)
                    record = dict(panel=scenario, kind=kind, l2=l2, weights_key=key,
                                  hit2=summary["all_hosts"]["hit@2"], tests2=float(metrics["tests@2"].mean()),
                                  validation_hosts=len(metrics), fit_rows=len(fit), iterations=model["iterations"],
                                  gradient_infinity_norm=model["gradient_infinity_norm"])
                    trials.append(record); options.append(record)
                    selection_predictions.append(scores.assign(panel=scenario, kind=kind, l2=l2, weights_key=key))
                    print(prefix, kind, "L2", l2, "validation Hit@2", record["hit2"], "cost@2", record["tests2"], flush=True)
                best = choose_trial(options)
                chosen[kind] = best["weights_key"]
            if scenario == "phylogroup:B2":
                encoder_key = "b2_outer_refit"
                encoders[encoder_key] = pooled_encoder(outer_train, raw, digests)
                transformed = matrices(raw, encoders[encoder_key])
                for kind, selection_key in list(chosen.items()):
                    l2 = all_models[selection_key]["l2"]
                    model = fit_pooled(outer_train, transformed, kind, l2, config["optimizer"])
                    key = f"b2_outer_{kind}"
                    store_model(model, key, encoder_key, outer_train)
                    chosen[kind] = key
            selected[scenario] = chosen
            panel = pd.MultiIndex.from_product([sorted(targets), phages], names=["bacterium", "phage"]).to_frame(index=False)
            stats = outer_train.groupby("phage").infection.agg(["sum", "count"])
            rates = (stats["sum"] + 1) / (stats["count"] + 2)
            scores_by_method = {"random": panel.assign(score=0.), "phage_prevalence": panel.assign(score=panel.phage.map(rates))}
            if scenario == "grouped_host":
                annotation_model = read_json_object(paths["annotation_weights"])
            else:
                annotation_model = fit_host_model(outer_train, annotations, annotation_columns, categorical, C=.1, max_iter=2000)
                per_phage_models["b2_annotation"] = annotation_model
                fitting_counts["per_phage_slots"] += len(annotation_model["phages"])
                fitting_counts["logistic"] += sum(m["kind"] == "logistic" for m in annotation_model["phages"].values())
            scores_by_method["annotation"] = panel.assign(score=predict_host_model(annotation_model, panel, annotations))
            host_model = fit_host_model(outer_train, fourmers, fourmers.columns.tolist(), [], C=.1, max_iter=2000)
            per_phage_models[prefix + "_host_4mer"] = host_model
            fitting_counts["per_phage_slots"] += len(host_model["phages"])
            fitting_counts["logistic"] += sum(m["kind"] == "logistic" for m in host_model["phages"].values())
            scores_by_method["host_4mer"] = panel.assign(score=predict_host_model(host_model, panel, fourmers))
            scores_by_method["genophi"] = frame(paths["genophi_primary_scores" if scenario == "grouped_host" else "genophi_b2_scores"])
            for kind, key in chosen.items():
                model = {**all_models[key], "theta": archive[key]}
                scores_by_method[kind] = score_pooled(model, transformed, targets, phages)
            for method in METHODS:
                check_scores(scores_by_method[method], targets, phages)
                partition = split.set_index("bacterium").split
                final_scores[scenario, method] = scores_by_method[method].assign(panel=scenario, partition=lambda x: x.bacterium.map(partition), method=method)

        if fitting_counts["pooled"] != plan["compute"]["max_pooled_fits"] or fitting_counts["per_phage_slots"] != plan["compute"]["max_per_phage_logistic_fits"]:
            raise ValueError("fit count differs from the frozen budget")
        # All selection, fitted weights and complete target score matrices now exist.
        # Outer-test labels enter metric evaluation only below this point.
        np.savez_compressed(stage / "weights.npz", **archive)
        host_paths = {}
        for name, model in per_phage_models.items():
            filename = f"{name}.json"
            atomic_write_json(model, stage / filename)
            host_paths[name] = filename
        model_info = dict(feature_columns=CANONICAL_WORDS, encoders=encoders, pooled_models=all_models, selected=selected,
                          host_models=host_paths, annotation_primary=plan["inputs"]["annotation_weights"],
                          weights_sha256=sha256_file(stage / "weights.npz"), score_semantics="uncalibrated source-MLC scores; pooled models use raw logits")
        atomic_write_json(model_info, stage / "models.json")
        save_frame(pd.DataFrame(trials), stage / "validation_trials.csv")
        save_frame(pd.concat(selection_predictions, ignore_index=True), stage / "selection_scores.csv.gz")
        complete_scores = pd.concat(final_scores.values(), ignore_index=True)
        save_frame(complete_scores, stage / "scores.csv.gz")

        per_host_frames, macros, contrasts, support, subgroups, phage_diagnostics, display = [], [], [], [], [], [], []
        for scenario in ("grouped_host", "phylogroup:B2"):
            split = splits[splits.scenario.eq(scenario)]
            groups = split.set_index("bacterium").host_group.to_dict()
            phylogroups = split.set_index("bacterium").phylogroup.to_dict()
            for partition in ("validation", "test"):
                ids = set(split.loc[split.split.eq(partition), "bacterium"])
                if not ids:
                    continue
                local_frames = []
                for method in METHODS:
                    scores = final_scores[scenario, method]
                    scores = scores[scores.partition.eq(partition)][["bacterium", "phage", "score"]]
                    summary, metrics = evaluate_scores(observed, scores, ids, phages)
                    context = dict(panel=scenario, partition=partition, method=method)
                    support.append(dict(**context, **summary))
                    macros += macro_rows(metrics, context)
                    metrics = metrics.assign(**context, host_group=lambda x: x.bacterium.map(groups), phylogroup=lambda x: x.bacterium.map(phylogroups))
                    local_frames.append(metrics)
                    for phylogroup, subset in metrics.groupby("phylogroup", sort=True):
                        subgroups.append(dict(**context, phylogroup=phylogroup, hosts=len(subset), eligible_hosts=int(subset.eligible.sum()),
                                              hit2=float(subset["hit@2"].mean()), tests2=float(subset["tests@2"].mean()), mrr=float(subset.mrr.mean())))
                    phage_diagnostics += per_phage_rows(scores, observed[observed.bacterium.isin(ids)], context)
                    top = scores.sort_values(["bacterium", "score", "phage"], ascending=[True, False, True], kind="stable").groupby("bacterium", sort=False).head(5).copy()
                    top["display_rank"] = top.groupby("bacterium", sort=False).cumcount() + 1
                    display.append(top.assign(**context))
                local = pd.concat(local_frames, ignore_index=True)
                per_host_frames.append(local)
                contrasts += contrast_rows(local, groups, dict(panel=scenario, partition=partition), plan)
        per_host = pd.concat(per_host_frames, ignore_index=True)
        comparison = pd.DataFrame(contrasts)
        save_frame(per_host, stage / "per_host_metrics.csv")
        save_frame(pd.DataFrame(macros), stage / "macro_metrics.csv")
        save_frame(comparison, stage / "contrasts.csv")
        save_frame(pd.DataFrame(subgroups), stage / "subgroups.csv")
        save_frame(pd.DataFrame(phage_diagnostics), stage / "per_phage_diagnostics.csv")
        save_frame(pd.concat(display, ignore_index=True), stage / "top_candidates.csv.gz")
        atomic_write_json(dict(panels=support, display_ties="Lexicographic phage ID for named top-five inspection only; expected metrics average exact ties.",
                               per_phage_scope="All scored hosts with observed labels/unknown slots explicit; observed-row diagnostics, not full-panel performance estimates."), stage / "support.json")
        primary = comparison[comparison.panel.eq("grouped_host") & comparison.partition.eq("test") & comparison.right.eq("phage_prevalence") & comparison.denominator.eq("all_complete")].set_index("metric")
        hit_gain, cost_gain = bool(primary.loc["hit@2", "lower"] > 0), bool(primary.loc["tests@2", "upper"] < 0)
        peak_kib = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        if sys.platform == "linux" and peak_kib > plan["compute"]["max_ram_gb"] * 1024**2:
            raise ValueError("peak RSS exceeds the frozen RAM budget")
        report = dict(study_id=plan["study_id"], status="completed_inspected_development_only", no_model_promotion=True,
                      decisions=dict(primary_hit_gain_supported=hit_gain, primary_cost_saving_supported=cost_gain, joint_hit_and_cost_benefit_supported=hit_gain and cost_gain),
                      primary_contrasts=primary.reset_index().to_dict("records"), macro_metrics=macros,
                      fit_counts=dict(fitting_counts), selected_l2={p: {kind: all_models[key]["l2"] for kind, key in kinds.items()} for p, kinds in selected.items()},
                      complete_score_rows=len(complete_scores), selection_score_rows=sum(len(x) for x in selection_predictions), per_host_metric_rows=len(per_host),
                      reused_feature_vectors=preflight["verified_sequence_vectors"], new_genome_downloads=0, new_genophi_fits=0, new_paid_provider_requests=0,
                      compute=dict(wall_seconds=time.monotonic()-started, peak_rss=peak_kib, peak_rss_unit="KiB on Linux", cpu_threads=1),
                      provenance=dict(protocol_sha256=PROTOCOL_SHA256, implementation_sha256=sha256_file(root / IMPLEMENTATION),
                                      code_ref=code_ref, code_ref_source="caller-supplied public pre-fit Git commit; local code bytes verified against the implementation manifest",
                                      files_sha256=implementation["files_sha256"], numerical_environment=versions, source_protocol_commit="f4ec41a7114d6bb1aa241cb34ceb46cf81871cd7"),
                      limitations=["Inspected development data; no independent release validation.", "Source MLC positivity is not productive infection or cocktail efficacy.",
                                   "GenoPHI genome vocabulary is transductive; inherited artifacts retain their source environments.", "All 94 source-named candidates remain, including exact assembly duplicates.",
                                   "Hit probability and capped testing consumption have separate criteria; secondary intervals are pointwise and unadjusted."],
                      artifacts_sha256={p.name: sha256_file(p) for p in sorted(stage.iterdir())})
        atomic_write_json(report, stage / "report.json")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--annotation-data-root", type=Path, default=ROOT / "artifacts/annotation_inputs")
    parser.add_argument("--fetch", action="store_true", help="fetch only the two hash-locked annotation tables if missing")
    parser.add_argument("--code-ref", required=True, help="public pre-fit implementation commit SHA")
    args = parser.parse_args()
    try:
        report = run(args.out, args.annotation_data_root, fetch=args.fetch, code_ref=args.code_ref)
    except Exception as error:
        failure = args.out.with_name(args.out.name + ".failure.json")
        if not args.out.exists() and not failure.exists():
            atomic_write_json(dict(status="failed_no_admitted_result", error_type=type(error).__name__, reason=str(error),
                                   protocol_sha256=PROTOCOL_SHA256, code_ref=args.code_ref), failure)
        print(f"study failed: {error}", file=sys.stderr)
        return 2
    print(json.dumps(dict(decisions=report["decisions"], selected_l2=report["selected_l2"], fit_counts=report["fit_counts"]), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
