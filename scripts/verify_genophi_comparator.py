"""Offline numerical QA for a completed frozen-split GenoPHI comparator run."""
from __future__ import annotations

import argparse
from importlib.metadata import version
from pathlib import Path
import time

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file, staged_directory
from phage_cocktail.cohort import curate_mlc
from phage_cocktail.data import load_interactions
from phage_cocktail.genophi_comparator import feature_columns, selected_occurrences
from phage_cocktail.genophi_proteomes import verify_proteomes
from phage_cocktail.genophi_seed_protocols import (
    load_study_protocol, verify_seed_baseline_reference, assert_baseline_model_split, BASE_COMPARATOR, CACHE_LOCK,
)
from phage_cocktail.host_study import verify_resources, host_features
from phage_cocktail.host_genome import host_genome_features
from phage_cocktail.host_baseline import predict_host_model, HostFeatureEncoder
from phage_cocktail.retrieval import evaluate_retrieval


def verify_artifact_hashes(root, records):
    root = Path(root).resolve()
    verified = set()
    for record in records:
        if record["status"] != "complete":
            raise ValueError("numeric QA requires completed stage receipts")
        for relative, digest in record["output_sha256"].items():
            path = (root / relative).resolve()
            if not path.is_relative_to(root) or not path.is_file() or sha256_file(path) != digest:
                raise ValueError(f"missing or changed receipted artifact: {relative}")
            verified.add(relative)
    return len(verified)


def verify_predicted_scores(model, candidates, saved, *, threads):
    """Reconstruct scores directly from the saved estimator, aligning pair IDs."""
    keys = ["bacterium", "phage"]
    candidate_keys = candidates[["strain", "phage"]].rename(columns={"strain": "bacterium"})
    expected = pd.MultiIndex.from_frame(candidate_keys)
    actual = pd.MultiIndex.from_frame(saved[keys])
    if actual.has_duplicates or set(actual) != set(expected):
        raise ValueError("saved model scores lack the complete candidate panel")
    columns = list(model.feature_names_)
    if not columns or any(not c.startswith(("sc_", "pc_")) for c in columns):
        raise ValueError("saved model uses non-protein features")
    if list(model.classes_) != [0, 1]:
        raise ValueError("saved model lacks binary class support")
    # Do not reuse the adapter prediction helper: call CatBoost independently.
    fresh = np.asarray(model.predict_proba(candidates[columns], thread_count=threads))[:, 1]
    recorded = saved.set_index(keys).score.reindex(expected).to_numpy(dtype=float)
    if not np.isfinite(fresh).all() or not np.array_equal(fresh, recorded):
        raise ValueError("saved candidate scores do not exactly reproduce model inference")
    return fresh


def verify_summary(expected, actual):
    """Compare every support field and metric, allowing only float arithmetic noise."""
    for key in actual:
        if key in {"all_hosts", "eligible_hosts_only"}:
            if set(expected[key]) != set(actual[key]):
                raise ValueError("retrieval metric fields differ")
            for metric, value in actual[key].items():
                previous = expected[key][metric]
                if value is None:
                    if previous is not None:
                        raise ValueError(f"changed undefined metric: {metric}")
                elif previous is None or not np.isclose(value, previous, rtol=0, atol=1e-14):
                    raise ValueError(f"recalculated retrieval metric differs: {key}/{metric}")
        elif expected[key] != actual[key]:
            raise ValueError(f"recalculated retrieval support differs: {key}")


def verify_run(repo, root, data_root, proteomes, proteome_receipt):
    repo, root = Path(repo).resolve(), Path(root).resolve()
    started = time.perf_counter()
    report_path = root / "evaluation/report.json"
    report = read_json_object(report_path)
    provenance, plan = report["provenance"], report["protocol"]
    if report["status"] != "completed_aligned_development_comparator":
        raise ValueError("the comparator has not completed")
    protocol_path, original_plan = load_study_protocol(repo, provenance.get("protocol_path", BASE_COMPARATOR), "comparator")
    if plan != original_plan or read_json_object(root / "protocol.json") != original_plan:
        raise ValueError("changed run protocol")
    if sha256_file(protocol_path) != provenance["protocol_sha256"]:
        raise ValueError("changed frozen comparator protocol")
    erratum_path = repo / "data/studies/genophi_comparator.v1.errata.json"
    if read_json_object(erratum_path)["protocol_sha256"] != sha256_file(repo / BASE_COMPARATOR):
        raise ValueError("metadata erratum belongs to a different protocol")
    if "split_seed" in plan:
        binding = verify_seed_baseline_reference(repo, plan)
        if binding != provenance["matching_seed_baseline"]:
            raise ValueError("matching-seed baseline provenance changed")
        if sha256_file(repo / CACHE_LOCK) != provenance["feature_artifact_lock_sha256"]:
            raise ValueError("changed feature reuse lock")
        imported = read_json_object(root / "checkpoints/reuse_features.json")
        if imported["result"]["new_clustering_search_assignment_calls"] != 0:
            raise ValueError("sensitivity unexpectedly rebuilt protein-family features")
    for relative, digest in {**plan["repo_files_sha256"], **provenance["code_sha256"]}.items():
        if sha256_file(repo / relative) != digest:
            raise ValueError(f"changed frozen input/code: {relative}")
    if version("catboost") != provenance["versions"]["catboost"]:
        raise ValueError("reload models with the recorded CatBoost version")
    lock = read_json_object(repo / plan["resource_lock"])
    resources = verify_resources(lock, data_root)
    if sha256_file(proteome_receipt) != provenance["proteome_receipt_sha256"]:
        raise ValueError("changed proteome receipt")
    admitted = verify_proteomes(read_json_object(proteome_receipt), proteomes, lock)
    observed, excluded, _ = curate_mlc(load_interactions(resources[("table", "subset")]),
        pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0))
    saved_observed = pd.read_csv(root / "cohort/observed_pairs.csv", dtype={"bacterium": str, "phage": str})
    pd.testing.assert_frame_equal(observed, saved_observed)
    pd.testing.assert_frame_equal(excluded, pd.read_csv(root / "cohort/excluded_pairs.csv", dtype=str).astype({"source_binary": int}))
    splits = pd.read_csv(root / "cohort/scenario_splits.csv", dtype=str)
    full = pd.read_csv(root / "full_feature_table.csv", dtype={"strain": str, "phage": str})
    full_features = feature_columns(full)
    host_features = pd.read_csv(root / "clustering/strain/features/feature_table.csv", dtype={"strain": str})
    phage_features = pd.read_csv(root / "clustering/phage/features/feature_table.csv", dtype={"phage": str})
    checkpoints = [read_json_object(p) for p in sorted((root / "checkpoints").glob("*.json"))]
    artifacts = verify_artifact_hashes(root, checkpoints)
    if not any(r["stage"] == "evaluation" for r in checkpoints):
        raise ValueError("missing completed evaluation checkpoint")
    panel = sorted(observed.phage.unique())
    models_verified, evaluated_panels, max_ensemble_error = 0, 0, 0.
    prevalence_verified = []
    identical_feature_repeats = {}
    for scenario in plan["scenarios"]:
        split = splits[splits.scenario.eq(scenario)]
        train_hosts = set(split.loc[split.split.eq("train"), "bacterium"])
        held_hosts = sorted(split.loc[split.split.isin(["validation", "test"]), "bacterium"])
        # Build an independent complete Cartesian product, not observed-only rows.
        candidates = pd.MultiIndex.from_product([held_hosts, panel], names=["strain", "phage"]).to_frame(index=False)
        candidates = candidates.merge(host_features, on="strain", validate="many_to_one").merge(
            phage_features, on="phage", validate="many_to_one")
        feature_columns(candidates)
        train_expected = full[full.strain.isin(train_hosts)].sort_values(["strain", "phage"]).reset_index(drop=True)
        all_values, repeat_values, repeat_feature_sets = [], [], []
        for repeat in range(plan["repeats"]):
            tag = f"{scenario.replace(':', '_')}_repeat_{repeat}"
            fit = root / "fits" / tag
            train = pd.read_csv(fit / "train_full_feature_table.csv", dtype={"strain": str, "phage": str})
            pd.testing.assert_frame_equal(train_expected, train.sort_values(["strain", "phage"]).reset_index(drop=True))
            if set(train.strain) != train_hosts or set(train.strain) & set(held_hosts):
                raise ValueError("outer-training host leakage")
            fs = plan["feature_selection"]
            occurrences = selected_occurrences(fit / "feature_selection", num_runs=fs["num_runs_fs"], allowed_features=full_features)
            chosen = set(occurrences.loc[occurrences.Occurrence.ge(fs["cutoff"]), "Feature"])
            selected = pd.read_csv(fit / "selected" / f"select_feature_table_cutoff_{fs['cutoff']}.csv",
                                   dtype={"strain": str, "phage": str})
            if set(feature_columns(selected)) != chosen or not fs["min_features"] <= len(chosen) < fs["max_features"]:
                raise ValueError("selected features differ from the fixed stratum")
            pd.testing.assert_frame_equal(selected.sort_values(["strain", "phage"]).reset_index(drop=True),
                train_expected[["strain", "phage", "interaction", *sorted(chosen)]], check_dtype=False)
            repeat_feature_sets.append(chosen)
            values = []
            for run in range(plan["modeling"]["num_runs_modeling"]):
                model_dir = fit / "models" / f"run_{run}"
                internal = pd.read_csv(model_dir / "best_model_predictions.csv", dtype={"strain": str, "phage": str})
                if not set(internal.strain).issubset(train_hosts):
                    raise ValueError("internal selection includes an outer held-out host")
                model = joblib.load(model_dir / "best_model.pkl")
                saved = pd.read_csv(model_dir / "candidate_scores.csv", dtype={"bacterium": str, "phage": str}, float_precision="round_trip")
                fresh = verify_predicted_scores(model, candidates, saved, threads=plan["compute"]["threads"])
                values.append(fresh)
                models_verified += 1
            # NumPy median is independent of the runner's pandas groupby median.
            repeat_values.append(np.median(np.stack(values), axis=0))
            all_values.extend(values)
        identical_feature_repeats[scenario] = all(s == repeat_feature_sets[0] for s in repeat_feature_sets)
        ensemble = np.median(np.stack(all_values), axis=0)
        ensemble_file = pd.read_csv(root / "evaluation" / f"{scenario.replace(':', '_')}_candidate_scores.csv",
                                   dtype={"bacterium": str, "phage": str}, float_precision="round_trip")
        keys = candidates[["strain", "phage"]].rename(columns={"strain": "bacterium"})
        expected_index = pd.MultiIndex.from_frame(keys)
        recorded_ensemble = ensemble_file.set_index(["bacterium", "phage"]).score.reindex(expected_index).to_numpy()
        if not np.isfinite(recorded_ensemble).all():
            raise ValueError("ensemble candidate coverage is incomplete")
        error = float(np.max(np.abs(ensemble - recorded_ensemble)))
        if error > 5e-16:
            raise ValueError("ensemble median differs from reloaded model scores")
        max_ensemble_error = max(max_ensemble_error, error)
        counts = observed[observed.bacterium.isin(train_hosts)].groupby("phage").infection.agg(["sum", "count"])
        rates = (counts["sum"] + 1) / (counts["count"] + 2)
        prevalence = keys.phage.map(rates).to_numpy(dtype=float)
        if "split_seed" in plan:
            saved_prevalence = pd.read_csv(root / "evaluation" / f"{scenario.replace(':', '_')}_prevalence_scores.csv",
                                          dtype={"bacterium": str, "phage": str}, float_precision="round_trip")
            np.testing.assert_array_equal(prevalence,
                saved_prevalence.set_index(["bacterium", "phage"]).score.reindex(expected_index).to_numpy())
            baseline_root = (repo / plan["reference_report"]).parent
            baseline_scores = pd.read_csv(baseline_root / "predictions.csv.gz", dtype={"bacterium": str, "phage": str},
                                          float_precision="round_trip")
            baseline_scores = baseline_scores[baseline_scores.scenario.eq(scenario) & baseline_scores.model.eq("phage_prevalence")]
            np.testing.assert_array_equal(prevalence,
                baseline_scores.set_index(["bacterium", "phage"]).score.reindex(expected_index).to_numpy())
            prevalence_verified.append({"scenario": scenario, "split_seed": plan["split_seed"],
                                         "candidate_pairs": len(prevalence), "same_train_only_rates_in_both_studies": True})
        for partition in plan["evaluation"]["partitions"]:
            hosts = sorted(split.loc[split.split.eq(partition), "bacterium"])
            if not hosts:
                continue
            mask = keys.bacterium.isin(hosts).to_numpy()
            held = observed[observed.bacterium.isin(hosts)]
            for model_name, probabilities in [("genophi_median", ensemble), ("phage_prevalence", prevalence),
                    *[(f"genophi_repeat_{r}", value) for r, value in enumerate(repeat_values)]]:
                scores = keys.loc[mask].assign(score=probabilities[mask])
                calculated, _ = evaluate_retrieval(held, scores, hosts=hosts, candidate_phages=panel, ks=plan["evaluation"]["k"])
                saved_summary = next(m for m in report["metrics"] if m["scenario"] == scenario
                                     and m["partition"] == partition and m["model"] == model_name)
                verify_summary(saved_summary, calculated)
                evaluated_panels += len(hosts)
    return {
        "status": "verified", "study_id": plan["study_id"], "report_sha256": sha256_file(report_path),
        "verification_script_sha256": sha256_file(__file__), "protocol_sha256": provenance["protocol_sha256"],
        "pre_fit_metadata_erratum_sha256": sha256_file(erratum_path),
        "source_resources_verified": len(resources), "proteomes_verified": len(admitted),
        "observed_labels_recurated": len(observed), "unknown_labels_excluded": len(excluded),
        "checkpoint_receipts_verified": len(checkpoints), "unique_checkpoint_artifacts_verified": artifacts,
        "models_reloaded": models_verified, "model_predictions_bit_exact": True,
        "numpy_median_max_abs_difference_from_saved_ensemble": max_ensemble_error,
        "retrieval_summaries_reproduced_from_reloaded_models": True,
        "host_method_panels_recalculated": evaluated_panels,
        "feature_sets_identical_across_complete_repeats": identical_feature_repeats,
        "prevalence_checks": prevalence_verified,
        "checkpointed_compute_wall_seconds_including_evaluation": sum(r["wall_seconds"] for r in checkpoints),
        "wall_time_clock": "Sum of time.perf_counter stage durations; excludes between-stage overhead and is independent of wall-clock adjustments",
        "verification_wall_seconds": time.perf_counter() - started,
    }


def publish_verified_results(root, destination, qa, *, runtime_receipt=None):
    """Publish verified derived results and receipts, excluding raw/model payloads."""
    root, destination = Path(root), Path(destination)
    if qa.get("status") != "verified" or sha256_file(root / "evaluation/report.json") != qa["report_sha256"]:
        raise ValueError("publication requires the verified, unchanged run report")
    if destination.exists():
        raise ValueError("published result directory already exists")
    report = read_json_object(root / "evaluation/report.json")
    paths = [(root / "evaluation/report.json", Path("report.json"))]
    for name, expected in report["output_sha256"].items():
        if Path(name).name != name or not name.endswith(".csv"):
            raise ValueError("unexpected evaluation artifact for publication")
        source = root / "evaluation" / name
        if sha256_file(source) != expected:
            raise ValueError("changed verified evaluation artifact")
        paths.append((source, Path(name)))
    for name in ("provenance.json", "protocol.json"):
        paths.append((root / name, Path(name)))
    for path in sorted((root / "checkpoints").glob("*.json")):
        paths.append((path, Path("checkpoints") / path.name))
    for path in sorted((root / "fits").glob("*/feature_occurrences.csv")):
        paths.append((path, Path("feature_occurrences") / f"{path.parent.name}.csv"))
    for name in ("host_manifest.csv", "scenario_splits.csv"):
        paths.append((root / "cohort" / name, Path(name)))
    if runtime_receipt is not None:
        runtime = read_json_object(runtime_receipt)
        if runtime.get("status") != "completed" or runtime.get("elapsed_wall_seconds", 0) <= 0:
            raise ValueError("runtime receipt does not describe a completed run")
        paths.append((Path(runtime_receipt), Path("runtime.json")))
    with staged_directory(destination) as stage:
        for source, relative in paths:
            target = stage / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
        result = {**qa, "published_output_sha256": {
            p.relative_to(stage).as_posix(): sha256_file(p) for p in sorted(stage.rglob("*")) if p.is_file()
        }}
        atomic_write_json(result, stage / "NUMERIC_QA.json")
    return result


def verify_seed_baseline(repo, root, data_root, protocol):
    """Reload original-parameter refits under their exact numerical environment."""
    repo, root = Path(repo).resolve(), Path(root).resolve()
    started = time.perf_counter()
    protocol_path, plan = load_study_protocol(repo, protocol, "genome")
    report = read_json_object(root / "report.json")
    if report["provenance"]["protocol_sha256"] != sha256_file(protocol_path):
        raise ValueError("baseline output belongs to another seed protocol")
    for relative, digest in report["output_sha256"].items():
        if sha256_file(root / relative) != digest:
            raise ValueError("changed baseline output")
    if set(m["scenario"] for m in report["metrics"]) != {"grouped_host"}:
        raise ValueError("baseline sensitivity refitted a seed-invariant scenario")
    resources = verify_resources(plan, data_root)
    observed, _, _ = curate_mlc(load_interactions(resources[("table", "subset")]),
        pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0))
    manifest = pd.read_csv(repo / plan["host_manifest"], dtype=str)
    hosts = sorted(manifest.bacterium)
    annotations, numeric = host_features(hosts,
        pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0),
        pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0), plan["categorical_features"])
    genome, _ = host_genome_features(manifest, resources)
    feature_sets = {
        "host_logistic": (annotations, numeric, plan["categorical_features"], "model_00.json"),
        "host_4mer_logistic": (genome, list(genome), [], "model_00_host_4mer_logistic.json"),
        "hybrid_4mer_logistic": (annotations.join(genome), numeric + list(genome),
                                   plan["categorical_features"], "model_00_hybrid_4mer_logistic.json"),
    }
    train_hosts = set(manifest.loc[manifest.split.eq("train"), "bacterium"])
    targets = sorted(manifest.loc[manifest.split.isin(["validation", "test"]), "bacterium"])
    panel = sorted(observed.phage.unique())
    candidates = pd.MultiIndex.from_product([targets, panel], names=["bacterium", "phage"]).to_frame(index=False)
    saved = pd.read_csv(root / "predictions.csv.gz", dtype={"bacterium": str, "phage": str}, float_precision="round_trip")
    fresh_scores = {}
    for name, (features, columns, categories, filename) in feature_sets.items():
        model = read_json_object(root / filename)
        assert_baseline_model_split(model, manifest)
        encoder = HostFeatureEncoder.fit(features.loc[sorted(train_hosts)], columns, categories)
        if encoder.to_dict() != model["preprocessing"]:
            raise ValueError("baseline preprocessing was not fitted on its exact training hosts")
        with threadpool_limits(limits=1):
            predicted = predict_host_model(model, candidates, features)
        recorded = saved[saved.model.eq(name)].set_index(["bacterium", "phage"]).score.reindex(
            pd.MultiIndex.from_frame(candidates)).to_numpy()
        np.testing.assert_array_equal(predicted, recorded)
        fresh_scores[name] = predicted
    counts = observed[observed.bacterium.isin(train_hosts)].groupby("phage").infection.agg(["sum", "count"])
    rates = (counts["sum"] + 1) / (counts["count"] + 2)
    fresh_scores["phage_prevalence"] = candidates.phage.map(rates).to_numpy(dtype=float)
    np.testing.assert_array_equal(fresh_scores["phage_prevalence"],
        saved[saved.model.eq("phage_prevalence")].set_index(["bacterium", "phage"]).score.reindex(
            pd.MultiIndex.from_frame(candidates)).to_numpy())
    for partition in ("validation", "test"):
        target_hosts = sorted(manifest.loc[manifest.split.eq(partition), "bacterium"])
        mask = candidates.bacterium.isin(target_hosts)
        for name, values in fresh_scores.items():
            summary, _ = evaluate_retrieval(observed[observed.bacterium.isin(target_hosts)],
                candidates.loc[mask].assign(score=values[mask]), hosts=target_hosts, candidate_phages=panel, ks=plan["k"])
            previous = next(m for m in report["metrics"] if m["partition"] == partition and m["model"] == name)
            verify_summary(previous, summary)
    return {
        "status": "verified", "study_id": plan["study_id"], "split_seed": plan["split_seed"],
        "report_sha256": sha256_file(root / "report.json"), "protocol_sha256": sha256_file(protocol_path),
        "verification_script_sha256": sha256_file(__file__), "baseline_models_reloaded": len(feature_sets),
        "baseline_predictions_bit_exact": True, "train_only_prevalence_bit_exact": True,
        "training_hosts": sorted(train_hosts), "target_hosts": targets,
        "numerical_versions": {p: version(p) for p in ("numpy", "pandas", "scikit-learn", "scipy", "joblib", "threadpoolctl")},
        "verification_wall_seconds": time.perf_counter() - started,
    }


def publish_seed_baseline(root, destination, qa):
    root, destination = Path(root), Path(destination)
    if destination.exists() or qa.get("status") != "verified" or sha256_file(root / "report.json") != qa["report_sha256"]:
        raise ValueError("baseline publication requires a new destination and verified report")
    names = ["report.json", "protocol.json", "host_splits.csv", "per_host_metrics.csv", "contrasts.csv",
             "predictions.csv.gz", "genome_audit.csv", "model_00.json",
             "model_00_host_4mer_logistic.json", "model_00_hybrid_4mer_logistic.json"]
    with staged_directory(destination) as stage:
        for name in names:
            (stage / name).write_bytes((root / name).read_bytes())
        result = {**qa, "published_output_sha256": {name: sha256_file(stage / name) for name in names}}
        atomic_write_json(result, stage / "NUMERIC_QA.json")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--proteomes", type=Path, default=Path("data/derived/genophi_proteomes"))
    parser.add_argument("--proteome-receipt", type=Path, default=Path("data/receipts/genophi_proteomes.v1.json"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--publish", type=Path, help="New directory for verified derived results and receipts")
    parser.add_argument("--runtime-receipt", type=Path, help="Optional external process-clock measurement to publish")
    parser.add_argument("--baseline-protocol", type=Path, help="Verify/publish the matching-seed baseline in its original numerical environment")
    args = parser.parse_args()
    if args.out.exists():
        parser.error("numeric QA receipt already exists; choose a new immutable path")
    if args.publish is not None and args.publish.exists():
        parser.error("published result directory already exists")
    if args.baseline_protocol is not None:
        result = verify_seed_baseline(Path(__file__).resolve().parents[1], args.run, args.data_root, args.baseline_protocol)
        if args.publish is not None:
            result = publish_seed_baseline(args.run, args.publish, result)
        atomic_write_json(result, args.out)
        print({k: result[k] for k in ("status", "split_seed", "baseline_models_reloaded", "baseline_predictions_bit_exact")})
        return
    result = verify_run(Path(__file__).resolve().parents[1], args.run, args.data_root, args.proteomes, args.proteome_receipt)
    if args.publish is not None:
        result = publish_verified_results(args.run, args.publish, result, runtime_receipt=args.runtime_receipt)
    atomic_write_json(result, args.out)
    print({k: result[k] for k in ("status", "models_reloaded", "model_predictions_bit_exact",
                                  "retrieval_summaries_reproduced_from_reloaded_models", "verification_wall_seconds")})


if __name__ == "__main__":
    main()
