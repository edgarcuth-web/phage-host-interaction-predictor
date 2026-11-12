"""Offline numerical QA for a completed frozen-split GenoPHI comparator run."""
from __future__ import annotations

import argparse
from importlib.metadata import version
from pathlib import Path
import time

import joblib
import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file, staged_directory
from phage_cocktail.cohort import curate_mlc
from phage_cocktail.data import load_interactions
from phage_cocktail.genophi_comparator import feature_columns, selected_occurrences
from phage_cocktail.genophi_proteomes import verify_proteomes
from phage_cocktail.host_study import verify_resources
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
    original_plan = read_json_object(repo / "data/studies/genophi_comparator.v1.json")
    if plan != original_plan or read_json_object(root / "protocol.json") != original_plan:
        raise ValueError("changed run protocol")
    if sha256_file(repo / "data/studies/genophi_comparator.v1.json") != provenance["protocol_sha256"]:
        raise ValueError("changed frozen comparator protocol")
    erratum_path = repo / "data/studies/genophi_comparator.v1.errata.json"
    if read_json_object(erratum_path)["protocol_sha256"] != provenance["protocol_sha256"]:
        raise ValueError("metadata erratum belongs to a different protocol")
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
        for partition in plan["evaluation"]["partitions"]:
            hosts = sorted(split.loc[split.split.eq(partition), "bacterium"])
            if not hosts:
                continue
            mask = keys.bacterium.isin(hosts).to_numpy()
            held = observed[observed.bacterium.isin(hosts)]
            for model_name, probabilities in [("genophi_median", ensemble),
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--proteomes", type=Path, default=Path("data/derived/genophi_proteomes"))
    parser.add_argument("--proteome-receipt", type=Path, default=Path("data/receipts/genophi_proteomes.v1.json"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--publish", type=Path, help="New directory for verified derived results and receipts")
    parser.add_argument("--runtime-receipt", type=Path, help="Optional external process-clock measurement to publish")
    args = parser.parse_args()
    if args.out.exists():
        parser.error("numeric QA receipt already exists; choose a new immutable path")
    if args.publish is not None and args.publish.exists():
        parser.error("published result directory already exists")
    result = verify_run(Path(__file__).resolve().parents[1], args.run, args.data_root, args.proteomes, args.proteome_receipt)
    if args.publish is not None:
        result = publish_verified_results(args.run, args.publish, result, runtime_receipt=args.runtime_receipt)
    atomic_write_json(result, args.out)
    print({k: result[k] for k in ("status", "models_reloaded", "model_predictions_bit_exact",
                                  "retrieval_summaries_reproduced_from_reloaded_models", "verification_wall_seconds")})


if __name__ == "__main__":
    main()
