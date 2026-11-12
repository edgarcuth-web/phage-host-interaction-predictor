"""Train-only GenoPHI adapter for the frozen host splits and retrieval estimand.

GenoPHI remains optional. Its public Python helpers are invoked separately so
expensive stages can be checkpointed without its automatic cutoff selection.
"""
from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from importlib.metadata import version
from importlib.util import find_spec
import logging
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from unittest.mock import patch

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .artifacts import (
    atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory,
)
from .cohort import canonical_assembly_digest, curate_mlc, group_hosts
from .data import load_interactions
from .genophi_proteomes import verify_proteomes
from .host_study import verify_resources
from .registry import verify_registered_dataset
from .retrieval import evaluate_retrieval, paired_retrieval_bootstrap
from .retrieval_study import retrieval_splits
from .split import split_by_host_group

PROTOCOL = "data/studies/genophi_comparator.v1.json"
IDS = ["strain", "phage"]
TARGET = "interaction"


class ComparatorBlocked(ValueError):
    """A prespecified comparison cannot produce a complete, fair score panel."""


def select_best_hits_streaming(assigned_tsv, best_hits_tsv, clusters_tsv, *, chunksize=250_000):
    """Bound memory without changing upstream select_best_hits semantics.

    Upstream materializes and sorts the entire all-versus-all alignment table.
    Reduce chunks by the same float32 identity/e-value keys, retaining the first
    input row on exact ties, then merge the unchanged target-to-family mapping.
    """
    names = ["Query", "Target", "Score", "SeqIdentity", "E-value",
             "qStartPos", "qEndPos", "qLen", "tStartPos", "tEndPos", "tLen"]
    dtypes = {**dict.fromkeys(names[:2], "string"),
              **dict.fromkeys(names[2:5], "float32"), **dict.fromkeys(names[5:], "int32")}
    best, alignments = {}, 0
    for chunk in pd.read_csv(assigned_tsv, sep="\t", header=None, names=names,
                             dtype=dtypes, chunksize=chunksize):
        alignments += len(chunk)
        if chunk[names[:5]].isna().any().any():
            raise ComparatorBlocked("missing alignment identity, target or ranking key")
        reduced = chunk.sort_values(["Query", "SeqIdentity", "E-value"],
                                    ascending=[True, False, True]).drop_duplicates("Query")
        for query, target, identity, evalue in reduced[["Query", "Target", "SeqIdentity", "E-value"]].itertuples(index=False, name=None):
            previous = best.get(query)
            if previous is None or identity > previous[1] or (identity == previous[1] and evalue < previous[2]):
                best[query] = (target, identity, evalue)
    if not best:
        raise ComparatorBlocked("empty upstream assignment table")
    clusters = pd.read_csv(clusters_tsv, sep="\t", header=None, names=["Cluster", "Contig"], dtype=str)
    if clusters.Contig.duplicated().any():
        raise ComparatorBlocked("protein assigned to multiple original clusters")
    selected = pd.DataFrame([(q, best[q][0]) for q in sorted(best)], columns=["Query", "Target"])
    mapped = selected.merge(clusters, left_on="Target", right_on="Contig", how="left", validate="many_to_one")
    if mapped.Cluster.isna().any():
        raise ComparatorBlocked("assignment target missing from original clusters")
    mapped[["Query", "Cluster"]].to_csv(best_hits_tsv, sep="\t", header=False, index=False)
    logging.info("Reduced %d alignments to %d best hits with bounded memory", alignments, len(mapped))


def feature_columns(frame):
    columns = [c for c in frame if c.startswith(("sc_", "pc_"))]
    if (not columns or frame.columns.duplicated().any()
            or set(frame) - set([*IDS, TARGET, *columns])):
        raise ValueError("GenoPHI inputs must contain only identifiers, target and sc_/pc_ features")
    values = frame[columns].to_numpy(dtype=float)
    if not np.isfinite(values).all() or not np.isin(values, [0, 1]).all():
        raise ValueError("protein-family features must be finite binary values")
    return columns


def assemble_feature_table(phenotype, host_features, phage_features):
    """Join by exact resource IDs, with no silent row loss or outcome metadata."""
    if list(phenotype.columns) != [*IDS, TARGET] or phenotype.duplicated(IDS).any():
        raise ValueError("phenotype must have unique strain,phage,interaction rows")
    result = phenotype.copy()
    for features, identity, prefix in [(host_features, "strain", "sc_"),
                                        (phage_features, "phage", "pc_")]:
        if (identity not in features or features[identity].duplicated().any()
                or set(features[identity]) != set(phenotype[identity])
                or any(not c.startswith(prefix) for c in features if c != identity)):
            raise ValueError("feature-table IDs/columns differ from the frozen genome sets")
        result = result.merge(features, on=identity, how="left", validate="many_to_one")
    feature_columns(result)
    if result[TARGET].isna().any() or not result[TARGET].isin([0, 1]).all():
        raise ValueError("training labels must be observed binary values")
    return result


def training_table(full, split):
    """Admit only outer-training rows before any supervised GenoPHI function."""
    feature_columns(full)
    if split.bacterium.duplicated().any() or split.groupby("host_group").split.nunique().max() != 1:
        raise ValueError("outer host identities/groups cross partitions")
    if set(full.strain) != set(split.bacterium):
        raise ValueError("feature table and outer host manifest differ")
    hosts = set(split.loc[split.split.eq("train"), "bacterium"])
    held = set(split.loc[split.split.isin(["validation", "test", "purged"]), "bacterium"])
    train = full[full.strain.isin(hosts)].copy().reset_index(drop=True)
    if not hosts or set(train.strain) != hosts or set(train.strain) & held:
        raise ValueError("empty training support or held-out host leakage")
    if set(train.phage) != set(full.phage):
        raise ComparatorBlocked("training rows lack a fixed-panel phage")
    return train


def candidate_features(host_features, phage_features, hosts, panel):
    """Construct the complete label-free product, including unknown-label pairs."""
    ids = pd.MultiIndex.from_product([sorted(hosts), sorted(panel)], names=IDS).to_frame(index=False)
    if host_features.strain.duplicated().any() or phage_features.phage.duplicated().any():
        raise ValueError("genome feature rows must be unique")
    result = ids.merge(host_features, on="strain", how="left", validate="many_to_one").merge(
        phage_features, on="phage", how="left", validate="many_to_one")
    if TARGET in result:
        raise ValueError("prediction inputs must not contain held-out labels")
    feature_columns(result)
    return result


def predict_fitted_model(model, candidates, *, threads):
    """Use GenoPHI's positive-class confidence, rejecting partial model support."""
    available = feature_columns(candidates)
    required = list(model.feature_names_)
    if TARGET in candidates or not required or set(required) - set(available):
        raise ValueError("model requires missing or non-protein features")
    if list(model.classes_) != [0, 1]:
        raise ComparatorBlocked("model does not support both binary classes")
    values = np.asarray(model.predict_proba(candidates[required], thread_count=threads), dtype=float)
    if (values.shape != (len(candidates), 2) or not np.isfinite(values).all()
            or (values < 0).any() or (values > 1).any()):
        raise ComparatorBlocked("invalid or incomplete GenoPHI predictions")
    return candidates[IDS].rename(columns={"strain": "bacterium"}).assign(score=values[:, 1])


def selected_occurrences(directory, *, num_runs, allowed_features):
    """Count each feature at most once per completed FS run, including resumes."""
    counts = {}
    for i in range(num_runs):
        path = Path(directory) / f"run_{i}" / "feature_importances.csv"
        if not path.is_file():
            raise ComparatorBlocked(f"missing feature-selection run: {i}")
        frame = pd.read_csv(path)
        if ("Feature" not in frame or frame.Feature.duplicated().any()
                or set(frame.Feature) - set(allowed_features)):
            raise ComparatorBlocked("invalid or non-protein selected feature")
        for feature in frame.Feature:
            counts[feature] = counts.get(feature, 0) + 1
    return pd.DataFrame(sorted(counts.items()), columns=["Feature", "Occurrence"]).sort_values(
        ["Occurrence", "Feature"], ascending=[False, True]).reset_index(drop=True)


@contextmanager
def _capture_output(stream):
    """Capture Python output and native MMseqs/CatBoost child-process output."""
    sys.stdout.flush()
    sys.stderr.flush()
    stdout, stderr = os.dup(1), os.dup(2)
    try:
        os.dup2(stream.fileno(), 1)
        os.dup2(stream.fileno(), 2)
        with redirect_stdout(stream), redirect_stderr(stream):
            yield
            stream.flush()
    finally:
        os.dup2(stdout, 1)
        os.dup2(stderr, 2)
        os.close(stdout)
        os.close(stderr)


def checkpoint(root, name, action):
    """Publish a stage receipt only on success; hash-check every reused artifact."""
    root = Path(root)
    receipts = root / "checkpoints"
    receipts.mkdir(parents=True, exist_ok=True)
    path = receipts / f"{name}.json"
    if path.exists():
        record = read_json_object(path)
        for relative, digest in record["output_sha256"].items():
            payload = (root / relative).resolve()
            if not payload.is_relative_to(root.resolve()) or not payload.is_file() or sha256_file(payload) != digest:
                raise ValueError(f"changed checkpoint artifact: {relative}")
        print(f"Verified checkpoint {name}", flush=True)
        return record
    started = time.perf_counter()
    utc = datetime.now(timezone.utc).isoformat()
    logs = root / "logs"
    logs.mkdir(exist_ok=True)
    existing = list(logs.glob(f"{name}.*.log"))
    log_path = logs / f"{name}.{len(existing):02d}.log"
    print(f"Starting {name}; log: {log_path}", flush=True)
    with log_path.open("x", encoding="utf-8", buffering=1) as stream:
        handler = logging.StreamHandler(stream)
        logging.getLogger().addHandler(handler)
        try:
            with _capture_output(stream):
                outputs, result = action()
        except Exception as exc:
            atomic_write_json({
                "stage": name, "status": "blocked" if isinstance(exc, ComparatorBlocked) else "failed",
                "error_type": type(exc).__name__, "error": str(exc), "started_utc": utc,
                "wall_seconds": time.perf_counter() - started,
            }, logs / f"{name}.{len(existing):02d}.failure.json")
            raise
        finally:
            logging.getLogger().removeHandler(handler)
    record = {
        "stage": name, "status": "complete", "started_utc": utc,
        "wall_seconds": time.perf_counter() - started, "result": result,
        "output_sha256": {Path(p).relative_to(root).as_posix(): sha256_file(p) for p in outputs},
    }
    atomic_write_json(record, path)
    print(f"Completed {name}: {record['wall_seconds']:.1f} s", flush=True)
    return record


def runtime_identity():
    executable = shutil.which("mmseqs")
    if executable is None:
        raise ValueError("MMseqs2 is required on PATH")
    spec = find_spec("genophi")
    if spec is None:
        raise ValueError("genophi==1.0.3 is required")
    package = Path(spec.origin).parent
    return {
        "versions": {"python": platform.python_version(),
                     **{p: version(p) for p in ("genophi", "pyrodigal", "numpy", "pandas", "scikit-learn",
                                                "scipy", "catboost", "biopython", "joblib", "threadpoolctl",
                                                "matplotlib", "seaborn", "plotnine", "shap", "hdbscan")},
                     "mmseqs2": subprocess.check_output([executable, "version"], text=True).strip()},
        "genophi_source_sha256": {p.relative_to(package).as_posix(): sha256_file(p)
                                  for p in sorted(package.rglob("*.py"))},
        "mmseqs_executable_sha256": sha256_file(executable),
        "platform": platform.platform(), "machine": platform.machine(),
    }


def _freeze_cohort(repo, data_root, root, plan, lock):
    resources = verify_resources(lock, data_root)
    parent = read_json_object(repo / plan["parent_protocol"])
    verify_registered_dataset(repo / parent["subset_manifest"], resources[("table", "subset")])
    observed, excluded, aliases = curate_mlc(
        load_interactions(resources[("table", "subset")]),
        pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0))
    hosts, panel = sorted(observed.bacterium.unique()), sorted(observed.phage.unique())
    digests = {h: canonical_assembly_digest(resources[("host_genome", h)]) for h in hosts}
    distances = pd.read_csv(resources[("table", "core_distances")], sep="\t", index_col=0)
    groups = group_hosts(hosts, distances, digests, threshold=parent["grouping"]["core_distance_threshold"])
    assigned = split_by_host_group(observed, groups, seed=parent["split"]["seed"], fractions=parent["split"]["fractions"])
    manifest = assigned[["bacterium", "host_group", "split"]].drop_duplicates().sort_values("bacterium").reset_index(drop=True)
    manifest["source_bacterium"] = manifest.bacterium.map(aliases)
    manifest["assembly_sha256"] = manifest.bacterium.map(digests)
    frozen = pd.read_csv(repo / plan["host_manifest"], dtype=str).sort_values("bacterium").reset_index(drop=True)
    if not manifest.equals(frozen):
        raise ComparatorBlocked("recomputed host split differs from the frozen study")
    audit = {"hosts": len(hosts), "phages": len(panel), "observed_pairs": len(observed),
             "excluded_pairs": len(excluded), "host_groups": len(set(groups.values()))}
    if audit != plan["expected"]:
        raise ComparatorBlocked("curated cohort differs from prespecified support")
    scenarios = retrieval_splits(manifest,
        pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0),
        seed=parent["split"]["seed"], fractions=parent["split"]["fractions"])
    with staged_directory(root / "cohort") as stage:
        observed.to_csv(stage / "observed_pairs.csv", index=False)
        excluded.to_csv(stage / "excluded_pairs.csv", index=False)
        manifest.to_csv(stage / "host_manifest.csv", index=False)
        phenotype = observed.rename(columns={"bacterium": "strain", "infection": TARGET})[[*IDS, TARGET]]
        phenotype.to_csv(stage / "phenotype.csv", index=False)
        pd.concat([scenarios[s].assign(scenario=s) for s in plan["scenarios"]], ignore_index=True).to_csv(stage / "scenario_splits.csv", index=False)
    return list((root / "cohort").iterdir()), audit


def _cluster_set(root, proteomes, genome_type, plan):
    from genophi import mmseqs2_clustering as upstream

    directory = "strain_AAs" if genome_type == "strain" else "phage_AAs"
    output, tmp = root / "clustering" / genome_type, root / "tmp" / genome_type
    # Only the memory reduction is adapted; no search, clustering or feature
    # parameter is changed, and equivalence is checked against the upstream reducer.
    with patch.object(upstream, "select_best_hits", select_best_hits_streaming):
        upstream.run_clustering_workflow(
            str(proteomes / directory), str(output), str(tmp),
            **{k: plan["clustering"][k] for k in ("min_seq_id", "coverage", "sensitivity", "suffix")},
            threads=plan["compute"]["threads"], strain_column=genome_type)
    matrix = output / "presence_absence_matrix.csv"
    # Check the best-hit reduction independently on a deterministic alignment
    # prefix, including at least one chunk boundary, before feature use.
    assigned = output / "assigned_clusters.tsv"
    if assigned.is_file():
        sample = output / "reducer_sample.tsv"
        with assigned.open("rb") as source, sample.open("wb") as target:
            for i, line in enumerate(source):
                if i == 50_000:
                    break
                target.write(line)
        upstream.select_best_hits(str(sample), str(output / "reducer_upstream.tsv"), str(output / "clusters.tsv"))
        select_best_hits_streaming(sample, output / "reducer_streaming.tsv", output / "clusters.tsv", chunksize=997)
        if (output / "reducer_upstream.tsv").read_bytes() != (output / "reducer_streaming.tsv").read_bytes():
            raise ComparatorBlocked("streaming assignment differs from upstream best-hit semantics")
    shape = pd.read_csv(matrix).shape
    files = [matrix, output / "best_hits.tsv", output / "clusters.tsv"]
    files.extend(output.glob("reducer_*.tsv"))
    return files, {"genomes": shape[0], "protein_families": shape[1] - 1,
                   "best_hit_reducer": "streaming; deterministic prefix byte-matches upstream"}


def _assign_set(root, genome_type, plan):
    from genophi.mmseqs2_clustering import run_feature_assignment

    output = root / "clustering" / genome_type
    features = output / "features"
    run_feature_assignment(str(output / "presence_absence_matrix.csv"), str(features),
                           source=genome_type, select_column=genome_type,
                           max_ram=plan["compute"]["max_ram_gb"], threads=plan["compute"]["threads"])
    table = pd.read_csv(features / "feature_table.csv")
    return [p for p in features.iterdir() if p.is_file()], {"genomes": len(table), "features": table.shape[1] - 1}


def _merge_features(root):
    phenotype = pd.read_csv(root / "cohort/phenotype.csv", dtype={"strain": str, "phage": str})
    hosts = pd.read_csv(root / "clustering/strain/features/feature_table.csv", dtype={"strain": str})
    phages = pd.read_csv(root / "clustering/phage/features/feature_table.csv", dtype={"phage": str})
    full = assemble_feature_table(phenotype, hosts, phages)
    output = root / "full_feature_table.csv"
    full.to_csv(output, index=False)
    return [output], {"observed_rows": len(full), "features": len(feature_columns(full))}


def _fit_repeat(root, scenario, repeat, plan):
    import joblib
    from genophi.feature_selection import run_feature_selection_iterations, generate_feature_tables
    from genophi.select_feature_modeling import model_testing_select_MCC

    split = pd.read_csv(root / "cohort/scenario_splits.csv", dtype=str)
    split = split[split.scenario.eq(scenario)]
    full = pd.read_csv(root / "full_feature_table.csv", dtype={"strain": str, "phage": str})
    train = training_table(full, split)
    held_hosts = sorted(split.loc[split.split.isin(["validation", "test"]), "bacterium"])
    hosts = pd.read_csv(root / "clustering/strain/features/feature_table.csv", dtype={"strain": str})
    phages = pd.read_csv(root / "clustering/phage/features/feature_table.csv", dtype={"phage": str})
    candidates = candidate_features(hosts, phages, held_hosts, sorted(full.phage.unique()))
    tag = f"{scenario.replace(':', '_')}_repeat_{repeat}"
    fit = root / "fits" / tag
    fit.mkdir(parents=True, exist_ok=True)
    train_path = fit / "train_full_feature_table.csv"
    if train_path.exists():
        if not pd.read_csv(train_path, dtype={"strain": str, "phage": str}).equals(train):
            raise ValueError("changed outer-training feature table")
    else:
        train.to_csv(train_path, index=False)
    fs, compute = plan["feature_selection"], plan["compute"]
    selection = fit / "feature_selection"
    common = dict(task_type="classification", sample_column="strain", phenotype_column=TARGET,
                  phage_column="phage", use_dynamic_weights=False, use_clustering=True,
                  cluster_method="hierarchical", n_clusters=fs["n_clusters"], max_ram=compute["max_ram_gb"])
    for run in range(fs["num_runs_fs"]):
        def select(run=run):
            with threadpool_limits(limits=compute["threads"]):
                run_feature_selection_iterations(str(train_path), str(selection), threads=compute["threads"],
                    num_features=fs["num_features"], filter_type="strain", num_runs=run + 1,
                    method=fs["method"], check_feature_presence=False, filter_by_cluster_presence=False, **common)
            directory = selection / f"run_{run}"
            if not (directory / "feature_importances.csv").is_file():
                raise ComparatorBlocked("upstream feature selection produced no supported model")
            return list(directory.glob("*.csv")) + [train_path], {"outer_train_hosts": train.strain.nunique(), "outer_train_pairs": len(train)}
        checkpoint(root, f"{tag}_fs_{run}", select)
    tables = fit / "selected"
    selected = tables / f"select_feature_table_cutoff_{fs['cutoff']}.csv"

    def finalize_selection():
        occurrences = selected_occurrences(selection, num_runs=fs["num_runs_fs"], allowed_features=feature_columns(train))
        chosen = occurrences.loc[occurrences.Occurrence.ge(fs["cutoff"]), "Feature"].tolist()
        if not fs["min_features"] <= len(chosen) < fs["max_features"]:
            raise ComparatorBlocked(f"fixed cutoff {fs['cutoff']} has {len(chosen)} features outside the prespecified window")
        generate_feature_tables(str(selection), str(train_path), str(tables), phenotype_column=TARGET,
                                sample_column="strain", cut_offs=[fs["cutoff"]], binary_data=False,
                                max_features=fs["max_features"], min_features=fs["min_features"], filter_type="strain")
        selected_frame = pd.read_csv(selected, dtype={"strain": str, "phage": str})
        if set(feature_columns(selected_frame)) != set(chosen):
            raise ComparatorBlocked("upstream feature table differs from the fixed occurrence rule")
        # The upstream pivot orders identifier columns via set(). Canonicalize
        # once so a fresh interpreter cannot change a checkpoint's row bytes.
        selected_frame[[*IDS, TARGET, *sorted(chosen)]].sort_values(IDS).to_csv(selected, index=False)
        occurrences.to_csv(fit / "feature_occurrences.csv", index=False)
        return [selected, fit / "feature_occurrences.csv"], {"features": len(chosen), "cutoff": fs["cutoff"]}

    selection_record = checkpoint(root, f"{tag}_selection", finalize_selection)
    selected_count = selection_record["result"]["features"]
    for run in range(plan["modeling"]["num_runs_modeling"]):
        def model_run(run=run):
            directory = fit / "models" / f"run_{run}"
            directory.mkdir(parents=True, exist_ok=True)
            with threadpool_limits(limits=compute["threads"]):
                model_testing_select_MCC(str(selected), str(directory), compute["threads"], run,
                    set_filter="strain", binary_data=False, use_shap=False, **common)
            model_path = directory / "best_model.pkl"
            if not model_path.is_file():
                raise ComparatorBlocked("upstream modeling produced no model; no run may be dropped")
            internal = pd.read_csv(directory / "best_model_predictions.csv", dtype={"strain": str})
            if not set(internal.strain).issubset(set(train.strain)):
                raise ComparatorBlocked("outer held-out host entered internal model selection")
            scores = predict_fitted_model(joblib.load(model_path), candidates, threads=compute["threads"])
            score_path = directory / "candidate_scores.csv"
            scores.to_csv(score_path, index=False)
            return [model_path, score_path, selected, train_path, *directory.glob("*.csv")], {
                "outer_train_hosts": train.strain.nunique(), "outer_train_pairs": len(train),
                "prediction_hosts": len(held_hosts), "predictions": len(scores), "features": selected_count}
        checkpoint(root, f"{tag}_model_{run}", model_run)
    return fit


def evaluate_score_runs(observed, splits, score_runs, plan):
    """Evaluate complete runs and the prespecified median with the project metrics."""
    summaries, per_host, run_rows, contrasts = [], [], [], []
    panel = sorted(observed.phage.unique())
    for scenario in plan["scenarios"]:
        split = splits[splits.scenario.eq(scenario)]
        training = observed[observed.bacterium.isin(split.loc[split.split.eq("train"), "bacterium"])]
        stats = training.groupby("phage").infection.agg(["sum", "count"])
        rates = (stats["sum"] + 1) / (stats["count"] + 2)
        groups = split.set_index("bacterium").host_group.to_dict()
        for partition in plan["evaluation"]["partitions"]:
            hosts = sorted(split.loc[split.split.eq(partition), "bacterium"])
            if not hosts:
                continue
            held = observed[observed.bacterium.isin(hosts)]
            context = {"scenario": scenario, "partition": partition}
            individual = []
            for repeat in range(plan["repeats"]):
                repeat_scores = []
                for run in range(plan["modeling"]["num_runs_modeling"]):
                    scores = score_runs[(scenario, repeat, run)]
                    scores = scores[scores.bacterium.isin(hosts)].copy()
                    summary, _ = evaluate_retrieval(held, scores, hosts=hosts, candidate_phages=panel, ks=plan["evaluation"]["k"])
                    run_rows.append({**context, "repeat": repeat, "run": run,
                                     **{f"all_{k}": v for k, v in summary["all_hosts"].items()},
                                     **{f"eligible_{k}": v for k, v in summary["eligible_hosts_only"].items()}})
                    repeat_scores.append(scores)
                combined = pd.concat(repeat_scores, ignore_index=True)
                individual.append(combined)
                median = combined.groupby(["bacterium", "phage"], as_index=False).score.median()
                summary, metrics = evaluate_retrieval(held, median, hosts=hosts, candidate_phages=panel, ks=plan["evaluation"]["k"])
                model = f"genophi_repeat_{repeat}"
                summaries.append({**context, "model": model, **summary})
                per_host.append(metrics.assign(**context, model=model))
            ensemble = pd.concat(individual).groupby(["bacterium", "phage"], as_index=False).score.median()
            prevalence = ensemble[["bacterium", "phage"]].assign(score=lambda f: f.phage.map(rates))
            host_tables = {}
            for model, scores in [("genophi_median", ensemble), ("phage_prevalence", prevalence)]:
                summary, metrics = evaluate_retrieval(held, scores, hosts=hosts, candidate_phages=panel, ks=plan["evaluation"]["k"])
                summaries.append({**context, "model": model, **summary})
                per_host.append(metrics.assign(**context, model=model))
                host_tables[model] = metrics
            for eligible in (False, True):
                contrasts.append({**context, "comparison": "genophi_median - phage_prevalence",
                    **paired_retrieval_bootstrap(host_tables["genophi_median"], host_tables["phage_prevalence"], groups,
                        eligible_only=eligible, **plan["evaluation"]["bootstrap"])})
    return summaries, pd.concat(per_host, ignore_index=True), pd.DataFrame(run_rows), contrasts


def _evaluate(repo, root, plan):
    if "split_seed" in plan:
        from .genophi_seed_protocols import verify_seed_baseline_reference
        verify_seed_baseline_reference(repo, plan)
    observed = pd.read_csv(root / "cohort/observed_pairs.csv", dtype={"bacterium": str, "phage": str})
    splits = pd.read_csv(root / "cohort/scenario_splits.csv", dtype=str)
    scores = {}
    for scenario in plan["scenarios"]:
        for repeat in range(plan["repeats"]):
            tag = f"{scenario.replace(':', '_')}_repeat_{repeat}"
            for run in range(plan["modeling"]["num_runs_modeling"]):
                scores[(scenario, repeat, run)] = pd.read_csv(
                    root / "fits" / tag / "models" / f"run_{run}" / "candidate_scores.csv",
                    dtype={"bacterium": str, "phage": str})
    summaries, per_host, run_metrics, contrasts = evaluate_score_runs(observed, splits, scores, plan)
    source_report = read_json_object(repo / plan["reference_report"])
    source_metrics = repo / plan["reference_per_host"]
    if sha256_file(source_metrics) != source_report["output_sha256"]["per_host_metrics.csv"]:
        raise ValueError("changed frozen in-house per-host result")
    reference = pd.read_csv(source_metrics)
    for scenario in plan["scenarios"]:
        for partition in plan["evaluation"]["partitions"]:
            ours = per_host[(per_host.scenario.eq(scenario)) & (per_host.partition.eq(partition)) & per_host.model.eq("genophi_median")]
            if ours.empty:
                continue
            previous = reference[reference.scenario.eq(scenario) & reference.partition.eq(partition)]
            for model in ("phage_prevalence", "host_4mer_logistic", "hybrid_4mer_logistic"):
                values = previous[previous.model.eq(model)]
                if set(values.bacterium) != set(ours.bacterium):
                    raise ComparatorBlocked("reference result has different complete-panel host support")
                if model == "phage_prevalence":
                    current = per_host[per_host.scenario.eq(scenario) & per_host.partition.eq(partition) & per_host.model.eq(model)]
                    a, b = current.set_index("bacterium").sort_index(), values.set_index("bacterium").sort_index()
                    np.testing.assert_allclose(a["hit@2"], b["hit@2"], rtol=0, atol=1e-14)
                else:
                    groups = splits[splits.scenario.eq(scenario)].set_index("bacterium").host_group.to_dict()
                    for eligible in (False, True):
                        contrasts.append({"scenario": scenario, "partition": partition,
                            "comparison": f"genophi_median - frozen_{model}",
                            **paired_retrieval_bootstrap(ours, values, groups, eligible_only=eligible,
                                **plan["evaluation"]["bootstrap"])})
    checkpoints = [read_json_object(p) for p in sorted((root / "checkpoints").glob("*.json"))]
    variability = []
    for (scenario, partition), frame in run_metrics.groupby(["scenario", "partition"]):
        medians = [m["all_hosts"]["hit@2"] for m in summaries if m["scenario"] == scenario
                   and m["partition"] == partition and m["model"].startswith("genophi_repeat_")]
        variability.append({
            "scenario": scenario, "partition": partition, "individual_models": len(frame),
            "individual_model_hit@2_mean": float(frame["all_hit@2"].mean()),
            "individual_model_hit@2_std": float(frame["all_hit@2"].std(ddof=1)),
            "individual_model_hit@2_range": [float(frame["all_hit@2"].min()), float(frame["all_hit@2"].max())],
            "complete_repeat_median_hit@2": medians,
            "complete_repeat_hit@2_std": float(np.std(medians, ddof=1)),
            "scope": plan["variance_scope"],
        })
    report = {
        "study_id": plan["study_id"], "status": "completed_aligned_development_comparator",
        "provenance": read_json_object(root / "provenance.json"), "protocol": plan,
        "metrics": summaries, "contrasts": contrasts, "variability": variability,
        "frozen_in_house_metrics": [m for m in source_report["metrics"] if m["scenario"] in plan["scenarios"]
                                    and m["model"] in {"host_4mer_logistic", "hybrid_4mer_logistic", "phage_prevalence"}],
        "stage_timings": [{k: r[k] for k in ("stage", "wall_seconds", "result")} for r in checkpoints],
        "compute_wall_seconds": sum(r["wall_seconds"] for r in checkpoints),
        "limitations": plan["limitations"],
    }
    with staged_directory(root / "evaluation") as stage:
        per_host.to_csv(stage / "per_host_metrics.csv", index=False)
        run_metrics.to_csv(stage / "individual_model_metrics.csv", index=False)
        pd.DataFrame(contrasts).to_csv(stage / "contrasts.csv", index=False)
        for scenario in plan["scenarios"]:
            combined = pd.concat([s.assign(repeat=r, run=i) for (sc, r, i), s in scores.items() if sc == scenario])
            median = combined.groupby(["bacterium", "phage"], as_index=False).score.median()
            median.to_csv(stage / f"{scenario.replace(':', '_')}_candidate_scores.csv", index=False)
            if "split_seed" in plan:
                split = splits[splits.scenario.eq(scenario)]
                train = observed[observed.bacterium.isin(split.loc[split.split.eq("train"), "bacterium"])]
                counts = train.groupby("phage").infection.agg(["sum", "count"])
                rates = (counts["sum"] + 1) / (counts["count"] + 2)
                median[["bacterium", "phage"]].assign(score=lambda f: f.phage.map(rates)).to_csv(
                    stage / f"{scenario.replace(':', '_')}_prevalence_scores.csv", index=False)
        report["output_sha256"] = {p.name: sha256_file(p) for p in sorted(stage.iterdir())}
        atomic_write_json(report, stage / "report.json")
    return list((root / "evaluation").iterdir()), {"metrics": len(summaries), "individual_models": len(run_metrics)}


def run_comparator(repo, data_root, proteomes, proteome_receipt, output_dir, *, stage="all", scenario=None, repeat=None,
                   protocol=None, feature_cache=None):
    """Execute or resume immutable stages; completed stages are always reverified."""
    repo, root, proteomes = Path(repo).resolve(), Path(output_dir).resolve(), Path(proteomes).resolve()
    from .genophi_seed_protocols import load_study_protocol, import_feature_cache, verify_seed_baseline_reference, CACHE_LOCK

    protocol_path, plan = load_study_protocol(repo, protocol or PROTOCOL, "comparator")
    sensitivity = "split_seed" in plan
    if sensitivity and feature_cache is None:
        raise ValueError("registered sensitivity runs require the verified v1 feature cache")
    if feature_cache is not None and not sensitivity:
        raise ValueError("feature import requires a registered sensitivity protocol")
    for relative, digest in plan["repo_files_sha256"].items():
        if sha256_file(repo / relative) != digest:
            raise ValueError(f"changed frozen input (preserve LF bytes): {relative}")
    lock = read_json_object(repo / plan["resource_lock"])
    receipt = read_json_object(proteome_receipt)
    if receipt["resource_lock_sha256"] != sha256_file(repo / plan["resource_lock"]):
        raise ValueError("proteome receipt uses a different resource lock")
    verify_proteomes(receipt, proteomes, lock)
    runtime = runtime_identity()
    if any(runtime["versions"][k] != v for k, v in plan["versions"].items()):
        raise ValueError("external tool version differs from the frozen comparator protocol")
    reference_binding = None
    if sensitivity:
        cache_lock = read_json_object(repo / CACHE_LOCK)
        if runtime["versions"] != cache_lock["versions"]:
            raise ValueError("sensitivity runtime differs from the frozen GenoPHI numerical environment")
        reference_binding = verify_seed_baseline_reference(repo, plan)
    provenance = {
        **runtime, **code_provenance(repo, repo / "scripts/run_genophi_comparator.py"),
        "protocol_sha256": sha256_file(protocol_path), "proteome_receipt_sha256": sha256_file(proteome_receipt),
        "reference_report_sha256": sha256_file(repo / plan["reference_report"]),
    }
    if sensitivity:
        provenance.update(protocol_path=protocol_path.relative_to(repo).as_posix(),
                          matching_seed_baseline=reference_binding,
                          feature_artifact_lock_sha256=sha256_file(repo / CACHE_LOCK))
    if root.exists():
        if not (root / "provenance.json").is_file() or read_json_object(root / "provenance.json") != provenance:
            raise ValueError("output exists with a different or missing run binding")
    else:
        root.mkdir(parents=True)
        atomic_write_json(provenance, root / "provenance.json")
        atomic_write_json(plan, root / "protocol.json")
    logging.getLogger().setLevel(logging.INFO)
    checkpoint(root, "cohort", lambda: _freeze_cohort(repo, data_root, root, plan, lock))
    if stage == "prepare":
        return root
    if sensitivity:
        checkpoint(root, "reuse_features", lambda: import_feature_cache(
            repo, feature_cache, root, sha256_file(proteome_receipt)))
    else:
        for genome_type in ("strain", "phage"):
            checkpoint(root, f"cluster_{genome_type}", lambda g=genome_type: _cluster_set(root, proteomes, g, plan))
            checkpoint(root, f"features_{genome_type}", lambda g=genome_type: _assign_set(root, g, plan))
        checkpoint(root, "merge_features", lambda: _merge_features(root))
    if stage == "cluster":
        return root
    chosen_scenarios = [scenario] if scenario is not None else plan["scenarios"]
    chosen_repeats = [repeat] if repeat is not None else range(plan["repeats"])
    if set(chosen_scenarios) - set(plan["scenarios"]) or any(r not in range(plan["repeats"]) for r in chosen_repeats):
        raise ValueError("scenario or repeat outside the prespecified protocol")
    for current in chosen_scenarios:
        for r in chosen_repeats:
            _fit_repeat(root, current, r, plan)
    if stage == "fit" or scenario is not None or repeat is not None:
        return root
    checkpoint(root, "evaluation", lambda: _evaluate(repo, root, plan))
    return root
