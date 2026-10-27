"""Locked, table-only known-phage retrieval development study."""
from __future__ import annotations

import gzip
import importlib.metadata
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import sklearn
from threadpoolctl import threadpool_limits

from .artifacts import atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory
from .cohort import curate_mlc, host_aliases, host_key
from .data import load_interactions
from .host_baseline import fit_host_model, predict_host_model
from .host_genome import host_genome_features
from .host_sampling import sampling_features
from .evo2_host import load_cache
from .host_similarity import align_core_distances, neighbour_scores, prevalence_shift
from .host_study import host_features, verify_resources
from .registry import verify_registered_dataset
from .retrieval import evaluate_retrieval, paired_retrieval_bootstrap
from .split import assign_partitions


def retrieval_splits(host_manifest, metadata, *, seed, fractions):
    """Freeze all scenarios without consulting any interaction outcomes.

    Broad phylogroups are source metadata, not newly inferred tree clades. Any
    near-clonal component touching a held-out phylogroup is purged from fitting.
    """
    table = host_manifest.copy()
    required = {"bacterium", "host_group", "split"}
    if not required.issubset(table) or table[list(required)].isna().any().any():
        raise ValueError("host manifest requires complete identities, groups and splits")
    if table.bacterium.duplicated().any() or set(table.split) != {"train", "validation", "test"}:
        raise ValueError("host manifest requires unique hosts and all three partitions")
    if table.groupby("host_group").split.nunique().max() != 1:
        raise ValueError("frozen host groups cross partitions")
    table = table.sort_values("bacterium").reset_index(drop=True)
    lookup = host_aliases(metadata.index)
    table["phylogroup"] = [metadata.loc[lookup[host_key(h)], "Clermont_Phylo"] for h in table.bacterium]
    if table.phylogroup.isna().any():
        raise ValueError("every cohort host needs published phylogroup metadata")
    scenarios = {"grouped_host": table.copy()}
    random = table.copy()
    random["split"] = random.bacterium.map(assign_partitions(table.bacterium, seed=seed, fractions=fractions))
    scenarios["random_host"] = random
    for phylogroup in sorted(table.phylogroup.unique()):
        held = table.phylogroup.eq(phylogroup)
        held_groups = set(table.loc[held, "host_group"])
        split = table.copy()
        split["split"] = np.where(held, "test", np.where(table.host_group.isin(held_groups), "purged", "train"))
        scenarios[f"phylogroup:{phylogroup}"] = split
    return scenarios


def score_known_panel(train, candidates, features, numeric, categorical, estimator):
    """Fit only on train; scoring accepts identities/features and no test labels."""
    if set(train.bacterium) & set(candidates.bacterium):
        raise ValueError("retrieval training and target hosts must be disjoint")
    if set(candidates.phage) - set(train.phage):
        raise ValueError("fixed library includes a phage without training support")
    with threadpool_limits(limits=1):
        model = fit_host_model(train, features, numeric, categorical,
                               C=estimator["C"], max_iter=estimator["max_iter"])
        host_scores = predict_host_model(model, candidates, features)
    stats = train.groupby("phage").infection.agg(["sum", "count"])
    rates = (stats["sum"] + 1) / (stats["count"] + 2)
    return model, {
        # Uniform random order and constant scores have the same exact expectation.
        "random": np.zeros(len(candidates)),
        "global_prevalence": np.full(len(candidates), float(train.infection.mean())),
        "phage_prevalence": candidates.phage.map(rates).to_numpy(dtype=float),
        "host_logistic": host_scores,
    }


def run_retrieval_study(repo, data_root, output_dir, *, fetch=False, genome=False,
                        similarity=False, embedding_cache=None, sampling=False):
    """Verify pinned tables/manifests, freeze splits, refit and publish atomically."""
    repo, output = Path(repo).resolve(), Path(output_dir).resolve()
    if output.exists():
        raise ValueError("retrieval output directory must be new")
    if sampling:
        if genome or similarity or embedding_cache is not None:
            raise ValueError("sampling control must use its own protocol")
        genome = True
    if embedding_cache is not None:
        similarity = True
    if genome and similarity:
        raise ValueError("choose a single study protocol")
    protocol_path = repo / ("data/studies/host_window_sampling.v1.json" if sampling else
                            "data/studies/evo2_host_representation.v1.json" if embedding_cache is not None else
                            "data/studies/host_similarity.v1.json" if similarity else
                            "data/studies/genome_host_retrieval.v1.json" if genome else
                            "data/studies/host_retrieval.v1.json")
    plan = read_json_object(protocol_path)
    if sampling:
        amendment = plan
        parent = repo / amendment["parent_protocol"]
        if sha256_file(parent) != amendment["parent_sha256"]:
            raise ValueError("changed frozen sampling parent protocol")
        plan = read_json_object(parent)
        plan.update(amendment["overrides"])
        plan["sampling_control"] = amendment["sampling_control"]
        plan["repo_files_sha256"][amendment["parent_protocol"]] = amendment["parent_sha256"]
    for path, expected in plan["repo_files_sha256"].items():
        if sha256_file(repo / path) != expected:
            raise ValueError(f"changed frozen retrieval input: {path}")
    resources = verify_resources(plan, data_root, fetch=fetch)
    subset = resources[("table", "subset")]
    verify_registered_dataset(repo / plan["subset_manifest"], subset)
    observed, excluded, _ = curate_mlc(
        load_interactions(subset), pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0)
    )
    manifest = pd.read_csv(repo / plan["host_manifest"], dtype=str)
    if set(manifest.bacterium) != set(observed.bacterium):
        raise ValueError("host manifest and admitted cohort differ")
    metadata = pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0)
    defenses = pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0)
    scenarios = retrieval_splits(manifest, metadata, **plan["random_split"])
    panel = sorted(observed.phage.unique())
    if len(panel) != plan["expected_candidate_phages"]:
        raise ValueError("candidate library differs from the fixed panel")
    features, numeric = host_features(sorted(observed.bacterium.unique()), metadata, defenses,
                                      plan["categorical_features"])
    extra_features, genome_audit = {}, None
    if similarity:
        distances = align_core_distances(
            pd.read_csv(resources[("table", "core_distances")], sep="\t", index_col=0),
            manifest.bacterium)
    if genome:
        sequence_features, genome_audit = host_genome_features(manifest, resources)
        extra_features = {
            "host_4mer_logistic": (sequence_features, sequence_features.columns.tolist(), []),
            "hybrid_4mer_logistic": (features.join(sequence_features),
                                      numeric + sequence_features.columns.tolist(),
                                      plan["categorical_features"]),
        }
    if sampling:
        sampled_features, sampling_audit, sampling_windows = sampling_features(
            manifest, resources, plan["sampling_control"])
        for budget, frame in sampled_features.items():
            extra_features[f"host_window{budget}_4mer_logistic"] = (frame, frame.columns.tolist(), [])
    if embedding_cache is not None:
        embedding_features, sequence_features, genome_audit, cache_sha = load_cache(
            embedding_cache, repo, manifest)
        extra_features = {
            "host_4mer_logistic": (sequence_features, sequence_features.columns.tolist(), []),
            "hybrid_4mer_logistic": (features.join(sequence_features),
                                      numeric + sequence_features.columns.tolist(), plan["categorical_features"]),
            "host_evo2_logistic": (embedding_features, embedding_features.columns.tolist(), []),
            "hybrid_evo2_logistic": (features.join(embedding_features),
                                      numeric + embedding_features.columns.tolist(), plan["categorical_features"]),
        }
    provenance = {
        **code_provenance(repo, repo / ("scripts/run_host_window_sampling.py" if sampling else
                                       "scripts/run_evo2_host_comparison.py" if embedding_cache is not None else
                                       "scripts/run_host_similarity.py" if similarity else
                                       "scripts/run_genome_host_retrieval.py" if genome else
                                       "scripts/run_host_retrieval.py")),
        "protocol_sha256": sha256_file(protocol_path),
        "repo_files_sha256": plan["repo_files_sha256"],
        "resources_sha256": {entry["path"]: entry["sha256"] for entry in plan["resources"]},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scikit_learn": sklearn.__version__},
    }
    if embedding_cache is not None:
        provenance["embedding_cache_sha256"] = cache_sha
    if genome:
        provenance["versions"].update({name: importlib.metadata.version(name)
                                       for name in ("scipy", "joblib", "threadpoolctl")})
    summaries, contrasts, per_host_tables, prediction_tables = [], [], [], []
    distance_tables, shift_tables = [], []
    groups = manifest.set_index("bacterium").host_group.to_dict()
    with staged_directory(output) as stage:
        # Publish every scenario's partition before the first fit.
        pd.concat([frame.assign(scenario=name) for name, frame in scenarios.items()],
                  ignore_index=True).to_csv(stage / "host_splits.csv", index=False)
        excluded.to_csv(stage / "excluded_pairs.csv", index=False)
        observed.to_csv(stage / "observed_pairs.csv", index=False)
        atomic_write_json(plan, stage / "protocol.json")
        if sampling:
            (stage / "sampling_protocol.json").write_bytes(protocol_path.read_bytes())
            sampling_audit.to_csv(stage / "sampling_audit.csv", index=False)
            sampling_windows.to_csv(stage / "sampling_windows.csv", index=False)
            for budget, frame in sampled_features.items():
                frame.to_csv(stage / f"host_window{budget}_features.csv", index_label="bacterium")
        if genome or embedding_cache is not None:
            genome_audit.to_csv(stage / "genome_audit.csv", index=False)
            sequence_features.to_csv(stage / "host_genome_features.csv", index_label="bacterium")
        for number, (scenario, split) in enumerate(scenarios.items()):
            train_hosts = set(split.loc[split.split.eq("train"), "bacterium"])
            train = observed[observed.bacterium.isin(train_hosts)].copy()
            target_hosts = sorted(split.loc[split.split.isin(["validation", "test"]), "bacterium"])
            candidates = pd.MultiIndex.from_product(
                [target_hosts, panel], names=["bacterium", "phage"]
            ).to_frame(index=False)
            model, score_arrays = score_known_panel(
                train, candidates, features, numeric, plan["categorical_features"], plan["estimator"]
            )
            model["provenance"] = provenance
            atomic_write_json(model, stage / f"model_{number:02d}.json")
            if similarity:
                values, distance_diagnostics = neighbour_scores(
                    train, candidates, distances, k=plan["neighbours"]["k"])
                score_arrays["core_neighbour"] = values
                distance_tables.append(distance_diagnostics.merge(
                    split[["bacterium", "split", "phylogroup", "host_group"]],
                    on="bacterium", validate="one_to_one").assign(scenario=scenario))
            for name, (frame, columns, categories) in extra_features.items():
                fitted, values = score_known_panel(
                    train, candidates, frame, columns, categories, plan["estimator"]
                )
                fitted["provenance"] = provenance
                atomic_write_json(fitted, stage / f"model_{number:02d}_{name}.json")
                score_arrays[name] = values["host_logistic"]
            for partition in ("validation", "test"):
                hosts = sorted(split.loc[split.split.eq(partition), "bacterium"])
                if not hosts:
                    continue
                mask = candidates.bacterium.isin(hosts).to_numpy()
                ids = candidates.loc[mask].copy()
                held = observed[observed.bacterium.isin(hosts)]
                if similarity:
                    shift_tables.append(prevalence_shift(train, held, panel).assign(
                        scenario=scenario, partition=partition))
                metrics = {}
                overlap = sorted(set(split.loc[split.bacterium.isin(train_hosts), "host_group"])
                                 & set(split.loc[split.bacterium.isin(hosts), "host_group"]))
                if scenario != "random_host" and overlap:
                    raise ValueError("group leakage in a group-separated scenario")
                for method, values in score_arrays.items():
                    scored = ids.assign(score=values[mask])
                    summary, host_metrics = evaluate_retrieval(
                        held, scored, hosts=hosts, candidate_phages=panel, ks=plan["k"]
                    )
                    context = {"scenario": scenario, "partition": partition, "model": method}
                    summaries.append({**context, **summary, "train_hosts": len(train_hosts),
                                      "train_pairs": len(train), "overlapping_train_host_groups": overlap,
                                      "purged_hosts": int(split.split.eq("purged").sum())})
                    metrics[method] = host_metrics
                    per_host_tables.append(host_metrics.assign(**context))
                    scored["rank_min"] = scored.groupby("bacterium").score.rank(method="min", ascending=False)
                    scored["rank_max"] = scored.groupby("bacterium").score.rank(method="max", ascending=False)
                    prediction_tables.append(scored.assign(**context))
                comparisons = [("host_logistic", "phage_prevalence")]
                if similarity:
                    comparisons.extend([("core_neighbour", "phage_prevalence"),
                                        ("core_neighbour", "host_logistic")])
                for name in extra_features:
                    comparisons.extend([(name, "phage_prevalence"), (name, "host_logistic")])
                if genome or embedding_cache is not None:
                    comparisons.append(("hybrid_4mer_logistic", "host_4mer_logistic"))
                if sampling:
                    budgets = plan["sampling_control"]["budgets"]
                    comparisons.extend((f"host_window{k}_4mer_logistic", "host_4mer_logistic") for k in budgets)
                    comparisons.extend((f"host_window{b}_4mer_logistic", f"host_window{a}_4mer_logistic")
                                       for a, b in zip(budgets, budgets[1:]))
                if embedding_cache is not None:
                    comparisons.extend([
                        ("host_evo2_logistic", "host_4mer_logistic"),
                        ("host_evo2_logistic", "core_neighbour"),
                        ("hybrid_evo2_logistic", "host_evo2_logistic"),
                        ("hybrid_evo2_logistic", "hybrid_4mer_logistic"),
                    ])
                for left, right in comparisons:
                    for eligible_only in (False, True):
                        contrast = paired_retrieval_bootstrap(
                            metrics[left], metrics[right], groups,
                            eligible_only=eligible_only, **plan["bootstrap"]
                        )
                        contrasts.append({"scenario": scenario, "partition": partition,
                                          "comparison": f"{left} - {right}", **contrast})
            print(f"Completed {scenario}: {len(train_hosts)} fit hosts", flush=True)
        pd.DataFrame(contrasts).to_csv(stage / "contrasts.csv", index=False)
        per_host = pd.concat(per_host_tables, ignore_index=True)
        per_host.to_csv(stage / "per_host_metrics.csv", index=False)
        if similarity:
            distance_table = pd.concat(distance_tables, ignore_index=True)
            distance_table.to_csv(stage / "host_distances.csv", index=False)
            pd.concat(shift_tables, ignore_index=True).to_csv(stage / "phage_prevalence_shift.csv", index=False)
            diagnostics = per_host.merge(
                distance_table.rename(columns={"split": "partition"}),
                on=["scenario", "partition", "bacterium"], validate="many_to_one")
            reference = per_host.loc[per_host.model.eq("phage_prevalence"),
                                    ["scenario", "partition", "bacterium", "hit@2"]]
            diagnostics = diagnostics.merge(reference.rename(columns={"hit@2": "prevalence_hit@2"}),
                                            on=["scenario", "partition", "bacterium"], validate="many_to_one")
            diagnostics["hit@2_delta_vs_prevalence"] = diagnostics["hit@2"] - diagnostics["prevalence_hit@2"]
            diagnostics.to_csv(stage / "host_diagnostics.csv", index=False)
        payload = pd.concat(prediction_tables, ignore_index=True).to_csv(index=False, lineterminator="\n").encode()
        (stage / "predictions.csv.gz").write_bytes(gzip.compress(payload, mtime=0))
        report = {
            "study_id": plan["study_id"], "provenance": provenance,
            "audit": {"observed_pairs": len(observed), "unknown_pairs_excluded": len(excluded),
                      "hosts": len(manifest), "host_groups": manifest.host_group.nunique(),
                      "candidate_phages": len(panel)},
            "metrics": summaries, "contrasts": contrasts,
            "not_run": ({k: v for k, v in plan["not_run"].items() if k != "Evo2"}
                        if embedding_cache is not None else plan["not_run"]), "limitations": plan["limitations"],
            "output_sha256": {p.name: sha256_file(p) for p in sorted(stage.iterdir())},
        }
        if genome or embedding_cache is not None:
            report["audit"].update(host_genomes_verified=len(genome_audit),
                                    genome_features=sequence_features.shape[1])
        if sampling:
            report["audit"]["sampling_control_rows"] = len(sampling_audit)
        atomic_write_json(report, stage / "report.json")
    return report
