"""Independently check real similarity scores, model reloads and historical metrics."""
import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.cohort import host_aliases, host_key
from phage_cocktail.host_baseline import predict_host_model
from phage_cocktail.host_study import host_features, verify_resources


def verify(repo, source, run, repeat=None):
    report = read_json_object(run / "report.json")
    for name, digest in report["output_sha256"].items():
        if sha256_file(run / name) != digest:
            raise ValueError(f"changed output: {name}")
    plan = read_json_object(run / "protocol.json")
    resources = verify_resources(plan, source)
    observed = pd.read_csv(run / "observed_pairs.csv")
    splits = pd.read_csv(run / "host_splits.csv")
    predictions = pd.read_csv(run / "predictions.csv.gz", float_precision="round_trip")
    metrics = pd.read_csv(run / "per_host_metrics.csv", float_precision="round_trip")
    historical = pd.read_csv(repo / "data/studies/results/host_retrieval.v1/per_host_metrics.csv",
                             float_precision="round_trip")
    retained = metrics[metrics.model.isin(historical.model.unique())].reset_index(drop=True)
    pd.testing.assert_frame_equal(retained, historical, check_exact=True)
    metadata = pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0)
    defenses = pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0)
    features, _ = host_features(sorted(observed.bacterium.unique()), metadata, defenses,
                                plan["categorical_features"])
    source_distance = pd.read_csv(resources[("table", "core_distances")], sep="\t", index_col=0)
    rows, cols = host_aliases(source_distance.index), host_aliases(source_distance.columns)
    checked_scores, reloads = 0, 0
    for number, scenario in enumerate(splits.scenario.drop_duplicates()):
        target = predictions[predictions.scenario.eq(scenario)]
        annotation = target[target.model.eq("host_logistic")]
        model = read_json_object(run / f"model_{number:02d}.json")
        np.testing.assert_array_equal(predict_host_model(model, annotation, features), annotation.score)
        reloads += 1
        fit_hosts = sorted(splits.loc[splits.scenario.eq(scenario) & splits.split.eq("train"), "bacterium"])
        train = observed[observed.bacterium.isin(fit_hosts)]
        # Scalar implementation independent of the vectorized transfer scorer.
        for host, scores in target[target.model.eq("core_neighbour")].groupby("bacterium"):
            d = [float(source_distance.loc[rows[host_key(host)], cols[host_key(h)]]) for h in fit_hosts]
            boundary = sorted(d)[min(plan["neighbours"]["k"], len(d)) - 1]
            remaining = min(plan["neighbours"]["k"], len(d)) - sum(v < boundary for v in d)
            tied = sum(v == boundary for v in d)
            weights = {h: (1. if v < boundary else remaining/tied if v == boundary else 0.)
                       for h, v in zip(fit_hosts, d)}
            for r in scores.itertuples(index=False):
                labels = train[train.phage.eq(r.phage)]
                total = sum(weights[x.bacterium] for x in labels.itertuples(index=False))
                positive = sum(weights[x.bacterium] * x.infection for x in labels.itertuples(index=False))
                expected = ((positive + 1) / (total + 2) if total else
                            (labels.infection.sum() + 1) / (len(labels) + 2))
                np.testing.assert_allclose(r.score, expected, rtol=0, atol=1e-15)
                checked_scores += 1
    # Enumerate top-2 tie choices directly, without the metric implementation.
    panels = predictions[predictions.scenario.eq("grouped_host") & predictions.partition.eq("test")]
    checked_panels = 0
    for (method, host), scores in panels.groupby(["model", "bacterium"]):
        label = observed[observed.bacterium.eq(host)].set_index("phage").infection
        if len(label) != plan["expected_candidate_phages"]:
            continue
        boundary = scores.score.nlargest(2).iloc[-1]
        higher = scores.loc[scores.score.gt(boundary), "phage"].tolist()
        tied = scores.loc[scores.score.eq(boundary), "phage"].tolist()
        trials = [any(label.loc[p] for p in higher + list(chosen))
                  for chosen in combinations(tied, 2-len(higher))]
        expected = sum(trials)/len(trials)
        actual = metrics.loc[metrics.scenario.eq("grouped_host") & metrics.partition.eq("test")
                             & metrics.model.eq(method) & metrics.bacterium.eq(host), "hit@2"].item()
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-15)
        checked_panels += 1
    same = None
    if repeat is not None:
        other = read_json_object(repeat / "report.json")
        same = (report == other and all(sha256_file(repeat / n) == d
                                       for n, d in report["output_sha256"].items()))
        if not same:
            raise ValueError("repeat outputs differ")
    return {"status": "passed", "historical_per_host_metrics_exact": True,
            "reloaded_annotation_models_exact": reloads,
            "independent_scalar_neighbour_scores": checked_scores,
            "independent_top2_tie_panels": checked_panels,
            "second_run_all_report_and_output_hashes_equal": same,
            "report_sha256": sha256_file(run / "report.json")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--repeat", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = verify(Path(__file__).resolve().parents[1], args.data_root, args.run, args.repeat)
    atomic_write_json(result, args.out)
    print(result)


if __name__ == "__main__":
    main()
