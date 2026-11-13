"""Independently verify budget metrics from first-positive rank probabilities."""
import argparse
from fractions import Fraction
from itertools import groupby
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
import run_screening_budget_audit as audit


def first_hit_metrics(labels, scores, budgets):
    """Construct the exact first-positive PMF; do not reuse retrieval metrics."""
    ordered = sorted(zip(scores, labels), key=lambda item: -item[0])
    offset, distribution = 0, {}
    for _, block in groupby(ordered, key=lambda item: item[0]):
        values = list(block)
        n, p = len(values), sum(label for _, label in values)
        if p:
            distribution = {offset+j: Fraction(comb(n-j, p-1), comb(n, p))
                            for j in range(1, n-p+2)}
            break
        offset += n
    result = {"mrr": float(sum((v/rank for rank, v in distribution.items()), Fraction()))}
    for k in budgets:
        result[f"hit@{k}"] = float(sum((v for rank, v in distribution.items() if rank <= k), Fraction()))
        result[f"tests@{k}"] = (float(sum((v*min(rank, k) for rank, v in distribution.items()), Fraction()))
                                 if distribution else float(k))
    return result


def require_close(actual, expected, label):
    if not np.isfinite(actual) or not np.isclose(actual, expected, rtol=0, atol=1e-12):
        raise ValueError(f"independent numeric check failed: {label}")


def verify(run, repo, data_root):
    run, repo = Path(run), Path(repo)
    audit.verify_implementation(repo)
    report = read_json_object(run / "report.json")
    if (report["protocol_sha256"] != sha256_file(repo / audit.PROTOCOL)
            or report["implementation_sha256"] != sha256_file(repo / audit.IMPLEMENTATION)):
        raise ValueError("report freeze identities differ")
    audit.verify_hashes(run, report["output_sha256"])
    if report["new_learned_model_fits"] != 0 or report["new_provider_requests"] != 0:
        raise ValueError("unexpected execution accounting")
    plan, phages, panels = audit.load_study(repo, data_root)
    per_host = audit.read_frame(run / "per_host_metrics.csv")
    curves = audit.read_frame(run / "curves.csv")
    contrasts = audit.read_frame(run / "contrasts.csv")
    availability = audit.read_frame(run / "availability.csv")
    if (per_host.duplicated(["panel", "partition", "model", "bacterium"]).any()
            or len(per_host) != report["per_host_rows"] or len(curves) != report["curve_rows"]
            or len(contrasts) != report["contrast_rows"]
            or contrasts.duplicated(["panel", "partition", "model", "denominator", "metric"]).any()):
        raise ValueError("duplicate or incomplete audit output")
    support = read_json_object(run / "support.json")["panels"]
    checked, intervals_checked = 0, 0
    for panel in panels:
        section = per_host[per_host.panel.eq(panel["id"]) & per_host.partition.eq(panel["partition"])]
        labels = {(r.bacterium, r.phage): int(r.infection) for r in panel["observed"].itertuples()}
        complete = [h for h in panel["hosts"] if all((h, p) in labels for p in phages)]
        expected_methods = set(panel["references"])
        if set(section.model) != expected_methods:
            raise ValueError("method support differs")
        counts = [s for s in support if s["panel"] == panel["id"] and s["partition"] == panel["partition"]]
        if len(counts) != 1 or counts[0]["incomplete_hosts_excluded"] != sorted(set(panel["hosts"])-set(complete)):
            raise ValueError("incomplete-host ledger differs")
        frames = {}
        for method in sorted(expected_methods):
            values = section[section.model.eq(method)].set_index("bacterium").sort_index()
            if set(values.index) != set(complete):
                raise ValueError("complete-host support differs")
            frames[method] = values
            source = panel["scores"].get(method)
            reference = panel["references"][method].set_index("bacterium")
            for host, row in values.iterrows():
                y = [labels[(host, p)] for p in phages]
                if row.positive_candidates != sum(y) or bool(row.eligible) != bool(sum(y)) or row.host_group != panel["groups"][host]:
                    raise ValueError("label or group identity differs")
                if source is not None:
                    s = source[source.bacterium.eq(host)].set_index("phage").score
                    expected = first_hit_metrics(y, [float(s[p]) for p in phages], range(1, max(plan["budgets"])+1))
                    expected = {k: v for k, v in expected.items()
                                if not k.startswith("tests@") or int(k.split("@")[1]) in plan["budgets"]}
                    for metric, value in expected.items():
                        require_close(row[metric], value, f"{method}/{host}/{metric}")
                    if row.source_kind != "full_scores":
                        raise ValueError("full-score availability mislabeled")
                else:
                    for metric in audit.HISTORICAL:
                        require_close(row[metric], float(reference.loc[host, metric]), metric)
                    for k in range(1, max(plan["budgets"])+1):
                        if k not in plan["historical_hit_budgets"] and pd.notna(row[f"hit@{k}"]):
                            raise ValueError("missing historical probabilities were filled")
                    for k in plan["budgets"]:
                        earlier = [j for j in range(1, k) if j not in plan["historical_hit_budgets"]]
                        if earlier:
                            if pd.notna(row[f"tests@{k}"]):
                                raise ValueError("unsupported capped cost was filled")
                        else:
                            expected = sum([1., *[1-float(reference.loc[host, f"hit@{j}"]) for j in range(1, k)]])
                            require_close(row[f"tests@{k}"], expected, "historical capped consumption")
                    if row.source_kind != "published_metrics":
                        raise ValueError("historical availability mislabeled")
                checked += 1
            for k in plan["budgets"]:
                a = availability[availability.panel.eq(panel["id"]) & availability.partition.eq(panel["partition"])
                    & availability.model.eq(method) & availability.budget.eq(k)]
                if len(a) != 1 or bool(a.hit.iloc[0]) != bool(values[f"hit@{k}"].notna().all()) or bool(a.tests.iloc[0]) != bool(values[f"tests@{k}"].notna().all()):
                    raise ValueError("availability ledger differs")
            for denominator in ("all_complete_hosts", "eligible_complete_hosts"):
                selected = values[values.eligible] if denominator.startswith("eligible") else values
                c = curves[curves.panel.eq(panel["id"]) & curves.partition.eq(panel["partition"])
                           & curves.model.eq(method) & curves.denominator.eq(denominator)]
                if set(c.budget) != set(plan["budgets"]) or len(c) != len(plan["budgets"]):
                    raise ValueError("curve support differs")
                for row in c.itertuples():
                    if row.hosts != len(selected):
                        raise ValueError("curve denominator differs")
                    for output, metric in (("hit_probability", f"hit@{row.budget}"),
                                           ("expected_tests_capped", f"tests@{row.budget}"), ("mrr", "mrr")):
                        mean = selected[metric].mean()
                        if pd.isna(mean):
                            if pd.notna(getattr(row, output)):
                                raise ValueError("unavailable curve became a number")
                        else:
                            require_close(getattr(row, output), mean, output)
        # Independently expand every group draw to host indices. This does not
        # use the runner's group-sum/multiplicity matrix implementation.
        hosts = sorted(complete)
        names = sorted({panel["groups"][h] for h in hosts})
        indices = [np.array([i for i, h in enumerate(hosts) if panel["groups"][h] == g]) for g in names]
        for denominator in ("all_complete_hosts", "eligible_complete_hosts"):
            rows = contrasts[contrasts.panel.eq(panel["id"]) & contrasts.partition.eq(panel["partition"])
                             & contrasts.denominator.eq(denominator)]
            metrics = [*[f"hit@{k}" for k in plan["budgets"]],
                       *[f"tests@{k}" for k in plan["budgets"]], "mrr"]
            expected_keys = {(method, metric) for method, values in frames.items()
                             if method != plan["baseline"] for metric in metrics
                             if values[metric].notna().all()}
            if set(zip(rows.model, rows.metric)) != expected_keys:
                raise ValueError("paired contrast coverage differs")
            reference = frames[plan["baseline"]].loc[hosts]
            columns = [frames[row.model].loc[hosts, row.metric].to_numpy(float)
                       - reference[row.metric].to_numpy(float) for row in rows.itertuples()]
            differences = np.column_stack(columns)
            keep = reference.eligible.to_numpy(bool) if denominator.startswith("eligible") else np.ones(len(hosts), bool)
            settings = plan["statistics"]["bootstrap"]
            rng = np.random.default_rng(settings["seed"])
            samples = []
            for _ in range(settings["resamples"]):
                positions = np.concatenate([indices[i] for i in rng.integers(len(names), size=len(names))])
                positions = positions[keep[positions]]
                if len(positions):
                    samples.append(differences[positions].mean(axis=0))
            alpha = (1-settings["confidence"])/2
            bounds = np.quantile(samples, [alpha, 1-alpha], axis=0)
            points = differences[keep].mean(axis=0)
            for i, row in enumerate(rows.itertuples()):
                require_close(row.delta_model_minus_reference, points[i], "paired point")
                require_close(row.ci_low, bounds[0, i], "paired lower interval")
                require_close(row.ci_high, bounds[1, i], "paired upper interval")
                if row.valid_resamples != len(samples):
                    raise ValueError("bootstrap resample support differs")
                intervals_checked += 1
    return {"status": "verified", "report_sha256": sha256_file(run / "report.json"),
            "protocol_sha256": report["protocol_sha256"],
            "verification_script_sha256": sha256_file(Path(__file__)),
            "per_host_rows_independently_checked": checked, "contrasts_independently_checked": intervals_checked,
            "arithmetic": "Exact Fraction first-positive PMF; independently expanded paired group bootstrap",
            "new_learned_model_fits": 0, "new_provider_requests": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/screening_budget"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("verification output already exists")
    result = verify(args.run, Path(__file__).resolve().parents[1], args.data_root)
    atomic_write_json(result, args.out)
    print(result)


if __name__ == "__main__":
    main()
