"""Offline, source-endpoint-specific comparison of two candidate discovery routes.

All scoring files must cover the *same complete* host x phage panel. Missing
observations are unknown, not negatives. Existing retrieval metrics use exact
tie expectations; membership/rescue use a separately declared lexicographic tie
break. This evaluator does not train a model or infer productive infection.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .artifacts import (
    atomic_write_json, code_provenance, read_json_object, sha256_file, staged_directory,
)
from .data import validate_interactions
from .retrieval import (
    DEFAULT_K, evaluate_retrieval, paired_retrieval_bootstrap,
)

SHA_RE = re.compile(r"^[0-9a-f]{64}$")
METHODS_REQUIRED = frozenset({"evidence", "learned"})
COLUMN_ORDER = ("bacterium", "phage", "score")


def _ids(values, label):
    if not isinstance(values, list) or not values or any(
        not isinstance(v, str) or not v.strip() for v in values
    ) or len(set(values)) != len(values):
        raise ValueError(f"{label} must be a nonempty list of unique string IDs")
    return sorted(values)


def _panel(hosts, phages):
    return pd.MultiIndex.from_product(
        [hosts, phages], names=["bacterium", "phage"],
    )


def _check_scores(data, hosts, phages, name):
    if not isinstance(data, pd.DataFrame) or any(k not in data for k in COLUMN_ORDER):
        raise ValueError(f"{name} must contain bacterium, phage and score columns")
    expected = _panel(hosts, phages)
    actual = pd.MultiIndex.from_frame(data[["bacterium", "phage"]])
    if actual.has_duplicates or set(actual) != set(expected):
        raise ValueError(f"{name} must cover exactly the complete fixed panel")
    try:
        numbers = pd.to_numeric(data["score"], errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} scores must be numeric") from exc
    if not np.isfinite(numbers).all():
        raise ValueError(f"{name} scores must be finite")
    clean = data[["bacterium", "phage"]].copy()
    clean["score"] = numbers
    return clean.set_index(["bacterium", "phage"]).reindex(expected)["score"]


def _union_order(rankings, schedule, phages):
    """Cycle a precommitted route schedule, skipping already emitted identities.

    Each turn takes the highest as-yet-unselected phage from that route. A full
    union ordering has one position per library phage; no extra screening slots.
    """
    selected, seen = [], set()
    pointers = {name: 0 for name in METHODS_REQUIRED}
    while len(selected) < len(phages):
        route = schedule[len(selected) % len(schedule)]
        order = rankings[route]
        index = pointers[route]
        while index < len(order) and order[index] in seen:
            index += 1
        if index >= len(order):
            # The branch is exhausted only when every phage is already selected.
            raise ValueError("union schedule could not fill the fixed panel")
        phage = order[index]
        pointers[route] = index + 1
        selected.append((phage, route))
        seen.add(phage)
    return selected


def compare_discovery_routes(
    observed, score_tables, *, hosts, phages, host_groups,
    ks=DEFAULT_K, union_schedule=("evidence", "learned"),
    bootstrap_resamples=2000, bootstrap_seed=20261024,
):
    """Compare exact tie-expected metrics and deterministic candidate rescue.

    Only complete-label hosts are scored. Unknown hosts are retained in the
    ranking audit, but earn no partial rescue credit.
    """
    hosts, phages = _ids(list(hosts), "hosts"), _ids(list(phages), "phages")
    if not isinstance(host_groups, dict) or set(host_groups) != set(hosts) or any(
        not isinstance(g, str) or not g.strip() for g in host_groups.values()
    ):
        raise ValueError("host_groups must cover all evaluation hosts")
    ks = tuple(ks)
    if not ks or any(type(k) is not int or k < 1 or k > len(phages) for k in ks) or len(set(ks)) != len(ks):
        raise ValueError("K values must be distinct integers within the phage library size")
    ks = tuple(sorted(ks))
    schedule = tuple(union_schedule)
    if not schedule or set(schedule) != METHODS_REQUIRED:
        raise ValueError("union_schedule must use both evidence and learned routes")
    if not isinstance(score_tables, dict) or not METHODS_REQUIRED.issubset(score_tables) or "union" in score_tables:
        raise ValueError("score_tables must include evidence and learned, but not union")
    if type(bootstrap_resamples) is not int or bootstrap_resamples < 2 or type(bootstrap_seed) is not int:
        raise ValueError("invalid bootstrap settings")
    validate_interactions(observed)
    expected = _panel(hosts, phages)
    known = pd.MultiIndex.from_frame(observed[["bacterium", "phage"]])
    if not set(known).issubset(set(expected)):
        raise ValueError("observations refer to pairs outside the fixed panel")
    labels = observed.set_index(["bacterium", "phage"])["infection"].to_dict()
    complete = {
        h: all((h, p) in labels for p in phages)
        for h in hosts
    }
    ordered_methods = sorted(score_tables)
    indexed = {
        name: _check_scores(frame, hosts, phages, name)
        for name, frame in score_tables.items()
    }
    ranked, ranking_rows = {}, []
    for host in hosts:
        rankings = {
            name: sorted(phages, key=lambda p: (-float(series.loc[(host, p)]), p))
            for name, series in indexed.items()
        }
        union = _union_order(rankings, schedule, phages)
        rankings["union"] = [p for p, _ in union]
        routes = {p: source for p, source in union}
        for name, ordered in sorted(rankings.items()):
            ranked[host, name] = ordered
            for position, phage in enumerate(ordered, 1):
                source_score = None if name == "union" else float(indexed[name].loc[(host, phage)])
                ranking_rows.append({
                    "bacterium": host, "method": name, "phage": phage,
                    "rank": position, "score": source_score,
                    "tie_policy": "deterministic_phage_id_ascending",
                    "selected_by_route": routes[phage] if name == "union" else name,
                    "observed_label": labels.get((host, phage)),
                    "label_status": "observed" if (host, phage) in labels else "unknown",
                })
    union_scores = pd.DataFrame(
        [(h, p, len(phages) - i) for h in hosts
         for i, p in enumerate(ranked[h, "union"])],
        columns=COLUMN_ORDER,
    )
    all_scores = {
        **{name: frame[["bacterium", "phage", "score"]].copy()
           for name, frame in score_tables.items()},
        "union": union_scores,
    }
    summaries, by_method, metrics_rows = {}, {}, []
    for name in [*ordered_methods, "union"]:
        summary, rows = evaluate_retrieval(
            observed, all_scores[name], hosts=hosts, candidate_phages=phages, ks=ks,
        )
        summaries[name] = summary
        by_method[name] = rows
        for row in rows.to_dict(orient="records"):
            metrics_rows.append({"method": name, **row})
    rescue_rows = []
    for host in hosts:
        for k in ks:
            evidence = set(ranked[host, "evidence"][:k])
            learned = set(ranked[host, "learned"][:k])
            union = set(ranked[host, "union"][:k])
            valid = complete[host]
            positives = {p for p in phages if labels.get((host, p)) == 1} if valid else set()
            rescue_rows.append({
                "bacterium": host, "k": k, "complete_panel": valid,
                "positive_candidates": len(positives) if valid else None,
                "evidence_selected": sorted(evidence),
                "learned_selected": sorted(learned),
                "union_selected": ranked[host, "union"][:k],
                "learned_only_positive_phages": sorted((learned - evidence) & positives) if valid else None,
                "union_only_positive_phages": sorted((union - evidence) & positives) if valid else None,
                "learned_only_positive_count": len((learned - evidence) & positives) if valid else None,
                "union_only_positive_count": len((union - evidence) & positives) if valid else None,
            })
    rescue = pd.DataFrame(rescue_rows)
    primary_k = 2 if 2 in ks else ks[0]
    support = rescue.loc[rescue.k.eq(primary_k) & rescue.complete_panel]
    positives = int(support.positive_candidates.sum())
    intervals = {
        name: paired_retrieval_bootstrap(
            by_method[name], by_method["evidence"], host_groups,
            metric=f"hit@{primary_k}",
            resamples=bootstrap_resamples, seed=bootstrap_seed,
        )
        for name in [*ordered_methods, "union"] if name != "evidence"
    }
    report = {
        "research_status": "offline_development_evaluation_only",
        "scope": "named_source_endpoint_not_mechanistic_infection",
        "hosts": len(hosts), "phages": len(phages),
        "ks": list(ks), "primary_k": primary_k,
        "membership_tie_policy": "score_descending_then_phage_id_ascending",
        "metric_tie_policy": "exact expectation over uniform within-score-tie permutations (union ordering deterministic)",
        "union_schedule": list(schedule),
        "observed_pairs": len(observed),
        "unknown_pairs": len(expected) - len(observed),
        "excluded_incomplete_hosts": sorted(h for h in hosts if not complete[h]),
        "method_summaries": summaries,
        "paired_hit_bootstrap_vs_evidence": intervals,
        "primary_k_rescue": {
            "complete_hosts": len(support), "observed_positive_phages": positives,
            "learned_only_positive_instances": int(support.learned_only_positive_count.sum()),
            "union_only_positive_instances": int(support.union_only_positive_count.sum()),
            "learned_only_fraction_of_observed_positives": (
                float(support.learned_only_positive_count.sum() / positives) if positives else None
            ),
            "union_only_fraction_of_observed_positives": (
                float(support.union_only_positive_count.sum() / positives) if positives else None
            ),
        },
        "limitations": [
            "Results apply only to source-defined observed assay labels in the declared library.",
            "Unmeasured/unknown pairs are never converted to negative labels.",
            "Unknown-label hosts do not contribute to rescue or metric estimates.",
            "Tie-expected retrieval metrics and deterministic membership/rescue are distinct estimands.",
            "A learned-only positive is not evidence of a previously unknown causal gene.",
            "No cocktail efficacy, productive-infection probability, clinical outcome or unseen-phage generalization is inferred.",
            "Declared training_hosts checks cannot prove absence of training leakage.",
        ],
    }
    return report, pd.DataFrame(ranking_rows), rescue, pd.DataFrame(metrics_rows)


def _input_path(root, entry, endpoint, label):
    if not isinstance(entry, dict) or not all(k in entry for k in ("path", "sha256", "endpoint")):
        raise ValueError(f"{label} must supply path, sha256 and endpoint")
    if entry["endpoint"] != endpoint:
        raise ValueError(f"{label} declares a different source endpoint")
    raw = entry["path"]
    if not isinstance(raw, str) or not raw or Path(raw).is_absolute() or ".." in Path(raw).parts:
        raise ValueError(f"{label} path must be relative and confined to manifest directory")
    source = (root / raw).resolve()
    if not source.is_relative_to(root.resolve()) or not source.is_file():
        raise ValueError(f"{label} source is outside manifest directory or missing")
    if not isinstance(entry["sha256"], str) or not SHA_RE.fullmatch(entry["sha256"]):
        raise ValueError(f"{label} requires exact lowercase SHA-256")
    if sha256_file(source) != entry["sha256"]:
        raise ValueError(f"{label} changed frozen source bytes")
    return source


def run_discovery_audit(manifest_path: str | Path, outdir: str | Path):
    """Run locked file-based evaluation with staged, non-overwriting outputs."""
    manifest_path, outdir = Path(manifest_path).resolve(), Path(outdir)
    if outdir.exists():
        raise ValueError("output directory must be new")
    manifest = read_json_object(manifest_path, label="retrieval benchmark manifest")
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("study_id"), str) or not manifest["study_id"].strip():
        raise ValueError("unsupported manifest schema or study_id")
    endpoint = manifest.get("endpoint")
    if not isinstance(endpoint, str) or not endpoint.strip():
        raise ValueError("manifest must name an exact source endpoint")
    if manifest.get("research_status") != "development_only":
        raise ValueError("only explicitly scoped development audits are supported")
    hosts, phages = _ids(manifest.get("hosts"), "hosts"), _ids(manifest.get("phages"), "phages")
    groups = manifest.get("host_groups")
    if not isinstance(groups, dict) or set(groups) != set(hosts) or any(
        not isinstance(g, str) or not g.strip() for g in groups.values()
    ):
        raise ValueError("host_groups must cover the exact host manifest")
    sources = manifest.get("score_sources")
    if not isinstance(sources, dict) or not METHODS_REQUIRED.issubset(sources) or "union" in sources:
        raise ValueError("score_sources must include evidence and learned without union")
    root = manifest_path.parent
    labels_path = _input_path(root, manifest.get("observations"), endpoint, "observations")
    observed = pd.read_csv(labels_path)
    score_tables = {}
    for name, entry in sources.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(entry, dict):
            raise ValueError("invalid score_source method")
        producer = entry.get("producer")
        if not isinstance(producer, str) or not producer.strip():
            raise ValueError(f"{name} requires named producer provenance")
        trained = entry.get("training_hosts")
        if not isinstance(trained, list) or any(not isinstance(h, str) for h in trained):
            raise ValueError(f"{name} must declare training_hosts, possibly empty")
        if set(trained).intersection(hosts):
            raise ValueError(f"{name} training_hosts overlap evaluation hosts")
        score_tables[name] = pd.read_csv(_input_path(root, entry, endpoint, name))
    ks = manifest.get("ks", list(DEFAULT_K))
    schedule = manifest.get("union_schedule")
    if not isinstance(ks, list) or not isinstance(schedule, list):
        raise ValueError("ks and union_schedule must be declared as lists")
    bootstrap = manifest.get("bootstrap", {"resamples": 2000, "seed": 20261024})
    if not isinstance(bootstrap, dict) or set(bootstrap) != {"resamples", "seed"}:
        raise ValueError("bootstrap must declare resamples and seed")
    report, rankings, rescue, metrics = compare_discovery_routes(
        observed, score_tables, hosts=hosts, phages=phages, host_groups=groups,
        ks=ks, union_schedule=schedule,
        bootstrap_resamples=bootstrap["resamples"],
        bootstrap_seed=bootstrap["seed"],
    )
    report.update({
        "study_id": manifest["study_id"], "endpoint": endpoint,
        "manifest_sha256": sha256_file(manifest_path),
        "source_sha256": {
            "observations": manifest["observations"]["sha256"],
            **{name: entry["sha256"] for name, entry in sorted(sources.items())},
        },
        "score_producers": {name: entry["producer"] for name, entry in sorted(sources.items())},
    })
    repo = Path(__file__).resolve().parents[2]
    report["code_provenance"] = code_provenance(repo, repo / "scripts/run_retrieval_rescue.py")
    with staged_directory(outdir) as stage:
        rankings.to_csv(stage / "ranked_candidates.csv", index=False)
        rescue.to_json(stage / "rescue_by_host.json", orient="records", indent=2)
        metrics.to_csv(stage / "per_host_metrics.csv", index=False)
        report["output_sha256"] = {
            name: sha256_file(stage / name)
            for name in ("ranked_candidates.csv", "rescue_by_host.json", "per_host_metrics.csv")
        }
        atomic_write_json(report, stage / "report.json")
    return report
