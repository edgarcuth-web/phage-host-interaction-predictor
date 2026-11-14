"""The frozen Kaneko grade-3, fixed-library research comparison."""
from __future__ import annotations

from collections import Counter
import csv
import gzip
import hashlib
import importlib.metadata
import io
from itertools import product
import json
import math
from pathlib import Path
import resource
import struct
import sys
import tempfile

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .cohort import canonical_assembly_digest
from .host_baseline import fit_host_model, predict_host_model
from .host_genome import CANONICAL_WORDS, host_fourmers
from .retrieval import evaluate_retrieval, paired_retrieval_bootstrap

PROTOCOL_PATH = "data/studies/results/kaneko_baseline_freeze.v1/protocol.json"
PROTOCOL_SHA256 = "33f5b7e22b722fb2507f33d6646a84b58be7f8cc97b6d36902b796b454cd8f9a"
ENDPOINT = "kaneko-2025:complete-spot-clearing-grade3:v1"
SECONDARY_ENDPOINT = "kaneko-2025:visible-spot-clearing-grade-ge2:v1"
SCORE_SEMANTICS = "uncalibrated_source_grade3_complete_spot_clearing"
FEATURE_COLUMNS = ["host_4mer:" + word for word in CANONICAL_WORDS]
METHODS = ("prevalence", "host_4mer")
KS = (1, 2, 3, 4, 5)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_rows(path, rows):
    require(bool(rows), "Cannot emit an empty result table")
    with Path(path).open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(clean(rows))


def binding(path):
    data = Path(path).read_bytes()
    return {"sha256": digest(data), "size_bytes": len(data)}


def check_binding(path, expected):
    actual = binding(path)
    require(actual == {k: expected[k] for k in actual}, f"Changed bound file: {path}")


def environment():
    return {"python": sys.version.split()[0], **{name: importlib.metadata.version(name)
            for name in ("numpy", "pandas", "scipy", "scikit-learn", "threadpoolctl")}}


def check_memory():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = int(peak if sys.platform == "darwin" else peak * 1024)
    require(peak_bytes <= 4 * 1024**3, "Frozen 4 GiB peak-RSS budget exceeded")
    return peak_bytes


def load_inputs(root, grades_path):
    root = Path(root)
    protocol_path = root / PROTOCOL_PATH
    require(digest(protocol_path.read_bytes()) == PROTOCOL_SHA256, "Frozen resolved protocol changed")
    plan = json.loads(protocol_path.read_text())
    require(environment() == plan["future_fit_environment"], "Frozen numerical environment differs")
    design_path = root / plan["design_path"]
    require(digest(design_path.read_bytes()) == plan["design_sha256"], "Frozen pre-run design changed")
    design = json.loads(design_path.read_text())
    for name, sha in plan["frozen_math_code_git_blobs"].items():
        data = (root/name).read_bytes()
        require(hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest() == sha,
                f"Frozen mathematical helper changed: {name}")
    resolved = {}
    for role, item in design["inputs"].items():
        check_binding(root/item["path"], item)
        resolved[role] = json.loads((root/item["path"]).read_text())
    for name, item in plan["frozen_outputs"].items():
        check_binding(protocol_path.parent/name, item)
    hosts = read_rows(protocol_path.parent/"hosts.csv")
    host_map = {r["host_id"]: r for r in hosts}
    require(len(hosts) == len(host_map) == 143 and "127" not in host_map, "Host scope changed")
    assembly_map = {r["source_label"]: r["assembly_sha256"] for r in resolved["sequence_qa"]["qc"]
                    if r["kind"] == "host_genome"}
    require(all(r["assembly_sha256"] == assembly_map.get(r["host_id"]) for r in hosts),
            "Frozen source assembly ledger differs")
    require(len(plan["library"]) == len(set(plan["library"])) == 42 and plan["endpoint_id"] == ENDPOINT,
            "Candidate library or endpoint changed")
    group_partitions = {}
    for row in hosts:
        group_partitions.setdefault(row["host_group"], set()).add(row["partition"])
    require(len(group_partitions) == 27 and all(len(v) == 1 for v in group_partitions.values()),
            "Group leakage or group scope changed")
    require(digest(Path(grades_path).read_bytes()) == design["source_grades"]["sha256"],
            "Original grade export changed")
    raw = read_rows(grades_path)
    require(len(raw) == 6480, "Original source row count changed")
    seen, selected = set(), []
    for row in raw:
        key = row["host_label"], row["phage_label"]
        require(key not in seen and row["source_grade"] in {"0", "1", "2", "3"},
                "Duplicate, invalid or unknown source cell")
        seen.add(key)
        if key[0] not in host_map or key[1] not in plan["library"]:
            continue
        require(row["phage_role"] == "study" and row["source_sheet"] == "SF2", "Wrong source cell role")
        selected.append({"bacterium": key[0], "phage": key[1], "source_grade": int(row["source_grade"]),
                         "host_group": host_map[key[0]]["host_group"], "partition": host_map[key[0]]["partition"]})
    labels = pd.DataFrame(selected).sort_values(["bacterium", "phage"]).reset_index(drop=True)
    require(len(labels) == 6006 and set(zip(labels.bacterium, labels.phage)) ==
            set(product(host_map, plan["library"])), "Incomplete 143-by-42 source support")
    labels["infection"] = (labels.source_grade == 3).astype(int)  # Compatibility column, grade-3 target only.
    for part, expected in plan["split_support"].items():
        local = labels[labels.partition == part]
        require(local.bacterium.nunique() == expected["hosts"] and local.host_group.nunique() == expected["groups"]
                and len(local) == expected["source_pairs"] and int(local.infection.sum()) == expected["grade3_targets"],
                f"Frozen label/split support differs: {part}")
    return plan, design, resolved, hosts, labels


def fit_grade3(train, features, library):
    require((train.infection == (train.source_grade == 3).astype(int)).all(), "Training target is not source grade 3")
    require(set(train.phage) == set(library), "Incomplete training library")
    with threadpool_limits(limits=1):
        model = fit_host_model(train, features, FEATURE_COLUMNS, [], C=0.1, max_iter=2000)
    model["score_semantics"] = SCORE_SEMANTICS
    model["endpoint_id"] = ENDPOINT
    model["library"] = list(library)
    require(set(model["phages"]) == set(library), "Fit dropped a candidate")
    return model


def predict_grade3(model, pairs, features, library):
    require(model.get("endpoint_id") == ENDPOINT and model.get("score_semantics") == SCORE_SEMANTICS,
            "Wrong model endpoint or score semantics")
    require(model.get("library") == list(library) and set(model.get("phages", {})) == set(library),
            "Wrong fitted candidate library")
    with threadpool_limits(limits=1):
        scores = predict_host_model(model, pairs, features)
    require(np.isfinite(scores).all() and ((scores >= 0) & (scores <= 1)).all(), "Invalid candidate scores")
    return scores


def selection_mass(scores, k):
    values = np.asarray(scores, dtype=float)
    require(values.ndim == 1 and len(values) > 0 and np.isfinite(values).all() and k >= 1, "Bad selection input")
    k = min(k, len(values))
    boundary = np.sort(values)[-k]
    above, tied = values > boundary, values == boundary
    return above.astype(float) + tied * ((k - int(above.sum())) / int(tied.sum()))


def metric_contrast(left, right, groups, metric, uncertainty):
    eligible_only = metric.startswith("recall@")
    a, b = left.copy(), right.copy()
    if eligible_only:
        # Keep all groups in the sampling universe. Undefined ineligible recalls are
        # finite carriers only; paired_retrieval_bootstrap masks them after each draw.
        for frame in (a, b):
            require(frame.loc[frame.eligible, metric].notna().all(), "Undefined eligible recall")
            frame.loc[~frame.eligible.astype(bool), metric] = 0.0
    return paired_retrieval_bootstrap(a, b, groups, metric=metric, eligible_only=eligible_only,
              resamples=uncertainty["resamples"], seed=uncertainty["seed"], confidence=uncertainty["confidence"])


def check_implementation(root):
    path = Path(root)/"data/studies/kaneko_baseline.implementation.v1.json"
    receipt = json.loads(path.read_text())
    require(receipt["protocol_sha256"] == PROTOCOL_SHA256, "Implementation references a different protocol")
    for name, expected in receipt["code_sha256"].items():
        require(digest((Path(root)/name).read_bytes()) == expected, f"Implementation changed: {name}")
    return path, receipt


def fit_study(root, sequence_root, grades_path, out, implementation_commit):
    root, out, sequence_root = Path(root), Path(out), Path(sequence_root)
    require(not out.exists(), "Fit output directory must be new")
    require(len(implementation_commit) == 40, "Implementation commit is required")
    implementation_path, _ = check_implementation(root)
    plan, design, sources, hosts, labels = load_inputs(root, grades_path)
    out.mkdir(parents=True)
    vectors, feature_rows = {}, []
    qa = {r["source_label"]: r for r in sources["sequence_qa"]["qc"] if r["kind"] == "host_genome"}
    with tempfile.TemporaryDirectory(prefix="fasta_", dir=out) as temporary:
        fasta = Path(temporary)/"source.fasta"
        for index, source in enumerate(sources["sequence_intake"]["host_assemblies"], 1):
            item = source["fasta"]
            path = (sequence_root/item["filename"]).resolve()
            require(path.is_relative_to(sequence_root.resolve()), "Source path escapes its root")
            check_binding(path, item)
            compressed = path.read_bytes()
            require(hashlib.md5(compressed).hexdigest() == source["compressed_md5"], "Archive MD5 mismatch")
            payload = gzip.decompress(compressed)
            fasta.write_bytes(payload)
            identity = canonical_assembly_digest(fasta)
            host = source["source_label"]
            require(identity == qa[host]["assembly_sha256"], "Genome differs from frozen QC")
            vector, qc = host_fourmers(payload)
            require(len(vector) == 136 and np.isfinite(vector).all() and abs(float(vector.sum())-1) < 1e-12,
                    "Invalid host representation")
            vectors[host] = vector
            feature_rows.append({"host_id": host, "assembly_sha256": identity,
                "source_archive_sha256": item["sha256"],
                "feature_values_sha256": digest(struct.pack("<"+"d"*136, *vector)), **qc})
            check_memory()
            if index % 20 == 0:
                print(json.dumps({"stage": "features", "completed_hosts": index, "total_hosts": 143}), flush=True)
    require(set(vectors) == {r["host_id"] for r in hosts}, "Feature host scope differs")
    features = pd.DataFrame.from_dict(vectors, orient="index", columns=FEATURE_COLUMNS).sort_index()
    write_rows(out/"host_features.csv", [{"host_id": host, **dict(zip(CANONICAL_WORDS, vector))}
                                        for host, vector in features.iterrows()])
    write_rows(out/"feature_manifest.csv", sorted(feature_rows, key=lambda r: r["host_id"]))
    train = labels[labels.partition == "train"].reset_index(drop=True)
    model = fit_grade3(train, features, plan["library"])
    model["protocol_sha256"] = PROTOCOL_SHA256
    model["hyperparameters"] = {"C": 0.1, "solver": "lbfgs", "max_iter": 2000, "class_weight": None}
    model["calibrated"] = False
    require(model["training_hosts"] == sorted(train.bacterium.unique()), "Fitted host scope differs")
    stats = train.groupby("phage").infection.agg(["sum", "count"])
    prevalence = {"schema_version": 1, "endpoint_id": ENDPOINT, "score_semantics": SCORE_SEMANTICS,
                  "protocol_sha256": PROTOCOL_SHA256, "calibrated": False,
                  "library": plan["library"], "training_hosts": model["training_hosts"],
                  "phages": {p: {"positive_count": int(stats.loc[p,"sum"]),
                                 "observed_count": int(stats.loc[p,"count"]),
                                 "score": float((stats.loc[p,"sum"]+1)/(stats.loc[p,"count"]+2))}
                             for p in plan["library"]}}
    write_json(out/"host_4mer_model.json", model)
    write_json(out/"prevalence_model.json", prevalence)
    kinds = Counter(r["kind"] for r in model["phages"].values())
    require(sum(kinds.values()) == 42, "Fitted candidate count differs")
    fitted_files = {p.name: binding(p) for p in sorted(out.iterdir()) if p.is_file()}
    receipt = {"schema_version": 1, "study_id": design["study_id"], "stage": "fitted_no_held_out_scoring",
       "implementation_commit": implementation_commit, "implementation_receipt": binding(implementation_path),
       "protocol_sha256": PROTOCOL_SHA256, "environment": environment(),
       "training_hosts": 123, "training_groups": 17, "training_pairs": 5166,
       "feature_hosts": 143, "feature_dimensions": 136, "candidate_count": 42,
       "fitted_candidate_kinds": dict(kinds), "logistic_fits": kinds["logistic"],
       "single_class_fallbacks": kinds["single_class_prevalence"], "convergence_warnings": 0,
       "held_out_score_runs": 0, "files": fitted_files, "peak_rss_bytes": check_memory(),
       "cpu_threads": 1, "gpu_requests": 0, "paid_provider_requests": 0,
       "independent_release_validation_established": False}
    write_json(out/"FIT.json", receipt)
    print(json.dumps({"stage": "fit_complete", "candidate_kinds": dict(kinds), "held_out_score_runs": 0}), flush=True)
    return receipt


def read_fit(root, fit_dir, plan):
    fit_dir = Path(fit_dir)
    receipt = json.loads((fit_dir/"FIT.json").read_text())
    require(receipt["protocol_sha256"] == PROTOCOL_SHA256 and receipt["environment"] == environment(),
            "Fit provenance differs")
    implementation_path, _ = check_implementation(root)
    check_binding(implementation_path, receipt["implementation_receipt"])
    for name, item in receipt["files"].items():
        check_binding(fit_dir/name, item)
    rows = read_rows(fit_dir/"host_features.csv")
    require(len(rows) == 143 and len({r["host_id"] for r in rows}) == 143, "Feature cache scope differs")
    features = pd.DataFrame([[float(r[w]) for w in CANONICAL_WORDS] for r in rows],
                  index=[r["host_id"] for r in rows], columns=FEATURE_COLUMNS).sort_index()
    manifest = {r["host_id"]: r for r in read_rows(fit_dir/"feature_manifest.csv")}
    for host, row in features.iterrows():
        require(digest(struct.pack("<"+"d"*136, *row)) == manifest[host]["feature_values_sha256"],
                "Feature cache did not round-trip")
    model = json.loads((fit_dir/"host_4mer_model.json").read_text())
    prevalence = json.loads((fit_dir/"prevalence_model.json").read_text())
    for artifact in (model, prevalence):
        require(artifact["endpoint_id"] == ENDPOINT and artifact["score_semantics"] == SCORE_SEMANTICS
                and artifact["library"] == plan["library"] and artifact["protocol_sha256"] == PROTOCOL_SHA256,
                "Wrong artifact task")
    return receipt, features, model, prevalence


def evaluate_study(root, fit_dir, grades_path, out, fit_commit):
    root, out = Path(root), Path(out)
    require(not out.exists() and len(fit_commit) == 40, "New output and frozen fit commit required")
    plan, design, _, hosts, labels = load_inputs(root, grades_path)
    fit_receipt, features, model, prevalence = read_fit(root, fit_dir, plan)
    train_ids = sorted({r["host_id"] for r in hosts if r["partition"] == "train"})
    require(model["training_hosts"] == prevalence["training_hosts"] == train_ids, "Training host identity differs")
    held = labels[labels.partition.isin(["validation", "test"])].reset_index(drop=True)
    require(len(held) == 840, "Held-out scoring scope differs")
    score_vectors = {"prevalence": np.array([prevalence["phages"][p]["score"] for p in held.phage]),
                     "host_4mer": predict_grade3(model, held, features, plan["library"])}
    score_rows = []
    for method, values in score_vectors.items():
        require(len(values) == 840 and np.isfinite(values).all(), "Incomplete scores")
        for row, value in zip(held.itertuples(index=False), values):
            score_rows.append({"partition": row.partition, "host_id": row.bacterium,
                               "phage_id": row.phage, "method": method, "score_endpoint_id": ENDPOINT,
                               "score": float(value)})
    groups = {r["host_id"]: r["host_group"] for r in hosts}
    summaries, metric_rows, contrast_rows, phage_rows = {}, [], [], []
    for evaluation_endpoint, minimum in ((ENDPOINT, 3), (SECONDARY_ENDPOINT, 2)):
        summaries[evaluation_endpoint] = {}
        for part in ("validation", "test"):
            observed = held[held.partition == part].copy()
            positions = np.flatnonzero((held.partition == part).to_numpy())
            observed["infection"] = (observed.source_grade >= minimum).astype(int)
            per_method, summaries[evaluation_endpoint][part] = {}, {}
            for method, vector in score_vectors.items():
                scores = observed[["bacterium", "phage"]].copy()
                scores["score"] = vector[positions]
                summary, per_host = evaluate_retrieval(observed, scores,
                    hosts=sorted(observed.bacterium.unique()), candidate_phages=plan["library"], ks=KS)
                require(summary["unknown_pairs"] == 0 and not summary["incomplete_hosts_excluded"], "Unexpected missing labels")
                for k in (2, 5):
                    name = f"capped_tests@{k}"
                    per_host[name] = k-sum(per_host[f"hit@{j}"] for j in range(1,k))
                    summary["all_hosts"][name] = float(per_host[name].mean())
                per_method[method] = per_host
                summaries[evaluation_endpoint][part][method] = summary
                for row in per_host.to_dict("records"):
                    metric_rows.append({"evaluation_endpoint_id": evaluation_endpoint, "score_endpoint_id": ENDPOINT,
                        "partition": part, "method": method, "host_group": groups[row["bacterium"]], **row})
                selected = []
                for host in sorted(observed.bacterium.unique()):
                    local = scores[scores.bacterium == host]
                    mass = selection_mass(local.score.to_numpy(), 2)
                    local_y = observed[observed.bacterium == host].set_index("phage").infection
                    selected.extend({"phage": p, "selection_mass": float(m), "positive_mass": float(m)*int(local_y[p])}
                                    for p, m in zip(local.phage, mass))
                selection = pd.DataFrame(selected).groupby("phage").sum()
                counts = observed.groupby("phage").infection.agg(["sum", "count"])
                for phage in plan["library"]:
                    phage_rows.append({"evaluation_endpoint_id": evaluation_endpoint, "partition": part,
                        "method": method, "phage_id": phage, "observed_hosts": int(counts.loc[phage,"count"]),
                        "positive_hosts": int(counts.loc[phage,"sum"]),
                        "expected_top2_selections": float(selection.loc[phage,"selection_mass"]),
                        "expected_positive_top2_selections": float(selection.loc[phage,"positive_mass"])})
            metrics = [f"{kind}@{k}" for kind in ("hit", "recall") for k in (1,2,3,5)] + ["mrr", "capped_tests@2", "capped_tests@5"]
            for metric in metrics:
                contrast = metric_contrast(per_method["host_4mer"], per_method["prevalence"],
                                            groups, metric, design["uncertainty"])
                contrast_rows.append({"evaluation_endpoint_id": evaluation_endpoint, "partition": part,
                                     "left": "host_4mer", "right": "prevalence", **contrast})
    primary = [r for r in contrast_rows if r["evaluation_endpoint_id"] == ENDPOINT and r["partition"] == "test"]
    hit = next(r for r in primary if r["metric"] == "hit@2")
    cost = next(r for r in primary if r["metric"] == "capped_tests@2")
    decisions = {"hit_gain_supported": hit["interval"] is not None and hit["interval"][0] > 0,
                 "cost_reduction_supported": cost["interval"] is not None and cost["interval"][1] < 0}
    decisions["joint_benefit_supported"] = all(decisions.values())
    out.mkdir(parents=True)
    write_rows(out/"scores.csv", score_rows)
    write_rows(out/"per_host_metrics.csv", metric_rows)
    write_rows(out/"per_phage_diagnostics.csv", phage_rows)
    write_json(out/"contrasts.json", contrast_rows)
    report = {"schema_version": 1, "study_id": design["study_id"], "protocol_sha256": PROTOCOL_SHA256,
              "implementation_commit": fit_receipt["implementation_commit"], "frozen_fit_commit": fit_commit,
              "fit_receipt": binding(Path(fit_dir)/"FIT.json"), "environment": environment(),
              "score_rows": len(score_rows), "score_endpoints": [ENDPOINT], "evaluation_endpoints": [ENDPOINT,SECONDARY_ENDPOINT],
              "summaries": summaries, "primary_hit_contrast": hit, "primary_cost_contrast": cost,
              "decisions": decisions, "fit_summary": fit_receipt,
              "fit_calls_during_evaluation": 0, "peak_rss_bytes": check_memory(),
              "limitations": [plan["test_precision_limit"], "Inspected within-source research panel; no independent release validation.",
                              "Source-target scores are uncalibrated and do not estimate productive infection or treatment benefit."],
              "independent_release_validation_established": False}
    write_json(out/"report.json", report)
    write_json(out/"OUTPUTS.json", {p.name: binding(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps({"stage": "evaluation_complete", "score_rows": len(score_rows), "decisions": decisions,
                      "hit_contrast": clean(hit), "cost_contrast": clean(cost)}), flush=True)
    return report
