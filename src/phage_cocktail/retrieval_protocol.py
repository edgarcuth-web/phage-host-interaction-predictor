"""Read-only integrity and compatibility checks for the next retrieval study."""
from __future__ import annotations

from collections import Counter
import csv
import gzip
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import struct

PROTOCOL = "data/studies/new_host_retrieval.v1.json"
PROTOCOL_SHA256 = "d75017e4831bad99dbb4fabd5e8f0a1c0c99b7b1eeb20ff3162d558c983ed289"


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _path(root, name):
    root = Path(root).resolve()
    path = (root / name).resolve()
    _require(path.is_relative_to(root), f"input escapes repository: {name}")
    return path


def _rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _unique(rows, keys, name):
    index = {tuple(row[k] for k in keys): row for row in rows}
    _require(len(index) == len(rows), f"duplicate {name} identity")
    return index


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def validate_inputs(root, plan, *, annotation_data_root=None):
    """Validate bound bytes and scientific scopes without fitting or scoring."""
    root = Path(root).resolve()
    paths = {}
    for role, item in plan["inputs"].items():
        path = _path(root, item["path"])
        data = path.read_bytes()
        _require(len(data) == item["size_bytes"] and sha256(data).hexdigest() == item["sha256"],
                 f"changed input: {role} ({item['path']})")
        paths[role] = path
    for name, digest in plan["frozen_code_sha256"].items():
        _require(sha256(_path(root, name).read_bytes()).hexdigest() == digest,
                 f"changed frozen scorer/evaluator: {name}")

    expected = plan["dataset"]
    hosts = _unique(_rows(paths["host_manifest"]), ("bacterium",), "host")
    library = _unique(_rows(paths["library"]), ("phage",), "library")
    host_ids = {key[0] for key in hosts}
    phage_ids = {key[0] for key in library}
    _require(len(hosts) == expected["hosts"] and len(library) == expected["phages"],
             "host/library roster differs")
    _require(len({r["host_group"] for r in hosts.values()}) == expected["host_groups"],
             "host group count differs")
    _require(len({r["assembly_sha256"] for r in library.values()}) == expected["unique_phage_assemblies"],
             "unique phage assembly count differs")

    observed = _rows(paths["observations"])
    obs = _unique(observed, ("bacterium", "phage"), "observation")
    unknown = _unique(_rows(paths["unknown_exclusions"]), ("bacterium", "phage"), "unknown")
    universe = set(product(host_ids, phage_ids))
    _require(not (set(obs) & set(unknown)) and set(obs) | set(unknown) == universe,
             "observed/unknown labels do not partition the fixed library")
    _require(len(obs) == expected["observed_pairs"] and len(unknown) == expected["unknown_pairs"],
             "observed/unknown counts differ")
    for pair, row in obs.items():
        mlc = float(row["source_mlc"])
        _require(math.isfinite(mlc) and row["infection"] in ("0", "1")
                 and row["source_binary"] == row["infection"]
                 and int(row["infection"]) == int(mlc > 0), f"source MLC label differs: {pair}")
    _require(sum(int(r["infection"]) for r in observed) == expected["positives"], "positive count differs")
    _require(all(r["reason"] == "source_mlc_missing" for r in unknown.values()),
             "unknown exclusion reason differs")

    split_rows = _rows(paths["scenario_splits"])
    _unique(split_rows, ("scenario", "bacterium"), "scenario host")
    _require({r["scenario"] for r in split_rows} == set(plan["panels"]), "scenario roster differs")
    panels, support = {}, {}
    for scenario, spec in plan["panels"].items():
        rows = [r for r in split_rows if r["scenario"] == scenario]
        _require({r["bacterium"] for r in rows} == host_ids, f"host scope differs: {scenario}")
        _require(dict(Counter(r["split"] for r in rows)) == spec["hosts"], f"split counts differ: {scenario}")
        partitions = {}
        for row in rows:
            original = hosts[row["bacterium"],]
            _require(all(row[k] == original[k] for k in ("host_group", "assembly_sha256")),
                     f"host identity differs: {scenario}/{row['bacterium']}")
            if scenario == "grouped_host":
                _require(row["split"] == original["split"], "primary host split differs")
            partitions.setdefault(row["host_group"], set()).add(row["split"])
        _require(all(len(values) == 1 for values in partitions.values()), f"host group leakage: {scenario}")
        groups = {part: len({r["host_group"] for r in rows if r["split"] == part}) for part in spec["hosts"]}
        _require(groups == spec["groups"], f"split group counts differ: {scenario}")
        panels[scenario] = {part: {r["bacterium"] for r in rows if r["split"] == part} for part in spec["hosts"]}
        train = [r for r in observed if r["bacterium"] in panels[scenario]["train"]]
        _require(len(train) == plan["panel_support"][scenario]["train_rows"]
                 and {r["phage"] for r in train} == phage_ids, f"training label support differs: {scenario}")
        held = {}
        for part, want in plan["panel_support"][scenario].items():
            if part == "train_rows":
                continue
            ids = panels[scenario][part]
            complete = {h for h in ids if all((h, p) in obs for p in phage_ids)}
            eligible = {h for h in complete if any(int(obs[h, p]["infection"]) for p in phage_ids)}
            actual = dict(scored_hosts=len(ids), complete_hosts=len(complete), eligible_hosts=len(eligible),
                          unknown_pairs=sum((h, p) in unknown for h in ids for p in phage_ids))
            _require(all(actual[k] == want[k] for k in actual), f"evaluation support differs: {scenario}/{part}")
            _require(sorted(ids - complete) == want.get("excluded_hosts", []), f"incomplete-host exclusions differ: {scenario}/{part}")
            held[part] = {**actual, "excluded_hosts": sorted(ids - complete)}
        support[scenario] = dict(train_hosts=len(panels[scenario]["train"]), train_rows=len(train), held=held)

    inner = _rows(paths["b2_selection"])
    _unique(inner, ("bacterium",), "B2 selection host")
    b2_train = panels["phylogroup:B2"]["train"]
    _require({r["bacterium"] for r in inner} == b2_train, "B2 selection must contain only all outer-training hosts")
    groups = sorted({hosts[h,]["host_group"] for h in b2_train},
                    key=lambda g: (sha256(f"20261113:{g}".encode()).hexdigest(), g))
    val_groups = set(groups[:len(groups) // 5])
    for row in inner:
        group = hosts[row["bacterium"],]["host_group"]
        _require(row["host_group"] == group and row["selection_split"] == ("validation" if group in val_groups else "fit"),
                 "B2 inner split/group differs from frozen label-free rule")
    _require(dict(Counter(r["selection_split"] for r in inner)) == plan["b2_inner_selection"]["hosts"],
             "B2 selection host counts differ")
    inner_support = {}
    for part, want in plan["b2_inner_selection"]["support"].items():
        ids = {r["bacterium"] for r in inner if r["selection_split"] == part}
        complete = {h for h in ids if all((h, p) in obs for p in phage_ids)}
        eligible = {h for h in complete if any(int(obs[h, p]["infection"]) for p in phage_ids)}
        actual = dict(hosts=len(ids), rows=sum(r["bacterium"] in ids for r in observed),
                      complete_hosts=len(complete), eligible_hosts=len(eligible),
                      unknown_pairs=sum((h, p) in unknown for h in ids for p in phage_ids),
                      excluded_hosts=sorted(ids - complete))
        _require(actual == want, f"B2 inner label support differs: {part}")
        inner_support[part] = actual

    resources = json.loads(paths["resource_lock"].read_text())["resources"]
    resources = {(r["kind"], r["id"]): r for r in resources}
    vectors = _rows(paths["sequence_vectors"])
    vector_map = _unique(vectors, ("entity_kind", "entity_id"), "sequence vector")
    manifest = _unique(_rows(paths["sequence_manifest"]), ("entity_kind", "entity_id"), "sequence receipt")
    rc = str.maketrans("ACGT", "TGCA")
    words = sorted({min(w, w.translate(rc)[::-1]) for w in map("".join, product("ACGT", repeat=4))})
    _require(len(words) == expected["feature_dimensions"], "feature dimension differs")
    _require(list(vectors[0]) == ["entity_kind", "entity_id", *words], "feature order/schema differs")
    expected_vectors = {("host_genome", h) for h in host_ids} | {("phage_genome", p) for p in phage_ids}
    _require(set(vector_map) == set(manifest) == expected_vectors, "sequence vector/receipt scope differs")
    for key, row in vector_map.items():
        receipt = manifest[key]
        identity = hosts[key[1],] if key[0] == "host_genome" else library[key[1],]
        assembly, representation = identity["assembly_sha256"], expected["representation_sha256"]
        feature_key = sha256(_canonical(dict(assembly_sha256=assembly, representation_sha256=representation)).encode()).hexdigest()
        _require(receipt["assembly_sha256"] == assembly and receipt["representation_sha256"] == representation
                 and receipt["feature_key"] == feature_key and receipt["raw_sha256"] == resources[key]["sha256"],
                 f"sequence receipt identity differs: {key}")
        values = [float(row[w]) for w in words]
        _require(all(math.isfinite(v) and v >= 0 for v in values) and abs(math.fsum(values) - 1) <= 1e-12,
                 f"invalid sequence frequencies: {key}")
        digest = sha256(struct.pack("<" + "d" * len(values), *values)).hexdigest()
        _require(digest == receipt["feature_values_sha256"], f"sequence vector values differ: {key}")

    weights = json.loads(paths["annotation_weights"].read_text())
    _require(weights["provenance"]["protocol_sha256"] == plan["inputs"]["host_protocol"]["sha256"]
             and weights["provenance"]["resource_lock_sha256"] == plan["inputs"]["resource_lock"]["sha256"],
             "annotation weight source protocol/resources differ")
    train_hosts = panels["grouped_host"]["train"]
    _require(set(weights["training_hosts"]) == train_hosts and len(weights["training_hosts"]) == len(train_hosts),
             "annotation weights have incompatible training hosts")
    _require(set(weights["phages"]) == phage_ids and weights["score_semantics"] == "uncalibrated_source_mlc_positive",
             "annotation weights have incompatible phages/endpoint")
    counts = Counter(r["phage"] for r in observed if r["bacterium"] in train_hosts)
    _require(all(weights["phages"][p]["training_rows"] == n for p, n in counts.items()),
             "annotation weights have incompatible training label counts")

    score_support = {}
    for scenario, role in (("grouped_host", "genophi_primary_scores"), ("phylogroup:B2", "genophi_b2_scores")):
        rows = _rows(paths[role])
        scores = _unique(rows, ("bacterium", "phage"), "GenoPHI score")
        targets = set().union(*(ids for part, ids in panels[scenario].items() if part in ("validation", "test")))
        _require(set(scores) == set(product(targets, phage_ids)), f"incomplete or train-host GenoPHI score support: {scenario}")
        _require(all(math.isfinite(float(r["score"])) and 0 <= float(r["score"]) <= 1 for r in rows),
                 f"invalid GenoPHI score: {scenario}")
        score_support[scenario] = dict(hosts=len(targets), phages=len(phage_ids), scores=len(rows))

    prior_fit = [r for r in _rows(paths["paired_split"]) if r["partition"] == "fit"]
    old_phages = {r["phage"] for r in prior_fit}
    _require(len(prior_fit) == 6834 and len(old_phages) == 65
             and len({library[p,]["assembly_sha256"] for p in old_phages}) == 61,
             "prior crossed-fit scope differs; revisit compatibility decision")
    prior_models = json.loads(paths["paired_models"].read_text())
    _require(prior_models["encoders"]["development"]["phage"]["unique_fit_assemblies"] == 61
             and prior_models["encoders"]["full_refit"]["host"]["unique_fit_assemblies"] == 177,
             "prior model scaler scopes differ; revisit compatibility decision")
    ledger = _rows(paths["compatibility"])
    _unique(ledger, ("panel", "method"), "compatibility decision")
    required_methods = {"random", "phage_prevalence", "annotation", "host_4mer", "genophi", "phage_only", "additive", "bilinear"}
    _require({(r["panel"], r["method"]) for r in ledger} == set(product(plan["panels"], required_methods)),
             "compatibility ledger scope differs")

    annotation_checked = 0
    if annotation_data_root is not None:
        for item in plan["annotation_resources"]:
            data = _path(annotation_data_root, item["path"]).read_bytes()
            _require(len(data) == item["size_bytes"] and sha256(data).hexdigest() == item["sha256"],
                     f"changed annotation source: {item['id']}")
            annotation_checked += 1
    return dict(status="verified_inputs_and_compatibility_only", source_revision=plan["source_revision"],
                bound_inputs=len(paths), frozen_code_files=len(plan["frozen_code_sha256"]),
                cohort={k: expected[k] for k in ("hosts", "host_groups", "phages", "unique_phage_assemblies", "observed_pairs", "positives", "unknown_pairs")},
                verified_sequence_vectors=len(vectors), feature_dimensions=len(words), panels=support,
                b2_inner_hosts=dict(Counter(r["selection_split"] for r in inner)),
                b2_inner_support=inner_support,
                reusable_annotation_weights=dict(hosts=len(train_hosts), phages=len(counts), training_rows=sum(counts.values())),
                existing_genophi_scores=score_support,
                prior_crossed_fit=dict(rows=len(prior_fit), named_phages=len(old_phages), unique_phage_assemblies=61, reusable_for_matched_fit=False),
                compatibility_rows=len(ledger), annotation_source_files_checked=annotation_checked,
                new_model_fits=0, new_candidate_scores=0, new_genome_downloads=0,
                independent_release_validation=False, training_scoring_runner="pending_step_2")


def check_protocol(root, *, annotation_data_root=None):
    root = Path(root).resolve()
    payload = _path(root, PROTOCOL).read_bytes()
    digest = sha256(payload).hexdigest()
    _require(digest == PROTOCOL_SHA256, "protocol differs from its pre-fit freeze")
    report = validate_inputs(root, json.loads(payload), annotation_data_root=annotation_data_root)
    return dict(protocol_sha256=digest, checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), **report)
