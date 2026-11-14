"""Compute label-free host groups and freeze the 42-candidate baseline protocol."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import re
import resource
import subprocess
import time

from phage_cocktail.cohort import DisjointSet, canonical_assembly_digest

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "data/studies/kaneko_baseline.design.v1.json"
DESIGN_SHA256 = "90731654b84a91466672b49fc50948eaef6d8d7679653380a35eb9e45e30fe54"
DESIGN_PRE_RUN_COMMIT = "2da042742fca9042f3053d1ae390ac30b6c92028"
FROZEN_CODE_GIT_BLOBS = {
    "src/phage_cocktail/cohort.py": "781a43883474851a47cfba4dd480924abbe888e7",
    "src/phage_cocktail/host_genome.py": "744f65f9048891cfa72e6727a5bd03e35a3ecd63",
    "src/phage_cocktail/host_baseline.py": "593622089294f278a75c55aaf5a99c68c7ebd7f4",
    "src/phage_cocktail/retrieval.py": "213942caa4b96c342db59204badd3b6c29c66d3a",
}
FUTURE_FIT_ENVIRONMENT = {
    "python": "3.12.14", "numpy": "2.3.5", "pandas": "2.2.3",
    "scipy": "1.17.0", "scikit-learn": "1.8.0", "threadpoolctl": "3.6.0",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def bound_json(binding):
    data = (ROOT / binding["path"]).read_bytes()
    require(len(data) == binding["size_bytes"] and sha256(data) == binding["sha256"],
            f"Input changed: {binding['path']}")
    return json.loads(data)


def primary_target(grade):
    require(type(grade) is str and grade in {"0", "1", "2", "3"},
            "Unknown or invalid source grade; cannot create a target")
    return int(grade == "3")


def read_pairs(path, host_ids):
    """Require every unordered pair, not just neighbours reported by a tool."""
    pairs = {}
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t" if path.suffix == ".tsv" else ",")
        raw = {"Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"}
        normalized = {"host_left", "host_right", "ani_percent", "af_left_percent", "af_right_percent"}
        fields = set(reader.fieldnames or [])
        require(fields >= raw or fields == normalized, "Unexpected skani edge schema")
        is_raw = fields >= raw
        for row in reader:
            names = []
            for column in (("Ref_file", "Query_file") if is_raw else ("host_left", "host_right")):
                if is_raw:
                    match = re.fullmatch(r"host_([0-9]{3})\.fna", Path(row[column]).name)
                    require(match is not None, "Unexpected skani genome name")
                    name = match[1]
                else:
                    name = row[column]
                require(name in host_ids, "Unexpected skani host ID")
                names.append(name)
            require(names[0] != names[1], "Unexpected diagonal in lower-triangle output")
            values = [float(row[k]) for k in (("ANI", "Align_fraction_ref", "Align_fraction_query")
                       if is_raw else ("ani_percent", "af_left_percent", "af_right_percent"))]
            require(all(math.isfinite(x) and 0 <= x <= 100 for x in values), "Invalid ANI/AF values")
            key = tuple(sorted(names))
            require(key not in pairs, "Duplicate unordered ANI pair")
            ani, af_left, af_right = values
            if names != list(key):
                af_left, af_right = af_right, af_left
            pairs[key] = {"host_left": key[0], "host_right": key[1], "ani_percent": ani,
                          "af_left_percent": af_left, "af_right_percent": af_right}
    require(set(pairs) == set(combinations(sorted(host_ids), 2)),
            "Missing ANI comparisons; skipped pairs cannot be treated as distant")
    return [pairs[key] for key in sorted(pairs)]


def make_groups(host_ids, pairs, digests, ani, minimum_af):
    components = DisjointSet(host_ids)
    for left, right in combinations(sorted(host_ids), 2):
        if digests[left] == digests[right]:
            components.union(left, right)
    edges = []
    for row in pairs:
        if row["ani_percent"] >= ani and min(row["af_left_percent"], row["af_right_percent"]) >= minimum_af:
            components.union(row["host_left"], row["host_right"])
            edges.append(row)
    return components.labels(), edges


def assign_splits(groups, seed, minimum_groups):
    ids = sorted(set(groups.values()),
                 key=lambda group: (sha256(f"{seed}:{group}".encode()), group))
    require(len(ids) >= minimum_groups, "Frozen design has insufficient genomic groups")
    size = len(ids) // 5
    test, validation = set(ids[:size]), set(ids[size:2*size])
    return {host: "test" if group in test else "validation" if group in validation else "train"
            for host, group in groups.items()}


def curate_grades(path, host_ids, library, expected):
    data = path.read_bytes()
    require(sha256(data) == expected["sha256"], "Source grade export changed")
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == expected["rows"], "Source grade row count changed")
    observed = {}
    for row in rows:
        target = primary_target(row["source_grade"])
        key = row["host_label"], row["phage_label"]
        require(key not in observed, "Duplicate source pair")
        observed[key] = (row, target)
    selected = {}
    for key in product(sorted(host_ids), library):
        require(key in observed, "Missing required source cell")
        row, target = observed[key]
        require(row["phage_role"] == "study" and row["source_sheet"] == "SF2",
                "Control or wrong source endpoint entered the library")
        selected[key] = {"source_grade": int(row["source_grade"]), "target": target}
    require(len(selected) == 6006, "Research pair scope changed")
    return selected


def write_csv(path, rows):
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare_genomes(sequence, qa, data_root, work):
    expected = {r["source_label"]: r for r in qa["qc"] if r["kind"] == "host_genome"}
    require(len(expected) == len(sequence["host_assemblies"]) == 143 and "127" not in expected,
            "Available-genome cohort changed")
    entries = []
    for source in sequence["host_assemblies"]:
        item = source["fasta"]
        data = (data_root / item["filename"]).read_bytes()
        require(len(data) == item["size_bytes"] and sha256(data) == item["sha256"],
                f"Host archive changed: {source['source_label']}")
        path = work / f"host_{source['source_label']}.fna"
        with path.open("xb") as stream:
            stream.write(gzip.decompress(data))
        digest = canonical_assembly_digest(path)
        require(digest == expected[source["source_label"]]["assembly_sha256"], "Host assembly identity changed")
        entries.append({"host_id": source["source_label"], "assembly_sha256": digest,
                        "sample_accession": source["sample_accession"],
                        "source_run": source["run_accession"], "path": str(path.resolve())})
    return sorted(entries, key=lambda r: r["host_id"])


def run_similarity(executable, roster, work, design):
    require(sha256(executable.read_bytes()) == design["grouping"]["binary_sha256"], "skani binary differs")
    version = subprocess.run([str(executable.resolve()), "--version"], check=True,
                             capture_output=True, text=True).stdout.strip()
    require(version == "skani 0.4.0", "skani version differs")
    listing = work / "genomes.txt"
    listing.write_text("\n".join(r["path"] for r in roster) + "\n")
    output = work / "skani_edges.tsv"
    command = [str(executable.resolve()), *design["grouping"]["command_options"],
               "-l", str(listing.resolve()), "-o", str(output.resolve())]
    def limits():
        cap = design["compute"]["max_ram_gb"] * 1024**3
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    started = time.monotonic()
    with (work / "skani.stdout").open("w") as stdout, (work / "skani.stderr").open("w") as stderr:
        subprocess.run(command, check=True, stdout=stdout, stderr=stderr, timeout=1200, preexec_fn=limits)
    return output, {"version": version, "binary_sha256": design["grouping"]["binary_sha256"],
                    "options": design["grouping"]["command_options"],
                    "wall_seconds": time.monotonic() - started,
                    "raw_edges_sha256": sha256(output.read_bytes())}


def freeze(design, sequence, qa, roster, pairs, grades, out, execution, design_commit):
    host_ids = {r["host_id"] for r in roster}
    digests = {r["host_id"]: r["assembly_sha256"] for r in roster}
    settings = design["grouping"]
    groups, edges = make_groups(host_ids, pairs, digests, settings["primary_ani_percent"],
                                settings["minimum_both_aligned_fraction_percent"])
    partitions = assign_splits(groups, design["split"]["seed"], design["split"]["minimum_total_groups"])
    require(all(partitions[r["host_left"]] == partitions[r["host_right"]] for r in edges),
            "Qualifying near-neighbour edge crosses partitions")
    hosts = [{k: row[k] for k in ("host_id", "assembly_sha256", "sample_accession", "source_run")} |
             {"host_group": groups[row["host_id"]], "partition": partitions[row["host_id"]]}
             for row in roster]
    sensitivity = {}
    for threshold in settings["sensitivity_ani_percent"]:
        labels, _ = make_groups(host_ids, pairs, digests, threshold,
                                settings["minimum_both_aligned_fraction_percent"])
        sensitivity[str(threshold)] = len(set(labels.values()))
    support = {}
    for part in ("train", "validation", "test"):
        ids = {h for h, p in partitions.items() if p == part}
        cells = [grades[h, phage] for h in ids for phage in design["dataset"]["library"]]
        support[part] = {"hosts": len(ids), "groups": len({groups[h] for h in ids}),
                         "source_pairs": len(cells), "grade3_targets": sum(x["target"] for x in cells),
                         "raw_grade_counts": {str(i): sum(x["source_grade"] == i for x in cells) for i in range(4)},
                         "all_primary_negative_hosts": sum(not any(grades[h, p]["target"]
                            for p in design["dataset"]["library"]) for h in ids)}
    out.mkdir(parents=True)
    write_csv(out / "host_pairs.csv", pairs)
    write_csv(out / "hosts.csv", hosts)
    bound_outputs = {name: {"sha256": sha256((out/name).read_bytes()), "size_bytes": (out/name).stat().st_size}
                     for name in ("host_pairs.csv", "hosts.csv")}
    protocol = {"schema_version": 1, "study_id": design["study_id"],
                "status": "frozen_research_protocol_no_susceptibility_fit_or_scoring",
                "design_path": "data/studies/kaneko_baseline.design.v1.json",
                "design_sha256": DESIGN_SHA256, "design_pre_run_commit": design_commit,
                "frozen_math_code_git_blobs": FROZEN_CODE_GIT_BLOBS,
                "future_fit_environment": FUTURE_FIT_ENVIRONMENT,
                "test_precision_limit": "Only seven test hosts in five genomic groups; report conditional research uncertainty and do not infer stable superiority.",
                "frozen_outputs": bound_outputs, "grouping": settings,
                "host_groups": len(set(groups.values())), "split_support": support,
                "library": design["dataset"]["library"], "endpoint_id": design["endpoint"]["id"],
                "host127_excluded": True, "model_fits": 0, "susceptibility_score_runs": 0,
                "independent_release_validation_established": False}
    (out/"protocol.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2)+"\n")
    receipt = {"schema_version": 1, "study_id": design["study_id"],
               "design_sha256": DESIGN_SHA256, "design_pre_run_commit": design_commit,
               "execution": execution, "observed_host_pairs": len(pairs), "required_host_pairs": 10153,
               "primary_edges": len(edges), "primary_groups": len(set(groups.values())),
               "largest_group_hosts": max(Counter(groups.values()).values()),
               "sensitivity_group_counts_only": sensitivity, "split_support": support,
               "cross_partition_primary_edges": 0, "source_pairs_scoped": len(grades),
               "candidate_count_every_host": 42, "unknowns_filled_as_zero": 0,
               "model_fits": 0, "susceptibility_score_runs": 0,
               "independent_release_validation_established": False}
    (out/"NUMERIC_QA.json").write_text(json.dumps(receipt, indent=2)+"\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-root", type=Path, required=True)
    parser.add_argument("--grades", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--skani", type=Path)
    parser.add_argument("--pairs", type=Path, help="Replay raw skani TSV or the committed normalized pair CSV")
    parser.add_argument("--design-pre-run-commit", default=DESIGN_PRE_RUN_COMMIT)
    args = parser.parse_args()
    require(not args.out_dir.exists() and not args.work_dir.exists(), "Output/work directories must be new")
    require(bool(args.skani) != bool(args.pairs), "Choose exactly one of --skani or --pairs")
    require(args.design_pre_run_commit == DESIGN_PRE_RUN_COMMIT, "Pre-run design commit differs")
    data = DESIGN.read_bytes()
    require(sha256(data) == DESIGN_SHA256, "Pre-run design changed")
    design = json.loads(data)
    for name, expected in FROZEN_CODE_GIT_BLOBS.items():
        code = (ROOT/name).read_bytes()
        digest = hashlib.sha1(b"blob " + str(len(code)).encode() + b"\0" + code).hexdigest()
        require(digest == expected, f"Frozen mathematical implementation changed: {name}")
    sequence = bound_json(design["inputs"]["sequence_intake"])
    qa = bound_json(design["inputs"]["sequence_qa"])
    bound_json(design["inputs"]["workbook_intake"])
    args.work_dir.mkdir(parents=True)
    roster = prepare_genomes(sequence, qa, args.sequence_root, args.work_dir)
    if args.skani:
        pair_path, execution = run_similarity(args.skani, roster, args.work_dir, design)
    else:
        pair_path, execution = args.pairs, {"mode": "offline_replay",
                                          "raw_edges_sha256": sha256(args.pairs.read_bytes())}
    pairs = read_pairs(pair_path, {r["host_id"] for r in roster})
    grades = curate_grades(args.grades, {r["host_id"] for r in roster},
                            design["dataset"]["library"], design["source_grades"])
    receipt = freeze(design, sequence, qa, roster, pairs, grades, args.out_dir,
                     execution, args.design_pre_run_commit)
    print(json.dumps({k: receipt[k] for k in ("observed_host_pairs", "primary_groups",
                      "largest_group_hosts", "cross_partition_primary_edges", "model_fits")}))


if __name__ == "__main__":
    main()
