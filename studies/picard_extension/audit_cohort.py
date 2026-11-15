#!/usr/bin/env python3
"""Audit source identity, missingness and genomic overlap before panel scoring."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "studies/picard_extension"
SPLITS = "data/studies/results/genophi_comparator.v1/scenario_splits.csv"
RESOURCES = "data/studies/genophi_host_baseline.resources.v1.json"
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def check(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(name):
    return str(name).upper().replace("-", "")


def unique(records, field):
    result = {}
    for row in records:
        k = key(row[field])
        check(k and k not in result, f"ambiguous alias: {row[field]}")
        result[k] = row
    return result


def read(path, delimiter=","):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(path, records):
    check(bool(records), "empty result without an output schema")
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


class Components:
    def __init__(self, names):
        self.parent = {h: h for h in names}

    def find(self, h):
        while self.parent[h] != h:
            self.parent[h] = self.parent[self.parent[h]]
            h = self.parent[h]
        return h

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        self.parent[max(a, b)] = min(a, b)


def group_records(names, distances, digests, threshold):
    """Transitive core-distance and exact-assembly components, including bridges."""
    names = sorted(names)
    check(0 < threshold and math.isfinite(threshold), "invalid core-distance threshold")
    check(set(digests) <= set(names), "assembly identity outside the grouping scope")
    result, edges = Components(names), []
    for a in names:
        for b in names:
            check((a, b) in distances and (b, a) in distances, "missing core-distance cell")
            value, reverse = distances[a, b], distances[b, a]
            check(math.isfinite(value) and value >= 0 and abs(value - reverse) <= 1e-10, "invalid asymmetric core distance")
            check(a != b or abs(value) <= 1e-10, "nonzero core-distance diagonal")
            if a < b and value < threshold:
                result.union(a, b)
                edges.append({"left": a, "right": b, "kind": "source_core_distance", "distance": value})
    duplicates = defaultdict(list)
    for h, digest in sorted(digests.items()):
        check(bool(digest), "empty assembly digest")
        duplicates[digest].append(h)
    for members in duplicates.values():
        for other in members[1:]:
            result.union(members[0], other)
            edges.append({"left": members[0], "right": other, "kind": "exact_canonical_assembly", "distance": ""})
    return result, edges


def fasta_profile(path, source):
    data = Path(path).read_bytes()
    check(len(data) == source["size"] and hashlib.md5(data).hexdigest() == source["computed_md5"],
          f"FASTA differs from Figshare: {source['name']}")
    contigs, pieces = [], None
    for line in data.decode("ascii").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if pieces is not None:
                contigs.append("".join(pieces).upper())
            pieces = []
        else:
            check(pieces is not None, "FASTA sequence precedes its header")
            pieces.append(line)
    check(pieces is not None, "FASTA is empty")
    contigs.append("".join(pieces).upper())
    check(all(c and not (set(c) - set("ACGTRYSWKMBDHVN")) for c in contigs), "invalid IUPAC DNA")
    per_contig = [hashlib.sha256(min(c, c.translate(COMPLEMENT)[::-1]).encode("ascii")).hexdigest() for c in contigs]
    canonical = hashlib.sha256("\n".join(sorted(per_contig)).encode("ascii")).hexdigest()
    lengths = sorted(map(len, contigs), reverse=True)
    total, cumulative, n50 = sum(lengths), 0, None
    for length in lengths:
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    ambiguous = sum(len(c.translate(str.maketrans("", "", "ACGT"))) for c in contigs)
    return {"raw_sha256": hashlib.sha256(data).hexdigest(), "assembly_sha256": canonical,
            "contigs": len(contigs), "total_bases": total, "n50": n50,
            "ambiguous_bases": ambiguous}


def load_sources(source_root, fetch=False):
    implementation = json.loads((STUDY / "implementation.json").read_text())
    for path, digest in implementation["files_sha256"].items():
        check(sha(ROOT / path) == digest, f"changed frozen implementation: {path}")
    plan = json.loads((STUDY / "protocol.json").read_text())
    lock = json.loads((STUDY / "input_manifest.json").read_text())
    for path, digest in lock["files_sha256"].items():
        check((ROOT / path).is_file() and sha(ROOT / path) == digest, f"changed frozen input: {path}")
    check(sha(STUDY / "genome_sources.json") == lock["genome_sources_sha256"], "genome source roster changed")
    source_root.mkdir(parents=True, exist_ok=True)
    names = {"source_metadata": "source_metadata.csv", "source_mlc": "source_mlc.csv", "core_distances": "core_distances.tsv"}
    for resource, filename in names.items():
        path, spec = source_root / filename, plan["resources"][resource]
        if not path.exists() and fetch:
            request = urllib.request.Request(spec["url"], headers={"User-Agent": "picard-extension-v1"})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read(spec["size_bytes"] + 1)
            check(len(data) == spec["size_bytes"] and hashlib.sha256(data).hexdigest() == spec["sha256"], "downloaded source differs")
            path.write_bytes(data)
        check(path.is_file() and path.stat().st_size == spec["size_bytes"] and sha(path) == spec["sha256"], f"missing/changed source: {filename}")
    return plan


def audit(source_root, genomes, out, fetch=False):
    plan = load_sources(source_root, fetch)
    matrix = read(source_root / "source_mlc.csv", ";")
    metadata = unique(read(source_root / "source_metadata.csv", ";"), "bacteria")
    source_hosts = unique(matrix, "bacteria")
    check(len(matrix) == 402 and len(metadata) == 403, "source roster differs")
    old = {}
    for row in read(ROOT / SPLITS):
        name = row["source_bacterium"]
        if name in old:
            check(all(old[name][f] == row[f] for f in ("bacterium", "assembly_sha256", "host_group")), "inconsistent frozen host identity")
        else:
            old[name] = row
    check(len(old) == 177 and set(old) <= {r["bacteria"] for r in matrix}, "frozen crosswalk differs")
    additional = sorted(r["bacteria"] for r in matrix if r["bacteria"] not in old)
    check(len(additional) == 225, "additional host roster differs")
    files = json.loads((STUDY / "genome_sources.json").read_text())["files"]
    files = unique([{**f, "host": f["name"].removesuffix(".fasta")} for f in files], "host")
    check(len(files) == 403, "Figshare roster differs")
    old_resources = {r["id"]: r for r in json.loads((ROOT / RESOURCES).read_text())["resources"] if r["kind"] == "host_genome"}
    for name, row in old.items():
        a, b = files[key(name)], old_resources[row["bacterium"]]
        check(a["id"] == b["source_file_id"] and a["size"] == b["size_bytes"] and a["computed_md5"] == b["md5"], "old Figshare identity changed")
    profiles = {}
    for index, name in enumerate(additional, 1):
        source = files[key(name)]
        path = genomes / f"{source['id']}.fasta"
        if not path.exists() and fetch:
            from fetch_genomes import retrieve
            genomes.mkdir(parents=True, exist_ok=True)
            retrieve(source, genomes)
        check(path.is_file(), f"missing genome: {name}; run fetch_genomes.py first")
        profiles[name] = fasta_profile(path, source)
        if index % 50 == 0:
            print(f"{index}/225 new assemblies profiled", flush=True)
    core_rows = read(source_root / "core_distances.tsv", "\t")
    core = unique(core_rows, "bacteria")
    check(len(core) == 404, "core-distance axis differs")
    columns = {key(c): c for c in core_rows[0] if c != "bacteria"}
    check(set(core) == set(columns) and set(source_hosts) <= set(core), "core axes do not cover source hosts")
    names = {k: source_hosts[k]["bacteria"] if k in source_hosts else core[k]["bacteria"] for k in core}
    distances = {(names[a], names[b]): float(core[a][columns[b]]) for a in core for b in core}
    threshold = plan["admission"]["core_distance_threshold"]
    digests = {h: r["assembly_sha256"] for h, r in old.items()}
    digests.update({h: r["assembly_sha256"] for h, r in profiles.items()})
    full, edges = group_records(names.values(), distances, digests, threshold)
    old_only, _ = group_records(old, distances, {h: digests[h] for h in old}, threshold)
    for a in old:
        for b in old:
            check((old_only.find(a) == old_only.find(b)) == (old[a]["host_group"] == old[b]["host_group"]), "original 147-group relation does not reproduce")
    groups = defaultdict(list)
    for h in sorted(names.values()):
        groups[full.find(h)].append(h)
    phages = sorted(p for p in matrix[0] if p not in ("bacteria", "T4LD", "T7_Portugal"))
    check(len(phages) == 94, "candidate library differs")
    ledger, unknowns, admitted = [], [], []
    for h in additional:
        meta, file = metadata[key(h)], files[key(h)]
        row = source_hosts[key(h)]
        missing = [p for p in phages if row[p].strip().lower() in ("", "nan", "na")]
        check(all(row[p].strip().lower() in ("", "nan", "na") or float(row[p]) in (0, 1, 2, 3, 4) for p in phages), "invalid source MLC")
        group = full.find(h)
        old_members = sorted(x for x in groups[group] if x in old)
        scope = meta["Clermont_Phylo"] in plan["admission"]["taxon_scope"]
        reasons = []
        if not scope:
            reasons.append("source_taxon_out_of_scope")
        if missing:
            reasons.append("incomplete_94_candidate_panel")
        if old_members:
            reasons.append("component_contains_existing_host")
        nearest = min(distances[h, x] for x in old)
        record = {"source_host": h, "gembase": meta["Gembase"], "source_phylogroup": meta["Clermont_Phylo"],
            "taxon_evidence": "source_metadata_only", "in_taxon_scope": scope,
            "figshare_file_id": file["id"], "figshare_filename": file["name"], "figshare_url": file["download_url"],
            "file_bytes": file["size"], "verified_md5": file["computed_md5"], **profiles[h],
            "observed_candidates": 94 - len(missing), "unknown_candidates": len(missing), "complete_panel": not missing,
            "component": group, "component_size_404": len(groups[group]),
            "existing_component_hosts": ";".join(old_members),
            "existing_frozen_groups": ";".join(sorted({old[x]["host_group"] for x in old_members})),
            "nearest_existing_core_distance": nearest,
            "nearest_existing_host_ids": ";".join(sorted(x for x in old if distances[h, x] == nearest)),
            "exact_duplicate_of_existing": ";".join(sorted(x for x in old if digests[x] == digests[h])),
            "eligible": not reasons, "blocking_reasons": ";".join(reasons)}
        ledger.append(record)
        if not reasons:
            admitted.append({"source_host": h, "component": group, "source_phylogroup": meta["Clermont_Phylo"], "assembly_sha256": digests[h]})
        unknowns.extend({"source_host": h, "phage": p, "status": "unknown", "source_value": row[p]} for p in missing)
    members_rows = [{"source_host": h, "component": g, "existing_host": h in old,
        "matrix_host": key(h) in source_hosts, "assembly_sha256": digests.get(h, "")}
        for g in sorted(groups) for h in groups[g]]
    write_csv(out / "host_admission.csv", ledger)
    write_csv(out / "admitted_hosts.csv", admitted)
    write_csv(out / "component_members.csv", members_rows)
    write_csv(out / "grouping_edges.csv", sorted(edges, key=lambda r: (r["kind"], r["left"], r["right"])))
    write_csv(out / "unknown_pairs.csv", unknowns)
    summary = {"status": "admission_completed_before_panel_scoring", "source_matrix_hosts": 402,
        "old_hosts": 177, "additional_hosts": 225, "new_fastas_verified": len(profiles),
        "new_fasta_bytes": sum(r["file_bytes"] for r in ledger), "old_canonical_digests_reused": 177,
        "old_groups_reproduced": len({old_only.find(h) for h in old}), "core_axis_hosts": 404,
        "distance_only_bridge_records": sorted(h for h in names.values() if key(h) not in source_hosts),
        "additional_source_phylogroups": dict(sorted(Counter(r["source_phylogroup"] for r in ledger).items())),
        "complete_additional_panels": sum(r["complete_panel"] for r in ledger),
        "unknown_pairs": len(unknowns), "in_taxon_scope": sum(r["in_taxon_scope"] for r in ledger),
        "overlaps_existing_component": sum(bool(r["existing_component_hosts"]) for r in ledger),
        "exact_duplicates_of_existing": sum(bool(r["exact_duplicate_of_existing"]) for r in ledger),
        "eligible_hosts": len(admitted), "eligible_groups": len({r["component"] for r in admitted}),
        "eligible_source_phylogroups": dict(sorted(Counter(r["source_phylogroup"] for r in admitted).items())),
        "blocking_reason_counts_nonexclusive": dict(sorted(Counter(reason for r in ledger for reason in r["blocking_reasons"].split(";") if reason).items())),
        "protocol_sha256": sha(STUDY / "protocol.json"), "input_manifest_sha256": sha(STUDY / "input_manifest.json"),
        "genome_sources_sha256": sha(STUDY / "genome_sources.json"), "code_sha256": sha(Path(__file__)),
        "output_sha256": {p.name: sha(p) for p in sorted(out.iterdir())}}
    write_json(out / "ADMISSION.json", summary)
    print(json.dumps({k: summary[k] for k in ("additional_hosts", "new_fastas_verified", "in_taxon_scope", "overlaps_existing_component", "eligible_hosts", "eligible_groups")}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, default=ROOT / "data/raw/picard_extension_sources")
    p.add_argument("--genomes", type=Path, default=ROOT / "data/raw/picard_extension_genomes")
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    check(not args.out.exists(), "output must be a new directory")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".admission-", dir=args.out.parent))
    try:
        audit(args.sources, args.genomes, stage, args.fetch)
        stage.rename(args.out)
    except BaseException:
        shutil.rmtree(stage)
        raise


if __name__ == "__main__":
    main()
