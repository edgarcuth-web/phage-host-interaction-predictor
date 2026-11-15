#!/usr/bin/env python3
"""Independently reconstruct assembly identities, graph components and admission."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "studies/picard_extension"


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, sep=","):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=sep))


def norm(value):
    return value.replace("-", "").upper()


def canonical_profile(path):
    raw = Path(path).read_bytes()
    # Split on FASTA headers instead of using the runner's line-state parser.
    require(raw.startswith(b">"), "invalid FASTA header")
    contigs = []
    for block in raw.decode("ascii").split(">"):
        if not block:
            continue
        header, *lines = block.splitlines()
        seq = "".join(line.strip() for line in lines if line.strip()).upper()
        require(header and seq and not set(seq) - set("ACGTRYSWKMBDHVN"), "invalid DNA record")
        contigs.append(seq)
    complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    canonical_contigs = sorted(hashlib.sha256(min(s, s[::-1].translate(complement)).encode()).hexdigest() for s in contigs)
    assembly = hashlib.sha256("\n".join(canonical_contigs).encode()).hexdigest()
    lengths = sorted((len(s) for s in contigs), reverse=True)
    halfway = sum(lengths) / 2
    accumulated = 0
    for n50 in lengths:
        accumulated += n50
        if accumulated >= halfway:
            break
    return {"verified_md5": hashlib.md5(raw).hexdigest(), "raw_sha256": hashlib.sha256(raw).hexdigest(),
        "assembly_sha256": assembly, "contigs": len(contigs), "total_bases": sum(lengths), "n50": n50,
        "ambiguous_bases": sum(len(s) - sum(s.count(b) for b in "ACGT") for s in contigs), "file_bytes": len(raw)}


def components(nodes, adjacency):
    """Connected components by graph traversal, with no runner union-find."""
    result = {}
    for start in sorted(nodes):
        if start in result:
            continue
        pending, reached = [start], {start}
        while pending:
            h = pending.pop()
            for other in adjacency[h] - reached:
                reached.add(other)
                pending.append(other)
        label = min(reached)
        result.update({h: label for h in reached})
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, default=ROOT / "data/raw/picard_extension_sources")
    p.add_argument("--genomes", type=Path, default=ROOT / "data/raw/picard_extension_genomes")
    p.add_argument("--admission", type=Path, required=True)
    args = p.parse_args()
    plan = json.loads((STUDY / "protocol.json").read_text())
    lock = json.loads((STUDY / "input_manifest.json").read_text())
    implementation = json.loads((STUDY / "implementation.json").read_text())
    for path, digest in {**lock["files_sha256"], **implementation["files_sha256"]}.items():
        require(sha(ROOT / path) == digest, f"frozen input/code changed: {path}")
    for name, spec in plan["resources"].items():
        if "sha256" in spec:
            filename = name + (".tsv" if name == "core_distances" else ".csv")
            require(sha(args.sources / filename) == spec["sha256"], "source bytes changed")
    require(sha(STUDY / "genome_sources.json") == lock["genome_sources_sha256"], "source roster changed")
    receipt = json.loads((args.admission / "ADMISSION.json").read_text())
    for filename, digest in receipt["output_sha256"].items():
        require(sha(args.admission / filename) == digest, "admission output changed")
    matrix = read(args.sources / "source_mlc.csv", ";")
    metadata = {norm(r["bacteria"]): r for r in read(args.sources / "source_metadata.csv", ";")}
    source = {r["bacteria"]: r for r in matrix}
    splits = read(ROOT / "data/studies/results/genophi_comparator.v1/scenario_splits.csv")
    old = {r["source_bacterium"]: r for r in splits}
    additional = sorted(set(source) - set(old))
    require(len(source) == 402 and len(old) == 177 and len(additional) == 225, "roster differs")
    ledger = read(args.admission / "host_admission.csv")
    require(len(ledger) == 225 and [r["source_host"] for r in ledger] == additional, "ledger support differs")
    files = {norm(f["name"].removesuffix(".fasta")): f for f in json.loads((STUDY / "genome_sources.json").read_text())["files"]}
    digests = {h: r["assembly_sha256"] for h, r in old.items()}
    for r in ledger:
        host = r["source_host"]
        file = files[norm(host)]
        profile = canonical_profile(args.genomes / f"{file['id']}.fasta")
        require(profile["verified_md5"] == file["computed_md5"] and profile["file_bytes"] == file["size"], "Figshare bytes differ")
        for field, value in profile.items():
            require(str(value) == r[field], f"genome profile differs: {host}/{field}")
        require(str(file["id"]) == r["figshare_file_id"] and file["name"] == r["figshare_filename"] and file["download_url"] == r["figshare_url"], "genome route differs")
        digests[host] = profile["assembly_sha256"]
    raw_core = read(args.sources / "core_distances.tsv", "\t")
    core_lookup = {norm(r["bacteria"]): r for r in raw_core}
    columns = {norm(c): c for c in raw_core[0] if c != "bacteria"}
    canonical_names = {norm(h): h for h in source}
    names = {k: canonical_names.get(k, r["bacteria"]) for k, r in core_lookup.items()}
    require(len(names) == 404 and set(names) == set(columns), "core axes differ")
    adjacency, edge_set = defaultdict(set), set()
    distances = {}
    for a, row in core_lookup.items():
        for b in columns:
            value = float(row[columns[b]])
            require(math.isfinite(value) and value >= 0 and abs(value - float(core_lookup[b][columns[a]])) <= 1e-10, "invalid core-distance matrix")
            require(a != b or abs(value) <= 1e-10, "nonzero core diagonal")
            distances[names[a], names[b]] = value
            if a != b and value < plan["admission"]["core_distance_threshold"]:
                adjacency[names[a]].add(names[b])
                left, right = sorted((names[a], names[b]))
                edge_set.add((left, right, "source_core_distance"))
    same = defaultdict(list)
    for h in sorted(digests):
        same[digests[h]].append(h)
    for members in same.values():
        for other in members[1:]:
            adjacency[members[0]].add(other)
            adjacency[other].add(members[0])
            edge_set.add((members[0], other, "exact_canonical_assembly"))
    labels = components(names.values(), adjacency)
    old_labels = components(old, {h: adjacency[h] & set(old) for h in old})
    require(len(set(old_labels.values())) == 147, "old grouping does not reproduce")
    for a in old:
        for b in old:
            require((old_labels[a] == old_labels[b]) == (old[a]["host_group"] == old[b]["host_group"]), "old group relation differs")
    members = defaultdict(set)
    for h, group in labels.items():
        members[group].add(h)
    reported = read(args.admission / "component_members.csv")
    require(len(reported) == len({r["source_host"] for r in reported}) == 404, "component ledger support differs")
    for r in reported:
        require(labels[r["source_host"]] == r["component"], "component membership differs")
        require(r["assembly_sha256"] == digests.get(r["source_host"], ""), "component assembly differs")
    retained_edges = read(args.admission / "grouping_edges.csv")
    require(len(retained_edges) == len(edge_set) and {(r["left"], r["right"], r["kind"]) for r in retained_edges} == edge_set, "group edge support differs")
    phages = [p for p in matrix[0] if p not in ("bacteria", "T4LD", "T7_Portugal")]
    admitted, missing_expected = [], []
    for r in ledger:
        host = r["source_host"]
        meta = metadata[norm(host)]
        require(r["gembase"] == meta["Gembase"] and r["source_phylogroup"] == meta["Clermont_Phylo"], "metadata identity differs")
        unknowns = sorted(p for p in phages if source[host][p].strip().lower() in ("", "nan", "na"))
        old_members = sorted(members[labels[host]] & set(old))
        scoped = meta["Clermont_Phylo"] in plan["admission"]["taxon_scope"]
        reasons = ([] if scoped else ["source_taxon_out_of_scope"]) + ([] if not unknowns else ["incomplete_94_candidate_panel"]) + ([] if not old_members else ["component_contains_existing_host"])
        require(r["blocking_reasons"] == ";".join(reasons) and (r["eligible"] == "True") == (not reasons), "admission rule differs")
        require(r["existing_component_hosts"] == ";".join(old_members), "existing-host overlap differs")
        require(int(r["unknown_candidates"]) == len(unknowns) and int(r["observed_candidates"]) == 94 - len(unknowns), "unknown denominator differs")
        require(r["component"] == labels[host], "admission group differs")
        require(abs(float(r["nearest_existing_core_distance"]) - min(distances[host, x] for x in old)) < 1e-15, "nearest existing distance differs")
        if not reasons:
            admitted.append({"source_host": host, "component": labels[host], "source_phylogroup": meta["Clermont_Phylo"], "assembly_sha256": digests[host]})
        missing_expected.extend((host, p) for p in unknowns)
    require(read(args.admission / "admitted_hosts.csv") == admitted, "frozen admitted roster differs")
    unknown_rows = read(args.admission / "unknown_pairs.csv")
    require(len(unknown_rows) == len(missing_expected) and {(r["source_host"], r["phage"]) for r in unknown_rows} == set(missing_expected), "unknown-pair support differs")
    require(int(receipt["eligible_hosts"]) == len(admitted) and int(receipt["eligible_groups"]) == len({r["component"] for r in admitted}), "receipt denominator differs")
    result = {"status": "verified", "imports_runner": False, "new_assemblies_rehashed": 225,
        "old_canonical_digests_reused": 177, "graph_nodes_reconstructed": 404, "graph_edges_checked": len(edge_set),
        "admission_records_checked": 225, "unknown_pairs_checked": len(unknown_rows),
        "eligible_hosts": len(admitted), "eligible_groups": len({r["component"] for r in admitted}),
        "verifier_sha256": sha(Path(__file__)), "admission_receipt_sha256": sha(args.admission / "ADMISSION.json")}
    (args.admission / "INDEPENDENT_QA.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
