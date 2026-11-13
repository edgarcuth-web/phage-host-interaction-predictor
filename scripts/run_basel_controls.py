#!/usr/bin/env python3
"""Run the separately frozen BASEL v2 composition and random controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd
import sklearn

from phage_cocktail.artifacts import canonical_json, code_provenance, read_json_object, sha256_file, staged_directory
from phage_cocktail.basel_controls import cosine_panel, evaluate_control, sequence_feature, validate_feature
from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock
from phage_cocktail.host_genome import CANONICAL_WORDS
from phage_cocktail.host_study import verify_resources

PROTOCOL = "data/studies/basel_ecor_evaluation_protocol.v2.json"
FREEZE = "c0c01dae8d558643ec279e8cc45387b0d8780ecf"
PROTOCOL_SHA256 = "6e2f370fd24f7f3997285a8ca59ecf123424202e748dfae819bae7b6aac4a182"
AUDIT = "data/studies/results/basel_ecor_sequence_audit.v1"


def run(data_root, out, *, fetch=False):
    out = Path(out).resolve()
    if out.exists():
        raise ValueError("output directory must be new")
    # A digest check also works after a squash merge or shallow checkout.
    if sha256_file(ROOT / PROTOCOL) != PROTOCOL_SHA256:
        raise ValueError("protocol differs from its pre-scoring freeze commit")
    protocol = read_json_object(ROOT / PROTOCOL)
    for relative, digest in protocol["files_sha256"].items():
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT) or sha256_file(path) != digest:
            raise ValueError(f"changed frozen evidence/implementation: {relative}")
    representation = protocol["representation"]
    if representation["ordered_features"] != CANONICAL_WORDS or representation["dimensions"] != len(CANONICAL_WORDS):
        raise ValueError("representation columns differ from the frozen encoding")
    spec_hash = sha256(canonical_json(representation).encode()).hexdigest()
    lock = read_json_object(ROOT / "data/studies/basel_ecor_sequences.resources.v1.json")
    training = read_json_object(ROOT / "data/studies/genophi_host_baseline.resources.v1.json")
    selected = [r for r in lock["resources"] if r["kind"] in
                {"basel_phage_genome", "published_host_genome", "phenotype_source"}]
    selected += [dict(r, path="training/" + r["path"]) for r in training["resources"]
                 if r["kind"] == "host_genome" and r["id"] in lock["host_ids"]]
    paths = verify_resources({"resources": selected}, data_root, fetch=fetch)
    source_path = paths[("phenotype_source", "BASEL_ECOR")]
    source = inspect_basel_ecor(source_path.read_bytes(), load_lock(ROOT / "data/studies/basel_ecor_intake.v1.json"))
    phages = pd.read_csv(ROOT / AUDIT / "basel_accessions.csv").set_index("phage")
    hosts = pd.read_csv(ROOT / AUDIT / "host_assemblies.csv").set_index("host")
    if set(phages.index) != set(source["phage_ids"]) or set(hosts.index) != set(source["host_ids"]):
        raise ValueError("sequence evidence and fixed source roster differ")
    audit = read_json_object(ROOT / AUDIT / "report.json")
    near = next(s["candidates_in_training_components"] for s in audit["phage_overlap"] if s["identity"] == .90)
    vectors, receipts, features = {}, [], []
    for entry in selected:
        kind, name = entry["kind"], entry["id"]
        if kind == "phenotype_source":
            continue
        raw = paths[(kind, name)].read_bytes()
        if kind == "basel_phage_genome":
            role, build, payload = "phage", "basel_accession", raw
            expected = phages.loc[name, "assembly_sha256"]
            if entry["accession_version"] != phages.loc[name, "accession_version"]:
                raise ValueError("phage accession.version differs from the identity receipt")
        elif kind == "published_host_genome":
            role, build, payload = "host", "published_wgs_primary", gzip.decompress(raw)
            expected = hosts.loc[name, "published_assembly_sha256"]
        else:
            role, build, payload = "host", "training_assembly_sensitivity", raw
            expected = hosts.loc[name, "training_assembly_sha256"]
        vector, receipt = sequence_feature(payload, expected, spec_hash)
        validate_feature(vector, receipt, expected, spec_hash)
        identity = dict(entity_id=name, entity_kind=role, source_build=build)
        receipts.append(dict(**identity, source_file_sha256=sha256(raw).hexdigest(),
                             fasta_payload_sha256=sha256(payload).hexdigest(), **receipt))
        features.append(dict(**identity, **dict(zip(CANONICAL_WORDS, vector))))
        vectors[role, build, name] = vector
    phage_vectors = {p: vectors["phage", "basel_accession", p] for p in source["phage_ids"]}
    panels = {}
    for build in ("published_wgs_primary", "training_assembly_sensitivity"):
        host_vectors = {h: vectors["host", build, h] for h in source["host_ids"]}
        panels["cosine_" + build] = cosine_panel(host_vectors, phage_vectors)
    panels["uniform_random"] = panels["cosine_published_wgs_primary"].assign(score=0.)
    # Phenotypes enter only evaluation; the complete score matrices already exist.
    observed = pd.read_csv(source_path).rename(columns={"strain": "bacterium", "interaction": "infection"})
    summaries, per_host, rankings = {}, [], []
    ks = (protocol["target"]["primary_K"], *protocol["target"]["secondary_K"])
    for method, scores in panels.items():
        summary, frame = evaluate_control(observed, scores, source["host_ids"], source["phage_ids"],
                                          ks=ks, near_training=near)
        summaries[method] = summary
        per_host.append(frame.assign(method=method))
        for host, group in scores.groupby("bacterium", sort=True):
            names = group.sort_values(["score", "phage"], ascending=[False, True]).phage.head(2).tolist()
            rankings.append(dict(method=method, bacterium=host, first=names[0], second=names[1]))
    random = summaries["uniform_random"]["macro"]
    for method, summary in summaries.items():
        if method == "uniform_random":
            continue
        summary["hit_minus_random_conservative_outer_bounds"] = {
            k: dict(lower=v["hit_lower"] - random[k]["hit_upper"],
                    upper=v["hit_upper"] - random[k]["hit_lower"])
            for k, v in summary["macro"].items()
        }
    sets = {(r["method"], r["bacterium"]): {r["first"], r["second"]} for r in rankings}
    changed = [h for h in source["host_ids"] if sets["cosine_published_wgs_primary", h]
               != sets["cosine_training_assembly_sensitivity", h]]
    with staged_directory(out) as stage:
        pd.DataFrame(receipts).to_csv(stage / "feature_manifest.csv", index=False)
        pd.DataFrame(features).to_csv(stage / "sequence_features.csv", index=False, float_format="%.17g")
        for method, scores in panels.items():
            scores.to_csv(stage / (method + "_scores.csv"), index=False, float_format="%.17g")
        pd.concat(per_host, ignore_index=True).to_csv(stage / "per_host_bounds.csv", index=False, float_format="%.17g")
        pd.DataFrame(rankings).to_csv(stage / "deterministic_top2.csv", index=False)
        report = dict(
            schema_version=1, protocol_id=protocol["protocol_id"],
            status="two_frozen_controls_complete_inspected_reference_only", learned_model_trained=False,
            primary_host_build="published_wgs_primary", source_counts=source["source_summary"],
            feature_vectors=len(receipts), feature_dimensions=len(CANONICAL_WORDS),
            representation_sha256=spec_hash, methods=summaries,
            assembly_sensitivity=dict(changed_deterministic_top2_hosts=changed,
                                      changed_host_count=len(changed), total_hosts=len(hosts)),
            source_decisions=protocol["source_adjudication"],
            training_neighbours_i90_retained=near,
            independent_release_validation_established=False,
            limitations=["Fixed composition similarity is not a learned interaction model or infection probability.",
                         "All bounds describe missing-label assignments, not sampling confidence.",
                         "Cosine-minus-random intervals are conservative outer bounds, not sharp paired bounds.",
                         "Host identifiers are known; distinct assemblies do not establish cold hosts.",
                         "No claim to reconstruct supplementary-only phenotypes or explain their absence.",
                         "This inspected published reference cannot alone validate deployment."],
            provenance=dict(**code_provenance(ROOT, Path(__file__)), protocol_freeze_commit=FREEZE,
                            protocol_sha256=sha256_file(ROOT / PROTOCOL),
                            python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                            sklearn=sklearn.__version__),
            artifacts_sha256={p.name: sha256_file(p) for p in sorted(stage.iterdir())},
        )
        (stage / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    report = run(args.data_root, args.out, fetch=args.fetch)
    for method, result in report["methods"].items():
        values = result["macro"]["2"]
        print(method, "macro Hit@2 identification bounds:", values["hit_lower"], values["hit_upper"])
    print("Changed top-2 hosts:", report["assembly_sensitivity"]["changed_host_count"], "/ 25")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
