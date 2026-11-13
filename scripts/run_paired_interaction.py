#!/usr/bin/env python3
"""Train the frozen bilinear comparison and evaluate the inspected BASEL panel."""
from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import canonical_json, code_provenance, read_json_object, sha256_file, staged_directory
from phage_cocktail.basel_controls import evaluate_control, sequence_feature, validate_feature
from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock
from phage_cocktail.cohort import curate_mlc
from phage_cocktail.data import load_interactions
from phage_cocktail.host_genome import CANONICAL_WORDS
from phage_cocktail.host_study import verify_resources
from phage_cocktail.paired_baseline import grouped_ap_contrast, sequence_pair_weights, weighted_ap
from phage_cocktail.paired_interaction import (
    KINDS, classification_report, fit_entity_scaler, fit_interaction, predict_logits, transform_entities,
)
from phage_cocktail.registry import verify_registered_dataset
from scripts.run_paired_baseline import validate_manifest

PROTOCOL = "data/studies/paired_interaction.v1.json"
PROTOCOL_SHA256 = "bf871a3ac0e39e9e42660c88bc120499f99156b5f805bf92b6cadb1be8e62fc5"
FREEZE = "1b5c7cb314845067ef944231259cd3c36e659e4b"
BASEL = "data/studies/results/basel_ecor_fixed_controls.v1"
SPLIT = "data/studies/results/paired_split_freeze.v1/paired_split_i95_s42.csv.gz"


def check_files(hashes):
    for name, digest in hashes.items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or sha256_file(path) != digest:
            raise ValueError(f"changed frozen input: {name}")


def read_frame(path):
    return pd.read_csv(path, dtype={"bacterium": str, "phage": str, "entity_id": str}, float_precision="round_trip")


def write_gzip(frame, path):
    path.write_bytes(gzip.compress(frame.to_csv(index=False, float_format="%.17g").encode(), mtime=0))


def make_encoder(rows, hosts, phages, hd, pdig):
    return dict(host=fit_entity_scaler(hosts, rows.bacterium, hd),
                phage=fit_entity_scaler(phages, rows.phage, pdig))


def fit_one(rows, hosts, phages, kind, l2, optimizer):
    hids, pids = sorted(set(rows.bacterium)), sorted(set(rows.phage))
    return fit_interaction(hosts.loc[hids].to_numpy(), phages.loc[pids].to_numpy(),
                           pd.Index(hids).get_indexer(rows.bacterium),
                           pd.Index(pids).get_indexer(rows.phage), rows.infection.to_numpy(),
                           rows.weight.to_numpy(), kind=kind, l2=l2, optimizer=optimizer)


def row_logits(model, rows, hosts, phages):
    scores = predict_logits(model, hosts.to_numpy(), phages.to_numpy())
    return scores[pd.Index(hosts.index).get_indexer(rows.bacterium),
                  pd.Index(phages.index).get_indexer(rows.phage)]


def model_metadata(model, key):
    return {**{k: v for k, v in model.items() if k != "theta"}, "weights_key": key}


def run(data_root, out, *, fetch=False):
    out = Path(out).resolve()
    if out.exists():
        raise ValueError("output directory must be new")
    if sha256_file(ROOT / PROTOCOL) != PROTOCOL_SHA256:
        raise ValueError("interaction protocol differs from its pre-fit freeze")
    plan = read_json_object(ROOT / PROTOCOL)
    check_files(plan["files_sha256"])
    for name, version in plan["numerical_environment"].items():
        actual = sys.version.split()[0] if name == "python" else importlib.metadata.version(name)
        if actual != version:
            raise ValueError(f"numerical environment differs: {name} {actual}")
    base_protocol = read_json_object(ROOT / "data/studies/basel_ecor_evaluation_protocol.v2.json")
    check_files(base_protocol["files_sha256"])
    spec_hash = sha256(canonical_json(base_protocol["representation"]).encode()).hexdigest()
    resource_lock = read_json_object(ROOT / "data/studies/genophi_host_baseline.resources.v1.json")
    entries = [r for r in resource_lock["resources"] if r["kind"] in {"host_genome", "phage_genome"}
               or (r["kind"] == "table" and r["id"] in {"subset", "mlc"})]
    basel_lock = read_json_object(ROOT / "data/studies/basel_ecor_sequences.resources.v1.json")
    entries += [dict(r, path="basel_source.csv") for r in basel_lock["resources"] if r["kind"] == "phenotype_source"]
    paths = verify_resources({"resources": entries}, data_root, fetch=fetch)
    verify_registered_dataset(ROOT / "data/registry/genophi_ecoli_subset.v1.json", paths["table", "subset"])
    observed, excluded, _ = curate_mlc(load_interactions(paths["table", "subset"]),
                                      pd.read_csv(paths["table", "mlc"], sep=";", index_col=0))
    if len(observed) != 16597 or len(excluded) != 41:
        raise ValueError("source observations/unknown exclusions changed")
    host_audit = read_frame(ROOT / "data/studies/results/genome_host_retrieval.v1/genome_audit.csv").set_index("bacterium")
    phage_audit = read_frame(ROOT / "data/studies/results/paired_split_freeze.v1/phage_groups_i95.csv").set_index("phage")
    hd, pdig = host_audit.assembly_sha256.to_dict(), phage_audit.assembly_sha256.to_dict()
    features, receipts = [], []
    for kind, ids, digests in (("host_genome", sorted(set(observed.bacterium)), hd),
                               ("phage_genome", sorted(set(observed.phage)), pdig)):
        for name in ids:
            raw = paths[kind, name].read_bytes()
            vector, receipt = sequence_feature(raw, digests[name], spec_hash)
            validate_feature(vector, receipt, digests[name], spec_hash)
            identity = dict(entity_kind=kind, entity_id=name)
            features.append(dict(**identity, **dict(zip(CANONICAL_WORDS, vector))))
            receipts.append(dict(**identity, raw_sha256=sha256(raw).hexdigest(), **receipt))
    feature_frame = pd.DataFrame(features)
    hosts = feature_frame[feature_frame.entity_kind.eq("host_genome")].set_index("entity_id")[CANONICAL_WORDS].sort_index()
    phages = feature_frame[feature_frame.entity_kind.eq("phage_genome")].set_index("entity_id")[CANONICAL_WORDS].sort_index()
    if len(hosts) != 177 or len(phages) != 94:
        raise ValueError("sequence cohort differs")
    rows = validate_manifest(read_frame(ROOT / SPLIT), observed, hd, pdig)
    fit = rows[rows.partition.eq("fit")]
    val = rows[rows.partition.eq("validation") & rows.regime.eq("cold_phage")]
    if len(fit) != 6834 or len(val) != 1746:
        raise ValueError("frozen fit/selection support differs")
    encoders = dict(development=make_encoder(fit, hosts, phages, hd, pdig),
                    full_refit=make_encoder(rows, hosts, phages, hd, pdig))
    matrices = {stage: (transform_entities(hosts, enc["host"]), transform_entities(phages, enc["phage"]))
                for stage, enc in encoders.items()}
    weights = dict(observation_equal=np.ones(len(rows)), sequence_pair_equal=sequence_pair_weights(rows, hd, pdig))
    metadata, archive, trials, metrics, predictions, contrasts, final_models = {}, {}, [], [], [], [], {}
    with staged_directory(out) as stage, threadpool_limits(limits=1):
        write_gzip(feature_frame, stage / "training_features.csv.gz")
        pd.DataFrame(receipts).to_csv(stage / "training_feature_manifest.csv", index=False)
        write_gzip(observed, stage / "observed_pairs.csv.gz")
        excluded.to_csv(stage / "excluded_unknown_pairs.csv", index=False)
        for weighting, values in weights.items():
            current = rows.assign(weight=values)
            train = current[current.partition.eq("fit")]
            selection = current[current.partition.eq("validation") & current.regime.eq("cold_phage")]
            H, P = matrices["development"]
            selected, metadata[weighting] = {}, {}
            for kind in KINDS:
                best, best_ap = None, -np.inf
                for l2 in plan["training"]["l2_grid_in_tie_preference_order"]:
                    model = fit_one(train, H, P, kind, l2, plan["training"]["optimizer"])
                    ap = weighted_ap(selection.infection.to_numpy(), row_logits(model, selection, H, P), selection.weight.to_numpy())
                    trials.append(dict(weighting=weighting, model=kind, l2=l2, validation_AP=ap,
                                       **{k: v for k, v in model.items() if k not in {"kind", "l2", "theta"}}))
                    if ap > best_ap + plan["training"]["selection_tie_tolerance"]:
                        best, best_ap = model, ap
                selected[kind] = best
                key = "development_" + weighting + "_" + kind
                archive[key] = best["theta"]
                metadata[weighting][kind] = dict(development=model_metadata(best, key), validation_AP=best_ap)
                print(weighting, kind, "selected L2", best["l2"], "validation AP", best_ap, flush=True)
            for (regime, partition), held in current[current.partition.isin(["validation", "test"])].groupby(["regime", "partition"], sort=True):
                scores = {kind: row_logits(model, held, H, P) for kind, model in selected.items()}
                for kind, z in scores.items():
                    metrics.append(dict(weighting=weighting, model=kind, regime=regime, partition=partition,
                                        **classification_report(held, z, held.weight.to_numpy())))
                    predictions.append(held[["bacterium", "phage", "regime", "partition"]].assign(
                        weighting=weighting, model=kind, logit=z))
                if regime == "cold_phage" and partition == "test":
                    config = plan["development_evaluation"]["bootstrap"]
                    contrast = grouped_ap_contrast(held, scores["bilinear"], scores["additive"],
                                                   held.weight.to_numpy(), "cold_phage", resamples=config["resamples"],
                                                   seed=config["seed"], minimum_groups=config["minimum_groups"])
                    contrasts.append(dict(weighting=weighting, **contrast))
            Hfull, Pfull = matrices["full_refit"]
            for kind in KINDS:
                model = fit_one(current, Hfull, Pfull, kind, selected[kind]["l2"], plan["training"]["optimizer"])
                key = "full_refit_" + weighting + "_" + kind
                archive[key], final_models[weighting, kind] = model["theta"], model
                metadata[weighting][kind]["full_refit"] = model_metadata(model, key)
        # All hyperparameters and full-refit parameters now exist without BASEL labels.
        np.savez_compressed(stage / "weights.npz", **archive)
        artifact = dict(encoders=encoders, models=metadata, representation_sha256=spec_hash,
                        weights_sha256=sha256_file(stage / "weights.npz"), feature_columns=CANONICAL_WORDS,
                        score_semantics="uncalibrated_logit_for_original_observed_MLC_gt_0")
        (stage / "models.json").write_text(json.dumps(artifact, indent=2) + "\n")
        base_features = read_frame(ROOT / BASEL / "sequence_features.csv").set_index(["entity_kind", "source_build", "entity_id"], verify_integrity=True)
        base_manifest = read_frame(ROOT / BASEL / "feature_manifest.csv").set_index(["entity_kind", "source_build", "entity_id"], verify_integrity=True)
        for identity, row in base_manifest.iterrows():
            validate_feature(base_features.loc[identity].to_numpy(float), row.to_dict(), row.assembly_sha256, spec_hash)
        hids = sorted(base_features.loc[("host", "published_wgs_primary")].index)
        pids = sorted(base_features.loc[("phage", "basel_accession")].index)
        panels = {}
        for (weighting, kind), model in final_models.items():
            enc = encoders["full_refit"]
            Pbase = transform_entities(base_features.loc[("phage", "basel_accession")].loc[pids], enc["phage"])
            for build in ("published_wgs_primary", "training_assembly_sensitivity"):
                Hbase = transform_entities(base_features.loc[("host", build)].loc[hids], enc["host"])
                z = predict_logits(model, Hbase.to_numpy(), Pbase.to_numpy())
                panel = pd.MultiIndex.from_product([hids, pids], names=["bacterium", "phage"]).to_frame(index=False)
                panels[weighting, kind, build] = panel.assign(score=z.ravel())
        source = inspect_basel_ecor(paths["phenotype_source", "BASEL_ECOR"].read_bytes(), load_lock(ROOT / "data/studies/basel_ecor_intake.v1.json"))
        if hids != source["host_ids"] or pids != source["phage_ids"]:
            raise ValueError("BASEL scores do not cover the frozen panel")
        truth = pd.read_csv(paths["phenotype_source", "BASEL_ECOR"]).rename(columns={"strain": "bacterium", "interaction": "infection"})
        base_summary, per_host, base_predictions = {}, [], []
        near = read_json_object(ROOT / BASEL / "report.json")["training_neighbours_i90_retained"]
        for (weighting, kind, build), panel in panels.items():
            summary, host_rows = evaluate_control(truth, panel, hids, pids, near_training=near)
            method = weighting + ":" + kind + ":" + build
            base_summary[method] = summary
            per_host.append(host_rows.assign(weighting=weighting, model=kind, host_build=build))
            base_predictions.append(panel.assign(weighting=weighting, model=kind, host_build=build))
        for weighting in weights:
            for build in ("published_wgs_primary", "training_assembly_sensitivity"):
                a = base_summary[weighting + ":additive:" + build]["macro"]
                b = base_summary[weighting + ":bilinear:" + build]
                b["hit_minus_additive_conservative_outer_bounds"] = {
                    k: dict(lower=v["hit_lower"] - a[k]["hit_upper"], upper=v["hit_upper"] - a[k]["hit_lower"])
                    for k, v in b["macro"].items()}
        pd.DataFrame(trials).to_csv(stage / "validation_trials.csv", index=False, float_format="%.17g")
        pd.DataFrame(metrics).to_csv(stage / "development_metrics.csv", index=False, float_format="%.17g")
        write_gzip(pd.concat(predictions, ignore_index=True), stage / "development_predictions.csv.gz")
        write_gzip(pd.concat(base_predictions, ignore_index=True), stage / "basel_scores.csv.gz")
        pd.concat(per_host, ignore_index=True).to_csv(stage / "basel_per_host_bounds.csv", index=False, float_format="%.17g")
        report = dict(study_id=plan["study_id"], status="completed_inspected_development_and_published_reference",
                      observed_training_pairs=len(observed), excluded_unknown_training_pairs=len(excluded),
                      training_feature_vectors=len(features), feature_dimensions=136, development_contrasts=contrasts,
                      basel_source_counts=source["source_summary"], basel=base_summary,
                      frozen_fixed_comparators=read_json_object(ROOT / BASEL / "report.json")["methods"],
                      no_independent_release_validation=True,
                      provenance=dict(**code_provenance(ROOT, Path(__file__)), protocol_sha256=PROTOCOL_SHA256,
                                      protocol_freeze_commit=FREEZE, numerical_environment=plan["numerical_environment"]),
                      artifacts_sha256={p.name: sha256_file(p) for p in sorted(stage.iterdir())})
        (stage / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    try:
        report = run(args.data_root, args.out, fetch=args.fetch)
    except Exception as error:
        failure = args.out.with_name(args.out.name + ".failure.json")
        failure.parent.mkdir(parents=True, exist_ok=True)
        with failure.open("x") as stream:
            json.dump(dict(study_id="sequence-paired-interaction:v1", status="failed_no_admitted_result",
                           error_type=type(error).__name__, reason=str(error),
                           protocol_sha256=PROTOCOL_SHA256, **code_provenance(ROOT, Path(__file__))), stream, indent=2)
            stream.write("\n")
        raise
    for name, values in report["basel"].items():
        print(name, "Hit@2", values["macro"]["2"]["hit_lower"], values["macro"]["2"]["hit_upper"])


if __name__ == "__main__":
    main()
