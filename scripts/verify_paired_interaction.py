#!/usr/bin/env python3
"""Verify serialized interaction weights, independent matrix scores and metrics."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import read_json_object, sha256_file
from phage_cocktail.basel_controls import validate_feature
from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock
from phage_cocktail.cohort import curate_mlc
from phage_cocktail.data import load_interactions
from phage_cocktail.paired_baseline import grouped_ap_contrast
from scripts.verify_basel_controls import independent_metrics


def close(a, b, atol=1e-10):
    if not np.isfinite(np.asarray(a)).all() or not np.allclose(a, b, rtol=0, atol=atol):
        raise ValueError("independent numerical check failed")


def frame(path):
    return pd.read_csv(path, dtype={"bacterium": str, "phage": str, "entity_id": str}, float_precision="round_trip")


def unpack_independent(theta, kind, dh, dp):
    position, a = 1, np.zeros(dh)
    if kind != "phage_only":
        a, position = theta[position:position + dh], position + dh
    b, position = theta[position:position + dp], position + dp
    W = theta[position:].reshape(dh, dp) if kind == "bilinear" else np.zeros((dh, dp))
    if len(theta) != 1 + dp + (dh if kind != "phage_only" else 0) + (dh*dp if kind == "bilinear" else 0):
        raise ValueError("weight archive has wrong parameter dimensions")
    return theta[0], a, b, W


def independent_scores(theta, kind, H, P):
    c, a, b, W = unpack_independent(theta, kind, H.shape[1], P.shape[1])
    return c + (H @ a)[:, None] + (P @ b)[None, :] + H @ W @ P.T


def verify(bundle, data_root):
    bundle, data_root = Path(bundle), Path(data_root)
    report = read_json_object(bundle / "report.json")
    protocol_path = ROOT / "data/studies/paired_interaction.v1.json"
    plan = read_json_object(protocol_path)
    if sha256_file(protocol_path) != report["provenance"]["protocol_sha256"]:
        raise ValueError("protocol hash differs")
    for name, digest in plan["files_sha256"].items():
        if sha256_file(ROOT / name) != digest: raise ValueError(f"frozen input changed: {name}")
    for name, digest in report["artifacts_sha256"].items():
        path = (bundle / name).resolve()
        if not path.is_relative_to(bundle.resolve()) or sha256_file(path) != digest:
            raise ValueError(f"artifact changed: {name}")
    for name, digest in report["provenance"]["code_sha256"].items():
        if sha256_file(ROOT / name) != digest: raise ValueError(f"execution code changed: {name}")
    lock = read_json_object(ROOT / "data/studies/genophi_host_baseline.resources.v1.json")
    resources = {(r["kind"], r["id"]): r for r in lock["resources"]}
    for name in ("subset", "mlc"):
        entry = resources["table", name]
        if sha256_file(data_root / entry["path"]) != entry["sha256"]:
            raise ValueError("primary phenotype source changed")
    truth, excluded, _ = curate_mlc(load_interactions(data_root / resources["table", "subset"]["path"]),
                                    pd.read_csv(data_root / resources["table", "mlc"]["path"], sep=";", index_col=0))
    saved = frame(bundle / "observed_pairs.csv.gz")
    if len(truth) != 16597 or len(excluded) != 41 or len(saved) != 16597:
        raise ValueError("training label support differs")
    joined = truth.merge(saved, on=["bacterium", "phage"], validate="one_to_one", suffixes=("_source", "_saved"))
    if len(joined) != 16597 or not joined.infection_source.equals(joined.infection_saved):
        raise ValueError("saved observations differ from reconciled primary MLC")
    split = frame(ROOT / "data/studies/results/paired_split_freeze.v1/paired_split_i95_s42.csv.gz")
    rows = split.merge(truth, on=["bacterium", "phage"], validate="one_to_one")
    features = frame(bundle / "training_features.csv.gz").set_index(["entity_kind", "entity_id"], verify_integrity=True)
    manifest = frame(bundle / "training_feature_manifest.csv").set_index(["entity_kind", "entity_id"], verify_integrity=True)
    model_info = read_json_object(bundle / "models.json")
    columns, spec = model_info["feature_columns"], model_info["representation_sha256"]
    ha = frame(ROOT / "data/studies/results/genome_host_retrieval.v1/genome_audit.csv").set_index("bacterium")
    pa = frame(ROOT / "data/studies/results/paired_split_freeze.v1/phage_groups_i95.csv").set_index("phage")
    if len(features) != 271 or set(features.index) != set(manifest.index) or features.columns.tolist() != columns:
        raise ValueError("training feature panel differs")
    for identity, entry in manifest.iterrows():
        kind, name = identity
        expected = ha.loc[name, "assembly_sha256"] if kind == "host_genome" else pa.loc[name, "assembly_sha256"]
        validate_feature(features.loc[identity].to_numpy(float), entry.to_dict(), expected, spec)
        if entry.raw_sha256 != resources[kind, name]["sha256"]:
            raise ValueError("feature input bytes differ from source lock")
    raw = dict(host=features.loc["host_genome"].sort_index(), phage=features.loc["phage_genome"].sort_index())
    matrices = {}
    for stage, enc in model_info["encoders"].items():
        admitted = rows[rows.partition.eq("fit")] if stage == "development" else rows
        matrices[stage] = {}
        for role, idcolumn, audit in (("host", "bacterium", ha), ("phage", "phage", pa)):
            names = sorted(set(admitted[idcolumn]))
            unique = {}
            for name in names: unique.setdefault(audit.loc[name, "assembly_sha256"], name)
            values = raw[role].loc[list(unique.values())].to_numpy(float)
            means, scales = values.mean(axis=0), values.std(axis=0)
            scales[scales == 0] = 1.
            close(enc[role]["means"], means, 1e-15); close(enc[role]["scales"], scales, 1e-15)
            if enc[role]["unique_fit_assemblies"] != len(unique): raise ValueError("fit-only entity scope differs")
            matrices[stage][role] = pd.DataFrame((raw[role].to_numpy(float)-means)/scales,
                                                 index=raw[role].index, columns=columns)
    key_counts = Counter((ha.loc[r.bacterium, "assembly_sha256"], pa.loc[r.phage, "assembly_sha256"]) for r in rows.itertuples())
    rows["sequence_weight"] = [1/key_counts[ha.loc[r.bacterium, "assembly_sha256"], pa.loc[r.phage, "assembly_sha256"]] for r in rows.itertuples()]
    archive = np.load(bundle / "weights.npz", allow_pickle=False)
    if sha256_file(bundle / "weights.npz") != model_info["weights_sha256"]:
        raise ValueError("model archive identity differs")
    predictions = frame(bundle / "development_predictions.csv.gz")
    metric_table, trials = frame(bundle / "development_metrics.csv"), frame(bundle / "validation_trials.csv")
    if len(predictions) != 53334 or len(trials) != 18 or len(metric_table) != 48:
        raise ValueError("development trial or score support differs")
    model_count, development_checks = 0, 0
    for weighting, models in model_info["models"].items():
        weights = np.ones(len(rows)) if weighting == "observation_equal" else rows.sequence_weight.to_numpy()
        current = rows.assign(weight=weights)
        for kind, metadata in models.items():
            t = trials[trials.weighting.eq(weighting) & trials.model.eq(kind)].sort_values("l2", ascending=False)
            best, best_ap = None, -np.inf
            for trial in t.itertuples():
                if trial.validation_AP > best_ap + plan["training"]["selection_tie_tolerance"]:
                    best, best_ap = trial.l2, trial.validation_AP
            close(metadata["development"]["l2"], best)
            close(metadata["validation_AP"], best_ap)
            for stage in ("development", "full_refit"):
                saved_model = metadata[stage]
                theta = archive[saved_model["weights_key"]]
                H, P = matrices[stage]["host"], matrices[stage]["phage"]
                Z = independent_scores(theta, kind, H.to_numpy(), P.to_numpy())
                fit = current[current.partition.eq("fit")] if stage == "development" else current
                hi, pi = H.index.get_indexer(fit.bacterium), P.index.get_indexer(fit.phage)
                z, y, w = Z[hi, pi], fit.infection.to_numpy(), fit.weight.to_numpy()
                residual = w*(expit(z)-y)/w.sum()
                # Observation-indexed derivatives are independent of the fitted cell-scatter implementation.
                h, p = H.to_numpy()[hi], P.to_numpy()[pi]
                blocks = [np.array([residual.sum()])]
                if kind != "phage_only": blocks.append(h.T @ residual)
                blocks.append(p.T @ residual)
                if kind == "bilinear": blocks.append((h.T @ (residual[:, None]*p)).ravel())
                gradient = np.concatenate(blocks); gradient[1:] += saved_model["l2"]*theta[1:]
                objective = np.average(np.logaddexp(0.,z)-y*z,weights=w)+.5*saved_model["l2"]*np.sum(theta[1:]**2)
                close(objective, saved_model["objective"])
                if np.max(np.abs(gradient)) > plan["training"]["optimizer"]["maximum_admitted_gradient_infinity_norm"]:
                    raise ValueError("serialized weights fail independent stationarity check")
                model_count += 1
                if stage != "development": continue
                saved_scores = predictions[predictions.weighting.eq(weighting) & predictions.model.eq(kind)]
                expected_cells = set(map(tuple, current[current.partition.isin(["validation","test"])][["bacterium","phage"]].to_numpy()))
                if saved_scores.duplicated(["bacterium","phage"]).any() or set(map(tuple,saved_scores[["bacterium","phage"]].to_numpy())) != expected_cells:
                    raise ValueError("development predictions differ from the frozen held-out rows")
                close(saved_scores.logit.to_numpy(), Z[H.index.get_indexer(saved_scores.bacterium),P.index.get_indexer(saved_scores.phage)])
                development_checks += len(saved_scores)
                checked = current.merge(saved_scores[["bacterium","phage","logit"]],on=["bacterium","phage"],validate="one_to_one")
                for (regime, partition), group in checked.groupby(["regime","partition"]):
                    m = metric_table[metric_table.weighting.eq(weighting)&metric_table.model.eq(kind)&metric_table.regime.eq(regime)&metric_table.partition.eq(partition)].iloc[0]
                    close(m.average_precision,average_precision_score(group.infection,group.logit,sample_weight=group.weight))
                    close(m.auroc,roc_auc_score(group.infection,group.logit,sample_weight=group.weight))
                    if regime=="cold_phage" and partition=="validation": close(m.average_precision,best_ap)
    for contrast in report["development_contrasts"]:
        weighting = contrast["weighting"]
        held = rows[rows.regime.eq("cold_phage") & rows.partition.eq("test")].copy()
        values = {}
        for kind in ("bilinear", "additive"):
            lookup = predictions[predictions.weighting.eq(weighting) & predictions.model.eq(kind)].set_index(["bacterium", "phage"]).logit
            values[kind] = lookup.reindex(pd.MultiIndex.from_frame(held[["bacterium", "phage"]])).to_numpy()
        w = np.ones(len(held)) if weighting == "observation_equal" else held.sequence_weight.to_numpy()
        config = plan["development_evaluation"]["bootstrap"]
        replay = grouped_ap_contrast(held, values["bilinear"], values["additive"], w, "cold_phage",
                                    resamples=config["resamples"], seed=config["seed"], minimum_groups=config["minimum_groups"])
        close(contrast["delta_average_precision"], replay["delta_average_precision"])
        close(contrast["interval"], replay["interval"])
        if contrast["improvement_supported"] != replay["improvement_supported"]:
            raise ValueError("primary acceptance decision differs")
    base_features = frame(ROOT / "data/studies/results/basel_ecor_fixed_controls.v1/sequence_features.csv").set_index(["entity_kind","source_build","entity_id"]).sort_index()
    source_path = data_root / "basel_source.csv"
    source = inspect_basel_ecor(source_path.read_bytes(),load_lock(ROOT / "data/studies/basel_ecor_intake.v1.json"))
    truth_base = pd.read_csv(source_path)
    labels = {h:dict(zip(g.phage,g.interaction.astype(int))) for h,g in truth_base.groupby("strain")}
    base_scores, base_rows = frame(bundle / "basel_scores.csv.gz"), frame(bundle / "basel_per_host_bounds.csv")
    if len(base_scores)!=15600 or len(base_rows)!=1200: raise ValueError("BASEL panel support differs")
    near = set(read_json_object(ROOT / "data/studies/results/basel_ecor_fixed_controls.v1/report.json")["training_neighbours_i90_retained"])
    base_checks = 0
    for (weighting,kind,build), group in base_scores.groupby(["weighting","model","host_build"]):
        enc=model_info["encoders"]["full_refit"]
        H=(base_features.loc[("host",build)].loc[source["host_ids"]].to_numpy(float)-enc["host"]["means"])/enc["host"]["scales"]
        P=(base_features.loc[("phage","basel_accession")].loc[source["phage_ids"]].to_numpy(float)-enc["phage"]["means"])/enc["phage"]["scales"]
        theta=archive[model_info["models"][weighting][kind]["full_refit"]["weights_key"]]
        Z=independent_scores(theta,kind,H,P)
        if group.duplicated(["bacterium","phage"]).any() or len(group)!=1300: raise ValueError("incomplete BASEL scores")
        close(group.score.to_numpy(),Z[pd.Index(source["host_ids"]).get_indexer(group.bacterium),pd.Index(source["phage_ids"]).get_indexer(group.phage)])
        base_checks+=len(group)
        independent=[]
        for host, panel in group.groupby("bacterium"):
            scores=dict(zip(panel.phage,panel.score)); unknown=set(source["phage_ids"])-set(labels[host])
            hits={j:[independent_metrics(scores,labels[host],unknown,j,v)[0] for v in (0,1)] for j in range(1,6)}
            for k in (1,2,3,5):
                expected=dict(hit_lower=hits[k][0],hit_upper=hits[k][1],
                              tests_lower=k-sum(hits[j][1] for j in range(1,k)),tests_upper=k-sum(hits[j][0] for j in range(1,k)),
                              unknown_slot_fraction=independent_metrics(scores,labels[host],unknown,k,0)[1],
                              near_training_i90_slot_fraction=independent_metrics(scores,labels[host],near,k,0)[1])
                saved=base_rows[base_rows.weighting.eq(weighting)&base_rows.model.eq(kind)&base_rows.host_build.eq(build)&base_rows.bacterium.eq(host)&base_rows.k.eq(k)]
                if len(saved)!=1: raise ValueError("BASEL host metric duplicated or missing")
                for metric,value in expected.items():close(saved.iloc[0][metric],value)
                independent.append(dict(k=k,**expected))
        summary=report["basel"][weighting+":"+kind+":"+build]["macro"]
        for k, metric_rows in pd.DataFrame(independent).groupby("k"):
            for metric,value in metric_rows.drop(columns="k").mean().items():close(summary[str(k)][metric],value)
    for weighting in model_info["models"]:
        for build in ("published_wgs_primary", "training_assembly_sensitivity"):
            a = report["basel"][weighting+":additive:"+build]["macro"]
            b = report["basel"][weighting+":bilinear:"+build]
            for k, values in b["macro"].items():
                delta = b["hit_minus_additive_conservative_outer_bounds"][k]
                close(delta["lower"], values["hit_lower"]-a[k]["hit_upper"])
                close(delta["upper"], values["hit_upper"]-a[k]["hit_lower"])
    return dict(status="verified", serialized_models=model_count, training_features=271,
                independently_checked_development_scores=development_checks,independently_checked_BASEL_scores=base_checks,
                independently_checked_BASEL_host_K_rows=len(base_rows),
                checks=["frozen input/code/output identities", "primary MLC reconciliation and 41 unknown exclusions",
                        "fit-only unique-assembly scalers", "18 validation-only trials and deterministic selection",
                        "12 serialized models: independent observation-level gradients and objectives",
                        "independent matrix predictions and weighted sklearn AP/AUROC",
                        "replayed primary phage-group bootstrap and acceptance decision",
                        "separate hypergeometric BASEL Hit/cost/exposure bounds"],
                verifier_sha256=sha256_file(Path(__file__)),
                published_output_sha256={p.name:sha256_file(p) for p in sorted(bundle.iterdir()) if p.is_file() and p.name!="NUMERIC_QA.json"})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle",type=Path,required=True);parser.add_argument("--data-root",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True);args=parser.parse_args()
    with threadpool_limits(limits=1):result=verify(args.bundle,args.data_root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("x") as stream:json.dump(result,stream,indent=2);stream.write("\n")
    print("Verified",result["serialized_models"],"models,",result["independently_checked_development_scores"],
          "development and",result["independently_checked_BASEL_scores"],"BASEL scores")


if __name__=="__main__":main()
