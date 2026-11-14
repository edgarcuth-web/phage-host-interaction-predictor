"""Recheck the pinned public-data intake; this never fits or scores a model."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import urllib.request

from phage_cocktail.external_overlap import IntakeRow, audit_overlap

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/studies/public_new_host_intake.v1.json"
SOURCE_STUDY = "10.1038/s41564-024-01832-5"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def checked(data: bytes, lock: dict) -> bytes:
    if len(data) != lock["size_bytes"] or digest(data) != lock["sha256"]:
        raise ValueError(f"Source bytes differ from lock: {lock.get('id', lock.get('path'))}")
    if "git_blob_sha1" in lock:
        blob = b"blob " + str(len(data)).encode() + b"\0" + data
        if hashlib.sha1(blob).hexdigest() != lock["git_blob_sha1"]:
            raise ValueError(f"Git blob differs from lock: {lock['id']}")
    return data


def table(data: bytes, delimiter: str = ";") -> tuple[list[str], list[dict]]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter=delimiter)
    if not reader.fieldnames:
        raise ValueError("Missing source header")
    rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError("Malformed source row")
    return reader.fieldnames, rows


def value_counts(rows: list[dict], columns: list[str]) -> dict:
    return dict(sorted(Counter(row[col] for row in rows for col in columns).items()))


def audit(data: dict[str, bytes], registry: dict, observations: list[dict]) -> tuple[dict, dict]:
    phages = {x["phage_id"]: x for x in registry["library"]["candidates"]}
    hosts = {x["host_id"]: x for x in registry["training"]["known_cohort_members"]}
    columns, matrix = table(data["gaborieau_matrix"])
    source_hosts = {row["bacteria"]: row for row in matrix}
    if len(source_hosts) != len(matrix) or not set(phages) <= set(columns[1:]):
        raise ValueError("Duplicate source host or missing fixed-library candidate")
    crosswalk: dict[str, set[str]] = {}
    same_mlc = same_binary = 0
    for row in observations:
        crosswalk.setdefault(row["bacterium"], set()).add(row["source_bacterium"])
        raw = source_hosts[row["source_bacterium"]][row["phage"]]
        if not raw:
            raise ValueError("Frozen observed cell is missing in original source")
        same_mlc += float(raw) == float(row["source_mlc"])
        same_binary += int(float(raw) > 0) == int(row["source_binary"])
    if set(crosswalk) != set(hosts) or any(len(names) != 1 for names in crosswalk.values()):
        raise ValueError("Frozen development-to-source crosswalk is incomplete or ambiguous")

    _, metadata = table(data["gaborieau_hosts"])
    _, phage_metadata = table(data["gaborieau_phages"])
    _, prospective = table(data["gaborieau_test_hosts"])
    _, selected = table(data["gaborieau_selected"])
    _, composition = table(data["gaborieau_composition"])
    _, steps = table(data["gaborieau_selection_steps"])
    _, cocktail_readouts = table(data["gaborieau_cocktail_readouts"])
    prospective_ids = {row["bacteria"] for row in prospective}
    selected_ids = {row["bacteria"] for row in selected}
    if any({row["bacteria"] for row in rows} != selected_ids
           for rows in (composition, cocktail_readouts)):
        raise ValueError("Measured readout and composition host rosters disagree")
    if {row["bacteria"] for row in steps} != prospective_ids:
        raise ValueError("Selection-step and metadata host rosters disagree")
    pairs = [(row["bacteria"], row[col]) for row in selected for col in ("p1", "p2", "p3")]
    selected_phages = {phage for _, phage in pairs}
    if not selected_phages <= set(phages):
        raise ValueError("Selected source phage absent from pinned library")
    step_pairs = {(row["bacteria"], row["phage"]) for row in steps}
    for row in composition:
        other = next(x for x in selected if x["bacteria"] == row["bacteria"])
        if any(row[col] != other[col] for col in ("p1", "p2", "p3")):
            raise ValueError("Composition and measured singles disagree")

    phl_columns, phl_rows = table(data["phl_interactions"], delimiter=",")
    phl_values = value_counts(phl_rows, phl_columns[1:])
    if set(phl_values) != {"", "0.0", "1.0"}:
        raise ValueError("Unexpected PhageHostLearn label encoding")
    counts = value_counts(matrix, columns[1:])
    fixed_counts = value_counts(matrix, sorted(phages))
    if set(counts) != {"", "0.0", "1.0", "2.0", "3.0", "4.0"}:
        raise ValueError("Unexpected Gaborieau MLC encoding")
    metrics = {
        "gaborieau_full": {
            "matrix_hosts": len(matrix), "matrix_phages": len(columns) - 1,
            "cells": len(matrix) * (len(columns) - 1), "value_counts": counts,
            "metadata_hosts": len(metadata), "metadata_phages": len(phage_metadata),
            "metadata_only_hosts": sorted({r["bacteria"] for r in metadata} - set(source_hosts)),
            "fixed_library_candidates": len(phages), "fixed_library_value_counts": fixed_counts,
            "source_only_phages": sorted(set(columns[1:]) - set(phages)),
            "literal_development_name_overlap": len(set(hosts) & set(source_hosts)),
            "frozen_crosswalk_hosts": len(crosswalk),
            "inherited_observed_cells": len(observations),
            "identical_source_mlc_cells": same_mlc, "identical_source_binary_cells": same_binary,
        },
        "gaborieau_prospective": {
            "metadata_hosts": len(prospective), "unique_host_ids": len(prospective_ids),
            "unique_gembase_ids": len({r["Gembase"] for r in prospective}),
            "selected_single_slots": len(pairs), "unique_selected_pairs": len(set(pairs)),
            "unique_selected_phages": len(selected_phages),
            "metadata_hosts_with_same_readout_name": len(prospective_ids & selected_ids),
            "metadata_only_host_names": sorted(prospective_ids - selected_ids),
            "readout_only_host_names": sorted(selected_ids - prospective_ids),
            "selected_pairs_with_exact_step_support": len(set(pairs) & step_pairs),
            "missing_selected_single_readouts": sum(not r[col].strip() for r in selected
                                                     for col in ("p1_MLC", "p2_MLC", "p3_MLC")),
            "selected_single_value_counts": value_counts(selected, ["p1_MLC", "p2_MLC", "p3_MLC"]),
            "absent_from_selected_single_table": len(prospective_ids) * len(phages) - len(set(pairs)),
            "literal_development_name_overlap": len(prospective_ids & set(hosts)),
            "literal_picard_name_overlap": sorted(prospective_ids & set(source_hosts)),
            "picard_gembase_overlap": len({r["Gembase"] for r in prospective}
                                          & {r["Gembase"] for r in metadata}),
        },
        "phagehostlearn_raw_export": {
            "hosts": len(phl_rows), "phages": len(phl_columns) - 1,
            "cells": len(phl_rows) * (len(phl_columns) - 1), "value_counts": phl_values,
            "literal_fixed_library_phage_id_overlap": len(set(phl_columns[1:]) & set(phages)),
        },
    }

    # Known library identities come from the same pinned source FASTAs. Candidate
    # host assemblies remain unverified: never attach DNA hashes by host name.
    expected_prefix = ("https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/"
                       "99656d8c020fcd1dbbc30a02ee10274a07655562/data/genomics/phages/FNA/")
    if any(p["source_url"] != expected_prefix + name + ".fna" for name, p in phages.items()):
        raise ValueError("Library phage provenance does not match source revision")
    development = [IntakeRow(
        row["bacterium"], row["phage"], hosts[row["bacterium"]]["assembly_sha256"],
        phages[row["phage"]]["assembly_sha256"], SOURCE_STUDY, "",
        "genophi-ecoli:source-mlc-positive:v1",
        "observed_positive" if int(row["source_binary"]) else "observed_negative",
    ) for row in observations]
    candidate = [IntakeRow(
        row["bacteria"], name, "", phage["assembly_sha256"], SOURCE_STUDY, "",
        "gaborieau-picard:exported-mlc-positive:v1",
        "unknown" if not row[name] else
        ("observed_positive" if float(row[name]) > 0 else "observed_negative"),
    ) for row in matrix for name, phage in sorted(phages.items())]
    gate = audit_overlap(development, candidate, regime="cold_host")
    if gate["status"] != "blocked" or "shared_primary_source_lineage" not in gate["blocking_reasons"]:
        raise ValueError("Existing overlap gate failed to block same-source validation")
    return metrics, gate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Fetch only missing, checksum-pinned source tables")
    parser.add_argument("--out", type=Path, required=True, help="New software QA report; will not overwrite")
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    data = {}
    for lock in manifest["resources"]:
        path = args.data_root / lock["cache_path"]
        if not path.exists() and args.fetch:
            with urllib.request.urlopen(lock["url"], timeout=30) as response:
                payload = checked(response.read(lock["size_bytes"] + 1), lock)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        data[lock["id"]] = checked(path.read_bytes(), lock)
    registry = json.loads(checked((ROOT / manifest["registry"]["path"]).read_bytes(), manifest["registry"]))
    observed = checked((ROOT / manifest["observations"]["path"]).read_bytes(), manifest["observations"])
    _, observations = table(gzip.decompress(observed), delimiter=",")
    checked((ROOT / manifest["overlap_gate"]["path"]).read_bytes(), manifest["overlap_gate"])
    metrics, gate = audit(data, registry, observations)
    if metrics != manifest["expected_numeric_checks"]:
        raise ValueError("Numeric intake checks differ from the frozen review")
    report = {
        "schema_version": 1, "intake_id": manifest["intake_id"], "status": "verified",
        "manifest_sha256": digest(manifest_bytes), "verifier_sha256": digest(Path(__file__).read_bytes()),
        "source_files_verified": len(data), "numeric_checks": metrics,
        "same_source_cold_host_gate": gate,
        "decisions": {x["panel_id"]: x["role"] for x in manifest["panels"]},
        "release_cohorts_admitted": 0, "model_fits": 0, "model_score_runs": 0,
        "independent_release_validation_established": False,
        "limits": ["No near-duplicate or new-host genome audit was performed.",
                   "Candidate host hashes stay blank; zero exact-host matches does not establish disjointness.",
                   "A software receipt is not validation of model performance."],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(report, output, indent=2, sort_keys=True)
        output.write("\n")
    print(f"verified {len(data)} source files; same-source validation blocked; admitted release cohorts: 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
