#!/usr/bin/env python3
"""Independent cutoff-based checks. Does not import the audit builder."""
import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path


def read(path, delimiter=","):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    base, out = args.input_root, args.out
    primary = read(base / "sources/mlc.csv", ";")
    source = {r["bacteria"]: (i, r) for i, r in enumerate(primary, 2)}
    labels = {}
    checked = 0
    source_rows = read(out / "source_evidence.csv")
    for r in source_rows:
        line, original = source[r["source_bacterium"]]
        assert line == int(r["source_matrix_row"])
        assert list(original)[int(r["source_matrix_column"]) - 1] == r["phage"]
        assert original[r["phage"]] == r["source_mlc_original"]
        checked += 1
    observed = read(base / "repo/data/studies/results/paired_interaction.v1/observed_pairs.csv.gz")
    for r in observed:
        original = source[r["source_bacterium"]][1][r["phage"]]
        assert float(original) == float(r["source_mlc"])
        labels[r["bacterium"], r["phage"]] = int(float(original) > 0)
    excluded = read(base / "repo/data/studies/results/paired_interaction.v1/excluded_unknown_pairs.csv")
    for r in excluded:
        assert source[r["source_bacterium"]][1][r["phage"]] == ""
        labels[r["bacterium"], r["phage"]] = None
    saved = read(base / "repo/data/studies/results/new_host_retrieval.v1/scores.csv.gz")
    scores = defaultdict(dict)
    for r in saved:
        scores[r["panel"], r["partition"], r["method"], r["bacterium"]][r["phage"]] = float(r["score"])
    ledger = read(out / "all_host_method_decisions.csv")
    total_checks = 0
    for r in ledger:
        h = r["bacterium"]
        values = scores[r["panel"], r["partition"], r["method"], h]
        cutoff = sorted(values.values(), reverse=True)[1]
        above = [phage for phage in values if values[phage] > cutoff]
        tied = [phage for phage in values if values[phage] == cutoff]
        remaining = 2 - len(above)
        probabilities = []
        for optimistic in (False, True):
            positive = lambda phage: labels[h, phage] == 1 or (optimistic and labels[h, phage] is None)
            if any(positive(phage) for phage in above):
                probability = 1.0
            else:
                n = sum(not positive(phage) for phage in tied)
                probability = 1 - (math.comb(n, remaining) / math.comb(len(tied), remaining) if n >= remaining else 0)
            probabilities.append(probability)
        assert abs(probabilities[0] - float(r["hit2_lower"])) < 1e-12
        assert abs(probabilities[1] - float(r["hit2_upper"])) < 1e-12
        assert set(r["top2_possible_candidates"].split("; ")) == set(above + tied)
        ys = [labels[h, phage] for phage in values]
        pos, unk = sum(v == 1 for v in ys), sum(v is None for v in ys)
        assert pos == int(r["positive_candidates"]) and unk == int(r["unknown_candidates"])
        if r["diagnosis_code"] == "ranking_miss":
            assert pos > 0 and unk == 0 and probabilities == [0.0, 0.0]
        if r["diagnosis_code"] == "no_source_positive":
            assert pos == 0 and unk == 0
        if r["diagnosis_code"] == "unknown":
            assert unk > 0
        total_checks += 1
    ranked = read(out / "ranked_candidates.csv.gz")
    for r in ranked:
        values = scores[r["panel"], r["partition"], r["method"], r["bacterium"]]
        value = values[r["phage"]]
        assert float(r["score"]) == value
        assert int(r["rank_min"]) == 1 + sum(x > value for x in values.values())
        assert int(r["rank_max"]) == sum(x >= value for x in values.values())
    raw_export = read(out / "raw_annotations.csv.gz")
    wanted = {int(r["source_raw_row"]): r for r in raw_export}
    assert len(wanted) == len(raw_export)
    matched = 0
    raw_counts = defaultdict(lambda: [0, 0, 0, 0])
    with (base / "sources/raw_interactions.csv").open(encoding="utf-8", newline="") as f:
        for line, original in enumerate(csv.DictReader(f, delimiter=";"), 2):
            if line in wanted:
                exported = wanted[line]
                assert all(exported[k] == v for k, v in original.items())
                matched += 1
                count = raw_counts[exported["evidence_id"]]
                count[0] += 1
                count[{"1": 1, "0": 2, "n": 3}[original["score"]]] += 1
    assert matched == len(raw_export)
    for r in source_rows:
        assert raw_counts[r["evidence_id"]] == [int(r[c]) for c in ("raw_annotation_count", "raw_score1_count", "raw_score0_count", "raw_score_n_count")]
    result = {
        "status": "passed", "independent_of_builder": True,
        "source_matrix_cells_checked": checked, "observed_source_labels_checked": len(observed),
        "unknown_source_cells_checked": len(excluded), "host_cutoff_checks": total_checks,
        "candidate_rank_intervals_checked": len(ranked), "raw_annotation_rows_checked": matched,
    }
    (out / "INDEPENDENT_QA.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
