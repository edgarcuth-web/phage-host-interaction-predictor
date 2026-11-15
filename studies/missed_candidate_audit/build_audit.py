#!/usr/bin/env python3
"""Audit saved 94-candidate rankings against pinned source readings; no fitting."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import urllib.request
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

SNAPSHOT = "5962f5caed0221455ec34b542e400a6438bb49e7"
SOURCE_REF = "99656d8c020fcd1dbbc30a02ee10274a07655562"
SOURCE_URL = f"https://github.com/mdmparis/coli_phage_interactions_2023/blob/{SOURCE_REF}/data/interactions/interaction_matrix.csv"
RAW_URL = f"https://github.com/mdmparis/coli_phage_interactions_2023/blob/{SOURCE_REF}/data/interactions/raw/raw_interactions.csv"
REPO_URL = f"https://github.com/edgarcuth-web/phage-host-interaction-predictor/blob/{SNAPSHOT}"
FOCUS = ("NILS53", "NILS78", "AL505")
METHODS = ("phage_prevalence", "host_4mer", "genophi", "annotation", "phage_only", "additive", "bilinear")
STATUS = {
    "hit": "已命中",
    "ranking_miss": "有阳性但排序漏掉",
    "no_source_positive": "库里没有来源阳性",
    "unknown": "数据未知（面板不完整）",
    "tie_risk": "并列导致部分漏检风险",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path, delimiter=","):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def dump_csv(path, records):
    if not records:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=list(records[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    data = text.getvalue().encode("utf-8-sig")
    if path.suffix == ".gz":
        data = gzip.compress(data, mtime=0)
    path.write_bytes(data)


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def no_hit_probability(blocks, labels, k, unknown_positive=False):
    """Exact expectation over uniform permutations within exact-score blocks."""
    result = Fraction(1)
    remaining = k
    for _, members in blocks:
        if not remaining:
            break
        take = min(remaining, len(members))
        positive = sum(labels[p] == 1 or (labels[p] is None and unknown_positive) for p in members)
        negative = len(members) - positive
        if negative < take:
            return Fraction(0)
        result *= Fraction(math.comb(negative, take), math.comb(len(members), take))
        remaining -= take
    assert remaining == 0
    return result


def rank_candidates(local):
    ordered = sorted(local, key=lambda row: (-float(row["score"]), row["phage"]))
    score_blocks = defaultdict(list)
    for row in ordered:
        score_blocks[float(row["score"])].append(row["phage"])
    blocks = list(score_blocks.items())
    ranks = {}
    offset = 0
    for score, names in blocks:
        for ordinal, phage in enumerate(names, offset + 1):
            ranks[phage] = {
                "score": score,
                "rank_min": offset + 1,
                "rank_max": offset + len(names),
                "display_rank": ordinal,
                "top2_inclusion_probability": min(max(2 - offset, 0), len(names)) / len(names),
            }
        offset += len(names)
    return blocks, ranks


def source_evidence(root):
    repo = root / "repo"
    observed = rows(repo / "data/studies/results/paired_interaction.v1/observed_pairs.csv.gz")
    excluded = rows(repo / "data/studies/results/paired_interaction.v1/excluded_unknown_pairs.csv")
    subset = rows(root / "sources/subset.csv")
    matrix = rows(root / "sources/mlc.csv", delimiter=";")
    assert len(observed) == 16597 and len(excluded) == 41 and len(subset) == 16638
    aliases = {}
    for row in observed + excluded:
        host, source_host = row["bacterium"], row["source_bacterium"]
        assert aliases.get(host, source_host) == source_host
        aliases[host] = source_host
    library = sorted({r["phage"] for r in subset})
    assert len(aliases) == 177 and len(library) == 94
    assert len({(r["strain"], r["phage"]) for r in subset}) == 16638
    source = {r["bacteria"]: (line, r) for line, r in enumerate(matrix, 2)}
    columns = {name: index for index, name in enumerate(matrix[0], 1)}
    observed_map = {(r["bacterium"], r["phage"]): r for r in observed}
    unknown_map = {(r["bacterium"], r["phage"]): r for r in excluded}
    evidence = {}
    for subset_line, row in enumerate(subset, 2):
        host, phage = row["strain"], row["phage"]
        key = (host, phage)
        source_host = aliases[host]
        matrix_line, matrix_row = source[source_host]
        original = matrix_row[phage]
        missing = original.strip().lower() in {"", "nan", "na", "none"}
        label = None if missing else int(float(original) > 0)
        if label is None:
            assert key in unknown_map and key not in observed_map
            assert int(row["interaction"]) == 0
        else:
            assert key in observed_map and key not in unknown_map
            expected = observed_map[key]
            assert float(original) == float(expected["source_mlc"])
            assert label == int(expected["infection"]) == int(row["interaction"])
        evidence[key] = {
            "evidence_id": host + "::" + phage,
            "bacterium": host,
            "source_bacterium": source_host,
            "phage": phage,
            "source_mlc_original": original,
            "source_mlc_numeric": None if missing else float(original),
            "source_label": label,
            "label_status": "未知" if label is None else ("阳性" if label else "阴性"),
            "geno_subset_original_binary": int(row["interaction"]),
            "source_matrix_row": matrix_line,
            "source_matrix_column": columns[phage],
            "source_matrix_header": phage,
            "subset_row": subset_line,
            "matrix_url": SOURCE_URL + f"#L{matrix_line}",
            "raw_annotation_count": 0,
            "raw_score1_count": 0,
            "raw_score0_count": 0,
            "raw_score_n_count": 0,
            "raw_row_numbers": [],
            "raw_images": set(),
        }
    assert set(evidence) == {(h, p) for h in aliases for p in library}
    return aliases, library, evidence


def load_raw_annotations(root, aliases, evidence):
    inverse = {source: host for host, source in aliases.items()}
    assert len(inverse) == len(aliases)
    raw_by_pair = defaultdict(list)
    with (root / "sources/raw_interactions.csv").open(newline="", encoding="utf-8") as handle:
        for line, row in enumerate(csv.DictReader(handle, delimiter=";"), 2):
            host = inverse.get(row["bacteria"])
            key = (host, row["phage"])
            if key not in evidence:
                continue
            entry = evidence[key]
            entry["raw_annotation_count"] += 1
            assert row["score"] in {"0", "1", "n"}
            count_column = "raw_score_n_count" if row["score"] == "n" else f"raw_score{row['score']}_count"
            entry[count_column] += 1
            entry["raw_row_numbers"].append(line)
            entry["raw_images"].add(row["image"])
            raw_by_pair[key].append({
                "evidence_id": entry["evidence_id"],
                "source_raw_row": line,
                **row,
            })
    for entry in evidence.values():
        entry["raw_row_numbers"] = ",".join(map(str, entry["raw_row_numbers"]))
        entry["raw_images"] = "; ".join(sorted(entry["raw_images"]))
        entry["raw_evidence_status"] = "逐次原始注释可追溯" if entry["raw_annotation_count"] else "公开原始注释表未收录（仅汇总矩阵）"
    return raw_by_pair


def audit(root, out):
    out.mkdir(parents=True, exist_ok=False)
    repo = root / "repo"
    aliases, library, evidence = source_evidence(root)
    raw_by_pair = load_raw_annotations(root, aliases, evidence)
    source_hosts = []
    for host in sorted(aliases):
        local = [evidence[host, p] for p in library]
        pos = sum(r["source_label"] == 1 for r in local)
        neg = sum(r["source_label"] == 0 for r in local)
        unknown = sum(r["source_label"] is None for r in local)
        source_hosts.append({
            "bacterium": host, "source_bacterium": aliases[host],
            "library_size": 94, "positive_candidates": pos,
            "negative_candidates": neg, "unknown_candidates": unknown,
            "source_status": "数据未知（面板不完整）" if unknown else ("库里没有来源阳性" if not pos else "库内有来源阳性"),
            "raw_annotation_count": sum(r["raw_annotation_count"] for r in local),
        })
        assert pos + neg + unknown == 94
    host_source = {r["bacterium"]: r for r in source_hosts}
    split_rows = rows(repo / "data/studies/results/genophi_comparator.v1/scenario_splits.csv")
    split_map = {(r["scenario"], r["bacterium"]): r for r in split_rows}
    training = defaultdict(list)
    for row in split_rows:
        assert aliases[row["bacterium"]] == row["source_bacterium"]
        if row["split"] == "train":
            training[row["scenario"]].append(row["bacterium"])
    support = {}
    for panel, hosts in training.items():
        assert len(hosts) == (117 if panel == "grouped_host" else 115)
        for phage in library:
            values = [evidence[h, phage]["source_label"] for h in hosts]
            support[panel, phage] = (sum(v == 1 for v in values), sum(v is not None for v in values))
    score_file = repo / "data/studies/results/new_host_retrieval.v1/scores.csv.gz"
    saved_scores = rows(score_file)
    assert len(saved_scores) == 91744
    groups = defaultdict(list)
    for row in saved_scores:
        assert math.isfinite(float(row["score"]))
        groups[row["panel"], row["partition"], row["method"], row["bacterium"]].append(row)
    original_metrics = {(r["panel"], r["partition"], r["method"], r["bacterium"]): r
                        for r in rows(repo / "data/studies/results/new_host_retrieval.v1/per_host_metrics.csv")}
    assert len(original_metrics) == 968
    host_ledger, random_ledger, ranked_rows = [], [], []
    candidate_rows, top2_rows, focus_rows = [], [], []
    evidence_export_keys = set()
    metric_checks = 0
    max_error = 0.0
    for key in sorted(groups):
        panel, partition, method, host = key
        local = groups[key]
        assert len(local) == 94 and {r["phage"] for r in local} == set(library)
        context = split_map[panel, host]
        assert context["split"] == partition and partition in {"validation", "test"}
        labels = {p: evidence[host, p]["source_label"] for p in library}
        blocks, ranks = rank_candidates(local)
        hits = {k: 1 - no_hit_probability(blocks, labels, k) for k in (1, 2, 3, 4, 5)}
        upper = 1 - no_hit_probability(blocks, labels, 2, unknown_positive=True)
        src = host_source[host]
        if src["unknown_candidates"]:
            status_key = "unknown"
            assert key not in original_metrics
        elif not src["positive_candidates"]:
            status_key = "no_source_positive"
        elif hits[2] == 0:
            status_key = "ranking_miss"
        elif hits[2] < 1:
            status_key = "tie_risk"
        else:
            status_key = "hit"
        if not src["unknown_candidates"]:
            saved = original_metrics[key]
            assert int(saved["positive_candidates"]) == src["positive_candidates"]
            for k in (1, 2, 3, 4, 5):
                error = abs(float(hits[k]) - float(saved[f"hit@{k}"]))
                max_error = max(max_error, error)
                assert error < 1e-12, (key, k, error)
                metric_checks += 1
            assert abs(float(2 - hits[1]) - float(saved["tests@2"])) < 1e-12
        ordered = [p for _, members in blocks for p in members]
        positive_order = [p for p in ordered if labels[p] == 1]
        possible_top2 = [p for p in ordered if ranks[p]["rank_min"] <= 2]
        certain_top2 = [p for p in ordered if ranks[p]["rank_max"] <= 2]
        top2_text = "; ".join(p + " (MLC=" + ("未知" if labels[p] is None else evidence[host, p]["source_mlc_original"]) + ")" for p in possible_top2)
        best = positive_order[0] if positive_order else ""
        ledger = {
            "panel": panel, "partition": partition, "method": method,
            "bacterium": host, "source_bacterium": aliases[host],
            "host_group": context["host_group"], "phylogroup": context["phylogroup"],
            "diagnosis": STATUS[status_key], "diagnosis_code": status_key,
            "positive_candidates": src["positive_candidates"],
            "negative_candidates": src["negative_candidates"],
            "unknown_candidates": src["unknown_candidates"],
            "hit2_lower": float(hits[2]), "hit2_upper": float(upper),
            "top2_outcome": ("已知读数已命中" if hits[2] == 1 else
                             "已知读数未命中" if upper == 0 else
                             "缺失读数影响判定" if hits[2] != upper else "并列带来部分命中概率"),
            "official_complete_panel_eligible": not bool(src["unknown_candidates"]),
            "top2_possible_candidates": "; ".join(possible_top2),
            "top2_guaranteed_candidates": "; ".join(certain_top2),
            "top2_source_readings": top2_text,
            "top2_boundary_tie": len(possible_top2) > 2,
            "best_positive_candidate": best,
            "best_positive_rank_min": ranks[best]["rank_min"] if best else None,
            "best_positive_rank_max": ranks[best]["rank_max"] if best else None,
            "focus_case": host in FOCUS,
        }
        if method == "random":
            random_ledger.append(ledger)
            continue
        assert method in METHODS
        host_ledger.append(ledger)
        is_miss = status_key in {"ranking_miss", "tie_risk"}
        for phage in ordered:
            entry, rank = evidence[host, phage], ranks[phage]
            train_pos, train_observed = support[panel, phage]
            pair = {
                "panel": panel, "partition": partition, "method": method,
                "bacterium": host, "phage": phage,
                "host_diagnosis": STATUS[status_key],
                "candidate_role": ("漏检宿主的已知阳性" if is_miss else "重点宿主阳性对照") if entry["source_label"] == 1 else "前两名候选对照",
                "source_label": entry["source_label"],
                "source_mlc_original": entry["source_mlc_original"],
                **rank,
                "training_positive_hosts": train_pos,
                "training_observed_hosts": train_observed,
                "training_prevalence_beta11": (train_pos + 1) / (train_observed + 2),
                "evidence_id": entry["evidence_id"],
                "source_matrix_row": entry["source_matrix_row"],
                "source_matrix_column": entry["source_matrix_column"],
                "raw_annotation_count": entry["raw_annotation_count"],
                "raw_evidence_status": entry["raw_evidence_status"],
                "focus_case": host in FOCUS,
            }
            if method == "phage_prevalence":
                assert abs(pair["score"] - pair["training_prevalence_beta11"]) < 1e-15
            ranked_rows.append(pair)
            if rank["rank_min"] <= 2:
                top2_rows.append(pair)
                evidence_export_keys.add((host, phage))
            if entry["source_label"] == 1 and (is_miss or host in FOCUS):
                candidate_rows.append(pair)
                evidence_export_keys.add((host, phage))
            if host in FOCUS and (entry["source_label"] == 1 or rank["rank_min"] <= 2):
                focus_rows.append(pair)
        if src["unknown_candidates"]:
            evidence_export_keys.update((host, p) for p in library if labels[p] is None)
    assert len(host_ledger) == 854 and len(random_ledger) == 122 and len(ranked_rows) == 80276
    assert metric_checks == 4840
    misses = [r for r in host_ledger if r["diagnosis_code"] != "hit"]
    summary = []
    for panel, part, method in sorted({(r["panel"], r["partition"], r["method"]) for r in host_ledger}):
        local = [r for r in host_ledger if (r["panel"], r["partition"], r["method"]) == (panel, part, method)]
        counts = Counter(r["diagnosis_code"] for r in local)
        summary.append({
            "panel": panel, "partition": part, "method": method,
            "scored_hosts": len(local), "complete_hosts": sum(r["official_complete_panel_eligible"] for r in local),
            "hit_hosts": counts["hit"], "ranking_miss_hosts": counts["ranking_miss"],
            "no_source_positive_hosts": counts["no_source_positive"],
            "unknown_hosts": counts["unknown"], "tie_risk_hosts": counts["tie_risk"],
            "complete_host_hit2": sum(r["hit2_lower"] for r in local if r["official_complete_panel_eligible"]) / sum(r["official_complete_panel_eligible"] for r in local),
        })
    for k, entry in evidence.items():
        if entry["source_label"] is None:
            evidence_export_keys.add(k)
    exported_evidence = [evidence[k] for k in sorted(evidence_export_keys)]
    unknown = [evidence[k] for k in sorted(evidence) if evidence[k]["source_label"] is None]
    raw_export = [r for k in sorted(evidence_export_keys) for r in raw_by_pair[k]]
    focus_lookup = {(r["panel"], r["method"], r["bacterium"], r["phage"]): r for r in focus_rows if r["source_label"] == 1}
    focus_comparison = []
    for host in FOCUS:
        for phage in library:
            source = evidence[host, phage]
            if source["source_label"] != 1:
                continue
            comparison = {"bacterium": host, "phage": phage, "source_mlc": source["source_mlc_numeric"]}
            for panel, prefix in (("grouped_host", "primary"), ("phylogroup:B2", "b2")):
                for method in ("phage_prevalence", "host_4mer", "genophi"):
                    record = focus_lookup[panel, method, host, phage]
                    for bound in ("min", "max"):
                        comparison[f"{prefix}_{method}_rank_{bound}"] = record[f"rank_{bound}"]
            comparison.update({k: source[k] for k in ("evidence_id", "raw_annotation_count", "raw_evidence_status", "source_matrix_row", "source_matrix_column", "matrix_url")})
            focus_comparison.append(comparison)
    assert len(focus_comparison) == 59
    for name, records in [
        ("source_host_status.csv", source_hosts),
        ("all_host_method_decisions.csv", host_ledger),
        ("missed_hosts.csv", misses),
        ("summary.csv", summary),
        ("rescue_candidates.csv", candidate_rows),
        ("top2_candidates.csv", top2_rows),
        ("focus_cases.csv", focus_rows),
        ("focus_positive_candidates.csv", focus_comparison),
        ("source_evidence.csv", exported_evidence),
        ("unknown_pairs.csv", unknown),
        ("ranked_candidates.csv.gz", ranked_rows),
        ("raw_annotations.csv.gz", raw_export),
        ("random_reference.csv", random_ledger),
    ]:
        dump_csv(out / name, records)
    qa = {
        "status": "passed", "snapshot": SNAPSHOT,
        "source_hosts": 177, "candidate_phages": 94,
        "observed_pairs": 16597, "source_positive_pairs": sum(r["source_label"] == 1 for r in evidence.values()),
        "unknown_pairs": 41, "scored_host_panel_occurrences": 122,
        "unique_scored_hosts": len({r["bacterium"] for r in host_ledger}),
        "deterministic_host_method_records": len(host_ledger),
        "random_reference_records": len(random_ledger),
        "ranked_candidate_rows": len(ranked_rows),
        "original_metric_comparisons": metric_checks,
        "maximum_hit_metric_error": max_error,
        "miss_or_unknown_records": len(misses),
        "confirmed_ranking_miss_records": sum(r["diagnosis_code"] == "ranking_miss" for r in host_ledger),
        "tie_risk_records": sum(r["diagnosis_code"] == "tie_risk" for r in host_ledger),
        "exported_source_pairs": len(exported_evidence),
        "exported_raw_annotations": len(raw_export),
        "no_source_positive_hosts_all177": [r["bacterium"] for r in source_hosts if r["source_status"] == "库里没有来源阳性"],
        "incomplete_hosts_all177": [{k: r[k] for k in ("bacterium", "positive_candidates", "unknown_candidates")} for r in source_hosts if r["unknown_candidates"]],
        "focus_positive_counts": {h: host_source[h]["positive_candidates"] for h in FOCUS},
        "focus_raw_annotation_counts": {h: host_source[h]["raw_annotation_count"] for h in FOCUS},
        "new_fits": 0, "new_model_inference": 0,
    }
    assert qa["source_positive_pairs"] == 3813
    assert qa["focus_positive_counts"] == {"NILS53": 11, "NILS78": 7, "AL505": 41}
    dump_json(out / "QA.json", qa)
    dump_json(out / "workbook_data.json", {
        "summary": summary, "ledger": host_ledger,
        "candidates": candidate_rows, "evidence": exported_evidence, "unknown": unknown,
        "focus": focus_rows, "focus_comparison": focus_comparison, "qa": qa,
    })
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("data/raw/missed_candidate_audit"))
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("input_manifest.json"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    checked = []
    for item in manifest["inputs"]:
        target = args.input_root / item["path"]
        if not target.exists() and args.fetch:
            target.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(item["url"], timeout=60) as response:
                target.write_bytes(response.read())
        assert target.is_file(), f"Missing input: {target}. Use --fetch to obtain pinned public inputs."
        assert target.stat().st_size == item["size_bytes"], f"Size mismatch: {target}"
        assert sha256(target) == item["sha256"], f"SHA-256 mismatch: {target}"
        if item.get("git_blob_sha"):
            data = target.read_bytes()
            digest = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            assert digest == item["git_blob_sha"], f"Git blob mismatch: {target}"
        checked.append(item["path"])
    qa = audit(args.input_root, args.out)
    qa["verified_input_files"] = len(checked)
    dump_json(args.out / "QA.json", qa)
    print(json.dumps(qa, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
