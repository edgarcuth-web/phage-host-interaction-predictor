"""Audit pinned Kaneko supplements; export source grades outside training data.

Uses only Python's standard library. XLSX files are read, never modified.
No binary infection labels, genome identities, model fits or predictions result.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import posixpath
import re
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data/studies/kaneko_workbook_intake.v1.json"
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
ENDPOINT = "kaneko-2025:overnight-spot-clearing-grade-0-3:v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def canonical_digest(rows):
    """Order-independent digest of all fields, not just accession names."""
    normalized = sorted(rows, key=lambda r: json.dumps(r, sort_keys=True))
    return sha256(json.dumps(normalized, sort_keys=True, separators=(",", ":"),
                             ensure_ascii=False).encode())


def column(index):
    result = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(65 + remainder) + result
    return result


def read_xlsx(path, sheet_name):
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ["".join(e.itertext()) for e in
                      ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook.findall("s:sheets/s:sheet", NS)
        require(len(sheets) == 1 and sheets[0].get("name") == sheet_name,
                f"Unexpected sheets in {path.name}")
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_id = sheets[0].get(f"{{{REL}}}id")
        target = next(e.get("Target") for e in rels if e.get("Id") == rel_id)
        member = (target.lstrip("/") if target.startswith("/") else
                  posixpath.normpath(posixpath.join("xl", target)))
        worksheet = ET.fromstring(archive.read(member))
        result = {}
        for cell in worksheet.findall("s:sheetData/s:row/s:c", NS):
            address = cell.get("r")
            require(address not in result, f"Duplicate cell {address}")
            require(cell.find("s:f", NS) is None, f"Formula at {address}")
            value = cell.find("s:v", NS)
            kind = cell.get("t")
            if kind == "inlineStr":
                value = "".join(cell.find("s:is", NS).itertext())
            elif value is None or value.text is None:
                value = None
            elif kind == "s":
                value = shared[int(value.text)]
            elif kind in (None, "n"):
                raw = value.text
                value = int(raw) if re.fullmatch(r"-?\d+", raw) else float(raw)
            else:
                raise ValueError(f"Unsupported cell type {kind} at {address}")
            result[address] = value
        return result


def parse_matrix(cells):
    hosts = [cells.get(f"{column(c)}2") for c in range(2, 146)]
    require(hosts == [f"{i:03d}" for i in range(1, 145)], "Host axis changed")
    phages = [cells.get(f"A{r}") for r in range(3, 48)]
    require(len(set(phages)) == 45 and phages[:3] == ["T2", "T4", "T7"],
            "Missing, duplicate or reordered control phages")
    require(all(isinstance(p, str) and re.fullmatch(r"ΦWec\d+", p)
                for p in phages[3:]), "Study phage axis changed")
    rows = []
    for r, phage in enumerate(phages, 3):
        for c, host in enumerate(hosts, 2):
            address = f"{column(c)}{r}"
            value = cells.get(address)
            require(type(value) is int and value in (0, 1, 2, 3),
                    f"Missing or invalid spot grade at SF2!{address}")
            rows.append({"host_label": host, "phage_label": phage,
                         "phage_role": "study" if phage.startswith("Φ") else "control",
                         "source_grade": value, "endpoint_id": ENDPOINT,
                         "source_sheet": "SF2", "source_cell": address})
    return hosts, phages, rows


def reconcile_s5(cells, rows):
    require([cells.get(f"{c}2") for c in "BCD"] ==
            ["score > 1", "score > 2", "score > 3"], "S5 headings changed")
    by_phage = {}
    for row in rows:
        by_phage.setdefault(row["phage_label"], []).append(row["source_grade"])
    names = [cells.get(f"A{r}") for r in range(3, 48)]
    require(len(set(names)) == 45 and set(names) == set(by_phage), "S5 phage mismatch")
    for r, name in enumerate(names, 3):
        for c, threshold in zip("BCD", (1, 2, 3)):
            require(cells.get(f"{c}{r}") == sum(x >= threshold for x in by_phage[name]),
                    f"S5 aggregate mismatch: {name}, threshold {threshold}")
    return {"matching_phages": 45, "matching_counts": 135,
            "printed_operator": ">", "operator_reproducing_all_counts": ">=",
            "binary_threshold_selected": False}


def parse_metadata(data, resource):
    reader = csv.DictReader(io.StringIO(data.decode()), delimiter="\t")
    require(len(reader.fieldnames or []) == len(set(reader.fieldnames or [])),
            "Duplicate metadata headers")
    rows = list(reader)
    require(len(rows) == resource["rows"], "Metadata row count changed")
    require(canonical_digest(rows) == resource["canonical_rows_sha256"],
            f"Metadata content changed: {resource['id']}")
    return rows


def join_metadata(hosts, phages, metadata):
    host_join = []
    for row in metadata["host"]:
        match = re.fullmatch(r"JMU_No_([1-9]\d*)_read", row["sample_title"])
        require(match is not None, "Unexpected host title")
        host_join.append({"host_label": f"{int(match[1]):03d}", **row,
                          "join_evidence": "paper host number to ENA sample title"})
    require(len({r["host_label"] for r in host_join}) == len(hosts) and
            {r["host_label"] for r in host_join} == set(hosts), "Host title join failed")
    require(all(r["study_accession"] == "PRJDB18063" and r["fastq_ftp"]
                for r in host_join), "Host project/read links changed")
    require({r["experiment_accession"] for r in host_join} ==
            {f"DRX{i}" for i in range(534018, 534162)}, "Host experiment range changed")
    phage_join = []
    for row in metadata["phage"]:
        match = re.fullmatch(r"phi(Wec\d+)_sequence_read", row["sample_title"])
        require(match is not None and row["study_accession"] == "PRJDB18044" and
                row["fastq_ftp"], "Unexpected phage metadata")
        phage_join.append({"phage_label": "Φ" + match[1], "sequence_kind": "raw_reads",
                           "accession": row["experiment_accession"], **row})
    require({r["experiment_accession"] for r in phage_join} ==
            {f"DRX{i}" for i in range(534162, 534200)}, "Phage experiment range changed")
    for row in metadata["legacy"]:
        match = re.fullmatch(r"Escherichia phage phi(Wec\d+) DNA, complete genome\.",
                             row["description"])
        require(match is not None, "Unexpected legacy sequence description")
        phage_join.append({"phage_label": "Φ" + match[1],
                           "sequence_kind": "complete_genome_record_metadata", **row})
    study_phages = set(phages) - {"T2", "T4", "T7"}
    require(len(phage_join) == 42 and len({r["phage_label"] for r in phage_join}) == 42
            and {r["phage_label"] for r in phage_join} == study_phages,
            "Study phage metadata coverage is not exactly 42")
    for key in ("sample_accession", "experiment_accession", "run_accession"):
        for name in ("host", "phage"):
            require(len({r[key] for r in metadata[name]}) == len(metadata[name]),
                    f"Duplicate {key} in {name}")
    return sorted(host_join, key=lambda r: r["host_label"]), sorted(
        phage_join, key=lambda r: r["phage_label"])


def audit(xlsx_dir, metadata_dir, manifest):
    candidates = {}
    for path in xlsx_dir.glob("*.xlsx"):
        candidates.setdefault(sha256(path.read_bytes()), []).append(path)
    sheets, file_checks = {}, []
    for resource in manifest["workbooks"]:
        matches = candidates.get(resource["sha256"], [])
        require(len(matches) == 1, f"Need exactly one pinned {resource['sheet']} workbook")
        path = matches[0]
        require(path.stat().st_size == resource["size_bytes"], "Workbook size changed")
        sheets[resource["sheet"]] = read_xlsx(path, resource["sheet"])
        file_checks.append({"id": resource["id"], "sha256": resource["sha256"],
                            "bytes_match": True})
    hosts, phages, rows = parse_matrix(sheets["SF2"])
    require([sheets["ST1"].get(f"A{r}") for r in range(4, 148)] ==
            [f"No.{i}" for i in range(1, 145)], "ST1 host identity mismatch")
    all_phages = [sheets["ST2"].get(f"A{r}") for r in range(3, 49)]
    require(len(set(all_phages)) == 46 and set(all_phages) == set(phages) | {"T6"},
            "ST2 phage inventory mismatch")
    for name in ("ST3", "ST4"):
        require([sheets[name].get(f"{column(c)}2") for c in range(2, 48)] == all_phages
                and [sheets[name].get(f"A{r}") for r in range(3, 49)] == all_phages,
                f"{name} ANI axes mismatch")
        require(all(sheets[name].get(f"{column(i+2)}{i+3}") == 100
                    for i in range(46)), f"{name} ANI diagonal mismatch")
    for address, expected in manifest["captions_sha256"].items():
        sheet, cell = address.split("!")
        require(sha256(sheets[sheet][cell].encode()) == expected, "Assay caption changed")
    aggregates = reconcile_s5(sheets["ST5"], rows)
    metadata = {}
    metadata_checks = []
    for resource in manifest["metadata"]:
        data = (metadata_dir / resource["filename"]).read_bytes()
        metadata[resource["id"]] = parse_metadata(data, resource)
        metadata_checks.append({"id": resource["id"], "sha256": sha256(data),
                                "snapshot_bytes_match": sha256(data) == resource["sha256"],
                                "canonical_rows_match": True})
    host_join, phage_join = join_metadata(hosts, phages, metadata)
    def counts(selected):
        counted = Counter(r["source_grade"] for r in selected)
        return {str(i): counted[i] for i in range(4)}
    summary = {"schema_version": 1, "review_id": manifest["review_id"],
               "status": "audited_source_ordinal_data_not_admitted_for_training",
               "workbook_checks": file_checks, "metadata_checks": metadata_checks,
               "hosts": 144, "study_phages": 42, "measured_control_phages": ["T2", "T4", "T7"],
               "inventory_only_control_not_measured_in_SF2": "T6",
               "observed_source_cells": len(rows), "study_source_cells": 6048,
               "control_source_cells": 432, "missing_matrix_cells": 0,
               "grade_counts_all": counts(rows),
               "grade_counts_study": counts([r for r in rows if r["phage_role"] == "study"]),
               "s5_reconciliation": aggregates, "host_number_to_read_record_joins": len(host_join),
               "phage_name_to_read_record_joins": 38, "phage_name_to_legacy_genome_metadata_joins": 4,
               "endpoint_id": ENDPOINT, "binary_threshold_selected": False,
               "assembly_and_stock_identity_verified": False,
               "training_pairs_admitted": 0, "independent_release_validation_established": False,
               "model_fits": 0, "model_score_runs": 0}
    require(summary["grade_counts_study"] == manifest["expected_grade_counts_study"],
            "Study grade counts changed")
    return summary, rows, host_join, phage_join


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xlsx-dir", type=Path, required=True)
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--fetch-metadata", action="store_true")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    require(not args.out_dir.exists(), "Output directory must be new")
    manifest = json.loads(args.manifest.read_text())
    if args.fetch_metadata:
        args.metadata_dir.mkdir(parents=True, exist_ok=False)
        for resource in manifest["metadata"]:
            with urllib.request.urlopen(resource["url"], timeout=30) as response:
                data = response.read()
            (args.metadata_dir / resource["filename"]).write_bytes(data)
    summary, rows, hosts, phages = audit(args.xlsx_dir, args.metadata_dir, manifest)
    args.out_dir.mkdir(parents=True)
    source = next(r for r in manifest["workbooks"] if r["sheet"] == "SF2")
    for row in rows:
        row.update({"source_study_doi": manifest["doi"],
                    "source_filename": source["filename"], "source_sha256": source["sha256"]})
    for filename, output in (("source_spot_grades.csv", rows),
                             ("host_read_crosswalk.csv", hosts),
                             ("phage_sequence_routes.csv", phages)):
        write_csv(args.out_dir / filename, output)
    summary["output_sha256"] = {p.name: sha256(p.read_bytes())
                                 for p in sorted(args.out_dir.glob("*.csv"))}
    (args.out_dir / "AUDIT.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("status", "observed_source_cells",
                                             "study_source_cells", "training_pairs_admitted")}))


if __name__ == "__main__":
    main()
