"""Verify source-linked Kaneko assemblies and report sequence QC, without training."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tempfile
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from phage_cocktail.cohort import canonical_assembly_digest, parse_dna_fasta

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data/studies/kaneko_sequence_intake.v1.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def canonical_digest(rows):
    rows = sorted(rows, key=lambda r: json.dumps(r, sort_keys=True))
    return sha256(json.dumps(rows, sort_keys=True, separators=(",", ":"),
                             ensure_ascii=False).encode())


def read_metadata(data, resource):
    reader = csv.DictReader(io.StringIO(data.decode()), delimiter="\t")
    require(len(reader.fieldnames or []) == len(set(reader.fieldnames or [])),
            "Duplicate metadata headers")
    require(set(reader.fieldnames or []) == set(resource["columns"]),
            f"Metadata schema changed: {resource['filename']}")
    rows = list(reader)
    require(len(rows) == resource["rows"] and
            canonical_digest(rows) == resource["canonical_rows_sha256"],
            f"Metadata content changed: {resource['filename']}")
    return rows


def read_file(path, resource):
    data = path.read_bytes()
    require(len(data) == resource["size_bytes"] and sha256(data) == resource["sha256"],
            f"Source bytes changed: {path.name}")
    return data


def paired_read_files(row):
    fields = [row[key].split(";") for key in ("fastq_ftp", "fastq_bytes", "fastq_md5")]
    require(row["library_layout"] == "PAIRED" and {len(x) for x in fields} == {2},
            f"Incomplete paired-file metadata for {row['run_accession']}")
    output = []
    for location, size, checksum in zip(*fields):
        require(location.startswith("ftp.sra.ebi.ac.uk/") and
                location.endswith(".fastq.gz"), "Unexpected read-file location")
        require(size.isdigit() and int(size) > 0 and re.fullmatch(r"[0-9a-f]{32}", checksum),
                "Invalid read size or checksum")
        output.append({"url": "https://" + location, "size_bytes": int(size), "md5": checksum})
    require(len({x["url"] for x in output}) == 2, "Duplicate paired-file location")
    return output


def verify_analysis_xml(data, source, run):
    root = ET.fromstring(data)
    analyses = root.findall("ANALYSIS")
    require(len(analyses) == 1 and analyses[0].get("accession") == source["analysis_accession"],
            "Analysis accession mismatch")
    require([x.get("accession") for x in root.findall(".//RUN_REF")] == [run],
            "Analysis does not reference the expected original run")
    biosamples = [x.text for x in root.findall(".//SAMPLE_REF/IDENTIFIERS/EXTERNAL_ID")
                  if x.get("namespace") == "BioSample"]
    require(biosamples == [source["sample_accession"]], "Analysis BioSample mismatch")
    files = root.findall(".//FILES/FILE")
    require(len(files) == 1 and files[0].get("checksum_method") == "MD5" and
            files[0].get("checksum") == source["compressed_md5"], "Analysis file checksum mismatch")
    return {"program": root.findtext(".//PROGRAM"),
            "declared_coverage": root.findtext(".//COVERAGE"),
            "source_run_match": True}


def fasta_qc(payload, work_path):
    text = payload.decode("ascii")
    contigs = parse_dna_fasta(text, "source assembly")
    lengths = sorted((len(x) for x in contigs), reverse=True)
    total = sum(lengths)
    cumulative = 0
    for size in lengths:
        cumulative += size
        if 2 * cumulative >= total:
            n50 = size
            break
    counts = Counter("".join(contigs))
    called = sum(counts[x] for x in "ACGT")
    work_path.write_bytes(payload)
    return {"assembly_sha256": canonical_assembly_digest(work_path),
            "contigs": len(contigs), "bases": total, "n50": n50,
            "ambiguous_bases": total - called,
            "gc_fraction_called_acgt": (counts["G"] + counts["C"]) / called if called else None}


def fetch_sources(data_root, manifest):
    require(not data_root.exists(), "Fetch directory must be new")
    data_root.mkdir(parents=True)
    resources = list(manifest["metadata"])
    for source in manifest["host_assemblies"]:
        resources.extend([source["fasta"], source["xml"]])
    resources.extend(source["fasta"] for source in manifest["legacy_phages"])
    def fetch(resource):
        path = data_root / resource["filename"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(resource["url"], timeout=60) as response:
            data = response.read()
        path.write_bytes(data)
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(fetch, resources))


def audit(data_root, manifest):
    reference = manifest["development_reference"]
    reference_bytes = (ROOT / reference["path"]).read_bytes()
    require(sha256(reference_bytes) == reference["sha256"], "Development reference changed")
    development = list(csv.DictReader(io.StringIO(reference_bytes.decode())))
    require(Counter(x["entity_kind"] for x in development) ==
            {"host_genome": 177, "phage_genome": 94}, "Development entity scope changed")
    identities = {kind: {} for kind in ("host_genome", "phage_genome")}
    for row in development:
        identities[row["entity_kind"]].setdefault(row["assembly_sha256"], []).append(row["entity_id"])
    prior = manifest["workbook_intake_reference"]
    prior_bytes = (ROOT / prior["path"]).read_bytes()
    require(sha256(prior_bytes) == prior["sha256"], "Workbook intake reference changed")
    prior_manifest = json.loads(prior_bytes)
    metadata, metadata_checks = {}, []
    for resource in manifest["metadata"]:
        data = (data_root / resource["filename"]).read_bytes()
        metadata[resource["filename"]] = read_metadata(data, resource)
        metadata_checks.append({"filename": resource["filename"],
                                "canonical_rows_match": True,
                                "snapshot_bytes_match": sha256(data) == resource["sha256"]})
    runs_by_sample = {}
    read_summary = {}
    for role, expected_count in (("host", 144), ("phage", 38)):
        rows = metadata[f"{role}_runs.tsv"]
        previous = next(r for r in prior_manifest["metadata"] if r["id"] == role)
        core = urllib.parse.parse_qs(urllib.parse.urlparse(previous["url"]).query)["fields"][0].split(",")
        require(canonical_digest([{k: r[k] for k in core} for r in rows]) ==
                previous["canonical_rows_sha256"], "Run identities differ from workbook audit")
        require(len(rows) == expected_count and len({r["sample_accession"] for r in rows}) == expected_count,
                "Run sample scope changed")
        files = [file for row in rows for file in paired_read_files(row)]
        runs_by_sample.update({r["sample_accession"]: r for r in rows})
        read_summary[role] = {"samples": len(rows), "fastq_files": len(files),
                              "compressed_fastq_bytes": sum(x["size_bytes"] for x in files),
                              "instrument_models": dict(Counter(r["instrument_model"] for r in rows))}
    analysis_rows = [r for name, rows in metadata.items()
                     if name.startswith("host_analysis_samples_") for r in rows]
    require(len(analysis_rows) == 143 and len({r["sample_accession"] for r in analysis_rows}) == 143,
            "Host analysis sample scope changed")
    analyses = {r["sample_accession"]: r for r in analysis_rows}
    host_samples = {r["sample_accession"] for r in metadata["host_runs.tsv"]}
    require(host_samples - set(analyses) == {"SAMD00774523"}, "Missing host set changed")
    for name, rows in metadata.items():
        if "_project.tsv" in name or name.startswith(("phage_assembly_samples_",
                "phage_analysis_samples_", "phage_sequence_samples_", "host_sequence_samples_")):
            require(not rows, f"New assembled source located: {name}; review it before reuse")
    qc = []
    with tempfile.TemporaryDirectory(prefix="kaneko_qc_", dir=data_root) as temporary:
        work = Path(temporary) / "sequence.fasta"
        for source in manifest["host_assemblies"]:
            sample = source["sample_accession"]
            row = analyses[sample]
            require(source["analysis_accession"] == row["analysis_accession"] and
                    source["compressed_md5"] == row["submitted_md5"] and
                    source["fasta"]["url"] == "https://" + row["submitted_ftp"], "Host source join changed")
            run = runs_by_sample[sample]
            match = re.fullmatch(r"JMU_No_(\d+)_read", run["sample_title"])
            require(match is not None and source["source_label"] == "{:03d}".format(int(match[1])),
                    "Host phenotype number mismatch")
            provenance = verify_analysis_xml(read_file(data_root / source["xml"]["filename"], source["xml"]),
                                             source, run["run_accession"])
            compressed = read_file(data_root / source["fasta"]["filename"], source["fasta"])
            require(hashlib.md5(compressed).hexdigest() == source["compressed_md5"], "Host MD5 mismatch")
            metrics = fasta_qc(gzip.decompress(compressed), work)
            qc.append({"kind": "host_genome", "source_label": source["source_label"],
                       "sample_accession": sample, "analysis_accession": source["analysis_accession"],
                       **metrics, **provenance})
        legacy = {r["accession"]: r for r in metadata["legacy_sequences.tsv"]}
        for source in manifest["legacy_phages"]:
            record = legacy[source["accession"]]
            require(record["accession"] + "." + record["sequence_version"] == source["accession_version"],
                    "Legacy accession version changed")
            require(record["scientific_name"] == "Escherichia phage phi" + source["source_label"][1:],
                    "Legacy phage name mismatch")
            payload = read_file(data_root / source["fasta"]["filename"], source["fasta"])
            contigs = parse_dna_fasta(payload.decode())
            require(len(contigs) == 1 and source["accession_version"] in payload.decode().splitlines()[0],
                    "Unexpected legacy FASTA record")
            require(hashlib.md5(contigs[0].encode()).hexdigest() == record["sequence_md5"] and
                    len(contigs[0]) == int(record["base_count"]), "Legacy sequence MD5/length mismatch")
            qc.append({"kind": "phage_genome", "source_label": source["source_label"],
                       "accession_version": source["accession_version"], **fasta_qc(payload, work)})
    require(len(qc) == 147 and len({(r["kind"], r["source_label"]) for r in qc}) == 147,
            "Sequence scope or source labels changed")
    for row in qc:
        row["development_exact_matches"] = sorted(identities[row["kind"]].get(row["assembly_sha256"], []))
    duplicates = []
    for kind in ("host_genome", "phage_genome"):
        grouped = {}
        for row in qc:
            if row["kind"] == kind:
                grouped.setdefault(row["assembly_sha256"], []).append(row["source_label"])
        duplicates.extend({"kind": kind, "assembly_sha256": digest, "source_labels": sorted(labels)}
                          for digest, labels in grouped.items() if len(labels) > 1)
    summary = {"schema_version": 1, "review_id": manifest["review_id"],
               "metadata_checks": metadata_checks, "read_inventory": read_summary,
               "host_assemblies_verified": 143, "host_original_run_refs_verified": 143,
               "legacy_phage_sequences_verified": 4,
               "host_missing_assembly": {"source_label": "127", "sample_accession": "SAMD00774523",
                                         "run_accession": "DRR550410"},
               "phages_with_only_read_routes": 38,
               "compressed_host_assembly_bytes": sum(x["fasta"]["size_bytes"] for x in manifest["host_assemblies"]),
               "development_exact_assembly_matches": [r for r in qc if r["development_exact_matches"]],
               "within_panel_exact_assembly_duplicates": duplicates,
               "qc": qc, "status": "sequence_source_and_basic_qc_audit_only",
               "near_duplicate_audit_completed": False, "contamination_and_completeness_assessed": False,
               "training_pairs_admitted": 0, "model_fits": 0, "model_score_runs": 0,
               "independent_release_validation_established": False}
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), "Output file must be new")
    manifest = json.loads(args.manifest.read_text())
    if args.fetch:
        fetch_sources(args.data_root, manifest)
    report = audit(args.data_root, manifest)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        stream.write(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: report[k] for k in ("host_assemblies_verified",
          "legacy_phage_sequences_verified", "training_pairs_admitted")}))


if __name__ == "__main__":
    main()
