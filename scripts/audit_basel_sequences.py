#!/usr/bin/env python3
"""Verify published BASEL accessions, host builds and the joint 94+52 phage audit.

No labels are used to select sequences or thresholds; no prediction is generated.
Exit 2 means the evidence was written but scientific execution remains blocked.
"""
from __future__ import annotations

import argparse
from hashlib import sha1, sha256
import gzip
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from phage_cocktail.artifacts import code_provenance, read_json_object, sha256_file, staged_directory
from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock
from phage_cocktail.basel_sequences import (
    assembly_digest, overlap_summary, read_s5, read_supplement_tables, verify_basel_fasta,
)
from phage_cocktail.cohort import parse_dna_fasta
from phage_cocktail.host_study import verify_resources
from phage_cocktail.paired_splits import parse_paf
import pandas as pd


def run(data_root, out, minimap2, *, fetch=False, pdftotext="pdftotext"):
    out, minimap2 = Path(out).resolve(), Path(minimap2).resolve()
    if out.exists():
        raise ValueError("output directory must be new")
    lock_path = ROOT / "data/studies/basel_ecor_sequences.resources.v1.json"
    lock = read_json_object(lock_path)
    training_path = ROOT / "data/studies/genophi_host_baseline.resources.v1.json"
    training_lock = read_json_object(training_path)
    if sha256_file(training_path) != lock["training_resource_lock_sha256"]:
        raise ValueError("training resource lock changed")
    intake_path = ROOT / "data/studies/basel_ecor_intake.v1.json"
    intake = load_lock(intake_path)
    host_ids = set(lock["host_ids"])
    training_entries = [dict(r, path="training/" + r["path"])
                        for r in training_lock["resources"] if r["kind"] == "phage_genome"
                        or (r["kind"] == "host_genome" and r["id"] in host_ids)]
    paths = verify_resources({"resources": [*lock["resources"], *training_entries]}, data_root, fetch=fetch)
    source = inspect_basel_ecor(paths[("phenotype_source", "BASEL_ECOR")].read_bytes(), intake)
    candidates = [r for r in lock["resources"] if r["kind"] == "basel_phage_genome"]
    published_hosts = [r for r in lock["resources"] if r["kind"] == "published_host_genome"]
    if (set(source["phage_ids"]) != {r["id"] for r in candidates}
            or set(source["host_ids"]) != host_ids
            or host_ids != {r["id"] for r in published_hosts}):
        raise ValueError("sequence rosters differ from pinned phenotype source")
    s5 = read_s5(paths[("supplement", "maffei_s5")].read_bytes())
    pdf = paths[("supplement", "genophi_supplement")]
    extracted = subprocess.check_output([pdftotext, "-layout", str(pdf), "-"], text=True)
    table8, table9 = read_supplement_tables(extracted)
    if set(table9) != host_ids:
        raise ValueError("Supplementary Table 9 host roster differs")
    grouping_path = ROOT / "data/studies/phage_grouping.v1.json"
    if sha256_file(grouping_path) != lock["grouping_lock_sha256"]:
        raise ValueError("predeclared grouping lock changed")
    grouping = read_json_object(grouping_path)
    version = subprocess.check_output([str(minimap2), "--version"], text=True).strip()
    if version != grouping["tool"]["version"]:
        raise ValueError("minimap2 version differs from the frozen grouping rule")

    sequences, digests, phage_rows, conflicts = {}, {}, [], []
    for entry in candidates:
        name = entry["id"]
        row = s5[name]
        if table8[name]["accession"] != row["accession"]:
            raise ValueError(f"published phage accessions conflict: {name}")
        raw = paths[(entry["kind"], name)].read_bytes()
        sequence, receipt = verify_basel_fasta(raw, entry, row)
        sequences[name], digests[name] = sequence, receipt["assembly_sha256"]
        receipt["supplement_table8_genome_size_bp"] = table8[name]["genome_size_bp"]
        if len(sequence) != table8[name]["genome_size_bp"]:
            conflicts.append(dict(phage=name, s5_and_sequence_bp=len(sequence),
                                  genophi_table8_bp=table8[name]["genome_size_bp"]))
        phage_rows.append(receipt)
    training_ids = set()
    for entry in training_entries:
        if entry["kind"] != "phage_genome":
            continue
        name = entry["id"]
        raw = paths[(entry["kind"], name)].read_bytes()
        blob = sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        contigs = parse_dna_fasta(raw.decode("ascii"))
        if blob != entry["git_blob_sha"] or len(contigs) != 1 or raw.splitlines()[0] != b">" + name.encode():
            raise ValueError(f"training phage not the locked single-contig source: {name}")
        training_ids.add(name)
        sequences[name], digests[name] = contigs[0], assembly_digest(contigs)
    if len(training_ids) != 94 or len(candidates) != 52:
        raise ValueError("joint audit must contain exactly 94+52 phages")

    host_rows = []
    for entry in published_hosts:
        name = entry["id"]
        raw = paths[(entry["kind"], name)].read_bytes()
        decoded = gzip.decompress(raw).decode("ascii")
        contigs = parse_dna_fasta(decoded)
        headers = [line for line in decoded.splitlines() if line.startswith(">")]
        master = table9[name]["wgs_master_accession"]
        pattern = r">ENA\|(" + master[:4] + r"\d{8})\|\1\.\d+ .*\b" + re.escape(name) + r"\b"
        if (len(contigs) != table9[name]["contigs"] or len(headers) != len(contigs)
                or not all(re.match(pattern, h) for h in headers)):
            raise ValueError(f"published host FASTA does not match Table 9: {name}")
        previous = paths[("host_genome", name)].read_bytes()
        training_contigs = parse_dna_fasta(previous.decode("ascii"))
        previous_digest, published_digest = assembly_digest(training_contigs), assembly_digest(contigs)
        host_rows.append(dict(
            host=name, wgs_master_accession=master,
            training_raw_sha256=sha256(previous).hexdigest(),
            published_raw_gzip_sha256=sha256(raw).hexdigest(),
            training_assembly_sha256=previous_digest, published_assembly_sha256=published_digest,
            training_contigs=len(training_contigs), published_contigs=len(contigs),
            training_bases=sum(map(len, training_contigs)), published_bases=sum(map(len, contigs)),
            normalized_assemblies_equal=previous_digest == published_digest,
            contig_headers_sha256=sha256("\n".join(headers).encode()).hexdigest(),
        ))

    differences = dict(table8_listed_basel_ids_absent_from_csv=sorted(set(table8) - set(source["phage_ids"])),
                       genome_length_conflicts=conflicts)
    if differences != lock["expected_source_discrepancies"]:
        raise ValueError("source discrepancies changed; review rather than silently remap")
    with staged_directory(out) as stage:
        fasta, paf = stage / "phages.fasta", stage / "alignments.paf"
        fasta.write_text("".join(f">{p}\n{sequences[p]}\n" for p in sorted(sequences)))
        command = [str(minimap2), *grouping["arguments"], str(fasta), str(fasta)]
        with paf.open("w") as output, (stage / "minimap2.log").open("w") as errors:
            subprocess.run(command, stdout=output, stderr=errors, check=True)
        paf.write_text("\n".join(sorted(paf.read_text().splitlines())) + "\n")
        blocks = parse_paf(paf, {p: len(s) for p, s in sequences.items()})
        summaries = []
        for identity in [grouping["primary_identity"], *grouping["sensitivity_identities"]]:
            summary, groups, edges = overlap_summary(
                digests, blocks, source["phage_ids"], training_ids, identity=identity,
                coverage=grouping["minimum_bilateral_coverage"],
            )
            tag = round(100 * identity)
            pd.DataFrame([dict(phage=p, panel="BASEL" if p in source["phage_ids"] else "training",
                               assembly_sha256=digests[p], group=groups[p])
                          for p in sorted(groups)]).to_csv(stage / f"phage_groups_i{tag}.csv", index=False)
            pd.DataFrame(edges, columns=["left", "right"]).to_csv(stage / f"phage_edges_i{tag}.csv", index=False)
            summaries.append(summary)
        pd.DataFrame(phage_rows).to_csv(stage / "basel_accessions.csv", index=False)
        pd.DataFrame(host_rows).to_csv(stage / "host_assemblies.csv", index=False)
        report = dict(
            schema_version=1, status="sequence_audit_complete_execution_blocked_pending_source_and_host_build_review",
            model_scoring_allowed=False, independent_release_validation_established=False,
            full_panel_cold_phage_claim_allowed=False,
            source_counts=source["source_summary"], phages_accession_version_and_sequence_verified=52,
            training_phages_verified=94, published_host_assemblies_verified=25,
            matching_normalized_host_assemblies=sum(r["normalized_assemblies_equal"] for r in host_rows),
            source_discrepancies=differences, phage_overlap=summaries,
            normalization="existing_contig_multiset_order_header_case_and_strand_invariant_not_rotation_invariant",
            limitations=["Host assembly mismatch does not establish a different biological strain.",
                         "No edge under a heuristic genomic threshold is not proof of unrelatedness.",
                         "Canonical exact hashes do not normalize genome rotation or contig segmentation.",
                         "The 52-candidate CSV panel and all missing labels remain unchanged.",
                         "No model is fitted or scored; execution needs host-build and source-discrepancy adjudication.",
                         "Published inspected data cannot justify an untouched release claim."],
            provenance=dict(
                **code_provenance(ROOT, Path(__file__)),
                sequence_resource_lock_sha256=sha256_file(lock_path),
                training_resource_lock_sha256=sha256_file(training_path),
                grouping_lock_sha256=sha256_file(grouping_path), intake_lock_sha256=sha256_file(intake_path),
                minimap2_version=version, minimap2_binary_sha256=sha256_file(minimap2),
                minimap2_source_revision=grouping["tool"]["source_revision"],
                command=["minimap2", *grouping["arguments"], "phages.fasta", "phages.fasta"],
                pdftotext_version=subprocess.run([pdftotext, "-v"], capture_output=True, text=True, check=True).stderr.splitlines()[0],
                python=platform.python_version(),
            ),
        )
        # Complete raw FASTA/PAF remain local; receipts and compact group tables may be committed.
        report["artifacts_sha256"] = {p.name: sha256_file(p) for p in sorted(stage.iterdir())
                                      if p.name != "minimap2.log"}
        (stage / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--minimap2", type=Path, required=True)
    parser.add_argument("--pdftotext", default="pdftotext")
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    report = run(args.data_root, args.out, args.minimap2, fetch=args.fetch, pdftotext=args.pdftotext)
    print("Accessions/sequences:", report["phages_accession_version_and_sequence_verified"], "/ 52")
    print("Same host assembly:", report["matching_normalized_host_assemblies"], "/ 25")
    print("BASEL phages in training-neighbour components:",
          len(report["phage_overlap"][0]["candidates_in_training_components"]))
    print("Status:", report["status"])
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
