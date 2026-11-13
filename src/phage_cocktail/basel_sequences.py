"""Accession reconciliation and sequence-overlap evidence, without model scoring."""
from __future__ import annotations

from hashlib import sha256
from io import BytesIO
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

from .cohort import parse_dna_fasta, reverse_complement
from .paired_splits import neighbour_groups

_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def read_s5(payload: bytes) -> dict:
    """Read the original S5 accession cells, preserving their worksheet row numbers.

    Uses the pinned workbook's existing sheet/columns, not numerical accession
    offsets or a new spreadsheet dependency. Raw bytes are verified by the runner.
    """
    with ZipFile(BytesIO(payload)) as archive:
        strings = ["".join(si.itertext()) for si in
                   ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("s:si", _NS)]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = {}
    for row in sheet.findall("s:sheetData/s:row", _NS):
        cells = {}
        for cell in row.findall("s:c", _NS):
            value = cell.find("s:v", _NS)
            text = "" if value is None else value.text
            if cell.get("t") == "s":
                text = strings[int(text)]
            elif cell.get("t") == "inlineStr":
                text = "".join(cell.find("s:is", _NS).itertext())
            cells[re.sub(r"\d+$", "", cell.attrib["r"])] = text
        name = cells.get("A", "")
        if not re.fullmatch(r"Bas\d{2}", name):
            continue
        accession = cells.get("O", "")
        if name in rows:
            raise ValueError("duplicate S5 phage")
        # S5 also lists the T3(K12) chimera Bas66 with accession n.a.
        # Preserve that cell; only a requested candidate must resolve to an accession.
        rows[name] = dict(s5_row=int(row.attrib["r"]), accession=accession,
                          name=cells["J"], genome_size_bp=int(cells["P"]))
    if not rows:
        raise ValueError("no phage accession rows in S5")
    return rows


def read_supplement_tables(text: str) -> tuple[dict, dict]:
    """Extract accession/length and host assembly fields from pdftotext -layout.

    Matrix/RB-TnSeq column marks are deliberately not inferred from whitespace.
    Their separate manual review is recorded in the source discrepancy notes.
    """
    phage_text = text.split("Supplementary Table 8.", 1)[1].split("Supplementary Table 9.", 1)[0]
    host_text = text.split("Supplementary Table 9.", 1)[1].split("Supplementary Figures", 1)[0]
    phages, hosts = {}, {}
    pattern = r"\b(Bas\d{2})\b(?:(?!\bBas\d{2}\b).)*?\b(MZ\d{6})\b\s+(\d+)"
    for name, accession, size in re.findall(pattern, phage_text, re.S):
        if name in phages:
            raise ValueError("duplicate Supplementary Table 8 phage")
        phages[name] = dict(accession=accession, genome_size_bp=int(size))
    for name, accession, count in re.findall(
        r"\b(ECOR\d+)\s+(Q[A-Z]{3}00000000)\s+\d\.\d Mb\s+(\d+)", host_text
    ):
        if name in hosts:
            raise ValueError("duplicate Supplementary Table 9 host")
        hosts[name] = dict(wgs_master_accession=accession, contigs=int(count))
    if not phages or not hosts:
        raise ValueError("supplement accession tables were not extracted")
    return phages, hosts


def assembly_digest(contigs: list[str]) -> str:
    """Exactly the existing canonical_assembly_digest convention, on parsed DNA.

    Ignores headers, case, contig order and strand. Retains contig multiplicity,
    segmentation and ambiguous bases; does not normalize circular rotation.
    """
    digests = [sha256(min(s, reverse_complement(s)).encode("ascii")).hexdigest() for s in contigs]
    return sha256("\n".join(sorted(digests)).encode("ascii")).hexdigest()


def verify_basel_fasta(payload: bytes, entry: dict, source_row: dict) -> tuple[str, dict]:
    """Require accession.version, source name or explicit strain alias, and length."""
    name, version = entry["id"], entry["accession_version"]
    if (not re.fullmatch(r"MZ\d{6}", source_row["accession"])
            or not re.fullmatch(re.escape(source_row["accession"]) + r"\.\d+", version)):
        raise ValueError(f"accession/version differs from S5: {name}")
    text = payload.decode("ascii")
    contigs = parse_dna_fasta(text)
    headers = [line for line in text.splitlines() if line.startswith(">")]
    if len(contigs) != 1 or len(headers) != 1:
        raise ValueError(f"BASEL phage must have one complete genome record: {name}")
    header = headers[0]
    if not header.startswith(f">ENA|{source_row['accession']}|{version} "):
        raise ValueError(f"FASTA header has wrong accession.version: {name}")
    source_name = source_row["name"] + ", complete genome."
    strain_alias = re.search(r"\bstrain " + re.escape(name) + r", complete genome\.$", header)
    if not (header.endswith(source_name) or strain_alias):
        raise ValueError(f"FASTA phage name/strain alias disagrees: {name}")
    if len(contigs[0]) != source_row["genome_size_bp"]:
        raise ValueError(f"sequence length differs from original S5: {name}")
    return contigs[0], dict(
        phage=name, accession_version=version, s5_row=source_row["s5_row"],
        published_name=source_row["name"], retrieved_header=header,
        genome_size_bp=len(contigs[0]), raw_fasta_sha256=sha256(payload).hexdigest(),
        assembly_sha256=assembly_digest(contigs),
        source_name_in_header=header.endswith(source_name),
    )


def overlap_summary(digests, blocks, candidates, training, *, identity, coverage):
    """Report direct edges and transitive overlap without silently dropping phages."""
    candidates, training = set(candidates), set(training)
    if not candidates or not training or candidates & training or set(digests) != candidates | training:
        raise ValueError("candidate/training rosters must be disjoint and cover the genomes")
    groups, edges = neighbour_groups(digests, blocks, identity=identity, coverage=coverage)
    training_groups = {groups[p] for p in training}
    shared = sorted(p for p in candidates if groups[p] in training_groups)
    cross_edges = [edge for edge in edges if bool(edge[0] in training) != bool(edge[1] in training)]
    exact = sorted((p, t) for p in candidates for t in training if digests[p] == digests[t])
    summary = dict(identity=identity, minimum_bilateral_coverage=coverage,
                   group_count=len(set(groups.values())), edge_count=len(edges),
                   exact_cross_panel_pairs=exact, direct_cross_panel_edges=cross_edges,
                   candidates_in_training_components=shared,
                   candidates_outside_training_components=sorted(candidates - set(shared)))
    return summary, groups, edges
