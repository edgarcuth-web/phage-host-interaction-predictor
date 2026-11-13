"""Identity and overlap controls that do not use experimental outcome labels."""
from hashlib import sha256
from io import BytesIO
from zipfile import ZipFile

import pytest

from phage_cocktail.basel_sequences import (
    assembly_digest, overlap_summary, read_s5, read_supplement_tables, verify_basel_fasta,
)
from phage_cocktail.cohort import canonical_assembly_digest, parse_dna_fasta


def workbook(duplicate=False, unavailable=False):
    # Bas02 maps to an explicitly supplied accession, not a numerical offset.
    row = ('<row r="9"><c r="A9" t="s"><v>0</v></c>'
           '<c r="J9" t="s"><v>1</v></c><c r="O9" t="s"><v>2</v></c>'
           '<c r="P9"><v>4</v></c></row>')
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    out = BytesIO()
    with ZipFile(out, "w") as archive:
        name, accession = ("Bas66", "n.a.") if unavailable else ("Bas02", "MZ501099")
        archive.writestr("xl/sharedStrings.xml", f'<sst {ns}><si><t>{name}</t></si>'
                         f'<si><t>Escherichia phage Example</t></si><si><t>{accession}</t></si></sst>')
        archive.writestr("xl/worksheets/sheet1.xml", f'<worksheet {ns}><sheetData>'
                         + row * (2 if duplicate else 1) + '</sheetData></worksheet>')
    return out.getvalue()


def test_accession_is_read_from_source_cells_not_bas_number():
    assert read_s5(workbook())["Bas02"] == dict(
        s5_row=9, accession="MZ501099", name="Escherichia phage Example", genome_size_bp=4,
    )
    with pytest.raises(ValueError, match="duplicate"):
        read_s5(workbook(duplicate=True))


def test_unavailable_nonpanel_accession_is_preserved_but_cannot_admit_a_sequence():
    # The real S5 includes Bas66 with accession n.a.; it is outside our frozen panel.
    row = read_s5(workbook(unavailable=True))["Bas66"]
    assert row["accession"] == "n.a."
    with pytest.raises(ValueError, match="accession/version"):
        verify_basel_fasta(b">anything\nACGT\n", dict(id="Bas66", accession_version="MZ501099.1"), row)


def test_pdf_tables_preserve_length_conflicts_and_host_builds():
    text = ("Supplementary Table 8. Phages\nBas02 Escherichia phage Example MZ501099 3 x\n"
            "Supplementary Table 9. Assemblies\nECOR12 QOWX00000000 5.1 Mb 218 128.9 kb\n"
            "Supplementary Figures")
    phages, hosts = read_supplement_tables(text)
    assert phages["Bas02"]["genome_size_bp"] == 3  # does not silently replace with S5's 4
    assert hosts["ECOR12"] == dict(wgs_master_accession="QOWX00000000", contigs=218)


@pytest.mark.parametrize("title", ["Escherichia phage Example", "Renamedvirus example strain Bas02"])
def test_versioned_sequence_with_original_name_or_explicit_strain_alias(title):
    row = read_s5(workbook())["Bas02"]
    raw = f">ENA|MZ501099|MZ501099.1 {title}, complete genome.\nACGT\n".encode()
    sequence, receipt = verify_basel_fasta(raw, dict(id="Bas02", accession_version="MZ501099.1"), row)
    assert sequence == "ACGT"
    assert receipt["raw_fasta_sha256"] == sha256(raw).hexdigest()
    assert receipt["source_name_in_header"] == (title == row["name"])


@pytest.mark.parametrize("change,match", [
    ("version", "accession.version"), ("alias", "name/strain"),
    ("length", "length"), ("contig", "one complete"), ("iupac", "invalid"),
])
def test_inconsistent_sequence_identity_fails_closed(change, match):
    raw = b">ENA|MZ501099|MZ501099.1 Escherichia phage Example, complete genome.\nACGT\n"
    if change == "version":
        raw = raw.replace(b"MZ501099.1", b"MZ501099.2")
    elif change == "alias":
        raw = raw.replace(b"Escherichia phage Example", b"Other virus strain Bas03")
    elif change == "length":
        raw = raw.replace(b"ACGT", b"ACG")
    elif change == "contig":
        raw += b">second\nA\n"
    else:
        raw = raw.replace(b"ACGT", b"ACGZ")
    with pytest.raises(ValueError, match=match):
        verify_basel_fasta(raw, dict(id="Bas02", accession_version="MZ501099.1"), read_s5(workbook())["Bas02"])


def test_normalization_matches_existing_convention_and_retains_contig_multiplicity(tmp_path):
    first = tmp_path / "first.fasta"
    first.write_text(">a\nacgtt\n>b\nNNAC\n")
    normalized = assembly_digest(parse_dna_fasta(first.read_text()))
    assert normalized == canonical_assembly_digest(first)
    assert normalized == assembly_digest(["GTNN", "AACGT"])
    assert normalized != assembly_digest(["GTNN", "AACGT", "AACGT"])
    assert normalized != assembly_digest(["ACGTTNNAC"])


def test_training_overlap_uses_transitive_components_and_reports_exact_duplicates():
    digests = {"train": "a", "Bas01": "b", "Bas02": "c", "Bas03": "a"}
    blocks = [dict(query=a, target=b, qlen=100, tlen=100, qs=0, qe=90,
                   ts=0, te=90, matches=90, block=90, strand="+")
              for a, b in [("train", "Bas01"), ("Bas01", "Bas02")]]
    result, groups, edges = overlap_summary(
        digests, blocks, ["Bas01", "Bas02", "Bas03"], ["train"], identity=.95, coverage=.85,
    )
    assert result["candidates_in_training_components"] == ["Bas01", "Bas02", "Bas03"]
    assert result["exact_cross_panel_pairs"] == [("Bas03", "train")]
    assert len(result["direct_cross_panel_edges"]) == 2  # Bas02 is connected indirectly
    assert len(set(groups.values())) == 1


def test_incomplete_or_colliding_roster_cannot_clear_overlap():
    with pytest.raises(ValueError, match="rosters"):
        overlap_summary({"Bas01": "a", "train": "b"}, [], ["Bas01"], ["Bas01"], identity=.95, coverage=.85)
