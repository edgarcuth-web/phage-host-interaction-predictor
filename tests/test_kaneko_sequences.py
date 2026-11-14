"""Checks for source-run mismatches and sequence/read integrity semantics."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "kaneko_sequences", Path(__file__).resolve().parents[1] / "scripts/audit_kaneko_sequences.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class SequenceAuditTests(unittest.TestCase):
    def reads(self):
        return {"run_accession": "RUN1", "library_layout": "PAIRED",
                "fastq_ftp": "ftp.sra.ebi.ac.uk/a_1.fastq.gz;ftp.sra.ebi.ac.uk/a_2.fastq.gz",
                "fastq_bytes": "100;200", "fastq_md5": "a" * 32 + ";" + "b" * 32}

    def test_paired_file_budget_counts_both_mates(self):
        rows = audit.paired_read_files(self.reads())
        self.assertEqual(sum(r["size_bytes"] for r in rows), 300)
        self.assertEqual(len(rows), 2)

    def test_incomplete_or_bad_file_metadata_is_rejected(self):
        for field, value in (("fastq_bytes", "100"), ("fastq_md5", "unknown;unknown"),
                             ("fastq_bytes", "100;0"), ("library_layout", "SINGLE")):
            with self.subTest(field=field, value=value):
                row = self.reads()
                row[field] = value
                with self.assertRaises(ValueError):
                    audit.paired_read_files(row)

    def test_duplicate_mate_cannot_count_as_two_files(self):
        row = self.reads()
        row["fastq_ftp"] = ";".join(["ftp.sra.ebi.ac.uk/a.fastq.gz"] * 2)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            audit.paired_read_files(row)

    def test_fasta_metrics_retain_ambiguous_bases(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as temporary:
            qc = audit.fasta_qc(b">one\nACGTNN\n>two\nACGT\n", Path(temporary) / "work.fa")
        self.assertEqual((qc["bases"], qc["contigs"], qc["ambiguous_bases"], qc["n50"]),
                         (10, 2, 2, 6))
        self.assertEqual(qc["gc_fraction_called_acgt"], 0.5)

    def test_source_run_biosample_and_checksum_are_all_required(self):
        source = {"analysis_accession": "ERZ1", "sample_accession": "SAMPLE1",
                  "compressed_md5": "a" * 32}
        data = ('<ANALYSIS_SET><ANALYSIS accession="ERZ1"><RUN_REF accession="RUN1"/>'
                '<SAMPLE_REF><IDENTIFIERS><EXTERNAL_ID namespace="BioSample">SAMPLE1</EXTERNAL_ID>'
                '</IDENTIFIERS></SAMPLE_REF><FILES><FILE checksum_method="MD5" checksum="' +
                "a" * 32 + '"/></FILES></ANALYSIS></ANALYSIS_SET>').encode()
        self.assertTrue(audit.verify_analysis_xml(data, source, "RUN1")["source_run_match"])
        for old, new in ((b"RUN1", b"RUN2"), (b"SAMPLE1", b"SAMPLE2"),
                         (("a" * 32).encode(), ("b" * 32).encode())):
            with self.subTest(change=old):
                with self.assertRaises(ValueError):
                    audit.verify_analysis_xml(data.replace(old, new), source, "RUN1")

    def test_metadata_row_reordering_allowed_but_content_change_rejected(self):
        rows = [{"sample": "A", "run": "X"}, {"sample": "B", "run": "Y"}]
        resource = {"filename": "synthetic", "rows": 2,
                    "columns": ["sample", "run"],
                    "canonical_rows_sha256": audit.canonical_digest(rows)}
        data = b"sample\trun\nB\tY\nA\tX\n"
        self.assertEqual(len(audit.read_metadata(data, resource)), 2)
        with self.assertRaises(ValueError):
            audit.read_metadata(data.replace(b"X", b"Z"), resource)

    def test_error_text_cannot_be_treated_as_an_empty_search(self):
        resource = {"filename": "empty_search", "rows": 0,
                    "columns": ["accession", "sample"],
                    "canonical_rows_sha256": audit.canonical_digest([])}
        self.assertEqual(audit.read_metadata(b"accession\tsample\n", resource), [])
        with self.assertRaisesRegex(ValueError, "schema changed"):
            audit.read_metadata(b"Service unavailable\n", resource)


if __name__ == "__main__":
    unittest.main()
