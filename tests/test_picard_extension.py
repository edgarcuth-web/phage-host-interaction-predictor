"""Identity and leakage-boundary checks for the same-source extension."""
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "studies/picard_extension/audit_cohort.py"
SPEC = importlib.util.spec_from_file_location("picard_audit", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def distances(names, default=.1):
    return {(a, b): 0. if a == b else default for a in names for b in names}


class PicardExtensionTests(unittest.TestCase):
    def profile(self, text):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "host.fasta"
            p.write_text(text)
            source = {"size": p.stat().st_size, "computed_md5": hashlib.md5(p.read_bytes()).hexdigest(), "name": "host.fasta"}
            return audit.fasta_profile(p, source)

    def test_rejects_alias_collisions(self):
        with self.assertRaises(ValueError):
            audit.unique([{"host": "ECOR-1"}, {"host": "ecor1"}], "host")

    def test_canonical_content_ignores_header_order_and_strand(self):
        a = self.profile(">one\nACGTR\n>two\nGGA\n")
        b = self.profile(">other\nTCC\n>renamed\nYACGT\n")
        self.assertEqual(a["assembly_sha256"], b["assembly_sha256"])
        self.assertNotEqual(a["raw_sha256"], b["raw_sha256"])
        self.assertEqual(a["total_bases"], 8)
        self.assertEqual(a["ambiguous_bases"], 1)

    def test_different_fragmentation_remains_distinct(self):
        self.assertNotEqual(self.profile(">a\nACGT\n")["assembly_sha256"],
                            self.profile(">a\nAC\n>b\nGT\n")["assembly_sha256"])

    def test_rejects_malformed_or_non_dna_fasta(self):
        for value in ("ACGT\n", ">empty\n", ">protein\nACEF\n"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.profile(value)

    def test_genome_digest_is_checked_before_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "host.fasta"
            p.write_text(">a\nACGT\n")
            with self.assertRaises(ValueError):
                audit.fasta_profile(p, {"size": p.stat().st_size, "computed_md5": "0" * 32, "name": "host.fasta"})

    def test_distance_only_bridge_prevents_false_new_group(self):
        d = distances(["old", "new", "unmeasured"])
        d["old", "unmeasured"] = d["unmeasured", "old"] = .00005
        d["new", "unmeasured"] = d["unmeasured", "new"] = .00005
        components, _ = audit.group_records(["old", "new", "unmeasured"], d,
                                            {"old": "a", "new": "b"}, .0001)
        self.assertEqual(components.find("new"), components.find("old"))

    def test_threshold_is_strict(self):
        d = distances(["a", "b"], .0001)
        components, _ = audit.group_records(["a", "b"], d, {"a": "x", "b": "y"}, .0001)
        self.assertNotEqual(components.find("a"), components.find("b"))

    def test_exact_assembly_match_connects_distant_core_records(self):
        components, edges = audit.group_records(["a", "b"], distances(["a", "b"]),
                                               {"a": "identical", "b": "identical"}, .0001)
        self.assertEqual(components.find("a"), components.find("b"))
        self.assertEqual(edges[0]["kind"], "exact_canonical_assembly")

    def test_invalid_core_distances_are_rejected(self):
        for wrong in (float("nan"), -.1, .2):
            d = distances(["a", "b"])
            d["a", "b"] = wrong
            with self.subTest(wrong=wrong), self.assertRaises(ValueError):
                audit.group_records(["a", "b"], d, {}, .0001)


if __name__ == "__main__":
    unittest.main()
