"""Synthetic checks for pair completeness, transitive grouping and source targets."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "kaneko_freeze", Path(__file__).resolve().parents[1] / "scripts/freeze_kaneko_baseline.py")
freeze = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(freeze)


class KanekoFreezeTests(unittest.TestCase):
    def pair(self, left, right, ani=99.5, af_left=90.0, af_right=90.0):
        return {"host_left": left, "host_right": right, "ani_percent": ani,
                "af_left_percent": af_left, "af_right_percent": af_right}

    def test_transitive_groups_keep_chained_neighbours_together(self):
        hosts = {"001", "002", "003"}
        pairs = [self.pair("001", "002"), self.pair("002", "003"),
                 self.pair("001", "003", ani=98.0)]
        labels, edges = freeze.make_groups(hosts, pairs, {h: h for h in hosts}, 99.5, 90.0)
        self.assertEqual(set(labels.values()), {"001"})
        self.assertEqual(len(edges), 2)

    def test_high_ani_requires_both_aligned_fractions(self):
        hosts = {"001", "002"}
        labels, _ = freeze.make_groups(hosts, [self.pair("001", "002", ani=100, af_right=89.99)],
                                       {h: h for h in hosts}, 99.5, 90.0)
        self.assertEqual(len(set(labels.values())), 2)

    def test_exact_assembly_match_is_never_split_by_low_alignment(self):
        hosts = {"001", "002"}
        labels, _ = freeze.make_groups(hosts, [self.pair("001", "002", ani=98, af_right=50)],
                                       {h: "same" for h in hosts}, 99.5, 90.0)
        self.assertEqual(labels["001"], labels["002"])

    def test_group_split_is_deterministic_and_never_splits_members(self):
        groups = {f"{i:03d}": f"{i:03d}" for i in range(1, 11)}
        groups["011"] = "001"
        a = freeze.assign_splits(groups, 20261114, 10)
        b = freeze.assign_splits(dict(reversed(list(groups.items()))), 20261114, 10)
        self.assertEqual(a, b)
        self.assertEqual(a["001"], a["011"])
        self.assertEqual(set(a.values()), {"train", "validation", "test"})
        with self.assertRaises(ValueError):
            freeze.assign_splits({"001": "001"}, 20261114, 10)

    def test_partial_clearing_is_not_primary_complete_clearing(self):
        self.assertEqual([freeze.primary_target(str(i)) for i in range(4)], [0, 0, 0, 1])
        for value in ("", "unknown", "4", "0.0", None, 0, True):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    freeze.primary_target(value)

    def test_missing_pair_or_bad_numeric_output_fails_closed(self):
        header = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"]
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as temporary:
            path = Path(temporary) / "edges.tsv"
            def write(ani):
                with path.open("w", newline="") as stream:
                    writer = csv.writer(stream, delimiter="\t")
                    writer.writerow(header)
                    writer.writerow(["host_002.fna", "host_001.fna", ani, 95, 93])
            write(99.6)
            row = freeze.read_pairs(path, {"001", "002"})[0]
            self.assertEqual((row["af_left_percent"], row["af_right_percent"]), (93, 95))
            normalized = Path(temporary) / "pairs.csv"
            freeze.write_csv(normalized, [row])
            self.assertEqual(freeze.read_pairs(normalized, {"001", "002"}), [row])
            with self.assertRaisesRegex(ValueError, "Missing ANI"):
                freeze.read_pairs(path, {"001", "002", "003"})
            write("nan")
            with self.assertRaisesRegex(ValueError, "Invalid ANI"):
                freeze.read_pairs(path, {"001", "002"})

    def test_missing_source_cell_is_not_filled_as_zero(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as temporary:
            path = Path(temporary) / "grades.csv"
            path.write_text("host_label,phage_label,source_grade,phage_role,source_sheet\n001,p,3,study,SF2\n")
            expected = {"sha256": freeze.sha256(path.read_bytes()), "rows": 1}
            with self.assertRaisesRegex(ValueError, "Missing required source cell"):
                freeze.curate_grades(path, {"001"}, ["p", "q"], expected)

    def test_committed_frozen_outputs_match_their_byte_bindings(self):
        root = Path(__file__).resolve().parents[1]
        out = root / "data/studies/results/kaneko_baseline_freeze.v1"
        protocol = json.loads((out/"protocol.json").read_text())
        for name, item in protocol["frozen_outputs"].items():
            data = (out/name).read_bytes()
            self.assertEqual(len(data), item["size_bytes"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), item["sha256"])
            self.assertNotIn(b"\r", data)
        design = root/protocol["design_path"]
        self.assertEqual(hashlib.sha256(design.read_bytes()).hexdigest(), protocol["design_sha256"])


if __name__ == "__main__":
    unittest.main()
