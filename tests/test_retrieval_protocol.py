"""Exercise study admission failures without fitting or generating predictions."""
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from phage_cocktail.retrieval_protocol import PROTOCOL, check_protocol, validate_inputs

ROOT = Path(__file__).resolve().parents[1]


class RetrievalProtocolTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.plan = json.loads((ROOT / PROTOCOL).read_text())
        names = {v["path"] for v in self.plan["inputs"].values()} | set(self.plan["frozen_code_sha256"]) | {PROTOCOL}
        for name in names:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, path)

    def change_bytes(self, role, data):
        item = self.plan["inputs"][role]
        (self.root / item["path"]).write_bytes(data)
        item.update(sha256=sha256(data).hexdigest(), size_bytes=len(data))

    def table(self, role):
        return list(csv.DictReader(io.StringIO((self.root / self.plan["inputs"][role]["path"]).read_text())))

    def change_csv(self, role, change):
        path = self.root / self.plan["inputs"][role]["path"]
        data = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
        reader = csv.DictReader(io.StringIO(data.decode()))
        rows = list(reader)
        change(rows)
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, reader.fieldnames, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
        data = stream.getvalue().encode()
        self.change_bytes(role, gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)

    def test_real_artifact_admission_without_training(self):
        report = check_protocol(self.root)
        self.assertEqual(report["verified_sequence_vectors"], 271)
        self.assertEqual(report["existing_genophi_scores"]["grouped_host"]["scores"], 5640)
        self.assertEqual(report["cohort"]["phages"], 94)
        self.assertEqual(report["cohort"]["unique_phage_assemblies"], 89)
        self.assertEqual(report["new_model_fits"], 0)
        self.assertEqual(report["new_candidate_scores"], 0)

    def test_protocol_change_is_rejected_before_reading_inputs(self):
        path = self.root / PROTOCOL
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "pre-fit freeze"):
            check_protocol(self.root)

    def test_changed_bound_bytes_are_rejected(self):
        path = self.root / self.plan["inputs"]["library"]["path"]
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "changed input: library"):
            validate_inputs(self.root, self.plan)

    def test_input_path_escape_is_rejected(self):
        self.plan["inputs"]["library"]["path"] = "../outside.csv"
        with self.assertRaisesRegex(ValueError, "escapes repository"):
            validate_inputs(self.root, self.plan)

    def test_host_group_crossing_partitions_is_rejected(self):
        original = self.table("host_manifest")
        groups = {}
        for row in original:
            groups.setdefault(row["host_group"], []).append(row)
        train = next(rows[0]["bacterium"] for rows in groups.values() if len(rows) > 1 and rows[0]["split"] == "train")
        valid = next(r["bacterium"] for r in original if r["split"] == "validation")
        def swap(rows):
            for row in rows:
                if row.get("scenario", "grouped_host") != "grouped_host":
                    continue
                if row["bacterium"] == train:
                    row["split"] = "validation"
                elif row["bacterium"] == valid:
                    row["split"] = "train"
        self.change_csv("host_manifest", swap)
        self.change_csv("scenario_splits", swap)
        with self.assertRaisesRegex(ValueError, "host group leakage"):
            validate_inputs(self.root, self.plan)

    def test_unknown_label_cannot_be_added_as_a_negative(self):
        unknown = self.table("unknown_exclusions")[0]
        new_row = {k: unknown[k] for k in ("bacterium", "source_bacterium", "phage", "source_binary")}
        new_row.update(source_mlc="0", infection="0")
        self.change_csv("observations", lambda rows: rows.append(new_row))
        with self.assertRaisesRegex(ValueError, "partition the fixed library"):
            validate_inputs(self.root, self.plan)

    def test_incomplete_genophi_candidate_panel_is_rejected(self):
        self.change_csv("genophi_primary_scores", lambda rows: rows.pop())
        with self.assertRaisesRegex(ValueError, "GenoPHI score support"):
            validate_inputs(self.root, self.plan)

    def test_nonfinite_genophi_score_is_rejected(self):
        self.change_csv("genophi_b2_scores", lambda rows: rows[0].update(score="nan"))
        with self.assertRaisesRegex(ValueError, "invalid GenoPHI score"):
            validate_inputs(self.root, self.plan)

    def test_wrong_assembly_receipt_is_rejected(self):
        self.change_csv("sequence_manifest", lambda rows: rows[0].update(assembly_sha256="0" * 64))
        with self.assertRaisesRegex(ValueError, "sequence receipt identity"):
            validate_inputs(self.root, self.plan)

    def test_changed_normalized_vector_is_rejected(self):
        def edit(rows):
            rows[0]["AAAA"] = str(float(rows[0]["AAAA"]) + 1e-5)
            rows[0]["AAAC"] = str(float(rows[0]["AAAC"]) - 1e-5)
        self.change_csv("sequence_vectors", edit)
        with self.assertRaisesRegex(ValueError, "sequence vector values"):
            validate_inputs(self.root, self.plan)

    def test_annotation_weights_for_different_hosts_are_rejected(self):
        path = self.root / self.plan["inputs"]["annotation_weights"]["path"]
        weights = json.loads(path.read_text())
        table = self.table("host_manifest")
        weights["training_hosts"][0] = next(r["bacterium"] for r in table if r["split"] == "test")
        self.change_bytes("annotation_weights", (json.dumps(weights) + "\n").encode())
        with self.assertRaisesRegex(ValueError, "incompatible training hosts"):
            validate_inputs(self.root, self.plan)

    def test_b2_test_host_cannot_enter_inner_tuning(self):
        outer = self.table("scenario_splits")
        host = next(r for r in outer if r["scenario"] == "phylogroup:B2" and r["split"] == "test")
        self.change_csv("b2_selection", lambda rows: rows[0].update(bacterium=host["bacterium"], host_group=host["host_group"]))
        with self.assertRaisesRegex(ValueError, "only all outer-training hosts"):
            validate_inputs(self.root, self.plan)

    def test_annotation_weights_from_another_protocol_are_rejected(self):
        path = self.root / self.plan["inputs"]["annotation_weights"]["path"]
        weights = json.loads(path.read_text())
        weights["provenance"]["protocol_sha256"] = "0" * 64
        self.change_bytes("annotation_weights", (json.dumps(weights) + "\n").encode())
        with self.assertRaisesRegex(ValueError, "source protocol/resources differ"):
            validate_inputs(self.root, self.plan)

    def test_cli_refuses_to_overwrite_a_run_record(self):
        output = self.root / "report.json"
        output.write_text("existing record\n")
        result = subprocess.run([sys.executable, str(ROOT / "scripts/check_new_host_retrieval.py"),
                                 "--root", str(self.root), "--out", str(output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("output must be new", result.stderr)
        self.assertEqual(output.read_text(), "existing record\n")


if __name__ == "__main__":
    unittest.main()
