"""Synthetic regression checks for source semantics and metadata identity."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "kaneko_audit", Path(__file__).resolve().parents[1] / "scripts/audit_kaneko_workbooks.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def matrix():
    cells = {f"{audit.column(c)}2": f"{c-1:03d}" for c in range(2, 146)}
    for r in range(3, 48):
        cells[f"A{r}"] = ["T2", "T4", "T7"][r-3] if r < 6 else f"ΦWec{100+r}"
        for c in range(2, 146):
            cells[f"{audit.column(c)}{r}"] = 0
    cells["AA6"] = 1
    cells["EO47"] = 3
    return cells


class KanekoAuditTests(unittest.TestCase):
    def test_matrix_keeps_zero_and_intermediate_grades_and_columns_past_z(self):
        hosts, phages, rows = audit.parse_matrix(matrix())
        assert len(rows) == 6480 and len(hosts) == 144 and len(phages) == 45
        assert "T6" not in phages
        assert next(r for r in rows if r["source_cell"] == "AA6")["source_grade"] == 1
        assert next(r for r in rows if r["source_cell"] == "EO47")["source_grade"] == 3
        assert rows[0]["source_grade"] == 0
        assert "infection" not in rows[0] and "binary_label" not in rows[0]


    def test_invalid_or_missing_grade_never_becomes_negative(self):
        for replacement in [None, -1, 4, "0", 1.0, True]:
            with self.subTest(replacement=replacement):
                cells = matrix()
                cells["AA6"] = replacement
                with self.assertRaisesRegex(ValueError, "invalid spot grade"):
                    audit.parse_matrix(cells)


    def test_duplicate_host_is_rejected(self):
        cells = matrix()
        cells["C2"] = "001"
        with self.assertRaisesRegex(ValueError, "Host axis"):
            audit.parse_matrix(cells)


    def test_s5_checks_inclusive_thresholds_without_selecting_a_binary_endpoint(self):
        _, phages, rows = audit.parse_matrix(matrix())
        s5 = {"B2": "score > 1", "C2": "score > 2", "D2": "score > 3"}
        for r, name in enumerate(phages, 3):
            values = [x["source_grade"] for x in rows if x["phage_label"] == name]
            s5[f"A{r}"] = name
            for col, threshold in zip("BCD", (1, 2, 3)):
                s5[f"{col}{r}"] = sum(x >= threshold for x in values)
        assert audit.reconcile_s5(s5, rows)["binary_threshold_selected"] is False
        s5["D47"] += 1
        with self.assertRaisesRegex(ValueError, "aggregate mismatch"):
            audit.reconcile_s5(s5, rows)


    def test_metadata_row_order_allowed_but_identity_changes_fail(self):
        rows = [{"accession": "A", "title": "first"}, {"accession": "B", "title": "second"}]
        resource = {"id": "synthetic", "rows": 2, "canonical_rows_sha256": audit.canonical_digest(rows)}
        data = b"title\taccession\nsecond\tB\nfirst\tA\n"
        assert len(audit.parse_metadata(data, resource)) == 2
        with self.assertRaisesRegex(ValueError, "content changed"):
            audit.parse_metadata(data.replace(b"first", b"unrelated"), resource)


    def test_duplicate_metadata_headers_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate metadata headers"):
            audit.parse_metadata(b"accession\taccession\nA\tA\n", {})

if __name__ == "__main__":
    unittest.main()
