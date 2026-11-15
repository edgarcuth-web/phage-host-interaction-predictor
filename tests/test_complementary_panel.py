"""Safety and ranking semantics for the exploratory fixed-panel comparison."""
import importlib.util
import itertools
import unittest
from fractions import Fraction
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "studies/complementary_panel/run_comparison.py"
SPEC = importlib.util.spec_from_file_location("complementary_panel", SCRIPT)
study = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(study)


class ComplementaryPanelTests(unittest.TestCase):
    def test_selects_complement_instead_of_redundant_second(self):
        labels = {(h, p): int(h in positives) for p, positives in
                  {"a": "123", "b": "123", "c": "4"}.items() for h in "1234"}
        panel, trace, _ = study.select_panel(list("1234"), list("abc"), labels)
        self.assertEqual(panel, ["a", "c"])
        self.assertEqual([r["marginal_training_hosts"] for r in trace], [3, 1])

    def test_selector_cannot_read_target_labels(self):
        class TrainingOnly(dict):
            def __getitem__(self, key):
                if key[0] == "unseen":
                    raise AssertionError("target label accessed")
                return super().__getitem__(key)
        labels = TrainingOnly({("train", "a"): 1, ("train", "b"): 0,
                               ("unseen", "a"): 0, ("unseen", "b"): 1})
        self.assertEqual(study.select_panel(["train"], ["a", "b"], labels)[0], ["a", "b"])

    def test_unknowns_not_added_to_prevalence_denominator(self):
        labels = {("h1", "a"): 1, ("h2", "a"): None,
                  ("h1", "b"): 0, ("h2", "b"): 1}
        stats = study.candidate_stats(["h1", "h2"], ["a", "b"], labels)
        self.assertEqual(stats["a"]["rate"], Fraction(2, 3))
        self.assertEqual(stats["a"]["observed"], 1)
        self.assertEqual(stats["b"]["rate"], Fraction(1, 2))

    def test_training_group_multiplicities_affect_selection(self):
        labels = {(h, p): int((h == "h1") == (p == "a"))
                  for h in ("h1", "h2") for p in ("a", "b")}
        panel, trace, _ = study.select_panel(["h1", "h2"], ["a", "b"], labels,
                                           weights={"h1": 1, "h2": 3})
        self.assertEqual(panel, ["b", "a"])
        self.assertEqual([r["marginal_training_hosts"] for r in trace], [3, 1])
        with self.assertRaises(ValueError):
            study.select_panel(["h1", "h2"], ["a", "b"], labels, weights={"target": 1})

    def test_exact_ties_match_exhaustive_permutations(self):
        scores, labels = {"a": 1., "b": 1., "c": 1., "d": 0.}, {"a": 1, "b": 0, "c": 0, "d": 1}
        for k in (1, 2):
            expected = Fraction(sum(any(labels[p] for p in order[:k])
                                    for order in itertools.permutations("abc")), 6)
            self.assertEqual(study.hit_probability(scores, labels, k), expected)

    def test_incomplete_outcome_bounds_and_known_hit(self):
        scores, labels = {"a": 2., "b": 1., "c": 0.}, {"a": 0, "b": None, "c": 1}
        self.assertEqual(study.hit_probability(scores, labels, 2), 0)
        self.assertEqual(study.hit_probability(scores, labels, 2, unknown_positive=True), 1)
        labels["a"] = 1
        self.assertEqual(study.hit_probability(scores, labels, 2), 1)
        self.assertEqual(study.hit_probability(scores, labels, 2, unknown_positive=True), 1)

    def test_tie_changes_test_consumption_without_changing_budget(self):
        scores, labels = {"a": 1., "b": 1., "c": 0.}, {"a": 1, "b": 0, "c": 0}
        self.assertEqual(study.hit_probability(scores, labels, 2), 1)
        self.assertEqual(2 - study.hit_probability(scores, labels, 1), Fraction(3, 2))


if __name__ == "__main__":
    unittest.main()
