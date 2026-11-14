"""Synthetic checks for the frozen comparison; no real-data model fitting."""
from itertools import permutations
import math
import unittest

import numpy as np
import pandas as pd

from scripts.run_new_host_retrieval import (
    assert_selection_scope, choose_trial, contrast_rows, evaluate_scores, per_phage_rows,
)
from scripts.verify_new_host_retrieval import independent_metrics


class NewHostRetrievalTests(unittest.TestCase):
    def test_hit_gain_has_priority_over_cost(self):
        a = dict(hit2=.6, tests2=1., l2=.1)
        b = dict(hit2=.8, tests2=1.9, l2=.001)
        self.assertIs(choose_trial([a, b]), b)

    def test_cost_breaks_hit_tie(self):
        a = dict(hit2=.8, tests2=1.8, l2=.1)
        b = dict(hit2=.8, tests2=1.4, l2=.001)
        self.assertIs(choose_trial([a, b]), b)

    def test_stronger_regularization_breaks_metric_tie(self):
        a = dict(hit2=.8, tests2=1.4, l2=.1)
        b = dict(hit2=.8+1e-14, tests2=1.4-1e-14, l2=.001)
        self.assertIs(choose_trial([b, a]), a)

    def test_nonfinite_trial_cannot_be_selected(self):
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            choose_trial([dict(hit2=np.nan, tests2=1., l2=.1)])

    def test_outer_test_hosts_cannot_enter_selection(self):
        with self.assertRaisesRegex(ValueError, "outer-test"):
            assert_selection_scope({"a", "b"}, {"c"}, {"b"})
        with self.assertRaisesRegex(ValueError, "overlap"):
            assert_selection_scope({"a", "b"}, {"b"}, {"c"})
        assert_selection_scope({"a", "b"}, {"c"}, {"d"})

    def panel(self, incomplete=False):
        phages = [f"p{i}" for i in range(5)]
        hosts = ["positive", "negative"] + (["unknown"] if incomplete else [])
        scores = pd.MultiIndex.from_product([hosts, phages], names=["bacterium", "phage"]).to_frame(index=False).assign(score=0.)
        observed = scores.drop(columns="score").assign(infection=lambda x: (x.bacterium.eq("positive") & x.phage.eq("p0")).astype(int))
        if incomplete:
            observed = observed[~(observed.bacterium.eq("unknown") & observed.phage.eq("p4"))]
        return phages, hosts, observed, scores

    def test_capped_five_test_cost_uses_intermediate_hit_four(self):
        phages, hosts, observed, scores = self.panel()
        summary, per_host = evaluate_scores(observed, scores, hosts, phages)
        positive = per_host.set_index("bacterium").loc["positive"]
        negative = per_host.set_index("bacterium").loc["negative"]
        self.assertAlmostEqual(positive["hit@4"], .8)
        self.assertAlmostEqual(positive["tests@5"], 3.)
        self.assertAlmostEqual(negative["tests@5"], 5.)
        self.assertAlmostEqual(summary["all_hosts"]["hit@2"], .2)

    def test_unknown_host_is_excluded_whole_after_full_scoring(self):
        phages, hosts, observed, scores = self.panel(incomplete=True)
        summary, per_host = evaluate_scores(observed, scores, hosts, phages)
        self.assertEqual(len(scores), 15)
        self.assertEqual(summary["incomplete_hosts_excluded"], ["unknown"])
        self.assertEqual(len(per_host), 2)
        self.assertIn("negative", per_host.bacterium.to_list())

    def test_partial_candidate_scores_fail_before_metrics(self):
        phages, hosts, observed, scores = self.panel()
        with self.assertRaisesRegex(ValueError, "full fixed"):
            evaluate_scores(observed, scores.iloc[:-1], hosts, phages)

    def test_slot_diagnostics_keep_unknown_selected_outcomes(self):
        phages, hosts, observed, scores = self.panel(incomplete=True)
        rows = per_phage_rows(scores, observed, dict(panel="toy", method="random", partition="test"))
        self.assertAlmostEqual(sum(r["expected_top2_slots"] for r in rows), 6.)
        self.assertAlmostEqual(sum(r["unknown_top2_slots"] for r in rows), .4)
        self.assertEqual(sum(r["unknown"] for r in rows), 1)

    def test_bootstrap_masks_undefined_recall_after_group_sampling(self):
        phages, hosts, observed, scores = self.panel()
        _, metrics = evaluate_scores(observed, scores, hosts, phages)
        names = ("bilinear", "phage_prevalence", "host_4mer", "genophi", "annotation", "additive", "phage_only")
        table = pd.concat([metrics.assign(method=m) for m in names], ignore_index=True)
        plan = dict(uncertainty=dict(resamples=100, seed=20261113, confidence=.95))
        rows = contrast_rows(table, dict(positive="g1", negative="g2"), dict(panel="toy", partition="test"), plan)
        eligible = [r for r in rows if r["denominator"] == "eligible" and r["metric"] == "recall@2"]
        self.assertTrue(all(0 < r["valid_resamples"] < 100 for r in eligible))
        self.assertTrue(all(r["delta"] == r["lower"] == r["upper"] == 0 for r in eligible))
        self.assertTrue(math.isnan(metrics.set_index("bacterium").loc["negative", "recall@2"]))

    def test_fraction_metrics_agree_with_all_tie_permutations(self):
        labels = dict(a=1, b=0, c=1, d=0, e=0)
        scores = dict(a=2., b=2., c=1., d=1., e=1.)
        orders = [list(x)+list(y) for x in permutations("ab") for y in permutations("cde")]
        expected = independent_metrics(labels, scores)
        for k in range(1, 6):
            hit = np.mean([any(labels[p] for p in order[:k]) for order in orders])
            cost = np.mean([min(next(i+1 for i,p in enumerate(order) if labels[p]), k) for order in orders])
            self.assertAlmostEqual(expected[f"hit@{k}"], hit)
            self.assertAlmostEqual(expected[f"tests@{k}"], cost)


if __name__ == "__main__":
    unittest.main()
