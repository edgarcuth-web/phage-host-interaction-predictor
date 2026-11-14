"""Endpoint, preprocessing, fallback and grouped-metric contracts for the runner."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from phage_cocktail.kaneko_baseline import (
    ENDPOINT, FEATURE_COLUMNS, PROTOCOL_PATH, PROTOCOL_SHA256, SCORE_SEMANTICS,
    check_implementation, digest, fit_grade3, metric_contrast, predict_grade3,
    selection_mass,
)


class KanekoBaselineTests(unittest.TestCase):
    def features(self):
        data = np.zeros((3,136))
        data[:,0] = [0.1,0.3,999.0]
        return pd.DataFrame(data,index=["a","b","held"],columns=FEATURE_COLUMNS)

    def training(self, single_class=False):
        rows=[]
        for host in ["a","b"]:
            for phage in ["p","q"]:
                grade = (0 if phage=="p" else 3) if single_class else (3 if (host=="a")== (phage=="p") else 0)
                rows.append({"bacterium":host,"phage":phage,"source_grade":grade,"infection":int(grade==3)})
        return pd.DataFrame(rows)

    def test_scaling_uses_only_training_hosts_and_model_round_trips(self):
        features=self.features()
        model=fit_grade3(self.training(),features,["p","q"])
        self.assertAlmostEqual(model["preprocessing"]["means"][0],0.2)
        self.assertAlmostEqual(model["preprocessing"]["scales"][0],0.1)
        self.assertEqual(model["training_hosts"],["a","b"])
        self.assertEqual(model["endpoint_id"],ENDPOINT)
        self.assertEqual(model["score_semantics"],SCORE_SEMANTICS)
        pairs=pd.DataFrame({"bacterium":["held","held"],"phage":["p","q"]})
        left=predict_grade3(model,pairs,features,["p","q"])
        right=predict_grade3(json.loads(json.dumps(model)),pairs,features,["p","q"])
        np.testing.assert_array_equal(left,right)

    def test_single_class_fallback_retains_both_candidates(self):
        model=fit_grade3(self.training(True),self.features(),["p","q"])
        self.assertEqual(set(model["phages"]),{"p","q"})
        self.assertEqual(model["phages"]["p"]["score"],0.25)
        self.assertEqual(model["phages"]["q"]["score"],0.75)

    def test_wrong_target_and_missing_training_candidate_fail(self):
        train=self.training()
        train.loc[0,"source_grade"]=2
        with self.assertRaisesRegex(ValueError,"not source grade 3"):
            fit_grade3(train,self.features(),["p","q"])
        with self.assertRaisesRegex(ValueError,"Incomplete training library"):
            fit_grade3(self.training().query('phage == "p"'),self.features(),["p","q"])

    def test_historical_endpoint_and_incomplete_model_are_rejected(self):
        features=self.features()
        model=fit_grade3(self.training(True),features,["p","q"])
        pairs=pd.DataFrame({"bacterium":["held","held"],"phage":["p","q"]})
        legacy=deepcopy(model)
        legacy["score_semantics"]="uncalibrated_source_mlc_positive"
        with self.assertRaisesRegex(ValueError,"Wrong model endpoint"):
            predict_grade3(legacy,pairs,features,["p","q"])
        incomplete=deepcopy(model)
        del incomplete["phages"]["q"]
        with self.assertRaisesRegex(ValueError,"Wrong fitted candidate library"):
            predict_grade3(incomplete,pairs,features,["p","q"])

    def test_exact_tie_selection_mass_is_not_lexicographic(self):
        np.testing.assert_array_equal(selection_mass([1,1,1,1],2),[0.5]*4)
        np.testing.assert_array_equal(selection_mass([3,2,2,0],2),[1,0.5,0.5,0])
        self.assertEqual(float(selection_mass([3,2,2,0],2).sum()),2.0)

    def test_ineligible_recall_is_masked_after_group_sampling(self):
        left=pd.DataFrame({"bacterium":["a","b"],"positive_candidates":[1,0],
              "eligible":[True,False],"recall@2":[0.5,np.nan]})
        right=left.copy()
        right.loc[0,"recall@2"]=0.25
        result=metric_contrast(left,right,{"a":"g1","b":"g2"},"recall@2",
                 {"resamples":100,"seed":123,"confidence":0.95})
        self.assertEqual(result["groups"],2)
        self.assertGreater(result["undefined_resamples"],0)
        self.assertEqual(result["interval"],[0.25,0.25])
        self.assertTrue(pd.isna(left.loc[1,"recall@2"]))

    def test_implementation_code_mutation_fails(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as temporary:
            root=Path(temporary)
            (root/"src").mkdir()
            path=root/"src/runner.py"
            path.write_text("frozen\n")
            receipt=root/"data/studies/kaneko_baseline.implementation.v1.json"
            receipt.parent.mkdir(parents=True)
            receipt.write_text(json.dumps({"protocol_sha256":PROTOCOL_SHA256,
                 "code_sha256":{"src/runner.py":digest(path.read_bytes())}}))
            check_implementation(root)
            path.write_text("changed\n")
            with self.assertRaisesRegex(ValueError,"Implementation changed"):
                check_implementation(root)

    def test_resolved_protocol_bytes_are_frozen(self):
        root=Path(__file__).resolve().parents[1]
        self.assertEqual(digest((root/PROTOCOL_PATH).read_bytes()),PROTOCOL_SHA256)


if __name__ == "__main__":
    unittest.main()
