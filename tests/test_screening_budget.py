"""Budget stopping, ties, eligibility, grouped uncertainty and immutable inputs."""
from hashlib import sha256
import importlib.util
from itertools import permutations
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.retrieval import evaluate_retrieval, paired_retrieval_bootstrap

REPO = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


audit = load_script("run_screening_budget_audit", "run_screening_budget_audit.py")
qa = load_script("verify_screening_budget_audit", "verify_screening_budget_audit.py")


def test_exact_rank_distribution_matches_every_permutation_inside_a_tie():
    labels = [1, 0, 1, 0]
    ranks = [next(i+1 for i, p in enumerate(order) if labels[p])
             for order in permutations(range(4))]
    result = qa.first_hit_metrics(labels, [1.]*4, [1, 2, 3, 4])
    assert result["mrr"] == pytest.approx(np.mean([1/r for r in ranks]))
    for k in range(1, 5):
        assert result[f"hit@{k}"] == pytest.approx(np.mean([r <= k for r in ranks]))
        assert result[f"tests@{k}"] == pytest.approx(np.mean([min(r, k) for r in ranks]))
    negative = qa.first_hit_metrics([0]*4, [1.]*4, [1, 3])
    assert negative == {"mrr": 0., "hit@1": 0., "tests@1": 1., "hit@3": 0., "tests@3": 3.}


def test_stopping_cost_keeps_negative_hosts_and_excludes_incomplete_hosts_whole():
    phages, hosts = list("PQRS"), list("ABC")
    observed = pd.DataFrame([(h, p, int(h == "A" and p == "P"))
        for h in hosts for p in phages if (h, p) != ("C", "S")],
        columns=["bacterium", "phage", "infection"])
    scores = pd.DataFrame([(h, p, float(p in "PQ")) for h in hosts for p in phages],
                          columns=["bacterium", "phage", "score"])
    summary, frame = evaluate_retrieval(observed, scores, hosts=hosts, candidate_phages=phages,
                                        ks=(1, 2, 3, 4))
    frame = audit.budget_metrics(frame, [1, 2, 3]).set_index("bacterium")
    assert summary["incomplete_hosts_excluded"] == ["C"]
    assert set(frame.index) == {"A", "B"}
    assert frame.loc["A", "tests@2"] == 1.5
    assert frame.loc["B", "tests@3"] == 3
    assert frame["tests@3"].mean() == 2.25
    assert frame[frame.eligible]["tests@3"].mean() == 1.5
    assert frame.loc["B", "hit@3"] == 0
    with pytest.raises(ValueError, match="exact finite"):
        audit.check_scores(scores.iloc[:-1], hosts, phages)


def test_sparse_historical_points_never_invent_intermediate_probabilities():
    points = pd.DataFrame({"hit@1": [.2], "hit@2": [.4], "hit@3": [.6], "hit@5": [.8], "mrr": [.4]})
    result = audit.budget_metrics(points, [1, 2, 3, 5, 10])
    assert result["tests@1"].iloc[0] == 1
    assert result["tests@3"].iloc[0] == pytest.approx(2.4)
    assert result["tests@5"].isna().all() and result["tests@10"].isna().all()
    assert result["hit@10"].isna().all() and result["hit@5"].iloc[0] == .8


@pytest.mark.parametrize("eligible_only", [False, True])
def test_batched_intervals_match_existing_host_group_bootstrap(eligible_only):
    left = pd.DataFrame({"bacterium": list("ABCDEF"), "eligible": [True, True, True, False, True, True],
                        "positive_candidates": [2, 3, 1, 0, 2, 1],
                        "hit@2": [1., .5, 0., 0., 1., .25], "tests@3": [1., 1.5, 3., 3., 2., 2.5]})
    right = left.assign(**{"hit@2": [.5, 0., 1., 0., 1., .5], "tests@3": [2., 3., 1., 3., 1., 2.]})
    groups = dict(zip("ABCDEF", ["g1", "g1", "g2", "g3", "g4", "g4"]))
    kwargs = {"resamples": 97, "seed": 17, "confidence": .95, "eligible_only": eligible_only}
    batch = audit.paired_intervals(left, right, groups, ["hit@2", "tests@3"], **kwargs)
    for metric in batch:
        original = paired_retrieval_bootstrap(left, right, groups, metric=metric, **kwargs)
        assert batch[metric]["delta"] == pytest.approx(original["delta"], abs=1e-14)
        np.testing.assert_allclose(batch[metric]["interval"], original["interval"], rtol=0, atol=1e-14)
        assert batch[metric]["valid_resamples"] == original["valid_resamples"]
    with pytest.raises(ValueError, match="label support"):
        audit.paired_intervals(left, right.assign(eligible=True), groups, ["hit@2"], **kwargs)


def test_frozen_inputs_and_historical_support_fail_on_drift(tmp_path):
    path = tmp_path / "input.txt"
    path.write_bytes(b"frozen")
    hashes = {"input.txt": sha256(b"frozen").hexdigest()}
    audit.verify_hashes(tmp_path, hashes)
    path.write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed frozen"):
        audit.verify_hashes(tmp_path, hashes)
    with pytest.raises(ValueError, match="frozen input"):
        audit.verify_hashes(tmp_path, {"../outside": "0"*64})
    reference = pd.DataFrame({"bacterium": ["A"], "positive_candidates": [1], "eligible": [True],
                              **{m: [.5] for m in audit.HISTORICAL}})
    audit.check_reference(reference, reference.copy())
    with pytest.raises(ValueError, match="metric differs"):
        audit.check_reference(reference.assign(**{"hit@2": .6}), reference)
    with pytest.raises(ValueError, match="support"):
        audit.check_reference(reference.iloc[:0], reference)


def test_committed_budget_bundle_is_bound_to_independent_qa():
    root = REPO / "data/studies/results/screening_budget.v1"
    if not root.exists():
        pytest.skip("execution follows implementation freeze")
    report = json.loads((root / "report.json").read_text())
    verification = json.loads((root / "NUMERIC_QA.json").read_text())
    assert verification["status"] == "verified"
    assert sha256((root / "report.json").read_bytes()).hexdigest() == verification["report_sha256"]
    assert verification["per_host_rows_independently_checked"] == report["per_host_rows"]
    assert verification["contrasts_independently_checked"] == report["contrast_rows"]
    assert report["new_learned_model_fits"] == report["new_provider_requests"] == 0
    for name, digest in report["output_sha256"].items():
        assert sha256((root / name).read_bytes()).hexdigest() == digest
    freeze = json.loads((REPO / audit.IMPLEMENTATION).read_text())
    assert sha256((REPO / audit.PROTOCOL).read_bytes()).hexdigest() == freeze["protocol_sha256"] == report["protocol_sha256"]
    assert sha256((REPO / audit.IMPLEMENTATION).read_bytes()).hexdigest() == report["implementation_sha256"]
    assert verification["verification_script_sha256"] == freeze["files_sha256"]["scripts/verify_screening_budget_audit.py"]
