# Kaneko prevalence versus host 4-mer results

The frozen comparison **did not establish an improvement from host 4-mer
ranking**. Both methods reached test Hit@2 of 42.86%; host 4-mer used more capped
tests at the point estimate. Keep this result and the original split.

## Primary endpoint and results

The target is **source grade 3, complete spot clearing**. Observed grades 0–2
mean not complete clearing for this target; they do not establish no infection.
All 42 study phages remain candidates. The fixed genomic-group split contains
123 training hosts, 13 validation hosts and seven test hosts in five test groups.
The two test hosts with no grade-3 candidate remain in the primary denominator.

| Method | Validation Hit@2 | Test Hit@2 | Test capped tests at K=2 |
| --- | ---: | ---: | ---: |
| Training prevalence, Beta(1,1) | 76.92% | 42.86% | 1.571 |
| Host 4-mer logistic | 76.92% | 42.86% | 1.714 |

Hit@2 is the chance of finding at least one source-positive candidate in the
first two tests. Capped cost is tests consumed until the first positive or the
budget of two: its expectation is `2 - Hit@1`. Both metrics use exact expectations
for score ties, with uniform permutations inside each tied block.

The precommitted paired whole-group bootstrap uses 2,000 draws and seed 20261114.
Contrasts below are host 4-mer minus prevalence on the primary test endpoint:

| Contrast | Difference | Pointwise 95% interval | Frozen criterion met? |
| --- | ---: | ---: | --- |
| Hit@2 | 0.00 percentage points | −27.27 to +60.00 percentage points | No; lower bound must exceed zero |
| Capped tests at K=2 | +0.143 tests | 0.000 to +0.273 tests | No; upper bound must be below zero |

Neither improvement criterion nor their joint criterion passed. Seven hosts in
five genomic groups give imprecise conditional research intervals. This
within-source comparison does not establish stable superiority, validate the
earlier 94-output weights, or support a release promotion. Scores are
uncalibrated source-target scores, not productive-infection or treatment-benefit
probabilities.

The prespecified secondary grade ≥2 readout reused the same primary-target
scores without refitting: test Hit@2 was 57.14% for both methods, with capped
costs of 1.429 for prevalence and 1.571 for host 4-mer. The complete report also
retains K=1,3,5, positive-eligible recall, MRR and per-phage selection diagnostics.
These descriptive readouts do not replace the primary endpoint.

## What two tests recovered on the inspected splits

A lab that runs both tests can ask how many grade-3 host–phage pairs it would
recover. A lab that stops at the first positive can ask how many hosts it would
cover and how many tests it would consume. The same frozen rankings give
different answers to these questions.

| Partition | Method | Hosts with a hit in the first two | Positive pairs if both tests run | Tests consumed if stopping at first hit |
| --- | --- | ---: | ---: | ---: |
| Validation (13 hosts) | Prevalence | 10/13 | 10/26 | 16 |
| Validation (13 hosts) | Host 4-mer | 10/13 | 17/26 | 18 |
| Test (7 hosts) | Prevalence | 3/7 | 3/14 | 11 |
| Test (7 hosts) | Host 4-mer | 3/7 | 5/14 | 12 |

The test-set tie in host coverage hides one gain and one loss. Host 125 gained
a hit on its second test, while host 082 lost its first-test hit. Hosts 025 and
126 each contributed one extra positive pair when both tests ran. These two
hosts and host 082 share genomic group 025. Host 094 had five grade-3 candidates
in the source panel and both methods missed all five at this budget. Hosts 092
and 135 had no grade-3 candidate among the 42 and remain in the denominator.

To reproduce the pair counts, filter the saved per-host metrics to the primary
grade-3 endpoint. Within each partition and method, sum `positive_candidates`
times `recall@2`, assigning zero contributions to hosts with zero positives.
The totals agree with an independent sum of `expected_positive_top2_selections`
in the per-phage diagnostics. Both calculations retain the frozen score-tie
expectation. Pair counts refer to host–phage pairs, so one phage can count on
more than one host.

This post hoc readout preserves the original acceptance decision. It motivates
a prospective distinction for the next comparison: if the question is how many
positives a fixed screening budget buys, prespecify mean positive pairs per host
after both tests as the primary measure. Report host coverage and first-hit cost
alongside it. A claim about learned representations would also require a matched
annotation baseline and a cohort reserved before model selection.

## Frozen implementation and saved artifacts

The [protocol](kaneko-baseline-protocol.md) and its bound design, input hashes,
groups, seed, endpoint and hyperparameters remain unchanged. The
[runner implementation](../data/studies/kaneko_baseline.implementation.v1.json)
was [committed before any real study fit](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/374c6abba3172488817772169ee9962f2a0cfe53).
The weights and feature cache were then
[committed before held-out scoring](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/f6c796b725a98618e2f8a4bb003bbc317ecc4213).

Training used the 123 training hosts only, including standardization of the 136
canonical host 4-mer features. Each of the 42 fixed `C=0.1`, `lbfgs`,
`max_iter=2000` logistic fits converged: zero convergence warnings and zero
single-class fallbacks. Validation selected nothing; evaluation performed no
fit. One CPU thread was used, with about 174 MiB peak training RSS and no GPU,
paid-provider or new sequence-download requests.

| Artifact | Contents |
| --- | --- |
| [Fit receipt](../data/studies/results/kaneko_baselines.v1/fit/FIT.json) | Training scope, environment, resource use and artifact byte hashes |
| [Host model](../data/studies/results/kaneko_baselines.v1/fit/host_4mer_model.json) and [prevalence model](../data/studies/results/kaneko_baselines.v1/fit/prevalence_model.json) | Frozen weights, train-only preprocessing and endpoint semantics |
| [Feature cache](../data/studies/results/kaneko_baselines.v1/fit/host_features.csv) and [manifest](../data/studies/results/kaneko_baselines.v1/fit/feature_manifest.csv) | All 143 host vectors and source bindings |
| [Report](../data/studies/results/kaneko_baselines.v1/evaluation/report.json) and [contrasts](../data/studies/results/kaneko_baselines.v1/evaluation/contrasts.json) | Both endpoints/partitions, summaries, intervals and acceptance decisions |
| [Candidate scores](../data/studies/results/kaneko_baselines.v1/evaluation/scores.csv) | All 1,680 validation/test candidate scores for both methods |
| [Per-host metrics](../data/studies/results/kaneko_baselines.v1/evaluation/per_host_metrics.csv) and [per-phage diagnostics](../data/studies/results/kaneko_baselines.v1/evaluation/per_phage_diagnostics.csv) | Ranking, cost, eligibility and expected top-two selections |
| [Output bindings](../data/studies/results/kaneko_baselines.v1/evaluation/OUTPUTS.json) | SHA-256 and byte sizes of evaluation outputs |
| [Independent numerical QA](../data/studies/results/kaneko_baselines.v1/NUMERIC_QA.json) | Recomputed features, scores, metrics, intervals and decisions |

The independent verifier imports no project fitting, scoring or metric helpers.
It recomputes all 19,448 feature values directly from source FASTAs, all 1,680
scores from saved coefficients, 1,040 per-host metric values, 44 grouped intervals
and 336 per-phage diagnostic rows. Maximum score disagreement is
`4.44e-16`; no refit is performed. Aggregate report metrics were also checked
against the per-host export. Eight runner contract tests cover train-only
scaling, endpoint/library rejection, single-class support, tie selection,
recall eligibility after group sampling and frozen code/protocol bytes.

## Reproduce

Use the byte-verified [sequence sources](kaneko-sequence-intake.md) and the
`source_spot_grades.csv` from the [workbook audit](kaneko-workbook-intake.md).
Install the exact environment in `fit/FIT.json` (Python 3.12.14, NumPy 2.3.5,
pandas 2.2.3, SciPy 1.17.0, scikit-learn 1.8.0, threadpoolctl 3.6.0).
The supplied fit directory can be evaluated without fitting again:

```bash
PYTHONPATH=src python scripts/run_kaneko_baselines.py evaluate \
  --fit-dir data/studies/results/kaneko_baselines.v1/fit \
  --grades /path/to/source_spot_grades.csv \
  --out-dir artifacts/kaneko_evaluation_new \
  --frozen-fit-commit f6c796b725a98618e2f8a4bb003bbc317ecc4213

python scripts/verify_kaneko_baselines.py \
  --fit-dir data/studies/results/kaneko_baselines.v1/fit \
  --evaluation-dir artifacts/kaneko_evaluation_new \
  --sequence-root /path/to/kaneko_sequences \
  --grades /path/to/source_spot_grades.csv \
  --out artifacts/kaneko_numeric_qa_new.json
```

To reproduce training separately, use a new output directory:

```bash
PYTHONPATH=src python scripts/run_kaneko_baselines.py fit \
  --sequence-root /path/to/kaneko_sequences \
  --grades /path/to/source_spot_grades.csv \
  --out-dir artifacts/kaneko_fit_new \
  --implementation-commit 374c6abba3172488817772169ee9962f2a0cfe53
```

Resource receipts may vary by machine; saved artifact byte bindings are checked
before evaluation. Raw genomes and bulk phenotype grades remain outside Git.
Next, use the saved failures to select one concrete data or representation
hypothesis. Any new fit or selection requires a separate protocol; this test set
has now been inspected.
