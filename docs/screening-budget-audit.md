# Fixed-library screening-budget audit

Status: completed and independently verified.
This is a descriptive reanalysis of inspected development scores. There is no
new model fit, extraction, provider inference or independent release validation.

Protocol freeze: `146fd27ba106cdbfa61db94b0139c69f0094cf2e`.
Executable/verifier freeze: `ef0729017428f65a3529908a9f8b231358d60458`.

## Question

At a fixed test budget, how often does a ranking find a source-positive phage,
and how many tests does it consume when testing stops at the first positive?
The [protocol](../data/studies/screening_budget.v1.json) binds the inputs,
availability rules, metrics and paired host-group uncertainty.

The four source panels are grouped-host seeds 42, 43 and 44, plus B2 once.
Retain validation and test partitions where available. GenoPHI scores cover
all panels; full annotation/4-mer/hybrid exports exist for seeds 43/44.
Other methods contribute their published per-host metric points. Missing
curves remain unavailable rather than being reconstructed by a new fit.

## Metrics and support

- Freeze K = 1, 2, 3, 5, 10 and the same 94 candidates for every query.
- Report exact-tie Hit@K and expected tests consumed together. With first-hit
  rank T, capped consumption is `E[min(T,K)] = K - sum(Hit@j, j=1..K-1)`.
  An all-negative host spends K tests. MRR is a separate metric.
- Primary summaries include every complete-label host. Eligible-only summaries
  are secondary. Exclude an entire incomplete-label host, preserve its ledger,
  and keep its full candidate panel in the score-support check.
- Published Hit@1/2/3/5 points support costs at K=1/2/3, but not K=5/10:
  intermediate probabilities are missing. No interpolation is permitted.
- Use paired frozen-host-group bootstrap intervals: 2,000 resamples,
  seed 20261024, pointwise 95% percentiles. Compare against train-only
  prevalence on the same hosts. Negative cost differences mean fewer tests.

One candidate test is one equal-cost unit. This sequential, perfect-readout
proxy excludes batching, assay noise, turnaround and clinical costs. A useful
ranking needs both adequate hit probability and acceptable test consumption.
An interval describes this fixed development comparison; no new statistical
promotion or median-reduction gate is introduced.

Both denominators draw from the complete-host group universe; the secondary
eligible-only mask is applied after each draw, preserving the existing bootstrap
convention. Reported `groups` counts contributing groups for the selected
denominator. All intervals are pointwise and unadjusted across budgets/methods.

## Observed results

At a cap of three tests, **GenoPHI saves tests in the three grouped-host test
panels, but consumes more on B2**. Its seed-43 Hit@3 is also below prevalence,
despite lower mean test consumption. Cost and hit probability answer different
questions and must be considered together.

All rows include complete all-negative hosts. `Saved` is prevalence consumption
minus GenoPHI consumption; negative values mean extra tests. Intervals are
conditional 95% paired host-group bootstrap intervals.

| Panel | Partition | Complete hosts | Hit@3 prevalence / GenoPHI | Tests prevalence / GenoPHI | Saved tests, 95% interval |
| --- | --- | ---: | ---: | ---: | ---: |
| Seed 42 | Validation | 26 | 88.46% / 80.77% | 1.423 / 1.538 | -0.115 [-0.465, +0.217] |
| Seed 42 | Test | 34 | 85.29% / 91.18% | 1.706 / 1.265 | +0.441 [+0.125, +0.767] |
| B2 | Test | 61 | 90.16% / 90.16% | 1.352 / 1.508 | -0.156 [-0.313, -0.018] |
| Seed 43 | Validation | 24 | 83.33% / 87.50% | 1.542 / 1.500 | +0.042 [-0.273, +0.360] |
| Seed 43 | Test | 27 | 81.48% / 77.78% | 2.000 / 1.648 | +0.352 [+0.043, +0.680] |
| Seed 44 | Validation | 24 | 83.33% / 83.33% | 1.681 / 1.792 | -0.111 [-0.500, +0.304] |
| Seed 44 | Test | 22 | 72.73% / 81.82% | 1.879 / 1.545 | +0.333 [+0.101, +0.587] |

The full [curves](../data/studies/results/screening_budget.v1/curves.csv) retain
K=1/2/3/5/10, both denominators and every available method. The
[contrasts](../data/studies/results/screening_budget.v1/contrasts.csv) retain all
764 comparisons, including unfavorable and inconclusive estimates. Larger
budgets do not uniformly resolve the failures; validation cost differences at
K=3 all have intervals crossing zero.

This budget view does not replace the original
[Hit@2 seed-sensitivity decision](genophi-seed-sensitivity.md), which remains
`not_established`. The new cost intervals are descriptive findings on inspected
data, not a newly established robustness claim. They support evaluating ranking
and test cost together, while preserving the cross-lineage generalization gap.
No new residual model, fallback policy or architecture is selected by this audit.

## Availability and verification

The [availability ledger](../data/studies/results/screening_budget.v1/availability.csv)
distinguishes full scores from published metric points. Full 4-mer, annotation
and hybrid curves exist only at seeds 43/44; seed-42/B2 references retain their
reported points and derivable K=1/2/3 costs. Neighbour controls likewise retain
their published points. Missing quantities remain blank.

The [support ledger](../data/studies/results/screening_budget.v1/support.json)
retains incomplete-host exclusions: ECOR53 on B2, 370D on seed-43 validation and
seed-44 test. Full score exports are checked against the entire 94-phage library
before exclusions; published references retain their original panel contract.

The [run record](../data/studies/results/screening_budget.v1/report.json) binds
25 existing repository inputs and two small source-label tables. Independent
[numeric QA](../data/studies/results/screening_budget.v1/NUMERIC_QA.json) reconstructs
the exact first-positive rank distribution using Fraction arithmetic and expands
group draws directly, independently of the runner's batched implementation.
It verifies all 1,429 host-method rows, 450 curve rows and 764 paired contrasts.
No learned model was refitted and no provider request was made.

## Reproduce

Use the frozen code/settings above. CPU project dependencies suffice; no GenoPHI
package, GPU, model weights or source genomes are needed. `--fetch` downloads only
the two SHA-256/size-locked label tables (about 442 KB). Output paths must be new.

```bash
python scripts/run_screening_budget_audit.py --fetch --outdir artifacts/screening_budget_v1
python scripts/verify_screening_budget_audit.py \
  --run artifacts/screening_budget_v1 --out artifacts/screening_budget_QA.json
```
