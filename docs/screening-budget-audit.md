# Fixed-library screening-budget audit

Status: definition frozen before the new budget curves; execution pending.
This is a descriptive reanalysis of inspected development scores. There is no
new model fit, extraction, provider inference or independent release validation.

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

The implementation and independent verifier will be frozen before execution.
Completed results and reproduction commands will be recorded here.
