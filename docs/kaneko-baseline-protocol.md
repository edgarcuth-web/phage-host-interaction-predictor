# Kaneko 42-candidate baseline protocol

Status: the host grouping, cohort and comparison protocol are frozen. No
susceptibility model was fit or scored. This is an inspected, single-source
research comparison, not validation of the current 94-output weights.

The [design](../data/studies/kaneko_baseline.design.v1.json) was
[committed before genomic similarity analysis](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/2da042742fca9042f3053d1ae390ac30b6c92028).
The [resolved protocol](../data/studies/results/kaneko_baseline_freeze.v1/protocol.json),
[host ledger](../data/studies/results/kaneko_baseline_freeze.v1/hosts.csv) and
[numeric QA](../data/studies/results/kaneko_baseline_freeze.v1/NUMERIC_QA.json)
bind the actual groups and partitions.

## Host groups and frozen partitions

All 10,153 unordered pairs among the 143 available host assemblies were computed
with [skani 0.4.0](https://github.com/bluenote-1577/skani), using the pinned binary
and options in the design. The complete
[pair ledger](../data/studies/results/kaneko_baseline_freeze.v1/host_pairs.csv)
retains ANI and aligned fractions for both genomes. Missing comparisons fail
verification rather than becoming evidence of dissimilarity.

The primary graph connects pairs with reported ANI ≥99.5% and aligned fractions
≥90% on **both** genomes, plus exact assembly matches. Connected components stay
together. These are operational split groups, not a biological strain definition;
transitive membership does not mean every within-group pair exceeds the cutoff.

The graph has **27 components**, including one with **90 hosts**. Fixed-threshold
descriptive checks give 26 groups at 99.0% and 49 at 99.9%, retaining AF ≥90%.
Those alternatives do not change the primary split or trigger model fits.

| Partition | Hosts | Genomic groups | Observed study-phage cells |
| --- | ---: | ---: | ---: |
| Train | 123 | 17 | 5,166 |
| Validation | 13 | 5 | 546 |
| Test | 7 | 5 | 294 |

Groups were sorted by SHA-256 of `20261114:<group_id>`, then group ID. The first
five groups are test, the next five validation and the rest train. No seed,
threshold or allocation was retried for phenotype balance. All 3,779 qualifying
edges stay within partitions.

**The seven-host, five-group test set has limited precision.** The largest
component falls in training under the precommitted rule. Retain this imbalance
and report conditional research uncertainty; it cannot establish stable superiority.
Near-duplicate comparisons against the prior 177 hosts and contamination/
completeness assessment remain outside this completed grouping step.

## Cohort, endpoint and methods

- Use the 143 source-run-linked assemblies and all 42 named study phages:
  **6,006 measured source cells**. Host 127 remains an explicit genome-availability
  exclusion from every partition; its phenotype rows are preserved in the source.
  T2/T4/T7 controls are outside this candidate library. T6 has no SF2 row.
- Primary target: **source grade 3, complete spot clearing**. Grades 0, 1 and 2
  are observed “not complete clearing” for this target. They are not renamed
  infection-negative results. Retain original grades; reject unknown/invalid cells.
  The workbook captions define this target despite the unreconciled abstract
  percentage and S5 heading discrepancy.
- Secondary readout: grade ≥2, evaluated on the same primary-target scores as a
  descriptive sensitivity. It selects neither a different target nor a refit.
- Compare Beta(1,1)-smoothed training prevalence with one host 4-mer logistic
  model per named phage: 136 existing strand-canonical bins, training-host-only
  standardization, C=0.1, `lbfgs`, `max_iter=2000`, no class weighting.
  Single-class candidates use training prevalence; all 42 remain eligible.
- Fit once on the 123 training hosts. Validation reports the fixed methods;
  it does not select C, features, thresholds or seeds. No train-plus-validation
  refit. The future runner must bind new artifacts to the grade-3 endpoint;
  the historical fitter's hardcoded MLC descriptor is not appropriate here.

The resolved protocol pins the existing mathematical helper blobs and the future
fit environment. Phage genomes are not inputs to these two known-library baselines.
All source scores remain uncalibrated; productive infection and treatment benefit
are different endpoints.

## Evaluation and completion rules

Primary contrast: host 4-mer minus prevalence, all-complete-host macro **Hit@2**
on the frozen test set. Keep the two test hosts without any grade-3 candidate in
the primary denominator. Report K=1,3,5, positive-eligible Recall@K, MRR and capped
test costs as secondary outcomes. Ties use exact uniform-permutation expectations.

Resample whole genomic groups, paired across methods, with 2,000 draws and seed
20261114 for pointwise 95% intervals. Research hit gain needs a lower contrast
bound above zero; cost reduction at K=2 needs an upper difference bound below
zero. Joint benefit needs both. Keep failures and the small-test-set limitation;
no result here promotes a model. Pooled pair AP cannot replace the ranking target.

## Reproduce the freeze

Use the verified [sequence-source directory](kaneko-sequence-intake.md) and the
`source_spot_grades.csv` from the [workbook audit](kaneko-workbook-intake.md):

```bash
python scripts/freeze_kaneko_baseline.py \
  --sequence-root /path/to/kaneko_sequences \
  --grades /path/to/source_spot_grades.csv \
  --work-dir artifacts/kaneko_freeze_work_new \
  --out-dir artifacts/kaneko_freeze_new \
  --pairs data/studies/results/kaneko_baseline_freeze.v1/host_pairs.csv
```

This offline replay verifies source assemblies and reconstructs the frozen
groups, split, targets and protocol. Work/output directories must be new.
To recompute similarities, replace `--pairs` with `--skani /path/to/skani-0.4.0`;
the binary must match the recorded SHA-256. The official download URL is rolling,
so a different binary fails verification. No source genomes or bulk grades are
redistributed in Git.

Next: implement and run the two fixed baselines, verify endpoint-aware artifacts
and full 42-candidate support, then report validation/test results under this protocol.
