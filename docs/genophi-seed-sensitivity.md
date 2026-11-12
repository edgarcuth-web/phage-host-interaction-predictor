# GenoPHI grouped-host split sensitivity: seeds 43 and 44

**The registered robustness criterion is not met (`not_established`).** GenoPHI's
test Hit@2 point advantage over train-only phage prevalence stays positive in
both new seeds, but both paired group-bootstrap intervals cross zero and both
Holm-adjusted group sign-flip p-values exceed 0.05. This is failure to establish
the registered statistical robustness claim, not a directional reversal or
proof that the methods are equal.

Prevalence's validation-to-test fall is also present under both sensitivities:
**-16.20 points at seed 43 and -9.85 at seed 44**, beside seed 42's -11.09.
Seed 42 is therefore not the only split with a substantial prevalence swing.
The exact registered “at least 10 points in both seeds” rule is nevertheless
not met: seed 44 is 0.15 percentage points short. The registered classification
is **`mixed_or_different`**, with the threshold preserved.

## Registration and scope

The [decision rule](../data/studies/genophi_seed_sensitivity.v1.json) was committed
and pushed in `81c4a969ce5f6c2c74a924839c7c8e8c3663d6c1` before sensitivity fitting.
The manifests/protocols were frozen in `9024086b7a948796c8feadf52a997cb8f3ee7304`;
execution/statistical code was frozen in
`c215d86a497a3135322722676b80a63ebb7225b2`. The original seed-42 protocols and
result bundles retain their original bytes.

Two clarifications to the task brief were settled before running:

- Seeds 43/44 are taken from `paired_benchmark.v2.json`'s prespecified
  `splitting.sensitivity_seeds`. GenoPHI's direct host-study parent specified
  only seed 42, so this is a newly registered sensitivity study using those
  existing seed choices.
- A host-level sign-flip test reduces to exact McNemar for binary pairs.
  **Group-level** flips have that equivalence when discordant hosts occur in
  distinct contributing groups; binary outcomes alone do not guarantee it.

The [request audit](../data/studies/genophi_seed_sensitivity.audit.v1.json)
reproduces every quoted seed-42 count, median, fractional-value count and
McNemar contrast. GenoPHI changes +14.2534 points from validation to test while
prevalence changes -11.0860; the *contrast* changes 25.3394 points.

Only `grouped_host` is refitted. Reading `retrieval_splits` and comparing its B2
output under all three manifests verifies that B2 depends on phylogroup and
frozen host groups, not the seed. **B2 is not refitted by either study.**
Clustering, search, assignment and the full protein-family feature table are
split-independent: their v1 bytes are imported under the
[feature reuse lock](../data/studies/genophi_features_reuse.v1.json), rather than
recomputed. Both new runs record zero new clustering/search/assignment calls.

### Registered decision

Both seeds must complete and pass QA. For the primary all-complete-host **test
Hit@2, GenoPHI minus prevalence**, *each* new seed must have:

1. a strictly positive point difference;
2. a strictly positive lower endpoint of the original 95% paired host-group
   bootstrap interval; and
3. two-sided group sign-flip p <= 0.05 after Holm adjustment across the two
   new primary contrasts.

Any nonpositive new point difference would be a directional counterexample.
Positive differences without the required interval/p-value support are
`not_established`. Any incomplete or unverified seed is `blocked`; no seed may
be replaced or silently omitted. No seed, cutoff, repeat or model is selected
after inspecting its outer result. The two new panels overlap and are not
treated as independent samples for pooling or p-value combination.

The separate prevalence-swing rule fixes test-minus-validation <= -0.10 as a
large drop and absolute swing <= 0.05 as stable. Both new seeds must satisfy a
category to call it reproduced or stable; intermediate/mixed patterns are
reported as such.

## Frozen splits and matching refits

The same grouping code recomputes 147 components using core distance <1e-4 or
canonical assembly identity, and `split_by_host_group` applies fractions
`[0.7,0.15,0.15]`. `_freeze_cohort` retains its exact manifest equality check
and original `ComparatorBlocked` exception. Host-parent variants carry the
new seeds; separate comparator and genome-retrieval variants bind the new
manifest and matching baseline result.

| Seed | Train / validation / test hosts | Train / validation / test groups | Groups moved relative to 42 | Complete / eligible validation hosts | Complete / eligible test hosts |
| --- | ---: | ---: | ---: | ---: | ---: |
| 42, original | 117 / 26 / 34 | 103 / 22 / 22 | 0 | 26 / 26 | 34 / 33 |
| 43 | 125 / 25 / 27 | 103 / 22 / 22 | 80 | 24 / 24 | 27 / 25 |
| 44 | 130 / 24 / 23 | 103 / 22 / 22 | 87 | 24 / 24 | 22 / 21 |

The 16,597 observed MLC labels and 41 unknown exclusions are unchanged.
`370D` is scored over all 94 candidates but excluded whole from seed-43
validation and seed-44 test metrics because of its 39 unknown labels. Eligible
denominators also exclude complete all-negative hosts. The
[input freeze](../data/studies/genophi_seed_sensitivity.inputs.v1.json) records
every changed group, partition count and input digest.

All cited learned baselines are refitted for each sensitivity seed. Their
representations, training-only preprocessing, per-phage logistic estimator,
C=0.1, lbfgs, max_iter=2000 and single-class fallback are the v1 settings.
The original numerical lock is used in a separate **Python 3.12.14** environment:
NumPy 2.3.5, pandas 2.2.3, scikit-learn 1.8.0, SciPy 1.17.0, joblib 1.5.3 and
threadpoolctl 3.6.0. Annotation, genome-only and hybrid models are all refitted;
the comparison below reports genome-only and hybrid, with prevalence and GenoPHI.

Every baseline model's serialized `training_hosts` must equal the same seed's
frozen training hosts and be disjoint from its targets. Protocol, manifest,
split table, model and score hashes are checked before the comparator trains.
A regression test explicitly rejects weights trained under another seed.

GenoPHI retains all original decisions: distribution 1.0.3, Pyrodigal 3.7.1,
MMseqs2 15-6f452+ds-2, original Python 3.12.3 numerical environment, cutoff 3,
`20 <= n < 501`, RFE 200 features x five runs, ten retained models x three
complete repeats, and the all-30-model median confidence rank score. The
unchanged grid, eight threads, outer-training-only supervised input and disabled
dynamic weights/SHAP are preserved. Cutoff 3 retains 143 features for seed 43
and 144 for seed 44, identically across their complete repeats.

## Results

### Test Hit@2: all complete hosts

| Seed | GenoPHI | Prevalence | Host 4-mer | Hybrid | GenoPHI - prevalence, pp | Paired 95% interval, pp | Group p | Holm p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 42, observed reference | 91.18% | 73.53% | 88.24% | 88.24% | +17.65 | [+6.00, +35.49] | 0.03125 | Descriptive reference |
| 43 | **70.37%** | **62.96%** | 74.07% | 74.07% | **+7.41** | **[-10.00, +24.16]** | **0.6875** | **0.6875** |
| 44 | **77.27%** | **65.15%** | 72.73% | 63.64% | **+12.12** | **[-1.45, +28.57]** | **0.1875** | **0.3750** |

These are point advantages against prevalence on all three test panels, so the
directional advantage is not confined to seed 42. Neither sensitivity interval
has a positive lower endpoint, and neither adjusted test passes: the registered
robustness claim remains **not established**. Against 4-mer, GenoPHI is one host
ahead at seed 42, one behind at 43, and one ahead at 44. This small margin changes
sign and does not establish representation superiority.

### Hit@2 under both denominators

Cells are **all-complete-host / eligible-complete-host percentages**.

| Seed | Panel | Complete / eligible | GenoPHI | Prevalence | Host 4-mer | Hybrid |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 42 | Validation | 26 / 26 | 76.92 / 76.92 | 84.62 / 84.62 | 69.23 / 69.23 | 84.62 / 84.62 |
| 42 | Test | 34 / 33 | 91.18 / 93.94 | 73.53 / 75.76 | 88.24 / 90.91 | 88.24 / 90.91 |
| 43 | Validation | 24 / 24 | 83.33 / 83.33 | 79.17 / 79.17 | 70.83 / 70.83 | 66.67 / 66.67 |
| 43 | Test | 27 / 25 | 70.37 / 76.00 | 62.96 / 68.00 | 74.07 / 80.00 | 74.07 / 80.00 |
| 44 | Validation | 24 / 24 | 66.67 / 66.67 | 75.00 / 75.00 | 75.00 / 75.00 | 66.67 / 66.67 |
| 44 | Test | 22 / 21 | 77.27 / 80.95 | 65.15 / 68.25 | 72.73 / 76.19 | 63.64 / 66.67 |

The [full metric table](../data/studies/results/genophi_seed_sensitivity.v1/metrics.csv)
reports **Hit@1/2/3/5 and MRR under both denominators for every seed, method and
validation/test panel**, plus eligible Recall@K, raw/complete/eligible host counts
and positive-candidate medians. Values are retained at full floating precision
in CSV/JSON; display rounding here is not used by the statistics.

### Prevalence validation-to-test swing

| Seed | Prevalence validation Hit@2 | Prevalence test Hit@2 | Test - validation, pp | GenoPHI test - validation, pp | Positive-candidate median, validation / test |
| --- | ---: | ---: | ---: | ---: | ---: |
| 42 | 84.62% | 73.53% | **-11.09** | +14.25 | 17.5 / 25.0 |
| 43 | 79.17% | 62.96% | **-16.20** | -12.96 | 13.0 / 15.0 |
| 44 | 75.00% | 65.15% | **-9.85** | +10.61 | 10.0 / 12.5 |

The seed-42 prevalence test result is actually the **highest** of these three
test point estimates, not a uniquely low one. Both sensitivities have drops
larger than the registered five-point “stable” bound. Seed 43 meets the
ten-point large-drop rule, while 44 narrowly misses it. We therefore report the
literal **`mixed_or_different`** classification, not “both reproduced,” and do
not move the threshold after seeing 44. The proposed pattern “43/44 are stable
and 42 is the isolated unusual drop” is not supported by these measurements.
Three overlapping development splits still do not identify the cause of a swing
or support a population outlier test.

## Statistics and fractional hosts

The [paired contrast table](../data/studies/results/genophi_seed_sensitivity.v1/paired_contrasts.csv)
contains **360** contrasts: three seeds x two panels x six method pairs x five
metrics x two denominators. Every row includes the paired 2,000-draw host-group
bootstrap interval (seed 20261024), exact two-sided group sign-flip p-value,
both-binary discordances and the non-binary host count, including zero.
Only the two preregistered new-seed primary tests receive the declared Holm
adjustment; other contrasts are descriptive unadjusted results.

The sign-flip statistic sums per-host differences within each frozen group,
then flips one sign per group. It preserves the host-macro weighting. Stored
binary64 values are converted to exact common-denominator integers; a
meet-in-the-middle calculation counts all tail assignments, including ties.
No fractional host is dropped, rounded or thresholded. Interpretation assumes
group-wise sign exchangeability/symmetry, conditional on the fitted models and
the fixed phage panel.

For example, the new seed-44 GenoPHI/prevalence test Hit@2 contrast has **22
hosts: 17 both-binary and five non-binary**. Its binary-subset discordance is
2:0, but the reported test uses all 22 hosts and gives **p=0.1875**. Applying
McNemar only to the binary subset would give 0.5 and would discard information.
The prevalence sum of expected Hit@2 values is 14 1/3 over 22 hosts; it is not a
binary count to round to 14 or 15. Seed-44 validation has six non-binary hosts
in that contrast as well. Seed-43 test has zero non-binary hosts and 4:2
discordances in six groups.

The original seed-42 primary and B2 GenoPHI/prevalence contrasts remain 6:0 and
0:7, respectively, with zero non-binary hosts and one discordant host in each
contributing group. Their group tests equal exact McNemar, 0.03125 and 0.015625.
The original primary GenoPHI/4-mer contrast is 1:0, p=1.0. Its nonnegative
bootstrap interval is structurally nonnegative because 4-mer never wins alone,
as now explained in the [v1 report](genophi-frozen-comparator.md).

## Execution, verification and provenance

| Seed | Baseline fit, s | Baseline QA, s | GenoPHI fit/scoring, s | GenoPHI QA process, s |
| --- | ---: | ---: | ---: | ---: |
| 43 | 69.65 | 66.46 | 3,368.8 | 147.9 |
| 44 | 44.9 | 40.0 | 3,710.4 | 169.1 |

The execution ledger records **7,617.38 s (2.12 h)** of completed initial stages
plus monotonic continued execution. It excludes the launcher interruption.
The original launcher hit a 120-second tool limit after seed-43 baseline fit
and QA, but **before the GenoPHI run directory, binding or any fit existed**.
The unchanged baseline outputs were hash-verified and reused; only originally
scheduled unfinished commands continued. The prior ledger is preserved beside
the final one. There were no additional completed-model refits, seed omissions,
replacements, extra repeats, or scientific blocked runs.

Both new comparator bundles retain 51 complete checkpoints and independently
reload all 30 models. Saved per-model scores reproduce bit-for-bit; independent
NumPy median aggregation reproduces the retrieval summaries. Each baseline QA
independently reconstructs the original feature representations, checks exact
training-only preprocessing, reloads all three learned models and reproduces
their candidate scores bit-for-bit under its original numerical environment.
Per-seed Beta(1,1) prevalence rates are independently recomputed from training
observations and agree bit-for-bit with both studies' full score panels.

The aggregate [numeric QA](../data/studies/results/genophi_seed_sensitivity.v1/NUMERIC_QA.json)
also independently brute-forces the six real primary Hit@2 p-values (three
seeds, both denominators) using `Fraction` arithmetic rather than the production
sign-flip engine, and reproduces the registered decisions. See the
[complete report](../data/studies/results/genophi_seed_sensitivity.v1/report.json)
for protocol/input/output hashes, source versions, accounting, all metrics and
the unchanged seed-42 reference.

### Reproduction

Use the frozen executable commit `c215d86` in a fresh checkout/worktree, before
the published sensitivity result directories exist. Keep published bundles
immutable. Prepare the original v1 raw inputs/proteomes using the
[v1 instructions](genophi-frozen-comparator.md), supply the **original v1
feature-cache artifact bound by the reuse lock**, and install the declared
baseline Python/package lock in a separate environment. This sensitivity's
cache import checks both payload and original source-provenance bytes; a new
cache derivation with different provenance is a new versioned input, rather
than an implicit replacement. From that checkout:

```bash
python scripts/run_genophi_seed_sensitivity.py \
  --outdir artifacts/genophi_seed_sensitivity_v1 \
  --baseline-python /path/to/python-3.12.14-baseline-env/bin/python \
  --genophi-python /path/to/original-genophi-env/bin/python \
  --feature-cache /path/to/verified/genophi_comparator_v1
```

The pure reporting scripts were added after the executable freeze. Place their
published versions in that checkout's `scripts/` directory after fitting; they
can then aggregate and verify its result bundles without new fits:

```bash
python scripts/summarize_genophi_seed_sensitivity.py \
  --execution artifacts/genophi_seed_sensitivity_v1 \
  --outdir /path/to/new-sensitivity-summary
python scripts/verify_genophi_seed_sensitivity.py \
  --run /path/to/new-sensitivity-summary --out /path/to/new-qa.json
```

## Boundaries

The all-genome protein-family vocabulary remains **transductive**, while the
4-mer representation is **inductive**. Changing split seeds does not remove
that asymmetry; the sweep does not isolate protein families versus 4-mers.
Nor does it make this previously inspected source cohort an independent release
holdout. The outcome concerns robustness on these two prespecified group splits
and source-defined MLC positivity, not calibration or clinical efficacy.
