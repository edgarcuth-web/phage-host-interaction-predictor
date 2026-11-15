# Training-only complementary screening panels

**Complementary selection improves on prevalence in the primary test panel,
but loses to prevalence on the B2 holdout. It does not establish a stable
cross-group benefit or the proposed low-prevalence rescue mechanism.**

This completes the first comparison following the
[missed-candidate audit](missed-candidate-audit.md). All data were already
inspected. The held-out hosts were excluded from candidate selection, but are
not new independent validation data. No model or screening policy is promoted.

## Same budget, training labels only

The [protocol](../studies/complementary_panel/protocol.json) was committed as
[7b12a95](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/7b12a954a911bf9e50a84dca107178a8f8060ed8),
and the runner and its checks as
[0f4d78c](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/0f4d78ce9b93f74db39a6d3ae37b16a85460dcff),
before this full comparison ran. This fixes an exploratory analysis; it is not
a preregistration on unseen data.

For each existing training split, select an ordered pair from the same 94 named
phages. Each step maximizes the number of training hosts with a known positive
not yet covered. Ties prefer higher smoothed training prevalence, then ascending
phage name. Use the pair unchanged for every target host. Neither validation nor
test labels select a candidate. The runner also checks that no host group crosses
the training/target boundary.

Known positives from incomplete training hosts can contribute; missing readings
never enter the prevalence denominator or become negatives. Primary summaries
exclude whole incomplete target panels, retain all-negative hosts, and use the
same source MLC > 0 endpoint and saved comparator scores. Exact-score ties are
averaged for prevalence, host 4-mer and GenoPHI. A test costs one unit; stop after
the first source-positive result, or after two unsuccessful tests. Capped
consumption is `2 - Hit@1`.

## Test results

| Method | Primary Hit@2, 34 hosts | Primary capped tests | B2 Hit@2, 61 complete hosts | B2 capped tests |
| --- | ---: | ---: | ---: | ---: |
| Fixed complementary pair | 85.29% (29/34) | 1.441 | 83.61% (51/61) | 1.279 |
| Train-only prevalence | 73.53% (25/34) | 1.441 | 90.16% (55/61) | 1.254 |
| Host 4-mer | 88.24% (30/34) | 1.235 | 75.41% (46/61) | 1.410 |
| GenoPHI | 91.18% (31/34) | 1.176 | 78.69% (48/61) | 1.295 |

B2 scores all 62 hosts; ECOR53 has two unknown readings and remains outside the
complete-panel denominator. Its top-two hit is known under all four methods;
the separate lower/upper bounds preserve that information. Primary validation
results for 26 hosts are retained separately, without selecting a new pair.

| Complementary pair versus comparator | Primary rescued / newly missed | B2 rescued / newly missed |
| --- | ---: | ---: |
| Prevalence | 4 / 0 | 0 / 4 |
| Host 4-mer | 0 / 1 | 8 / 3 |
| GenoPHI | 0 / 2 | 5 / 2 |

Primary prevalence rescues are **NILS36, NILS44, NILS46 and NILS56**, four
different frozen host groups. B2's newly missed hosts versus prevalence are
**381A, DEC1A, ED1A and NILS70**, also four different groups. Favorable comparisons
against the weaker B2 genome models do not erase the loss against prevalence.

Relative to prevalence, the primary Hit@2 difference is +11.76 percentage points
[+2.56, +25.00]; B2 is −6.56 points [−13.85, −1.41]. Capped-consumption differences
are 0.000 [0.000, 0.000] and +0.025 [−0.076, +0.117], respectively. Brackets are
pointwise 95% paired host-group percentile intervals from 2,000 draws, conditional
on these fixed pairs and existing models. They do not account for prior data
inspection or candidate-selection uncertainty. Do not pool panels or methods
as independent samples or interpret these intervals as confirmatory tests.

## Named candidates and stability

| Training split | Ordered fixed pair | Second candidate's training positives / observed | Extra training hosts beyond first | Extra test hosts beyond first |
| --- | --- | ---: | ---: | ---: |
| Primary, 117 hosts | DIJ07_P2 → LF73_P1 | 70/116 | 17 | 10 in 6 groups |
| Non-B2, 115 hosts | 536_P7 → LF110_P2 | 53/114 | 21 | 7 in 5 groups |

The last column measures the second candidate's marginal coverage relative to
the first candidate alone. It is different from rescue relative to the
prevalence **pair**: LF73_P1 adds ten primary-test hits beyond DIJ07_P2 alone,
but only four beyond the two-candidate prevalence policy. LF110_P2 adds seven
B2 hits beyond 536_P7 alone; prevalence's pair covers those hosts too.

LF73_P1 is third in training prevalence (60.17% after smoothing). LF110_P2 has
prevalence rank interval 9–10 (46.55%). Both exceed their split's median rate
(17.80% and 14.53%). Thus this run does not support a claim that low-prevalence
phages supplied the observed improvement. The full candidate table retains all
94 candidates per target partition; no target contribution changed selection.

Reselecting from 2,000 training host-group bootstrap samples puts LF73_P1 in the
two-candidate panel in **29.45%** of draws, and LF110_P2 in **17.20%**. These are
selection frequencies, not success probabilities or independent replications.
They indicate sensitivity to training composition. In the B2 training set,
536_P7 and DIJ07_P2 also tie on initial coverage and prevalence; the predeclared
name rule picks 536_P7. All maximum-gain ties are retained in the selection trace.

## The original difficult cases remain visible

| Host | Primary fixed pair | B2 fixed pair | Source evidence |
| --- | --- | --- | --- |
| NILS53 | Both MLC 0; missed | Both MLC 0; missed | Still 11 positive library candidates; nine raw annotations per audited library pair |
| NILS78 | Both MLC 0; missed | Both MLC 0; missed | Still seven positive library candidates; nine raw annotations per audited library pair |
| AL505 | DIJ07_P2 MLC 4; LF73_P1 MLC 3 | 536_P7 MLC 4; LF110_P2 MLC 3 | Aggregate matrix only; individual raw annotations remain unavailable |

AL505 is rescued relative to B2 host 4-mer, but prevalence already finds it.
Neither fixed pair solves NILS53 or NILS78. Their full positive-candidate lists
remain in the earlier audit. No stable general rescue candidate has emerged, so
this study does not advance a receptor or defense-system explanation.

For a laboratory using DIJ07_P2 plus 536_P7, LF73_P1 is a concrete candidate to
check against the four named primary-test gaps, using the source readings and
matched strain identities. This evidence does not justify replacing that panel
across host populations. A future independent cohort can test the frozen
candidates; further tuning on these inspected cases would remain exploratory.

## Evidence and reproduction

The [result directory](../data/studies/results/complementary_panel.v1/) retains:

- `selected_panels.csv`: chosen names, ordered gains, training host IDs and ties.
- `candidate_marginal_coverage.csv`: all candidates' train-only rates and named
  target contributions, with unknown-anchor cases kept unresolved.
- `summary.csv`, `per_host_metrics.csv`, `paired_comparisons.csv`: complete-panel
  denominators, exact-tie outcomes, consumption, differences and intervals.
- `host_transitions.csv`, `changed_host_evidence.csv`: every paired outcome and
  the original MLC evidence for changed hosts, including newly missed controls.
- `focus_source_evidence.csv`, `top2_source_evidence.csv`: all score blocks
  intersecting the two-test budget, original row/column links and prior raw-audit
  status where available. Unreviewed raw annotations are not claimed as verified.
- `training_bootstrap_selection.csv`, `QA.json`, `INDEPENDENT_QA.json`,
  `PUBLICATION.json`: stability counts, input/code hashes and verification.

Python 3.11+ standard library suffices. Run from the repository root; output
directories must be new. `--fetch` downloads only the pinned 158 KB source matrix,
which is kept locally rather than republished. Its primary archive reuse review
and attribution remain in `data/studies/paired_data_reuse_review.v1.json`.

```bash
python -m unittest discover -s tests -p test_complementary_panel.py -v
python studies/complementary_panel/run_comparison.py --fetch \
  --out artifacts/complementary_panel
python studies/complementary_panel/verify_comparison.py \
  --bundle artifacts/complementary_panel
```

The independent verifier imports no runner helpers. It checks all 16,638
source pairs, 488 target metric rows including incomplete bounds, 363 paired
transitions, 978 top-two evidence records, 4,000 training resamples and 36,000
paired interval draws. A second full run is byte-identical for all eleven runner
outputs. Seven targeted tests cover complement selection, target isolation,
unknown denominators, group multiplicity, exact ties and capped consumption.

GenoPHI retains its inherited transductive genome-vocabulary limitation. Source
MLC positivity does not establish productive infection, cocktail efficacy or
treatment benefit.
