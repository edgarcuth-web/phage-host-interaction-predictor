# Known-phage retrieval on unseen hosts

The host-feature model reaches **28/34 test hosts at Hit@2 (82.35%)**, versus
**25/34 (73.53%)** for train-only phage prevalence. The paired host-group 95%
interval for the gain is **[-3.23, +23.53] percentage points**. Improvement remains
inconclusive. Broader phylogroup holdouts reveal failures that the random-host
split conceals; these scores do not justify model promotion.

## Question and fixed protocol

Given an unseen host and a fixed collection of existing phages, does the model
place at least one source-positive phage in its first K recommendations?
This measures individual-phage retrieval, without asserting mixture efficacy.

The [v1 protocol](../data/studies/host_retrieval.v1.json) was written before
calculating these retrieval metrics. It reuses an already inspected development
cohort, so this is a new analysis, not an independent confirmatory benchmark.
The protocol SHA-256 is
`32d75d4b54717d01f47a792e208685e3789671741071a2e5fa97734e92fe25b0`.

The runner admits the four checksum-locked source tables from the existing host
study, reconciles MLC labels, and verifies the registered subset and frozen host
manifest. It uses 16,597 observed pairs, 177 host IDs and all 94 phage IDs.
The 41 unknown MLC pairs remain unknown. It reuses the frozen manifest's
147 genomic host groups; it does not download or re-audit the approximately
920 MB of genomes/core-distance resources.

All models receive the same complete 94-phage candidate panel. Every candidate
must have a finite score. A host missing any panel label is excluded as a whole
from primary metrics, with its identity and unknown-pair count reported. The
evaluator never shrinks a host's ranked library to its observed or positive pairs.
Both incomplete hosts happen to be in the primary training partition, so all
34 primary test hosts have complete labels.

## Metrics and denominators

- Hit@1/2/3/5: host-macro fraction with at least one positive among the first K.
  This is the top-K host-coverage question; it is not Recall@K.
- Recall@K: fraction of a host's positive panel members retrieved, averaged over
  eligible hosts. It is undefined for an all-negative host.
- MRR: reciprocal rank of the first positive, zero for an all-negative host.
- Complete-panel host fraction: hosts with every candidate label / held-out hosts.
- Eligible-host fraction: complete-panel hosts with at least one positive /
  complete-panel hosts. Hit/MRR are reported with both this eligible denominator
  and the denominator including all complete-panel hosts.

Exact score ties are averaged analytically over uniform permutations within the
tied block. Thus a constant global-prevalence score has the same expected ranking
performance as uniform random ordering. Neither candidate IDs nor source CSV
ordering break ties. A score export records minimum/maximum rank for every tie;
these expectations are not one realized ordered recommendation list.

## Fits and splits

The unchanged per-phage host-feature logistic specification uses C=0.1 and 2,000
maximum iterations. Its 137 defense counts and four surface categories use
train-host-only preprocessing. Single-class training phages use the existing
Beta(1,1) fallback. Each scenario refits independently; frozen all-host weights
are never reused in a new holdout. Baselines are uniform random ordering,
train-only global prevalence and Beta(1,1) per-phage prevalence. No hyperparameter
search or evaluation-driven feature selection was performed.

1. **Primary: frozen grouped-host split.** 117/26/34 train/validation/test hosts,
   with no genomic host component shared across partitions.
2. **Random-host sensitivity.** The existing deterministic FNV assignment of host
   IDs, seed 42 and 70/15/15 group-count allocation, gives 124/27/26 hosts.
   Five host components occur in both train and test, making this split permissive.
3. **Published phylogroup sensitivity.** Leave each observed `Clermont_Phylo`
   category out, fit on the rest, and report every category. If a near-clonal
   component touches the test group, remove its other members from fitting.
   The actual cohort needs no such collateral purges. These are published
   metadata categories, not newly reconstructed phylogenetic clades.

All split manifests are written before the first fit. The paired Hit@2 bootstrap
draws whole frozen host groups, preserving all their hosts and pairing the two
models. It uses 2,000 draws, seed 20261024 and percentile 95% intervals. The
primary contrast includes all complete-panel hosts; eligible-only contrasts are
also saved. Intervals condition on the fitted models, library and chosen split.
One-group panels have no estimated interval. No favorable sensitivity is selected
as a new primary result, and no multiplicity-adjusted claim is made.

## Primary results

Test panel: 34 complete hosts in 22 genomic groups; 33 eligible hosts and one
all-negative host. Metrics below include all 34 hosts.

| Method | Hit@1 | Hit@2 | Hit@3 | Hit@5 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| Uniform random expectation | 0.2656 | 0.4371 | 0.5520 | 0.6892 | 0.4482 |
| Global prevalence | 0.2656 | 0.4371 | 0.5520 | 0.6892 | 0.4482 |
| Phage prevalence | 0.5588 | 0.7353 | 0.8529 | 0.8529 | 0.7021 |
| Host-feature logistic | 0.6471 | 0.8235 | 0.9118 | 0.9412 | 0.7770 |

Eligible-only Hit@2 is 28/33 (84.85%) for host logistic and 25/33 (75.76%) for
phage prevalence. Eligible-host fraction is 33/34 (97.06%). Eligible-only Recall@2
is 0.0687 versus 0.0516. Recall is lower than Hit because each host can have many
positive phages, while Hit asks only for the first success.

Validation has 26 complete, eligible hosts. Both models have Hit@2 of 22/26
(84.62%). Test improvement is +8.82 percentage points, with the interval crossing
zero. This does not overturn the earlier inconclusive host-feature study.

## Sensitivity results

Every row below uses **all complete-panel hosts**, including all-negative hosts.

| Holdout | Held-out / complete / eligible hosts | Phage prevalence Hit@2 | Host-feature Hit@2 |
| --- | ---: | ---: | ---: |
| Random host | 26 / 26 / 25 | 0.7308 | 0.8846 |
| Phylogroup A | 37 / 37 / 36 | 0.6216 | 0.6216 |
| Phylogroup B1 | 22 / 22 / 21 | 0.6364 | 0.5909 |
| Phylogroup B2 | 62 / 61 / 60 | 0.9016 | 0.8033 |
| Phylogroup C | 7 / 7 / 7 | 1.0000 | 0.8571 |
| Clade I | 4 / 3 / 3 | 0.6667 | 0.6667 |
| Clade IV | 1 / 1 / 1 | 1.0000 | 0.0000 |
| Phylogroup D | 17 / 17 / 17 | 0.8824 | 0.8235 |
| Phylogroup E | 15 / 15 / 15 | 0.9333 | 0.9333 |
| Phylogroup F | 8 / 8 / 8 | 0.3750 | 0.6250 |
| Phylogroup G | 4 / 4 / 4 | 0.5000 | 0.7500 |

The random-host gain has an interval above zero, but related-host train/test
overlap limits its interpretation. B2 holdout instead favors phage prevalence:
gain -9.84 points, descriptive interval [-17.74, -3.28]. Small groups yield unstable
results; the one-host Clade IV panel has no group-bootstrap interval.
`ECOR53` (two unknown labels) and `370D` (39 unknown labels) are excluded from
complete-panel B2 and Clade I metrics, respectively.

## Reproduce and inspect

From an installed checkout, using a new output directory:

```bash
pip install -e ".[dev]" -c data/studies/host_retrieval_environment.v1.txt
python scripts/run_host_retrieval.py --fetch --outdir artifacts/host_retrieval_v1
```

The four source tables total approximately 600 kB. Existing host-study caches can
be reused with `--data-root data/raw/genophi_host_study`; every required file is
still hash-checked. No GPU, embedding generation or genome download is required.

The full run publishes the protocol, curated/excluded labels, all host splits,
12 refitted models, all candidate scores and tie-rank intervals, per-host metrics,
contrasts and a report with code/input/output hashes. Existing outputs are
rejected and failed runs do not publish partial output directories. The committed
[results](../data/studies/results/host_retrieval.v1/) contain the aggregate report,
split/per-host tables, contrasts and numerical QA; raw tables, full scores and
model weights are reproducible local outputs.

Numerical QA verified all 19 report-listed output hashes and current source
hashes. An independent repeat produced identical scores, splits, host metrics and
contrasts; all 12 reloaded models reproduced their host scores bit for bit.
Direct enumeration of choices at the top-2 tie boundary independently matched
the primary four-method results. The full suite passes 89 tests plus eight
subtests, including exhaustive small-panel rank permutations and failure cases
for unknown labels, score support, group leakage and partial run publication.

The reusable `evaluate_retrieval` API in `phage_cocktail.retrieval` accepts separate
observed labels and complete-panel score tables. Future GenoPHI or genome-model
scores must use the same admitted hosts, library and partitions before comparison.
The official 25x25 GenoPHI demo and the species-level Evo2 retrieval study use
different tasks/cohorts and are explicitly **not run as comparators here**.

This study adopts the unseen-host/known-library recommendation question discussed
in relation to Omnii. It does not reproduce Omnii's models, panel, eligibility
rules or reported values. Source MLC positivity, individual productive infection,
cocktail effects and patient benefit remain separate endpoints. Named library
members also include sequence duplicates; no outcome-based merging is performed.
The next model comparison needs a separately frozen protocol and the same
denominators, especially under phylogroup shift.
