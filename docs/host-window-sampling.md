# Host-window sampling controls before a larger Evo2 run

This CPU experiment measures how a fixed composition model changes when it sees
1, 4 or up to 32 sampled windows instead of whole host assemblies. The purpose is
to quantify a sampling control and its prospective hosted transfer cost. It does
not measure Evo2 quality or identify an optimal embedding sampling budget.

The [amendment](../data/studies/host_window_sampling.v1.json) pins the complete
[genome study](../data/studies/genome_host_retrieval.v1.json) by SHA-256 rather
than duplicating its 177-genome resource lock. Cohort identities, all 94 named
candidate phages, MLC unknown exclusions, complete-panel/eligible denominators,
exact score-tie expectations and all 12 split scenarios remain unchanged.
Every numeric encoder is fitted over unique training hosts only; C=0.1 and
maximum 2,000 iterations remain fixed. All sampling contrasts are exploratory
on an already inspected development cohort, without noninferiority, model
selection, release promotion or an independent confirmation claim.

## Frozen sampling

Canonicalize contig orientation/order exactly as in the hosted engineering pilot.
Generate its 32 systematically spaced windows within contigs, then retain only
full 8,192-base windows. Short tails are never padded or joined. From this same
pool, budget 1 takes the first full window, budget 4 takes evenly spaced floor
indices, and budget 32 takes all full members. The three selections are nested.
If fewer windows exist, use all available members and report actual support.
Reject hosts without a usable full sampled window; do not silently drop them.

The first three hosts' one-window inputs must match the actual hosted pilot's
input hashes. This experiment excludes short sampled tails, so its up-to-32
control is distinct from the original local Evo2 protocol that includes short
tails. Selection can omit short contigs and useful loci; genome composition is
not a guarantee of receptor/defense coverage or an embedding's information needs.

Each sampled window is treated as a separate contig for canonical 4-mer counting.
Skip ambiguous windows in place, with no artificial cross-window adjacency.
Normalize pooled counts over valid within-window 4-mers. Reverse complement has
identical canonical 4-mer counts, so a separate strand count is unnecessary for
this control. That invariance does not apply to contextual Evo2 representations.
QC and coverage are exported but are not additional predictive features.

| Method | Inputs |
| --- | --- |
| Random / global / phage prevalence | Existing train-only controls |
| Annotation model | Existing defense/surface features |
| Whole-genome 4-mer | Same 136 canonical frequencies as the parent |
| Whole-genome hybrid | Parent genome + annotation control |
| 1 / 4 / up-to-32-window 4-mer | Same 136 frequencies computed only within selected windows |

## Reproduce and verify

```bash
pip install -e '.[dev]' -c data/studies/host_retrieval_environment.v1.txt
python scripts/run_host_window_sampling.py --fetch --outdir artifacts/host_window_sampling_v1
python scripts/run_host_window_sampling.py --outdir artifacts/host_window_sampling_v1_repeat

python scripts/verify_host_window_sampling.py \
  --run artifacts/host_window_sampling_v1 \
  --repeat artifacts/host_window_sampling_v1_repeat \
  --out artifacts/host_window_sampling_qa.json \
  --budget-out artifacts/host_window_sampling_budget.csv
```

`--fetch` obtains approximately 911 MB of pinned host FASTAs and tables, with
exact byte/size verification. Existing verified files are reused. Figshare 503s
during this run required retries using the same file IDs via a fresh download
query or its API download route; all admitted bytes match the original lock.
No NVIDIA request, key, GPU inference or model fitting through a provider is
needed. Use a new output directory on each run; failed studies publish no
partial experiment.

Artifacts include the effective inherited protocol and source amendment, full
genome and sampled feature tables, per-budget assembly coverage, exact sampled
window coordinates/input hashes, all 72 fitted models, complete scores, per-host
metrics and paired group contrasts. All historical six methods are retained so
their per-host metrics can be checked exactly against the parent experiment.

Offline QA reloads all saved models, enumerates primary top-2 ties independently,
and checks scalar 4-mer counts for the first three frozen hosts at every budget.
The budget projection uses one actual real-host 8,192-base API response and the
observed full-window counts in this cohort. It is not a measured full extraction,
invoice or schedule. Different windows can compress and run differently. The
projection gives no credit for cached pilot inputs and does not authorize an
automatic provider run.

## Results

Protocol and implementation were published at `851c4fc042315309ebc67873d9d1dfa62f062109`
before new metrics. All 177 original FASTAs and four tables matched their frozen
byte locks. Every host had usable full-window support; no source host was dropped.
The first three one-window hashes match the actual hosted pilot exactly.

The grouped test contains 34 complete-panel hosts (33 eligible). Validation has
26 complete/eligible hosts. Both columns below use every complete-panel host.

| Method | Test Hit@2 | Validation Hit@2 |
| --- | ---: | ---: |
| Phage prevalence | 0.7353 | 0.8462 |
| Annotations | 0.8235 | 0.8462 |
| Whole genome 4-mer | 0.8824 | 0.6923 |
| Whole genome + annotations | 0.8824 | 0.8462 |
| 1 full window | 0.6765 | 0.8077 |
| 4 full windows | 0.7059 | 0.6923 |
| Up to 32 full windows | 0.6765 | 0.8462 |

The three sampled genome controls achieve 23/34, 24/34 and 23/34 test Hit@2,
versus whole-genome composition at 30/34 and prevalence at 25/34. Their paired
conditional group-bootstrap contrasts against whole-genome composition are:

| Sampling control | Hit@2 difference, percentage points | Conditional 95% interval |
| --- | ---: | ---: |
| 1 full window | -20.59 | [-36.37, -7.41] |
| 4 full windows | -17.65 | [-33.33, -2.94] |
| Up to 32 full windows | -20.59 | [-37.78, -2.50] |

These exploratory intervals condition on the fitted models, split and fixed
library. They do not establish an independent biological discovery. More sampled
windows do not consistently improve this model: validation and several
phylogroups move differently from the main test. The comparison cannot identify
why sampling changes performance, or determine how Evo2 would behave.

| Test holdout | Complete / eligible | Prevalence | Whole genome | 1 window | 4 windows | Up to 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| random_host | 26 / 25 | 0.7308 | 0.8846 | 0.7308 | 0.7308 | 0.6923 |
| phylogroup:A | 37 / 36 | 0.6216 | 0.5405 | 0.4595 | 0.3514 | 0.5676 |
| phylogroup:B1 | 22 / 21 | 0.6364 | 0.6364 | 0.6818 | 0.6364 | 0.5909 |
| phylogroup:B2 | 61 / 60 | 0.9016 | 0.7541 | 0.7541 | 0.7869 | 0.8033 |
| phylogroup:C | 7 / 7 | 1.0000 | 0.7143 | 0.7143 | 0.7143 | 0.7143 |
| phylogroup:Clade I | 3 / 3 | 0.6667 | 0.6667 | 0.6667 | 1.0000 | 0.6667 |
| phylogroup:Clade IV | 1 / 1 | 1.0000 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| phylogroup:D | 17 / 17 | 0.8824 | 0.8824 | 0.6471 | 0.7647 | 0.9412 |
| phylogroup:E | 15 / 15 | 0.9333 | 0.8667 | 0.8667 | 0.7333 | 0.8667 |
| phylogroup:F | 8 / 8 | 0.3750 | 1.0000 | 0.7500 | 0.5000 | 0.5000 |
| phylogroup:G | 4 / 4 | 0.5000 | 0.7500 | 0.2500 | 0.7500 | 0.7500 |

B2 remains difficult: whole genome is 75.41%, sampled controls are 75.41%,
78.69% and 80.33%, and prevalence is 90.16%. This does not solve cross-lineage
generalization. No model or sampling strategy is promoted.

## Coverage and transfer projection

| Full-window cap | Median selected windows | Median assembly coverage | Projected JSON download |
| --- | ---: | ---: | ---: |
| 1 | 1 | 0.161% | 63.24 GB |
| 4 | 4 | 0.643% | 252.97 GB |
| 32 | 21 | 3.423% | 1308.10 GB |

For cap 32, actual full-window support ranges from 9 to 31, with median 21.
Short tails explain much of the gap between a nominal 32-window count and
its effective input coverage. The original Evo2 protocol retains those short
windows; it remains a distinct representation from this full-window control.

Using actual original-window lengths, its linear transfer estimate is **1.425 TB**,
rather than the earlier approximately 2.024 TB estimate that treated every
window as a full 8,192 bases. Both are projections of intermediate API responses,
not raw FASTA sizes, model weights or retained host-vector cache sizes.
If one float32 8,192-dimensional vector is retained per host, all 177 vectors
contain only 5,799,936 numerical payload bytes (about 5.8 MB). API downloads remain
large because token-level matrices must reach the client before local pooling.

The one-window projection is about 63.24 GB, four windows about 252.97 GB, and
up to 32 full windows about 1.308 TB. These give no credit for reuse. Service-side
pooling would avoid transferring full token matrices, but that capability has
not been verified on the hosted endpoint. No provider extraction was launched.

## Numeric verification

All six historical methods reproduce every committed parent per-host metric
exactly. All 72 saved models reload with exact scores; sampled frames restore
their original C-order layout when loaded from CSV. Scalar reference counts
verify nine first-three-host feature vectors. Independent top-2 tie enumeration
checks 306 grouped-test host/method panels. Two complete runs reproduce the
entire report and every listed output hash. The suite passes 125 tests and eight
subtests. The [committed record](../data/studies/results/host_window_sampling.v1/report.json)
and adjacent per-host, coverage, coordinate, contrast, cost and QA tables retain
all results. Large models and full candidate-score exports remain reproducible
outputs; source genomes are not redistributed. No NVIDIA inference call or
credential was used in this study.
