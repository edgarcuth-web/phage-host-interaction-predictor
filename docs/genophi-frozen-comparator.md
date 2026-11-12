# GenoPHI on the frozen 177-host x 94-phage cohort

GenoPHI now runs as an aligned external protein-family comparator on this
project's frozen split. Primary test Hit@2 is **31/34 = 91.1765%**, beside the
existing host 4-mer model's **30/34 = 88.2353%** and prevalence's **25/34 =
73.5294%**. On B2 holdout it falls to **48/61 = 78.6885%**, below prevalence's
**55/61 = 90.1639%**. The external comparator therefore retains the important
cross-phylogroup failure rather than resolving it.

Protocol and implementation were committed at
`40ffe55a1bbefc5ddbca749e9585c99493479d84` before comparator fitting. This is a
development comparison on the already inspected source-MLC cohort. The
[verified result bundle](../data/studies/results/genophi_comparator.v1/report.json)
and [independent numerical QA](../data/studies/results/genophi_comparator.v1/NUMERIC_QA.json)
retain the full results, support, paired contrasts and provenance.

## Inputs and gene calling

The 276 resources in `genophi_host_baseline.resources.v1.json` have been rebuilt
and verified: 177 host assemblies, 94 phage assemblies and five source tables.
Receipts are committed; nucleotide assemblies, source tables, proteomes and
intermediate alignment databases are local generated inputs.

- [`host_fastas.v1.json`](../data/receipts/host_fastas.v1.json): 177 matching host
  files, 910,399,604 bytes. This final verification reuses the downloaded files.
- [`genophi_inputs.v1.json`](../data/receipts/genophi_inputs.v1.json): every
  source file's SHA-256 and byte size against the original resource lock.
- [`genophi_proteomes.v1.json`](../data/receipts/genophi_proteomes.v1.json): pinned
  caller/settings, input/output digests and counts for all 271 proteomes.

Pyrodigal **3.7.1**, genetic code **11**, is used for both sets. Each bacterial
assembly is trained independently in normal/single mode, using all its contigs;
genes are then called separately on each contig. Every phage uses meta/anonymous
mode, without short-genome self-training. The
[current Pyrodigal API documentation](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html)
distinguishes trained single mode from pretrained meta mode; Prodigal's
[input-type guidance](https://github.com/hyattpd/Prodigal/wiki/Advice-By-Input-Type#plasmids-phages-viruses-and-other-short-sequences)
recommends anonymous mode for isolated sequences below approximately 100 kb.

Importantly, meta mode uses **only built-in bins whose training information has
translation table 11**. The default bin collection also includes three code-4
models. Merely requesting code 11 when translating already predicted genes
would not constrain the preceding gene prediction.

Ends remain open, masking is disabled, minimum gene/edge lengths are 90/60 nt,
maximum overlap is 60 nt, and terminal stops are omitted from protein FASTAs.
Contig order is preserved. A single counter across all contigs in one assembly
generates `>370D_1`, `>370D_2`, etc.; the file is exactly `370D.faa`. Neither a
contig-name mapping nor a new genome-ID normalization is introduced.

The measured gene-calling run used eight workers and **52.90 s** wall time:

| Set | Genomes | Proteins | Protein FASTA bytes |
| --- | ---: | ---: | ---: |
| Bacteria, single mode | 177 | 872,342 | 278,579,663 |
| Phages, meta mode | 94 | 9,841 | 2,187,439 |

## Prespecified comparator

The machine-readable [protocol](../data/studies/genophi_comparator.v1.json)
fixes these decisions before fitting:

- GenoPHI distribution **1.0.3**, MMseqs2 **15-6f452+ds-2**, Pyrodigal **3.7.1**.
  Runtime receipts also hash all installed GenoPHI Python sources and the
  MMseqs executable. `genophi --version` still says 1.0.2; distribution metadata
  is the version check. Numerical dependencies are pinned in
  [`genophi_comparator_environment.v1.txt`](../data/studies/genophi_comparator_environment.v1.txt).
- One host protein-family vocabulary over all 177 host genomes and one phage
  vocabulary over all 94 phages, shared across scenarios/repeats. Identity 0.4,
  coverage 0.8, sensitivity 7.5; upstream perfect-co-occurrence collapse and no
  singleton-family or sample-cluster presence filter. The original protocol's
  descriptive phrase "singleton removal" is corrected in the
  [pre-fit metadata erratum](../data/studies/genophi_comparator.v1.errata.json):
  GenoPHI runs that filter only when a genome-selection list is supplied, and
  this adapter uses `select='none'` over the complete locked genome sets.
- `curate_mlc` admits **16,597** observed labels; **41** missing MLC readings
  remain unknown. The phenotype export contains only `strain,phage,interaction`.
- Recompute `group_hosts` at threshold `1e-4`, then `split_by_host_group` at seed
  42 and fractions `[0.7,0.15,0.15]`. Require exact equality with the frozen host
  manifest. The primary split has **117/26/34** train/validation/test hosts.
  B2 uses the existing phylogroup holdout with touching host-group purge.
- Subset full feature-table rows to outer-training hosts **before** calling
  any supervised GenoPHI function. Its internal feature selection, MCC grid
  search and early stopping operate entirely within this subset. Candidate
  prediction receives IDs and protein features, with no held-out labels.
- RFE selects **200** features in each of **5** feature-selection runs. Fix
  occurrence **cutoff 3 a priori**, with the explicit window **20 <= n < 501**.
  Occurrence counts selection runs, not genomes. At most `5*200/3` features can
  pass, so the upper bound is explicit and nonbinding. An unsupported fixed
  stratum is a blocked run, not permission to change the cutoff.
- Retain **10** internally selected models per complete workflow and repeat the
  entire feature-selection/modeling stage **three** times. Internal filtering
  uses strain protein-content hierarchical clustering with 20 clusters, as in
  the upstream machinery. Dynamic weighting and SHAP export are disabled.
  Grid: iterations 500/1000, learning rate 0.05/0.1, depth 4/6, Logloss, and
  early stopping after 100 rounds. Upstream models are not refitted on all
  outer-training rows after their internal selection.
- The primary rank score is the **median positive-class confidence over all 30
  models**. Report complete-repeat medians and individual-model spread as well;
  never select a run using an outer result. Fixed upstream seeds mean these
  replays can agree exactly. Repeatability is not independent sampling uncertainty.

The protein vocabulary is transductive in genome sequences, as requested; the
supervised boundary is outer-train-only. This distinction is explicit because
the earlier fixed 4-mer representation is inductive. This bounded configuration
is an aligned comparator, not a reconstruction of the authors' paper benchmark
or the unrecorded upstream reference CSV.

### Memory-bounded assignment

The installed upstream implementation searches its protein database against
itself, then reads and sorts the entire all-versus-all alignment table to choose
one best hit per query. The adapter keeps the same search and family mappings,
but reduces alignment chunks before merging: float32 sequence identity
descending, float32 e-value ascending, first input row on exact ties. No search
threshold, family definition, or label-dependent rule changes.

Tests compare against the installed upstream reducer on randomized ties, extreme
e-values and chunk boundaries. Each real clustering run also compares the first
50,000 alignments byte-for-byte against the upstream reducer before feature use.
This limits the reduction's memory cost without silently changing the method.

## Execution and wall time

The run used eight threads on WSL2 Ubuntu 24.04, with approximately 15 GiB RAM
available to the VM and GenoPHI's advisory `max_ram=8` setting. It produced
45,538 host and 2,177 phage families, collapsed to **11,455 host + 346 phage =
11,801** protein-family features. The fixed occurrence cutoff retained **147**
features on the grouped-host split and **133** on B2, identically in all three
complete repeats. Both are inside the prespecified `20 <= n < 501` window.

| Measurement | Wall seconds |
| --- | ---: |
| Gene calling, outside the comparator run | 52.90 |
| Host clustering, all-versus-all search and bounded assignment | 8,504.80 |
| Phage clustering/search/assignment | 26.81 |
| Complete comparator, including all six complete fits and scoring | **14,883.32 (4.13 h)** |
| Independent post-run verification | 242.88 |

The [runtime receipt](../data/studies/results/genophi_comparator.v1/runtime.json)
uses the Linux boot-time clock and process start ticks (one-second polling).
Completed stages sum to 14,453.88 seconds using `time.perf_counter`; the total
also includes work between stage checkpoints. These monotonic measurements are
used because the host clock changed during the run, affecting upstream elapsed
time strings and UTC timestamp differences. The runtime receipt separately
identifies its memory-sampling start and scope.

The strain assignment reduced **136,946,512** alignments to **872,342** best
hits. Its first 50,000 alignment rows and the corresponding phage sample match
the unmodified upstream reducer byte-for-byte. All **103** stage receipts are
complete, including 30 feature-selection runs and 60 retained modeling runs.

## Reproduce and resume

On Ubuntu 24.04 / WSL2, install MMseqs2 (the measured apt package is
`mmseqs2=15-6f452+ds-2`) and create/activate a Python 3.12 environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev,reference]' -c data/studies/genophi_comparator_environment.v1.txt
python -c "from importlib.metadata import version; print(version('genophi'))"
mmseqs version
```

Rebuild the nucleotide inputs with the commands in
[`data/README.md`](../data/README.md#rebuilding-the-locked-genome-inputs), then:

```bash
python scripts/call_genophi_genes.py \
  --out data/receipts/genophi_proteomes.local.json \
  --source-receipt data/receipts/genophi_inputs.local.json --workers 8

python scripts/run_genophi_comparator.py \
  --proteome-receipt data/receipts/genophi_proteomes.local.json \
  --outdir artifacts/genophi_comparator_v1

python scripts/verify_genophi_comparator.py \
  --run artifacts/genophi_comparator_v1 \
  --proteome-receipt data/receipts/genophi_proteomes.local.json \
  --out artifacts/genophi_comparator_v1.NUMERIC_QA.local.json
```

Use a **new receipt filename**: committed receipts are immutable historical
records, not files to overwrite. Per-genome checkpoints under
`data/derived/genophi_proteomes/checkpoints/` allow verified gene-calling reuse.
The same comparator command resumes an interrupted run if tool versions,
protocol, code and completed artifact hashes still match. Named stages
`--stage prepare`, `--stage cluster`, and `--stage fit` stop at boundaries;
`--scenario grouped_host --repeat 0` restricts fitting to that frozen unit.

Each feature-selection/modeling run has its own UTF-8 log and hash receipt.
Failures record the exception and elapsed time, and do not publish a completed
stage or comparison. Altered inputs, missing models, partial candidate panels,
and a missing fixed feature stratum are rejected. `data/derived/` and
`artifacts/` are ignored by Git.

Byte-hashed tracked files must retain LF line endings. `.gitattributes` enforces
this for fresh Windows/WSL checkouts; no frozen protocol or source hash is
redefined to accommodate CRLF conversion.

## Metric contract and results

All 94 candidates are scored for each held-out host. Metrics use
`phage_cocktail.retrieval.evaluate_retrieval`: exact tie expectations,
Hit@1/2/3/5, eligible-only Recall@K and MRR, with both all-complete-host and
positive-eligible-host denominators. A host with any unknown panel label is
excluded whole, and an all-negative complete host contributes Hit/MRR zero.
Group bootstrap contrasts use the existing paired retrieval implementation.

The frozen [host-genome result](genome-host-retrieval.md) is the in-house
reference. The earlier models below are its committed results; the prevalence
scores were also independently recomputed on each scenario's training rows and
their per-host Hit@2 agrees with that reference.

### Hit@2 beside the existing host-genome result

| Panel | Complete / eligible hosts | Prevalence, all hosts | Host 4-mer, all hosts | GenoPHI, all hosts | GenoPHI, eligible only |
| --- | ---: | ---: | ---: | ---: | ---: |
| Grouped-host validation | 26 / 26 | 22/26 = 84.6154% | 18/26 = 69.2308% | 20/26 = 76.9231% | 20/26 = 76.9231% |
| Grouped-host test, primary | 34 / 33 | 25/34 = 73.5294% | 30/34 = 88.2353% | **31/34 = 91.1765%** | **31/33 = 93.9394%** |
| B2 test | 61 / 60 | 55/61 = 90.1639% | 46/61 = 75.4098% | **48/61 = 78.6885%** | **48/60 = 80.0000%** |

The primary fit receives 117 training hosts and 10,957 observed pairs. B2 uses
115 training hosts and 10,771 observed pairs. All 62 B2 target hosts receive
all 94 scores (5,828 pairs); `ECOR53`, with two unknown labels, is excluded whole
from metrics, leaving 61 complete hosts and 60 eligible hosts. The primary test
and B2 each retain one all-negative complete host in their all-host denominator.
The existing hybrid's test Hit@2 is also 30/34 primary and 46/61 B2.

### All-complete-host retrieval

| Panel | Method | Hit@1 | Hit@2 | Hit@3 | Hit@5 | MRR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Primary test | Prevalence | 0.5588 | 0.7353 | 0.8529 | 0.8529 | 0.7021 |
| Primary test | Host 4-mer | 0.7647 | 0.8824 | 0.8824 | 0.9118 | 0.8339 |
| Primary test | GenoPHI | **0.8235** | **0.9118** | **0.9118** | **0.9706** | **0.8809** |
| B2 test | Prevalence | 0.7459 | 0.9016 | 0.9016 | 0.9672 | 0.8395 |
| B2 test | Host 4-mer | 0.5902 | 0.7541 | 0.7541 | 0.8033 | 0.7002 |
| B2 test | GenoPHI | 0.7049 | 0.7869 | 0.9016 | 0.9508 | 0.7999 |

### GenoPHI positive-eligible-host retrieval

| Panel | Hosts | Hit@1 | Hit@2 | Hit@3 | Hit@5 | MRR | Recall@1 | Recall@2 | Recall@3 | Recall@5 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Validation | 26 | 0.6923 | 0.7692 | 0.8077 | 0.9231 | 0.7773 | 0.0463 | 0.0927 | 0.1326 | 0.2039 |
| Primary test | 33 | 0.8485 | 0.9394 | 0.9394 | 1.0000 | 0.9076 | 0.0576 | 0.1046 | 0.1276 | 0.1965 |
| B2 test | 60 | 0.7167 | 0.8000 | 0.9167 | 0.9667 | 0.8132 | 0.0327 | 0.0562 | 0.0973 | 0.1569 |

### Paired uncertainty and repeatability

All-complete-host Hit@2 contrasts use the existing paired host-group bootstrap
(2,000 resamples, seed 20261024). These intervals condition on fitted models,
this development split and the fixed phage library.

| Test panel | Contrast | Difference, percentage points | Paired 95% interval |
| --- | --- | ---: | ---: |
| Primary | GenoPHI - prevalence | +17.65 | [+6.00, +35.49] |
| Primary | GenoPHI - host 4-mer | +2.94 | [0.00, +11.54] |
| B2 | GenoPHI - prevalence | -11.48 | [-20.00, -4.84] |
| B2 | GenoPHI - host 4-mer | +3.28 | [-9.72, +15.69] |

The primary comparison with 4-mer is one additional successful host, and its
interval touches zero. The B2 comparison with 4-mer is also inconclusive. The
external comparator supplies a useful reference point, but neither the small
primary difference nor the B2 failure supports model promotion.

Three complete same-machine replays produced identical selected feature sets
and retrieval metrics; complete-repeat median Hit@2 SD is **0** for every panel.
Individual models differ across the ten internal run splits:

| Panel | Individual-model mean Hit@2 | SD, percentage points | Range |
| --- | ---: | ---: | ---: |
| Validation | 79.04% | 4.67 | 71.15–84.62% |
| Primary test | 88.82% | 2.24 | 85.29–91.18% |
| B2 test | 79.51% | 2.71 | 75.41–83.61% |

These descriptive SDs use all 30 recorded models per panel; they contain ten
internal run splits repeated three times, rather than 30 independent draws.
No individual model or repeat was selected for reporting on an outer metric.

## Independent verification and saved artifacts

The [numeric QA receipt](../data/studies/results/genophi_comparator.v1/NUMERIC_QA.json)
records fresh verification of all 276 resources and 271 proteomes, re-curation
of 16,597 labels, and hash verification of **103 receipts / 528 unique local
artifacts**. All **60 models** were independently reloaded; their full candidate
score arrays reproduce saved per-model scores **bit-for-bit**. NumPy median
aggregation agrees with the saved ensembles within `2.22e-16`, and every
retrieval summary reproduces from reloaded models, including exact tie handling
and unknown-host exclusion.

The committed bundle includes aggregate candidate scores, per-host metrics,
individual-model metrics, feature occurrence counts, host splits, all stage
receipts, runtime measurements, provenance and the full machine-readable report.
Large alignment, training-feature and model payloads are reproducible local
artifacts. The metadata erratum's public pre-fit record is linked from its JSON
and its hash is included in numeric QA.

The earlier [host-baseline report](host-baseline-results.md) describes source
reconciliation and the classification experiment; its AP is not Hit@K.
GenoPHI's internal AUC/MCC outputs are training-fold diagnostics, not the
aligned retrieval comparison. Scores still concern source-coded MLC positivity,
not productive infection, treatment response or cocktail efficacy.
