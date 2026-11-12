# Phage–Host Interaction Predictor

Research infrastructure for phage-host interaction prediction and phage-cocktail
baselines. Dataset provenance, label semantics and evaluation boundaries are
explicit; this repository does not establish clinical efficacy.

## Next research direction: host-conditioned phage discovery

The proposed next-generation objective is **input a bacterial genome and retrieve
activity-supported existing phages from a versioned library**, including candidates
that known receptor/defense annotations cannot confidently identify. Model
retrieval and evidence/mechanism retrieval should be **complementary discovery
routes**, not an irreversible "mechanism gate, then model reranking" cascade.
Unknown annotations must remain unknown, not hard evidence of non-infection.

This does **not** imply the learned retriever exists, that Evo2 embeddings alone
predict infection, or that an AI-retrieved candidate has a verified novel
mechanism. The predictor owns representation/learned retrieval and evaluation;
SDL owns candidate-pool provenance, mechanism evidence, union/decision policies
and experimental feedback. A pair scorer may also be used to rank an entire
fixed library; a separate dual encoder is optional and must earn its complexity
through controlled evaluation.

See [host-conditioned retrieval roadmap](docs/host-conditioned-phage-retrieval-roadmap.md)
for non-goals, frozen-study boundaries, proposed data handoff, fair retrieval
benchmarks, and what it would take to claim discovery beyond known annotations.
Current implementation status and the original source-pair endpoint remain as
described below.

## Offline two-route retrieval rescue audit

A **provider-free, source-endpoint-specific evaluation utility** now accepts
immutable, hash-verified candidate-score panels for independent `evidence`
and `learned` retrieval routes. It evaluates each complete host x phage
panel with the existing exact-tie metrics, builds a deterministic **equal-K**
union candidate ordering, and reports known-positive phages rescued from
evidence-only top-K, unknown-label exclusions and grouped uncertainty.

This is **software for comparing externally supplied score files**, not a
trained learned retriever, mechanism detector, Evo2 inference, or clinical
prediction. It neither filters the phage library based on missing mechanism
annotations nor converts source-MLC labels to productive infection. See
[offline retrieval-rescue evaluator](docs/offline-retrieval-rescue-evaluator.md)
for the manifest, CLI and output contract.

```bash
pytest -q tests/test_retrieval_rescue.py
python scripts/run_retrieval_rescue.py --manifest /path/to/manifest.json --outdir /path/to/new-audit-dir
```

## First task

Predict a named source-defined binary outcome for one bacterial strain and one
existing phage, then evaluate on unseen strains with possibly known phages.
The task, label evidence, generalization limits and acceptance criteria are in
[`docs/first-model-task.md`](docs/first-model-task.md).

Current delivery: a working file-based model boundary, a real-data label/genome
audit, frozen host-feature weights, a completed additive host–phage benchmark
and fixed-library unseen-host retrieval comparisons using annotations and host genome composition.
The first host model did not establish improvement over phage prevalence. The new
4-mer pair baseline also did not meet its cold-phage improvement criterion across
the prespecified sensitivities. Both remain uncalibrated research scores; see the
[host study](docs/host-baseline-results.md) and
[paired study](docs/paired-baseline-results.md).
The [retrieval study](docs/host-retrieval-benchmark.md) finds test Hit@2 of 82.4%
versus 73.5% for phage prevalence, with an inconclusive paired gain interval and
weaker performance on several phylogroup holdouts.
The [host-genome comparison](docs/genome-host-retrieval.md) reaches 88.2% Hit@2
on the primary split, but falls to 75.4% on B2 holdout versus prevalence's 90.2%.
This inspected development result does not justify model promotion.
The [aligned GenoPHI comparator](docs/genophi-frozen-comparator.md) now reaches
91.2% primary Hit@2 and 78.7% on B2; its B2 result also trails prevalence's 90.2%.
Its [preregistered split-seed sensitivity](docs/genophi-seed-sensitivity.md)
retains positive test point differences against prevalence at seeds 43/44,
but does not meet its statistical robustness criterion. Prevalence's large
validation-to-test fall is not confined to seed 42.
Pairwise matrices and offline union coverage do not establish cocktail efficacy.

| Capability | Current status |
| --- | --- |
| File exchange | Bundle/resource validation, result validation and explicit abstention |
| Trained host baseline | Host features plus known phage identity; uncalibrated source-MLC scores; no demonstrated improvement |
| Baseline serving through PairPrediction | Not implemented; weights and the generic exchange are separate paths |
| Mechanistic predictions / unseen phages | Not established by the current baseline |
| Host genome composition | Fixed 4-mer and annotation-hybrid retrieval comparison; cross-phylogroup failures retained |
| External genome-aware comparator | GenoPHI on the same frozen host split and 94-phage panel; primary and B2 results verified |
| Pretrained genomic models | Evo2 published species-score reproduction only; no aligned strain-level embeddings |

## Repository layout

```text
data/
  raw/              # downloaded source data; not committed
  processed/        # cleaned/model-ready data; not committed
  registry/         # versioned dataset manifests and label semantics
  endpoints/        # explicit source-label endpoint descriptors
  studies/          # frozen protocols, resource locks and audit reports
models/             # frozen development JSON weights and model card
notebooks/          # exploration and experiments
src/phage_cocktail/
  artifacts.py        # canonical JSON hashes, file digests, code provenance, staged outputs
  model_interface.py  # exchange validation and CLI
  input_bundle.py     # frozen input/resource consumer and PairPrediction
  analysis.py         # interaction-table summary statistics
  data.py             # observed interaction-table loading/validation
  registry.py         # dataset checksum and summary verification
  split.py            # deterministic strain- and group-level splits
  benchmark.py        # host-agnostic comparator baselines and metrics
  cohort.py           # MLC reconciliation, FASTA digests and genomic host grouping
  host_baseline.py    # train-only preprocessing and known-phage logistic baseline
  host_study.py       # locked real-data audit/training/evaluation runner
  pair_cohort.py      # exact-sequence identity audit for paired inputs
  paired_splits.py    # phage-neighbour grouping and crossed paired partitions
  paired_baseline.py  # 4-mer features, additive pair models and grouped metrics
  retrieval.py        # fixed-panel Hit@K/Recall@K/MRR, exact ties and paired intervals
  retrieval_study.py  # locked known-phage retrieval fits and host/phylogroup holdouts
  genophi_reference.py # optional GenoPHI reference workflow helpers
  genophi_proteomes.py # pinned, receipted single/meta code-11 gene calling
  genophi_comparator.py # train-only GenoPHI adapter and checkpointed retrieval study
  schemas/            # ModelInputBundle / PairPrediction v1 wire schemas
tests/              # automated tests
```

## Dataset and labels

The first dataset is the pinned GenoPHI E. coli subset. The fetch script verifies
the registered SHA-256, and analysis scripts verify the checksum and summary
statistics before use.

The canonical `infection` column represents a source-coded binary label. A later
row-level audit found **41 zeros whose primary MLC readings are missing**. The
new study excludes these, retaining 16,597 observed labels that exactly match
`MLC > 0`. The target combines source spotting-assay categories; a mapping to
adsorption, productive infection or resistance emergence remains unverified.
Spotting/clearing, plaque formation, productive infection and growth suppression
are not interchangeable labels; [endpoint semantics](docs/endpoint-semantics.md)
defines their review boundaries without changing the frozen source endpoint.
See [the audit](docs/host-baseline-results.md),
[`docs/dataset_contract.md`](docs/dataset_contract.md) and
[`data/README.md`](data/README.md).

Genomes are not committed either. The 177 host genomes (910.4 MB) and 94 phage
genomes (6.8 MB) are pinned by URL, SHA-256 and byte size in the resource locks,
so they are rebuilt by verification rather than re-derived.
[`data/README.md`](data/README.md#rebuilding-the-locked-genome-inputs) gives the
commands, the expected sizes and which download hosts each half needs; the phage
half and every table need only GitHub.

## Model exchange

SDL exports **ModelInputBundle v1** with complete pair snapshots, an explicit
source endpoint/context, scientific FeatureArtifact provenance and optional local
FASTA payloads. The consumer verifies all supplied bytes and rejects missing
model-required sequences. **PairPrediction v1** binds the exact bundle and
separates unresolved, uncalibrated-score and probability outputs.

```bash
phage-model-interface validate-bundle tests/fixtures/input_bundle/bundle.json --require-sequence host
```

The fixture is synthetic and tests software only. See
[`docs/model-interface.md`](docs/model-interface.md) for the complete round trip.

## Setup

Requires Python 3.11+.

```bash
python -m venv .venv
```

Activate the environment, then install:

```bash
pip install -e ".[dev]"
pytest
```

## Step 1: fetch and inspect the registered dataset

```bash
python scripts/fetch_genophi.py
python scripts/explore_interactions.py
```

For interactive exploration, open `notebooks/01_explore_interactions.ipynb`.

## First observed-label host-feature study

```bash
python scripts/run_host_baseline.py --fetch --outdir artifacts/host_baseline_v1
```

This verifies 177 host genomes, 94 phage genomes and all annotation/source tables
(approximately 920 MB), excludes unknown MLC readings, freezes host groups and
fits the [fixed protocol](docs/host-baseline-protocol.md). Existing run directories
and tampered input files are rejected. Weights, metrics, row predictions and
provenance are saved together. Frozen development weights are in `models/`.

Test AP is 0.5413 versus 0.4952 for per-phage prevalence; the paired 95% AP-gain
interval is [-0.0229, 0.1092]. This is inconclusive, and Brier/log loss are worse.
See [the full report](docs/host-baseline-results.md) for the split and limitations.

## Paired benchmark protocol

The next paired-model study has a separate
[development protocol and executable input preflight](docs/paired-benchmark-protocol.md).
It audits sequence/label identity and duplicate-genome disagreements without
training. [The follow-up source review and split freeze](docs/paired-split-freeze.md)
now fixes 69 operational phage-neighbour groups for the primary scheme and all
nine development split manifests. It traces source annotation differences without
relabeling them. The [completed fixed baseline study](docs/paired-baseline-results.md)
runs host-only, phage-only and additive pair models, plus applicable comparators,
under both weighting policies. Cold-phage gains remain inconclusive. The already
inspected cohort cannot serve as independent release validation.

## Unseen-host candidate retrieval

Evaluate how often the first 1/2/3/5 recommended phages include a source-positive
candidate, on a fixed 94-phage panel:

```bash
python scripts/run_host_retrieval.py --fetch --outdir artifacts/host_retrieval_v1
```

This table-only study refits the fixed host baseline on the frozen grouped-host
split, a random-host sensitivity and every published phylogroup holdout. It
reports all-host and eligible-host denominators, exact expectations for ties,
Recall@K/MRR and paired host-group Hit@2 intervals. Missing panel labels exclude
whole hosts from primary metrics. No GPU is required. See the
[protocol, results and limits](docs/host-retrieval-benchmark.md).

## Host genome and annotation comparison

```bash
python scripts/run_genome_host_retrieval.py --fetch --outdir artifacts/genome_host_retrieval_v1
```

This separately frozen experiment verifies 177 host FASTAs (approximately 911 MB
including tables) and compares host 4-mer composition, annotations and their
concatenation under the same splits, 94-candidate panel and fixed per-known-phage
logistic specification. It needs no GPU and reports every phylogroup failure.
See [the protocol, full results and verification](docs/genome-host-retrieval.md).

## Host similarity and Evo2 comparison protocol

The next host comparison adds a fixed five-neighbour control, per-host distance
and prevalence-shift diagnostics, and a separately specified Evo2-40B cache path:

```bash
python scripts/run_host_similarity.py --fetch --outdir artifacts/host_similarity_v1
```

See [similarity results, the Evo2 protocol and GPU handoff](docs/host-similarity-evo2.md).
The local-checkpoint 40B extraction and GPU hardware pilot have not run; the deployment lock must
be created on the existing GPU installation before extraction.

Authenticated NVIDIA-hosted 40B forward access has now been verified on short
controls, including the fixed layer's 8,192-dimensional outputs. This is separate
from the local GPU/genome pilot. See [API access, transfer measurements and the
bounded probe](docs/nvidia-evo2-api.md). Actual selected lengths refine the full
plan's JSON-transfer projection to about 1.425 TB, versus the earlier full-window
estimate of about 2 TB. Full hosted extraction has not been started.

The [real-host feasibility and bounded engineering pilot](docs/nvidia-evo2-host-pilot.md)
use a separately frozen hosted protocol, per-orientation checkpoints and an
aggregate download budget. They measure access, pooling and cache recovery;
susceptibility prediction and the full representation comparison remain separate.

The [CPU window-sampling control](docs/host-window-sampling.md) compares 1, 4
and up to 32 full windows with whole-genome 4-mer composition under the same
fixed-panel retrieval protocol. It reports sampling loss, actual coverage and
transfer projections without new NVIDIA calls. Its results do not establish
Evo2 quality or a sufficient embedding sampling budget.

The [resumable hosted executor](docs/nvidia-evo2-extraction.md) adds local/S3
immutable checkpoints, frozen call/download budgets, exact pilot-cache imports
and offline verification. Its defaults make no new NVIDIA requests. A real-vector
replay checks interrupted and completed-cache recovery; a live cloud deployment
and hosted prediction-study adapter remain separate work.

## Optional GenoPHI genome-aware reference

GenoPHI's official 25-strain x 25-phage protein-family demo can verify the path
from protein FASTA files through MMseqs2-derived features to strain-separated
prediction.

On Linux/macOS, or inside WSL2 on Windows:

```bash
conda env create -f environment.yml
conda activate phage-cocktail
python scripts/bootstrap_genophi_demo.py
python scripts/run_genophi_reference.py
python scripts/summarize_genophi_reference.py
```

GenoPHI is pinned to version 1.0.3 and its demo data to commit
`f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f`. See
[`docs/step3_genophi_reference.md`](docs/step3_genophi_reference.md).

This demo now runs end to end here, taking 1,893 s on 2 cores. It verifies the
pipeline only. The reference output shipped in the upstream repository is not a
valid comparison target: its `cutoff_8` to `cutoff_10` rows are unreachable at
the documented demo setting of `--num_runs_fs 5`, because feature occurrence is
capped by the number of feature-selection runs, and its generating parameters
are recorded nowhere upstream. The upstream code at the commit that ships that
output also cannot complete the documented demo at all. Measured curves,
the excluded explanations and the consequences for using GenoPHI as a
comparator are in
[`docs/genophi-demo-reproduction.md`](docs/genophi-demo-reproduction.md).

## GenoPHI comparator on the frozen host split

GenoPHI 1.0.3 now provides an aligned external comparator on the same reconciled
16,597 labels, grouped-host split and fixed 94-phage library as the in-house
models. Protein-family vocabularies are built once over the locked genomes;
supervised feature selection and modeling use outer-training hosts only.
Occurrence cutoff 3 is fixed a priori, and scoring uses this project's exact-tie
Hit@K implementation. The whole comparator took 4.13 hours on eight threads.

| Test panel | GenoPHI Hit@2 | Host 4-mer Hit@2 | Phage prevalence Hit@2 |
| --- | ---: | ---: | ---: |
| Primary, all 34 complete hosts | **31/34 = 91.18%** | 30/34 = 88.24% | 25/34 = 73.53% |
| B2, all 61 complete hosts | **48/61 = 78.69%** | 46/61 = 75.41% | 55/61 = 90.16% |

Eligible-only GenoPHI Hit@2 is 31/33 (93.94%) primary and 48/60 (80.00%) B2.
The primary gain over the host 4-mer model is one host, with a paired interval
touching zero. B2 still exposes a substantial failure against prevalence.
Three complete replays agree, and all 60 retained models reproduce their saved
scores bit-for-bit on independent reload. See the
[protocol, versions, receipts, full metrics and verification](docs/genophi-frozen-comparator.md)
for reproducible gene calling and the resumable comparator command.

The [seed-43/44 sensitivity](docs/genophi-seed-sensitivity.md) refits matching
baselines and reuses the same protein vocabulary. Test GenoPHI/prevalence
Hit@2 is **70.37%/62.96%** at seed 43 and **77.27%/65.15%** at seed 44. Both
group-bootstrap gain intervals cross zero; Holm-adjusted group sign-flip
p-values are 0.6875 and 0.375, so the preregistered robustness claim is
**not established**. Prevalence falls 16.20 and 9.85 points from validation to
test, rather than remaining stable outside seed 42. Fractional tie expectations
are preserved in every test; B2 is seed-invariant and was not refitted.

## Evo2 retrieval reference

Reproduce the upstream species-level retrieval metrics without a GPU:

```bash
python scripts/run_evo2_reference.py --fetch
```

This downloads approximately 49 MB of pinned upstream scripts, mapping and
precomputed score matrices, verifies every file's SHA-256, and runs the original
five single-method evaluations and prespecified four-way RRF. All six methods
match the published rounded reference metrics on 4,490 phages. This is a
**precomputed-score reproduction**; Evo2 embedding generation and the 40B
comparison have not been run. See the
[results and next experiment](docs/evo2-retrieval-reference.md).

## Milestones

- [x] Create a reproducible Python project skeleton
- [x] Import a pinned public phage-host interaction dataset
- [x] Register source revision, checksum, expected statistics and label caveats
- [x] Define deterministic held-out-strain splitting and basic baselines
- [x] Guard against train/test strain overlap and missing-pair coverage bias
- [x] Wire in the official GenoPHI genome-aware reference workflow
- [x] Define the first pairwise task and preserve unverified assay semantics
- [x] Consume frozen input bundles and return exact source-endpoint results
- [x] Reconcile source binary labels with published MLC readings and exclude unknowns
- [x] Verify all cohort genomes and supplied biological annotations
- [x] Freeze genomic host groups and run a fixed host-feature baseline
- [x] Report paired host-group uncertainty, including an inconclusive result
- [x] Evaluate fixed-panel unseen-host Hit@K/Recall@K/MRR with explicit eligibility and phylogroup sensitivities
- [ ] Finish row-level assay-context and separate dataset-license review
- [x] Run the GenoPHI reference workflow end to end and establish that its shipped reference output is not a reproduction target
- [x] Evaluate GenoPHI on the frozen grouped-host and B2 splits with a priori cutoff selection and aligned retrieval metrics
- [x] Run preregistered grouped-host seed sensitivities with matching refits, exact group-level fractional-metric tests and an inconclusive robustness result
- [ ] Define target-specific label mapping before training mechanistic targets
- [ ] Demonstrate improvement and validate calibration before model promotion
- [x] Evaluate unseen-phage regimes with frozen sequence groups and report inconclusive additive-baseline gains
- [ ] Add repeated group evaluation and an independent release holdout
- [ ] Freeze a separate host-conditioned retrieval study with a fixed phage-library snapshot and labels
- [ ] Evaluate genome-only, evidence-only, fusion and union retrieval on the same candidate pool
- [ ] Measure rescued source-positive candidates, recall at equal K, and cold-host / cold-phage failures
- [ ] Define a reviewed, versioned learned-retrieval result and SDL admission contract before serving
- [ ] Compare explicit host–phage pair features with available pretrained representations under the same protocol
- [x] Reproduce Evo2 species-level retrieval metrics from pinned published score matrices
- [ ] Regenerate the Evo2-7B embeddings and compare 40B under a separately frozen protocol

## Scope

This repository is a research prototype for computational prediction. It does
not establish clinical efficacy or substitute for experimental validation.
