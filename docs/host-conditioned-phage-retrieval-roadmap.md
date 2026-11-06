# Host-conditioned phage discovery and retrieval — proposed next task

**Status: architecture and evaluation proposal only.** This document does not
define a trained model, change existing endpoints, publish model weights, or
authorize a provider/GPU run. The historical source-pair model task remains
frozen in [first-model-task.md](first-model-task.md); the full-cohort Evo2
representation comparison in issue #24 retains its unchanged protocol.

## Next milestone: proposed offline discovery benchmark

A concrete but **not yet frozen** evaluation design now lives in
[retrieval-rescue benchmark v0](retrieval-rescue-benchmark-v0.md). It
predefines equal-budget learned/evidence candidate discovery comparisons,
source-positive rescue accounting, tie and missing-label semantics, and the
separate cold-phage generalization requirements. Review its checklist and
freeze a new versioned protocol before fitting any new retrieval model. The
already-inspected 177-host/94-phage development cohort remains exploratory.

## Why a retrieval-first objective

A mechanism-first workflow can lose true activity-positive phages before a
learned model sees them. An unannotated receptor, defense feature or phage
counter-defense is **missing information**, not proof of incompatibility.
Consequently a learned model restricted to an already-filtered candidate pool
cannot recover positives filtered out upstream.

**Research objective:** Given one admitted bacterial genome and a fixed,
versioned library of **existing phage genomes**, learn to discover and prioritize
phages with activity under a **named, verified source endpoint**, including
candidates unsupported by the currently annotated mechanism features.
This is candidate discovery for screening, not a declaration of productive
infection, mechanistic causality, cocktail benefit or therapy.

We want both high coverage and an auditable explanation of **where candidates
came from**, not an unqualified top-N score.

## Proposed scientific architecture

```text
Host genome / HostProfile / time-indexed HostState
                   |
           representation(s)
                   |
        +----------+-----------+
        |                      |
        v                      v
  Learned retrieval       Evidence retrieval
  host-conditioned        pinned host annotations,
  phage genome library    recorded host range, PHIST,
  scores + unresolved     receptor/defense evidence
        |                      |
        +----------+-----------+
                   v
        versioned UNION of candidates
        + provenance and missingness
                   |
          policy and evidence review
                   |
       screened candidate prioritization
                   |
       later pair-specific evaluation
                   |
          SDL CandidateSet / tests
```

Both branches must query a defined phage library or be explicit about their
separate coverage. A model-only hit can enter the candidate pool without passing
a prior mechanism rule. A mechanism-only hit is not silently dropped because a
model is uncertain. Hard exclusions need independent, explicit applicability
and safety/provenance evidence; unknown annotations are not hard negatives.
Any score fusion or quota allocation is versioned and measured, not assumed.

The predictor owns **learned representations, model-specific phage-library
encoding/index, candidate scores, training and evaluation**. SDL owns phage
source snapshots, HostProfile/HostState identity, mechanism evidence, candidate
pool union, decision/review policies, assay/execution context and feedback.
An indexed approximate-nearest-neighbour service is an optional future
optimization; exhaustive scoring is appropriate for the current 94-phage panel.
Neither repository should import the other's private implementation.

## Modeling options, not an architectural commitment

1. **Low-cost baseline:** known-library source-activity ranking from host
   features, training-only phage prevalence, or host 4-mer features. This is
   already studied on the fixed panel. Known-phage identity models do *not*
   demonstrate generalization to unseen phages.
2. **Sequence-conditioned pair scorer:** host and phage genome representations
   (with optional annotated features) predict the declared source activity
   score for each existing phage. Batching the scorer gives host-conditioned
   retrieval without a separate two-tower architecture.
3. **Two-tower retriever (hypothesis):** learn host/phage encoders in a shared
   interaction-supervised space. Precompute phage vectors, then retrieve by
   host-vector score. Plain cosine similarity of untrained Evo2 embeddings is
   *not* an infection model. Validate against the pair scorer before promotion.
4. **Optional reranker:** assess pair-specific evidence, uncertainty and
   applicability on a diverse, recall-preserving short list.

Evo2 host-embedding feasibility is an independent input study. It must neither
be represented as already providing learned retrieval nor be used to change
the frozen 177-host study in issue #24. Phage-side representation work is a
**separate subsequent protocol** with its own budget and data identity.

## Benchmark and acceptance design

Freeze **before fitting** the library identity/sequence hashes, host and phage
split manifests, training labels, missingness, assay context, feature availability,
deduplication, model versions, scoring budgets, top-K tie rules and evaluation
code. Do not tune a frozen development test after seeing its results.

Compare under the **same candidate pool and evaluation opportunity**:

- Training-only prevalence and random/uniform controls.
- Explicitly supported mechanism/annotation baselines (only when the
  available annotations and rule definitions justify them).
- Genome-only learned retrieval.
- Mechanism + genome fusion, with missing features represented as missing.
- Provenance-preserving union of the two discovery routes.

Primary research questions:

- At identical K and under an explicitly reported evaluated-label denominator,
  does learned retrieval improve **Recall@K** and/or **Hit@K** over the eligible
  evidence baseline? Report MRR and candidate coverage as secondary metrics.
- What fraction of source-positive phages **missed by the evidence route** are
  recovered by the learned route ("rescued positives")? This is a *conditional
  label-and-library measurement*, not proof of a novel biological mechanism.
- Do improvements persist under grouped-host holdout, frozen phage-family
  holdout and **double-cold** splits? Report subgroup failures and uncertainty,
  not only average performance.
- Does explicit mechanism information add independent value? Use ablations.
  Report costs, latency, candidate-pool recall and evidence missingness.

Important study distinctions:

- The pinned 177-host x 94-phage GenoPHI cohort is a **fixed-library
  development benchmark**. Its 16,597 audited observed readings support
  historical source-MLC positivity, with 41 missing-primary-readout entries
  excluded. Unobserved pairs are not negatives.
- PHIST shared-k-mer results and INPHARED depositor host annotations do not
  constitute observed infection or mechanistic ground truth.
- The current 177-host cohort does not independently prove cross-library,
  cross-species, or unseen-phage discovery. Cold-phage claims require separate,
  leak-safe data and evaluation.
- We cannot claim "previously unknown infectivity genes" merely because a
  learned model rescues a source-positive pair. Mechanistic attribution requires
  independent experimental investigation and controlled causal evidence.
- The already-inspected development dataset cannot serve as an independent
  release holdout. Model calibration, transportability and promotion are
  separate gates; source activity never becomes `productive_infection` merely
  by renaming the field.

## Proposed handoff (not yet implemented)

A future versioned **RetrievalResult** should record at least:

- exact target HostState/HostProfile and genome identity; timestamp/context;
- immutable library snapshot and eligible/unavailable phage records;
- per-candidate identity, phage genome identity, provenance and discovery route
  (`learned`, `evidence`, or both), plus any missingness;
- endpoint-specific ranking score vs explicit `not_estimated`, model and
  preprocessing digests, applicability limits, optional calibrated uncertainty;
- complete considered-pool counts, filtering reasons, K/ties and policy digest.

This is not `InteractionPrediction`, `PairPrediction v1`, or a probability
of `productive_infection`. Any schema/admission into SDL requires a separate
reviewed contract and a matching receiver policy. The current pair-exchange
contracts remain unchanged. Results must not be retrofitted into the existing
mechanistic CandidateSet policy.

## Ordered engineering milestones

1. Finish (or report infeasibility of) the **unchanged**, independently scoped
   Evo2 host-representation study under issue #24, with no unauthorized API
   calls and no favorable-result cherry-picking.
2. Specify and freeze a **new** retrieval benchmark, end-to-end candidate-pool
   coverage and a suitable independent release holdout. Review source labels
   and permitted uses before any new training.
3. Implement low-cost exhaustive library scoring and a separately tested
   host+phage sequence baseline. Measure recall gains and subgroup failures.
4. Compare the pair-scorer and two-tower alternatives if data volume, costs
   and cold-phage results justify complexity. Keep explicit uncertainty/abstention.
5. Define and test the retrieval-results-to-SDL admission contract, union
   candidate policy and reviewer; only then consider optimizer input and
   active-learning acquisition from physical QC-admitted observations.

**Non-goals:** new repository, changed source-label semantics, hypothetical
engineered phage generation, implicit clinical recommendations, automatic model
promotion, or modification of the existing PHIST baseline.
