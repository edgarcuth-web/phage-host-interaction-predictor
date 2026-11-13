# Host-conditioned phage retrieval roadmap

**Objective:** given an admitted bacterial genome and a versioned library of
existing phages, rank candidates with activity under a named source endpoint.
Candidates with missing receptor/defense annotations remain eligible: unknown
mechanism evidence is not evidence of non-infection.

[Reference cases and project decisions](reference-cases.md) connects the
literature to the priorities below, with separate evidence and applicability
limits for each case.

## Implemented and pending

| Component | Status |
| --- | --- |
| Annotation, host 4-mer, additive pair and GenoPHI comparisons | Completed development studies; see the [results index](README.md#results-and-reproduction) |
| Two-route score/rescue audit | [Offline evaluator implemented](offline-retrieval-rescue-evaluator.md); consumes supplied scores |
| Fixed-library screening budgets | [Existing-score audit complete](screening-budget-audit.md); grouped-host test cost savings coexist with B2 failure |
| Evo2 host representation | Paused optional branch; bounded pilots and extraction/comparison software retained; full-cohort comparison unexecuted |
| Sequence-conditioned learned retrieval | Not trained or validated |
| Retrieval benchmark v0 | [Proposed design](retrieval-rescue-benchmark-v0.md); scientific protocol not frozen |
| RetrievalResult-to-SDL admission | Proposed; no serving/admission contract implemented |

The 177-host/94-phage cohort is an already-inspected development benchmark.
Its 16,597 observed MLC labels exclude 41 missing readings. Known-library results
cannot establish unseen-phage, cross-species or release performance. The latest
[GenoPHI sensitivity](genophi-seed-sensitivity.md) did not establish its registered
robustness claim.

## Candidate discovery and ownership

Evidence retrieval and learned retrieval query a defined library independently.
Their versioned union retains route provenance and missingness, so a model can
nominate candidates unsupported by available annotations. A source-positive pair
rescued by a model shows ranking complementarity; causal mechanism discovery
requires separate experimental evidence.

| Predictor owns | SDL owns |
| --- | --- |
| Model-specific representations and library encoding | Source library snapshots and HostProfile/HostState identity |
| Candidate scores, training and evaluation | Mechanism evidence and provenance |
| Model preprocessing, applicability and calibration | Union/decision policies, assay context, experiments and feedback |

Exhaustive pair scoring is sufficient for the current 94-phage library. A separate
two-tower encoder, approximate-neighbour index or reranker must justify its added
complexity through comparison with that scorer. Untrained Evo2 cosine similarity
is not a validated infection score.

## Paused Evo2 branch

Decision, 2026-11-13: pause new Evo2 development and full-cohort extraction.
The completed engineering pilots establish access and cache integrity, with no
strain-level ranking result. The current hosted plan projects about 1.425 TB of
intermediate responses for roughly 5.8 MB of final host-vector payload; see the
[source-verified preflight](https://github.com/edgarcuth-web/phage-host-interaction-predictor/issues/24#issuecomment-6488966629).
Issue #24 is closed as `not_planned`.

Retain protocols, pilots, scores and extraction code for reproducibility. Shared
input validation and evaluation remain available to the main retrieval work.
This is a resource-priority decision, with Evo2 predictive value still untested
on this cohort. Any revival needs a separately frozen, bounded experiment and
prespecified acceptance against inexpensive baselines. It cannot block the work
below or reuse a zero-budget plan as authorization for inference.

## Next work

1. Use the [completed budget audit](screening-budget-audit.md) to specify a
   concrete generalization/data hypothesis worth a new study. Cost savings on
   inspected grouped-host tests coexist with B2 failure and mixed validation
   results; an explicit residual model remains optional.
2. Before fitting a new retrieval model,
   freeze library/sequence hashes, labels, missingness, feature availability,
   grouped splits, model versions, K/ties and candidate budgets. Reserve
   independent release data; the inspected cohort remains development data.
3. Compare prevalence, annotation proxy, genome-only pair scoring, fusion and
   two-route union at equal final K. Report Hit/Recall/MRR, rescued positives,
   coverage, costs, grouped uncertainty and subgroup failures. Cold-phage and
   double-cold claims require phage-sequence groups and separate held-out labels.
   Test alternative architectures when the evidence justifies them, and keep
   unsupported or unavailable predictions explicit.
4. Define a versioned retrieval result and matching SDL receiver policy before
   serving or passing scores to candidate decisions.

A retrieval result must bind the host/genome/context, library snapshot,
per-candidate identity and route, missingness, endpoint score or `not_estimated`,
model/preprocessing hashes, applicability, considered-pool counts, filtering
reasons and K/tie policy. This contract is still proposed. The existing
[PairPrediction exchange](model-interface.md) remains a separate interface.

Source activity does not become `productive_infection` by renaming a field.
Mechanistic targets, engineered phage generation, cocktail benefit and clinical
efficacy are outside the evidence supplied by these retrieval studies. See
[endpoint semantics](endpoint-semantics.md) and the
[benchmark design](retrieval-rescue-benchmark-v0.md) for the detailed boundaries.
