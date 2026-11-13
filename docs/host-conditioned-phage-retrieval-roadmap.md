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
| Published BASEL/ECOR reference | [Sequence audit](basel-ecor-sequence-audit.md), [fixed controls](basel-ecor-controls.md) and [bilinear comparison](paired-interaction-results.md) complete. Bilinear primary Hit@2 bounds 48%–60%; not a release holdout |
| Evo2 host representation | Paused optional branch; bounded pilots and extraction/comparison software retained; full-cohort comparison unexecuted |
| Learned host–phage interaction ranking | [Bilinear sequence model trained and compared](paired-interaction-results.md); primary development AP improvement not established. No validated serving/retriever release |
| New-host fixed-library comparison | [Compatibility audit and protocol freeze complete](new-host-retrieval-protocol.md); matched training/scoring is pending |
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

Plan updated 2026-11-13. [Step 1 is complete](new-host-retrieval-protocol.md):
the compatibility ledger and new study protocol are frozen before new fitting
or scoring. **Next: implement and run step 2 under that protocol.** Public-data
intake can proceed alongside it. Model changes depend on the comparison and its
failure analysis.

### 1. Align the new-host, fixed-library comparison

- [x] Reuse the completed budget audit and availability ledger; record only the
  remaining gaps for prevalence, annotation, host 4-mer, GenoPHI and bilinear
  scoring. Check full candidate support, endpoint, sequence/feature identity,
  training hosts/phages, preprocessing and validation selection. Separate
  reusable scores, score-only regeneration, required refits and unavailable
  methods. Sparse paired-test exports cannot stand in for a full library.
- [x] Freeze a new study protocol for unseen E. coli host groups and the existing
  94-phage library, with all library phages available to training. Reuse the
  existing primary grouped-host split and retain B2 as a prespecified stress
  test. Bind inputs, source-MLC labels/unknowns, exact sequences, feature
  availability, training/validation scopes, method versions and selection rules
  before any new fit or scoring.
- [x] Reuse a fitted artifact only when those scopes match. The current crossed
  bilinear fit uses a restricted training phage set; comparing it directly with
  models trained on all 94 is not this matched experiment. Freeze necessary
  refits and resource limits, or report the method unavailable.

**Done when:** a compact compatibility ledger and a hash-bound executable study
protocol specify the comparison. Historical protocols and acceptance decisions
remain unchanged; the inspected cohort remains development data.

### 2. Run the matched ranking comparison and diagnose failures

- [ ] Score every candidate using the existing retrieval and budget evaluators.
  Keep all-negative complete hosts. Preserve the frozen incomplete-label-host
  exclusions and all 94 candidates; no test-label-based candidate filtering.
- [ ] Report host-macro Hit@2 as primary; K=1,3,5, Recall/MRR, capped tests consumed,
  coverage/abstention and subgroup failures as secondary. Freeze contrasts,
  acceptance/tradeoff rules and paired host-group intervals in step 1. Pair AP
  remains a diagnostic, not a replacement for within-host ranking.
- [ ] Retain validation, primary and B2 outcomes, including failures against
  prevalence and GenoPHI. Inspect existing score/feature errors to choose one
  concrete data or representation hypothesis, without another test-seed search.

**Done when:** one reproducible report contains matched methods, per-host scores,
costs, exclusions and a recorded next decision. Reanalysis of inspected data
cannot promote a model, even if a new development contrast is favorable.

### 3. Deliver a research FASTA-to-ranking command

- [ ] Load declared frozen preprocessing and weights, verify host FASTA and the
  library snapshot, then emit scores/ranks, endpoint, model/sequence hashes,
  applicability, missingness and explicit unavailable results.
- [ ] Reproduce the step-2 scorer on admitted inputs; check wrong inputs, unknown
  candidates and incomplete score support. A software demo can use a declared
  baseline even if the learned comparison fails; it must report that limitation.

**Done when:** one documented command produces a reproducible research candidate
list. Source scores remain uncalibrated. This does not implement SDL decision
admission; serving into decisions needs the retrieval contract and receiver
policy below.

### Public-data intake alongside steps 1–3

- [ ] Identify public new-host datasets with linked genome identities, measured
  positive/negative outcomes, source assay definitions and explicit unknowns.
  Audit source lineage, host/phage overlap and reuse terms before admission.
- [ ] Reserve suitable evaluation data before model selection. Distinguish
  additional same-study hosts from independent laboratory/assay validation.
  BASEL has known host IDs and stays an inspected cold-phage reference.

**Done when:** each proposed source has an intake decision and a frozen role,
or a concrete rejection reason. If no suitable release cohort is available,
research work continues and release validation remains pending.

### Conditional follow-up

- [ ] If failure analysis supports a change, freeze one bounded feature/model
  experiment and compare it with the same controls. Select using validation
  only; retain all held-out failures. Revise later tasks from that evidence.

Evo2 full-cohort extraction remains paused. New retrieval infrastructure,
two-route fusion/union, engineered variants and cocktail efficacy are later
studies, not prerequisites for the matched baseline or research command.

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
