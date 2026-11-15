# Host-conditioned phage retrieval roadmap

**Objective:** given an admitted bacterial genome and a versioned library of
existing phages, rank candidates under a named source endpoint. Missing
receptor/defense annotations do not exclude a candidate.

## Implemented and pending

| Component | Status and evidence |
| --- | --- |
| 94-phage development comparisons | [Matched comparison](new-host-retrieval-results.md) and [budget audit](screening-budget-audit.md) complete; bilinear improvement not established |
| Fixed panels and Picard extension | [Initial comparison](complementary-panel-results.md); [138-host extension](picard-extension-results.md) gives LF73_P1 eight net extra hits; LF110_P2's negative result is retained; same-source evidence only |
| Research FASTA-to-ranking | [Command available](research-fasta-ranking.md) for frozen host 4-mer and prevalence models; no model promotion |
| Kaneko 42-candidate comparison | [Complete](kaneko-baseline-results.md): both methods reach test Hit@2 of 42.86%; host 4-mer has higher point-estimate cost |
| BASEL/ECOR reference | [Fixed controls](basel-ecor-controls.md) and [bilinear comparison](paired-interaction-results.md) complete; remains an inspected cold-phage reference |
| Public-data intake | [Source review](external-evaluation-intake.md), [Kaneko workbooks](kaneko-workbook-intake.md) and [sequence audit](kaneko-sequence-intake.md) recorded; independent release validation pending |
| SDL integration | [PairPrediction exchange](model-interface.md) validates files and abstention; trained serving and retrieval-result admission remain pending |
| Two-route retrieval | [Supplied-score evaluator](offline-retrieval-rescue-evaluator.md) implemented; [benchmark design](retrieval-rescue-benchmark-v0.md) remains unfrozen |

The [results index](README.md#results-and-reproduction) links earlier annotation,
host 4-mer, additive, GenoPHI and evidence-coverage studies. Their protocols,
acceptance decisions, scores, weights and negative results remain frozen.

## Next work

Updated 2026-11-15. The matched 94-candidate comparison, research ranking command,
first Kaneko comparison and fixed complementary-panel extension are complete.

1. **Prepare a matched check of LF73_P1's marginal coverage.** The
   [11 rescued and three newly missed extension hosts](picard-extension-results.md#all-primary-changed-cases)
   give a concrete strain/candidate comparison with DIJ07_P2 and 536_P7 controls.
   Independent measurements or another laboratory/assay are still needed;
   retain LF110_P2's unfavorable result and source MLC's endpoint limits.
   Kaneko remains a separate 42-candidate task; its seven test hosts in five
   groups do not validate the named 94-library candidates.
2. **Finish source checks.** Resolve phenotype reuse/stock questions, assess
   contamination/completeness and compare near-duplicate hosts across sources
   before making broader claims. Host 127 and 38 phage assemblies remain
   unavailable; the latter are needed for later paired-sequence comparisons.
   Vitt's pair-level table remains unlocated.
3. **Reserve suitable validation inputs.** The reviewed panels do not supply an
   eligible independent cohort for the current 94 stocks. The 100 prospective
   source files remain quarantined pending identity, endpoint, reuse and
   partial-label review. Additional same-study hosts and independent
   laboratory/assay validation have different roles.
4. **Freeze a follow-up only when justified.** Specify inputs, comparison,
   resource limits and acceptance before fitting; select using validation only
   and retain all held-out failures. No new experiment is authorized by this
   roadmap alone.

Kaneko has its own grade-3 complete-clearing endpoint and 42 named candidates.
Its within-source results do not validate the earlier 94-output models.
The earlier 177-host cohort is also inspected development data: 16,597 observed
MLC labels, with 41 missing readings excluded. Neither cohort establishes
unseen-phage, cross-species or release performance.

## Candidate discovery and ownership

| Predictor owns | SDL owns |
| --- | --- |
| Model-specific representations, library encoding and preprocessing | Source library snapshots and HostProfile/HostState identity |
| Candidate scores, training, evaluation and calibration | Mechanism evidence, experimental context and provenance |
| Applicability and model artifacts | Union/decision policies, experiments and feedback |

Exhaustive pair scoring is sufficient for the current 94-phage library. New
encoders, neighbour indexes or rerankers need evidence of added value. The
[reference cases](reference-cases.md) record the rationale.

Later two-route retrieval would preserve the versioned union of evidence and
model nominations, including route provenance and missingness. A rescued
source-positive pair demonstrates ranking complementarity; mechanism discovery
requires separate evidence. Fusion, engineered variants and cocktail studies
are follow-up work.

Retrieval-to-SDL admission still needs a versioned result contract binding the
host/genome/context, library snapshot, per-candidate identity and route,
missingness, endpoint score or `not_estimated`, model/preprocessing hashes,
applicability, considered-pool counts, filtering reasons and K/tie policy.
The [existing pair exchange](model-interface.md) remains separate. Source
activity cannot be relabelled as `productive_infection`; see
[endpoint semantics](endpoint-semantics.md).

## Paused Evo2 branch

Decision, 2026-11-13: pause new Evo2 development and full-cohort extraction.
Bounded pilots establish access and cache integrity, with no strain-level
ranking result. The [source-verified preflight](https://github.com/edgarcuth-web/phage-host-interaction-predictor/issues/24#issuecomment-6488966629)
projects about 1.425 TB of intermediate responses for roughly 5.8 MB of final
host vectors. Issue #24 is closed as `not_planned`.

Protocols, pilots, scores and extraction code remain available for reproduction.
The dedicated preflight Actions workflow has been removed; the
[local runbook](nvidia-evo2-extraction.md) retains source verification and
zero-request planning. Revival requires a separately frozen, bounded experiment
and acceptance against inexpensive baselines. A zero-budget plan does not
authorize inference.
