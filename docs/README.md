# Documentation

Start with the [project README](../README.md) for setup and current results,
or the [roadmap](host-conditioned-phage-retrieval-roadmap.md#next-work) for remaining work.
Frozen protocols and verified outputs live in [`data/studies/`](../data/studies/).

## Current commands and data

| Need | Documentation |
| --- | --- |
| Rank a host FASTA | [Research ranking command](research-fasta-ranking.md) |
| Exchange SDL files | [ModelInputBundle / PairPrediction](model-interface.md) |
| Fetch inputs or load weights | [Data guide](../data/README.md), [model cards](../models/README.md) |
| Interpret labels and missingness | [Endpoint semantics](endpoint-semantics.md), [dataset contract](dataset_contract.md) |
| Audit proposed external data | [Intake decisions and overlap gate](external-evaluation-intake.md) |
| Reproduce Kaneko intake | [Workbook audit](kaneko-workbook-intake.md), [sequence sources](kaneko-sequence-intake.md) |

## Results and reproduction

| Study | Documentation |
| --- | --- |
| Kaneko prevalence versus host 4-mer, 42 candidates | [Frozen protocol](kaneko-baseline-protocol.md), [results and reproduction](kaneko-baseline-results.md) |
| Matched new-host comparison, 94 candidates | [Protocol](new-host-retrieval-protocol.md), [results](new-host-retrieval-results.md) |
| Missed candidates and source evidence, 94 candidates | [Complete miss/exclusion ledger and reproduction](missed-candidate-audit.md) |
| Training-only complementary panels, two-test budget | [Results, named rescues and newly missed hosts](complementary-panel-results.md) |
| GenoPHI | [Demo and runbook](genophi-demo-reproduction.md), [aligned comparator](genophi-frozen-comparator.md), [seed sensitivity](genophi-seed-sensitivity.md) |
| Screening budgets from existing scores | [Hit/cost audit](screening-budget-audit.md) |
| Known-library retrieval | [Annotation baseline](host-retrieval-benchmark.md), [host 4-mer comparison](genome-host-retrieval.md) |
| First annotation baseline | [Frozen task](first-model-task.md), [protocol](host-baseline-protocol.md), [results](host-baseline-results.md) |
| Additive host–phage baseline | [Protocol](paired-benchmark-protocol.md), [split freeze](paired-split-freeze.md), [results](paired-baseline-results.md) |
| Bilinear sequence interaction | [Results, weights and reproduction](paired-interaction-results.md) |
| BASEL/ECOR reference | [Source intake](basel-ecor-published-reference.md), [identity and protocol freeze](basel-ecor-identity-protocol-freeze.md), [sequence audit](basel-ecor-sequence-audit.md), [fixed controls](basel-ecor-controls.md) |
| Host controls | [Similarity](host-similarity-evo2.md), [window sampling](host-window-sampling.md) |
| Evidence coverage | [Findings](annotation-recall-audit.md), [scripts](../studies/annotation_recall_audit/README.md), [estimand reconciliation](estimand-reconciliation.md) |

## Design references

- [Reference cases](reference-cases.md): source evidence and project decisions.
- [Offline retrieval-rescue evaluator](offline-retrieval-rescue-evaluator.md):
  implemented comparison of supplied score files.
- [Retrieval-rescue benchmark v0](retrieval-rescue-benchmark-v0.md): unfrozen
  scientific design; the evaluator does not make it a completed study.

## Paused Evo2 work

The [pause decision](host-conditioned-phage-retrieval-roadmap.md#paused-evo2-branch)
keeps these records available for reproduction. Full-cohort comparison remains
unexecuted; this work is not a prerequisite for the current ranking command.

- [Published-score reproduction](evo2-retrieval-reference.md).
- [API probe](nvidia-evo2-api.md) and [bounded host pilot](nvidia-evo2-host-pilot.md).
- [Extraction runbook](nvidia-evo2-extraction.md) and
  [unexecuted host-representation protocol](host-similarity-evo2.md).
