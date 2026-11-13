# Documentation

For current progress and setup, start with the [project README](../README.md).
Study reports below retain methods, results and limitations; machine-readable
protocols and verified outputs live in [`data/studies/`](../data/studies/).

## Results and reproduction

| Study | Documentation |
| --- | --- |
| First observed-label annotation baseline | [Protocol](host-baseline-protocol.md), [results](host-baseline-results.md) |
| Additive host–phage baseline | [Protocol](paired-benchmark-protocol.md), [split freeze](paired-split-freeze.md), [results](paired-baseline-results.md) |
| Known-library retrieval | [Annotation baseline](host-retrieval-benchmark.md), [host 4-mer comparison](genome-host-retrieval.md) |
| GenoPHI | [Demo and runbook](genophi-demo-reproduction.md), [aligned comparator](genophi-frozen-comparator.md), [seed sensitivity](genophi-seed-sensitivity.md) |
| Host controls | [Similarity](host-similarity-evo2.md), [window sampling](host-window-sampling.md) |
| Evo2 published score reproduction (paused branch) | [Reference results](evo2-retrieval-reference.md) |
| Hosted Evo2 (paused branch) | [API probe](nvidia-evo2-api.md), [bounded host pilot](nvidia-evo2-host-pilot.md), [extraction runbook](nvidia-evo2-extraction.md) |
| Evidence coverage audit | [Findings](annotation-recall-audit.md), [reproduction scripts](../studies/annotation_recall_audit/README.md), [estimand reconciliation](estimand-reconciliation.md) |

## Interfaces and data

| Topic | Documentation |
| --- | --- |
| Task, labels and acceptance | [First model task](first-model-task.md), [dataset contract](dataset_contract.md), [endpoint semantics](endpoint-semantics.md) |
| SDL file exchange | [ModelInputBundle / PairPrediction](model-interface.md) |
| Data and development weights | [Inputs and rebuilds](../data/README.md), [model cards](../models/README.md) |
| Supplied-score comparison | [Implemented offline retrieval-rescue evaluator](offline-retrieval-rescue-evaluator.md) |

## Design references

[Reference cases and project decisions](reference-cases.md) records seven source
cases, their evidence status and applicability, and the budget-audit-first
priority. Detailed results remain in the study reports above.

## Proposed work

[Host-conditioned retrieval roadmap](host-conditioned-phage-retrieval-roadmap.md)
records the objective and remaining engineering work.
[Retrieval-rescue benchmark v0](retrieval-rescue-benchmark-v0.md) is an unfrozen
scientific design; the implemented evaluator does not make it a completed study.
The Evo2 host-representation extension is [paused](host-conditioned-phage-retrieval-roadmap.md#paused-evo2-branch)
and is no prerequisite for retrieval development. Its unexecuted protocol and
hardware handoff remain in [host similarity and Evo2](host-similarity-evo2.md).
