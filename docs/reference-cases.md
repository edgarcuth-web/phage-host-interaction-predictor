# Reference cases and project decisions

Reviewed 2026-11-13. This note connects selected primary sources to engineering
choices. The project decisions below are our interpretation of those sources.
Performance comparisons belong in the versioned study reports, with matched
labels, candidate libraries and denominators.

## Cases

| Case and evidence | What to borrow | Applicability and project decision |
| --- | --- | --- |
| **GenoPHI**, Noonan et al., *Nature Microbiology* (2026). [Paper](https://doi.org/10.1038/s41564-026-02482-5), [code](https://github.com/Noonanav/GenoPHI). Peer reviewed. | Protein-family features from both genomes, feature selection and CatBoost; evaluate unseen-host, unseen-phage and double-cold tasks separately. | Retain the aligned comparator. Its candidate-set coverage analysis does not establish mixed-cocktail efficacy. Our [comparator](genophi-frozen-comparator.md) and [seed sensitivity](genophi-seed-sensitivity.md) determine the local evidence, including failures. |
| **Escherichia predictor and recommender**, Gaborieau et al., *Nature Microbiology* (2024). [Paper](https://doi.org/10.1038/s41564-024-01832-5), [code](https://github.com/mdmparis/coli_phage_interactions_2023). Peer reviewed. | Known-library recommendations for new isolates, strain-grouped evaluation and experimental testing of recommended combinations. | Use the deployment question and prospective-isolate validation as references. Preserve the source assay endpoint and distinguish E. coli-specific annotations from generally available features. |
| **PhageHostLearn**, Boeckaerts et al., *Nature Communications* (2024). [Paper](https://doi.org/10.1038/s41467-024-48675-6), [code](https://github.com/dimiboeckaerts/PhageHostLearn). Peer reviewed. | Pair phage RBP and host K-locus protein representations with a supervised model; evaluate ranked candidates experimentally. | Reference for focused biological representations. Its capsule-centered Klebsiella setting requires a separate feature/data assessment before transfer to E. coli. Missing annotations remain unknown. |
| **Klebsiella interaction atlas**, Selvakumar et al. (2026). [Preprint v1](https://www.biorxiv.org/content/10.64898/2026.08.12.744533v1). Not peer reviewed. | Connect an interaction atlas and model-selected mixtures to separate in-vitro and animal measurements. | Keep ranking, mixture suppression and in-vivo outcomes as distinct endpoints. The preprint reports that in-vitro suppression did not consistently translate into reduced colonization. |
| **Locus high-throughput cocktail platform**, Penke et al., *Nature Communications* (2026). [Paper](https://doi.org/10.1038/s41467-026-68684-x). Peer reviewed. | Standardized automated phenotypic measurements, run metadata and iterative assessment of candidate cocktails. | Apply measurement/provenance and feedback patterns in SDL. The workflow description alone does not specify a reproducible Bayesian or active-learning acquisition policy. |
| **T7 multiobjective modeling**, Novy et al., *Cell Systems* (2026). [Paper](https://doi.org/10.1016/j.cels.2026.101712). Peer reviewed. | Model variant-specific measurements across host contexts and validate new variant predictions experimentally. | Later-stage reference requiring variant-level sequence/fitness data. The current source-pair matrix is a different training target. |
| **Omnii-1 phage**, Radical Numerics (2026-10-20). [Company technical disclosure](https://www.radicalnumerics.ai/blog/omnii-phage). | Link compatibility prediction, candidate prioritization and experimental feedback; track recovery and host activity separately. | Longer-term system reference. Administration under compassionate use is company-reported; patient outcomes are undisclosed. The article supplies no transferable validation of this repository's model or Evo2 representation. |

## Decisions for this repository

1. **Screening budgets using existing outputs (B): completed.** The
   [separately frozen audit](screening-budget-audit.md) reuses full scores and
   published points under explicit availability rules, without a new fit.
   Grouped-host test cost savings coexist with B2 failure and mixed validation
   results. Unsupported curves remain unavailable.
2. **Then: diagnose the remaining modeling gap (A, if needed).** Per-phage host
   models already condition on host features. An explicit residual model is an
   optional hypothesis, not a prerequisite for host-conditioned prediction.
   Adding the same host-only scalar to every candidate cannot change that
   host's phage ordering.
3. **Later: study candidate-set diversity/coverage (C).** Use training-only
   interaction profiles or declared genomic features for diversity policies.
   Coverage and multiple-positive selection do not establish independent
   mechanisms, resistance prevention or mixture efficacy. SDL owns experimental
   measurements and feedback; see the [ownership boundary](host-conditioned-phage-retrieval-roadmap.md#candidate-discovery-and-ownership).

The [Evo2 pause](host-conditioned-phage-retrieval-roadmap.md#paused-evo2-branch)
remains in effect. These cases motivate task-specific data and evaluation;
they provide no local evidence for restarting full-cohort extraction.

## Rules for the budget audit

- Preserve the fixed 94-phage library and source endpoint. Incomplete-label
  hosts remain excluded whole from primary metrics; complete all-negative hosts
  remain included. Report eligible-host summaries separately.
- Report hit probability and test consumption together. Under a stop-at-first-hit
  policy capped at K, an unsuccessful host consumes K tests. Preserve the
  existing exact-tie metric convention; MRR and mean first-hit rank are distinct.
- A first-hit median can already be one test. Freeze a useful cost/coverage
  contrast rather than requiring an impossible one-test median reduction.
- Reusing inspected seed-42/43/44 and B2 panels remains development analysis.
  Grouped uncertainty does not create independent release validation, and an
  inconclusive gain does not establish absence of useful host information.

Keep detailed definitions in [endpoint semantics](endpoint-semantics.md),
[retrieval metrics](host-retrieval-benchmark.md#metrics-and-denominators) and
[estimand reconciliation](estimand-reconciliation.md). Any new study needs its
own versioned protocol; this literature note does not amend frozen studies.
