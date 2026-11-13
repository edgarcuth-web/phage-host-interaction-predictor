# BASEL x ECOR published validation source: frozen intake v1

**Status:** source-verified *published, already inspected reference only*. This is **not** a new untouched release holdout, a model fit, an infection-mechanism target, or a candidate-ranking performance result.

## Why this source, and what it measures

Noonan et al., [GenoPHI (2026)](https://doi.org/10.1038/s41564-026-02482-5), describe a separate-laboratory evaluation with **52 BASEL phages** and **25 ECOR E. coli hosts**. The publication reports **1,300 attempted spotting interactions**, including 60 ambiguous phenotypes. This source is particularly relevant as a *published cold-phage reference*, not as a prospective evaluation created by this repository. The assay's binary clearance is **not equivalent** to measured productive infection.

The [pinned upstream CSV](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/data/experimental_validation/BASEL_ECOR_interaction_matrix.csv) provides **1,240 explicit binary rows**: 302 recorded positives and 938 negatives; 25 unique host identifiers and 52 unique phage identifiers. The 60 other matrix positions are **missing from the CSV** and must remain unknown (the paper reports 60 ambiguous results; we have not recovered all underlying source readouts).

Crucial evaluation limit: **only 3 of the 25 hosts have complete labels for all 52 phages**; 22 hosts are missing at least one label. The current fixed-panel evaluator excludes incomplete hosts. Reusing its original all-complete-host Hit@K primary metric without review would leave just **3** hosts, not 25. Do not fill unreported pairs with 0, subset candidates based on measured test labels, or silently change the denominator.

This reference was already published and the GenoPHI approach was considered in our earlier development work. The 52 BASEL phages are distinct study candidates, but **cold-phage accession/sequence disjointness from our fitted library has not yet been audited**. The 25 ECOR hosts may overlap host training context; no cold-host claim is made.

## Machine-readable receipt and CLI

[data/studies/basel_ecor_intake.v1.json](../data/studies/basel_ecor_intake.v1.json) pins the original repository commit, Git blob SHA-1, exact byte size, endpoint ID, expected source-label counts, missingness and validation status. It does **not** register the CSV as a new training dataset or assume a data redistribution license.

Run from the repository root after installing the package:

```bash
python scripts/inspect_basel_ecor.py --fetch --out artifacts/basel_ecor_intake_v1.json
# Or audit a locally obtained, unmodified source file:
python scripts/inspect_basel_ecor.py --source /path/to/BASEL_ECOR_interaction_matrix.csv --out artifacts/basel_ecor_local.json
```

The read-only runner verifies the original Git blob and length; checks exact schema, duplicate pairs and label coding; derives **all 60 missing pair identities**; records complete-host coverage and a SHA-256 of the bytes it read. It writes a fresh JSON report, but does not commit, republish or automatically relabel the source CSV. Local execution requires a network connection only for `--fetch`; synthetic tests run offline.

## Scientific gates before model evaluation

1. **Resolve source identity:** pin GenBank/NCBI accession tables for each BASEL phage and ECOR host, obtain and hash actual FASTA, and reconcile identifiers against existing Gaborieau/ECOR resources. Exact and near-duplicate checks must precede any cold-phage/cold-host claim. This intake includes **no** genome downloads.
2. **Review dataset rights and assay equivalence:** publisher code license does not automatically license all experimental data. Retain original spotting criteria and ambiguous readouts; do not map clearance to productive infection.
3. **Freeze partial-label evaluation:** the original full-panel evaluator admits only 3 complete hosts. Predeclare how known, unknown and unscorable labels affect top-K/coverage and uncertainty **before producing model scores**. An observed-only metric is not an unbiased full-panel Hit@K; it must not be presented as such.
4. **Implement source-compatible inference:** the existing additive host/phage baseline is an exploratory model, not evidence of generalizable pair-specific compatibility. No trained candidate scorer for this BASEL library is delivered here.
5. **Reserve separate independent release validation:** these known, published measurements are useful for development stress testing and an honest retrospective comparison. They are **not** a previously untouched final release cohort. A genuinely independent evaluation still requires its own source/lineage audit, protocol freeze, and retained failure results.

See [external evaluation intake](external-evaluation-intake.md), [endpoint semantics](endpoint-semantics.md), and [current retrieval roadmap](host-conditioned-phage-retrieval-roadmap.md).
