# Retrieval-rescue benchmark v0 — design for review

**Status: proposed benchmark specification, NOT a frozen experiment or reported
result.** This does not implement or train a new model. Its purpose is to make
the next implementation testable without silently reusing the already-inspected
development results as independent confirmation. Review and version the protocol,
source licensing, split manifests and hashes **before** any new model fit.

## Question and hypothesis

Given a query bacterial genome and a fixed, versioned collection of existing
phage genomes, does **sequence-conditioned learned retrieval** find more
source-assay-positive phages at a fixed top-K budget than the available
annotation/evidence retrieval route?

The key additional question is whether the learned route **rescues source-positive
phages the evidence route fails to rank in its top K**. This is *ranking
complementarity*, not evidence that a previously unknown causal receptor or gene
has been discovered.

Reject the assumption that unknown receptor, defense, RBP, or counter-defense
annotations are hard-negative infection evidence. A screened-out true positive
cannot be recovered by a downstream reranker unless a parallel discovery route
can nominate it.

## Stage 0: existing development benchmark (no new claims)

- Library: the *same* 94 named phages in the pinned E. coli GenoPHI subset.
  Confirm each ID and any exact-sequence duplicates; keep the historical panel
  for comparability and do not deduplicate by label.
- Hosts: the existing 177-host manifest. The frozen observed readout has 16,597
  audited MLC > 0 / non-positive labels and 41 excluded missing-primary-readout
  pairs. Preserve exact row-level unknowns and source identifiers.
- Primary partition: the existing genomic-group split with 117 train,
  26 validation and 34 test hosts. All 34 primary test hosts have complete
  94-member label panels; 33 have one or more source positives.
- Status: *already-inspected development comparison*. No independently
  confirmatory gains or model promotion can be inferred from it.
- Label: **source-defined MLC positivity**. Do not convert it to adsorption,
  productive infection, growth suppression, resistance risk or treatment success.
- Transfer observations from existing representation studies are not evidence of
  full-cohort Evo2 feature completion. Issue #24 remains a separately frozen
  host-only study and no additional NVIDIA calls are authorized by this document.

The existing `src/phage_cocktail/retrieval.py` implements exact expected Hit@K,
Recall@K and MRR over within-score ties; it fails closed on incomplete score
support and excludes an entire host with missing labels. Reuse rather than
replace this evaluator in the first comparable analysis.

## Baselines and method labels (avoid overclaiming)

| ID | Candidate ordering | Interpretation |
| --- | --- | --- |
| P | Train-only per-phage prevalence | Known-library non-genomic control |
| A | Frozen annotation-based host features (surface categories and defense counts) | **Annotation proxy**; not a validated mechanistic exclusion rule |
| G | Genome-only learned host+phage sequence compatibility model | Proposed; *not yet trained* |
| AG | Combined annotation+genome features under a trained model | Proposed; evaluate incremental information |
| U | Predeclared union of the A and G top-ranked candidates | Proposed; complementarity at an *identical final K* |

A host 4-mer known-phage logistic baseline is also required as a strong,
inexpensive control on the existing panel. Do **not** label the annotation model
"mechanism-only retrieval" or present PHIST shared k-mers as a verified
receptor/defense screen. A genuine mechanistic-rule baseline needs a separately
documented, source-supported receptor/RBP/defense compatibility policy and
measured coverage of those inputs.

G must use **phage sequence information** if it is to support a claim of
generalization beyond phage identity. A learned host-only representation plus
per-phage identity coefficients is useful as a known-library comparator but is
NOT evidence for cold-phage retrieval. Untrained Evo2 vector cosine similarity
is not a validated compatibility score.

## Equal candidate-budget comparison

For each frozen host, score **every one of the 94 library phages**, including
unannotated or source-negative members, unless the entire host is excluded
because of missing source readout. Unknown labels must never change the
candidate library at inference time.

Predeclare primary K=2, secondary K=1,3,5. Report each method's:

- expected complete-panel host-macro **Hit@K**, eligible-host **Recall@K**,
  MRR, and both complete/eligible denominators;
- candidate pool size, unscorable/abstained count, duplicate accession/assembly
  identities and feature missingness;
- per-host, per-phage score plus the exact source endpoint and score producer;
- source-positive phage IDs selected only by G at K, only by A at K, and by both.

For the last "rescued positive" tabulation, a concrete top-K set requires a
**separate deterministic tie-break** (e.g. lexicographic phage ID) solely to
define membership, reported alongside ties. Do not silently mix those
deterministic set counts with the existing *analytical tie-averaged* Hit/Recall
numbers. Alternatively implement expected set-rescue under an explicitly
reviewed joint tie policy. If there are no pairwise labels for a host/library
combination, report unavailable, not zero.

For union U, a naive concatenation of A top K plus G top K uses **2K** screening
slots and is not an equal-budget comparison. The first candidate rule for a
reviewed v1 protocol should explicitly freeze A/G slot allocation, duplicate
handling, fill order and any score normalization *without looking at test
labels*. Include U's exact final K and compare against A, G and a strong
prevalence baseline. For K=1, a two-branch union cannot reserve a slot for each;
test a separately precommitted selection policy or report N/A.

The primary estimand should be specified and signed off before implementation:
a paired difference in **Hit@2** on complete-panel grouped-test hosts, alongside
the predeclared rescued-positive rate conditional on observed positives.
Paired confidence intervals must resample **host genomic groups**, not each
host–phage row independently; include all negative-only complete hosts in Hit.
Show all phylogroup failures and sensitivities; no new primary metric selection
after reading test results.

## What cannot be learned from this panel

The full 94-phage library includes the same named phages in training and
grouped-test retrieval. It is thus an **unseen-host, known-phage-library** test.
That panel cannot establish new phage-family or double-cold generalization.
A separate freeze must group by exact phage genome identity and phage family/
neighbourhood, and ensure there is no cross-partition host-genomic or phage-
sequence leakage. It needs both train coverage and held-out labels for each
claim, plus a completely untouched future release holdout.

The current dataset also cannot validate discovery of causally novel infection
genes. "Unexplained by current annotations" is an evidence state. Distinguish
missing annotation, low-confidence annotation, explicit supported negative
mechanism evidence, and independent wet-lab confirmation.

## Stage 1 implementation contract (next PR, not implemented here)

Build an **offline, provider-free evaluation adapter** which accepts:

- frozen host/library manifests with source hashes, endpoint, time and split;
- one observed-label table, with unknowns preserved;
- complete-panel per-pair score files from A, G, AG and P (or explicit absence);
- model+preprocessing hashes, applicability context and candidate provenance.

Outputs must include a validated per-host pairwise ranking table, aggregate
metrics, candidate membership and rescue audit, unresolved/unknown exclusion
ledger, group-bootstrap contrasts, hash-bound report and explicit scientific
limits. Fail on duplicate pair IDs, inconsistent library/split, nonfinite score,
mixed endpoints, changed manifests, test-derived training features and invalid
tie policy. A model output is an **endpoint score**, not a mechanistic
`InteractionPrediction`.

Suggested acceptance tests: zero rescue for identical deterministic orderings;
positive rescue when two orderings find disjoint observed positives; no rescue
credit for unknown labels; an all-negative host affects Hit but not eligible
Recall; collision-safe tied rankings; U never exceeds the specified K; no
model endpoint is recast as productive infection.

Only after the offline adapter and retrieval benchmark are accepted should a
subsequent PR define a versioned SDL learned-retrieval result, union policy and
scientific reviewer. No new repository, external-model calls or hardware
execution are necessary for Stage 1.

## Review checklist

- [ ] Approve exact input cohort/library manifests and source permissions
- [ ] Approve endpoint, missing-value treatment and applicability limits
- [ ] Approve baseline A as an **annotation proxy**, not a mechanistic oracle
- [ ] Approve equal-K union policy, tie conventions and rescue denominator
- [ ] Approve grouped-host and separate cold-phage/double-cold claim scopes
- [ ] Approve held-out-data policy and prefit hashes
- [ ] Only then version the experiment as a **frozen study** and implement
