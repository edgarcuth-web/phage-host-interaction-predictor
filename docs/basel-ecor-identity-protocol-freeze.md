# BASEL–ECOR identity audit and frozen retrospective scoring protocol

**Status (2026-11-13): protocol frozen; sequence-identity gate NOT cleared; no model scored under it.**

This is a **published, already-inspected reference** from Noonan et al., [GenoPHI](https://doi.org/10.1038/s41564-026-02482-5), **not** an untouched release-validation set. The [source intake](basel-ecor-published-reference.md) pins 25 host IDs, 52 BASEL phage IDs, 1,240 observed binary spot-clearance rows (302 positives/938 negatives) and 60 unreported/ambiguous pairs. We do not turn unknown pairs into negatives or source clearance into productive infection.

## Identity check already completed

The [machine-readable identity ledger](../data/studies/basel_ecor_identity_ledger.v1.json) was derived directly from the pinned upstream GenoPHI CSV (Git blob SHA-1 `dcbbca4b4e4b1dd0bc543c5e2348189813e55bb0`) and the local 177-host / 94-phage [development resource lock](../data/studies/genophi_host_baseline.resources.v1.json) (Git blob SHA-1 `58ba0786f7cf836ed6146316f865d5c4fc322a9e`).

| Question | Evidence | Status |
| --- | --- | --- |
| Is this a cold-host evaluation? | **25/25** BASEL–ECOR host identifiers match host IDs in the earlier training resource lock; host file hashes and source URLs are recorded | **No**; these are previously represented ECOR hosts |
| Are the 52 phages new by *name*? | **0/52** BASEL phage IDs exactly match the 94 existing library phage IDs | **Name-only** difference; does not prove sequence novelty |
| Have we established the exact host assembly matches the paper's accession? | Training resources contain raw FASTA file SHA-256, but the paper points to ECOR assemblies under BioProject PRJNA230969 | **Not verified**; do not confuse file hashes and normalized sequence identity |
| Can we claim cold-phage sequence independence? | Per-Bas phage GenBank accessions, accession versions, complete genome sequences and joint similarity edges are not present in this intake | **Blocked**, no such claim yet |

The GenoPHI article explicitly calls this validation a study of new phages on hosts **present in training**. It points to GenBank accessions in Supplementary Table 8 and ECOR assemblies in Supplementary Table 9. The original [BASEL study](https://doi.org/10.1371/journal.pbio.3001424) provides an accession list (original BASEL genome range MZ501046–MZ501113) in S5 Table. **Those accession ranges do not establish a specific BasXX-to-accession mapping**. Do not fill accession fields from numerical order or guessed offsets.

Run the reproducible identity check using the exact published input:

```bash
python scripts/audit_basel_identity.py \
  --fetch --out artifacts/basel_ecor_identity_local.json
# Alternatively use --source /path/to/BASEL_ECOR_interaction_matrix.csv
```

The runner verifies the original CSV bytes, reads the locked existing training manifest, writes a new non-overwriting receipt and exits with status **2** by design: the phage-sequence gate is still blocked. The checked-in identity ledger describes the reviewed source **identifier join**; it is not a validated FASTA comparison.

## Frozen evaluation decisions

The immutable [v1 protocol](../data/studies/basel_ecor_evaluation_protocol.v1.json) was committed separately in **`ea23ba73fea847760908cc8db4c5b7e469a6b24e`**, before the subsequent identity-audit implementation and without fitting, scoring or selecting a model for this source.

- **Regime:** potential *known-host / cold-phage* retrospective evaluation after identity checks. Not cold-host, double-cold, independent release, mechanism prediction, or treatment-effect validation.
- **Panel:** all 25 host IDs × all 52 BASEL phage IDs. All candidate IDs must receive a source-compatible finite score or an explicitly reported abstention. The panel is never reduced to pairs with available labels.
- **Primary estimand:** equal-host-weight macro **Hit@2** as a *partial identification interval* over the **25 hosts**. Secondary K=1,3,5. All-negative or partly-unknown hosts remain in the denominator.
- **Missingness:** for each *fixed* ranking, lower Hit@K treats all unreported outcomes as negative **solely to calculate the mathematical bound**; upper Hit@K treats them as positive **solely for the upper bound**. Neither assignment changes the source data or asserts an observed label. Exact score ties are averaged over all tie permutations using the repository's existing retrieval convention.
- **Interpretation:** an interval [0.6, 0.8] is **not a 95% confidence interval**. It describes what the published, incomplete label panel can and cannot identify. Do not report the midpoint as measured performance.
- **Membership:** when naming concrete top-K candidates, rank by score descending, then lexicographic phage ID; do not mix this deterministic set membership with analytical tie-averaged metrics.
- **Comparators:** uniform-random ranking as a minimal control; a sequence-only control after confirmed sequences. If original GenoPHI scores for this exact readout and panel are retrievable, compare on identical pairs/denominators. Do **not** apply train-only *phage-specific* prevalence to wholly unseen phages.
- **Statistics and deployment:** record score coverage, source unknown exposure, intervals and unfavorable results. Bootstrap/group uncertainty is conditional on documented genomic groups; if those groups cannot be verified, report no confidence interval. **No model promotion from this inspected reference alone.**

The pure function `phage_cocktail.partial_hit_bounds.bounded_hit_at_k` implements these exact-tie lower/upper identification bounds for synthetic tests; it is not a trained predictor and does not fetch source labels, generate scores or claim an empirical Hit@K outcome.

## Remaining hard gate before running predictions

1. Manually reconcile all 52 `BasXX` names to the **published** phage accession and accession version, using Supplementary Table 8 and/or the BASEL original study's S5 Table.
2. Retrieve/verify those genome bytes and 25 ECOR assemblies, keeping raw FASTA SHA-256 **separate** from normalized sequence identity. Record accession updates and conflicting variants instead of silently remapping them.
3. Compare all 52 sequences jointly with the 94 locked training phages for exact duplicates and genomically near-identical groups under the **existing predeclared minimap2 2.31 rule** (95% identity and 85% aligned coverage of **both** genomes, with predeclared 90/99% sensitivity). Record ambiguous/fragmentary assemblies and exclusions without looking at test labels.
4. Review phage sourcing and data licenses; distinguish access to source data from permission to redistribute genomes or raw phenotype tables.
5. Only after these identity checks have been reviewed, fit or load a frozen model **without** using the BASEL outcomes for selection. Evaluate on the full 25 × 52 score matrix with the above bounds and failures retained.

If mapping or genomic review fails, keep this protocol frozen and mark execution blocked. A revised question needs a newly versioned protocol, not an edit to v1.
