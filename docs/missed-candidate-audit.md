# Missed-positive candidates and source evidence

The saved 94-candidate rankings contain **179 confirmed ranking misses across
47 distinct hosts**. A miss means at least one observed source-positive candidate
exists in the library, but the first two recommendations find none. These are
overlapping, inspected development comparisons, not independent observations or
a newly validated improvement.

## Scope and classifications

The audit retains the primary validation/test partitions (26/34 hosts) and B2
test partition (62 scored hosts). These are 122 host-panel occurrences involving
93 distinct hosts. Seven deterministic methods give 854 decisions; random-score
expectations remain separate. All 177 source hosts and 94 candidates are checked
against the primary MLC table: 16,597 observed readings and 41 unknowns.

| Classification | Evidence and scope |
| --- | --- |
| Ranking miss | 179 complete-panel host-method records; observed positive exists, exact-tie Hit@2 is zero |
| No source-positive candidate | NILS24 and DEC2A account for 14 scored records. NILS22 also has 94 observed zeros but has no target scores in this comparison |
| Incomplete source panel | 370D has 39 missing MLC readings; ECOR53 has two. Missing MLC remains unknown despite the subset's original zero coding |

ECOR53 has 16 known positives and a known top-two hit under all seven methods.
It remains excluded from the original complete-panel metrics; an incomplete
panel does not make its top-two outcome unknown. No complete-panel record has
a partial Hit@2 probability here, but all exact-score ties and rank intervals
are retained.

## Requested cases

| Host | Source positives | Concrete missed positive |
| --- | ---: | --- |
| NILS53 | 11 | NIC06_P2, MLC 1.0, ranks fourth in primary GenoPHI. All seven methods miss at two in both panels |
| NILS78 | 7 | LF82_P9, MLC 1.0, ranks fifth in primary GenoPHI. B2 phage-only/additive successfully select LF82_P8 second; these controls remain |
| AL505 | 41 | LF82_P2, MLC 4.0, ranks fourth in B2 host 4-mer; its first two candidates both have MLC zero |

[Focus candidates](../data/studies/results/missed_candidate_audit.v1/focus_positive_candidates.csv)
lists all 59 positive pairs and the prevalence/4-mer/GenoPHI rank intervals in
both panels. For example, NILS53's NIC06_P2 has prevalence ranks 16–17, not a
unique rank implied by lexicographic display order.

## Evidence and retained outputs

The [result directory](../data/studies/results/missed_candidate_audit.v1/) contains
every deterministic host decision, all miss/exclusion records, successful focus
controls, positive-candidate ranks and training-only support, complete candidate
rankings, all 41 unknowns, and 1,478 deduplicated source-evidence records. Join
tables by `evidence_id = host::phage`. Source records bind original host names,
matrix row/column, MLC, pinned URL, raw-annotation line numbers and image names.

NILS53 and NILS78 have nine public annotation rows per library pair. AL505 has
no rows in the public raw-annotation table, so its evidence stops at the aggregate
matrix. Six missing aggregate MLCs have at least one raw `score=1` annotation;
these remain unknown rather than being relabelled. Images were not independently
read. Source MLC positivity is not a productive-infection or treatment endpoint.

Full source tables and raw-annotation extracts are local reproduction inputs or
outputs, outside Git. The
[source-reuse review](../data/studies/paired_data_reuse_review.v1.json) records the
primary archive's CC BY 4.0 declaration and MLC byte correspondence; it does not
establish separate terms for the GenoPHI subset or the raw annotations. The
[input manifest](../studies/missed_candidate_audit/input_manifest.json) pins URLs,
sizes, SHA-256 values and Git blob identities. Primary data attribution is to
Gaborieau, Vaysset, Tesson et al. (DOI 10.1038/s41564-024-01832-5; Zenodo
10.5281/zenodo.13831957), with cohort membership from Noonan et al., GenoPHI.

## Reproduce and verify

Python 3.11+ standard library suffices. Run from the repository root with a new
output directory. `--fetch` obtains approximately 16 MB of pinned public tables
and saved outputs. Without it, missing inputs fail. No model is fitted or inferred.

```bash
python studies/missed_candidate_audit/build_audit.py --fetch \
  --input-root data/raw/missed_candidate_audit \
  --out artifacts/missed_candidate_audit
python studies/missed_candidate_audit/verify_audit.py \
  --input-root data/raw/missed_candidate_audit \
  --out artifacts/missed_candidate_audit
```

The runner checks 11 inputs and reproduces 4,840 saved Hit@1–5 values (maximum
error 1.11e-16). The independent cutoff-based verifier imports no runner helpers;
it checks 854 decisions, 80,276 candidate rank intervals and 11,979 locally
exported raw annotation rows. QA receipts and the
[publication manifest](../data/studies/results/missed_candidate_audit.v1/PUBLICATION.json)
bind the checked code and retained outputs. The next coverage-selection
hypothesis remains untested; this audit neither changes the frozen models nor
establishes a biological cause of their misses.
