# Picard extension: fixed panels on new host groups

**DIJ07_P2 → LF73_P1 finds a source-positive candidate for 108/138 admitted
hosts, versus 100/138 for prevalence: 11 rescues, three new misses, and unchanged
mean two-test consumption.** This supports replication within the same primary
study. It is not independent or blind validation. The secondary
536_P7 → LF110_P2 panel loses nine net hits.

No model was fitted, no candidate was reselected, and no genome model was scored.
This tests the named candidates from the
[earlier comparison](complementary-panel-results.md). It does not restore the
unsupported low-prevalence hypothesis or promote a serving default.

## Admission before scoring

The [protocol](../studies/picard_extension/protocol.json) and genome roster were
committed as [ec5c8c3](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/ec5c8c3a8af23cce7d5f60309f7f402668923029).
The admission implementation was committed as
[bd41417](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/bd414170b26cb6ae564c86920d6be67fdd7bfa5d)
before the full audit. After independent admission checks, the cohort and
evaluation code were committed as
[11cdac4](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/11cdac4b0135d473ca3b8ea794f3783bd34cabd3)
before extension scoring. The public matrix had already been inspected; these
freezes do not make it a blind holdout.

| Sequential admission step | Remaining records | Removed at this step |
| --- | ---: | ---: |
| Original 402 matrix hosts minus the frozen 177-host crosswalk | 225 | — |
| Source-labelled phylogroups A, B1, B2, C, D, E, F, G | 186 | 39 |
| All 94 candidate readings observed | 179 | 7 |
| Component contains no original-cohort host | **138 in 114 groups** | 41 |

Across all 225 records, 13 panels contain 115 unknown readings, 43 share a
component with the old cohort, and 39 are outside the taxon scope. These reasons
overlap. Before taxon filtering, 212/225 panels are complete. The
[admission ledger](../data/studies/results/picard_extension.v1/admission/host_admission.csv)
retains every record, source identity, sequence digest and blocking reason.

All 225 additional FASTAs from Figshare article 25941691 version 1 were checked
against file IDs, names, sizes and MD5s: 1,144,096,161 bytes in total. DNA parsing,
raw SHA-256, canonical assembly digests, contig counts, N50 and ambiguity counts
were verified. Canonical hashing ignores headers, contig order and strand, but
not fragmentation. The old 177 file identities still match their resource lock;
their previously verified canonical digests were reused without redownloading
those genomes. No new assembly exactly duplicates an old canonical assembly.

Grouping uses the existing strict core-distance threshold `< 0.0001`, exact
canonical matches and transitive components. All 404 distance-table nodes are
kept, including `H27` and `H1-005-0065-L-P`, which lack matrix rows and never enter
evaluation. This retains possible unmeasured bridges. The old 147-group relation
reproduces exactly, and those groups are not merged by the expansion.

Names join only by unique uppercase/hyphen-removal aliases. Taxon scope is from
source metadata, not an independent sequence-based species assignment. Generic
FASTA contig headers do not independently bind Gembase IDs; that mapping relies
on the published metadata/filename crosswalk. Source core distances were checked
and reused, not recalculated from all assemblies. N50 and ambiguity do not
certify completeness, contamination freedom or physical stock identity.

## Same two-test budget

The primary pair remains **DIJ07_P2 → LF73_P1**, originally selected on 117
training hosts. Its prevalence control is DIJ07_P2 then 536_P7. The secondary
pair remains **536_P7 → LF110_P2**, selected on 115 non-B2 training hosts.
Prevalence is reconstructed only from the respective original training labels,
checked against all 16,597 old observations. All 94 candidates remain ranked;
exact-score ties are averaged.

All four source-all-negative admitted hosts remain in the denominator: FN-B4,
FN-B7, H1-002-0060-C-T and ROAR205. Unknown readings elsewhere remain unknown.

| Policy on the same 138 hosts | Hit@2 | Mean capped tests |
| --- | ---: | ---: |
| Primary prevalence: DIJ07_P2 then 536_P7 | 72.46% (100/138) | 1.333 |
| Primary fixed pair: DIJ07_P2 then LF73_P1 | **78.26% (108/138)** | **1.333** |
| Secondary prevalence: tied 536_P7 / DIJ07_P2 | 72.46% (100/138) | 1.402 |
| Secondary fixed pair: 536_P7 then LF110_P2 | 65.94% (91/138) | 1.471 |

The prevalence policies select the same unordered pair at two tests; their
original training rates give different first-test order/tie expectations and
thus different costs. One test is one equal-cost unit. Stop at the first
source-positive result or after two unsuccessful tests: capped cost is
`2 - Hit@1`.

The primary net gain is **8/138 = +5.80 percentage points**, conditional 95%
interval **[+0.70, +11.51]**, meeting the protocol's supportive same-source
replication rule. Consumption is unchanged because both primary policies start
with DIJ07_P2. Eleven rescues occupy ten genomic groups; three new misses occupy
three groups.

| Prespecified stratum: primary pair versus prevalence | Hosts / groups | Rescued / newly missed | Hit@2 difference, percentage points |
| --- | ---: | ---: | ---: |
| All admitted | 138 / 114 | 11 / 3 | +5.80 [+0.70, +11.51] |
| Source B2 | 32 / 27 | 2 / 0 | +6.25 [0.00, +16.67] |
| Source non-B2 | 106 / 87 | 9 / 3 | +5.66 [−0.98, +12.82] |

Both subgroup intervals include zero; the pooled gain does not establish benefit
in every phylogroup. The primary pair's training included B2 hosts, so its B2
extension is new-host-group evaluation, not unseen-phylogroup training.
Intervals use 2,000 paired whole-component draws, seed 20261115, and are
pointwise and conditional on the fixed policies/source. Prior inspection,
assay dependence, grouping limitations and candidate-selection uncertainty
remain. Do not pool panels, methods or strata as independent replications.

The secondary pair rescues eight but newly misses 17, net **−9/138**. Its Hit@2
difference is −6.52 points [−15.00, +1.39], and its mean cost difference is
+0.069 tests [+0.026, +0.115]. Its B2 subset rescues none and newly misses two.
These failures remain in the result tables.

## All primary changed cases

These are original aggregate MLCs, not independently verified productive
infection. No raw image or replicate annotation was rescored in this extension.

| Source host | Source phylogroup | DIJ07_P2 MLC | 536_P7 MLC | LF73_P1 MLC | Change |
| --- | --- | ---: | ---: | ---: | --- |
| ECOR-03 | A | 0 | 0 | 1 | Rescued |
| ECOR-08 | A | 0 | 0 | 2 | Rescued |
| IAI22 | A | 0 | 0 | 3 | Rescued |
| IAI28 | A | 0 | 0 | 2 | Rescued |
| IAI29 | A | 0 | 0 | 3 | Rescued |
| IAI37 | A | 0 | 0 | 3 | Rescued |
| IAI40 | A | 0 | 0 | 4 | Rescued |
| IAI53 | B2 | 0 | 0 | 4 | Rescued |
| IAI67 | B2 | 0 | 0 | 4 | Rescued |
| PDP351 | C | 0 | 0 | 3 | Rescued |
| ROAR419 | B1 | 0 | 0 | 1 | Rescued |
| ECOR-07 | A | 0 | 2 | 0 | Newly missed |
| IAI10 | A | 0 | 1 | 0 | Newly missed |
| LM08 | D | 0 | 3 | 0 | Newly missed |

[Changed-host evidence](../data/studies/results/picard_extension.v1/evaluation/changed_host_evidence.csv)
includes original row/column links, candidate rank intervals and successful
controls for both policies. These cases support prioritizing LF73_P1 for a
matched follow-up check; the three lost hosts preclude a universal replacement
claim. LF73_P1 is relatively prevalent, and a receptor/defense explanation has
not been established.

## Artifacts and reproduction

The [study directory](../studies/picard_extension/) contains the protocol, public
genome-source roster, input/implementation locks and
[cohort seal](../studies/picard_extension/cohort_freeze.json). The
[result directory](../data/studies/results/picard_extension.v1/) retains:

- `admission/`: all 225 eligibility records, the 138-host roster, 404 graph
  memberships, 183 edges, 115 unknown pairs and independent QA.
- `evaluation/`: 552 host-policy rows, 276 transitions, 1,104 source-evidence
  records, train-only priors, marginal coverage, subgroup comparisons and QA.
- `PUBLICATION.json`: code/result identities and reproduction checks.

Python 3.11+ standard library suffices. Fresh fetching obtains about 1.14 GB of
new assemblies plus three source tables. Changed caches are rejected; raw
genomes and the full source matrix stay outside Git. Output directories must
be new:

```bash
python -m unittest discover -s tests -p test_picard_extension.py -v
python studies/picard_extension/audit_cohort.py --fetch \
  --out artifacts/picard_extension_admission
python studies/picard_extension/verify_admission.py \
  --admission artifacts/picard_extension_admission
python studies/picard_extension/evaluate_panels.py \
  --admission artifacts/picard_extension_admission \
  --out artifacts/picard_extension_evaluation
python studies/picard_extension/verify_evaluation.py \
  --admission artifacts/picard_extension_admission \
  --evaluation artifacts/picard_extension_evaluation
```

Independent verifiers import neither runner. They rehash 225 assemblies,
reconstruct components by graph traversal, check all admission rules and verify
12,972 target cells, 16,597 old observations, all metrics/evidence and 24,000
paired interval draws. All six admission-runner outputs and nine evaluation
outputs reproduce byte for byte. A changed host outcome is rejected; existing
results cannot be overwritten.

Attribution: Gaborieau, Vaysset, Tesson et al., DOI 10.1038/s41564-024-01832-5;
primary archive Zenodo 10.5281/zenodo.13831957; assemblies Figshare
10.6084/m9.figshare.25941691.v1. The
[reuse review](../data/studies/paired_data_reuse_review.v1.json) records archive
byte correspondence and CC BY 4.0 declarations. This does not clear unrelated
sources or establish independent validation, cocktail efficacy or treatment
benefit.
