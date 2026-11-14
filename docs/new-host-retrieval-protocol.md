# New-host fixed-library study: compatibility and freeze

Step 1 is complete. The [protocol](../data/studies/new_host_retrieval.v1.json)
was frozen before new fitting or scoring. The [step-2 comparison is complete](new-host-retrieval-results.md)
and retains its negative result. This is a development study on inspected data,
with no model promotion.

## What can be reused

The [compatibility ledger](../data/studies/new_host_retrieval.compatibility.v1.csv)
records both panels and every method. It extends the completed budget audit only
where the new comparison needs additional scores or training.

| Method | Primary grouped-host panel | B2 panel |
| --- | --- | --- |
| Train-only phage prevalence | Reconstruct from training observations | Reconstruct from training observations |
| Annotation proxy | Regenerate scores from the existing 117-host weights | Refit fixed C=0.1; matching weights/scores are not committed |
| Host 4-mer | Refit fixed C=0.1; matching weights/full scores are not committed | Same |
| GenoPHI | Reuse all 5,640 scores for 60 held-out hosts | Reuse all 5,828 scores for 62 held-out hosts |
| Phage-only, additive and bilinear sequence models | Refit under the host-only split with all 94 phages | Refit under the B2 split |

All **271 actual 136-dimensional sequence vectors** pass their assembly,
representation and vector-digest bindings. Reuse these raw features, then fit
scalers on the permitted training assemblies. No new genome download or feature
extraction is needed. These checks verify derived artifacts against their locked
receipts; they do not reread all source DNA.

The existing crossed sequence fit uses **65 named phages / 61 unique assemblies**
and 6,834 rows. Its full-cohort refit includes all 177 hosts. Neither supplies
matched weights for this new evaluation. Historical results remain unchanged.

## Frozen comparison

| Panel | Training | Parameter selection | Evaluation support |
| --- | --- | --- | --- |
| Primary | 117 hosts, all 94 phages, 10,957 observed rows | Existing 26 outer-validation hosts | 34 test hosts, all complete; include the all-negative host |
| B2 | 115 outer-training hosts, all 94 phages, 10,771 observed rows | Frozen inner split: 94 fit / 21 validation hosts; select L2 then refit all 115 | Score all 62 test hosts; 61 have complete labels |

The [B2 selection manifest](../data/studies/new_host_retrieval.b2_selection.v1.csv)
is derived from host-group identities only: 84 fit groups and 20 validation
groups. B2 test labels and primary-panel selected parameters do not select B2
parameters. There is no new seed search.

- Primary: all-complete-host macro Hit@2, bilinear minus train-only prevalence.
  Report annotation, host 4-mer, GenoPHI and additive/phage-only controls alongside it.
- Secondary: K=1,3,5, Recall/MRR, capped tests consumed, coverage and failures.
  Compute intermediate Hit@1–4 for the K=5 cost calculation. Pair AP is diagnostic.
- Pooled sequence models use the fixed L2 grid 0.1/0.01/0.001. Selection prefers
  validation Hit@2, then lower two-test consumption, then stronger L2. All new
  fits use the numerical environment recorded in the protocol.
- Keep all 94 named candidates and 41 source-unknown exclusions. B2 excludes
  ECOR53 whole from complete-panel metrics; it remains in score-support checks.
  Missing labels are never negative. Exact-tie expectations determine metrics.
- Use 2,000 paired host-group bootstrap draws. Record supported development
  hit gain and cost saving separately; claim joint benefit only if both frozen
  criteria hold. Retain unfavorable B2 and subgroup results.

GenoPHI's inherited protein-family vocabulary uses all cohort genomes without
labels. Its supervised fitting uses training labels only; retain the transductive
vocabulary limitation. Inherited artifacts keep their original environments;
new baseline refits are new-task artifacts, not bitwise historical reproductions.

## Check inputs without training

```bash
python scripts/check_new_host_retrieval.py
```

Optionally verify the already pinned annotation tables and save a new receipt:

```bash
python scripts/check_new_host_retrieval.py \
  --annotation-data-root artifacts/annotation_inputs \
  --out artifacts/new_host_retrieval_preflight.json
```

The checker uses the Python standard library. It rejects altered protocol/input
bytes, incompatible scopes, missing candidate scores and changed feature
identities. The [committed preflight](../data/studies/results/new_host_retrieval_preflight.v1.json)
also verifies the two annotation source tables. It performs zero fits and emits
zero new candidate scores. Step 2 is bounded to 21 pooled fits and at most 282
per-phage logistic fits, one CPU thread and no paid inference or GenoPHI retraining.
