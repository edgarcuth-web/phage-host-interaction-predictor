# Phage–Host Interaction Predictor

Given a bacterial genome and a versioned library of existing phages, the research
goal is to rank candidates for experimental screening. Current code provides
source-label baselines, fixed-library retrieval studies and an SDL file interface.
It does not yet provide a validated genome-to-candidate service or cocktail model.

## Current progress

- **Data:** 177 E. coli hosts, 94 phages and 16,597 observed source-MLC labels;
  41 missing primary readings are excluded. Missing observations stay unknown.
- **Evaluation:** completed annotation, host 4-mer, additive host–phage and
  GenoPHI comparisons, with grouped splits, exact-tie retrieval metrics and
  uncertainty estimates. No model has established the evidence needed for promotion.
- **Evo2 (paused):** published-score reproduction and bounded hosted pilots are
  retained as optional research records. Full-cohort extraction and further Evo2
  infrastructure are outside the current work plan; retrieval work can proceed
  independently. See the [decision](docs/host-conditioned-phage-retrieval-roadmap.md#paused-evo2-branch).
- **Integration:** ModelInputBundle / PairPrediction validation and abstention
  work. Trained baseline serving through that interface is pending. The offline
  two-route evaluator audits supplied score files; it does not train a retriever.

Test Hit@2 on the frozen 94-phage panel, including all complete-label hosts:

| Panel | GenoPHI | Host 4-mer | Train-only phage prevalence |
| --- | ---: | ---: | ---: |
| Primary seed 42, 34 hosts | 91.18% | 88.24% | 73.53% |
| B2 holdout, 61 hosts | 78.69% | 75.41% | 90.16% |

The [seed-43/44 sensitivity](docs/genophi-seed-sensitivity.md) retains positive
GenoPHI/prevalence point differences, but both gain intervals cross zero and the
registered robustness criterion is **not established**. These inspected development
results do not establish superiority on new data. See the
[full comparator report](docs/genophi-frozen-comparator.md).

The [existing-score budget audit](docs/screening-budget-audit.md) finds test-cost
savings on the grouped-host splits, with B2 failure and mixed validation results.

The [external-evaluation intake gate](docs/external-evaluation-intake.md) can audit
source-lineage and exact assembly overlap for proposed new data. This is a
conservative software check, **not** an independently validated holdout.

A [published BASEL/ECOR source intake](docs/basel-ecor-published-reference.md)
pins a 25-host / 52-phage spotting assay matrix. Only 3 hosts have complete
52-phage labels, so this is a development reference rather than a directly
usable independent full-panel Hit@K benchmark or a promoted model.
Its [frozen evaluation protocol](docs/basel-ecor-identity-protocol-freeze.md)
specifies Hit@2 missing-label bounds before model scoring. The subsequent
[sequence audit](docs/basel-ecor-sequence-audit.md) resolves all 52 accessions:
no training overlap at the declared 95% identity boundary, but eight candidates
share training components at 90%. All 25 host IDs are known, while their cited
assemblies differ from the training FASTAs. The subsequent
[v2 assembly policy and fixed controls](docs/basel-ecor-controls.md) use explicitly
bound sequence features: cosine Hit@2 is 4.00% versus random identification bounds
of 37.32%–44.06%. This negative composition control is not a validated recommender;
learned interaction models remain to be compared.

## Quick start

Use Python 3.11+ in an activated virtual environment, from the repository root:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -e ".[dev]"
pytest -q
```

Optional extras: `.[cloud]` for S3, `.[proteomes]` for pinned gene calling,
and `.[reference]` for GenoPHI. Its MMseqs2 dependency and conda setup are in the
[demo runbook](docs/genophi-demo-reproduction.md#run-the-pinned-demo).
The Python import path remains `phage_cocktail` for compatibility.

Fetch and inspect the checksum-verified source table:

```bash
python scripts/fetch_genophi.py
python scripts/explore_interactions.py
```

This table contains 16,638 source-coded entries; the studies reconcile primary
MLC readings before retaining the 16,597 observed labels. Genome downloads and
receipts are documented in [data/README.md](data/README.md#rebuilding-the-locked-genome-inputs).

Validate the synthetic SDL exchange fixture:

```bash
phage-model-interface validate-bundle tests/fixtures/input_bundle/bundle.json --require-sequence host
```

See the [model-interface runbook](docs/model-interface.md) for a complete round trip.
The fixture verifies software, not prediction quality.

## Research boundary and next work

The target is the named source assay label. It is not automatically a measure of
adsorption, productive infection, resistance emergence or treatment success;
see [endpoint semantics](docs/endpoint-semantics.md).
Known-library retrieval does not establish unseen-phage performance, and pairwise
scores do not establish cocktail efficacy.

The next sequence-conditioned retrieval study needs a separately frozen library,
labels and evaluation protocol. Compare the evidence and learned routes at equal
candidate budgets, retain missing annotations as unknown, and reserve independent
data for release validation. The [roadmap](docs/host-conditioned-phage-retrieval-roadmap.md)
sets out the predictor/SDL boundary and the pending work.

## Where to look

- [Documentation index](docs/README.md): results, runbooks and proposed studies.
- [Data and protocols](data/README.md): resource locks, labels and rebuild commands.
- [Research weights](models/README.md): frozen development artifacts and limitations.
- `src/phage_cocktail/`, `scripts/`, `tests/`: implementation, commands and checks.
- `studies/annotation_recall_audit/`: separate evidence-coverage audit and outputs.
