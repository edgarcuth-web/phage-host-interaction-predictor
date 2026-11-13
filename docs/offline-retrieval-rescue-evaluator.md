# Offline retrieval-rescue evaluator v1 (software-only)

Status: **implemented offline scoring-audit utility, not a trained retriever,
not an approved/frozen scientific benchmark and not a validated mechanism
screen**. No provider API, Evo2 inference, GPU, experiment, or clinical claim
is involved. The proposed scientific comparison is in
[retrieval-rescue benchmark v0](retrieval-rescue-benchmark-v0.md). The utility
implements score-file auditing; the scientific comparison remains proposed.

> **Related, different estimand.** `docs/annotation-recall-audit.md` measures the
> evidence route's own recall loss on different data, with a global rather than
> per-host budget and without tie-expected metrics. Its numbers are not values of
> this evaluator's `primary_k_rescue`. See `docs/estimand-reconciliation.md`.

## Question it can answer

Given **the same exact full host x phage panel** and existing source-endpoint
scores from two independent retrieval routes, count source-positive candidates
the learned route finds in top-K but the evidence route does not. The tool also
reports whether an **equal-size union** of the two routes rescues positives.

- Route `evidence`: source-mechanism annotations or other evidence-derived
  scores. A phage can be unannotated; it must still receive a finite score
  in the comparable full panel. Today's GenoPHI host feature baseline is only
  an **annotation proxy**, not verified host-receptor compatibility.
- Route `learned`: scores from an independently evaluated learned genome
  compatibility model, or synthetic scores *for software tests*. The tool does
  **not** compute learned representations or train the model.
- Other score methods (e.g. `prevalence`, `genome_4mer`, `fusion`) may be
  supplied for comparable full-panel ranking and paired comparisons. The
  `union` identifier is reserved.

## Inputs and source identity

The required local JSON manifest supplies all IDs and immutable byte hashes.
File references must be relative to the manifest directory, inside that
directory, with no `..` traversal, and every referenced file must already
exist and pass its exact SHA-256.

```json
{
  "schema_version": 1,
  "research_status": "development_only",
  "study_id": "reviewed-local-run-name",
  "endpoint": "source-specific-endpoint-id:v1",
  "hosts": ["H1", "H2"],
  "phages": ["P1", "P2", "P3"],
  "host_groups": {"H1": "G1", "H2": "G2"},
  "observations": {
    "path": "inputs/observed.csv",
    "sha256": "64-lowercase-hex-digest",
    "endpoint": "source-specific-endpoint-id:v1"
  },
  "score_sources": {
    "evidence": {
      "path": "inputs/evidence.csv",
      "sha256": "64-lowercase-hex-digest",
      "endpoint": "source-specific-endpoint-id:v1",
      "producer": "registered-evidence-model:v1",
      "training_hosts": ["T1"]
    },
    "learned": {
      "path": "inputs/learned.csv",
      "sha256": "64-lowercase-hex-digest",
      "endpoint": "source-specific-endpoint-id:v1",
      "producer": "registered-sequence-model:v1",
      "training_hosts": ["T1"]
    }
  },
  "ks": [1, 2, 3],
  "union_schedule": ["evidence", "learned"],
  "bootstrap": {"resamples": 2000, "seed": 20261024}
}
```

The values above are placeholders and **not hashes or model evidence**.
Use real SHA-256 digests computed from the exact on-disk CSV bytes.

`observed.csv` columns: `bacterium,phage,infection`; each row is an
*actually observed* source-coded binary 0/1 pair. Missing pairs are **omitted**,
not filled with 0. A host lacking any one of its phage labels is excluded
entirely from primary metrics and rescue counts, while remaining present in
the ranking/missingness audit. No outcome-derived subsetting is allowed for
score tables.

Every score file has exactly `bacterium,phage,score` (additional columns
ignored by the v1 evaluator), with one finite score per pair of the full
host x phage Cartesian product. Higher scores rank first; scores are **not
probabilities**. All methods must declare the same exact endpoint. The code
checks the manifest-declared training-host IDs against evaluation hosts, but
cannot independently establish that training, feature fitting, normalization,
model selection or underlying pretraining were leakage-free. Evidence/model
card review remains necessary.

The `union_schedule` is a nonempty repeating list containing both route
names. At each turn it takes the highest-ranked *not yet selected* phage
from the named route. This creates a single deterministic full-library order
in which top K has **exactly K unique phages**. For example
`["evidence", "learned"]` alternates routes, skipping duplicates. It is a
fixed **policy hypothesis**, not a measured optimal fusion algorithm. Do not
tune its schedule on the held-out test labels.

## Metrics versus deterministic candidate membership

Per-route Hit@K, eligible-host Recall@K and MRR reuse
`phage_cocktail.retrieval.evaluate_retrieval`, which **analytically averages
all orderings inside equal-score tie blocks**. A union output is a
deterministically assembled ordered set and has no score ties.

The exact phages selected and rescued require concrete membership, so
the audit uses a *separate* fixed tie break: descending score, then ascending
phage ID. The rescued-positive count is the number of source-positive
phage/host pairs appearing in the learned (or union) top-K but **not** the
evidence top-K. The primary reported rescue uses K=2 when available, else
the smallest declared K. The fraction is over **all observed positive pairs
among complete-panel hosts**, not over missing labels, and a zero denominator
is returned as null. The audit intentionally does not call these counts
tie-expected rescue estimates.

Each method is also compared against `evidence` with paired host-group
bootstrap for primary Hit@K; these are **conditional development intervals**.
The frozen host panel and evaluation labels must not be used to fit/tune the
incoming model scores.

## Run locally (no network)

```bash
pip install -e ".[dev]"
pytest -q tests/test_retrieval_rescue.py
python scripts/run_retrieval_rescue.py \
  --manifest /absolute/path/to/manifest.json \
  --outdir /absolute/path/to/new-audit-dir
```

The `--outdir` path must not exist. The process only reads the local manifest
and its referenced CSVs. It writes:

- `report.json`: summary metrics, grouped intervals, source SHA-256,
  manifest SHA-256, code hashes, support and limitations;
- `ranked_candidates.csv`: complete deterministic rankings for every route
  and the union, including missing label status and selected union route;
- `rescue_by_host.json`: each complete host/K selected sets and recovered
  positives; incomplete hosts carry null rescue, not zero;
- `per_host_metrics.csv`: existing evaluator per-host output, one method per row.

All three data outputs have SHA-256 digests inside `report.json`, and the
publication uses an atomic staged directory. There is no output if a
validated preflight fails. Training and model deployment are entirely
outside this script.

## Scientific gates still outstanding

1. Freeze and review the actual 177-host / 94-phage score inputs, source
   semantics, phage identity/duplicate policy and data-reuse permissions.
2. Establish a genuinely independent **mechanism-based** comparator;
   a currently available annotation model is not automatically that.
3. Fit a learned host+phage representation model **without using evaluated
   host/labels**; there is no such model file included in this PR.
4. Run the new scorer against the same frozen library, inspect all subgroup
   failures, and reserve an untouched release holdout.
5. Seek biological evidence before claiming a newly identified infection
   mechanism or any therapeutic/cocktail effect.
