# Phage–Host Interaction Predictor

Rank existing phages for experimental screening from an assembled bacterial
genome. The current research command supports E. coli and a frozen 94-phage
library, using either host 4-mer features or training-only phage prevalence.
No model has established the evidence needed for promotion.

## Current results

Hit@2 measures finding at least one source-positive candidate in the first two
tests. These are inspected development results; each panel retains its own
endpoint, candidate library and frozen split.

| Test panel | Candidates | Prevalence Hit@2 | Host 4-mer Hit@2 | GenoPHI Hit@2 |
| --- | ---: | ---: | ---: | ---: |
| Primary seed 42, 34 hosts | 94 | 73.53% | 88.24% | 91.18% |
| B2 holdout, 61 hosts | 94 | 90.16% | 75.41% | 78.69% |
| Kaneko, 7 hosts in 5 genomic groups | 42 | 42.86% | 42.86% | Not evaluated |

- **94-phage studies:** 177 hosts and 16,597 observed source-MLC labels; 41
  missing primary readings remain unknown. The [matched comparison](docs/new-host-retrieval-results.md)
  retains the bilinear model's failure against prevalence on both test panels.
  [GenoPHI seed sensitivity](docs/genophi-seed-sensitivity.md) did not establish
  the registered robustness criterion.
  The [source-evidence audit](docs/missed-candidate-audit.md) separates ranking
  misses, source-all-negative hosts and missing readings, with named candidates
  and source pointers.
  The [training-only complementary panel comparison](docs/complementary-panel-results.md)
  rescues four primary-test hosts versus prevalence at two tests, but newly
  misses four on B2. A stable benefit is not established.
- **Kaneko:** a separate 143-host, 42-candidate study targeting grade-3 complete
  spot clearing. Both methods tie on test Hit@2; capped tests consumed are 1.571
  for prevalence and 1.714 for host 4-mer. [Neither improvement criterion passed](docs/kaneko-baseline-results.md).
  The seven-host test has limited precision and does not validate the 94-output models.

Scores are uncalibrated source-endpoint scores. They do not establish productive
infection, cocktail efficacy or treatment benefit. Missing observations stay
unknown; see [endpoint semantics](docs/endpoint-semantics.md).

## Quick start

Use Python 3.11+ from the repository root:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -e ".[dev]"
python -m pytest -q
```

Rank an assembled host FASTA with an explicitly selected method:

```bash
phage-research-rank --fasta host.fasta --method host-4mer --top-k 5 \
  --host-species "Escherichia coli" --out artifacts/host_ranking.json
```

Use `--method prevalence` for the non-genomic control. The command reuses frozen
weights, keeps all 94 candidates, and records score ties, provenance,
applicability and unavailable estimates. The [runbook](docs/research-fasta-ranking.md)
covers inputs and verification. Its local report is separate from SDL admission.

Optional extras: `.[cloud]` for S3, `.[proteomes]` for pinned gene calling,
and `.[reference]` for GenoPHI; the [GenoPHI runbook](docs/genophi-demo-reproduction.md#run-the-pinned-demo)
includes MMseqs2 setup. CI tests install `.[dev,cloud,proteomes]`.
The Python import path remains `phage_cocktail` for compatibility.

## Next work

Inspect the saved Kaneko failures before choosing one data or representation
hypothesis for a separately frozen experiment. Continue auditing suitable
validation inputs; independent release validation remains pending.
See the [roadmap](docs/host-conditioned-phage-retrieval-roadmap.md#next-work).
Evo2 full-cohort extraction remains [paused](docs/host-conditioned-phage-retrieval-roadmap.md#paused-evo2-branch).

## Where to look

| Need | Entry point |
| --- | --- |
| Study results, protocols and reproduction | [Documentation index](docs/README.md) |
| Source data, resource locks and downloads | [Data guide](data/README.md) |
| Frozen weights and loading instructions | [Model cards](models/README.md) |
| SDL bundle validation and unresolved predictions | [Model interface](docs/model-interface.md) |
| Implementation and checks | `src/phage_cocktail/`, `scripts/`, `tests/` |
