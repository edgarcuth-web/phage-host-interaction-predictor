# Research FASTA-to-ranking command

`phage-research-rank` reads assembled DNA, scores the frozen 94-phage library
with an explicitly selected baseline and writes a hash-bound JSON report.
It reuses the [matched-study artifacts](new-host-retrieval-results.md); no model
is trained or promoted by this command.

## Run

From a source checkout with an activated environment:

```bash
pip install -e .
phage-research-rank \
  --fasta host.fasta --method host-4mer --top-k 5 \
  --host-species "Escherichia coli" \
  --out artifacts/host_ranking.json
```

The output must be a new file. Input is ASCII IUPAC DNA FASTA, at most 64 MiB;
FASTQ, protein sequences and empty contigs are rejected. Completeness and
taxonomy are not inferred. `--host-species` is an optional caller declaration;
an explicitly supplied other species produces `not_estimated` for every candidate.

| Required method choice | Scorer | Genome conditioning |
| --- | --- | --- |
| `host-4mer` | Frozen primary 117-host, per-phage logistic weights | Actual within-contig, strand-canonical DNA 4-mer frequencies |
| `prevalence` | Beta(1,1)-smoothed counts from the same 10,957 training observations | None; the FASTA supplies identity and QC only |

GenoPHI's stored comparison scores do not supply a new-FASTA model adapter here.
The failed bilinear model is not selected as a default. Both supported methods
remain research baselines for the source MLC endpoint.

## Report and applicability

The report contains all 94 candidates, scores and exact-tie rank bounds, an
inspection top-K list, the library/assembly identities, raw FASTA and feature
digests, model/registry hashes, runtime and scorer-code provenance. The library
has 94 source IDs and 89 unique assemblies; duplicates are retained. There is
no annotation-based filtering or implicit change of the candidate pool.

Named top-K entries use lexicographic phage IDs within exact ties. `rank_min`,
`rank_max` and `expected_top_k_slot` preserve the ambiguity; the named list is
not a cocktail efficacy estimate. Scores are not calibrated infection probabilities.
Query uncertainty is explicitly `not_estimated`.

FASTA headers never attach a source host identity. Membership is checked only
against canonical assembled sequence hashes. Exact training and inspected
development-cohort matches are reported. Length and ambiguity checks describe
the training envelope; they neither certify assembly quality nor validate a query.
An unlabelled new genome has unverified taxonomic applicability.

A valid FASTA without usable 4-mer windows yields `not_estimated` for the genomic
method, with null scores/ranks for all candidates. The non-genomic prevalence
control remains available. Corrupted registries/weights, unsupported library
identities, incomplete score support and invalid arguments fail before output.

This local `research_candidate_ranking` v1 report is separate from
ModelInputBundle/PairPrediction and the proposed SDL retrieval receiver policy.
It does not admit source scores into productive-infection decisions.

## Pinned artifacts and verification

The [registry](../models/research_ranker.v1.json) binds the candidate snapshot,
source receipts, training counts, cohort identities and exact host-model bytes.
It loads JSON data only. For a relocated or non-editable installation, use
`--model-root` with a directory preserving the registry and its referenced
weight path. The prevalence method needs only the registry; host 4-mer also
needs `data/studies/results/new_host_retrieval.v1/primary_host_4mer.json`.

The [software conformance receipt](../data/studies/results/research_ranker.v1/NUMERIC_QA.json)
verifies 11,280 scores over all 60 primary held-out hosts and both methods.
One exact source AL505 FASTA also reproduces 188 scores end to end; its maximum
difference from saved study scores is about 2.2e-16. The installed console command
produces the same reports as the Python API. These checks establish software
conformance on inspected data, not independent biological validation.

- [AL505 host 4-mer report](../data/studies/results/research_ranker.v1/AL505.host_4mer.json)
- [AL505 prevalence report](../data/studies/results/research_ranker.v1/AL505.prevalence.json)

Run the cached-feature check without downloading genomes:

```bash
python scripts/verify_research_ranker.py --out artifacts/research_ranker_check.json
```

To repeat the full-FASTA check, obtain the single source file from its pinned
Figshare file ID, then pass it to the verifier. Source genomes remain outside Git.

```bash
mkdir -p artifacts/research_ranker_inputs
curl --fail --location https://ndownloader.figshare.com/files/46741495 \
  --output artifacts/research_ranker_inputs/AL505.fasta
python scripts/verify_research_ranker.py \
  --fasta artifacts/research_ranker_inputs/AL505.fasta \
  --out artifacts/research_ranker_full_fasta_check.json
```
