# Evo2 retrieval: precomputed-score reproduction

The six prespecified evaluations match the rounded upstream reference metrics.
This establishes agreement between the released scores and upstream evaluation,
not an end-to-end Evo2 reproduction or validation of our deployed 40B model.

## Run

Install this repository's normal dependencies, then run from the repository root:

```bash
pip install -e .
python scripts/run_evo2_reference.py --fetch
```

No Evo2, PyTorch, CUDA or MMseqs2 installation is required for this step.
Approximately 49 MB are downloaded to `external/Evo2-Host-Prediction`. Only missing
files are fetched. Existing files with changed bytes fail verification rather than
being overwritten. Omit `--fetch` for an offline repeat after download.

The output directory must be new. To retain a second run:

```bash
python scripts/run_evo2_reference.py --outdir artifacts/evo2_reference_repeat
```

Each run saves five upstream metric JSON files, the fusion CSV, six logs and a
`report.json` with resource/code hashes, versions, actual commands, denominators
and deviations from the reference. Commands record the temporary paths actually
used before successful results are moved together to the output directory.
Metric mismatch saves the report but exits unsuccessfully; execution failure
does not publish partial results. No prediction is sent to SDL.

## Frozen reference and observed results

Source: Wang and Arsuaga,
[Evo2-Host-Prediction](https://github.com/Arsuaga-Vazquez-Lab/Evo2-Host-Prediction/tree/ec5b77384609dd8a11e3476b5212e893faaa6b27),
commit `ec5b77384609dd8a11e3476b5212e893faaa6b27`.
The [protocol and resource lock](../data/studies/evo2_retrieval_reference.v1.json)
bind the original evaluator, fusion module, README, license, mapping and matrices
by exact SHA-256 and byte count. The wrapper invokes the upstream CLI unchanged;
it does not implement its own rank weighting or tie handling.

The Gram-negative reference contains 4,490 evaluated phages and 308 candidate host
species. Scores are phage rows by host TaxID columns. Metrics use the upstream
host balancing and stable tie handling, including the upstream treatment of
phages with multiple recorded hosts. RRF uses only BLASTN, VirHostMatcher, PHIST
and Evo2, equal weights, `k=60`, and the intersection policy. No method subset or
parameter is selected using these results.

| Method | MRR | Hit@1 | Hit@5 | Hit@10 |
| --- | ---: | ---: | ---: | ---: |
| BLASTN | 0.303015 | 0.222067 | 0.380163 | 0.460248 |
| VirHostMatcher | 0.320724 | 0.231539 | 0.390200 | 0.492867 |
| PHIST | 0.267648 | 0.220878 | 0.308275 | 0.338145 |
| WIsH | 0.247588 | 0.169210 | 0.302036 | 0.424552 |
| Evo2 published 7B scores | 0.310563 | 0.193885 | 0.436781 | 0.554052 |
| Four-way RRF | 0.368147 | 0.268762 | 0.468856 | 0.584731 |

Every method evaluates 4,490 phages. The comparison tolerance is 0.000051:
upstream README metrics are rounded to four decimals, and the fusion CLI writes
six decimals. This tolerance checks numerical reproduction, not model quality.
The saved [aggregate run record](../data/studies/results/evo2_retrieval_reference.v1.json)
retains the observed precision and hashes. Large source matrices are fetched,
not redistributed in this repository. Upstream's Apache-2.0 declaration is
retained; it is not a separate audit of all underlying database reuse terms.

## Relation to the existing pairwise benchmark

This task ranks host **species** for a phage among a fixed candidate panel.
The existing GenoPHI-derived task predicts observed source `MLC > 0` for a
**strain–phage pair** using experimentally measured labels. These endpoints,
candidate sets and denominators are different and must remain separate.

Transposing a score matrix makes reverse queries computationally possible. It
does not validate strain-level phage recommendation. The recorded associations
are incomplete; absent links are unknown rather than measured negatives. A
reverse-retrieval study needs its own relevance definition, candidate panel,
query balancing and coverage assessment. Likewise, individual-pair scores or
union coverage do not demonstrate cocktail efficacy.

## Next: embedding generation and the 40B comparison

The present command does not require or access the deployed GPU server. Before
running it there, freeze a separate embedding experiment with these inputs:

| Item | Required record |
| --- | --- |
| Cohort | Exact genome accessions/versions, FASTA hashes and TaxID mapping; preserve the upstream development/test cohorts |
| Normalization | The disjoint 1,500-phage reference bank and its genome hashes; never fit normalization on the test cohort |
| Models | Exact 7B and 40B checkpoint revisions/digests, Evo2/Vortex versions, precision and GPU layout |
| Extraction | Layer names and selection rule, window size, stride, contig handling, ambiguous bases and pooling |
| Budget | Same admitted sequence inputs and a prespecified compute policy; report runtime, failures and support |
| Evaluation | The same upstream evaluator and fixed fusion components; retain the published-score reference separately |

The upstream example uses `blocks.24`, 8,192-token windows and stride 6,144.
Its DNA-only mode deletes non-ACGT characters, and its pooling avoids counting
overlapping tokens twice. These are implementation choices to document, not
silently replace during strict reproduction. The released repository provides
TaxID lists and scores; this delivery has not resolved or frozen the FASTA files
used to generate those scores.

`blocks.24` in a 40B model is not automatically the comparable representation
depth. Any layer search must be prespecified and limited to development data.
Replacing the model name alone is not yet an admitted 40B experiment, and 40B
performance cannot be assumed to exceed 7B.

After that comparison, test embeddings as representations on the repository's
existing frozen strain/phage groups and measured labels. Preserve its cold-host,
cold-phage and double-cold evaluation and all comparators. A new architecture or
LoRA fit is not needed to complete this first comparison. The independently
reserved release holdout remains a separate requirement.
