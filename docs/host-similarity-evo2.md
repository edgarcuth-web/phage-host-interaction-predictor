# Host similarity diagnostics and the Evo2 host comparison

The immediate question is whether a representation improves recommendation for
unseen bacterial hosts, including distant phylogroups. This is development on an
already inspected cohort, with no release promotion or clinical interpretation.

## Similarity control

The [similarity protocol](../data/studies/host_similarity.v1.json) preserves all
94 named candidates, unknown labels, complete-panel/eligible denominators,
host-macro metrics, exact score-tie expectations and all 12 split scenarios.
The implementation and protocol were committed before computing the new results.

`core_neighbour` transfers observed labels from five training hosts with the
smallest published core-genome distances. Neighbours are selected once per target,
without looking at phage labels. A tie at the fifth distance shares remaining
neighbour mass uniformly; there is no search over k. For each phage the score is
`(weighted positives + 1) / (weighted observed labels + 2)`. Unknown neighbour
labels contribute to neither sum. If all selected labels for a phage are unknown,
use its train-only Beta(1,1) prevalence. These fractional kernel scores are not an
expectation over a nonlinear missing-label estimator.

The source distance matrix was constructed jointly across the source collection.
It is label-free but transductive: this is a diagnostic control, not an online
distance implementation or an independent prospective benchmark. Raw source
distances are not an ANI estimate.

```bash
pip install -e '.[dev]' -c data/studies/host_retrieval_environment.v1.txt
python scripts/run_host_similarity.py --fetch --outdir artifacts/host_similarity_v1
```

Only five source tables are needed, about 2.5 MB. Every byte hash is verified.
All historical four methods are refitted without changing their parameters.
Outputs retain full candidate scores, models, per-host metrics, paired group
intervals and the following descriptive diagnostics:

- `host_distances.csv`: nearest train distance/IDs, selected neighbours, boundary ties and missing-label fallback counts, including incomplete-panel hosts.
- `phage_prevalence_shift.csv`: train/held observed counts and positive fractions for every phage. Held labels are inspected only after scoring.
- `host_diagnostics.csv`: complete-panel host metrics, distance and phylogroup, plus each model's Hit@2 difference from prevalence. Negative deltas expose failures without attributing a biological mechanism.

## Frozen Evo2 comparison

The [Evo2 protocol](../data/studies/evo2_host_representation.v1.json) fixes the
representation and comparison rules. **The actual deployment is not yet frozen:**
no local 40B weight hash, installed runtime, GPU access or real cache is available
in this session. A separate deployment lock must precede even the pilot.

| Setting | Fixed choice |
| --- | --- |
| Model | `evo2_40b`, existing local merged checkpoint; no download or model substitution |
| Layer | `blocks.28.mlp.l3`, expected width 8,192; fail if unavailable |
| Input | Canonical contig strand and lexicographic order, keeping duplicates and boundaries |
| Windows | 8,192 bases, non-overlapping, include short ends; select at most 32 systematically spaced windows per host |
| Ambiguity | Map non-ACGT IUPAC bases to N, preserving positions |
| Orientation | Encode each selected window and its reverse complement independently |
| Pooling | FP32 token mean, mean both strands, base-weighted host mean accumulated in FP64, store FP32 |
| Supervised model | Existing per-known-phage logistic specification, C=0.1, at most 2,000 iterations; preprocessing fits unique train hosts only |
| Compared inputs | Annotations; full-genome 4-mer; 4-mer + annotations; sampled Evo2; Evo2 + annotations; prevalence/random and core-neighbour controls |

No BOS/EOS, padding, species tags, context across windows, PCA, layer sweep or
fine-tuning. Contigs are not rotated. This is a sampled-window representation,
not whole-genome context encoding: selected bases and assembly coverage are
reported. Full-genome 4-mer is calculated from the same verified FASTAs during
extraction, not from the sampled windows. Fixed C does not equalize feature-set
capacity or effective regularization.

The upstream [embedding API](https://github.com/ArcInstitute/evo2/blob/53f195997257c56c00e5ef8d33a54f5baad143a6/README.md#embeddings)
demonstrates this module on 7B. Its use on 40B here is a prespecified engineering
choice, not a published optimum for host susceptibility. The 40B
[configuration](https://github.com/ArcInstitute/evo2/blob/53f195997257c56c00e5ef8d33a54f5baad143a6/evo2/configs/evo2-40b-1m.yml)
has 50 blocks and width 8,192. The pilot must confirm the module's actual output.
No alternative is selected after looking at results.

Run these commands in the existing **local Evo2 Python GPU environment**. This
script does not assume a remote generation endpoint supports embeddings. Use the
same repository revision on the GPU server and CPU evaluation machine.

```bash
# Pin checkpoint bytes, configuration, software and hardware before any pilot.
python scripts/extract_evo2_hosts.py \
  --checkpoint /absolute/path/to/evo2_40b.pt \
  --deployment-lock artifacts/evo2_deployment.json --freeze-deployment

# First three sorted hosts only: time, memory, shape, finiteness and coverage QC.
python scripts/extract_evo2_hosts.py \
  --checkpoint /absolute/path/to/evo2_40b.pt \
  --deployment-lock artifacts/evo2_deployment.json \
  --data-root data/raw/genophi_host_study --pilot --outdir artifacts/evo2_pilot

# After hardware QC, keep the same lock/settings for all 177 hosts.
python scripts/extract_evo2_hosts.py \
  --checkpoint /absolute/path/to/evo2_40b.pt \
  --deployment-lock artifacts/evo2_deployment.json \
  --data-root data/raw/genophi_host_study --outdir artifacts/evo2_hosts

# Transfer the entire cache, then evaluate with the frozen five source tables.
python scripts/run_evo2_host_comparison.py --fetch \
  --cache artifacts/evo2_hosts --outdir artifacts/evo2_host_comparison_v1
```

The existing 177 locked FASTAs must already be present on the extraction machine;
the prior genome study's `--fetch` command can obtain and verify them. Extracting
does not fetch weights or genomes. A deployment change needs a versioned amendment
before any scoring, rather than reusing the same cache identity. Pilot interaction
outcomes must not guide extraction settings.

The cache stores pooled embeddings, full-genome 4-mer frequencies, input assembly
and FASTA hashes, sampled-window coordinates/input hashes, code hashes, timing,
GPU peak memory and the deployment/protocol snapshots. The evaluator rejects
pilot/partial caches, mismatched hosts or assemblies, nonfinite/wrong-width
embeddings and changed payload hashes. It retains signed embedding dimensions
while preserving nonnegative validation for counts and frequencies. Cache hashes
establish integrity, not independent attestation that an external encoder was run.

No real 40B extraction or GPU pilot has been executed here. CPU fixtures exercise
the cache and comparison contract only. The full GPU extraction path is prepared
but awaits validation against the actual deployed runtime. The current extraction
publishes complete caches atomically; it does not resume failed full extractions.
This limits the first full run and is a reason to complete pilot QC first.

The primary development contrast is Evo2 minus prevalence on grouped-test Hit@2.
Every split, denominator and paired contrast remains visible, including B2 and A.
Unknown pretraining overlap remains a limitation. Independent locked release data
and calibration are separate requirements.
