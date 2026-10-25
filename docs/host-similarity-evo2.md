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

## Similarity results

The protocol and scorer were locally committed at `f1ea565` before computing
new results. The published implementation snapshot is `4843cfba55393efdf98c02e8a84a456dcdf27702`;
the reported reproducibility runs use that reachable revision. The
[committed run record](../data/studies/results/host_similarity.v1/report.json)
and adjacent diagnostics preserve all scenarios. No Evo2 result is included.

Grouped-host primary test: 34 complete-panel hosts, including 33 eligible hosts.
Hit@2 uses all 34. Values are exact expectations over score ties, so the neighbour
score is not an integer number of successful hosts.

| Method | Hit@1 | Hit@2 | Hit@3 | Hit@5 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| random | 0.2656 | 0.4371 | 0.5520 | 0.6892 | 0.4482 |
| phage_prevalence | 0.5588 | 0.7353 | 0.8529 | 0.8529 | 0.7021 |
| host_logistic | 0.6471 | 0.8235 | 0.9118 | 0.9412 | 0.7770 |
| core_neighbour | 0.7985 | 0.8678 | 0.9207 | 0.9345 | 0.8575 |

Neighbour-minus-prevalence Hit@2 is +13.25 percentage points, conditional 95%
group-bootstrap interval [+3.98, +26.49]. Eligible-only neighbour Hit@2 is
89.41% on 33 hosts. Validation neighbour Hit@2 is 81.39%, below prevalence and
annotations at 84.62%. The same cohort has already informed development; the
positive primary interval is not independent confirmation.

| Test holdout | Complete / eligible hosts | Prevalence Hit@2 | Annotation Hit@2 | Neighbour Hit@2 |
| --- | ---: | ---: | ---: | ---: |
| random_host | 26 / 25 | 0.7308 | 0.8846 | 0.8858 |
| phylogroup:A | 37 / 36 | 0.6216 | 0.6216 | 0.4865 |
| phylogroup:B1 | 22 / 21 | 0.6364 | 0.5909 | 0.5909 |
| phylogroup:B2 | 61 / 60 | 0.9016 | 0.8033 | 0.5574 |
| phylogroup:C | 7 / 7 | 1.0000 | 0.8571 | 0.5714 |
| phylogroup:Clade I | 3 / 3 | 0.6667 | 0.6667 | 0.9000 |
| phylogroup:Clade IV | 1 / 1 | 1.0000 | 0.0000 | 1.0000 |
| phylogroup:D | 17 / 17 | 0.8824 | 0.8235 | 0.6275 |
| phylogroup:E | 15 / 15 | 0.9333 | 0.9333 | 0.6238 |
| phylogroup:F | 8 / 8 | 0.3750 | 0.6250 | 0.6375 |
| phylogroup:G | 4 / 4 | 0.5000 | 0.7500 | 0.6250 |

B2 neighbour-minus-prevalence is -34.43 points, conditional interval
[-50.00, -20.00]. Its nearest-training-distance median is 0.018775 versus
0.001107 in grouped-test, about 17 times larger in source distance units.
This associates distant evaluation with transfer failure but does not establish
causality or attribute the earlier 4-mer failure to one mechanism. Broad
phylogroups remain distinct from a reconstructed clade holdout.

All historical four baselines reproduce their committed per-host metrics
exactly. Numerical QA independently recomputes 27,260 neighbour scores, reloads
all 12 annotation models exactly and enumerates top-2 ties in 170 primary
host/method panels. Two complete runs reproduce the report and every listed
output hash. The suite passes 107 tests and eight subtests. Large model and
candidate-score files remain reproducible run outputs; source genomes are not
redistributed.

```bash
python scripts/verify_host_similarity.py \
  --run artifacts/host_similarity_v1 \
  --repeat artifacts/host_similarity_v1_repeat \
  --out artifacts/host_similarity_numeric_qa.json
```

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
