# Host genome composition on the fixed phage library

This development experiment asks whether host genome composition improves
unseen-host recommendation over train-only phage prevalence and the existing
annotation model. It does not use Evo2, infer mechanisms or predict unseen phages.

## Frozen comparison

The [protocol](../data/studies/genome_host_retrieval.v1.json) inherits the exact
cohort, 94 named candidates, 41 unknown-pair exclusions, host groups, splits,
tie expectations and denominators of the
[table-only study](host-retrieval-benchmark.md). All 12 scenarios are refitted:
grouped hosts, random hosts and all ten published phylogroups. Broad phylogroups
are not equivalent to Omnii's clade holdout. The already inspected development
cohort is not independent validation.

Three fixed feature sets feed the same per-phage logistic estimator (C=0.1,
maximum 2,000 iterations, no parameter search):

| Representation | Host inputs | Phage input |
| --- | --- | --- |
| Existing annotation model | 137 defense counts and four surface categories | Known phage identity |
| Genome composition | 136 strand-canonical 4-mer frequencies | Known phage identity |
| Hybrid composition | Both host feature sets concatenated | Known phage identity |

Each host FASTA must match both its locked byte hash and the frozen manifest's
canonical assembly digest. Counts never cross contig boundaries. Any window
containing an ambiguous IUPAC base is skipped, rather than deleting the base and
inventing adjacency. Frequencies are normalized over valid windows. The complete
vocabulary is fixed without looking at data. Numeric means/scales and category
vocabularies use unique training hosts only. Test genomes are encoded separately;
they do not influence any other host's features or the fitted preprocessing.

The primary contrast is grouped-test all-complete-host Hit@2 for genome composition
minus train-only phage prevalence. Annotation/hybrid contrasts and all sensitivity
results are descriptive. Intervals resample complete host groups, conditional on
each fit and the fixed library. Adding hybrid features changes effective
regularization even with the same C. A result cannot isolate a causal biological
mechanism or justify selecting a representation for release.

## Reproduce

```bash
pip install -e . -c data/studies/host_retrieval_environment.v1.txt
python scripts/run_genome_host_retrieval.py --fetch --outdir artifacts/genome_host_retrieval_v1
```

The four tables and 177 host FASTAs total approximately 911 MB. Existing verified
files in `data/raw/genophi_host_study` are reused; omit `--fetch` for an offline
repeat. No phage genome, GPU, Evo2 installation or lab execution is required.
Use a new output directory on every run. Failed runs publish no partial study.

Outputs include the frozen protocol, feature and genome-QC tables, all splits,
36 models, complete candidate scores and tie-rank intervals, per-host metrics,
all six paired contrasts and input/code/output hashes. The historical table-only
command remains available. Published Evo2 species-level scores and GenoPHI's
25x25 demo remain different tasks; neither is an aligned comparator here.

## Results

The protocol and implementation were committed at
`88d031bafe759f0bf4a79faa2b5b7394358dc719` before computing the new metrics.
All 177 assemblies matched their frozen identities. All six methods scored the
same 94 candidates for every held-out host; no incomplete-label host was rescued
by dropping candidates.

Primary test: 34 complete hosts in 22 genomic groups, including 33 eligible hosts
and one all-negative host. The table includes all 34 hosts.

| Method | Hit@1 | Hit@2 | Hit@3 | Hit@5 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| Uniform random expectation / global prevalence | 0.2656 | 0.4371 | 0.5520 | 0.6892 | 0.4482 |
| Phage prevalence | 0.5588 | 0.7353 | 0.8529 | 0.8529 | 0.7021 |
| Annotations | 0.6471 | 0.8235 | 0.9118 | 0.9412 | 0.7770 |
| Host genome 4-mer | 0.7647 | 0.8824 | 0.8824 | 0.9118 | 0.8339 |
| Hybrid 4-mer + annotations | 0.7353 | 0.8824 | 0.9412 | 0.9412 | 0.8299 |

Genome-only and hybrid Hit@2 are both 30/34, versus 25/34 for prevalence and
28/34 for annotations. Eligible-only Hit@2 is 30/33 (90.91%) for both new models;
eligibility is 33/34, not a different library or an independently selected cohort.

The prespecified genome-minus-prevalence Hit@2 gain is **+14.71 percentage points**,
with a paired 95% group interval **[+4.00, +30.43]**. The hybrid-minus-prevalence
interval is [+3.45, +31.03]. These conditional development intervals exclude zero
on this split. Genome-minus-annotations remains inconclusive at +5.88 points
[-5.88, +19.35]; hybrid-minus-annotations is +5.88 [0.00, +16.13]. The identical
aggregate Hit@2 of the new models hides different host-level successes; their
paired difference interval is [-9.38, +8.82]. No representation is promoted.

On the 26-host validation panel, genome-only Hit@2 is 18/26 (69.23%), versus
22/26 (84.62%) for annotations, hybrid and prevalence. The validation and test
results therefore do not establish a consistently better representation.

## Every sensitivity

All figures use complete-panel hosts, including all-negative hosts. ECOR53 and
370D remain excluded only where held out because their library labels are
incomplete, exactly as in the parent study.

| Test holdout | Complete / eligible hosts | Prevalence Hit@2 | Annotation Hit@2 | Genome Hit@2 | Hybrid Hit@2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Random host | 26 / 25 | 0.7308 | 0.8846 | 0.8846 | 0.8846 |
| Phylogroup A | 37 / 36 | 0.6216 | 0.6216 | 0.5405 | 0.5676 |
| Phylogroup B1 | 22 / 21 | 0.6364 | 0.5909 | 0.6364 | 0.5909 |
| Phylogroup B2 | 61 / 60 | 0.9016 | 0.8033 | 0.7541 | 0.7541 |
| Phylogroup C | 7 / 7 | 1.0000 | 0.8571 | 0.7143 | 0.7143 |
| Clade I | 3 / 3 | 0.6667 | 0.6667 | 0.6667 | 0.6667 |
| Clade IV | 1 / 1 | 1.0000 | 0.0000 | 1.0000 | 1.0000 |
| Phylogroup D | 17 / 17 | 0.8824 | 0.8235 | 0.8824 | 0.8824 |
| Phylogroup E | 15 / 15 | 0.9333 | 0.9333 | 0.8667 | 0.9333 |
| Phylogroup F | 8 / 8 | 0.3750 | 0.6250 | 1.0000 | 0.7500 |
| Phylogroup G | 4 / 4 | 0.5000 | 0.7500 | 0.7500 | 0.7500 |

B2 holdout exposes a substantial failure: both genome-based models score 46/61
versus prevalence's 55/61. Genome-minus-prevalence is -14.75 points, descriptive
interval [-25.87, -3.77]; hybrid-minus-prevalence is -14.75 [-24.62, -5.76].
Genome composition can reflect lineage; these measurements do not identify why a
particular host succeeds or fails. The small F panel's 8/8 genome result and the
single Clade IV host cannot override the larger shift failures. Random-host
results still have five genomic components shared between train and test.

## Verification and next comparison

The committed [run record](../data/studies/results/genome_host_retrieval.v1/report.json)
and adjacent split, per-host, contrast and genome-QC tables retain every result.
Large full candidate-score exports, model files and extracted features are local
reproducible outputs rather than redistributed source genomes.

Numerical QA checks every output and input hash, reloads all 36 models and
reconstructs their scores exactly from the saved feature table. Direct enumeration
at the top-2 tie boundary checks all 204 primary host/method panels. The old four
methods' per-host metrics reproduce the parent study exactly. The QA record also
reports whether a second full run reproduced every report-listed output hash.
The suite passes 96 tests plus eight subtests; Python 3.11 and 3.12 CI pass.

The useful next question is whether aligned Evo2 host representations improve
these same frozen strain-level scenarios, especially B2 and A. Published
species-level scores cannot answer it. Any embedding extraction needs a new
protocol fixing checkpoint, layer, windows, contig/ambiguity handling and pooling;
any model or layer selection must stay outside test outcomes. An independent
release holdout and calibration remain separate requirements. No SDL prediction,
automatic retraining, engineering proposal or clinical inference is emitted.
