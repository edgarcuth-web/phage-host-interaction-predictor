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

The protocol and implementation were frozen before computing this experiment's
metrics. Results will be recorded after the locked run and numerical verification.
