# GenoPHI on the frozen 177-host x 94-phage cohort

This study adds an aligned external protein-family comparator to the existing
host-genome retrieval experiment. The protocol and implementation are frozen
before comparator fitting; measured results will be recorded after execution.
It is a development comparison on the already inspected source-MLC cohort.

## Inputs and gene calling

The 276 resources in `genophi_host_baseline.resources.v1.json` have been rebuilt
and verified: 177 host assemblies, 94 phage assemblies and five source tables.
Receipts are committed; nucleotide assemblies, source tables, proteomes and
intermediate alignment databases are local generated inputs.

- [`host_fastas.v1.json`](../data/receipts/host_fastas.v1.json): 177 matching host
  files, 910,399,604 bytes. This final verification reuses the downloaded files.
- [`genophi_inputs.v1.json`](../data/receipts/genophi_inputs.v1.json): every
  source file's SHA-256 and byte size against the original resource lock.
- [`genophi_proteomes.v1.json`](../data/receipts/genophi_proteomes.v1.json): pinned
  caller/settings, input/output digests and counts for all 271 proteomes.

Pyrodigal **3.7.1**, genetic code **11**, is used for both sets. Each bacterial
assembly is trained independently in normal/single mode, using all its contigs;
genes are then called separately on each contig. Every phage uses meta/anonymous
mode, without short-genome self-training. The
[current Pyrodigal API documentation](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html)
distinguishes trained single mode from pretrained meta mode; Prodigal's
[input-type guidance](https://github.com/hyattpd/Prodigal/wiki/Advice-By-Input-Type#plasmids-phages-viruses-and-other-short-sequences)
recommends anonymous mode for isolated sequences below approximately 100 kb.

Importantly, meta mode uses **only built-in bins whose training information has
translation table 11**. The default bin collection also includes three code-4
models. Merely requesting code 11 when translating already predicted genes
would not constrain the preceding gene prediction.

Ends remain open, masking is disabled, minimum gene/edge lengths are 90/60 nt,
maximum overlap is 60 nt, and terminal stops are omitted from protein FASTAs.
Contig order is preserved. A single counter across all contigs in one assembly
generates `>370D_1`, `>370D_2`, etc.; the file is exactly `370D.faa`. Neither a
contig-name mapping nor a new genome-ID normalization is introduced.

The measured gene-calling run used eight workers and **52.90 s** wall time:

| Set | Genomes | Proteins | Protein FASTA bytes |
| --- | ---: | ---: | ---: |
| Bacteria, single mode | 177 | 872,342 | 278,579,663 |
| Phages, meta mode | 94 | 9,841 | 2,187,439 |

## Prespecified comparator

The machine-readable [protocol](../data/studies/genophi_comparator.v1.json)
fixes these decisions before fitting:

- GenoPHI distribution **1.0.3**, MMseqs2 **15-6f452+ds-2**, Pyrodigal **3.7.1**.
  Runtime receipts also hash all installed GenoPHI Python sources and the
  MMseqs executable. `genophi --version` still says 1.0.2; distribution metadata
  is the version check. Numerical dependencies are pinned in
  [`genophi_comparator_environment.v1.txt`](../data/studies/genophi_comparator_environment.v1.txt).
- One host protein-family vocabulary over all 177 host genomes and one phage
  vocabulary over all 94 phages, shared across scenarios/repeats. Identity 0.4,
  coverage 0.8, sensitivity 7.5; upstream singleton removal and co-occurrence
  collapse, with no additional sample-cluster presence filter.
- `curate_mlc` admits **16,597** observed labels; **41** missing MLC readings
  remain unknown. The phenotype export contains only `strain,phage,interaction`.
- Recompute `group_hosts` at threshold `1e-4`, then `split_by_host_group` at seed
  42 and fractions `[0.7,0.15,0.15]`. Require exact equality with the frozen host
  manifest. The primary split has **117/26/34** train/validation/test hosts.
  B2 uses the existing phylogroup holdout with touching host-group purge.
- Subset full feature-table rows to outer-training hosts **before** calling
  any supervised GenoPHI function. Its internal feature selection, MCC grid
  search and early stopping operate entirely within this subset. Candidate
  prediction receives IDs and protein features, with no held-out labels.
- RFE selects **200** features in each of **5** feature-selection runs. Fix
  occurrence **cutoff 3 a priori**, with the explicit window **20 <= n < 501**.
  Occurrence counts selection runs, not genomes. At most `5*200/3` features can
  pass, so the upper bound is explicit and nonbinding. An unsupported fixed
  stratum is a blocked run, not permission to change the cutoff.
- Retain **10** internally selected models per complete workflow and repeat the
  entire feature-selection/modeling stage **three** times. Internal filtering
  uses strain protein-content hierarchical clustering with 20 clusters, as in
  the upstream machinery. Dynamic weighting and SHAP export are disabled.
  Grid: iterations 500/1000, learning rate 0.05/0.1, depth 4/6, Logloss, and
  early stopping after 100 rounds. Upstream models are not refitted on all
  outer-training rows after their internal selection.
- The primary rank score is the **median positive-class confidence over all 30
  models**. Report complete-repeat medians and individual-model spread as well;
  never select a run using an outer result. Fixed upstream seeds mean these
  replays can agree exactly. Repeatability is not independent sampling uncertainty.

The protein vocabulary is transductive in genome sequences, as requested; the
supervised boundary is outer-train-only. This distinction is explicit because
the earlier fixed 4-mer representation is inductive. This bounded configuration
is an aligned comparator, not a reconstruction of the authors' paper benchmark
or the unrecorded upstream reference CSV.

### Memory-bounded assignment

The installed upstream implementation searches its protein database against
itself, then reads and sorts the entire all-versus-all alignment table to choose
one best hit per query. The adapter keeps the same search and family mappings,
but reduces alignment chunks before merging: float32 sequence identity
descending, float32 e-value ascending, first input row on exact ties. No search
threshold, family definition, or label-dependent rule changes.

Tests compare against the installed upstream reducer on randomized ties, extreme
e-values and chunk boundaries. Each real clustering run also compares the first
50,000 alignments byte-for-byte against the upstream reducer before feature use.
This limits the reduction's memory cost without silently changing the method.

## Reproduce and resume

On Ubuntu 24.04 / WSL2, install MMseqs2 (the measured apt package is
`mmseqs2=15-6f452+ds-2`) and create/activate a Python 3.12 environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev,reference]' -c data/studies/genophi_comparator_environment.v1.txt
python -c "from importlib.metadata import version; print(version('genophi'))"
mmseqs version
```

Rebuild the nucleotide inputs with the commands in
[`data/README.md`](../data/README.md#rebuilding-the-locked-genome-inputs), then:

```bash
python scripts/call_genophi_genes.py \
  --out data/receipts/genophi_proteomes.local.json \
  --source-receipt data/receipts/genophi_inputs.local.json --workers 8

python scripts/run_genophi_comparator.py \
  --proteome-receipt data/receipts/genophi_proteomes.local.json \
  --outdir artifacts/genophi_comparator_v1
```

Use a **new receipt filename**: committed receipts are immutable historical
records, not files to overwrite. Per-genome checkpoints under
`data/derived/genophi_proteomes/checkpoints/` allow verified gene-calling reuse.
The same comparator command resumes an interrupted run if tool versions,
protocol, code and completed artifact hashes still match. Named stages
`--stage prepare`, `--stage cluster`, and `--stage fit` stop at boundaries;
`--scenario grouped_host --repeat 0` restricts fitting to that frozen unit.

Each feature-selection/modeling run has its own UTF-8 log and hash receipt.
Failures record the exception and elapsed time, and do not publish a completed
stage or comparison. Altered inputs, missing models, partial candidate panels,
and a missing fixed feature stratum are rejected. `data/derived/` and
`artifacts/` are ignored by Git.

Byte-hashed tracked files must retain LF line endings. `.gitattributes` enforces
this for fresh Windows/WSL checkouts; no frozen protocol or source hash is
redefined to accommodate CRLF conversion.

## Metric contract and existing results

All 94 candidates are scored for each held-out host. Metrics use
`phage_cocktail.retrieval.evaluate_retrieval`: exact tie expectations,
Hit@1/2/3/5, eligible-only Recall@K and MRR, with both all-complete-host and
positive-eligible-host denominators. A host with any unknown panel label is
excluded whole, and an all-negative complete host contributes Hit/MRR zero.
Group bootstrap contrasts use the existing paired retrieval implementation.

The frozen [host-genome result](genome-host-retrieval.md) is the reference:

| Test panel | Complete / eligible hosts | Phage prevalence Hit@2 | Host 4-mer Hit@2 |
| --- | ---: | ---: | ---: |
| Grouped-host primary | 34 / 33 | 25/34 = 73.5294% | 30/34 = 88.2353% |
| B2 holdout | 61 / 60 | 55/61 = 90.1639% | 46/61 = 75.4098% |

The earlier [host-baseline report](host-baseline-results.md) describes source
reconciliation and the classification experiment; its AP is not Hit@K.
GenoPHI's internal AUC/MCC outputs are training-fold diagnostics, not the
aligned retrieval comparison. Scores still concern source-coded MLC positivity,
not productive infection, treatment response or cocktail efficacy.
