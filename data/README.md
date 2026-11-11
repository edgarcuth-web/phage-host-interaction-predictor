# Data provenance

Raw and processed datasets are intentionally not committed. The pinned source,
checksum, label semantics and expected summary for the first dataset are recorded
in the machine-readable registry at
[`registry/genophi_ecoli_subset.v1.json`](registry/genophi_ecoli_subset.v1.json).

## Step 1 dataset: GenoPHI E. coli subset

Source project: GenoPHI by Noonan et al.

Pinned upstream commit:

    f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f

Upstream file:

    data/interaction_matrices/ecoli_interaction_matrix_subset.csv

Download it reproducibly with:

    python scripts/fetch_genophi.py

The download script verifies SHA-256 before writing the file. Analysis scripts
verify the checksum and expected summary again before using the dataset.

The upstream schema is:

    strain,phage,interaction

The project loader normalizes it to:

    bacterium,phage,infection

This is a column-name normalization, not a claim that the source label directly
measures adsorption, productive infection, or resistance emergence. The
[later MLC audit](../docs/host-baseline-results.md) establishes exact binary
agreement on known primary readings and identifies 41 source zeros with missing
MLC. Mapping the target to mechanistic outcomes remains unverified.

Expected statistics for this pinned file:

- 16,638 source-coded rows
- 177 bacterial strains
- 94 phages
- 3,813 source-positive labels
- 12,825 source-negative labels
- 22.9% positive rate
- 100% of the 177 x 94 strain-phage matrix syntactically populated
- no duplicate strain-phage pairs

**Missing-pair rule:** a pair absent from a dataset is unknown/unobserved, not
negative. Do not fill missing pairs with zero. The source CSV is syntactically
complete, but 41 primary readings are missing. The new study admits 16,597
observed rows with 12,784 explicit negatives.

License status has not yet been independently verified. Do not redistribute the
downloaded data until the upstream data terms have been checked.

The later [paired-study review](studies/paired_data_reuse_review.v1.json) separately
records the primary Zenodo/Figshare CC BY 4.0 declarations and verifies archive
byte correspondence for 94 phage genomes and four analysis tables, plus selected
host-file metadata. This does not change the frozen GenoPHI registry's status or
treat a software license as a separate dataset license. Paired-study targets are
reconciled primary MLC observations.

Publication:

Noonan AJC et al. "Phylogeny-agnostic strain-level prediction of phage-host
interactions from genomes using machine learning." Nature Microbiology (2026).

The processed publication data are also archived on Zenodo:
DOI 10.5281/zenodo.21696332.

## Rebuilding the locked genome inputs

The 177 host genomes and 94 phage genomes are not committed either, and they are
large enough that losing the machine that held them is a realistic event. They
are fully reconstructible: every file is pinned by URL, SHA-256 and exact byte
size in a resource lock under [`studies/`](studies/), so a rebuild is a
verification, not a re-derivation.

Which lock covers what:

| Lock | Contents | Download hosts |
| --- | --- | --- |
| `studies/genophi_host_baseline.resources.v1.json` | 177 host genomes (910.4 MB), 94 phage genomes (6.8 MB), 5 tables (2.5 MB) | figshare, raw.githubusercontent |
| `studies/genome_host_retrieval.v1.json` | the same 177 host genomes, 4 tables | figshare, raw.githubusercontent |
| `studies/host_retrieval.v1.json`, `host_similarity.v1.json`, `evo2_host_representation.v1.json`, `evo2_retrieval_reference.v1.json`, `paired_source_review.resources.v1.json` | tables and upstream files only | raw.githubusercontent |

The 177 `host_genome` entries in the first two locks have been verified identical
in id, URL, SHA-256, size and path, so either lock rebuilds the same host cohort.

### Host genomes

    python scripts/fetch_pinned_host_fastas.py --out data/receipts/host_fastas.v1.json

Reads the `host_genome` entries of `studies/genome_host_retrieval.v1.json`,
writes under `--data-root` (default `data/raw/genophi_host_study`), and expects
910.4 MB across 177 files. It verifies SHA-256 and size per file, replaces
atomically, retries 429/500/502/503/504 with exponential backoff, and skips any
file already present and matching, so an interrupted rebuild can simply be rerun.
`--out` must not already exist; the receipt it writes is immutable.

### Phage genomes and tables

`verify_resources` handles every resource kind, so no separate phage fetcher is
needed. Filtering the lock fetches the non-figshare half on its own:

    python -c "
    from pathlib import Path
    from phage_cocktail.artifacts import read_json_object
    from phage_cocktail.host_study import verify_resources

    lock = read_json_object(Path('data/studies/genophi_host_baseline.resources.v1.json'))
    lock = dict(lock, resources=[r for r in lock['resources'] if r['kind'] != 'host_genome'])
    verify_resources(lock, Path('data/raw/genophi_host_study'), fetch=True)
    "

That retrieves 94 phage genomes and 5 tables, 9.1 MB on disk, in about 9 seconds,
each verified against its locked SHA-256. Dropping the filter line verifies the
whole lock including the host genomes.

### Fetching as part of a study

The study runners materialise whatever their lock requires when given `--fetch`,
for example `python scripts/run_paired_baseline.py --fetch --outdir <new dir>`.
Those runs also perform the study, and require a fresh `--outdir`, so they are
not a way to prepare inputs on their own.

### Network requirements

Host genomes come from `ndownloader.figshare.com`; phage genomes and tables come
from `raw.githubusercontent.com`. An environment that reaches GitHub but not
figshare can still rebuild the phage half and every table, which is enough for
work that does not need host sequence. Figshare declares CC BY 4.0 for the host
assemblies (article 25941691); see the
[paired-study review](studies/paired_data_reuse_review.v1.json) for the recorded
primary declarations.

### What fetching does not do

Fetching is opt-in everywhere. Without it, a missing resource is an error rather
than a silent download, and files already on disk are re-verified against the
lock on every run, so a modified or truncated input fails loudly instead of
being reused.
