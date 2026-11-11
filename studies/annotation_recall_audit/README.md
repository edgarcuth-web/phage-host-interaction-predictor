# Annotation-route recall audit

Quantifies how many experimentally positive phage-host pairs a mechanism/annotation
screening route loses, on the public Gaborieau *Escherichia* interaction data.

Results and limitations: [`docs/annotation-recall-audit.md`](../../docs/annotation-recall-audit.md).

This study is **separate from this repository's GenoPHI-based infrastructure**. It
does not import `src/phage_cocktail`, does not touch the frozen source-pair endpoint,
and establishes no clinical efficacy. It trains no released model; sections 4 and 6
of the report are in-sample upper bounds by construction.

## Data

Not redistributed. Clone the authors' repository and point the scripts at the
staged tables:

```bash
git clone --depth 1 https://github.com/mdmparis/coli_phage_interactions_2023
```

Then copy these files into one directory, renaming where noted:

| File in the authors' repository | Name expected here |
|---|---|
| `data/interactions/interaction_matrix.csv` | `interaction_matrix.csv` |
| `data/genomics/bacteria/picard_collection.csv` | `picard_collection.csv` |
| `data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv` | `defense_subtypes.csv` |
| `data/genomics/phages/guelin_collection.csv` | `guelin_collection.csv` |
| `dev/cocktails/data/picard+test_collection_phylogenetic_distances.tsv` | `phylo_dist.tsv` |
| `dev/predictions/per_phage_perf.csv` | `per_phage_perf.csv` |

All tables are semicolon-separated except `phylo_dist.tsv` (tab) and
`per_phage_perf.csv` (comma). The authors' repository carries no licence file;
check its terms before reuse.

The ablation figures the report compares against come from
`github.com/viti1igo/coli_phage_interactions_2023_replica` (also no licence file).

Section 9 uses the interaction matrices shipped in `github.com/Noonanav/GenoPHI`
under `data/interaction_matrices/` and `data/experimental_validation/`
(`step4_cross_dataset.py` takes a directory of those long-format CSVs).

Section 10 uses the PhageHostLearn archive, Zenodo 10.5281/zenodo.11061100:
`phage_host_interactions.csv`, `Locibase.json`, `all_loci_score_matrix.txt`,
`esm2_embeddings_loci.csv` and `RBPbase.csv`, all in one directory. The
genome archives in that record are not needed.

## Running

Python 3.13 with pandas, numpy, scipy, scikit-learn. No GPU. Each script takes a
data directory as its first argument.

```bash
python -I step1_ceiling.py   DATA_DIR        # scores, power, phylogenetic structuring
python -I step1b.py          DATA_DIR        # annotation/phylogeny redundancy; recall artefact check
python -I step2_rules.py     DATA_DIR 1      # single-feature rules, transfer, ABC decomposition
python -I step2c_recall_at_k.py DATA_DIR 1   # recall@k, precision@k, hit@k against the oracle
python -I step3a_rule_robust.py DATA_DIR     # all thresholds, paired test, chance level, bootstrap CI
python -I step3b_recall_robust.py DATA_DIR 1 5   # recall@k over 5 cross-validation partitions
python -I step4_cross_dataset.py GENOPHI_MATRIX_DIR   # which datasets can support the analysis
python -I phl_inspect.py PHL_DIR                 # PhageHostLearn shapes and identifier overlap
python -I step5_klebsiella_coverage.py PHL_DIR   # K-locus identity coverage ceiling
python -I step5b_rescue.py PHL_DIR               # centred ESM-2 versus K-locus identity
python -I step5c_phage_side.py PHL_DIR           # RBP similarity on the phage side
```

`step5b_rescue.py` supersedes the embedding comparison inside
`step5_klebsiella_coverage.py`: raw ESM-2 cosine over K-locus embeddings is
degenerate (median 0.999 between strains) and must be mean-centred first.

The second argument of the `step2*`/`step3b` scripts is the binarisation threshold
on the ordinal 0-4 plaque score. `step3a` sweeps all four itself.
`step3b`'s third argument is the number of partitions; it is the slow one
(roughly five minutes per partition on one core).

## Outputs

`outputs/` holds the raw stdout of every script, regenerated in one pass, plus
`00_provenance.md`: interpreter and package versions, and sha256 of every input
table and every script. Every figure quoted in the report traces to a file there.
File numbering follows report order, not run order.

## What supersedes what

`step3a_rule_robust.py` supersedes the transfer table in `step2_rules.py`: it
compares a transferred rule against the target phage's own rule on identical
labels, which the earlier version did not. `step3b_recall_robust.py` supersedes
the single-partition numbers in `step2c_recall_at_k.py`. The earlier scripts are
kept because they also produce the rule inventory, the ABC decomposition and the
precision@k / hit@k tables.
