# GenoPHI demo reproduction: outcome and why the shipped reference is not a target

The GenoPHI protein-family workflow now runs end to end in this project on the
authors' 25-strain demo. The shipped reference output in the upstream repository
is **not** a valid comparison target for the pinned commit, so the original
milestone wording ("reproduce GenoPHI demo metrics locally") cannot be satisfied
as written.

This document records what ran, what came out, and the evidence for that claim.

## What ran

Pinned GenoPHI package version `1.0.3`, demo data at upstream commit
`f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f`, run with the settings in
`build_demo_command` (`--num_features 50`, `--num_runs_fs 5`,
`--num_runs_modeling 10`, `--method rfe`, `--filter_type strain`), which are the
settings documented in the upstream README's Demo section.

Observed cost on 2 cores: 1,893 s wall clock, 0.65 GB peak RAM, 14,853 protein
families, 2,217 features in the final table. Upstream README quotes roughly
25 minutes on a MacBook Pro M2.

The expensive stage is feature assignment, not clustering: clustering finished in
about 100 s, the remaining ~26 minutes is the search of all proteins against the
cluster representatives.

## The reproduction result

Our run produced two rows, at `cutoff_3` and `cutoff_4`. The upstream reference
output (`data/test_data/output/model_performance_metrics.csv`) has seven rows,
`cutoff_4` through `cutoff_10`.

On the single overlapping stratum, `cutoff_4`:

- ours: AUC 0.8138, accuracy 0.8741, precision 0.8276, recall 0.2637, F1 0.4000, MCC 0.4224
- reference: AUC 0.8944, accuracy 0.8891, precision 0.7963, recall 0.4624, F1 0.5850, MCC 0.5521

## Why the row sets differ

`Occurrence` in `generate_feature_tables` counts **how many feature-selection
runs selected a feature**, not how many genomes contain it. Its ceiling is
therefore `--num_runs_fs`. With the documented demo value of 5, every cutoff
above 5 is structurally empty, so rows at `cutoff_8`, `cutoff_9` and `cutoff_10`
cannot be produced by the demo settings at any data scale.

Which of the reachable cutoffs get modelled is set by the window
`min_features <= n < max_features`. The automatic values are
`min_features = 5 if interaction_count < 500 else 20` and
`max_features = interaction_count / 10 if interaction_count < 500 else interaction_count / 20`.
The demo has 625 interactions, so `min_features = 20`.

Two entry points resolve the upper bound differently on identical input:

- `protein-family-workflow` overrides it at the call site with `--num_features`, so 50 here
- `select-and-train` leaves it at the automatic value, 31.25 here

Their `cut_offs` lists also differ.

Measured occurrence curves on this demo, and the cutoffs selected under the
window `20 <= n < 50`:

| `--num_runs_fs` | occurrence curve (cutoff→features) | modelled |
| --- | --- | --- |
| 5 | 3→39, 4→21, 5→13, 6+→0 | 3, 4 |
| 10 | 3→64, 4→47, 5→37, 6→30, 7→24, 8→18, 9→13, 10→8, 11+→0 | 4, 5, 6, 7 |
| 25 | 3→88, 4→76, 5→68, 6→59, 7→53, 8→50, 9→44, 10→41, 12→33, 15→24, 17→20, 20→12, 22→6, 25→0 | 9 … 17 |

Occurrence scales with the number of runs, so the band of cutoffs whose feature
count falls inside the window slides upward as `--num_runs_fs` grows. None of the
three settings yields the reference's `{4, 5, 6, 7, 8, 9, 10}`; the documented
demo value (5) and the documented production value (25) both miss it.

## The reference output predates the pinned code and its parameters are unrecorded

The reference CSV was last committed in `3d2c97a` (2025-12-10). The pinned demo
data commit is `f7e30ab` (2026-04-22). The demo *inputs* are byte-identical
between those two commits, so input drift is excluded, but the code changed:
`b6df116` lowered the `min_seq_id` default from 0.6 to 0.4, `84909f8` fixed
metadata column leakage in feature tables, and `3343d9b` added the
`min_features` parameter.

Two of those turn out not to matter here, and are recorded so they are not
re-investigated:

- The clustering default change is immaterial on this data: `min_seq_id` 0.6
  gives 14,843 protein families and 2,219 final features, against 14,853 and
  2,217 at 0.4.
- The metadata leakage fix cannot affect this demo. The fix filters the
  phenotype matrix to `{sample_column, phage, phenotype_column}` before merging,
  and the demo's `ecoli_test_interaction_matrix.csv` has exactly those three
  columns (`strain`, `phage`, `interaction`). There is no metadata to leak.

Running the code at `3d2c97a` itself does not resolve the discrepancy, because
**that code cannot complete the documented demo**. At
`feature_selection.py:1606` it evaluates `min_features < num_features < max_features`
while the `--max_features` CLI default is the string `'none'`, with no type
coercion, so the first cutoff raises
`TypeError: '<' not supported between instances of 'int' and 'str'` and the run
dies after feature selection, before any modelling. Reaching the modelling stage
on that commit requires passing `--max_features` explicitly.

That same version also hard-sets `min_features = 5 if interaction_count < 500 else 20`
unconditionally, with no `if min_features is None` guard, so a supplied
`--min_features` is silently ignored and the window's lower edge is fixed at 20
for this demo. The guard appears only in the later code.

Under that fixed lower edge, the reference row set requires the occurrence curve
to cross 20 between cutoff 10 and cutoff 12 (12 is the next value in
`protein-family-workflow`'s list and is absent from the reference) while also
crossing `max_features` between cutoff 3 and cutoff 4 (`cutoff_3` is absent and
`cutoff_4` is present). Interpolating the three measured curves puts both
crossings in the right place at roughly `--num_runs_fs 15` with
`max_features = 50`. That is an untested inference, offered only to bound the
search: the generating parameters are not recorded anywhere in the upstream
repository, so the reference output is unusable as a reproduction target.

## A second, independent discrepancy

Recall in all three of our runs declines monotonically as the cutoff rises: 0.462
down to 0.264 at 5 runs, 0.253 down to 0.138 at 25 runs. The reference holds
recall between 0.440 and 0.482 across all seven of its rows while its AUC falls
from 0.894 to 0.682. AUC degrading that far without recall moving did not occur
in any of our runs. The cause is not established here.

## Consequences for how GenoPHI is used in this repository

- The demo verifies the pipeline. It does not reproduce the paper's benchmark,
  and it cannot reproduce the shipped reference numbers.
- The reference CSV should not be cited as GenoPHI's performance. It comes from
  an unrecorded configuration whose recall behaviour differs systematically from
  the documented demo settings. This matters wherever those numbers are used as
  a published baseline, including third-party comparisons against GenoPHI.
- When GenoPHI is used as a comparator on this project's frozen split, the
  cutoff must be fixed in advance. The workflow sorts its performance table by
  MCC descending and feeds the top row's cutoff downstream, so the reported
  stratum is selected on the same data it is evaluated on.
- `genophi --version` reports `1.0.2` while the distribution metadata reports
  `1.0.3`; verify the pinned version with `importlib.metadata.version("genophi")`
  instead.
- `scripts/run_genophi_reference.py` resolves the executable with
  `shutil.which("genophi")`, so it requires the environment to be activated
  rather than merely invoking its interpreter by path.
- Reproducibility is better than the absent `--seed` flag suggests: the
  implementation uses fixed incrementing split seeds (0, 1000, … ) and logs
  them, so same-machine reruns are expected to agree. Feature-selection
  stochasticity is still not exposed as a parameter.

## Remaining gap for a comparator on the 177 x 94 cohort

GenoPHI consumes per-genome amino-acid FASTAs. This project's locked resources
are nucleotide assemblies: 177 host `.fasta` files from figshare article
25941691 (CC BY 4.0, MD5-pinned) and the phage `.fna` set from the Gaborieau
revision. Neither GenoPHI nor the Gaborieau repository ships proteomes for
anything but the 25-strain demo, so gene calling is required, and its output
becomes derived data needing its own hashes and a pinned caller version.

`genophi train` exposes only `--set_filter {none,strain,phage}` and `--num_runs`
and performs its own splitting, so forcing it onto the frozen grouped-host split
means clustering once over all genomes, subsetting the feature table rows to the
training hosts, training, then predicting the held-out hosts.

This handoff has subsequently been completed in the
[frozen-split GenoPHI comparator](genophi-frozen-comparator.md): the 271
proteomes are receipted, supervised fitting uses outer-training hosts only,
cutoff 3 is fixed a priori, and primary/B2 predictions are evaluated through
the project's exact-tie retrieval code. The demo and shipped upstream reference
remain separate historical reproduction records.
