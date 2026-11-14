# Kaneko workbook audit

Status, 2026-11-14: the seven supplied publisher workbooks contain a complete
source-grade matrix. They support a **separate E. coli panel candidate**. No
training data, independent release cohort or current-94-stock validation is
admitted. The [receipt](../data/studies/kaneko_workbook_intake.v1.json) pins the
files and metadata; [numeric QA](../data/studies/results/kaneko_workbook_intake.v1/NUMERIC_QA.json)
records the reproduced results. The prior
[follow-up receipt](../data/studies/public_panel_followup.v1.json) remains a record
of what was known before these workbooks were supplied.

Source: [Kaneko et al., PHAGE, 2025](https://doi.org/10.1089/phage.2024.0023).

## What is actually in the files

The `sj-xlsx-N` publisher suffix identifies the workbook, not the table number.

| N | Sheet | Verified content |
| --- | --- | --- |
| 1 | SF1 | Illustrated score definitions; machine-readable assay caption in A6 |
| 2 | SF2 | 45 phages × 144 hosts, numeric source grades in B3:EO47; assay caption in A49 |
| 3 | ST1 | 144 numbered host records, genome lengths, sequence types and AMR metadata |
| 4 | ST2 | Inventory of 42 study phages and four references, including T6 |
| 5 | ST3 | 46 × 46 phage ANI matrix; not host-range observations |
| 6 | ST4 | 46 × 46 coverage matrix associated with ANI; not host-range observations |
| 7 | ST5 | Three aggregate host counts for each of the 45 measured phages |

SF2 has **6,480 observed source cells**, comprising **6,048 study-phage cells**
and **432 T2/T4/T7 control cells**. There are no missing cells within that
rectangle. Host names remain strings `001`–`144`. **T6 is absent from SF2**;
its presence in the inventory/ANI files does not supply 144 negative labels.

| Source grade | Meaning from SF1/SF2 captions | Study-phage cells | All cells |
| --- | --- | ---: | ---: |
| 0 | No visible change from surrounding bacterial lawn | 4,160 | 4,530 |
| 1 | Discoloration without clearance | 547 | 565 |
| 2 | Turbid spot with clearance visible from below | 367 | 375 |
| 3 | Complete clearing | 974 | 1,010 |

These are high-titre, overnight **spot-clearing grades**. A zero is an observed
no-clearing result under this assay; grades 1 and 2 are observed intermediate
results. They are not missing values. Neither a binary threshold nor a productive
infection endpoint is inferred. Replicate-level measurements are not supplied.
Full article methods remain behind the publisher's access wall; the available
workbook captions establish this endpoint but do not resolve every analysis rule.

## Aggregate and publication discrepancies

ST5 B2:D2 literally prints `>1`, `>2` and `>3`. All **135 aggregate counts**
instead reproduce SF2 at **≥1, ≥2 and ≥3**, respectively. The extractor checks
this agreement and retains the original headings as a recorded discrepancy.
It does not alter the source or select a training threshold.

For ΦWec277, SF2 and ST5 agree on 135/144 at ≥1, 129/144 at ≥2, and 105/144
at grade 3. Those fractions are 93.75%, 89.58% and 72.92%. None reproduces the
abstract's 78.5%. An uninspected full-text analysis definition/subset or author
clarification is needed; this audit neither corrects the abstract nor changes
the matrix to match it.

## Sequence routes, with identity limits

Public ENA metadata provide one-to-one **number/name joins**:

- All 144 ST1/SF2 host numbers match `JMU_No_N_read` sample titles in PRJDB18063,
  with distinct samples, experiments and runs at DRX534018–DRX534161.
- Thirty-eight study-phage names match `phiWecN_sequence_read` in PRJDB18044,
  with distinct samples, experiments and runs at DRX534162–DRX534199.
- Four remaining study-phage names match pre-existing complete-genome records:

| Workbook phage | ENA sequence accession |
| --- | --- |
| ΦWec189 | [LC739538](https://www.ebi.ac.uk/ena/browser/view/LC739538) |
| ΦWec191 | [LC739540](https://www.ebi.ac.uk/ena/browser/view/LC739540) |
| ΦWec193 | [LC739541](https://www.ebi.ac.uk/ena/browser/view/LC739541) |
| ΦWec196 | [LC739542](https://www.ebi.ac.uk/ena/browser/view/LC739542) |

Thus all 42 names have a sequence route, rather than four unaccounted-for
study-phage names. The paper's DRX534200–DRX534203 range error still stands:
those records belong to unrelated projects. The legacy records are a separate
route, not guessed replacement DRX IDs.

The receipt supplies exact ENA requests and response hashes. Metadata row ordering
can change, so replay requires an order-independent digest over **all returned
fields**, reports the actual byte hash, and fails if field contents change.

These joins do not establish exact experimental-stock identity, assembled
sequence content/quality or independence from prior development strains. No reads
or genome sequences were downloaded, and no assembly hashes or accession versions
were verified. T2/T4/T7 names must not be silently joined to other stocks, including
stocks in the current 94-candidate library.

## Reproduce locally

Place the seven original workbooks in one directory. Local download suffixes such
as `(3)` are accepted because files are selected by SHA-256, not by their names.
The command uses only Python's standard library and does not modify workbooks.

```bash
python scripts/audit_kaneko_workbooks.py \
  --xlsx-dir /path/to/seven_workbooks \
  --fetch-metadata --metadata-dir artifacts/kaneko_metadata_new \
  --out-dir artifacts/kaneko_audit_new
```

For offline replay, omit `--fetch-metadata` and supply the existing metadata
directory. Output directories must be new. The command emits:

- `source_spot_grades.csv`: the 6,480 original ordinal grades, with study/control
  roles, host labels, endpoint, source sheet/cell and workbook hash.
- `host_read_crosswalk.csv` and `phage_sequence_routes.csv`: number/name joins
  to metadata, preserving sequence-route differences and unresolved stock identity.
- `AUDIT.json`: counts, source checks, discrepancies and export hashes.

The publisher gives a non-OA supplement reuse notice; an open phenotype-data
license has not been verified. Source workbooks, bulk grades and raw metadata
remain outside Git. The exports are source-audit intermediates, not inputs to
the existing training cohort or overlap gate's binary endpoint ledger.

## Next decision

The [sequence follow-up](kaneko-sequence-intake.md) has now retrieved and checked
143 third-party host assemblies linked to the original runs and four versioned
phage genomes. Host 127 and 38 phages still have only read routes in the inspected
sources. The workbook receipt remains unchanged.

Resolve reuse terms and the endpoint/publication discrepancy, then obtain or
assemble quality-checked genomes and verify the phenotype/sequence/stock joins.
Audit exact and near-duplicate overlap against development sources. A separate
panel protocol must freeze endpoint, controls and held-out host groups before
model fitting or scoring. The inspected source is not a blind release set and
cannot validate the current 94-output weights by itself.
