# External evaluation intake: identity and provenance gate

Status: the row-level gate is **software-only**. The bounded public-source review
below is complete; it admits no independent validation cohort or new training
data. Model validation remains pending.

## Why

The existing GenoPHI E. coli subset is derived from Gaborieau et al.'s source observations. A renamed or resplit subset does **not** become a new independent release cohort. The already-inspected 177-host / 94-phage development panel, B2, and sensitivity seeds remain development evidence. See the [source reconciliation](host-baseline-results.md), [paired-source review](paired-split-freeze.md), and [roadmap](host-conditioned-phage-retrieval-roadmap.md).

The CLI compares two researcher-curated **row-level CSV ledgers**. It checks source-study lineage, original observation ID, exact host/phage assembly SHA-256, and exact sequence pairs. It accounts for unknown results without treating them as negatives. It does not fetch new datasets, compare near-identical sequences, infer assay semantics, or certify independence.

## Input contract

Each CSV contains exactly:

```csv
host_id,phage_id,host_assembly_sha256,phage_assembly_sha256,source_study_id,source_observation_id,endpoint_id,observation_status
```

- `host_id` and `phage_id` preserve source names. Names alone are not sequence identities.
- `host_assembly_sha256` / `phage_assembly_sha256` are verified **sequence-content** hashes, using the *same documented canonicalization* across ledgers. Do not mix raw FASTA byte hashes with normalized assembly hashes. If not verified, leave blank. This tool verifies SHA format, not the canonicalization procedure itself.
- `source_study_id` is the reviewed **underlying primary collection lineage**. Original and derivative subsets must share it. Distinct text IDs do not prove different sources.
- `source_observation_id` tracks an original measurement, where available. If unknown, leave blank rather than fabricating an ID.
- `endpoint_id` must state the source-specific assay/readout and version. `observation_status` must be `observed_positive`, `observed_negative` or `unknown`. Unknown is not a negative.

Synthetic input example (letters stand in for digests, not actual sequences):

```csv
host_id,phage_id,host_assembly_sha256,phage_assembly_sha256,source_study_id,source_observation_id,endpoint_id,observation_status
H1,P1,aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa,bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb,study-a,row-1,source-assay:mlc-positive:v1,observed_positive
```

After `pip install -e '.[dev]'`:

```bash
python scripts/audit_external_overlap.py \
  --development path/to/locked_development_ledger.csv \
  --candidate path/to/proposed_external_ledger.csv \
  --regime cold_host \
  --out artifacts/external_overlap_v1.json
```

Regimes: `cold_host` (a known phage library is allowed), `cold_phage` (known host genomes allowed), `double_cold` (neither allowed). All regimes block known common primary source/observation lineage, exact sequence-pair overlap, and missing assembly identity evidence; the requested cold axis must also lack exact assembly overlap.

Exit code **2** is a blocked claim; **0** means only `requires_manual_review`, never scientific approval. Invalid input fails nonzero; the report cannot overwrite an existing output path. The output records source readout IDs, missingness, exact overlap counts and blocking reasons without redistributing sequence data.

## Public-source review, 2026-11-13

The [frozen intake](../data/studies/public_new_host_intake.v1.json) records five
proposed panels, exact resource hashes, endpoint limits, reuse terms and roles.
This is a bounded review, not an exhaustive survey. The
[verification receipt](../data/studies/results/public_new_host_intake.v1/NUMERIC_QA.json)
rechecks ten original files and exercises the existing overlap gate. No model
was fit or scored on these proposals.

| Proposed panel | Verified evidence | Role and admission decision |
| --- | --- | --- |
| [Gaborieau full Picard export](https://github.com/mdmparis/coli_phage_interactions_2023/tree/99656d8c020fcd1dbbc30a02ee10274a07655562) | 402 hosts × 96 phages; 157 blank cells. All 177 development hosts join through the existing `source_bacterium` crosswalk; all 16,597 observed MLC and binary labels agree. | Same-study development extension only. Shared source lineage and reused observations block independent validation. New host assemblies and phenotype-table reuse terms remain unverified. |
| Gaborieau prospective 100 | 100 host metadata records; 300 selected-single measurements over 25 known phages, plus separate mixture/baseline readouts. Singles use source levels 0–6. | Quarantine the pinned files for a possible same-study follow-up. The singles table lacks the other 9,100 cells; they remain unknown. Selected measurements cannot establish full-library ranking performance. |
| [GenoPHI BASEL × ECOR](basel-ecor-sequence-audit.md) | Retain the prior 25-host/52-phage audit: 1,240 binary spot readouts and 60 unknowns. Host strain IDs are already represented in development. | Keep the inspected known-host/cold-phage reference. Different assemblies do not create new host strains; this is a different library from the fixed 94. |
| [PhageHostLearn data v2](https://doi.org/10.5281/zenodo.11061100) | The raw matrix has 215 hosts × 132 phages: 465 positive, 12,505 negative and 15,410 blank cells. The paper's model cohort instead has 200 hosts, 105 phages and 10,006 tested pairs. | Defer to a separately scoped Klebsiella task. Reconstruct cohort/feature selection before using the raw export. Visible spotting is not necessarily productive infection. Record data license CC BY 4.0 separately from MIT code. |
| [Berryhill enteric E. coli study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13410980/) | Publication screen: 54 isolates × 14 phages with no/turbid/clear plaque categories; raw reads at PRJNA1028583. | No verified match to the 94 library stocks, assembled-host join or admitted pair export. Defer to separate library/assay intake. Article CC BY 4.0 does not clear the supplement's distinct reuse notice. |

Keep source discrepancies visible. The Picard metadata has 403 hosts, with
`LF110` absent from the matrix. Its extra phages `T4LD` and `T7_Portugal` do not
extend the current 94. Prospective metadata/selection records use `BCH857-72`,
while readouts/compositions use `BCH857`: only 99 host names and 297 selected
pairs join literally. No alias is repaired here. `AN03` also occurs in Picard,
but with a different Gembase ID; names alone cannot establish strain identity.

The 100-host source files are held out of further model selection under the
intake's quarantine role. Their labels and published selection records were
inspected during this review, so they are neither a blind nor an independent
release set. Any follow-up needs assembly/grouping checks, resolved reuse terms,
a separate protocol and an estimand for selected, incomplete labels. Unknowns
and mixture outcomes cannot be substituted for unmeasured single-phage results.

**Reservation outcome:** no eligible independent new-host cohort for the exact
94-stock library was found among these panels. A release cohort and evaluation
protocol remain pending; research baselines can continue under their existing
limits. A different taxon or library requires a separate study.

Recheck the source bytes and numeric findings after installing the project:

```bash
python scripts/verify_public_data_intake.py \
  --fetch --data-root artifacts/public_intake_sources \
  --out artifacts/public_intake_review_new.json
```

Omit `--fetch` for an offline replay. The output path must be new. Exit 0 means
the intake findings reproduced, including the **blocked** same-source claim;
it does not admit a dataset. Source tables, reads and genomes stay outside Git.

## Public-panel follow-up, 2026-11-14

The [follow-up source receipt](../data/studies/public_panel_followup.v1.json)
records an inspection of Vitt's actual supplementary files and Kaneko's public
sequence metadata. These are candidates for a separate E. coli panel, not
admitted training data or a test of the current 94-output weights. The original
five-panel intake and its numeric receipt remain unchanged.

### Vitt: sequence metadata available, pair table not located

[Vitt et al. (2024)](https://doi.org/10.1016/j.isci.2024.108826) reports a
198-host/28-phage panel. The public Europe PMC
[supplement package](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10844046/supplementaryFiles)
was retrieved and its two supplementary documents were read:

- `mmc2.xlsx`, Table S1: 198 unique host IDs with 198 unique, nonmissing sample
  accessions. It links 78 samples to PRJEB22091, 80 to PRJEB14641 and 40 to
  PRJEB14086. These are sample identities, not verified assembled FASTAs.
- `mmc1.pdf`, Table S2, PDF page 3: 28 phage-level rows summarizing infected-host
  counts, phylogroups and sequence types. It is **not** a 198-by-28 pair table.
  The other tables concern experimental reagents; no full numerical pair matrix
  was found in these inspected supplements.

Figure 4 presents the host-range heatmap. The source uses a spot screen followed
by plaque assays for lysis-positive spots; its titre bins and no-plaques category
must retain that two-stage meaning. The nominal 5,544 pairs are not 5,544 verified
importable labels. No heatmap cells were digitized, and aggregate counts were
not converted into pair labels. A pair-level source table and sequence/stock
crosswalk are still needed.

### Kaneko: workbook access blocked and accession range inconsistent

[Kaneko et al. (2025)](https://doi.org/10.1089/phage.2024.0023) reports
144 hosts and 42 phages. The publisher lists seven XLSX supplements. Downloading
the first workbook returned HTTP 403 with a Cloudflare challenge, not workbook
bytes. None of the seven workbooks' contents has been inspected, so the measured
pair count, assay definitions and positive/negative/intermediate/unknown coding
remain unverified. The nominal 6,048 combinations must not be imported as labels.

The paper's experiment-accession ranges were checked against ENA read-run
metadata. The receipt contains exact request URLs and response hashes:

| Reported range | ENA evidence at review | Decision |
| --- | --- | --- |
| DRX534018–DRX534161, hosts | PRJDB18063 has 144 distinct experiments, samples and runs, exactly matching the range; all have FASTQ links | Raw-read metadata confirmed; phenotype joins and assembled genomes are not verified |
| DRX534162–DRX534203, phages | PRJDB18044 has 38 distinct experiments, samples and runs, DRX534162–DRX534199; all have FASTQ links | Does not resolve the sequence identities of all 42 paper phages |
| DRX534200 | PRJDB17934; Candidatus Magnetaquiglobus chichijimensis | Other-project record; do not import as a paper phage |
| DRX534201 and DRX534202 | PRJDB18047; Homo sapiens | Other-project records; do not import as paper phages |
| DRX534203 | PRJDB18042; Lentinula edodes | Other-project record; do not import as a paper phage |

This discrepancy is not a missing-read result: the last four reported accessions
resolve to other projects and organisms. It does not identify the correct
accessions for the four unresolved paper phages. Do not repair the range by
guessing adjacent IDs or silently reduce the study's 42 candidates to 38.
Sample titles and read availability do not establish stock identity, label joins
or new-host independence. No sequencing reads or genomes were downloaded.

**Next:** obtain the Kaneko workbooks and assay definitions, reconcile all 42
phage identities with original source metadata, and locate Vitt's pair table.
Endpoint/reuse review and exact/near-duplicate identity checks remain necessary
before freezing any separate-panel study or exact-current-library subset.
Neither source supplies an admitted independent release cohort at this stage.

## Mandatory manual steps even if this gate finds no overlap

1. Resolve study aliases, source-data derivations, shared experiments and source licenses; verify original accessions and dataset checksums.
2. Audit **near-duplicate genomes** with a prespecified grouped sequence-similarity or phylogenetic review. Exact SHA-256 disjointness alone is insufficient.
3. Review assay type, threshold, positive/negative/unknown semantics, and alignment of source endpoints. Source MLC/spotting positivity cannot silently become productive infection.
4. Audit training and model-selection exposure, feature preprocessing, reused test data and relevant pretrained-model contamination.
5. Freeze a truly independent cohort, allowed claims, subgroup definitions and evaluation protocol **before running predictions**.

A clean report never sets `independent_release_validation_established` to true. This tool blocks some invalid claims; it does not replace independent scientific review or experimental validation.

