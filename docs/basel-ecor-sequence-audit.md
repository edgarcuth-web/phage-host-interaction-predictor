# BASEL–ECOR accession and sequence audit

The 52 phage accession joins and joint sequence audit are complete. The [separately frozen v2 assembly policy](basel-ecor-controls.md) resolves the host-build/source choices for two fixed controls and records their negative composition result. Learned model fitting needs its own training protocol. This is the published, inspected reference from the [frozen v1 protocol](basel-ecor-identity-protocol-freeze.md), not an untouched release set.

## Verified inputs

The [new resource lock](../data/studies/basel_ecor_sequences.resources.v1.json) pins 80 downloaded files: 52 phage FASTAs, 25 published ECOR WGS sets, two original supplements and the existing source CSV. The runner also re-verifies the 94 training phage FASTAs and 25 training host FASTAs against the existing resource lock. Original locks, protocol v1, labels and weights are unchanged.

All 52 CSV phage IDs join to the same accession in both [Maffei et al. S5 Table](https://doi.org/10.1371/journal.pbio.3001424.s005) and [Noonan et al. Supplementary Table 8](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-026-02482-5/MediaObjects/41564_2026_2482_MOESM1_ESM.pdf). The retrieved versions are all **`.1`**. Each FASTA is one complete record with the expected accession.version and original S5 length. Bas02 and Bas09 have updated taxonomy in ENA headers; their explicit `strain BasXX` aliases are retained rather than treated as different phages.

See [all 52 accession receipts](../data/studies/results/basel_ecor_sequence_audit.v1/basel_accessions.csv) and the [machine-readable report](../data/studies/results/basel_ecor_sequence_audit.v1/report.json). Raw genomes, phenotype CSV and supplements are not copied into Git.

## Joint 94 + 52 phage comparison

The existing minimap2 **2.31-r1302**, `asm20 -c -P --seed 42 -t 2` rule was applied to all 146 genomes together. Alignment blocks, bilateral non-overlapping coverage and transitive components use the existing implementation. Thresholds were not selected from phenotype labels.

| Identity; coverage at least 85% of both genomes | All-panel groups | Direct training–BASEL edges | BASEL phages in a component containing training phages |
| --- | ---: | ---: | ---: |
| 95% primary | 120 | 0 | 0 / 52 |
| 90% sensitivity | 88 | 7 | 8 / 52 |
| 99% sensitivity | 136 | 0 | 0 / 52 |

There are **0 cross-panel canonical exact-assembly matches**. At 90%, the eight affected candidates are **Bas35, Bas36, Bas37, Bas39, Bas41, Bas42, Bas60 and Bas61**. Only Bas37, Bas60 and Bas61 have direct cross-panel edges; the others connect transitively. All 52 remain in the candidate panel.

These findings support separation at the declared 95% operational boundary, not absence of shared genes, distant ancestry or functional similarity. Canonical hashes ignore FASTA headers, case, contig order and strand, but retain segmentation and ambiguous bases and do not normalize circular rotation. Missing minimap2 edges are not proofs of unrelatedness.

## Host names match; assembled sequences do not

All **25/25 host IDs** remain represented in the training roster. However, **0/25 canonical host assembly hashes match** the WGS assemblies cited in Supplementary Table 9. For example, ECOR12 has 245 training contigs versus 218 in the cited QOWX assembly. Each published WGS set has the expected strain identifiers, accession prefix and scaffold count. The [host receipt](../data/studies/results/basel_ecor_sequence_audit.v1/host_assemblies.csv) records both raw-byte and canonical identities, counts and lengths.

This establishes different assembled sequence inputs; it does **not** establish a different biological strain, a host-name error or cold-host generalization. Existing host annotations or cached features cannot be attached to the published assemblies by name alone. A later scoring run must declare which assembly supplies each feature and regenerate or verify those features against that exact content.

## Source discrepancies retained

| Source issue | Evidence and treatment |
| --- | --- |
| Supplementary matrix roster differs from the pinned CSV | Table 8 lists 57 Bas IDs. Bas38, Bas44, Bas46 and Bas48 are marked for the interaction matrix but absent from the 52-phage CSV roster. Bas18 is also absent but is marked only for RB-TnSeq. Keep the locked 52-candidate CSV panel; do not infer extra observations or silently expand it. |
| Bas64 genome length differs | Original S5 and accession MZ501081.1 contain **39,842 bp**; GenoPHI Table 8 prints **39,451 bp**. Preserve both numbers and use the verified source genome in the identity receipt. No phenotype is reassigned. |

The CSV is still byte-verified as 1,240 observed pairs, 302 positive, 938 negative and 60 unknown combinations across 25 × 52. This audit does not explain why the four supplementary matrix phages are absent, infer their labels or replace unknown outcomes with negatives.

## Reproduce

Install the project dependencies, Poppler `pdftotext`, and the exact minimap2 revision from [the grouping lock](../data/studies/phage_grouping.v1.json). Then run:

```bash
python scripts/audit_basel_sequences.py \
  --fetch --data-root artifacts/basel_identity_inputs \
  --minimap2 /path/to/minimap2-2.31/minimap2 \
  --out artifacts/basel_sequence_audit_new
```

Omit `--fetch` to use verified local inputs. Output must be a new directory. Exit **2** means the audit completed but model scoring is still blocked. The report binds input locks, code hashes, binary identity, tool versions and output hashes. FASTA and sorted PAF outputs stay local; the committed report and CSV receipts are sufficient to inspect the joins and resulting groups. Source supplements and genomes require no new Python spreadsheet/PDF package.

The **known-host assembly policy**, 52-phage CSV scope and Bas64 sequence choice are now specified for the two controls in [protocol v2](basel-ecor-controls.md). The reason for the supplementary matrix discrepancy remains unknown. Neither this audit nor published-source model results alone can promote a model to release validation.
