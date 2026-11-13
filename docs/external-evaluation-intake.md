# External evaluation intake: identity and provenance gate

Status: **software-only intake gate**, not a new dataset, a completed benchmark, or independent model validation. No external source is admitted by this change.

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

## Mandatory manual steps even if this gate finds no overlap

1. Resolve study aliases, source-data derivations, shared experiments and source licenses; verify original accessions and dataset checksums.
2. Audit **near-duplicate genomes** with a prespecified grouped sequence-similarity or phylogenetic review. Exact SHA-256 disjointness alone is insufficient.
3. Review assay type, threshold, positive/negative/unknown semantics, and alignment of source endpoints. Source MLC/spotting positivity cannot silently become productive infection.
4. Audit training and model-selection exposure, feature preprocessing, reused test data and relevant pretrained-model contamination.
5. Freeze a truly independent cohort, allowed claims, subgroup definitions and evaluation protocol **before running predictions**.

A clean report never sets `independent_release_validation_established` to true. This tool blocks some invalid claims; it does not replace independent scientific review or experimental validation.
