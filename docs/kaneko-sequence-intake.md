# Kaneko sequence sources and basic QC

The [workbook audit](kaneko-workbook-intake.md) now has downloaded sequence
sources for **143 host assemblies and four phage genomes**. This follow-up checks
source links and sequence content. It does not admit a training/evaluation cohort.
The [source manifest](../data/studies/kaneko_sequence_intake.v1.json) pins the
requests and files; [numeric QA](../data/studies/results/kaneko_sequence_intake.v1/NUMERIC_QA.json)
records the measured metrics and exact-hash comparisons.

## Sources actually retrieved

| Source | Retrieved and checked | Remaining gap |
| --- | --- | --- |
| 144 original host BioSamples | 143 compressed assemblies from EMBL-EBI third-party project PRJEB126751; every analysis XML references the original BioSample and DRR run | Host `127`, SAMD00774523 / DRR550410, has paired reads but no assembly in the inspected indexes |
| 42 study-phage names | Four versioned public genomes, with sequence MD5 and length checks | The 38 BioSamples in PRJDB18044 have read routes but no assembled records in the inspected indexes |

Searching only PRJDB18063/PRJDB18044 returned no assemblies. Searching the same
182 BioSamples across projects found 143 host assembly and analysis records.
The selected EMBL-EBI FASTAs are third-party reassemblies, not the paper's original
assemblies. Their analysis XMLs identify Shovill v1.1.0 with SPAdes v3.15.5 and
explicitly link the original sequencing runs. Files and XML provenance are pinned
separately. Project/sample searches are bounded evidence, not proof that no other
assemblies exist anywhere.

The 143 compressed host assemblies total **209,011,797 bytes**. Original FASTQ
metadata provide this separate acquisition inventory:

| Reads | Paired samples | FASTQ files | Compressed bytes |
| --- | ---: | ---: | ---: |
| Hosts | 144 | 288 | 102,478,826,136 |
| Phages | 38 | 76 | 27,241,110,587 |
| Total | 182 | 364 | 129,719,936,723 |

No FASTQ was downloaded or assembled in this follow-up. Host 127 is the one
DNBSEQ-G400 host record; the other 143 host runs are Illumina HiSeq X. This is an
observed coverage gap, not evidence that sequencing platform caused the absence.
Do not silently exclude host 127 or invent its genome.

## Sequence checks

The audit verifies downloaded size/hash, archive MD5, gzip decoding, nonempty
IUPAC DNA, source BioSample/run links, contig count, total length, N50, ambiguity
and GC. It preserves all contigs and ambiguous bases, with no quality-based
filter chosen from these findings.

The 143 host assemblies span 4,799,870–5,758,484 bases and 52–472 contigs.
Median N50 is 191,251 bases (range 43,787–872,287); none contains non-ACGT bases.
All 147 acquired assembly digests are distinct, with no exact match to the frozen
177-host/94-phage reference. These are exact-assembly findings, not evidence of
strain independence.

| Phage | Versioned record | Bases | Non-ACGT bases |
| --- | --- | ---: | ---: |
| ΦWec189 | LC739538.1 | 135,501 | 0 |
| ΦWec191 | LC739540.1 | 137,785 | 10 |
| ΦWec193 | LC739541.1 | 138,673 | 0 |
| ΦWec196 | LC739542.1 | 138,673 | 0 |

For exact overlap, the command reuses `cohort.canonical_assembly_digest` and the
frozen 177-host/94-phage feature manifest. That digest normalizes FASTA headers,
case, contig order and strand. It does not establish distinct strains when
contig boundaries, circular starting positions or assemblies differ. Near-duplicate
relationships, contamination and completeness still require separate assessment.
Metadata/source links also do not replace verification of the experimental stocks.

## Reproduce

After installing the project, fetch the pinned metadata, provenance XMLs and
approximately 210 MB of compressed assemblies/public FASTAs:

```bash
python scripts/audit_kaneko_sequences.py \
  --fetch --data-root artifacts/kaneko_sequences_new \
  --out artifacts/kaneko_sequence_qa_new.json
```

For an offline replay, omit `--fetch`. The output file must be new; fetching also
requires a new data directory. Metadata replay checks both the column set and all
canonical row contents, allowing ordering differences. An HTTP/error response
cannot stand in for an empty assembly search. Genome and XML byte changes fail
verification. Source sequences, archives, XMLs and bulk metadata remain outside
Git; the report retains hashes and QC findings.

## Next

The [42-candidate protocol](kaneko-baseline-protocol.md) now freezes all 143
available hosts, 27 genomic groups and the grade-3 complete-clearing readout.
Next is the fixed prevalence/host 4-mer comparison, with the seven-host/five-group
test limitation retained. Complete broader stock, contamination/completeness,
reuse and cross-source near-duplicate reviews for any wider claim. Host 127 and
the 38 phage assemblies remain acquisition tasks for later scopes. No
susceptibility model fitting or scoring occurred here.
