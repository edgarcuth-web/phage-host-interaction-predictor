# Real-host hosted Evo2 feasibility and engineering pilot

This work tests whether the authenticated hosted API can process real frozen
host windows, pool outputs and recover cached results within an explicit budget.
No interaction labels, phenotype model, ranking metric or model promotion are
part of this pilot.

## One real 8192-base window

The [one-window protocol](../data/studies/nvidia_evo2_window_feasibility.v1.json)
and client were published at `ed417ea59d9b6982829ea716d52defb94b898d52` before the
request. The first lexicographically sorted frozen host is 370D. Its exact FASTA
bytes and canonical assembly identity were verified. The window is the first
full-length member of the existing canonical-contig 32-window sample, without
consulting outcomes. Its 316 ambiguous positions remain N.

The [record](../data/studies/results/nvidia_evo2_window_feasibility.v1.json)
reports HTTP 200 for `blocks.28.mlp.l3`: an 8,192 × 8,192 finite float64 tensor,
50,698 ms server time, 66.95 seconds client wall time and 178,653,288 JSON bytes.
The token mean was accumulated in float64 and stored as an 8,192-element float32
vector; maximum casting error was 1.36e-19. Full token outputs were released.
The socket I/O timeout is not an absolute end-to-end deadline or a memory bound.

Scaling this actual window to the original 32-window/two-strand plan gives about
34.30 GB for three hosts and 2.024 TB for 177 hosts. Other sequences can compress
differently. This confirms that long requests work, while the unmodified full
transfer remains large. It does not validate a full host representation.

## A separate engineering-only pilot

The [hosted pilot protocol](../data/studies/nvidia_evo2_host_pilot.v1.json) is
distinct from the local-checkpoint representation benchmark:

| Item | Frozen choice |
| --- | --- |
| Hosts | First three sorted cohort hosts: 370D, 381A, 536 |
| Sampling | One full 8,192-base window per host, first in the existing 32-window selection |
| Orientations | Forward and reverse complement, independently encoded |
| Layer | The same `blocks.28.mlp.l3`, width 8,192 |
| Pooling | Float64 token mean → float32 orientation vector; mean the two stored orientations in float64 → float32 host vector |
| Reuse | Import the exact verified 370D forward pooled result from the feasibility run |
| New inference budget | At most five POST attempts, no retries or status polling |
| New JSON transfer budget | At most 1 GiB cumulatively, at most 256 MiB per response |
| Identity | Hosted checkpoint/runtime hashes remain unknown; never manufacture a local deployment lock |

The selected-window assembly fractions are reported. One window is a small sample
of a bacterial genome; these vectors are not the original 32-window representation
and cannot justify an Evo2-versus-4-mer performance claim. 370D's incomplete
phenotype labels do not affect this label-free engineering exercise.

## Run and verify

All three pinned host FASTAs must already be present locally. The key is supplied
through a hidden prompt or a process-injected `NVIDIA_API_KEY`, never a command-line
argument or saved configuration.

```bash
# Restore the exact parent without another inference call.
mkdir -p artifacts/nvidia_evo2_window_feasibility_v1
cp data/studies/results/nvidia_evo2_window_feasibility.v1.json \
  artifacts/nvidia_evo2_window_feasibility_v1/report.json
cp data/studies/results/nvidia_evo2_host_pilot.v1/job_00.npy \
  artifacts/nvidia_evo2_window_feasibility_v1/pooled_mean.npy

python scripts/run_nvidia_host_pilot.py \
  --feasibility-cache artifacts/nvidia_evo2_window_feasibility_v1 \
  --outdir artifacts/nvidia_evo2_host_pilot_v1

# This path needs neither a key nor a network call for a completed cache.
python scripts/run_nvidia_host_pilot.py --verify-only \
  --feasibility-cache artifacts/nvidia_evo2_window_feasibility_v1 \
  --outdir artifacts/nvidia_evo2_host_pilot_v1

python scripts/verify_nvidia_host_pilot.py \
  --cache artifacts/nvidia_evo2_host_pilot_v1 \
  --out artifacts/nvidia_evo2_host_pilot_qa.json
```

The standalone `scripts/probe_nvidia_host_window.py --outdir NEW_DIRECTORY`
command makes a fresh feasibility inference request. Such a new run will have a
new report identity and needs an amended pilot protocol. The procedure above
instead restores the exact parent vector/report pinned by this pilot.

The journal marks a job attempted before entering the API call. Each completed
orientation has an atomically written vector, input identity, response metadata
and vector hash. Resume uses explicit `--resume`, unchanged code/protocol and
verified cached outputs. An unresolved attempt blocks further execution and
never triggers an automatic second inference. Fresh, unattempted jobs can resume
after verified completions; a fully completed cache can be verified without
any API request. Only one full token tensor is kept at a time.
Run a single writer process; concurrent launches are not coordinated by a lock.

An API error or download/decode limit stops the run immediately. The journal
retains evidence, and the aggregate budget counts successful newly downloaded
responses; an unresolved/error attempt blocks all subsequent requests. No promise
is made that an interrupted remote job has been cancelled.

## Measured pilot results

The hosted pilot protocol, implementation and pinned feasibility parent were
published at `1c0981f71b1e416032d3fb1f8342441e0c37d7a4` before the new pilot
inferences. The [complete cache](../data/studies/results/nvidia_evo2_host_pilot.v1/report.json)
retains all six pooled orientation vectors, the three pooled host vectors,
per-job metadata, journal and input/code/output hashes. Full token tensors and
source FASTAs are not redistributed.

| Host | Sampled assembly | Forward / reverse server seconds | Forward / reverse JSON MB |
| --- | ---: | ---: | ---: |
| 370D | 0.1552% | 50.70 / 50.56 | 178.65 / 178.67 |
| 381A | 0.1599% | 50.54 / 50.72 | 178.70 / 178.68 |
| 536 | 0.1684% | 50.80 / 50.72 | 178.72 / 178.72 |

All six oriented results returned HTTP 200, finite 8,192 × 8,192 float64
token matrices, and valid float32 pooled outputs. The first forward result was
reused exactly. Five new POSTs downloaded 893,489,254 JSON bytes (about 852.10 MiB),
within the 1 GiB new-download budget. Including the earlier reused response,
all six responses total 1,072,142,542 bytes. Seven retained vector files occupy
295,808 bytes in total. The three host vectors are distinct.

Independent offline QA checks every payload and recorded source-code hash,
verifies all six orientation vectors, reconstructs the three host vectors
exactly, validates the parent reuse and reconciles request/download accounting.
Completed-cache replay with a network callback that raises on any call makes
zero API calls and preserves every artifact hash. The QA record is committed as
`NUMERIC_QA.json`. All 121 tests and eight subtests pass; the implementation
snapshot also passed Python 3.11 and 3.12 CI.

Each host is represented here by only one window, covering roughly 0.16% of its
assembly. This is sufficient to validate the transport/pooling/cache path,
but does not test full-genome representations, unseen-host Hit@K or Evo2 superiority.
The unchanged 32-window comparison would still require a much larger hosted
budget or pooling at a controllable inference server. No full-cohort extraction,
benchmark fitting, calibration or model promotion was initiated.
