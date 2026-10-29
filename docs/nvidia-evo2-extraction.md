# Resumable hosted Evo2 extraction

The executor receives NVIDIA forward responses, pools each orientation, and
publishes small immutable checkpoints to a local directory or an S3 bucket.
It runs on CPU; the model remains on NVIDIA's service. A storage bucket needs a
separate Python worker. This implementation does not deploy a worker, create a
bucket, measure a live S3 service, or start the full 177-host extraction.

The default is the same first three sorted hosts and first full sampled window
as the measured engineering pilot. Both strands are retained. An explicit
`original_32_including_short_tails` selection implements the original systematic
32-window sampling, including partial tails; it is not whole-genome coverage.
No labels, fitting, ranking, or prediction metrics enter this executor.

## Running Phase 1 without a preexisting local FASTA folder

The manual [Evo2 full-host offline preflight](../.github/workflows/evo2-host-offline-preflight.yml)
GitHub Actions workflow reconstructs the missing input folder directly from the
177 published and SHA-256-pinned host FASTA URLs. It has no NVIDIA credentials,
performs no inference, and is not a budget approval. The downloaded raw genomes
stay on the ephemeral Actions runner and are **not uploaded**.

After the workflow has been merged into the default branch, open **Actions →
Evo2 full-host offline preflight → Run workflow**. A PR that adds or changes the
workflow can also exercise it. Its run uploads a seven-day artifact named
`full-host-evo2-offline-preflight`, containing only:

- `host_fasta_source_receipt.json` — SHA-256-locked source download/verification receipt;
- `full_hosted_plan_zero_budget.json` — immutable zero-call, zero-response-budget plan;
- `full_hosted_plan_preflight.json` — verified host coverage, request count and projected transfer.

This step fetches about 910 MB of public source FASTAs on the Actions runner,
but transfers **zero token-level Evo2 outputs** and makes **zero NVIDIA calls**.
If a public source is unavailable or differs from its pinned SHA-256/size, the
workflow fails closed; it never substitutes a sequence or runs provider inference.
The source fetch helper uses limited retries only for transient HTTP/network
errors and will not overwrite an existing mismatched resource.

The output plan is bound to the exact code and immutable manifest of that run.
Do not edit its budgets or assume it can be replayed under changed code. A
separate reviewed execution plan would be necessary for any future paid run.

## Full-cohort budget review (no provider calls)

A full hosted `original_32_including_short_tails` plan for the 177 previously
verified host assemblies requires an unusually large intermediate tensor
transfer. The existing independent CPU sampling audit predicts **11,328
orientation requests** and **about 1.425 TB of JSON body downloads**. It
retains only about 5.8 MB of final float32 host-vector data; the large transfer
occurs before local token pooling.

The current executor reserves the response cap plus one sentinel byte for
each live attempt. With the default 256 MiB cap, the **no-reuse full-cap
reservation** for 11,328 requests would be **3,040,836,856,896 bytes (3.04 TB)**.
This is a conservative *budget reservation*, not a predicted actual download,
invoice, provider quota or authorization. The independent published estimate
and the exact-plan estimate must match before an execution plan is reviewed.

Create a **zero-request** full-panel plan on a machine with the verified
GenoPHI FASTAs already downloaded:

```bash
python scripts/run_nvidia_host_extraction.py plan \
  --all-hosts --selection original_32_including_short_tails \
  --plan artifacts/full_hosted_plan_zero_budget.json \
  --max-new-calls 0 --max-new-response-bytes 0

python scripts/preflight_nvidia_host_extraction.py \
  --plan artifacts/full_hosted_plan_zero_budget.json \
  --out artifacts/full_hosted_plan_preflight.json \
  --verify-sources
```

The preflight cross-checks the actual plan job lengths and calls against the
independently committed sampling projection. It reports per-host coverage,
estimated download size, conservative full-cap reservation and the frozen
budget's sufficiency. `--verify-sources` also verifies each local FASTA,
the source resource lock and executor code. Without this flag the output
explicitly says source FASTA verification was not performed. Neither command
reads an NVIDIA key or issues a provider POST.

The zero-budget plan cannot later be edited to increase limits. Any proposed
paid run must receive a **new immutable plan**, separately reviewed against the
preflight and with appropriate external quota/cost approval. Do not run with
`--allow-live` solely because a preflight was produced. The independent
sampling experiment also found reduced 4-mer retrieval performance on
sampled windows; that result is a cost/information warning, not an Evo2
performance result.

## Hosted evaluation boundary

The hosted executor and the local-checkpoint Evo2 v1 cache intentionally use
different pooling arithmetic. Hosted outputs therefore **must not** be copied or
renamed into the local cache layout. A completed full-panel hosted run is
evaluated through the explicit hosted adapter, which verifies the extraction
plan, immutable receipts, final host matrix and exact source FASTAs before
reusing the existing host splits, candidate panel, estimator and retrieval
metrics.

The hosted study uses distinct method names
`host_hosted_evo2_logistic` and `hybrid_hosted_evo2_logistic`. Its output
protocol records the hosted extraction plan as the representation definition and
retains `local_v1_cache_compatible: false`.

After a full `original_32_including_short_tails` extraction is complete:

```bash
python scripts/run_hosted_evo2_host_comparison.py \
  --plan artifacts/full_hosted_extraction_plan.json \
  --store s3://YOUR_BUCKET/evo2/full-run-001 \
  --outdir artifacts/hosted_evo2_host_comparison_v1
```

Use a local checkpoint directory instead of the S3 URI when applicable. The
comparison command performs no NVIDIA inference. It refuses incomplete runs,
partial host panels, the one-window engineering selection, changed source
FASTA/code identity, corrupted final vectors, or a plan that differs from the
stored immutable plan. Hosted checkpoint and server-runtime hashes remain
unavailable from the service and are reported as a limitation rather than
invented.

## Offline acceptance using real existing results

Use the verified source FASTAs from the preceding study and the current checkout:

```bash
pip install -e '.[dev]'
python scripts/verify_nvidia_host_extraction.py \
  --outdir artifacts/nvidia_extraction_replay \
  --qa-out artifacts/nvidia_extraction_replay_qa.json
```

This command imports the six committed real orientation vectors, pauses after
two saved jobs, resumes the remaining four, reconstructs the three host vectors,
and compares every value with the original pilot. A final completed-run replay
must preserve every stored object hash. Request/key callbacks raise if accessed.
All request budgets are zero. No NVIDIA key or S3 account is used. Synthetic CPU
fixtures in the test suite separately exercise failures, budgets and claim races;
they are not biological results.

The [recorded replay](../data/studies/results/nvidia_host_extraction_replay.v1.json)
uses published implementation `db062825be52b9fa66cfb4f72d4a79b9b7ff5680`.
It passes the two-job pause, six-job resume, exact three-host reconstruction,
key-free verification and hash-preserving completed replay. Its 15 checkpoint
objects total 380,370 bytes, including the plan, base64 receipts and final matrix.
New NVIDIA calls are zero; no live S3 service was tested. The full local suite
passes 145 tests and eight subtests, including the Boto3 stubbed transport test.

## Plan first, then run

```bash
python scripts/run_nvidia_host_extraction.py plan \
  --plan artifacts/extraction_plan.json

python scripts/run_nvidia_host_extraction.py run \
  --plan artifacts/extraction_plan.json \
  --store artifacts/extraction_objects \
  --pilot-cache data/studies/results/nvidia_evo2_host_pilot.v1 \
  --max-jobs 2

# Same arguments without --max-jobs complete the remaining cached imports.
python scripts/run_nvidia_host_extraction.py run \
  --plan artifacts/extraction_plan.json \
  --store artifacts/extraction_objects \
  --pilot-cache data/studies/results/nvidia_evo2_host_pilot.v1

python scripts/run_nvidia_host_extraction.py verify \
  --plan artifacts/extraction_plan.json \
  --store artifacts/extraction_objects
```

`plan` verifies selected raw FASTAs and assembly identities before writing a new
manifest. It freezes exact input hashes, host/window/strand identities, source
file hashes, executor/package code hashes, endpoint, layer, pooling and budgets.
It contains no DNA strings or credentials. `run` reconstructs and verifies every
input before reading a key or making a request. Keep the same code checkout,
manifest and verified source files on the worker. `verify` needs only the plan
and stored checkpoint objects, without source FASTAs, NVIDIA key or GPU.

The default plan budgets permit **zero** new requests. Live execution requires
both nonzero budgets frozen in a new plan and `run --allow-live`. Only then is
`NVIDIA_API_KEY` read, or a hidden interactive prompt shown, when an uncached job
actually needs a key. Neither plans nor checkpoints store that key.

For example, this prepares a one-new-request engineering budget for one host's
original sampled windows; it does not submit any inference:

```bash
python scripts/run_nvidia_host_extraction.py plan \
  --plan artifacts/one_new_request_plan.json \
  --hosts 370D --selection original_32_including_short_tails \
  --max-new-calls 1 --max-new-response-bytes 268435457
```

Use `--hosts ID1,ID2` for an explicit panel or `--all-hosts` to select all 177.
Review the generated job count and budget before live execution. Unspecified
selection stays at one full window even with `--all-hosts`. `--max-jobs` pauses
after that many newly saved imports/live results without changing the plan.
An incomplete run returns a stop reason: live requests disabled, per-invocation
job limit, or frozen budget exhausted. It never silently increases a budget.

## Cloud storage

```bash
pip install -e '.[cloud]'
python scripts/run_nvidia_host_extraction.py run \
  --plan artifacts/extraction_plan.json \
  --store s3://YOUR_BUCKET/evo2/run-001 \
  --pilot-cache data/studies/results/nvidia_evo2_host_pilot.v1
```

This example still has zero NVIDIA request budget and no `--allow-live`.
Supply bucket credentials through the worker's normal Boto3 credential chain,
preferably its assigned role. The run prefix needs `GetObject` and `PutObject`;
the executor neither lists nor deletes objects, and does not create buckets.
Use a dedicated prefix, keep its objects immutable, and do not delete attempt
claims to make a failed job run again. An optional `--s3-endpoint-url` accepts an
HTTPS endpoint without embedded credentials for a compatible service.

The backend uses `PutObject(IfNoneMatch="*")`, SigV4, and one SDK attempt per
operation. At execution startup a small persistent probe checks that a second
conditional create is refused and the first value remains readable. A service
that rejects or ignores these conditions cannot be used. Strong read-after-write
consistency is also required. Compatibility with a specific non-AWS service must
be checked on that service; a successful CPU/stubbed test is not a live cloud test.

See official [conditional writes](https://docs.aws.amazon.com/AmazonS3/latest/userguide/conditional-writes.html),
[PutObject](https://docs.aws.amazon.com/boto3/latest/reference/services/s3/client/put_object.html),
[SDK retries](https://docs.aws.amazon.com/boto3/latest/guide/retries.html), and
[credentials](https://docs.aws.amazon.com/boto3/latest/guide/credentials.html).

## Durable publication and recovery

For every ordered job the store holds:

- `jobs/NNNNN/claim.json`: input identity, plan hash, origin and a unique live
  attempt ID/download reservation, or the verified pilot report hash for imports.
- `jobs/NNNNN/result.json`: claim hash, input identity, response metadata, vector
  SHA-256 and the base64-encoded float32 NPY vector. Keeping them in one object
  prevents a completed receipt from being published without its vector.

The plan is `plan.json`. After every job completes, `host_vectors.npy` and
`report.json` are derived from the persisted receipts. Verification reconstructs
the matrix, checks every receipt and reconciles the report and budgets.
Final export publication can be retried from complete receipts without inference.

A live claim is conditionally created **before** POST. A competing worker that
loses that create stops without submitting. A receipt appears only after the
response is decoded, checked, pooled and durably saved. All known jobs are checked
before further execution. A live claim without a valid result stops the entire
run, including after HTTP 202, timeout, decode failure, process termination or
storage failure after inference. There is no automatic NVIDIA retry or polling,
and no automatic claim removal. Provider status/manual reconciliation is needed
for ambiguous attempts; this version supplies no resubmit/reconciliation command.

An interrupted pilot import can finish using the same verified source cache,
because it does not submit inference. A completed cache runs again without that
pilot cache and without a key. Only the exact committed engineering pilot is an
automatic import source; arbitrary cross-plan cache migration is not implemented.
Use one worker per run. Conditional claims guard duplicate launches; this is not
a parallel scheduler or a distributed lease/recovery system. Do not edit a stored
plan to enlarge its budget. A changed protocol requires a new reviewed plan, with
existing and ambiguous attempts accounted for before any new extraction.

Local publication uses flushed/fsynced temporary files and atomic hard links,
then fsyncs the containing directory. Use a filesystem that supports these local
POSIX semantics. Unpublished temporary files left by a killed process are not
completion records.

## Budget and numerical semantics

Every live claim reserves its response cap **plus one byte** for the bounded
client's oversize sentinel read. Reservations stay charged across resumes and
are not refunded when a response is smaller. Thus a 268,435,457-byte budget permits
one request with a 256-MiB response cap. The journal reports both conservative
reserved bytes and actual bytes of successful new responses. Imported vectors
use no new NVIDIA budget. These limits bound application body reads/attempts,
not TLS/header traffic, cloud storage requests, provider billing, peak memory,
or an absolute wall-clock deadline. S3 transfers are separate from NVIDIA budgets.

The 180-second timeout is an I/O timeout and the server wait is 60 seconds.
Each decoded member is capped at 640 MiB; one 8192-base tensor was observed as
512 MiB float64. JSON, base64 and decoding intermediates add memory overhead.
Measure peak memory on the actual worker before choosing concurrency/capacity.
No full response or position-level tensor is stored by this executor.

Each orientation uses float64 token means, saved as float32. Two orientations
are averaged in float64 and saved as a float32 window vector; window vectors are
weighted by original base length and accumulated in float64, then the host vector
is saved as float32. This matches the measured one-window pilot. It differs from
the earlier local-checkpoint protocol's token-mean arithmetic and has its own
identity. Mean pooling loses position-level information. All 177 final vectors
have 5.8 MB of numerical payload; base64 receipts and metadata add overhead to
the approximately 371 MB of orientation-vector numerical payload in the original
32-window plan. Large NVIDIA-to-worker responses still have to be transferred.

Hosted checkpoint/runtime hashes remain unknown. This engineering cache is
explicitly incompatible with the local-checkpoint v1 evaluation cache. A hosted
prediction study still needs a separately frozen protocol and evaluator adapter;
no model performance or clinical result follows from completing this executor.
