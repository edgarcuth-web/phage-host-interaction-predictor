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
