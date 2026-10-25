# NVIDIA-hosted Evo2-40B: verified forward access

The authenticated hosted forward API works with `blocks.28.mlp.l3`. This is a
connectivity and transfer-feasibility probe, not the three-host genome pilot or
the host-prediction benchmark. No generated genome, phage design or prediction
is emitted.

The fixed endpoint is:

```text
https://health.api.nvidia.com/v1/biology/arc/evo2-40b/forward
```

NVIDIA documents [forward tensor capture](https://docs.nvidia.com/nim/bionemo/evo2/1.0.2/endpoints.html#run-model-forward-pass-and-save-layers-outputs)
with `sequence` and `output_layers`, returning base64-encoded NPZ data. The hosted
service accepted the existing protocol's exact layer name; no layer substitution
was made. Newer self-hosted NIM documentation uses different module names, so
do not infer that a renamed module has been validated against this service.

## Observed results

The [run record](../data/studies/results/nvidia_evo2_api_probe.v1.json) retains
request IDs, input/output hashes, code hashes, shapes and response sizes. Four
recorded requests use two fixed 16-base control strings, repeat the first string,
and repeat its motif to make a 256-base input. Preliminary connection probes are
not included in this four-request record. No host genomes were submitted.

| Request | HTTP | Tensor shape before removing batch | Dtype | Server elapsed | JSON bytes |
| --- | ---: | --- | --- | ---: | ---: |
| Control A, 16 bases | 200 | 1 × 16 × 8,192 | float64 | 301 ms | 349,720 |
| Control B, 16 bases | 200 | 1 × 16 × 8,192 | float64 | 303 ms | 349,928 |
| Repeat A, 16 bases | 200 | 1 × 16 × 8,192 | float64 | 296 ms | 349,720 |
| 256 bases | 200 | 1 × 256 × 8,192 | float64 | 1,716 ms | 5,589,613 |

Every tensor is finite and varies along the token axis. Repeating control A
reproduced its tensor exactly. Different controls changed the mean-pooled vector
by about 2.07% relative L2; their last-base vectors also differed despite both
strings ending in T. Absolute output values are small (roughly 1e-11); this alone
does not establish an error, and these values must not be downcast to float16.
These controls establish input sensitivity, not biological predictive quality
or agreement with local checkpoint inference.

## Transfer is the current bottleneck

The 256-base response is 5.59 MB. Linear extrapolation to the unchanged plan
(32 windows × 8,192 bases × both strands) gives approximately:

| Scope | Estimated JSON download |
| --- | ---: |
| Three-host pilot | 34.3 GB |
| Full 177-host cohort | 2.03 TB |

These are rough estimates from short repeated controls, not measured 8,192-base
responses. Compression, sequence content and service limits can change them.
One uncompressed 8,192-base float64 tensor alone would occupy 512 MiB. Pooling
after downloading reduces retained cache size but does not reduce network
traffic. The documented API does not expose a pooling parameter; server-side
pooling cannot be assumed.

The client therefore limits each download to 32 MiB and decoded arrays to
128 MiB by default. It deliberately does not start the original three-host/full
extraction. A hosted experiment needs its own versioned provider/protocol
amendment, explicit transfer budget and cache format. If window coverage changes,
all affected representation comparisons need the amended protocol before any
outcome analysis. The local-checkpoint v1 protocol remains intact.

Hosted weights and runtime hashes are not exposed by this endpoint. The client
records those identities as unknown; it does not invent a deployment lock or
admit hosted results as a local-checkpoint v1 cache. No claim is made about hosted
model revision pinning, hardware memory use, throughput on full genomes or
8,192-base request support.

## Run a bounded probe

```bash
python scripts/probe_nvidia_evo2.py --outdir artifacts/nvidia_evo2_api_probe
```

The script prompts for the key without echo. A process-injected `NVIDIA_API_KEY`
is also supported. The key is never a command-line argument, artifact field or
repository setting. Do not put credentials in committed files. The four requests
are authenticated inference calls; there are no automatic retries or asynchronous
polling. Errors report HTTP status without arbitrary response bodies that could
echo a credential. An existing output directory is rejected.

The small client is usable separately:

```python
import getpass
from phage_cocktail.nvidia_evo2 import forward

result = forward("ACGTTGCATGTCGCAT", getpass.getpass("NVIDIA API key: "))
tokens = result.tokens  # [16, 8192]; observed float64
```

The response decoder handles the observed MIME base64 newlines, refuses pickle
loading and requires the exact layer, batch/token layout, width and finite
floating-point output. Contract tests cover these boundaries, bounded downloads,
one-request behavior, credential-free error messages and metadata.
