"""Small NVIDIA-hosted forward client, separate from local-checkpoint caches."""
from __future__ import annotations

import base64
from dataclasses import dataclass
from hashlib import sha256
import io
import json
import urllib.error
import urllib.request
import zipfile

import numpy as np

FORWARD_URL = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/forward"
LAYER = "blocks.28.mlp.l3"
DIMENSION = 8192


class NvidiaForwardError(RuntimeError):
    """Errors exclude request headers, credentials and arbitrary server bodies."""


@dataclass
class ForwardResult:
    tokens: np.ndarray
    metadata: dict


def decode_layer(document, sequence_length, *, max_tensor_bytes):
    """Decode MIME-wrapped base64 NPZ; require the observed hosted output layout."""
    value = document.get("data")
    if not isinstance(value, str):
        raise NvidiaForwardError("NVIDIA response lacks base64 tensor data")
    try:
        payload = base64.b64decode("".join(value.split()), validate=True)
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            members = archive.infolist()
            if (len(members) != 1 or members[0].filename != LAYER + ".output.npy"
                    or members[0].file_size > max_tensor_bytes):
                raise NvidiaForwardError("unexpected or oversized NVIDIA tensor archive")
        with np.load(io.BytesIO(payload), allow_pickle=False) as arrays:
            values = arrays[LAYER + ".output"]
            if values.shape != (1, sequence_length, DIMENSION):
                raise NvidiaForwardError("NVIDIA layer token count or hidden dimension changed")
            if values.dtype.kind != "f" or not np.isfinite(values).all():
                raise NvidiaForwardError("NVIDIA layer must contain finite floating-point values")
            return values[0]
    except NvidiaForwardError:
        raise
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
        raise NvidiaForwardError("cannot decode NVIDIA tensor payload") from None


def forward(sequence, api_key, *, timeout=55, max_response_bytes=32 * 1024 * 1024,
            max_tensor_bytes=128 * 1024 * 1024):
    """One authenticated request, with bounded downloads and no automatic retries.

    This client does not claim a hosted service exposes local weight/runtime hashes.
    It cannot produce the local-checkpoint deployment lock from the v1 experiment.
    """
    if not isinstance(sequence, str) or not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("sequence must be non-empty uppercase ACGTN")
    if not isinstance(api_key, str) or not api_key.strip() or any(c.isspace() for c in api_key):
        raise ValueError("a non-empty NVIDIA API key without whitespace is required")
    if max_response_bytes < 1 or max_tensor_bytes < 1 or timeout <= 0:
        raise ValueError("request limits must be positive")
    body = json.dumps({"sequence": sequence, "output_layers": [LAYER]}).encode("ascii")
    request = urllib.request.Request(
        FORWARD_URL, data=body, method="POST",
        headers={"Authorization": "Bearer " + api_key, "Content-Type": "application/json",
                 "Accept": "application/json", "NVCF-POLL-SECONDS": "10"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            status = response.status
            request_id = response.headers.get("NVCF-REQID")
            if status != 200:
                raise NvidiaForwardError(f"NVIDIA forward returned HTTP {status}; no retry or polling was submitted")
            if "application/json" not in response.headers.get("Content-Type", ""):
                raise NvidiaForwardError("NVIDIA forward returned an unexpected content type")
            length = response.headers.get("Content-Length")
            if length is not None and int(length) > max_response_bytes:
                raise NvidiaForwardError("NVIDIA response exceeds the configured download limit")
            raw = response.read(max_response_bytes + 1)
            if len(raw) > max_response_bytes:
                raise NvidiaForwardError("NVIDIA response exceeds the configured download limit")
    except NvidiaForwardError:
        raise
    except urllib.error.HTTPError as exc:
        # Never include arbitrary remote messages: they can echo credentials.
        raise NvidiaForwardError(f"NVIDIA forward returned HTTP {exc.code}; no retry was submitted") from None
    except (OSError, urllib.error.URLError, ValueError):
        raise NvidiaForwardError("NVIDIA forward connection or response failed; no retry was submitted") from None
    try:
        document = json.loads(raw)
    except (ValueError, UnicodeError):
        raise NvidiaForwardError("NVIDIA forward returned invalid JSON") from None
    if not isinstance(document, dict):
        raise NvidiaForwardError("NVIDIA response must be an object")
    values = decode_layer(document, len(sequence), max_tensor_bytes=max_tensor_bytes)
    return ForwardResult(values, {
        "endpoint": FORWARD_URL, "layer": LAYER, "request_id": request_id,
        "http_status": status, "elapsed_ms": document.get("elapsed_ms"),
        "sequence_length": len(sequence), "input_sha256": sha256(sequence.encode("ascii")).hexdigest(),
        "response_bytes": len(raw), "response_sha256": sha256(raw).hexdigest(),
        "shape": list(values.shape), "returned_dtype": str(values.dtype),
        "finite": True, "weights_sha256": None, "server_runtime_version": None,
    })
