"""Hosted forward contract checks, with no real credentials or network calls."""
import base64
import io
import json
import urllib.error

import numpy as np
import pytest

from phage_cocktail.nvidia_evo2 import DIMENSION, FORWARD_URL, LAYER, NvidiaForwardError, decode_layer, forward


def document(values, *, layer=LAYER):
    stream = io.BytesIO()
    np.savez_compressed(stream, **{layer + ".output": values})
    # Real hosted responses contain MIME-style base64 newlines.
    return {"data": base64.encodebytes(stream.getvalue()).decode("ascii"), "elapsed_ms": 123}


def test_mime_base64_npz_retains_real_hosted_dtype_and_token_axis():
    expected = np.arange(2 * DIMENSION, dtype=np.float64).reshape(1, 2, DIMENSION)
    result = decode_layer(document(expected), 2, max_tensor_bytes=1024 * 1024)
    np.testing.assert_array_equal(result, expected[0])
    assert result.dtype == np.float64


@pytest.mark.parametrize("values", [np.zeros((2, 1, DIMENSION)), np.zeros((1, 2, DIMENSION - 1)),
                                    np.full((1, 2, DIMENSION), np.nan)])
def test_wrong_layout_dimension_or_nonfinite_tensors_are_rejected(values):
    with pytest.raises(NvidiaForwardError):
        decode_layer(document(values), 2, max_tensor_bytes=1024 * 1024)


def test_other_layer_and_oversized_archive_are_rejected():
    values = np.zeros((1, 2, DIMENSION))
    with pytest.raises(NvidiaForwardError, match="archive"):
        decode_layer(document(values, layer="embedding_layer"), 2, max_tensor_bytes=1024 * 1024)
    with pytest.raises(NvidiaForwardError, match="archive"):
        decode_layer(document(values), 2, max_tensor_bytes=16)
    with pytest.raises(NvidiaForwardError, match="decode"):
        decode_layer({"data": "***bad***"}, 2, max_tensor_bytes=1024)


class Response:
    def __init__(self, body, *, status=200, headers=None):
        self.body = body
        self.status = status
        self.headers = headers or {"Content-Type": "application/json", "NVCF-REQID": "fixture-request"}
        self.read_sizes = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self, size):
        self.read_sizes.append(size)
        return self.body[:size]


def test_one_authenticated_request_without_key_in_metadata(monkeypatch):
    response = Response(json.dumps(document(np.zeros((1, 2, DIMENSION)))).encode())
    requests = []
    def open_request(request, *, timeout):
        requests.append(request)
        return response
    monkeypatch.setattr("urllib.request.urlopen", open_request)
    result = forward("AC", "fixture-secret")
    assert len(requests) == 1 and requests[0].full_url == FORWARD_URL
    assert requests[0].get_header("Authorization") == "Bearer fixture-secret"
    assert json.loads(requests[0].data) == {"sequence": "AC", "output_layers": [LAYER]}
    assert "fixture-secret" not in json.dumps(result.metadata)
    assert result.metadata["weights_sha256"] is None
    assert result.metadata["shape"] == [2, DIMENSION]


def test_http_errors_never_echo_remote_body_or_retry(monkeypatch):
    calls = []
    def fail(*args, **kwargs):
        calls.append(1)
        raise urllib.error.HTTPError(FORWARD_URL, 401, "fixture-secret", {}, io.BytesIO(b"fixture-secret"))
    monkeypatch.setattr("urllib.request.urlopen", fail)
    with pytest.raises(NvidiaForwardError, match="HTTP 401") as error:
        forward("AC", "fixture-secret")
    assert "fixture-secret" not in str(error.value) and calls == [1]


def test_download_limit_and_async_response_do_not_trigger_more_calls(monkeypatch):
    response = Response(b"x" * 33)
    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: response)
    with pytest.raises(NvidiaForwardError, match="download limit"):
        forward("AC", "fixture-secret", max_response_bytes=32)
    assert response.read_sizes == [33]
    response.status = 202
    with pytest.raises(NvidiaForwardError, match="HTTP 202") as error:
        forward("AC", "fixture-secret")
    assert error.value.http_status == 202
    assert error.value.request_id == "fixture-request"


def test_server_wait_is_explicit_and_bounded(monkeypatch):
    requests = []
    response = Response(json.dumps(document(np.zeros((1, 2, DIMENSION)))).encode())
    def open_request(request, *, timeout):
        requests.append((request, timeout))
        return response
    monkeypatch.setattr("urllib.request.urlopen", open_request)
    forward("AC", "fixture-secret", poll_seconds=60, timeout=180)
    assert requests[0][0].get_header("Nvcf-poll-seconds") == "60"
    assert requests[0][1] == 180
    with pytest.raises(ValueError, match="between 1 and 60"):
        forward("AC", "fixture-secret", poll_seconds=61)
    assert len(requests) == 1


def test_invalid_input_is_rejected_before_network(monkeypatch):
    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: pytest.fail("network must not be called"))
    with pytest.raises(ValueError, match="ACGTN"):
        forward("hello", "fixture-secret")
    with pytest.raises(ValueError, match="API key"):
        forward("AC", "")
