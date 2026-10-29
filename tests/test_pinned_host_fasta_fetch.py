"""Safe source-fetch tests; no public URL is requested in tests."""
from hashlib import sha256
import importlib.util
import io
from pathlib import Path
import urllib.error

import pytest

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail import nvidia_host_extraction as extraction
from test_nvidia_host_extraction import extraction as source_fixture


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/fetch_pinned_host_fastas.py"
SPEC = importlib.util.spec_from_file_location("fetch_pinned_host_fastas", SCRIPT)
fetcher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fetcher)


def prepare(extraction):
    repo, root, store = extraction
    lock_path = repo / extraction_module_path()
    data = read_json_object(lock_path)
    for item in data["resources"]:
        original = root / item["path"]
        final = root / "host" / Path(item["path"]).name
        final.parent.mkdir(parents=True, exist_ok=True)
        original.rename(final)
        item["path"] = final.relative_to(root).as_posix()
        item["url"] = "https://ndownloader.figshare.com/files/" + str(len(item["id"]))
    atomic_write_json(data, lock_path)
    return repo, root, store


def extraction_module_path():
    return extraction.GENOMES


def test_existing_locked_fastas_are_checked_without_fetch(source_fixture):
    repo, root, _ = prepare(source_fixture)
    def no_network(*args, **kwargs):
        raise AssertionError("no network should be accessed for verified files")
    report = fetcher.fetch_locked_fastas(repo, root, opener=no_network)
    assert report["verified_hosts"] == 3
    assert report["reused_hosts"] == 3
    assert report["downloaded_hosts"] == 0
    assert report["new_nvidia_requests"] == 0
    assert report["total_locked_bytes"] == sum(
        item["size_bytes"] for item in read_json_object(repo / extraction.GENOMES)["resources"]
    )


def test_fetch_missing_exact_resource_and_reuse_without_network(source_fixture):
    repo, root, _ = prepare(source_fixture)
    payload = (root / "host/A.fasta").read_bytes()
    (root / "host/A.fasta").unlink()
    urls = []
    def response(request, timeout):
        urls.append(request.full_url)
        return io.BytesIO(payload)
    report = fetcher.fetch_locked_fastas(repo, root, opener=response)
    assert len(urls) == 1
    assert report["verified_hosts"] == 3
    assert report["downloaded_hosts"] == 1
    assert report["reused_hosts"] == 2
    assert (root / "host/A.fasta").read_bytes() == payload
    assert not (root / "host/A.fasta.download-partial").exists()


def test_retries_transient_http_error_only(source_fixture):
    repo, root, _ = prepare(source_fixture)
    payload = (root / "host/A.fasta").read_bytes()
    (root / "host/A.fasta").unlink()
    attempts, backoffs = [], []
    def response(request, timeout):
        attempts.append(1)
        if len(attempts) == 1:
            raise urllib.error.HTTPError(request.full_url, 503, "Service Unavailable", {}, None)
        return io.BytesIO(payload)
    report = fetcher.fetch_locked_fastas(
        repo, root, opener=response, sleep=backoffs.append, attempts=2
    )
    assert len(attempts) == 2
    assert backoffs == [2]
    assert report["downloaded_hosts"] == 1


def test_bad_source_digest_rejected_without_admission(source_fixture):
    repo, root, _ = prepare(source_fixture)
    original = (root / "host/A.fasta").read_bytes()
    (root / "host/A.fasta").unlink()
    changed = b"X" + original[1:]
    with pytest.raises(ValueError, match="digest/size mismatch"):
        fetcher.fetch_locked_fastas(repo, root, opener=lambda *a, **kw: io.BytesIO(changed))
    assert not (root / "host/A.fasta").exists()
    assert not (root / "host/A.fasta.download-partial").exists()


def test_existing_corrupt_source_fails_without_overwrite(source_fixture):
    repo, root, _ = prepare(source_fixture)
    (root / "host/B.fasta").write_bytes(b"corrupt")
    def no_network(*args, **kwargs):
        raise AssertionError("corrupt existing source must never be overwritten")
    with pytest.raises(ValueError, match="existing host resource differs"):
        fetcher.fetch_locked_fastas(repo, root, opener=no_network)


def test_duplicate_or_untrusted_resource_is_refused(source_fixture):
    repo, root, _ = prepare(source_fixture)
    lock = read_json_object(repo / extraction.GENOMES)
    lock["resources"][0]["url"] = "http://insecure.test/file"
    atomic_write_json(lock, repo / extraction.GENOMES)
    with pytest.raises(ValueError, match="untrusted host"):
        fetcher.fetch_locked_fastas(repo, root)
    lock["resources"][0]["url"] = "https://ndownloader.figshare.com/files/123"
    lock["resources"].append(dict(lock["resources"][0]))
    atomic_write_json(lock, repo / extraction.GENOMES)
    with pytest.raises(ValueError, match="duplicate"):
        fetcher.fetch_locked_fastas(repo, root)
