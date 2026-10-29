"""Fetch only the SHA-256-locked host FASTAs for zero-request extraction preflight.

No NVIDIA requests, genome prediction, or model fitting are performed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import time
import urllib.error
import urllib.request

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file


LOCK = "data/studies/genome_host_retrieval.v1.json"
RETRYABLE = {429, 500, 502, 503, 504}


def _download_one(root, entry, *, opener, sleep, attempts):
    relative = Path(entry["path"])
    if (
        relative.is_absolute() or ".." in relative.parts or not relative.parts
        or relative.parts[0] != "host"
        or not entry["url"].startswith("https://ndownloader.figshare.com/files/")
    ):
        raise ValueError("untrusted host resource path or URL")
    target = (root / relative).resolve()
    if not target.is_relative_to(root):
        raise ValueError("host resource escapes data root")
    if target.exists():
        if (not target.is_file() or target.stat().st_size != entry["size_bytes"]
                or sha256_file(target) != entry["sha256"]):
            raise ValueError(f"existing host resource differs: {relative}")
        return {"id": entry["id"], "path": relative.as_posix(), "status": "verified_existing"}

    target.parent.mkdir(parents=True, exist_ok=True)
    temp = target.with_name(target.name + ".download-partial")
    if temp.exists():
        raise ValueError(f"stale download partial requires review: {temp.name}")
    for attempt in range(attempts):
        digest = hashlib.sha256()
        downloaded = 0
        try:
            request = urllib.request.Request(
                entry["url"], headers={"User-Agent": "phage-host-locked-preflight/1"}
            )
            with opener(request, timeout=90) as response, temp.open("xb") as out:
                while chunk := response.read(1024 * 1024):
                    downloaded += len(chunk)
                    if downloaded > entry["size_bytes"]:
                        raise ValueError(f"host download exceeds locked size: {entry['id']}")
                    digest.update(chunk)
                    out.write(chunk)
            if downloaded != entry["size_bytes"] or digest.hexdigest() != entry["sha256"]:
                raise ValueError(f"host download digest/size mismatch: {entry['id']}")
            temp.replace(target)
            return {"id": entry["id"], "path": relative.as_posix(), "status": "downloaded_verified"}
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            status = exc.code if isinstance(exc, urllib.error.HTTPError) else None
            is_retryable = (status in RETRYABLE if status is not None else
                            isinstance(exc, (urllib.error.URLError, TimeoutError)))
            if not is_retryable or attempt == attempts - 1:
                raise ValueError(f"cannot fetch locked host {entry['id']}: {type(exc).__name__}") from None
            sleep(min(2 ** (attempt + 1), 16))
        finally:
            temp.unlink(missing_ok=True)
    raise AssertionError("unreachable")


def fetch_locked_fastas(repo, data_root, *, opener=urllib.request.urlopen, sleep=time.sleep,
                        attempts=5, workers=4):
    repo, root = Path(repo), Path(data_root).resolve()
    if attempts < 1 or workers < 1:
        raise ValueError("attempts and workers must be positive")
    lock_path = repo / LOCK
    lock = read_json_object(lock_path)
    entries = [r for r in lock["resources"] if r["kind"] == "host_genome"]
    paths = [r["path"] for r in entries]
    ids = [r["id"] for r in entries]
    if not entries or len(set(paths)) != len(paths) or len(set(ids)) != len(ids):
        raise ValueError("host resource lock has missing or duplicate identities")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(
            lambda item: _download_one(root, item, opener=opener, sleep=sleep, attempts=attempts),
            entries
        ))
    return {
        "status": "source_fastas_verified",
        "lock_sha256": sha256_file(lock_path),
        "verified_hosts": len(results),
        "total_locked_bytes": sum(r["size_bytes"] for r in entries),
        "downloaded_hosts": sum(r["status"] == "downloaded_verified" for r in results),
        "reused_hosts": sum(r["status"] == "verified_existing" for r in results),
        "new_nvidia_requests": 0,
        "resources": sorted(results, key=lambda r: r["id"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("output exists; do not overwrite immutable receipt")
    repo = Path(__file__).resolve().parents[1]
    receipt = fetch_locked_fastas(repo, args.data_root)
    atomic_write_json(receipt, args.out)
    print({key: receipt[key] for key in (
        "status", "verified_hosts", "total_locked_bytes", "downloaded_hosts",
        "reused_hosts", "new_nvidia_requests",
    )})


if __name__ == "__main__":
    main()
