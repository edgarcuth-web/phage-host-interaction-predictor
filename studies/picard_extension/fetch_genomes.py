#!/usr/bin/env python3
"""Retrieve the 225 additional assemblies from the frozen Figshare roster."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def key(name):
    return name.upper().replace("-", "")


def unique_lookup(records, field):
    lookup = {}
    for r in records:
        k = key(r[field])
        if k in lookup:
            raise ValueError(f"ambiguous source alias: {r[field]}")
        lookup[k] = r
    return lookup


def retrieve(item, out):
    target = out / f"{item['id']}.fasta"
    if target.exists():
        if target.stat().st_size == item["size"] and hashlib.md5(target.read_bytes()).hexdigest() == item["computed_md5"]:
            return item["name"]
        raise ValueError(f"changed cached genome: {item['name']}")
    request = urllib.request.Request(item["download_url"], headers={"User-Agent": "picard-extension-v1"})
    partial = target.with_suffix(".partial")
    for attempt in range(3):
        try:
            md5, size = hashlib.md5(), 0
            with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as handle:
                while chunk := response.read(1024 * 1024):
                    size += len(chunk)
                    if size > item["size"]:
                        raise ValueError(f"oversized genome: {item['name']}")
                    md5.update(chunk)
                    handle.write(chunk)
            if size != item["size"] or md5.hexdigest() != item["computed_md5"]:
                raise ValueError(f"genome size or MD5 mismatch: {item['name']}")
            os.replace(partial, target)
            return item["name"]
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == 2:
                raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, default=ROOT / "data/raw/picard_extension_sources")
    p.add_argument("--out", type=Path, default=ROOT / "data/raw/picard_extension_genomes")
    p.add_argument("--workers", type=int, default=6)
    args = p.parse_args()
    if not 1 <= args.workers <= 8:
        p.error("workers must be between one and eight")
    with (args.sources / "source_mlc.csv").open() as f:
        matrix = list(csv.DictReader(f, delimiter=";"))
    with (ROOT / "data/studies/results/genophi_comparator.v1/scenario_splits.csv").open() as f:
        old = {r["source_bacterium"] for r in csv.DictReader(f)}
    roster = sorted(r["bacteria"] for r in matrix if r["bacteria"] not in old)
    files = json.loads((ROOT / "studies/picard_extension/genome_sources.json").read_text())["files"]
    lookup = unique_lookup([{**r, "host": r["name"].removesuffix(".fasta")} for r in files], "host")
    if len(roster) != 225:
        raise ValueError("additional source roster differs")
    items = [lookup[key(h)] for h in roster]
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"Fetching {len(items)} assemblies; {sum(x['size'] for x in items):,} bytes", flush=True)
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(retrieve, item, args.out): item["name"] for item in items}
        for completed, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                future.result()
            except Exception as exc:
                failures.append({"file": futures[future], "error": str(exc)})
            if completed % 10 == 0 or completed == len(items):
                print(f"{completed}/{len(items)} checked; failures={len(failures)}", flush=True)
    if failures:
        raise RuntimeError(json.dumps(failures))
    print("All 225 additional source assemblies verified", flush=True)


if __name__ == "__main__":
    main()
