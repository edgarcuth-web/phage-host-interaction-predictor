"""Verify a pinned BASEL/ECOR source interaction matrix without model inference."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.request import Request, urlopen

from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock

REPO = Path(__file__).resolve().parents[1]
LOCK = REPO / "data/studies/basel_ecor_intake.v1.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--source", type=Path, help="already downloaded unmodified upstream CSV")
    inputs.add_argument("--fetch", action="store_true", help="download pinned GitHub raw CSV (18,625 bytes)")
    parser.add_argument("--out", type=Path, required=True, help="new JSON summary; no overwrite")
    args = parser.parse_args()
    lock = load_lock(LOCK)
    if args.fetch:
        request = Request(lock["source"]["url"], headers={"User-Agent": "phage-source-intake/1"})
        with urlopen(request, timeout=30) as response:
            payload = response.read(lock["source"]["byte_size"] + 1)
    else:
        payload = args.source.read_bytes()
    report = inspect_basel_ecor(payload, lock)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(report, output, sort_keys=True, indent=2)
        output.write("\n")
    print("Verified source-only intake: ", report["source_summary"])
    print("Not an independent release validation: ", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
