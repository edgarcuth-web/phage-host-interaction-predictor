"""Check published BASEL/ECOR IDs against frozen training resource identities."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from urllib.request import Request, urlopen

from phage_cocktail.basel_ecor import inspect_basel_ecor, load_lock
from phage_cocktail.basel_identity import audit_basel_identity, write_report_new

ROOT = Path(__file__).resolve().parents[1]
SOURCE_LOCK = ROOT / "data/studies/basel_ecor_intake.v1.json"
RESOURCE_LOCK = ROOT / "data/studies/genophi_host_baseline.resources.v1.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    input_mode = parser.add_mutually_exclusive_group(required=True)
    input_mode.add_argument("--source", type=Path, help="untouched, locally stored pinned BASEL CSV")
    input_mode.add_argument("--fetch", action="store_true", help="fetch pinned 18,625-byte upstream source CSV")
    parser.add_argument("--resource-lock", type=Path, default=RESOURCE_LOCK)
    parser.add_argument("--out", type=Path, required=True, help="new non-overwriting JSON receipt")
    args = parser.parse_args()

    source_lock = load_lock(SOURCE_LOCK)
    if args.fetch:
        req = Request(source_lock["source"]["url"], headers={"User-Agent": "phage-source-intake/1"})
        with urlopen(req, timeout=30) as response:
            payload = response.read(source_lock["source"]["byte_size"] + 1)
    else:
        payload = args.source.read_bytes()
    source_report = inspect_basel_ecor(payload, source_lock)
    resource_bytes = args.resource_lock.read_bytes()
    resource_lock = json.loads(resource_bytes)
    report = audit_basel_identity(source_report, resource_lock)
    report["development_resource_lock_sha256"] = sha256(resource_bytes).hexdigest()
    report["source_input_sha256"] = source_report["local_input_sha256"]
    write_report_new(report, args.out)
    print("Host-ID overlap:", report["host_name_overlap_count"], "/", source_report["source_summary"]["hosts"])
    print("Phage accession map:", report["phage_genome_accessions_resolved"], "/", source_report["source_summary"]["phages"])
    print("Status:", report["status"])
    return 2  # an identity gate that cannot yet admit cold-phage comparisons


if __name__ == "__main__":
    raise SystemExit(main())
