"""Plan, run or verify bounded hosted extraction using local/S3 checkpoints."""
import argparse
import getpass
import os
from pathlib import Path

import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, read_json_object
from phage_cocktail.extraction_store import open_store
from phage_cocktail.nvidia_host_extraction import GENOMES, make_plan, run_extraction, verify_extraction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "run", "verify"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--hosts", help="Comma-separated frozen host IDs; default first three")
    group.add_argument("--all-hosts", action="store_true")
    parser.add_argument("--selection", choices=("first_full_window", "original_32_including_short_tails"),
                        default="first_full_window")
    parser.add_argument("--max-new-calls", type=int, default=0)
    parser.add_argument("--max-new-response-bytes", type=int, default=0)
    parser.add_argument("--store", help="Local directory or s3://bucket/dedicated-run-prefix")
    parser.add_argument("--s3-endpoint-url", help="Optional HTTPS endpoint of a compatible object store")
    parser.add_argument("--pilot-cache", type=Path)
    parser.add_argument("--allow-live", action="store_true", help="Permit requests within the frozen plan budget")
    parser.add_argument("--max-jobs", type=int, help="Pause after this many newly saved jobs")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    if args.command == "plan":
        if args.plan.exists():
            parser.error("plan output already exists; use a new versioned path")
        if args.store or args.allow_live or args.pilot_cache:
            parser.error("plan only writes a local manifest; it cannot execute jobs")
        hosts = sorted(args.hosts.split(",")) if args.hosts else None
        if args.all_hosts:
            source = read_json_object(repo / GENOMES)
            hosts = sorted(pd.read_csv(repo / source["host_manifest"], dtype=str).bacterium)
        plan = make_plan(repo, args.data_root, hosts=hosts, selection=args.selection,
                         max_new_calls=args.max_new_calls, max_new_response_bytes=args.max_new_response_bytes)
        atomic_write_json(plan, args.plan)
        print({"status": "planned", "hosts": plan["hosts"], "jobs": len(plan["jobs"]), "budget": plan["budget"]})
        return
    if not args.store:
        parser.error("run and verify require --store")
    plan = read_json_object(args.plan)
    store = open_store(args.store, endpoint_url=args.s3_endpoint_url)
    if args.command == "verify":
        print(verify_extraction(plan, store))
        return
    def key_provider():
        return os.environ.get("NVIDIA_API_KEY") or getpass.getpass("NVIDIA API key (hidden): ")
    print(run_extraction(repo, args.data_root, plan, store, pilot_cache=args.pilot_cache,
                         allow_live=args.allow_live, key_provider=key_provider, max_jobs=args.max_jobs))


if __name__ == "__main__":
    main()
