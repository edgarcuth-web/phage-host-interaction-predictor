"""Compare a complete hosted Evo2 host extraction under frozen retrieval splits."""
import argparse
from pathlib import Path

from phage_cocktail.extraction_store import open_store
from phage_cocktail.retrieval_study import run_retrieval_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--store", required=True, help="Local directory or s3://bucket/dedicated-run-prefix")
    parser.add_argument("--s3-endpoint-url", help="Optional HTTPS endpoint of a compatible object store")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Fetch pinned tables only; never model weights")
    args = parser.parse_args()

    store = open_store(args.store, endpoint_url=args.s3_endpoint_url)
    run_retrieval_study(
        Path(__file__).resolve().parents[1],
        args.data_root,
        args.outdir,
        fetch=args.fetch,
        hosted_embedding=(args.plan, store),
    )


if __name__ == "__main__":
    main()
