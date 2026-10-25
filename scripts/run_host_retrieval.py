"""Evaluate unseen-host retrieval on a fixed known-phage library."""
import argparse
from pathlib import Path

from phage_cocktail.retrieval_study import run_retrieval_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True, help="A new run directory")
    parser.add_argument("--fetch", action="store_true", help="Download four locked source tables")
    args = parser.parse_args()
    report = run_retrieval_study(Path(__file__).resolve().parents[1], args.data_root,
                                 args.outdir, fetch=args.fetch)
    for row in report["metrics"]:
        if row["scenario"] == "grouped_host" and row["partition"] == "test":
            print(row["model"], row["all_hosts"], flush=True)


if __name__ == "__main__":
    main()
