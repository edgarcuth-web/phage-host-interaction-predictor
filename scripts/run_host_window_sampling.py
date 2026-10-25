"""Compare 1/4/32-window host 4-mer composition against whole-genome controls."""
import argparse
from pathlib import Path

from phage_cocktail.retrieval_study import run_retrieval_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    run_retrieval_study(Path(__file__).resolve().parents[1], args.data_root,
                        args.outdir, fetch=args.fetch, sampling=True)


if __name__ == "__main__":
    main()
